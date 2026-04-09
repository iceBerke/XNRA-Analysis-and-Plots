"""
XNRA Uncertainty Calculator v2
==============================
Calculates statistical uncertainties for IBA data from SIMNRA .xnra files.

CHANGELOG v2:
- FIXED: Consistent use of fitted range (handles partial spectrum fitting)
- ADDED: Automatic detection of fitted channel range from simulated spectrum
- ADDED: fit_range_start/fit_range_end in output for transparency
- IMPROVED: Chi-squared and counting statistics now use same fitted region

Authors: Based on work by Berke Santos & Giuseppe Legrottaglie
Extended for uncertainty analysis
"""

import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re


# =============================================================================
# PARSING FUNCTIONS
# =============================================================================

def parse_xnra_full(filepath: str) -> Dict:
    """
    Parse .xnra file and extract all relevant data including layers.
    
    Returns dict with:
        - metadata (filename, beam, geometry)
        - layers (thickness, element concentrations)
        - spectra (raw counts, simulated counts)
        - calibration
        - fit_range (detected from simulation extent)
    """
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    ns = {
        'idf': 'http://idf.schemas.itn.pt',
        'simnra': 'http://www.simnra.com/simnra'
    }
    
    data = {
        'filepath': str(filepath),
        'filename': Path(filepath).stem,
        'layers': [],
        'elements': [],
        'raw_channels': None,
        'raw_counts': None,
        'simulated_channels': None,
        'simulated_counts': None,
        'beam_energy': None,
        'beam_particle': None,
        'scattering_angle': None,
        'cal_offset': 0,
        'cal_gain': 1,
        # NEW: Fit range info
        'fit_range_start': None,
        'fit_range_end': None,
    }
    
    # ===== Extract Metadata =====
    beam_particle = root.find('.//idf:beamparticle', ns)
    if beam_particle is not None:
        data['beam_particle'] = beam_particle.text
    
    beam_energy = root.find('.//idf:beamenergy', ns)
    if beam_energy is not None:
        data['beam_energy'] = float(beam_energy.text)
    
    scattering_angle = root.find('.//idf:scatteringangle', ns)
    if scattering_angle is not None:
        data['scattering_angle'] = float(scattering_angle.text)
    
    # ===== Extract Calibration =====
    cal_params = root.findall('.//idf:energycalibration[1]/idf:calibrationparameters/idf:calibrationparameter', ns)
    if len(cal_params) >= 2:
        data['cal_offset'] = float(cal_params[0].text)
        data['cal_gain'] = float(cal_params[1].text)
    
    # ===== Extract Elements List =====
    elements = root.findall('.//idf:elementsandmolecules/idf:elements/idf:element/idf:n', ns)
    data['elements'] = [el.text for el in elements]
    
    # ===== Extract Layer Data =====
    layers = root.findall('.//idf:layeredstructure/idf:layers/idf:layer', ns)
    
    for i, layer in enumerate(layers):
        layer_data = {
            'layer_num': i + 1,
            'thickness': None,
            'concentrations': {}
        }
        
        # Thickness
        thickness_el = layer.find('idf:layerthickness', ns)
        if thickness_el is not None:
            layer_data['thickness'] = float(thickness_el.text)
        
        # Element concentrations
        layer_elements = layer.findall('.//idf:layerelement', ns)
        for el in layer_elements:
            # Try both 'name' and 'n' tags (different SIMNRA versions)
            name_el = el.find('idf:name', ns)
            if name_el is None:
                name_el = el.find('idf:n', ns)
            conc_el = el.find('idf:concentration', ns)
            
            if name_el is not None and conc_el is not None:
                element_name = name_el.text
                concentration = float(conc_el.text)
                
                # Normalize isotope names (e.g., "40Ca" -> "Ca")
                base_element = re.sub(r'^\d+', '', element_name)
                
                # Sum isotopes for the same element
                if base_element in layer_data['concentrations']:
                    layer_data['concentrations'][base_element] += concentration
                else:
                    layer_data['concentrations'][base_element] = concentration
        
        data['layers'].append(layer_data)
    
    # ===== Extract Raw Spectrum Data =====
    raw_data = root.find('.//idf:data/idf:simpledata', ns)
    if raw_data is not None:
        x_el = raw_data.find('idf:x', ns)
        y_el = raw_data.find('idf:y', ns)
        if x_el is not None and y_el is not None:
            data['raw_channels'] = np.array([float(x) for x in x_el.text.split()])
            data['raw_counts'] = np.array([float(y) for y in y_el.text.split()])
    
    # ===== Extract Simulated Spectrum Data =====
    sim_data = root.find('.//idf:simulation[1]/idf:simpledata', ns)
    if sim_data is None:
        sim_data = root.find('.//simnra:calculateddata/idf:simpledata', ns)
    
    if sim_data is not None:
        x_el = sim_data.find('idf:x', ns)
        y_el = sim_data.find('idf:y', ns)
        if x_el is not None and y_el is not None:
            data['simulated_channels'] = np.array([float(x) for x in x_el.text.split()])
            data['simulated_counts'] = np.array([float(y) for y in y_el.text.split()])
    
    # ===== Detect Fit Range =====
    data['fit_range_start'], data['fit_range_end'] = detect_fit_range(data)
    
    return data


def detect_fit_range(data: Dict) -> Tuple[int, int]:
    """
    Detect the channel range that was actually fitted.
    
    Strategy:
    1. Use simulated spectrum extent (most reliable)
    2. Look for non-zero regions if simulated channels are different from raw
    3. Fall back to full raw spectrum range
    
    Returns:
        (start_channel, end_channel) - inclusive range
    """
    raw_counts = data.get('raw_counts')
    sim_counts = data.get('simulated_counts')
    sim_channels = data.get('simulated_channels')
    raw_channels = data.get('raw_channels')
    
    # Default: use full raw spectrum
    if raw_counts is None:
        return 0, 0
    
    start = 0
    end = len(raw_counts)
    
    if sim_counts is not None:
        # Method 1: Simulated spectrum is shorter than raw
        # This is the most common case when fitting partial spectrum
        if len(sim_counts) < len(raw_counts):
            end = len(sim_counts)
        
        # Method 2: Check if simulated channels indicate a subset
        elif sim_channels is not None and raw_channels is not None:
            # Simulated might start/end at different channels
            if len(sim_channels) > 0:
                sim_start = int(sim_channels[0])
                sim_end = int(sim_channels[-1]) + 1
                
                # Only use if it's a true subset
                if sim_start > 0 or sim_end < len(raw_counts):
                    start = sim_start
                    end = sim_end
        
        # Method 3: Look for non-zero simulation region
        # (handles case where sim array is same length but padded with zeros)
        if len(sim_counts) == len(raw_counts):
            non_zero = np.where(sim_counts > 0)[0]
            if len(non_zero) > 0:
                # Use the extent of non-zero simulated data
                # Add small buffer to avoid edge effects
                start = max(0, non_zero[0])
                end = min(len(raw_counts), non_zero[-1] + 1)
    
    return start, end


# =============================================================================
# UNCERTAINTY CALCULATIONS
# =============================================================================

def calculate_chi_squared_reduced(raw_counts: np.ndarray, 
                                   simulated_counts: np.ndarray,
                                   fit_start: int = 0,
                                   fit_end: int = None,
                                   min_counts: float = 4.0) -> Tuple[float, int]:
    """
    Calculate reduced chi-squared for spectrum fit.
    
    Uses the FITTED REGION only (as specified by fit_start/fit_end).
    
    Parameters:
        raw_counts: Experimental counts array
        simulated_counts: Simulated/fitted counts array
        fit_start: Start channel of fitted region
        fit_end: End channel of fitted region (None = auto from simulation length)
        min_counts: Minimum counts for uncertainty (avoid div by zero)
    
    Returns:
        (chi2_reduced, degrees_of_freedom)
    """
    if fit_end is None:
        fit_end = min(len(raw_counts), len(simulated_counts))
    
    # Ensure we don't exceed array bounds
    fit_end = min(fit_end, len(raw_counts), len(simulated_counts))
    fit_start = max(0, fit_start)
    
    if fit_start >= fit_end:
        return np.nan, 0
    
    # Extract fitted region
    raw_roi = raw_counts[fit_start:fit_end]
    sim_roi = simulated_counts[fit_start:fit_end]
    
    if len(raw_roi) == 0:
        return np.nan, 0
    
    # Calculate uncertainty (Poisson: sqrt(N), minimum of min_counts)
    sigma = np.sqrt(np.maximum(raw_roi, min_counts))
    
    # Chi-squared
    residuals = (raw_roi - sim_roi) / sigma
    chi2 = np.sum(residuals**2)
    
    # Degrees of freedom (data points - 1)
    dof = len(raw_roi) - 1
    
    chi2_reduced = chi2 / dof if dof > 0 else np.inf
    
    return chi2_reduced, dof


def calculate_counting_statistics_error(raw_counts: np.ndarray,
                                         concentration: float,
                                         fit_start: int = 0,
                                         fit_end: int = None) -> Tuple[float, int]:
    """
    Calculate statistical uncertainty on concentration from counting statistics.
    
    The relative uncertainty on concentration equals the relative uncertainty
    on the total counts in the fitted region.
    
    σ_C / C = σ_N / N = √N / N = 1 / √N
    
    Parameters:
        raw_counts: Experimental counts array
        concentration: Fitted concentration (fraction, not %)
        fit_start: Start channel of fitted region
        fit_end: End channel of fitted region
    
    Returns:
        (sigma_concentration, total_counts_in_region)
    """
    if fit_end is None:
        fit_end = len(raw_counts)
    
    # Ensure bounds
    fit_start = max(0, fit_start)
    fit_end = min(fit_end, len(raw_counts))
    
    # Total counts in FITTED region only
    raw_roi = raw_counts[fit_start:fit_end]
    N_total = np.sum(raw_roi)
    
    if N_total <= 0:
        return np.nan, 0
    
    # Relative uncertainty
    relative_error = 1.0 / np.sqrt(N_total)
    
    # Absolute uncertainty on concentration
    sigma_concentration = concentration * relative_error
    
    return sigma_concentration, int(N_total)


def estimate_systematic_uncertainty(technique: str) -> float:
    """
    Return typical systematic uncertainty for IBA techniques.
    
    These are literature values for stopping power and cross-section uncertainties.
    
    Parameters:
        technique: 'RBS', 'ERDA', or 'NRA'
    
    Returns:
        Relative systematic uncertainty (e.g., 0.03 = 3%)
    """
    systematic = {
        'RBS': 0.03,   # 3% - stopping powers
        'ERDA': 0.05,  # 5% - stopping powers + cross-sections
        'NRA': 0.05,   # 5% - reaction cross-sections
    }
    return systematic.get(technique.upper(), 0.05)


def combine_uncertainties(sigma_stat: float, sigma_sys: float) -> float:
    """
    Combine statistical and systematic uncertainties in quadrature.
    
    σ_total = √(σ_stat² + σ_sys²)
    """
    return np.sqrt(sigma_stat**2 + sigma_sys**2)


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def detect_technique(filepath: str) -> str:
    """
    Detect IBA technique from filename or file content.
    """
    filename = Path(filepath).stem.upper()
    
    if 'ERDA' in filename:
        return 'ERDA'
    elif 'NRA' in filename:
        return 'NRA'
    elif 'RBS' in filename:
        return 'RBS'
    else:
        # Try to detect from beam particle
        data = parse_xnra_full(filepath)
        if data['beam_particle']:
            particle = data['beam_particle'].upper()
            if 'H' in particle or 'D' in particle:
                return 'NRA'  # Likely NRA with protons/deuterons
        return 'RBS'  # Default


def process_single_file(filepath: str, 
                        elements_of_interest: List[str] = None,
                        technique: str = None,
                        force_fit_range: Tuple[int, int] = None) -> List[Dict]:
    """
    Process a single .xnra file and calculate uncertainties.
    
    Parameters:
        filepath: Path to .xnra file
        elements_of_interest: List of elements to analyze (None = all)
        technique: 'RBS', 'ERDA', 'NRA' (None = auto-detect)
        force_fit_range: Override auto-detected fit range with (start, end)
    
    Returns:
        List of dicts with results per element per layer
    """
    data = parse_xnra_full(filepath)
    
    if technique is None:
        technique = detect_technique(filepath)
    
    # Get fit range (auto-detected or forced)
    if force_fit_range is not None:
        fit_start, fit_end = force_fit_range
    else:
        fit_start = data['fit_range_start']
        fit_end = data['fit_range_end']
    
    results = []
    
    # Calculate chi-squared using FITTED REGION ONLY
    if data['raw_counts'] is not None and data['simulated_counts'] is not None:
        chi2_r, dof = calculate_chi_squared_reduced(
            data['raw_counts'], 
            data['simulated_counts'],
            fit_start=fit_start,
            fit_end=fit_end
        )
    else:
        chi2_r, dof = np.nan, 0
    
    # Total counts in FITTED REGION ONLY (not full spectrum!)
    if data['raw_counts'] is not None:
        fitted_counts = np.sum(data['raw_counts'][fit_start:fit_end])
        total_raw_counts = np.sum(data['raw_counts'])  # For reference
    else:
        fitted_counts = 0
        total_raw_counts = 0
    
    # Process each layer
    for layer in data['layers']:
        for element, concentration in layer['concentrations'].items():
            
            # Skip if not in elements of interest
            if elements_of_interest and element not in elements_of_interest:
                continue
            
            # Skip zero concentrations
            if concentration <= 0:
                continue
            
            # Calculate statistical uncertainty using FITTED REGION
            if fitted_counts > 0:
                # Element contributes proportionally to its concentration
                # in the fitted region
                effective_counts = fitted_counts * concentration
                sigma_stat_rel = 1.0 / np.sqrt(max(effective_counts, 1))
                sigma_stat = concentration * sigma_stat_rel
            else:
                sigma_stat = np.nan
            
            # Systematic uncertainty
            sigma_sys_rel = estimate_systematic_uncertainty(technique)
            sigma_sys = concentration * sigma_sys_rel
            
            # Combined uncertainty
            if not np.isnan(sigma_stat):
                sigma_total = combine_uncertainties(sigma_stat, sigma_sys)
            else:
                sigma_total = sigma_sys
            
            results.append({
                'filename': data['filename'],
                'filepath': data['filepath'],
                'technique': technique,
                'beam_energy_keV': data['beam_energy'],
                'layer': layer['layer_num'],
                'thickness_1e15_at_cm2': layer['thickness'],
                'element': element,
                'concentration_fraction': concentration,
                'concentration_percent': concentration * 100,
                'sigma_stat_fraction': sigma_stat,
                'sigma_stat_percent': sigma_stat * 100 if not np.isnan(sigma_stat) else np.nan,
                'sigma_sys_fraction': sigma_sys,
                'sigma_sys_percent': sigma_sys * 100,
                'sigma_total_fraction': sigma_total,
                'sigma_total_percent': sigma_total * 100,
                'chi2_reduced': chi2_r,
                # NEW: Fit range info for transparency
                'fit_range_start': fit_start,
                'fit_range_end': fit_end,
                'fitted_counts': fitted_counts,
                'total_raw_counts': total_raw_counts,
                'fit_fraction': fitted_counts / total_raw_counts if total_raw_counts > 0 else np.nan,
            })
    
    return results


def process_batch(filepaths: List[str],
                  elements_of_interest: List[str] = None,
                  output_csv: str = None,
                  output_excel: str = None,
                  force_fit_range: Tuple[int, int] = None) -> pd.DataFrame:
    """
    Process multiple .xnra files and compile results.
    
    Parameters:
        filepaths: List of paths to .xnra files
        elements_of_interest: Elements to analyze (None = all)
        output_csv: Path for CSV output (optional)
        output_excel: Path for Excel output (optional)
        force_fit_range: Override auto-detected fit range for ALL files
    
    Returns:
        DataFrame with all results
    """
    all_results = []
    
    for filepath in filepaths:
        try:
            results = process_single_file(
                filepath, 
                elements_of_interest,
                force_fit_range=force_fit_range
            )
            all_results.extend(results)
            
            # Show fit range info
            if results:
                r = results[0]
                print(f"✓ {Path(filepath).name}: channels {r['fit_range_start']}-{r['fit_range_end']} "
                      f"({r['fit_fraction']*100:.1f}% of spectrum)")
        except Exception as e:
            print(f"✗ Error processing {filepath}: {e}")
    
    df = pd.DataFrame(all_results)
    
    # Sort by filename, layer, element
    if not df.empty:
        df = df.sort_values(['filename', 'layer', 'element'])
    
    # Save outputs
    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"\n✓ Saved CSV: {output_csv}")
    
    if output_excel:
        df.to_excel(output_excel, index=False, sheet_name='Uncertainties')
        print(f"✓ Saved Excel: {output_excel}")
    
    return df


def find_xnra_files(directory: str, pattern: str = "*.xnra") -> List[str]:
    """
    Find all .xnra files in a directory.
    """
    path = Path(directory)
    return sorted([str(f) for f in path.glob(pattern)])


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

def create_depth_profile_table(df: pd.DataFrame,
                                elements: List[str] = None,
                                sample_name: str = None) -> pd.DataFrame:
    """
    Create a depth profile table for a single sample.
    
    Rows = layers (depth)
    Columns = element concentrations and uncertainties
    
    This is the format you want for ion depletion studies!
    """
    if elements is None:
        elements = sorted(df['element'].unique())
    
    # Filter to single sample if specified
    if sample_name:
        df = df[df['filename'] == sample_name]
    
    # Get unique layers sorted
    layers = sorted(df['layer'].unique())
    
    rows = []
    for layer_num in layers:
        layer_df = df[df['layer'] == layer_num]
        
        # Get thickness (should be same for all elements in layer)
        thickness = layer_df['thickness_1e15_at_cm2'].iloc[0] if len(layer_df) > 0 else None
        
        row = {
            'Layer': layer_num,
            'Thickness (10¹⁵ at/cm²)': thickness,
        }
        
        for element in elements:
            el_row = layer_df[layer_df['element'] == element]
            
            if len(el_row) > 0:
                conc = el_row['concentration_percent'].iloc[0]
                unc = el_row['sigma_total_percent'].iloc[0]
                row[f'{element} (%)'] = conc
                row[f'±{element}'] = unc
            else:
                row[f'{element} (%)'] = None
                row[f'±{element}'] = None
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def create_all_depth_profiles(df: pd.DataFrame,
                               elements: List[str] = None,
                               output_excel: str = None) -> Dict[str, pd.DataFrame]:
    """
    Create depth profile tables for ALL samples in the dataset.
    
    Returns dict of {sample_name: DataFrame}
    Optionally saves to Excel with one sheet per sample.
    """
    samples = df['filename'].unique()
    profiles = {}
    
    for sample in samples:
        sample_df = df[df['filename'] == sample]
        profile = create_depth_profile_table(sample_df, elements)
        profiles[sample] = profile
    
    if output_excel:
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            for sample, profile in profiles.items():
                # Sheet names max 31 chars
                sheet_name = sample[:31]
                profile.to_excel(writer, sheet_name=sheet_name, index=False)
        print(f"✓ Saved depth profiles: {output_excel}")
    
    return profiles


def create_summary_by_sample(df: pd.DataFrame, 
                              elements: List[str] = None) -> pd.DataFrame:
    """
    Create summary table with one row per sample, columns for each element.
    
    Shows average concentration ± uncertainty across layers.
    """
    if elements is None:
        elements = df['element'].unique()
    
    summary_data = []
    
    for filename in df['filename'].unique():
        sample_df = df[df['filename'] == filename]
        
        row = {'sample': filename}
        
        # Add fit range info
        if 'fit_range_start' in sample_df.columns:
            row['fit_range'] = f"{sample_df['fit_range_start'].iloc[0]}-{sample_df['fit_range_end'].iloc[0]}"
        
        for element in elements:
            el_df = sample_df[sample_df['element'] == element]
            
            if len(el_df) > 0:
                # Weighted average by layer thickness
                weights = el_df['thickness_1e15_at_cm2'].values
                concentrations = el_df['concentration_percent'].values
                uncertainties = el_df['sigma_total_percent'].values
                
                if weights.sum() > 0:
                    avg_conc = np.average(concentrations, weights=weights)
                    # Propagate uncertainty
                    avg_unc = np.sqrt(np.sum((weights * uncertainties)**2)) / weights.sum()
                else:
                    avg_conc = np.mean(concentrations)
                    avg_unc = np.mean(uncertainties)
                
                row[f'{element}_conc_%'] = avg_conc
                row[f'{element}_unc_%'] = avg_unc
            else:
                row[f'{element}_conc_%'] = np.nan
                row[f'{element}_unc_%'] = np.nan
        
        summary_data.append(row)
    
    return pd.DataFrame(summary_data)


# =============================================================================
# MAIN / CLI
# =============================================================================

def main():
    """
    Command-line interface for batch processing.
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Calculate uncertainties for SIMNRA .xnra files (v2 - consistent fit range)'
    )
    parser.add_argument('input', nargs='+', 
                        help='Input .xnra file(s) or directory')
    parser.add_argument('-e', '--elements', nargs='+',
                        help='Elements to analyze (default: all)')
    parser.add_argument('-o', '--output', 
                        help='Output file (CSV or Excel based on extension)')
    parser.add_argument('--summary', action='store_true',
                        help='Also create summary by sample')
    parser.add_argument('--fit-range', nargs=2, type=int, metavar=('START', 'END'),
                        help='Force fit range (channels), overrides auto-detection')
    
    args = parser.parse_args()
    
    # Collect input files
    filepaths = []
    for inp in args.input:
        p = Path(inp)
        if p.is_dir():
            filepaths.extend(find_xnra_files(inp))
        elif p.is_file() and p.suffix.lower() == '.xnra':
            filepaths.append(str(p))
    
    if not filepaths:
        print("No .xnra files found!")
        return
    
    print(f"Processing {len(filepaths)} files...")
    
    # Parse fit range override
    force_fit_range = tuple(args.fit_range) if args.fit_range else None
    if force_fit_range:
        print(f"Using forced fit range: channels {force_fit_range[0]}-{force_fit_range[1]}")
    
    # Determine output format
    output_csv = None
    output_excel = None
    if args.output:
        if args.output.endswith('.xlsx'):
            output_excel = args.output
        else:
            output_csv = args.output if args.output.endswith('.csv') else args.output + '.csv'
    
    # Process
    df = process_batch(filepaths, args.elements, output_csv, output_excel, force_fit_range)
    
    print(f"\nProcessed {len(df)} element-layer combinations")
    
    if args.summary and not df.empty:
        summary = create_summary_by_sample(df, args.elements)
        summary_path = args.output.replace('.csv', '_summary.csv').replace('.xlsx', '_summary.xlsx') if args.output else 'summary.csv'
        summary.to_csv(summary_path, index=False)
        print(f"✓ Saved summary: {summary_path}")


if __name__ == '__main__':
    main()
