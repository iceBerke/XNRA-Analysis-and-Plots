"""
XNRA Uncertainty Calculator v5
==============================
Calculates statistical uncertainties for IBA data from SIMNRA .xnra files.

This module implements ROI-based (Region of Interest) uncertainty analysis,
which is physically correct because different elements produce signals in
different channel ranges of the spectrum.

=============================================================================
THEORETICAL BACKGROUND
=============================================================================

1. COUNTING STATISTICS (Statistical Uncertainty)
-------------------------------------------------
Ion beam analysis is fundamentally a counting experiment. The number of 
detected particles follows Poisson statistics, where the standard deviation
equals the square root of the number of counts:

    σ_N = √N

For a relative uncertainty:

    σ_N / N = 1 / √N

For an element with concentration C in a region with total counts N_roi,
the element contributes approximately (C × N_roi) counts. Therefore:

    σ_stat(C) = C × (1 / √(C × N_roi)) = √(C / N_roi)

Or equivalently:

    σ_stat(C) / C = 1 / √(C × N_roi)

This means:
- Higher counts → lower statistical uncertainty
- Higher concentration → lower relative uncertainty (more signal)
- The ROI should be chosen where the element's signal is strongest


2. CHI-SQUARED (χ²) - FIT QUALITY METRIC
----------------------------------------
Chi-squared quantifies how well the simulation matches the experimental data:

    χ² = Σ [(O_i - E_i)² / σ_i²]

Where:
- O_i = Observed counts in channel i (raw experimental data)
- E_i = Expected counts in channel i (simulated/fitted data)
- σ_i = Uncertainty in channel i = √(max(O_i, min_counts))
  (We use min_counts=4 to avoid division issues at low statistics)

The reduced chi-squared normalizes by degrees of freedom:

    χ²/ν = χ² / (n - p)

Where:
- n = number of data points (channels in ROI)
- p = number of fitted parameters (we use p=1 as approximation)

INTERPRETATION:
- χ²/ν ≈ 1.0: Excellent fit, uncertainties are reliable
- χ²/ν < 1.0: Possibly overestimated uncertainties or overfitting
- χ²/ν = 2-5: Acceptable fit, uncertainties may be underestimated
- χ²/ν > 10: Poor fit in this region, uncertainties are questionable


3. SYSTEMATIC UNCERTAINTIES
---------------------------
These arise from fundamental limitations in the physics models:

RBS (Rutherford Backscattering Spectrometry): ~3%
- Dominated by stopping power uncertainties
- Well-established Coulomb scattering cross-sections

ERDA (Elastic Recoil Detection Analysis): ~5%
- Stopping power uncertainties for both beam and recoils
- Recoil cross-section uncertainties

NRA (Nuclear Reaction Analysis): ~5-10%
- Nuclear reaction cross-section uncertainties
- These can vary significantly by reaction

The systematic uncertainty on concentration is:

    σ_sys(C) = C × σ_sys_relative


4. TOTAL UNCERTAINTY (Combined)
-------------------------------
Statistical and systematic uncertainties are independent, so they combine
in quadrature:

    σ_total = √(σ_stat² + σ_sys²)


5. ROI-BASED ANALYSIS RATIONALE
-------------------------------
Different elements produce peaks at different energies (channels) in the
spectrum due to kinematics:

- Heavy elements (Ba, Ca, Fe) → High channels (high backscatter energy)
- Medium elements (Si, Al, Na) → Medium channels
- Light elements (O, C, H) → Low channels (for RBS) or detected via ERDA/NRA

Using one global chi² for all elements is physically incorrect because:
- The fit quality may vary across the spectrum
- An element's uncertainty should depend on counts where ITS signal appears

Therefore, we calculate separate chi² and counting statistics for each ROI.

=============================================================================
AUTHORS
=============================================================================
Based on work by Berke Santos & Giuseppe Legrottaglie
Extended for ROI-based uncertainty analysis

=============================================================================
"""

import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re


# =============================================================================
# DEFAULT SYSTEMATIC UNCERTAINTIES (can be overridden)
# =============================================================================

DEFAULT_SYSTEMATIC_UNCERTAINTIES = {
    'RBS': 0.03,    # 3% - stopping powers
    'ERDA': 0.05,   # 5% - stopping powers + recoil cross-sections
    'NRA': 0.05,    # 5% - nuclear reaction cross-sections (can be higher)
}


# =============================================================================
# PARSING FUNCTIONS
# =============================================================================

def parse_xnra_full(filepath: str) -> Dict:
    """
    Parse .xnra file and extract all relevant data.
    
    The .xnra format is XML-based, following the IDF (Ion beam Data Format)
    schema with SIMNRA-specific extensions.
    
    Parameters
    ----------
    filepath : str
        Path to the .xnra file
    
    Returns
    -------
    dict
        Dictionary containing:
        - filepath, filename: File identification
        - layers: List of layer dicts with thickness and concentrations
        - raw_channels, raw_counts: Experimental spectrum
        - simulated_channels, simulated_counts: Fitted/simulated spectrum
        - beam_energy, beam_particle, scattering_angle: Experimental setup
        - cal_offset, cal_gain: Energy calibration (E = offset + gain × channel)
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
        
        # Thickness (in 10^15 atoms/cm²)
        thickness_el = layer.find('idf:layerthickness', ns)
        if thickness_el is not None:
            layer_data['thickness'] = float(thickness_el.text)
        
        # Element concentrations
        layer_elements = layer.findall('.//idf:layerelement', ns)
        for el in layer_elements:
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
    
    return data


# =============================================================================
# UNCERTAINTY CALCULATIONS
# =============================================================================

def calculate_chi_squared_for_roi(raw_counts: np.ndarray, 
                                   simulated_counts: np.ndarray,
                                   roi_start: int,
                                   roi_end: int,
                                   min_counts: float = 4.0) -> Tuple[float, int, int]:
    """
    Calculate reduced chi-squared for a specific ROI (Region of Interest).
    
    The chi-squared statistic measures the goodness-of-fit between the
    experimental (raw) and simulated spectra within the specified channel range.
    
    Formula:
        χ² = Σ [(raw_i - sim_i)² / σ_i²]
        
    Where σ_i = √(max(raw_i, min_counts)) is the Poisson uncertainty.
    
    Parameters
    ----------
    raw_counts : np.ndarray
        Experimental spectrum counts
    simulated_counts : np.ndarray
        Simulated/fitted spectrum counts
    roi_start : int
        Start channel of ROI (inclusive)
    roi_end : int
        End channel of ROI (exclusive)
    min_counts : float, default=4.0
        Minimum counts for uncertainty calculation (avoids division by zero
        and stabilizes chi² at low statistics)
    
    Returns
    -------
    tuple
        (chi2_reduced, degrees_of_freedom, total_counts_in_roi)
        
    Notes
    -----
    - χ²/ν ≈ 1 indicates a good fit
    - χ²/ν >> 1 indicates poor fit or underestimated uncertainties
    - χ²/ν << 1 indicates overestimated uncertainties
    """
    # Ensure bounds are valid
    roi_start = max(0, roi_start)
    roi_end = min(roi_end, len(raw_counts), len(simulated_counts))
    
    if roi_start >= roi_end:
        return np.nan, 0, 0
    
    # Extract ROI
    raw_roi = raw_counts[roi_start:roi_end]
    sim_roi = simulated_counts[roi_start:roi_end]
    
    if len(raw_roi) == 0:
        return np.nan, 0, 0
    
    # Total counts in ROI (for counting statistics)
    total_counts = int(np.sum(raw_roi))
    
    # Poisson uncertainty: σ = √N, with minimum to avoid div/0
    sigma = np.sqrt(np.maximum(raw_roi, min_counts))
    
    # Chi-squared: Σ[(O-E)²/σ²]
    residuals = (raw_roi - sim_roi) / sigma
    chi2 = np.sum(residuals**2)
    
    # Degrees of freedom: n_points - n_parameters
    # We use n-1 as a simple approximation
    dof = len(raw_roi) - 1
    
    chi2_reduced = chi2 / dof if dof > 0 else np.inf
    
    return chi2_reduced, dof, total_counts


def calculate_uncertainty_for_element(concentration: float,
                                       roi_counts: int,
                                       technique: str = 'RBS',
                                       systematic_uncertainty: float = None) -> Dict[str, float]:
    """
    Calculate statistical and systematic uncertainties for an element.
    
    Statistical Uncertainty (Counting Statistics)
    ---------------------------------------------
    For an element with concentration C in a region with N_roi total counts,
    the element contributes approximately (C × N_roi) counts.
    
    The relative statistical uncertainty is:
        σ_stat / C = 1 / √(C × N_roi)
    
    Therefore:
        σ_stat = C / √(C × N_roi) = √(C / N_roi)
    
    Systematic Uncertainty
    ----------------------
    Arises from stopping power and cross-section uncertainties:
        σ_sys = C × systematic_uncertainty_relative
    
    Total Uncertainty
    -----------------
    Combined in quadrature (independent errors):
        σ_total = √(σ_stat² + σ_sys²)
    
    Parameters
    ----------
    concentration : float
        Element concentration as a fraction (0 to 1, not percentage)
    roi_counts : int
        Total counts in the ROI where this element's signal appears
    technique : str, default='RBS'
        IBA technique: 'RBS', 'ERDA', or 'NRA'
    systematic_uncertainty : float, optional
        Override the default systematic uncertainty (as fraction, e.g., 0.03 for 3%)
    
    Returns
    -------
    dict
        {
            'sigma_stat': statistical uncertainty (absolute, same units as C),
            'sigma_sys': systematic uncertainty (absolute),
            'sigma_total': total combined uncertainty (absolute)
        }
    """
    # Get systematic uncertainty
    if systematic_uncertainty is None:
        sigma_sys_rel = DEFAULT_SYSTEMATIC_UNCERTAINTIES.get(technique.upper(), 0.05)
    else:
        sigma_sys_rel = systematic_uncertainty
    
    # Statistical uncertainty from counting statistics
    # σ_stat = √(C / N_roi) = C / √(C × N_roi)
    if roi_counts > 0 and concentration > 0:
        effective_counts = roi_counts * concentration
        sigma_stat_rel = 1.0 / np.sqrt(max(effective_counts, 1))
        sigma_stat = concentration * sigma_stat_rel
    else:
        sigma_stat = np.nan
    
    # Systematic uncertainty
    sigma_sys = concentration * sigma_sys_rel
    
    # Combined uncertainty (quadrature sum)
    if not np.isnan(sigma_stat):
        sigma_total = np.sqrt(sigma_stat**2 + sigma_sys**2)
    else:
        sigma_total = sigma_sys
    
    return {
        'sigma_stat': sigma_stat,
        'sigma_sys': sigma_sys,
        'sigma_total': sigma_total,
    }


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def process_single_file_with_rois(filepath: str,
                                   regions_of_interest: List[Dict],
                                   technique: str = 'RBS',
                                   systematic_uncertainty: float = None) -> List[Dict]:
    """
    Process a single .xnra file using ROI-based analysis.
    
    Parameters
    ----------
    filepath : str
        Path to .xnra file
    regions_of_interest : list of dict
        Each ROI dict must contain:
        - 'name': str, ROI identifier
        - 'channels': tuple (start, end), channel range
        - 'elements': list of str, element symbols in this ROI
    technique : str
        'RBS', 'ERDA', or 'NRA'
    systematic_uncertainty : float, optional
        Override default systematic uncertainty
    
    Returns
    -------
    list of dict
        Results for each element in each layer in each ROI
    """
    data = parse_xnra_full(filepath)
    
    results = []
    
    for roi in regions_of_interest:
        roi_name = roi.get('name', f"ch{roi['channels'][0]}-{roi['channels'][1]}")
        roi_start, roi_end = roi['channels']
        roi_elements = roi['elements']
        
        # Calculate chi-squared for this ROI
        if data['raw_counts'] is not None and data['simulated_counts'] is not None:
            chi2_r, dof, roi_counts = calculate_chi_squared_for_roi(
                data['raw_counts'],
                data['simulated_counts'],
                roi_start,
                roi_end
            )
        else:
            chi2_r, dof, roi_counts = np.nan, 0, 0
        
        # Process each layer for elements in this ROI
        for layer in data['layers']:
            for element in roi_elements:
                concentration = layer['concentrations'].get(element, 0)
                
                if concentration <= 0:
                    continue
                
                unc = calculate_uncertainty_for_element(
                    concentration,
                    roi_counts,
                    technique,
                    systematic_uncertainty
                )
                
                results.append({
                    # File info
                    'filename': data['filename'],
                    'technique': technique,
                    # ROI info
                    'roi_name': roi_name,
                    'roi_start': roi_start,
                    'roi_end': roi_end,
                    'roi_counts': roi_counts,
                    'chi2_reduced': chi2_r,
                    'chi2_dof': dof,
                    # Layer info
                    'layer': layer['layer_num'],
                    'thickness_1e15_at_cm2': layer['thickness'],
                    # Element info
                    'element': element,
                    'concentration_fraction': concentration,
                    'concentration_percent': concentration * 100,
                    # Uncertainties (as fractions)
                    'sigma_stat_fraction': unc['sigma_stat'],
                    'sigma_sys_fraction': unc['sigma_sys'],
                    'sigma_total_fraction': unc['sigma_total'],
                    # Uncertainties (as percentages)
                    'sigma_stat_percent': unc['sigma_stat'] * 100 if not np.isnan(unc['sigma_stat']) else np.nan,
                    'sigma_sys_percent': unc['sigma_sys'] * 100,
                    'sigma_total_percent': unc['sigma_total'] * 100,
                })
    
    return results


def process_batch_with_rois(filepaths: List[str],
                             regions_of_interest: List[Dict],
                             technique: str = 'RBS',
                             systematic_uncertainty: float = None,
                             output_csv: str = None,
                             output_excel: str = None) -> pd.DataFrame:
    """
    Process multiple .xnra files using ROI-based analysis.
    
    Parameters
    ----------
    filepaths : list of str
        Paths to .xnra files
    regions_of_interest : list of dict
        ROI definitions
    technique : str
        'RBS', 'ERDA', or 'NRA'
    systematic_uncertainty : float, optional
        Override default systematic uncertainty
    output_csv : str, optional
        Path for CSV output
    output_excel : str, optional
        Path for Excel output
    
    Returns
    -------
    pd.DataFrame
        All results combined
    """
    all_results = []
    
    print(f"\nRegions of Interest:")
    for roi in regions_of_interest:
        name = roi.get('name', f"ch{roi['channels'][0]}-{roi['channels'][1]}")
        print(f"  • {name}: channels {roi['channels'][0]}-{roi['channels'][1]} → {roi['elements']}")
    print()
    
    for filepath in filepaths:
        try:
            results = process_single_file_with_rois(
                filepath, regions_of_interest, technique, systematic_uncertainty
            )
            all_results.extend(results)
            
            if results:
                print(f"✓ {Path(filepath).name}:")
                rois_seen = {}
                for r in results:
                    roi_name = r['roi_name']
                    if roi_name not in rois_seen:
                        rois_seen[roi_name] = r
                        print(f"    {roi_name}: χ²/ν = {r['chi2_reduced']:.2f}, counts = {r['roi_counts']}")
        except Exception as e:
            print(f"✗ Error processing {filepath}: {e}")
    
    df = pd.DataFrame(all_results)
    
    if not df.empty:
        df = df.sort_values(['filename', 'roi_name', 'layer', 'element'])
    
    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"\n✓ Saved CSV: {output_csv}")
    
    if output_excel:
        df.to_excel(output_excel, index=False, sheet_name='Uncertainties')
        print(f"✓ Saved Excel: {output_excel}")
    
    return df


# =============================================================================
# DEPTH PROFILE GENERATION
# =============================================================================

def create_depth_profile(df: pd.DataFrame,
                          sample_name: str = None,
                          elements_order: List[str] = None) -> pd.DataFrame:
    """
    Create a depth profile table from uncertainty results.
    
    The depth profile shows concentrations and uncertainties for each element
    at each depth (layer), formatted for publication or reporting.
    
    Parameters
    ----------
    df : pd.DataFrame
        Results DataFrame from process_batch_with_rois or calculate_uncertainties
    sample_name : str, optional
        Filter to specific sample (by filename)
    elements_order : list of str, optional
        Preferred order of elements in columns
    
    Returns
    -------
    pd.DataFrame
        Depth profile with columns:
        Layer, Thickness, Element1 (%), ±Element1, Element2 (%), ±Element2, ...
    """
    if df.empty:
        return pd.DataFrame()
    
    # Filter by sample if specified
    if sample_name:
        df = df[df['filename'] == sample_name]
    
    if df.empty:
        return pd.DataFrame()
    
    # Get unique layers and elements
    layers = sorted(df['layer'].unique())
    
    if elements_order is None:
        elements_order = sorted(df['element'].unique())
    else:
        # Only include elements that exist in data
        elements_order = [e for e in elements_order if e in df['element'].unique()]
    
    # Determine thickness column name (handle both conventions)
    thickness_col = None
    for col_name in ['thickness_1e15_at_cm2', 'thickness']:
        if col_name in df.columns:
            thickness_col = col_name
            break
    
    rows = []
    for layer_num in layers:
        layer_df = df[df['layer'] == layer_num]
        
        # Get thickness (should be same for all elements in layer)
        thickness = None
        if thickness_col and len(layer_df) > 0:
            thickness = layer_df[thickness_col].iloc[0]
        
        row = {
            'Layer': layer_num,
            'Thickness (10¹⁵ at/cm²)': thickness,
        }
        
        for element in elements_order:
            el_df = layer_df[layer_df['element'] == element]
            
            if len(el_df) > 0:
                # If multiple ROIs have this element, pick the one with best chi²
                best_idx = el_df['chi2_reduced'].idxmin()
                best = el_df.loc[best_idx]
                
                row[f'{element} (%)'] = round(best['concentration_percent'], 2)
                row[f'±{element} (%)'] = round(best['sigma_total_percent'], 3)
            else:
                row[f'{element} (%)'] = None
                row[f'±{element} (%)'] = None
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def create_depth_profiles_excel(df: pd.DataFrame,
                                 output_path: str,
                                 elements_order: List[str] = None) -> Dict[str, pd.DataFrame]:
    """
    Create Excel file with depth profiles for all samples.
    
    Parameters
    ----------
    df : pd.DataFrame
        Results DataFrame
    output_path : str
        Output Excel file path
    elements_order : list of str, optional
        Preferred element column order
    
    Returns
    -------
    dict
        {sample_name: depth_profile_dataframe}
    """
    profiles = {}
    
    # Create profile for each unique sample
    for sample in df['filename'].unique():
        profile = create_depth_profile(df, sample_name=sample, elements_order=elements_order)
        profiles[sample] = profile
    
    # Save to Excel with one sheet per sample
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        # Summary sheet
        summary_rows = []
        for sample, profile in profiles.items():
            sample_df = df[df['filename'] == sample]
            summary_rows.append({
                'Sample': sample,
                'Technique': sample_df['technique'].iloc[0] if len(sample_df) > 0 else '',
                'Layers': len(profile),
                'Elements': ', '.join(sorted(sample_df['element'].unique())),
            })
        
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Summary', index=False)
        
        # Individual sample sheets
        for sample, profile in profiles.items():
            sheet_name = sample[:31]  # Excel sheet name limit
            profile.to_excel(writer, sheet_name=sheet_name, index=False)
    
    print(f"✓ Saved depth profiles: {output_path}")
    
    return profiles


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def find_xnra_files(directory: str, pattern: str = "*.xnra") -> List[str]:
    """Find all .xnra files in a directory."""
    path = Path(directory)
    return sorted([str(f) for f in path.glob(pattern)])


def auto_detect_technique(filepath: str) -> str:
    """Auto-detect IBA technique from filename."""
    filename = Path(filepath).stem.upper()
    
    if 'ERDA' in filename:
        return 'ERDA'
    elif 'NRA' in filename:
        return 'NRA'
    else:
        return 'RBS'


# =============================================================================
# MAIN / CLI
# =============================================================================

def main():
    """Command-line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Calculate uncertainties for SIMNRA .xnra files (v5 - ROI-based)'
    )
    parser.add_argument('input', nargs='+', 
                        help='Input .xnra file(s) or directory')
    parser.add_argument('-r', '--roi', action='append', nargs=3,
                        metavar=('START', 'END', 'ELEMENTS'),
                        help='Define ROI: start_ch end_ch "El1,El2,..." (can repeat)')
    parser.add_argument('-t', '--technique', default='RBS',
                        choices=['RBS', 'ERDA', 'NRA'],
                        help='IBA technique (default: RBS)')
    parser.add_argument('-s', '--systematic', type=float,
                        help='Systematic uncertainty as fraction (e.g., 0.03 for 3%%)')
    parser.add_argument('-o', '--output', 
                        help='Output file (CSV or Excel based on extension)')
    parser.add_argument('--profiles', 
                        help='Output Excel for depth profiles')
    
    args = parser.parse_args()
    
    # Parse ROI definitions
    if not args.roi:
        print("Error: At least one ROI must be defined with -r START END ELEMENTS")
        print("Example: -r 380 420 'Si,Ca,O'")
        return
    
    regions_of_interest = []
    for start, end, elements in args.roi:
        regions_of_interest.append({
            'name': f"ch{start}-{end}",
            'channels': (int(start), int(end)),
            'elements': [e.strip() for e in elements.split(',')]
        })
    
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
    
    print(f"Processing {len(filepaths)} files with {args.technique}...")
    
    # Determine output format
    output_csv = None
    output_excel = None
    if args.output:
        if args.output.endswith('.xlsx'):
            output_excel = args.output
        else:
            output_csv = args.output if args.output.endswith('.csv') else args.output + '.csv'
    
    # Process
    df = process_batch_with_rois(
        filepaths, 
        regions_of_interest,
        technique=args.technique,
        systematic_uncertainty=args.systematic,
        output_csv=output_csv,
        output_excel=output_excel
    )
    
    print(f"\nProcessed {len(df)} element-layer-ROI combinations")
    
    # Create depth profiles if requested
    if args.profiles and not df.empty:
        create_depth_profiles_excel(df, args.profiles)


if __name__ == '__main__':
    main()
