"""
XNRA Uncertainty Calculator v3
==============================
Calculates statistical uncertainties for IBA data from SIMNRA .xnra files.

CHANGELOG v3:
- ROI-based analysis: define channel ranges and associated elements
- Each ROI gets its own chi² and uncertainty calculation
- Physically correct: only elements whose peaks are in the ROI get reported

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
# DATA STRUCTURES
# =============================================================================

# Example ROI configuration:
# REGIONS_OF_INTEREST = [
#     {'name': 'Heavy', 'channels': (700, 900), 'elements': ['Ba', 'Ca']},
#     {'name': 'Medium', 'channels': (400, 550), 'elements': ['Si', 'Al']},
#     {'name': 'Light', 'channels': (100, 250), 'elements': ['O', 'C']},
# ]


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
    Calculate reduced chi-squared for a specific ROI.
    
    Parameters:
        raw_counts: Experimental counts array
        simulated_counts: Simulated/fitted counts array
        roi_start: Start channel of ROI
        roi_end: End channel of ROI (exclusive)
        min_counts: Minimum counts for uncertainty (avoid div by zero)
    
    Returns:
        (chi2_reduced, degrees_of_freedom, total_counts_in_roi)
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
    
    # Total counts in ROI
    total_counts = int(np.sum(raw_roi))
    
    # Calculate uncertainty (Poisson: sqrt(N), minimum of min_counts)
    sigma = np.sqrt(np.maximum(raw_roi, min_counts))
    
    # Chi-squared
    residuals = (raw_roi - sim_roi) / sigma
    chi2 = np.sum(residuals**2)
    
    # Degrees of freedom (data points - 1)
    dof = len(raw_roi) - 1
    
    chi2_reduced = chi2 / dof if dof > 0 else np.inf
    
    return chi2_reduced, dof, total_counts


def calculate_uncertainty_for_element(concentration: float,
                                       roi_counts: int,
                                       technique: str = 'RBS') -> Dict[str, float]:
    """
    Calculate statistical and systematic uncertainties for an element.
    
    Parameters:
        concentration: Fitted concentration (fraction, not %)
        roi_counts: Total counts in the ROI for this element
        technique: 'RBS', 'ERDA', or 'NRA' (for systematic uncertainty)
    
    Returns:
        Dict with sigma_stat, sigma_sys, sigma_total (all as fractions)
    """
    # Statistical uncertainty from counting statistics in ROI
    # Element contributes proportionally to its concentration
    if roi_counts > 0 and concentration > 0:
        effective_counts = roi_counts * concentration
        sigma_stat_rel = 1.0 / np.sqrt(max(effective_counts, 1))
        sigma_stat = concentration * sigma_stat_rel
    else:
        sigma_stat = np.nan
    
    # Systematic uncertainty (literature values)
    systematic = {
        'RBS': 0.03,   # 3% - stopping powers
        'ERDA': 0.05,  # 5% - stopping powers + cross-sections
        'NRA': 0.05,   # 5% - reaction cross-sections
    }
    sigma_sys_rel = systematic.get(technique.upper(), 0.05)
    sigma_sys = concentration * sigma_sys_rel
    
    # Combined uncertainty
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
        return 'RBS'  # Default


def process_single_file_with_rois(filepath: str,
                                   regions_of_interest: List[Dict],
                                   technique: str = None) -> List[Dict]:
    """
    Process a single .xnra file using ROI-based analysis.
    
    Parameters:
        filepath: Path to .xnra file
        regions_of_interest: List of ROI definitions, each with:
            - 'name': ROI name (optional)
            - 'channels': (start, end) tuple
            - 'elements': list of element symbols
        technique: 'RBS', 'ERDA', 'NRA' (None = auto-detect)
    
    Returns:
        List of dicts with results per element per layer per ROI
    """
    data = parse_xnra_full(filepath)
    
    if technique is None:
        technique = detect_technique(filepath)
    
    results = []
    
    # Process each ROI
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
        
        # Process each layer, but only for elements in this ROI
        for layer in data['layers']:
            for element in roi_elements:
                # Get concentration for this element
                concentration = layer['concentrations'].get(element, 0)
                
                # Skip if element not in this layer
                if concentration <= 0:
                    continue
                
                # Calculate uncertainties using this ROI's counts
                unc = calculate_uncertainty_for_element(
                    concentration,
                    roi_counts,
                    technique
                )
                
                results.append({
                    'filename': data['filename'],
                    'filepath': data['filepath'],
                    'technique': technique,
                    'beam_energy_keV': data['beam_energy'],
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
                    # Uncertainties
                    'sigma_stat_fraction': unc['sigma_stat'],
                    'sigma_stat_percent': unc['sigma_stat'] * 100 if not np.isnan(unc['sigma_stat']) else np.nan,
                    'sigma_sys_fraction': unc['sigma_sys'],
                    'sigma_sys_percent': unc['sigma_sys'] * 100,
                    'sigma_total_fraction': unc['sigma_total'],
                    'sigma_total_percent': unc['sigma_total'] * 100,
                })
    
    return results


def process_batch_with_rois(filepaths: List[str],
                             regions_of_interest: List[Dict],
                             output_csv: str = None,
                             output_excel: str = None) -> pd.DataFrame:
    """
    Process multiple .xnra files using ROI-based analysis.
    
    Parameters:
        filepaths: List of paths to .xnra files
        regions_of_interest: List of ROI definitions
        output_csv: Path for CSV output (optional)
        output_excel: Path for Excel output (optional)
    
    Returns:
        DataFrame with all results
    """
    all_results = []
    
    print(f"\nRegions of Interest:")
    for roi in regions_of_interest:
        name = roi.get('name', f"ch{roi['channels'][0]}-{roi['channels'][1]}")
        print(f"  • {name}: channels {roi['channels'][0]}-{roi['channels'][1]} → {roi['elements']}")
    print()
    
    for filepath in filepaths:
        try:
            results = process_single_file_with_rois(filepath, regions_of_interest)
            all_results.extend(results)
            
            # Show summary per ROI
            if results:
                print(f"✓ {Path(filepath).name}:")
                # Group by ROI
                rois_seen = {}
                for r in results:
                    roi_name = r['roi_name']
                    if roi_name not in rois_seen:
                        rois_seen[roi_name] = r
                        print(f"    {roi_name}: χ²/ν = {r['chi2_reduced']:.2f}, counts = {r['roi_counts']}")
        except Exception as e:
            print(f"✗ Error processing {filepath}: {e}")
    
    df = pd.DataFrame(all_results)
    
    # Sort by filename, ROI, layer, element
    if not df.empty:
        df = df.sort_values(['filename', 'roi_name', 'layer', 'element'])
    
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
# DEPTH PROFILE TABLES
# =============================================================================

def create_depth_profile_table(df: pd.DataFrame,
                                roi_name: str = None,
                                sample_name: str = None) -> pd.DataFrame:
    """
    Create a depth profile table for a single sample/ROI combination.
    
    Rows = layers (depth)
    Columns = element concentrations and uncertainties
    """
    # Filter to specific sample and/or ROI
    filtered = df.copy()
    if sample_name:
        filtered = filtered[filtered['filename'] == sample_name]
    if roi_name:
        filtered = filtered[filtered['roi_name'] == roi_name]
    
    if filtered.empty:
        return pd.DataFrame()
    
    # Get elements in this ROI
    elements = sorted(filtered['element'].unique())
    
    # Get unique layers sorted
    layers = sorted(filtered['layer'].unique())
    
    # Get ROI info for header
    roi_info = filtered.iloc[0]
    
    rows = []
    for layer_num in layers:
        layer_df = filtered[filtered['layer'] == layer_num]
        
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
                               output_excel: str = None) -> Dict[str, pd.DataFrame]:
    """
    Create depth profile tables for all sample/ROI combinations.
    
    Returns dict of {sample_roi: DataFrame}
    Optionally saves to Excel with one sheet per combination.
    """
    profiles = {}
    
    for sample in df['filename'].unique():
        sample_df = df[df['filename'] == sample]
        
        for roi_name in sample_df['roi_name'].unique():
            key = f"{sample}_{roi_name}"
            profile = create_depth_profile_table(sample_df, roi_name=roi_name)
            profiles[key] = profile
    
    if output_excel:
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            for key, profile in profiles.items():
                # Sheet names max 31 chars
                sheet_name = key[:31]
                profile.to_excel(writer, sheet_name=sheet_name, index=False)
        print(f"✓ Saved depth profiles: {output_excel}")
    
    return profiles


# =============================================================================
# MAIN / CLI
# =============================================================================

def main():
    """
    Command-line interface for batch processing.
    """
    import argparse
    import json
    
    parser = argparse.ArgumentParser(
        description='Calculate uncertainties for SIMNRA .xnra files (v3 - ROI-based)'
    )
    parser.add_argument('input', nargs='+', 
                        help='Input .xnra file(s) or directory')
    parser.add_argument('-r', '--roi', action='append', nargs=3,
                        metavar=('START', 'END', 'ELEMENTS'),
                        help='Define ROI: start_ch end_ch "El1,El2,..." (can repeat)')
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
    
    print(f"Processing {len(filepaths)} files...")
    
    # Determine output format
    output_csv = None
    output_excel = None
    if args.output:
        if args.output.endswith('.xlsx'):
            output_excel = args.output
        else:
            output_csv = args.output if args.output.endswith('.csv') else args.output + '.csv'
    
    # Process
    df = process_batch_with_rois(filepaths, regions_of_interest, output_csv, output_excel)
    
    print(f"\nProcessed {len(df)} element-layer-ROI combinations")
    
    # Create depth profiles if requested
    if args.profiles and not df.empty:
        create_all_depth_profiles(df, args.profiles)


if __name__ == '__main__':
    main()
