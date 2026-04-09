"""
XNRA Uncertainty Calculator v4
==============================
Multi-technique IBA uncertainty analysis from SIMNRA .xnra files.

Combines RBS, ERDA, and NRA spectra to get complete elemental uncertainties:
  - RBS:  Heavy/medium elements (Ca, Ba, Si, Na, Fe, etc.)
  - ERDA: Light elements, especially H
  - NRA:  Light elements via nuclear reactions (C, O, N)

Each technique has its own ROIs and systematic uncertainties.

Authors: Berke Santos & Giuseppe Legrottaglie
"""

import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import re


# =============================================================================
# SYSTEMATIC UNCERTAINTIES BY TECHNIQUE
# =============================================================================

SYSTEMATIC_UNCERTAINTIES = {
    'RBS': 0.03,    # 3% - stopping powers
    'ERDA': 0.05,   # 5% - stopping powers + recoil cross-sections
    'NRA': 0.05,    # 5% - nuclear reaction cross-sections
}


# =============================================================================
# PARSING FUNCTIONS
# =============================================================================

def parse_xnra_full(filepath: str) -> Dict:
    """
    Parse .xnra file and extract all relevant data including layers.
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
            name_el = el.find('idf:name', ns)
            if name_el is None:
                name_el = el.find('idf:n', ns)
            conc_el = el.find('idf:concentration', ns)
            
            if name_el is not None and conc_el is not None:
                element_name = name_el.text
                concentration = float(conc_el.text)
                
                # Normalize isotope names (e.g., "40Ca" -> "Ca")
                base_element = re.sub(r'^\d+', '', element_name)
                
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
    
    Returns:
        (chi2_reduced, degrees_of_freedom, total_counts_in_roi)
    """
    roi_start = max(0, roi_start)
    roi_end = min(roi_end, len(raw_counts), len(simulated_counts))
    
    if roi_start >= roi_end:
        return np.nan, 0, 0
    
    raw_roi = raw_counts[roi_start:roi_end]
    sim_roi = simulated_counts[roi_start:roi_end]
    
    if len(raw_roi) == 0:
        return np.nan, 0, 0
    
    total_counts = int(np.sum(raw_roi))
    sigma = np.sqrt(np.maximum(raw_roi, min_counts))
    residuals = (raw_roi - sim_roi) / sigma
    chi2 = np.sum(residuals**2)
    dof = len(raw_roi) - 1
    chi2_reduced = chi2 / dof if dof > 0 else np.inf
    
    return chi2_reduced, dof, total_counts


def calculate_uncertainty_for_element(concentration: float,
                                       roi_counts: int,
                                       technique: str = 'RBS') -> Dict[str, float]:
    """
    Calculate statistical and systematic uncertainties for an element.
    """
    # Statistical uncertainty
    if roi_counts > 0 and concentration > 0:
        effective_counts = roi_counts * concentration
        sigma_stat_rel = 1.0 / np.sqrt(max(effective_counts, 1))
        sigma_stat = concentration * sigma_stat_rel
    else:
        sigma_stat = np.nan
    
    # Systematic uncertainty
    sigma_sys_rel = SYSTEMATIC_UNCERTAINTIES.get(technique.upper(), 0.05)
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
# MULTI-TECHNIQUE PROCESSING
# =============================================================================

def process_technique_file(filepath: str,
                           technique: str,
                           rois: List[Dict],
                           sample_name: str = None) -> List[Dict]:
    """
    Process a single .xnra file for a specific technique.
    
    Parameters:
        filepath: Path to .xnra file
        technique: 'RBS', 'ERDA', or 'NRA'
        rois: List of ROI definitions for this technique
        sample_name: Override sample name (for grouping multiple techniques)
    
    Returns:
        List of result dicts
    """
    data = parse_xnra_full(filepath)
    
    # Use provided sample name or derive from filename
    if sample_name is None:
        sample_name = data['filename']
    
    results = []
    
    for roi in rois:
        roi_name = roi.get('name', f"{technique}_{roi['channels'][0]}-{roi['channels'][1]}")
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
                    technique
                )
                
                results.append({
                    # Sample identification
                    'sample': sample_name,
                    'source_file': data['filename'],
                    # Technique info
                    'technique': technique,
                    'beam_energy_keV': data['beam_energy'],
                    'beam_particle': data['beam_particle'],
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


def process_multi_technique(analysis_config: List[Dict],
                            output_csv: str = None,
                            output_excel: str = None) -> pd.DataFrame:
    """
    Process multiple techniques for comprehensive uncertainty analysis.
    
    Parameters:
        analysis_config: List of technique configurations, each with:
            - 'file': path to .xnra file
            - 'technique': 'RBS', 'ERDA', or 'NRA'
            - 'sample': sample name (for grouping)
            - 'rois': list of ROI definitions
        output_csv: Path for CSV output
        output_excel: Path for Excel output
    
    Returns:
        DataFrame with all results
    """
    all_results = []
    
    print("\nAnalysis Configuration:")
    print("="*60)
    
    for config in analysis_config:
        filepath = config['file']
        technique = config['technique']
        sample = config.get('sample', Path(filepath).stem)
        rois = config['rois']
        
        print(f"\n{technique} - {Path(filepath).name} (sample: {sample})")
        for roi in rois:
            name = roi.get('name', f"ch{roi['channels'][0]}-{roi['channels'][1]}")
            print(f"  • {name}: ch {roi['channels'][0]}-{roi['channels'][1]} → {roi['elements']}")
        
        try:
            results = process_technique_file(filepath, technique, rois, sample)
            all_results.extend(results)
            
            # Show chi² per ROI
            rois_seen = {}
            for r in results:
                if r['roi_name'] not in rois_seen:
                    rois_seen[r['roi_name']] = r
                    status = "✓" if r['chi2_reduced'] < 5 else "⚠" if r['chi2_reduced'] < 20 else "✗"
                    print(f"    {status} {r['roi_name']}: χ²/ν = {r['chi2_reduced']:.2f}")
                    
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    df = pd.DataFrame(all_results)
    
    if not df.empty:
        df = df.sort_values(['sample', 'layer', 'element', 'technique'])
    
    # Save outputs
    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"\n✓ Saved CSV: {output_csv}")
    
    if output_excel:
        df.to_excel(output_excel, index=False, sheet_name='All_Results')
        print(f"✓ Saved Excel: {output_excel}")
    
    return df


# =============================================================================
# COMBINED DEPTH PROFILES
# =============================================================================

def create_combined_depth_profile(df: pd.DataFrame,
                                   sample_name: str) -> pd.DataFrame:
    """
    Create a combined depth profile for a sample using all techniques.
    
    Each element appears once with its best measurement (from appropriate technique).
    """
    sample_df = df[df['sample'] == sample_name]
    
    if sample_df.empty:
        return pd.DataFrame()
    
    # Get all layers and elements
    layers = sorted(sample_df['layer'].unique())
    elements = sorted(sample_df['element'].unique())
    
    rows = []
    for layer_num in layers:
        layer_df = sample_df[sample_df['layer'] == layer_num]
        
        # Get thickness (should be consistent across techniques)
        thickness = layer_df['thickness_1e15_at_cm2'].iloc[0]
        
        row = {
            'Layer': layer_num,
            'Thickness (10¹⁵ at/cm²)': thickness,
        }
        
        for element in elements:
            el_df = layer_df[layer_df['element'] == element]
            
            if len(el_df) > 0:
                # If multiple techniques measured this element, pick best chi²
                best_idx = el_df['chi2_reduced'].idxmin()
                best = el_df.loc[best_idx]
                
                row[f'{element} (%)'] = best['concentration_percent']
                row[f'±{element}'] = best['sigma_total_percent']
                row[f'{element}_technique'] = best['technique']
                row[f'{element}_chi2'] = best['chi2_reduced']
            else:
                row[f'{element} (%)'] = None
                row[f'±{element}'] = None
                row[f'{element}_technique'] = None
                row[f'{element}_chi2'] = None
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def create_all_combined_profiles(df: pd.DataFrame,
                                  output_excel: str = None) -> Dict[str, pd.DataFrame]:
    """
    Create combined depth profiles for all samples.
    """
    profiles = {}
    
    for sample in df['sample'].unique():
        profile = create_combined_depth_profile(df, sample)
        profiles[sample] = profile
    
    if output_excel:
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            # Summary sheet
            summary_rows = []
            for sample, profile in profiles.items():
                sample_df = df[df['sample'] == sample]
                techniques = sample_df['technique'].unique()
                elements = sample_df['element'].unique()
                summary_rows.append({
                    'Sample': sample,
                    'Techniques': ', '.join(techniques),
                    'Elements': ', '.join(sorted(elements)),
                    'Layers': len(profile),
                })
            pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Summary', index=False)
            
            # Individual sample sheets
            for sample, profile in profiles.items():
                sheet_name = sample[:31]
                profile.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print(f"✓ Saved combined profiles: {output_excel}")
    
    return profiles


# =============================================================================
# SIMPLIFIED DEPTH PROFILE (concentrations + uncertainties only)
# =============================================================================

def create_simple_depth_profile(df: pd.DataFrame,
                                 sample_name: str,
                                 elements_order: List[str] = None) -> pd.DataFrame:
    """
    Create a clean depth profile with just concentrations and uncertainties.
    
    Parameters:
        df: Full results DataFrame
        sample_name: Sample to extract
        elements_order: Preferred order of elements (optional)
    """
    sample_df = df[df['sample'] == sample_name]
    
    if sample_df.empty:
        return pd.DataFrame()
    
    layers = sorted(sample_df['layer'].unique())
    
    # Determine element order
    if elements_order is None:
        elements_order = sorted(sample_df['element'].unique())
    
    rows = []
    for layer_num in layers:
        layer_df = sample_df[sample_df['layer'] == layer_num]
        thickness = layer_df['thickness_1e15_at_cm2'].iloc[0]
        
        row = {
            'Layer': layer_num,
            'Thickness': thickness,
        }
        
        for element in elements_order:
            el_df = layer_df[layer_df['element'] == element]
            
            if len(el_df) > 0:
                # Pick measurement with best chi²
                best = el_df.loc[el_df['chi2_reduced'].idxmin()]
                row[f'{element}'] = best['concentration_percent']
                row[f'±{element}'] = best['sigma_total_percent']
            else:
                row[f'{element}'] = None
                row[f'±{element}'] = None
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def create_all_simple_profiles(df: pd.DataFrame,
                                elements_order: List[str] = None,
                                output_excel: str = None) -> Dict[str, pd.DataFrame]:
    """
    Create simple depth profiles for all samples.
    """
    profiles = {}
    
    for sample in df['sample'].unique():
        profile = create_simple_depth_profile(df, sample, elements_order)
        profiles[sample] = profile
    
    if output_excel:
        with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
            for sample, profile in profiles.items():
                sheet_name = sample[:31]
                profile.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print(f"✓ Saved simple profiles: {output_excel}")
    
    return profiles


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def find_xnra_files(directory: str, pattern: str = "*.xnra") -> List[str]:
    """Find all .xnra files in a directory."""
    path = Path(directory)
    return sorted([str(f) for f in path.glob(pattern)])


def auto_detect_technique(filepath: str) -> str:
    """Auto-detect technique from filename."""
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
        description='Multi-technique IBA uncertainty analysis (v4)'
    )
    parser.add_argument('--config', '-c',
                        help='JSON config file with analysis setup')
    parser.add_argument('-o', '--output',
                        help='Output Excel file')
    
    args = parser.parse_args()
    
    if args.config:
        import json
        with open(args.config) as f:
            config = json.load(f)
        
        df = process_multi_technique(config['analysis'], output_excel=args.output)
        
        if 'profiles_output' in config:
            create_all_simple_profiles(df, output_excel=config['profiles_output'])
    else:
        print("Usage: python xnra_uncertainty_v4.py -c config.json -o results.xlsx")
        print("\nSee xnra_process_uncertainties_v4.py for Python API usage.")


if __name__ == '__main__':
    main()
