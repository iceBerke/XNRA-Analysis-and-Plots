"""
XNRA Uncertainty Calculator v9
==============================
Calculates statistical uncertainties for IBA data from SIMNRA .xnra files.

This version implements the uncertainty formulas from the French IBA expert,
with proper separation of uncertainty components and optional χ²-adjustment.

VERSION 9 UPDATES:
- Added isotope-specific uncertainty calculations
- Clarified that Q·Ω and σ uncertainties are RELATIVE (constant across layers/elements)
- Isotopes can have different cross-section uncertainties based on their reactions
- Element spectrum parsing now preserves isotope data separately

=============================================================================
THEORETICAL BACKGROUND (from IBA-uncertainties-v3.pdf)
=============================================================================

1. DETECTED COUNTS RELATION
---------------------------
The number of detected counts for an element A can be expressed as:

    A_A = (Q · Ω · ∫[0→x] f · N_A(x) · σ_A(x) dx) / cos(α)

where:
- Q = collected charge
- Ω = detector solid angle  
- f = normalization factor
- N_A(x) = atomic concentration profile of element A
- σ_A(x) = reaction cross section
- α = incidence angle


2. STANDARD UNCERTAINTY ON CONCENTRATION (Equation 3)
-----------------------------------------------------
The uncertainty on concentration C_A combines independent uncertainties:

    ΔC_A = √[(ΔQΩ)² + (Δσ_A)² + (ΔA_A)²]

All terms are RELATIVE uncertainties (fractions, not absolute values).

IMPORTANT: ΔQΩ and Δσ are CONSTANT for all elements/layers in the same measurement.
Only ΔA_A (statistical) varies because it depends on element-specific counts.


3. χ²-ADJUSTED UNCERTAINTY (Equation 6)
---------------------------------------
When χ²/ν > 1, the statistical uncertainty may be underestimated. The
adjusted formula rescales the statistical term:

    ΔC_A^adj = √[(ΔQΩ)² + (Δσ_A)² + (√χ²/ν × ΔA_A)²]

This correction is only meaningful when χ²/ν significantly exceeds unity.
Otherwise, the standard expression should be used.


4. UNCERTAINTY COMPONENTS
-------------------------

a) ΔQΩ - Charge × Solid Angle Uncertainty (CONSTANT for all elements)
   - Estimated as standard deviation of ratios between simulated and 
     measured Particle·sr values during normalization
   - Valid for samples from the same analysis session
   - Typically 1-5% depending on setup stability
   - THIS IS A DETECTOR/SETUP PROPERTY - SAME FOR ALL ELEMENTS

b) Δσ_A - Cross-Section Uncertainty (varies by REACTION, not element)
   - Depends on the specific nuclear reaction:
     • 0%  for Mx(α,α)Mx (pure Rutherford backscattering)
     • ±2% for ¹⁶O(d,α₀)¹⁴N (deuteron NRA for oxygen)
     • ±2% for ¹H(⁴He,¹H)⁴He (forward recoil/ERDA for hydrogen)
     • ±5% for ¹²C(d,p₀)¹³C (deuteron NRA for carbon)
   - For isotope-specific calculations, each isotope may have different Δσ

c) ΔA_A - Statistical Uncertainty (Poisson) - VARIES by element
   - Given by: ΔA_A = √A_A / A_A = 1/√A_A
   - Where A_A = detected counts for element/isotope A at its reference peak
   - Calculated within a user-defined ROI (Region of Interest)


5. REDUCED CHI-SQUARED (χ²/ν)
-----------------------------
Defined as:

    χ²/ν = (1/ν) × Σ[(y_exp - y_sim)² / σ²]

where:
- y_exp, y_sim = experimental and simulated values
- σ = statistical uncertainty (√y_exp)
- ν = N - p = degrees of freedom (N = data points, p = fitted parameters)


=============================================================================
AUTHORS
=============================================================================
Based on work by Berke Santos & Giuseppe Legrottaglie
Formulas from French IBA expert documentation (v3)
Version 9: Isotope support + clarified constant Q·Ω/σ uncertainties

=============================================================================
"""

import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import re

# Import element spectrum functions from dedicated module
from xnra_element_spectra_v1 import (
    parse_element_spectrum,
    get_available_elements,
    get_available_isotopes,
    get_element_spectrum_summary,
    calculate_element_counts_from_spectrum,
    plot_element_spectra,
    ELEMENT_COLORS,
    COLOR_NAMES,
)


# =============================================================================
# CROSS-SECTION UNCERTAINTIES (from IBA-uncertainties-v3.pdf)
# =============================================================================

CROSS_SECTION_UNCERTAINTIES = {
    # Reaction name: (uncertainty as fraction, description)
    # From IBA-uncertainties-v3.pdf (French IBA expert)
    'RBS: Mx(α,α)Mx': (0.00, 'Rutherford backscattering (pure Coulomb)'),
    'NRA: ¹⁶O(d,α₀)¹⁴N': (0.02, 'Deuteron NRA for oxygen'),
    'NRA: ¹²C(d,p₀)¹³C': (0.05, 'Deuteron NRA for carbon'),
    'NRA: ¹⁵N(p,αγ)¹²C': (0.03, 'Proton NRA for nitrogen'),
    'ERDA: ¹H(⁴He,¹H)⁴He': (0.02, 'Forward recoil for hydrogen'),
    'ERDA: ²H(⁴He,²H)⁴He': (0.03, 'Forward recoil for deuterium'),
}

# Default cross-section uncertainties by isotope (for common IBA scenarios)
# Users can override these in the GUI
ISOTOPE_CROSS_SECTION_DEFAULTS = {
    # Format: 'isotope': (default_uncertainty, default_reaction)
    # Hydrogen isotopes (typically ERDA)
    '1H': (0.02, 'ERDA: ¹H(⁴He,¹H)⁴He'),
    '2H': (0.03, 'ERDA: ²H(⁴He,²H)⁴He'),
    'H': (0.02, 'ERDA: ¹H(⁴He,¹H)⁴He'),
    
    # Carbon isotopes (NRA or RBS depending on setup)
    '12C': (0.05, 'NRA: ¹²C(d,p₀)¹³C'),
    '13C': (0.05, 'NRA: ¹²C(d,p₀)¹³C'),  # Same reaction uncertainty
    'C': (0.05, 'NRA: ¹²C(d,p₀)¹³C'),
    
    # Nitrogen isotopes
    '14N': (0.03, 'NRA: ¹⁵N(p,αγ)¹²C'),
    '15N': (0.03, 'NRA: ¹⁵N(p,αγ)¹²C'),
    'N': (0.03, 'NRA: ¹⁵N(p,αγ)¹²C'),
    
    # Oxygen isotopes (NRA or RBS)
    '16O': (0.02, 'NRA: ¹⁶O(d,α₀)¹⁴N'),
    '17O': (0.02, 'NRA: ¹⁶O(d,α₀)¹⁴N'),
    '18O': (0.02, 'NRA: ¹⁶O(d,α₀)¹⁴N'),
    'O': (0.02, 'NRA: ¹⁶O(d,α₀)¹⁴N'),
    
    # Heavy elements (typically RBS - Rutherford)
    # These have 0% cross-section uncertainty for pure Coulomb scattering
}

# Legacy mapping for backward compatibility
TECHNIQUE_TO_REACTION = {
    'RBS': 'RBS: Mx(α,α)Mx',
    'ERDA': 'ERDA: ¹H(⁴He,¹H)⁴He',
    'NRA': 'NRA: ¹²C(d,p₀)¹³C',  # Default NRA
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
            'concentrations': {},
            'isotope_concentrations': {},  # NEW: Keep isotope-specific data
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
                
                # Store isotope-specific concentration
                layer_data['isotope_concentrations'][element_name] = concentration
                
                # Also store summed element concentration (backward compatible)
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
# UNCERTAINTY CALCULATIONS (v9 - French Expert Formulas + Isotope Support)
# =============================================================================

def calculate_chi_squared_for_roi(raw_counts: np.ndarray, 
                                   simulated_counts: np.ndarray,
                                   roi_start: int,
                                   roi_end: int,
                                   min_counts: float = 4.0,
                                   n_params: int = 1) -> Tuple[float, int, int, float]:
    """
    Calculate reduced chi-squared and AUC ratio for a specific ROI.
    
    Formula (Equation 5 from IBA-uncertainties-v3.pdf):
    
        χ²/ν = (1/ν) × Σ[(y_exp - y_sim)² / σ²]
    
    where:
    - y_exp, y_sim = experimental and simulated values
    - σ = √(max(y_exp, min_counts)) = Poisson uncertainty
    - ν = N - p = degrees of freedom
    
    Also calculates AUC ratio = Σ(y_sim) / Σ(y_exp) for normalization check.
    
    Parameters
    ----------
    raw_counts : np.ndarray
        Experimental spectrum counts (y_exp)
    simulated_counts : np.ndarray
        Simulated/fitted spectrum counts (y_sim)
    roi_start : int
        Start channel of ROI (inclusive)
    roi_end : int
        End channel of ROI (exclusive)
    min_counts : float, default=4.0
        Minimum counts for uncertainty calculation (avoids division by zero)
    n_params : int, default=1
        Number of fitted parameters (p) for degrees of freedom calculation
    
    Returns
    -------
    tuple
        (chi2_reduced, degrees_of_freedom, total_counts_in_roi, auc_ratio)
    """
    # Ensure bounds are valid
    roi_start = max(0, roi_start)
    roi_end = min(roi_end, len(raw_counts), len(simulated_counts))
    
    if roi_start >= roi_end:
        return np.nan, 0, 0, np.nan
    
    # Extract ROI
    raw_roi = raw_counts[roi_start:roi_end]
    sim_roi = simulated_counts[roi_start:roi_end]
    
    if len(raw_roi) == 0:
        return np.nan, 0, 0, np.nan
    
    # Total counts in ROI (for counting statistics)
    total_counts = int(np.sum(raw_roi))
    total_sim_counts = np.sum(sim_roi)
    
    # AUC ratio: simple normalization check
    auc_ratio = total_sim_counts / total_counts if total_counts > 0 else np.nan
    
    # Poisson uncertainty: σ = √N, with minimum to avoid div/0
    sigma = np.sqrt(np.maximum(raw_roi, min_counts))
    
    # Chi-squared: Σ[(y_exp - y_sim)² / σ²]
    chi2 = np.sum((raw_roi - sim_roi)**2 / sigma**2)
    
    # Degrees of freedom: ν = N - p
    n_points = len(raw_roi)
    dof = n_points - n_params
    
    chi2_reduced = chi2 / dof if dof > 0 else np.inf
    
    return chi2_reduced, dof, total_counts, auc_ratio


def interpret_chi2(chi2_reduced: float) -> Tuple[str, str]:
    """
    Interpret the reduced chi-squared value.
    
    Parameters
    ----------
    chi2_reduced : float
        The reduced chi-squared value
    
    Returns
    -------
    tuple
        (status_symbol, description)
    """
    if np.isnan(chi2_reduced):
        return "?", "N/A"
    elif 0.90 <= chi2_reduced <= 1.10:
        return "✓✓", "Perfect"
    elif 1.10 < chi2_reduced <= 3.0:
        return "✓", "Great"
    elif 3.0 < chi2_reduced <= 10.0:
        return "⚠", "Not OK"
    elif chi2_reduced > 10.0:
        return "✗", "Very bad"
    elif chi2_reduced < 0.90:
        return "⚠", "Too low"
    else:
        return "?", "Unknown"


# =============================================================================
# NOTE: Element spectrum functions moved to xnra_element_spectra.py
# Import: parse_element_spectrum, get_element_spectrum_summary, 
#         plot_element_spectra, calculate_element_counts_from_spectrum, etc.
# =============================================================================


def get_isotope_cross_section_uncertainty(isotope: str, 
                                           custom_overrides: Dict[str, float] = None) -> Tuple[float, str]:
    """
    Get the cross-section uncertainty for a specific isotope.
    
    Parameters
    ----------
    isotope : str
        Isotope name (e.g., '12C', '16O', 'H')
    custom_overrides : dict, optional
        Dict mapping isotope name → uncertainty (as fraction)
    
    Returns
    -------
    tuple
        (uncertainty as fraction, reaction description)
    """
    # Check custom overrides first
    if custom_overrides and isotope in custom_overrides:
        return custom_overrides[isotope], 'Custom'
    
    # Check default isotope uncertainties
    if isotope in ISOTOPE_CROSS_SECTION_DEFAULTS:
        unc, reaction = ISOTOPE_CROSS_SECTION_DEFAULTS[isotope]
        return unc, reaction
    
    # Try stripping mass number to find element default
    base_element = re.sub(r'^\d+', '', isotope)
    if base_element in ISOTOPE_CROSS_SECTION_DEFAULTS:
        unc, reaction = ISOTOPE_CROSS_SECTION_DEFAULTS[base_element]
        return unc, reaction
    
    # Default to Rutherford (0% uncertainty for heavy elements)
    return 0.0, 'RBS: Mx(α,α)Mx'


def calculate_uncertainty_v9(concentration: float,
                              element_counts: int,
                              qomega_uncertainty: float,
                              cross_section_uncertainty: float,
                              chi2_reduced: float = 1.0,
                              auc_ratio: float = 1.0,
                              adjustment_method: str = 'none') -> Dict[str, float]:
    """
    Calculate uncertainties using the French expert's formulas (v9).
    
    Standard Formula (Equation 3):
    ------------------------------
        ΔC_A = √[(ΔQΩ)² + (Δσ_A)² + (ΔA_A)²]
    
    IMPORTANT: ΔQΩ and Δσ_A are RELATIVE uncertainties that are CONSTANT
    for all elements/layers. Only ΔA_A varies based on counting statistics.
    
    Parameters
    ----------
    concentration : float
        Element concentration as a fraction (0 to 1, not percentage)
    element_counts : int
        Number of detected counts for element A in the ROI (A_A)
    qomega_uncertainty : float
        Q·Ω uncertainty as fraction (e.g., 0.03 for 3%) - CONSTANT for all elements
    cross_section_uncertainty : float
        Cross-section uncertainty as fraction (e.g., 0.02 for 2%) - per reaction
    chi2_reduced : float, default=1.0
        Reduced chi-squared value for the ROI
    auc_ratio : float, default=1.0
        AUC ratio (Σsim/Σexp) for the ROI
    adjustment_method : str, default='none'
        Which adjustment to apply: 'none', 'chi2', or 'auc'
    
    Returns
    -------
    dict
        Contains standard and all adjusted uncertainties for comparison.
        
        KEY DISTINCTION:
        - sigma_*_rel: RELATIVE uncertainties (fractions) - use for display
        - sigma_*: ABSOLUTE uncertainties (fractions of concentration)
    """
    # Statistical uncertainty (Poisson): ΔA_A = 1/√A_A
    if element_counts > 0:
        sigma_stat_rel = 1.0 / np.sqrt(element_counts)
    else:
        sigma_stat_rel = np.nan
    
    # Store relative uncertainties - THESE ARE CONSTANT (except stat)
    sigma_qomega_rel = qomega_uncertainty   # CONSTANT - detector property
    sigma_xsec_rel = cross_section_uncertainty  # CONSTANT per reaction type
    
    # χ² adjustment factor: √(max(χ²/ν, 1))
    if not np.isnan(chi2_reduced):
        chi2_factor = np.sqrt(max(chi2_reduced, 1.0))
    else:
        chi2_factor = 1.0
    
    # AUC adjustment factor: max(AUC, 1/AUC)
    if not np.isnan(auc_ratio) and auc_ratio > 0:
        auc_factor = max(auc_ratio, 1.0 / auc_ratio)
    else:
        auc_factor = 1.0
    
    # Calculate all versions of statistical uncertainty (relative)
    if not np.isnan(sigma_stat_rel):
        sigma_stat_chi2_rel = chi2_factor * sigma_stat_rel
        sigma_stat_auc_rel = auc_factor * sigma_stat_rel
    else:
        sigma_stat_chi2_rel = np.nan
        sigma_stat_auc_rel = np.nan
    
    # Standard combined RELATIVE uncertainty (Equation 3)
    if not np.isnan(sigma_stat_rel):
        sigma_total_rel = np.sqrt(
            sigma_qomega_rel**2 + 
            sigma_xsec_rel**2 + 
            sigma_stat_rel**2
        )
    else:
        sigma_total_rel = np.sqrt(sigma_qomega_rel**2 + sigma_xsec_rel**2)
    
    # χ²-adjusted combined uncertainty (Equation 6)
    if not np.isnan(sigma_stat_chi2_rel):
        sigma_total_chi2_rel = np.sqrt(
            sigma_qomega_rel**2 + 
            sigma_xsec_rel**2 + 
            sigma_stat_chi2_rel**2
        )
    else:
        sigma_total_chi2_rel = np.sqrt(sigma_qomega_rel**2 + sigma_xsec_rel**2)
    
    # AUC-adjusted combined uncertainty
    if not np.isnan(sigma_stat_auc_rel):
        sigma_total_auc_rel = np.sqrt(
            sigma_qomega_rel**2 + 
            sigma_xsec_rel**2 + 
            sigma_stat_auc_rel**2
        )
    else:
        sigma_total_auc_rel = np.sqrt(sigma_qomega_rel**2 + sigma_xsec_rel**2)
    
    # Convert to absolute uncertainties (scaled by concentration)
    # These represent the actual uncertainty in concentration percentage
    sigma_stat_abs = concentration * sigma_stat_rel if not np.isnan(sigma_stat_rel) else np.nan
    sigma_stat_chi2_abs = concentration * sigma_stat_chi2_rel if not np.isnan(sigma_stat_chi2_rel) else np.nan
    sigma_stat_auc_abs = concentration * sigma_stat_auc_rel if not np.isnan(sigma_stat_auc_rel) else np.nan
    sigma_qomega_abs = concentration * sigma_qomega_rel
    sigma_xsec_abs = concentration * sigma_xsec_rel
    sigma_total_abs = concentration * sigma_total_rel
    sigma_total_chi2_abs = concentration * sigma_total_chi2_rel
    sigma_total_auc_abs = concentration * sigma_total_auc_rel
    
    # Choose which to report based on adjustment method
    if adjustment_method == 'chi2':
        sigma_stat_final_abs = sigma_stat_chi2_abs
        sigma_stat_final_rel = sigma_stat_chi2_rel
        sigma_total_final_abs = sigma_total_chi2_abs
        sigma_total_final_rel = sigma_total_chi2_rel
        adjustment_factor = chi2_factor
    elif adjustment_method == 'auc':
        sigma_stat_final_abs = sigma_stat_auc_abs
        sigma_stat_final_rel = sigma_stat_auc_rel
        sigma_total_final_abs = sigma_total_auc_abs
        sigma_total_final_rel = sigma_total_auc_rel
        adjustment_factor = auc_factor
    else:  # 'none' or default
        sigma_stat_final_abs = sigma_stat_abs
        sigma_stat_final_rel = sigma_stat_rel
        sigma_total_final_abs = sigma_total_abs
        sigma_total_final_rel = sigma_total_rel
        adjustment_factor = 1.0
    
    return {
        # ============================================
        # PRIMARY UNCERTAINTIES (based on selected method)
        # ============================================
        # Absolute (concentration-scaled)
        'sigma_stat': sigma_stat_final_abs,
        'sigma_qomega': sigma_qomega_abs,
        'sigma_xsec': sigma_xsec_abs,
        'sigma_total': sigma_total_final_abs,
        
        # Relative (CONSTANT for Q·Ω and σ) - USE THESE FOR DISPLAY
        'sigma_stat_rel': sigma_stat_final_rel,
        'sigma_qomega_rel': sigma_qomega_rel,  # CONSTANT
        'sigma_xsec_rel': sigma_xsec_rel,       # CONSTANT
        'sigma_total_rel': sigma_total_final_rel,
        
        # ============================================
        # ALL VERSIONS FOR COMPARISON (absolute)
        # ============================================
        'sigma_stat_standard': sigma_stat_abs,
        'sigma_stat_chi2_adjusted': sigma_stat_chi2_abs,
        'sigma_stat_auc_adjusted': sigma_stat_auc_abs,
        'sigma_total_standard': sigma_total_abs,
        'sigma_total_chi2_adjusted': sigma_total_chi2_abs,
        'sigma_total_auc_adjusted': sigma_total_auc_abs,
        
        # All versions (relative)
        'sigma_stat_standard_rel': sigma_stat_rel,
        'sigma_stat_chi2_adjusted_rel': sigma_stat_chi2_rel,
        'sigma_stat_auc_adjusted_rel': sigma_stat_auc_rel,
        'sigma_total_standard_rel': sigma_total_rel,
        'sigma_total_chi2_adjusted_rel': sigma_total_chi2_rel,
        'sigma_total_auc_adjusted_rel': sigma_total_auc_rel,
        
        # ============================================
        # ADJUSTMENT FACTORS
        # ============================================
        'chi2_factor': chi2_factor,
        'auc_factor': auc_factor,
        'adjustment_factor': adjustment_factor,
        'adjustment_method': adjustment_method,
    }


# Alias for backward compatibility
calculate_uncertainty_v6 = calculate_uncertainty_v9


def get_cross_section_uncertainty(reaction: str) -> float:
    """
    Get the cross-section uncertainty for a given reaction.
    
    Parameters
    ----------
    reaction : str
        Reaction name (must match keys in CROSS_SECTION_UNCERTAINTIES)
        or a technique name ('RBS', 'ERDA', 'NRA')
    
    Returns
    -------
    float
        Cross-section uncertainty as a fraction
    """
    # Direct lookup
    if reaction in CROSS_SECTION_UNCERTAINTIES:
        return CROSS_SECTION_UNCERTAINTIES[reaction][0]
    
    # Technique mapping
    if reaction.upper() in TECHNIQUE_TO_REACTION:
        mapped_reaction = TECHNIQUE_TO_REACTION[reaction.upper()]
        return CROSS_SECTION_UNCERTAINTIES[mapped_reaction][0]
    
    # Default fallback
    return 0.05


# =============================================================================
# BATCH PROCESSING
# =============================================================================

def process_single_file_v9(filepath: str,
                            regions_of_interest: List[Dict],
                            qomega_uncertainty: float,
                            reaction: str = 'RBS: Mx(α,α)Mx',
                            cross_section_override: float = None,
                            adjustment_method: str = 'none',
                            element_spectrum_data: Dict = None,
                            isotope_mode: bool = False,
                            isotope_xsec_overrides: Dict[str, float] = None) -> List[Dict]:
    """
    Process a single .xnra file using v9 uncertainty formulas.
    
    Parameters
    ----------
    filepath : str
        Path to .xnra file
    regions_of_interest : list of dict
        Each ROI dict must contain:
        - 'name': str, ROI identifier
        - 'channels': tuple (start, end), channel range
        - 'elements': list of str, element symbols OR isotope names
    qomega_uncertainty : float
        Q·Ω uncertainty as fraction (e.g., 0.03 for 3%) - CONSTANT
    reaction : str
        Default reaction name for cross-section lookup
    cross_section_override : float, optional
        Override cross-section uncertainty (as fraction)
    adjustment_method : str, default='none'
        Adjustment method: 'none', 'chi2', or 'auc'
    element_spectrum_data : dict, optional
        Parsed element spectrum data for accurate A_A calculation
    isotope_mode : bool, default=False
        If True, process isotopes separately instead of summing to elements
    isotope_xsec_overrides : dict, optional
        Per-isotope cross-section overrides {isotope: uncertainty}
    
    Returns
    -------
    list of dict
        Results for each element/isotope in each layer in each ROI
    """
    data = parse_xnra_full(filepath)
    
    # Default cross-section uncertainty (may be overridden per isotope)
    if cross_section_override is not None:
        default_xsec_unc = cross_section_override
    else:
        default_xsec_unc = get_cross_section_uncertainty(reaction)
    
    results = []
    
    for roi in regions_of_interest:
        roi_name = roi.get('name', f"ch{roi['channels'][0]}-{roi['channels'][1]}")
        roi_start, roi_end = roi['channels']
        roi_elements = roi['elements']
        
        # Calculate chi-squared for this ROI
        if data['raw_counts'] is not None and data['simulated_counts'] is not None:
            chi2_r, dof, roi_counts, auc_ratio = calculate_chi_squared_for_roi(
                data['raw_counts'],
                data['simulated_counts'],
                roi_start,
                roi_end
            )
            chi2_symbol, chi2_status = interpret_chi2(chi2_r)
        else:
            chi2_r, dof, roi_counts, auc_ratio = np.nan, 0, 0, np.nan
            chi2_symbol, chi2_status = "?", "N/A"
        
        # Process each layer for elements/isotopes in this ROI
        for layer in data['layers']:
            for element_or_isotope in roi_elements:
                # Determine if this is an isotope or element
                is_isotope = element_or_isotope[0].isdigit() if element_or_isotope else False
                
                if isotope_mode and is_isotope:
                    # Use isotope-specific concentration
                    concentration = layer['isotope_concentrations'].get(element_or_isotope, 0)
                    # Get isotope-specific cross-section uncertainty
                    xsec_unc, xsec_reaction = get_isotope_cross_section_uncertainty(
                        element_or_isotope, isotope_xsec_overrides
                    )
                else:
                    # Use summed element concentration
                    base_element = re.sub(r'^\d+', '', element_or_isotope)
                    concentration = layer['concentrations'].get(base_element, 0)
                    xsec_unc = default_xsec_unc
                    xsec_reaction = reaction
                
                if concentration <= 0:
                    continue
                
                # Calculate element counts
                if element_spectrum_data is not None and data['raw_counts'] is not None:
                    element_counts, element_fraction, total_exp = calculate_element_counts_from_spectrum(
                        element_spectrum_data,
                        data['raw_counts'],
                        element_or_isotope,
                        roi_start,
                        roi_end,
                        use_isotope=(isotope_mode and is_isotope)
                    )
                    aa_method = 'spectrum'
                else:
                    element_counts = int(roi_counts * concentration)
                    element_fraction = concentration
                    aa_method = 'approximate'
                
                unc = calculate_uncertainty_v9(
                    concentration,
                    element_counts,
                    qomega_uncertainty,
                    xsec_unc,
                    chi2_reduced=chi2_r,
                    auc_ratio=auc_ratio,
                    adjustment_method=adjustment_method
                )
                
                results.append({
                    # File info
                    'filename': data['filename'],
                    'reaction': xsec_reaction,
                    'is_isotope': is_isotope,
                    # ROI info
                    'roi_name': roi_name,
                    'roi_start': roi_start,
                    'roi_end': roi_end,
                    'roi_counts': roi_counts,
                    'element_counts': element_counts,
                    'element_fraction': element_fraction,
                    'aa_method': aa_method,
                    # Fit quality metrics
                    'chi2_reduced': chi2_r,
                    'chi2_dof': dof,
                    'chi2_symbol': chi2_symbol,
                    'chi2_status': chi2_status,
                    'auc_ratio': auc_ratio,
                    # Adjustment info
                    'adjustment_method': adjustment_method,
                    'chi2_factor': unc['chi2_factor'],
                    'auc_factor': unc['auc_factor'],
                    'adjustment_factor': unc['adjustment_factor'],
                    # Layer info
                    'layer': layer['layer_num'],
                    'thickness_1e15_at_cm2': layer['thickness'],
                    # Element/isotope info
                    'element': element_or_isotope,
                    'concentration_fraction': concentration,
                    'concentration_percent': concentration * 100,
                    # Input uncertainties (CONSTANT)
                    'qomega_input_percent': qomega_uncertainty * 100,
                    'xsec_input_percent': xsec_unc * 100,
                    # Primary RELATIVE uncertainties (Q·Ω and σ are CONSTANT)
                    'sigma_stat_rel_percent': unc['sigma_stat_rel'] * 100 if not np.isnan(unc['sigma_stat_rel']) else np.nan,
                    'sigma_qomega_rel_percent': unc['sigma_qomega_rel'] * 100,  # CONSTANT
                    'sigma_xsec_rel_percent': unc['sigma_xsec_rel'] * 100,      # CONSTANT
                    'sigma_total_rel_percent': unc['sigma_total_rel'] * 100,
                    # Primary ABSOLUTE uncertainties (concentration-scaled)
                    'sigma_stat_percent': unc['sigma_stat'] * 100 if not np.isnan(unc['sigma_stat']) else np.nan,
                    'sigma_qomega_percent': unc['sigma_qomega'] * 100,
                    'sigma_xsec_percent': unc['sigma_xsec'] * 100,
                    'sigma_total_percent': unc['sigma_total'] * 100,
                    # All versions for comparison
                    'sigma_stat_standard_percent': unc['sigma_stat_standard'] * 100 if not np.isnan(unc['sigma_stat_standard']) else np.nan,
                    'sigma_stat_chi2_adjusted_percent': unc['sigma_stat_chi2_adjusted'] * 100 if not np.isnan(unc['sigma_stat_chi2_adjusted']) else np.nan,
                    'sigma_stat_auc_adjusted_percent': unc['sigma_stat_auc_adjusted'] * 100 if not np.isnan(unc['sigma_stat_auc_adjusted']) else np.nan,
                    'sigma_total_standard_percent': unc['sigma_total_standard'] * 100,
                    'sigma_total_chi2_adjusted_percent': unc['sigma_total_chi2_adjusted'] * 100,
                    'sigma_total_auc_adjusted_percent': unc['sigma_total_auc_adjusted'] * 100,
                })
    
    return results


# Alias for backward compatibility
process_single_file_v6 = process_single_file_v9


def process_batch_v9(filepaths: List[str],
                      regions_of_interest: List[Dict],
                      qomega_uncertainty: float,
                      reaction: str = 'RBS: Mx(α,α)Mx',
                      cross_section_override: float = None,
                      adjustment_method: str = 'none',
                      element_spectrum_data: Dict = None,
                      isotope_mode: bool = False,
                      isotope_xsec_overrides: Dict[str, float] = None,
                      output_csv: str = None,
                      output_excel: str = None) -> pd.DataFrame:
    """
    Process multiple .xnra files using v9 uncertainty formulas.
    """
    all_results = []
    
    default_xsec_unc = cross_section_override if cross_section_override else get_cross_section_uncertainty(reaction)
    
    method_names = {
        'none': 'Standard (no adjustment)',
        'chi2': 'χ²-adjusted',
        'auc': 'AUC-adjusted'
    }
    method_formulas = {
        'none': 'ΔC = √[(ΔQ·Ω)² + (Δσ)² + (1/√N)²]',
        'chi2': 'ΔC = √[(ΔQ·Ω)² + (Δσ)² + (√χ²/ν × 1/√N)²]',
        'auc': 'ΔC = √[(ΔQ·Ω)² + (Δσ)² + (max(AUC,1/AUC) × 1/√N)²]'
    }
    
    print(f"\n{'='*60}")
    print(f"IBA Uncertainty Analysis (v9 - Isotope Support)")
    print(f"{'='*60}")
    print(f"\nInput Parameters:")
    print(f"  • Q·Ω uncertainty:     {qomega_uncertainty*100:.2f}% (CONSTANT)")
    print(f"  • Default reaction:    {reaction}")
    print(f"  • Default σ:           {default_xsec_unc*100:.2f}%")
    print(f"  • Isotope mode:        {'Enabled' if isotope_mode else 'Disabled'}")
    print(f"  • Adjustment method:   {method_names.get(adjustment_method, adjustment_method)}")
    print(f"\nFormula: {method_formulas.get(adjustment_method, method_formulas['none'])}")
    
    if isotope_xsec_overrides:
        print(f"\nIsotope-specific cross-sections:")
        for iso, unc in isotope_xsec_overrides.items():
            print(f"  • {iso}: {unc*100:.2f}%")
    
    print(f"\nRegions of Interest:")
    for roi in regions_of_interest:
        name = roi.get('name', f"ch{roi['channels'][0]}-{roi['channels'][1]}")
        print(f"  • {name}: channels {roi['channels'][0]}-{roi['channels'][1]} → {roi['elements']}")
    print()
    
    for filepath in filepaths:
        try:
            results = process_single_file_v9(
                filepath,
                regions_of_interest,
                qomega_uncertainty,
                reaction,
                cross_section_override,
                adjustment_method,
                element_spectrum_data,
                isotope_mode,
                isotope_xsec_overrides
            )
            all_results.extend(results)
            
            if results:
                print(f"✓ {Path(filepath).name}:")
                rois_seen = {}
                for r in results:
                    roi_name = r['roi_name']
                    if roi_name not in rois_seen:
                        rois_seen[roi_name] = r
                        chi2_str = f"{r['chi2_reduced']:.2f}" if not np.isnan(r['chi2_reduced']) else "N/A"
                        auc_str = f"{r['auc_ratio']:.3f}" if not np.isnan(r['auc_ratio']) else "N/A"
                        print(f"    {roi_name}: χ²/ν={chi2_str}, AUC={auc_str}")
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


# Aliases for backward compatibility
process_batch_v6 = process_batch_v9


# =============================================================================
# DEPTH PROFILE GENERATION
# =============================================================================

def create_depth_profile(df: pd.DataFrame,
                          sample_name: str = None,
                          elements_order: List[str] = None,
                          adjustment_method: str = 'none') -> pd.DataFrame:
    """
    Create a depth profile table from uncertainty results.
    """
    if df.empty:
        return pd.DataFrame()
    
    if sample_name:
        df = df[df['filename'] == sample_name]
    
    if df.empty:
        return pd.DataFrame()
    
    layers = sorted(df['layer'].unique())
    
    if elements_order is None:
        elements_order = sorted(df['element'].unique())
    else:
        elements_order = [e for e in elements_order if e in df['element'].unique()]
    
    unc_col_map = {
        'none': 'sigma_total_standard_percent',
        'chi2': 'sigma_total_chi2_adjusted_percent',
        'auc': 'sigma_total_auc_adjusted_percent'
    }
    unc_col = unc_col_map.get(adjustment_method, 'sigma_total_percent')
    if unc_col not in df.columns:
        unc_col = 'sigma_total_percent'
    
    thickness_col = None
    for col_name in ['thickness_1e15_at_cm2', 'thickness']:
        if col_name in df.columns:
            thickness_col = col_name
            break
    
    rows = []
    for layer_num in layers:
        layer_df = df[df['layer'] == layer_num]
        
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
                best_idx = el_df['chi2_reduced'].idxmin()
                best = el_df.loc[best_idx]
                
                row[f'{element} (%)'] = round(best['concentration_percent'], 2)
                row[f'±{element} (%)'] = round(best[unc_col], 3)
            else:
                row[f'{element} (%)'] = None
                row[f'±{element} (%)'] = None
        
        rows.append(row)
    
    return pd.DataFrame(rows)


def create_depth_profiles_excel(df: pd.DataFrame,
                                 output_path: str,
                                 elements_order: List[str] = None,
                                 adjustment_method: str = 'none') -> Dict[str, pd.DataFrame]:
    """
    Create Excel file with depth profiles for all samples.
    """
    profiles = {}
    
    method_names = {
        'none': 'Standard (none)',
        'chi2': 'χ²-adjusted',
        'auc': 'AUC-adjusted'
    }
    
    for sample in df['filename'].unique():
        profile = create_depth_profile(df, sample_name=sample,
                                        elements_order=elements_order,
                                        adjustment_method=adjustment_method)
        profiles[sample] = profile
    
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        summary_rows = []
        for sample, profile in profiles.items():
            sample_df = df[df['filename'] == sample]
            summary_rows.append({
                'Sample': sample,
                'Reaction': sample_df['reaction'].iloc[0] if 'reaction' in sample_df.columns else '',
                'Adjustment': method_names.get(adjustment_method, adjustment_method),
                'Layers': len(profile),
                'Elements': ', '.join(sorted(sample_df['element'].unique())),
            })
        
        pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Summary', index=False)
        
        for sample, profile in profiles.items():
            sheet_name = sample[:31]
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


def list_available_reactions() -> None:
    """Print available reactions and their uncertainties."""
    print("\nAvailable Reactions and Cross-Section Uncertainties:")
    print("=" * 55)
    for reaction, (unc, desc) in CROSS_SECTION_UNCERTAINTIES.items():
        print(f"  {reaction:<30} {unc*100:>5.1f}%  ({desc})")
    print()


def list_isotope_defaults() -> None:
    """Print default isotope cross-section uncertainties."""
    print("\nDefault Isotope Cross-Section Uncertainties:")
    print("=" * 55)
    for isotope, (unc, reaction) in ISOTOPE_CROSS_SECTION_DEFAULTS.items():
        print(f"  {isotope:<6} {unc*100:>5.1f}%  ({reaction})")
    print()


# =============================================================================
# BACKWARD COMPATIBILITY
# =============================================================================

def calculate_uncertainty_for_element(concentration: float,
                                       roi_counts: int,
                                       technique: str = 'RBS',
                                       systematic_uncertainty: float = None) -> Dict[str, float]:
    """Backward-compatible function matching v5 API."""
    xsec_unc = get_cross_section_uncertainty(technique)
    qomega_unc = systematic_uncertainty if systematic_uncertainty else 0.03
    element_counts = int(roi_counts * concentration) if roi_counts > 0 else 0
    
    v9_result = calculate_uncertainty_v9(concentration, element_counts, qomega_unc, xsec_unc)
    
    return {
        'sigma_stat': v9_result['sigma_stat'],
        'sigma_sys': np.sqrt(v9_result['sigma_qomega']**2 + v9_result['sigma_xsec']**2),
        'sigma_total': v9_result['sigma_total'],
    }


DEFAULT_SYSTEMATIC_UNCERTAINTIES = {
    'RBS': 0.03,
    'ERDA': 0.05,
    'NRA': 0.05,
}


# =============================================================================
# MAIN / CLI
# =============================================================================

def main():
    """Command-line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Calculate IBA uncertainties (v9 - Isotope Support)'
    )
    parser.add_argument('input', nargs='+',
                        help='Input .xnra file(s) or directory')
    parser.add_argument('-r', '--roi', action='append', nargs=3,
                        metavar=('START', 'END', 'ELEMENTS'),
                        help='Define ROI: start_ch end_ch "El1,El2,..." (can repeat)')
    parser.add_argument('-q', '--qomega', type=float, default=0.03,
                        help='Q·Ω uncertainty as fraction (default: 0.03 = 3%%)')
    parser.add_argument('--reaction', default='RBS: Mx(α,α)Mx',
                        help='Default reaction type for cross-section lookup')
    parser.add_argument('--xsec', type=float,
                        help='Override cross-section uncertainty (fraction)')
    parser.add_argument('--adjustment', choices=['none', 'chi2', 'auc'], default='none',
                        help='Adjustment method')
    parser.add_argument('--isotope-mode', action='store_true',
                        help='Enable isotope-specific calculations')
    parser.add_argument('--element-spectrum',
                        help='Path to element spectrum file for accurate A_A')
    parser.add_argument('-o', '--output',
                        help='Output file (CSV or Excel based on extension)')
    parser.add_argument('--profiles',
                        help='Output Excel for depth profiles')
    parser.add_argument('--list-reactions', action='store_true',
                        help='List available reactions and exit')
    parser.add_argument('--list-isotopes', action='store_true',
                        help='List default isotope uncertainties and exit')
    
    args = parser.parse_args()
    
    if args.list_reactions:
        list_available_reactions()
        return
    
    if args.list_isotopes:
        list_isotope_defaults()
        return
    
    if not args.roi:
        print("Error: At least one ROI must be defined with -r START END ELEMENTS")
        return
    
    regions_of_interest = []
    for start, end, elements in args.roi:
        regions_of_interest.append({
            'name': f"ch{start}-{end}",
            'channels': (int(start), int(end)),
            'elements': [e.strip() for e in elements.split(',')]
        })
    
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
    
    # Load element spectrum if provided
    element_spectrum_data = None
    if args.element_spectrum:
        try:
            element_spectrum_data = parse_element_spectrum(args.element_spectrum)
            print(f"✓ Loaded element spectrum: {args.element_spectrum}")
        except Exception as e:
            print(f"⚠ Failed to load element spectrum: {e}")
    
    print(f"Processing {len(filepaths)} files...")
    
    output_csv = None
    output_excel = None
    if args.output:
        if args.output.endswith('.xlsx'):
            output_excel = args.output
        else:
            output_csv = args.output if args.output.endswith('.csv') else args.output + '.csv'
    
    df = process_batch_v9(
        filepaths,
        regions_of_interest,
        qomega_uncertainty=args.qomega,
        reaction=args.reaction,
        cross_section_override=args.xsec,
        adjustment_method=args.adjustment,
        element_spectrum_data=element_spectrum_data,
        isotope_mode=args.isotope_mode,
        output_csv=output_csv,
        output_excel=output_excel
    )
    
    print(f"\nProcessed {len(df)} element-layer-ROI combinations")
    
    if args.profiles and not df.empty:
        create_depth_profiles_excel(df, args.profiles, adjustment_method=args.adjustment)


if __name__ == '__main__':
    main()
