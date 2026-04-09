"""
XNRA Multi-Technique Uncertainty Processor v4
==============================================
Combine RBS, ERDA, and NRA spectra for complete elemental analysis.

Typical IBA workflow:
  - RBS:  Heavy/medium elements (Ca, Ba, Si, Na, Fe, Al, Mg, etc.)
  - ERDA: Hydrogen (H) - and sometimes D, He
  - NRA:  Light elements via nuclear reactions (C, O, N)

Each technique has:
  - Different systematic uncertainties
  - Different optimal channel ranges for different elements
  - Its own chi² evaluation

This script combines them into unified depth profiles.

Usage:
    1. Edit the ANALYSIS_CONFIG below with your files and ROIs
    2. Run: python xnra_process_uncertainties_v4.py

Authors: Berke Santos & Giuseppe Legrottaglie
"""

from pathlib import Path
from xnra_uncertainty_v4 import (
    process_multi_technique,
    create_all_simple_profiles,
    create_all_combined_profiles
)

# =============================================================================
# SETTINGS
# =============================================================================

# Output files
OUTPUT_FULL = "uncertainty_results_v4.xlsx"           # All raw data
OUTPUT_PROFILES = "depth_profiles_v4.xlsx"            # Clean depth profiles
OUTPUT_DETAILED = "depth_profiles_detailed_v4.xlsx"   # Profiles with technique info

# Preferred element order in output (left to right)
ELEMENTS_ORDER = ['H', 'C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'S', 'K', 'Ca', 'Fe', 'Ba']

# =============================================================================
# ANALYSIS CONFIGURATION - THIS IS THE KEY PART!
# =============================================================================
#
# Define your analysis here. Each entry has:
#   - 'file': path to the .xnra file
#   - 'technique': 'RBS', 'ERDA', or 'NRA'
#   - 'sample': sample name (use SAME name for files from same sample!)
#   - 'rois': list of ROI definitions for this technique
#
# EXAMPLE for a sample measured with all three techniques:
#
# ANALYSIS_CONFIG = [
#     # RBS spectrum - heavy and medium elements
#     {
#         'file': r"C:\Data\Sample1_RBS.xnra",
#         'technique': 'RBS',
#         'sample': 'Sample1',
#         'rois': [
#             {'name': 'Heavy', 'channels': (700, 900), 'elements': ['Ba', 'Ca', 'Fe']},
#             {'name': 'Medium', 'channels': (400, 550), 'elements': ['Si', 'Al', 'Na', 'Mg']},
#         ]
#     },
#     # ERDA spectrum - hydrogen
#     {
#         'file': r"C:\Data\Sample1_ERDA.xnra",
#         'technique': 'ERDA',
#         'sample': 'Sample1',
#         'rois': [
#             {'name': 'Hydrogen', 'channels': (200, 400), 'elements': ['H']},
#         ]
#     },
#     # NRA spectrum - light elements (C, O via nuclear reactions)
#     {
#         'file': r"C:\Data\Sample1_NRA.xnra",
#         'technique': 'NRA',
#         'sample': 'Sample1',
#         'rois': [
#             {'name': 'Carbon', 'channels': (100, 200), 'elements': ['C']},
#             {'name': 'Oxygen', 'channels': (300, 500), 'elements': ['O']},
#         ]
#     },
# ]

# YOUR CONFIGURATION HERE:
ANALYSIS_CONFIG = [
    # Example with just RBS (replace with your files)
    {
        'file': r"RBS_300.xnra",
        'technique': 'RBS',
        'sample': 'Sample300',
        'rois': [
            {'name': 'Heavy', 'channels': (700, 900), 'elements': ['Ba', 'Ca', 'Fe']},
            {'name': 'Medium', 'channels': (380, 550), 'elements': ['Si', 'Al', 'Na', 'Mg']},
        ]
    },
    # Add ERDA file for H
    # {
    #     'file': r"ERDA_300.xnra",
    #     'technique': 'ERDA',
    #     'sample': 'Sample300',  # Same sample name!
    #     'rois': [
    #         {'name': 'Hydrogen', 'channels': (200, 400), 'elements': ['H']},
    #     ]
    # },
    # Add NRA file for C, O
    # {
    #     'file': r"NRA_300.xnra",
    #     'technique': 'NRA',
    #     'sample': 'Sample300',  # Same sample name!
    #     'rois': [
    #         {'name': 'Carbon', 'channels': (100, 200), 'elements': ['C']},
    #         {'name': 'Oxygen', 'channels': (300, 500), 'elements': ['O']},
    #     ]
    # },
]

# =============================================================================
# PROCESSING
# =============================================================================

def main():
    print("="*70)
    print("XNRA MULTI-TECHNIQUE UNCERTAINTY PROCESSOR v4")
    print("="*70)
    
    # Check files exist
    missing = []
    for config in ANALYSIS_CONFIG:
        if not Path(config['file']).exists():
            missing.append(config['file'])
    
    if missing:
        print("\n⚠ WARNING: Some files not found:")
        for f in missing:
            print(f"  - {f}")
        print("\nPlease check the file paths in ANALYSIS_CONFIG.")
        # Continue anyway in case of relative paths
    
    # Process all techniques
    df = process_multi_technique(
        ANALYSIS_CONFIG,
        output_excel=OUTPUT_FULL
    )
    
    if df.empty:
        print("\nNo results generated. Check your configuration.")
        return
    
    # Create depth profiles
    print("\n" + "="*70)
    print("CREATING DEPTH PROFILES")
    print("="*70)
    
    # Simple profiles (just concentrations and uncertainties)
    simple_profiles = create_all_simple_profiles(
        df,
        elements_order=ELEMENTS_ORDER,
        output_excel=OUTPUT_PROFILES
    )
    
    # Detailed profiles (includes which technique was used)
    detailed_profiles = create_all_combined_profiles(
        df,
        output_excel=OUTPUT_DETAILED
    )
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    samples = df['sample'].unique()
    print(f"\nSamples analyzed: {len(samples)}")
    
    for sample in samples:
        sample_df = df[df['sample'] == sample]
        techniques = sample_df['technique'].unique()
        elements = sorted(sample_df['element'].unique())
        
        print(f"\n{sample}:")
        print(f"  Techniques: {', '.join(techniques)}")
        print(f"  Elements: {', '.join(elements)}")
        
        # Chi² summary per technique
        for tech in techniques:
            tech_df = sample_df[sample_df['technique'] == tech]
            avg_chi2 = tech_df.groupby('roi_name')['chi2_reduced'].first().mean()
            print(f"  {tech} avg χ²/ν: {avg_chi2:.2f}")
    
    # Preview first sample
    if simple_profiles:
        first_sample = list(simple_profiles.keys())[0]
        print(f"\nPreview of '{first_sample}' depth profile:")
        print("-"*70)
        preview = simple_profiles[first_sample]
        # Show first few columns
        cols_to_show = list(preview.columns)[:8]
        print(preview[cols_to_show].head(5).to_string(index=False))
        print("...")
    
    print(f"\n✓ Done!")
    print(f"  Full results: {OUTPUT_FULL}")
    print(f"  Depth profiles: {OUTPUT_PROFILES}")
    print(f"  Detailed profiles: {OUTPUT_DETAILED}")


if __name__ == '__main__':
    main()
