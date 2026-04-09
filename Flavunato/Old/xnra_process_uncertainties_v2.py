"""
XNRA Batch Uncertainty Processor v2
===================================
Quick script to process all your .xnra files and generate uncertainty reports.

CHANGELOG v2:
- Uses consistent fit range (handles partial spectrum fitting correctly!)
- Shows fit range info in output
- New option to manually override fit range

Usage:
    1. Put all your .xnra files in a folder (or subfolders)
    2. Edit the SETTINGS section below
    3. Run: python xnra_process_uncertainties_v2.py

Output:
    - Depth profile Excel (one sheet per sample) - THIS IS WHAT YOU WANT
    - Full results CSV/Excel with all layers (raw data)

Authors: Berke Santos & Giuseppe Legrottaglie
"""

from pathlib import Path
from xnra_uncertainty_v2 import (
    process_batch, 
    find_xnra_files, 
    create_depth_profile_table,
    create_all_depth_profiles
)

# =============================================================================
# SETTINGS - EDIT THESE
# =============================================================================

# Folder containing your .xnra files (can include subfolders)
INPUT_FOLDER = r"C:\Users\berke.santos\Documents\Giu\MultiSimnra-test\xnra"  # Current directory, or specify path like r"C:\Data\IBA"

# Elements you want to analyze (set to None for ALL elements)
ELEMENTS_OF_INTEREST = ['H', 'C', 'O', 'Na', 'Si', 'Ca', 'Ba']

# Output filenames
OUTPUT_DEPTH_PROFILES = "depth_profiles.xlsx"      # One sheet per sample - MAIN OUTPUT
OUTPUT_FULL = "uncertainty_results_full.xlsx"      # All raw data

# =============================================================================
# FIT RANGE SETTINGS (NEW in v2!)
# =============================================================================

# Set to None for AUTO-DETECTION (recommended - detects from simulation extent)
# Or specify (start_channel, end_channel) to force a specific range for ALL files
# Example: FIT_RANGE = (0, 1000)  # Only use channels 0-1000
FIT_RANGE = None  # Auto-detect from simulation

# =============================================================================
# PROCESSING - DON'T EDIT BELOW
# =============================================================================

def main():
    print("="*60)
    print("XNRA UNCERTAINTY BATCH PROCESSOR v2")
    print("="*60)
    
    # Find all .xnra files
    input_path = Path(INPUT_FOLDER)
    
    # Search recursively
    xnra_files = []
    for pattern in ['*.xnra', '**/*.xnra']:
        xnra_files.extend(input_path.glob(pattern))
    
    # Remove duplicates and sort
    xnra_files = sorted(set(str(f) for f in xnra_files))
    
    if not xnra_files:
        print(f"No .xnra files found in: {input_path.absolute()}")
        print("Please check the INPUT_FOLDER setting.")
        return
    
    print(f"\nFound {len(xnra_files)} .xnra files")
    print(f"Elements to analyze: {ELEMENTS_OF_INTEREST or 'ALL'}")
    
    if FIT_RANGE:
        print(f"Fit range: channels {FIT_RANGE[0]}-{FIT_RANGE[1]} (manual override)")
    else:
        print(f"Fit range: AUTO-DETECT from simulation extent")
    print()
    
    # Process all files
    df = process_batch(
        xnra_files,
        elements_of_interest=ELEMENTS_OF_INTEREST,
        output_excel=OUTPUT_FULL,
        force_fit_range=FIT_RANGE
    )
    
    if df.empty:
        print("No results generated. Check your files.")
        return
    
    # Create depth profile tables (ONE SHEET PER SAMPLE)
    print("\nCreating depth profile tables...")
    profiles = create_all_depth_profiles(
        df, 
        elements=ELEMENTS_OF_INTEREST,
        output_excel=OUTPUT_DEPTH_PROFILES
    )
    
    # Print quick stats
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    print(f"Files processed: {len(xnra_files)}")
    print(f"Samples with depth profiles: {len(profiles)}")
    
    # Show fit range stats
    if 'fit_fraction' in df.columns:
        avg_fit_frac = df.groupby('filename')['fit_fraction'].first().mean() * 100
        print(f"Average spectrum coverage: {avg_fit_frac:.1f}%")
    
    # Show first sample's profile as preview
    first_sample = list(profiles.keys())[0]
    print(f"\nPreview of '{first_sample}' depth profile:")
    print("-"*60)
    preview = profiles[first_sample][['Layer', 'Thickness (10¹⁵ at/cm²)']].copy()
    for el in ELEMENTS_OF_INTEREST[:4]:  # First 4 elements
        col_name = f'{el} (%)'
        if col_name in profiles[first_sample].columns:
            preview[col_name] = profiles[first_sample][col_name]
    print(preview.head(5).to_string(index=False))
    print("...")
    
    print(f"\n✓ Done! Your depth profiles are in: {OUTPUT_DEPTH_PROFILES}")
    print(f"  (Full data with fit range info in: {OUTPUT_FULL})")


if __name__ == '__main__':
    main()
