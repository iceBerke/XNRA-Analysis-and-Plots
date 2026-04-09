"""
XNRA Batch Uncertainty Processor v3
===================================
Process .xnra files using ROI-based analysis.

Each ROI defines:
  - A channel range to analyze
  - Which elements have peaks in that range

This is physically correct: chi² and uncertainties are calculated
separately for each ROI, using only the counts in that region.

Usage:
    1. Put all your .xnra files in a folder
    2. Edit the SETTINGS section below (especially REGIONS_OF_INTEREST)
    3. Run: python xnra_process_uncertainties_v3.py

Authors: Berke Santos & Giuseppe Legrottaglie
"""

from pathlib import Path
from xnra_uncertainty_v3 import (
    process_batch_with_rois,
    find_xnra_files,
    create_all_depth_profiles
)

# =============================================================================
# SETTINGS - EDIT THESE
# =============================================================================

# Folder containing your .xnra files (can include subfolders)
INPUT_FOLDER = r"C:\Users\berke.santos\Documents\Giu\MultiSimnra-test\xnra"  # Current directory, or specify path

# Output filenames
OUTPUT_DEPTH_PROFILES = "depth_profiles_v3.xlsx"    # One sheet per sample/ROI
OUTPUT_FULL = "uncertainty_results_v3.xlsx"          # All raw data

# =============================================================================
# REGIONS OF INTEREST - THIS IS THE KEY PART!
# =============================================================================
# 
# Define your ROIs here. Each ROI has:
#   - 'name': A descriptive name (optional, will auto-generate if missing)
#   - 'channels': (start, end) - the channel range to analyze
#   - 'elements': List of elements whose peaks fall in this range
#
# The script will:
#   1. Calculate chi² separately for each ROI
#   2. Use counts from that ROI to calculate uncertainties
#   3. Only report uncertainties for elements assigned to that ROI
#
# TIPS:
#   - Look at your spectrum in SIMNRA to identify channel ranges
#   - Heavy elements (Ba, Ca) → higher channels
#   - Light elements (O, C, H) → lower channels
#   - ROIs can overlap if needed

REGIONS_OF_INTEREST = [
    {
        'name': 'Heavy',
        'channels': (700, 900),
        'elements': ['Ba', 'Ca', 'Fe']
    },
    {
        'name': 'Medium', 
        'channels': (400, 550),
        'elements': ['Si', 'Al', 'Na', 'Mg']
    }
]

# =============================================================================
# PROCESSING - DON'T EDIT BELOW
# =============================================================================

def main():
    print("="*60)
    print("XNRA UNCERTAINTY BATCH PROCESSOR v3 (ROI-based)")
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
    
    # Process all files
    df = process_batch_with_rois(
        xnra_files,
        regions_of_interest=REGIONS_OF_INTEREST,
        output_excel=OUTPUT_FULL
    )
    
    if df.empty:
        print("No results generated. Check your files and ROI definitions.")
        return
    
    # Create depth profile tables
    print("\nCreating depth profile tables...")
    profiles = create_all_depth_profiles(df, output_excel=OUTPUT_DEPTH_PROFILES)
    
    # Print summary
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    print(f"Files processed: {len(xnra_files)}")
    print(f"ROIs defined: {len(REGIONS_OF_INTEREST)}")
    print(f"Depth profiles created: {len(profiles)}")
    
    # Chi² summary per ROI
    print("\nChi² summary by ROI:")
    for roi in REGIONS_OF_INTEREST:
        roi_name = roi.get('name', f"ch{roi['channels'][0]}-{roi['channels'][1]}")
        roi_df = df[df['roi_name'] == roi_name]
        if not roi_df.empty:
            avg_chi2 = roi_df.groupby('filename')['chi2_reduced'].first().mean()
            print(f"  {roi_name}: avg χ²/ν = {avg_chi2:.2f}")
    
    # Preview first profile
    if profiles:
        first_key = list(profiles.keys())[0]
        print(f"\nPreview of '{first_key}':")
        print("-"*60)
        print(profiles[first_key].head(5).to_string(index=False))
        print("...")
    
    print(f"\n✓ Done!")
    print(f"  Depth profiles: {OUTPUT_DEPTH_PROFILES}")
    print(f"  Full data: {OUTPUT_FULL}")


if __name__ == '__main__':
    main()
