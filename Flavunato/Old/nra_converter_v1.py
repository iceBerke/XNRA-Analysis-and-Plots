"""
XNRA Spectrum Viewer
A simple viewer for SIMNRA .xnra files with matplotlib popup display

HOW TO USE:
1. Change the FILE_PATH below to your .xnra file location
2. Adjust visual parameters in the CONFIGURATION section if desired
3. Run this script in VSCode (F5 or right-click > Run Python File)
4. The plot will appear in a popup window

Author: Your Name
Date: 2025-12-11
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================================
# CONFIGURATION - CHANGE THESE VALUES TO CUSTOMIZE THE PLOT
# ============================================================================

# FILE PATH - Change this to your .xnra file location
FILE_PATH = r"C:\Users\berke\Desktop\Giuseppe\input\NRA Re 550 RIC 166°.xnra"
#FILE_PATH = r"C:\Users\berke\Desktop\Giuseppe\input\NRA\NRA_250_20mg_NA_150°.xnra"


# PLOT APPEARANCE
FIGURE_SIZE = (14, 8)           # Width, height in inches
DPI = 100                       # Resolution (higher = sharper but slower)

# COLORS
RAW_DATA_COLOR = 'blue'         # Color for raw data points
SMOOTHED_DATA_COLOR = 'red'     # Color for smoothed line
GRID_COLOR = 'gray'             # Grid line color
BACKGROUND_COLOR = 'white'      # Plot background

# LINE STYLES
RAW_MARKER_SIZE = 1.5           # Size of raw data points
RAW_ALPHA = 0.4                 # Transparency of raw data (0=invisible, 1=solid)
SMOOTHED_LINE_WIDTH = 2.0       # Width of smoothed line
SMOOTHED_ALPHA = 1.0            # Transparency of smoothed line

# TEXT SIZES
TITLE_SIZE = 16                 # Main title font size
AXIS_LABEL_SIZE = 13            # X and Y axis label size
TICK_LABEL_SIZE = 11            # Numbers on axes
LEGEND_SIZE = 11                # Legend text size
INFO_BOX_SIZE = 10              # Info box text size

# PLOT LIMITS
ENERGY_MIN = 0                  # Minimum energy to display (keV)
ENERGY_MAX = 6500               # Maximum energy to display (keV)
AUTO_Y_SCALE = True             # Automatically scale Y axis?

# DISPLAY OPTIONS
SHOW_RAW_DATA = True            # Show raw data points?
SHOW_SMOOTHED_DATA = True       # Show smoothed line?
SHOW_GRID = True                # Show grid lines?
SHOW_LEGEND = True              # Show legend?
SHOW_INFO_BOX = True            # Show info box with parameters?

# ============================================================================
# END CONFIGURATION - Don't change anything below unless you know what you're doing
# ============================================================================

def parse_xnra_file(filepath):
    """
    Parse XNRA file and extract spectrum data and metadata
    
    Args:
        filepath: Path to .xnra file
        
    Returns:
        dict: Dictionary containing all extracted data
    """
    print(f"Reading file: {filepath}")
    
    # Parse XML
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    # Define namespaces
    ns = {
        'idf': 'http://idf.schemas.itn.pt',
        'simnra': 'http://www.simnra.com/simnra'
    }
    
    # Initialize data dictionary
    data = {
        'filename': None,
        'raw_channels': None,
        'raw_counts': None,
        'smoothed_channels': None,
        'smoothed_counts': None,
        'beam_particle': None,
        'beam_energy': None,
        'scattering_angle': None,
        'calibration_offset': None,
        'calibration_gain': None,
        'detector_type': None,
        'resolution': None
    }
    
    # Extract filename
    filename_node = root.find('.//idf:filename', ns)
    if filename_node is not None:
        data['filename'] = filename_node.text
    
    # Extract beam parameters
    beam_particle_node = root.find('.//idf:beamparticle', ns)
    if beam_particle_node is not None:
        data['beam_particle'] = beam_particle_node.text
    
    beam_energy_node = root.find('.//idf:beamenergy', ns)
    if beam_energy_node is not None:
        data['beam_energy'] = float(beam_energy_node.text)
    
    # Extract geometry
    scattering_angle_node = root.find('.//idf:scatteringangle', ns)
    if scattering_angle_node is not None:
        data['scattering_angle'] = float(scattering_angle_node.text)
    
    # Extract calibration
    calibration_params = root.findall('.//idf:energycalibration[1]/idf:calibrationparameters/idf:calibrationparameter', ns)
    if len(calibration_params) >= 2:
        data['calibration_offset'] = float(calibration_params[0].text)
        data['calibration_gain'] = float(calibration_params[1].text)
    
    # Extract detector info
    detector_node = root.find('.//idf:detectortype', ns)
    if detector_node is not None:
        data['detector_type'] = detector_node.text
    
    resolution_node = root.find('.//idf:detectorresolution[1]/idf:resolutionparameters/idf:resolutionparameter', ns)
    if resolution_node is not None:
        data['resolution'] = float(resolution_node.text)
    
    # Extract raw data
    raw_data = root.find('.//idf:data/idf:simpledata', ns)
    if raw_data is not None:
        x_text = raw_data.find('idf:x', ns).text
        y_text = raw_data.find('idf:y', ns).text
        data['raw_channels'] = np.array([float(x) for x in x_text.split()])
        data['raw_counts'] = np.array([float(y) for y in y_text.split()])
        print(f"  ✓ Raw data: {len(data['raw_channels'])} channels")
    
    # Extract smoothed data
    smoothed_data = root.find('.//simnra:smootheddata/idf:simpledata', ns)
    if smoothed_data is not None:
        x_text = smoothed_data.find('idf:x', ns).text
        y_text = smoothed_data.find('idf:y', ns).text
        data['smoothed_channels'] = np.array([float(x) for x in x_text.split()])
        data['smoothed_counts'] = np.array([float(y) for y in y_text.split()])
        print(f"  ✓ Smoothed data: {len(data['smoothed_channels'])} channels")
    
    return data


def plot_spectrum(data):
    """
    Create matplotlib plot of spectrum data
    
    Args:
        data: Dictionary containing spectrum data and metadata
    """
    print("\nCreating plot...")
    
    # Calculate energies
    if data['calibration_offset'] is not None and data['calibration_gain'] is not None:
        raw_energy = data['calibration_offset'] + data['calibration_gain'] * data['raw_channels']
        if data['smoothed_channels'] is not None:
            smoothed_energy = data['calibration_offset'] + data['calibration_gain'] * data['smoothed_channels']
        use_energy = True
        x_label = 'Energy (keV)'
        print(f"  ✓ Using energy calibration: E = {data['calibration_offset']} + {data['calibration_gain']} × Channel")
    else:
        raw_energy = data['raw_channels']
        if data['smoothed_channels'] is not None:
            smoothed_energy = data['smoothed_channels']
        use_energy = False
        x_label = 'Channel'
        print("  ⚠ No calibration found, using channel numbers")
    
    # Create figure
    fig, ax = plt.subplots(figsize=FIGURE_SIZE, dpi=DPI)
    fig.patch.set_facecolor(BACKGROUND_COLOR)
    ax.set_facecolor(BACKGROUND_COLOR)
    
    # Plot raw data
    if SHOW_RAW_DATA and data['raw_counts'] is not None:
        ax.plot(raw_energy, data['raw_counts'], 
                'o', 
                color=RAW_DATA_COLOR,
                markersize=RAW_MARKER_SIZE,
                alpha=RAW_ALPHA,
                label='Raw data')
        print("  ✓ Plotted raw data")
    
    # Plot smoothed data
    if SHOW_SMOOTHED_DATA and data['smoothed_counts'] is not None:
        ax.plot(smoothed_energy, data['smoothed_counts'],
                '-',
                color=SMOOTHED_DATA_COLOR,
                linewidth=SMOOTHED_LINE_WIDTH,
                alpha=SMOOTHED_ALPHA,
                label='Smoothed data')
        print("  ✓ Plotted smoothed data")
    
    # Set labels and title
    ax.set_xlabel(x_label, fontsize=AXIS_LABEL_SIZE, fontweight='bold')
    ax.set_ylabel('Yield (counts)', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
    
    # Create title
    title = 'XNRA Spectrum'
    if data['filename']:
        title = data['filename']
    ax.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=20)
    
    # Set tick label sizes
    ax.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)
    
    # Set x-axis limits
    if use_energy:
        ax.set_xlim(ENERGY_MIN, ENERGY_MAX)
    
    # Add grid
    if SHOW_GRID:
        ax.grid(True, alpha=0.3, color=GRID_COLOR, linestyle='--', linewidth=0.5)
    
    # Add legend
    if SHOW_LEGEND:
        ax.legend(loc='upper right', fontsize=LEGEND_SIZE, framealpha=0.9)
    
    # Add info box
    if SHOW_INFO_BOX:
        info_lines = []
        if data['beam_particle'] and data['beam_energy']:
            info_lines.append(f"Beam: {data['beam_particle']}, {data['beam_energy']:.0f} keV")
        if data['scattering_angle']:
            info_lines.append(f"Angle: {data['scattering_angle']:.1f}°")
        if data['detector_type']:
            info_lines.append(f"Detector: {data['detector_type']}")
        if data['resolution']:
            info_lines.append(f"Resolution: {data['resolution']:.1f} keV FWHM")
        
        if info_lines:
            info_text = '\n'.join(info_lines)
            ax.text(0.02, 0.98, info_text,
                   transform=ax.transAxes,
                   fontsize=INFO_BOX_SIZE,
                   verticalalignment='top',
                   bbox=dict(boxstyle='round',
                           facecolor='wheat',
                           alpha=0.8,
                           edgecolor='black',
                           linewidth=1))
    
    # Tight layout
    plt.tight_layout()
    
    print("  ✓ Plot created successfully!")
    print("\n" + "="*60)
    print("PLOT READY - Close the window to exit")
    print("="*60)
    
    # Show plot (this will block until window is closed)
    plt.show()


def main():
    """Main function to run the XNRA viewer"""
    
    print("="*60)
    print("XNRA SPECTRUM VIEWER")
    print("="*60)
    
    # Check if file exists
    filepath = Path(FILE_PATH)
    if not filepath.exists():
        print(f"\n❌ ERROR: File not found!")
        print(f"   Looked for: {filepath.absolute()}")
        print(f"\n💡 TIP: Update the FILE_PATH variable at the top of this script")
        print(f"   Current value: {FILE_PATH}")
        return
    
    try:
        # Parse file
        data = parse_xnra_file(filepath)
        
        # Check if we have data
        if data['raw_counts'] is None and data['smoothed_counts'] is None:
            print("\n❌ ERROR: No spectrum data found in file!")
            return
        
        # Plot
        plot_spectrum(data)
        
    except ET.ParseError as e:
        print(f"\n❌ ERROR: Invalid XML file!")
        print(f"   {e}")
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()