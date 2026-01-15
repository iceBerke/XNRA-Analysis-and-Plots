"""
XNRA Spectrum Viewer - SIMNRA Style with Dual X-Axes
Displays Channel on bottom axis and Energy on top axis, just like SIMNRA

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

# PLOT APPEARANCE
FIGURE_SIZE = (14, 6)           # Width, height in inches
DPI = 100                       # Resolution (higher = sharper but slower)

# COLORS
RAW_DATA_COLOR = 'blue'         # Color for raw data points
SMOOTHED_DATA_COLOR = 'red'     # Color for smoothed line
GRID_COLOR = 'gray'             # Grid line color
BACKGROUND_COLOR = 'white'      # Plot background

# LINE STYLES
RAW_MARKER_SIZE = 1.5           # Size of raw data points
RAW_ALPHA = 0.4                 # Transparency of raw data (0=invisible, 1=solid)
SMOOTHED_LINE_WIDTH = 1.0       # Width of smoothed line
SMOOTHED_ALPHA = 1.0            # Transparency of smoothed line

# TEXT SIZES
TITLE_SIZE = 16                 # Main title font size
AXIS_LABEL_SIZE = 13            # X and Y axis label size
TICK_LABEL_SIZE = 11            # Numbers on axes
LEGEND_SIZE = 11                # Legend text size
INFO_BOX_SIZE = 10              # Info box text size

# PLOT LIMITS
CHANNEL_MIN = 0                 # Minimum channel to display
CHANNEL_MAX = 1600              # Maximum channel to display (auto if None)
AUTO_Y_SCALE = True             # Automatically scale Y axis?

# DISPLAY OPTIONS
SHOW_RAW_DATA = True            # Show raw data points?
SHOW_SMOOTHED_DATA = True       # Show smoothed line?
SHOW_GRID = True                # Show grid lines?
SHOW_LEGEND = True              # Show legend?
SHOW_INFO_BOX = True            # Show info box with parameters?

# AXIS STYLE (SIMNRA-like)
DUAL_X_AXES = True              # Show both Channel and Energy axes?
ENERGY_ON_TOP = True            # Energy on top, Channel on bottom (SIMNRA style)

# ============================================================================
# END CONFIGURATION
# ============================================================================

def parse_xnra_file(filepath):
    """Parse XNRA file and extract spectrum data and metadata"""
    print(f"Reading file: {filepath}")
    
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    ns = {
        'idf': 'http://idf.schemas.itn.pt',
        'simnra': 'http://www.simnra.com/simnra'
    }
    
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
    
    # Extract metadata
    filename_node = root.find('.//idf:filename', ns)
    if filename_node is not None:
        data['filename'] = filename_node.text
    
    beam_particle_node = root.find('.//idf:beamparticle', ns)
    if beam_particle_node is not None:
        data['beam_particle'] = beam_particle_node.text
    
    beam_energy_node = root.find('.//idf:beamenergy', ns)
    if beam_energy_node is not None:
        data['beam_energy'] = float(beam_energy_node.text)
    
    scattering_angle_node = root.find('.//idf:scatteringangle', ns)
    if scattering_angle_node is not None:
        data['scattering_angle'] = float(scattering_angle_node.text)
    
    # Extract calibration
    calibration_params = root.findall('.//idf:energycalibration[1]/idf:calibrationparameters/idf:calibrationparameter', ns)
    if len(calibration_params) >= 2:
        data['calibration_offset'] = float(calibration_params[0].text)
        data['calibration_gain'] = float(calibration_params[1].text)
    
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


def channel_to_energy(channel, offset, gain):
    """Convert channel number to energy"""
    return offset + gain * channel


def energy_to_channel(energy, offset, gain):
    """Convert energy to channel number"""
    return (energy - offset) / gain


def plot_spectrum_dual_axes(data):
    """Create matplotlib plot with dual x-axes (Channel + Energy)"""
    print("\nCreating plot with dual x-axes...")
    
    # Calculate energies
    cal_offset = data['calibration_offset'] if data['calibration_offset'] is not None else 0
    cal_gain = data['calibration_gain'] if data['calibration_gain'] is not None else 1
    
    raw_channels = data['raw_channels']
    raw_counts = data['raw_counts']
    raw_energy = channel_to_energy(raw_channels, cal_offset, cal_gain)
    
    if data['smoothed_channels'] is not None:
        smoothed_channels = data['smoothed_channels']
        smoothed_counts = data['smoothed_counts']
        smoothed_energy = channel_to_energy(smoothed_channels, cal_offset, cal_gain)
    
    # Create figure
    fig, ax_bottom = plt.subplots(figsize=FIGURE_SIZE, dpi=DPI)
    fig.patch.set_facecolor(BACKGROUND_COLOR)
    ax_bottom.set_facecolor(BACKGROUND_COLOR)
    
    # Plot data using CHANNELS on the primary (bottom) axis
    if SHOW_RAW_DATA and raw_counts is not None:
        ax_bottom.plot(raw_channels, raw_counts, 
                       'o', 
                       color=RAW_DATA_COLOR,
                       markersize=RAW_MARKER_SIZE,
                       alpha=RAW_ALPHA,
                       label='Raw data')
        print("  ✓ Plotted raw data")
    
    if SHOW_SMOOTHED_DATA and data['smoothed_counts'] is not None:
        ax_bottom.plot(smoothed_channels, smoothed_counts,
                       '-',
                       color=SMOOTHED_DATA_COLOR,
                       linewidth=SMOOTHED_LINE_WIDTH,
                       alpha=SMOOTHED_ALPHA,
                       label='Smoothed data')
        print("  ✓ Plotted smoothed data")
    
    # Set bottom axis (Channel)
    ax_bottom.set_xlabel('Channel', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
    ax_bottom.set_ylabel('Counts', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
    ax_bottom.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)
    
    # Set channel limits
    if CHANNEL_MAX is not None:
        ax_bottom.set_xlim(CHANNEL_MIN, CHANNEL_MAX)
    else:
        ax_bottom.set_xlim(CHANNEL_MIN, raw_channels.max())
    
    # Create top axis (Energy) using secondary_xaxis
    if DUAL_X_AXES:
        ax_top = ax_bottom.secondary_xaxis('top', 
                                           functions=(
                                               lambda ch: channel_to_energy(ch, cal_offset, cal_gain),
                                               lambda en: energy_to_channel(en, cal_offset, cal_gain)
                                           ))
        ax_top.set_xlabel('Energy [keV]', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
        ax_top.tick_params(axis='x', labelsize=TICK_LABEL_SIZE)
        print("  ✓ Added dual x-axes (Channel + Energy)")
    
    # Title
    title = 'XNRA Spectrum'
    if data['filename']:
        title = data['filename']
    ax_bottom.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=40)  # More padding for top axis
    
    # Grid
    if SHOW_GRID:
        ax_bottom.grid(True, alpha=0.3, color=GRID_COLOR, linestyle='--', linewidth=0.5)
    
    # Legend
    if SHOW_LEGEND:
        ax_bottom.legend(loc='upper right', fontsize=LEGEND_SIZE, framealpha=0.9)
    
    # Info box
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
        if cal_offset is not None and cal_gain is not None:
            info_lines.append(f"Calibration: E = {cal_offset:.1f} + {cal_gain:.3f}×Ch")
        
        if info_lines:
            info_text = '\n'.join(info_lines)
            ax_bottom.text(0.02, 0.98, info_text,
                          transform=ax_bottom.transAxes,
                          fontsize=INFO_BOX_SIZE,
                          verticalalignment='top',
                          bbox=dict(boxstyle='round',
                                  facecolor='wheat',
                                  alpha=0.8,
                                  edgecolor='black',
                                  linewidth=1))
    
    plt.tight_layout()
    
    print("  ✓ Plot created successfully!")
    print("\n" + "="*60)
    print("PLOT READY - Close the window to exit")
    print("="*60)
    
    plt.show()


def main():
    """Main function to run the XNRA viewer"""
    
    print("="*60)
    print("XNRA SPECTRUM VIEWER - SIMNRA Style (Dual X-Axes)")
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
        
        # Plot with dual axes
        plot_spectrum_dual_axes(data)
        
    except ET.ParseError as e:
        print(f"\n❌ ERROR: Invalid XML file!")
        print(f"   {e}")
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()