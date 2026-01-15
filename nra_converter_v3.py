"""
XNRA Spectrum Viewer - SIMNRA Style with Dual X-Axes
Enhanced version with customizable markers and line connections

Features:
- Dual x-axes: Channel (bottom) and Energy (top)
- Customizable raw data markers (shape, fill, size)
- Customizable smoothed data markers (shape, fill, size)
- Optional line connecting data points
- Marker range control for both raw and smoothed data
- Full control over colors, sizes, and appearance

HOW TO USE:
1. Change FILE_PATH to your .xnra file location
2. Customize appearance in CONFIGURATION section
3. Run in VSCode (F5)
4. Plot appears in popup window

Author: Your Name
Date: 2025-12-12
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================================
# CONFIGURATION - CUSTOMIZE YOUR PLOT HERE
# ============================================================================

# ==================== FILE PATH ====================
#FILE_PATH = r"C:\Users\berke\Desktop\Giuseppe\input\NRA\NRA Re 550 RIC 166°.xnra"
FILE_PATH = r"C:\Users\berke\Desktop\Giuseppe\input\RBS\RBS_NT_166°_NORM_v2.xnra"

# ==================== WINDOW SETTINGS ====================
FIGURE_SIZE = (14, 8)           # Width x Height in inches
DPI = 100                       # Resolution (higher = sharper)

# ==================== RAW DATA APPEARANCE ====================
# --- Markers ---
SHOW_RAW_MARKERS = True                # Show raw data markers?
RAW_MARKER_SHAPE = 'square'            # Options: 'circle', 'square', 'triangle'
RAW_MARKER_FILL = 'hollow'             # Options: 'filled', 'hollow'
RAW_MARKER_SIZE = 2                    # Marker size
RAW_MARKER_COLOR = 'black'              # Marker color
RAW_MARKER_EDGE_WIDTH = 1              # Border thickness for hollow markers

# --- Marker Range ---
RAW_MARKER_RANGE_ENABLED = True       # Limit markers to specific channel range?
RAW_MARKER_RANGE_MIN = 0               # Start channel for markers
RAW_MARKER_RANGE_MAX = 200            # End channel for markers

# --- Line connecting raw data points ---
SHOW_RAW_LINE = True                   # Connect raw data with line?
RAW_LINE_WIDTH = 1.0                   # Line thickness
RAW_LINE_COLOR = 'blue'                # Line color
RAW_LINE_ALPHA = 1.0                   # Transparency (applies to both line and markers)

# ==================== SMOOTHED DATA APPEARANCE ====================
SHOW_SMOOTHED_DATA = True              # Show smoothed data?

# --- Line ---
SMOOTHED_LINE_WIDTH = 1.0              # Line thickness
SMOOTHED_LINE_COLOR = 'red'            # Line color
SMOOTHED_LINE_ALPHA = 1.0              # Transparency

# --- Markers ---
SHOW_SMOOTHED_MARKERS = False           # Show smoothed data markers?
SMOOTHED_MARKER_SHAPE = 'circle'       # Options: 'circle', 'square', 'triangle'
SMOOTHED_MARKER_FILL = 'hollow'        # Options: 'filled', 'hollow'
SMOOTHED_MARKER_SIZE = 3               # Marker size
SMOOTHED_MARKER_COLOR = 'orange'          # Marker color
SMOOTHED_MARKER_EDGE_WIDTH = 1         # Border thickness for hollow markers

# --- Marker Range ---
SMOOTHED_MARKER_RANGE_ENABLED = True  # Limit markers to specific channel range?
SMOOTHED_MARKER_RANGE_MIN = 200          # Start channel for markers
SMOOTHED_MARKER_RANGE_MAX = 400       # End channel for markers

# ==================== PLOT STYLING ====================
BACKGROUND_COLOR = 'white'             # Plot background
GRID_COLOR = 'gray'                    # Grid line color
SHOW_GRID = True                       # Show grid?

# ==================== TEXT SIZES ====================
TITLE_SIZE = 16                        # Main title
AXIS_LABEL_SIZE = 13                   # Axis labels
TICK_LABEL_SIZE = 11                   # Tick numbers
LEGEND_SIZE = 11                       # Legend text
INFO_BOX_SIZE = 10                     # Info box text

# ==================== AXIS LIMITS ====================
CHANNEL_MIN = 0                        # Minimum channel
CHANNEL_MAX = 1600                     # Maximum channel (None = auto)

# ==================== DISPLAY OPTIONS ====================
SHOW_LEGEND = True                     # Show legend?
SHOW_INFO_BOX = True                   # Show info box?

# ============================================================================
# END CONFIGURATION
# ============================================================================


def parse_xnra_file(filepath):
    """
    Parse XNRA file and extract all data and metadata.
    
    Args:
        filepath: Path to .xnra file
        
    Returns:
        dict: Dictionary containing spectrum data and metadata
    """
    print(f"📂 Reading file: {filepath}")
    
    tree = ET.parse(filepath)
    root = tree.getroot()
    
    # XML namespaces
    ns = {
        'idf': 'http://idf.schemas.itn.pt',
        'simnra': 'http://www.simnra.com/simnra'
    }
    
    data = {}
    
    # ========== Extract Metadata ==========
    filename_node = root.find('.//idf:filename', ns)
    data['filename'] = filename_node.text if filename_node is not None else "XNRA Spectrum"
    
    beam_particle_node = root.find('.//idf:beamparticle', ns)
    data['beam_particle'] = beam_particle_node.text if beam_particle_node is not None else None
    
    beam_energy_node = root.find('.//idf:beamenergy', ns)
    data['beam_energy'] = float(beam_energy_node.text) if beam_energy_node is not None else None
    
    scattering_angle_node = root.find('.//idf:scatteringangle', ns)
    data['scattering_angle'] = float(scattering_angle_node.text) if scattering_angle_node is not None else None
    
    detector_node = root.find('.//idf:detectortype', ns)
    data['detector_type'] = detector_node.text if detector_node is not None else None
    
    resolution_node = root.find('.//idf:detectorresolution[1]/idf:resolutionparameters/idf:resolutionparameter', ns)
    data['resolution'] = float(resolution_node.text) if resolution_node is not None else None
    
    # ========== Extract Calibration ==========
    cal_params = root.findall('.//idf:energycalibration[1]/idf:calibrationparameters/idf:calibrationparameter', ns)
    if len(cal_params) >= 2:
        data['cal_offset'] = float(cal_params[0].text)
        data['cal_gain'] = float(cal_params[1].text)
    else:
        data['cal_offset'] = 0
        data['cal_gain'] = 1
        print("⚠️  No calibration found, using default (E = Ch)")
    
    # ========== Extract Raw Data ==========
    raw_data = root.find('.//idf:data/idf:simpledata', ns)
    if raw_data is not None:
        x_text = raw_data.find('idf:x', ns).text
        y_text = raw_data.find('idf:y', ns).text
        data['raw_channels'] = np.array([float(x) for x in x_text.split()])
        data['raw_counts'] = np.array([float(y) for y in y_text.split()])
        print(f"✓ Raw data: {len(data['raw_channels'])} channels")
    else:
        data['raw_channels'] = None
        data['raw_counts'] = None
    
    # ========== Extract Smoothed Data ==========
    smoothed_data = root.find('.//simnra:smootheddata/idf:simpledata', ns)
    if smoothed_data is not None:
        x_text = smoothed_data.find('idf:x', ns).text
        y_text = smoothed_data.find('idf:y', ns).text
        data['smoothed_channels'] = np.array([float(x) for x in x_text.split()])
        data['smoothed_counts'] = np.array([float(y) for y in y_text.split()])
        print(f"✓ Smoothed data: {len(data['smoothed_channels'])} channels")
    else:
        data['smoothed_channels'] = None
        data['smoothed_counts'] = None
    
    return data


def get_marker_style(shape, fill):
    """
    Get matplotlib marker style based on shape and fill options.
    
    Args:
        shape: 'circle', 'square', or 'triangle'
        fill: 'filled' or 'hollow'
        
    Returns:
        tuple: (marker, fillstyle)
    """
    # Map shape names to matplotlib markers
    shape_map = {
        'circle': 'o',
        'square': 's',
        'triangle': '^'
    }
    
    marker = shape_map.get(shape.lower(), 'o')
    
    # Determine fill style
    if fill.lower() == 'hollow':
        fillstyle = 'none'
    else:
        fillstyle = 'full'
    
    return marker, fillstyle


def apply_marker_range(channels, counts, range_enabled, range_min, range_max):
    """
    Filter data to only include points within the specified marker range.
    
    Args:
        channels: Array of channel values
        counts: Array of count values
        range_enabled: Whether to apply range filtering
        range_min: Minimum channel value
        range_max: Maximum channel value
        
    Returns:
        tuple: (filtered_channels, filtered_counts) or (None, None) if no points in range
    """
    if not range_enabled:
        return channels, counts
    
    # Find indices within range
    mask = (channels >= range_min) & (channels <= range_max)
    
    if not np.any(mask):
        return None, None
    
    return channels[mask], counts[mask]


def plot_spectrum(data):
    """
    Create matplotlib plot with dual x-axes and customizable appearance.
    
    Args:
        data: Dictionary containing spectrum data and metadata
    """
    print("\n🎨 Creating plot...")
    
    # ========== Setup Figure ==========
    fig, ax = plt.subplots(figsize=FIGURE_SIZE, dpi=DPI)
    fig.patch.set_facecolor(BACKGROUND_COLOR)
    ax.set_facecolor(BACKGROUND_COLOR)
    
    # ========== Get Data ==========
    raw_channels = data['raw_channels']
    raw_counts = data['raw_counts']
    smoothed_channels = data['smoothed_channels']
    smoothed_counts = data['smoothed_counts']
    
    cal_offset = data['cal_offset']
    cal_gain = data['cal_gain']
    
    # ========== Plot Raw Data ==========
    if raw_channels is not None and raw_counts is not None:
        # Get marker style
        marker, fillstyle = get_marker_style(RAW_MARKER_SHAPE, RAW_MARKER_FILL)
        
        # Determine line style
        if SHOW_RAW_LINE:
            linestyle = '-'
            linewidth = RAW_LINE_WIDTH
        else:
            linestyle = ''  # No line
            linewidth = 0
        
        # Plot the line (full data)
        if SHOW_RAW_LINE:
            ax.plot(raw_channels, raw_counts,
                   linestyle=linestyle,
                   linewidth=linewidth,
                   color=RAW_LINE_COLOR,
                   alpha=RAW_LINE_ALPHA,
                   label='Raw data',
                   zorder=1)
            print(f"✓ Plotted raw data line")
        
        # Plot markers (possibly with range restriction)
        if SHOW_RAW_MARKERS:
            marker_channels, marker_counts = apply_marker_range(
                raw_channels, raw_counts,
                RAW_MARKER_RANGE_ENABLED,
                RAW_MARKER_RANGE_MIN,
                RAW_MARKER_RANGE_MAX
            )
            
            if marker_channels is not None:
                # Handle hollow markers
                if RAW_MARKER_FILL.lower() == 'hollow':
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=RAW_MARKER_SIZE,
                           markerfacecolor='none',
                           markeredgecolor=RAW_MARKER_COLOR,
                           markeredgewidth=RAW_MARKER_EDGE_WIDTH,
                           alpha=RAW_LINE_ALPHA,
                           label='Raw data' if not SHOW_RAW_LINE else '',
                           zorder=3)
                else:
                    # Filled markers
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=RAW_MARKER_SIZE,
                           markerfacecolor=RAW_MARKER_COLOR,
                           markeredgecolor=RAW_MARKER_COLOR,
                           alpha=RAW_LINE_ALPHA,
                           label='Raw data' if not SHOW_RAW_LINE else '',
                           zorder=3)
                
                range_info = f" (channels {RAW_MARKER_RANGE_MIN}-{RAW_MARKER_RANGE_MAX})" if RAW_MARKER_RANGE_ENABLED else ""
                print(f"✓ Plotted raw data markers: {RAW_MARKER_SHAPE}, {RAW_MARKER_FILL}{range_info}")
            else:
                print(f"⚠️  No raw data points in marker range {RAW_MARKER_RANGE_MIN}-{RAW_MARKER_RANGE_MAX}")
    
    # ========== Plot Smoothed Data ==========
    if SHOW_SMOOTHED_DATA and smoothed_channels is not None and smoothed_counts is not None:
        # Get marker style for smoothed data
        marker, fillstyle = get_marker_style(SMOOTHED_MARKER_SHAPE, SMOOTHED_MARKER_FILL)
        
        # Plot the line (full data)
        ax.plot(smoothed_channels, smoothed_counts,
               '-',
               color=SMOOTHED_LINE_COLOR,
               linewidth=SMOOTHED_LINE_WIDTH,
               alpha=SMOOTHED_LINE_ALPHA,
               label='Smoothed data',
               zorder=2)
        print("✓ Plotted smoothed data line")
        
        # Plot markers (possibly with range restriction)
        if SHOW_SMOOTHED_MARKERS:
            marker_channels, marker_counts = apply_marker_range(
                smoothed_channels, smoothed_counts,
                SMOOTHED_MARKER_RANGE_ENABLED,
                SMOOTHED_MARKER_RANGE_MIN,
                SMOOTHED_MARKER_RANGE_MAX
            )
            
            if marker_channels is not None:
                # Handle hollow markers
                if SMOOTHED_MARKER_FILL.lower() == 'hollow':
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=SMOOTHED_MARKER_SIZE,
                           markerfacecolor='none',
                           markeredgecolor=SMOOTHED_MARKER_COLOR,
                           markeredgewidth=SMOOTHED_MARKER_EDGE_WIDTH,
                           alpha=SMOOTHED_LINE_ALPHA,
                           zorder=4)
                else:
                    # Filled markers
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=SMOOTHED_MARKER_SIZE,
                           markerfacecolor=SMOOTHED_MARKER_COLOR,
                           markeredgecolor=SMOOTHED_MARKER_COLOR,
                           alpha=SMOOTHED_LINE_ALPHA,
                           zorder=4)
                
                range_info = f" (channels {SMOOTHED_MARKER_RANGE_MIN}-{SMOOTHED_MARKER_RANGE_MAX})" if SMOOTHED_MARKER_RANGE_ENABLED else ""
                print(f"✓ Plotted smoothed data markers: {SMOOTHED_MARKER_SHAPE}, {SMOOTHED_MARKER_FILL}{range_info}")
            else:
                print(f"⚠️  No smoothed data points in marker range {SMOOTHED_MARKER_RANGE_MIN}-{SMOOTHED_MARKER_RANGE_MAX}")
    
    # ========== Setup Axes ==========
    # Bottom axis: Channel
    ax.set_xlabel('Channel', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
    ax.set_ylabel('Counts', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
    ax.tick_params(axis='both', labelsize=TICK_LABEL_SIZE)
    
    # Set channel limits
    if CHANNEL_MAX is not None:
        ax.set_xlim(CHANNEL_MIN, CHANNEL_MAX)
    elif raw_channels is not None:
        ax.set_xlim(CHANNEL_MIN, raw_channels.max())
    
    # Top axis: Energy (dual x-axis)
    ax_top = ax.secondary_xaxis('top', functions=(
        lambda ch: cal_offset + cal_gain * ch,      # channel → energy
        lambda en: (en - cal_offset) / cal_gain     # energy → channel
    ))
    ax_top.set_xlabel('Energy [keV]', fontsize=AXIS_LABEL_SIZE, fontweight='bold')
    ax_top.tick_params(axis='x', labelsize=TICK_LABEL_SIZE)
    
    print("✓ Added dual x-axes (Channel + Energy)")
    
    # ========== Title ==========
    ax.set_title(data['filename'], fontsize=TITLE_SIZE, fontweight='bold', pad=40)
    
    # ========== Grid ==========
    if SHOW_GRID:
        ax.grid(True, alpha=0.3, color=GRID_COLOR, linestyle='--', linewidth=0.5)
    
    # ========== Legend ==========
    if SHOW_LEGEND:
        # Get unique labels to avoid duplicates
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), 
                 loc='upper right', fontsize=LEGEND_SIZE, framealpha=0.9)
    
    # ========== Info Box ==========
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
        
        info_lines.append(f"Cal: E = {cal_offset:.1f} + {cal_gain:.3f}×Ch")
        
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
    
    # ========== Finalize ==========
    plt.tight_layout()
    
    print("✓ Plot created successfully!")
    print("\n" + "="*60)
    print("PLOT READY - Close the window to exit")
    print("="*60)
    
    plt.show()


def main():
    """Main function - Run the XNRA viewer"""
    
    print("="*60)
    print("XNRA SPECTRUM VIEWER - Enhanced with Marker Control")
    print("="*60)
    
    # Check if file exists
    filepath = Path(FILE_PATH)
    if not filepath.exists():
        print(f"\n❌ ERROR: File not found!")
        print(f"   Looked for: {filepath.absolute()}")
        print(f"\n💡 TIP: Update FILE_PATH at the top of this script")
        return
    
    try:
        # Parse file
        data = parse_xnra_file(filepath)
        
        # Check if we have data
        if data['raw_counts'] is None and data['smoothed_counts'] is None:
            print("\n❌ ERROR: No spectrum data found in file!")
            return
        
        # Create plot
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
