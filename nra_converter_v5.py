### XNRA Spectrum Viewer for NRA, RBS, and ERDA data ###
# Single file version with no UI

## Create a full comment section detailing every part of the code or PDF

# Script developed by Berke Santos & Giuseppe Legrottaglie
# Last updated: 13/12/2025

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import traceback

# ============================================================================
# USER INPUTS + CONFIGURATION 
# ============================================================================

# File path (xnra complete filepath)
#FILE_PATH = r"C:\Users\berke\Desktop\Giuseppe\input\NRA\NRA Re 550 RIC 166°.xnra"
FILE_PATH = r"C:\Users\berke\Desktop\Giuseppe\input\RBS\RBS_NT_166°_NORM_v2.xnra"

# Visualization window settings 
FIGURE_SIZE = (12, 6)                  # Width x Height in inches 
DPI = 100                              # Resolution (higher = sharper)
# Stands for Dots Per Inch (dots as in pixels)

SHOW_LEGEND = True                     # True or False
SHOW_INFO_BOX = True                   # True or False
INFO_BOX_POSITION = 'left'            # 'left' or 'right'

# Raw/Experimental data settings
SHOW_RAW_DATA = True                   # True or False

# --- Markers ---
SHOW_RAW_MARKERS = True                # True or False
RAW_MARKER_SHAPE = 'square'            # 'circle', 'square', or 'triangle'
RAW_MARKER_FILL = 'hollow'             # 'filled' or 'hollow'
RAW_MARKER_COLOR = 'black'             # set colour
RAW_MARKER_SIZE = 2                    # set size (linear dimension, as in the diameter - in points)
# 1 point corresponds to 1/72 inches
RAW_MARKER_EDGE_WIDTH = 1              # set border thickness of the hollow markers

# --- Marker Range ---
RAW_MARKER_RANGE_ENABLED = True        # Limit markers to specific channel/energy range?
RAW_MARKER_RANGE_MODE = 'channel'      # 'channel' or 'energy'
RAW_MARKER_RANGE_MIN = 0               # set min value of channel/energy range
RAW_MARKER_RANGE_MAX = 200             # set max value of channel/energy range

# --- Line ---
SHOW_RAW_LINE = False                   # Connect raw data with line?
RAW_LINE_COLOR = 'blue'                # set line colour
RAW_LINE_WIDTH = 1.0                   # set line thickness
RAW_LINE_ALPHA = 1.0                   # set line and marker transparency

# Smoothed/Simulated data settings
SHOW_SMOOTHED_DATA = True               # True or False

# --- Markers ---
SHOW_SMOOTHED_MARKERS = True           # True or False
SMOOTHED_MARKER_SHAPE = 'triangle'        # 'circle', 'square', or 'triangle'
SMOOTHED_MARKER_FILL = 'hollow'         # 'filled' or 'hollow'
SMOOTHED_MARKER_COLOR = 'orange'        # set colour
SMOOTHED_MARKER_SIZE = 3                # set size (linear dimension, as in the diameter - in points)
# 1 point corresponds to 1/72 inches
SMOOTHED_MARKER_EDGE_WIDTH = 1          # set border thickness of the hollow markers

# --- Marker Range ---
SMOOTHED_MARKER_RANGE_ENABLED = True    # Limit markers to specific channel/energy range?
SMOOTHED_MARKER_RANGE_MODE = 'channel'  # 'channel' or 'energy'
SMOOTHED_MARKER_RANGE_MIN = 200         # set min value of channel/energy range
SMOOTHED_MARKER_RANGE_MAX = 400         # set max value of channel/energy range

# --- Line ---
SHOW_SMOOTHED_LINE = True              # True or False
SMOOTHED_LINE_WIDTH = 1.0              # set line thickness
SMOOTHED_LINE_COLOR = 'red'            # set line colour
SMOOTHED_LINE_ALPHA = 1.0              # set line and marker transparency

# Plot settings
#--- Bkg ---
BACKGROUND_COLOR = 'white'              # set plot background colour

#--- Grid ---
SHOW_GRID = True                        # True or False
GRID_COLOR = 'white'                    # set grid line colour 

#--- Text sizes ---
TITLE_SIZE = 16                         # plot title
AXIS_LABEL_SIZE = 12                    # axis labels
TICK_LABEL_SIZE = 10                    # tick numbers
LEGEND_SIZE = 10                        # legend text (for the lines)
INFO_BOX_SIZE = 10                      # info box text

#--- Text colours ---
TITLE_COLOR = 'black'                   # plot title colour
AXIS_LABEL_COLOR = 'black'              # axis labels colour
TICK_LABEL_COLOR = 'violet'              # tick numbers colour
LEGEND_TEXT_COLOR = 'black'             # legend text colour
INFO_BOX_TEXT_COLOR = 'black'           # info box text colour

#--- Tick marks (the small lines on axes) ---
TICK_LENGTH = 6                         # length of tick marks (in points)
TICK_WIDTH = 1.0                        # thickness of tick marks (in points)
TICK_COLOR = 'orange'                    # colour of tick marks
TICK_DIRECTION = 'in'                   # 'in', 'out', or 'inout'

#--- X-axis settings ---
X_AXIS_MODE = 'energy'                  # 'channel' or 'energy'

# Channel-based limits 
#(None = full display)
CHANNEL_MIN = None                      # set minimum channel displayed
CHANNEL_MAX = None                      # set maximum channel displayed
# Energy-based limits 
ENERGY_MIN = None                       # set minimum energy [keV] displayed 
ENERGY_MAX = None                       # set maximum energy [keV] displayed

#--- Y-axis limits ---
#(None = full display)
COUNTS_MIN = None                      # set minimum counts 
COUNTS_MAX = None                      # set maxixum counts 

# ============================================================================
# FUNCTIONS
# ============================================================================

def parse_xnra_file(filepath):

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
        print("Calibration parameters found !")

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


def apply_marker_range(channels, counts, range_enabled, range_mode, range_min, range_max, cal_offset, cal_gain):
    
    if not range_enabled:
        return channels, counts
    
    # Convert energy range to channel range if needed
    if range_mode.lower() == 'energy':
        # Energy → Channel conversion: Ch = (E - offset) / gain
        ch_min = (range_min - cal_offset) / cal_gain
        ch_max = (range_max - cal_offset) / cal_gain
    else:  # channel mode
        ch_min = range_min
        ch_max = range_max
    
    # Find indices within range (inclusive)
    mask = (channels >= ch_min) & (channels <= ch_max)
    
    # Check if any points are in range
    if not np.any(mask):
        return None, None
    
    return channels[mask], counts[mask]


def plot_spectrum(data):
    
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
    if SHOW_RAW_DATA and raw_channels is not None and raw_counts is not None:
        # Get marker style
        marker, fillstyle = get_marker_style(RAW_MARKER_SHAPE, RAW_MARKER_FILL)
        
        # Determine line style
        if SHOW_RAW_LINE:
            linestyle = '-'
            linewidth = RAW_LINE_WIDTH

            ax.plot(raw_channels, raw_counts,
                   linestyle=linestyle,
                   linewidth=linewidth,
                   color=RAW_LINE_COLOR,
                   alpha=RAW_LINE_ALPHA,
                   label='Raw data',
                   zorder=1)
            print("✓ Plotted raw data line")
                        
        # Plot markers 
        if SHOW_RAW_MARKERS:
            # Apply marker range filter
            marker_channels, marker_counts = apply_marker_range(
                raw_channels, raw_counts,
                RAW_MARKER_RANGE_ENABLED,
                RAW_MARKER_RANGE_MODE,
                RAW_MARKER_RANGE_MIN,
                RAW_MARKER_RANGE_MAX,
                cal_offset,
                cal_gain
            )
            
            # Check if we have valid marker data
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
                
                # Print status with appropriate units
                if RAW_MARKER_RANGE_ENABLED:
                    unit = 'keV' if RAW_MARKER_RANGE_MODE.lower() == 'energy' else 'channels'
                    range_info = f" ({RAW_MARKER_RANGE_MIN}-{RAW_MARKER_RANGE_MAX} {unit})"
                else:
                    range_info = ""
                    
                print(f"✓ Plotted raw data markers: {RAW_MARKER_SHAPE}, {RAW_MARKER_FILL}{range_info}")
            else:
                # No markers in range - print warning
                unit = 'keV' if RAW_MARKER_RANGE_MODE.lower() == 'energy' else 'channels'
                print(f"⚠️  No raw data points in marker range {RAW_MARKER_RANGE_MIN}-{RAW_MARKER_RANGE_MAX} {unit}")
                
    # ========== Plot Smoothed Data ==========
    if SHOW_SMOOTHED_DATA and smoothed_channels is not None and smoothed_counts is not None:
        # Get marker style for smoothed data
        marker, fillstyle = get_marker_style(SMOOTHED_MARKER_SHAPE, SMOOTHED_MARKER_FILL)
        
        # Determine line style
        if SHOW_SMOOTHED_LINE:
            linestyle = '-'
            linewidth = SMOOTHED_LINE_WIDTH

            # Plot the line (full data)
            ax.plot(smoothed_channels, smoothed_counts, 
                    linestyle=linestyle,
                    linewidth=linewidth,
                    color=SMOOTHED_LINE_COLOR,
                    alpha=SMOOTHED_LINE_ALPHA,
                    label='Smoothed data',
                zorder=1)
            print("✓ Plotted smoothed data line")
        
        # Plot markers 
        if SHOW_SMOOTHED_MARKERS:
            # Apply marker range filter
            marker_channels, marker_counts = apply_marker_range(
                smoothed_channels, smoothed_counts,
                SMOOTHED_MARKER_RANGE_ENABLED,
                SMOOTHED_MARKER_RANGE_MODE,
                SMOOTHED_MARKER_RANGE_MIN,
                SMOOTHED_MARKER_RANGE_MAX,
                cal_offset,
                cal_gain
            )
            
            # Check if we have valid marker data
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
                           zorder=3)
                else:
                    # Filled markers
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=SMOOTHED_MARKER_SIZE,
                           markerfacecolor=SMOOTHED_MARKER_COLOR,
                           markeredgecolor=SMOOTHED_MARKER_COLOR,
                           alpha=SMOOTHED_LINE_ALPHA,
                           zorder=3)
                
                # Print status with appropriate units
                if SMOOTHED_MARKER_RANGE_ENABLED:
                    unit = 'keV' if SMOOTHED_MARKER_RANGE_MODE.lower() == 'energy' else 'channels'
                    range_info = f" ({SMOOTHED_MARKER_RANGE_MIN}-{SMOOTHED_MARKER_RANGE_MAX} {unit})"
                else:
                    range_info = ""
                    
                print(f"✓ Plotted smoothed data markers: {SMOOTHED_MARKER_SHAPE}, {SMOOTHED_MARKER_FILL}{range_info}")
            else:
                # No markers in range - print warning
                unit = 'keV' if SMOOTHED_MARKER_RANGE_MODE.lower() == 'energy' else 'channels'
                print(f"⚠️  No smoothed data points in marker range {SMOOTHED_MARKER_RANGE_MIN}-{SMOOTHED_MARKER_RANGE_MAX} {unit}")

    # ========== Setup Axes ==========
    # Bottom axis: Channel
    ax.set_xlabel('Channel', fontsize=AXIS_LABEL_SIZE, fontweight='bold', color=AXIS_LABEL_COLOR)
    ax.set_ylabel('Counts', fontsize=AXIS_LABEL_SIZE, fontweight='bold', color=AXIS_LABEL_COLOR)
    ax.tick_params(axis='both', 
                   labelsize=TICK_LABEL_SIZE,
                   labelcolor=TICK_LABEL_COLOR,
                   length=TICK_LENGTH,
                   width=TICK_WIDTH,
                   color=TICK_COLOR,
                   direction=TICK_DIRECTION)
    
    # Set X-axis limits based on mode (channel or energy)
    if X_AXIS_MODE.lower() == 'energy':
        # Convert energy limits to channel limits
        if ENERGY_MIN is not None and ENERGY_MAX is not None:
            ch_min = (ENERGY_MIN - cal_offset) / cal_gain
            ch_max = (ENERGY_MAX - cal_offset) / cal_gain
            ax.set_xlim(ch_min, ch_max)

        elif ENERGY_MIN is not None:
            ch_min = (ENERGY_MIN - cal_offset) / cal_gain
            ax.set_xlim(left=ch_min)

        elif ENERGY_MAX is not None:
            ch_max = (ENERGY_MAX - cal_offset) / cal_gain
            ax.set_xlim(right=ch_max)

        else:
            if raw_channels is not None:
                ax.set_xlim(raw_channels.min(), raw_channels.max())
                
    else:  # channel mode
        
        if CHANNEL_MIN is not None and CHANNEL_MAX is not None:
            ax.set_xlim(CHANNEL_MIN, CHANNEL_MAX)

        elif CHANNEL_MIN is not None:
            ax.set_xlim(left=CHANNEL_MIN)

        elif CHANNEL_MAX is not None:
            ax.set_xlim(right=CHANNEL_MAX)

        else:
            if raw_channels is not None:
                ax.set_xlim(raw_channels.min(), raw_channels.max())
    
    # Set counts limits (Y-axis)
    if COUNTS_MIN is not None or COUNTS_MAX is not None:

        current_ylim = ax.get_ylim()
        y_min = COUNTS_MIN if COUNTS_MIN is not None else current_ylim[0]
        y_max = COUNTS_MAX if COUNTS_MAX is not None else current_ylim[1]
        ax.set_ylim(y_min, y_max)
    
    # Top axis: Energy (dual x-axis)
    ax_top = ax.secondary_xaxis('top', functions=(
        lambda ch: cal_offset + cal_gain * ch,      # channel → energy
        lambda en: (en - cal_offset) / cal_gain     # energy → channel
    ))
    ax_top.set_xlabel('Energy [keV]', fontsize=AXIS_LABEL_SIZE, fontweight='bold', color=AXIS_LABEL_COLOR)
    ax_top.tick_params(axis='x', 
                       labelsize=TICK_LABEL_SIZE,
                       labelcolor=TICK_LABEL_COLOR,
                       length=TICK_LENGTH,
                       width=TICK_WIDTH,
                       color=TICK_COLOR,
                       direction=TICK_DIRECTION)
        
    # ========== Title ==========
    ax.set_title(data['filename'], fontsize=TITLE_SIZE, fontweight='bold', color=TITLE_COLOR, pad=40)
    
    # ========== Grid ==========
    if SHOW_GRID:
        ax.grid(True, alpha=0.3, color=GRID_COLOR, linestyle='--', linewidth=0.5)
    
    # ========== Legend ==========
    # Legend always on the right
    if SHOW_LEGEND:
        # Get unique labels to avoid duplicates
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), 
                 loc='upper right', fontsize=LEGEND_SIZE, framealpha=0.9)
        
        for text in legend.get_texts():
            text.set_color(LEGEND_TEXT_COLOR)
    
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
            
            # Position info box based on configuration
            if INFO_BOX_POSITION.lower() == 'right':
                # On the right side, below the legend
                x_pos = 0.99
                y_pos = 0.85 if SHOW_LEGEND else 0.98  # Lower if legend is shown
                h_align = 'right'
            else:  # 'left'
                # On the left side
                x_pos = 0.01
                y_pos = 0.98
                h_align = 'left'
            
            ax.text(x_pos, y_pos, info_text,
                   transform=ax.transAxes,
                   fontsize=INFO_BOX_SIZE,
                   color=INFO_BOX_TEXT_COLOR,
                   verticalalignment='top',
                   horizontalalignment=h_align,
                   bbox=dict(boxstyle='round',
                           facecolor='wheat',
                           alpha=0.8,
                           edgecolor='black',
                           linewidth=1))
    
    # ========== Finalize ==========
    plt.tight_layout()
    
    print("✓ Plot created successfully!")
    print("="*60)
    
    plt.show()


def main():
    
    print("XNRA SPECTRUM VIEWER - Enhanced with Marker Control")
    print("="*60)
    
    # Check if file exists
    filepath = Path(FILE_PATH)
    if not filepath.exists():
        print(f"\n ERROR: File not found!")
        return
    
    try:
        # Parse file
        data = parse_xnra_file(filepath)
        
        # Check if we have data
        if data['raw_counts'] is None and data['smoothed_counts'] is None:
            print("\n ERROR: No spectrum data found in file!")
            return
        
        # Create plot
        plot_spectrum(data)
        
    except ET.ParseError as e:
        print(f"\n  ERROR: Invalid XML file!")
        print(f"   {e}")
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
