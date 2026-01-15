### XNRA Spectrum Viewer for NRA, RBS, and ERDA data ###
# Core plotting module - designed to be imported by GUI

# Difference from v3: changed defaults for info and legend box positions

# Script developed by Berke Santos & Giuseppe Legrottaglie
# Last updated: 17/12/2025

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================================

# HELPER FUNCTIONS

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
        print("✓ Calibration parameters found")
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
    
    shape_map = {
        'circle': 'o',
        'square': 's',
        'triangle': '^'
    }
    
    marker = shape_map.get(shape.lower(), 'o')
    fillstyle = 'none' if fill.lower() == 'hollow' else 'full'
    
    return marker, fillstyle


def apply_marker_range(channels, counts, range_enabled, range_mode, range_min, range_max, cal_offset, cal_gain):
    
    if not range_enabled:
        return channels, counts
    
    # Convert energy range to channel range if needed
    if range_mode.lower() == 'energy':
        ch_min = (range_min - cal_offset) / cal_gain
        ch_max = (range_max - cal_offset) / cal_gain
    else:
        ch_min = range_min
        ch_max = range_max
    
    # Find indices within range
    mask = (channels >= ch_min) & (channels <= ch_max)
    
    if not np.any(mask):
        return None, None
    
    return channels[mask], counts[mask]


# ============================================================================
# MAIN PLOTTING FUNCTION

def plot_xnra_spectrum(
    filepath,
    title=None,
    # Figure settings
    figure_size=(12, 6),
    dpi=100,
    show_legend=True,
    show_info_box=True,
    info_box_position='left',
    info_box_x_offset=0.05,  
    info_box_y_offset=0.95,  
    # Raw data settings
    show_raw_data=True,
    show_raw_markers=True,
    raw_marker_shape='square',
    raw_marker_fill='hollow',
    raw_marker_color='black',
    raw_marker_size=2,
    raw_marker_edge_width=1,
    raw_marker_range_enabled=True,
    raw_marker_range_mode='channel',
    raw_marker_range_min=0,
    raw_marker_range_max=200,
    show_raw_line=True,
    raw_line_color='blue',
    raw_line_width=1.0,
    raw_line_alpha=1.0,
    # Smoothed data settings
    show_smoothed_data=True,
    show_smoothed_markers=True,
    smoothed_marker_shape='triangle',
    smoothed_marker_fill='hollow',
    smoothed_marker_color='orange',
    smoothed_marker_size=3,
    smoothed_marker_edge_width=1,
    smoothed_marker_range_enabled=True,
    smoothed_marker_range_mode='channel',
    smoothed_marker_range_min=200,
    smoothed_marker_range_max=400,
    show_smoothed_line=True,
    smoothed_line_color='red',
    smoothed_line_width=1.0,
    smoothed_line_alpha=1.0,
    # Plot appearance
    background_color='white',
    show_grid=True,
    grid_color='white',
    title_size=16,
    axis_label_size=12,
    tick_label_size=10,
    legend_size=10,
    info_box_size=10,
    title_color='black',
    axis_label_color='black',
    tick_label_color='black',
    legend_text_color='black',
    info_box_text_color='black',
    tick_length=6,
    tick_width=1.0,
    tick_color='black',
    tick_direction='in',
    # Axis settings
    x_axis_mode='energy',
    channel_min=None,
    channel_max=None,
    energy_min=None,
    energy_max=None,
    counts_min=None,
    counts_max=None,
    # Display mode
    show_plot=True
):
    
    # Parse the data file
    data = parse_xnra_file(filepath)
    
    # Check if we have data
    if data['raw_counts'] is None and data['smoothed_counts'] is None:
        raise ValueError("No spectrum data found in file!")
    
    # ========== Setup Figure ==========
    #fig, ax = plt.subplots(figsize=figure_size, dpi=dpi)
    fig, ax = plt.subplots(figsize=figure_size, dpi=dpi, constrained_layout=True)

    fig.get_layout_engine().set(h_pad=4/72, w_pad=4/72, 
                            hspace=0.05, wspace=0.05,
                            rect=[0, 0, 1, 0.95])  
    
    fig.patch.set_facecolor(background_color)
    ax.set_facecolor(background_color)
    
    # ========== Get Data ==========
    raw_channels = data['raw_channels']
    raw_counts = data['raw_counts']
    smoothed_channels = data['smoothed_channels']
    smoothed_counts = data['smoothed_counts']
    
    cal_offset = data['cal_offset']
    cal_gain = data['cal_gain']
    
    # ========== Plot Raw Data ==========
    if show_raw_data and raw_channels is not None and raw_counts is not None:
        marker, fillstyle = get_marker_style(raw_marker_shape, raw_marker_fill)
        
        if show_raw_line:
            ax.plot(raw_channels, raw_counts,
                   linestyle='-',
                   linewidth=raw_line_width,
                   color=raw_line_color,
                   alpha=raw_line_alpha,
                   label='Raw data',
                   zorder=1)
            print("✓ Plotted raw data line")
                        
        if show_raw_markers:
            marker_channels, marker_counts = apply_marker_range(
                raw_channels, raw_counts,
                raw_marker_range_enabled,
                raw_marker_range_mode,
                raw_marker_range_min,
                raw_marker_range_max,
                cal_offset,
                cal_gain
            )
            
            if marker_channels is not None:
                if raw_marker_fill.lower() == 'hollow':
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=raw_marker_size,
                           markerfacecolor='none',
                           markeredgecolor=raw_marker_color,
                           markeredgewidth=raw_marker_edge_width,
                           alpha=raw_line_alpha,
                           label='Raw data' if not show_raw_line else '',
                           zorder=3)
                else:
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=raw_marker_size,
                           markerfacecolor=raw_marker_color,
                           markeredgecolor=raw_marker_color,
                           alpha=raw_line_alpha,
                           label='Raw data' if not show_raw_line else '',
                           zorder=3)
                
                if raw_marker_range_enabled:
                    unit = 'keV' if raw_marker_range_mode.lower() == 'energy' else 'channels'
                    print(f"✓ Plotted raw markers: {raw_marker_shape}, {raw_marker_fill} ({raw_marker_range_min}-{raw_marker_range_max} {unit})")
            else:
                unit = 'keV' if raw_marker_range_mode.lower() == 'energy' else 'channels'
                print(f"⚠️  No raw data points in range {raw_marker_range_min}-{raw_marker_range_max} {unit}")
                
    # ========== Plot Smoothed Data ==========
    if show_smoothed_data and smoothed_channels is not None and smoothed_counts is not None:
        marker, fillstyle = get_marker_style(smoothed_marker_shape, smoothed_marker_fill)
        
        if show_smoothed_line:
            ax.plot(smoothed_channels, smoothed_counts, 
                    linestyle='-',
                    linewidth=smoothed_line_width,
                    color=smoothed_line_color,
                    alpha=smoothed_line_alpha,
                    label='Smoothed data',
                    zorder=1)
            print("✓ Plotted smoothed data line")
        
        if show_smoothed_markers:
            marker_channels, marker_counts = apply_marker_range(
                smoothed_channels, smoothed_counts,
                smoothed_marker_range_enabled,
                smoothed_marker_range_mode,
                smoothed_marker_range_min,
                smoothed_marker_range_max,
                cal_offset,
                cal_gain
            )
            
            if marker_channels is not None:
                if smoothed_marker_fill.lower() == 'hollow':
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=smoothed_marker_size,
                           markerfacecolor='none',
                           markeredgecolor=smoothed_marker_color,
                           markeredgewidth=smoothed_marker_edge_width,
                           alpha=smoothed_line_alpha,
                           zorder=3)
                else:
                    ax.plot(marker_channels, marker_counts,
                           linestyle='',
                           marker=marker,
                           markersize=smoothed_marker_size,
                           markerfacecolor=smoothed_marker_color,
                           markeredgecolor=smoothed_marker_color,
                           alpha=smoothed_line_alpha,
                           zorder=3)
                
                if smoothed_marker_range_enabled:
                    unit = 'keV' if smoothed_marker_range_mode.lower() == 'energy' else 'channels'
                    print(f"✓ Plotted smoothed markers: {smoothed_marker_shape}, {smoothed_marker_fill} ({smoothed_marker_range_min}-{smoothed_marker_range_max} {unit})")
            else:
                unit = 'keV' if smoothed_marker_range_mode.lower() == 'energy' else 'channels'
                print(f"⚠️  No smoothed data points in range {smoothed_marker_range_min}-{smoothed_marker_range_max} {unit}")

    # ========== Setup Axes ==========
    ax.set_xlabel('Channel', fontsize=axis_label_size, fontweight='bold', color=axis_label_color)
    ax.set_ylabel('Counts', fontsize=axis_label_size, fontweight='bold', color=axis_label_color)
    ax.tick_params(axis='both', 
                   labelsize=tick_label_size,
                   labelcolor=tick_label_color,
                   length=tick_length,
                   width=tick_width,
                   color=tick_color,
                   direction=tick_direction)
    
    # Set X-axis limits
    if x_axis_mode.lower() == 'energy':
        if energy_min is not None and energy_max is not None:
            ch_min = (energy_min - cal_offset) / cal_gain
            ch_max = (energy_max - cal_offset) / cal_gain
            ax.set_xlim(ch_min, ch_max)
        elif energy_min is not None:
            ch_min = (energy_min - cal_offset) / cal_gain
            ax.set_xlim(left=ch_min)
        elif energy_max is not None:
            ch_max = (energy_max - cal_offset) / cal_gain
            ax.set_xlim(right=ch_max)
        else:
            if raw_channels is not None:
                ax.set_xlim(raw_channels.min(), raw_channels.max())
    else:  # channel mode
        if channel_min is not None and channel_max is not None:
            ax.set_xlim(channel_min, channel_max)
        elif channel_min is not None:
            ax.set_xlim(left=channel_min)
        elif channel_max is not None:
            ax.set_xlim(right=channel_max)
        else:
            if raw_channels is not None:
                ax.set_xlim(raw_channels.min(), raw_channels.max())
    
    # Set Y-axis limits
    if counts_min is not None or counts_max is not None:
        current_ylim = ax.get_ylim()
        y_min = counts_min if counts_min is not None else current_ylim[0]
        y_max = counts_max if counts_max is not None else current_ylim[1]
        ax.set_ylim(y_min, y_max)
    
    # Top axis: Energy
    ax_top = ax.secondary_xaxis('top', functions=(
        lambda ch: cal_offset + cal_gain * ch,
        lambda en: (en - cal_offset) / cal_gain
    ))
    ax_top.set_xlabel('Energy [keV]', fontsize=axis_label_size, fontweight='bold', color=axis_label_color)
    ax_top.tick_params(axis='x', 
                       labelsize=tick_label_size,
                       labelcolor=tick_label_color,
                       length=tick_length,
                       width=tick_width,
                       color=tick_color,
                       direction=tick_direction)
        
    # ========== Title ==========
    plot_title = title if title is not None else data['filename']
    title_pad = max(10, title_size * 0.8)  # scales with fontsize
    ax.set_title(plot_title, fontsize=title_size, fontweight='bold', color=title_color, pad=title_pad)
    
    # ========== Grid ==========
    if show_grid:
        ax.grid(True, alpha=0.3, color=grid_color, linestyle='--', linewidth=0.5)
    
    # ========== Legend ==========
    if show_legend:
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(),
                      loc='upper right',
                      bbox_to_anchor=(0.95, 0.85),  # Fine-tuned to avoid top axis
                      fontsize=legend_size,
                      framealpha=0.95,
                      edgecolor='black',
                      fancybox=True,
                      shadow=True)
        
        for text in legend.get_texts():
            text.set_color(legend_text_color)
    
    # ========== Info Box ==========
    if show_info_box:
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

            x_pos = info_box_x_offset
            y_pos = info_box_y_offset
            
            if info_box_position.lower() == 'right':
                h_align = 'right'
            else:
                h_align = 'left'
            
            ax.text(x_pos, y_pos, info_text,
                   transform=ax.transAxes,
                   fontsize=info_box_size,
                   color=info_box_text_color,
                   verticalalignment='top',
                   horizontalalignment=h_align,
                   bbox=dict(boxstyle='round',
                           facecolor='wheat',
                           alpha=0.8,
                           edgecolor='black',
                           linewidth=1))
    
    # ========== Finalize ==========
    #plt.tight_layout()
    
    print("✓ Plot created successfully!")
    print("="*60)
    
    if show_plot:
        plt.show()
    
    return fig
