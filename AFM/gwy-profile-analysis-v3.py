# AFM Profile Analysis Tool (raw data from Gwyddion)
# version 2 is the latest update of version 1

# PURPOSE:
# Analyzes AFM (Atomic Force Microscopy) profile data exported from Gwyddion.
# Calculates vertical and horizontal amplitude metrics for surface features.

# METRICS CALCULATED:
# 1. Vertical Amplitude Left (V_left):
#    - Left local max (nearest to global min on left side) minus global min
# 2. Vertical Amplitude Right (V_right):
#    - Right local max (nearest to global min on right side) minus global min
# 3. Horizontal Amplitude (H):
#    - Distance between local maxima on either side of the global minimum
#    - Tie-breaking: When multiple points share the max value (plateaus),
#      selects the point closest to the minimum to avoid overestimation
# 4. Width (W):
#    - Builds a cubic spline through the profile data
#    - Identifies the smaller of the two local maxima
#    - Finds where the spline crosses z = z_ref on the opposite side of
#      the global minimum (z_ref = smaller max's z-value)
#    - Returns the horizontal distance |Δx| between the reference point
#      and the spline intersection (perfectly horizontal by construction)
#    - Falls back to nearest discrete point if intersection not found

# INPUT:
# - Gwyddion-exported files with x and y columns (need to be converted to .txt)
# - Files should be in a single folder (user provides the directory to this folder)

# OUTPUT:
# 1. Individual annotated plots for each profile (saved in the input-directory/plots/ directory)
#    - Plotted in user-selected unit (meters or nanometers)
#    - Color-coded markers (red=global max, blue=global min, orange=local max,
#      cyan=width reference point)
#    - Visual amplitude annotations
# 2. CSV file with per-file results and summary statistics
#    - Includes V_High and V_Low (max/min of V_left and V_right per profile)

# Author: Berke Santos & Giuseppe Legrottaglie
# Developed with the help of Claude.AI and ChatGPT version 5.2
# Created: 16/02/2026

import os
import re
import glob
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


def load_profile_txt(path: str):

    xs, zs = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#") or s.startswith(";") or s.lower().startswith("gwyddion"):
                continue

            nums = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
            if len(nums) >= 2:
                xs.append(float(nums[0]))
                zs.append(float(nums[1]))

    if len(xs) < 3:
        raise ValueError("Not enough numeric data lines found.")

    x = np.array(xs)
    z = np.array(zs)

    # sort by x just in case
    order = np.argsort(x)
    return x[order], z[order]


def find_local_maxima_around_min(x: np.ndarray, y: np.ndarray):
    
    i_min = int(np.argmin(y))

    # Left local maximum (closest to i_min if plateau)
    if i_min > 0:
        left_section = y[: i_min + 1]
        max_val_left = np.max(left_section)
        candidates_left = np.where(left_section == max_val_left)[0]
        i_left = int(candidates_left[-1])  # rightmost = closest to min
    else:
        i_left = 0

    # Right local maximum (closest to i_min if plateau)
    if i_min < len(y) - 1:
        right_section = y[i_min:]
        max_val_right = np.max(right_section)
        candidates_right = np.where(right_section == max_val_right)[0]
        i_right = int(i_min + candidates_right[0])  # leftmost = closest to min
    else:
        i_right = len(y) - 1

    return i_left, i_min, i_right


def horizontal_amplitude_around_min(x: np.ndarray, y: np.ndarray):

    i_left, i_min, i_right = find_local_maxima_around_min(x, y)
    width = float(abs(x[i_right] - x[i_left]))

    return width, (i_left, i_min, i_right)


def vertical_amplitudes_left_right(x, z, idxs):
    
    i_left, i_min, i_right = idxs
    v_left = float(z[i_left] - z[i_min])
    v_right = float(z[i_right] - z[i_min])

    return v_left, v_right


def width_metric(x, z, idxs):
    
    i_left, i_min, i_right = idxs

    z_left = z[i_left]
    z_right = z[i_right]

    # Determine reference side (smaller max) and search side (opposite)
    if z_left <= z_right:
        ref_idx = i_left
        ref_z = float(z_left)
        search_side = 'right'
    else:
        ref_idx = i_right
        ref_z = float(z_right)
        search_side = 'left'

    # Build cubic spline
    cs = CubicSpline(x, z)

    match_x = None

    # Search for zero-crossing of  spline(x) - ref_z  on the correct side
    if search_side == 'right':
        x_lo, x_hi = float(x[i_min]), float(x[-1])
    else:
        x_lo, x_hi = float(x[0]), float(x[i_min])

    # Dense evaluation to find bracket(s) where spline crosses ref_z
    n_eval = max(1000, 10 * len(x))
    x_dense = np.linspace(x_lo, x_hi, n_eval)
    z_dense = cs(x_dense) - ref_z

    # Find sign changes (zero crossings)
    sign_changes = np.where(np.diff(np.sign(z_dense)))[0]

    if len(sign_changes) > 0:
        # Pick the crossing closest to the global minimum
        best = sign_changes[0] if search_side == 'right' else sign_changes[-1]
        # Linear interpolation between the two bracketing points for precision
        x_a, x_b = x_dense[best], x_dense[best + 1]
        z_a, z_b = z_dense[best], z_dense[best + 1]
        if abs(z_b - z_a) > 0:
            match_x = float(x_a - z_a * (x_b - x_a) / (z_b - z_a))
        else:
            match_x = float((x_a + x_b) / 2.0)

    # Fallback: nearest discrete point
    if match_x is None:
        if search_side == 'right':
            search_slice = z[i_min:]
            search_offset = i_min
        else:
            search_slice = z[: i_min + 1]
            search_offset = 0
        diffs = np.abs(search_slice - ref_z)
        best_local = int(np.argmin(diffs))
        fallback_idx = search_offset + best_local
        match_x = float(x[fallback_idx])

    # Width = horizontal distance (the line is horizontal by construction)
    match_z = ref_z
    w = float(abs(match_x - x[ref_idx]))

    return w, ref_idx, (match_x, match_z), cs


def mean_sd(values):

    v = np.array(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan"), 0

    mean = float(np.mean(v))
    sd = float(np.std(v, ddof=1)) if v.size > 1 else 0.0

    return mean, sd, int(v.size)


def safe_stem(filename: str):

    stem = os.path.splitext(filename)[0]
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem)


def format_number(val, decimal_sep='.', precision=':.4g'):
    """Format number with chosen decimal separator"""
    if not np.isfinite(val):
        return "NaN"
    # Extract the format spec from precision string (remove the leading ':')
    fmt_spec = precision.lstrip(':')
    formatted = f"{val:{fmt_spec}}"
    if decimal_sep == ',':
        formatted = formatted.replace('.', ',')
    return formatted


def plot_profile(path, x, z, out_dir, idxs, width_ref_idx,
                 width_match_point, spline, scale=1e9, unit='nm', decimal_sep='.'):

    name = os.path.basename(path)
    stem = safe_stem(name)

    i_left, i_min, i_right = idxs
    w_match_x, w_match_z = width_match_point

    # Convert to display units
    S = scale
    x_d = x * S
    z_d = z * S
    w_match_x_d = w_match_x * S
    w_match_z_d = w_match_z * S

    # Global extrema
    i_max = int(np.argmax(z_d))
    i_min_global = int(np.argmin(z_d))

    x_max, z_max = float(x_d[i_max]), float(z_d[i_max])
    x_min, z_min = float(x_d[i_min_global]), float(z_d[i_min_global])

    v_left = float(z_d[i_left] - z_d[i_min_global])
    v_right = float(z_d[i_right] - z_d[i_min_global])
    h_amp = float(abs(x_d[i_right] - x_d[i_left]))

    # Width metric value (horizontal by construction), in nm
    w_val = float(abs(w_match_x_d - x_d[width_ref_idx]))

    # Calculate plot ranges
    x_span = float(np.max(x_d) - np.min(x_d))
    z_span = float(np.max(z_d) - np.min(z_d))
    if x_span == 0:
        x_span = 1.0
    if z_span == 0:
        z_span = 1.0

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))

    # Plot spline interpolation (smooth curve through data)
    if spline is not None:
        x_fine = np.linspace(float(np.min(x)), float(np.max(x)), 500)
        z_fine = spline(x_fine)
        ax.plot(x_fine * S, z_fine * S, 'r-', linewidth=1.0, alpha=0.5,
                label='Cubic spline', zorder=2)

    # Plot the profile data points
    ax.plot(x_d, z_d, 'k-', linewidth=1.5, label='Profile', zorder=3)

    # Mark GLOBAL extrema
    ax.scatter([x_max], [z_max], s=120, c='red', marker='o',
               edgecolors='darkred', linewidths=2, label='Global Max', zorder=5)
    ax.scatter([x_min], [z_min], s=120, c='blue', marker='o',
               edgecolors='darkblue', linewidths=2, label='Global Min', zorder=5)

    # Mark LOCAL maxima (used for H metric) — only if different from global max
    local_markers_x = []
    local_markers_z = []

    if i_left != i_max:
        local_markers_x.append(x_d[i_left])
        local_markers_z.append(z_d[i_left])
    if i_right != i_max:
        local_markers_x.append(x_d[i_right])
        local_markers_z.append(z_d[i_right])

    if local_markers_x:
        ax.scatter(local_markers_x, local_markers_z, s=80, c='orange', marker='s',
                   edgecolors='darkorange', linewidths=1.5,
                   label='Local Max (H metric)', zorder=4)

    # Mark WIDTH match point (from spline intersection)
    ax.scatter([w_match_x_d], [w_match_z_d], s=80, c='cyan',
               marker='D', edgecolors='darkcyan', linewidths=1.5,
               label='Width match pt (spline)', zorder=6)

    # ---- V_left annotation (left local max → global min) ----
    x_vl = float(x_d[i_left] - 0.06 * x_span)
    ax.annotate(
        "",
        xy=(x_vl, z_d[i_left]),
        xytext=(x_vl, z_min),
        arrowprops=dict(arrowstyle="<->", lw=1.5, color='#CC3366'),
        annotation_clip=False,
    )
    ax.text(
        x_vl - 0.02 * x_span,
        (float(z_d[i_left]) + z_min) / 2.0,
        f"V_L = {format_number(v_left, decimal_sep)} {unit}",
        va="center", ha="right", fontsize=9, fontweight="bold", color='#CC3366',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                  edgecolor='#CC3366', alpha=0.8)
    )

    # ---- V_right annotation (right local max → global min) ----
    x_vr = float(x_d[i_right] + 0.06 * x_span)
    ax.annotate(
        "",
        xy=(x_vr, z_d[i_right]),
        xytext=(x_vr, z_min),
        arrowprops=dict(arrowstyle="<->", lw=1.5, color='#3366CC'),
        annotation_clip=False,
    )
    ax.text(
        x_vr + 0.02 * x_span,
        (float(z_d[i_right]) + z_min) / 2.0,
        f"V_R = {format_number(v_right, decimal_sep)} {unit}",
        va="center", ha="left", fontsize=9, fontweight="bold", color='#3366CC',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                  edgecolor='#3366CC', alpha=0.8)
    )

    # ---- Horizontal amplitude annotation (between local maxima) ----
    y_arrow = float(np.min(z_d) - 0.14 * z_span)
    ax.annotate(
        "",
        xy=(x_d[i_left], y_arrow),
        xytext=(x_d[i_right], y_arrow),
        arrowprops=dict(arrowstyle="<->", lw=2, color='green'),
        annotation_clip=False,
    )
    ax.text(
        (x_d[i_left] + x_d[i_right]) / 2.0,
        y_arrow - 0.05 * z_span,
        f"H = {format_number(h_amp, decimal_sep)} {unit}",
        va="top", ha="center", fontsize=10, fontweight="bold", color='green',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='green', alpha=0.8)
    )

    # Dashed guide lines for H metric
    ax.plot([x_d[i_left], x_d[i_left]], [z_d[i_left], y_arrow],
            linestyle='--', linewidth=1, color='gray', alpha=0.6)
    ax.plot([x_d[i_right], x_d[i_right]], [z_d[i_right], y_arrow],
            linestyle='--', linewidth=1, color='gray', alpha=0.6)

    # ---- Width annotation (line between ref point and spline match point) ----
    ax.plot([x_d[width_ref_idx], w_match_x_d],
            [z_d[width_ref_idx], w_match_z_d],
            linestyle='-', linewidth=2, color='darkcyan', alpha=0.8, zorder=3)
    mid_x_w = (float(x_d[width_ref_idx]) + w_match_x_d) / 2.0
    mid_z_w = (float(z_d[width_ref_idx]) + w_match_z_d) / 2.0
    ax.text(
        mid_x_w,
        mid_z_w + 0.06 * z_span,
        f"W = {format_number(w_val, decimal_sep)} {unit}",
        va="bottom", ha="center", fontsize=9, fontweight="bold", color='darkcyan',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                  edgecolor='darkcyan', alpha=0.8)
    )

    # Labels and title
    ax.set_title(name, fontsize=13, fontweight='bold')
    ax.set_xlabel(f'x [{unit}]', fontsize=11)
    ax.set_ylabel(f'z [{unit}]', fontsize=11)

    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.2), fontsize=9, framealpha=0.9)

    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

    # Expand limits to keep annotations visible
    x_left_pad = max(0.10 * x_span, abs(x_vl - np.min(x_d)) + 0.08 * x_span)
    x_right_pad = max(0.12 * x_span, abs(x_vr - np.max(x_d)) + 0.14 * x_span)
    ax.set_xlim(float(np.min(x_d) - x_left_pad),
                float(np.max(x_d) + x_right_pad))
    ax.set_ylim(float(y_arrow - 0.10 * z_span),
                float(np.max(z_d) + 0.10 * z_span))

    plt.tight_layout()

    out_path = os.path.join(out_dir, f"{stem}_profile.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    return out_path


def save_results_to_csv(per_file_results, summary_stats, output_path, unit='m', 
                        decimal_sep='.', csv_delimiter=';'):

    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile, delimiter=csv_delimiter)

        # Header
        writer.writerow([f'=== PER-FILE RESULTS (values in {unit}) ==='])
        writer.writerow([
            'Filename',
            'V_Left',
            'V_Right',
            'V_High',
            'V_Low',
            'H_Amplitude',
            'Width'
        ])

        for name, vl, vr, vh, vlo, h, w in per_file_results:
            def fmt(val):
                return format_number(val, decimal_sep, precision=':.10g')
            writer.writerow([name, fmt(vl), fmt(vr), fmt(vh), fmt(vlo), fmt(h), fmt(w)])

        writer.writerow([])
        writer.writerow(['=== SUMMARY STATISTICS ==='])
        writer.writerow(['Metric', 'Mean', 'SD', 'N'])
        for key, label in [
            ('vl',   'V_Left'),
            ('vr',   'V_Right'),
            ('vh',   'V_High'),
            ('vlo',  'V_Low'),
            ('h',    'H_Amplitude'),
            ('w',    'Width'),
        ]:
            writer.writerow([
                label,
                format_number(summary_stats[f'{key}_mean'], decimal_sep, precision=':.10g'),
                format_number(summary_stats[f'{key}_sd'], decimal_sep, precision=':.10g'),
                summary_stats[f'{key}_n']
            ])


def main():

    print("=" * 70)
    print("=== AFM Profile Analysis Tool (Enhanced Version) ===")
    print("=" * 70)
    print("\nThis tool calculates:")
    print("  • Vertical amplitude left (V_L): Left local max − Global min")
    print("  • Vertical amplitude right (V_R): Right local max − Global min")
    print("  • V_High / V_Low: Highest and lowest of V_L and V_R per profile")
    print("  • Horizontal amplitude (H): Distance between local maxima")
    print("    on either side of the global minimum")
    print("  • Width (W): Horizontal distance between the smaller local max")
    print("    and the cubic-spline intersection on the opposite side")
    print("    (perfectly horizontal line by construction)")
    print("=" * 70)
    print("\n  NOTE: This script assumes input data is in METERS (SI units),")
    print("  which is the default export format of Gwyddion.\n")

    folder = input("Paste the folder path containing .txt profiles: ").strip().strip('"').strip("'")

    if not os.path.isdir(folder):
        print(f"\nERROR: Not a folder: {folder}")
        return

    files = sorted(glob.glob(os.path.join(folder, "*.txt")))
    if not files:
        print(f"\nNo .txt files found in: {folder}")
        return

    plots_dir = os.path.join(folder, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Decimal delimiter prompt
    print("\n" + "=" * 70)
    print("DECIMAL DELIMITER:")
    print("  A: USA/Canada standard (e.g., 3.14)")
    print("  B: European standard (e.g., 3,14)")
    print("=" * 70)
    delimiter_choice = input("Choose decimal delimiter (A/B) [default: A]: ").strip().upper()

    if delimiter_choice == 'B':
        decimal_sep = ','
        csv_delimiter = ';'
        print("  → Using European format (comma as decimal separator)\n")
    else:
        decimal_sep = '.'
        csv_delimiter = ','
        print("  → Using USA/Canada format (period as decimal separator)\n")

    # Conversion prompt
    convert_input = input("  Convert output to nanometers (check if the input is in meters)? (y/n) [default: y]: ").strip().lower()
    if convert_input in ('n', 'no'):
        scale = 1.0
        unit = 'm'
        print("  → Output will be in METERS.\n")
    else:
        scale = 1e9
        unit = 'nm'
        print("  → Output will be in NANOMETERS.\n")

    per_file = []
    vlamps, vramps, vhamps, vloamps, hamps, wamps = [], [], [], [], [], []
    plot_paths = []

    print(f"\nFound {len(files)} .txt file(s). Processing...\n")

    for i, path in enumerate(files, 1):
        name = os.path.basename(path)
        print(f"  [{i}/{len(files)}] Processing: {name}...", end=" ")
        try:
            x, z = load_profile_txt(path)

            h_width, idxs = horizontal_amplitude_around_min(x, z)
            vl, vr = vertical_amplitudes_left_right(x, z, idxs)
            w, w_ref, w_match_pt, spline = width_metric(x, z, idxs)

            # V_High and V_Low (computed in original meter scale)
            v_high = max(vl, vr)
            v_low = min(vl, vr)

            out_plot = plot_profile(path, x, z, plots_dir, idxs,
                                   w_ref, w_match_pt, spline,
                                   scale=scale, unit=unit, decimal_sep=decimal_sep)

            # Apply scale for CSV and display
            vl_s = vl * scale
            vr_s = vr * scale
            vh_s = v_high * scale
            vlo_s = v_low * scale
            h_s = h_width * scale
            w_s = w * scale

            per_file.append((name, vl_s, vr_s, vh_s, vlo_s, h_s, w_s))
            vlamps.append(vl_s)
            vramps.append(vr_s)
            vhamps.append(vh_s)
            vloamps.append(vlo_s)
            hamps.append(h_s)
            wamps.append(w_s)
            plot_paths.append(out_plot)

            print("-> SUCCESSFUL !")

        except Exception as e:
            per_file.append((name, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))
            print(f"-> ERROR: {e}")

    # Print per-file results
    print(f"\n" + "=" * 120)
    print(f"PER-FILE METRICS ({unit}):")
    print("=" * 120)
    print(f"{'Filename':35s}  {'V_Left':>12s}  {'V_Right':>12s}  {'V_High':>12s}  {'V_Low':>12s}  {'H_amp':>12s}  {'Width':>12s}")
    print("-" * 120)

    for name, vl, vr, vh, vlo, h, w in per_file:
        def fmt(val):
            return format_number(val, decimal_sep, precision=':.6g')
        print(f"{name[:35]:35s}  {fmt(vl):>12s}  {fmt(vr):>12s}  {fmt(vh):>12s}  {fmt(vlo):>12s}  {fmt(h):>12s}  {fmt(w):>12s}")

    # Summary stats
    stats = {}
    for key, vals in [('vl', vlamps), ('vr', vramps),
                       ('vh', vhamps), ('vlo', vloamps),
                       ('h', hamps), ('w', wamps)]:
        m, s, n = mean_sd(vals)
        stats[f'{key}_mean'] = m
        stats[f'{key}_sd'] = s
        stats[f'{key}_n'] = n

    print(f"\n" + "=" * 78)
    print(f"SUMMARY STATISTICS ({unit}, finite values only):")
    print("=" * 78)

    for key, label in [('vl', 'Vertical amplitude left'),
                        ('vr', 'Vertical amplitude right'),
                        ('vh', 'Vertical amplitude HIGH'),
                        ('vlo', 'Vertical amplitude LOW'),
                        ('h', 'Horizontal amplitude'),
                        ('w', 'Width (spline)')]:
        print(f"{label}:")
        print(f"  Mean = {format_number(stats[f'{key}_mean'], decimal_sep, precision=':.6g')}, "
              f"SD = {format_number(stats[f'{key}_sd'], decimal_sep, precision=':.6g')}, "
              f"N = {stats[f'{key}_n']}")

    # Save CSV
    csv_path = os.path.join(folder, "afm_analysis_results.csv")
    save_results_to_csv(per_file, stats, csv_path, unit=unit, 
                        decimal_sep=decimal_sep, csv_delimiter=csv_delimiter)

    print("\n" + "=" * 78)
    print("OUTPUT FILES:")
    print("=" * 78)
    print(f"  • Individual plots: {plots_dir}/")
    print(f"    ({len(plot_paths)} profile plot(s) created)")
    print(f"  • Results CSV: {csv_path}")
    print(f"  • Units: {unit}")
    print(f"  • Decimal separator: '{decimal_sep}'")
    print(f"  • CSV field delimiter: '{csv_delimiter}'")
    print("=" * 78)
    print("\n--> Analysis complete !")


if __name__ == "__main__":
    main()
    