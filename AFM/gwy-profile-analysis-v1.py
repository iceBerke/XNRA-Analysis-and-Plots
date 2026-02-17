# AFM Profile Analysis Tool (raw data from Gwyddion)

# PURPOSE:
# Analyzes AFM (Atomic Force Microscopy) profile data exported from Gwyddion.
# Calculates vertical and horizontal amplitude metrics for surface features.

# METRICS CALCULATED:
# 1. Vertical Amplitude (V):
# - Peak-to-peak distance (global maximum - global minimum)
# 2. Horizontal Amplitude (H):
# - Distance between local maxima on either side of the global minimum
# - Tie-breaking: When multiple points share the max value (plateaus),
# selects the point closest to the minimum to avoid overestimation

# INPUT:
# - Gwyddion-exported files with x and y columns (need to be converted to .txt)
# - Files should be in a single folder (user provides the directory to this folder)

# OUTPUT:
# 1. Individual annotated plots for each profile (saved in the input-directory/plots/ directory)
# - Color-coded markers (red=global max, blue=global min, orange=local max)
# - Visual amplitude annotations
# 2. CSV file with per-file results and summary statistics
# 3. Two separate histogram distributions:
# - Vertical amplitude distribution
# - Horizontal amplitude distribution
# (Both include mean and ±1 SD markers)

# Author: Berke Santos & Giuseppe Legrotaglie
# Developed with the help of Claude.AI and ChatGPT version 5.2
# Created: 16/02/2026

import os
import re
import glob
import csv
import numpy as np
import matplotlib.pyplot as plt


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
    #print(x[order])
    #print(z[order])
    return x[order], z[order]


def vertical_amplitude_ptp(z: np.ndarray):
    
    return abs(float(np.max(z) - np.min(z)))


def horizontal_amplitude_around_min(x: np.ndarray, y: np.ndarray):
    
    # Find global minimum
    i_min = int(np.argmin(y))
    #print(i_min)

    # Find local maximum on the left side (use closest to i_min if there are ties)
    if i_min > 0:
        left_section = y[: i_min + 1]
        max_val_left = np.max(left_section)
        # Find all indices with this maximum value
        candidates_left = np.where(left_section == max_val_left)[0]
        # Pick the rightmost (closest to i_min)
        i_left = int(candidates_left[-1])
    else:
        i_left = 0

    # Find local maximum on the right side (use closest to i_min if there are ties)
    if i_min < len(y) - 1:
        right_section = y[i_min:]
        max_val_right = np.max(right_section)
        # Find all indices with this maximum value (relative to right_section)
        candidates_right = np.where(right_section == max_val_right)[0]
        # Pick the leftmost (closest to i_min), then offset by i_min
        i_right = int(i_min + candidates_right[0])
    else:
        i_right = len(y) - 1

    width = float(abs(x[i_right] - x[i_left]))
    return width, (i_left, i_min, i_right)


def mean_sd(values):
    
    v = np.array(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan"), 0
    
    mean = float(np.mean(v))
    sd = float(np.std(v, ddof=1)) if v.size > 1 else 0.0

    return mean, sd, int(v.size)


def plot_histogram(values, stats, metric_name, color, output_path):
    
    # Filter out NaN values
    finite_vals = np.array([v for v in values if np.isfinite(v)])
    
    if len(finite_vals) == 0:
        print(f"  Warning: No finite values for {metric_name}")
        return None
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Use 'auto' binning, but ensure at least 5 bins if we have enough data
    n_bins = 'auto' if len(finite_vals) > 10 else max(5, len(finite_vals) // 2)
    
    n, bins, patches = ax.hist(finite_vals, bins=n_bins, color=color, 
                                edgecolor='black', alpha=0.7, linewidth=1.2)
    
    # Add mean line
    ax.axvline(stats['mean'], color='red', linestyle='--', linewidth=2.5, 
                label=f"Mean = {stats['mean']:.6g}", zorder=5)
    
    # Add ±1 SD lines if SD exists
    if stats['sd'] > 0:
        ax.axvline(stats['mean'] - stats['sd'], color='orange', 
                   linestyle=':', linewidth=2, alpha=0.8,
                   label=f"±1 SD = {stats['sd']:.6g}", zorder=5)
        ax.axvline(stats['mean'] + stats['sd'], color='orange', 
                   linestyle=':', linewidth=2, alpha=0.8, zorder=5)
    
    ax.set_xlabel(f'{metric_name} [m]', fontsize=12, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax.set_title(f'{metric_name} Distribution (N={stats["n"]})', 
                 fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    return output_path


def safe_stem(filename: str):
    
    stem = os.path.splitext(filename)[0]
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem)


def plot_profile(path: str, x: np.ndarray, z: np.ndarray, out_dir: str, idxs):
    
    name = os.path.basename(path)
    stem = safe_stem(name)

    i_left, i_min, i_right = idxs

    # Global extrema
    i_max = int(np.argmax(z))
    i_min_global = int(np.argmin(z))

    x_max, z_max = float(x[i_max]), float(z[i_max])
    x_min, z_min = float(x[i_min_global]), float(z[i_min_global])

    v_amp = z_max - z_min
    h_amp = float(abs(x[i_right] - x[i_left]))

    # Calculate plot ranges
    x_span = float(np.max(x) - np.min(x))
    z_span = float(np.max(z) - np.min(z))
    if x_span == 0:
        x_span = 1.0
    if z_span == 0:
        z_span = 1.0

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot the profile
    ax.plot(x, z, 'k-', linewidth=1.5, label='Profile', zorder=1)

    # Mark GLOBAL extrema with larger markers
    ax.scatter([x_max], [z_max], s=120, c='red', marker='o', 
               edgecolors='darkred', linewidths=2, label='Global Max', zorder=5)
    ax.scatter([x_min], [z_min], s=120, c='blue', marker='o', 
               edgecolors='darkblue', linewidths=2, label='Global Min', zorder=5)

    # Mark LOCAL maxima (used for horizontal metric) with smaller, different color
    # Only plot them if they're different from global max
    local_markers_x = []
    local_markers_z = []
    
    if i_left != i_max:
        local_markers_x.append(x[i_left])
        local_markers_z.append(z[i_left])
    if i_right != i_max:
        local_markers_x.append(x[i_right])
        local_markers_z.append(z[i_right])
    
    if local_markers_x:
        ax.scatter(local_markers_x, local_markers_z, s=80, c='orange', marker='s',
                   edgecolors='darkorange', linewidths=1.5, label='Local Max (H metric)', zorder=4)

    # --- Vertical amplitude annotation (global max-min) ---
    x_bracket = float(np.max(x) + 0.10 * x_span)
    ax.annotate(
        "",
        xy=(x_bracket, z_max),
        xytext=(x_bracket, z_min),
        arrowprops=dict(arrowstyle="<->", lw=2, color='purple'),
        annotation_clip=False,
    )
    ax.text(
        x_bracket + 0.02 * x_span,
        (z_max + z_min) / 2.0,
        f"V = {v_amp:.6g}",
        va="center",
        ha="left",
        fontsize=11,
        fontweight="bold",
        color='purple',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='purple', alpha=0.8)
    )

    # --- Horizontal amplitude annotation (between local maxima) ---
    y_arrow = float(np.min(z) - 0.14 * z_span)
    ax.annotate(
        "",
        xy=(x[i_left], y_arrow),
        xytext=(x[i_right], y_arrow),
        arrowprops=dict(arrowstyle="<->", lw=2, color='green'),
        annotation_clip=False,
    )
    ax.text(
        (x[i_left] + x[i_right]) / 2.0,
        y_arrow - 0.05 * z_span,
        f"H = {h_amp:.6g}",
        va="top",
        ha="center",
        fontsize=11,
        fontweight="bold",
        color='green',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='green', alpha=0.8)
    )

    # Dashed guide lines for horizontal metric
    ax.plot([x[i_left], x[i_left]], [z[i_left], y_arrow], 
            linestyle='--', linewidth=1, color='gray', alpha=0.6)
    ax.plot([x[i_right], x[i_right]], [z[i_right], y_arrow], 
            linestyle='--', linewidth=1, color='gray', alpha=0.6)

    # Labels and title
    ax.set_title(name, fontsize=13, fontweight='bold')
    ax.set_xlabel('x [m]', fontsize=11)
    ax.set_ylabel('z [m]', fontsize=11)
    
    # Legend positioned outside the plot area (to the right)
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.2), fontsize=9, framealpha=0.9)
    
    # Grid for easier reading
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)

    # Expand limits to keep annotations visible
    ax.set_xlim(float(np.min(x) - 0.02 * x_span), float(x_bracket + 0.12 * x_span))
    ax.set_ylim(float(y_arrow - 0.10 * z_span), float(np.max(z) + 0.10 * z_span))

    plt.tight_layout()

    out_path = os.path.join(out_dir, f"{stem}_profile.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    return v_amp, h_amp, out_path


def save_results_to_csv(per_file_results, summary_stats, output_path):
    
    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        writer.writerow(['=== PER-FILE RESULTS ==='])
        writer.writerow(['Filename', 'Vertical_Amplitude_ptp', 'Horizontal_Amplitude_around_min'])
        
        # Write per-file data
        for name, v, h in per_file_results:
            v_str = f"{v:.10g}" if np.isfinite(v) else "NaN"
            h_str = f"{h:.10g}" if np.isfinite(h) else "NaN"
            writer.writerow([name, v_str, h_str])
        
        # Write summary
        writer.writerow([])
        writer.writerow(['=== SUMMARY STATISTICS ==='])
        writer.writerow(['Metric', 'Mean', 'SD', 'N'])
        writer.writerow(['Vertical_Amplitude_ptp', 
                        f"{summary_stats['v_mean']:.10g}", 
                        f"{summary_stats['v_sd']:.10g}", 
                        summary_stats['v_n']])
        writer.writerow(['Horizontal_Amplitude_around_min', 
                        f"{summary_stats['h_mean']:.10g}", 
                        f"{summary_stats['h_sd']:.10g}", 
                        summary_stats['h_n']])


def main():
    print("=" * 70)
    print("=== AFM Profile Analysis Tool (Enhanced Version) ===")
    print("=" * 70)
    print("\nThis tool calculates:")
    print("  • Vertical amplitude: Global max - Global min")
    print("  • Horizontal amplitude: Distance between local maxima")
    print("    on either side of the global minimum")
    print("=" * 70)
    
    folder = input("\nPaste the folder path containing .txt profiles: ").strip().strip('"').strip("'")

    if not os.path.isdir(folder):
        print(f"\nERROR: Not a folder: {folder}")
        return

    files = sorted(glob.glob(os.path.join(folder, "*.txt")))
    if not files:
        print(f"\nNo .txt files found in: {folder}")
        return

    plots_dir = os.path.join(folder, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    per_file = []
    vamps = []
    hamps = []
    plot_paths = []

    print(f"\nFound {len(files)} .txt file(s). Processing...\n")

    for i, path in enumerate(files, 1):
        name = os.path.basename(path)
        print(f"  [{i}/{len(files)}] Processing: {name}...", end=" ")
        try:
            x, z = load_profile_txt(path)

            v = vertical_amplitude_ptp(z)
            h_width, idxs = horizontal_amplitude_around_min(x, z)

            # Save annotated plot
            v2, h2, out_plot = plot_profile(path, x, z, plots_dir, idxs)

            per_file.append((name, v, h_width))
            vamps.append(v)
            hamps.append(h_width)
            plot_paths.append(out_plot)
            
            print("✓")

        except Exception as e:
            per_file.append((name, np.nan, np.nan))
            print(f"✗ ERROR: {e}")

    # Print per-file results
    print("\n" + "=" * 78)
    print("PER-FILE METRICS:")
    print("=" * 78)
    print(f"{'Filename':40s}  {'V_amp (ptp)':>16s}  {'H_amp (local max)':>18s}")
    print("-" * 78)
    for name, v, h in per_file:
        v_str = f"{v:.6g}" if np.isfinite(v) else "NaN"
        h_str = f"{h:.6g}" if np.isfinite(h) else "NaN"
        print(f"{name[:40]:40s}  {v_str:>16s}  {h_str:>18s}")

    # Summary stats
    v_mean, v_sd, v_n = mean_sd(vamps)
    h_mean, h_sd, h_n = mean_sd(hamps)

    summary = {
        'v_mean': v_mean, 'v_sd': v_sd, 'v_n': v_n,
        'h_mean': h_mean, 'h_sd': h_sd, 'h_n': h_n
    }

    print("\n" + "=" * 78)
    print("SUMMARY STATISTICS (finite values only):")
    print("=" * 78)
    print(f"Vertical amplitude (ptp):")
    print(f"  Mean = {v_mean:.6g}, SD = {v_sd:.6g}, N = {v_n}")
    print(f"\nHorizontal amplitude (local max distance):")
    print(f"  Mean = {h_mean:.6g}, SD = {h_sd:.6g}, N = {h_n}")

    # Generate histograms (separate files for each metric)
    print("\nGenerating distribution histograms...")
    v_stats = {'mean': v_mean, 'sd': v_sd, 'n': v_n}
    h_stats = {'mean': h_mean, 'sd': h_sd, 'n': h_n}
    
    v_hist_path = os.path.join(folder, "vertical_amplitude_distribution.png")
    h_hist_path = os.path.join(folder, "horizontal_amplitude_distribution.png")
    
    v_hist_saved = plot_histogram(vamps, v_stats, 'Vertical Amplitude (ptp)', 
                                   'steelblue', v_hist_path)
    h_hist_saved = plot_histogram(hamps, h_stats, 'Horizontal Amplitude (local max)', 
                                   'seagreen', h_hist_path)

    # Save results to CSV
    csv_path = os.path.join(folder, "afm_analysis_results.csv")
    save_results_to_csv(per_file, summary, csv_path)
    
    print("\n" + "=" * 78)
    print("OUTPUT FILES:")
    print("=" * 78)
    print(f"  • Individual plots: {plots_dir}/")
    print(f"    ({len(plot_paths)} profile plot(s) created)")
    print(f"  • Distribution histograms:")
    if v_hist_saved:
        print(f"    - Vertical: {v_hist_path}")
    if h_hist_saved:
        print(f"    - Horizontal: {h_hist_path}")
    print(f"  • Results CSV: {csv_path}")
    print("=" * 78)
    print("\n✓ Analysis complete!")


if __name__ == "__main__":
    main()
