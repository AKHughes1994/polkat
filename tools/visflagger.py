#!/usr/bin/env python3
"""
Visibility data flagger for radio interferometry observations.
Performs outlier detection on amplitude vs baseline distributions.
"""

import warnings
import os

# Suppress all warnings before any other imports
warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'

import numpy as np
import matplotlib.pyplot as plt
import datashader as dsh
import datashader.transfer_functions as tf
import pandas as pd
import dask
from daskms import xds_from_ms, xds_from_table
from datetime import datetime
import logging
import bottleneck as bn
import argparse
import psutil
import sys

# Setup file logging (will be reconfigured in main() with output suffix)
# logging.basicConfig(
#     filename='baseline-flagger.log',
#     level=logging.INFO,
#     format='[%(asctime)s] %(message)s',
#     datefmt='%Y-%m-%d %H:%M:%S',
#     filemode='w'  # Overwrite log file each run
# )

# Configuration defaults (overridden by command-line arguments)
def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Visibility data flagger for radio interferometry observations.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Mandatory arguments
    parser.add_argument('ms_file', type=str,
                        help='Path to measurement set file')
    parser.add_argument('field', type=str,
                        help='Field ID (integer) or field name (string) to process')
    
    # Optional arguments with defaults
    parser.add_argument('--data-column', type=str, default='CORRECTED_DATA',
                        help='Data column to read')
    parser.add_argument('--channel-group-size', type=int, default=32,
                        help='Number of channels per chunk')
    parser.add_argument('--freq-target', type=float, default=1.0e9,
                        help='Target frequency in Hz (starting point for chunk search)')
    parser.add_argument('--save-dir', type=str, default='VIS_plots/',
                        help='Directory to save plots')
    parser.add_argument('--output-suffix', type=str, default='',
                        help='Suffix for output files (e.g., "_test" -> baseline_flags_test.txt)')
    parser.add_argument('--sigma-threshold', type=float, default=4.0,
                        help='Sigma threshold for outlier detection')
    parser.add_argument('--outlier-mode', type=str, default='high',
                        choices=['high', 'low', 'both', 'mixed'],
                        help='Outlier detection mode: high, low, both, or mixed (both for XX/YY, high for XY/YX)')
    parser.add_argument('--outlier-fraction-threshold', type=float, default=0.1,
                        help='Flag baselines with >X fraction outliers')
    parser.add_argument('--correlation-products', type=str, default='XX,YY,XY,YX',
                        help='Correlation products to process (comma-separated or list)')
    parser.add_argument('--antenna-flag-threshold', type=float, default=0.90,
                        help='Flag entire antenna if >X fraction of its baselines are bad')
    parser.add_argument('--antenna-relative-sigma', type=float, default=30.0,
                        help='Sigma threshold for relative antenna outlier detection')
    parser.add_argument('--antenna-flag-cap', type=int, default=3,
                        help='Maximum number of antennas to flag globally')
    parser.add_argument('--num-chunks', type=int, default=None,
                        help='Number of chunks to process (None = all chunks)')
    parser.add_argument('--no-printing', action='store_true',
                        help='Disable console printing (logging still active)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    
    return parser.parse_args()


def check_system_resources():
    """Check available system resources (non-swap memory and CPU cores)."""
    mem = psutil.virtual_memory()
    
    resources = {
        'available_memory_gb': mem.available / (1024**3),
        'total_memory_gb': mem.total / (1024**3),
        'memory_percent_used': mem.percent,
        'physical_cores': psutil.cpu_count(logical=False),
        'logical_cores': psutil.cpu_count(logical=True),
    }
    
    return resources


def estimate_chunk_memory(n_baselines, n_time, n_freq_per_chunk, n_corr):
    """
    Estimate total memory required to process one chunk.
    
    Returns
    -------
    float
        Estimated memory in GB
    """
    # All baselines loaded into memory
    vis_size = n_baselines * n_time * n_freq_per_chunk * n_corr * 8
    weight_size = n_baselines * n_time * n_corr * 4
    flag_size = n_baselines * n_time * n_freq_per_chunk * n_corr * 1
    
    # DataFrame construction overhead (~2.5x the raw data)
    raw_data_gb = (vis_size + weight_size + flag_size) / (1024**3)
    with_overhead_gb = raw_data_gb * 2.5
    
    return with_overhead_gb


def check_chunk_memory_safe(chunk_memory_gb, available_memory_gb):
    """
    Check if chunk will fit safely in memory (85% limit).
    
    Returns
    -------
    tuple
        (is_safe, max_safe_gb, suggested_chunk_size)
    """
    max_safe_gb = available_memory_gb * 0.85
    is_safe = chunk_memory_gb <= max_safe_gb
    
    # Calculate suggested chunk size if not safe
    if not is_safe:
        reduction_factor = chunk_memory_gb / max_safe_gb
        # This gets filled in by caller
        suggested_chunk_size = None
    else:
        suggested_chunk_size = None
    
    return is_safe, max_safe_gb, suggested_chunk_size


def check_memory_critical():
    """
    Check if memory usage is critical. Raise error if >95%.
    
    Returns
    -------
    float
        Current memory usage percentage
    """
    mem = psutil.virtual_memory()
    percent_used = mem.percent
    
    if percent_used > 95:
        raise MemoryError(
            f"CRITICAL: Memory usage at {percent_used:.1f}% - stopping to prevent crash. "
            f"Available: {mem.available / (1024**3):.2f} GB"
        )
    
    if percent_used > 85:
        log(f"WARNING: Memory usage at {percent_used:.1f}% - approaching limit")
    
    return percent_used


# Global variable for controlling console output (set from args in main())
ENABLE_PRINTING = True


def log(message=""):
    """Print message with timestamp if ENABLE_PRINTING is True, always write to log file."""
    # Always write to log file
    if message == "":
        logging.info("")
    else:
        logging.info(message)
    
    # Conditionally print to console
    if ENABLE_PRINTING:
        if message == "":
            print()
        else:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{timestamp}] {message}")


def parse_correlation_products(corr_input):
    """
    Parse correlation product input and return list of (name, index) tuples.
    
    Mapping: XX=0, XY=1, YX=2, YY=3
    
    Parameters
    ----------
    corr_input : str or list
        Either 'XX,YY,XY,YX' or ['XX', 'YY'] or ['XX,YY']
    
    Returns
    -------
    list of tuples
        [(name, index), ...] e.g., [('XX', 0), ('YY', 3)]
    """
    corr_map = {'XX': 0, 'XY': 1, 'YX': 2, 'YY': 3}
    
    # Parse input
    if isinstance(corr_input, str):
        corr_names = [c.strip().upper() for c in corr_input.split(',')]
    elif isinstance(corr_input, list):
        # Flatten in case of ['XX,YY'] format
        corr_names = []
        for item in corr_input:
            if ',' in item:
                corr_names.extend([c.strip().upper() for c in item.split(',')])
            else:
                corr_names.append(item.strip().upper())
    else:
        raise ValueError(f"CORRELATION_PRODUCTS must be string or list, got {type(corr_input)}")
    
    # Convert to indices
    result = []
    for name in corr_names:
        if name not in corr_map:
            raise ValueError(f"Unknown correlation product: {name}. Must be one of {list(corr_map.keys())}")
        result.append((name, corr_map[name]))
    
    return result


def weighted_median(values, weights):
    """
    Calculate weighted median.
    
    Parameters
    ----------
    values : array-like
        Data values
    weights : array-like
        Weights for each value
        
    Returns
    -------
    float
        Weighted median
    """
    # Remove NaN values
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values_clean = values[mask]
    weights_clean = weights[mask]
    
    if len(values_clean) == 0:
        return np.nan
    
    # Sort by values
    sorted_indices = np.argsort(values_clean)
    sorted_values = values_clean[sorted_indices]
    sorted_weights = weights_clean[sorted_indices]
    
    # Cumulative sum of weights
    cumsum_weights = np.cumsum(sorted_weights)
    total_weight = cumsum_weights[-1]
    
    # Find where cumulative weight crosses 50%
    midpoint = total_weight / 2.0
    
    # Find the weighted median
    idx = np.searchsorted(cumsum_weights, midpoint)
    
    return sorted_values[idx]


def load_antenna_names(ms_file):
    """Load antenna names from the ANTENNA subtable."""
    ant_table = xds_from_table(f"{ms_file}::ANTENNA", columns=["NAME"])[0]
    names = ant_table.NAME.compute().values
    # Handle bytes vs unicode
    return np.array([n.decode() if isinstance(n, (bytes, np.bytes_)) else str(n) for n in names])


def load_chan_freq(ms_file):
    """Load channel frequencies from SPECTRAL_WINDOW subtable in Hz."""
    spw_table = xds_from_table(f"{ms_file}::SPECTRAL_WINDOW", columns=["CHAN_FREQ"])[0]
    chan_freq = spw_table.CHAN_FREQ.compute().values
    # Return first SPW (flatten if multiple SPWs)
    if chan_freq.ndim > 1:
        return chan_freq[0]
    return chan_freq


def load_field_names(ms_file):
    """Load field names from FIELD subtable."""
    field_table = xds_from_table(f"{ms_file}::FIELD", columns=["NAME"])[0]
    names = field_table.NAME.compute().values
    return names


def resolve_field(ms_file, field_arg):
    """
    Resolve field argument to field ID.
    
    Parameters
    ----------
    ms_file : str
        Path to measurement set
    field_arg : str
        Field argument (integer ID or field name)
    
    Returns
    -------
    int
        Field ID
    str
        Field name
    
    Raises
    ------
    ValueError
        If field not found
    """
    field_names = load_field_names(ms_file)
    
    # Try to parse as integer first
    try:
        field_id = int(field_arg)
        # Validate that field ID exists
        if field_id < 0 or field_id >= len(field_names):
            raise ValueError(
                f"Field ID {field_id} out of range. Available fields (0-{len(field_names)-1}): "
                f"{', '.join([f'{i}={name}' for i, name in enumerate(field_names)])}"
            )
        return field_id, field_names[field_id]
    except ValueError:
        # It's a string, look up field name (case-insensitive)
        field_arg_lower = field_arg.lower().strip()
        for i, name in enumerate(field_names):
            if name.lower().strip() == field_arg_lower:
                return i, name
        # Not found
        raise ValueError(
            f"Field '{field_arg}' not found. Available fields: "
            f"{', '.join([f'{i}={name}' for i, name in enumerate(field_names)])}"
        )


def load_grouped_ms(ms_file, field_id, data_column):
    """
    Load MS grouped by baseline (ANTENNA1, ANTENNA2, DATA_DESC_ID).
    Returns a list of datasets, one per baseline.
    """
    cols = [data_column, "FLAG", "WEIGHT", "ANTENNA1", "ANTENNA2"]
    taql = f"FIELD_ID=={field_id}"
    
    ds_list = xds_from_ms(
        ms_file,
        columns=cols,
        taql_where=taql,
        group_cols=["ANTENNA1", "ANTENNA2", "DATA_DESC_ID"],
    )
    
    # Rename data column to DATA for consistent access
    ds_list = [ds.rename({data_column: "DATA"}) for ds in ds_list]
    
    return ds_list


def find_first_unflagged_chunk(ds_list, chunk_size, start_channel=0):
    """
    Find the first channel chunk with unflagged data across all baselines.
    Returns start channel index.
    Computes a per-channel flag profile once, then analyzes chunks.
    """
    n_chan = ds_list[0].dims['chan']
    
    # Build lazy tasks: sum flags per channel across all baselines
    flag_per_chan_tasks = []
    total_elements_per_chan = 0
    
    for ds in ds_list:
        # Sum over row and corr dimensions, keeping chan dimension
        flag_per_chan_tasks.append(ds.FLAG.sum(dim=['row', 'corr']).data)
        # Count total elements per channel for this baseline
        total_elements_per_chan += ds.FLAG.shape[0] * ds.FLAG.shape[2]  # n_row * n_corr
    
    # Single compute call
    log("Computing channel flag profile...")
    flag_per_chan_arrays = dask.compute(*flag_per_chan_tasks)
    
    # Sum across all baselines to get total flags per channel
    flag_profile = sum(flag_per_chan_arrays)  # Shape: (n_chan,)
    
    # Find first unflagged chunk starting from start_channel
    for start_chan in range(start_channel, n_chan, chunk_size):
        end_chan = min(start_chan + chunk_size, n_chan)
        
        # Sum flags in this chunk
        chunk_flags = flag_profile[start_chan:end_chan].sum()
        chunk_total = total_elements_per_chan * (end_chan - start_chan)
        flag_fraction = chunk_flags / chunk_total if chunk_total > 0 else 1.0
        
        # DEBUG: Print flag fraction for this chunk
        log(f"  Channels {start_chan}-{end_chan}: flag fraction = {flag_fraction:.4f}")
        
        # If not fully flagged, use this chunk
        if chunk_flags < chunk_total:
            return start_chan
    
    # If no unflagged chunks found, return start_channel
    log(f"  Warning: No unflagged chunks found from channel {start_channel}, using that channel anyway")
    return start_channel


def create_histogram_plot(df, median_amp, mad_amp, threshold_mad, chunk_num, corr_name, 
                         outlier_mode, save_dir):
    """
    Create histogram of amplitude distribution showing median, MAD, and thresholds.
    """
    log("Generating histogram diagnostic plot...")
    
    amp_values = df['amplitude'].values
    amp_values_clean = amp_values[np.isfinite(amp_values)]
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Create histogram with 200 bins
    counts, bins, patches = ax.hist(amp_values_clean, bins=200, color='skyblue', 
                                      edgecolor='black', linewidth=0.5, alpha=0.7)
    
    # Use log scale for y-axis to see tails
    ax.set_yscale('log')
    
    # Add weighted median line
    ax.axvline(median_amp, color='red', linestyle='-', linewidth=2.5, 
               label=f'Weighted Median: {median_amp:.4f}', alpha=0.9, zorder=5)
    
    # Add MAD region
    ax.axvspan(median_amp - mad_amp, median_amp + mad_amp, 
               color='red', alpha=0.15, label=f'MAD region: ±{mad_amp:.4f}', zorder=3)
    
    # Add threshold lines based on outlier_mode
    if outlier_mode == 'high':
        ax.axvline(median_amp + threshold_mad, color='orange', linestyle='--', 
                   linewidth=2, label=f'High threshold: {median_amp + threshold_mad:.4f}', 
                   alpha=0.8, zorder=4)
    elif outlier_mode == 'low':
        ax.axvline(median_amp - threshold_mad, color='blue', linestyle='--', 
                   linewidth=2, label=f'Low threshold: {median_amp - threshold_mad:.4f}', 
                   alpha=0.8, zorder=4)
    else:  # 'both'
        ax.axvline(median_amp + threshold_mad, color='orange', linestyle='--', 
                   linewidth=2, label=f'Upper threshold: {median_amp + threshold_mad:.4f}', 
                   alpha=0.8, zorder=4)
        ax.axvline(median_amp - threshold_mad, color='blue', linestyle='--', 
                   linewidth=2, label=f'Lower threshold: {median_amp - threshold_mad:.4f}', 
                   alpha=0.8, zorder=4)
    
    ax.set_xlabel('Amplitude', fontsize=12)
    ax.set_ylabel('Count (log scale)', fontsize=12)
    ax.set_title(f'Amplitude Distribution - Chunk {chunk_num} - {corr_name}', fontsize=14)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save histogram
    os.makedirs(save_dir, exist_ok=True)
    output_path = os.path.join(save_dir, f'amp_histogram_chunk{chunk_num}_{corr_name}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    log(f"Histogram saved to {output_path}")
    plt.close()


def create_plot(df, baseline_names, median_amp, mad_amp, threshold_mad, 
                outlier_baselines, start_chan, end_chan, chunk_num, corr_name,
                outlier_fraction_threshold, outlier_mode, save_dir, field):
    """
    Create amplitude vs baseline plot with datashader and matplotlib.
    Two subplots vertically: upper colored by antenna, lower colored by weight.
    Marks baselines exceeding outlier_fraction_threshold with red stars.
    """
    log("Generating plot with datashader...")
    
    n_baselines = len(baseline_names)
    n_antennas = len(df['ant1'].cat.categories)
    
    # Generate distinct colors for all antennas using a colormap
    cmap = plt.cm.get_cmap('tab20')
    if n_antennas > 20:
        cmap = plt.cm.get_cmap('hsv')
    
    # Create color list with enough colors
    color_list = [cmap(i / n_antennas) for i in range(n_antennas)]
    color_key_ant = ['#%02x%02x%02x' % tuple(int(c * 255) for c in color[:3]) for color in color_list]
    
    # Use divergent colormap for weights (blue = high weight, red = low weight)
    weight_cmap = plt.cm.get_cmap('RdBu')
    
    # Create figure with two subplots vertically sharing x-axis
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), sharex=True)
    
    # Larger points: increase canvas resolution significantly
    canvas_width = 2000
    canvas_height = 800
    
    # Get amplitude range
    amp_min = df['amplitude'].min()
    amp_max = df['amplitude'].max()
    
    # ========== UPPER PLOT: Colored by antenna ==========
    canvas1 = dsh.Canvas(plot_width=canvas_width, plot_height=canvas_height,
                        x_range=(-0.5, n_baselines - 0.5),
                        y_range=(amp_min, amp_max))
    
    agg1 = canvas1.points(df, 'baseline_idx', 'amplitude', dsh.count_cat('ant1'))
    img1 = tf.shade(agg1, color_key=color_key_ant, how='linear', min_alpha=60)
    img1 = tf.dynspread(img1, threshold=0.5, max_px=4)
    img1 = tf.set_background(img1, 'white')
    
    ax1.imshow(img1.to_pil(), extent=[-0.5, n_baselines - 0.5, amp_min, amp_max],
               aspect='auto', origin='upper')
    
    # Identify baselines exceeding outlier fraction threshold
    bad_baselines = outlier_baselines[
        outlier_baselines['outlier_fraction'] > outlier_fraction_threshold
    ]
    
    if len(bad_baselines) > 0:
        n_total_baselines = len(baseline_names)
        baseline_loss_fraction = len(bad_baselines) / n_total_baselines
        log(f"Marking {len(bad_baselines)} baselines exceeding {outlier_fraction_threshold*100:.0f}% outlier threshold")
        log(f"Flagging these baselines would remove {baseline_loss_fraction:.4f} ({len(bad_baselines)}/{n_total_baselines}) of all baselines")
        
        # Mark 99.9th percentile max for each bad baseline with red stars (top outliers only)
        for bl_idx, row in bad_baselines.iterrows():
            bl_data = df[df['baseline_idx'] == bl_idx]
            if len(bl_data) > 0:
                amp_max = bl_data['amplitude'].quantile(0.999)
                ax1.plot(bl_idx, amp_max, 'r*', markersize=8, markeredgecolor='black', 
                        markeredgewidth=0.5)
                ax2.plot(bl_idx, amp_max, 'r*', markersize=8, markeredgecolor='black', 
                        markeredgewidth=0.5)
    
    # Add median line and MAD region to upper plot
    ax1.axhline(median_amp, color='red', linestyle='-', linewidth=2, 
                label=f'Weighted Median: {median_amp:.4f}', alpha=0.8)
    ax1.axhspan(median_amp - mad_amp, median_amp + mad_amp, 
                color='red', alpha=0.2, label=f'MAD: {mad_amp:.4f}')
    
    # Add threshold lines to upper plot
    if outlier_mode == 'high':
        ax1.axhline(median_amp + threshold_mad, color='orange', linestyle='--', 
                    linewidth=2, label=f'High threshold', alpha=0.8)
    elif outlier_mode == 'low':
        ax1.axhline(median_amp - threshold_mad, color='blue', linestyle='--', 
                    linewidth=2, label=f'Low threshold', alpha=0.8)
    else:  # 'both'
        ax1.axhline(median_amp + threshold_mad, color='orange', linestyle='--', 
                    linewidth=2, label=f'Upper threshold', alpha=0.8)
        ax1.axhline(median_amp - threshold_mad, color='blue', linestyle='--', 
                    linewidth=2, label=f'Lower threshold', alpha=0.8)
    
    ax1.set_ylabel('Amplitude', fontsize=12)
    ax1.set_title(f'Amplitude vs Baseline (colored by antenna) - Field {field}, Channels {start_chan}-{end_chan}, {corr_name}', 
                  fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=9)
    
    # ========== LOWER PLOT: Colored by weight ==========
    canvas2 = dsh.Canvas(plot_width=canvas_width, plot_height=canvas_height,
                        x_range=(-0.5, n_baselines - 0.5),
                        y_range=(df['amplitude'].min(), df['amplitude'].max()))
    
    # Use mean aggregation for continuous weight values
    agg2 = canvas2.points(df, 'baseline_idx', 'amplitude', dsh.mean('weight'))
    
    # Use span to set color range for weights
    weight_span = (df['weight'].quantile(0.01), df['weight'].quantile(0.99))
    img2 = tf.shade(agg2, cmap=weight_cmap, how='linear', span=weight_span, min_alpha=60)
    img2 = tf.dynspread(img2, threshold=0.5, max_px=4)
    img2 = tf.set_background(img2, 'white')
    
    ax2.imshow(img2.to_pil(), extent=[-0.5, n_baselines - 0.5, 
                                       df['amplitude'].min(), df['amplitude'].max()],
               aspect='auto', origin='upper')
    
    # Add median line and MAD region to lower plot
    ax2.axhline(median_amp, color='red', linestyle='-', linewidth=2, alpha=0.8)
    ax2.axhspan(median_amp - mad_amp, median_amp + mad_amp, 
                color='red', alpha=0.2)
    
    # Add threshold lines to lower plot
    if outlier_mode == 'high':
        ax2.axhline(median_amp + threshold_mad, color='orange', linestyle='--', 
                    linewidth=2, alpha=0.8)
    elif outlier_mode == 'low':
        ax2.axhline(median_amp - threshold_mad, color='blue', linestyle='--', 
                    linewidth=2, alpha=0.8)
    else:  # 'both'
        ax2.axhline(median_amp + threshold_mad, color='orange', linestyle='--', 
                    linewidth=2, alpha=0.8)
        ax2.axhline(median_amp - threshold_mad, color='blue', linestyle='--', 
                    linewidth=2, alpha=0.8)
    
    ax2.set_xlabel('Baseline Index', fontsize=12)
    ax2.set_ylabel('Amplitude', fontsize=12)
    ax2.set_title(f'Amplitude vs Baseline (colored by weight - blue is high weight)', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Create output directory and save
    os.makedirs(save_dir, exist_ok=True)
    output_path = os.path.join(save_dir, f'amp_vs_baseline_chunk{chunk_num}_{corr_name}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    log(f"Plot saved to {output_path}")
    plt.close()


def main():
    """Main processing routine."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Setup logging with output suffix
    log_filename = f'baseline-flagger{args.output_suffix}.log'
    logging.basicConfig(
        filename=log_filename,
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filemode='w',  # Overwrite log file each run
        force=True  # Reconfigure if already configured
    )
    
    # Set global variable for console output
    global ENABLE_PRINTING
    ENABLE_PRINTING = not args.no_printing
    
    # Resolve field argument to field ID and name
    field_id, field_name = resolve_field(args.ms_file, args.field)
    
    # =================================================================
    # CONFIGURATION SUMMARY
    # =================================================================
    log("="*70)
    log("VISFLAGGER CONFIGURATION")
    log("="*70)
    log("")
    
    log("INPUT DATA:")
    log(f"  MS file: {args.ms_file}")
    log(f"  Data column: {args.data_column}")
    log(f"  Field: {field_name} (ID: {field_id})")
    log("")
    
    log("PROCESSING PARAMETERS:")
    log(f"  Channel group size: {args.channel_group_size}")
    log(f"  Number of chunks: {args.num_chunks if args.num_chunks is not None else 'All (process entire band)'}")
    log(f"  Target frequency: {args.freq_target/1e9:.4f} GHz" if args.freq_target else "  Target frequency: Not specified (start from channel 0)")
    log(f"  Correlation products: {args.correlation_products}")
    log("")
    
    log("OUTLIER DETECTION:")
    log(f"  Sigma threshold: {args.sigma_threshold}σ")
    log(f"  Outlier mode: {args.outlier_mode.upper()}")
    if args.outlier_mode == 'high':
        log(f"    → Flag data points > median + {args.sigma_threshold}σ")
    elif args.outlier_mode == 'low':
        log(f"    → Flag data points < median - {args.sigma_threshold}σ")
    elif args.outlier_mode == 'both':
        log(f"    → Flag data points > median + {args.sigma_threshold}σ OR < median - {args.sigma_threshold}σ")
    elif args.outlier_mode == 'mixed':
        log(f"    → Parallel-hand (XX, YY): Flag high OR low outliers (both directions)")
        log(f"    → Cross-hand (XY, YX): Flag high outliers only")
    log(f"  Baseline flagging threshold: {args.outlier_fraction_threshold*100:.0f}% outliers")
    log("")
    
    log("ANTENNA FLAGGING:")
    log(f"  Method 1 - Absolute threshold: {args.antenna_flag_threshold*100:.0f}%")
    log(f"    → Flag antenna if >{args.antenna_flag_threshold*100:.0f}% of its baselines are bad")
    if args.antenna_relative_sigma <= 0:
        log(f"  Method 2 - Relative outlier detection: DISABLED")
    else:
        log(f"  Method 2 - Relative outlier detection: {args.antenna_relative_sigma:.1f}σ")
        log(f"    → Flag antenna if it's >{args.antenna_relative_sigma:.1f}σ above population median")
    log(f"  Maximum antennas to flag globally: {args.antenna_flag_cap}")
    log("")
    
    log("OUTPUT:")
    flag_filename = f'baseline_flags{args.output_suffix}.txt'
    log(f"  Flags file: {flag_filename}")
    log(f"  Log file: {log_filename}")
    log(f"  Plots directory: {args.save_dir}")
    log(f"  Generate plots: {'Yes' if args.num_chunks is not None and args.num_chunks < 3 else 'No (num_chunks >= 3 or processing all chunks)'}")
    if args.output_suffix:
        log(f"  Output suffix: '{args.output_suffix}'")
    log("")
    
    log("="*70)
    log("")
    
    # =================================================================
    # LOAD MS AND DISPLAY DATA DIMENSIONS
    # =================================================================
    log("Loading measurement set...")
    
    # Load antenna names
    antenna_names = load_antenna_names(args.ms_file)
    
    # Load channel frequencies
    chan_freq = load_chan_freq(args.ms_file)
    
    # Determine start channel based on args.num_chunks
    if args.num_chunks is None:
        # Process all chunks starting from channel 0
        start_search_channel = 0
        log("args.num_chunks=None: starting from channel 0")
    elif args.freq_target is not None:
        # Find closest channel to target frequency
        freq_diffs = np.abs(chan_freq - args.freq_target)
        start_search_channel = np.argmin(freq_diffs)
        log(f"Target frequency: {args.freq_target/1e9:.4f} GHz -> starting search at channel {start_search_channel}")
    else:
        start_search_channel = 0
        log("No args.freq_target specified: starting from channel 0")
    
    # Load MS grouped by baseline
    log("Loading MS data (grouped by baseline)...")
    ds_list = load_grouped_ms(args.ms_file, field_id, args.data_column)
    
    if len(ds_list) == 0:
        raise RuntimeError("No data found for specified field")
    
    # Get data dimensions from first baseline
    n_chan = ds_list[0].dims['chan']
    n_corr = ds_list[0].dims['corr']
    n_time_samples = ds_list[0].dims['row']
    
    # Parse correlation products
    correlations = parse_correlation_products(args.correlation_products)
    
    # Display data dimensions
    log("")
    log("="*70)
    log("DATA DIMENSIONS")
    log("="*70)
    log(f"  Antennas: {len(antenna_names)}")
    log(f"  Baselines: {len(ds_list)}")
    log(f"  Time samples: {n_time_samples:,}")
    log(f"  Frequency channels: {n_chan}")
    log(f"  Frequency range: {chan_freq[0]/1e9:.4f} - {chan_freq[-1]/1e9:.4f} GHz")
    log(f"  Total correlations in MS: {n_corr} (XX, XY, YX, YY)")
    log(f"  Processing correlations: {', '.join([name for name, _ in correlations])}")
    log("")
    
    # Calculate total data size
    total_vis = len(ds_list) * n_time_samples * n_chan * n_corr
    log(f"  Total visibility points: {total_vis:,}")
    log(f"  Estimated data size: {total_vis * 8 / (1024**3):.2f} GB (complex64)")
    log("="*70)
    log("")
    
    # Check system resources and validate chunk size
    log("")
    log("="*60)
    log("SYSTEM RESOURCE CHECK")
    log("="*60)
    resources = check_system_resources()
    log(f"Available RAM: {resources['available_memory_gb']:.2f} GB / {resources['total_memory_gb']:.2f} GB "
        f"({100 - resources['memory_percent_used']:.1f}% free)")
    log(f"CPU cores: {resources['physical_cores']} physical, {resources['logical_cores']} logical")
    log(f"Dask will automatically use all available cores for parallel processing")
    
    # Estimate chunk memory and check if safe
    # NOTE: Now loading all requested correlations together, so multiply by len(correlations)
    n_baselines = len(ds_list)
    n_corr_processing = len(correlations)
    chunk_memory_gb = estimate_chunk_memory(n_baselines, n_time_samples, 
                                            args.channel_group_size, n_corr_processing)
    is_safe, max_safe_gb, _ = check_chunk_memory_safe(chunk_memory_gb, resources['available_memory_gb'])
    
    log("")
    log(f"Estimated memory per chunk (loading {n_corr_processing} correlation(s) together): {chunk_memory_gb:.2f} GB")
    log(f"Memory safety limit (85%): {max_safe_gb:.2f} GB")
    
    if not is_safe:
        reduction_factor = chunk_memory_gb / max_safe_gb
        suggested_chunk_size = max(8, int(args.channel_group_size / reduction_factor))
        log("")
        log(f"WARNING: Chunk memory ({chunk_memory_gb:.2f} GB) exceeds safe limit ({max_safe_gb:.2f} GB)!")
        log(f"  Recommendation: Reduce --channel-group-size to {suggested_chunk_size}")
        log(f"  Current setting: {args.channel_group_size}")
        log(f"  Risk: May cause slowdown or out-of-memory errors")
        log("")
        response = input("Continue anyway? (yes/no): ").strip().lower()
        if response not in ['yes', 'y']:
            log("Exiting as requested.")
            sys.exit(0)
    else:
        log(f"Chunk size is safe for available memory ✓")
    
    log("="*60)
    log("")
    
    # Calculate total number of chunks
    total_chunks = (n_chan - start_search_channel + args.channel_group_size - 1) // args.channel_group_size
    log(f"Total chunks available: {total_chunks}")
    
    # Calculate existing flags and total data size
    log("Calculating existing flag statistics...")
    flag_tasks = [ds.FLAG.sum().data for ds in ds_list]
    total_size_tasks = [ds.FLAG.size for ds in ds_list]
    flag_sums = dask.compute(*flag_tasks)
    old_flags = sum(flag_sums)
    total_data = sum(total_size_tasks)
    old_flag_fraction = old_flags / total_data
    log(f"Existing flags: {old_flags:,} / {total_data:,} ({old_flag_fraction:.6f})")
    
    # Determine which chunks to process based on args.num_chunks
    if args.num_chunks is None:
        log(f"args.num_chunks=None: processing all {total_chunks} chunks...")
        chunks_to_process = list(range(start_search_channel, n_chan, args.channel_group_size))
    else:
        log(f"args.num_chunks={args.num_chunks}: processing {min(args.num_chunks, total_chunks)} chunks...")
        # Find first unflagged chunk
        log("Finding first unflagged channel chunk...")
        first_start = find_first_unflagged_chunk(ds_list, args.channel_group_size, start_search_channel)
        # Generate list of args.num_chunks consecutive chunk starts
        chunks_to_process = []
        for i in range(args.num_chunks):
            chunk_start = first_start + i * args.channel_group_size
            if chunk_start < n_chan:
                chunks_to_process.append(chunk_start)
    
    # Initialize flag file with suffix
    flag_file_path = f'baseline_flags{args.output_suffix}.txt'
    
    # If file exists, rename it with timestamp
    if os.path.exists(flag_file_path):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = f'baseline_flags{args.output_suffix}_{timestamp}.txt'
        os.rename(flag_file_path, backup_path)
        log(f"Renamed existing flag file to: {backup_path}")
    
    # Create new empty file
    with open(flag_file_path, 'w') as f:
        pass  # Create empty file
    log(f"Initialized flag file: {flag_file_path}")
    
    # Track cumulative new flags
    cumulative_new_flags = 0
    
    # Track antenna statistics across all chunks
    global_antenna_stats = {}  # {ant_name: {'total': N, 'flagged': M, 'chunks': [...]}}
    
    # Display processing plan
    log("")
    log("="*70)
    log("PROCESSING PLAN")
    log("="*70)
    log(f"  Total chunks to process: {len(chunks_to_process)}")
    if len(chunks_to_process) > 0:
        log(f"  Starting channel: {chunks_to_process[0]}")
        log(f"  Ending channel: {min(chunks_to_process[-1] + args.channel_group_size, n_chan)}")
        log(f"  Channel range: {chunks_to_process[0]} - {min(chunks_to_process[-1] + args.channel_group_size - 1, n_chan - 1)}")
    log(f"  Correlations per chunk: {len(correlations)} ({', '.join([name for name, _ in correlations])})")
    log(f"  Baselines per correlation: {len(ds_list)}")
    log(f"  Total operations: {len(chunks_to_process)} chunks × {len(correlations)} correlations = {len(chunks_to_process) * len(correlations)} correlation-chunks")
    log("")
    log(f"  Output files:")
    log(f"    - Flags: {flag_file_path}")
    log(f"    - Log: {log_filename}")
    if args.num_chunks is not None and args.num_chunks < 3:
        log(f"    - Plots: {args.save_dir}amp_vs_baseline_chunk{{N}}_{{corr}}.png")
        log(f"    - Plots: {args.save_dir}amp_histogram_chunk{{N}}_{{corr}}.png")
    log("="*70)
    log("")
    log("Starting processing...")
    log("")
    
    # Build baseline info once (doesn't change across chunks)
    log("Building baseline information...")
    baseline_info = []  # Store (a1, a2, bl_name, ant1_name, ds)
    for ds in ds_list:
        # Get baseline identity from attrs
        a1 = int(ds.attrs.get("ANTENNA1", -1))
        a2 = int(ds.attrs.get("ANTENNA2", -1))
        
        if a1 < 0 or a2 < 0 or a1 >= len(antenna_names) or a2 >= len(antenna_names):
            continue
        
        bl_name = f"{antenna_names[a1]}&{antenna_names[a2]}"
        ant1_name = antenna_names[a1]
        
        baseline_info.append((a1, a2, bl_name, ant1_name, ds))
    
    # Sort baselines by ANTENNA1, then ANTENNA2 for clear color blocks
    baseline_info.sort(key=lambda x: (x[0], x[1]))
    log(f"Processed {len(baseline_info)} baselines")
    
    # Process each chunk
    for chunk_idx, start_chan in enumerate(chunks_to_process):
        end_chan = min(start_chan + args.channel_group_size, n_chan)
        chunk_num = chunk_idx + 1
        total_to_process = len(chunks_to_process)
        
        log("")
        log(f"{'='*60}")
        log(f"Processing chunk {chunk_num}/{total_to_process}: channels {start_chan} to {end_chan}")
        log(f"{'='*60}")
        
        # Track if we hit memory error (initialize before any checks)
        memory_error_stop = False
        
        # Check memory before processing chunk (raises error if >95%)
        try:
            check_memory_critical()
        except MemoryError as e:
            log(str(e))
            log("Stopping processing to prevent system crash")
            memory_error_stop = True
            break
        
        # PRELIMINARY CHECK: Is this chunk fully flagged?
        # Check first correlation across all baselines to save time
        log("Checking if chunk is fully flagged...")
        chunk_flag_check_tasks = []
        for _, _, _, _, ds in baseline_info:
            # Just check first correlation (0) as representative
            flag_slice = ds.FLAG[:, start_chan:end_chan, 0]
            chunk_flag_check_tasks.append(flag_slice.data)
        
        # Compute flag checks
        chunk_flags_computed = dask.compute(*chunk_flag_check_tasks)
        
        # Check if all data in chunk is flagged
        total_unflagged = sum((~flag_arr).sum() for flag_arr in chunk_flags_computed)
        
        if total_unflagged == 0:
            log("")
            log(f"WARNING: Chunk {chunk_num} is fully flagged (checked correlation 0) - skipping all correlations")
            continue
        
        log(f"Chunk has {total_unflagged:,} unflagged data points - proceeding")
        
        # Track flagged baselines across all correlations for this chunk
        chunk_flagged_baselines = set()  # Union across all correlations
        
        # =================================================================
        # STEP 1: Build tasks for ALL correlations (lazy, no computation yet)
        # =================================================================
        log("")
        log(f"Preparing data extraction for {len(correlations)} correlation(s)...")
        
        all_correlation_tasks = []  # List of (corr_name, corr_idx, amp_tasks)
        baseline_names = [x[2] for x in baseline_info]  # Same for all correlations
        ant1_names = [x[3] for x in baseline_info]
        
        for corr_name, corr_idx in correlations:
            amp_tasks = []
            
            # Extract amplitude and weight data in sorted order (LAZY)
            for i, (a1, a2, bl_name, ant1_name, ds) in enumerate(baseline_info):
                # Slice and compute amplitude (still lazy until compute)
                vis_slice = ds.DATA[:, start_chan:end_chan, corr_idx]
                flag_slice = ds.FLAG[:, start_chan:end_chan, corr_idx]
                weight_slice = ds.WEIGHT[:, corr_idx]  # Weight for this correlation, shape: (row,)
            
                # Compute amplitude and apply flags
                amp = abs(vis_slice)
                amp_flagged = amp.where(~flag_slice, other=np.nan)
            
                # Broadcast weight to match amplitude shape
                weight_broadcasted, _ = dask.array.broadcast_arrays(weight_slice.data[:, None], amp.data)
                weight_flagged = amp.copy(data=weight_broadcasted)
                weight_flagged = weight_flagged.where(~flag_slice, other=np.nan)
            
                amp_tasks.append(amp_flagged.data)
                amp_tasks.append(weight_flagged.data)  # Append weight right after corresponding amplitude
            
            all_correlation_tasks.append((corr_name, corr_idx, amp_tasks))
        
        # =================================================================
        # STEP 2: Single compute for ALL correlations
        # =================================================================
        log(f"Computing amplitudes and weights for {len(baseline_info)} baselines × {len(correlations)} correlation(s)...")
        
        # Flatten all tasks into one list
        flat_tasks = []
        for _, _, tasks in all_correlation_tasks:
            flat_tasks.extend(tasks)
        
        # Single dask.compute() for all correlations
        all_computed_data = dask.compute(*flat_tasks)
        
        # Check memory after compute (large memory spike)
        try:
            mem_percent = check_memory_critical()
        except MemoryError as e:
            log(str(e))
            log("Stopping processing to prevent system crash")
            memory_error_stop = True
            break
        
        # =================================================================
        # STEP 3: Parse results by correlation and analyze sequentially
        # =================================================================
        result_idx = 0
        n_tasks_per_corr = len(baseline_info) * 2  # amp + weight per baseline
        
        for corr_name, corr_idx, _ in all_correlation_tasks:
            log("")
            log(f"Analyzing correlation: {corr_name} (index {corr_idx})")
            
            # Extract this correlation's data from computed results
            computed_data = all_computed_data[result_idx:result_idx + n_tasks_per_corr]
            result_idx += n_tasks_per_corr
    
            # Separate amplitudes and weights from interleaved results
            amp_arrays = [computed_data[i*2] for i in range(len(baseline_info))]
            weight_arrays = [computed_data[i*2 + 1] for i in range(len(baseline_info))]
            
            # SECONDARY CHECK: if no valid amplitude data for this correlation, skip
            total_valid = sum(np.isfinite(amp).sum() for amp in amp_arrays)
            if total_valid == 0:
                log("")
                log(f"WARNING: Chunk {chunk_num}, {corr_name} has no valid amplitude data - skipping correlation")
                continue
    
            # Flatten data for analysis (used for outlier detection and optional plotting)
            log("Preparing data for analysis...")
            
            # Collect data in lists (faster than building many DataFrames)
            baseline_idx_list = []
            amplitude_list = []
            weight_list = []
            ant1_list = []
    
            for i, (amp_arr, weight_arr, bl_name, ant1_name) in enumerate(zip(amp_arrays, weight_arrays, baseline_names, ant1_names)):
                # Flatten and filter valid data
                amp_flat = amp_arr.flatten()
                weight_flat = weight_arr.flatten()
                valid_mask = np.isfinite(amp_flat) & np.isfinite(weight_flat)
            
                if not np.any(valid_mask):
                    continue
            
                n_points = valid_mask.sum()
                
                # Append to lists (fast)
                baseline_idx_list.append(np.full(n_points, i, dtype=np.int32))
                amplitude_list.append(amp_flat[valid_mask])
                weight_list.append(weight_flat[valid_mask])
                ant1_list.append(np.full(n_points, ant1_name, dtype=object))
    
            # Concatenate numpy arrays (fast) then build single DataFrame (fast)
            df = pd.DataFrame({
                'baseline_idx': np.concatenate(baseline_idx_list),
                'amplitude': np.concatenate(amplitude_list),
                'weight': np.concatenate(weight_list),
                'ant1': np.concatenate(ant1_list),
            })
            log(f"Total data points: {len(df):,}")
    
            # Convert ant1 to categorical for datashader compatibility
            df['ant1'] = pd.Categorical(df['ant1'])
    
            # Calculate statistics using weighted median for better center estimation
            amp_values = df['amplitude'].values
            weight_values = df['weight'].values
            median_amp = weighted_median(amp_values, weight_values)
            # MAD calculation using weighted median as center
            mad_amp = bn.nanmedian(np.abs(amp_values - median_amp))
            log("")
            log("Amplitude statistics:")
            log(f"  Weighted Median: {median_amp:.4f}")
            log(f"  MAD: {mad_amp:.4f}")
    
            # Outlier detection per baseline
            log("")
            log(f"Outlier detection (threshold: {args.sigma_threshold} sigma):")
            threshold_mad = args.sigma_threshold * 1.4826 * mad_amp
            log(f"  MAD threshold: {threshold_mad:.4f}")
    
            # Determine effective outlier mode based on correlation type
            if args.outlier_mode == 'mixed':
                # Parallel-hand (XX, YY) -> both, Cross-hand (XY, YX) -> high
                if corr_name in ['XX', 'YY']:
                    effective_mode = 'both'
                    log(f"  Mode: MIXED → Using 'both' for parallel-hand correlation {corr_name}")
                else:  # XY, YX
                    effective_mode = 'high'
                    log(f"  Mode: MIXED → Using 'high' for cross-hand correlation {corr_name}")
            else:
                effective_mode = args.outlier_mode
                log(f"  Mode: {effective_mode.upper()}")
    
            # Mark outliers based on effective_mode
            if effective_mode == 'high':
                df['is_outlier'] = (df['amplitude'] - median_amp) > threshold_mad
            elif effective_mode == 'low':
                df['is_outlier'] = (median_amp - df['amplitude']) > threshold_mad
            else:  # 'both'
                df['is_outlier'] = np.abs(df['amplitude'] - median_amp) > threshold_mad
    
            # Per-baseline outlier statistics using groupby
            baseline_stats = df.groupby('baseline_idx').agg({
                'is_outlier': ['sum', 'count']
            })
            baseline_stats.columns = ['outlier_count', 'total_count']
            baseline_stats['outlier_fraction'] = baseline_stats['outlier_count'] / baseline_stats['total_count']
    
            # Filter to baselines with outliers and sort by fraction
            outlier_baselines = baseline_stats[baseline_stats['outlier_count'] > 0].sort_values(
                'outlier_fraction', ascending=False
            )
    
            if len(outlier_baselines) > 0:
                log("")
                log(f"Baselines with outliers ({len(outlier_baselines)} total):")
        
                # Truncate output if not verbose and more than 10 outlier baselines
                if not args.verbose and len(outlier_baselines) > 10:
                    # Print first 5
                    for i, (bl_idx, row) in enumerate(outlier_baselines.head(5).iterrows()):
                        bl_name = baseline_names[bl_idx]
                        log(f"  {bl_name}: {row['outlier_fraction']:.4f} "
                              f"({int(row['outlier_count'])}/{int(row['total_count'])} outliers)")
            
                    log("  ...")
            
                    # Print last 5
                    for i, (bl_idx, row) in enumerate(outlier_baselines.tail(5).iterrows()):
                        bl_name = baseline_names[bl_idx]
                        log(f"  {bl_name}: {row['outlier_fraction']:.4f} "
                              f"({int(row['outlier_count'])}/{int(row['total_count'])} outliers)")
                else:
                    # Print all
                    for bl_idx, row in outlier_baselines.iterrows():
                        bl_name = baseline_names[bl_idx]
                        log(f"  {bl_name}: {row['outlier_fraction']:.4f} "
                              f"({int(row['outlier_count'])}/{int(row['total_count'])} outliers)")
            else:
                log("")
                log("No baselines with outliers detected.")
    
            # Global outlier fraction
            global_outlier_frac = df['is_outlier'].sum() / len(df)
            log("")
            log(f"Global outlier fraction: {global_outlier_frac:.4f}")
        
            # Write bad baselines to flag file and calculate new flags
            bad_baselines_chunk = outlier_baselines[
                outlier_baselines['outlier_fraction'] > args.outlier_fraction_threshold
            ].copy()
            
            # Track flagged baselines for this correlation
            for bl_idx in bad_baselines_chunk.index:
                chunk_flagged_baselines.add(bl_idx)
        
            if len(bad_baselines_chunk) > 0:
                log("")
                log(f"Flagging {len(bad_baselines_chunk)} baselines for {corr_name}")
                
                # Just track for now, will write flags after all correlations processed
            else:
                log("")
                log(f"No baselines flagged for {corr_name}")
            
            # Create plot for this chunk/correlation if args.num_chunks < 3
            if args.num_chunks is not None and args.num_chunks < 3:
                log("")
                log(f"Generating plots for chunk {chunk_num}, {corr_name}...")
                create_plot(df, baseline_names, median_amp, mad_amp, threshold_mad, 
                            outlier_baselines, start_chan, end_chan, chunk_num, corr_name,
                            args.outlier_fraction_threshold, effective_mode, args.save_dir, field_name)
                create_histogram_plot(df, median_amp, mad_amp, threshold_mad, chunk_num, corr_name,
                                     effective_mode, args.save_dir)
        
        # End of correlation loop - now process antenna-level flagging
        
        # Check if we stopped due to memory error
        if memory_error_stop:
            log("Exiting chunk processing loop due to memory constraints")
            break
        
        log("")
        log(f"Analyzing antenna-level flagging for chunk {chunk_num}...")
        
        # Count how many baselines each antenna participates in (total and flagged)
        antenna_stats = {}  # {ant_name: {'total': N, 'flagged': M}}
        
        for bl_idx, (a1, a2, bl_name, ant1_name, ds) in enumerate(baseline_info):
            ant1 = antenna_names[a1]
            ant2 = antenna_names[a2]
            
            # Initialize if not seen
            if ant1 not in antenna_stats:
                antenna_stats[ant1] = {'total': 0, 'flagged': 0}
            if ant2 not in antenna_stats:
                antenna_stats[ant2] = {'total': 0, 'flagged': 0}
            
            # Count baseline for both antennas
            antenna_stats[ant1]['total'] += 1
            antenna_stats[ant2]['total'] += 1
            
            # If this baseline is flagged, count for both antennas
            if bl_idx in chunk_flagged_baselines:
                antenna_stats[ant1]['flagged'] += 1
                antenna_stats[ant2]['flagged'] += 1
        
        # Identify antennas exceeding threshold
        flagged_antennas = []
        for ant_name, stats in antenna_stats.items():
            if stats['total'] > 0:
                fraction = stats['flagged'] / stats['total']
                if fraction > args.antenna_flag_threshold:
                    flagged_antennas.append((ant_name, fraction, stats['flagged'], stats['total']))
        
        # Update global antenna statistics
        for ant_name, stats in antenna_stats.items():
            if ant_name not in global_antenna_stats:
                global_antenna_stats[ant_name] = {'total': 0, 'flagged': 0, 'chunks': []}
            global_antenna_stats[ant_name]['total'] += stats['total']
            global_antenna_stats[ant_name]['flagged'] += stats['flagged']
            if stats['flagged'] > 0:
                global_antenna_stats[ant_name]['chunks'].append(chunk_num)
        
        if len(flagged_antennas) > 0:
            log("")
            log(f"Antennas exceeding {args.antenna_flag_threshold*100:.0f}% threshold:")
            for ant_name, frac, flagged, total in sorted(flagged_antennas, key=lambda x: x[1], reverse=True):
                log(f"  {ant_name}: {frac:.2%} ({flagged}/{total} baselines flagged) - FLAGGING ENTIRE ANTENNA")
        else:
            log("No antennas exceed the antenna flagging threshold")
        
        # Write flags to file (baselines + antennas)
        if len(chunk_flagged_baselines) > 0 or len(flagged_antennas) > 0:
            log("")
            log("Writing flags to baseline_flags.txt...")
            
            # Prepare baseline flags (sorted by ANTENNA1, ANTENNA2)
            bad_baseline_data = []
            for bl_idx in chunk_flagged_baselines:
                a1, a2, bl_name, _, _ = baseline_info[bl_idx]
                bad_baseline_data.append((a1, a2, bl_idx, bl_name))
            bad_baseline_data.sort(key=lambda x: (x[0], x[1]))
            
            # Write baseline flags (all correlations since any bad = all bad)
            with open(flag_file_path, 'a') as f:
                for base_num, (a1, a2, bl_idx, bl_name) in enumerate(bad_baseline_data):
                    reason = f"BASED_FLAGGER_CHUNK_{chunk_num}_BASE{base_num}"
                    flag_line = f"mode='manual' antenna='{bl_name}' spw='0:{start_chan}~{end_chan-1}' reason='{reason}'\n"
                    f.write(flag_line)
            
            # Write antenna flags
            with open(flag_file_path, 'a') as f:
                for ant_idx, (ant_name, frac, flagged, total) in enumerate(flagged_antennas):
                    reason = f"BASED_FLAGGER_CHUNK_{chunk_num}_ANT{ant_idx}"
                    flag_line = f"mode='manual' antenna='{ant_name}' spw='0:{start_chan}~{end_chan-1}' reason='{reason}'\n"
                    f.write(flag_line)
            
            # Calculate new flags for this chunk
            n_freq_in_chunk = end_chan - start_chan
            n_corr_flagged = len(correlations)
            n_bad_baselines = len(bad_baseline_data)
            
            # Baseline flags
            chunk_new_flags = n_bad_baselines * n_time_samples * n_freq_in_chunk * n_corr_flagged
            # Antenna flags (each antenna participates in many baselines)
            for ant_name, _, _, _ in flagged_antennas:
                # Count baselines involving this antenna
                ant_baseline_count = sum(1 for a1, a2, _, _, _ in baseline_info 
                                        if antenna_names[a1] == ant_name or antenna_names[a2] == ant_name)
                chunk_new_flags += ant_baseline_count * n_time_samples * n_freq_in_chunk * n_corr_flagged
            
            cumulative_new_flags += chunk_new_flags
            log(f"New flags this chunk: {chunk_new_flags:,}")
            log(f"Cumulative new flags: {cumulative_new_flags:,}")
        
        # Report running flag tally
        total_flags_with_new = old_flags + cumulative_new_flags
        total_flag_fraction = total_flags_with_new / total_data
        log(f"Total flags (old + new): {total_flags_with_new:,} / {total_data:,} ({total_flag_fraction:.6f})")
    
    # Global antenna flagging across all chunks
    log("")
    log("="*60)
    log("GLOBAL ANTENNA FLAGGING (ALL CHUNKS)")
    log("="*60)
    
    # Build list of all antennas with their flagging fractions
    antenna_fractions = []
    for ant_name, stats in global_antenna_stats.items():
        if stats['total'] > 0:
            fraction = stats['flagged'] / stats['total']
            antenna_fractions.append((ant_name, fraction, stats['flagged'], stats['total']))
    
    if len(antenna_fractions) == 0:
        log("No antenna statistics available")
    else:
        log(f"Analyzing {len(antenna_fractions)} antennas...")
        log("")
        
        # METHOD 1: Absolute threshold
        absolute_flagged = []
        for ant_name, frac, flagged, total in antenna_fractions:
            if frac > args.antenna_flag_threshold:
                absolute_flagged.append((ant_name, frac, flagged, total, 'absolute'))
        
        if len(absolute_flagged) > 0:
            log(f"Method 1 - Absolute threshold (>{args.antenna_flag_threshold*100:.0f}%):")
            for ant_name, frac, flagged, total, _ in sorted(absolute_flagged, key=lambda x: x[1], reverse=True):
                log(f"  {ant_name}: {frac:.2%} ({flagged}/{total} baselines)")
        else:
            log(f"Method 1 - Absolute threshold: No antennas exceed {args.antenna_flag_threshold*100:.0f}%")
        
        log("")
        
        # METHOD 2: Relative outlier detection
        relative_flagged = []
        
        if args.antenna_relative_sigma <= 0:
            log(f"Method 2 - Relative outlier detection: DISABLED (sigma={args.antenna_relative_sigma})")
        else:
            fracs = np.array([frac for _, frac, _, _ in antenna_fractions])
            median_frac = np.median(fracs)
            mad_frac = np.median(np.abs(fracs - median_frac))
            threshold_frac = median_frac + args.antenna_relative_sigma * 1.4826 * mad_frac
            
            log(f"Method 2 - Relative outlier detection ({args.antenna_relative_sigma:.1f}σ):")
            log(f"  Population median: {median_frac:.2%}")
            log(f"  Population MAD: {mad_frac:.2%}")
            log(f"  Outlier threshold: {threshold_frac:.2%}")
            
            for ant_name, frac, flagged, total in antenna_fractions:
                if frac > threshold_frac:
                    relative_flagged.append((ant_name, frac, flagged, total, 'relative'))
            
            if len(relative_flagged) > 0:
                log(f"  Found {len(relative_flagged)} outlier(s):")
                for ant_name, frac, flagged, total, _ in sorted(relative_flagged, key=lambda x: x[1], reverse=True):
                    sigma_above = (frac - median_frac) / (mad_frac * 1.4826) if mad_frac > 0 else 0
                    log(f"    {ant_name}: {frac:.2%} ({sigma_above:.1f}σ above median)")
            else:
                log(f"  No antennas are {args.antenna_relative_sigma:.1f}σ outliers")
        
        log("")
        
        # COMBINE: Union of both methods
        # Track which method flagged each antenna
        flagged_dict = {}
        for ant_name, frac, flagged, total, method in absolute_flagged:
            flagged_dict[ant_name] = {'frac': frac, 'flagged': flagged, 'total': total, 'methods': [method]}
        for ant_name, frac, flagged, total, method in relative_flagged:
            if ant_name in flagged_dict:
                flagged_dict[ant_name]['methods'].append(method)
            else:
                flagged_dict[ant_name] = {'frac': frac, 'flagged': flagged, 'total': total, 'methods': [method]}
        
        # Sort by flagging fraction
        combined_flagged = sorted(
            [(ant, info['frac'], info['flagged'], info['total'], info['methods']) 
             for ant, info in flagged_dict.items()],
            key=lambda x: x[1],
            reverse=True
        )
        
        # Apply cap
        if len(combined_flagged) > args.antenna_flag_cap:
            log(f"Combined: {len(combined_flagged)} antennas flagged by at least one method")
            log(f"Applying cap: Flagging worst {args.antenna_flag_cap} antennas only")
            capped_flagged = combined_flagged[:args.antenna_flag_cap]
            excluded = combined_flagged[args.antenna_flag_cap:]
            log("")
            log(f"Antennas being flagged:")
            for ant_name, frac, flagged, total, methods in capped_flagged:
                methods_str = '+'.join(methods)
                log(f"  {ant_name}: {frac:.2%} ({flagged}/{total} baselines) [{methods_str}]")
            log("")
            log(f"Antennas excluded by cap (not flagged):")
            for ant_name, frac, flagged, total, methods in excluded[:5]:  # Show first 5
                methods_str = '+'.join(methods)
                log(f"  {ant_name}: {frac:.2%} ({flagged}/{total} baselines) [{methods_str}]")
            if len(excluded) > 5:
                log(f"  ... and {len(excluded)-5} more")
        else:
            capped_flagged = combined_flagged
            if len(capped_flagged) > 0:
                log(f"Flagging {len(capped_flagged)} antenna(s):")
                for ant_name, frac, flagged, total, methods in capped_flagged:
                    methods_str = '+'.join(methods)
                    log(f"  {ant_name}: {frac:.2%} ({flagged}/{total} baselines) [{methods_str}]")
            else:
                log("No antennas meet flagging criteria")
        
        # Write global antenna flags to file
        if len(capped_flagged) > 0 and len(chunks_to_process) > 0:
            min_chan = min(chunks_to_process)
            max_chan = max(chunks_to_process) + args.channel_group_size
            max_chan = min(max_chan, n_chan)  # Don't exceed total channels
            
            log("")
            log(f"Writing {len(capped_flagged)} global antenna flag(s) for channels {min_chan}~{max_chan-1}...")
            
            with open(flag_file_path, 'a') as f:
                for ant_idx, (ant_name, frac, flagged, total, methods) in enumerate(capped_flagged):
                    reason = f"BASED_FLAGGER_GLOBAL_ANT{ant_idx}"
                    flag_line = f"mode='manual' antenna='{ant_name}' spw='0:{min_chan}~{max_chan-1}' reason='{reason}'\n"
                    f.write(flag_line)
            
            log(f"✓ Wrote {len(capped_flagged)} global antenna flags to {flag_file_path}")
    
    # Print top 5 antennas with highest flagging fraction across all chunks
    log("")
    log("="*60)
    log("TOP 5 ANTENNAS WITH HIGHEST FLAGGING FRACTION (ALL CHUNKS)")
    log("="*60)
    
    # Calculate fractions and sort
    antenna_fractions = []
    for ant_name, stats in global_antenna_stats.items():
        if stats['total'] > 0:
            fraction = stats['flagged'] / stats['total']
            antenna_fractions.append((ant_name, fraction, stats['flagged'], stats['total'], stats['chunks']))
    
    antenna_fractions.sort(key=lambda x: x[1], reverse=True)
    
    # Print top 5
    top_5 = antenna_fractions[:5]
    if len(top_5) > 0:
        for rank, (ant_name, frac, flagged, total, chunks) in enumerate(top_5, start=1):
            chunks_str = ','.join(map(str, sorted(set(chunks))))
            log(f"{rank}. {ant_name}: {frac:.2%} ({flagged}/{total} baselines) - chunks: {chunks_str}")
    else:
        log("No antenna flagging data available")
    log("")


if __name__ == "__main__":
    main()
