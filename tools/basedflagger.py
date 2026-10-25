#!/usr/bin/env python3
"""
Visibility data flagger for radio interferometry observations.
Performs outlier detection on amplitude vs baseline distributions, and with
--phase on phase vs baseline as well (parallel hands only).
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
from datetime import datetime, timedelta
import json
import logging
import argparse
import psutil
from statistics import NormalDist
import re
import sys
import shutil
import atexit

# Correlations the phase pass is allowed to run on. Cross-hand phase carries
# the polarisation position angle, which sweeps with parallactic angle over a
# track, so testing it for outliers would flag real signal.
PHASE_CORRELATIONS = ('XX', 'YY')


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
    parser.add_argument('--flag-fields', type=str, default='',
                        help='Comma-separated field IDs/names that flags found on '
                             '`field` are passed on to: a flag on a scan of `field` also '
                             'goes to every scan of these fields between the previous and '
                             'the next scan of `field` (or the start/end of the '
                             'observation where there is none)')

    # Optional arguments with defaults
    parser.add_argument('--data-column', type=str, default='CORRECTED_DATA',
                        help='Data column to read')
    parser.add_argument('--channel-group-size', type=int, default=16,
                        help='Number of channels per chunk')
    parser.add_argument('--freq-target', type=float, default=1.0e9,
                        help='Target frequency in Hz (starting point for chunk search)')
    parser.add_argument('--save-dir', type=str, default='VISPLOTS/',
                        help='Directory to save plots')
    parser.add_argument('--output-suffix', type=str, default='',
                        help='Suffix for output files (e.g., "_test" -> baseline_flags_test.txt)')
    parser.add_argument('--sigma-threshold', type=float, default=5.0,
                        help='Sigma threshold for outlier detection')
    parser.add_argument('--outlier-mode', type=str, default='mixed',
                        choices=['high', 'low', 'both', 'mixed'],
                        help='Outlier detection mode: high, low, both, or mixed (both for XX/YY, high for XY/YX)')
    parser.add_argument('--outlier-fraction-threshold', type=float, default=0.05,
                        help='Flag baselines with >X fraction outliers. Together with '
                             '--sigma-threshold this sets how much broadening trips a '
                             'flag: a baseline broadened by factor b is caught when '
                             'b > sigma_threshold/z(X), for z the normal quantile. '
                             'The resulting sensitivity is reported in the run log')
    parser.add_argument('--point-sigma-threshold', type=float, default=7.0,
                        help='Stricter sigma threshold for individual points. A baseline '
                             'is flagged outright the moment >= --point-outlier-min-count '
                             'points cross this, independent of --outlier-fraction-threshold '
                             '-- catches a short, severe event that a fraction test would '
                             'dilute away in an otherwise-good baseline. <=0 disables it')
    parser.add_argument('--point-outlier-min-count', type=int, default=1,
                        help='Number of --point-sigma-threshold outliers required to flag '
                             'a baseline outright; <=0 disables this point-level override')
    parser.add_argument('--scan-window-count', type=int, default=4,
                        help="Number of nearest scans in time used to build a scan's local "
                             'reference amplitude (each neighbour weighted by 1/MAD^2, so a '
                             "noisier neighbour counts for less). A scan whose own weighted "
                             'median deviates from that reference by more than '
                             '--scan-window-sigma has its centre REPLACED by the reference '
                             'before baselines in it are tested -- otherwise a scan-wide '
                             "shift (bad solution, gain drop) moves every baseline's residual "
                             "together, so nothing in that scan looks anomalous against its "
                             'own (equally shifted) median. 0 disables this entirely')
    parser.add_argument('--scan-window-sigma', type=float, default=5.0,
                        help='Sigma threshold for --scan-window-count')
    parser.add_argument('--band-reference', type=str, default='poly',
                        choices=['none', 'poly'],
                        help="'poly' runs a statistics-only pass over every chunk first and "
                             'fits a smooth polynomial across the band to each '
                             "(correlation, scan)'s amplitude centre and MAD (in log space, "
                             'iteratively dropping chunks that sit above the fit). In the '
                             "flagging pass, a scan whose own MAD or centre in a chunk is "
                             'inflated relative to that fit -- broadband RFI filling most of '
                             'a chunk raises both, which otherwise lifts the threshold with '
                             'it -- is tested against the fitted values instead. Amplitude '
                             'passes only; the phase pass is unaffected')
    parser.add_argument('--band-poly-degree', type=int, default=3,
                        help='Polynomial degree for --band-reference poly')
    parser.add_argument('--band-clip-sigma', type=float, default=3.0,
                        help='Chunks lying more than this (MAD-scaled) sigma above the band '
                             'fit are dropped from it, iteratively, for --band-reference poly')
    parser.add_argument('--band-mad-cap', type=float, default=1.5,
                        help="A scan's own MAD in a chunk above this factor times the band "
                             'fit marks it contaminated (--band-reference poly)')
    parser.add_argument('--band-centre-sigma', type=float, default=5.0,
                        help="A scan's own centre in a chunk more than this many fitted-MAD "
                             'sigma above the band fit marks it contaminated '
                             '(--band-reference poly)')
    parser.add_argument('--log-dir', type=str, default='LOGS/',
                        help='Directory for the run log')
    parser.add_argument('--cache-dir', type=str, default='.',
                        help='Directory under which the temporary per-chunk data cache is '
                             'written (float32 amplitudes, plus phases for --phase; '
                             'removed when the run ends). Needs roughly '
                             'n_baselines x n_rows x n_channels x 4 bytes per pass')
    parser.add_argument('--correlation-products', type=str, default='XX,YY,XY,YX',
                        help='Correlation products to process (comma-separated or list)')
    parser.add_argument('--phase', action='store_true',
                        help='Additionally run the outlier tests on visibility phase. '
                             'The amplitude pass is unaffected and still runs on every '
                             f'selected correlation; the phase pass runs on '
                             f'{"/".join(PHASE_CORRELATIONS)} only and is never applied to the '
                             'cross hands, whose phase tracks the polarisation position '
                             'angle and sweeps with parallactic angle. Requires calibrated '
                             'data (a --data-column with gains applied), since on '
                             'uncalibrated data every baseline carries its own delay slope')
    parser.add_argument('--antenna-flag-threshold', type=float, default=0.90,
                        help='Flag an antenna in a scan if >X fraction of its baselines are bad '
                             'in that scan (per chunk, and over all chunks for the global pass)')
    parser.add_argument('--antenna-relative-sigma', type=float, default=5.0,
                        help='Sigma threshold for relative antenna outlier detection '
                             '(<=0 disables it; off by default because the threshold '
                             'degenerates to the population median when the MAD is zero)')
    parser.add_argument('--antenna-flag-cap', type=int, default=3,
                        help='Maximum number of antennas to flag globally in each scan')
    parser.add_argument('--time-chunk', type=float, default=300.0,
                        help='Targets only: flag in time windows of this many seconds '
                             'instead of whole scans (0 = whole scans). '
                             'Not used on calibrators (runs with --flag-fields).')
    parser.add_argument('--protect-antenna', type=str, default='auto',
                        help="Antenna that is never flagged as a whole, normally the reference "
                             "antenna. 'auto' = first ref_ant in project_info.json, 'none' = off.")
    parser.add_argument('--num-chunks', type=int, default=None,
                        help='Number of chunks to process (None = all chunks)')
    parser.add_argument('--no-printing', action='store_true',
                        help='Disable console printing (logging still active)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    
    return parser.parse_args()


def _parse_slurm_mem(value):
    """SLURM memory string -> bytes. Bare numbers are MB; M/G/T suffixes honoured."""
    suffix = value[-1].lower()
    if suffix == 'g':
        multiplier = 1024 ** 3
    elif suffix == 't':
        multiplier = 1024 ** 4
    else:
        multiplier = 1024 ** 2
    return int(value.rstrip('MmGgTt')) * multiplier


def get_slurm_memory_bytes():
    """
    Memory allocated to this job by SLURM, in bytes, or (None, None).

    Returns (bytes, source_string) so the caller can report which variable
    the number came from.
    """
    mem_per_node = os.environ.get('SLURM_MEM_PER_NODE')
    mem_per_cpu = os.environ.get('SLURM_MEM_PER_CPU')

    try:
        if mem_per_node:
            return _parse_slurm_mem(mem_per_node), f'SLURM_MEM_PER_NODE={mem_per_node}'

        if mem_per_cpu:
            per_cpu_bytes = _parse_slurm_mem(mem_per_cpu)
            cpus = os.environ.get('SLURM_CPUS_PER_TASK')
            ntasks = os.environ.get('SLURM_NTASKS')
            if cpus and ntasks:
                n_cpus = int(cpus) * int(ntasks)
            elif cpus:
                n_cpus = int(cpus)
            elif ntasks:
                n_cpus = int(ntasks)
            else:
                n_cpus = 1
            return (per_cpu_bytes * n_cpus,
                    f'SLURM_MEM_PER_CPU={mem_per_cpu} x {n_cpus} CPUs')
    except ValueError:
        pass

    return None, None


def get_allocated_cpus():
    """CPUs this job may actually use, honouring SLURM. Returns (count, source)."""
    cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    ntasks = os.environ.get('SLURM_NTASKS')
    job_cpus = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    job_cpus_match = re.match(r'(\d+)', job_cpus)

    try:
        if cpus_per_task and ntasks:
            return (int(cpus_per_task) * int(ntasks),
                    f'SLURM_CPUS_PER_TASK={cpus_per_task} x SLURM_NTASKS={ntasks}')
        if cpus_per_task:
            return int(cpus_per_task), f'SLURM_CPUS_PER_TASK={cpus_per_task}'
        if job_cpus_match:
            n = int(job_cpus_match.group(1))
            return n, f'SLURM_JOB_CPUS_PER_NODE={job_cpus} -> {n}'
        if ntasks:
            return int(ntasks), f'SLURM_NTASKS={ntasks}'
    except ValueError:
        pass

    return os.cpu_count() or 1, 'os.cpu_count() (no SLURM allocation detected)'


def check_system_resources():
    """
    Check available memory and CPU cores, honouring any SLURM allocation.

    Inside a batch job psutil describes the node, not the job: on a shared
    node it reports memory and cores this job is not entitled to, so sizing
    chunks from it invites an OOM kill. The reported figures are therefore the
    more restrictive of the SLURM allocation and what the node actually has
    free. Matches the resource handling in tools/pbcor_katbeam_loop.py and
    oxkat/3GC_split_model_images.py.

    Returns
    -------
    dict
        Memory in GB, node core counts, CPUs allocated to this job, and a
        source string for each limit so the log records where it came from.
    """
    mem = psutil.virtual_memory()
    slurm_mem, slurm_src = get_slurm_memory_bytes()
    n_cpus, cpu_src = get_allocated_cpus()

    if slurm_mem is not None:
        total = min(slurm_mem, mem.total)
        available = min(slurm_mem, mem.available)
        mem_src = f'min(SLURM, available)  [{slurm_src}]'
    else:
        total = mem.total
        available = mem.available
        mem_src = 'node available (no SLURM allocation detected)'

    percent_used = 100.0 * (1.0 - available / total) if total > 0 else 100.0

    return {
        'available_memory_gb': available / (1024**3),
        'total_memory_gb': total / (1024**3),
        'memory_percent_used': percent_used,
        'physical_cores': psutil.cpu_count(logical=False),
        'logical_cores': psutil.cpu_count(logical=True),
        'allocated_cpus': n_cpus,
        'memory_source': mem_src,
        'cpu_source': cpu_src,
    }


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
        (is_safe, max_safe_gb)
    """
    max_safe_gb = available_memory_gb * 0.85
    return chunk_memory_gb <= max_safe_gb, max_safe_gb


def check_memory_critical():
    """
    Check if memory usage is critical. Raise error if >95%.
    
    Returns
    -------
    float
        Current memory usage percentage
    """
    slurm_mem, slurm_src = get_slurm_memory_bytes()

    if slurm_mem is not None:
        # Inside an allocation the node's usage is the wrong measure in both
        # directions: a busy shared node would abort a job sitting well within
        # its own limit, while a job overrunning its cgroup would look fine on
        # a large node right up until it is OOM-killed. Measure this process
        # against what it was actually given.
        used = psutil.Process().memory_info().rss
        percent_used = 100.0 * used / slurm_mem
        detail = (f"{used / (1024**3):.2f} GB of "
                  f"{slurm_mem / (1024**3):.2f} GB allocated [{slurm_src}]")
    else:
        mem = psutil.virtual_memory()
        percent_used = mem.percent
        detail = f"Available: {mem.available / (1024**3):.2f} GB"

    if percent_used > 95:
        raise MemoryError(
            f"CRITICAL: Memory usage at {percent_used:.1f}% - stopping to prevent crash. "
            f"{detail}"
        )

    if percent_used > 85:
        log(f"WARNING: Memory usage at {percent_used:.1f}% - approaching limit ({detail})")

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
    
    # Conditionally print to console. Flushed every line: piped into a SLURM
    # log, stdout is otherwise block-buffered and output arrives late.
    if ENABLE_PRINTING:
        if message == "":
            print(flush=True)
        else:
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"[{timestamp}] {message}", flush=True)


def broadening_sensitivity(sigma_threshold, outlier_fraction, two_tailed):
    """
    Smallest broadening factor this configuration will flag.

    A baseline whose scatter is b times the array norm puts a fraction of its
    samples beyond the sigma cut; it is flagged once that fraction exceeds
    --outlier-fraction-threshold. Inverting that gives the b at which the
    criterion starts to bite, so the settings can be read as a sensitivity
    rather than as two abstract numbers.
    """
    tail = outlier_fraction / 2.0 if two_tailed else outlier_fraction
    if not 0.0 < tail < 0.5:
        return float('nan')
    z = NormalDist().inv_cdf(1.0 - tail)
    return sigma_threshold / z if z > 0 else float('nan')


def wrap_angle(delta):
    """
    Signed angular difference(s) folded into (-pi, pi].

    Phase differences are only meaningful modulo 2*pi: a point at +179 deg and
    a centre at -179 deg are 2 deg apart, not 358. Folding here keeps that
    separation, so the outlier tests measure the real excursion.
    """
    return np.arctan2(np.sin(delta), np.cos(delta))


def circular_weighted_mean(values, weights):
    """Weighted mean direction of angles, in (-pi, pi]."""
    return float(np.arctan2(np.sum(weights * np.sin(values)),
                            np.sum(weights * np.cos(values))))


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


def scrub_name(name):
    """Make a field name safe to embed in a filename."""
    return re.sub(r'[^A-Za-z0-9._+-]', '_', str(name))


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


def nearest_scans(scan_ids, idx, window_count):
    """
    The `window_count` nearest scans to scan_ids[idx], by rank in the
    time-ordered scan_ids list (not raw scan number, which need not be evenly
    spaced), excluding scan_ids[idx] itself. Expands outward one side at a
    time so an edge scan still gets up to `window_count` neighbours from
    whichever side has them, rather than being limited by the shorter side.
    """
    n = len(scan_ids)
    lo, hi = idx - 1, idx + 1
    neighbours = []
    while len(neighbours) < window_count and (lo >= 0 or hi < n):
        if lo >= 0:
            neighbours.append(scan_ids[lo])
            lo -= 1
        if len(neighbours) >= window_count:
            break
        if hi < n:
            neighbours.append(scan_ids[hi])
            hi += 1
    return neighbours


def local_scan_reference(scan_centers, scan_mads, ensemble_mad, window_count,
                          sigma_threshold, effective_mode, circular=False):
    """
    For each scan, a reference amplitude built from its `window_count`
    nearest scans in time, each neighbour weighted by 1/MAD^2 so a noisier
    neighbour counts for less. A scan whose own weighted median deviates from
    that reference by more than `sigma_threshold` MAD-scaled sigma is not a
    reliable centre to test its own baselines against -- a scan-wide shift
    (bad solution, gain drop) moves every baseline's residual together, so
    nothing in that scan looks anomalous against its own equally-shifted
    median.

    With `circular` the centres are angles: the reference is a mean direction
    and the deviations are folded into (-pi, pi].

    Returns
    -------
    dict
        {scan_id: local_reference_amplitude}, one entry per scan whose own
        centre should be replaced by this reference. Scans not replaced are
        absent from the returned dict.
    """
    scan_ids = sorted(scan_centers.keys())
    replacements = {}

    if window_count <= 0 or len(scan_ids) < 2:
        return replacements

    for idx, scan_id in enumerate(scan_ids):
        neighbours = nearest_scans(scan_ids, idx, window_count)

        weights, values = [], []
        for nb in neighbours:
            nb_mad = scan_mads.get(nb, np.nan)
            if not (np.isfinite(nb_mad) and nb_mad > 0):
                nb_mad = ensemble_mad
            if not (np.isfinite(nb_mad) and nb_mad > 0):
                continue
            weights.append(1.0 / (nb_mad ** 2))
            values.append(scan_centers[nb])

        if not values:
            continue

        weights = np.array(weights)
        values = np.array(values)
        if circular:
            local_mean = circular_weighted_mean(values, weights)
            deviations = np.abs(wrap_angle(values - local_mean))
        else:
            local_mean = float(np.sum(weights * values) / np.sum(weights))
            deviations = np.abs(values - local_mean)

        local_mad = float(weighted_median(deviations, weights))
        if not (np.isfinite(local_mad) and local_mad > 0):
            local_mad = ensemble_mad
        if not (np.isfinite(local_mad) and local_mad > 0):
            continue

        threshold = sigma_threshold * 1.4826 * local_mad
        this_value = scan_centers[scan_id]
        offset = wrap_angle(this_value - local_mean) if circular else this_value - local_mean
        if effective_mode == 'high':
            hit = offset > threshold
        elif effective_mode == 'low':
            hit = -offset > threshold
        else:
            hit = abs(offset) > threshold

        if hit:
            replacements[scan_id] = local_mean

    return replacements


class ChunkCache:
    """
    The field's amplitudes (and phases, for the phase pass), read from the MS
    once and laid out on disk by channel chunk, so every pass over a chunk
    reads only that chunk. Flagged samples are NaN. Weights (one per row) and
    scan numbers are held in memory.
    """

    def __init__(self, cache_dir, chunk_starts, chunk_size, n_chan, corr_names,
                 memmaps, weights, scans, times, rows_per_bl, unflagged):
        self.cache_dir = cache_dir
        self.chunk_index = {s: k for k, s in enumerate(chunk_starts)}
        self.chunk_width = {s: min(s + chunk_size, n_chan) - s for s in chunk_starts}
        self.corr_slot = {name: j for j, name in enumerate(corr_names)}
        self.memmaps = memmaps
        self.weights = weights
        self.rows_per_bl = rows_per_bl
        self.unflagged = unflagged
        self.scan_arrays = [scans[i, :rows_per_bl[i]] for i in range(len(rows_per_bl))]
        self.time_arrays = [times[i, :rows_per_bl[i]] for i in range(len(rows_per_bl))]

    def has_data(self, start_chan):
        return self.unflagged[self.chunk_index[start_chan]] > 0

    def pass_arrays(self, corr_name, quantity, start_chan):
        """
        Per-baseline (row, chan) value arrays for one pass over one chunk, and
        the matching weights broadcast across the channel axis.
        """
        k = self.chunk_index[start_chan]
        w = self.chunk_width[start_chan]
        j = self.corr_slot[corr_name]
        mm = self.memmaps[(corr_name, quantity)]
        values, weights = [], []
        for i, rows in enumerate(self.rows_per_bl):
            values.append(np.asarray(mm[k, i, :rows, :w]))
            weights.append(np.broadcast_to(self.weights[i, :rows, j][:, None], (rows, w)))
        return values, weights

    def cleanup(self):
        self.memmaps = {}
        shutil.rmtree(self.cache_dir, ignore_errors=True)


def build_chunk_cache(ms_file, field_id, data_column, baseline_info, correlations,
                      phase_corr_names, chunk_starts, chunk_size, n_chan, cache_dir,
                      memory_budget_bytes):
    """
    Read the field once, in contiguous row blocks (the MS's own row order, so
    each storage tile is read once rather than once per baseline), sort the
    rows into baselines, and write each pass's values into per-chunk
    on-disk arrays. Within a baseline, rows keep the MS's time order.

    Returns
    -------
    ChunkCache
    """
    ds_all = xds_from_ms(
        ms_file,
        columns=[data_column, "FLAG", "WEIGHT", "ANTENNA1", "ANTENNA2", "SCAN_NUMBER", "TIME"],
        taql_where=f"FIELD_ID=={field_id}",
        group_cols=["DATA_DESC_ID"],
    )
    if len(ds_all) != 1:
        raise RuntimeError(f"Expected one spectral window for this field, found {len(ds_all)}")
    ds = ds_all[0].rename({data_column: "DATA"})

    n_bl = len(baseline_info)
    ant1, ant2 = dask.compute(ds.ANTENNA1.data, ds.ANTENNA2.data)
    n_ant = int(max(ant1.max(), ant2.max(), max(max(b[0], b[1]) for b in baseline_info))) + 1
    lookup = np.full((n_ant, n_ant), -1, dtype=np.int64)
    for i, (a1, a2, _, _, _) in enumerate(baseline_info):
        lookup[a1, a2] = i
    row_bl = lookup[ant1, ant2]
    rows_per_bl = np.bincount(row_bl[row_bl >= 0], minlength=n_bl)
    max_rows = int(rows_per_bl.max())

    lo = min(chunk_starts)
    hi = min(max(chunk_starts) + chunk_size, n_chan)
    n_chunks = len(chunk_starts)
    corr_names = [name for name, _ in correlations]
    corr_ms_idx = [idx for _, idx in correlations]
    quantities = [(name, 'amplitude') for name in corr_names]
    quantities += [(name, 'phase') for name in corr_names if name in phase_corr_names]

    need_bytes = len(quantities) * n_chunks * n_bl * max_rows * chunk_size * 4
    os.makedirs(cache_dir, exist_ok=True)
    free_bytes = shutil.disk_usage(cache_dir).free
    log(f"Chunk cache: {len(quantities)} pass(es) x {n_chunks} chunk(s), "
        f"{need_bytes / 1024**3:.2f} GB on disk in {cache_dir} "
        f"({free_bytes / 1024**3:.1f} GB free)")
    if need_bytes > 0.9 * free_bytes:
        raise RuntimeError(f"Not enough disk space for the chunk cache in {cache_dir}: "
                           f"need {need_bytes / 1024**3:.2f} GB, "
                           f"{free_bytes / 1024**3:.1f} GB free (use --cache-dir)")

    memmaps = {
        q: np.memmap(os.path.join(cache_dir, f"{q[0]}_{q[1]}.f32"), dtype=np.float32,
                     mode='w+', shape=(n_chunks, n_bl, max_rows, chunk_size))
        for q in quantities}
    weights = np.full((n_bl, max_rows, len(correlations)), np.nan, dtype=np.float32)
    scans = np.zeros((n_bl, max_rows), dtype=np.int64)
    times = np.zeros((n_bl, max_rows), dtype=np.float64)
    unflagged = np.zeros(n_chunks, dtype=np.int64)
    counters = np.zeros(n_bl, dtype=np.int64)

    # Several working copies of a row block are alive at once (complex data,
    # flags, one float32 pass), so the budget is divided generously.
    n_corr_ms = ds.dims['corr']
    bytes_per_row = (hi - lo) * n_corr_ms * (8 + 1) * 3
    row_block = int(np.clip(memory_budget_bytes / bytes_per_row, 500, 200000))
    n_rows = len(row_bl)
    n_blocks = (n_rows + row_block - 1) // row_block
    log(f"Reading {n_rows:,} rows, channels {lo}-{hi - 1}, in {n_blocks} block(s) "
        f"of up to {row_block:,} rows")

    t_start = datetime.now()
    for b, r0 in enumerate(range(0, n_rows, row_block)):
        r1 = min(r0 + row_block, n_rows)
        blk_bl = row_bl[r0:r1]
        keep = blk_bl >= 0
        if not np.any(keep):
            continue

        sub = ds.isel(row=slice(r0, r1))
        data, flag, wt, scan, time = dask.compute(sub.DATA.data[:, lo:hi, :],
                                                  sub.FLAG.data[:, lo:hi, :],
                                                  sub.WEIGHT.data, sub.SCAN_NUMBER.data,
                                                  sub.TIME.data)
        data, flag, wt, scan, time, blk_bl = (data[keep], flag[keep], wt[keep], scan[keep],
                                              time[keep], blk_bl[keep])

        # Position of each row within its baseline: rows already placed for
        # that baseline, plus its rank among this block's rows of the same
        # baseline (the stable sort keeps them in MS order).
        order = np.argsort(blk_bl, kind='stable')
        sorted_bl = blk_bl[order]
        rank = np.empty_like(order)
        rank[order] = np.arange(len(order)) - np.searchsorted(sorted_bl, sorted_bl, side='left')
        row_pos = counters[blk_bl] + rank
        counters += np.bincount(blk_bl, minlength=n_bl)

        scans[blk_bl, row_pos] = scan
        times[blk_bl, row_pos] = time
        weights[blk_bl, row_pos, :] = wt[:, corr_ms_idx]

        for (name, quantity), mm in memmaps.items():
            ms_idx = corr_ms_idx[corr_names.index(name)]
            vis = data[:, :, ms_idx]
            values = (np.abs(vis) if quantity == 'amplitude'
                      else np.arctan2(vis.imag, vis.real)).astype(np.float32)
            values[flag[:, :, ms_idx]] = np.nan
            for k, s in enumerate(chunk_starts):
                e = min(s + chunk_size, n_chan)
                mm[k][blk_bl, row_pos, :e - s] = values[:, s - lo:e - lo]

        sel_flag = flag[:, :, corr_ms_idx]
        for k, s in enumerate(chunk_starts):
            e = min(s + chunk_size, n_chan)
            unflagged[k] += int((~sel_flag[:, s - lo:e - lo, :]).sum())

        elapsed = (datetime.now() - t_start).total_seconds()
        log(f"  Block {b + 1}/{n_blocks} read ({r1:,}/{n_rows:,} rows, {elapsed / 60:.1f} min)")

    for mm in memmaps.values():
        mm.flush()

    return ChunkCache(cache_dir, chunk_starts, chunk_size, n_chan, corr_names,
                      memmaps, weights, scans, times, rows_per_bl, unflagged)


def flatten_valid(value_arrays, weight_arrays, *row_arrays):
    """
    Flatten per-baseline (row, chan) value and weight arrays into 1-D arrays
    of the samples where both are finite, alongside each sample's baseline
    index and its value in each per-row array (scan number, time window).

    Returns
    -------
    tuple of np.ndarray
        (baseline_idx, values, weights, *row values); all empty if nothing is valid.
    """
    bl_list, val_list, wt_list = [], [], []
    row_lists = [[] for _ in row_arrays]
    for i, (val_arr, wt_arr) in enumerate(zip(value_arrays, weight_arrays)):
        val_flat = val_arr.flatten()
        wt_flat = wt_arr.flatten()
        valid = np.isfinite(val_flat) & np.isfinite(wt_flat)
        if not np.any(valid):
            continue

        # (row, chan) flattens row-major, so repeating each row's value
        # n_chan times lines the row axis up exactly.
        n_chan_this = val_arr.shape[1] if val_arr.ndim > 1 else 1
        for row_list, arrays in zip(row_lists, row_arrays):
            row_list.append(np.repeat(arrays[i], n_chan_this)[valid])

        bl_list.append(np.full(valid.sum(), i, dtype=np.int32))
        val_list.append(val_flat[valid])
        wt_list.append(wt_flat[valid])

    if not val_list:
        empty = np.array([])
        return (np.array([], dtype=np.int32), empty, empty) + tuple(empty for _ in row_arrays)
    return ((np.concatenate(bl_list), np.concatenate(val_list), np.concatenate(wt_list))
            + tuple(np.concatenate(row_list) for row_list in row_lists))


def casa_time(mjd_seconds):
    """MS TIME (MJD seconds, UTC) in flagdata's timerange format."""
    t = datetime(1858, 11, 17) + timedelta(seconds=float(mjd_seconds))
    return t.strftime('%Y/%m/%d/%H:%M:%S.') + f'{t.microsecond // 1000:03d}'


def time_windows(scan_arrays, time_arrays, time_chunk):
    """
    Split every scan into equal time windows of roughly `time_chunk` seconds.

    A scan is cut into the whole number of windows closest to `time_chunk`
    long, so a 12-minute scan becomes two 6-minute windows rather than 5 + 5 + 2.
    Window edges fall between integrations, so each row is in exactly one window.

    Returns
    -------
    window_arrays : list of np.ndarray
        Window number of every row, per baseline (all 0 when time_chunk <= 0).
    window_ranges : dict
        {(scan, window): flagdata timerange}; empty when time_chunk <= 0.
    """
    window_arrays = [np.zeros(len(s), dtype=np.int64) for s in scan_arrays]
    if time_chunk <= 0 or not any(len(s) for s in scan_arrays):
        return window_arrays, {}

    scans = np.concatenate(scan_arrays)
    times = np.concatenate(time_arrays)
    unique_times = np.unique(times)
    dt = float(np.median(np.diff(unique_times))) if len(unique_times) > 1 else 0.0

    edges = {}  # {scan: window edge times}
    for scan in np.unique(scans):
        t = times[scans == scan]
        start, length = t.min() - dt / 2, t.max() - t.min() + dt
        n = max(1, int(round(length / time_chunk)))
        steps = np.arange(n + 1) * length / n
        edges[int(scan)] = start + (np.round(steps / dt) * dt if dt > 0 else steps)

    for windows, scan_arr, time_arr in zip(window_arrays, scan_arrays, time_arrays):
        for scan, e in edges.items():
            in_scan = scan_arr == scan
            windows[in_scan] = np.clip(np.searchsorted(e, time_arr[in_scan], side='right') - 1,
                                       0, len(e) - 2)

    window_ranges = {(scan, k): f"{casa_time(e[k])}~{casa_time(e[k + 1])}"
                     for scan, e in edges.items() for k in range(len(e) - 1)}
    return window_arrays, window_ranges


def resolve_protected_antenna(value, antenna_names):
    """
    Name of the antenna that is never flagged as a whole, or None.

    'auto' takes the first antenna of ref_ant in ./project_info.json (an index
    or a name); 'none' turns protection off.
    """
    if value == 'none':
        return None
    if value == 'auto':
        try:
            with open('project_info.json') as f:
                value = str(json.load(f)['ref_ant']).split(',')[0].strip()
        except (OSError, KeyError, ValueError):
            return None
    return antenna_names[int(value)] if value.isdigit() else value


def compute_scan_stats(values, weights, scans, circular=False):
    """
    Per-scan weighted centre and weighted MAD about that centre.

    With `circular` the values are angles: the centre is the median offset
    from a provisional mean direction, folded back into (-pi, pi], and the
    deviations are wrapped the same way.

    Returns
    -------
    tuple of dict
        ({scan: centre}, {scan: mad})
    """
    scan_centers = {}
    scan_mads = {}
    for scan_id in np.unique(scans):
        in_scan = scans == scan_id
        if circular:
            provisional = circular_weighted_mean(values[in_scan], weights[in_scan])
            offset = weighted_median(wrap_angle(values[in_scan] - provisional),
                                     weights[in_scan])
            centre = float(wrap_angle(provisional + offset))
            deviations = np.abs(wrap_angle(values[in_scan] - centre))
        else:
            centre = weighted_median(values[in_scan], weights[in_scan])
            deviations = np.abs(values[in_scan] - centre)
        scan_centers[scan_id] = centre
        scan_mads[scan_id] = weighted_median(deviations, weights[in_scan])
    return scan_centers, scan_mads


def robust_poly_fit(x, y, degree, clip_sigma, max_iter=10):
    """
    Polynomial fit of y(x) that iteratively drops points lying more than
    `clip_sigma` (MAD-scaled) ABOVE the current fit. The clipping is
    one-sided because contamination (RFI) only ever pushes the statistics
    being fitted up; a genuinely low point is kept.

    The degree is reduced when there are too few points to constrain it.

    Returns
    -------
    tuple or None
        (np.polynomial.Polynomial, used_mask), or None with fewer than 3 points.
    """
    if len(x) < 3:
        return None
    mask = np.ones(len(x), dtype=bool)
    poly = None
    for _ in range(max_iter):
        deg = min(degree, int(mask.sum()) - 2)
        if deg < 0:
            break
        poly = np.polynomial.Polynomial.fit(x[mask], y[mask], deg)
        residuals = y - poly(x)
        res_centre = np.median(residuals[mask])
        res_mad = np.median(np.abs(residuals[mask] - res_centre))
        if not res_mad > 0:
            break
        new_mask = (residuals - res_centre) <= clip_sigma * 1.4826 * res_mad
        if np.array_equal(new_mask, mask) or new_mask.sum() < 3:
            break
        mask = new_mask
    if poly is None:
        return None
    return poly, mask


def collect_band_statistics(cache, chunks_to_process, correlations):
    """
    Statistics-only pass over every chunk: the per-scan amplitude centre and
    MAD for each correlation, computed exactly as the flagging pass computes
    them. No flags are decided here.

    Returns
    -------
    dict
        {corr_name: {scan: {start_chan: (centre, mad)}}}
    """
    band_stats = {corr_name: {} for corr_name, _ in correlations}
    for start_chan in chunks_to_process:
        if not cache.has_data(start_chan):
            continue
        for corr_name, _ in correlations:
            value_arrays, weight_arrays = cache.pass_arrays(corr_name, 'amplitude', start_chan)
            _, values, weights, scans = flatten_valid(value_arrays, weight_arrays,
                                                      cache.scan_arrays)
            if values.size == 0:
                continue
            scan_centers, scan_mads = compute_scan_stats(values, weights, scans)
            for scan_id in scan_centers:
                band_stats[corr_name].setdefault(scan_id, {})[start_chan] = (
                    scan_centers[scan_id], scan_mads[scan_id])
    return band_stats


def fit_band_reference(band_stats, chunk_freqs, degree, clip_sigma):
    """
    For each (correlation, scan), fit log(centre) and log(MAD) against chunk
    frequency with robust_poly_fit, and evaluate both fits at every chunk.

    Returns
    -------
    dict
        {corr_name: {scan: {start_chan: (fitted_centre, fitted_mad)}}}; a
        (correlation, scan) whose fit could not be made is absent.
    """
    band_ref = {}
    all_starts = sorted(chunk_freqs)
    all_x = np.array([chunk_freqs[s] for s in all_starts])

    for corr_name, per_scan in band_stats.items():
        for scan_id, per_chunk in sorted(per_scan.items()):
            starts = [s for s in sorted(per_chunk)
                      if all(np.isfinite(v) and v > 0 for v in per_chunk[s])]
            if len(starts) < 3:
                log(f"  {corr_name} scan {int(scan_id)}: only {len(starts)} usable chunk(s) "
                    f"- no band fit")
                continue
            x = np.array([chunk_freqs[s] for s in starts])
            log_centre = np.log([per_chunk[s][0] for s in starts])
            log_mad = np.log([per_chunk[s][1] for s in starts])

            centre_fit = robust_poly_fit(x, log_centre, degree, clip_sigma)
            mad_fit = robust_poly_fit(x, log_mad, degree, clip_sigma)
            if centre_fit is None or mad_fit is None:
                log(f"  {corr_name} scan {int(scan_id)}: band fit failed")
                continue

            centre_poly, centre_used = centre_fit
            mad_poly, mad_used = mad_fit
            fitted_centre = np.exp(centre_poly(all_x))
            fitted_mad = np.exp(mad_poly(all_x))
            band_ref.setdefault(corr_name, {})[scan_id] = {
                s: (float(c), float(m)) for s, c, m in zip(all_starts, fitted_centre, fitted_mad)}

            dropped = sorted(set(np.array(starts)[~centre_used]) | set(np.array(starts)[~mad_used]))
            log(f"  {corr_name} scan {int(scan_id)}: centre fit on {centre_used.sum()}/{len(starts)}, "
                f"MAD fit on {mad_used.sum()}/{len(starts)} chunk(s); dropped chunk start "
                f"channel(s): {', '.join(str(int(s)) for s in dropped) or 'none'}")

    return band_ref


def apply_band_reference(scan_centers, scan_mads, chunk_ref, mad_cap, centre_sigma):
    """
    Compare each scan's own centre and MAD in this chunk with the band fit.
    A scan is contaminated when its MAD exceeds `mad_cap` times the fitted
    MAD, or its centre lies more than `centre_sigma` fitted-MAD sigma above
    the fitted centre. Only an excess counts: RFI only adds power.

    Returns
    -------
    dict
        {scan: (fitted_centre, fitted_mad, reason)} for each contaminated
        scan, whose statistics are then replaced by the fitted values.
    """
    replacements = {}
    for scan_id, centre in scan_centers.items():
        ref = chunk_ref.get(scan_id)
        if ref is None:
            continue
        fit_centre, fit_mad = ref
        own_mad = scan_mads.get(scan_id, np.nan)

        reasons = []
        if np.isfinite(own_mad) and own_mad > mad_cap * fit_mad:
            reasons.append(f"MAD {own_mad:.4f} > {mad_cap}x fit {fit_mad:.4f}")
        if np.isfinite(centre) and centre - fit_centre > centre_sigma * 1.4826 * fit_mad:
            reasons.append(f"centre {centre:.4f} > fit {fit_centre:.4f} + "
                           f"{centre_sigma}sigma")
        if reasons:
            replacements[scan_id] = (fit_centre, fit_mad, '; '.join(reasons))
    return replacements


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


def resolve_field_list(ms_file, fields_arg):
    """Resolve a comma-separated field list to a list of (field_id, field_name)."""
    return [resolve_field(ms_file, f.strip()) for f in fields_arg.split(',') if f.strip()]


def load_scan_order(ms_file):
    """
    Every scan in the MS in time order, with the field it belongs to.

    Returns
    -------
    scan_ids, field_ids : ndarray
        Parallel arrays, ordered by each scan's first timestamp.
    """
    main = xds_from_table(ms_file, columns=["SCAN_NUMBER", "FIELD_ID", "TIME"])[0]
    scan, field, time = dask.compute(main.SCAN_NUMBER.data, main.FIELD_ID.data,
                                     main.TIME.data)
    per_scan = (pd.DataFrame({'scan': scan, 'field': field, 'time': time})
                .groupby('scan')
                .agg(field=('field', 'first'), start=('time', 'min'))
                .sort_values('start'))
    return per_scan.index.values, per_scan['field'].values


def propagation_scans(scan_ids, field_ids, check_field_id, flag_field_ids):
    """
    For each scan of the check field, the other scans its flags are passed to.

    A check-field scan's flags go to every scan of a field in `flag_field_ids`
    lying between the previous check-field scan and the next one, in time
    across the whole MS. With no previous (next) check-field scan the span runs
    to the start (end) of the observation, so a check field with a single scan
    passes its flags to every scan of those fields. For example, with scans
    P1 T2 T3 P4 T5, a flag on P1 goes to T2 and T3, and a flag on P4 goes to
    T2, T3 and T5. Scans of any other field are never included.

    Returns
    -------
    dict
        {check_scan: sorted list of other scan numbers}
    """
    flag_field_ids = set(flag_field_ids) - {check_field_id}
    check_idx = [i for i, f in enumerate(field_ids) if f == check_field_id]
    bounds = [-1] + check_idx + [len(scan_ids)]

    extras = {}
    for k, idx in enumerate(check_idx, start=1):
        lo, hi = bounds[k - 1], bounds[k + 1]
        extras[int(scan_ids[idx])] = sorted(
            int(scan_ids[j]) for j in range(lo + 1, hi)
            if field_ids[j] in flag_field_ids)
    return extras


def scan_selection(scans):
    """flagdata scan selection string for a collection of scan numbers."""
    return ','.join(str(s) for s in sorted({int(s) for s in scans}))


def load_grouped_ms(ms_file, field_id, data_column):
    """
    Load MS grouped by baseline (ANTENNA1, ANTENNA2, DATA_DESC_ID).
    Returns a list of datasets, one per baseline.
    """
    # SCAN_NUMBER is a column rather than a group column, so each row can be
    # attributed to its scan for the per-scan statistics while every baseline
    # stays whole. Grouping on it instead would split each baseline per scan
    # and change what a flag applies to.
    cols = [data_column, "FLAG", "WEIGHT", "ANTENNA1", "ANTENNA2", "SCAN_NUMBER"]
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
                         outlier_mode, save_dir, run_tag, quantity='amplitude'):
    """
    Create histogram of the tested quantity showing median, MAD, and thresholds.
    """
    log("Generating histogram diagnostic plot...")

    axis_label = 'Phase (rad)' if quantity == 'phase' else 'Amplitude'

    amp_values = df['value'].values
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
    
    ax.set_xlabel(axis_label, fontsize=12)
    ax.set_ylabel('Count (log scale)', fontsize=12)
    ax.set_title(f'{quantity.capitalize()} Distribution - Chunk {chunk_num} - {corr_name}',
                 fontsize=14)
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save histogram
    os.makedirs(save_dir, exist_ok=True)
    prefix = 'phase' if quantity == 'phase' else 'amp'
    output_path = os.path.join(save_dir,
                               f'{prefix}_histogram_{run_tag}_chunk{chunk_num}_{corr_name}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    log(f"Histogram saved to {output_path}")
    plt.close()


def create_plot(df, baseline_names, median_amp, mad_amp, threshold_mad,
                flagged_baselines, start_chan, end_chan, chunk_num, corr_name,
                outlier_mode, save_dir, field, run_tag, quantity='amplitude'):
    """
    Create a quantity vs baseline plot with datashader and matplotlib.
    Two subplots vertically: upper colored by antenna, lower colored by weight.
    `flagged_baselines` is the caller's final per-(baseline, scan) decision
    (fraction threshold and/or point-level override), indexed by
    (baseline_idx, scan); every baseline with >=1 flagged scan gets a red star.
    """
    log("Generating plot with datashader...")

    axis_label = 'Phase (rad)' if quantity == 'phase' else 'Amplitude'
    title_label = quantity.capitalize()

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
    amp_min = df['value'].min()
    amp_max = df['value'].max()
    
    # ========== UPPER PLOT: Colored by antenna ==========
    canvas1 = dsh.Canvas(plot_width=canvas_width, plot_height=canvas_height,
                        x_range=(-0.5, n_baselines - 0.5),
                        y_range=(amp_min, amp_max))
    
    agg1 = canvas1.points(df, 'baseline_idx', 'value', dsh.count_cat('ant1'))
    img1 = tf.shade(agg1, color_key=color_key_ant, how='linear', min_alpha=60)
    img1 = tf.dynspread(img1, threshold=0.5, max_px=4)
    img1 = tf.set_background(img1, 'white')
    
    ax1.imshow(img1.to_pil(), extent=[-0.5, n_baselines - 0.5, amp_min, amp_max],
               aspect='auto', origin='upper')
    
    # This plot has no scan axis, so collapse to unique baseline_idx: a
    # baseline with several bad scans gets one star.
    flagged_bl_idxs = (sorted(set(flagged_baselines.index.get_level_values('baseline_idx')))
                        if len(flagged_baselines) > 0 else [])

    if len(flagged_bl_idxs) > 0:
        n_total_baselines = len(baseline_names)
        baseline_loss_fraction = len(flagged_bl_idxs) / n_total_baselines
        log(f"Marking {len(flagged_bl_idxs)} flagged baselines "
            f"({len(flagged_baselines)} baseline-scan flag(s))")
        log(f"Flagging these baselines would remove {baseline_loss_fraction:.4f} ({len(flagged_bl_idxs)}/{n_total_baselines}) of all baselines")

        # Mark 99.9th percentile max for each bad baseline with red stars (top outliers only)
        for bl_idx in flagged_bl_idxs:
            bl_data = df[df['baseline_idx'] == bl_idx]
            if len(bl_data) > 0:
                amp_max = bl_data['value'].quantile(0.999)
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
    
    ax1.set_ylabel(axis_label, fontsize=12)
    ax1.set_title(f'{title_label} vs Baseline (colored by antenna) - Field {field}, Channels {start_chan}-{end_chan}, {corr_name}',
                  fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=9)
    
    # ========== LOWER PLOT: Colored by weight ==========
    canvas2 = dsh.Canvas(plot_width=canvas_width, plot_height=canvas_height,
                        x_range=(-0.5, n_baselines - 0.5),
                        y_range=(df['value'].min(), df['value'].max()))
    
    # Use mean aggregation for continuous weight values
    agg2 = canvas2.points(df, 'baseline_idx', 'value', dsh.mean('weight'))
    
    # Use span to set color range for weights
    weight_span = (df['weight'].quantile(0.01), df['weight'].quantile(0.99))
    img2 = tf.shade(agg2, cmap=weight_cmap, how='linear', span=weight_span, min_alpha=60)
    img2 = tf.dynspread(img2, threshold=0.5, max_px=4)
    img2 = tf.set_background(img2, 'white')
    
    ax2.imshow(img2.to_pil(), extent=[-0.5, n_baselines - 0.5, 
                                       df['value'].min(), df['value'].max()],
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
    ax2.set_ylabel(axis_label, fontsize=12)
    ax2.set_title(f'{title_label} vs Baseline (colored by weight - blue is high weight)', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Create output directory and save
    os.makedirs(save_dir, exist_ok=True)
    prefix = 'phase' if quantity == 'phase' else 'amp'
    output_path = os.path.join(save_dir,
                               f'{prefix}_vs_baseline_{run_tag}_chunk{chunk_num}_{corr_name}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    log(f"Plot saved to {output_path}")
    plt.close()


def main():
    """Main processing routine."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Set global variable for console output
    global ENABLE_PRINTING
    ENABLE_PRINTING = not args.no_printing

    # Resolve the field before logging is configured, so the log can be named
    # after it. resolve_field() raises rather than logging, so nothing is lost.
    field_id, field_name = resolve_field(args.ms_file, args.field)
    flag_fields = resolve_field_list(args.ms_file, args.flag_fields)

    # Every run gets its own field-named, timestamped log and plots, so repeated
    # runs -- and the five separate call sites in 1GC_05_casa_refcal.py -- never
    # overwrite each other's output. The flag file is the deliberate exception:
    # it keeps a predictable name because flagdata() is handed it by bare
    # filename, and an existing one is rotated to a timestamped backup instead.
    run_stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_tag = f'{scrub_name(field_name)}_{run_stamp}'

    os.makedirs(args.log_dir, exist_ok=True)
    log_filename = os.path.join(args.log_dir,
                                f'baseline-flagger_{run_tag}{args.output_suffix}.log')
    logging.basicConfig(
        filename=log_filename,
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filemode='w',
        force=True  # Reconfigure if already configured
    )

    # =================================================================
    # CONFIGURATION SUMMARY
    # =================================================================
    log("="*70)
    log("BASEDFLAGGER CONFIGURATION")
    log("="*70)
    log("")
    
    log("INPUT DATA:")
    log(f"  MS file: {args.ms_file}")
    log(f"  Data column: {args.data_column}")
    log(f"  Field: {field_name} (ID: {field_id})")
    log(f"  Flags passed on to fields: "
        f"{', '.join(name for _, name in flag_fields) or 'none'}")
    log("")
    
    log("PROCESSING PARAMETERS:")
    log(f"  Channel group size: {args.channel_group_size}")
    log(f"  Number of chunks: {args.num_chunks if args.num_chunks is not None else 'All (process entire band)'}")
    log(f"  Target frequency: {args.freq_target/1e9:.4f} GHz" if args.freq_target else "  Target frequency: Not specified (start from channel 0)")
    log(f"  Correlation products: {args.correlation_products}")
    log(f"  Phase pass: {'ENABLED' if args.phase else 'DISABLED'}")
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
    point_override_enabled = args.point_outlier_min_count > 0 and args.point_sigma_threshold > 0
    if point_override_enabled:
        log(f"  Point-level override: flag baseline outright if >= "
            f"{args.point_outlier_min_count} point(s) exceed {args.point_sigma_threshold}σ")
    else:
        log(f"  Point-level override: DISABLED")
    if args.scan_window_count > 0:
        log(f"  Scan-window centre correction: replace a scan's own centre with its "
            f"{args.scan_window_count}-nearest-scan reference if it's >"
            f"{args.scan_window_sigma}σ off")
    else:
        log(f"  Scan-window centre correction: DISABLED")
    if args.band_reference == 'poly':
        log(f"  Band reference: degree-{args.band_poly_degree} polynomial across the band "
            f"(clip {args.band_clip_sigma}σ above fit); a scan in a chunk is tested against "
            f"the fit if its MAD > {args.band_mad_cap}x fit or its centre > fit + "
            f"{args.band_centre_sigma}σ (amplitude only)")
    else:
        log(f"  Band reference: DISABLED")
    log("")

    log("ANTENNA FLAGGING:")
    log(f"  Method 1 - Absolute threshold: {args.antenna_flag_threshold*100:.0f}%")
    log(f"    → Flag antenna in a scan if >{args.antenna_flag_threshold*100:.0f}% of its baselines are bad in that scan")
    if args.antenna_relative_sigma <= 0:
        log(f"  Method 2 - Relative outlier detection: DISABLED")
    else:
        log(f"  Method 2 - Relative outlier detection: {args.antenna_relative_sigma:.1f}σ")
        log(f"    → Flag antenna if it's >{args.antenna_relative_sigma:.1f}σ above population median")
    log(f"  Maximum antennas to flag globally per scan: {args.antenna_flag_cap}")
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

    # The reference antenna is never flagged as a whole (its baselines still can be)
    protected_antenna = resolve_protected_antenna(args.protect_antenna, antenna_names)
    log(f"Protected antenna (never flagged as a whole): {protected_antenna or 'none'}")
    protected_hits = []  # where it would have been flagged, for the warning at the end

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

    # The phase pass covers the parallel hands only; anything else selected is
    # tested on amplitude alone.
    if args.phase:
        phase_corrs = [name for name, _ in correlations if name in PHASE_CORRELATIONS]
        skipped_corrs = [name for name, _ in correlations if name not in PHASE_CORRELATIONS]
        if not phase_corrs:
            log("")
            log(f"WARNING: --phase given but none of {args.correlation_products} are "
                f"{'/'.join(PHASE_CORRELATIONS)} - no phase pass will run")
        else:
            log("")
            log(f"Phase pass enabled for: {', '.join(phase_corrs)}")
        if skipped_corrs:
            log(f"Phase pass NOT run for: {', '.join(skipped_corrs)} "
                f"(cross-hand phase tracks the polarisation angle and sweeps with "
                f"parallactic angle; amplitude is still tested)")
    
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
    log(f"Memory basis: {resources['memory_source']}")
    log(f"CPU cores: {resources['physical_cores']} physical, {resources['logical_cores']} logical on this node")
    log(f"CPUs allocated to this job: {resources['allocated_cpus']}  [{resources['cpu_source']}]")

    # Cap dask to the allocation. Its threaded scheduler otherwise defaults to
    # os.cpu_count(), which oversubscribes a SLURM job on a shared node.
    dask.config.set(num_workers=resources['allocated_cpus'])
    log(f"Dask capped to {resources['allocated_cpus']} thread(s) for parallel processing")
    
    # Estimate chunk memory and check if safe
    # All requested correlations are loaded together, so the estimate covers
    # every one of them
    n_baselines = len(ds_list)
    n_corr_processing = len(correlations)
    chunk_memory_gb = estimate_chunk_memory(n_baselines, n_time_samples,
                                            args.channel_group_size, n_corr_processing)
    is_safe, max_safe_gb = check_chunk_memory_safe(chunk_memory_gb, resources['available_memory_gb'])
    
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
    global_antenna_stats = {}  # {ant_name: {'scans': {scan: {'total': N, 'flagged': M}}, 'chunks': [...]}}
    
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
    n_passes = len(correlations) + (len([n for n, _ in correlations if n in PHASE_CORRELATIONS])
                                     if args.phase else 0)
    log(f"  Passes per chunk: {n_passes} (amplitude on every correlation"
        f"{', phase on the parallel hands' if args.phase else ''})")
    log(f"  Total operations: {len(chunks_to_process)} chunks × {n_passes} passes = {len(chunks_to_process) * n_passes} pass-chunks")
    log("")
    log(f"  Output files:")
    log(f"    - Flags: {flag_file_path}")
    log(f"    - Log: {log_filename}")
    if args.num_chunks is not None and args.num_chunks < 3:
        log(f"    - Plots: {args.save_dir}amp_vs_baseline_{run_tag}_chunk{{N}}_{{corr}}.png")
        log(f"    - Plots: {args.save_dir}amp_histogram_{run_tag}_chunk{{N}}_{{corr}}.png")
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

    # Read the field once into the per-chunk cache; every pass below reads
    # from it. Removed at exit (an abnormal kill can leave it behind).
    log("")
    log("="*60)
    log("BUILDING CHUNK CACHE")
    log("="*60)
    cache_dir = os.path.join(args.cache_dir,
                             f'basedflagger_cache_{run_tag}{args.output_suffix}')
    cache = build_chunk_cache(
        args.ms_file, field_id, args.data_column, baseline_info, correlations,
        PHASE_CORRELATIONS if args.phase else (), chunks_to_process,
        args.channel_group_size, n_chan, cache_dir,
        0.25 * resources['available_memory_gb'] * 1024**3)
    atexit.register(cache.cleanup)
    log("="*60)

    # Per-row scan numbers, in the cache's row order: SCAN_NUMBER has no
    # channel axis, so it is identical for every chunk and every correlation.
    scan_arrays = cache.scan_arrays
    all_scans = np.unique(np.concatenate([s for s in scan_arrays if s.size]))
    log(f"Found {len(all_scans)} scan(s) for this field: "
        f"{', '.join(str(int(s)) for s in all_scans)}")

    # Which other scans each of this field's flags is passed on to. Every flag
    # line is written with an explicit scan selection built from this, so a
    # flag lands only on this field's own scan(s) plus those other scans.
    ms_scan_ids, ms_field_ids = load_scan_order(args.ms_file)
    extra_scans = propagation_scans(ms_scan_ids, ms_field_ids, field_id,
                                    [fid for fid, _ in flag_fields])
    for scan in all_scans:
        extra_scans.setdefault(int(scan), [])
    if flag_fields:
        log("Flag propagation (scan of this field -> other scans flagged with it):")
        for scan in sorted(extra_scans):
            log(f"  {scan} -> {scan_selection(extra_scans[scan]) or 'none'}")

    # Targets are flagged in time windows of about --time-chunk seconds.
    # Calibrators keep whole scans, since their flags are copied onto the
    # target scans around them.
    time_chunk = 0 if flag_fields else args.time_chunk
    window_arrays, window_ranges = time_windows(scan_arrays, cache.time_arrays, time_chunk)
    all_windows = [tuple(int(v) for v in sw) for sw in np.unique(
        np.stack([np.concatenate(scan_arrays), np.concatenate(window_arrays)], axis=1), axis=0)]
    if window_ranges:
        log(f"Flagging in time windows of about {time_chunk:.0f} s: "
            f"{len(all_windows)} window(s) over {len(all_scans)} scan(s)")
    else:
        log("Flagging whole scans" + (" (calibrator: flags are copied to other fields)"
                                      if flag_fields else ""))

    # The protected antenna's baselines, and the rows of each newly covered
    # by a flag in each channel chunk, for its before/after report at the end
    protected_bls = {i for i, (a1, a2, _, _, _) in enumerate(baseline_info)
                     if protected_antenna in (antenna_names[a1], antenna_names[a2])}
    protected_new_rows = {}  # {(start_chan, bl_idx): bool array over that baseline's rows}

    def mark_protected(start_chans, bl_idx, scan, window):
        """Record rows of a protected-antenna baseline covered by a new flag."""
        rows = (scan_arrays[bl_idx] == scan) & (window_arrays[bl_idx] == window)
        for start_chan in start_chans:
            mask = protected_new_rows.setdefault(
                (start_chan, bl_idx), np.zeros(len(scan_arrays[bl_idx]), dtype=bool))
            mask |= rows

    def flag_selection(scan, window):
        """flagdata selection for a flag: the scan (or its time window) plus any propagated scans."""
        if window_ranges:
            return f"scan='{scan}' timerange='{window_ranges[(scan, window)]}'"
        return f"scan='{scan_selection([scan] + extra_scans.get(scan, []))}'"

    def window_label(scan, window):
        """Scan, plus the window number when flagging in time windows."""
        return f"{scan}" + (f"_WIN{window}" if window_ranges else "")

    # The scan-window reference compares a scan against its neighbours in time,
    # so a single-scan field has nothing to compare against. The rest of the
    # detection still works on one scan: its centre and MAD are simply that
    # scan's own, which is an ordinary single-population outlier test.
    if args.scan_window_count > 0 and len(all_scans) < 2:
        log(f"WARNING: --scan-window-count={args.scan_window_count} needs at least "
            f"2 scans to build a local reference; this field has {len(all_scans)}. "
            f"Scan-window check disabled for this run.")
        args.scan_window_count = 0

    # Band reference: a statistics-only pass over every chunk, then a smooth
    # fit across frequency per (correlation, scan). The flagging pass below
    # then knows the band-wide shape when it reaches each chunk.
    band_ref = None
    band_empty_chunks = {s for s in chunks_to_process if not cache.has_data(s)}
    if args.band_reference == 'poly':
        if len(chunks_to_process) < args.band_poly_degree + 3:
            log(f"WARNING: --band-reference poly needs at least {args.band_poly_degree + 3} "
                f"chunks to fit a degree-{args.band_poly_degree} polynomial; "
                f"{len(chunks_to_process)} to process. Band reference disabled for this run.")
        else:
            log("")
            log("="*60)
            log("BAND STATISTICS PASS")
            log("="*60)
            band_stats = collect_band_statistics(cache, chunks_to_process, correlations)
            chunk_freqs = {
                s: float(np.mean(chan_freq[s:min(s + args.channel_group_size, n_chan)])) / 1e9
                for s in chunks_to_process}
            log("")
            log("Band fits (log centre and log MAD vs frequency):")
            band_ref = fit_band_reference(band_stats, chunk_freqs,
                                          args.band_poly_degree, args.band_clip_sigma)

            # Preview of what the flagging pass will do with the fit:
            # which chunks have at least one scan tested against it.
            chunk_numbers = {s: i + 1 for i, s in enumerate(chunks_to_process)}
            log("")
            log("BAND REFERENCE SUMMARY (chunks where >=1 scan will be tested "
                "against the band fit):")
            for corr_name, _ in correlations:
                per_scan_ref = band_ref.get(corr_name, {})
                hit_chunks = {}
                for start in chunks_to_process:
                    own_centres = {}
                    own_mads = {}
                    for scan_id, per_chunk in band_stats.get(corr_name, {}).items():
                        if start in per_chunk:
                            own_centres[scan_id], own_mads[scan_id] = per_chunk[start]
                    chunk_ref = {scan_id: ref[start] for scan_id, ref in per_scan_ref.items()
                                 if start in ref}
                    replaced = apply_band_reference(own_centres, own_mads, chunk_ref,
                                                    args.band_mad_cap, args.band_centre_sigma)
                    if replaced:
                        hit_chunks[start] = sorted(int(s) for s in replaced)
                n_with_data = len(chunks_to_process) - len(band_empty_chunks)
                log(f"  {corr_name}: {len(hit_chunks)}/{n_with_data} chunk(s) with data")
                for start, scans_hit in sorted(hit_chunks.items()):
                    log(f"    chunk {chunk_numbers[start]} ({chunk_freqs[start]:.4f} GHz): "
                        f"scan(s) {', '.join(str(s) for s in scans_hit)}")
            log(f"Fully flagged chunks found: {len(band_empty_chunks)}")
            log("="*60)

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
        
        # PRELIMINARY CHECK: Is this chunk fully flagged? Counted, over every
        # selected correlation, when the cache was built.
        if start_chan in band_empty_chunks:
            log("")
            log(f"WARNING: Chunk {chunk_num} is fully flagged - skipping all correlations")
            continue
        log(f"Chunk has {cache.unflagged[cache.chunk_index[start_chan]]:,} unflagged "
            f"data points - proceeding")

        # Track flagged (baseline, scan) pairs across all correlations for this
        # chunk. {(bl_idx, scan): {trigger_str, ...}} -- a dict so the flag file
        # /log can still say which criterion(a) (FRACTION / SEVERE) caught each
        # one, even though membership is unioned across every correlation.
        chunk_flagged_baseline_scans = {}

        baseline_names = [x[2] for x in baseline_info]  # Same for all correlations
        ant1_names = [x[3] for x in baseline_info]

        # Amplitude on every correlation, phase on the parallel hands with --phase
        chunk_passes = []
        for corr_name, corr_idx in correlations:
            chunk_passes.append((corr_name, corr_idx, 'amplitude'))
            if args.phase and corr_name in PHASE_CORRELATIONS:
                chunk_passes.append((corr_name, corr_idx, 'phase'))

        for corr_name, corr_idx, quantity in chunk_passes:
            is_phase = quantity == 'phase'
            pass_label = f'{corr_name} {quantity}'
            log("")
            log(f"Analyzing correlation: {corr_name} (index {corr_idx}), {quantity}")

            # This pass's per-baseline values and weights, from the cache
            amp_arrays, weight_arrays = cache.pass_arrays(corr_name, quantity, start_chan)

            # SECONDARY CHECK: if no valid data for this pass, skip
            total_valid = sum(np.isfinite(amp).sum() for amp in amp_arrays)
            if total_valid == 0:
                log("")
                log(f"WARNING: Chunk {chunk_num}, {pass_label} has no valid data - skipping")
                continue
    
            # Flatten data for analysis (used for outlier detection and optional plotting)
            log("Preparing data for analysis...")
            
            # Build a single DataFrame. 'value' holds whichever quantity this
            # pass tests -- amplitude, or phase in radians.
            bl_idx_flat, values_flat, weights_flat, scans_flat, windows_flat = flatten_valid(
                amp_arrays, weight_arrays, scan_arrays, window_arrays)
            df = pd.DataFrame({
                'baseline_idx': bl_idx_flat,
                'value': values_flat,
                'weight': weights_flat,
                'ant1': np.asarray(ant1_names, dtype=object)[bl_idx_flat],
                'scan': scans_flat,
                'window': windows_flat,
            })
            log(f"Total data points: {len(df):,}")
    
            # Convert ant1 to categorical for datashader compatibility
            df['ant1'] = pd.Categorical(df['ant1'])
    
            # =========================================================
            # Amplitude statistics: per-scan centre, ensemble spread
            #
            # Centre is taken per scan. A source's amplitude is not
            # constant across a track -- on a polarised calibrator the
            # cross-hand amplitude |XY| = |p.I.cos(2X + phi)| sweeps
            # through zero with parallactic angle X -- so a single median
            # pooled over all scans sits mid-sweep and represents none of
            # them. Comparing against it flags whole scans for behaving
            # exactly as they should.
            #
            # Spread is taken across scans, not within one. A scan allowed
            # to set its own width normalises away its own broadening and
            # escapes detection; the median of the per-scan MADs holds the
            # width at the array norm and tolerates a minority of bad
            # scans. Broadened baselines are then measured against how the
            # array normally behaves rather than against their own noise.
            # =========================================================
            amp_values = df['value'].values
            weight_values = df['weight'].values
            scan_values = df['scan'].values

            raw_scan_centers, raw_scan_mads = compute_scan_stats(
                amp_values, weight_values, scan_values, circular=is_phase)

            # A scan whose statistics in this chunk are inflated relative to
            # the band fit (RFI filling most of the chunk raises both centre
            # and MAD, lifting the threshold with them) is tested against the
            # fitted values instead.
            scan_centers = dict(raw_scan_centers)
            scan_mads = dict(raw_scan_mads)
            band_replacements = {}
            if band_ref is not None and not is_phase:
                chunk_ref = {scan_id: per_chunk[start_chan]
                             for scan_id, per_chunk in band_ref.get(corr_name, {}).items()
                             if start_chan in per_chunk}
                band_replacements = apply_band_reference(
                    scan_centers, scan_mads, chunk_ref,
                    args.band_mad_cap, args.band_centre_sigma)
                for scan_id, (fit_centre, fit_mad, _) in band_replacements.items():
                    scan_centers[scan_id] = fit_centre
                    scan_mads[scan_id] = fit_mad

            finite_mads = [m for m in scan_mads.values() if np.isfinite(m) and m > 0]
            ensemble_mad = float(np.median(finite_mads)) if finite_mads else np.nan

            if not np.isfinite(ensemble_mad) or ensemble_mad <= 0:
                log("")
                log(f"WARNING: Chunk {chunk_num}, {pass_label} has a zero/non-finite "
                    f"ensemble MAD - no usable spread, skipping")
                continue

            # Phase works in radians throughout; the log carries degrees alongside
            # the headline numbers, which is how they are easier to judge.
            def in_degrees(value):
                return f' ({np.degrees(value):.2f} deg)' if is_phase else ''

            # Representative values, used for plotting and for the log summary
            if is_phase:
                median_amp = float(wrap_angle(circular_weighted_mean(
                    np.array(list(scan_centers.values())),
                    np.ones(len(scan_centers)))))
            else:
                median_amp = float(np.median(list(scan_centers.values())))
            mad_amp = ensemble_mad

            log("")
            log(f"{quantity.capitalize()} statistics:")
            log(f"  Scans in chunk: {len(scan_centers)}")
            log(f"  Per-scan weighted medians: "
                f"{', '.join(f'{int(s)}:{c:.4f}' for s, c in sorted(raw_scan_centers.items()))}")
            log(f"  Per-scan weighted MADs: "
                f"{', '.join(f'{int(s)}:{m:.4f}' for s, m in sorted(raw_scan_mads.items()))}")
            if band_replacements:
                log(f"  Scan statistics replaced with band fit (inflated relative to the band):")
                for scan_id, (fit_centre, fit_mad, reason) in sorted(band_replacements.items()):
                    log(f"    Scan {int(scan_id)}: centre {fit_centre:.4f}, MAD {fit_mad:.4f} "
                        f"({reason})")
            log(f"  Median of per-scan centres: {median_amp:.4f}{in_degrees(median_amp)}")
            log(f"  Ensemble MAD (median of per-scan weighted MADs): "
                f"{ensemble_mad:.4f}{in_degrees(ensemble_mad)}")

            # Outlier detection per baseline
            point_override_enabled = (args.point_outlier_min_count > 0
                                       and args.point_sigma_threshold > 0)
            point_override_str = (f"{args.point_sigma_threshold} sigma" if point_override_enabled
                                   else "DISABLED")
            log("")
            log(f"Outlier detection (threshold: {args.sigma_threshold} sigma, "
                f"point override: {point_override_str}):")
            threshold_mad = args.sigma_threshold * 1.4826 * ensemble_mad
            # Only meaningful when the override is enabled; left at +inf otherwise
            # so is_severe_outlier below is trivially all-False without a branch.
            point_threshold_mad = (args.point_sigma_threshold * 1.4826 * ensemble_mad
                                    if point_override_enabled else float('inf'))
            log(f"  MAD threshold: {threshold_mad:.4f}{in_degrees(threshold_mad)}")
            if point_override_enabled:
                log(f"  Point-override MAD threshold: {point_threshold_mad:.4f}"
                    f"{in_degrees(point_threshold_mad)}")

            # Determine effective outlier mode based on correlation type
            if is_phase:
                # A phase excursion is bad in either direction, so the
                # directional modes have nothing to say about it.
                effective_mode = 'both'
                log(f"  Mode: BOTH (enforced for phase)")
            elif args.outlier_mode == 'mixed':
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
    
            sensitivity = broadening_sensitivity(
                args.sigma_threshold, args.outlier_fraction_threshold,
                two_tailed=(effective_mode == 'both'))
            if np.isfinite(sensitivity):
                log(f"  Sensitivity: flags a baseline once its scatter reaches "
                    f"{sensitivity:.2f}x the array norm "
                    f"(sigma={args.sigma_threshold}, "
                    f"outlier fraction={args.outlier_fraction_threshold}, "
                    f"{'two' if effective_mode == 'both' else 'one'}-tailed)")

            # A scan whose own weighted median isn't representative (a
            # scan-wide shift) has its centre replaced by a reference built
            # from its nearest scans in time, before any baseline in it is
            # tested. Otherwise the shift moves every baseline's residual
            # together and nothing in that scan looks anomalous against its
            # own equally-shifted median.
            scan_replacements = local_scan_reference(
                scan_centers, scan_mads, ensemble_mad,
                args.scan_window_count, args.scan_window_sigma, effective_mode,
                circular=is_phase)
            corrected_scan_centers = dict(scan_centers)
            corrected_scan_centers.update(scan_replacements)
            if scan_replacements:
                log(f"  Scan centre(s) replaced with local reference "
                    f"(own median not representative):")
                for scan_id, ref in sorted(scan_replacements.items()):
                    raw = scan_centers[scan_id]
                    shift = wrap_angle(ref - raw) if is_phase else ref - raw
                    log(f"    Scan {int(scan_id)}: raw median {raw:.4f} -> "
                        f"local reference {ref:.4f} ({shift:+.4f})")

            centres = df['scan'].map(corrected_scan_centers).values
            offsets = wrap_angle(amp_values - centres) if is_phase else amp_values - centres
            residuals = np.abs(offsets)

            # Mark outliers against each point's own scan centre. is_severe_outlier
            # uses the same directional logic as is_outlier (cross-hand tests high
            # excursions only) but the stricter --point-sigma-threshold, so a
            # single severe point can flag a baseline that a fraction test alone
            # would let through.
            if effective_mode == 'high':
                df['is_outlier'] = offsets > threshold_mad
                df['is_severe_outlier'] = offsets > point_threshold_mad
            elif effective_mode == 'low':
                df['is_outlier'] = -offsets > threshold_mad
                df['is_severe_outlier'] = -offsets > point_threshold_mad
            else:  # 'both'
                df['is_outlier'] = residuals > threshold_mad
                df['is_severe_outlier'] = residuals > point_threshold_mad

            # Per-baseline outlier statistics using groupby
            baseline_stats = df.groupby('baseline_idx').agg({
                'is_outlier': ['sum', 'count'],
                'is_severe_outlier': ['sum'],
            })
            baseline_stats.columns = ['outlier_count', 'total_count', 'severe_outlier_count']
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
        
            # Outlier statistics per (baseline, scan): a scan sweeping through
            # a polarised calibrator's cross-hand null, or simply carrying its
            # own noise level, shouldn't have its fate decided by every other
            # scan on the same baseline. On a target each scan is further split
            # into its time windows, so a bad few minutes only flags those
            # minutes.
            bl_scan_stats = df.groupby(['baseline_idx', 'scan', 'window']).agg({
                'is_outlier': ['sum', 'count'],
                'is_severe_outlier': ['sum'],
            })
            bl_scan_stats.columns = ['outlier_count', 'total_count', 'severe_outlier_count']
            bl_scan_stats['outlier_fraction'] = bl_scan_stats['outlier_count'] / bl_scan_stats['total_count']

            # A (baseline, scan) is flagged outright if EITHER trigger fires:
            # its outlier fraction within that scan crosses
            # --outlier-fraction-threshold, or it has >= --point-outlier-min-count
            # points past the stricter --point-sigma-threshold.
            fraction_hit = bl_scan_stats['outlier_fraction'] > args.outlier_fraction_threshold
            if point_override_enabled:
                severe_hit = bl_scan_stats['severe_outlier_count'] >= args.point_outlier_min_count
            else:
                severe_hit = pd.Series(False, index=bl_scan_stats.index)

            bad_baseline_scans = bl_scan_stats[fraction_hit | severe_hit].copy()
            # The phase pass tags its triggers so the flag file and logs say
            # which quantity caught each baseline-scan.
            trigger_prefix = 'PHASE_' if is_phase else ''
            bad_baseline_scans['trigger'] = np.select(
                [fraction_hit.loc[bad_baseline_scans.index] & severe_hit.loc[bad_baseline_scans.index],
                 fraction_hit.loc[bad_baseline_scans.index]],
                [f'{trigger_prefix}FRACTION+SEVERE', f'{trigger_prefix}FRACTION'],
                default=f'{trigger_prefix}SEVERE'
            )

            # Track flagged (baseline, scan) pairs for this correlation,
            # remembering why so the flag file and logs can say which
            # criterion caught each one.
            for (bl_idx, scan, window), row in bad_baseline_scans.iterrows():
                chunk_flagged_baseline_scans.setdefault((bl_idx, int(scan), int(window)),
                                                        set()).add(row['trigger'])

            if len(bad_baseline_scans) > 0:
                n_severe_only = (bad_baseline_scans['trigger'] == f'{trigger_prefix}SEVERE').sum()
                n_distinct_bl = bad_baseline_scans.index.get_level_values('baseline_idx').nunique()
                log("")
                log(f"Flagging {len(bad_baseline_scans)} baseline-scan pair(s) "
                    f"({n_distinct_bl} distinct baseline(s)) for {pass_label} "
                    f"({n_severe_only} by point-level override alone)")

                # Just track for now, will write flags after all correlations processed
            else:
                log("")
                log(f"No baselines flagged for {pass_label}")
            
            # Create plot for this chunk/correlation if args.num_chunks < 3
            if args.num_chunks is not None and args.num_chunks < 3:
                log("")
                log(f"Generating plots for chunk {chunk_num}, {pass_label}...")
                create_plot(df, baseline_names, median_amp, mad_amp, threshold_mad,
                            bad_baseline_scans, start_chan, end_chan, chunk_num, corr_name,
                            effective_mode, args.save_dir,
                            field_name, run_tag, quantity=quantity)
                create_histogram_plot(df, median_amp, mad_amp, threshold_mad, chunk_num, corr_name,
                                     effective_mode, args.save_dir, run_tag, quantity=quantity)
        
        # End of correlation loop - now process antenna-level flagging
        
        # Check if we stopped due to memory error
        if memory_error_stop:
            log("Exiting chunk processing loop due to memory constraints")
            break
        
        log("")
        log(f"Analyzing antenna-level flagging for chunk {chunk_num}...")

        # Antenna flagging is decided scan by scan (window by window on a
        # target): an antenna is flagged when more than --antenna-flag-threshold
        # of its baselines failed there, and the flag covers only that scan
        # (or window). The protected antenna is never flagged.
        flagged_bls_by_window = {}  # {(scan, window): {bl_idx, ...}}
        for bl_idx, scan, window in chunk_flagged_baseline_scans:
            flagged_bls_by_window.setdefault((scan, window), set()).add(bl_idx)

        flagged_antennas = []  # [(ant_name, scan, window, fraction, flagged, total)]
        for scan, window in all_windows:
            window_flagged_bls = flagged_bls_by_window.get((scan, window), set())

            # Count how many baselines each antenna participates in (total and
            # flagged in this scan or window)
            antenna_stats = {}  # {ant_name: {'total': N, 'flagged': M}}
            for bl_idx, (a1, a2, bl_name, ant1_name, ds) in enumerate(baseline_info):
                for ant in (antenna_names[a1], antenna_names[a2]):
                    stats = antenna_stats.setdefault(ant, {'total': 0, 'flagged': 0})
                    stats['total'] += 1
                    if bl_idx in window_flagged_bls:
                        stats['flagged'] += 1

            for ant_name, stats in antenna_stats.items():
                if stats['total'] > 0:
                    fraction = stats['flagged'] / stats['total']
                    if fraction > args.antenna_flag_threshold:
                        if ant_name == protected_antenna:
                            protected_hits.append(f"chunk {chunk_num}, scan {window_label(scan, window)}: "
                                                  f"{fraction:.0%} of its baselines bad")
                        else:
                            flagged_antennas.append((ant_name, scan, window, fraction,
                                                     stats['flagged'], stats['total']))

                # Global statistics are kept per antenna per scan (per window
                # on a target), in baseline counts summed over chunks
                ant_global = global_antenna_stats.setdefault(ant_name, {'scans': {}, 'chunks': []})
                scan_global = ant_global['scans'].setdefault((scan, window), {'total': 0, 'flagged': 0})
                scan_global['total'] += stats['total']
                scan_global['flagged'] += stats['flagged']
                if stats['flagged'] > 0 and chunk_num not in ant_global['chunks']:
                    ant_global['chunks'].append(chunk_num)

        if len(flagged_antennas) > 0:
            log("")
            log(f"Antenna-scans exceeding {args.antenna_flag_threshold*100:.0f}% threshold:")
            for ant_name, scan, window, frac, flagged, total in sorted(
                    flagged_antennas, key=lambda x: (x[1], x[2], -x[3])):
                log(f"  {ant_name} scan {window_label(scan, window)}: {frac:.2%} "
                    f"({flagged}/{total} baselines flagged) - FLAGGING ANTENNA HERE")
        else:
            log("No antenna-scans exceed the antenna flagging threshold")

        # Write flags to file (baseline-scans + antennas)
        if len(chunk_flagged_baseline_scans) > 0 or len(flagged_antennas) > 0:
            log("")
            log(f"Writing flags to {flag_file_path}...")

            # Prepare baseline-scan flags (sorted by ANTENNA1, ANTENNA2, scan, window)
            bad_baseline_scan_data = []
            for (bl_idx, scan, window), triggers in chunk_flagged_baseline_scans.items():
                a1, a2, bl_name, _, _ = baseline_info[bl_idx]
                trigger_str = '+'.join(sorted(triggers))
                bad_baseline_scan_data.append((a1, a2, bl_idx, bl_name, scan, window, trigger_str))
            bad_baseline_scan_data.sort(key=lambda x: (x[0], x[1], x[4], x[5]))

            # Each line flags one scan (or window) on one baseline, plus the
            # scans it propagates to, across all correlations
            with open(flag_file_path, 'a') as f:
                for base_num, (a1, a2, bl_idx, bl_name, scan, window, trigger_str) in enumerate(bad_baseline_scan_data):
                    reason = (f"BASED_FLAGGER_CHUNK_{chunk_num}_BASE{base_num}_"
                              f"SCAN{window_label(scan, window)}_{trigger_str}")
                    flag_line = (f"mode='manual' antenna='{bl_name}' spw='0:{start_chan}~{end_chan-1}' "
                                 f"{flag_selection(scan, window)} reason='{reason}'\n")
                    f.write(flag_line)

            # Antenna flags cover that antenna's baselines in the one scan (or
            # window) it failed, plus the scans that scan propagates to
            with open(flag_file_path, 'a') as f:
                for ant_idx, (ant_name, scan, window, frac, flagged, total) in enumerate(flagged_antennas):
                    reason = f"BASED_FLAGGER_CHUNK_{chunk_num}_ANT{ant_idx}_SCAN{window_label(scan, window)}"
                    flag_line = (f"mode='manual' antenna='{ant_name}' spw='0:{start_chan}~{end_chan-1}' "
                                 f"{flag_selection(scan, window)} reason='{reason}'\n")
                    f.write(flag_line)

            # Calculate new flags for this chunk: the rows of every
            # (baseline, scan, window) covered by a baseline flag or by an
            # antenna flag, each counted once.
            n_freq_in_chunk = end_chan - start_chan
            n_corr_flagged = len(correlations)

            flagged_ants_by_window = {}  # {(scan, window): {ant_name, ...}}
            for ant_name, scan, window, _, _, _ in flagged_antennas:
                flagged_ants_by_window.setdefault((scan, window), set()).add(ant_name)

            covered_baseline_scans = set(chunk_flagged_baseline_scans)
            for (scan, window), ant_names in flagged_ants_by_window.items():
                for bl_idx, (a1, a2, _, _, _) in enumerate(baseline_info):
                    if antenna_names[a1] in ant_names or antenna_names[a2] in ant_names:
                        covered_baseline_scans.add((bl_idx, scan, window))

            flagged_rows = sum(
                int(np.sum((scan_arrays[bl_idx] == scan) & (window_arrays[bl_idx] == window)))
                for (bl_idx, scan, window) in covered_baseline_scans
            )
            for bl_idx, scan, window in covered_baseline_scans:
                if bl_idx in protected_bls:
                    mark_protected([start_chan], bl_idx, scan, window)

            chunk_new_flags = flagged_rows * n_freq_in_chunk * n_corr_flagged

            cumulative_new_flags += chunk_new_flags
            n_distinct_bl_flagged = len({key[0] for key in covered_baseline_scans})
            log(f"New flags this chunk: {chunk_new_flags:,} "
                f"({len(chunk_flagged_baseline_scans)} baseline-scan flag(s), "
                f"{n_distinct_bl_flagged} distinct baseline(s))")
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
    
    # Antenna statistics summed over all chunks, scan by scan (window by window
    # on a target): an antenna is flagged across the band in a scan or window
    # it fails, and in that one only. The cap applies to each scan or window.
    global_scans = sorted({key for stats in global_antenna_stats.values() for key in stats['scans']})

    capped_flagged = []  # [(ant_name, scan, window, frac, flagged, total, methods)]

    if len(global_scans) == 0:
        log("No antenna statistics available")

    for scan, window in global_scans:
        antenna_fractions = []
        for ant_name, stats in global_antenna_stats.items():
            scan_stats = stats['scans'].get((scan, window))
            if scan_stats and scan_stats['total'] > 0:
                fraction = scan_stats['flagged'] / scan_stats['total']
                antenna_fractions.append((ant_name, fraction, scan_stats['flagged'], scan_stats['total']))
        if len(antenna_fractions) == 0:
            continue

        log("")
        log(f"--- Scan {window_label(scan, window)}: analyzing {len(antenna_fractions)} antennas ---")
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

        # The protected antenna is never flagged, and doesn't use up a place under the cap
        if protected_antenna in flagged_dict:
            info = flagged_dict.pop(protected_antenna)
            log(f"Not flagging protected antenna {protected_antenna} ({info['frac']:.2%} of its baselines bad)")
            protected_hits.append(f"global pass, scan {window_label(scan, window)}: "
                                  f"{info['frac']:.0%} of its baselines bad")

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
            scan_capped = combined_flagged[:args.antenna_flag_cap]
            excluded = combined_flagged[args.antenna_flag_cap:]
            log("")
            log(f"Antennas being flagged:")
            for ant_name, frac, flagged, total, methods in scan_capped:
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
            scan_capped = combined_flagged
            if len(scan_capped) > 0:
                log(f"Flagging {len(scan_capped)} antenna(s):")
                for ant_name, frac, flagged, total, methods in scan_capped:
                    methods_str = '+'.join(methods)
                    log(f"  {ant_name}: {frac:.2%} ({flagged}/{total} baselines) [{methods_str}]")
            else:
                log("No antennas meet flagging criteria")
        
        # Write this scan's (or window's) global antenna flags to file, plus
        # the scans it propagates to
        if len(scan_capped) > 0 and len(chunks_to_process) > 0:
            min_chan = min(chunks_to_process)
            max_chan = max(chunks_to_process) + args.channel_group_size
            max_chan = min(max_chan, n_chan)  # Don't exceed total channels

            log("")
            log(f"Writing {len(scan_capped)} global antenna flag(s) for scan "
                f"{window_label(scan, window)}, channels {min_chan}~{max_chan-1}...")

            with open(flag_file_path, 'a') as f:
                for ant_idx, (ant_name, frac, flagged, total, methods) in enumerate(scan_capped):
                    reason = f"BASED_FLAGGER_GLOBAL_ANT{ant_idx}_SCAN{window_label(scan, window)}"
                    flag_line = (f"mode='manual' antenna='{ant_name}' spw='0:{min_chan}~{max_chan-1}' "
                                 f"{flag_selection(scan, window)} reason='{reason}'\n")
                    f.write(flag_line)

            log(f"✓ Wrote {len(scan_capped)} global antenna flags to {flag_file_path}")

        capped_flagged += [(ant_name, scan, window, frac, flagged, total, methods)
                           for ant_name, frac, flagged, total, methods in scan_capped]
    
    # Print top 5 antennas with highest flagging fraction across all chunks
    log("")
    log("="*60)
    log("TOP 5 ANTENNAS WITH HIGHEST FLAGGING FRACTION (ALL CHUNKS)")
    log("="*60)
    
    # Calculate fractions and sort
    antenna_fractions = []
    for ant_name, stats in global_antenna_stats.items():
        total = sum(s['total'] for s in stats['scans'].values())
        flagged = sum(s['flagged'] for s in stats['scans'].values())
        if total > 0:
            antenna_fractions.append((ant_name, flagged / total, flagged, total, stats['chunks']))
    
    antenna_fractions.sort(key=lambda x: x[1], reverse=True)
    
    # Print top 5
    top_5 = antenna_fractions[:5]
    if len(top_5) > 0:
        for rank, (ant_name, frac, flagged, total, chunks) in enumerate(top_5, start=1):
            chunks_str = ','.join(map(str, sorted(set(chunks))))
            log(f"{rank}. {ant_name}: {frac:.2%} ({flagged}/{total} baseline-scans) - chunks: {chunks_str}")
    else:
        log("No antenna flagging data available")
    log("")

    # Global antenna flags span every chunk processed, so they are counted here
    # rather than in the per-chunk tally, which never saw them.
    if capped_flagged and chunks_to_process:
        span_chan = min(max(chunks_to_process) + args.channel_group_size, n_chan) \
                    - min(chunks_to_process)
        for ant_name, scan, window, _, _, _, _ in capped_flagged:
            ant_bls = [bl_idx for bl_idx, (a1, a2, _, _, _) in enumerate(baseline_info)
                       if ant_name in (antenna_names[a1], antenna_names[a2])]
            ant_scan_rows = sum(int(np.sum((scan_arrays[bl_idx] == scan)
                                           & (window_arrays[bl_idx] == window)))
                                for bl_idx in ant_bls)
            cumulative_new_flags += ant_scan_rows * span_chan * len(correlations)
            # A global antenna flag spans every processed chunk
            for bl_idx in protected_bls.intersection(ant_bls):
                mark_protected(chunks_to_process, bl_idx, scan, window)

    # The protected antenna met the antenna-flagging criteria but was kept
    if protected_hits:
        log("")
        log(f"WARNING: protected antenna {protected_antenna} would have been flagged "
            f"{len(protected_hits)} time(s) but was kept (its baselines can still be flagged):")
        for hit in protected_hits:
            log(f"  {hit}")

    # =================================================================
    # FINAL SUMMARY
    # =================================================================
    new_fraction = cumulative_new_flags / total_data if total_data else 0.0
    total_with_new = old_flags + cumulative_new_flags
    total_fraction = total_with_new / total_data if total_data else 0.0

    log("")
    log("="*70)
    log(f"SUMMARY [{field_name}]: flagged {old_flag_fraction:.4%} -> {total_fraction:.4%} "
        f"(+{new_fraction:.4%}, {cumulative_new_flags:,} new of {total_data:,}) | "
        f"log: {log_filename} | flags: {flag_file_path}")
    log("="*70)
    log("")
    log("NOTE: the new-flag count is this field's ESTIMATE, not a measurement -- "
        "basedflagger only writes the flag commands. It covers this field only: flags "
        "passed on to the --flag-fields fields are not counted. The true before/after is reported by the flagdata "
        "mode='summary' pair around the apply in 1GC_05_casa_refcal.py.")
    log("")

    # Protected antenna: how much of its baselines' data was flagged before this
    # run, and would be after this run's flags (estimated, like the summary)
    if protected_antenna is None:
        log("PROTECTED ANTENNA: none")
    else:
        total = before = after = 0
        for start_chan in chunks_to_process:
            k, w = cache.chunk_index[start_chan], cache.chunk_width[start_chan]
            for bl_idx in protected_bls:
                rows = cache.rows_per_bl[bl_idx]
                new = protected_new_rows.get((start_chan, bl_idx), np.zeros(rows, dtype=bool))
                for corr_name, _ in correlations:
                    flagged = ~np.isfinite(cache.memmaps[(corr_name, 'amplitude')][k, bl_idx, :rows, :w])
                    total += flagged.size
                    before += int(flagged.sum())
                    after += int((flagged | new[:, None]).sum())
        before_frac = before / total if total else 0.0
        after_frac = after / total if total else 0.0
        log(f"PROTECTED ANTENNA {protected_antenna} ({len(protected_bls)} baselines): flagged "
            f"{before_frac:.4%} before -> {after_frac:.4%} after "
            f"(+{after_frac - before_frac:.4%}; never flagged as a whole, "
            f"{len(protected_hits)} time(s) it would have been)")
    log("")


if __name__ == "__main__":
    main()
