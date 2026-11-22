#!/usr/bin/env python3
"""
Per-antenna full-Stokes diagnostics of a calibrator field.

Reads a calibrated column (CORRECTED_DATA by default) of one field with
dask-ms, averages it in time over each scan and in frequency into a few
channel bins, and caches the result to disk. From the cache it forms Stokes
I, Q, U, V per baseline, per antenna (every baseline to that antenna) and for
the whole array, and asks whether a single antenna drives the time evolution.

For every antenna a these series are compared over the scans

    A(t)      the array, all baselines
    S_a(t)    the baselines to antenna a
    A_-a(t)   the array without antenna a

and three numbers are reported. All statistics are measured against the
array, using medians and MADs so that outlying antennas do not skew them.
The formal error bars (scatter across baselines) are only drawn, never used.

    R(t)      per-scan median of S_a(t) over antennas
    D_a(t)    S_a(t) - R(t), antenna a's deviation from the array

    offset_z  STATIC offset. O_a = median over scans of D_a(t), turned into
              a robust z-score across antennas,
                  offset_z = (O_a - median_a O) / (1.4826 MAD_a O)
              |offset_z| > OFFSET_SIGMA marks an antenna wrong by a constant
    chi2r     TIME EVOLUTION. With the static offset removed,
                  E_a(t) = D_a(t) - O_a
                  sigma(t) = 1.4826 MAD over antennas of E_a(t)   (per scan)
                  z_a(t) = E_a(t) / sigma(t)
                  chi2r = sum_t z_a(t)^2 / (n_scan - 1)
              so a typical antenna sits near 1 by construction, whatever the
              noise or systematics level, and an antenna that evolves
              differently from the rest is large. chi2r > CHI2_THRESHOLD is
              highlighted. z_a(t) is what the heatmap shows
    loo_drop  1 - rms_t[A_-a] / rms_t[A]; positive means removing antenna a
              flattens the array's time evolution

The per-scan MAD needs a reasonable number of antennas (>~10) to be stable.

Plots (in SAVE_DIR) and their colours

    *_stokes_vs_<xaxis>.png   IQUV against time, parallactic angle or scan
        black squares, thick  the array, all baselines
        blue, orange, green   top HIGHLIGHT antennas by chi2r, in that order,
                              among those with chi2r > CHI2_THRESHOLD (ranked separately per Stokes; the
                              legend names them). Fewer, or none, if no
                              antenna stands out
        red, thin             static-offset outliers, |offset_z| > OFFSET_SIGMA
        light grey, thin      every other antenna
    *_grid_<S>_vs_<xaxis>.png one panel per antenna for Stokes S
        grey line             the reference R(t), per-scan antenna median
        points                that antenna; blue/orange/green for the top
                              chi2r antennas, red for offset outliers, dark
                              grey otherwise. Highlighted panels have a bold
                              title, which carries chi2r and offset_z
    *_ranking.png             bars of chi2r (log, dotted at 1, dashed at
                              CHI2_THRESHOLD), offset_z
                              (dashed lines at +-OFFSET_SIGMA) and loo_drop,
                              per antenna, same colours as above
    *_antenna_scan_heatmap.png  z_a(t) per antenna and scan on a red-white-
                              blue scale clipped at +-5; red above, blue
                              below that antenna's usual offset, in units of
                              that scan's robust spread across antennas.
                              Each row's static offset is removed, so this
                              shows time evolution only (offset_z covers the
                              rest)
    *_ranking.txt             the three numbers per antenna and Stokes

Linear feeds, Stokes as

    I = (XX + YY)/2   Q = (XX - YY)/2   U = (XY + YX)/2   V = (XY - YX)/(2i)

The sign of V depends on the feed and correlator convention (FLIP_V).

The read is the expensive step and is done once. A second run on the same MS,
field, column, weights and channel binning loads the cache and only replots,
so the plotting settings (COMPONENT, FRACTIONAL, XAXIS, PLOT_BINS, PER_BIN)
can be changed freely. The cache is invalidated automatically when
the MS is modified.

All inputs are globals in the CONFIGURATION block below. Edit them and run
    python akhtools_check_antenna_time_diagnostic.py
"""

import warnings
import os

warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'

from types import SimpleNamespace
import hashlib
import json
import logging
import re
from datetime import datetime, timezone

import numpy as np
import psutil
import dask
from daskms import xds_from_ms, xds_from_table

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

STOKES = ('I', 'Q', 'U', 'V')
CORR_TYPE_NAMES = {9: 'XX', 10: 'XY', 11: 'YX', 12: 'YY'}
CACHE_VERSION = 1

# Fixed categorical order for highlighted antennas (never cycled)
HIGHLIGHT_COLOURS = ['#2a6fdb', '#e0702a', '#1e9e6a', '#5b4bc4', '#8a6d1e', '#0f8b9c']
OFFSET_COLOUR = '#c62828'      # static-offset outliers, never used for anything else
ARRAY_COLOUR = '#1a1a1a'
ORDINARY_COLOUR = '#555555'    # unhighlighted antennas in the grid
BACKGROUND_COLOUR = '#b8b8b8'


# =====================================================================
# CONFIGURATION (edit here)
# =====================================================================

# Input
MS_FILE = '/path/to/calibrator.ms'
FIELD = 'J1331+3030'              # field ID (int) or name

# What is read (these define the cache)
DATA_COLUMN = 'CORRECTED_DATA'
WEIGHT_COLUMN = 'WEIGHT'          # 'WEIGHT', 'WEIGHT_SPECTRUM' or 'none'
CHAN_BINS = 8                     # equal-width channel bins across CHAN_RANGE
CHAN_RANGE = None                 # (lo, hi) with hi exclusive, or None for all

# Cache
CACHE_DIR = 'CACHE/'
OVERWRITE_CACHE = False           # re-read the MS even if a matching cache exists
MEMORY_FRACTION = 0.3             # fraction of available memory one read block may use

# What is plotted (free to change between runs, the cache is reused)
COMPONENT = 'real'                # 'real', 'imag' or 'amp' of the vector average
FRACTIONAL = False                # plot Q/I, U/I, V/I instead of Q, U, V
FLIP_V = False                    # multiply Stokes V by -1
XAXIS = 'time'                    # 'time', 'parang' or 'scan'
PLOT_BINS = 'all'                 # list of bin indices combined for the main plots, or 'all'
PER_BIN = False                   # also make the full plot set for every bin
MIN_BASELINES = 3                 # minimum unflagged baselines for an antenna/scan point
HIGHLIGHT = 3                     # top-ranked antennas by chi2r highlighted per Stokes
CHI2_THRESHOLD = 6.0              # ... but only if their robust chi2r exceeds this (typical ~1)
OFFSET_SIGMA = 5.0                # |offset_z| above which an antenna is a static-offset outlier
MAKE_GRIDS = True                 # per-antenna small-multiple grids

# Output
SAVE_DIR = 'PLOTS/'
LOG_DIR = 'LOGS/'
OUTPUT_SUFFIX = ''
PRINTING = True                   # False logs to file only


def get_config():
    """The configuration globals, under the names the rest of the script uses."""
    return SimpleNamespace(
        ms_file=MS_FILE, field=str(FIELD),
        data_column=DATA_COLUMN, weight_column=WEIGHT_COLUMN,
        chan_bins=CHAN_BINS, chan_range=CHAN_RANGE,
        cache_dir=CACHE_DIR, overwrite_cache=OVERWRITE_CACHE,
        memory_fraction=MEMORY_FRACTION,
        component=COMPONENT, fractional=FRACTIONAL, flip_v=FLIP_V, xaxis=XAXIS,
        plot_bins=PLOT_BINS, per_bin=PER_BIN, min_baselines=MIN_BASELINES,
        highlight=HIGHLIGHT, chi2_threshold=CHI2_THRESHOLD, offset_sigma=OFFSET_SIGMA, no_grid=not MAKE_GRIDS,
        save_dir=SAVE_DIR, log_dir=LOG_DIR, output_suffix=OUTPUT_SUFFIX,
        no_printing=not PRINTING)


# =====================================================================
# RESOURCES AND LOGGING (as in basedflagger)
# =====================================================================

def _parse_slurm_mem(value):
    """SLURM memory string -> bytes. Bare numbers are MB; M/G/T suffixes honoured."""
    suffix = value[-1].lower()
    multiplier = {'g': 1024 ** 3, 't': 1024 ** 4}.get(suffix, 1024 ** 2)
    return int(value.rstrip('MmGgTt')) * multiplier


def get_slurm_memory_bytes():
    mem_per_node = os.environ.get('SLURM_MEM_PER_NODE')
    mem_per_cpu = os.environ.get('SLURM_MEM_PER_CPU')
    try:
        if mem_per_node:
            return _parse_slurm_mem(mem_per_node), f'SLURM_MEM_PER_NODE={mem_per_node}'
        if mem_per_cpu:
            cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
            ntasks = int(os.environ.get('SLURM_NTASKS', 1))
            n = cpus * ntasks
            return _parse_slurm_mem(mem_per_cpu) * n, f'SLURM_MEM_PER_CPU={mem_per_cpu} x {n} CPUs'
    except ValueError:
        pass
    return None, None


def get_allocated_cpus():
    cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    ntasks = os.environ.get('SLURM_NTASKS')
    job_cpus = re.match(r'(\d+)', os.environ.get('SLURM_JOB_CPUS_PER_NODE', ''))
    try:
        if cpus_per_task and ntasks:
            return int(cpus_per_task) * int(ntasks), 'SLURM_CPUS_PER_TASK x SLURM_NTASKS'
        if cpus_per_task:
            return int(cpus_per_task), 'SLURM_CPUS_PER_TASK'
        if job_cpus:
            return int(job_cpus.group(1)), 'SLURM_JOB_CPUS_PER_NODE'
        if ntasks:
            return int(ntasks), 'SLURM_NTASKS'
    except ValueError:
        pass
    return os.cpu_count() or 1, 'os.cpu_count() (no SLURM allocation detected)'


def available_memory_bytes():
    """The more restrictive of the SLURM allocation and the node's free memory."""
    slurm_mem, src = get_slurm_memory_bytes()
    avail = psutil.virtual_memory().available
    if slurm_mem is not None:
        used = psutil.Process().memory_info().rss
        return min(slurm_mem - used, avail), f'min(SLURM, available) [{src}]'
    return avail, 'node available (no SLURM allocation detected)'


ENABLE_PRINTING = True


def log(message=""):
    logging.info(message)
    if ENABLE_PRINTING:
        if message == "":
            print(flush=True)
        else:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def scrub_name(name):
    return re.sub(r'[^A-Za-z0-9_.+-]', '_', str(name))


# =====================================================================
# MS METADATA
# =====================================================================

def _decode(values):
    # String columns come back from casacore as either str or bytes depending
    # on the python-casacore version, so normalise to plain str.
    return np.array([v.decode() if isinstance(v, (bytes, np.bytes_)) else str(v) for v in values])


# ---------------------------------------------------------------------
# Subtable reads
#
# An MS is a directory of casacore tables. The main table holds the
# visibilities, the subtables (ANTENNA, FIELD, SPECTRAL_WINDOW,
# POLARIZATION, ...) hold the metadata. dask-ms addresses a subtable with
# the "<ms>::<SUBTABLE>" syntax.
#
# xds_from_table() returns a LIST of xarray Datasets, one per group. With no
# group_cols the whole table is a single group, hence the [0]. Each column
# becomes a lazy dask array (ds.COLUMN.data); nothing is read from disk until
# dask.compute() (or .compute()) is called, and passing several arrays to one
# dask.compute() reads them in one pass over the table.
# ---------------------------------------------------------------------

def load_antenna_table(ms_file):
    # ANTENNA has one row per antenna, the row index is the antenna ID used by
    # ANTENNA1/ANTENNA2 in the main table. POSITION is ITRF (x, y, z) in metres.
    t = xds_from_table(f"{ms_file}::ANTENNA", columns=["NAME", "POSITION"])[0]
    names, pos = dask.compute(t.NAME.data, t.POSITION.data)
    return _decode(names), np.asarray(pos)


def load_field_table(ms_file):
    # FIELD has one row per field, the row index is FIELD_ID. PHASE_DIR has
    # shape (field, poly, 2), the polynomial term 0 being the (RA, Dec) in rad.
    t = xds_from_table(f"{ms_file}::FIELD", columns=["NAME", "PHASE_DIR"])[0]
    names, phase_dir = dask.compute(t.NAME.data, t.PHASE_DIR.data)
    return _decode(names), np.asarray(phase_dir).reshape(len(names), -1, 2)[:, 0, :]


def load_chan_freq(ms_file):
    # SPECTRAL_WINDOW has one row per SPW, CHAN_FREQ shape (spw, chan) in Hz.
    # MeerKAT data have a single SPW, so row 0 is taken.
    t = xds_from_table(f"{ms_file}::SPECTRAL_WINDOW", columns=["CHAN_FREQ"])[0]
    freq = np.asarray(t.CHAN_FREQ.data.compute())
    return freq[0] if freq.ndim > 1 else freq


def load_corr_names(ms_file):
    # POLARIZATION.CORR_TYPE gives the Stokes enum of each slot along the
    # main table's corr axis (9=XX, 10=XY, 11=YX, 12=YY for linear feeds).
    # Reading it rather than assuming the order means the Stokes formation
    # below indexes the right slots whatever order the MS was written in.
    t = xds_from_table(f"{ms_file}::POLARIZATION", columns=["CORR_TYPE"])[0]
    corr = np.asarray(t.CORR_TYPE.data.compute())
    corr = corr[0] if corr.ndim > 1 else corr
    names = [CORR_TYPE_NAMES.get(int(c), f'C{int(c)}') for c in corr]
    missing = {'XX', 'XY', 'YX', 'YY'} - set(names)
    if missing:
        raise RuntimeError(f"MS correlations {names} lack {sorted(missing)}; "
                           f"full-Stokes needs XX, XY, YX, YY")
    return names


def resolve_field(field_names, field_arg):
    listing = ', '.join(f'{i}={n}' for i, n in enumerate(field_names))
    try:
        fid = int(field_arg)
    except ValueError:
        for i, name in enumerate(field_names):
            if name.lower().strip() == field_arg.lower().strip():
                return i, name
        raise ValueError(f"Field '{field_arg}' not found. Available fields: {listing}")
    if not 0 <= fid < len(field_names):
        raise ValueError(f"Field ID {fid} out of range. Available fields: {listing}")
    return fid, field_names[fid]


def ms_modification_time(ms_file):
    """Latest mtime of the main table's files, so any column rewrite invalidates the cache."""
    latest = 0.0
    for entry in os.scandir(ms_file):
        # table.lock is rewritten by every open, including this script's own read
        if entry.is_file() and entry.name != 'table.lock':
            latest = max(latest, entry.stat().st_mtime)
    return latest


def parse_chan_range(chan_range, n_chan):
    if chan_range is None:
        return 0, n_chan
    lo, hi = (int(v) for v in chan_range)
    lo, hi = max(lo, 0), min(hi, n_chan)
    if hi <= lo:
        raise ValueError(f"Empty channel range {chan_range} for {n_chan} channels")
    return lo, hi


def scan_parallactic_angles(ant_pos, phase_dir, times):
    """Parallactic angle (deg) at the array centre for each scan's mean time, or None."""
    try:
        from astropy.utils import iers
        iers.conf.auto_download = False
        from astropy.coordinates import EarthLocation
        from astropy.time import Time
        import astropy.units as u

        centre = ant_pos.mean(axis=0)
        loc = EarthLocation.from_geocentric(*centre, unit=u.m)
        ra, dec = phase_dir
        t = Time(times / 86400.0, format='mjd', scale='utc')
        lst = t.sidereal_time('apparent', longitude=loc.lon).rad
        ha = lst - ra
        lat = loc.lat.rad
        q = np.arctan2(np.sin(ha), np.tan(lat) * np.cos(dec) - np.sin(dec) * np.cos(ha))
        return np.degrees(q)
    except Exception as exc:
        log(f"WARNING: parallactic angle not computed ({exc})")
        return None


# =====================================================================
# READ AND CACHE
# =====================================================================

def cache_path(args, field_id, field_name, chan_lo, chan_hi):
    """
    Cache file CACHE_DIR/<obsid>_<field>_<hash>.npz. The obs ID is the leading
    digits of the MS name and the field is by name, so different fields of one
    MS never share a file. The hash covers everything else that changes what
    is read, and the full key is stored inside the file and checked on load.
    """
    key = json.dumps({
        'version': CACHE_VERSION,
        'ms': os.path.abspath(args.ms_file).rstrip('/'),
        'mtime': ms_modification_time(args.ms_file),
        'field_id': field_id,
        'field_name': str(field_name),
        'data_column': args.data_column,
        'weight_column': args.weight_column,
        'chan_bins': args.chan_bins,
        'chan_range': [chan_lo, chan_hi],
    }, sort_keys=True)
    digest = hashlib.sha1(key.encode()).hexdigest()[:8]
    ms_base = os.path.basename(os.path.abspath(args.ms_file).rstrip('/'))
    match = re.match(r'\d+', ms_base)
    obsid = match.group(0) if match else scrub_name(ms_base)[:16]
    return os.path.join(args.cache_dir, f'{obsid}_{scrub_name(field_name)}_{digest}.npz'), key


def build_cache(args, field_id, corr_names, chan_freq, chan_lo, chan_hi):
    """
    Read the field once in contiguous row blocks (the MS's own row order) and
    accumulate, per (scan, baseline, channel bin, correlation), the weighted
    sum of the visibilities, the sum of weights and the count of unflagged
    samples. Autocorrelations are excluded.
    """
    # -----------------------------------------------------------------
    # 1. Open the main table lazily with dask-ms
    # -----------------------------------------------------------------
    # Only the columns listed are exposed. Everything else in the MS (DATA,
    # MODEL_DATA, UVW, ...) is never touched, which matters because the
    # visibility columns are by far the largest thing on disk.
    cols = [args.data_column, "FLAG", "ANTENNA1", "ANTENNA2", "SCAN_NUMBER", "TIME"]
    if args.weight_column != 'none':
        cols.append(args.weight_column)
    try:
        ds_all = xds_from_ms(
            args.ms_file, columns=cols,
            # taql_where is a TaQL selection applied by casacore before
            # dask-ms sees the rows, so the dataset only contains this field's
            # cross-correlations. Autocorrelations carry no Stokes information
            # for a point source and would dominate any average by amplitude.
            taql_where=f"FIELD_ID=={field_id} AND ANTENNA1!=ANTENNA2",
            # group_cols splits the selected rows into one Dataset per unique
            # value of these columns. Grouping only on DATA_DESC_ID (i.e. SPW
            # plus polarisation setup) keeps every baseline and scan in ONE
            # dataset, so rows stay in the MS's on-disk order and each storage
            # tile is read once. Grouping on ANTENNA1/ANTENNA2 instead (as
            # basedflagger's load_grouped_ms does) gives one dataset per
            # baseline, which scatters the reads across the whole file.
            group_cols=["DATA_DESC_ID"])
    except Exception as exc:
        raise RuntimeError(f"Could not open {args.ms_file} with columns {cols}: {exc}")
    if len(ds_all) != 1:
        raise RuntimeError(f"Expected one spectral window for this field, found {len(ds_all)}")

    # ds is an xarray Dataset with dimensions
    #   row   one per (time, baseline) sample of this field
    #   chan  frequency channels
    #   corr  correlation slots, ordered as POLARIZATION.CORR_TYPE says
    # so the data column has shape (row, chan, corr) complex64, FLAG the same
    # in bool, WEIGHT (row, corr) float32, and ANTENNA1/2, SCAN_NUMBER, TIME
    # are (row,).
    #
    # Without xarray installed, dask-ms returns its own minimal Dataset, which
    # has attribute access to variables (ds.FLAG) and .data, but none of the
    # xarray methods (rename, isel, ...). So the underlying dask arrays are
    # pulled out here and everything below slices those directly, which works
    # with or without xarray.
    ds = ds_all[0]
    data_col = getattr(ds, args.data_column).data
    flag_col = ds.FLAG.data
    weight_col = getattr(ds, args.weight_column).data if args.weight_column != 'none' else None

    ant1, ant2, scan, time = dask.compute(ds.ANTENNA1.data, ds.ANTENNA2.data,
                                          ds.SCAN_NUMBER.data, ds.TIME.data)
    n_rows = len(ant1)
    if n_rows == 0:
        raise RuntimeError("No cross-correlation rows for this field")

    # -----------------------------------------------------------------
    # 3. Map each row to an output cell (scan, baseline)
    # -----------------------------------------------------------------
    # np.unique(..., return_inverse=True) gives the sorted unique values and,
    # for every row, the index of its value in that list. That index is used
    # directly as the accumulator position, so no per-baseline loop is needed.
    # Antenna pairs are put in (low, high) order in case the MS stores any
    # baseline both ways round.
    pairs = np.stack([np.minimum(ant1, ant2), np.maximum(ant1, ant2)], axis=1)
    baselines, row_bl = np.unique(pairs, axis=0, return_inverse=True)
    row_bl = row_bl.ravel()
    scans, row_scan = np.unique(scan, return_inverse=True)
    row_scan = row_scan.ravel()
    n_bl, n_scan, n_corr = len(baselines), len(scans), len(corr_names)

    # Channel bins, as offsets within the read range. np.add.reduceat below
    # sums the channel axis between consecutive start indices, so only the
    # starts are needed.
    edges = np.unique(np.linspace(0, chan_hi - chan_lo, args.chan_bins + 1).round().astype(int))
    bin_starts = edges[:-1]
    n_bin = len(bin_starts)
    bin_freq = np.array([chan_freq[chan_lo + a:chan_lo + b].mean()
                         for a, b in zip(edges[:-1], edges[1:])])

    log(f"Reading {n_rows:,} rows: {n_scan} scans x {n_bl} baselines, "
        f"channels {chan_lo}-{chan_hi - 1} into {n_bin} bins")

    # Accumulators, flattened over (scan, baseline) into one "cell" axis so a
    # single np.bincount per (bin, corr) column can add a whole block at once.
    # float64 throughout, since millions of float32 samples are summed.
    n_cell = n_scan * n_bl
    vsum_re = np.zeros((n_cell, n_bin, n_corr))
    vsum_im = np.zeros((n_cell, n_bin, n_corr))
    wsum = np.zeros((n_cell, n_bin, n_corr))
    nsum = np.zeros((n_cell, n_bin, n_corr))
    row_cell = row_scan * n_bl + row_bl

    # -----------------------------------------------------------------
    # 4. Choose the row block size from the memory available
    # -----------------------------------------------------------------
    # Per row in a block, several (chan, corr) arrays are alive at once: the
    # complex64 data, the bool flags, the float32 weights, the float32
    # weighted copy and its complex product, plus float64 reduceat outputs.
    # 64 bytes per (chan, corr) element is a deliberately generous bound.
    n_chan_read = chan_hi - chan_lo
    bytes_per_row = n_chan_read * n_corr * 64
    budget, mem_src = available_memory_bytes()
    budget *= args.memory_fraction
    row_block = int(np.clip(budget / bytes_per_row, 500, 500000))
    n_blocks = (n_rows + row_block - 1) // row_block
    log(f"Memory budget per block {budget / 1024**3:.2f} GB [{mem_src}]; "
        f"{n_blocks} block(s) of up to {row_block:,} rows")

    # -----------------------------------------------------------------
    # 5. Stream the visibilities block by block
    # -----------------------------------------------------------------
    t_start = datetime.now()
    for b, r0 in enumerate(range(0, n_rows, row_block)):
        r1 = min(r0 + row_block, n_rows)

        # Slicing a dask array is still lazy: it narrows the graph to these
        # rows and to CHAN_RANGE without reading anything. Because rows are
        # contiguous in the MS, casacore reads them as a few large sequential
        # chunks.
        tasks = [data_col[r0:r1, chan_lo:chan_hi, :], flag_col[r0:r1, chan_lo:chan_hi, :]]
        if args.weight_column == 'WEIGHT_SPECTRUM':
            tasks.append(weight_col[r0:r1, chan_lo:chan_hi, :])   # (row, chan, corr)
        elif args.weight_column == 'WEIGHT':
            tasks.append(weight_col[r0:r1])                       # (row, corr)

        # This is the only expensive I/O. dask.compute() executes the graph
        # with the threaded scheduler (capped to the SLURM allocation in
        # main) and returns plain numpy arrays.
        out = dask.compute(*tasks)
        data, flag = out[0], out[1]

        # A sample contributes only if unflagged, finite, and with a finite
        # positive weight. WEIGHT (one per row and corr) is broadcast across
        # channels without copying.
        good = ~flag & np.isfinite(data)
        if args.weight_column == 'WEIGHT_SPECTRUM':
            w = out[2]
        elif args.weight_column == 'WEIGHT':
            w = np.broadcast_to(out[2][:, None, :], data.shape)
        else:
            w = np.ones(data.shape, dtype=np.float32)
        wg = np.where(good & np.isfinite(w), w, 0).astype(np.float32)
        good &= wg > 0
        dw = np.where(good, data, 0) * wg
        del data, flag, w

        # Reduce (row, chan, corr) -> (row, bin, corr) with reduceat, then
        # add every row into its (scan, baseline) cell with bincount. The
        # time average over each scan and the frequency average within each
        # bin are therefore both weighted sums, divided by wsum later.
        cells = row_cell[r0:r1]
        for acc, arr in ((vsum_re, dw.real), (vsum_im, dw.imag), (wsum, wg), (nsum, good)):
            binned = np.add.reduceat(arr, bin_starts, axis=1, dtype=np.float64)
            binned = binned.reshape(len(cells), n_bin * n_corr)
            flat = acc.reshape(n_cell, n_bin * n_corr)
            for j in range(n_bin * n_corr):
                flat[:, j] += np.bincount(cells, weights=binned[:, j], minlength=n_cell)
        del dw, wg, good

        elapsed = (datetime.now() - t_start).total_seconds()
        log(f"  Block {b + 1}/{n_blocks} read ({r1:,}/{n_rows:,} rows, {elapsed / 60:.1f} min)")

    # Mean MS TIME (MJD seconds, UTC) of each scan's rows, used for the x-axis
    # and the parallactic angle.
    scan_time = (np.bincount(row_scan, weights=time, minlength=n_scan)
                 / np.bincount(row_scan, minlength=n_scan))
    shape = (n_scan, n_bl, n_bin, n_corr)
    return {
        'vsum': (vsum_re + 1j * vsum_im).reshape(shape),
        'wsum': wsum.reshape(shape),
        'nsum': nsum.reshape(shape),
        'baselines': baselines.astype(np.int32),
        'scans': scans.astype(np.int32),
        'scan_time': scan_time,
        'bin_edges': edges + chan_lo,
        'bin_freq': bin_freq,
        'corr_names': np.array(corr_names),
    }


def save_cache(path, key, cache):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    tmp = path + '.tmp.npz'
    np.savez(tmp, key=np.array(key), **cache)
    os.replace(tmp, path)


def load_cache(path, key):
    with np.load(path, allow_pickle=False) as f:
        if str(f['key']) != key:
            return None
        return {k: f[k] for k in f.files if k != 'key'}


# =====================================================================
# STOKES AND AGGREGATION
# =====================================================================

def corr_means(cache, bins):
    """Weighted time- and frequency-averaged correlations over the chosen bins."""
    vsum = cache['vsum'][:, :, bins, :].sum(axis=2)
    wsum = cache['wsum'][:, :, bins, :].sum(axis=2)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = np.where(wsum > 0, vsum / wsum, np.nan)
    return mean, wsum


def stokes_per_baseline(mean, wsum, corr_names, flip_v):
    """
    Complex Stokes (scan, baseline, 4) and their weights. The weight of a
    half-sum or half-difference of two correlations is 4 / (1/w1 + 1/w2).
    """
    c = {name: k for k, name in enumerate(corr_names)}
    xx, xy, yx, yy = (mean[..., c[n]] for n in ('XX', 'XY', 'YX', 'YY'))
    wxx, wxy, wyx, wyy = (wsum[..., c[n]] for n in ('XX', 'XY', 'YX', 'YY'))
    v_sign = -1.0 if flip_v else 1.0
    S = np.stack([(xx + yy) / 2, (xx - yy) / 2, (xy + yx) / 2,
                  v_sign * (xy - yx) / 2j], axis=-1)
    with np.errstate(divide='ignore', invalid='ignore'):
        w_par = np.where((wxx > 0) & (wyy > 0), 4.0 / (1.0 / wxx + 1.0 / wyy), 0.0)
        w_crs = np.where((wxy > 0) & (wyx > 0), 4.0 / (1.0 / wxy + 1.0 / wyx), 0.0)
    W = np.stack([w_par, w_par, w_crs, w_crs], axis=-1)
    W = np.where(np.isfinite(S), W, 0.0)
    return S, W


def aggregate(S, W, sel, min_baselines):
    """
    Vector-average Stokes over the selected baselines, per scan. The error is
    the standard error of one real component, from the weighted scatter
    across baselines.

    Returns mean (scan, 4) complex, err (scan, 4), nbl (scan, 4).
    """
    w = np.where(sel[None, :, None], W, 0.0)
    Sz = np.where(w > 0, S, 0)
    wtot = w.sum(axis=1)
    nbl = (w > 0).sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = (w * Sz).sum(axis=1) / wtot
        var = (w * np.abs(Sz - mean[:, None, :]) ** 2).sum(axis=1) / wtot
        neff = wtot ** 2 / (w ** 2).sum(axis=1)
        err = np.sqrt(var / 2.0 / neff)
    bad = nbl < min_baselines
    mean[bad] = np.nan
    err[bad] = np.nan
    return mean, err, nbl


def view(mean, err, component, fractional):
    """Real-valued series to plot from complex aggregated Stokes."""
    m = mean.copy()
    e = err.copy()
    if fractional:
        with np.errstate(invalid='ignore', divide='ignore'):
            m[..., 1:] = mean[..., 1:] / mean[..., :1]
            e[..., 1:] = err[..., 1:] / np.abs(mean[..., :1])
    if component == 'real':
        return m.real, e
    if component == 'imag':
        return m.imag, e
    return np.abs(m), e


def weighted_rms_about_mean(x, e):
    """Inverse-variance weighted rms of a series about its weighted mean."""
    ok = np.isfinite(x) & np.isfinite(e) & (e > 0)
    if ok.sum() < 2:
        return np.nan
    w = 1.0 / e[ok] ** 2
    m = np.sum(w * x[ok]) / w.sum()
    return np.sqrt(np.sum(w * (x[ok] - m) ** 2) / w.sum())


def analyse(cache, bins, args, n_ant):
    """Array, per-antenna and leave-one-out series plus the driver metrics."""
    mean, wsum = corr_means(cache, bins)
    S, W = stokes_per_baseline(mean, wsum, list(cache['corr_names']), args.flip_v)
    bl = cache['baselines']
    all_bl = np.ones(len(bl), dtype=bool)

    arr_m, arr_e, _ = aggregate(S, W, all_bl, args.min_baselines)
    arr_x, arr_ex = view(arr_m, arr_e, args.component, args.fractional)

    n_scan = S.shape[0]
    ant_x = np.full((n_ant, n_scan, 4), np.nan)
    ant_ex = np.full_like(ant_x, np.nan)
    loo_x = np.full_like(ant_x, np.nan)
    loo_ex = np.full_like(ant_x, np.nan)
    chi2r = np.full((n_ant, 4), np.nan)
    loo_drop = np.full((n_ant, 4), np.nan)
    present = np.zeros(n_ant, dtype=bool)

    arr_rms = np.array([weighted_rms_about_mean(arr_x[:, s], arr_ex[:, s]) for s in range(4)])

    for a in range(n_ant):
        sel = (bl[:, 0] == a) | (bl[:, 1] == a)
        if not sel.any():
            continue
        m, e, _ = aggregate(S, W, sel, args.min_baselines)
        if not np.isfinite(m).any():
            continue
        present[a] = True
        ant_x[a], ant_ex[a] = view(m, e, args.component, args.fractional)
        m, e, _ = aggregate(S, W, ~sel, args.min_baselines)
        loo_x[a], loo_ex[a] = view(m, e, args.component, args.fractional)
        for s in range(4):
            rms_loo = weighted_rms_about_mean(loo_x[a, :, s], loo_ex[a, :, s])
            loo_drop[a, s] = 1.0 - rms_loo / arr_rms[s] if arr_rms[s] > 0 else np.nan

    # Everything below is relative to the array and robust: medians and
    # MADs over antennas, so a few bad antennas cannot drag the reference or
    # inflate the scale. MAD * 1.4826 is the Gaussian-equivalent sigma.
    def robust(x, axis):
        med = np.nanmedian(x, axis=axis, keepdims=True)
        mad = 1.4826 * np.nanmedian(np.abs(x - med), axis=axis, keepdims=True)
        return med, mad

    P = np.where(present)[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # R(t), per scan and Stokes, median over antennas
        ref_x = np.nanmedian(ant_x[P], axis=0)                       # (scan, 4)
        dev = ant_x - ref_x[None]                                    # D_a(t)

        # Static offset O_a, median over scans, then robust z across antennas
        offset = np.nanmedian(dev, axis=1)                           # (ant, 4)
        med_o, mad_o = robust(offset[P], axis=0)
        offset_z = (offset - med_o) / mad_o

        # Time evolution, offset removed, scaled by each scan's robust spread
        resid = dev - offset[:, None, :]                             # E_a(t)
        _, sig_t = robust(resid[P], axis=0)                          # (1, scan, 4)
        z = resid / sig_t
        n_ok = np.isfinite(z).sum(axis=1)
        chi2r = np.nansum(z ** 2, axis=1) / np.maximum(n_ok - 1, 1)
        chi2r[n_ok < 2] = np.nan

    for arr in (offset, offset_z, chi2r):
        arr[~present] = np.nan
    z[~present] = np.nan
    dev[~present] = np.nan

    return {
        'arr_x': arr_x, 'arr_ex': arr_ex, 'arr_rms': arr_rms,
        'ant_x': ant_x, 'ant_ex': ant_ex, 'loo_x': loo_x, 'loo_ex': loo_ex,
        'ref_x': ref_x, 'dev': dev, 'z': z,
        'chi2r': chi2r, 'loo_drop': loo_drop, 'present': present,
        'offset': offset, 'offset_z': offset_z,
    }


# =====================================================================
# PLOTS
# =====================================================================

def stokes_labels(args):
    unit = '' if args.fractional else ' [Jy]'
    comp = {'real': 'Re', 'imag': 'Im', 'amp': '|·|'}[args.component]
    labels = []
    for s in STOKES:
        q = s if (s == 'I' or not args.fractional) else f'{s}/I'
        labels.append(f'{comp} {q}{unit if (s == "I" or not args.fractional) else ""}')
    return labels


def top_antennas(res, s, n, threshold):
    """Up to n antennas with robust chi2r > threshold for Stokes index s, largest first."""
    c = res['chi2r'][:, s]
    idx = np.where(res['present'] & np.isfinite(c) & (c > threshold))[0]
    return list(idx[np.argsort(-c[idx])][:n])


def offset_outliers(res, s, sigma, exclude=()):
    """Antennas whose static offset exceeds sigma, largest first, for Stokes index s."""
    z = res['offset_z'][:, s]
    idx = np.where(res['present'] & np.isfinite(z) & (np.abs(z) > sigma))[0]
    return [a for a in idx[np.argsort(-np.abs(z[idx]))] if a not in exclude]


def antenna_colour(a, hl, off):
    if a in hl:
        return HIGHLIGHT_COLOURS[hl.index(a) % len(HIGHLIGHT_COLOURS)]
    if a in off:
        return OFFSET_COLOUR
    return None


def style_axis(ax):
    ax.grid(True, color='#e6e6e6', lw=0.6)
    ax.set_axisbelow(True)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)


def plot_stokes_vs_x(res, x, xlabel, ant_names, args, title, path):
    labels = stokes_labels(args)
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    for s, ax in enumerate(axes):
        style_axis(ax)
        hl = top_antennas(res, s, args.highlight, args.chi2_threshold)
        off = offset_outliers(res, s, args.offset_sigma, exclude=hl)
        for a in np.where(res['present'])[0]:
            if a in hl or a in off:
                continue
            ax.plot(x, res['ant_x'][a, :, s], '-', color=BACKGROUND_COLOUR, lw=0.7,
                    alpha=0.6, zorder=1)
        handles = [Line2D([], [], color=BACKGROUND_COLOUR, lw=1, label='other antennas')]
        for k, a in enumerate(hl):
            col = HIGHLIGHT_COLOURS[k % len(HIGHLIGHT_COLOURS)]
            ax.errorbar(x, res['ant_x'][a, :, s], yerr=res['ant_ex'][a, :, s], fmt='o-',
                        color=col, ms=4, lw=1.4, elinewidth=0.8, capsize=0, zorder=3)
            handles.append(Line2D([], [], color=col, marker='o', lw=1.4,
                                  label=f"{ant_names[a]}  χ²ᵣ={res['chi2r'][a, s]:.1f}"))
        for a in off:
            ax.plot(x, res['ant_x'][a, :, s], '.-', color=OFFSET_COLOUR, lw=0.9, ms=4, zorder=2)
            handles.append(Line2D([], [], color=OFFSET_COLOUR, marker='.', lw=0.9,
                                  label=f"{ant_names[a]}  offset {res['offset_z'][a, s]:+.1f}σ"))
        ax.errorbar(x, res['arr_x'][:, s], yerr=res['arr_ex'][:, s], fmt='s-',
                    color=ARRAY_COLOUR, ms=4, lw=2, elinewidth=1, capsize=0, zorder=4)
        handles.append(Line2D([], [], color=ARRAY_COLOUR, marker='s', lw=2, label='array'))
        ax.set_ylabel(labels[s])
        ax.legend(handles=handles, fontsize=8, loc='upper left', bbox_to_anchor=(1.01, 1),
                  frameon=False)
    axes[-1].set_xlabel(xlabel)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def plot_heatmap(res, scans, ant_names, args, title, path):
    """
    z_a(t), each antenna's deviation from the per-scan antenna median with
    its static offset removed, in units of that scan's robust spread across
    antennas. One row per antenna.
    """
    labels = stokes_labels(args)
    ants = np.where(res['present'])[0]
    fig, axes = plt.subplots(1, 4, figsize=(16, max(6, 0.17 * len(ants) + 2)), sharey=True)
    for s, ax in enumerate(axes):
        z = res['z'][ants, :, s]
        im = ax.imshow(z, aspect='auto', cmap='RdBu_r', vmin=-5, vmax=5,
                       interpolation='nearest')
        ax.set_title(labels[s], fontsize=10)
        ax.set_xticks(range(len(scans)))
        ax.set_xticklabels(scans, rotation=90, fontsize=6)
        ax.set_xlabel('scan')
    axes[0].set_yticks(range(len(ants)))
    axes[0].set_yticklabels(ant_names[ants], fontsize=6)
    cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.01)
    cb.set_label('z = (D − median_t D) / (1.4826 MAD over antennas),  D = antenna − array median')
    fig.suptitle(title, fontsize=11)
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def plot_ranking(res, ant_names, args, title, path):
    """chi2r and leave-one-out drop in temporal rms, per antenna, per Stokes."""
    labels = stokes_labels(args)
    ants = np.where(res['present'])[0]
    xpos = np.arange(len(ants))
    fig, axes = plt.subplots(4, 3, figsize=(20, 12), sharex=True)
    for s in range(4):
        hl = top_antennas(res, s, args.highlight, args.chi2_threshold)
        off = offset_outliers(res, s, args.offset_sigma, exclude=hl)
        cols = [antenna_colour(a, hl, off) or BACKGROUND_COLOUR for a in ants]
        for col, (ax, key, ylab) in enumerate((
                (axes[s, 0], 'chi2r', 'robust χ²ᵣ (time evolution)'),
                (axes[s, 1], 'offset_z', 'static offset [robust σ]'),
                (axes[s, 2], 'loo_drop', '1 − rms(without) / rms(all)'))):
            style_axis(ax)
            ax.bar(xpos, res[key][ants, s], color=cols, width=0.8, edgecolor='white', lw=0.5)
            ax.set_ylabel(f'{labels[s]}\n{ylab}', fontsize=9)
            if key == 'chi2r':
                ax.axhline(1, color=ARRAY_COLOUR, lw=0.8, ls=':')
                ax.axhline(args.chi2_threshold, color=ARRAY_COLOUR, lw=0.8, ls='--')
                ax.set_yscale('log')
            elif key == 'offset_z':
                for sign in (-1, 1):
                    ax.axhline(sign * args.offset_sigma, color=ARRAY_COLOUR, lw=0.8, ls='--')
                ax.axhline(0, color=ARRAY_COLOUR, lw=0.8)
            else:
                ax.axhline(0, color=ARRAY_COLOUR, lw=0.8)
    for ax in axes[-1]:
        ax.set_xticks(xpos)
        ax.set_xticklabels(ant_names[ants], rotation=90, fontsize=6)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def plot_antenna_grid(res, x, xlabel, ant_names, args, s, title, path):
    """Small multiples, one panel per antenna, against the per-scan antenna median."""
    labels = stokes_labels(args)
    ants = np.where(res['present'])[0]
    ncol = int(np.ceil(np.sqrt(len(ants))))
    nrow = int(np.ceil(len(ants) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.2 * ncol, 1.7 * nrow),
                             sharex=True, sharey=True, squeeze=False)
    hl = top_antennas(res, s, args.highlight, args.chi2_threshold)
    off = offset_outliers(res, s, args.offset_sigma, exclude=hl)
    for k, ax in enumerate(axes.flat):
        if k >= len(ants):
            ax.axis('off')
            continue
        a = ants[k]
        col = antenna_colour(a, hl, off) or ORDINARY_COLOUR
        ax.plot(x, res['ref_x'][:, s], '-', color=BACKGROUND_COLOUR, lw=1.5)
        ax.errorbar(x, res['ant_x'][a, :, s], yerr=res['ant_ex'][a, :, s], fmt='o',
                    color=col, ms=2.5, elinewidth=0.6, capsize=0)
        ax.set_title(f"{ant_names[a]}  χ²ᵣ={res['chi2r'][a, s]:.1f}  "
                     f"off={res['offset_z'][a, s]:+.1f}σ", fontsize=7,
                     fontweight='bold' if (a in hl or a in off) else 'normal')
        ax.tick_params(labelsize=6)
        ax.grid(True, color='#eeeeee', lw=0.5)
    fig.supxlabel(xlabel, fontsize=9)
    fig.supylabel(f'{labels[s]}  (points antenna, grey antenna median)', fontsize=9)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def write_ranking(res, ant_names, path):
    ants = np.where(res['present'])[0]
    with open(path, 'w') as f:
        f.write('# antenna  ' + '  '.join(f'chi2r_{s:<3s} offz_{s:<4s} drop_{s:<4s}'
                                          for s in STOKES) + '\n')
        f.write('# sorted by the larger of max chi2r and max |offset_z| over Stokes\n')
        score = np.fmax(np.nanmax(res['chi2r'][ants], axis=1),
                        np.nanmax(np.abs(res['offset_z'][ants]), axis=1))
        order = ants[np.argsort(-score)]
        for a in order:
            vals = '  '.join(f"{res['chi2r'][a, s]:9.2f} {res['offset_z'][a, s]:+9.2f} "
                             f"{res['loo_drop'][a, s]:+9.3f}" for s in range(4))
            f.write(f'{ant_names[a]:<9s}  {vals}\n')


def log_ranking(res, ant_names, n, threshold, offset_sigma):
    for s, name in enumerate(STOKES):
        order = top_antennas(res, s, n, threshold)
        summary = ', '.join(f"{ant_names[a]} (χ²ᵣ={res['chi2r'][a, s]:.1f}, "
                            f"drop={res['loo_drop'][a, s]:+.2f})" for a in order)
        log(f"  {name} time-variable: {summary or f'none with χ²ᵣ > {threshold:g}'}")
        off = offset_outliers(res, s, offset_sigma)
        if off:
            log(f"  {name} static offset > {offset_sigma:g}σ: "
                + ', '.join(f"{ant_names[a]} ({res['offset'][a, s]:+.4g}, "
                            f"{res['offset_z'][a, s]:+.1f}σ)" for a in off))


def make_plot_set(cache, bins, args, ant_names, x, xlabel, tag, label):
    res = analyse(cache, bins, args, len(ant_names))
    if not res['present'].any():
        log(f"  {label}: no antenna has data, skipped")
        return
    freq = cache['bin_freq'][bins]
    title = (f"{tag}   {label}   {freq.min() / 1e9:.3f}-{freq.max() / 1e9:.3f} GHz   "
             f"{args.data_column}")
    base = os.path.join(args.save_dir, f'{tag}_{label}{args.output_suffix}')

    log(f"{label}: array temporal rms "
        + ', '.join(f'{s}={r:.4g}' for s, r in zip(STOKES, res['arr_rms'])))
    log_ranking(res, ant_names, args.highlight, args.chi2_threshold, args.offset_sigma)

    plot_stokes_vs_x(res, x, xlabel, ant_names, args, title, f'{base}_stokes_vs_{args.xaxis}.png')
    plot_heatmap(res, cache['scans'], ant_names, args, title, f'{base}_antenna_scan_heatmap.png')
    plot_ranking(res, ant_names, args, title, f'{base}_ranking.png')
    write_ranking(res, ant_names, f'{base}_ranking.txt')
    if not args.no_grid:
        for s, name in enumerate(STOKES):
            plot_antenna_grid(res, x, xlabel, ant_names, args, s, title,
                              f'{base}_grid_{name}_vs_{args.xaxis}.png')
    log(f"  plots written to {base}_*")


# =====================================================================
# MAIN
# =====================================================================

def main():
    args = get_config()
    global ENABLE_PRINTING
    ENABLE_PRINTING = not args.no_printing

    field_names, phase_dirs = load_field_table(args.ms_file)
    field_id, field_name = resolve_field(field_names, args.field)
    tag = f'{scrub_name(field_name)}'

    os.makedirs(args.log_dir, exist_ok=True)
    os.makedirs(args.save_dir, exist_ok=True)
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(args.log_dir, f'stokes-diagnostics_{tag}_{stamp}{args.output_suffix}.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='[%(asctime)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', filemode='w', force=True)

    n_cpu, cpu_src = get_allocated_cpus()
    dask.config.set(num_workers=n_cpu)

    log("=" * 70)
    log("STOKES ANTENNA DIAGNOSTICS")
    log("=" * 70)
    log(f"  MS: {args.ms_file}")
    log(f"  Field: {field_name} (ID {field_id})")
    log(f"  Data column: {args.data_column}   weights: {args.weight_column}")
    log(f"  Component: {args.component}   fractional: {args.fractional}   "
        f"flip V: {args.flip_v}   x-axis: {args.xaxis}")
    log(f"  Dask threads: {n_cpu} [{cpu_src}]")
    log(f"  Log: {log_file}")
    log("")

    ant_names, ant_pos = load_antenna_table(args.ms_file)
    chan_freq = load_chan_freq(args.ms_file)
    corr_names = load_corr_names(args.ms_file)
    chan_lo, chan_hi = parse_chan_range(args.chan_range, len(chan_freq))

    path, key = cache_path(args, field_id, field_name, chan_lo, chan_hi)
    cache = None
    if os.path.exists(path) and not args.overwrite_cache:
        cache = load_cache(path, key)
        if cache is not None:
            log(f"Loaded cache {path}")
    if cache is None:
        log(f"Building cache {path}")
        cache = build_cache(args, field_id, corr_names, chan_freq, chan_lo, chan_hi)
        save_cache(path, key, cache)
        log(f"Cache written ({os.path.getsize(path) / 1024**2:.1f} MB)")

    n_scan, n_bl, n_bin, _ = cache['vsum'].shape
    log(f"{n_scan} scans, {n_bl} baselines, {n_bin} channel bins "
        f"({cache['bin_freq'][0] / 1e9:.3f}-{cache['bin_freq'][-1] / 1e9:.3f} GHz centres)")
    log("")

    # Horizontal axis
    scan_time = cache['scan_time']
    xaxis = args.xaxis
    if xaxis == 'parang':
        x = scan_parallactic_angles(ant_pos, phase_dirs[field_id], scan_time)
        xlabel = 'parallactic angle [deg]'
        if x is None:
            log("Falling back to XAXIS = 'time'")
            xaxis = 'time'
    if xaxis == 'time':
        x = (scan_time - scan_time.min()) / 3600.0
        t0 = datetime.fromtimestamp(scan_time.min() - 3506716800.0, tz=timezone.utc)
        xlabel = f'hours since {t0:%Y-%m-%d %H:%M} UTC'
    elif xaxis == 'scan':
        x = cache['scans'].astype(float)
        xlabel = 'scan number'
    args.xaxis = xaxis

    if isinstance(args.plot_bins, str) and args.plot_bins.strip().lower() == 'all':
        bins = list(range(n_bin))
    else:
        bins = sorted({int(b) for b in np.atleast_1d(args.plot_bins)})
        if not bins or min(bins) < 0 or max(bins) >= n_bin:
            raise ValueError(f"PLOT_BINS {args.plot_bins} outside 0-{n_bin - 1}")
    label = 'band' if len(bins) == n_bin else 'bins' + '-'.join(map(str, bins))
    make_plot_set(cache, bins, args, ant_names, x, xlabel, tag, label)

    if args.per_bin and n_bin > 1:
        for b in range(n_bin):
            log("")
            make_plot_set(cache, [b], args, ant_names, x, xlabel, tag, f'bin{b:02d}')

    log("")
    log("Done")


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        log(f"ERROR: {exc}")
        raise
