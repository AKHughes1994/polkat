#!/usr/bin/env python3
"""
Median XX and YY amplitude per time bin for one field, and a three-panel plot
of it: XX and YY, the Stokes I proxy (XX+YY)/2, and the feed-frame Stokes Q
proxy (XX-YY)/2.

The median in each time bin is taken over every unflagged baseline and
channel in the frequency range. Its error is the standard error on a median,
1.2533 * 1.4826 * MAD / sqrt(N).

The medians are saved to a text file. If that file already exists the
processing is skipped and only the plot is made; delete the file to
reprocess.
"""

import os
import re
import sys
from datetime import datetime, timedelta

import numpy as np
import dask
import psutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import AutoMinorLocator
from daskms import xds_from_ms, xds_from_table


# ---------------------------------------------------------------------------- #
# ---------------------------------- INPUTS ---------------------------------- #
# ---------------------------------------------------------------------------- #

MS_FILE = '1780433330-sdp-l0_2026-06-18T14-28-58_j6E_1024ch.ms'
FIELD = 'SS433'                    # field name, or field index as an int or string
TIME_BIN_S = 16.0             # time bin width in seconds
FREQ_RANGE_GHZ = (1.1, 1.4)   # (min, max) in GHz; () or None for the whole band
SAVE_DIR = os.getcwd()
DATA_COLUMN = 'CORRECTED_DATA'

# Fraction of the available memory the processing may use. Half of it goes to
# the batch of time bins being reduced, half to the dask chunks being read.
MEMORY_SAFE_FRACTION = 0.85

# Bytes held per (row, channel, correlation) while a batch is processed:
# complex64 visibility and bool flag as read, then float32 amplitude and bool
# mask for the selected correlations, with headroom for dask's copies.
MEMORY_OVERHEAD = 2.5

# MS CORR_TYPE codes for the linear feeds
CORR_CODES = {9: 'XX', 10: 'XY', 11: 'YX', 12: 'YY'}

# Scale factors for the median's standard error from the MAD
MAD_TO_SIGMA = 1.4826
MEDIAN_SE_FACTOR = 1.2533

# MAD-based outlier clipping applied to the finished per-bin medians and
# errors before they are saved or plotted (see mad_clip_bins)
CLIP_SIGMA = 10.0
CLIP_ON = 'error'   # 'data' (med_xx/med_yy), 'error' (err_xx/err_yy), or 'both'

# MS TIME is in MJD seconds, counted from this epoch
MJD_EPOCH = datetime(1858, 11, 17)


# ---------------------------------------------------------------------------- #
# ------------------------------ RESOURCES ----------------------------------- #
# ---------------------------------------------------------------------------- #

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
    """Memory allocated to this job by SLURM, in bytes, or None."""
    mem_per_node = os.environ.get('SLURM_MEM_PER_NODE')
    mem_per_cpu = os.environ.get('SLURM_MEM_PER_CPU')
    try:
        if mem_per_node:
            return _parse_slurm_mem(mem_per_node)
        if mem_per_cpu:
            cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 1))
            ntasks = int(os.environ.get('SLURM_NTASKS', 1))
            return _parse_slurm_mem(mem_per_cpu) * cpus * ntasks
    except ValueError:
        pass
    return None


def get_allocated_cpus():
    """CPUs this job may actually use, honouring SLURM."""
    cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    ntasks = os.environ.get('SLURM_NTASKS')
    job_cpus = re.match(r'(\d+)', os.environ.get('SLURM_JOB_CPUS_PER_NODE', ''))
    try:
        if cpus_per_task:
            return int(cpus_per_task) * int(ntasks or 1)
        if job_cpus:
            return int(job_cpus.group(1))
        if ntasks:
            return int(ntasks)
    except ValueError:
        pass
    return os.cpu_count() or 1


def available_memory_bytes():
    """
    Memory this job can use: the more restrictive of the SLURM allocation and
    what the node has free. Inside a batch job psutil describes the node, not
    the job, so on a shared node it alone would overstate what is available.
    """
    node_available = psutil.virtual_memory().available
    slurm_mem = get_slurm_memory_bytes()
    return min(slurm_mem, node_available) if slurm_mem else node_available


def check_memory_critical():
    """
    Raise MemoryError if this job is above 95% of its memory. Inside a SLURM
    allocation this process is measured against what it was given; outside
    one, the node's usage is used.
    """
    slurm_mem = get_slurm_memory_bytes()
    if slurm_mem:
        percent_used = 100.0 * psutil.Process().memory_info().rss / slurm_mem
    else:
        percent_used = psutil.virtual_memory().percent

    if percent_used > 95:
        raise MemoryError(f'CRITICAL: memory usage at {percent_used:.1f}%')
    if percent_used > 85:
        print(f'WARNING: memory usage at {percent_used:.1f}%')


# ---------------------------------------------------------------------------- #
# ------------------------------- MS HELPERS --------------------------------- #
# ---------------------------------------------------------------------------- #

def resolve_field(ms_file, field):
    """Field name or index -> (field_id, field_name)."""
    names = xds_from_table(f'{ms_file}::FIELD', columns=['NAME'])[0].NAME.values
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in names]

    try:
        field_id = int(field)
    except ValueError:
        lowered = [n.lower().strip() for n in names]
        if str(field).lower().strip() not in lowered:
            raise ValueError(f"Field '{field}' not found. Available: "
                             f"{', '.join(f'{i}={n}' for i, n in enumerate(names))}")
        field_id = lowered.index(str(field).lower().strip())

    if not 0 <= field_id < len(names):
        raise ValueError(f'Field ID {field_id} out of range. Available: '
                         f"{', '.join(f'{i}={n}' for i, n in enumerate(names))}")
    return field_id, names[field_id]


def select_channels(ms_file, freq_range_ghz):
    """
    Channel range [c0, c1) of the first SPW inside `freq_range_ghz`, and the
    SPW's channel frequencies in Hz. An empty range selects the whole band.
    """
    chan_freq = xds_from_table(f'{ms_file}::SPECTRAL_WINDOW',
                               columns=['CHAN_FREQ'])[0].CHAN_FREQ.values[0]
    if not freq_range_ghz:
        return 0, len(chan_freq), chan_freq

    fmin, fmax = sorted(f * 1e9 for f in freq_range_ghz)
    inside = np.flatnonzero((chan_freq >= fmin) & (chan_freq <= fmax))
    if inside.size == 0:
        raise ValueError(f'No channels between {fmin/1e9:.4f} and {fmax/1e9:.4f} GHz; '
                         f'the band covers {chan_freq.min()/1e9:.4f} - '
                         f'{chan_freq.max()/1e9:.4f} GHz')
    return int(inside[0]), int(inside[-1]) + 1, chan_freq


def parallel_hand_indices(ms_file):
    """Correlation indices of XX and YY."""
    corr_types = xds_from_table(f'{ms_file}::POLARIZATION',
                                columns=['CORR_TYPE'])[0].CORR_TYPE.values[0]
    names = [CORR_CODES.get(int(c), str(c)) for c in corr_types]
    missing = [c for c in ('XX', 'YY') if c not in names]
    if missing:
        raise ValueError(f"Correlation(s) {', '.join(missing)} not in the MS "
                         f"(it has {', '.join(names)})")
    return names.index('XX'), names.index('YY'), len(names)


# ---------------------------------------------------------------------------- #
# ------------------------------- PROCESSING --------------------------------- #
# ---------------------------------------------------------------------------- #

def median_and_error(values):
    """Median of `values` and its standard error; NaN for an empty array."""
    if values.size == 0:
        return np.nan, np.nan
    med = np.median(values)
    mad = np.median(np.abs(values - med))
    return med, MEDIAN_SE_FACTOR * MAD_TO_SIGMA * mad / np.sqrt(values.size)


def mad_clip_bins(data, sigma=CLIP_SIGMA, mode=CLIP_ON):
    """
    Drop time bins whose XX/YY median amplitude and/or its error is a
    MAD-based outlier against the rest of the run.

    `mode` selects which columns are tested: 'data' tests med_xx/med_yy only,
    'error' tests err_xx/err_yy only, 'both' tests all four. Each tested
    column is checked independently against its own median and MAD; a bin
    failing any of them is dropped entirely, not just in that column, so XX,
    YY, and the derived I/Q panels stay on the same time axis. A bin with no
    unflagged data (NaN) is left alone: a NaN comparison is never True, so it
    can't trigger the clip. A column with zero MAD (no spread to measure) is
    skipped rather than flagging every value that isn't exactly the median.

    Returns
    -------
    tuple
        (kept rows, {column name: number of bins it flagged})
    """
    all_cols = {'med_xx': 2, 'err_xx': 3, 'med_yy': 5, 'err_yy': 6}
    if mode == 'data':
        col_idx = {k: v for k, v in all_cols.items() if k.startswith('med_')}
    elif mode == 'error':
        col_idx = {k: v for k, v in all_cols.items() if k.startswith('err_')}
    elif mode == 'both':
        col_idx = all_cols
    else:
        raise ValueError(f"CLIP_ON must be 'data', 'error', or 'both' (got: '{mode}')")

    bad = np.zeros(len(data), dtype=bool)
    counts = {}
    for name, idx in col_idx.items():
        values = data[:, idx]
        med = np.nanmedian(values)
        mad = np.nanmedian(np.abs(values - med))
        threshold = sigma * MAD_TO_SIGMA * mad
        outlier = (np.abs(values - med) > threshold) if threshold > 0 \
            else np.zeros(len(values), dtype=bool)
        counts[name] = int(outlier.sum())
        bad |= outlier
    return data[~bad], counts


def get_dump_time(ms_file, field_id):
    """Median integration (dump) time for this field's cross-correlations, in seconds."""
    ds_list = xds_from_ms(ms_file, columns=['INTERVAL'],
                          taql_where=f'FIELD_ID=={field_id} AND ANTENNA1!=ANTENNA2',
                          group_cols=['DATA_DESC_ID'])
    if not ds_list:
        raise RuntimeError('No cross-correlation rows for this field')
    interval = dask.compute(ds_list[0].INTERVAL.data)[0]
    return float(np.median(interval))


def snap_time_bin(bin_s, dump_s):
    """
    Round `bin_s` to the nearest whole number of dump times, with a minimum of
    one dump -- a bin narrower than one dump can't hold a whole dump, and one
    that isn't a whole multiple would mix a partial dump in at each edge.

    Returns
    -------
    tuple
        (snapped bin width in seconds, number of dumps per bin)
    """
    n_dumps = max(1, round(bin_s / dump_s))
    return n_dumps * dump_s, n_dumps


def time_bin_bounds(time, bin_s):
    """
    Split time-sorted rows into bins of `bin_s` seconds from the first
    timestamp. Returns the bin centres and the [start, stop) row range of each
    non-empty bin.
    """
    bin_idx = np.floor((time - time[0]) / bin_s).astype(np.int64)
    starts = np.concatenate([[0], np.flatnonzero(np.diff(bin_idx)) + 1])
    stops = np.concatenate([starts[1:], [len(time)]])
    centres = time[0] + (bin_idx[starts] + 0.5) * bin_s
    return centres, starts, stops


def batch_bins(starts, stops, bytes_per_row, budget_bytes, bin_s):
    """
    Group consecutive bins into batches whose rows fit in `budget_bytes`.
    Returns a list of (first_bin, last_bin_exclusive).
    """
    batches = []
    first = 0
    for k in range(len(starts)):
        bin_bytes = (stops[k] - starts[k]) * bytes_per_row
        if bin_bytes > budget_bytes:
            raise MemoryError(
                f'One {bin_s:.4f} s bin needs ~{bin_bytes/1024**3:.2f} GB, over the '
                f'{budget_bytes/1024**3:.2f} GB budget. Narrow FREQ_RANGE_GHZ or '
                f'TIME_BIN_S, or run with more memory.')
        if (stops[k] - starts[first]) * bytes_per_row > budget_bytes:
            batches.append((first, k))
            first = k
    batches.append((first, len(starts)))
    return batches


def compute_medians(ms_file, field_id, c0, c1, n_chan_total, bin_s):
    """
    Median XX and YY amplitude and error per time bin. Returns an array with
    columns time_mjd_s, n_xx, med_xx, err_xx, n_yy, med_yy, err_yy.
    """
    i_xx, i_yy, n_corr = parallel_hand_indices(ms_file)

    n_cpus = get_allocated_cpus()
    dask.config.set(scheduler='threads', num_workers=n_cpus)
    print(f'Dask capped to {n_cpus} thread(s)')

    # Every read holds all correlations of the selected channels for its rows.
    n_chan = c1 - c0
    bytes_per_row = n_chan * n_corr * (8 + 1) * MEMORY_OVERHEAD
    budget = available_memory_bytes() * MEMORY_SAFE_FRACTION
    batch_budget = budget / 2

    # One chunk per thread can be in flight, so the row chunk is sized for all
    # of them to fit in the read half of the budget.
    row_chunk = max(1, int((budget / 2) / (n_cpus * bytes_per_row)))

    # The selected channels are one chunk of their own, so no read touches a
    # channel outside the range.
    chan_chunks = tuple(n for n in (c0, n_chan, n_chan_total - c1) if n > 0)

    ds_list = xds_from_ms(ms_file,
                          columns=[DATA_COLUMN, 'FLAG', 'TIME'],
                          taql_where=f'FIELD_ID=={field_id} AND ANTENNA1!=ANTENNA2',
                          group_cols=['DATA_DESC_ID'],
                          index_cols=['TIME'],
                          chunks={'row': row_chunk, 'chan': chan_chunks})
    if not ds_list:
        raise RuntimeError('No cross-correlation rows for this field')
    if len(ds_list) > 1:
        print(f'WARNING: {len(ds_list)} DATA_DESC_IDs found; using the first only')
    ds = ds_list[0]

    print('Reading TIME...')
    time = ds.TIME.values
    if np.any(np.diff(time) < 0):
        raise RuntimeError('Rows are not in time order')

    centres, starts, stops = time_bin_bounds(time, bin_s)
    batches = batch_bins(starts, stops, bytes_per_row, batch_budget, bin_s)

    print(f'{len(time):,} rows in {len(starts)} time bin(s), {n_chan} channel(s)')
    print(f'Memory budget {budget/1024**3:.2f} GB: {len(batches)} batch(es), '
          f'{row_chunk:,}-row read chunks')

    # Sliced directly rather than via Dataset.isel(), which daskms's own
    # lightweight Dataset (unlike xarray's) does not implement.
    vis = ds[DATA_COLUMN].data[:, c0:c1, :]
    flag = ds.FLAG.data[:, c0:c1, :]
    rows = []

    for b, (k0, k1) in enumerate(batches, start=1):
        r0, r1 = starts[k0], stops[k1 - 1]
        print(f'Batch {b}/{len(batches)}: bins {k0}-{k1 - 1}, rows {r0:,}-{r1:,}')

        batch_vis, batch_flag = dask.compute(vis[r0:r1][..., [i_xx, i_yy]],
                                             flag[r0:r1][..., [i_xx, i_yy]])
        try:
            check_memory_critical()
        except MemoryError as e:
            print(f'{e} -- stopping; saving the {len(rows)} bin(s) done so far')
            break

        for k in range(k0, k1):
            s0, s1 = starts[k] - r0, stops[k] - r0
            row = [centres[k]]
            for p in (0, 1):
                amp = np.abs(batch_vis[s0:s1, :, p])
                good = ~batch_flag[s0:s1, :, p] & np.isfinite(amp) & (amp > 0)
                med, err = median_and_error(amp[good])
                row += [int(good.sum()), med, err]
            rows.append(row)

        del batch_vis, batch_flag

    return np.array(rows, dtype=float)


def save_medians(path, data, field_name, c0, c1, chan_freq, bin_s,
                 n_before, clip_counts):
    n_removed = n_before - len(data)
    clip_lines = '\n'.join(f'    {name}: {n}' for name, n in clip_counts.items())
    header = (f'MS: {MS_FILE}\n'
              f'Field: {field_name}\n'
              f'Data column: {DATA_COLUMN}\n'
              f'Time bin: {bin_s:.4f} s\n'
              f'Channels: {c0}-{c1 - 1} ({chan_freq[c0]/1e9:.4f} - '
              f'{chan_freq[c1 - 1]/1e9:.4f} GHz)\n'
              f'Outlier clipping: {CLIP_SIGMA:g}-sigma MAD, tested independently per column\n'
              f'  {n_removed} of {n_before} time bin(s) removed for being an outlier in at '
              f'least one column (counts below can overlap, so may sum to more than '
              f'{n_removed}):\n'
              f'{clip_lines}\n'
              'time_mjd_s n_xx med_xx err_xx n_yy med_yy err_yy')
    np.savetxt(path, data, header=header,
               fmt=['%.3f', '%d', '%.6e', '%.6e', '%d', '%.6e', '%.6e'])
    print(f'Saved: {path}')


# ---------------------------------------------------------------------------- #
# -------------------------------- PLOTTING ---------------------------------- #
# ---------------------------------------------------------------------------- #

def style_axis(ax):
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=11)
    ax.tick_params(which='major', length=7, width=1.5)
    ax.tick_params(which='minor', length=4, width=1.0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)


def plot_medians(txt_path, png_path, title):
    t, _, xx, exx, _, yy, eyy = np.loadtxt(txt_path, unpack=True, ndmin=2)
    times = [MJD_EPOCH + timedelta(seconds=s) for s in t]

    stokes_i = (xx + yy) / 2
    stokes_q = (xx - yy) / 2
    err_iq = np.sqrt(exx**2 + eyy**2) / 2

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True,
                             gridspec_kw={'hspace': 0})
    marker = dict(fmt='o', ms=3, capsize=0, elinewidth=1)

    axes[0].errorbar(times, xx, yerr=exx, color='tab:blue', label='XX', **marker)
    axes[0].errorbar(times, yy, yerr=eyy, color='tab:orange', label='YY', **marker)
    axes[0].set_ylabel('Median amplitude')
    axes[0].legend(frameon=False, fontsize=11)
    axes[0].set_title(title, fontsize=13)

    axes[1].errorbar(times, stokes_i, yerr=err_iq, color='black', **marker)
    axes[1].set_ylabel('(XX + YY) / 2')

    axes[2].errorbar(times, stokes_q, yerr=err_iq, color='tab:red', **marker)
    axes[2].axhline(0, color='grey', lw=1, ls='--')
    axes[2].set_ylabel('(XX − YY) / 2')
    axes[2].set_xlabel('Time (UTC)')
    axes[2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    for ax in axes:
        style_axis(ax)
        ax.yaxis.label.set_size(12)
    axes[2].xaxis.label.set_size(12)
    fig.align_ylabels(axes)

    fig.savefig(png_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {png_path}')


# ---------------------------------------------------------------------------- #
# ---------------------------------- MAIN ------------------------------------ #
# ---------------------------------------------------------------------------- #

def main():
    if CLIP_ON not in ('data', 'error', 'both'):
        sys.exit(f"CLIP_ON must be 'data', 'error', or 'both' (got: '{CLIP_ON}')")

    field_id, field_name = resolve_field(MS_FILE, FIELD)
    c0, c1, chan_freq = select_channels(MS_FILE, FREQ_RANGE_GHZ)

    dump_s = get_dump_time(MS_FILE, field_id)
    bin_s, n_dumps = snap_time_bin(TIME_BIN_S, dump_s)

    ms_base = os.path.basename(os.path.normpath(MS_FILE)).replace('.ms', '')
    safe_field = re.sub(r'[^A-Za-z0-9._+-]', '_', field_name)
    stem = os.path.join(SAVE_DIR, f'med_amp_{ms_base}_{safe_field}_{bin_s:g}s')
    txt_path, png_path = f'{stem}.txt', f'{stem}.png'

    print('=' * 60)
    print(f'  MS file:      {MS_FILE}')
    print(f'  Field:        {field_name} (ID {field_id})')
    print(f'  Data column:  {DATA_COLUMN}')
    print(f'  Dump time:    {dump_s:.4f} s')
    print(f'  Time bin:     {TIME_BIN_S:g} s requested -> {bin_s:.4f} s '
          f'({n_dumps} dump(s)/bin)')
    print(f'  Frequency:    {chan_freq[c0]/1e9:.4f} - {chan_freq[c1 - 1]/1e9:.4f} GHz '
          f'(channels {c0}-{c1 - 1}'
          f"{', whole band' if not FREQ_RANGE_GHZ else ''})")
    print(f'  Outlier clip: {CLIP_SIGMA:g}σ MAD on {CLIP_ON}')
    print(f'  Save dir:     {SAVE_DIR}')
    print(f'  Output:       {txt_path}')
    print(f'                {png_path}')
    print('=' * 60)

    try:
        answer = input('Proceed? [y/n]: ').strip().lower()
    except EOFError:
        answer = ''
    if answer not in ('y', 'yes'):
        sys.exit('Aborted.')

    os.makedirs(SAVE_DIR, exist_ok=True)

    if os.path.isfile(txt_path):
        print(f'{txt_path} already exists -- skipping to plotting. '
              f'Delete it if reprocessing is necessary.')
    else:
        data = compute_medians(MS_FILE, field_id, c0, c1, len(chan_freq), bin_s)
        if data.size == 0:
            sys.exit('No time bins processed; nothing to save.')

        # Clipped once here, against the full set of bins from every batch --
        # not per-batch, since the population median/MAD it clips against
        # would otherwise be built from whatever partial set happened to be
        # in memory at the time.
        n_before = len(data)
        data, clip_counts = mad_clip_bins(data)
        n_removed = n_before - len(data)
        if n_removed:
            print(f'MAD clip ({CLIP_SIGMA:g}σ): removed {n_removed} of {n_before} bin(s) -- '
                  + ', '.join(f'{name}: {n}' for name, n in clip_counts.items()))

        save_medians(txt_path, data, field_name, c0, c1, chan_freq, bin_s,
                    n_before, clip_counts)

    title = (f'{field_name}  |  {chan_freq[c0]/1e9:.3f}–{chan_freq[c1 - 1]/1e9:.3f} GHz'
             f'  |  {bin_s:.3g} s bins')
    plot_medians(txt_path, png_path, title)


if __name__ == '__main__':
    main()
