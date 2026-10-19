# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
RMSYNTH step 1: fit Stokes I, P, Q, U, V of every source in rmsynth_info.json
on its MFS and per-channel images, and write the inputs for RM synthesis.

For each prefix (one per time interval, or one for the whole observation):
  - MFS: fit Stokes I (the reference position), then P, Q and U (anchored to
    P) and V (anchored to I). A fit that drifts more than 1/3 beam from its
    reference is redone with positions fixed.
  - Channels: the same sequence for every channel, in a pool of workers.

Outputs (RESULTS/):
  <prefix>_polarization.json           every MFS and channel fit
  <prefix>_<component>_rmsynth.txt     freq, I, Q, U, dI, dQ, dU (Hz, Jy)
  <prefix>_<component>_IQUV.txt        freq, I, Q, U, V, dI, dQ, dU, dV
  <prefix>_<component>_rmsynth_UVfix.txt  U/V-merged spectrum, only without
                                       a polarization angle calibrator
  fitting_plots/                       IQUV spectra; MFS light curves for
                                       time-resolved runs

Runs under python-pycasa (casatasks/casatools are imported directly).
"""

import glob
import os
import re
import datetime
import sys
import json
import numpy as np
import os.path as o
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.spatial.distance import cdist
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

from casatasks import imfit, imstat, imhead
from casatools import table
tb = table()


# =============================================================================
# SETTINGS (defaults from oxkat/config.py)
# =============================================================================

FORCE_FIX_STOKES_V    = cfg.RMSYN_FORCE_FIX_STOKES_V     # always fix Stokes V positions to Stokes I
OVERWRITE             = cfg.RMSYN_OVERWRITE              # False: reuse an existing JSON, only re-plot
SPEC_PLOT_SNR_THRESH  = cfg.RMSYN_SPEC_PLOT_SNR_THRESH   # MFS S/N needed for a spectrum plot (time-resolved)
SPEC_INDEX_SNR_THRESH = cfg.RMSYN_SPEC_INDEX_SNR_THRESH  # MFS S/N needed for a spectral index fit
SPEC_INDEX_MAD_CLIP   = cfg.RMSYN_SPEC_INDEX_MAD_CLIP    # Stokes I outlier clip (robust sigma)
MAX_I_DRIFT_PIX       = cfg.RMSYN_MAX_I_DRIFT_PIX        # max MFS Stokes I drift from the input position (False = off)
USE_RESIDUAL_RMS      = cfg.RMSYN_USE_RESIDUAL_RMS       # measure noise on the matching residual image

# Plot y ranges cover points within this many robust sigmas of the median
# (+/- error bars), padded by this fraction of the range
PLOT_YLIM_MAD_SIGMA   = 10.0
PLOT_YLIM_PAD         = 0.25

# S/N at or above which a component's position is fitted freely rather than
# fixed (check_position): Q, U and V images, and Plin/Ptot images
FREE_FIT_SNR          = 5.0
P_FREE_FIT_SNR        = 7.0

# No polarization angle calibrator: a linear cross-hand phase fit merges U and
# V into one spectrum (see write_uvfix_rmsynth)
UVFIX_MAX_DELAY_NS    = 10.0  # half-width of the cross-hand delay search (ns)
UVFIX_MIN_SNR         = 6.0   # coherent S/N below which the fit is flagged unreliable
UVFIX_OBSERVATORY     = EarthLocation.from_geodetic(lon="21:26:35.736", lat="-30:42:44.838",
                                                    height=1059.662443 * u.m)  # MeerKAT

# Channel-fit workers leak a few MB per prefix (casacore caches), so the pool is
# restarted every N prefixes: N = spare memory per worker / growth rate, within [MIN, MAX]
ASSUMED_WORKER_GROWTH_MB_PER_PREFIX = 5.0
MIN_POOL_RECYCLE_INTERVAL = 10
MAX_POOL_RECYCLE_INTERVAL = 200


# =============================================================================
# UTILITIES
# =============================================================================

def msg(txt):
    """Print a timestamped log line."""
    print(time.strftime(' %Y-%m-%d %H:%M:%S | ') + txt, flush=True)


_MJD_EPOCH = datetime.datetime(1858, 11, 17)


def _iso_to_mjd(date_str):
    """ISO date string (YYYY-MM-DDTHH:MM:SS[.f]) to MJD."""
    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S'):
        try:
            return (datetime.datetime.strptime(date_str, fmt) - _MJD_EPOCH).total_seconds() / 86400.0
        except ValueError:
            continue
    raise ValueError(f'Cannot parse date string for MJD conversion: {date_str!r}')


def _mjd_to_isot(mjd):
    """MJD to an ISO date string."""
    return (_MJD_EPOCH + datetime.timedelta(days=mjd)).isoformat()


def _apply_label(path, src_name, label):
    """
    Replace src_name with label in the file name of path only, never in its
    directories. Move a path into RESULTS/ before calling this, since that
    move matches on the original directory name.
    """
    if label == src_name:
        return path
    dirname, basename = os.path.split(path)
    return os.path.join(dirname, basename.replace(src_name, label))


def case_insensitive(text):
    """Glob pattern matching text in any letter case."""
    return ''.join(f'[{c.lower()}{c.upper()}]' if c.isalpha() else glob.escape(c) for c in text)


# A declination is the sexagesimal value after the comma of an 'RA,Dec' pair,
# with or without a sign
_DEC_COLON_RE = re.compile(r'(?<=,)(\s*[+-]?\d{1,3}):(\d{1,2}):(\d{1,2}(?:\.\d+)?)')


def fix_dec_colons(pos_str, context_label=''):
    """
    Convert colon-separated declinations ('-63:42:45.6', '12:30:00') to the
    period form CASA needs ('-63.42.45.6'), in a position ('RA,Dec') or every
    coordinate pair of a region string. RA is left alone. Non-strings (e.g.
    an unset rms_region) are returned unchanged.
    """
    if not isinstance(pos_str, str) or not pos_str:
        return pos_str
    fixed = _DEC_COLON_RE.sub(lambda m: f'{m.group(1)}.{m.group(2)}.{m.group(3)}', pos_str)
    if fixed != pos_str:
        msg(f'NOTE: {context_label} declination uses colons -- converting to CASA period format: '
            f'"{pos_str}" -> "{fixed}"')
    return fixed


def parse_casa_position(position_str):
    """
    Position string to (ra_deg, dec_deg), or (None, None) if it cannot be
    parsed. Accepts 'HH:MM:SS.SS,+DD.MM.SS.SS' or decimal degrees
    ('260.43,-16.20', optionally with 'deg').
    """
    try:
        if ':' not in position_str:
            ra, dec = position_str.replace('deg', '').split(',')
            return float(ra), float(dec)

        ra_str, dec_str = position_str.split(',')
        ra_parts = ra_str.strip().split(':')
        if len(ra_parts) != 3:
            msg(f'WARNING: Expected 3 parts in RA "{ra_str}", got {len(ra_parts)}')
            return None, None
        ra_deg = 15.0 * (float(ra_parts[0]) + float(ra_parts[1]) / 60.0 + float(ra_parts[2]) / 3600.0)

        dec_str  = dec_str.strip()
        dec_sign = -1.0 if dec_str.startswith('-') else 1.0
        dec_parts = dec_str.lstrip('+-').split('.')
        if len(dec_parts) < 3:
            msg(f'WARNING: Expected at least 3 parts in Dec "{dec_str}", got {len(dec_parts)}')
            return None, None
        dec_arcsec = float(dec_parts[2] + ('.' + dec_parts[3] if len(dec_parts) > 3 else ''))
        dec_deg = dec_sign * (float(dec_parts[0]) + float(dec_parts[1]) / 60.0 + dec_arcsec / 3600.0)
        return ra_deg, dec_deg

    except Exception as e:
        msg(f'ERROR parsing CASA position "{position_str}": {e}')
        return None, None


# =============================================================================
# WORKER POOL (channel fits): memory, worker count, chunking
# =============================================================================

def _proc_rss_mb(pid):
    """Resident memory (MB) of a process, read from /proc."""
    try:
        with open(f'/proc/{pid}/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    return 0.0


def log_memory_usage(label=''):
    """Log the memory of this process and the pool workers, and the change
    since the last call (steady growth points to a casacore cache leak)."""
    try:
        import psutil
        proc = psutil.Process()
        self_mb = proc.memory_info().rss / 1024 ** 2
        children = [c for c in proc.children(recursive=True) if c.is_running()]
        worker_mbs = [c.memory_info().rss / 1024 ** 2 for c in children]
    except ImportError:
        self_mb = _proc_rss_mb(os.getpid())
        worker_pids = [p.pid for p in mp.active_children() if p.pid is not None]
        worker_mbs = [_proc_rss_mb(pid) for pid in worker_pids]

    n_workers  = len(worker_mbs)
    workers_mb = sum(worker_mbs)
    avg_mb     = workers_mb / n_workers if n_workers else 0.0
    min_mb     = min(worker_mbs) if worker_mbs else 0.0
    max_mb     = max(worker_mbs) if worker_mbs else 0.0
    total_mb   = self_mb + workers_mb

    prev_total = getattr(log_memory_usage, '_prev_total_mb', None)
    delta_str  = '' if prev_total is None else f'  delta_vs_prev={total_mb - prev_total:+.0f}MB'
    log_memory_usage._prev_total_mb = total_mb

    msg(f'  [mem] {label} main={self_mb:.0f}MB  '
        f'workers({n_workers}): sum={workers_mb:.0f}MB avg={avg_mb:.0f}MB '
        f'min={min_mb:.0f}MB max={max_mb:.0f}MB  '
        f'total={total_mb:.0f}MB{delta_str}')


def get_available_memory_bytes():
    """Memory available to this job: the SLURM allocation, else the node's free memory."""
    slurm_mem_node = os.environ.get('SLURM_MEM_PER_NODE')
    if slurm_mem_node:
        try:
            allocated = int(slurm_mem_node) * 1024 ** 2
            msg(f'SLURM memory detected: SLURM_MEM_PER_NODE={slurm_mem_node} MB -> {allocated / 1024**3:.1f} GB')
            return allocated
        except ValueError:
            pass

    slurm_mem_cpu = os.environ.get('SLURM_MEM_PER_CPU')
    cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    ntasks        = os.environ.get('SLURM_NTASKS')
    alloc_cpus    = (int(cpus_per_task) * int(ntasks) if cpus_per_task and ntasks
                     else int(cpus_per_task) if cpus_per_task
                     else int(ntasks) if ntasks
                     else None)
    if slurm_mem_cpu and alloc_cpus:
        try:
            allocated = int(slurm_mem_cpu) * alloc_cpus * 1024 ** 2
            msg(f'SLURM memory detected: SLURM_MEM_PER_CPU={slurm_mem_cpu} MB x {alloc_cpus} CPUs -> {allocated / 1024**3:.1f} GB')
            return allocated
        except ValueError:
            pass

    try:
        import psutil
        available = psutil.virtual_memory().available
        msg(f'No SLURM memory allocation detected -- using node available memory: {available / 1024**3:.1f} GB')
        return available
    except ImportError:
        pass

    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    available = int(line.split()[1]) * 1024
                    msg(f'No SLURM memory allocation detected -- using /proc/meminfo: {available / 1024**3:.1f} GB')
                    return available
    except OSError:
        pass

    msg('WARNING: Could not determine available memory. Assuming 4 GB.')
    return 4 * 1024 ** 3


def estimate_memory_per_worker_bytes(image, casa_worker_floor_mb=350):
    """Peak memory of one worker: 4x the image size, or the cost of a fresh
    casatasks import, whichever is larger."""
    with fits.open(image) as hdul:
        shape       = hdul[0].data.shape
        dtype       = hdul[0].data.dtype
        image_bytes = 4 * hdul[0].data.nbytes
    floor_bytes  = casa_worker_floor_mb * 1024 ** 2
    chosen_bytes = max(image_bytes, floor_bytes)
    msg(f'  [mem-est] sample={os.path.basename(image)}  shape={shape}  '
        f'dtype={dtype}  x4 copies = {image_bytes / 1024**2:.1f}MB  '
        f'vs  floor={floor_bytes / 1024**2:.0f}MB  '
        f'-> using {"image size" if image_bytes >= floor_bytes else "casatasks import floor"} '
        f'({chosen_bytes / 1024**2:.1f}MB/worker)')
    return chosen_bytes


def compute_max_workers(sample_images, n_jobs):
    """
    Returns (workers, recycle_interval). workers is the smallest of the
    memory cap (from sample_images), the allocated CPUs and n_jobs;
    recycle_interval is how many prefixes the pool runs before a restart.
    """
    available  = get_available_memory_bytes()
    per_worker = max(estimate_memory_per_worker_bytes(im) for im in sample_images)
    mem_cap    = max(1, int(available // per_worker))

    cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    ntasks        = os.environ.get('SLURM_NTASKS')
    job_cpus      = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    job_cpus_m    = re.match(r'(\d+)', job_cpus)

    if cpus_per_task and ntasks:
        cpu_cap  = int(cpus_per_task) * int(ntasks)
        cpu_src  = f'SLURM_CPUS_PER_TASK={cpus_per_task} x SLURM_NTASKS={ntasks}'
    elif cpus_per_task:
        cpu_cap  = int(cpus_per_task)
        cpu_src  = f'SLURM_CPUS_PER_TASK={cpus_per_task}'
    elif job_cpus_m:
        cpu_cap  = int(job_cpus_m.group(1))
        cpu_src  = f'SLURM_JOB_CPUS_PER_NODE={job_cpus} -> {cpu_cap}'
    elif ntasks:
        cpu_cap  = int(ntasks)
        cpu_src  = f'SLURM_NTASKS={ntasks}'
    else:
        cpu_cap  = os.cpu_count() or 1
        cpu_src  = 'os.cpu_count() (no SLURM allocation detected)'

    workers        = max(1, min(mem_cap, cpu_cap, n_jobs))
    projected_used = per_worker * workers
    msg(f'Memory available   : {available / 1024**3:.2f} GB')
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB ({per_worker / 1024**2:.0f} MB)')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, jobs={n_jobs})]')
    msg(f'Projected usage    : {workers} workers x {per_worker / 1024**2:.0f} MB '
        f'= {projected_used / 1024**3:.2f} GB of {available / 1024**3:.2f} GB available '
        f'({100.0 * projected_used / available:.0f}%)')

    headroom = max(0.0, available / workers - per_worker)
    raw_interval = int(headroom / (ASSUMED_WORKER_GROWTH_MB_PER_PREFIX * 1024 ** 2))
    recycle_interval = max(MIN_POOL_RECYCLE_INTERVAL, min(MAX_POOL_RECYCLE_INTERVAL, raw_interval))
    msg(f'Pool recycle interval : {recycle_interval} prefixes  '
        f'[{headroom / 1024**2:.0f} MB spare/worker / {ASSUMED_WORKER_GROWTH_MB_PER_PREFIX:.0f} MB/prefix '
        f'= {raw_interval}, clamped to [{MIN_POOL_RECYCLE_INTERVAL}, {MAX_POOL_RECYCLE_INTERVAL}]]')

    return workers, recycle_interval


def chunk_jobs(jobs, n_chunks):
    """Split jobs round-robin, so channels clustered in cost (e.g. by RFI) spread across workers."""
    n_chunks = max(1, min(n_chunks, len(jobs)))
    chunks = [[] for _ in range(n_chunks)]
    for i, job in enumerate(jobs):
        chunks[i % n_chunks].append(job)
    return [c for c in chunks if c]


# =============================================================================
# CASA FITTING
# =============================================================================

def return_max(im, region):
    """Peak or trough in region, whichever is larger in magnitude (Q, U, V can be negative)."""
    ims = imstat(im, region=region)
    fmax, fmin = ims['max'][0], ims['min'][0]
    return fmax if fmax > abs(fmin) else fmin


def get_imstat_values(image, xpix, ypix, manual_rms_region=False):
    """
    [peak, peak x, peak y, rms, rms region] near (xpix, ypix). The rms is
    measured in an annulus of ~500 beam areas around the source (or
    manual_rms_region), on the matching residual image when USE_RESIDUAL_RMS
    is set and one exists.
    """
    bmaj = imhead(image, mode='get', hdkey='BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey='BMIN')['value']

    r_in  = 5.0 * bmaj
    r_out = np.sqrt(500 * 0.25 * bmaj * bmin + r_in ** 2)
    src_region = f'circle[[{xpix}pix,{ypix}pix],3.0pix]'
    rms_region = manual_rms_region or f'annulus[[{xpix}pix,{ypix}pix],[{r_in}arcsec,{r_out}arcsec]]'

    ims  = imstat(image, region=src_region)
    flux = return_max(image, src_region)

    rms_image = image
    if USE_RESIDUAL_RMS:
        residual = re.sub(r'-image(\.homogenized)?\.fits$', r'-residual\1.fits', image)
        if os.path.exists(residual):
            rms_image = residual
        elif not getattr(get_imstat_values, 'warned', False):
            msg(f'  WARNING: no residual image for the RMS (e.g. {os.path.basename(residual)}); using the image itself')
            get_imstat_values.warned = True
    rms = imstat(rms_image, region=rms_region)['rms'][0]

    return [flux, ims['maxpos'][0], ims['maxpos'][1], rms, rms_region]


def make_estimate(fname, image, xpix, ypix, fix_var):
    """
    Write an imfit estimate file: one beam-shaped Gaussian per component,
    seeded with the local peak. fix_var holds one CASA fix string per
    component ('abp': position free, 'xyabp': position fixed).
    """
    bmaj = imhead(image, mode='get', hdkey='BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey='BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey='BPA')['value']
    with open(fname, 'w') as f:
        for x, y, fix in zip(xpix, ypix, fix_var):
            f.write(f'{get_imstat_values(image, x, y)[0]},{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix}\n')


def get_imfit_values(fname, image, xpix, ypix):
    """
    imfit result for image with estimates fname, in a circle centred on the
    first component: 5 beams across, or twice the largest component
    separation if that is bigger.
    """
    bmaj = imhead(image, mode='get', hdkey='BMAJ')['value']
    r = 5 * bmaj
    if len(xpix) > 1:
        pixel_asec = np.degrees(imhead(image, mode='get', hdkey='cdelt2')['value']) * 3600
        points = np.array([xpix, ypix]).T
        r = max(r, 2.0 * np.amax(cdist(points, points)) * pixel_asec)
    return imfit(image, estimates=fname, region=f'circle[[{xpix[0]}pix,{ypix[0]}pix],{r}arcsec]')


def check_position(fname, image, xpix, ypix, snr_thresh=None, P_image=False,
                   fix_additional_comps=False, manual_rms_region=False, force_fix_all=False):
    """
    Write the estimate file for image, choosing per component whether its
    position is fitted ('abp') or fixed ('xyabp'):
      force_fix_all                          fixed
      fix_additional_comps and component > 0 fixed (faint secondary components)
      S/N >= snr_thresh                      free
      otherwise                              fixed
    snr_thresh defaults to P_FREE_FIT_SNR for P images, FREE_FIT_SNR otherwise.
    For P images (Rice noise) the S/N uses the largest rms of the Stokes
    images P is built from, matching calculate_P0: Q and U for Plin, Q, U
    and V for Ptot.
    """
    bmaj = imhead(image, mode='get', hdkey='bmaj')['value']
    bmin = imhead(image, mode='get', hdkey='bmin')['value']
    bpa  = imhead(image, mode='get', hdkey='bpa')['value']
    name = os.path.basename(image)
    if snr_thresh is None:
        snr_thresh = P_FREE_FIT_SNR if P_image else FREE_FIT_SNR
    noise_stokes = 'QUV' if '-Ptot-' in name else 'QU'

    fix_var = []
    for k, (x, y) in enumerate(zip(xpix, ypix)):
        # Flux at the trial position: a zero-amplitude, fixed-position fit
        check_file = f'check_pos_{os.path.basename(fname)}_{k}.txt'
        with open(check_file, 'w') as f:
            f.write(f'0.0,{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, xyabp')
        test_flux = abs(imfit(image, region=f'circle[[{x}pix,{y}pix],{2 * bmaj}arcsec]',
                              estimates=check_file)['results']['component0']['peak']['value'])

        if P_image:
            rms = max(get_imstat_values(image.replace('-Plin-', f'-{s}-').replace('-Ptot-', f'-{s}-'),
                                        x, y, manual_rms_region)[3] for s in noise_stokes)
        else:
            rms = get_imstat_values(image, x, y, manual_rms_region)[3]
        snr = test_flux / rms

        if force_fix_all:
            fix, why = 'xyabp', 'FIXED (force_fix_all)'
        elif fix_additional_comps and k > 0:
            fix, why = 'xyabp', 'FIXED (secondary component)'
        elif snr >= snr_thresh:
            fix, why = 'abp', f'FREE (S/N >= {snr_thresh})'
        else:
            fix, why = 'xyabp', f'FIXED (S/N < {snr_thresh})'
        fix_var.append(fix)
        msg(f'  Component {k}: {why} | S/N = {snr:.1f} | {name}')

    make_estimate(fname, image, xpix, ypix, fix_var)


def check_fit_offset_and_refix(imfit_result, image, ref_xpix, ref_ypix, bmaj,
                               src_name, stokes_label, tag=''):
    """
    If any fitted component lies more than 1/3 beam from its reference
    position (normally Stokes I), refit the whole image with every position
    fixed to the reference. Returns (imfit_result, refitted).
    """
    cell_arcsec = np.degrees(abs(imhead(image, mode='get', hdkey='CDELT2')['value'])) * 3600.0
    threshold_arcsec = bmaj / 3.0

    exceeded = False
    for z, component in enumerate(comp_keys(imfit_result)):
        x, y = imfit_result['results'][component]['pixelcoords'][:2]
        offset_pix = np.hypot(x - ref_xpix[z], y - ref_ypix[z])
        offset_arcsec = offset_pix * cell_arcsec
        if offset_arcsec > threshold_arcsec:
            exceeded = True
            msg(f'  WARNING: Stokes {stokes_label} component {z} offset {offset_arcsec:.2f}" '
                f'({offset_pix:.1f} pix) > 1/3 BMAJ ({threshold_arcsec:.2f}") -- refitting with fixed positions')
        else:
            msg(f'  Position check OK: Stokes {stokes_label} component {z} offset {offset_arcsec:.2f}" '
                f'({offset_pix:.1f} pix)')

    if exceeded:
        fixed_est = f'estimate_fixed_{stokes_label}_{src_name}{tag}.txt'
        make_estimate(fixed_est, image, ref_xpix, ref_ypix, ['xyabp'] * len(ref_xpix))
        imfit_result = get_imfit_values(fixed_est, image, ref_xpix, ref_ypix)
    return imfit_result, exceeded


def fit_stokes(est_fname, image, seed_xpix, seed_ypix, ref_xpix, ref_ypix, bmaj,
               src_name, stokes_label, tag='', **check_kwargs):
    """
    Fit one Stokes image: choose free or fixed positions at the seed
    (check_position), fit around the reference position, and refit with
    positions fixed to the reference if any component drifts > 1/3 beam.
    """
    check_position(est_fname, image, seed_xpix, seed_ypix, **check_kwargs)
    imfit_result = get_imfit_values(est_fname, image, ref_xpix, ref_ypix)
    return check_fit_offset_and_refix(imfit_result, image, ref_xpix, ref_ypix, bmaj,
                                      src_name, stokes_label, tag=tag)[0]


def comp_keys(imfit_result):
    """Component keys ('component0', ...) of an imfit result."""
    return [key for key in imfit_result['results'] if 'component' in key]


def fitted_pix(imfit_result, components):
    """Fitted pixel positions of the components, as (x list, y list)."""
    return ([imfit_result['results'][c]['pixelcoords'][0] for c in components],
            [imfit_result['results'][c]['pixelcoords'][1] for c in components])


def peak(imfit_result, component):
    """Fitted peak and its error, in mJy."""
    p = imfit_result['results'][component]['peak']
    return p['value'] * 1e3, p['error'] * 1e3


def sky_deg(imfit_result, component):
    """Fitted RA and Dec, in degrees."""
    d = imfit_result['results'][component]['shape']['direction']
    return np.degrees(d['m0']['value']), np.degrees(d['m1']['value'])


def cleanup_estimate_files(src_name):
    """Remove the estimate and check-position scratch files written for src_name."""
    removed = 0
    for pattern in (f'estimate_*{src_name}*.txt', f'check_pos_*{src_name}*.txt'):
        for fpath in glob.glob(pattern):
            try:
                os.remove(fpath)
                removed += 1
            except OSError:
                pass
    if removed:
        msg(f'Cleaned up {removed} temporary estimate/check-position file(s) for {src_name}.')


# =============================================================================
# POLARIZATION QUANTITIES
# =============================================================================

def calculate_P0(flux_P, rms_Q, rms_U, rms_V, pol_flag, Aq=0.8):
    """
    De-biased linear polarization and its noise, (P0, rms_P).

    With a polarization angle calibrator: noise from Hales et al. (2012),
    de-biased when S/N >= 3 (Vaillancourt 2006). Without one, P includes V,
    so the noise is the largest of rms_Q, rms_U, rms_V and P is always
    de-biased with a factor of 2 (from Monte Carlo tests).
    """
    if pol_flag:
        A, B = (Aq, 1.0 - Aq) if rms_Q >= rms_U else (1.0 - Aq, Aq)
        rms_P = (A * rms_Q ** 2 + B * rms_U ** 2) ** 0.5
        flux_P0 = (flux_P ** 2 - rms_P ** 2) ** 0.5 if flux_P / rms_P >= 3 else flux_P
    else:
        rms_P = np.amax([rms_Q, rms_U, rms_V])
        flux_P0 = (flux_P ** 2 - 2.0 * rms_P ** 2) ** 0.5
    return flux_P0, rms_P


def component_results(imfits, images, component, manual_rms_region, pol_flag):
    """
    Results for one component from the imfit result and image of each Stokes
    ({'I': ...} or {'I', 'P', 'Q', 'U', 'V'}): fluxes and fit errors (mJy),
    rms (mJy, all measured at the component's Stokes I position), positions
    (deg) and, with polarization, P0, LP fraction (%) and EVPA (deg; only
    with a polarization angle calibrator).
    """
    x, y = imfits['I']['results'][component]['pixelcoords'][:2]
    rms = {s: get_imstat_values(images[s], x, y, manual_rms_region)[3] * 1e3
           for s in images if s != 'P'}

    out = {}
    for s in imfits:
        out[f'{s}_flux_mJy'], out[f'{s}_err_mJy'] = peak(imfits[s], component)
        if s != 'P':
            out[f'{s}_rms_mJy'] = rms[s]
    out['I_RA_deg'], out['I_DEC_deg'] = sky_deg(imfits['I'], component)
    if 'P' not in imfits:
        return out

    out['P_RA_deg'], out['P_DEC_deg'] = sky_deg(imfits['P'], component)
    flux_I, flux_Q, flux_U = out['I_flux_mJy'], out['Q_flux_mJy'], out['U_flux_mJy']
    flux_P0, rms_P = calculate_P0(out['P_flux_mJy'], rms['Q'], rms['U'], rms['V'], pol_flag)
    out['P0_flux_mJy'], out['P_rms_mJy'] = flux_P0, rms_P
    out['LP_frac']     = flux_P0 / flux_I * 100.0
    out['LP_frac_err'] = abs(out['LP_frac']) * np.sqrt((rms['I'] / flux_I) ** 2 + (rms_P / flux_P0) ** 2)
    if pol_flag:
        out['LP_EVPA']     = 0.5 * np.degrees(np.arctan2(flux_U, flux_Q))
        out['LP_EVPA_err'] = 0.5 * np.degrees(np.sqrt(flux_U ** 2 * rms['Q'] ** 2 + flux_Q ** 2 * rms['U'] ** 2)
                                              / (flux_U ** 2 + flux_Q ** 2))
    else:
        out['LP_EVPA'] = out['LP_EVPA_err'] = None
    return out


# =============================================================================
# U/V FIX (no polarization angle calibrator): merge U and V into one spectrum
# =============================================================================

def fit_linear_xy_phase(freq_hz, U, V, rms_U, rms_V, max_delay_ns, oversample=8):
    '''
    Fit a linear (delay + offset) cross-hand phase to uncalibrated U and V.

    Without a polarization angle calibrator the cross-hand phase is never
    solved, so if Stokes V is intrinsically zero the observed U + iV is the
    true U rotated by an unknown phase phi(nu):

        U + iV = U_true * exp(i * phi(nu)),  phi = phi0 + 2*pi*tau*(nu - nu_ref)

    Squaring removes the sign of U_true (including its Faraday rotation sign
    flips), leaving (U + iV)^2 = U_true^2 * exp(2i * phi): a non-negative
    amplitude times a phase that is purely instrumental. tau is found by a
    coherent delay search over that quantity, so the whole band is used at once
    and the fit works when individual channels are noise-dominated. Channels
    are weighted by S/N^2, and the noise bias of the square, rms_U^2 - rms_V^2,
    is subtracted.

    The sign of U_true cannot be recovered (phi0 -> phi0 + pi gives the same
    data), so phi0 is returned in (-pi/2, pi/2].

    Returns a dict with tau_s, phi0_rad, nu_ref_hz, snr (coherent S/N of the
    fit), n_chan (channels used) and at_edge (best delay within 5% of the search
    limit), or None if fewer than 8 channels are usable.
    '''

    freq  = np.asarray(freq_hz, dtype=float)
    U     = np.asarray(U, dtype=float)
    V     = np.asarray(V, dtype=float)
    rms_U = np.asarray(rms_U, dtype=float)
    rms_V = np.asarray(rms_V, dtype=float)

    good = (np.isfinite(freq) & np.isfinite(U) & np.isfinite(V) &
            np.isfinite(rms_U) & np.isfinite(rms_V) & (rms_U > 0) & (rms_V > 0))
    n_chan = int(good.sum())
    if n_chan < 8:
        return None

    nu_ref = np.mean(freq[good])
    dnu    = freq[good] - nu_ref

    c2 = (U[good] + 1j * V[good]) ** 2 - (rms_U[good] ** 2 - rms_V[good] ** 2)
    z  = c2 / (0.5 * (rms_U[good] ** 2 + rms_V[good] ** 2))

    # Phase of the squared quantity advances by 4*pi*tau*dnu, so one cycle
    # across the band is tau = 1/(2*bandwidth).
    bandwidth = dnu.max() - dnu.min()
    step    = 1.0 / (2.0 * bandwidth * oversample)
    tau_max = max_delay_ns * 1e-9
    tau     = np.arange(-tau_max, tau_max + step, step)
    Z       = np.exp(-4j * np.pi * np.outer(tau, dnu)) @ z

    # Refine the best grid point on a finer grid spanning its neighbours
    tau_fine = np.linspace(tau[np.argmax(np.abs(Z))] - step, tau[np.argmax(np.abs(Z))] + step, 41)
    Z_fine   = np.exp(-4j * np.pi * np.outer(tau_fine, dnu)) @ z
    best     = np.argmax(np.abs(Z_fine))

    # Noise floor of |Z| for the coherent S/N. The median |Z| across the delay
    # grid (Rayleigh median = sigma*sqrt(2*ln2)) measures it directly, so the
    # S/N does not depend on the quoted rms being right. That needs enough
    # independent delays (resolution is 1/(2*bandwidth)); for a narrower band
    # the floor is the analytic one, where each of Re(Z), Im(Z) has variance
    # 4*n_chan if the quoted rms are right.
    if 4.0 * tau_max * bandwidth >= 16:
        floor = np.median(np.abs(Z)) / np.sqrt(2.0 * np.log(2.0))
    else:
        floor = np.sqrt(4.0 * n_chan)

    return {
        'tau_s'       : float(tau_fine[best]),
        'phi0_rad'    : float(0.5 * np.angle(Z_fine[best])),
        'nu_ref_hz'   : float(nu_ref),
        'snr'         : float(np.abs(Z_fine[best]) / floor),
        'n_chan'      : n_chan,
        'at_edge'     : bool(abs(tau_fine[best]) >= 0.95 * tau_max),
        # For the diagnostic plot: the delay search, its noise floor, and the
        # squared phasor the fit was made to
        'floor'       : float(floor),
        'tau_grid_s'  : tau,
        'Z_abs'       : np.abs(Z),
        'freq_used_hz': freq[good],
        'z'           : z,
    }



def derotate_uv(freq_hz, U, V, rms_U, rms_V, fit):
    '''
    Rotate (U, V) by minus the fitted linear cross-hand phase. Returns the
    derotated U and residual V with their propagated errors. The rotation is
    linear, so noise stays zero-mean Gaussian (no Rice bias).
    '''

    phi = fit['phi0_rad'] + 2.0 * np.pi * fit['tau_s'] * (np.asarray(freq_hz, dtype=float) - fit['nu_ref_hz'])
    cos_phi, sin_phi = np.cos(phi), np.sin(phi)

    U_alt = U * cos_phi + V * sin_phi
    V_res = V * cos_phi - U * sin_phi
    rms_U_alt = np.sqrt((cos_phi * rms_U) ** 2 + (sin_phi * rms_V) ** 2)
    rms_V_res = np.sqrt((sin_phi * rms_U) ** 2 + (cos_phi * rms_V) ** 2)

    return U_alt, V_res, rms_U_alt, rms_V_res



def parallactic_angle_deg(ra_deg, dec_deg, mjd_days, observatory=UVFIX_OBSERVATORY):
    '''
    Parallactic angle in degrees, wrapped to (-180, 180], of (ra_deg, dec_deg) (ICRS) at mjd_days (UTC):
    the position angle, east of north, between the field direction and the local zenith, from astropy.
    '''

    field  = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    zenith = SkyCoord(AltAz(alt=90 * u.deg, az=0 * u.deg, obstime=Time(mjd_days, format='mjd', scale='utc'),
                            location=observatory)).transform_to('icrs')
    chi    = field.position_angle(zenith).wrap_at(180 * u.deg)

    return float(chi.degree)



def rotate_to_sky(Q, U, rms_Q, rms_U, chi_deg):
    '''Rotate (Q, U) from the feed frame to the sky frame by -2*chi, propagating the errors.'''

    rot = np.radians(-2.0 * chi_deg)
    cos_r, sin_r = np.cos(rot), np.sin(rot)

    Q_sky = Q * cos_r + U * sin_r
    U_sky = -Q * sin_r + U * cos_r
    rms_Q_sky = np.sqrt((cos_r * rms_Q) ** 2 + (sin_r * rms_U) ** 2)
    rms_U_sky = np.sqrt((sin_r * rms_Q) ** 2 + (cos_r * rms_U) ** 2)

    return Q_sky, U_sky, rms_Q_sky, rms_U_sky



def mad_rms_outlier_channels(rms_arrays, nsigma=5.0):
    '''
    True for channels where ANY of the noise arrays is a MAD outlier (more than
    nsigma robust sigmas from its median, or non-finite). This is the test
    plot_stokes_spectrum uses to mask pathological channels from I, Q, U, V and P
    alike.
    '''

    bad = np.zeros(len(rms_arrays[0]), dtype=bool)
    for rms in rms_arrays:
        rms = np.asarray(rms, dtype=float)
        med = np.nanmedian(rms)
        sig = 1.4826 * np.nanmedian(np.abs(rms - med))
        if sig > 0:
            bad |= ~(np.abs(rms - med) <= nsigma * sig)

    return bad



def plot_uvfix_diagnostics(out_png, freq, U, V, rms_U, rms_V, U_alt, V_res, rms_U_alt, rms_V_res, fit, summary):
    '''
    Diagnostics of the linear cross-hand phase fit, for fits that pass and for
    fits that do not (when nothing else is written):
        1. delay search: |Z(tau)| over the search range, with the best delay,
           the noise floor and the S/N threshold
        2. phase of (U + iV)^2 against frequency, with the fitted 2*phi(nu)
           (dots are more opaque the higher their S/N^2)
        3. observed U and V
        4. U and V after derotating by the fit; V should be consistent with
           zero if the assumption holds
    '''

    from matplotlib.colors import to_rgba

    ok = summary['status'] == 'ok'
    fghz = freq / 1e9
    fused = fit['freq_used_hz'] / 1e9

    try:
        fig = plt.figure(figsize=(9, 14), constrained_layout=True)
        ax_per, ax_ph, ax_obs, ax_fix = fig.subplots(4, 1)

        # 1. delay search
        ax_per.plot(fit['tau_grid_s'] * 1e9, fit['Z_abs'], lw=0.8, color='C0')
        ax_per.axvline(fit['tau_s'] * 1e9, color='C3', ls='--', label=f'best tau = {fit["tau_s"] * 1e9:.3f} ns')
        ax_per.axhline(fit['floor'], color='0.5', ls=':', label='noise floor')
        ax_per.axhline(UVFIX_MIN_SNR * fit['floor'], color='C1', ls='-.', label=f'threshold ({UVFIX_MIN_SNR:g} x floor)')
        ax_per.set_xlim(-UVFIX_MAX_DELAY_NS, UVFIX_MAX_DELAY_NS)
        ax_per.set_xlabel('Cross-hand delay tau (ns)')
        ax_per.set_ylabel('|Z(tau)|')
        ax_per.set_title('Delay search: coherent S/N = %.1f%s%s' % (
            fit['snr'],
            '' if ok else ' -- below threshold, derotation unreliable',
            ' -- best delay at the edge of the search range' if fit['at_edge'] else ''), fontsize=9)
        ax_per.legend(loc='best', fontsize=8)

        # 2. phase of the squared phasor against the fit
        weight = np.abs(fit['z'])
        alpha  = np.clip(weight / max(np.percentile(weight, 95), 1e-12), 0.05, 1.0)
        rgba   = np.tile(to_rgba('C0'), (len(weight), 1))
        rgba[:, 3] = alpha
        psi_model = 2.0 * fit['phi0_rad'] + 4.0 * np.pi * fit['tau_s'] * (fit['freq_used_hz'] - fit['nu_ref_hz'])
        ax_ph.scatter(fused, np.degrees(np.angle(fit['z'])), c=rgba, s=10, label='phase of (U+iV)^2')
        ax_ph.scatter(fused, np.degrees(np.angle(np.exp(1j * psi_model))), c='C3', s=2, label='fitted 2*phi(nu)')
        ax_ph.set_ylim(-190, 190)
        ax_ph.set_yticks(range(-180, 181, 90))
        ax_ph.set_ylabel('Phase (deg)')
        ax_ph.set_title('Cross-hand phase: %.1f turn(s) across the band' % summary['xy_phase_turns_across_band'], fontsize=9)
        ax_ph.legend(loc='upper right', fontsize=8)

        # 3 and 4. observed and derotated spectra
        ax_obs.errorbar(fghz, U * 1e3, rms_U * 1e3, fmt='.', ms=3, lw=0.5, label='U (observed)')
        ax_obs.errorbar(fghz, V * 1e3, rms_V * 1e3, fmt='.', ms=3, lw=0.5, label='V (observed)')
        ax_fix.errorbar(fghz, U_alt * 1e3, rms_U_alt * 1e3, fmt='.', ms=3, lw=0.5, label='U (derotated)')
        ax_fix.errorbar(fghz, V_res * 1e3, rms_V_res * 1e3, fmt='.', ms=3, lw=0.5, label='V (residual)')
        ax_fix.set_title('V residual chi2_red = %.2f%s' % (
            summary['V_residual_chi2_red'], '' if ok else ' (fit below threshold)'), fontsize=9)
        for ax in (ax_ph, ax_obs, ax_fix):
            ax.set_xlim(fghz.min(), fghz.max())
        for ax in (ax_obs, ax_fix):
            ax.axhline(0, color='k', lw=0.5)
            ax.set_ylabel('Flux density (mJy)')
            ax.legend(loc='best', fontsize=8)
        ax_ph.set_xlabel('Frequency (GHz)')
        ax_fix.set_xlabel('Frequency (GHz)')

        fig.suptitle('%s: %d channel(s) used, %d MAD-clipped' % (
            os.path.basename(out_png), fit['n_chan'], summary['n_mad_clipped']), fontsize=9)
        fig.savefig(out_png, dpi=100)
        plt.close(fig)
    except Exception as e:
        msg(f'  WARNING: could not write U/V fix diagnostics plot ({e})')



def write_uvfix_rmsynth(rmsynth_fname, rmsynth_arr, V, rms_V, ra_deg, dec_deg, mjd_days):
    '''
    No polarization angle calibrator: write a companion to an RM synthesis
    input file in which all of the linear polarization is assumed to be U,
    spread into V by a linear cross-hand phase (see fit_linear_xy_phase).

    rmsynth_arr is the (freq, I, Q, U, dI, dQ, dU) array in Hz and Jy;
    V and rms_V are the matching Stokes V spectrum and noise in Jy.

    Channels whose I, Q, U or V noise is a MAD outlier (mad_rms_outlier_channels)
    are excluded from the fit and dropped from the companion file, as are
    channels with non-finite U or V.

    Outputs, alongside rmsynth_fname (<base> = rmsynth_fname minus '.txt'):
        <base>_UVfix.txt       freq, I, Q, U, dI, dQ, dU, cross-hand phase and
                               parallactic angle corrected; written
                               whenever a fit could be made, with a warning if
                               its coherent S/N is below UVFIX_MIN_SNR
        <base>_UVfix_fit.json  fit parameters, significance and V residual
                               statistics (always written)
        <base>_UVfix.png       diagnostics of the fit (plot_uvfix_diagnostics),
                               written whenever the fit could be made
    '''

    base     = rmsynth_fname[:-len('.txt')] + '_UVfix'
    out_txt  = base + '.txt'
    out_json = base + '_fit.json'
    out_png  = base + '.png'

    rmsynth_arr = np.asarray(rmsynth_arr, dtype=float)
    V           = np.asarray(V, dtype=float)
    rms_V       = np.asarray(rms_V, dtype=float)

    mad_clipped = mad_rms_outlier_channels([rmsynth_arr[4], rmsynth_arr[5], rmsynth_arr[6], rms_V])
    use = ~mad_clipped & np.isfinite(rmsynth_arr[3]) & np.isfinite(V)

    freq, I, Q, U, rms_I, rms_Q, rms_U = rmsynth_arr[:, use]
    V, rms_V = V[use], rms_V[use]

    fit = fit_linear_xy_phase(freq, U, V, rms_U, rms_V, UVFIX_MAX_DELAY_NS)

    summary = {
        'max_delay_ns'   : UVFIX_MAX_DELAY_NS,
        'min_snr'        : UVFIX_MIN_SNR,
        'n_chan_total'   : int(len(use)),
        'n_mad_clipped'  : int(mad_clipped.sum()),
    }
    if mad_clipped.any():
        msg(f'  U/V fix: excluding {int(mad_clipped.sum())}/{len(use)} MAD-clipped channels (outlier RMS)')

    if fit is None:
        summary['status'] = 'too few usable channels'
    else:
        summary.update({
            'tau_ns'      : fit['tau_s'] * 1e9,
            'phi0_deg'    : float(np.degrees(fit['phi0_rad'])),
            'nu_ref_GHz'  : fit['nu_ref_hz'] / 1e9,
            'xy_phase_turns_across_band': float(abs(fit['tau_s']) * (freq.max() - freq.min())),
            'coherent_snr': fit['snr'],
            'noise_floor' : fit['floor'],
            'n_chan_used' : fit['n_chan'],
            'at_search_edge': fit['at_edge'],
            'status'      : 'ok' if fit['snr'] >= UVFIX_MIN_SNR else 'coherent S/N below threshold',
        })

        # Derotated even when the fit is below threshold, for the diagnostics
        U_alt, V_res, rms_U_alt, rms_V_res = derotate_uv(freq, U, V, rms_U, rms_V, fit)
        usable = np.isfinite(V_res) & np.isfinite(rms_V_res) & (rms_V_res > 0)
        summary['V_residual_chi2_red'] = float(np.mean((V_res[usable] / rms_V_res[usable]) ** 2))

        # Feed frame -> sky frame, after the cross-hand phase is removed
        chi_deg = parallactic_angle_deg(ra_deg, dec_deg, mjd_days)
        Q_sky, U_sky, rms_Q_sky, rms_U_sky = rotate_to_sky(Q, U_alt, rms_Q, rms_U_alt, chi_deg)
        summary['parang_deg'] = chi_deg
        summary['parang_mjd'] = float(mjd_days)
        summary['parang_ra_deg'] = float(ra_deg)
        summary['parang_dec_deg'] = float(dec_deg)
        msg(f'  U/V fix: parallactic angle chi = {chi_deg:+.3f} deg at MJD {mjd_days:.5f}; (Q, U) rotated by {-2.0 * chi_deg:+.3f} deg')

        msg(f'  U/V fix: tau = {summary["tau_ns"]:.3f} ns ({summary["xy_phase_turns_across_band"]:.1f} turn(s) across the band), '
            f'phi0 = {summary["phi0_deg"]:.1f} deg, coherent S/N = {fit["snr"]:.1f} (threshold {UVFIX_MIN_SNR:g}), '
            f'V residual chi2_red = {summary["V_residual_chi2_red"]:.2f}')
        if fit['at_edge']:
            msg(f'    WARNING: fitted delay is at the edge of the +/-{UVFIX_MAX_DELAY_NS:.1f} ns search range')

    if fit is not None:
        np.savetxt(out_txt, np.array([freq, I, Q_sky, U_sky, rms_I, rms_Q_sky, rms_U_sky]).T)
        msg(f'  U/V-merged RM synthesis file written: {out_txt}')
        if summary['status'] != 'ok':
            msg(f'  WARNING: coherent S/N = {fit["snr"]:.1f} is below the threshold ({UVFIX_MIN_SNR:g}) -- '
                f'file written anyway, treat the derotation as unreliable')
    else:
        # Drop an earlier run's file so a stale one is never picked up
        if os.path.exists(out_txt):
            os.remove(out_txt)
        msg(f'  U/V-merged RM synthesis file NOT written: {summary["status"]}')

    if fit is not None:
        plot_uvfix_diagnostics(out_png, freq, U, V, rms_U, rms_V, U_alt, V_res, rms_U_alt, rms_V_res, fit, summary)
        msg(f'  U/V fix diagnostics written: {out_png}')
    elif os.path.exists(out_png):
        os.remove(out_png)

    with open(out_json, 'w') as j:
        json.dump(summary, j, indent=4)


# =============================================================================
# IMAGE DISCOVERY
# =============================================================================

def scan_stokes_images(src_im_identifier, mfs_im_suffix, chan_im_suffix, pol_flag):
    '''
    Locate the MFS and per-channel images for one source, matching every image
    to its Stokes parameter by name, and decide whether full Stokes (I, P, Q,
    U, V) can be fitted or only Stokes I.

    Multi-Stokes images are named '{prefix}-MFS-{stokes}-{suffix}' and
    '{prefix}-{chan}-{stokes}-{suffix}'; Stokes I only imaging has no Stokes
    tag, '{prefix}-MFS-{suffix}' and '{prefix}-{chan}-{suffix}'. P is Plin
    when pol_flag (a polarization angle calibrator was used), Ptot otherwise.

    Full Stokes is used if ANY prefix has the complete MFS set (I, P, Q, U, V)
    and at least one complete channel set; otherwise every prefix is fitted as
    Stokes I only. In a full Stokes run a prefix (time interval) without a
    complete MFS set is dropped, since its I and polarization results could not
    line up with the other epochs, and a channel without a complete set of
    Stokes images (e.g. a flagged channel) is left out. A prefix with a complete
    MFS set but no complete channel set is kept for its MFS results only.

    Returns a dict:
        prefixes     prefixes to fit, sorted: those with an MFS Stokes I image,
                     or in a full Stokes run those with a complete MFS set
        mfs          {prefix: {'I'|'P'|'Q'|'U'|'V': path}}
        chan_groups  {prefix: [[paths], ...]}, one group per channel in
                     frequency order: [I, P, Q, U, V] or [I]
        full_stokes  True if full Stokes will be fitted
        notes        why Stokes I only was chosen (empty if full Stokes)
        n_chan_incomplete  {prefix: channels left out for missing Stokes}
    '''

    p_tag      = 'Plin' if pol_flag else 'Ptot'
    tag_to_key = {'I': 'I', 'Q': 'Q', 'U': 'U', 'V': 'V', p_tag: 'P'}
    full_keys  = ['I', 'P', 'Q', 'U', 'V']
    tag_re     = '|'.join(['I', 'Q', 'U', 'V', 'Plin', 'Ptot'])

    mfs_re  = re.compile(rf'(?P<prefix>.+?)-MFS-(?:(?P<tag>{tag_re})-)?{re.escape(mfs_im_suffix)}')
    chan_re = re.compile(rf'(?P<prefix>.+)-(?P<chan>\d+)-(?:(?P<tag>{tag_re})-)?{re.escape(chan_im_suffix)}')

    def add(entry, tag, path):
        # A tagged Stokes I image is part of the multi-Stokes set and takes
        # precedence over an untagged one; the unused P type is ignored.
        if tag is None:
            entry.setdefault('I', path)
        elif tag in tag_to_key:
            entry[tag_to_key[tag]] = path

    def count_tags(counts, tag):
        label = 'I (untagged)' if tag is None else tag
        counts[label] = counts.get(label, 0) + 1

    def tag_summary(counts):
        return ', '.join(f'{k}={counts[k]}' for k in sorted(counts))

    mfs_pattern = f'{src_im_identifier}*MFS*{mfs_im_suffix}'
    msg(f'Locating MFS images: {mfs_pattern}')
    mfs_files = glob.glob(mfs_pattern)
    mfs, mfs_counts, mfs_unmatched = {}, {}, []
    for f in mfs_files:
        m = mfs_re.fullmatch(f)
        if m is None:
            mfs_unmatched.append(f)
            continue
        count_tags(mfs_counts, m.group('tag'))
        add(mfs.setdefault(m.group('prefix'), {}), m.group('tag'), f)
    msg(f'  {len(mfs_files)} file(s) matched; by Stokes: {tag_summary(mfs_counts) or "none"}')
    if mfs_unmatched:
        msg(f'  {len(mfs_unmatched)} file(s) not recognised as MFS Stokes images and ignored, e.g. '
            f'{os.path.basename(sorted(mfs_unmatched)[0])}')
    unused_p = 'Ptot' if pol_flag else 'Plin'
    if unused_p in mfs_counts:
        msg(f'  Ignoring {mfs_counts[unused_p]} {unused_p} MFS image(s): '
            f'{p_tag} is used {"with" if pol_flag else "without"} a polarization angle calibrator')

    # One glob for the channel images of every prefix, partitioned in memory.
    # Filtered by a '-MFS-' substring check rather than a [!MFS] glob bracket,
    # since fnmatch's '*' backtracking can defeat that bracket.
    chan_pattern = f'{src_im_identifier}*-{chan_im_suffix}'
    msg(f'Locating channel images: {chan_pattern}')
    chan, chan_counts, chan_unmatched, n_chan_files = {}, {}, [], 0
    for f in glob.glob(chan_pattern):
        if '-MFS-' in f:
            continue
        n_chan_files += 1
        m = chan_re.fullmatch(f)
        if m is None:
            chan_unmatched.append(f)
            continue
        count_tags(chan_counts, m.group('tag'))
        add(chan.setdefault(m.group('prefix'), {}).setdefault(int(m.group('chan')), {}), m.group('tag'), f)
    msg(f'  {n_chan_files} file(s) matched; by Stokes: {tag_summary(chan_counts) or "none"}')
    if chan_unmatched:
        msg(f'  {len(chan_unmatched)} file(s) not recognised as channel Stokes images and ignored, e.g. '
            f'{os.path.basename(sorted(chan_unmatched)[0])}')
    if unused_p in chan_counts:
        msg(f'  Ignoring {chan_counts[unused_p]} {unused_p} channel image(s)')

    prefixes = sorted(p for p in mfs if 'I' in mfs[p])
    for p in sorted(set(mfs) - set(prefixes)):
        msg(f'  WARNING: {os.path.basename(p)} has MFS images but no Stokes I image -- skipped')
    msg(f'  {len(prefixes)} prefix(es) with an MFS Stokes I image')

    groups_full, groups_I, n_incomplete = {}, {}, {}
    for p in prefixes:
        channels = sorted(chan.get(p, {}).items())
        groups_I[p]     = [[ch['I']] for _, ch in channels if 'I' in ch]
        groups_full[p]  = [[ch[k] for k in full_keys] for _, ch in channels if all(k in ch for k in full_keys)]
        n_incomplete[p] = len(groups_I[p]) - len(groups_full[p])

    no_mfs   = {p: [k for k in full_keys if k not in mfs[p]] for p in prefixes}
    complete = [p for p in prefixes if not no_mfs[p]]
    usable   = [p for p in complete if groups_full[p]]
    full_stokes = bool(usable)

    def missing_text(p):
        return (f'{os.path.basename(p)}: MFS image(s) missing for {",".join(no_mfs[p])}'
                + (f' (is {p_tag} made by make_pol_images.py?)' if 'P' in no_mfs[p] else ''))

    notes = []
    if not full_stokes:
        if not prefixes:
            notes.append('no MFS Stokes I image found')
        elif not any(set(mfs[p]) - {'I'} for p in prefixes):
            notes.append('no polarization images found (Stokes I only imaging)')
        else:
            notes.append(f'no prefix has a complete set of {",".join(full_keys)} images in both MFS and at least one channel')
            notes += [missing_text(p) for p in [p for p in prefixes if no_mfs[p]][:5]]
            notes += [f'{os.path.basename(p)}: MFS set complete but no channel has a complete set'
                      for p in [p for p in complete if not groups_full[p]][:5]]

    kept = prefixes
    if full_stokes:
        kept = complete
        dropped = [p for p in prefixes if no_mfs[p]]
        msg(f'  Full Stokes: {len(usable)} prefix(es) with a complete MFS set and at least one complete channel set')
        for p in dropped[:10]:
            msg(f'  WARNING: {missing_text(p)} -- prefix dropped (flagged interval?)')
        if dropped:
            msg(f'  {len(dropped)} of {len(prefixes)} prefix(es) dropped for missing polarization images')
        for p in kept:
            if not groups_full[p]:
                msg(f'  {os.path.basename(p)}: MFS set complete but no complete channel set -- MFS results only')
            elif n_incomplete[p]:
                msg(f'  {os.path.basename(p)}: {len(groups_full[p])} complete channel set(s), '
                    f'{n_incomplete[p]} channel(s) missing a Stokes image')

    return {
        'prefixes'         : kept,
        'mfs'              : {p: mfs[p] for p in kept},
        'chan_groups'      : {p: (groups_full if full_stokes else groups_I)[p] for p in kept},
        'full_stokes'      : full_stokes,
        'notes'            : notes,
        'n_chan_incomplete': {p: (n_incomplete[p] if full_stokes else 0) for p in kept},
    }


# =============================================================================
# PER-CHANNEL FITS (run in pool workers)
# =============================================================================

def fit_channel_group(job):
    """
    Fit one channel's images of one prefix (job['CHAN_images']: [I] or
    [I, P, Q, U, V]) in the same order as the MFS fit: I seeded at the MFS I
    position, P at the MFS P position, Q and U at the channel P position, and
    V at the channel I position.

    Returns (chan_idx, result). result holds freq_GHz, the beam and
    'components' (component_results per component); it is None when the
    channel's beam is over 3x the scaled MFS beam (heavily flagged) or a fit
    fails.
    """
    chan_idx, images, components = job['chan_idx'], job['CHAN_images'], job['components']
    src_name, tag, rms_region = job['src_name'], job['tag'], job['manual_rms_region']

    try:
        msg(f'Fitting image: {images[0]}')
        bmaj     = imhead(images[0], mode='get', hdkey='bmaj')['value']
        bmin     = imhead(images[0], mode='get', hdkey='bmin')['value']
        bpa      = imhead(images[0], mode='get', hdkey='bpa')['value']
        freq_GHz = imhead(images[0], mode='get', hdkey='CRVAL3')['value'] / 1.0e9

        if bmaj > 3.0 * job['mfs_freq_ghz'] / freq_GHz * job['mfs_bmaj_asec']:
            msg(f'Skipping channel {chan_idx}: beam over 3x the scaled MFS beam (likely heavily flagged)')
            return chan_idx, None

        common = dict(bmaj=bmaj, src_name=src_name, tag=tag, manual_rms_region=rms_region)
        I_x, I_y = job['MFS_I_ra_pix'], job['MFS_I_dec_pix']
        imfits = {'I': fit_stokes(f'estimate_I_{src_name}_{tag}.txt', images[0],
                                  I_x, I_y, I_x, I_y, stokes_label='I_chan', **common)}

        if not job['only_intensity']:
            pol = dict(fix_additional_comps=job['fix_additional_comps'], **common)
            cI_x, cI_y = fitted_pix(imfits['I'], components)
            imfits['P'] = fit_stokes(f'estimate_P_{src_name}_{tag}.txt', images[1],
                                     job['MFS_P_ra_pix'], job['MFS_P_dec_pix'], cI_x, cI_y,
                                     stokes_label='P', P_image=True, **pol)
            cP_x, cP_y = fitted_pix(imfits['P'], components)
            for s, image in (('Q', images[2]), ('U', images[3])):
                imfits[s] = fit_stokes(f'estimate_{s}_{src_name}_{tag}.txt', image,
                                       cP_x, cP_y, cP_x, cP_y, stokes_label=s, **pol)
            imfits['V'] = fit_stokes(f'estimate_V_{src_name}_{tag}.txt', images[4],
                                     cI_x, cI_y, cI_x, cI_y, stokes_label='V',
                                     force_fix_all=FORCE_FIX_STOKES_V, **pol)

        stokes_images = dict(zip('IPQUV', images))
        return chan_idx, {
            'freq_GHz': freq_GHz, 'bmaj_asec': bmaj, 'bmin_asec': bmin, 'bpa_deg': bpa,
            'components': {c: component_results(imfits, stokes_images, c, rms_region, job['pol_flag'])
                           for c in components}}

    except Exception:
        msg(f'Fitting failed: channel {chan_idx} is likely flagged')
        return chan_idx, None


def fit_channel_group_batch(jobs):
    """Pool entry point: fit a worker's chunk of channels."""
    return [fit_channel_group(job) for job in jobs]


# =============================================================================
# EXTRACTION
# =============================================================================

def ms_timing(src_name):
    """
    (centre MJD, span in s) of src_name's observation, from the TIME column of
    its split MS (project_info target_ms) or, failing that, the working MS
    filtered by its FIELD_ID. None when neither MS can be read.
    """
    if not os.path.isfile('project_info.json'):
        msg('  [timing] project_info.json not found -- will use image header fallback')
        return None
    with open('project_info.json') as f:
        pinfo = json.load(f)

    names = [name.lower() for name in pinfo.get('target_names', [])]
    idx = names.index(src_name.lower()) if src_name.lower() in names else None
    target_ms = pinfo.get('target_ms', [])
    target_ids = pinfo.get('target_ids', [])

    ms_path, field_id = None, None
    if idx is not None and idx < len(target_ms) and target_ms[idx] and os.path.exists(target_ms[idx]):
        ms_path = target_ms[idx]
        msg(f'  [timing] using split MS: {os.path.basename(ms_path)}')
    elif pinfo.get('working_ms') and os.path.exists(pinfo['working_ms']):
        ms_path = pinfo['working_ms']
        if idx is not None and idx < len(target_ids):
            field_id = int(target_ids[idx])
        msg(f'  [timing] split MS absent; using working MS: {os.path.basename(ms_path)}'
            + (f' (FIELD_ID={field_id})' if field_id is not None else ''))
    else:
        msg(f'  [timing] no MS found for "{src_name}" -- will use image header fallback')
        return None

    try:
        tb.open(ms_path)
        times = tb.getcol('TIME')
        if field_id is not None:
            times = times[tb.getcol('FIELD_ID') == field_id]
        tb.close()
    except Exception as e:
        msg(f'  [timing] failed to read MS ({e}) -- will use image header fallback')
        return None
    if len(times) == 0:
        msg('  [timing] no timestamps found -- will use image header fallback')
        return None

    t_start, t_end = float(times.min()), float(times.max())
    ctr_mjd = (t_start + t_end) / 2.0 / 86400.0
    msg(f'  [timing] MS: centre = {ctr_mjd:.8f} MJD, span = {t_end - t_start:.0f} s')
    return ctr_mjd, t_end - t_start


def record(store, values, k=None, n_epochs=None):
    """
    Append each of values to store[key]: one entry per epoch (MFS), or, given
    epoch k, into epoch k's list of channel values (store[key] holds n_epochs
    lists).
    """
    for key, val in values.items():
        if k is None:
            store.setdefault(key, []).append(val)
        else:
            store.setdefault(key, [[] for _ in range(n_epochs)])[k].append(val)


def write_spectra(chan, k, base, pol_flag, time_ctr_mjd):
    """
    Write epoch k's channel spectra of one component: <base>_rmsynth.txt
    (freq, I, Q, U, dI, dQ, dU), <base>_IQUV.txt (adds V and dV) and, without
    a polarization angle calibrator, the U/V-merged companion. Hz and Jy, all
    channels.
    """
    freq = np.array(chan['freq_GHz'][k]) * 1e9
    jy = {f'{s}_{q}': np.array(chan[f'{s}_{q}_mJy'][k]) / 1e3 for s in 'IQUV' for q in ('flux', 'rms')}

    rmsynth_fname = base + '_rmsynth.txt'
    rmsynth_arr = np.array([freq] + [jy[f'{s}_{q}'] for q in ('flux', 'rms') for s in 'IQU'])
    np.savetxt(rmsynth_fname, rmsynth_arr.T)
    msg(f'  RM synthesis file written: {rmsynth_fname}')

    np.savetxt(base + '_IQUV.txt', np.column_stack([freq] + [jy[f'{s}_{q}'] for q in ('flux', 'rms') for s in 'IQUV']),
               header='freq_Hz I_Jy Q_Jy U_Jy V_Jy dI_Jy dQ_Jy dU_Jy dV_Jy')
    msg(f'  IQUV spectrum file written: {base}_IQUV.txt')

    if not pol_flag:
        try:
            write_uvfix_rmsynth(rmsynth_fname, rmsynth_arr, jy['V_flux'], jy['V_rms'],
                                float(np.nanmedian(chan['I_RA_deg'][k])),
                                float(np.nanmedian(chan['I_DEC_deg'][k])), time_ctr_mjd)
        except Exception as e:
            msg(f'  WARNING: U/V fix failed for {os.path.basename(base)} ({e}) -- continuing without it')


def make_plots(output_dictionary, src_name, prefix_arr, label, pol_flag):
    """IQUV spectrum plots (above SPEC_PLOT_SNR_THRESH when time-resolved) and, when time-resolved, the light curve."""
    for component in sorted(output_dictionary['MFS']):
        if len(prefix_arr) > 1:
            _plot_detected_spectra(output_dictionary, src_name, prefix_arr, component, label, pol_flag=pol_flag)
        else:
            plot_stokes_spectrum(output_dictionary, src_name, prefix_arr[0], component, 0,
                                 save_plot=True, label=label, pol_flag=pol_flag)
    if len(prefix_arr) > 1:
        plot_light_curve(output_dictionary, src_name, prefix_arr, label=label)


def extract_polarization_properties(src_name, src_im_identifier, src_im_suffix, src_ra, src_dec,
                                    src_ulims, pol_flag, manual_rms_region, image_directory,
                                    image_identifier, fix_additional_comps=False, out_label=''):
    """
    Fit every component of one source on its MFS and channel images, for each
    prefix (time interval), then write the JSON, the spectra files and plots.

    src_name             name used to find the images and the MS (for timing)
    src_im_identifier    glob prefix of the source's images
    src_im_suffix        'image.fits' or 'image.homogenized.fits'; the other
                         is used for any image type with none of these
    src_ra, src_dec      input positions (deg), one per component
    src_ulims            per component: True fixes its MFS Stokes I position
    pol_flag             True with a polarization angle calibrator (Plin, EVPA);
                         False fits Ptot
    manual_rms_region    CASA region for the noise, or False for an annulus
    image_directory      image directory (replaced by RESULTS in output paths)
    image_identifier     image identifier, e.g. 'pcalmask'
    fix_additional_comps fix secondary component positions in P, Q, U, V
    out_label            name for output files and plot titles (default src_name)
    """
    msg('')
    msg(f'{"="*70}')
    msg(f'  Starting polarization extraction for source: {src_name}')
    msg(f'  Image identifier : {src_im_identifier}')
    msg(f'  RMS region       : {"MANUAL -- " + manual_rms_region if manual_rms_region else "default annulus (~500 beam areas) centred on source"}')
    msg(f'  RMS image        : {"matching residual where it exists" if USE_RESIDUAL_RMS else "the image itself"}')
    msg(f'  Polarization angle calibrator present: {pol_flag}')
    msg(f'  Fix secondary component positions    : {fix_additional_comps}')
    msg(f'{"="*70}')
    msg('')

    label = out_label or src_name
    if label != src_name:
        msg(f'  Output label     : "{label}" (overrides source name "{src_name}" in output filenames)')
    output_dictionary = {'name': label, 'MFS': {}, 'CHAN': {}}

    # Use the requested suffix where images with it exist; homogenized images
    # stand in for a type with none (e.g. chunked channelisation)
    other_suffix = 'image.homogenized.fits' if src_im_suffix == 'image.fits' else 'image.fits'

    def pick_suffix(mfs):
        for suffix in (src_im_suffix, other_suffix):
            if any(('-MFS-' in f) == mfs for f in glob.glob(f'{src_im_identifier}*{suffix}')):
                return suffix
        return None

    mfs_found, chan_found = pick_suffix(True), pick_suffix(False)

    # Confirm at least some images exist before proceeding
    if mfs_found is None and chan_found is None:
        raise FileNotFoundError(
            f'No images found matching: {src_im_identifier}*{src_im_suffix} or *{other_suffix}\n'
            f'Check the image_directory and image_identifier fields in rmsynth_info.json.')

    mfs_im_suffix = mfs_found or src_im_suffix
    chan_im_suffix = chan_found or src_im_suffix
    for kind, found in (('MFS', mfs_found), ('channel', chan_found)):
        if found is None:
            msg(f'WARNING: No {kind} images found with {src_im_suffix} or {other_suffix}')
        elif found != src_im_suffix:
            msg(f'WARNING: No {kind} {src_im_suffix} images found; using {found}')
    src_im_suffix = chan_im_suffix
    msg(f'MFS image suffix: {mfs_im_suffix}; channel image suffix: {src_im_suffix}')

    # One directory scan matches every MFS and channel image to its Stokes and
    # decides between full Stokes and Stokes I only. Each prefix is one time
    # interval (or the whole observation).
    stokes_scan    = scan_stokes_images(src_im_identifier, mfs_im_suffix, src_im_suffix, pol_flag)
    prefix_arr     = stokes_scan['prefixes']
    only_intensity = not stokes_scan['full_stokes']
    n_epochs       = len(prefix_arr)

    if only_intensity:
        msg('Fitting Stokes I only:')
        for note in stokes_scan['notes']:
            msg(f'  {note}')
    else:
        msg(f'Full Stokes images found -- fitting I, {"Plin" if pol_flag else "Ptot"}, Q, U, V')
    if not prefix_arr:
        msg('WARNING: No prefixes found -- nothing to process. '
            'Check image_directory, image_identifier, and image_suffix in rmsynth_info.json.')
        return 0
    n_chan = sum(len(v) for v in stokes_scan['chan_groups'].values())
    msg(f'Found {n_epochs} prefix(es) (time interval(s)) with {n_chan} channel(s) in total')
    n_chan_skipped = sum(stokes_scan['n_chan_incomplete'].values())
    if n_chan_skipped:
        msg(f'WARNING: {n_chan_skipped} channel(s) left out because a Stokes image is missing')

    # Exposure centre (MJD) and span (s) per prefix. Time-resolved: each
    # header DATE-OBS plus half the median step, spanning first to last prefix.
    # Single epoch: from the MS, else (in the loop) DATE-OBS with a 15 min span.
    if n_epochs > 1:
        hdr_mjd = np.array([_iso_to_mjd(imhead(stokes_scan['mfs'][p]['I'], mode='get', hdkey='DATE-OBS')
                                        .replace('/', '-', 2).replace('/', 'T')) for p in prefix_arr])
        step_s = float(np.median(np.diff(hdr_mjd))) * 86400.0
        span_s = float(hdr_mjd[-1] - hdr_mjd[0]) * 86400.0
        timing = [(mjd + step_s / 2.0 / 86400.0, span_s) for mjd in hdr_mjd]
        msg(f'  [timing] time-resolved ({n_epochs} prefixes): '
            f'median step = {step_s:.1f} s, total span = {span_s:.1f} s')
    else:
        timing = [ms_timing(src_name)]

    # JSON name: the last prefix without its time index, moved into RESULTS/
    # before the label is applied (see _apply_label)
    split_id  = image_identifier.rstrip('-')
    json_base = prefix_arr[-1]
    if f'{split_id}-t' in json_base:
        json_base = json_base.split(f'{split_id}-t')[0] + split_id
    json_path = _apply_label(json_base.replace(image_directory, 'RESULTS'), src_name, label) + '_polarization.json'

    if not OVERWRITE and os.path.isfile(json_path):
        msg(f'OVERWRITE is False and {json_path} exists -- skipping fitting, loading for plots')
        with open(json_path) as j:
            make_plots(json.load(j), src_name, prefix_arr, label, pol_flag)
        return 0

    # Channel fits run in a pool of spawned workers, started at the first
    # prefix with channels and restarted every recycle_interval prefixes to
    # bound the workers' memory growth (casacore caches)
    chan_pool, n_workers, recycle_interval, since_recycle = None, None, None, 0

    for k, prefix in enumerate(prefix_arr):
        msg('')
        msg(f'--- Prefix {k+1}/{n_epochs}: {os.path.basename(prefix)} ---')

        # MFS images: I, then P (Plin with a polarization angle calibrator,
        # else Ptot), Q, U, V
        images = {s: stokes_scan['mfs'][prefix][s] for s in ('I' if only_intensity else 'IPQUV')}
        msg(f'MFS images selected ({len(images)}):')
        for image in images.values():
            msg(f'  {os.path.basename(image)}')

        image_I  = images['I']
        freq_GHz = imhead(image_I, mode='get', hdkey='CRVAL3')['value'] / 1.0e9
        date_obs = imhead(image_I, mode='get', hdkey='DATE-OBS').replace('/', '-', 2).replace('/', 'T')
        bmaj     = imhead(image_I, mode='get', hdkey='bmaj')['value']
        bmin     = imhead(image_I, mode='get', hdkey='bmin')['value']
        bpa      = imhead(image_I, mode='get', hdkey='bpa')['value']
        msg(f'Observation: date={date_obs}  freq={freq_GHz:.4f} GHz  '
            f'beam={bmaj:.2f}"x{bmin:.2f}" PA={bpa:.1f} deg')

        if timing[k] is not None:
            time_ctr_mjd, time_dt = timing[k]
            msg(f'  Timing ({"image header" if n_epochs > 1 else "MS"}): centre MJD = {time_ctr_mjd:.8f}, delta = {time_dt:.0f} s')
        else:
            time_ctr_mjd, time_dt = _iso_to_mjd(date_obs), 15.0 * 60.0
            msg(f'  Timing (header fallback): centre MJD = {time_ctr_mjd:.8f}, delta = 900 s (15 min default)')

        # Input positions to pixels in this image
        src_ra_pix, src_dec_pix = [], []
        for ra, dec in zip(src_ra, src_dec):
            maxpos = imstat(image_I, region=f'circle[[{ra}deg,{dec}deg],1.0pix]')['maxpos']
            src_ra_pix.append(maxpos[0])
            src_dec_pix.append(maxpos[1])

        # ---- MFS Stokes I: the reference position for everything else ----
        # Upper-limit components are fitted at a fixed position
        fix = ['xyabp' if ulim else 'abp' for ulim in src_ulims]
        for z, ulim in enumerate(src_ulims):
            msg(f'  Component {z}: {"upper limit -- position will be fixed" if ulim else "detected -- position will be fitted freely"}')
        est_I = f'estimate_I_{src_name}.txt'
        make_estimate(est_I, image_I, src_ra_pix, src_dec_pix, fix)
        imfits = {'I': get_imfit_values(est_I, image_I, src_ra_pix, src_dec_pix)}
        components = comp_keys(imfits['I'])
        I_x, I_y = fitted_pix(imfits['I'], components)

        # Refit with every position fixed if any component drifted more than
        # MAX_I_DRIFT_PIX from its input position
        if MAX_I_DRIFT_PIX is not False:
            drift = np.hypot(np.subtract(I_x, src_ra_pix), np.subtract(I_y, src_dec_pix))
            for z, d in enumerate(drift):
                msg(f'  MFS Stokes I component {z} drift {d:.1f} pix (threshold {MAX_I_DRIFT_PIX} pix)'
                    + (' -- WARNING: refitting with positions fixed to the input' if d > MAX_I_DRIFT_PIX else ''))
            if np.any(drift > MAX_I_DRIFT_PIX):
                make_estimate(est_I, image_I, src_ra_pix, src_dec_pix, ['xyabp'] * len(src_ra_pix))
                imfits['I'] = get_imfit_values(est_I, image_I, src_ra_pix, src_dec_pix)
                I_x, I_y = fitted_pix(imfits['I'], components)

        # ---- MFS P, Q, U, V ----
        if not only_intensity:
            # P is fitted first with every component free to move; if that
            # fails or a component lands > 3 pix from Stokes I, it is refitted
            # with secondary components fixed (fix_additional_comps)
            est_P = f'estimate_P_{src_name}.txt'
            try:
                check_position(est_P, images['P'], I_x, I_y, P_image=True, manual_rms_region=manual_rms_region)
                imfits['P'] = get_imfit_values(est_P, images['P'], I_x, I_y)
                P_x, P_y = fitted_pix(imfits['P'], components)
                offset = float(np.max(np.hypot(np.subtract(P_x, I_x), np.subtract(P_y, I_y))))
                refit = offset > 3.0
            except Exception as e:
                msg(f'  P free-fit failed (likely complex multi-component field): {e}')
                offset, refit = -1.0, True
            if refit:
                msg(f'  Stokes I--P offset = {offset:.1f} pix; re-fitting P with '
                    f'secondary component positions fixed to Stokes I')
                check_position(est_P, images['P'], I_x, I_y, P_image=True,
                               fix_additional_comps=fix_additional_comps, manual_rms_region=manual_rms_region)
                imfits['P'] = get_imfit_values(est_P, images['P'], I_x, I_y)
            imfits['P'] = check_fit_offset_and_refix(imfits['P'], images['P'], I_x, I_y, bmaj, src_name, 'P')[0]
            P_x, P_y = fitted_pix(imfits['P'], components)

            # Q and U are anchored to the P position, V to Stokes I (V is weak
            # and leakage-dominated; FORCE_FIX_STOKES_V fixes it outright)
            pol = dict(bmaj=bmaj, src_name=src_name, manual_rms_region=manual_rms_region,
                       fix_additional_comps=fix_additional_comps)
            for s in 'QU':
                imfits[s] = fit_stokes(f'estimate_{s}_{src_name}.txt', images[s],
                                       P_x, P_y, P_x, P_y, stokes_label=s, **pol)
            imfits['V'] = fit_stokes(f'estimate_V_{src_name}.txt', images['V'],
                                     I_x, I_y, I_x, I_y, stokes_label='V',
                                     force_fix_all=FORCE_FIX_STOKES_V, **pol)

        # ---- MFS results ----
        msg('')
        msg('  --- MFS fit results (Stokes I is the astrometric reference) ---')
        for z, component in enumerate(components):
            res = component_results(imfits, images, component, manual_rms_region, pol_flag)
            store = output_dictionary['MFS'].setdefault(component, {'upperlimit': [src_ulims[z]]})
            record(store, {'freq_GHz': freq_GHz, 'bmaj_asec': bmaj, 'bmin_asec': bmin, 'bpa_deg': bpa,
                           **res, 'time_ctr_mjd': time_ctr_mjd,
                           'time_ctr_isot': _mjd_to_isot(time_ctr_mjd), 'time_dt': time_dt})
            msg(f'  {component}: RA={res["I_RA_deg"]:.6f} deg  Dec={res["I_DEC_deg"]:.6f} deg  '
                f'({I_x[z]:.1f}, {I_y[z]:.1f}) pix')
            msg(f'    I = {res["I_flux_mJy"]:.3f} +/- {res["I_err_mJy"]:.3f} mJy  '
                f'rms = {res["I_rms_mJy"]:.3f} mJy  S/N = {res["I_flux_mJy"] / res["I_rms_mJy"]:.1f}')
            if not only_intensity:
                msg(f'    P0 = {res["P0_flux_mJy"]:.3f} mJy  LP = {res["LP_frac"]:.2f}%  '
                    f'V = {res["V_flux_mJy"]:.3f} mJy')
        msg('')

        # ---- Channels: one job per channel, fitted in the pool ----
        msg(f'Fitting CHAN image(s) for prefix {k} with {src_im_suffix}: {prefix}')
        jobs = [{'chan_idx': i, 'CHAN_images': [str(im) for im in group], 'components': components,
                 'MFS_I_ra_pix': I_x, 'MFS_I_dec_pix': I_y,
                 'MFS_P_ra_pix': None if only_intensity else P_x,
                 'MFS_P_dec_pix': None if only_intensity else P_y,
                 'only_intensity': only_intensity, 'pol_flag': pol_flag,
                 'fix_additional_comps': fix_additional_comps, 'manual_rms_region': manual_rms_region,
                 'src_name': src_name, 'tag': f'{k}_{i}',
                 'mfs_freq_ghz': freq_GHz, 'mfs_bmaj_asec': bmaj}
                for i, group in enumerate(stokes_scan['chan_groups'][prefix])]

        chan_results = [None] * len(jobs)
        if jobs:
            if chan_pool is None:
                n_workers, recycle_interval = compute_max_workers([jobs[0]['CHAN_images'][0]], len(jobs))
                chan_pool = ProcessPoolExecutor(max_workers=n_workers, mp_context=mp.get_context('spawn'))
                since_recycle = 0
                msg(f'CHAN-fitting pool started with {n_workers} worker(s), '
                    f'recycle interval={recycle_interval} prefixes.')
            futures = [chan_pool.submit(fit_channel_group_batch, chunk) for chunk in chunk_jobs(jobs, n_workers)]
            for future in as_completed(futures):
                for chan_idx, result in future.result():
                    chan_results[chan_idx] = result

        # Record results in channel order
        for chan_idx, result in enumerate(chan_results):
            if result is None:
                msg(f'Fitting failed or channel skipped: chan_idx={chan_idx} is likely flagged')
                continue
            beam = {key: result[key] for key in ('freq_GHz', 'bmaj_asec', 'bmin_asec', 'bpa_deg')}
            for component in components:
                record(output_dictionary['CHAN'].setdefault(component, {}),
                       {**beam, **result['components'][component]}, k=k, n_epochs=n_epochs)

        # ---- Outputs for this prefix ----
        for component in components:
            if component not in output_dictionary['CHAN']:
                msg(f'  WARNING: no channel fits for {component} yet -- no spectra written')
                continue
            if not only_intensity:
                base = _apply_label(f'{prefix}_{component}'.replace(image_directory, 'RESULTS'), src_name, label)
                write_spectra(output_dictionary['CHAN'][component], k, base, pol_flag, time_ctr_mjd)
            # Spectral index of every epoch (stored in the JSON); plots come
            # after the JSON is saved
            plot_stokes_spectrum(output_dictionary, src_name, prefix, component, k,
                                 save_plot=False, label=label, pol_flag=pol_flag)

        cleanup_estimate_files(src_name)
        log_memory_usage(label=f'prefix {k+1}/{n_epochs} ({os.path.basename(prefix)}) done --')

        if chan_pool is not None:
            since_recycle += 1
            if since_recycle >= recycle_interval:
                msg(f'Recycling CHAN-fitting pool after {since_recycle} prefixes to bound worker memory.')
                chan_pool.shutdown(wait=True)
                chan_pool = None

    if chan_pool is not None:
        chan_pool.shutdown(wait=True)
        msg('CHAN-fitting pool shut down.')

    # Save everything before plotting, so a plotting error loses nothing
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, 'w') as j:
        json.dump(output_dictionary, j, indent=4)
    msg(f'  Output JSON saved: {json_path}')

    make_plots(output_dictionary, src_name, prefix_arr, label, pol_flag)
    return 0


# =============================================================================
# SPECTRAL INDEX AND PLOTS
# =============================================================================

def mad_inliers(values, nsigma):
    """Values within nsigma robust sigmas (1.4826 x MAD) of the median; all of them when the MAD is 0."""
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    sig = 1.4826 * np.nanmedian(np.abs(values - med))
    return np.abs(values - med) <= nsigma * sig if sig > 0 else np.ones(len(values), dtype=bool)


def iterative_mad_clip(values, nsigma):
    """Clip values more than nsigma robust sigmas from the median, repeated until none remain; True = clipped."""
    clipped = np.zeros(len(values), dtype=bool)
    while np.sum(~clipped) >= 3:
        med = np.nanmedian(values[~clipped])
        sig = 1.4826 * np.nanmedian(np.abs(values[~clipped] - med))
        if not sig > 0:
            break
        new = ~clipped & (np.abs(values - med) > nsigma * sig)
        if not new.any():
            break
        clipped |= new
    return clipped


def fit_spectral_index(freq, flux, rms, ref_freq):
    """
    Power-law fit log S = alpha log(nu / ref_freq) + c, weighted by S/N.
    Returns (alpha, alpha_err, chi2, ndof, fit_nu, fit_S).
    """
    w = flux / rms
    log_nu, log_S = np.log(freq / ref_freq), np.log(flux)
    coeffs, cov = np.polyfit(log_nu, log_S, 1, w=w, cov=True)
    chi2 = float(np.sum(((log_S - np.polyval(coeffs, log_nu)) * w) ** 2))
    fit_nu = np.linspace(freq.min(), freq.max(), 200)
    return coeffs[0], np.sqrt(cov[0, 0]), chi2, len(freq) - 2, fit_nu, np.exp(coeffs[1]) * (fit_nu / ref_freq) ** coeffs[0]


def set_robust_ylim(ax, series, refs):
    """
    y range covering, for each (flux, rms) in series, its PLOT_YLIM_MAD_SIGMA
    inliers +/- error bars, and the reference values, padded by PLOT_YLIM_PAD.
    Each series is clipped on its own, so panels given the same series share
    one range.
    """
    vals = [np.asarray(refs, dtype=float)]
    for flux, rms in series:
        keep = mad_inliers(flux, PLOT_YLIM_MAD_SIGMA)
        vals += [flux[keep] + rms[keep], flux[keep] - rms[keep]]
    vals = np.concatenate(vals)
    lo, hi = np.nanmin(vals), np.nanmax(vals)
    if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
        ax.set_ylim(lo - PLOT_YLIM_PAD * (hi - lo), hi + PLOT_YLIM_PAD * (hi - lo))


def _plot_detected_spectra(output_dictionary, src_name, prefix_arr, component, label, pol_flag=False):
    """Spectrum plots for the time-resolved epochs whose MFS Stokes I S/N >= SPEC_PLOT_SNR_THRESH."""
    mfs = output_dictionary['MFS'][component]
    mfs_rms = np.array(mfs['I_rms_mJy'])
    snr = np.where(mfs_rms > 0, np.array(mfs['I_flux_mJy']) / mfs_rms, 0.0)
    msg(f'  {component}: {int(np.sum(snr >= SPEC_PLOT_SNR_THRESH))}/{len(snr)} epochs with '
        f'MFS S/N >= {SPEC_PLOT_SNR_THRESH} -- plotting spectra')
    for k, s in enumerate(snr):
        saved = s >= SPEC_PLOT_SNR_THRESH
        msg(f'    epoch {k} ({os.path.basename(prefix_arr[k])}): MFS S/N = {s:.1f} -- '
            f'spectrum plot {"SAVED" if saved else "SKIPPED"}')
        if saved:
            plot_stokes_spectrum(output_dictionary, src_name, prefix_arr[k], component, k,
                                 save_plot=True, label=label, pol_flag=pol_flag)


def plot_stokes_spectrum(output_dictionary, src_name, prefix, component, k,
                         save_plot=True, label=None, pol_flag=False):
    """
    Spectral index of one component in epoch k, stored in its MFS results
    (alpha, alpha_err, chi2, ndof, chi2_red; None when not fitted), and with
    save_plot its I, Q, U, V spectrum plot (RESULTS/fitting_plots/).

    Channels are dropped where the rms of any Stokes is a 5-sigma MAD
    outlier, then where Stokes I is a SPEC_INDEX_MAD_CLIP outlier (iterated;
    marked red on the plot). The power law is fitted to the remaining
    channels with I > 0 and S/N >= 3, when the MFS S/N exceeds
    SPEC_INDEX_SNR_THRESH. With a polarization angle calibrator, Plin is
    overlaid on the Q and U panels (Ptot includes V, so is not shown).
    """
    label = label or src_name
    mfs = output_dictionary['MFS'][component]
    fit = dict(alpha=None, alpha_err=None, chi2=None, ndof=None, chi2_red=None)

    def store_fit():
        for key, val in fit.items():
            values = mfs.setdefault(key, [])
            values.extend([None] * (k + 1 - len(values)))
            values[k] = val

    # This epoch's channel values: freq and every flux and rms
    chan = output_dictionary['CHAN'].get(component, {})
    spec = {key: np.array(vals[k], dtype=float) for key, vals in chan.items()
            if key == 'freq_GHz' or key.endswith(('_flux_mJy', '_rms_mJy'))}
    stokes = 'IQUV' if 'Q_flux_mJy' in spec else 'I'
    if not len(spec.get('freq_GHz', [])):
        msg(f'  plot_stokes_spectrum: no channel data for {component} epoch {k}, skipping plot')
        store_fit()
        return

    # Drop channels with an outlier rms in any Stokes
    good = np.logical_and.reduce([mad_inliers(spec[f'{s}_rms_mJy'], 5.0) for s in stokes])
    if not good.all():
        msg(f'  Masking {int(np.sum(~good))}/{len(good)} channels with outlier RMS values')
    spec = {key: v[good] for key, v in spec.items()}
    if not good.any():
        msg(f'  plot_stokes_spectrum: all channels masked for {component} epoch {k}, skipping plot')
        store_fit()
        return

    # Spectral index, after clipping Stokes I outliers
    freq, I, I_rms = spec['freq_GHz'], spec['I_flux_mJy'], spec['I_rms_mJy']
    clipped = iterative_mad_clip(I, SPEC_INDEX_MAD_CLIP)
    if clipped.any():
        msg(f'  MAD clip ({component}, epoch {k}): {int(np.sum(clipped))} Stokes I outlier channel(s) '
            f'beyond {SPEC_INDEX_MAD_CLIP} sigma')
    mfs_I, mfs_freq = mfs['I_flux_mJy'][k], mfs['freq_GHz'][k]
    mfs_snr = mfs_I / mfs['I_rms_mJy'][k] if mfs['I_rms_mJy'][k] > 0 else 0.0
    valid = ~clipped & (I > 0) & (I_rms > 0) & (I / I_rms >= 3.0)
    fit_curve = None
    if mfs_snr <= SPEC_INDEX_SNR_THRESH:
        msg(f'  Spectral index fit skipped ({component}, epoch {k}): '
            f'MFS S/N = {mfs_snr:.1f} (< {SPEC_INDEX_SNR_THRESH})')
    elif valid.sum() < 3:
        msg(f'  Spectral index fit skipped ({component}, epoch {k}): '
            f'only {int(valid.sum())} unclipped channels with S/N >= 3 (need >= 3)')
    else:
        try:
            alpha, alpha_err, chi2, ndof, fit_nu, fit_S = fit_spectral_index(freq[valid], I[valid], I_rms[valid], mfs_freq)
            fit = dict(alpha=alpha, alpha_err=alpha_err, chi2=chi2, ndof=ndof,
                       chi2_red=chi2 / ndof if ndof > 0 else None)
            fit_curve = (fit_nu, fit_S)
            msg(f'  Spectral index fit ({component}, epoch {k}): alpha = {alpha:.3f} +/- {alpha_err:.3f}  '
                f'chi2 = {chi2:.2f}  ndof = {ndof}  (MFS S/N = {mfs_snr:.1f}, N_chan = {int(valid.sum())} / {len(freq)})')
        except Exception as e:
            msg(f'  WARNING: Spectral index fit failed for {component} epoch {k}: {e}')
    store_fit()
    if not save_plot:
        return

    # ---- Plot: one panel per Stokes; flux-clipped channels left out, marked red ----
    flagged_freqs = freq[clipped]
    spec = {key: v[~clipped] for key, v in spec.items()}
    freq = spec['freq_GHz']
    colours = {'I': 'royalblue', 'Q': 'darkorange', 'U': 'forestgreen', 'V': 'mediumpurple'}
    # Q and U share one y range, which also covers Plin where it is overlaid
    show_P = pol_flag and 'P_flux_mJy' in spec
    qu_series = [(spec[f'{s}_flux_mJy'], spec[f'{s}_rms_mJy']) for s in 'QU' if s in stokes]
    if show_P:
        qu_series.append((spec['P_flux_mJy'], spec['P_rms_mJy']))
    qu_refs = [mfs[f'{s}_flux_mJy'][k] for s in 'QU' if s in stokes] + [0.0]
    fig, axes = plt.subplots(len(stokes), 1, figsize=(10, 3 * len(stokes)), sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, s in zip(axes, stokes):
        flux, rms, mfs_val = spec[f'{s}_flux_mJy'], spec[f'{s}_rms_mJy'], mfs[f'{s}_flux_mJy'][k]
        ax.errorbar(freq, flux, yerr=rms, fmt='o', markersize=4, capsize=3, color=colours[s],
                    label=f'Stokes {s}', zorder=3)
        if s == 'I':
            ax.axhline(mfs_val, color='k', linestyle='--', linewidth=1.5,
                       label=f'MFS: {mfs_val:.2f} mJy  (S/N = {mfs_snr:.1f})', zorder=2)
            if fit_curve is not None:
                ax.plot(*fit_curve, color='tomato', linewidth=1.8, zorder=4,
                        label=rf'$\alpha$ = {fit["alpha"]:.2f} $\pm$ {fit["alpha_err"]:.2f}')
            set_robust_ylim(ax, [(flux, rms)], [mfs_val])
        else:
            ax.axhline(mfs_val, color='k', linestyle='--', linewidth=1.5, label=f'MFS: {mfs_val:.2f} mJy', zorder=2)
            ax.axhline(0, color='grey', linestyle=':', linewidth=0.8, zorder=1)
            if s in 'QU':
                set_robust_ylim(ax, qu_series, qu_refs)
            else:
                set_robust_ylim(ax, [(flux, rms)], [mfs_val, 0.0])
            if show_P and s in 'QU':
                ax.errorbar(freq, spec['P_flux_mJy'], yerr=spec['P_rms_mJy'], fmt='o', markersize=3, capsize=2,
                            color='black', alpha=0.5, label='Stokes P (lin)', zorder=2.5)
        ax.set_ylabel(f'{s}  (mJy/beam)', fontsize=10)
        ax.grid(True, which='both', alpha=0.3, linestyle=':')
        ax.legend(fontsize=9, loc='best')

    if len(flagged_freqs) and len(freq) >= 2:
        width = float(np.min(np.diff(np.sort(freq))))
        for f in flagged_freqs:
            for ax in axes:
                ax.axvspan(f - width / 2, f + width / 2, color='red', alpha=0.2, zorder=0, linewidth=0)

    axes[-1].set_xlabel('Frequency (GHz)', fontsize=11)
    fig.suptitle(f'{label}  |  {component}  |  epoch {k}', fontsize=12, y=1.01)
    plt.tight_layout()

    plot_dir = os.path.join('RESULTS', 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)
    plot_fname = os.path.join(plot_dir, f'{_apply_label(os.path.basename(prefix), src_name, label)}'
                                        f'_{component}_IQUV_spectrum.png')
    plt.savefig(plot_fname, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  IQUV spectrum plot saved: {plot_fname}')


def plot_light_curve(output_dictionary, src_name, prefix_arr, label=None):
    """
    MFS light curve of a time-resolved run: Stokes I (and Q, U, V when
    fitted) against time, one column per component. Every epoch is plotted;
    epochs whose error bar is a SPEC_INDEX_MAD_CLIP MAD outlier are left out
    of the y range only.
    """
    label = label or src_name
    mfs = output_dictionary['MFS']
    components = sorted(mfs)
    stokes = 'IQUV' if 'Q_flux_mJy' in mfs[components[0]] else 'I'
    colours = ['royalblue', 'tomato', 'forestgreen', 'mediumpurple',
               'darkorange', 'teal', 'crimson', 'goldenrod']

    fig, axes = plt.subplots(len(stokes), len(components), squeeze=False, sharex=True,
                             figsize=(8 * len(components), 3.5 * len(stokes)))
    all_dates = []
    for j, component in enumerate(components):
        comp = mfs[component]
        dates = [datetime.datetime.fromisoformat(d) for d in comp['time_ctr_isot']]
        all_dates += dates
        axes[0, j].set_title(component, fontsize=10)
        for i, s in enumerate(stokes):
            if f'{s}_flux_mJy' not in comp:
                continue
            ax = axes[i, j]
            flux, rms = np.array(comp[f'{s}_flux_mJy']), np.array(comp[f'{s}_rms_mJy'])
            ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8)
            ax.errorbar(dates, flux, yerr=rms, fmt='o', markersize=5, capsize=3,
                        color=colours[j % len(colours)], label=component if s == 'I' else None)
            if s == 'I':
                ax.legend(fontsize=8, loc='best')
            ax.set_ylabel(f'Stokes {s} (mJy/beam)', fontsize=10)
            ax.grid(True, alpha=0.3, linestyle=':')

            keep = mad_inliers(rms, SPEC_INDEX_MAD_CLIP)
            if not keep.any():
                keep[:] = True
            lo, hi = np.nanmin(flux[keep] - rms[keep]), np.nanmax(flux[keep] + rms[keep])
            if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
                ax.set_ylim(lo - PLOT_YLIM_PAD * (hi - lo), hi + PLOT_YLIM_PAD * (hi - lo))

    # Time axis: HH:MM, with the date when the run spans more than a day
    multi_day = all_dates and (max(all_dates) - min(all_dates)).total_seconds() > 86400
    for ax in axes[-1, :]:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M' if multi_day else '%H:%M'))
        ax.set_xlabel('Time (UT)', fontsize=10)
        for tick in ax.get_xticklabels():
            tick.set_rotation(30)
            tick.set_ha('right')

    fig.suptitle(f'{label}  |  MFS light curve', fontsize=12, y=1.02)
    plt.tight_layout()

    # File name: the first prefix without its time index (e.g. -t0051)
    plot_dir = os.path.join('RESULTS', 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)
    lc_basename = _apply_label(re.sub(r'-t\d+$', '', os.path.basename(prefix_arr[0])), src_name, label)
    plot_fname = os.path.join(plot_dir, f'{lc_basename}_{label}_MFS_lightcurve.png')
    plt.savefig(plot_fname, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  MFS light curve plot saved: {plot_fname}')


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Extract every source listed in rmsynth_info.json (cfg.RMSYN_INFO_FILE)."""
    with open(cfg.RMSYN_INFO_FILE) as j:
        rmsynth_info = json.load(j)

    pol_flag = cfg.POLANG_NAME != ''   # a polarization angle calibrator was used
    labels = rmsynth_info.get('label', [])   # optional output names; '' = source name

    for k in range(len(rmsynth_info['image_directory'])):
        if not os.path.exists(rmsynth_info['image_directory'][k]):
            msg('Skipping {} as the directory does not exist'.format(rmsynth_info['image_directory'][k]))
            continue

        # Images: <image_directory>/*<source>*.ms_<identifier>- (then t<N>- when time-resolved)
        src_im_identifier = (f'{cfg.CWD}/{rmsynth_info["image_directory"][k]}/'
                             f'*{case_insensitive(rmsynth_info["source_name"][k])}*.ms_'
                             f'{rmsynth_info["image_identifier"][k]}-' + ('t' if rmsynth_info["image_timing"][k] else ''))

        # Source positions: 'HH:MM:SS.SS,+DD.MM.SS.SS' (colon declinations are
        # converted) or decimal degrees
        src_ra  = []
        src_dec = []
        for pos in rmsynth_info['source_pos'][k]:
            pos = fix_dec_colons(pos, 'source position')
            ra_deg, dec_deg = parse_casa_position(pos)
            if ra_deg is None:
                msg(f'ERROR: Could not parse position "{pos}" -- skipping source')
                continue
            msg(f'Parsed position: RA = {ra_deg:.6f} deg, Dec = {dec_deg:.6f} deg')
            src_ra.append(ra_deg)
            src_dec.append(dec_deg)

        src_ra, src_dec = np.array(src_ra), np.array(src_dec)

        msg(f'Starting property extraction for identifier: {src_im_identifier}')
        extract_polarization_properties(
            rmsynth_info['source_name'][k], src_im_identifier, rmsynth_info['image_suffix'][k],
            src_ra, src_dec, rmsynth_info['source_ulim'][k], pol_flag,
            fix_dec_colons(rmsynth_info['rms_region'][k], 'rms_region'),
            rmsynth_info['image_directory'][k], rmsynth_info['image_identifier'][k],
            fix_additional_comps=True, out_label=labels[k] if k < len(labels) else '')


if __name__ == "__main__":
    main()
