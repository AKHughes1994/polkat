# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

# Modular-CASA version: imfit/imstat/imhead/table are plain imports, not
# casashell-injected globals, so this runs under a normal python3 interpreter
# (not `casa -c`) as long as casatasks/casatools are installed.

import glob
import os
import re
import datetime
import subprocess
import sys
import json
import shutil
import numpy as np
import os.path as o
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.spatial.distance import cdist
from astropy.io import fits
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
# Global fitting, plotting, and spectral-index options
# (defaults pulled from config.py; override here if needed)
# =============================================================================

FORCE_FIX_STOKES_V    = cfg.RMSYN_FORCE_FIX_STOKES_V
OVERWRITE             = cfg.RMSYN_OVERWRITE
SPEC_PLOT_SNR_THRESH  = cfg.RMSYN_SPEC_PLOT_SNR_THRESH
SPEC_INDEX_SNR_THRESH = cfg.RMSYN_SPEC_INDEX_SNR_THRESH
SPEC_INDEX_MAD_CLIP   = cfg.RMSYN_SPEC_INDEX_MAD_CLIP
MAX_I_DRIFT_PIX       = cfg.RMSYN_MAX_I_DRIFT_PIX

# =============================================================================


def _iso_to_mjd(date_str):
    """Convert an ISO-format datetime string (YYYY-MM-DDTHH:MM:SS[.f]) to MJD."""
    _epoch = datetime.datetime(1858, 11, 17, 0, 0, 0)
    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S'):
        try:
            return (datetime.datetime.strptime(date_str, fmt) - _epoch).total_seconds() / 86400.0
        except ValueError:
            continue
    raise ValueError(f'Cannot parse date string for MJD conversion: {date_str!r}')


def _mjd_to_isot(mjd):
    """Convert MJD (float) to an ISOT-format datetime string (inverse of _iso_to_mjd)."""
    _epoch = datetime.datetime(1858, 11, 17, 0, 0, 0)
    return (_epoch + datetime.timedelta(days=mjd)).isoformat()


def _apply_label(path, src_name, label):
    """
    Substitute `label` for `src_name` in the filename component of `path`
    only -- never in directory segments.

    A plain path.replace(src_name, label) across the *whole* path is unsafe:
    if src_name also happens to appear in a directory segment (image_directory,
    CWD, etc.) that occurrence gets silently corrupted too. In particular, if
    this runs before an image_directory -> 'RESULTS' relocation, the corrupted
    directory segment no longer matches image_directory verbatim, that
    relocation silently no-ops, and the file is written under the original
    image directory instead of RESULTS -- with no error, just a "missing"
    output file. Always relocate to RESULTS (or otherwise fix up directories)
    using the untouched path *before* calling this.
    """
    if label == src_name:
        return path
    dirname, basename = os.path.split(path)
    return os.path.join(dirname, basename.replace(src_name, label))


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush=True)


# =============================================================================
# Channel-fit parallelization: worker count + job chunking.
# =============================================================================

def _proc_rss_mb(pid):
    """Read a process's resident memory (MB) straight from /proc -- no
    psutil dependency needed."""
    try:
        with open(f'/proc/{pid}/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return int(line.split()[1]) / 1024.0
    except OSError:
        pass
    return 0.0


def log_memory_usage(label=''):
    """Log main-process + CHAN-fitting pool worker memory (RSS): per-worker
    avg/min/max, and the delta since the last call so steady growth (e.g. a
    casacore table/image cache leak) shows up without diffing log lines.
    Prefers psutil; falls back to /proc if it isn't installed."""
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
    """Conservative peak estimate for one channel-fit worker: pixel data size
    (x4 for headroom) or a ~280 MB floor for a fresh casatasks/casatools
    import per spawned worker, whichever is larger."""
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


# Assumed per-worker memory growth rate (MB/prefix), used by compute_max_workers
# to size the pool recycle interval. Clamped to [MIN,MAX] so it's never too
# small (constant recycling overhead) or too large (no protection against leaks).
ASSUMED_WORKER_GROWTH_MB_PER_PREFIX = 5.0
MIN_POOL_RECYCLE_INTERVAL = 10
MAX_POOL_RECYCLE_INTERVAL = 200


def compute_max_workers(sample_images, n_jobs):
    """sample_images is only used to estimate per-worker memory (a few
    representative images is enough). n_jobs (e.g. channel groups) caps the
    worker count, not len(sample_images).

    Returns (workers, recycle_interval): recycle_interval scales with the
    memory headroom left per worker once its baseline footprint is
    subtracted, divided by the assumed growth rate -- more workers or less
    memory both mean less headroom, so a smaller interval."""
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

    # Logged step by step (not just the final number) so it's clear whether a
    # small interval is from tight memory or from many workers splitting it.
    mem_per_worker_budget = available / workers
    headroom_bytes        = max(0.0, mem_per_worker_budget - per_worker)
    growth_bytes_per_pfx  = ASSUMED_WORKER_GROWTH_MB_PER_PREFIX * 1024 ** 2
    raw_interval          = int(headroom_bytes / growth_bytes_per_pfx)
    recycle_interval      = max(MIN_POOL_RECYCLE_INTERVAL, min(MAX_POOL_RECYCLE_INTERVAL, raw_interval))
    msg(f'Recycle interval calc: budget/worker = available / workers '
        f'= {available / 1024**2:.0f}MB / {workers} = {mem_per_worker_budget / 1024**2:.0f}MB')
    msg(f'Recycle interval calc: headroom = budget/worker - est. memory/worker '
        f'= {mem_per_worker_budget / 1024**2:.0f}MB - {per_worker / 1024**2:.0f}MB '
        f'= {headroom_bytes / 1024**2:.0f}MB')
    msg(f'Recycle interval calc: raw = headroom / assumed growth '
        f'= {headroom_bytes / 1024**2:.0f}MB / {ASSUMED_WORKER_GROWTH_MB_PER_PREFIX:.0f}MB/prefix '
        f'= {raw_interval} prefixes')
    msg(f'Pool recycle interval : {recycle_interval} prefixes  '
        f'[clamped to [{MIN_POOL_RECYCLE_INTERVAL}, {MAX_POOL_RECYCLE_INTERVAL}]]')

    return workers, recycle_interval


def chunk_jobs(jobs, n_chunks):
    """Round-robin split so any systematic per-channel cost trend (e.g. RFI
    flagging clustered in frequency) gets spread evenly across workers."""
    n_chunks = max(1, min(n_chunks, len(jobs)))
    chunks = [[] for _ in range(n_chunks)]
    for i, job in enumerate(jobs):
        chunks[i % n_chunks].append(job)
    return [c for c in chunks if c]


_DEC_COLON_RE = re.compile(r'[+-]\d{1,3}:\d{1,2}:\d{1,2}(?:\.\d+)?')


def fix_dec_colons(pos_str, context_label=''):
    """
    Convert a colon-separated declination (HMS,DMS style, e.g. '-63:42:45.6')
    to CASA's period-separated style ('-63.42.45.6'), leaving RA untouched.

    Declination is identified as a signed sexagesimal token (only Dec carries
    a +/- sign in these position/region strings). If no such token is found,
    the string is returned unchanged.

    Parameters
    ----------
    pos_str       : str  -- position or region string to check
    context_label : str  -- what this string is, for the NOTE message (e.g.
                             'source position', 'rms_region')

    Returns
    -------
    str : the string with any colon-separated declination converted to periods,
          or pos_str itself unchanged if it isn't a (non-empty) string -- e.g.
          rms_region defaults to False/'' when no manual region is set.
    """
    if not isinstance(pos_str, str) or not pos_str:
        return pos_str
    fixed = _DEC_COLON_RE.sub(lambda m: m.group(0).replace(':', '.'), pos_str)
    if fixed != pos_str:
        msg(f'NOTE: {context_label} declination uses colons -- converting to CASA period format: '
            f'"{pos_str}" -> "{fixed}"')
    return fixed


def parse_casa_position(position_str):
    """
    Parse a CASA-format position string to decimal degrees.

    Accepts two formats:
      - CASA HMS/DMS : 'HH:MM:SS.SS,±DD.MM.SS.SS'   (colon RA, period Dec)
      - Decimal degrees: '260.430458,-16.204894'  or with 'deg' suffix

    Parameters
    ----------
    position_str : str
        Position string from rmsynth_info.json source_pos field.

    Returns
    -------
    ra_deg, dec_deg : float, float
        Position in decimal degrees, or (None, None) on parse failure.
    """
    try:
        # Decimal-degree format (no colons in the RA part)
        if ':' not in position_str:
            parts = position_str.replace('deg', '').split(',')
            return float(parts[0].strip()), float(parts[1].strip())

        # CASA HMS/DMS format: 'HH:MM:SS.SS,±DD.MM.SS.SS'
        ra_str, dec_str = position_str.split(',')

        # --- Parse RA (colon-separated: HH:MM:SS.SS) ---
        ra_parts = ra_str.strip().split(':')
        if len(ra_parts) != 3:
            msg(f'WARNING: Expected 3 parts in RA "{ra_str}", got {len(ra_parts)}')
            return None, None
        ra_deg = (float(ra_parts[0]) + float(ra_parts[1]) / 60.0
                  + float(ra_parts[2]) / 3600.0) * 15.0  # hours -> degrees

        # --- Parse Dec (period-separated: ±DD.MM.SS.SS) ---
        dec_str = dec_str.strip()
        dec_sign = -1.0 if dec_str.startswith('-') else 1.0
        dec_str  = dec_str.lstrip('+-')

        dec_parts = dec_str.split('.')
        if len(dec_parts) < 3:
            msg(f'WARNING: Expected at least 3 parts in Dec "{dec_str}", got {len(dec_parts)}')
            return None, None

        # Reassemble fractional arcseconds from any extra period-separated parts
        dec_arcsec = float(dec_parts[2] + ('.' + dec_parts[3] if len(dec_parts) > 3 else ''))
        dec_deg = dec_sign * (abs(float(dec_parts[0]))
                              + float(dec_parts[1]) / 60.0
                              + dec_arcsec / 3600.0)

        return ra_deg, dec_deg

    except Exception as e:
        msg(f'ERROR parsing CASA position "{position_str}": {e}')
        return None, None


def get_imfit_values(fname, image, xpix, ypix):

    '''
    Run imfit on a specific image
    Inputs: 
        fname = string containing estimate file name
        image = string containing path to image to be fit
        xpix  = pixel coordinate (RA)
        ypix  = pixel coordinate (Dec)
    Returns:
        imfit dictionary
    '''

    # Get the beam parameters
    bmaj        = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin        = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa         = imhead(image, mode='get', hdkey = 'BPA')['value']
    pixel_asec  = imhead(image, mode='get', hdkey='cdelt2')['value'] * 3600 * 180.0 / np.pi # pixel size in asecs


    # Get the number of components
    n_comp = len(xpix)
  
    # For single components
    if n_comp == 1: 
        x = xpix[0]
        y = ypix[0]
        r = 5 * bmaj
        src_region = f'circle[[{x}pix,{y}pix],{r}arcsec]'

    # For multi-component fits (times out if components are very far apart)
    else:
        point_array = np.array([xpix, ypix]).T
        max_dist = np.amax(cdist(point_array, point_array)) * pixel_asec

        # Bounding region radius: twice the max separation, or 5*bmaj, whichever is bigger
        x = xpix[0]
        y = ypix[0]
        r = np.amax((5 * bmaj, 2.0 * max_dist))
        src_region = f'circle[[{x}pix,{y}pix],{r}arcsec]'
   
    return imfit(image, estimates = fname, region = src_region)
    
    

def calculate_P0(flux_P, rms_Q, rms_U, rms_V, pol_flag, Aq = 0.8):
    '''
    Calculate the de-biased linearly polarized flux
    '''    

    # If polarization angle calibrator included
    if pol_flag:

        # Get the noise ratio coeffs, and calculate noise, from Hales 2012. https://arxiv.org/abs/1205.5310
        if rms_Q >= rms_U:
            A = Aq
            B = 1. - Aq 
        else: 
            B = Aq
            A = 1. - Aq 

        rms_P = (A * rms_Q ** 2 + B * rms_U ** 2) ** 0.5 

        # De-bias if SNR >= 3, following Vaillancourt 2006. https://arxiv.org/abs/astro-ph/0603110
        if flux_P / rms_P >= 3:
            flux_P0 = (flux_P ** 2 - rms_P ** 2) ** (0.5)

        else:
            flux_P0 = flux_P

    # No polarization angle calibrator: no reference for this case, so go
    # conservative -- adopt the max of Q/U/V as the noise, and always de-bias
    # (Monte Carlo experiments put the bias correction factor at 2 for
    # P^2 = Q^2 + U^2 + V^2).
    else:
        rms_P = np.amax([rms_Q, rms_U, rms_V])
        flux_P0 = (flux_P ** 2 - 2.0 * rms_P ** 2) ** 0.5

    return flux_P0, rms_P

        

def return_max(im, region):
    '''
    Return the value that has the higher absolute magnitude
    Necessary for fluxes that are non-positive definite
    '''

    ims = imstat(im, region=region)
    fmax = ims['max'][0]
    fmin = ims['min'][0]

    if fmax > abs(fmin):
        return fmax
    else:
        return fmin


def get_imstat_values(image, xpix, ypix, manual_rms_region = False):
    '''
    Take in an image and a position, 
    return the max, max pixel location(s), and rms
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']    

    # Define the regions of interest (rms is ~500 beam area)
    r_in  = 5.0 * bmaj
    r_out = np.sqrt(500 * 0.25 * bmaj * bmin + r_in ** 2)
    src_region = f'circle[[{xpix}pix,{ypix}pix],3.0pix]'
    rms_region = f'annulus[[{xpix}pix,{ypix}pix],[{r_in}arcsec,{r_out}arcsec]]'
    if manual_rms_region:
        rms_region = manual_rms_region
    #else:
    #    msg('Using default annular RMS region')

    # Values of interest -- Source
    ims = imstat(image, region = src_region)
    flux = return_max(image, src_region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]

    # Extract RMS
    rms = imstat(image, region = rms_region)['rms'][0]

    return [flux, xpix, ypix, rms, rms_region]
    
    
def check_position(fname, image, xpix, ypix, snr_thresh=5.0, P_image=False,
                   fix_additional_comps=False, manual_rms_region=False,
                   src_name='src', force_fix_all=False):
    '''
    Decide per component whether to fit position freely or fix it, then write
    the imfit estimate file.

    Decision logic:
      - force_fix_all=True                          --> all fixed ('xyabp'), no S/N check
      - k>0 component and fix_additional_comps=True  --> always fixed ('xyabp')
                                                          (faint ejecta/lobes: position is
                                                          unreliable off Stokes I, so anchor it)
      - k=0, S/N >= snr_thresh                       --> free ('abp')
      - k=0, S/N <  snr_thresh                       --> fixed ('xyabp')

    P images have non-Gaussian (Rice) noise, so the S/N check uses the larger
    of rms_Q/rms_U instead of measuring P directly.

    Parameters
    ----------
    fname               -- path for the output estimate file passed to imfit
    image               -- CASA image to be fitted
    xpix                -- list of RA pixel positions, one per component
    ypix                -- list of Dec pixel positions, one per component
    snr_thresh          -- minimum S/N to allow a free position fit (default 5)
    P_image             -- set True for polarized-intensity images so the noise
                           estimate is taken from Q/U rather than P directly
    fix_additional_comps-- set True to force secondary components (k>0) to fixed
                           positions regardless of S/N (recommended for faint ejecta)
    manual_rms_region   -- optional CASA region string to override the default
                           annular RMS region
    src_name            -- source name string; used to build unique temp filenames
    force_fix_all       -- set True to fix ALL components (including k=0) regardless
                           of S/N; overrides both snr_thresh and fix_additional_comps
    '''

    fix_var = []

    # Read beam shape once; it is the same for all components in this image
    bmaj = imhead(image, mode='get', hdkey='bmaj')['value']
    bmin = imhead(image, mode='get', hdkey='bmin')['value']
    bpa  = imhead(image, mode='get', hdkey='bpa')['value']

    # Also get the pixel scale so we can report offsets in arcsec
    cell_rad    = abs(imhead(image, mode='get', hdkey='CDELT2')['value'])
    cell_arcsec = np.degrees(cell_rad) * 3600.0

    for k, (x, y) in enumerate(zip(xpix, ypix)):

        # Unique per component/source name so concurrent Stokes passes don't clobber each other
        check_file = f'check_pos_{os.path.basename(fname)}_{k}.txt'
        with open(check_file, 'w') as f:
            # Zero amplitude, fixed position -- just want the peak flux at this pixel
            f.write(f'0.0,{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, xyabp')

        # Fit region: circle of radius 2*BMAJ centred on the trial position
        region    = f'circle[[{x}pix,{y}pix],{2 * bmaj}arcsec]'
        test_flux = abs(imfit(image, region=region,
                              estimates=check_file)['results']['component0']['peak']['value'])

        # P images have Rice-distributed noise -- use the larger of rms_Q/rms_U instead
        rms = get_imstat_values(image, x, y, manual_rms_region=manual_rms_region)[3]
        if P_image:
            image_Q = image.replace('-Plin-', '-Q-').replace('-Ptot-', '-Q-')
            image_U = image.replace('-Plin-', '-U-').replace('-Ptot-', '-U-')
            rms_Q   = get_imstat_values(image_Q, x, y, manual_rms_region=manual_rms_region)[3]
            rms_U   = get_imstat_values(image_U, x, y, manual_rms_region=manual_rms_region)[3]
            rms     = np.amax((rms_Q, rms_U))

        snr = test_flux / rms

        if force_fix_all:
            fix_var.append('xyabp')
            msg(f'  Component {k}: FIXED (force_fix_all=True) '
                f'| S/N = {snr:.1f} | {os.path.basename(image)}')

        elif fix_additional_comps and k > 0:
            fix_var.append('xyabp')
            msg(f'  Component {k}: FIXED (secondary component flag set) '
                f'| S/N = {snr:.1f} | {os.path.basename(image)}')

        elif test_flux > snr_thresh * rms:
            fix_var.append('abp')
            msg(f'  Component {k}: FREE position fit '
                f'| S/N = {snr:.1f} (>= {snr_thresh}) '
                f'| {os.path.basename(image)}')

        else:
            fix_var.append('xyabp')
            msg(f'  Component {k}: FIXED position (low S/N) '
                f'| S/N = {snr:.1f} (< {snr_thresh}) '
                f'| {os.path.basename(image)}')

    # Write the final estimate file that imfit will actually use
    make_estimate(fname, image, xpix, ypix, fix_var)


def make_estimate(fname, image, xpix, ypix, fix_var):
    '''
    Write a CASA imfit estimate file: one line per component, seeding imfit
    with a starting position and beam-sized Gaussian. fix_var controls which
    parameters are free ('abp' = amplitude/axes/PA only; 'xyabp' = also fixes x/y).

    Inputs:
        fname   -- output estimate file path
        image   -- CASA image from which to read the beam and peak flux seed
        xpix    -- list of RA pixel positions for each component
        ypix    -- list of Dec pixel positions for each component
        fix_var -- list of CASA fix strings, one per component (e.g. 'abp' or 'xyabp')
    '''

    # Read beam shape from the image header to initialise the Gaussian seed
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']
    
    # Write one line per component: amplitude, x, y, beam major, beam minor, PA, fix string
    f = open(fname, 'w')
    for x, y, fix in zip(xpix, ypix, fix_var):
        ims = get_imstat_values(image, x, y)
        f.write(f'{ims[0]},{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix}\n')
    f.close()            
    
    return 0    


def check_fit_offset_and_refix(imfit_result, image, ref_xpix, ref_ypix, bmaj,
                               src_name, stokes_label, manual_rms_region=False, tag=''):
    '''
    Check each fitted component hasn't drifted more than 1/3 beam BMAJ from
    its reference position (normally the Stokes I result); genuine sub-beam
    offsets between Stokes parameters should be much smaller than that, so
    anything larger is almost certainly a noise-driven imfit excursion. If
    any component exceeds it, the whole image is re-fit with all positions
    fixed to the reference, so one wandering component can't contaminate the
    others' fluxes.

    Parameters
    ----------
    imfit_result    : dict  -- result dictionary returned by get_imfit_values()
    image           : str   -- path to the CASA image (used to convert pixels to arcsec)
    ref_xpix        : list  -- reference RA pixel positions, one per component
                               (typically the Stokes I fitted positions)
    ref_ypix        : list  -- reference Dec pixel positions, one per component
    bmaj            : float -- beam major axis in arcsec; sets the offset threshold
    src_name        : str   -- source name, used to build a unique estimate filename
    stokes_label    : str   -- Stokes label for log messages (e.g. 'Q', 'Plin', 'V')
    manual_rms_region : str or False -- passed through to make_estimate() if a
                                        re-fit is needed
    tag             : str   -- appended to the re-fit estimate filename to keep
                               it unique across concurrently running channel-fit
                               workers on the same src_name (default '' for the
                               serial MFS-fitting call sites)

    Returns
    -------
    imfit_result : dict  -- original result, or updated result if a re-fit was done
    was_refitted : bool  -- True if at least one component exceeded the threshold
    '''

    # CDELT2 is in radians -- convert to arcsec, then BMAJ to pixels
    cell_rad    = abs(imhead(image, mode='get', hdkey='CDELT2')['value'])
    cell_arcsec = np.degrees(cell_rad) * 3600.0
    bmaj_pixels = bmaj / cell_arcsec

    # Threshold: flag fits that move more than 1/3 beam from the reference position
    threshold_fraction = 1.0 / 3.0
    threshold_pixels   = threshold_fraction * bmaj_pixels
    threshold_arcsec   = threshold_fraction * bmaj

    # imfit results also contain non-component metadata keys ('nelements' etc.) -- exclude them
    components = [key for key in imfit_result['results'].keys() if 'component' in key]
    any_exceeded = False

    for z, component in enumerate(components):
        fitted_xpix = imfit_result['results'][component]['pixelcoords'][0]
        fitted_ypix = imfit_result['results'][component]['pixelcoords'][1]

        # Euclidean pixel distance from the Stokes I reference for this component
        pixel_dist    = np.sqrt((fitted_xpix - ref_xpix[z])**2 +
                                (fitted_ypix - ref_ypix[z])**2)
        offset_arcsec = pixel_dist * cell_arcsec

        if pixel_dist > threshold_pixels:
            msg(f'  WARNING: Stokes {stokes_label} component {z} offset = '
                f'{offset_arcsec:.2f}" ({pixel_dist:.1f} pix) '
                f'exceeds 1/3 BMAJ ({threshold_arcsec:.2f}"). '
                f'Will re-fit entire image with all positions fixed to Stokes I reference.')
            any_exceeded = True
        else:
            msg(f'  Position check OK: Stokes {stokes_label} component {z} offset = '
                f'{offset_arcsec:.2f}" ({pixel_dist:.1f} pix) '
                f'< 1/3 BMAJ ({threshold_arcsec:.2f}")')

    if any_exceeded:
        # Refit with all positions fixed to reference. Filename is unique per
        # Stokes param so it doesn't clobber another still-in-use estimate file.
        fixed_est = f'estimate_fixed_{stokes_label}_{src_name}{tag}.txt'
        fix_var   = ['xyabp'] * len(ref_xpix)
        make_estimate(fixed_est, image, ref_xpix, ref_ypix, fix_var)
        imfit_result = get_imfit_values(fixed_est, image, ref_xpix, ref_ypix)
        msg(f'  Re-fit complete: all Stokes {stokes_label} positions fixed to reference.')

    return imfit_result, any_exceeded


def cleanup_estimate_files(src_name):
    '''
    Remove the temporary estimate/check_position scratch files written for
    src_name (estimate_I/P/Q/U/V_*, estimate_fixed_*, check_pos_*).

    Called once per prefix/timestep. CHAN-fit worker files are uniquely
    tagged per channel per timestep to avoid collisions between concurrent
    workers, so without this they'd pile up -- one file per channel per
    Stokes per timestep, across a run that can have hundreds of timesteps.
    '''
    patterns = [f'estimate_*{src_name}*.txt', f'check_pos_*{src_name}*.txt']
    removed = 0
    for pattern in patterns:
        for fpath in glob.glob(pattern):
            try:
                os.remove(fpath)
                removed += 1
            except OSError:
                pass
    if removed:
        msg(f'Cleaned up {removed} temporary estimate/check-position file(s) for {src_name}.')


def fit_channel_group(job):
    '''
    Fit Stokes I(/P/Q/U/V) for one channel-image group of one prefix/timestep.
    Self-contained (no output_dictionary access) so it can run inside a
    parallel worker: everything comes in via `job`, results go out as a plain
    dict.

    `job` keys: chan_idx, CHAN_images, components, MFS_I_ra_pix, MFS_I_dec_pix,
    MFS_P_ra_pix, MFS_P_dec_pix, only_intensity, pol_flag, fix_additional_comps,
    manual_rms_region, src_name, tag, mfs_freq_ghz, mfs_bmaj_asec.

    Returns (chan_idx, result) where result is None if the channel was skipped
    (elongated beam -- likely heavily flagged) or the fit raised, otherwise a
    dict with the per-channel scalars (freq_GHz/bmaj_asec/bmin_asec/bpa_deg)
    and a 'components' sub-dict keyed the same as `components`.
    '''
    chan_idx             = job['chan_idx']
    CHAN_images          = job['CHAN_images']
    components           = job['components']
    MFS_I_ra_pix         = job['MFS_I_ra_pix']
    MFS_I_dec_pix        = job['MFS_I_dec_pix']
    MFS_P_ra_pix         = job['MFS_P_ra_pix']
    MFS_P_dec_pix        = job['MFS_P_dec_pix']
    only_intensity       = job['only_intensity']
    pol_flag             = job['pol_flag']
    fix_additional_comps = job['fix_additional_comps']
    manual_rms_region    = job['manual_rms_region']
    src_name             = job['src_name']
    tag                  = job['tag']
    mfs_freq_ghz         = job['mfs_freq_ghz']
    mfs_bmaj_asec        = job['mfs_bmaj_asec']

    try:
        msg(f'Fitting image: {CHAN_images[0]}')

        bmaj     = imhead(CHAN_images[0], mode='get', hdkey='bmaj')['value']
        bmin     = imhead(CHAN_images[0], mode='get', hdkey='bmin')['value']
        bpa      = imhead(CHAN_images[0], mode='get', hdkey='bpa')['value']
        freq_GHz = imhead(CHAN_images[0], mode='get', hdkey='CRVAL3')['value'] / 1.0e9

        # Check if the beam is significantly elongated relative to the MFS
        # expectation (channel is probably very flagged and should be omitted)
        beam_scaled = mfs_freq_ghz / freq_GHz * mfs_bmaj_asec
        if bmaj > 3.0 * beam_scaled:
            msg('Skipping Channel: BMAJ is very different from expectation based on MFS image (likely high flagged)')
            return chan_idx, None

        est_I = f'estimate_I_{src_name}_{tag}.txt'
        check_position(est_I, CHAN_images[0], MFS_I_ra_pix, MFS_I_dec_pix,
                        manual_rms_region=manual_rms_region, src_name=src_name)
        CHAN_I_imfit    = get_imfit_values(est_I, CHAN_images[0], MFS_I_ra_pix, MFS_I_dec_pix)
        CHAN_I_imfit, _ = check_fit_offset_and_refix(
            CHAN_I_imfit, CHAN_images[0], MFS_I_ra_pix, MFS_I_dec_pix,
            bmaj, src_name, 'I_chan', manual_rms_region=manual_rms_region, tag=tag)
        CHAN_I_ra_pix  = [CHAN_I_imfit['results'][key]['pixelcoords'][0] for key in components]
        CHAN_I_dec_pix = [CHAN_I_imfit['results'][key]['pixelcoords'][1] for key in components]

        if not only_intensity:
            est_P = f'estimate_P_{src_name}_{tag}.txt'
            check_position(est_P, CHAN_images[1], MFS_P_ra_pix, MFS_P_dec_pix,
                            P_image=True, fix_additional_comps=fix_additional_comps,
                            manual_rms_region=manual_rms_region, src_name=src_name)
            CHAN_P_imfit    = get_imfit_values(est_P, CHAN_images[1], CHAN_I_ra_pix, CHAN_I_dec_pix)
            CHAN_P_imfit, _ = check_fit_offset_and_refix(
                CHAN_P_imfit, CHAN_images[1], CHAN_I_ra_pix, CHAN_I_dec_pix,
                bmaj, src_name, 'P', manual_rms_region=manual_rms_region, tag=tag)
            CHAN_P_ra_pix  = [CHAN_P_imfit['results'][key]['pixelcoords'][0] for key in components]
            CHAN_P_dec_pix = [CHAN_P_imfit['results'][key]['pixelcoords'][1] for key in components]

            est_Q = f'estimate_Q_{src_name}_{tag}.txt'
            check_position(est_Q, CHAN_images[2], CHAN_P_ra_pix, CHAN_P_dec_pix,
                            P_image=False, fix_additional_comps=fix_additional_comps,
                            manual_rms_region=manual_rms_region, src_name=src_name)
            CHAN_Q_imfit    = get_imfit_values(est_Q, CHAN_images[2], CHAN_P_ra_pix, CHAN_P_dec_pix)
            CHAN_Q_imfit, _ = check_fit_offset_and_refix(
                CHAN_Q_imfit, CHAN_images[2], CHAN_P_ra_pix, CHAN_P_dec_pix,
                bmaj, src_name, 'Q', manual_rms_region=manual_rms_region, tag=tag)

            est_U = f'estimate_U_{src_name}_{tag}.txt'
            check_position(est_U, CHAN_images[3], CHAN_P_ra_pix, CHAN_P_dec_pix,
                            P_image=False, fix_additional_comps=fix_additional_comps,
                            manual_rms_region=manual_rms_region, src_name=src_name)
            CHAN_U_imfit    = get_imfit_values(est_U, CHAN_images[3], CHAN_P_ra_pix, CHAN_P_dec_pix)
            CHAN_U_imfit, _ = check_fit_offset_and_refix(
                CHAN_U_imfit, CHAN_images[3], CHAN_P_ra_pix, CHAN_P_dec_pix,
                bmaj, src_name, 'U', manual_rms_region=manual_rms_region, tag=tag)

            # Fit Stokes V anchored to channel I. FORCE_FIX_STOKES_V=True (global
            # default) passes force_fix_all=True, bypassing the S/N check for
            # all components including k=0.
            est_V = f'estimate_V_{src_name}_{tag}.txt'
            check_position(est_V, CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix,
                            P_image=False, fix_additional_comps=fix_additional_comps,
                            manual_rms_region=manual_rms_region, src_name=src_name,
                            force_fix_all=FORCE_FIX_STOKES_V)
            CHAN_V_imfit, _ = check_fit_offset_and_refix(
                get_imfit_values(est_V, CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix),
                CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix,
                bmaj, src_name, 'V', manual_rms_region=manual_rms_region, tag=tag)

        result = {'freq_GHz': freq_GHz, 'bmaj_asec': bmaj, 'bmin_asec': bmin,
                  'bpa_deg': bpa, 'components': {}}

        for component in components:
            flux_I    = CHAN_I_imfit['results'][component]['peak']['value'] * 1e3
            err_I     = CHAN_I_imfit['results'][component]['peak']['error'] * 1e3
            RA_I      = CHAN_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_I     = CHAN_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            RA_pix_I  = CHAN_I_imfit['results'][component]['pixelcoords'][0]
            DEC_pix_I = CHAN_I_imfit['results'][component]['pixelcoords'][1]
            rms_I     = get_imstat_values(CHAN_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3

            comp_result = {
                'I_flux_mJy': flux_I, 'I_err_mJy': err_I, 'I_rms_mJy': rms_I,
                'I_RA_deg': RA_I, 'I_DEC_deg': DEC_I,
            }

            if not only_intensity:
                flux_P = CHAN_P_imfit['results'][component]['peak']['value'] * 1e3
                flux_Q = CHAN_Q_imfit['results'][component]['peak']['value'] * 1e3
                flux_U = CHAN_U_imfit['results'][component]['peak']['value'] * 1e3
                flux_V = CHAN_V_imfit['results'][component]['peak']['value'] * 1e3

                err_P = CHAN_P_imfit['results'][component]['peak']['error'] * 1e3
                err_Q = CHAN_Q_imfit['results'][component]['peak']['error'] * 1e3
                err_U = CHAN_U_imfit['results'][component]['peak']['error'] * 1e3
                err_V = CHAN_V_imfit['results'][component]['peak']['error'] * 1e3

                RA_P  = CHAN_P_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
                DEC_P = CHAN_P_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi

                rms_Q = get_imstat_values(CHAN_images[2], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
                rms_U = get_imstat_values(CHAN_images[3], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
                rms_V = get_imstat_values(CHAN_images[4], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3

                flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U, rms_V, pol_flag, Aq=0.8)

                LP_frac     = flux_P0 / flux_I * 100.0
                LP_frac_err = abs(LP_frac) * np.sqrt((rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2)

                if pol_flag:
                    LP_EVPA     = np.arctan2(flux_U, flux_Q) * 180.0 / np.pi * 0.5
                    LP_EVPA_err = np.sqrt(flux_U ** 2 * rms_Q ** 2 + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2 + flux_Q ** 2) * 180.0 / np.pi * 0.5
                else:
                    LP_EVPA     = None
                    LP_EVPA_err = None

                comp_result.update({
                    'Q_flux_mJy': flux_Q, 'U_flux_mJy': flux_U, 'V_flux_mJy': flux_V,
                    'P_flux_mJy': flux_P, 'P0_flux_mJy': flux_P0,
                    'Q_err_mJy': err_Q, 'U_err_mJy': err_U, 'V_err_mJy': err_V, 'P_err_mJy': err_P,
                    'Q_rms_mJy': rms_Q, 'U_rms_mJy': rms_U, 'V_rms_mJy': rms_V, 'P_rms_mJy': rms_P,
                    'P_RA_deg': RA_P, 'P_DEC_deg': DEC_P,
                    'LP_frac': LP_frac, 'LP_frac_err': LP_frac_err,
                    'LP_EVPA': LP_EVPA, 'LP_EVPA_err': LP_EVPA_err,
                })

            result['components'][component] = comp_result

        return chan_idx, result

    except Exception:
        msg(f'Fitting Failed: Channel {chan_idx} is likely flagged')
        return chan_idx, None


def fit_channel_group_batch(jobs):
    '''
    Pool entry point: processes one whole chunk of channel-groups in a
    spawned worker process. imfit/imstat/imhead/tb (top-level imports) resolve
    normally there since spawn re-executes the module's top-level code.
    '''
    return [fit_channel_group(job) for job in jobs]


def _try_get_ms_timing(src_name, n_intervals):
    """
    Derive per-interval MJD timing from the measurement set for src_name.

    Tries project_info['target_ms'] first (the per-target split MS).  Falls
    back to project_info['working_ms'] (filtered by FIELD_ID when possible).
    Returns None and lets the caller fall back to the DATE-OBS image header
    if neither MS is found or readable.

    Returns a list of n_intervals dicts:
        'time_ctr_mjd' : float  -- interval centre in MJD
        'time_dt'      : float  -- relevant time span in seconds
    """
    if not os.path.isfile('project_info.json'):
        msg('  [timing] project_info.json not found -- will use image header fallback')
        return None

    with open('project_info.json') as _f:
        _pinfo = json.load(_f)

    # -----------------------------------------------------------------------
    # Locate the MS: try the per-target split MS, then the full working MS
    # -----------------------------------------------------------------------
    _ms_to_use  = None
    _field_id   = None

    _target_names = _pinfo.get('target_names', [])
    _target_ms    = _pinfo.get('target_ms',    [])

    if src_name in _target_names:
        _idx = _target_names.index(src_name)
        _candidate = _target_ms[_idx] if _idx < len(_target_ms) else None
        if _candidate and os.path.exists(_candidate):
            _ms_to_use = _candidate
            msg(f'  [timing] using split MS: {os.path.basename(_ms_to_use)}')

    if _ms_to_use is None:
        _working = _pinfo.get('working_ms', '')
        if _working and os.path.exists(_working):
            _ms_to_use = _working
            _target_ids = _pinfo.get('target_ids', [])
            if src_name in _target_names:
                _idx = _target_names.index(src_name)
                if _idx < len(_target_ids):
                    _field_id = int(_target_ids[_idx])
            msg(f'  [timing] split MS absent; using working MS: {os.path.basename(_ms_to_use)}'
                + (f' (FIELD_ID={_field_id})' if _field_id is not None else ''))

    if _ms_to_use is None:
        msg(f'  [timing] no MS found for "{src_name}" -- will use image header fallback')
        return None

    # -----------------------------------------------------------------------
    # Read the TIME column
    # -----------------------------------------------------------------------
    try:
        tb.open(_ms_to_use)
        _all_times  = tb.getcol('TIME')
        _all_fields = tb.getcol('FIELD_ID') if _field_id is not None else None
        tb.close()
        if _field_id is not None:
            _times = np.unique(_all_times[_all_fields == _field_id])
        else:
            _times = np.unique(_all_times)
    except Exception as _e:
        msg(f'  [timing] failed to read MS ({_e}) -- will use image header fallback')
        return None

    if len(_times) == 0:
        msg('  [timing] no timestamps found -- will use image header fallback')
        return None

    _times   = np.sort(_times)
    _t_start = float(_times[0])
    _t_end   = float(_times[-1])
    _t_span  = _t_end - _t_start
    msg(f'  [timing] {len(_times)} unique timestamps, '
        f'span = {_t_span:.1f} s ({_t_span / 3600.0:.3f} h)')

    # -----------------------------------------------------------------------
    # Compute timing per interval
    # -----------------------------------------------------------------------
    if n_intervals == 1:
        _time_ctr_mjd = (_t_start + _t_end) / 2.0 / 86400.0
        _time_dt      = _t_span
        msg(f'  [timing] MFS: centre = {_time_ctr_mjd:.8f} MJD, delta = {_time_dt:.0f} s')
        return [{'time_ctr_mjd': _time_ctr_mjd, 'time_dt': _time_dt}]

    # Time-resolved: split into n_intervals evenly-sized chunks.
    # Each interval centre = start of chunk + half the median inter-chunk step.
    _n      = len(_times)
    _chunk  = max(1, _n // n_intervals)
    _starts = np.array([float(_times[min(i * _chunk, _n - 1)]) for i in range(n_intervals)])
    _diffs  = np.diff(_starts)
    _sep_s  = float(np.median(_diffs)) if len(_diffs) > 0 else _t_span
    _total  = float(_starts[-1] - _starts[0])

    msg(f'  [timing] {n_intervals} intervals, median step = {_sep_s:.1f} s, '
        f'total span = {_total:.1f} s')

    _result = []
    for _i, _t0 in enumerate(_starts):
        _cmjd = (_t0 + _sep_s / 2.0) / 86400.0
        _result.append({'time_ctr_mjd': _cmjd, 'time_dt': _total})
        msg(f'  [timing]   interval {_i:>3}: start = {_t0 / 86400.0:.8f} MJD, '
            f'centre = {_cmjd:.8f} MJD')

    return _result


def extract_polarization_properties(src_name,
    src_im_identifier,
    src_im_suffix,
    src_ra,
    src_dec,
    src_ulims,
    pol_flag,
    manual_rms_region,
    image_directory,
    image_identifier,
    fix_additional_comps = False,
    out_label = ''):

    '''
    Core extraction routine: fit Stokes I, Q, U, V (and polarized intensity P)
    for all source components across both the multi-frequency-synthesis (MFS)
    and per-channel images produced by WSCLEAN.

    For each prefix (one per time interval, or just one for a full-observation
    image), the function:
      1. Fits the MFS Stokes I image to establish the reference sky position.
      2. Fits MFS P, Q, U, V using that position as a prior, with post-fit
         offset checks that refix wandering components automatically.
      3. Iterates over per-channel image sets, repeating the same fitting
         sequence and accumulating results for RM synthesis output.

    Supports multi-component sources (e.g. partially resolved jet ejecta): set
    fix_additional_comps=True to anchor secondary components to the Stokes I
    position in non-intensity images, preventing noise-driven wandering.

    Parameters
    ----------
    src_name            -- source name string; drives image discovery (glob pattern
                           matching, MS lookup for timing) and is the default for
                           output filenames/titles when out_label is not set
    src_im_identifier   -- glob pattern identifying the image set, constructed
                           from rmsynth_info.json fields
    src_im_suffix       -- FITS suffix, e.g. 'image.fits' or
                           'image.homogenized.fits'; auto-detected if MFS
                           images are not found with the requested suffix
    src_ra              -- array of RA guesses in decimal degrees, one per component
    src_dec             -- array of Dec guesses in decimal degrees, one per component
    src_ulims           -- list of booleans, True = treat component as upper limit
                           (position held fixed, flux not de-biased)
    pol_flag            -- True if a polarization angle calibrator was used;
                           controls which P image type is fitted (Plin vs Ptot)
                           and which de-biasing formula is applied in calculate_P0
    manual_rms_region   -- CASA region string for the noise measurement box/annulus;
                           if False, a default annulus of ~500 beam areas is used
    image_directory     -- top-level image directory (e.g. 'IMAGES' or 'INTERVALS')
    image_identifier    -- image identifier string (e.g. 'pcalmask', 'restored')
    fix_additional_comps-- if True, secondary components (k>0) have their positions
                           frozen in all non-Stokes-I images; recommended for faint
                           partially unresolved jet ejecta
    out_label           -- optional display/output name override (rmsynth_info.json
                           'label' field). Empty string (default) means "use
                           src_name" -- current behaviour, unchanged. When set, it
                           replaces src_name in output_dictionary['name'], the
                           output JSON filename, and every plot filename/title
                           produced for this source. src_name itself keeps doing
                           image discovery and MS timing lookups either way, so
                           out_label never affects which input files are found.

    Returns
    -------
    output_dictionary : dict
        Nested dictionary with 'MFS' and 'CHAN' sub-dictionaries, each keyed by
        component name, containing fluxes, errors, RMS values, positions, beam
        parameters, LP fraction, EVPA, and observation metadata.
    '''

    # ------------------------------------------------------------------
    # Log key settings at the start so the output log is self-documenting
    # ------------------------------------------------------------------
    msg('')
    msg(f'{"="*70}')
    msg(f'  Starting polarization extraction for source: {src_name}')
    msg(f'  Image identifier : {src_im_identifier}')
    if manual_rms_region:
        msg(f'  RMS region       : MANUAL -- {manual_rms_region}')
    else:
        msg(f'  RMS region       : default annulus (~500 beam areas) centred on source')
    msg(f'  Polarization angle calibrator present: {pol_flag}')
    msg(f'  Fix secondary component positions    : {fix_additional_comps}')
    msg(f'{"="*70}')
    msg('')

    # Display/output name: out_label overrides src_name for everything written
    # out below (JSON filename+contents, plot filenames+titles). src_name keeps
    # doing image discovery/MS-timing lookups untouched -- see out_label in the
    # docstring above.
    label = out_label if out_label else src_name
    if label != src_name:
        msg(f'  Output label     : "{label}" (overrides source name "{src_name}" in output filenames)')

    # Initialize the output dictionary.  'MFS' holds one entry per time
    # interval; 'CHAN' holds per-frequency-channel data for RM synthesis.
    output_dictionary = {'name': label}
    output_dictionary['MFS']  = {}
    output_dictionary['CHAN'] = {}

    # ------------------------------------------------------------------
    # Locate images and determine the correct suffix to use
    # ------------------------------------------------------------------

    # Confirm at least some images exist before proceeding
    test_images = glob.glob(f'{src_im_identifier}*{src_im_suffix}')
    if not test_images:
        raise FileNotFoundError(
            f'No images found matching: {src_im_identifier}*{src_im_suffix}\n'
            f'Check the image_directory and image_identifier fields in rmsynth_info.json.')
    msg(f'Found images matching pattern (example): {test_images[0]}')

    # MFS images may use 'image.homogenized.fits' if per-channel splitting
    # has already consumed the standard image.fits products
    mfs_im_suffix = src_im_suffix[:]
    if not glob.glob(f'{src_im_identifier}*MFS*{src_im_suffix}'):
        msg(f'WARNING: No MFS {src_im_suffix} images found; '
            f'falling back to image.homogenized.fits')
        mfs_im_suffix = 'image.homogenized.fits'
    else:
        msg(f'MFS image suffix: {mfs_im_suffix}')

    # Check whether polarization images (Q, U, V, P) exist, or if this run
    # is Stokes I only (e.g. early-stage pipeline without polarization calibration)
    if not glob.glob(f'{src_im_identifier}*MFS-I-{mfs_im_suffix}'):
        only_intensity = True
        msg('No MFS polarization images found -- running in Stokes I only mode')
    else:
        only_intensity = False
        msg('MFS polarization images found -- fitting full Stokes IQUV')

    # Build the list of unique image prefixes.  Each prefix corresponds to one
    # time interval (or the full observation if image_timing=false in rmsynth_info.json).
    _mfs_glob_all = glob.glob(f'{src_im_identifier}*MFS*{mfs_im_suffix}')
    prefix_arr = sorted(list(set([x.split('-MFS')[0] for x in _mfs_glob_all])))
    is_time_resolved = len(prefix_arr) > 1
    msg(f'Found {len(prefix_arr)} prefix(es) (time interval(s)) to process')

    # Cache each prefix's MFS images from the single glob above, in memory,
    # instead of re-globbing the same (potentially huge) directory once per
    # prefix -- with hundreds of timesteps this turns hundreds of repeated
    # full-directory scans into zero additional filesystem calls.
    _mfs_images_by_prefix = {}
    for _f in _mfs_glob_all:
        _mfs_images_by_prefix.setdefault(_f.split('-MFS')[0], []).append(_f)
    msg(f'Cached MFS images for {len(_mfs_images_by_prefix)} prefix(es) from one directory scan.')

    # One glob for CHAN images across all prefixes, partitioned in memory by
    # prefix -- channel counts vary per timestep (flagging differs per
    # snapshot), so each prefix needs its own real result, not a suffix set
    # learned from prefix 0. Filtered by a '-MFS-' substring check rather than
    # a [!MFS] glob bracket, since fnmatch's '*' backtracking can defeat that
    # bracket and falsely match an MFS image.
    _all_chan_glob = glob.glob(f'{src_im_identifier}*-{src_im_suffix}')
    _chan_images_by_prefix = {}
    for _f in _all_chan_glob:
        if '-MFS-' in _f:
            continue
        # CHAN filenames are '{prefix}-{chan}-{stokes}-{suffix}' and chan/
        # stokes/suffix are all hyphen-free, so the prefix is always
        # everything before the last 3 hyphens.
        _prefix_key = _f.rsplit('-', 3)[0]
        _chan_images_by_prefix.setdefault(_prefix_key, []).append(_f)

    _n_chan_total = sum(len(v) for v in _chan_images_by_prefix.values())
    _n_stokes_per_chan = 1 if only_intensity else 5
    msg(f'Cached CHAN images for {len(_chan_images_by_prefix)} prefix(es) '
        f'({_n_chan_total} image(s) total, ~{_n_chan_total // max(1, len(_chan_images_by_prefix)) // max(1, _n_stokes_per_chan)} '
        f'channel(s)/prefix on average) from one directory scan.')

    # Pre-compute timing before the prefix loop.
    #   Time-resolved: read DATE-OBS from each prefix's MFS image header; the
    #     centre is the header time shifted forward by half the median inter-
    #     prefix step, and the span covers first-to-last prefix.
    #   MFS (single epoch): try the target MS for accurate start/end; fall
    #     back to DATE-OBS + 15-min default inside the loop if the MS is gone.
    if is_time_resolved:
        _ms_timings = None
        _hdr_mjds   = []
        for _pfx in prefix_arr:
            _pfx_imgs = sorted(_mfs_images_by_prefix.get(_pfx, []))
            _pfx_dobs = imhead(_pfx_imgs[0], mode='get', hdkey='DATE-OBS').replace('/','-',2).replace('/','T')
            _hdr_mjds.append(_iso_to_mjd(_pfx_dobs))
        _hdr_mjds   = np.array(_hdr_mjds)
        _hdr_diffs  = np.diff(_hdr_mjds) * 86400.0         # seconds between adjacent prefixes
        _hdr_step_s = float(np.median(_hdr_diffs)) if len(_hdr_diffs) > 0 else 0.0
        _hdr_span_s = float((_hdr_mjds[-1] - _hdr_mjds[0]) * 86400.0)
        msg(f'  [timing] time-resolved ({len(prefix_arr)} prefixes): '
            f'median step = {_hdr_step_s:.1f} s, total span = {_hdr_span_s:.1f} s')
    else:
        _ms_timings  = _try_get_ms_timing(src_name, 1)
        _hdr_mjds = _hdr_step_s = _hdr_span_s = None

    # Compute the output JSON filename from the prefix list so we can
    # check for an existing file when OVERWRITE is False.
    if len(prefix_arr) > 0:
        _last_prefix = prefix_arr[-1]
        _split_id    = image_identifier.rstrip('-')
        if len(_last_prefix.split(f'{_split_id}-t')) > 1:
            json_basename = _last_prefix.split(f'{_split_id}-t')[0] + _split_id
        else:
            json_basename = _last_prefix
        # Relocate into RESULTS/ first, while json_basename is still the
        # untouched on-disk path -- image_directory is only guaranteed to
        # appear verbatim in *that*. Applying the label swap first would risk
        # corrupting a directory segment that happens to contain src_name too
        # (e.g. image_directory or CWD named after the source), which then
        # silently breaks this relocation and writes the file under the
        # original image directory instead of RESULTS. See _apply_label().
        json_basename = json_basename.replace(image_directory, 'RESULTS')
        # Now the label override, restricted to the filename component only.
        json_basename = _apply_label(json_basename, src_name, label)
        json_path = '{}_polarization.json'.format(json_basename)
    else:
        json_path = None

    # --- If OVERWRITE is False and the JSON exists, skip fitting and just plot ---
    if not OVERWRITE and json_path is not None and os.path.isfile(json_path):
        msg(f'OVERWRITE is False and {json_path} exists -- skipping fitting, loading for plots')
        with open(json_path, 'r') as j:
            output_dictionary = json.load(j)

        components = sorted([key for key in output_dictionary['MFS'].keys()])

        if is_time_resolved:
            for component in components:
                _plot_detected_spectra(output_dictionary, src_name, prefix_arr, component, label,
                                       pol_flag=pol_flag)
            plot_light_curve(output_dictionary, src_name, prefix_arr, label=label)
        else:
            for component in components:
                plot_stokes_spectrum(output_dictionary, src_name, prefix_arr[0],
                                    component, 0, save_plot=True, label=label,
                                    pol_flag=pol_flag)
        return 0

    # Persistent pool for CHAN-image fitting, created lazily on first use (once
    # we've actually seen a channel image to size workers against) and reused
    # across every prefix/timestep -- MFS fitting stays serial in this process
    # since it's cheap; per-channel fitting is the part worth parallelizing.
    #
    # Recycled every chan_pool_recycle_interval prefixes (computed by
    # compute_max_workers, scaled to available headroom/worker -- see there):
    # casacore's table/image caching does not fully release memory between
    # imfit/imstat calls within a long-lived worker process, so per-worker
    # memory climbs several MB/prefix. Shutting down and recreating the pool
    # periodically bounds that growth instead of letting it accumulate for
    # the whole run; recreation cost is small (~1-3s, paid in parallel across
    # workers) relative to the risk of an OOM kill partway through a run with
    # hundreds of timesteps.
    chan_pool = None
    chan_pool_max_workers = None
    chan_pool_ctx = mp.get_context('spawn')
    chan_pool_prefixes_since_recycle = 0
    chan_pool_recycle_interval = 50   # placeholder until compute_max_workers sets it

    for k, prefix in enumerate(prefix_arr[:]):

        msg('')
        msg(f'--- Prefix {k+1}/{len(prefix_arr)}: {os.path.basename(prefix)} ---')

        # Collect MFS images for this prefix. sorted() puts them in the order
        # I, P(lin/tot), Q, U, V which is the order the rest of the code assumes.
        if only_intensity:
            MFS_images = [f'{prefix}-MFS-{mfs_im_suffix}']
        else:
            MFS_images = _mfs_images_by_prefix.get(prefix, [])
            # Exclude whichever polarized intensity type is not in use:
            # pol_flag=True  --> pol angle calibrator present --> use Plin (sqrt(Q^2+U^2))
            # pol_flag=False --> no pol angle calibrator       --> use Ptot (sqrt(Q^2+U^2+V^2))
            if pol_flag:
                MFS_images = sorted([im for im in MFS_images if '-Ptot-' not in im])
            else:
                MFS_images = sorted([im for im in MFS_images if '-Plin-' not in im])

        msg(f'MFS images selected ({len(MFS_images)}):')
        for MFS_image in MFS_images:
            msg(f'  {os.path.basename(MFS_image)}')

        # Read observation metadata from the Stokes I header
        freq_GHz = imhead(MFS_images[0], mode='get', hdkey='CRVAL3')['value'] / 1.0e9
        date_obs = imhead(MFS_images[0], mode='get', hdkey='DATE-OBS').replace('/','-',2).replace('/','T')
        bmaj     = imhead(MFS_images[0], mode='get', hdkey='bmaj')['value']
        bmin     = imhead(MFS_images[0], mode='get', hdkey='bmin')['value']
        bpa      = imhead(MFS_images[0], mode='get', hdkey='bpa')['value']
        msg(f'Observation: date={date_obs}  freq={freq_GHz:.4f} GHz  '
            f'beam={bmaj:.2f}"x{bmin:.2f}" PA={bpa:.1f} deg')

        # Per-prefix MJD timing -- time_ctr_mjd/time_ctr_isot/time_dt always
        # describe the centre of the relevant exposure, never its start.
        if is_time_resolved:
            # Centre = header DATE-OBS shifted by half the median inter-prefix step.
            _time_ctr_mjd = float(_hdr_mjds[k]) + _hdr_step_s / 2.0 / 86400.0
            _time_dt      = _hdr_span_s
            msg(f'  Timing (image header): centre MJD = {_time_ctr_mjd:.8f}, delta = {_time_dt:.0f} s')
        elif _ms_timings is not None:
            _time_ctr_mjd = _ms_timings[0]['time_ctr_mjd']
            _time_dt      = _ms_timings[0]['time_dt']
            msg(f'  Timing (MS): centre MJD = {_time_ctr_mjd:.8f}, delta = {_time_dt:.0f} s')
        else:
            # MFS fallback: DATE-OBS from image header + 15-min default width.
            _time_ctr_mjd = _iso_to_mjd(date_obs)
            _time_dt      = 15.0 * 60.0
            msg(f'  Timing (header fallback): centre MJD = {_time_ctr_mjd:.8f}, delta = 900 s (15 min default)')
        _time_ctr_isot = _mjd_to_isot(_time_ctr_mjd)

        # Convert the input RA/Dec guesses (decimal degrees) to pixel coordinates
        # in this specific image.  imstat finds the nearest pixel to the input position.
        src_ra_pix  = []
        src_dec_pix = []
        for z in range(src_ra.size):
            region = 'circle[[{}deg,{}deg],1.0pix]'.format(src_ra[z], src_dec[z])
            pixel_imstat = imstat(MFS_images[0], region=region)
            src_ra_pix.append(pixel_imstat['maxpos'][0])
            src_dec_pix.append(pixel_imstat['maxpos'][1])

        # -----------------------------------------------------------------
        # Fit MFS Stokes I -- establishes the reference position for all
        # subsequent polarization fits in this prefix epoch
        # -----------------------------------------------------------------
        # Upper-limit components have their position fixed (xyabp); detected
        # components are fitted freely for flux, size and PA (abp).
        fix = []
        for z, boolean in enumerate(src_ulims):
            if boolean:
                fix.append('xyabp')
                msg(f'  Component {z}: upper limit -- position will be fixed')
            else:
                fix.append('abp')
                msg(f'  Component {z}: detected -- position will be fitted freely')

        make_estimate(f'estimate_I_{src_name}.txt', MFS_images[0], src_ra_pix, src_dec_pix, fix)
        MFS_I_imfit = get_imfit_values(f'estimate_I_{src_name}.txt', MFS_images[0], src_ra_pix, src_dec_pix)

        # Extract the list of fitted component keys (e.g. ['component0', 'component1'])
        components    = [key for key in MFS_I_imfit['results'].keys() if 'component' in key]
        MFS_I_ra_pix  = [MFS_I_imfit['results'][key]['pixelcoords'][0] for key in components]
        MFS_I_dec_pix = [MFS_I_imfit['results'][key]['pixelcoords'][1] for key in components]

        # --- MFS Stokes I drift check vs input sky-coordinate seed -----------
        # If MAX_I_DRIFT_PIX is set and any component has moved further than
        # that threshold from the input reference, refix all positions and refit.
        if MAX_I_DRIFT_PIX is not False:
            mfs_i_drift_exceeded = False
            for z, component in enumerate(components):
                drift = np.sqrt((MFS_I_ra_pix[z]  - src_ra_pix[z])  ** 2 +
                                (MFS_I_dec_pix[z] - src_dec_pix[z]) ** 2)
                if drift > MAX_I_DRIFT_PIX:
                    msg(f'  WARNING: MFS Stokes I component {z} drifted {drift:.1f} pix '
                        f'from input reference (threshold {MAX_I_DRIFT_PIX} pix). '
                        f'Will re-fit with all positions fixed to input reference.')
                    mfs_i_drift_exceeded = True
                else:
                    msg(f'  MFS Stokes I component {z} drift OK: {drift:.1f} pix '
                        f'(threshold {MAX_I_DRIFT_PIX} pix)')
            if mfs_i_drift_exceeded:
                fix_all = ['xyabp'] * len(src_ra_pix)
                make_estimate(f'estimate_I_{src_name}.txt', MFS_images[0],
                              src_ra_pix, src_dec_pix, fix_all)
                MFS_I_imfit   = get_imfit_values(f'estimate_I_{src_name}.txt',
                                                  MFS_images[0], src_ra_pix, src_dec_pix)
                MFS_I_ra_pix  = [MFS_I_imfit['results'][key]['pixelcoords'][0] for key in components]
                MFS_I_dec_pix = [MFS_I_imfit['results'][key]['pixelcoords'][1] for key in components]
                msg(f'  MFS Stokes I re-fit complete: all positions fixed to input reference.')

        # Log the Stokes I fitted results -- this position is the astrometric
        # reference for all polarization fits below
        msg('')
        msg(f'  --- MFS Stokes I fit results (astrometric reference) ---')
        for z, component in enumerate(components):
            flux_I_log = MFS_I_imfit['results'][component]['peak']['value'] * 1e3
            err_I_log  = MFS_I_imfit['results'][component]['peak']['error'] * 1e3
            RA_I_log   = MFS_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_I_log  = MFS_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            _I_log_stats = get_imstat_values(MFS_images[0], MFS_I_ra_pix[z], MFS_I_dec_pix[z], manual_rms_region)
            rms_I_log  = _I_log_stats[3] * 1e3
            rms_region_I_log = _I_log_stats[4]
            msg(f'  {component}: RA={RA_I_log:.6f} deg  Dec={DEC_I_log:.6f} deg  '
                f'({MFS_I_ra_pix[z]:.1f}, {MFS_I_dec_pix[z]:.1f}) pix')
            msg(f'    Flux = {flux_I_log:.3f} +/- {err_I_log:.3f} mJy  '
                f'rms = {rms_I_log:.3f} mJy  S/N = {flux_I_log/rms_I_log:.1f}')
            msg(f'    RMS region: {rms_region_I_log}')
        msg('')

        # Populate the output dictionary with MFS Stokes I results
        for z, component in enumerate(components):

            flux_I   = MFS_I_imfit['results'][component]['peak']['value'] * 1e3
            err_I    = MFS_I_imfit['results'][component]['peak']['error'] * 1e3
            RA_I     = MFS_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_I    = MFS_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            RA_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][0]
            DEC_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][1]
            rms_I    = get_imstat_values(MFS_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3

            if component not in output_dictionary['MFS']:
                # First epoch: initialise lists
                output_dictionary['MFS'][component] = {}
                output_dictionary['MFS'][component]['upperlimit']  = [src_ulims[z]]
                output_dictionary['MFS'][component]['freq_GHz']    = [freq_GHz]
                output_dictionary['MFS'][component]['bmaj_asec']   = [bmaj]
                output_dictionary['MFS'][component]['bmin_asec']   = [bmin]
                output_dictionary['MFS'][component]['bpa_deg']     = [bpa]
                output_dictionary['MFS'][component]['I_flux_mJy']   = [flux_I]
                output_dictionary['MFS'][component]['I_err_mJy']    = [err_I]
                output_dictionary['MFS'][component]['I_rms_mJy']    = [rms_I]
                output_dictionary['MFS'][component]['I_RA_deg']     = [RA_I]
                output_dictionary['MFS'][component]['I_DEC_deg']    = [DEC_I]
                output_dictionary['MFS'][component]['time_ctr_mjd']  = [_time_ctr_mjd]
                output_dictionary['MFS'][component]['time_ctr_isot'] = [_time_ctr_isot]
                output_dictionary['MFS'][component]['time_dt']       = [_time_dt]
            else:
                # Subsequent epochs: append to existing lists
                output_dictionary['MFS'][component]['freq_GHz'].append(freq_GHz)
                output_dictionary['MFS'][component]['bmaj_asec'].append(bmaj)
                output_dictionary['MFS'][component]['bmin_asec'].append(bmin)
                output_dictionary['MFS'][component]['bpa_deg'].append(bpa)
                output_dictionary['MFS'][component]['I_flux_mJy'].append(flux_I)
                output_dictionary['MFS'][component]['I_err_mJy'].append(err_I)
                output_dictionary['MFS'][component]['I_rms_mJy'].append(rms_I)
                output_dictionary['MFS'][component]['I_RA_deg'].append(RA_I)
                output_dictionary['MFS'][component]['I_DEC_deg'].append(DEC_I)
                output_dictionary['MFS'][component]['time_ctr_mjd'].append(_time_ctr_mjd)
                output_dictionary['MFS'][component]['time_ctr_isot'].append(_time_ctr_isot)
                output_dictionary['MFS'][component]['time_dt'].append(_time_dt)
            
        
        # -----------------------------------------------------------------
        # Fit MFS polarization images (P, Q, U, V)
        # Each fit is followed by an offset check: if any component has
        # drifted more than 1/3 BMAJ from the Stokes I reference, the
        # fit is automatically re-run with all positions fixed.
        # -----------------------------------------------------------------
        if not only_intensity:
            
            # --- First pass: fit P freely to see if it lands near Stokes I ---
            # If it wanders too far (handled by check_fit_offset_and_refix later)
            # the second pass below will re-fit with positions fixed.
            try:
                check_position(f'estimate_P_{src_name}.txt', MFS_images[1],
                                MFS_I_ra_pix, MFS_I_dec_pix,
                                P_image=True, fix_additional_comps=False,
                                manual_rms_region=manual_rms_region, src_name=src_name)
                MFS_P_imfit   = get_imfit_values(f'estimate_P_{src_name}.txt',
                                                  MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
                MFS_P_ra_pix  = [MFS_P_imfit['results'][key]['pixelcoords'][0] for key in components]
                MFS_P_dec_pix = [MFS_P_imfit['results'][key]['pixelcoords'][1] for key in components]

                # Check whether any P component has drifted more than 1/3 BMAJ from Stokes I.
                # If so, fall through to the guarded re-fit below.
                dist_check = False
                dist_max   = 0.0
                for z in range(len(MFS_P_ra_pix)):
                    dist = np.sqrt((MFS_I_ra_pix[z] - MFS_P_ra_pix[z])**2 +
                                   (MFS_I_dec_pix[z] - MFS_P_dec_pix[z])**2)
                    if dist > 3.0:   # >3 pixel offset is unphysical for unresolved sources
                        dist_check = True
                        dist_max   = max(dist_max, dist)

            except Exception as e:
                msg(f'  P free-fit failed (likely complex multi-component field); '
                    f'will re-fit with secondary positions fixed to Stokes I: {e}')
                dist_max   = -1.0
                dist_check = True

            if dist_check:
                msg(f'  Stokes I--P offset = {dist_max:.1f} pix; re-fitting P with '
                    f'secondary component positions fixed to Stokes I')
                check_position(f'estimate_P_{src_name}.txt', MFS_images[1],
                                MFS_I_ra_pix, MFS_I_dec_pix,
                                P_image=True, fix_additional_comps=fix_additional_comps,
                                manual_rms_region=manual_rms_region, src_name=src_name)
                MFS_P_imfit   = get_imfit_values(f'estimate_P_{src_name}.txt',
                                                  MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
                MFS_P_ra_pix  = [MFS_P_imfit['results'][key]['pixelcoords'][0] for key in components]
                MFS_P_dec_pix = [MFS_P_imfit['results'][key]['pixelcoords'][1] for key in components]

            # Post-fit offset check for P -- refixes automatically if needed
            MFS_P_imfit, _ = check_fit_offset_and_refix(
                MFS_P_imfit, MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix,
                bmaj, src_name, 'P', manual_rms_region=manual_rms_region)
            MFS_P_ra_pix  = [MFS_P_imfit['results'][key]['pixelcoords'][0] for key in components]
            MFS_P_dec_pix = [MFS_P_imfit['results'][key]['pixelcoords'][1] for key in components]

            # --- Fit Q and U, referenced to the P position ------------------
            # Using the P position as the reference for Q/U is more accurate
            # than using Stokes I when the source is significantly polarized.
            check_position(f'estimate_Q_{src_name}.txt', MFS_images[2],
                            MFS_P_ra_pix, MFS_P_dec_pix,
                            P_image=False, fix_additional_comps=fix_additional_comps,
                            manual_rms_region=manual_rms_region, src_name=src_name)
            MFS_Q_imfit   = get_imfit_values(f'estimate_Q_{src_name}.txt',
                                              MFS_images[2], MFS_P_ra_pix, MFS_P_dec_pix)
            MFS_Q_imfit, _ = check_fit_offset_and_refix(
                MFS_Q_imfit, MFS_images[2], MFS_P_ra_pix, MFS_P_dec_pix,
                bmaj, src_name, 'Q', manual_rms_region=manual_rms_region)
            MFS_Q_ra_pix  = [MFS_Q_imfit['results'][key]['pixelcoords'][0] for key in components]
            MFS_Q_dec_pix = [MFS_Q_imfit['results'][key]['pixelcoords'][1] for key in components]

            check_position(f'estimate_U_{src_name}.txt', MFS_images[3],
                            MFS_P_ra_pix, MFS_P_dec_pix,
                            P_image=False, fix_additional_comps=fix_additional_comps,
                            manual_rms_region=manual_rms_region, src_name=src_name)
            MFS_U_imfit   = get_imfit_values(f'estimate_U_{src_name}.txt',
                                              MFS_images[3], MFS_P_ra_pix, MFS_P_dec_pix)
            MFS_U_imfit, _ = check_fit_offset_and_refix(
                MFS_U_imfit, MFS_images[3], MFS_P_ra_pix, MFS_P_dec_pix,
                bmaj, src_name, 'U', manual_rms_region=manual_rms_region)
            MFS_U_ra_pix  = [MFS_U_imfit['results'][key]['pixelcoords'][0] for key in components]
            MFS_U_dec_pix = [MFS_U_imfit['results'][key]['pixelcoords'][1] for key in components]

            # --- Fit Stokes V, referenced back to Stokes I ------------------
            # Circular polarization is intrinsically weak for synchrotron
            # sources and residual I->V leakage can dominate; anchoring to I
            # avoids noise-driven position bias.  When FORCE_FIX_STOKES_V
            # is True (the global default), ALL components are held fixed via
            # force_fix_all, bypassing the S/N check entirely.
            check_position(f'estimate_V_{src_name}.txt', MFS_images[4],
                            MFS_I_ra_pix, MFS_I_dec_pix,
                            P_image=False, fix_additional_comps=fix_additional_comps,
                            manual_rms_region=manual_rms_region, src_name=src_name,
                            force_fix_all=FORCE_FIX_STOKES_V)
            MFS_V_imfit   = get_imfit_values(f'estimate_V_{src_name}.txt',
                                              MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix)
            MFS_V_imfit, _ = check_fit_offset_and_refix(
                MFS_V_imfit, MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix,
                bmaj, src_name, 'V', manual_rms_region=manual_rms_region)
            MFS_V_ra_pix  = [MFS_V_imfit['results'][key]['pixelcoords'][0] for key in components]
            MFS_V_dec_pix = [MFS_V_imfit['results'][key]['pixelcoords'][1] for key in components]   
  
            # Once again initialize arrays that don't exist
            for component in components[:]:
            
                # For ease and readability separate out the parameters for calculations and rms extraction
                flux_I = MFS_I_imfit['results'][component]['peak']['value'] * 1e3
                flux_P = MFS_P_imfit['results'][component]['peak']['value'] * 1e3
                flux_Q = MFS_Q_imfit['results'][component]['peak']['value'] * 1e3
                flux_U = MFS_U_imfit['results'][component]['peak']['value'] * 1e3
                flux_V = MFS_V_imfit['results'][component]['peak']['value'] * 1e3
                
                err_P = MFS_P_imfit['results'][component]['peak']['error'] * 1e3
                err_Q = MFS_Q_imfit['results'][component]['peak']['error'] * 1e3
                err_U = MFS_U_imfit['results'][component]['peak']['error'] * 1e3
                err_V = MFS_V_imfit['results'][component]['peak']['error'] * 1e3

                RA_P   = MFS_P_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
                DEC_P = MFS_P_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
                
                rms_I = get_imstat_values(MFS_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
                rms_Q = get_imstat_values(MFS_images[2], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                rms_U = get_imstat_values(MFS_images[3], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                rms_V = get_imstat_values(MFS_images[4], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                

                # Calculate additional linear polarisation parameters
                flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U, rms_V, pol_flag, Aq = 0.8)
                
                # Calculate the other Polarisation parameters
                LP_frac     = flux_P0 / flux_I * 100.0
                LP_frac_err = abs(LP_frac) * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2 )

                if pol_flag:
                    LP_EVPA     = np.arctan2(flux_U, flux_Q) * 180.0 / np.pi * 0.5
                    LP_EVPA_err = np.sqrt(flux_U ** 2 * rms_Q ** 2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi * 0.5
   
                else:
                    LP_EVPA = None
                    LP_EVPA_err = None  

                # Append values to dictionary -- Initialize if it doesn't exist
                if 'P_flux_mJy' not in output_dictionary['MFS'][component]:

                    output_dictionary['MFS'][component]['Q_flux_mJy'] = [flux_Q]            
                    output_dictionary['MFS'][component]['U_flux_mJy'] = [flux_U]                
                    output_dictionary['MFS'][component]['V_flux_mJy'] = [flux_V]                
                    output_dictionary['MFS'][component]['P_flux_mJy'] = [flux_P]           
                    output_dictionary['MFS'][component]['P0_flux_mJy'] = [flux_P0]    

                    output_dictionary['MFS'][component]['Q_err_mJy'] = [err_Q]            
                    output_dictionary['MFS'][component]['U_err_mJy'] = [err_U]                
                    output_dictionary['MFS'][component]['V_err_mJy'] = [err_V]                
                    output_dictionary['MFS'][component]['P_err_mJy'] = [err_P]           

                    output_dictionary['MFS'][component]['Q_rms_mJy'] = [rms_Q]            
                    output_dictionary['MFS'][component]['U_rms_mJy'] = [rms_U]                
                    output_dictionary['MFS'][component]['V_rms_mJy'] = [rms_V]                
                    output_dictionary['MFS'][component]['P_rms_mJy'] = [rms_P]        

                    output_dictionary['MFS'][component]['P_RA_deg'] = [RA_P]
                    output_dictionary['MFS'][component]['P_DEC_deg'] = [DEC_P]    

                    output_dictionary['MFS'][component]['LP_frac'] = [LP_frac]
                    output_dictionary['MFS'][component]['LP_frac_err'] = [LP_frac_err]
                    output_dictionary['MFS'][component]['LP_EVPA'] = [LP_EVPA]
                    output_dictionary['MFS'][component]['LP_EVPA_err'] = [LP_EVPA_err]        
            
                else:

                    output_dictionary['MFS'][component]['Q_flux_mJy'].append(flux_Q)         
                    output_dictionary['MFS'][component]['U_flux_mJy'].append(flux_U)               
                    output_dictionary['MFS'][component]['V_flux_mJy'].append(flux_V)               
                    output_dictionary['MFS'][component]['P_flux_mJy'].append(flux_P)           
                    output_dictionary['MFS'][component]['P0_flux_mJy'].append(flux_P0)    

                    output_dictionary['MFS'][component]['Q_err_mJy'].append(err_Q)            
                    output_dictionary['MFS'][component]['U_err_mJy'].append(err_U)               
                    output_dictionary['MFS'][component]['V_err_mJy'].append(err_V)                
                    output_dictionary['MFS'][component]['P_err_mJy'].append(err_P)           

                    output_dictionary['MFS'][component]['Q_rms_mJy'].append(rms_Q)            
                    output_dictionary['MFS'][component]['U_rms_mJy'].append(rms_U)                
                    output_dictionary['MFS'][component]['V_rms_mJy'].append(rms_V)                
                    output_dictionary['MFS'][component]['P_rms_mJy'].append(rms_P)        

                    output_dictionary['MFS'][component]['P_RA_deg'].append(RA_P)
                    output_dictionary['MFS'][component]['P_DEC_deg'].append(DEC_P)    

                    output_dictionary['MFS'][component]['LP_frac'].append(LP_frac)
                    output_dictionary['MFS'][component]['LP_frac_err'].append(LP_frac_err)
                    output_dictionary['MFS'][component]['LP_EVPA'].append(LP_EVPA)
                    output_dictionary['MFS'][component]['LP_EVPA_err'].append(LP_EVPA_err)     
                

        ###########################
        # Extract the CHAN image parameters  #
        ###########################

        msg(f'Fitting CHAN image(s) for prefix {k} with {src_im_suffix}: {prefix}')

        # Look up this prefix's actual CHAN images from the one-time cached
        # scan above -- no glob here, and each prefix's own real channel
        # count/list is used (not assumed identical across timesteps).
        CHAN_images = sorted(_chan_images_by_prefix.get(prefix, []))

        # Reshape to so that each component is a set of Stokes parameters
        if only_intensity:
            CHAN_images_arr = np.array(CHAN_images).reshape(len(CHAN_images), 1)

        else:
            if pol_flag:
               CHAN_images = sorted([im for im in CHAN_images if '-Ptot-' not in im])
            else:
               CHAN_images = sorted([im for im in CHAN_images if '-Plin-' not in im])

            CHAN_images_arr = np.array(CHAN_images).reshape(int(len(CHAN_images) / 5), 5) # reshape to group in frequency for each set of Stokes parameters

        # Build one job per frequency channel, fully self-contained so it can
        # run in a parallel worker (see fit_channel_group). MFS_P_ra_pix/dec
        # are only defined when polarization images exist.
        mfs_freq_ghz  = output_dictionary['MFS']['component0']['freq_GHz'][k]
        mfs_bmaj_asec = output_dictionary['MFS']['component0']['bmaj_asec'][k]
        jobs = []
        for chan_idx, CHAN_images in enumerate(CHAN_images_arr):
            jobs.append({
                'chan_idx': chan_idx,
                'CHAN_images': [str(im) for im in CHAN_images],
                'components': components,
                'MFS_I_ra_pix': MFS_I_ra_pix, 'MFS_I_dec_pix': MFS_I_dec_pix,
                'MFS_P_ra_pix': MFS_P_ra_pix if not only_intensity else None,
                'MFS_P_dec_pix': MFS_P_dec_pix if not only_intensity else None,
                'only_intensity': only_intensity,
                'pol_flag': pol_flag,
                'fix_additional_comps': fix_additional_comps,
                'manual_rms_region': manual_rms_region,
                'src_name': src_name,
                'tag': f'{k}_{chan_idx}',
                'mfs_freq_ghz': mfs_freq_ghz,
                'mfs_bmaj_asec': mfs_bmaj_asec,
            })

        # Lazily size and start the persistent CHAN-fitting pool on first use
        # (and again each time it's recycled -- see below).
        if jobs and chan_pool is None:
            chan_pool_max_workers, chan_pool_recycle_interval = compute_max_workers(
                [jobs[0]['CHAN_images'][0]], len(jobs))
            chan_pool = ProcessPoolExecutor(max_workers=chan_pool_max_workers,
                                             mp_context=chan_pool_ctx)
            chan_pool_prefixes_since_recycle = 0
            msg(f'CHAN-fitting pool started with {chan_pool_max_workers} worker(s), '
                f'recycle interval={chan_pool_recycle_interval} prefixes.')

        # Dispatch in round-robin chunks (one task per worker) and collect
        # results indexed by chan_idx, so completion order doesn't matter --
        # we replay them below in the original channel order.
        chan_results = [None] * len(jobs)
        if jobs:
            chunks = chunk_jobs(jobs, chan_pool_max_workers)
            futures = [chan_pool.submit(fit_channel_group_batch, chunk) for chunk in chunks]
            for future in as_completed(futures):
                for chan_idx, result in future.result():
                    chan_results[chan_idx] = result

        # Replay results in original channel order and append to the output
        # dictionary exactly as the old serial loop did.
        for chan_idx, result in enumerate(chan_results):

            if result is None:
                msg(f'Fitting Failed or Channel Skipped: chan_idx={chan_idx} is likely flagged')
                continue

            freq_GHz = result['freq_GHz']
            bmaj     = result['bmaj_asec']
            bmin     = result['bmin_asec']
            bpa      = result['bpa_deg']

            for component in components:

                comp = result['components'][component]
                flux_I = comp['I_flux_mJy']
                err_I  = comp['I_err_mJy']
                rms_I  = comp['I_rms_mJy']
                RA_I   = comp['I_RA_deg']
                DEC_I  = comp['I_DEC_deg']

                # Append Stokes I values to the output dictionary
                if component not in output_dictionary['CHAN']:

                    output_dictionary['CHAN'][component] = {}

                    output_dictionary['CHAN'][component]['freq_GHz'] = [[freq_GHz]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['bmaj_asec'] = [[bmaj]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['bmin_asec'] = [[bmin]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['bpa_deg'] = [[bpa]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['I_flux_mJy'] = [[flux_I]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['I_err_mJy'] = [[err_I]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['I_rms_mJy'] = [[rms_I]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['I_RA_deg'] = [[RA_I]] + [[] for _ in range(len(prefix_arr) - 1)]
                    output_dictionary['CHAN'][component]['I_DEC_deg'] = [[DEC_I]] + [[] for _ in range(len(prefix_arr) - 1)]

                else:
                    output_dictionary['CHAN'][component]['freq_GHz'][k].append(freq_GHz)
                    output_dictionary['CHAN'][component]['bmaj_asec'][k].append(bmaj)
                    output_dictionary['CHAN'][component]['bmin_asec'][k].append(bmin)
                    output_dictionary['CHAN'][component]['bpa_deg'][k].append(bpa)
                    output_dictionary['CHAN'][component]['I_flux_mJy'][k].append(flux_I)
                    output_dictionary['CHAN'][component]['I_err_mJy'][k].append(err_I)
                    output_dictionary['CHAN'][component]['I_rms_mJy'][k].append(rms_I)
                    output_dictionary['CHAN'][component]['I_RA_deg'][k].append(RA_I)
                    output_dictionary['CHAN'][component]['I_DEC_deg'][k].append(DEC_I)

                # Append polarization values if applicable
                if not only_intensity:

                    flux_Q = comp['Q_flux_mJy']; flux_U = comp['U_flux_mJy']; flux_V = comp['V_flux_mJy']
                    flux_P = comp['P_flux_mJy']; flux_P0 = comp['P0_flux_mJy']
                    err_Q  = comp['Q_err_mJy'];  err_U  = comp['U_err_mJy'];  err_V  = comp['V_err_mJy']
                    err_P  = comp['P_err_mJy']
                    rms_Q  = comp['Q_rms_mJy'];  rms_U  = comp['U_rms_mJy'];  rms_V  = comp['V_rms_mJy']
                    rms_P  = comp['P_rms_mJy']
                    RA_P   = comp['P_RA_deg'];   DEC_P  = comp['P_DEC_deg']
                    LP_frac = comp['LP_frac'];   LP_frac_err = comp['LP_frac_err']
                    LP_EVPA = comp['LP_EVPA'];   LP_EVPA_err = comp['LP_EVPA_err']

                    # Append values to dictionary -- Initialize if it doesn't exist
                    if 'P_flux_mJy' not in output_dictionary['CHAN'][component]:

                        output_dictionary['CHAN'][component]['Q_flux_mJy'] = [[flux_Q]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['U_flux_mJy'] = [[flux_U]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['V_flux_mJy'] = [[flux_V]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['P_flux_mJy'] = [[flux_P]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['P0_flux_mJy'] = [[flux_P0]] + [[] for _ in range(len(prefix_arr) - 1)]

                        output_dictionary['CHAN'][component]['Q_err_mJy'] = [[err_Q]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['U_err_mJy'] = [[err_U]]  + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['V_err_mJy'] = [[err_V]]  + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['P_err_mJy'] = [[err_P]]  + [[] for _ in range(len(prefix_arr) - 1)]

                        output_dictionary['CHAN'][component]['Q_rms_mJy'] = [[rms_Q]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['U_rms_mJy'] = [[rms_U]]  + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['V_rms_mJy'] = [[rms_V]]  + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['P_rms_mJy'] = [[rms_P]]  + [[] for _ in range(len(prefix_arr) - 1)]

                        output_dictionary['CHAN'][component]['P_RA_deg'] = [[RA_P]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['P_DEC_deg'] = [[DEC_P]] + [[] for _ in range(len(prefix_arr) - 1)]

                        output_dictionary['CHAN'][component]['LP_frac'] = [[LP_frac]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['LP_frac_err'] = [[LP_frac_err]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['LP_EVPA'] = [[LP_EVPA]] + [[] for _ in range(len(prefix_arr) - 1)]
                        output_dictionary['CHAN'][component]['LP_EVPA_err'] = [[LP_EVPA_err]] + [[] for _ in range(len(prefix_arr) - 1)]

                    else:

                        output_dictionary['CHAN'][component]['Q_flux_mJy'][k].append(flux_Q)
                        output_dictionary['CHAN'][component]['U_flux_mJy'][k].append(flux_U)
                        output_dictionary['CHAN'][component]['V_flux_mJy'][k].append(flux_V)
                        output_dictionary['CHAN'][component]['P_flux_mJy'][k].append(flux_P)
                        output_dictionary['CHAN'][component]['P0_flux_mJy'][k].append(flux_P0)

                        output_dictionary['CHAN'][component]['Q_err_mJy'][k].append(err_Q)
                        output_dictionary['CHAN'][component]['U_err_mJy'][k].append(err_U)
                        output_dictionary['CHAN'][component]['V_err_mJy'][k].append(err_V)
                        output_dictionary['CHAN'][component]['P_err_mJy'][k].append(err_P)

                        output_dictionary['CHAN'][component]['Q_rms_mJy'][k].append(rms_Q)
                        output_dictionary['CHAN'][component]['U_rms_mJy'][k].append(rms_U)
                        output_dictionary['CHAN'][component]['V_rms_mJy'][k].append(rms_V)
                        output_dictionary['CHAN'][component]['P_rms_mJy'][k].append(rms_P)

                        output_dictionary['CHAN'][component]['P_RA_deg'][k].append(RA_P)
                        output_dictionary['CHAN'][component]['P_DEC_deg'][k].append(DEC_P)

                        output_dictionary['CHAN'][component]['LP_frac'][k].append(LP_frac)
                        output_dictionary['CHAN'][component]['LP_frac_err'][k].append(LP_frac_err)
                        output_dictionary['CHAN'][component]['LP_EVPA'][k].append(LP_EVPA)
                        output_dictionary['CHAN'][component]['LP_EVPA_err'][k].append(LP_EVPA_err)


        # Write RM Synthesis files (freq, I, Q, U, dI, dQ, dU) in Jy, one file
        # per component per epoch.  These are the direct inputs to RM-Tools.
        if not only_intensity:
            for component in components:
                rmsynth_arr = np.array([
                    np.array(output_dictionary['CHAN'][component]['freq_GHz'][k]) * 1e9,
                    np.array(output_dictionary['CHAN'][component]['I_flux_mJy'][k])  / 1e3,
                    np.array(output_dictionary['CHAN'][component]['Q_flux_mJy'][k])  / 1e3,
                    np.array(output_dictionary['CHAN'][component]['U_flux_mJy'][k])  / 1e3,
                    np.array(output_dictionary['CHAN'][component]['I_rms_mJy'][k])   / 1e3,
                    np.array(output_dictionary['CHAN'][component]['Q_rms_mJy'][k])   / 1e3,
                    np.array(output_dictionary['CHAN'][component]['U_rms_mJy'][k])   / 1e3])

                rmsynth_fname = '{}_{}_rmsynth.txt'.format(prefix, component).replace(image_directory, 'RESULTS')
                rmsynth_fname = _apply_label(rmsynth_fname, src_name, label)
                np.savetxt(rmsynth_fname, rmsynth_arr.T)
                msg(f'  RM synthesis file written: {rmsynth_fname}')

        # Compute the spectral index for every epoch unconditionally (stored in
        # JSON) -- this compute-only pass (save_plot=False) is NOT gated by
        # SPEC_PLOT_SNR_THRESH. That threshold only filters which epochs get an
        # actual PNG rendered, in the plotting pass below/after the JSON save.
        # Plotting is deferred until after the JSON is saved, so that data
        # is preserved even if a plotting call crashes.
        for component in components:
            plot_stokes_spectrum(output_dictionary, src_name, prefix, component, k,
                                save_plot=False, label=label, pol_flag=pol_flag)

        # Remove this timestep's estimate/check_position scratch files
        cleanup_estimate_files(src_name)
        log_memory_usage(label=f'prefix {k+1}/{len(prefix_arr)} ({os.path.basename(prefix)}) done --')

        if chan_pool is not None:
            chan_pool_prefixes_since_recycle += 1
            if chan_pool_prefixes_since_recycle >= chan_pool_recycle_interval:
                msg(f'Recycling CHAN-fitting pool after {chan_pool_prefixes_since_recycle} prefixes '
                    f'(interval={chan_pool_recycle_interval}) to bound per-worker memory growth.')
                chan_pool.shutdown(wait=True)
                chan_pool = None
                chan_pool_prefixes_since_recycle = 0

    if chan_pool is not None:
        chan_pool.shutdown(wait=True)
        msg('CHAN-fitting pool shut down.')

    # ------------------------------------------------------------------
    # Save all data files before any plotting
    # ------------------------------------------------------------------
    if len(prefix_arr) == 0:
        msg('WARNING: No prefixes found -- nothing to process. '
            'Check image_directory, image_identifier, and image_suffix in rmsynth_info.json.')
        return 0

    # Save the full dictionary
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, 'w') as j:
        json.dump(output_dictionary, j, indent = 4)
    msg(f'  Output JSON saved: {json_path}')

    # ------------------------------------------------------------------
    # Plotting (after all files are safely on disk)
    # ------------------------------------------------------------------
    if is_time_resolved:
        # Plot IQUV spectrum for all epochs meeting the MFS S/N threshold
        for component in components:
            _plot_detected_spectra(output_dictionary, src_name, prefix_arr, component, label,
                                   pol_flag=pol_flag)

        # Plot MFS light curve across all epochs
        plot_light_curve(output_dictionary, src_name, prefix_arr, label=label)
    else:
        # Single-epoch: plot spectrum for each component
        for component in components:
            plot_stokes_spectrum(output_dictionary, src_name, prefix_arr[0],
                                component, 0, save_plot=True, label=label,
                                pol_flag=pol_flag)

    return 0


def _plot_detected_spectra(output_dictionary, src_name, prefix_arr, component, label, pol_flag=False):
    '''
    Plot the IQUV spectrum for every time-resolved epoch of `component` whose
    MFS Stokes I S/N clears SPEC_PLOT_SNR_THRESH, and log a line per epoch
    either way, so it's visible from the log alone -- not just by checking
    which files landed on disk -- whether each epoch's spectrum plot was
    saved or skipped for being sub-threshold.
    '''
    mfs_I_arr   = np.array(output_dictionary['MFS'][component]['I_flux_mJy'])
    mfs_rms_arr = np.array(output_dictionary['MFS'][component]['I_rms_mJy'])
    mfs_snr_arr = np.where(mfs_rms_arr > 0, mfs_I_arr / mfs_rms_arr, 0.0)
    det_mask    = mfs_snr_arr >= SPEC_PLOT_SNR_THRESH
    n_det       = int(np.sum(det_mask))
    msg(f'  {component}: {n_det}/{len(det_mask)} epochs with MFS S/N >= {SPEC_PLOT_SNR_THRESH} -- plotting spectra')
    for k in range(len(det_mask)):
        epoch_desc = f'epoch {k} ({os.path.basename(prefix_arr[k])})'
        if det_mask[k]:
            msg(f'    {epoch_desc}: MFS S/N = {mfs_snr_arr[k]:.1f} >= {SPEC_PLOT_SNR_THRESH} -- spectrum plot SAVED')
            plot_stokes_spectrum(output_dictionary, src_name, prefix_arr[k],
                                component, int(k), save_plot=True, label=label,
                                pol_flag=pol_flag)
        else:
            msg(f'    {epoch_desc}: MFS S/N = {mfs_snr_arr[k]:.1f} <  {SPEC_PLOT_SNR_THRESH} -- spectrum plot SKIPPED')


def plot_stokes_spectrum(output_dictionary, src_name, prefix, component, k,
                         save_plot=True, label=None, pol_flag=False):
    '''
    Diagnostic plot of per-channel Stokes I, Q, U, V spectra for one source
    component in one time-interval epoch, stacked vertically with a shared
    frequency axis.

    Layout
    ------
    Four panels stacked vertically, sharing the same x-axis:
      - Panel 1 (top):  Stokes I  -- linear axes; MFS reference dashed line;
                                      spectral index power-law fit overlaid when
                                      MFS S/N > 30 and >= 3 valid channels.
      - Panels 2-4:     Stokes Q, U, V -- linear y-axis (fluxes can be negative);
                                          robust MAD-based y-axis range;
                                          MFS value as a dashed horizontal reference.

    All panels use RMS error bars (not the imfit formal errors, which are stored
    in the JSON for completeness but not plotted).

    Spectral index fitting (Stokes I panel only)
    --------------------------------------------
    Conditions: MFS S/N > 30, >= 3 channels with flux > 0 and per-channel S/N >= 3.
    Method: weighted log-log linear regression (np.polyfit), weights = S/N per channel.
    alpha and alpha_err are stored back into output_dictionary['MFS'][component]
    and therefore appear in the final JSON without a second write.  None is stored
    if the fit is skipped.

    Parameters
    ----------
    output_dictionary : dict  -- full extraction dictionary (modified in place)
    src_name          : str   -- true source name; prefix is a literal glob match
                                 on this string, so it's what label (if set)
                                 gets substituted for below
    prefix            : str   -- epoch path prefix (used for output filename)
    component         : str   -- component key, e.g. 'component0'
    k                 : int   -- epoch index within the prefix list
    save_plot         : bool  -- if False, compute and store spectral index but
                                 skip figure generation (default True)
    label             : str or None -- display/output name override (rmsynth_info.json
                                 'label' field, via extract_polarization_properties).
                                 None/unset means "use src_name" -- unchanged behaviour.
    pol_flag          : bool  -- True when a polarization angle calibrator was used, so
                                 P is Plin = sqrt(Q^2+U^2) rather than Ptot = sqrt(Q^2+U^2+V^2).
                                 Only then is P overlaid on the Q and U panels -- Ptot mixes
                                 in V and doesn't belong on Q/U axes.
    '''

    label = label or src_name

    # Ensure the output directory exists
    plot_dir = os.path.join('RESULTS', 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Retrieve per-channel data arrays for this epoch
    # ------------------------------------------------------------------
    chan  = output_dictionary['CHAN'][component]
    freq_arr  = np.array(chan['freq_GHz'][k])

    if len(freq_arr) == 0:
        msg(f'  plot_stokes_spectrum: no channel data for {component} epoch {k}, skipping plot')
        return

    I_flux = np.array(chan['I_flux_mJy'][k])
    I_rms  = np.array(chan['I_rms_mJy'][k])
    # err_arr (imfit formal errors) retrieved only to satisfy completeness;
    # RMS is used for all error bars -- see docstring.
    _ = np.array(chan['I_err_mJy'][k])

    # Q, U, V are only available when not running in Stokes-I-only mode
    only_intensity = 'Q_flux_mJy' not in chan
    if not only_intensity:
        Q_flux = np.array(chan['Q_flux_mJy'][k])
        Q_rms  = np.array(chan['Q_rms_mJy'][k])
        U_flux = np.array(chan['U_flux_mJy'][k])
        U_rms  = np.array(chan['U_rms_mJy'][k])
        V_flux = np.array(chan['V_flux_mJy'][k])
        V_rms  = np.array(chan['V_rms_mJy'][k])

    # P is only overlaid on the Q/U panels when it's Plin (pol_flag True) --
    # Ptot includes V and isn't meaningful on Q/U axes.
    show_P = pol_flag and not only_intensity and 'P_flux_mJy' in chan
    if show_P:
        P_flux = np.array(chan['P_flux_mJy'][k])
        P_rms  = np.array(chan['P_rms_mJy'][k])

    # ------------------------------------------------------------------
    # MAD-based outlier rejection on error bars: if ANY of I, Q, U, V has a
    # crazy RMS in a channel, mask that channel from I, Q, U, V, and P alike.
    # This prevents single pathological channels from distorting the plot.
    # ------------------------------------------------------------------
    all_rms = [I_rms]
    if not only_intensity:
        all_rms += [Q_rms, U_rms, V_rms]

    good = np.ones(len(freq_arr), dtype=bool)
    for rms_arr in all_rms:
        med  = np.nanmedian(rms_arr)
        mad  = np.nanmedian(np.abs(rms_arr - med))
        sig  = 1.4826 * mad
        if sig > 0:
            good &= np.abs(rms_arr - med) <= 5.0 * sig

    n_masked = int(np.sum(~good))
    if n_masked > 0:
        msg(f'  Masking {n_masked}/{len(good)} channels with outlier RMS values')
        freq_arr = freq_arr[good]
        I_flux   = I_flux[good]
        I_rms    = I_rms[good]
        if not only_intensity:
            Q_flux = Q_flux[good]
            Q_rms  = Q_rms[good]
            U_flux = U_flux[good]
            U_rms  = U_rms[good]
            V_flux = V_flux[good]
            V_rms  = V_rms[good]
        if show_P:
            P_flux = P_flux[good]
            P_rms  = P_rms[good]

    if len(freq_arr) == 0:
        msg(f'  plot_stokes_spectrum: all channels masked for {component} epoch {k}, skipping plot')
        return

    # MFS reference values for this epoch
    mfs       = output_dictionary['MFS'][component]
    mfs_I     = mfs['I_flux_mJy'][k]
    mfs_I_rms = mfs['I_rms_mJy'][k]
    mfs_freq  = mfs['freq_GHz'][k]
    mfs_snr   = mfs_I / mfs_I_rms if mfs_I_rms > 0 else 0.0
    if not only_intensity:
        mfs_Q = mfs['Q_flux_mJy'][k]
        mfs_U = mfs['U_flux_mJy'][k]
        mfs_V = mfs['V_flux_mJy'][k]

    # ------------------------------------------------------------------
    # Spectral index fit on Stokes I (log-log space)
    # ------------------------------------------------------------------
    alpha        = None
    alpha_err    = None
    chi2         = None
    ndof         = None
    chi2_red     = None
    fit_nu       = None
    fit_S        = None
    flagged_freqs = np.array([])

    # Quality mask: positive flux and per-channel S/N >= 3
    valid = (I_flux > 0) & (I_rms > 0) & (I_flux / I_rms >= 3.0)

    if mfs_snr > SPEC_INDEX_SNR_THRESH and np.sum(valid) >= 3:
        nu_v = freq_arr[valid]
        S_v  = I_flux[valid]
        w_v  = S_v / I_rms[valid]          # S/N weights (= 1/sigma_logS in log space)

        # Iterative MAD-based outlier rejection on flux values
        _iter = 0
        while len(nu_v) >= 3:
            _med  = np.median(S_v)
            _mad  = np.median(np.abs(S_v - _med))
            _sig  = 1.4826 * _mad
            if _sig <= 0:
                break
            _keep     = np.abs(S_v - _med) <= SPEC_INDEX_MAD_CLIP * _sig
            n_clipped = int(np.sum(~_keep))
            msg(f'  MAD clip iter {_iter + 1} ({component}, epoch {k}): '
                f'median = {_med:.3f} mJy  sigma = {_sig:.3f} mJy  '
                f'threshold = {SPEC_INDEX_MAD_CLIP * _sig:.3f} mJy  '
                f'clipped = {n_clipped}  remaining = {int(np.sum(_keep))}')
            if n_clipped == 0:
                break
            flagged_freqs = np.concatenate([flagged_freqs, nu_v[~_keep]])
            nu_v  = nu_v[_keep]
            S_v   = S_v[_keep]
            w_v   = w_v[_keep]
            _iter += 1

        if len(nu_v) < 3:
            msg(f'  Spectral index fit skipped ({component}, epoch {k}): '
                f'only {len(nu_v)} channels remain after MAD clipping (need >= 3)')
        else:
            # Fit in log space: log(S) = alpha*log(nu/nu_ref) + log(S_ref)
            log_nu = np.log(nu_v / mfs_freq)
            log_S  = np.log(S_v)
            try:
                coeffs, cov = np.polyfit(log_nu, log_S, 1, w=w_v, cov=True)
                alpha     = coeffs[0]
                alpha_err = np.sqrt(cov[0, 0])

                # Chi-squared: residuals in log space weighted by S/N
                # sigma_logS = 1/SNR, so chi2 = sum((residual * w_v)^2)
                log_S_model = coeffs[0] * log_nu + coeffs[1]
                residuals   = log_S - log_S_model
                chi2        = float(np.sum((residuals * w_v) ** 2))
                ndof        = len(nu_v) - 2
                chi2_red    = chi2 / ndof if ndof > 0 else None

                nu_plot = np.linspace(nu_v.min(), nu_v.max(), 200)
                fit_S   = np.exp(coeffs[1]) * (nu_plot / mfs_freq) ** alpha
                fit_nu  = nu_plot

                msg(f'  Spectral index fit ({component}, epoch {k}): '
                    f'alpha = {alpha:.3f} +/- {alpha_err:.3f}  '
                    f'chi2 = {chi2:.2f}  chi2_red = {chi2_red:.2f}  ndof = {ndof}  '
                    f'(MFS S/N = {mfs_snr:.1f} >= {SPEC_INDEX_SNR_THRESH}, '
                    f'N_chan = {len(nu_v)} / {int(np.sum(valid))})')
            except Exception as e:
                msg(f'  WARNING: Spectral index fit failed for {component} epoch {k}: {e}')
    else:
        if mfs_snr <= SPEC_INDEX_SNR_THRESH:
            msg(f'  Spectral index fit skipped ({component}, epoch {k}): '
                f'MFS S/N = {mfs_snr:.1f} (< {SPEC_INDEX_SNR_THRESH})')
        else:
            msg(f'  Spectral index fit skipped ({component}, epoch {k}): '
                f'only {int(np.sum(valid))} valid channels (need >= 3)')

    # Store alpha and chi2 in the MFS dict so they appear in the JSON.
    # Initialise as a list on the first epoch; append on subsequent epochs.
    # Guard against duplicates when the function is re-called for plotting only.
    if 'alpha' not in mfs:
        mfs['alpha']     = [alpha]
        mfs['alpha_err'] = [alpha_err]
        mfs['chi2']      = [chi2]
        mfs['ndof']      = [ndof]
        mfs['chi2_red']  = [chi2_red]
    elif len(mfs['alpha']) <= k:
        mfs['alpha'].append(alpha)
        mfs['alpha_err'].append(alpha_err)
        mfs['chi2'].append(chi2)
        mfs['ndof'].append(ndof)
        mfs['chi2_red'].append(chi2_red)

    # If save_plot is False, we only needed the spectral index computation above.
    if not save_plot:
        return

    # ------------------------------------------------------------------
    # Build figure: 4 panels stacked vertically, shared x-axis
    # (1 panel if Stokes-I-only mode)
    # ------------------------------------------------------------------
    n_panels = 1 if only_intensity else 4
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 3 * n_panels),
                             sharex=True)
    if n_panels == 1:
        axes = [axes]   # make iterable

    # --- Panel 0: Stokes I (linear axes) ---
    ax = axes[0]
    ax.errorbar(freq_arr, I_flux, yerr=I_rms,
                fmt='o', markersize=4, capsize=3, color='royalblue',
                label='Stokes I', zorder=3)

    # Robust y-axis range: MAD-based 5-sigma outlier rejection
    med_I   = np.nanmedian(I_flux)
    mad_I   = np.nanmedian(np.abs(I_flux - med_I))
    sigma_I = 1.4826 * mad_I
    if sigma_I > 0:
        inlier = np.abs(I_flux - med_I) <= 5.0 * sigma_I
    else:
        inlier = np.ones(len(I_flux), dtype=bool)

    range_vals = np.concatenate([
        I_flux[inlier] + I_rms[inlier],
        I_flux[inlier] - I_rms[inlier],
        [mfs_I]])
    lo = np.nanmin(range_vals)
    hi = np.nanmax(range_vals)
    if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
        span = hi - lo
        pad  = 0.15 * span
        ax.set_ylim(lo - pad, hi + pad)

    ax.axhline(y=mfs_I, color='k', linestyle='--', linewidth=1.5,
               label=f'MFS: {mfs_I:.2f} mJy  (S/N = {mfs_snr:.1f})', zorder=2)
    if fit_nu is not None:
        ax.plot(fit_nu, fit_S, color='tomato', linewidth=1.8,
                label=rf'$\alpha$ = {alpha:.2f} $\pm$ {alpha_err:.2f}', zorder=4)
    ax.set_ylabel('I  (mJy/beam)', fontsize=10)
    ax.grid(True, which='both', alpha=0.3, linestyle=':')
    ax.legend(fontsize=9, loc='best')

    # --- Panels 1-3: Stokes Q, U, V (linear y, actual data range) ---
    if not only_intensity:
        for ax, flux, rms, label, mfs_val, colour in zip(
                axes[1:],
                [Q_flux,   U_flux,   V_flux],
                [Q_rms,    U_rms,    V_rms],
                ['Q',      'U',      'V'],
                [mfs_Q,    mfs_U,    mfs_V],
                ['darkorange', 'forestgreen', 'mediumpurple']):

            ax.errorbar(freq_arr, flux, yerr=rms,
                        fmt='o', markersize=4, capsize=3, color=colour,
                        label=f'Stokes {label}', zorder=3)
            ax.axhline(y=mfs_val, color='k', linestyle='--', linewidth=1.5,
                       label=f'MFS: {mfs_val:.2f} mJy', zorder=2)
            ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8, zorder=1)

            # Robust y-axis range: MAD-based 5-sigma outlier rejection,
            # then use the actual min/max of inlier data +/- error bars
            # and the MFS reference value.
            med   = np.nanmedian(flux)
            mad   = np.nanmedian(np.abs(flux - med))
            sigma = 1.4826 * mad
            if sigma > 0:
                inlier = np.abs(flux - med) <= 5.0 * sigma
            else:
                inlier = np.ones(len(flux), dtype=bool)

            range_vals = np.concatenate([
                flux[inlier] + rms[inlier],
                flux[inlier] - rms[inlier],
                [mfs_val, 0.0]])
            lo = np.nanmin(range_vals)
            hi = np.nanmax(range_vals)
            if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
                span = hi - lo
                pad  = 0.15 * span
                ax.set_ylim(lo - pad, hi + pad)

            ax.set_ylabel(f'{label}  (mJy/beam)', fontsize=10)
            ax.grid(True, which='both', alpha=0.3, linestyle=':')
            ax.legend(fontsize=9, loc='best')

        # Overlay Plin on the Q and U panels only (axes[1], axes[2] -- not V)
        if show_P:
            for ax in axes[1:3]:
                ax.errorbar(freq_arr, P_flux, yerr=P_rms,
                            fmt='o', markersize=3, capsize=2, color='black',
                            alpha=0.5, label='Stokes P (lin)', zorder=2.5)
                ax.legend(fontsize=9, loc='best')

    # Overlay flagged channels on all panels as translucent red spans.
    # Width of each span = minimum channel separation in the frequency array.
    if len(flagged_freqs) > 0 and len(freq_arr) >= 2:
        _chan_width = float(np.min(np.diff(np.sort(freq_arr))))
        for _f in flagged_freqs:
            for _ax in axes:
                _ax.axvspan(_f - _chan_width / 2, _f + _chan_width / 2,
                            color='red', alpha=0.2, zorder=0, linewidth=0)

    # Shared x-axis label on the bottom panel only
    axes[-1].set_xlabel('Frequency (GHz)', fontsize=11)

    fig.suptitle(f'{label}  |  {component}  |  epoch {k}', fontsize=12, y=1.01)
    plt.tight_layout()

    # prefix is a literal glob match on src_name (see docstring); _apply_label
    # carries the label into the output filename only (never a directory).
    out_basename = _apply_label(os.path.basename(prefix), src_name, label)

    plot_fname = os.path.join(
        plot_dir,
        f'{out_basename}_{component}_IQUV_spectrum.png')
    plt.savefig(plot_fname, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  IQUV spectrum plot saved: {plot_fname}')


def plot_light_curve(output_dictionary, src_name, prefix_arr, label=None):
    '''
    Plot MFS IQUV light curves for time-resolved data. Always produced
    whenever the run is time-resolved, and always plots every epoch -- this
    is a diagnostic plot, not a detection-filtered one (per-channel spectrum
    plots from plot_stokes_spectrum() are the ones restricted to high-S/N
    epochs, see SPEC_PLOT_SNR_THRESH).

    One figure per source with one column per component and up to four rows:
      Row 1: Stokes I flux density (mJy/beam) vs time
      Row 2: Stokes Q flux density (mJy/beam) vs time
      Row 3: Stokes U flux density (mJy/beam) vs time
      Row 4: Stokes V flux density (mJy/beam) vs time

    Rows 2-4 are only shown when polarization data is available.

    Every epoch is always plotted. A MAD clip (the same threshold used for
    the spectral-index fit, SPEC_INDEX_MAD_CLIP) is applied to the error
    array of each row to identify epochs with anomalously large error bars;
    those epochs are excluded only when setting the y-axis range, so a
    single bad-error outlier can't blow out the scale for everything else.

    Parameters
    ----------
    output_dictionary : dict  -- full extraction dictionary
    src_name          : str   -- true source name; prefix_arr entries are literal
                                 glob matches on this string, so it's what label
                                 (if set) gets substituted for below
    prefix_arr        : list  -- list of prefix strings (one per epoch),
                                 used to build the output filename
    label             : str or None -- display/output name override (rmsynth_info.json
                                 'label' field, via extract_polarization_properties).
                                 None/unset means "use src_name" -- unchanged behaviour.
    '''

    label = label or src_name

    plot_dir = os.path.join('RESULTS', 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)

    mfs = output_dictionary['MFS']
    components = sorted(mfs.keys())

    # Determine whether polarization data exists
    has_pol = 'Q_flux_mJy' in mfs[components[0]]
    n_rows  = 4 if has_pol else 1
    n_cols  = len(components)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(8 * n_cols, 3.5 * n_rows),
                             squeeze=False, sharex=True)

    # Component colours
    colours = ['royalblue', 'tomato', 'forestgreen', 'mediumpurple',
               'darkorange', 'teal', 'crimson', 'goldenrod']

    def mad_clip_ylim(flux, rms):
        """MAD-clip the error array (a MAD clip, same threshold as the
        spectral-index fit) to find epochs with well-behaved error bars, and
        use only those to set the y-axis range. All epochs are still plotted
        regardless of this mask."""
        med  = np.nanmedian(rms)
        mad  = np.nanmedian(np.abs(rms - med))
        sig  = 1.4826 * mad
        if sig > 0:
            keep = np.abs(rms - med) <= SPEC_INDEX_MAD_CLIP * sig
        else:
            keep = np.ones(len(rms), dtype=bool)
        if not np.any(keep):
            keep = np.ones(len(rms), dtype=bool)
        lo = np.nanmin(flux[keep] - rms[keep])
        hi = np.nanmax(flux[keep] + rms[keep])
        if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
            span = hi - lo
            pad  = 0.15 * span if span > 0 else 0.5
            return lo - pad, hi + pad
        return None

    def plot_stokes_row(ax, dates, flux, rms, label, col, legend_label=None):
        """Plot every epoch for one Stokes parameter, unconditionally."""
        ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8)
        ax.set_ylabel(f'Stokes {label} (mJy/beam)', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle=':')
        ax.errorbar(dates, flux, yerr=rms, fmt='o', markersize=5, capsize=3,
                    color=col, label=legend_label)
        if legend_label is not None:
            ax.legend(fontsize=8, loc='best')
        ylim = mad_clip_ylim(flux, rms)
        if ylim is not None:
            ax.set_ylim(*ylim)

    for j, component in enumerate(components):
        comp = mfs[component]

        # Parse observation dates (centre of each exposure)
        dates = np.array([datetime.datetime.fromisoformat(d) for d in comp['time_ctr_isot']])
        col   = colours[j % len(colours)]

        I_flux = np.array(comp['I_flux_mJy'])
        I_rms  = np.array(comp['I_rms_mJy'])

        ax = axes[0, j]
        ax.set_title(component, fontsize=10)
        plot_stokes_row(ax, dates, I_flux, I_rms, 'I', col, legend_label=component)

        if has_pol and 'Q_flux_mJy' in comp:
            Q_flux = np.array(comp['Q_flux_mJy'])
            Q_rms  = np.array(comp['Q_rms_mJy'])
            plot_stokes_row(axes[1, j], dates, Q_flux, Q_rms, 'Q', col)

            U_flux = np.array(comp['U_flux_mJy'])
            U_rms  = np.array(comp['U_rms_mJy'])
            plot_stokes_row(axes[2, j], dates, U_flux, U_rms, 'U', col)

            V_flux = np.array(comp['V_flux_mJy'])
            V_rms  = np.array(comp['V_rms_mJy'])
            plot_stokes_row(axes[3, j], dates, V_flux, V_rms, 'V', col)

    # Format the time axis on the bottom row
    for ax in axes[-1, :]:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        ax.set_xlabel('Time (UT)', fontsize=10)
        for tick in ax.get_xticklabels():
            tick.set_rotation(30)
            tick.set_ha('right')

    # If the time span covers more than 1 day, use a date+time format instead
    all_dates = []
    for component in components:
        all_dates += [datetime.datetime.fromisoformat(d)
                      for d in mfs[component]['time_ctr_isot']]
    if all_dates:
        span = max(all_dates) - min(all_dates)
        if span.total_seconds() > 86400:
            for ax in axes[-1, :]:
                ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M'))

    fig.suptitle(f'{label}  |  MFS light curve', fontsize=12, y=1.02)
    plt.tight_layout()

    # Strip the time index (e.g. -t0051) from the prefix for the filename.
    # prefix_arr[0] is a literal glob match on src_name (see docstring); the
    # _apply_label call carries the label into the filename only, same
    # pattern as plot_stokes_spectrum().
    lc_basename = _apply_label(re.sub(r'-t\d+$', '', os.path.basename(prefix_arr[0])), src_name, label)
    plot_fname = os.path.join(
        plot_dir,
        f'{lc_basename}_{label}_MFS_lightcurve.png')
    plt.savefig(plot_fname, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  MFS light curve plot saved: {plot_fname}')


def main():
    
    # Load in the rmsynthesis data
    with open(cfg.RMSYN_INFO_FILE, 'r') as j:
        rmsynth_info = json.load(j)

    # Check to see if there is a Polarization angle calibrator -- pol_flag = Trye means that you do have
    pol_flag=False
    if cfg.POLANG_NAME != '':
        pol_flag = True

    # Optional per-entry output-name override. Missing key or empty string
    # means "use source_name" (unchanged default behaviour) -- see out_label
    # in extract_polarization_properties() docstring.
    labels = rmsynth_info.get('label', [])

    # Iterate through sources as specified in rmsynth_info.json
    for k in range(len(rmsynth_info['image_directory']))[:]:
   
        if not os.path.exists(rmsynth_info['image_directory'][k]):
            msg('Skipping {} as the directory does not exist'.format(rmsynth_info['image_directory'][k]))
            continue

        # Construct image identifier based on input options
        if rmsynth_info["image_timing"][k]:
            src_im_identifier = cfg.CWD +'/{}/*{}*.ms_{}-t'.format(rmsynth_info["image_directory"][k], 
                                                                                            rmsynth_info["source_name"][k], 
                                                                                            rmsynth_info["image_identifier"][k])
        else:
            src_im_identifier = cfg.CWD +'/{}/*{}*.ms_{}-'.format(rmsynth_info["image_directory"][k], 
                                                                                            rmsynth_info["source_name"][k], 
                                                                                           rmsynth_info["image_identifier"][k])                                                                                           
        # Parse source positions from rmsynth_info.json using parse_casa_position().
        # Accepts both CASA HMS/DMS format ('HH:MM:SS.SS,+-DD.MM.SS.SS') and
        # plain decimal degrees ('260.43,-16.20' or with 'deg' suffix).
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

        ###########################
        # Add the position prediction from here #
        ###########################
        
        msg(f'Starting property extraction for identifier: {src_im_identifier}')
        extract_polarization_properties(rmsynth_info["source_name"][k],
            src_im_identifier,
            rmsynth_info['image_suffix'][k], 
            src_ra, 
            src_dec, 
            rmsynth_info['source_ulim'][k],
            pol_flag,
            fix_dec_colons(rmsynth_info['rms_region'][k], 'rms_region'),
            rmsynth_info['image_directory'][k],
            rmsynth_info["image_identifier"][k],
            fix_additional_comps = True,
            out_label = labels[k] if k < len(labels) else '')



if __name__  == "__main__":
    main()
