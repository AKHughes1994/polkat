# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

# Modular-CASA version: imfit/imstat/imhead/table are plain imports from
# casatasks/casatools rather than casashell-injected globals, so this runs
# under a normal `python-pycasa` interpreter (not `casa -c`) as long as that
# interpreter has casatasks/casatools installed. This also means astropy is
# available, so estimate_memory_per_worker_bytes reads FITS headers directly
# instead of going through imhead.

import glob
import os
import datetime
import subprocess
import sys
import json
import shutil
import re
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
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

from casatasks import imfit, imstat, imhead
from casatools import table
tb = table()


# =============================================================================
# Global options
# =============================================================================

# Image/output identifier used in image matching and output filenames.
IDENTIFIER = 'pcalmask'

# --- Position referencing thresholds -----------------------------------------
MFS_P_SNR_THRESH   = 10.0   # MFS: minimum P S/N to reference Q/U against P position
CHAN_P_SNR_THRESH  = 10.0   # CHAN: minimum P S/N to reference Q/U against P position
MIN_FRAC_POL       = 0.003  # Minimum fractional polarisation (P/I > 0.3%) to allow P referencing
PLIN_FREE_FIT_SNR  = 8.0    # Plin: minimum S/N for free position fitting

# --- Calibrator position fixing ----------------------------------------------
# Set True to force all Stokes positions to the I position for non-target
# calibrators whose field name does not contain any EXCLUDED_CALIBRATOR_SUBSTRINGS.
USE_CALIBRATOR_POSITION_FIXING = True
EXCLUDED_CALIBRATOR_SUBSTRINGS = ['J1331', 'J1733', 'J1424']

# --- Stokes V position forcing -----------------------------------------------
# Sources whose field name does NOT contain any of these substrings will have
# their Stokes V position forced to the Stokes I reference position.
# Add a source substring here to allow V to be fit freely for that source.
V_FORCE_EXCLUDED_SUBSTRINGS = []

# --- Spectral index fit thresholds (defaults pulled from config.py) ----------
SPEC_INDEX_SNR_THRESH = cfg.RMSYN_SPEC_INDEX_SNR_THRESH  # Minimum MFS Stokes I S/N to attempt the fit
SPEC_INDEX_MAD_CLIP   = cfg.RMSYN_SPEC_INDEX_MAD_CLIP    # Iterative MAD outlier rejection threshold (sigma)

# =============================================================================


def _mjd_to_isot(mjd):
    """Convert MJD (float) to an ISOT-format datetime string."""
    _epoch = datetime.datetime(1858, 11, 17, 0, 0, 0)
    return (_epoch + datetime.timedelta(days=mjd)).isoformat()


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
    avg/min/max (to catch one runaway worker among many healthy ones) and the
    change since the previous call (to make steady growth across channels --
    e.g. from a leak in casacore's table/image caching inside long-lived
    worker processes -- visible without having to manually diff log lines).
    Prefers psutil when available; falls back to reading /proc directly
    (multiprocessing.active_children() for worker PIDs) if it isn't."""
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


def compute_max_workers(sample_images, n_jobs):
    """sample_images is used only to estimate per-worker memory (a handful of
    representative images is enough -- they're all roughly the same size).
    n_jobs is the actual number of parallel work units (channels), which is
    what should cap the worker count -- NOT len(sample_images).

    Returns (workers, _unused). The second value mirrors dev's pool-recycle
    interval, but tkat runs one field/epoch per script invocation (no
    multi-epoch loop to recycle the pool across), so it is computed for
    parity but not used here."""
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

    return workers, None


def chunk_jobs(jobs, n_chunks):
    """Round-robin split so any systematic per-channel cost trend (e.g. RFI
    flagging clustered in frequency) gets spread evenly across workers."""
    n_chunks = max(1, min(n_chunks, len(jobs)))
    chunks = [[] for _ in range(n_chunks)]
    for i, job in enumerate(jobs):
        chunks[i % n_chunks].append(job)
    return [c for c in chunks if c]


def get_source_from_ms(myms, source_name):
    """
    Query MS to find source by name and return position and intent.
    For per-scan sources (e.g., "J1744-5144_scan2"), strips the scan suffix
    to match the base field name in the MS.
    """
    msg(f"Querying MS for source: {source_name}")
    
    # Strip scan suffix if present (e.g., "J1744-5144_scan2" -> "J1744-5144")
    base_field_name = source_name
    scan_number = None
    scan_match = re.search(r'(.+)_scan(\d+)$', source_name)
    if scan_match:
        base_field_name = scan_match.group(1)
        scan_number = scan_match.group(2)
        msg(f"  Detected per-scan source: base={base_field_name}, scan={scan_number}")
    
    # Query FIELD table for source
    tb.open(myms + '/FIELD')
    field_names = tb.getcol('NAME')
    field_dirs = tb.getcol('PHASE_DIR')
    tb.close()
    
    # Find matching source (case-insensitive, using base field name)
    source_idx = None
    for i, fname in enumerate(field_names):
        if fname.lower() == base_field_name.lower():
            source_idx = i
            break

    if source_idx is None:
        msg(f"ERROR: Source '{base_field_name}' not found in MS FIELD table")
        msg(f"  Available fields: {', '.join(field_names)}")
        return {'name': source_name, 'ra_deg': None, 'dec_deg': None, 'intent': None, 'found': False}
    
    # Get position (convert radians to degrees)
    # PHASE_DIR shape is (2, 1, num_fields) where [0]=RA, [1]=DEC
    ra_rad = field_dirs[0][0][source_idx]
    dec_rad = field_dirs[1][0][source_idx]
    ra_deg = np.degrees(ra_rad)
    dec_deg = np.degrees(dec_rad)
    
    # Get intent from STATE table via main table
    tb.open(myms)
    field_col = tb.getcol('FIELD_ID')
    state_col = tb.getcol('STATE_ID')
    tb.close()
    
    # Find scans with this field
    field_mask = field_col == source_idx
    field_states = np.unique(state_col[field_mask])
    
    # Get intent from STATE table
    tb.open(myms + '/STATE')
    obs_modes = tb.getcol('OBS_MODE')
    tb.close()
    
    intent = obs_modes[field_states[0]] if len(field_states) > 0 else 'UNKNOWN'
    
    msg(f"Found source: {field_names[source_idx]}")
    msg(f"  Position: RA={ra_deg:.6f} deg, Dec={dec_deg:.6f} deg")
    msg(f"  Intent: {intent}")
    
    return {
        'name': field_names[source_idx],
        'ra_deg': ra_deg,
        'dec_deg': dec_deg,
        'intent': intent,
        'found': True
    }


def parse_casa_position(position_str):
    """
    Parse CASA format position string to decimal degrees.
    CASA format: RA is colon-separated (HH:MM:SS.SS), Dec is period-separated (±DD.MM.SS.SS)
    Example: '17:02:49.40,-48.47.23.40'
    
    Returns
    -------
    ra_deg, dec_deg : float, float
        Position in decimal degrees
    """
    try:
        # Split RA and Dec by comma
        ra_str, dec_str = position_str.split(',')
        
        # Parse RA (colon-separated: HH:MM:SS.SS)
        ra_parts = ra_str.strip().split(':')
        if len(ra_parts) != 3:
            msg(f"WARNING: Expected 3 parts in RA '{ra_str}', got {len(ra_parts)}")
            return None, None
        ra_hours = float(ra_parts[0])
        ra_mins = float(ra_parts[1])
        ra_secs = float(ra_parts[2])
        ra_deg = (ra_hours + ra_mins/60.0 + ra_secs/3600.0) * 15.0  # Convert hours to degrees
        
        # Parse Dec (period-separated: ±DD.MM.SS.SS)
        dec_str = dec_str.strip()
        dec_sign = -1 if dec_str.startswith('-') else 1
        dec_str = dec_str.lstrip('+-')  # Remove sign
        
        # Handle period-separated format
        dec_parts = dec_str.split('.')
        if len(dec_parts) < 3:
            msg(f"WARNING: Expected at least 3 parts in Dec '{dec_str}', got {len(dec_parts)}")
            return None, None
        dec_deg_part = float(dec_parts[0])
        dec_arcmin = float(dec_parts[1])
        # Reconstruct arcseconds from remaining parts
        dec_arcsec = float(dec_parts[2] + ('.' + dec_parts[3] if len(dec_parts) > 3 else ''))
        
        dec_deg = dec_sign * (abs(dec_deg_part) + dec_arcmin/60.0 + dec_arcsec/3600.0)
        
        return ra_deg, dec_deg
        
    except Exception as e:
        msg(f"ERROR parsing CASA position '{position_str}': {e}")
        return None, None


def read_position_file(position_file, source_name):
    """
    Read source position from text file for target sources.
    Format: FIELD_NAME  CASA_POSITION
    where CASA_POSITION is RA:colon-separated,Dec:period-separated
    """
    if not os.path.exists(position_file):
        msg(f"WARNING: Position file {position_file} not found")
        return None, None
    
    msg(f"Reading position from file: {position_file}")
    
    with open(position_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[0].lower() == source_name.lower():
                # Parse CASA format position
                position_str = parts[1]
                ra_deg, dec_deg = parse_casa_position(position_str)
                if ra_deg is not None:
                    msg(f"Found position in file: RA={ra_deg:.6f} deg, Dec={dec_deg:.6f} deg")
                    return ra_deg, dec_deg
    
    msg(f"WARNING: Source {source_name} not found in position file")
    return None, None


def read_rms_region_file(rms_region_file, source_name):
    """
    Read a manual RMS region override from text file, for sources whose
    default annulus (centred on the fitted position) is unsuitable -- e.g.
    known position issues or nearby confusing emission.
    Format: FIELD_NAME  CASA_REGION
    where CASA_REGION is any region string accepted by manual_rms_region
    (e.g. circle[[17:27:34.89,-16.13.19.84],75arcsec]).
    Checked before falling back to the default RMS annulus; returns None
    (not False) if the file is missing or the source has no entry, so the
    caller can decide the fallback value.
    """
    if not os.path.exists(rms_region_file):
        return None

    with open(rms_region_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split(None, 1)
            if len(parts) >= 2 and parts[0].lower() == source_name.lower():
                region = parts[1].strip()
                msg(f"Found manual RMS region in file: {region}")
                return region

    return None


def read_time_info(time_file, field_name, scan=None):
    """
    Read time information from RESULTS directory.
    Parses format:
    # MS_NAME  FIELD_NAME  SCAN  START_MJD  END_MJD
    
    Parameters
    ----------
    time_file : str
        Path to time info file
    field_name : str
        Field name to match
    scan : str, optional
        Scan number to match (for per-scan calibrators)
        
    Returns
    -------
    dict : Time information with start_mjd, end_mjd, middle_mjd, duration_hours
    """
    if not os.path.exists(time_file):
        msg(f"WARNING: Time file {time_file} not found")
        return None
    
    msg(f"Reading time info from: {time_file}")
    
    time_data = None
    with open(time_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split()
            if len(parts) >= 5:
                ms_file = parts[0]
                field_name_file = parts[1]
                scan_file = parts[2]
                start_mjd = float(parts[3])
                end_mjd = float(parts[4])
                
                # Match by field name and optionally by scan
                if field_name_file.lower() == field_name.lower():
                    if scan is None or scan == scan_file:
                        time_data = {
                            'ms_name': ms_file,
                            'field_name': field_name_file,
                            'scan': scan_file,
                            'start_mjd': start_mjd,
                            'end_mjd': end_mjd,
                            'middle_mjd': (start_mjd + end_mjd) / 2.0,
                            'duration_hours': (end_mjd - start_mjd) * 24.0
                        }
                        msg(f"Matched field: {field_name_file}, Scan: {scan_file}")
                        break
    
    if time_data:
        msg(f"  Start MJD: {time_data['start_mjd']:.10f}")
        msg(f"  End MJD: {time_data['end_mjd']:.10f}")
        msg(f"  Middle MJD: {time_data['middle_mjd']:.10f}")
        msg(f"  Duration: {time_data['duration_hours']:.4f} hours")
    else:
        msg(f"WARNING: No time info found for source {field_name}")
    
    return time_data


def get_imfit_values(fname, image, xpix, ypix):

    '''
    Run imfit on a specific image (single component only)
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

    # Single component fitting region
    x = xpix
    y = ypix
    r = 5 * bmaj
    src_region = f'circle[[{x}pix,{y}pix],{r}arcsec]'
   
    return imfit(image, estimates = fname, region = src_region)
    
    

def calculate_P0(flux_P, rms_Q, rms_U, rms_V, Aq = 0.8):
    '''
    Calculate the de-biased total polarized flux (Ptot = sqrt(Q^2 + U^2 + V^2))
    Always using total polarization (including V).
    '''    

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

    # Values of interest -- Source
    ims = imstat(image, region = src_region)
    flux = return_max(image, src_region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]

    # Extract RMS
    rms = imstat(image, region = rms_region)['rms'][0]

    return [flux, xpix, ypix, rms]
    
    
def check_position(fname, image, xpix, ypix, snr_thresh = 6.0, P_image = False, manual_rms_region = False, force_fix_to_I = False):

    '''
    Code to check whether there is sufficient flux at a position to allow 
    imfit to fit for position, or if said position should be frozen (single component)
    Inputs:
        fname = string containing name of output estimate file
        xpix = RA pixel position
        ypix = Dec pixel position
        image = name of the image to fit
        snr_thresh = 5.0
        P_image = Check if it's a P-image or not
        force_fix_to_I = If True, always fix position to Stokes I regardless of S/N
    Outputs:
        Nothing, but makes an estimate file
    '''
    
    # Get beam parameters
    bmaj = imhead(image, mode='get', hdkey='bmaj')['value']
    bmin = imhead(image, mode='get', hdkey='bmin')['value']
    bpa  = imhead(image, mode='get', hdkey='bpa')['value']

    region = f'circle[[{xpix}pix,{ypix}pix],{2 * bmaj}arcsec]'

    # Use same base name as fname for check_pos file to avoid conflicts
    check_file = fname.replace('estimate_', 'check_pos_')
    f = open(check_file, 'w')
    f.write(f'0.0,{xpix},{ypix},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, xyabp')
    f.close()

    # Get flux at test position
    test_flux = abs(imfit(image, region = region, estimates=check_file)['results']['component0']['peak']['value'])

    # Get the rms from the imstat rms calculation
    rms = get_imstat_values(image, xpix, ypix, manual_rms_region = manual_rms_region)[3]
    
    # If its a P-image don't use the image plane noise as the check criteria as it is (very) non-gaussian
    if P_image is True:
        image_Q = image.replace('-Plin-', '-Q-').replace('-Ptot-', '-Q-')
        image_U = image.replace('-Plin-', '-U-').replace('-Ptot-', '-U-')
        ims_Q = get_imstat_values(image_Q, xpix, ypix, manual_rms_region = manual_rms_region)
        ims_U = get_imstat_values(image_U, xpix, ypix, manual_rms_region = manual_rms_region)
        rms = np.amax((ims_Q[3],ims_U[3]))
        snr_thresh = PLIN_FREE_FIT_SNR  # Plin needs higher S/N for reliable free fitting
    
    # Determine if we should fit or fix position
    if force_fix_to_I:
        fix_var = 'xyabp'
        msg(f'Forcing position fix to Stokes I (non-target calibrator): {image}')
    elif test_flux > snr_thresh * rms:
        fix_var = 'abp'
        msg(f'Free fitting position with S/N = {test_flux/rms:.2f}: {image}')
    else:
        fix_var = 'xyabp'
        msg(f'Fixing position due to S/N = {test_flux/rms:.2f}: {image}')

    # Make the estimate file
    make_estimate(fname, image, xpix, ypix, fix_var, manual_rms_region = manual_rms_region)


def make_estimate(fname, image, xpix, ypix, fix_var, manual_rms_region = False):
    '''
    Take in imstat values and create an CASA imfit estimate file (single component)
    Inputs:
        fname  = string containing name of estimate file
        image  = string containing name of image to fit
        xpix   = Right ascension (pixel) estimate of source
        ypix   = Declination (pixel) estimate of source
        fix    = parameters to fix, default is assume a point source (abp) or fixing position (xyabp)
        manual_rms_region = option to specify the rms region
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']
    
    # Make estimate file
    f = open(fname, 'w')
    ims = get_imstat_values(image, xpix, ypix, manual_rms_region = manual_rms_region)
    f.write(f'{ims[0]},{xpix},{ypix},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix_var}\n')
    f.close()            
    
    return 0    


def check_fit_offset_and_refix(imfit_result, image, ref_ra_pix, ref_dec_pix, bmaj,
                               src_name, stokes_label, manual_rms_region=False,
                               is_homogenized=False):
    '''
    After an imfit call, check whether the fitted position has drifted more than
    a fraction of BMAJ (rough Euclidean pixel distance) from the reference position.
    If it has, re-fit with the position fixed to the reference and report the event.

    Parameters
    ----------
    imfit_result    : dict  - result from get_imfit_values
    image           : str   - image path (used for pixel-scale look-up via CDELT2)
    ref_ra_pix      : float - reference RA pixel (usually from Stokes I)
    ref_dec_pix     : float - reference Dec pixel (usually from Stokes I)
    bmaj            : float - beam major axis in arcsec (used as the distance threshold)
    src_name        : str   - source name (used for the temporary estimate filename)
    stokes_label    : str   - label for informative messages (e.g. "Q", "Ptot", "Plin")
    manual_rms_region : str/bool - passed through to make_estimate if a re-fit is needed
    is_homogenized  : bool  - if True, use tighter BMAJ/5 threshold (homogenized beams are uniform)

    Returns
    -------
    imfit_result : dict  - possibly updated if a re-fit was performed
    was_refitted : bool  - True if position was out of tolerance and re-fit was done
    '''
    component = 'component0'
    fitted_ra_pix  = imfit_result['results'][component]['pixelcoords'][0]
    fitted_dec_pix = imfit_result['results'][component]['pixelcoords'][1]

    # Euclidean pixel distance from the reference position
    pixel_dist = np.sqrt((fitted_ra_pix - ref_ra_pix)**2 +
                         (fitted_dec_pix - ref_dec_pix)**2)

    # Convert BMAJ (arcsec) to pixels so the threshold is physically meaningful
    # imhead returns CDELT2 in radians, so convert rad -> deg -> arcsec
    cell_rad    = abs(imhead(image, mode='get', hdkey='CDELT2')['value'])
    cell_arcsec = np.degrees(cell_rad) * 3600.0
    bmaj_pixels = bmaj / cell_arcsec
    offset_arcsec = pixel_dist * cell_arcsec

    # Homogenized images have uniform beams so allow less drift (BMAJ/5 vs BMAJ/3)
    if is_homogenized:
        threshold_fraction = 1.0 / 5.0
        threshold_label = '1/5'
    else:
        threshold_fraction = 1.0 / 3.0
        threshold_label = '1/3'

    if pixel_dist > threshold_fraction * bmaj_pixels:
        msg(f'  WARNING: Stokes {stokes_label} fitted position offset = {offset_arcsec:.2f}" '
            f'({pixel_dist:.1f} pix) exceeds {threshold_label} BMAJ ({threshold_fraction * bmaj:.2f}"). '
            f'Re-fitting with position fixed to reference.')
        fixed_est = f'estimate_fixed_{stokes_label}_{src_name}.txt'
        make_estimate(fixed_est, image, ref_ra_pix, ref_dec_pix, 'xyabp',
                      manual_rms_region=manual_rms_region)
        imfit_result = get_imfit_values(fixed_est, image, ref_ra_pix, ref_dec_pix)
        return imfit_result, True
    else:
        msg(f'  Position check OK: Stokes {stokes_label} offset = {offset_arcsec:.2f}" '
            f'({pixel_dist:.1f} pix) < {threshold_label} BMAJ ({threshold_fraction * bmaj:.2f}")')
        return imfit_result, False


def parse_channel_number(filename):
    """Extract the 4-digit channel number from a WSCLEAN filename."""
    match = re.search(r'-(\d{4})-(?:[IQUV]|Plin|Ptot)-', filename)
    if match:
        return match.group(1)
    # I-only images: prefix-NNNN-image.fits (no Stokes label)
    match = re.search(r'-(\d{4})-image', filename)
    if match:
        return match.group(1)
    return None


def build_channel_map(image_list):
    """Map channel-number string -> filepath for WSCLEAN channelised images."""
    cmap = {}
    for im in image_list:
        ch = parse_channel_number(im)
        if ch is not None:
            cmap[ch] = im
    return cmap


def fit_channel(job):
    '''
    Fit Stokes I(/P/Q/U/V) for one channel-image group. Runs the same fitting
    sequence as the old serial CHAN loop body, but is fully self-contained (no
    output_dictionary access) so it can run inside a parallel worker: all MFS
    reference positions/metadata come in via `job`, and results go out as a
    plain dict rather than being appended in place. Estimate-file names carry
    the channel number (unlike the old serial loop, which reused one fixed
    name per Stokes) since multiple channels may now run concurrently.

    `job` keys: ch_num, im_I, im_Q, im_U, im_V, im_P, stokes_I_only,
    pol_image_type, component, MFS_I_ra_pix, MFS_I_dec_pix, MFS_bmaj,
    MFS_bmin, MFS_freq_GHz, manual_rms_region, force_fix_to_I,
    force_fix_V_to_I, chan_is_homogenized, src_name.

    Returns (ch_num, result) where result is None if the channel was skipped
    (flagged/anomalous beam) or the fit raised, otherwise a dict of the
    per-channel scalars matching output_dictionary['CHAN']'s flat keys.
    '''
    ch_num               = job['ch_num']
    im_I                 = job['im_I']
    im_Q                 = job['im_Q']
    im_U                 = job['im_U']
    im_V                 = job['im_V']
    im_P                 = job['im_P']
    stokes_I_only        = job['stokes_I_only']
    pol_image_type       = job['pol_image_type']
    component            = job['component']
    MFS_I_ra_pix         = job['MFS_I_ra_pix']
    MFS_I_dec_pix        = job['MFS_I_dec_pix']
    MFS_bmaj             = job['MFS_bmaj']
    MFS_bmin             = job['MFS_bmin']
    MFS_freq_GHz         = job['MFS_freq_GHz']
    manual_rms_region    = job['manual_rms_region']
    force_fix_to_I       = job['force_fix_to_I']
    force_fix_V_to_I     = job['force_fix_V_to_I']
    chan_is_homogenized  = job['chan_is_homogenized']
    src_name             = job['src_name']

    try:
        # Get frequency and beam info
        freq_GHz = imhead(im_I, mode='get', hdkey='CRVAL3')['value'] / 1.0e9
        bmaj = imhead(im_I, mode='get', hdkey='bmaj')['value']
        bmin = imhead(im_I, mode='get', hdkey='bmin')['value']
        bpa = imhead(im_I, mode='get', hdkey='bpa')['value']

        msg('')
        msg(f'--- Channel {ch_num} | {freq_GHz:.4f} GHz ---')

        # Check if beam is valid (not flagged channel)
        if bmaj == 0.0 or bmin == 0.0:
            msg(f'Channel {ch_num} appears to be flagged (zero beam), skipping')
            return ch_num, None

        # Check beam scaling - skip if beam is >10x larger than expected from frequency scaling
        # Expected beam scales as freq_MFS / freq_channel
        expected_bmaj = MFS_bmaj * (MFS_freq_GHz / freq_GHz)
        expected_bmin = MFS_bmin * (MFS_freq_GHz / freq_GHz)

        if bmaj > 10.0 * expected_bmaj or bmin > 10.0 * expected_bmin:
            msg(f'Channel {ch_num} has anomalous beam '
                f'(bmaj={bmaj:.2f}, bmin={bmin:.2f} vs '
                f'expected {expected_bmaj:.2f}, {expected_bmin:.2f}), skipping')
            return ch_num, None

        # Fit Stokes I
        check_position(f'estimate_ch_I_{src_name}_{ch_num}.txt', im_I,
                      MFS_I_ra_pix, MFS_I_dec_pix, P_image=False,
                      manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
        ch_I_imfit = get_imfit_values(f'estimate_ch_I_{src_name}_{ch_num}.txt',
                                      im_I, MFS_I_ra_pix, MFS_I_dec_pix)

        flux_I = ch_I_imfit['results'][component]['peak']['value'] * 1e3
        err_I = ch_I_imfit['results'][component]['peak']['error'] * 1e3
        RA_I = ch_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
        DEC_I = ch_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
        RA_pix_ch = ch_I_imfit['results'][component]['pixelcoords'][0]
        DEC_pix_ch = ch_I_imfit['results'][component]['pixelcoords'][1]
        rms_I = get_imstat_values(im_I, RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3

        result = {'freq_GHz': freq_GHz, 'bmaj_asec': bmaj, 'bmin_asec': bmin, 'bpa_deg': bpa,
                  'I_flux_mJy': flux_I, 'I_err_mJy': err_I, 'I_rms_mJy': rms_I,
                  'I_RA_deg': RA_I, 'I_DEC_deg': DEC_I}

        if not stokes_I_only:
            # Fit P(lin/tot) first to check significance
            check_position(f'estimate_ch_P_{src_name}_{ch_num}.txt', im_P,
                          RA_pix_ch, DEC_pix_ch, P_image=True,
                          manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
            ch_P_imfit = get_imfit_values(f'estimate_ch_P_{src_name}_{ch_num}.txt',
                                          im_P, RA_pix_ch, DEC_pix_ch)
            ch_P_imfit, _ = check_fit_offset_and_refix(
                ch_P_imfit, im_P, RA_pix_ch, DEC_pix_ch,
                bmaj, src_name, pol_image_type, manual_rms_region=manual_rms_region,
                is_homogenized=chan_is_homogenized)
            flux_P = ch_P_imfit['results'][component]['peak']['value'] * 1e3
            err_P = ch_P_imfit['results'][component]['peak']['error'] * 1e3
            rms_P = get_imstat_values(im_P, RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3

            # Check if P is significant enough to use as Q/U position reference
            P_snr_ch = flux_P / rms_P
            frac_pol_ch = abs(flux_P / flux_I) if flux_I != 0 else 0.0
            if not force_fix_to_I and P_snr_ch > CHAN_P_SNR_THRESH and frac_pol_ch > MIN_FRAC_POL:
                ref_ra_pix_ch = ch_P_imfit['results'][component]['pixelcoords'][0]
                ref_dec_pix_ch = ch_P_imfit['results'][component]['pixelcoords'][1]
                msg(f'  Ch {ch_num}: {pol_image_type} SNR={P_snr_ch:.1f} > {CHAN_P_SNR_THRESH} '
                    f'and frac pol={frac_pol_ch*100:.2f}% > {MIN_FRAC_POL*100:.1f}% '
                    f'-> Q/U ref = {pol_image_type} position ({ref_ra_pix_ch:.2f}, {ref_dec_pix_ch:.2f}) pix')
            else:
                ref_ra_pix_ch = RA_pix_ch
                ref_dec_pix_ch = DEC_pix_ch
                if force_fix_to_I:
                    msg(f'  Ch {ch_num}: forced fix -> Q/U ref = Stokes I ({ref_ra_pix_ch:.2f}, {ref_dec_pix_ch:.2f}) pix')
                elif P_snr_ch <= CHAN_P_SNR_THRESH:
                    msg(f'  Ch {ch_num}: {pol_image_type} SNR={P_snr_ch:.1f} <= {CHAN_P_SNR_THRESH} '
                        f'-> Q/U ref = Stokes I ({ref_ra_pix_ch:.2f}, {ref_dec_pix_ch:.2f}) pix')
                else:
                    msg(f'  Ch {ch_num}: frac pol={frac_pol_ch*100:.2f}% <= {MIN_FRAC_POL*100:.1f}% '
                        f'-> Q/U ref = Stokes I ({ref_ra_pix_ch:.2f}, {ref_dec_pix_ch:.2f}) pix')

            # Fit Q
            check_position(f'estimate_ch_Q_{src_name}_{ch_num}.txt', im_Q,
                          ref_ra_pix_ch, ref_dec_pix_ch, P_image=False,
                          manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
            ch_Q_imfit = get_imfit_values(f'estimate_ch_Q_{src_name}_{ch_num}.txt',
                                          im_Q, ref_ra_pix_ch, ref_dec_pix_ch)
            ch_Q_imfit, _ = check_fit_offset_and_refix(
                ch_Q_imfit, im_Q, ref_ra_pix_ch, ref_dec_pix_ch,
                bmaj, src_name, 'Q', manual_rms_region=manual_rms_region,
                is_homogenized=chan_is_homogenized)
            flux_Q = ch_Q_imfit['results'][component]['peak']['value'] * 1e3
            err_Q = ch_Q_imfit['results'][component]['peak']['error'] * 1e3
            rms_Q = get_imstat_values(im_Q, RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3

            # Fit U
            check_position(f'estimate_ch_U_{src_name}_{ch_num}.txt', im_U,
                          ref_ra_pix_ch, ref_dec_pix_ch, P_image=False,
                          manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
            ch_U_imfit = get_imfit_values(f'estimate_ch_U_{src_name}_{ch_num}.txt',
                                          im_U, ref_ra_pix_ch, ref_dec_pix_ch)
            ch_U_imfit, _ = check_fit_offset_and_refix(
                ch_U_imfit, im_U, ref_ra_pix_ch, ref_dec_pix_ch,
                bmaj, src_name, 'U', manual_rms_region=manual_rms_region,
                is_homogenized=chan_is_homogenized)
            flux_U = ch_U_imfit['results'][component]['peak']['value'] * 1e3
            err_U = ch_U_imfit['results'][component]['peak']['error'] * 1e3
            rms_U = get_imstat_values(im_U, RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3

            # Fit V - ALWAYS referenced to Stokes I
            check_position(f'estimate_ch_V_{src_name}_{ch_num}.txt', im_V,
                          RA_pix_ch, DEC_pix_ch, P_image=False,
                          manual_rms_region=manual_rms_region,
                          force_fix_to_I=(force_fix_to_I or force_fix_V_to_I))
            ch_V_imfit = get_imfit_values(f'estimate_ch_V_{src_name}_{ch_num}.txt',
                                          im_V, RA_pix_ch, DEC_pix_ch)
            ch_V_imfit, _ = check_fit_offset_and_refix(
                ch_V_imfit, im_V, RA_pix_ch, DEC_pix_ch,
                bmaj, src_name, 'V', manual_rms_region=manual_rms_region,
                is_homogenized=chan_is_homogenized)
            flux_V = ch_V_imfit['results'][component]['peak']['value'] * 1e3
            err_V = ch_V_imfit['results'][component]['peak']['error'] * 1e3
            rms_V = get_imstat_values(im_V, RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3

            # Compute Euclidean distance of each fitted position from channelised Stokes I (arcsec)
            cell_rad_ch = abs(imhead(im_I, mode='get', hdkey='CDELT2')['value'])
            cell_asec_ch = np.degrees(cell_rad_ch) * 3600.0
            P_pix_ra = ch_P_imfit['results'][component]['pixelcoords'][0]
            P_pix_dec = ch_P_imfit['results'][component]['pixelcoords'][1]
            Q_pix_ra = ch_Q_imfit['results'][component]['pixelcoords'][0]
            Q_pix_dec = ch_Q_imfit['results'][component]['pixelcoords'][1]
            U_pix_ra = ch_U_imfit['results'][component]['pixelcoords'][0]
            U_pix_dec = ch_U_imfit['results'][component]['pixelcoords'][1]
            V_pix_ra = ch_V_imfit['results'][component]['pixelcoords'][0]
            V_pix_dec = ch_V_imfit['results'][component]['pixelcoords'][1]
            dist_P = np.sqrt((P_pix_ra - RA_pix_ch)**2 + (P_pix_dec - DEC_pix_ch)**2) * cell_asec_ch
            dist_Q = np.sqrt((Q_pix_ra - RA_pix_ch)**2 + (Q_pix_dec - DEC_pix_ch)**2) * cell_asec_ch
            dist_U = np.sqrt((U_pix_ra - RA_pix_ch)**2 + (U_pix_dec - DEC_pix_ch)**2) * cell_asec_ch
            dist_V = np.sqrt((V_pix_ra - RA_pix_ch)**2 + (V_pix_dec - DEC_pix_ch)**2) * cell_asec_ch

            msg(f"  Flux densities:")
            msg(f"    I:    {flux_I:+10.4f} +/- {rms_I:.4f} mJy  (S/N = {flux_I/rms_I:.1f})")
            msg(f"    {pol_image_type+':':5s} {flux_P:+10.4f} +/- {rms_P:.4f} mJy  (S/N = {flux_P/rms_P:.1f}, dI = {dist_P:.2f}\")")
            msg(f"    Q:    {flux_Q:+10.4f} +/- {rms_Q:.4f} mJy  (S/N = {flux_Q/rms_Q:.1f}, dI = {dist_Q:.2f}\")")
            msg(f"    U:    {flux_U:+10.4f} +/- {rms_U:.4f} mJy  (S/N = {flux_U/rms_U:.1f}, dI = {dist_U:.2f}\")")
            msg(f"    V:    {flux_V:+10.4f} +/- {rms_V:.4f} mJy  (S/N = {flux_V/rms_V:.1f}, dI = {dist_V:.2f}\")")

            result.update({
                'Q_flux_mJy': flux_Q, 'U_flux_mJy': flux_U, 'V_flux_mJy': flux_V,
                f'{pol_image_type}_flux_mJy': flux_P,
                'Q_err_mJy': err_Q, 'U_err_mJy': err_U, 'V_err_mJy': err_V,
                f'{pol_image_type}_err_mJy': err_P,
                'Q_rms_mJy': rms_Q, 'U_rms_mJy': rms_U, 'V_rms_mJy': rms_V,
                f'{pol_image_type}_rms_mJy': rms_P,
            })
        else:
            msg(f"  I: {flux_I:+10.4f} +/- {rms_I:.4f} mJy  (S/N = {flux_I/rms_I:.1f})")

        return ch_num, result

    except Exception as e:
        msg(f'Fitting failed for channel {ch_num}: {e}')
        msg('Channel is likely flagged, skipping')
        return ch_num, None


def fit_channel_batch(jobs):
    '''
    Pool entry point: runs in a freshly spawned worker process. imfit/imstat/
    imhead are plain top-level imports (see top of file), and
    multiprocessing's spawn context re-executes this module's top-level code
    in every child process same as any other Python module, so those names
    (and the CHAN_P_SNR_THRESH/MIN_FRAC_POL globals) resolve normally with no
    extra binding needed here. Processes its whole assigned chunk of channels
    in one call.
    '''
    return [fit_channel(job) for job in jobs]


def extract_polarization_properties(src_name,
    src_im_identifier,
    src_im_suffix, 
    src_ra, 
    src_dec, 
    src_ulims,
    manual_rms_region, 
    image_directory,
    image_identifier,
    time_info = None,
    force_fix_to_I = False,
    force_fix_V_to_I = False,
    use_plin = False):

    '''
    Fit the Stokes IQUV cube for single component in an image. This assumes
    that the Q, U, and V images follow the WSCLEAN naming
    convention.

    input parameters:
        src_name       = name of source
        src_im_identifier  = identifier for images that will have flux extraction
        src_im_suffix  = image suffix, included to differentiate between standard WSCLEAN image products (image.fits) and homogenized products (image.homogenized.fits)
        src_ra   = Estimated right ascension of source (degrees)
        src_dec  = Estimated declination of  source (degrees)
        manual_rms_region = option to specify the rms region, otherwise use an annulus centered on the source(s) with a ~100xPSF area
        time_info = dictionary with timing information from ms_time_info.txt
        force_fix_to_I = If True, always fix all Stokes positions to Stokes I (for non-target calibrators)
        force_fix_V_to_I = If True, force only Stokes V position to Stokes I regardless of force_fix_to_I
        use_plin = If True, use Plin (linear: sqrt(Q^2+U^2)) instead of Ptot (sqrt(Q^2+U^2+V^2)).
                   Should be set True when polang_name is present in project_info.

    Output parametres:
        output_dictionary = Dictionary containing all MFS and per-channel information for the IQUV fluxes
    ''' 

    # Initialize dictionary to contain the output parameters
    output_dictionary = {'name' : src_name}
    
    # Add time information if available
    if time_info:
        output_dictionary['timing'] = time_info
    
    output_dictionary['MFS'] = {}
    output_dictionary['CHAN'] = {}

    # Determine which polarized intensity image type to use.
    # Plin = sqrt(Q^2+U^2) is preferred when a polarization angle calibrator is present
    # because V is not used in angle calibration and Plin has lower noise.
    # Ptot = sqrt(Q^2+U^2+V^2) is the more general total polarized intensity.
    # If Ptot was requested but no Ptot images exist on disk (e.g. make_pol_images.py
    # skipped Ptot because Stokes V was unavailable), fall back to Plin instead of
    # failing downstream.
    plin_fallback = False
    if not use_plin and not glob.glob(f'{src_im_identifier}*-Ptot-*'):
        if glob.glob(f'{src_im_identifier}*-Plin-*'):
            msg('WARNING: Ptot requested but no Ptot images found on disk -- falling back to Plin')
            use_plin = True
            plin_fallback = True
        else:
            msg('WARNING: Neither Ptot nor Plin images found on disk for this source')

    pol_image_type    = 'Plin' if use_plin else 'Ptot'
    pol_image_exclude = '-Ptot-' if use_plin else '-Plin-'  # Filter out the non-selected type
    if plin_fallback:
        msg(f'Polarization image type: Plin (linear only) -- fallback, no Ptot images on disk')
    elif use_plin:
        msg(f'Polarization image type: Plin (linear only) -- polang_name is set, V excluded from P estimate')
    else:
        msg(f'Polarization image type: Ptot (total, including V)')
    output_dictionary['pol_image_type'] = pol_image_type

    # Check for MFS images: prioritize image.fits, fall back to image.homogenized.fits
    msg("Checking for MFS images...")
    mfs_standard = glob.glob(f'{src_im_identifier}MFS*image.fits')
    mfs_homogenized = glob.glob(f'{src_im_identifier}MFS*image.homogenized.fits')
    
    if mfs_standard:
        mfs_im_suffix = 'image.fits'
        msg(f"Using standard MFS images: {mfs_im_suffix}")
    elif mfs_homogenized:
        mfs_im_suffix = 'image.homogenized.fits'
        msg(f"Standard MFS images not found, using homogenized MFS images: {mfs_im_suffix}")
    else:
        msg(f"WARNING: No MFS images found with pattern: {src_im_identifier}MFS*")
        mfs_im_suffix = None
    
    mfs_is_homogenized = (mfs_im_suffix is not None and 'homogenized' in mfs_im_suffix)
    
    # Check for channelised images: prioritize image.fits, fall back to image.homogenized.fits
    # Auto-detect Stokes-I-only mode: if images have -I- in the name (e.g. -0042-I-image.fits)
    # it's IQUV. If they lack a Stokes label (e.g. -0042-image.fits) it's I-only.
    msg("Checking for channelised images...")
    stokes_I_only = False

    chan_standard = glob.glob(f'{src_im_identifier}*-I-image.fits')
    chan_standard = [im for im in chan_standard if '-MFS-' not in im]
    chan_homogenized = glob.glob(f'{src_im_identifier}*-I-image.homogenized.fits')
    chan_homogenized = [im for im in chan_homogenized if '-MFS-' not in im]
    
    if chan_standard:
        src_im_suffix = 'image.fits'
        msg(f"Using standard channelised IQUV images: {src_im_suffix}")
    elif chan_homogenized:
        src_im_suffix = 'image.homogenized.fits'
        msg(f"Using homogenized channelised IQUV images: {src_im_suffix}")
    else:
        # No -I- images found -- check for Stokes-I-only images (no Stokes label)
        chan_standard_ionly = glob.glob(f'{src_im_identifier}*-image.fits')
        chan_standard_ionly = [im for im in chan_standard_ionly
                              if '-MFS-' not in im and not re.search(r'-[IQUV]-image', im)
                              and '-Plin-' not in im and '-Ptot-' not in im]
        chan_homog_ionly = glob.glob(f'{src_im_identifier}*-image.homogenized.fits')
        chan_homog_ionly = [im for im in chan_homog_ionly
                           if '-MFS-' not in im and not re.search(r'-[IQUV]-image', im)
                           and '-Plin-' not in im and '-Ptot-' not in im]
        
        if chan_standard_ionly:
            src_im_suffix = 'image.fits'
            chan_standard = chan_standard_ionly
            stokes_I_only = True
            msg(f"Stokes-I-only images detected (no Stokes label): {src_im_suffix}")
        elif chan_homog_ionly:
            src_im_suffix = 'image.homogenized.fits'
            chan_homogenized = chan_homog_ionly
            stokes_I_only = True
            msg(f"Stokes-I-only homogenized images detected: {src_im_suffix}")
        else:
            msg(f"ERROR: No channelised images found")
            raise FileNotFoundError(
                f"No channelised images found with pattern: "
                f"{src_im_identifier}*-[I-]image[.homogenized].fits")
    
    # Now that src_im_suffix is finalised, determine if channelised images are homogenized
    chan_is_homogenized = 'homogenized' in src_im_suffix
    
    if stokes_I_only:
        msg('Stokes-I-only mode: skipping Q, U, V, P extraction')
        output_dictionary.pop('pol_image_type', None)
    
    # Get unique prefix from channelised images
    prefix_arr = chan_standard if chan_standard else chan_homogenized
    example_img = prefix_arr[0]
    if stokes_I_only:
        # Pattern: prefix-NNNN-image.fits
        match = re.search(r'(.+)-\d{4}-image', example_img)
    else:
        # Pattern: prefix-NNNN-I-image.fits
        match = re.search(r'(.+)-\d{4}-[IQUV]-', example_img)
    if match:
        prefix = match.group(1)
    else:
        prefix = example_img.rsplit('-', 3)[0]
    
    # Extract the MFS image parameters if they exist
    if mfs_im_suffix:
        msg(f'Fitting MFS image(s) for prefix: {prefix}')
        
        # Collect MFS images; order after sort will be I, P(lin/tot), Q, U, V.
        # Exclude whichever polarized intensity type is NOT being used (controlled by use_plin).
        MFS_images = glob.glob(f'{prefix}-MFS-*-{mfs_im_suffix}')
        if stokes_I_only:
            MFS_images = sorted([im for im in MFS_images if '-I-' in im])
            if not MFS_images:
                # I-only MFS: no Stokes label (prefix-MFS-image.fits)
                MFS_images = sorted(glob.glob(f'{prefix}-MFS-{mfs_im_suffix}'))
            msg(f'MFS images (Stokes I only): {[os.path.basename(im) for im in MFS_images]}')
        else:
            MFS_images = sorted([im for im in MFS_images if pol_image_exclude not in im])
            msg(f'MFS images selected (using {pol_image_type}, excluding {pol_image_exclude.strip("-")}): {[os.path.basename(im) for im in MFS_images]}')

        # Output image name
        for MFS_image in MFS_images:
            msg(f'Fitting MFS image name(s): {MFS_image}')

        # Get generalized properties from the first image (i.e., Stokes I header)
        MFS_freq_GHz = imhead(MFS_images[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
        MFS_bmaj = imhead(MFS_images[0], mode='get', hdkey = 'bmaj')['value']
        MFS_bmin = imhead(MFS_images[0], mode='get', hdkey = 'bmin')['value']
        MFS_bpa   = imhead(MFS_images[0], mode='get', hdkey = 'bpa')['value']

        # Convert the RA/DEC guess to pixel coordinates (single component)
        region = 'circle[[{}deg,{}deg],1.0pix]'.format(src_ra, src_dec)
        pixel_imstat = imstat(MFS_images[0], region = region)
        src_ra_pix = pixel_imstat['maxpos'][0]
        src_dec_pix = pixel_imstat['maxpos'][1]
            
        # First fit Stokes I -- use check_position to determine if position should be fixed
        check_position(f'estimate_I_{src_name}.txt', MFS_images[0], 
                      src_ra_pix, src_dec_pix, P_image=False, 
                      manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
        MFS_I_imfit = get_imfit_values(f'estimate_I_{src_name}.txt', 
                                       MFS_images[0], src_ra_pix, src_dec_pix)

        # Get component (single component only)
        component = 'component0'           
        MFS_I_ra_pix = MFS_I_imfit['results'][component]['pixelcoords'][0]
        MFS_I_dec_pix = MFS_I_imfit['results'][component]['pixelcoords'][1]
       
        # Extract values for single component
        flux_I = MFS_I_imfit['results'][component]['peak']['value'] * 1e3
        err_I = MFS_I_imfit['results'][component]['peak']['error'] * 1e3
        RA_I   = MFS_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
        DEC_I = MFS_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
        RA_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][0]
        DEC_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][1]
        rms_I = get_imstat_values(MFS_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
        
        msg('')
        msg(f"--- MFS Stokes I Fitted Position (Reference) ---")
        msg(f"RA:  {RA_I:.6f} deg  ({RA_pix_I:.2f} pix)")
        msg(f"Dec: {DEC_I:.6f} deg ({DEC_pix_I:.2f} pix)")
        msg(f"Flux: {flux_I:.3f} +/- {rms_I:.3f} mJy (S/N = {flux_I/rms_I:.1f})")
        msg('')
       
        # Initialize single-value MFS dictionary (not lists)
        output_dictionary['MFS'] = {}
        output_dictionary['MFS']['upperlimit'] = src_ulims[0]           
        output_dictionary['MFS']['freq_GHz'] = MFS_freq_GHz
        output_dictionary['MFS']['bmaj_asec'] = MFS_bmaj
        output_dictionary['MFS']['bmin_asec'] = MFS_bmin
        output_dictionary['MFS']['bpa_deg'] = MFS_bpa
        output_dictionary['MFS']['I_flux_mJy'] = flux_I
        output_dictionary['MFS']['I_err_mJy'] = err_I
        output_dictionary['MFS']['I_rms_mJy'] = rms_I
        output_dictionary['MFS']['I_RA_deg'] = RA_I
        output_dictionary['MFS']['I_DEC_deg'] = DEC_I
        
        # Next fit Q, U, V, and P(lin/tot) -- which P image is used is set by pol_image_type.
        if not stokes_I_only:
            # Fit P first (image 1 after I) - initially referenced to I
            check_position(f'estimate_P_{src_name}.txt', MFS_images[1], 
                          MFS_I_ra_pix, MFS_I_dec_pix, P_image=True, 
                          manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
            MFS_P_imfit  = get_imfit_values(f'estimate_P_{src_name}.txt', 
                                            MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
            # Check fit has not wandered more than 0.5 BMAJ from I reference
            MFS_P_imfit, _ = check_fit_offset_and_refix(
                MFS_P_imfit, MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix,
                MFS_bmaj, src_name, pol_image_type, manual_rms_region=manual_rms_region,
                is_homogenized=mfs_is_homogenized)
            
            # Get P flux and RMS to check significance
            flux_P = MFS_P_imfit['results'][component]['peak']['value'] * 1e3
            err_P = MFS_P_imfit['results'][component]['peak']['error'] * 1e3
            rms_P = get_imstat_values(MFS_images[1], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
            
            # Check if P is significant enough to use as Q/U position reference
            P_snr = flux_P / rms_P
            frac_pol = abs(flux_P / flux_I) if flux_I != 0 else 0.0
            if not force_fix_to_I and P_snr > MFS_P_SNR_THRESH and frac_pol > MIN_FRAC_POL:
                msg(f'{pol_image_type} has high SNR ({P_snr:.1f}) and frac pol ({frac_pol*100:.2f}%), '
                    f'referencing Q/U against {pol_image_type} position')
                # Get P position to use as reference for Q/U
                MFS_P_ra_pix = MFS_P_imfit['results'][component]['pixelcoords'][0]
                MFS_P_dec_pix = MFS_P_imfit['results'][component]['pixelcoords'][1]
                ref_ra_pix = MFS_P_ra_pix
                ref_dec_pix = MFS_P_dec_pix
            else:
                if force_fix_to_I:
                    msg(f'Referencing Q/U against Stokes I position (forced for non-target calibrator)')
                elif P_snr <= MFS_P_SNR_THRESH:
                    msg(f'{pol_image_type} has low SNR ({P_snr:.1f} <= {MFS_P_SNR_THRESH}), '
                        f'referencing Q/U against Stokes I position')
                else:
                    msg(f'{pol_image_type} frac pol too low ({frac_pol*100:.2f}% <= {MIN_FRAC_POL*100:.1f}%), '
                        f'referencing Q/U against Stokes I position')
                ref_ra_pix = MFS_I_ra_pix
                ref_dec_pix = MFS_I_dec_pix
            
            # Fit Q (image 2) - referenced to chosen position
            check_position(f'estimate_Q_{src_name}.txt', MFS_images[2], 
                          ref_ra_pix, ref_dec_pix, P_image=False, 
                          manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
            MFS_Q_imfit  = get_imfit_values(f'estimate_Q_{src_name}.txt', 
                                            MFS_images[2], ref_ra_pix, ref_dec_pix)
            # Check fit has not wandered more than 0.5 BMAJ from reference
            MFS_Q_imfit, _ = check_fit_offset_and_refix(
                MFS_Q_imfit, MFS_images[2], ref_ra_pix, ref_dec_pix,
                MFS_bmaj, src_name, 'Q', manual_rms_region=manual_rms_region,
                is_homogenized=mfs_is_homogenized)

            # Fit U (image 3) - referenced to chosen position
            check_position(f'estimate_U_{src_name}.txt', MFS_images[3], 
                          ref_ra_pix, ref_dec_pix, P_image=False, 
                          manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
            MFS_U_imfit  = get_imfit_values(f'estimate_U_{src_name}.txt', 
                                            MFS_images[3], ref_ra_pix, ref_dec_pix)
            # Check fit has not wandered more than 0.5 BMAJ from reference
            MFS_U_imfit, _ = check_fit_offset_and_refix(
                MFS_U_imfit, MFS_images[3], ref_ra_pix, ref_dec_pix,
                MFS_bmaj, src_name, 'U', manual_rms_region=manual_rms_region,
                is_homogenized=mfs_is_homogenized)

            # Fit V (image 4) - ALWAYS referenced to Stokes I regardless of P/Q/U reference choice.
            # V is not used in linear polarization analysis and is always checked against I.
            check_position(f'estimate_V_{src_name}.txt', MFS_images[4],
                          MFS_I_ra_pix, MFS_I_dec_pix, P_image=False,
                          manual_rms_region=manual_rms_region,
                          force_fix_to_I=(force_fix_to_I or force_fix_V_to_I))
            MFS_V_imfit  = get_imfit_values(f'estimate_V_{src_name}.txt', 
                                            MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix)
            # Check fit has not wandered more than 0.5 BMAJ from Stokes I reference
            MFS_V_imfit, _ = check_fit_offset_and_refix(
                MFS_V_imfit, MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix,
                MFS_bmaj, src_name, 'V', manual_rms_region=manual_rms_region,
                is_homogenized=mfs_is_homogenized)

            # Extract remaining flux values and errors
            flux_Q = MFS_Q_imfit['results'][component]['peak']['value'] * 1e3
            flux_U = MFS_U_imfit['results'][component]['peak']['value'] * 1e3
            flux_V = MFS_V_imfit['results'][component]['peak']['value'] * 1e3
            
            err_P = MFS_P_imfit['results'][component]['peak']['error'] * 1e3
            err_Q = MFS_Q_imfit['results'][component]['peak']['error'] * 1e3
            err_U = MFS_U_imfit['results'][component]['peak']['error'] * 1e3
            err_V = MFS_V_imfit['results'][component]['peak']['error'] * 1e3

            # Get remaining RMS values (P already computed above)
            rms_Q = get_imstat_values(MFS_images[2], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
            rms_U = get_imstat_values(MFS_images[3], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
            rms_V = get_imstat_values(MFS_images[4], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
            
            # Extract and display fitted positions for verification
            RA_Q = MFS_Q_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_Q = MFS_Q_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            RA_U = MFS_U_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_U = MFS_U_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            RA_V = MFS_V_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_V = MFS_V_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            
            msg('')
            msg(f"--- MFS Fitted Positions Summary ---")
            msg(f"Stokes Q: RA={RA_Q:.6f} deg, Dec={DEC_Q:.6f} deg (S/N={flux_Q/rms_Q:.1f})")
            msg(f"Stokes U: RA={RA_U:.6f} deg, Dec={DEC_U:.6f} deg (S/N={flux_U/rms_U:.1f})")
            msg(f"Stokes V: RA={RA_V:.6f} deg, Dec={DEC_V:.6f} deg (S/N={flux_V/rms_V:.1f})")
            
            # Calculate offsets from Stokes I to verify position fixing
            # Apply cos(dec) correction to RA term for proper angular separation
            cosdec = np.cos(np.radians(DEC_I))
            offset_Q = np.sqrt(((RA_Q - RA_I) * cosdec)**2 + (DEC_Q - DEC_I)**2) * 3600  # arcsec
            offset_U = np.sqrt(((RA_U - RA_I) * cosdec)**2 + (DEC_U - DEC_I)**2) * 3600  # arcsec
            offset_V = np.sqrt(((RA_V - RA_I) * cosdec)**2 + (DEC_V - DEC_I)**2) * 3600  # arcsec
            
            msg('')
            msg(f"Position offsets from Stokes I:")
            msg(f"Q offset: {offset_Q:.4f} arcsec")
            msg(f"U offset: {offset_U:.4f} arcsec")
            msg(f"V offset: {offset_V:.4f} arcsec")
            if force_fix_to_I:
                msg(f"(Should be ~0 since positions are fixed to Stokes I)")
            msg('')
            msg(f"--- MFS Flux Densities ---")
            msg(f"  I:    {flux_I:+10.4f} +/- {rms_I:.4f} mJy  (S/N = {flux_I/rms_I:.1f})")
            msg(f"  {pol_image_type+':':5s} {flux_P:+10.4f} +/- {rms_P:.4f} mJy  (S/N = {flux_P/rms_P:.1f})")
            msg(f"  Q:    {flux_Q:+10.4f} +/- {rms_Q:.4f} mJy  (S/N = {flux_Q/rms_Q:.1f})")
            msg(f"  U:    {flux_U:+10.4f} +/- {rms_U:.4f} mJy  (S/N = {flux_U/rms_U:.1f})")
            msg(f"  V:    {flux_V:+10.4f} +/- {rms_V:.4f} mJy  (S/N = {flux_V/rms_V:.1f})")
            msg('')

            # Store in output dictionary
            output_dictionary['MFS'][f'{pol_image_type}_flux_mJy'] = flux_P
            output_dictionary['MFS']['Q_flux_mJy'] = flux_Q
            output_dictionary['MFS']['U_flux_mJy'] = flux_U
            output_dictionary['MFS']['V_flux_mJy'] = flux_V
            
            output_dictionary['MFS'][f'{pol_image_type}_err_mJy'] = err_P
            output_dictionary['MFS']['Q_err_mJy'] = err_Q
            output_dictionary['MFS']['U_err_mJy'] = err_U
            output_dictionary['MFS']['V_err_mJy'] = err_V

            output_dictionary['MFS'][f'{pol_image_type}_rms_mJy'] = rms_P
            output_dictionary['MFS']['Q_rms_mJy'] = rms_Q
            output_dictionary['MFS']['U_rms_mJy'] = rms_U
            output_dictionary['MFS']['V_rms_mJy'] = rms_V
    else:
        msg('Skipping MFS image processing (no MFS images found)')

    # Now process channelized data (simplified for single component, 1D arrays)
    msg('')
    msg(f"{'='*80}")
    msg('Extracting channelized data')
    msg(f"{'='*80}")
    msg('')
    
    # Find channelized images (not MFS) and build channel maps
    if stokes_I_only:
        # I-only: images may or may not have -I- in name
        images_I_all = sorted([im for im in glob.glob(f'{prefix}-*-I-{src_im_suffix}') if '-MFS-' not in im])
        if not images_I_all:
            images_I_all = sorted([im for im in glob.glob(f'{prefix}-*-{src_im_suffix}') if '-MFS-' not in im
                                   and not re.search(r'-[IQUV]-', im) and '-Plin-' not in im and '-Ptot-' not in im])
        cmap_I = build_channel_map(images_I_all)
        valid_channels = sorted(cmap_I.keys())
        msg(f'Found {len(valid_channels)} channelized Stokes I images (I-only mode)')
    else:
        images_I_all = sorted([im for im in glob.glob(f'{prefix}-*-I-{src_im_suffix}') if '-MFS-' not in im])
        images_Q_all = sorted([im for im in glob.glob(f'{prefix}-*-Q-{src_im_suffix}') if '-MFS-' not in im])
        images_U_all = sorted([im for im in glob.glob(f'{prefix}-*-U-{src_im_suffix}') if '-MFS-' not in im])
        images_V_all = sorted([im for im in glob.glob(f'{prefix}-*-V-{src_im_suffix}') if '-MFS-' not in im])
        images_P_all = sorted([im for im in glob.glob(f'{prefix}-*-{pol_image_type}-{src_im_suffix}') if '-MFS-' not in im])

        cmap_I = build_channel_map(images_I_all)
        cmap_Q = build_channel_map(images_Q_all)
        cmap_U = build_channel_map(images_U_all)
        cmap_V = build_channel_map(images_V_all)
        cmap_P = build_channel_map(images_P_all)

        # Strict intersection: only channels present in ALL Stokes
        common = set(cmap_I) & set(cmap_Q) & set(cmap_U) & set(cmap_V) & set(cmap_P)
        valid_channels = sorted(common)

        n_I = len(cmap_I)
        if len(valid_channels) < n_I:
            dropped = sorted(set(cmap_I) - common)
            msg(f'WARNING: {n_I - len(valid_channels)}/{n_I} channels lack complete IQUV{pol_image_type} coverage, dropped: {dropped}')
        msg(f'Using {len(valid_channels)} channels with complete IQUV coverage')

    # Initialize 1D arrays for channelized data
    output_dictionary['CHAN'] = {}
    output_dictionary['CHAN']['freq_GHz'] = []
    output_dictionary['CHAN']['bmaj_asec'] = []
    output_dictionary['CHAN']['bmin_asec'] = []
    output_dictionary['CHAN']['bpa_deg'] = []
    output_dictionary['CHAN']['I_flux_mJy'] = []
    output_dictionary['CHAN']['I_err_mJy'] = []
    output_dictionary['CHAN']['I_rms_mJy'] = []
    output_dictionary['CHAN']['I_RA_deg'] = []
    output_dictionary['CHAN']['I_DEC_deg'] = []
    if not stokes_I_only:
        output_dictionary['CHAN']['Q_flux_mJy'] = []
        output_dictionary['CHAN']['U_flux_mJy'] = []
        output_dictionary['CHAN']['V_flux_mJy'] = []
        output_dictionary['CHAN'][f'{pol_image_type}_flux_mJy'] = []
        output_dictionary['CHAN']['Q_err_mJy'] = []
        output_dictionary['CHAN']['U_err_mJy'] = []
        output_dictionary['CHAN']['V_err_mJy'] = []
        output_dictionary['CHAN'][f'{pol_image_type}_err_mJy'] = []
        output_dictionary['CHAN']['Q_rms_mJy'] = []
        output_dictionary['CHAN']['U_rms_mJy'] = []
        output_dictionary['CHAN']['V_rms_mJy'] = []
        output_dictionary['CHAN'][f'{pol_image_type}_rms_mJy'] = []

    # Build one self-contained job per valid channel so it can run inside a
    # parallel worker (see fit_channel).
    jobs = []
    for ch_num in valid_channels:
        jobs.append({
            'ch_num': ch_num,
            'im_I': cmap_I[ch_num],
            'im_Q': cmap_Q[ch_num] if not stokes_I_only else None,
            'im_U': cmap_U[ch_num] if not stokes_I_only else None,
            'im_V': cmap_V[ch_num] if not stokes_I_only else None,
            'im_P': cmap_P[ch_num] if not stokes_I_only else None,
            'stokes_I_only': stokes_I_only,
            'pol_image_type': pol_image_type,
            'component': component,
            'MFS_I_ra_pix': MFS_I_ra_pix,
            'MFS_I_dec_pix': MFS_I_dec_pix,
            'MFS_bmaj': MFS_bmaj,
            'MFS_bmin': MFS_bmin,
            'MFS_freq_GHz': MFS_freq_GHz,
            'manual_rms_region': manual_rms_region,
            'force_fix_to_I': force_fix_to_I,
            'force_fix_V_to_I': force_fix_V_to_I,
            'chan_is_homogenized': chan_is_homogenized,
            'src_name': src_name,
        })

    # Size and run a one-shot pool for this field's channel fits. Unlike dev
    # (which loops over multiple epochs per script invocation and recycles
    # the pool periodically to bound worker memory growth), tkat runs one
    # field/epoch per script invocation, so no recycling is needed here.
    chan_results = {}
    if jobs:
        max_workers, _ = compute_max_workers([jobs[0]['im_I']], len(jobs))
        pool_ctx = mp.get_context('spawn')
        msg(f'CHAN-fitting pool starting with {max_workers} worker(s) for {len(jobs)} channel(s).')
        with ProcessPoolExecutor(max_workers=max_workers, mp_context=pool_ctx) as chan_pool:
            chunks = chunk_jobs(jobs, max_workers)
            futures = [chan_pool.submit(fit_channel_batch, chunk) for chunk in chunks]
            for future in as_completed(futures):
                for ch_num, result in future.result():
                    chan_results[ch_num] = result

    # Replay results in original channel order and append to the output
    # dictionary exactly as the old serial loop did.
    for ch_num in valid_channels:
        result = chan_results.get(ch_num)

        if result is None:
            msg(f'Fitting Failed or Channel Skipped: channel {ch_num} is likely flagged')
            continue

        output_dictionary['CHAN']['freq_GHz'].append(result['freq_GHz'])
        output_dictionary['CHAN']['bmaj_asec'].append(result['bmaj_asec'])
        output_dictionary['CHAN']['bmin_asec'].append(result['bmin_asec'])
        output_dictionary['CHAN']['bpa_deg'].append(result['bpa_deg'])
        output_dictionary['CHAN']['I_flux_mJy'].append(result['I_flux_mJy'])
        output_dictionary['CHAN']['I_err_mJy'].append(result['I_err_mJy'])
        output_dictionary['CHAN']['I_rms_mJy'].append(result['I_rms_mJy'])
        output_dictionary['CHAN']['I_RA_deg'].append(result['I_RA_deg'])
        output_dictionary['CHAN']['I_DEC_deg'].append(result['I_DEC_deg'])
        if not stokes_I_only:
            output_dictionary['CHAN']['Q_flux_mJy'].append(result['Q_flux_mJy'])
            output_dictionary['CHAN']['U_flux_mJy'].append(result['U_flux_mJy'])
            output_dictionary['CHAN']['V_flux_mJy'].append(result['V_flux_mJy'])
            output_dictionary['CHAN'][f'{pol_image_type}_flux_mJy'].append(result[f'{pol_image_type}_flux_mJy'])
            output_dictionary['CHAN']['Q_err_mJy'].append(result['Q_err_mJy'])
            output_dictionary['CHAN']['U_err_mJy'].append(result['U_err_mJy'])
            output_dictionary['CHAN']['V_err_mJy'].append(result['V_err_mJy'])
            output_dictionary['CHAN'][f'{pol_image_type}_err_mJy'].append(result[f'{pol_image_type}_err_mJy'])
            output_dictionary['CHAN']['Q_rms_mJy'].append(result['Q_rms_mJy'])
            output_dictionary['CHAN']['U_rms_mJy'].append(result['U_rms_mJy'])
            output_dictionary['CHAN']['V_rms_mJy'].append(result['V_rms_mJy'])
            output_dictionary['CHAN'][f'{pol_image_type}_rms_mJy'].append(result[f'{pol_image_type}_rms_mJy'])


    # Create timestamp prefix from time_info if available
    timestamp_prefix = ""
    if time_info and 'start_mjd' in time_info:
        # Convert MJD to datetime and format as YYYYMMDDTHH
        # MJD epoch is November 17, 1858 at 00:00:00 UTC
        mjd_epoch = datetime.datetime(1858, 11, 17, 0, 0, 0)
        dt = mjd_epoch + datetime.timedelta(days=time_info['start_mjd'])
        timestamp_prefix = dt.strftime('%Y%m%dT%H%M_')
        msg(f'Using timestamp prefix: {timestamp_prefix}')

    # Store epoch timing (MJD/ISOT centre, duration) in the MFS dict from the
    # already-loaded time_info (read_time_info()) -- no new loading/lookup.
    if time_info and 'middle_mjd' in time_info:
        output_dictionary['MFS']['time_ctr_mjd']  = time_info['middle_mjd']
        output_dictionary['MFS']['time_ctr_isot'] = _mjd_to_isot(time_info['middle_mjd'])
        output_dictionary['MFS']['time_dt']       = time_info['duration_hours'] * 3600.0

    # Spectral index fit (chi2-aware) on the MFS Stokes I reference, stored
    # here so it is captured in the JSON below.
    mfs = output_dictionary['MFS']
    if 'I_flux_mJy' in mfs and 'I_rms_mJy' in mfs and 'freq_GHz' in mfs:
        spec_fit = _fit_alpha_chi2(
            output_dictionary['CHAN']['freq_GHz'],
            output_dictionary['CHAN']['I_flux_mJy'],
            output_dictionary['CHAN']['I_rms_mJy'],
            mfs['I_flux_mJy'], mfs['I_rms_mJy'], mfs['freq_GHz'])
        mfs['alpha']     = spec_fit['alpha']
        mfs['alpha_err'] = spec_fit['alpha_err']
        mfs['chi2']      = spec_fit['chi2']
        mfs['ndof']      = spec_fit['ndof']
        mfs['chi2_red']  = spec_fit['chi2_red']

    # Save output with timestamp prefix
    output_file = cfg.RESULTS + f'/{timestamp_prefix}{src_name}_{IDENTIFIER}_polarization.json'
    with open(output_file, 'w') as j:
        json.dump(output_dictionary, j, indent=4)
    
    msg(f'Saved output to: {output_file}')
    
    # Also save text file with timestamp prefix
    if stokes_I_only:
        txt_file = cfg.RESULTS + f'/{timestamp_prefix}{src_name}_{IDENTIFIER}_stokesI.txt'
    else:
        txt_file = cfg.RESULTS + f'/{timestamp_prefix}{src_name}_{IDENTIFIER}_iquv.txt'
    with open(txt_file, 'w') as f:
        f.write(f"# Source: {src_name}\n")
        if 'I_RA_deg' in output_dictionary.get('MFS', {}):
            f.write(f"# Position: RA={output_dictionary['MFS']['I_RA_deg']:.6f} deg, "
                    f"Dec={output_dictionary['MFS']['I_DEC_deg']:.6f} deg\n")
        if time_info:
            f.write(f"# MS: {time_info.get('ms_name', 'N/A')}\n")
            f.write(f"# Scan: {time_info.get('scan', 'N/A')}\n")
            f.write(f"# Start MJD: {time_info.get('start_mjd', 0):.10f}\n")
            f.write(f"# End MJD: {time_info.get('end_mjd', 0):.10f}\n")
            f.write(f"# Middle MJD: {time_info.get('middle_mjd', 0):.10f}\n")
            f.write(f"# Duration: {time_info.get('duration_hours', 0):.4f} hours\n")
        if stokes_I_only:
            f.write(f"# Columns: Channel Freq[GHz] I[mJy] rms_I\n")
        else:
            f.write(f"# Columns: Channel Freq[GHz] I[mJy] Q[mJy] U[mJy] V[mJy] "
                    f"{pol_image_type}[mJy] rms_I rms_Q rms_U rms_V rms_{pol_image_type}\n")
        f.write("#\n")
        for i in range(len(output_dictionary['CHAN']['freq_GHz'])):
            if stokes_I_only:
                f.write(f"{i:4d} "
                       f"{output_dictionary['CHAN']['freq_GHz'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['I_flux_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['I_rms_mJy'][i]:10.4f}\n")
            else:
                f.write(f"{i:4d} "
                       f"{output_dictionary['CHAN']['freq_GHz'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['I_flux_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['Q_flux_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['U_flux_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['V_flux_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN'][f'{pol_image_type}_flux_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['I_rms_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['Q_rms_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['U_rms_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN']['V_rms_mJy'][i]:10.4f} "
                       f"{output_dictionary['CHAN'][f'{pol_image_type}_rms_mJy'][i]:10.4f}\n")
    msg(f'Saved text file to: {txt_file}')
    
    # Create additional raw data file for polang if applicable
    # This file contains: Freq(Hz), I, Q, U, dI, dQ, dU in Jy with no header
    return output_dictionary, timestamp_prefix


def create_polang_raw_file(project_info, src_name, output_dict, timestamp_prefix):
    """
    Create raw data file for polang analysis if target matches polang_name.
    Format: Freq(Hz) I Q U dI dQ dU (all fluxes in Jy, no header)
    """
    # Check if polang_name exists and is not empty
    if 'polang_name' not in project_info or not project_info['polang_name']:
        return
    
    # Polang requires Q/U
    if 'Q_flux_mJy' not in output_dict.get('CHAN', {}):
        msg('Skipping polang raw data file (Stokes-I-only mode)')
        return
    
    polang_name = project_info['polang_name']
    
    # Check if target includes polang_name
    if polang_name.lower() not in src_name.lower():
        msg(f'Source {src_name} does not include polang_name {polang_name}, skipping raw data file')
        return
    
    msg(f'Creating raw data file for polang analysis (polang_name={polang_name})')
    
    # Create raw data file with timestamp prefix
    raw_file = cfg.RESULTS + f'/{timestamp_prefix}{src_name}_{IDENTIFIER}_rmsynth.txt'
    
    with open(raw_file, 'w') as f:
        # No header - just data
        for i in range(len(output_dict['CHAN']['freq_GHz'])):
            freq_hz = output_dict['CHAN']['freq_GHz'][i] * 1e9  # Convert GHz to Hz
            I_jy = output_dict['CHAN']['I_flux_mJy'][i] / 1e3  # Convert mJy to Jy
            Q_jy = output_dict['CHAN']['Q_flux_mJy'][i] / 1e3
            U_jy = output_dict['CHAN']['U_flux_mJy'][i] / 1e3
            dI_jy = output_dict['CHAN']['I_rms_mJy'][i] / 1e3
            dQ_jy = output_dict['CHAN']['Q_rms_mJy'][i] / 1e3
            dU_jy = output_dict['CHAN']['U_rms_mJy'][i] / 1e3
            
            f.write(f"{freq_hz:.6e} {I_jy:.6e} {Q_jy:.6e} {U_jy:.6e} "
                   f"{dI_jy:.6e} {dQ_jy:.6e} {dU_jy:.6e}\n")
    
    msg(f'Saved polang raw data file to: {raw_file}')


def _mad_ylim(arr, err=None, k=10.0, pad=0.10):
    """
    Robust MAD-based y-axis limits.

    1. Exclude flux points outside median +/- k*sigma_flux.
    2. Also exclude points whose error exceeds median(err) + k*sigma_err.
    3. Set limits from min/max of survivors, then pad by `pad` fraction.
    Data drives the range entirely — no zero-floor is imposed.
    """
    arr = np.asarray(arr, dtype=float)
    err = np.zeros_like(arr) if err is None else np.asarray(err, dtype=float)

    med_f   = np.median(arr)
    mad_f   = np.median(np.abs(arr - med_f))
    sigma_f = mad_f * 1.4826
    if sigma_f == 0:
        sigma_f = np.std(arr) if len(arr) > 1 else 1.0
    flux_mask = (arr >= med_f - k * sigma_f) & (arr <= med_f + k * sigma_f)

    med_e   = np.median(err)
    mad_e   = np.median(np.abs(err - med_e))
    sigma_e = mad_e * 1.4826
    if sigma_e == 0:
        sigma_e = np.std(err) if len(err) > 1 else 1.0
    err_mask = err <= med_e + k * sigma_e

    mask = flux_mask & err_mask
    if not mask.any():
        mask = np.ones(len(arr), dtype=bool)

    lo = np.min(arr[mask])
    hi = np.max(arr[mask])
    if lo == hi:
        lo -= 1.0
        hi += 1.0
    margin = pad * (hi - lo)
    return (lo - margin, hi + margin)


def _fit_alpha_chi2(freq, flux, rms, mfs_flux, mfs_rms, mfs_freq):
    """
    Fit a power-law spectral index alpha to Stokes I channel data, gated on
    MFS S/N (SPEC_INDEX_SNR_THRESH) with iterative MAD-based outlier rejection
    on flux (SPEC_INDEX_MAD_CLIP). Fitting is done in log-log space
    (S = S_ref * (nu/nu_ref)^alpha) with nu_ref = mfs_freq.

    Returns
    -------
    dict with keys: alpha, alpha_err, chi2, ndof, chi2_red, fit_nu, fit_S
    All values are None if the fit could not be performed.
    """
    result = {'alpha': None, 'alpha_err': None, 'chi2': None, 'ndof': None,
              'chi2_red': None, 'fit_nu': None, 'fit_S': None}

    freq = np.asarray(freq, dtype=float)
    flux = np.asarray(flux, dtype=float)
    rms  = np.asarray(rms,  dtype=float)

    mfs_snr = mfs_flux / mfs_rms if mfs_rms else 0.0

    # Quality mask: positive flux and per-channel S/N >= 3
    valid = (flux > 0) & (rms > 0) & (flux / rms >= 3.0)

    if mfs_snr <= SPEC_INDEX_SNR_THRESH:
        msg(f'  Spectral index fit skipped: MFS S/N = {mfs_snr:.1f} (< {SPEC_INDEX_SNR_THRESH})')
        return result
    if np.sum(valid) < 3:
        msg(f'  Spectral index fit skipped: only {int(np.sum(valid))} valid channels (need >= 3)')
        return result

    nu_v = freq[valid]
    S_v  = flux[valid]
    w_v  = S_v / rms[valid]   # S/N weights (= 1/sigma_logS in log space)

    # Iterative MAD-based outlier rejection on flux values
    _iter = 0
    while len(nu_v) >= 3:
        _med = np.median(S_v)
        _mad = np.median(np.abs(S_v - _med))
        _sig = 1.4826 * _mad
        if _sig <= 0:
            break
        _keep     = np.abs(S_v - _med) <= SPEC_INDEX_MAD_CLIP * _sig
        n_clipped = int(np.sum(~_keep))
        msg(f'  MAD clip iter {_iter + 1}: median = {_med:.3f} mJy  sigma = {_sig:.3f} mJy  '
            f'threshold = {SPEC_INDEX_MAD_CLIP * _sig:.3f} mJy  clipped = {n_clipped}  '
            f'remaining = {int(np.sum(_keep))}')
        if n_clipped == 0:
            break
        nu_v  = nu_v[_keep]
        S_v   = S_v[_keep]
        w_v   = w_v[_keep]
        _iter += 1

    if len(nu_v) < 3:
        msg(f'  Spectral index fit skipped: only {len(nu_v)} channels remain after MAD clipping (need >= 3)')
        return result

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

        msg(f'  Spectral index fit: alpha = {alpha:.3f} +/- {alpha_err:.3f}  '
            f'chi2 = {chi2:.2f}  chi2_red = {chi2_red if chi2_red is None else round(chi2_red, 2)}  ndof = {ndof}  '
            f'(MFS S/N = {mfs_snr:.1f} >= {SPEC_INDEX_SNR_THRESH}, N_chan = {len(nu_v)} / {int(np.sum(valid))})')

        result.update(alpha=alpha, alpha_err=alpha_err, chi2=chi2, ndof=ndof,
                       chi2_red=chi2_red, fit_nu=nu_plot, fit_S=fit_S)
    except Exception as e:
        msg(f'  WARNING: Spectral index fit failed: {e}')

    return result


def plot_stokes_spectrum(output_dictionary, src_name, timestamp_prefix):
    """
    Plot Stokes I(QUV) spectra and save to RESULTS/fitting_plots.

    Layout  : 4x1 grid (IQUV) or 1x1 (I only), all panels share x-axis.
    Y-ranges: MAD-based, data-driven, no zero-floor on any panel.
    y=0 line: drawn on Q, U, V panels only.
    Alpha   : MFS S/N-gated, MAD-clipped power-law fit to Stokes I channels
              (see _fit_alpha_chi2); overplotted with alpha/chi2_red annotation.
    """
    chan = output_dictionary['CHAN']
    if len(chan['freq_GHz']) == 0:
        return

    # Output directory
    plot_dir = os.path.join(cfg.RESULTS, 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)

    freq     = np.asarray(chan['freq_GHz'], dtype=float)
    has_pol  = 'Q_flux_mJy' in chan
    mfs      = output_dictionary.get('MFS', {})

    stokes = [('I', chan['I_flux_mJy'], chan['I_rms_mJy'], 'tab:blue')]
    if has_pol:
        stokes += [
            ('Q', chan['Q_flux_mJy'], chan['Q_rms_mJy'], 'tab:orange'),
            ('U', chan['U_flux_mJy'], chan['U_rms_mJy'], 'tab:green'),
            ('V', chan['V_flux_mJy'], chan['V_rms_mJy'], 'tab:red'),
        ]

    n_panels  = len(stokes)
    fig, axes = plt.subplots(n_panels, 1, figsize=(12, 3.5 * n_panels),
                             sharex=True, squeeze=False)
    axes      = axes[:, 0]

    for ax, (label, flux, rms, colour) in zip(axes, stokes):
        flux_arr = np.asarray(flux, dtype=float)
        rms_arr  = np.asarray(rms,  dtype=float)

        ax.errorbar(freq, flux_arr, yerr=rms_arr,
                    fmt='o', markersize=3, capsize=2, color=colour,
                    elinewidth=0.8, linewidth=0, label=f'Stokes {label}')

        # MFS horizontal reference line
        mfs_key = f'{label}_flux_mJy'
        if mfs_key in mfs:
            ax.axhline(y=mfs[mfs_key], color='k', linestyle='--', linewidth=1.5,
                       label=f'MFS: {mfs[mfs_key]:.2f} mJy')

        # y=0 line on Q, U, V only
        if label in ('Q', 'U', 'V'):
            ax.axhline(0, color='grey', linestyle=':', linewidth=0.8)

        # Spectral index fit on Stokes I (same chi2-aware fit stored in the JSON;
        # recomputed here purely for the plot curve, so the two never disagree)
        if label == 'I' and 'I_flux_mJy' in mfs and 'I_rms_mJy' in mfs and 'freq_GHz' in mfs:
            fit = _fit_alpha_chi2(freq, flux_arr, rms_arr,
                                   mfs['I_flux_mJy'], mfs['I_rms_mJy'], mfs['freq_GHz'])
            if fit['alpha'] is not None:
                label_str = rf"$\alpha = {fit['alpha']:.2f} \pm {fit['alpha_err']:.2f}$"
                if fit['chi2_red'] is not None:
                    label_str += rf"  ($\chi^2_\nu = {fit['chi2_red']:.2f}$)"
                ax.plot(fit['fit_nu'], fit['fit_S'], color='tomato', linewidth=1.5,
                        linestyle='-', label=label_str)

        # MAD-based y-limits, data-driven
        ax.set_ylim(_mad_ylim(flux_arr, err=rms_arr))

        ax.set_ylabel(f'S (mJy) [Stokes {label}]', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.tick_params(direction='in', which='both', top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

    axes[-1].set_xlabel('Frequency (GHz)', fontsize=11)
    suffix = 'IQUV' if has_pol else 'I'
    axes[0].set_title(f'{src_name}  —  Stokes {suffix} Spectra', fontsize=12)

    plt.tight_layout()

    plot_file = os.path.join(
        plot_dir,
        f'{timestamp_prefix}{src_name}_{suffix}_spectrum.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'Saved {suffix} spectrum plot to: {plot_file}')


def main():
    
    # Parse command line arguments
    if len(sys.argv) < 2:
        msg("ERROR: Usage: RMSYNTH_01_extract_fluxes.py <source_name>")
        msg("  source_name: Name of source to process")
        msg("    - For targets: source_name")
        msg("    - For per-scan calibrators: cal_name_scanN (e.g., J1939-6342_scan3)")
        sys.exit(1)
    
    source_name = sys.argv[1]
    
    # Parse source_name to detect per-scan calibrators
    scan_match = re.search(r'(.+)_scan(\d+)$', source_name)
    if scan_match:
        # Per-scan calibrator
        field_name = scan_match.group(1)
        scan_num = scan_match.group(2)
        msg(f"Per-scan calibrator detected: field={field_name}, scan={scan_num}")
    else:
        # Regular target or calibrator
        field_name = source_name
        scan_num = None
        msg(f"Regular source: {field_name}")
    
    # Construct image directory from source name
    image_dir = f'IMAGES/{source_name}'
    msg(f"Starting IQUV extraction for source: {source_name}")
    msg(f"Image directory: {image_dir}")
    
    # Get MS path from project info
    with open('project_info.json') as f:
        project_info = json.load(f)
    myms = project_info['working_ms']
    
    # Query MS for source information using field_name (not full source_name)
    source_info = get_source_from_ms(myms, field_name)
    
    # Manually use incase you accidentally deleted 1024ch like my dumbass:
    # source_info = {
    #     'name': 'SwiftJ1727',
    #     'ra_deg': 0.0,
    #     'dec_deg': 0.0,
    #     'intent': 'TARGET',
    #     'found': True
    # }
    
    if not source_info['found']:
        msg("ERROR: Cannot proceed without source information")
        sys.exit(1)
    
    # Determine position based on intent
    is_target = 'TARGET' in source_info['intent'].upper()

    # Determine if we should force all Stokes positions to the Stokes I position.
    # Applies to non-targets not excluded by EXCLUDED_CALIBRATOR_SUBSTRINGS.
    force_fix_to_I = False
    if USE_CALIBRATOR_POSITION_FIXING and not is_target:
        contains_excluded = any(substring in field_name for substring in EXCLUDED_CALIBRATOR_SUBSTRINGS)
        if not contains_excluded:
            force_fix_to_I = True
            msg('')
            msg(f"{'='*80}")
            msg(f"NON-TARGET CALIBRATOR DETECTED: {field_name}")
            msg(f"All positions will be FIXED to Stokes I position (no free fitting)")
            msg(f"This behavior is disabled for sources containing: {', '.join(EXCLUDED_CALIBRATOR_SUBSTRINGS)}")
            msg(f"{'='*80}")
            msg('')
        else:
            msg('')
            msg(f"Non-target calibrator '{field_name}' contains excluded substring")
            msg(f"Using standard S/N-based position fitting logic")
            msg('')
    elif not USE_CALIBRATOR_POSITION_FIXING:
        msg('')
        msg(f"USE_CALIBRATOR_POSITION_FIXING is disabled - using standard S/N-based fitting for all sources")
        msg('')

    # Determine if Stokes V should be forced to the Stokes I position independently
    # of force_fix_to_I.  V is forced unless the field name contains a substring
    # from V_FORCE_EXCLUDED_SUBSTRINGS.
    force_fix_V_to_I = not any(substring in field_name for substring in V_FORCE_EXCLUDED_SUBSTRINGS)
    if force_fix_V_to_I and not force_fix_to_I:
        msg(f"Stokes V position will be forced to Stokes I reference for: {field_name}")
    
    if is_target:
        msg("Source is a TARGET - reading position from file")
        position_file = cfg.DATA + '/positions/XRB_pos_list.txt'
        ra_deg, dec_deg = read_position_file(position_file, field_name)
        if ra_deg is None:
            msg("WARNING: Position not found in file, using MS position")
            ra_deg = source_info['ra_deg']
            dec_deg = source_info['dec_deg']
        else:
            # Calculate and report position difference
            ra_ms = source_info['ra_deg']
            dec_ms = source_info['dec_deg']
            delta_ra_arcsec = (ra_deg - ra_ms) * 3600.0
            delta_dec_arcsec = (dec_deg - dec_ms) * 3600.0
            msg(f"Position difference (file - MS):")
            msg(f"  Delta RA:  {delta_ra_arcsec:+.2f} arcsec")
            msg(f"  Delta Dec: {delta_dec_arcsec:+.2f} arcsec")
            source_info['ra_deg'] = ra_deg
            source_info['dec_deg'] = dec_deg
    else:
        msg("Source is NOT a target - using MS position")
        ra_deg = source_info['ra_deg']
        dec_deg = source_info['dec_deg']
    
    # Read time information
    # Match by field name and scan number (if per-scan)
    # Time file is named: {working_ms}_time_info.txt in RESULTS directory
    time_file = cfg.RESULTS + '/' + myms + '_time_info.txt'
    time_info = read_time_info(time_file, field_name, scan=scan_num)
    
    # Construct image identifier - check for non-zoom first, then fall back to zoom
    src_im_identifier_nozoom = f'{image_dir}/*{field_name}*.ms_{IDENTIFIER}-'
    src_im_identifier_zoom = f'{image_dir}/*{field_name}*.ms_{IDENTIFIER}_zoom-'
    
    # Check if non-zoom images exist — match both plain channel numbers and time-interval prefixes
    test_images = glob.glob(f'{src_im_identifier_nozoom}[0t]*-I-image.fits')
    if test_images:
        src_im_identifier = src_im_identifier_nozoom
        msg(f'Using non-zoom images: {src_im_identifier}')
    else:
        src_im_identifier = src_im_identifier_zoom
        msg(f'Non-zoom images not found, using zoom images: {src_im_identifier}')

    src_im_suffix = 'image.fits'  # Will be checked and adjusted inside the function

    src_ulims = [False]  # Single component, not an upper limit

    # Sources with known position/confusion issues get a manual RMS region
    # from data/positions/RMS_region_list.txt instead of the default annulus
    # (mirrors read_position_file's target-position override above).
    rms_region_file = cfg.DATA + '/positions/RMS_region_list.txt'
    manual_rms_region = read_rms_region_file(rms_region_file, field_name)
    if not manual_rms_region:
        manual_rms_region = False

    # Log the RMS region choice once here before extraction begins
    if manual_rms_region:
        msg(f'RMS region: MANUAL -- {manual_rms_region}')
    else:
        msg('RMS region: default annulus (~500 beam areas) centred on source')

    # Use Plin instead of Ptot when a polarization angle calibrator is present,
    # as Plin (sqrt(Q^2+U^2)) is more suitable for angle calibration than Ptot (sqrt(Q^2+U^2+V^2)).
    use_plin = bool(project_info.get('polang_name', ''))
    if use_plin:
        msg(f'polang_name is set ({project_info["polang_name"]}): will use Plin images for polarized intensity')
    else:
        msg('polang_name is not set: will use Ptot images for polarized intensity')

    # Detect time-interval images (e.g. ..._pcalmask-t0000-0501-I-image.fits)
    msg(f'Checking for time-interval images: {src_im_identifier}t*-I-image.fits')
    time_interval_images = glob.glob(f'{src_im_identifier}t*-I-image.fits')
    if time_interval_images:
        time_suffixes = sorted(set(
            re.search(r'(t\d+)', o.basename(im)).group(1)
            for im in time_interval_images
            if re.search(r'(t\d+)', o.basename(im))
        ))
        msg(f'Time-interval mode: found {len(time_interval_images)} images across {len(time_suffixes)} intervals: {time_suffixes}')
        for tsuf in time_suffixes:
            timed_identifier = f'{src_im_identifier}{tsuf}-'
            timed_name = f'{source_name}_{tsuf}'
            msg(f'Starting property extraction for time interval {tsuf}: {timed_identifier}')
            output_dict, timestamp_prefix = extract_polarization_properties(
                timed_name,
                timed_identifier,
                src_im_suffix,
                ra_deg,
                dec_deg,
                src_ulims,
                manual_rms_region,
                image_dir,
                '',
                time_info = time_info,
                force_fix_to_I = force_fix_to_I,
                force_fix_V_to_I = force_fix_V_to_I,
                use_plin = use_plin)
            create_polang_raw_file(project_info, timed_name, output_dict, timestamp_prefix)
            plot_stokes_spectrum(output_dict, timed_name, timestamp_prefix)
            # Trailing '*' before .txt also catches the per-channel estimate/check_pos
            # files from fit_channel, which append _{ch_num} after the name.
            for temp_file in glob.glob(f'check_pos_*_{timed_name}*.txt'):
                os.remove(temp_file)
            for temp_file in glob.glob(f'estimate_*_{timed_name}*.txt'):
                os.remove(temp_file)
    else:
        msg(f'No time-interval images found — using time-averaged (single-call) mode')
        msg(f'Starting property extraction for identifier: {src_im_identifier}')
        output_dict, timestamp_prefix = extract_polarization_properties(
            source_name,  # Use full source_name (includes _scanN for per-scan)
            src_im_identifier,
            src_im_suffix,
            ra_deg,  # Single value, not array
            dec_deg,  # Single value, not array
            src_ulims,
            manual_rms_region,
            image_dir,
            '',  # image_identifier not needed
            time_info = time_info,
            force_fix_to_I = force_fix_to_I,
            force_fix_V_to_I = force_fix_V_to_I,
            use_plin = use_plin)
        create_polang_raw_file(project_info, source_name, output_dict, timestamp_prefix)
        plot_stokes_spectrum(output_dict, source_name, timestamp_prefix)
        # Trailing '*' before .txt also catches the per-channel estimate/check_pos
        # files from fit_channel, which append _{ch_num} after the name.
        for temp_file in glob.glob(f'check_pos_*_{source_name}*.txt'):
            os.remove(temp_file)
        for temp_file in glob.glob(f'estimate_*_{source_name}*.txt'):
            os.remove(temp_file)



if __name__  == "__main__":
    main()
