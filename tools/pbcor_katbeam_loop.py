#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# Primary beam correction for MeerKAT channel/time image cubes using katbeam.
#
# Loops over all images matching a glob pattern, evaluates the katbeam model
# at each image's native frequency, and writes PB-corrected FITS files.
# Parallelised with ProcessPoolExecutor, SLURM-aware.
#
# Requires:
#   https://github.com/ludwigschwardt/katbeam
#   https://pypi.org/project/scikit-ued/  (only if AZAVG = True)


import glob
import numpy as np
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from astropy.io import fits
from katbeam import JimBeam
from shutil import copyfile


# =============================================================================
# Configuration -- edit these before running
# =============================================================================

# Glob pattern for input images.  Should match all channel/time FITS files
# to be corrected, e.g. the WSClean per-channel outputs.
IMAGE_GLOB = 'INTERVALS/*restored*image.fits'

# MeerKAT band: 'L', 'UHF', or 'S'
BAND = 'L'

# FITS header axis that contains the frequency (3 for WSClean, 4 for DDFacet)
FREQ_AXIS = '3'

# Primary beam gain level below which to blank pixels
PB_CUT = 0.2

# Azimuthally average the beam pattern (recommended False for speed; the
# katbeam analytic model is already very close to circularly symmetric)
AZAVG = False

# Output products: toggle which files are written per input image
SAVE_PBCOR = True        # input / beam  (PB-corrected image)
SAVE_PB    = False        # beam image
SAVE_WT    = False        # beam^2 weight image

# Overwrite existing output files
OVERWRITE = True

# =============================================================================
# End configuration
# =============================================================================


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


# ---- SLURM-aware worker calculation (from restore_model.py) ----------------

def get_available_memory_bytes():
    # --- Check SLURM memory allocation first ---
    slurm_mem = None
    slurm_src = None

    mem_per_node = os.environ.get('SLURM_MEM_PER_NODE')
    mem_per_cpu  = os.environ.get('SLURM_MEM_PER_CPU')

    if mem_per_node:
        # SLURM_MEM_PER_NODE is in MB (or with suffix like '64G')
        val = mem_per_node.rstrip('MmGgTt')
        multiplier = 1024 ** 3 if mem_per_node[-1].lower() in ('g', 't') else 1024 ** 2
        if mem_per_node[-1].lower() == 't':
            multiplier = 1024 ** 4
        slurm_mem = int(val) * multiplier
        slurm_src = f'SLURM_MEM_PER_NODE={mem_per_node}'

    elif mem_per_cpu:
        val = mem_per_cpu.rstrip('MmGgTt')
        multiplier = 1024 ** 3 if mem_per_cpu[-1].lower() in ('g', 't') else 1024 ** 2
        if mem_per_cpu[-1].lower() == 't':
            multiplier = 1024 ** 4
        per_cpu_bytes = int(val) * multiplier

        # Multiply by allocated CPUs
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
        slurm_mem = per_cpu_bytes * n_cpus
        slurm_src = f'SLURM_MEM_PER_CPU={mem_per_cpu} x {n_cpus} CPUs'

    # --- System-reported memory (total and available) ---
    sys_total = None
    sys_avail = None
    try:
        import psutil
        vm = psutil.virtual_memory()
        sys_total = vm.total
        sys_avail = vm.available
    except ImportError:
        pass
    if sys_total is None or sys_avail is None:
        try:
            with open('/proc/meminfo') as f:
                for line in f:
                    if line.startswith('MemTotal:') and sys_total is None:
                        sys_total = int(line.split()[1]) * 1024
                    if line.startswith('MemAvailable:') and sys_avail is None:
                        sys_avail = int(line.split()[1]) * 1024
        except OSError:
            pass

    # --- Log and return the more restrictive of the two ---
    if sys_total is not None:
        msg(f'Memory (total)     : {sys_total / 1024**3:.2f} GB')
    if sys_avail is not None:
        msg(f'Memory (available) : {sys_avail / 1024**3:.2f} GB')
    if slurm_mem is not None:
        msg(f'Memory (SLURM)     : {slurm_mem / 1024**3:.2f} GB  [{slurm_src}]')

    if slurm_mem is not None and sys_avail is not None:
        effective = min(slurm_mem, sys_avail)
        msg(f'Memory (effective) : {effective / 1024**3:.2f} GB  [min(SLURM, available)]')
        return effective
    elif slurm_mem is not None:
        return slurm_mem
    elif sys_avail is not None:
        return sys_avail
    else:
        msg('WARNING: Could not determine available memory. Assuming 4 GB.')
        return 4 * 1024 ** 3


def estimate_memory_per_worker_bytes(fitsfile):
    """Conservative peak estimate: 4 array copies + 100 MB process overhead."""
    with fits.open(fitsfile) as hdul:
        return 4 * hdul[0].data.nbytes + 100 * 1024 ** 2


def compute_max_workers(image_list):
    available  = get_available_memory_bytes()
    per_worker = estimate_memory_per_worker_bytes(image_list[0])
    mem_cap    = max(1, int(available // per_worker))

    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _job_cpus      = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    _job_cpus_m    = re.match(r'(\d+)', _job_cpus)

    if _cpus_per_task and _ntasks:
        cpu_cap  = int(_cpus_per_task) * int(_ntasks)
        _cpu_src = f'SLURM_CPUS_PER_TASK={_cpus_per_task} x SLURM_NTASKS={_ntasks}'
    elif _cpus_per_task:
        cpu_cap  = int(_cpus_per_task)
        _cpu_src = f'SLURM_CPUS_PER_TASK={_cpus_per_task}'
    elif _job_cpus_m:
        cpu_cap  = int(_job_cpus_m.group(1))
        _cpu_src = f'SLURM_JOB_CPUS_PER_NODE={_job_cpus} -> {cpu_cap}'
    elif _ntasks:
        cpu_cap  = int(_ntasks)
        _cpu_src = f'SLURM_NTASKS={_ntasks}'
    else:
        cpu_cap  = os.cpu_count() or 1
        _cpu_src = 'os.cpu_count() (no SLURM allocation detected)'

    workers = min(mem_cap, cpu_cap, len(image_list))
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={len(image_list)})]')
    return workers


# ---- FITS I/O helpers -------------------------------------------------------

def get_header(fitsfile, freqaxis):
    inphdu = fits.open(fitsfile)
    inphdr = inphdu[0].header
    nx   = inphdr.get('NAXIS1')
    ny   = inphdr.get('NAXIS2')
    dx   = inphdr.get('CDELT1')
    dy   = inphdr.get('CDELT2')
    freq = inphdr.get('CRVAL' + freqaxis)
    inphdu.close()
    return nx, ny, dx, dy, freq


def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
        image = np.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        image = np.array(input_hdu.data[0, :, :])
    else:
        image = np.array(input_hdu.data[0, 0, :, :])
    return image


def flush_fits(newimage, fitsfile):
    f = fits.open(fitsfile, mode='update')
    input_hdu = f[0]
    if len(input_hdu.data.shape) == 2:
        input_hdu.data[:, :] = newimage
    elif len(input_hdu.data.shape) == 3:
        input_hdu.data[0, :, :] = newimage
    else:
        input_hdu.data[0, 0, :, :] = newimage
    f.flush()
    f.close()


# ---- Beam evaluation --------------------------------------------------------

def azimuthal_average_fast(image, cx, cy):
    """
    Vectorised azimuthal average using np.bincount.  Replaces the O(N^2)
    Python double-loop in the original pbcor_katbeam.py with a single-pass
    numpy operation.  ~100x faster for typical 4k x 4k images.
    """
    ny, nx = image.shape
    y, x = np.ogrid[:ny, :nx]
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2).astype(int)

    # Mask NaN pixels so they don't contaminate the average
    valid = np.isfinite(image)
    r_flat = r[valid].ravel()
    v_flat = image[valid].ravel()

    sum_per_bin   = np.bincount(r_flat, weights=v_flat)
    count_per_bin = np.bincount(r_flat)

    # Avoid division by zero for radial bins with no valid pixels
    avg = np.where(count_per_bin > 0, sum_per_bin / count_per_bin, np.nan)

    # Map back: replace each pixel with the average at its radius
    result = avg[r]
    return result


def evaluate_beam(beam, nx, ny, dx, freq_mhz, pbcut, azavg):
    """
    Evaluate the katbeam model at a given frequency and return the beam image.

    Parameters
    ----------
    beam     : JimBeam object
    nx, ny   : int, image dimensions (must be square)
    dx       : float, pixel scale in degrees (CDELT1, typically negative)
    freq_mhz : float, frequency in MHz
    pbcut    : float, gain level below which to blank
    azavg    : bool, apply azimuthal averaging

    Returns
    -------
    beam_image : 2D numpy array (ny, nx), NaN below pbcut
    """
    extent = nx * dx  # degrees (signed)
    interval = np.linspace(-extent / 2.0, extent / 2.0, nx)
    xx, yy = np.meshgrid(interval, interval)
    beam_image = beam.I(xx, yy, freq_mhz)

    # Blank below the PB cut level
    beam_image[beam_image < pbcut] = np.nan

    if azavg:
        cx = int(nx / 2)
        cy = int(ny / 2)
        beam_image = azimuthal_average_fast(beam_image, cx, cy)

    return beam_image


# ---- Per-image worker --------------------------------------------------------

def pbcor_one(input_fits, beam_model_name, freqaxis, pbcut, azavg,
              save_pbcor, save_pb, save_wt, overwrite, idx, n_total):
    """
    Apply primary beam correction to a single FITS image.
    Designed to be called from ProcessPoolExecutor.

    Returns a status string.
    """
    try:
        msg(f'[{idx}/{n_total}] Processing: {os.path.basename(input_fits)}')

        # Read header
        nx, ny, dx, dy, fitsfreq = get_header(input_fits, freqaxis)
        if nx != ny or abs(dx) != abs(dy):
            return f'SKIPPED (non-square image/pixels): {input_fits}'

        freq_mhz = fitsfreq / 1e6

        # Generate output filenames
        pbcor_fits = input_fits.replace('.fits', '.pbcor.fits')
        pb_fits    = input_fits.replace('.fits', '.pb.fits')
        wt_fits    = input_fits.replace('.fits', '.wt.fits')

        # Skip if outputs exist and overwrite is off
        if not overwrite:
            skip = False
            if save_pbcor and os.path.isfile(pbcor_fits):
                skip = True
            if save_pb and os.path.isfile(pb_fits):
                skip = True
            if save_wt and os.path.isfile(wt_fits):
                skip = True
            if skip:
                return f'SKIPPED (outputs exist): {input_fits}'

        # Instantiate beam model (per-worker to avoid pickling issues)
        beam = JimBeam(beam_model_name)

        # Evaluate beam at this channel's frequency
        beam_image = evaluate_beam(beam, nx, ny, dx, freq_mhz, pbcut, azavg)

        # Write requested output products
        if save_pbcor:
            input_image = get_image(input_fits)
            pbcor_image = input_image / beam_image
            copyfile(input_fits, pbcor_fits)
            flush_fits(pbcor_image, pbcor_fits)

        if save_pb:
            copyfile(input_fits, pb_fits)
            flush_fits(beam_image, pb_fits)

        if save_wt:
            copyfile(input_fits, wt_fits)
            flush_fits(beam_image ** 2.0, wt_fits)

        return f'OK ({freq_mhz:.2f} MHz): {os.path.basename(input_fits)}'

    except Exception as e:
        return f'ERROR: {input_fits}: {e}'


# ---- Main -------------------------------------------------------------------

def main():

    # Resolve band to katbeam model name
    band = BAND[0].lower()
    if band == 'l':
        beam_model_name = 'MKAT-AA-L-JIM-2020'
        band_label = 'L-band'
    elif band == 'u':
        beam_model_name = 'MKAT-AA-UHF-JIM-2020'
        band_label = 'UHF'
    elif band == 's':
        beam_model_name = 'MKAT-AA-S-JIM-2020'
        band_label = 'S-band'
    else:
        msg(f'ERROR: Unrecognised band "{BAND}". Use L, UHF, or S.')
        sys.exit(1)

    msg(f'Band         : {band_label}')
    msg(f'Beam model   : {beam_model_name}')
    msg(f'PB cut       : {PB_CUT}')
    msg(f'Azimuthal avg: {AZAVG}')
    msg(f'Overwrite    : {OVERWRITE}')
    msg(f'Save PBCOR   : {SAVE_PBCOR}')
    msg(f'Save PB      : {SAVE_PB}')
    msg(f'Save WT      : {SAVE_WT}')

    # Glob input images
    image_list = sorted(glob.glob(IMAGE_GLOB))
    if not image_list:
        msg(f'ERROR: No images found matching: {IMAGE_GLOB}')
        sys.exit(1)
    msg(f'Found {len(image_list)} images matching: {IMAGE_GLOB}')

    # Determine number of workers
    max_workers = compute_max_workers(image_list)

    # Dispatch
    msg('')
    msg(f'Starting PB correction with {max_workers} workers ...')
    msg('')

    n_ok   = 0
    n_skip = 0
    n_err  = 0

    n_total = len(image_list)

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(pbcor_one, im, beam_model_name, FREQ_AXIS, PB_CUT,
                        AZAVG, SAVE_PBCOR, SAVE_PB, SAVE_WT, OVERWRITE,
                        i + 1, n_total): im
            for i, im in enumerate(image_list)
        }
        for future in as_completed(futures):
            try:
                result = future.result()
                msg(result)
                if result.startswith('OK'):
                    n_ok += 1
                elif result.startswith('SKIPPED'):
                    n_skip += 1
                else:
                    n_err += 1
            except Exception as e:
                msg(f'ERROR processing {futures[future]}: {e}')
                n_err += 1

    msg('')
    msg(f'Done: {n_ok} corrected, {n_skip} skipped, {n_err} errors '
        f'(total {len(image_list)} images)')


if __name__ == '__main__':
    main()
