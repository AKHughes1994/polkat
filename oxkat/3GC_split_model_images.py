#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
#
# Split WSClean model images into direction-dependent components using DS9
# region files.  For each per-channel model FITS file, pixels inside the
# region(s) are extracted into a new "direction 1" model, and optionally
# the complement (everything outside the region) is written as a
# "subtracted" model for use in peeling / DD calibration.
#
# Supports circle and ellipse DS9 regions.  Ellipses are conservatively
# treated as circles using the larger semi-axis.
#
# Parallelised with ProcessPoolExecutor; SLURM-aware resource capping.
#
# Usage:
#   python 3GC_split_model_images.py --region target.reg --prefix img_prefix [--subtract]


import glob
import numpy
import os
import re
import shutil
import time
from astropy import wcs
from astropy.io import fits
from concurrent.futures import ProcessPoolExecutor, as_completed
from optparse import OptionParser


# ---------------------------------------------------------------------------------------
# Coordinate conversion utilities
# ---------------------------------------------------------------------------------------


def msg(txt):
    """Timestamped log message."""
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


def hms2deg(hms, delimiter=':'):
    """Convert a right ascension string (HH:MM:SS.SS) to decimal degrees."""
    h, m, s = hms.split(delimiter)
    deg = 15.0 * (float(h) + float(m) / 60.0 + float(s) / 3600.0)
    return deg


def dms2deg(dms, delimiter=':'):
    """Convert a declination string (±DD:MM:SS.SS) to decimal degrees."""
    d, m, s = dms.split(delimiter)
    decsign = -1.0 if d[0] == '-' else 1.0
    d = float(d.lstrip('+-'))
    deg = decsign * (d + float(m) / 60.0 + float(s) / 3600.0)
    return deg


def radius2deg(radius):
    """
    Convert a DS9 angular size string to decimal degrees.
    Accepts arcsec ("), arcmin ('), or bare float (assumed degrees).
    """
    if radius[-1] == '"':
        return float(radius[:-1]) / 3600.0
    elif radius[-1] == "'":
        return float(radius[:-1]) / 60.0
    else:
        return float(radius)


# ---------------------------------------------------------------------------------------
# DS9 region file parser
# ---------------------------------------------------------------------------------------


def process_region_file(region_file):
    """
    Parse a DS9 region file and return a list of (RA, Dec, radius) tuples
    in decimal degrees.

    Supported shapes:
      - circle(RA, Dec, radius)
      - ellipse(RA, Dec, semi_a, semi_b, angle)
        --> converted to a circle using max(semi_a, semi_b)

    RA/Dec can be in HMS/DMS (colon-delimited) or decimal degrees.
    Radius/semi-axes can have arcsec/arcmin suffixes.
    """
    circles = []

    f = open(region_file, 'r')
    line = f.readline()
    while line:
        # Identify the shape type from the start of the line
        shape = None
        if line.strip().startswith('circle'):
            shape = 'circle'
        elif line.strip().startswith('ellipse'):
            shape = 'ellipse'

        if shape is not None:
            # Strip whitespace and parentheses to isolate the parameter list
            line = line.replace(' ', '')
            line = line.rstrip('\n').replace('(', ' ').replace(')', ' ')
            params = line.split()[1].split(',')

            # Parse RA: HMS string or decimal degrees
            ra = params[0]
            ra = hms2deg(ra) if ':' in ra else float(ra)

            # Parse Dec: DMS string or decimal degrees
            dec = params[1]
            dec = dms2deg(dec) if ':' in dec else float(dec)

            if shape == 'circle':
                radius = radius2deg(params[2])
            elif shape == 'ellipse':
                # DS9 ellipse format: ellipse(RA, Dec, semi_a, semi_b, PA)
                # We take the larger semi-axis as a conservative circular mask
                semi_a = radius2deg(params[2])
                semi_b = radius2deg(params[3])
                radius = max(semi_a, semi_b)
                msg(f'Ellipse -> circle: semi_a={semi_a*3600:.2f}", '
                    f'semi_b={semi_b*3600:.2f}", using r={radius*3600:.2f}"')

            circles.append((ra, dec, radius))

        line = f.readline()
    f.close()

    return circles


# ---------------------------------------------------------------------------------------
# FITS I/O helpers
# ---------------------------------------------------------------------------------------


def get_image(fits_file):
    """
    Read the 2D image plane from a FITS file, handling 2D/3D/4D data cubes.
    WSClean models are typically 4D (Stokes, freq, dec, ra) with length-1
    leading axes.
    """
    input_hdu = fits.open(fits_file)[0]
    if len(input_hdu.data.shape) == 2:
        image = numpy.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        image = numpy.array(input_hdu.data[0, :, :])
    else:
        image = numpy.array(input_hdu.data[0, 0, :, :])
    return image


def flush_fits(image, fits_file):
    """
    Write a 2D numpy array back into an existing FITS file in-place,
    preserving the original header and data shape.
    """
    f = fits.open(fits_file, mode='update')
    input_hdu = f[0]
    if len(input_hdu.data.shape) == 2:
        input_hdu.data[:, :] = image
    elif len(input_hdu.data.shape) == 3:
        input_hdu.data[0, :, :] = image
    else:
        input_hdu.data[0, 0, :, :] = image
    f.flush()
    f.close()


# ---------------------------------------------------------------------------------------
# Masking
# ---------------------------------------------------------------------------------------


def apply_circle(image, xpix, ypix, rpix):
    """
    Set all pixels within a circle of radius rpix (in pixels) centred on
    (xpix, ypix) to 1.0 in the mask array.  Uses a bounding-box + distance
    check to avoid iterating over the full image.
    """
    # Define a bounding box around the circle to limit the pixel loop
    xg, yg = numpy.mgrid[int(xpix - rpix):int(xpix + rpix) + 1,
                          int(ypix - rpix):int(ypix + rpix) + 1]
    xg = xg.ravel()
    yg = yg.ravel()

    # Only set pixels whose Euclidean distance from centre < rpix
    for i, j in zip(xg, yg):
        sep = ((i - xpix)**2.0 + (j - ypix)**2.0)**0.5
        if sep < rpix:
            image[j, i] = 1.0

    return image


def fmt(xx):
    """Format a float to 5 decimal places for logging."""
    return str(round(xx, 5))


def spacer():
    """Print a horizontal separator line."""
    print('--------------|---------------------------------------------')


# ---------------------------------------------------------------------------------------
# SLURM-aware resource calculation
# ---------------------------------------------------------------------------------------


def get_available_memory_bytes():
    """
    Determine available memory in bytes.  Checks SLURM allocation variables
    first (SLURM_MEM_PER_NODE, SLURM_MEM_PER_CPU), then falls back to
    psutil or /proc/meminfo.  Returns the more restrictive of SLURM vs
    system-reported available memory.
    """
    slurm_mem = None
    slurm_src = None

    # --- Check SLURM memory allocation ---
    mem_per_node = os.environ.get('SLURM_MEM_PER_NODE')
    mem_per_cpu  = os.environ.get('SLURM_MEM_PER_CPU')

    if mem_per_node:
        # SLURM_MEM_PER_NODE: bare number = MB, or with G/T suffix
        val = mem_per_node.rstrip('MmGgTt')
        multiplier = 1024 ** 3 if mem_per_node[-1].lower() in ('g', 't') else 1024 ** 2
        if mem_per_node[-1].lower() == 't':
            multiplier = 1024 ** 4
        slurm_mem = int(val) * multiplier
        slurm_src = f'SLURM_MEM_PER_NODE={mem_per_node}'

    elif mem_per_cpu:
        # SLURM_MEM_PER_CPU: per-CPU allocation, multiply by number of CPUs
        val = mem_per_cpu.rstrip('MmGgTt')
        multiplier = 1024 ** 3 if mem_per_cpu[-1].lower() in ('g', 't') else 1024 ** 2
        if mem_per_cpu[-1].lower() == 't':
            multiplier = 1024 ** 4
        per_cpu_bytes = int(val) * multiplier

        # Determine total allocated CPUs from SLURM variables
        cpus   = os.environ.get('SLURM_CPUS_PER_TASK')
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
            with open('/proc/meminfo') as mf:
                for line in mf:
                    if line.startswith('MemTotal:') and sys_total is None:
                        sys_total = int(line.split()[1]) * 1024
                    if line.startswith('MemAvailable:') and sys_avail is None:
                        sys_avail = int(line.split()[1]) * 1024
        except OSError:
            pass

    # --- Log and return the more restrictive value ---
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
    """
    Conservative estimate of peak memory per worker process.
    Each worker holds: original image, mask, masked product, and possibly
    the subtracted product — ~4 copies, plus ~100 MB process overhead.
    """
    with fits.open(fitsfile) as hdul:
        return 4 * hdul[0].data.nbytes + 100 * 1024 ** 2


def compute_max_workers(image_list):
    """
    Determine the maximum number of parallel workers, constrained by:
      1. Available memory / estimated memory per worker
      2. Allocated CPU count (SLURM-aware)
      3. Number of images to process
    """
    available  = get_available_memory_bytes()
    per_worker = estimate_memory_per_worker_bytes(image_list[0])
    mem_cap    = max(1, int(available // per_worker))

    # --- SLURM-aware CPU cap ---
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

    # Take the most restrictive of all three caps
    workers = min(mem_cap, cpu_cap, len(image_list))
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={len(image_list)})]')
    return workers


# ---------------------------------------------------------------------------------------
# Per-image worker function (runs in subprocess)
# ---------------------------------------------------------------------------------------


def process_one(fits_file, circles, model_pattern, suffix, subtract, idx, n_total):
    """
    Process a single per-channel model FITS file:
      1. Read the model image
      2. Build a binary mask from the DS9 regions (1 inside, 0 outside)
      3. Write image * mask  as the direction-dependent model
      4. Optionally write image * (1 - mask)  as the subtracted model

    Each worker is fully self-contained: reads its own WCS, applies all
    regions, and writes output independently.  No shared state between workers.

    Parameters
    ----------
    fits_file     : str  -- path to the input model FITS file
    circles       : list -- list of (RA, Dec, radius) tuples in degrees
    model_pattern : str  -- wsclean prefix used to construct output filenames
    suffix        : str  -- region file basename (sans extension), appended to output
    subtract      : bool -- if True, also write the complementary (subtracted) model
    idx           : int  -- 1-based index for progress logging
    n_total       : int  -- total number of images being processed

    Returns
    -------
    status : str -- 'OK: filename' or 'ERROR: filename: message'
    """
    try:
        msg(f'[{idx}/{n_total}] Processing: {os.path.basename(fits_file)}')

        # Output filename: insert the region suffix into the model prefix
        dir1_fits = fits_file.replace(model_pattern, model_pattern + '-' + suffix)

        # Read the model image and initialise a zero-valued mask of the same shape
        img  = get_image(fits_file)
        mask = img * 0.0

        # Read WCS and pixel scale for sky-to-pixel coordinate conversion
        hdulist  = fits.open(fits_file)
        w        = wcs.WCS(hdulist[0].header)
        pixscale = hdulist[0].header['CDELT2']   # degrees per pixel
        hdulist.close()

        # Stamp each region into the binary mask
        for circle in circles:
            ra, dec, radius = circle
            # Convert sky coordinates to pixel coordinates via the WCS
            coord  = (ra, dec, 0, 0)   # 4-element tuple for 4D WCS
            pixels = w.wcs_world2pix([coord], 0)
            xpix   = pixels[0][0]
            ypix   = pixels[0][1]
            rpix   = radius / pixscale   # convert angular radius to pixels
            mask   = apply_circle(mask, xpix, ypix, rpix)

        # Direction 1: model components inside the region(s)
        dir1 = img * mask
        shutil.copyfile(fits_file, dir1_fits)   # copy to preserve header
        flush_fits(dir1, dir1_fits)              # overwrite pixel data

        # Subtracted model: everything outside the region(s)
        if subtract:
            subtract_fits = fits_file.replace(
                model_pattern, model_pattern + '-' + suffix + '-subtracted')
            subt = img * (1.0 - mask)
            shutil.copyfile(fits_file, subtract_fits)
            flush_fits(subt, subtract_fits)

        return f'OK: {os.path.basename(fits_file)}'

    except Exception as e:
        return f'ERROR: {fits_file}: {e}'


# ---------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------


def main():

    # --- Parse command-line arguments ---
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('--region', dest='region_file',
                      help='DS9 region file (circles and/or ellipses)')
    parser.add_option('--prefix', dest='model_pattern',
                      help='WSClean image prefix (e.g. img_1234_target)')
    parser.add_option('--subtract', dest='subtract', default=False,
                      action='store_true',
                      help='Also produce model with region components subtracted')
    (options, args) = parser.parse_args()

    region_file   = options.region_file
    model_pattern = options.model_pattern
    subtract      = options.subtract

    # --- Parse the DS9 region file ---
    circles = process_region_file(region_file)

    # Use the region filename (without extension) as the output suffix
    suffix = region_file.split('/')[-1].split('.')[0]

    spacer()
    msg('DS9 region    : ' + region_file)
    msg('Contains      : ' + str(len(circles)) + ' regions (circles + ellipses)')
    msg('Model suffix  : ' + suffix)
    spacer()

    # --- Glob for per-channel model images ---
    # WSClean names channels as -0000-, -0001-, etc. in the model files
    model_list = sorted(glob.glob(model_pattern + '-0*model*fits'))

    if not model_list:
        msg(f'ERROR: No model images found matching: {model_pattern}-0*model*fits')
        return

    msg(f'Found {len(model_list)} model images')

    # --- Log sky-to-pixel mapping for the first image as a sanity check ---
    msg('')
    msg('Region verification (first image):')
    hdulist  = fits.open(model_list[0])
    w        = wcs.WCS(hdulist[0].header)
    pixscale = hdulist[0].header['CDELT2']
    hdulist.close()
    for circle in circles:
        ra, dec, radius = circle
        coord  = (ra, dec, 0, 0)
        pixels = w.wcs_world2pix([coord], 0)
        xpix   = pixels[0][0]
        ypix   = pixels[0][1]
        rpix   = radius / pixscale
        msg(f'  sky {fmt(ra)} {fmt(dec)} {fmt(radius)}'
            f'  ->  pixel {fmt(xpix)} {fmt(ypix)} {fmt(rpix)}')
    msg('')

    # --- Determine parallelism (memory + CPU aware) ---
    max_workers = compute_max_workers(model_list)

    msg('')
    msg(f'Starting model splitting with {max_workers} workers ...')
    msg('')

    # --- Dispatch all images to the worker pool ---
    n_total = len(model_list)
    n_ok  = 0
    n_err = 0

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        # Submit one future per model image
        futures = {
            pool.submit(process_one, im, circles, model_pattern, suffix,
                        subtract, i + 1, n_total): im
            for i, im in enumerate(model_list)
        }
        # Collect results as they complete (order is non-deterministic)
        for future in as_completed(futures):
            try:
                result = future.result()
                msg(result)
                if result.startswith('OK'):
                    n_ok += 1
                else:
                    n_err += 1
            except Exception as e:
                msg(f'ERROR processing {futures[future]}: {e}')
                n_err += 1

    msg('')
    msg(f'Done: {n_ok} processed, {n_err} errors (total {n_total} images)')


if __name__ == '__main__':

    main()
