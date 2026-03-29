#!/usr/bin/env python
# Modify WSClean model images within a spatial region prior to self-calibration predict.
#
# Usage:
#   python mod_model_selfcal.py --identifier <prefix> --stokes <I|V|VI>
#                               --spatial <file.reg|file.fits> [--zero-all-neg-I]

import argparse
import glob
import os
import re
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import FITSFixedWarning

warnings.simplefilter('ignore', FITSFixedWarning)
warnings.filterwarnings('ignore', category=FITSFixedWarning)
warnings.filterwarnings('ignore', message='.*datfix.*')
warnings.filterwarnings('ignore', module='astropy.*')


# -----------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


# -----------------------------------------------------------------------
# FITS helpers
# -----------------------------------------------------------------------

def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
        return np.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        return np.array(input_hdu.data[0, :, :])
    else:
        return np.array(input_hdu.data[0, 0, :, :])


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


# -----------------------------------------------------------------------
# DS9 region parsing
# -----------------------------------------------------------------------

def parse_ds9_regions(regfile):
    """Parse a DS9 region file; return list of region dicts (circle/ellipse)."""
    with open(regfile) as f:
        lines = [l.strip() for l in f if l.strip() and not l.strip().startswith('#')]

    coord_type = 'image'
    sky_frames = {'fk5', 'fk4', 'icrs', 'galactic', 'j2000', 'b1950', 'wcs'}
    for line in lines:
        if line.lower() in sky_frames:
            coord_type = 'sky'
            break
        if line.lower() == 'image':
            coord_type = 'image'
            break

    def parse_coord(s):
        s = s.strip()
        if ':' in s:
            parts = s.split(':')
            val = abs(float(parts[0])) + float(parts[1]) / 60.0 + float(parts[2]) / 3600.0
            if s.startswith('-'):
                val = -val
            if not s.startswith('-') and float(parts[0]) < 24:
                val *= 15.0
            return val
        return float(s)

    def parse_radius(s):
        s = s.strip()
        if coord_type == 'sky':
            if s.endswith('"'):
                return float(s[:-1])
            if s.lower().endswith('deg') or s.endswith('d'):
                return float(re.sub(r'[^\d.]', '', s)) * 3600.0
            return float(s) * 3600.0
        else:
            return float(re.sub(r'[^\d.]', '', s))

    regions = []
    for line in lines:
        is_circle  = bool(re.match(r'circle\s*\(', line, re.IGNORECASE))
        is_ellipse = bool(re.match(r'ellipse\s*\(', line, re.IGNORECASE))
        if not (is_circle or is_ellipse):
            continue
        m = re.match(r'\w+\s*\(([^)]+)\)', line, re.IGNORECASE)
        if not m:
            raise ValueError(f'Could not parse region: {line}')
        parts = [p.strip() for p in m.group(1).split(',')]
        if is_circle:
            if len(parts) != 3:
                raise ValueError(f'Expected 3 values in circle(), got: {parts}')
            x = parse_coord(parts[0]) if coord_type == 'sky' else float(parts[0])
            y = parse_coord(parts[1]) if coord_type == 'sky' else float(parts[1])
            r = parse_radius(parts[2])
        else:
            if len(parts) < 4:
                raise ValueError(f'Expected at least 4 values in ellipse(), got: {parts}')
            x  = parse_coord(parts[0]) if coord_type == 'sky' else float(parts[0])
            y  = parse_coord(parts[1]) if coord_type == 'sky' else float(parts[1])
            r  = max(parse_radius(parts[2]), parse_radius(parts[3]))
        regions.append({'coord_type': coord_type, 'x': x, 'y': y, 'r': r,
                        'shape': 'circle' if is_circle else 'ellipse'})

    if not regions:
        raise ValueError(f'No circle() or ellipse() definitions found in {regfile}')
    return regions


# -----------------------------------------------------------------------
# Spatial input loader — auto-detects .reg vs .fits
# -----------------------------------------------------------------------

def load_spatial(path):
    """
    Returns (mode, data):
      mode == 'fits' : data is the path string; mask loaded per-worker
      mode == 'reg'  : data is a list of region dicts
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.fits', '.fit'):
        return ('fits', path)
    else:
        return ('reg', parse_ds9_regions(path))


# -----------------------------------------------------------------------
# Mask construction
# -----------------------------------------------------------------------

def build_pixel_mask_from_region(fitsfile, region):
    with fits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        ny, nx = hdr['NAXIS2'], hdr['NAXIS1']
    if region['coord_type'] == 'sky':
        with fits.open(fitsfile) as hdul:
            wcs = WCS(hdul[0].header).celestial
        cx, cy = wcs.all_world2pix(region['x'], region['y'], 0)
        pixscale_arcsec = abs(wcs.pixel_scale_matrix[0, 0]) * 3600.0
        r_pix = region['r'] / pixscale_arcsec
    else:
        cx, cy, r_pix = region['x'] - 1, region['y'] - 1, region['r']
    yy, xx = np.ogrid[:ny, :nx]
    return (xx - cx) ** 2 + (yy - cy) ** 2 <= r_pix ** 2


def build_combined_mask(fitsfile, spatial_mode, spatial_data):
    if spatial_mode == 'fits':
        mask_img  = get_image(spatial_data)
        model_img = get_image(fitsfile)
        if mask_img.shape != model_img.shape:
            raise ValueError(
                f'FITS mask shape {mask_img.shape} != model shape {model_img.shape}: {fitsfile}')
        return mask_img.astype(bool)
    else:
        with fits.open(fitsfile) as hdul:
            hdr = hdul[0].header
            ny, nx = hdr['NAXIS2'], hdr['NAXIS1']
        combined = np.zeros((ny, nx), dtype=bool)
        for region in spatial_data:
            combined |= build_pixel_mask_from_region(fitsfile, region)
        return combined


# -----------------------------------------------------------------------
# Memory-aware worker count
# -----------------------------------------------------------------------

def estimate_memory_per_worker_bytes(fitsfile):
    with fits.open(fitsfile) as hdul:
        return 4 * hdul[0].data.nbytes + 100 * 1024 ** 2


def get_available_memory_bytes():
    try:
        import psutil
        return psutil.virtual_memory().available
    except ImportError:
        pass
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    msg('WARNING: Could not determine available memory. Assuming 4 GB.')
    return 4 * 1024 ** 3


def compute_max_workers(fitslist):
    available  = get_available_memory_bytes()
    per_worker = estimate_memory_per_worker_bytes(fitslist[0])
    mem_cap    = max(1, int(available // per_worker))
    # --- SLURM-aware CPU cap ---
    # SLURM_CPUS_ON_NODE is intentionally excluded: it reports total node CPUs,
    # not the allocation. We multiply cpus_per_task x ntasks when both are set.
    import re as _re
    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _job_cpus      = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    _job_cpus_m    = _re.match(r'(\d+)', _job_cpus)

    if _cpus_per_task and _ntasks:
        cpu_cap    = int(_cpus_per_task) * int(_ntasks)
        _cpu_src   = f'SLURM_CPUS_PER_TASK={_cpus_per_task} x SLURM_NTASKS={_ntasks}'
    elif _cpus_per_task:
        cpu_cap    = int(_cpus_per_task)
        _cpu_src   = f'SLURM_CPUS_PER_TASK={_cpus_per_task}'
    elif _job_cpus_m:
        cpu_cap    = int(_job_cpus_m.group(1))
        _cpu_src   = f'SLURM_JOB_CPUS_PER_NODE={_job_cpus} -> {cpu_cap}'
    elif _ntasks:
        cpu_cap    = int(_ntasks)
        _cpu_src   = f'SLURM_NTASKS={_ntasks}'
    else:
        cpu_cap    = os.cpu_count() or 1
        _cpu_src   = 'os.cpu_count() (no SLURM allocation detected)'
    workers    = min(mem_cap, cpu_cap, len(fitslist))
    msg(f'Memory available   : {available / 1024**3:.2f} GB  [{"SLURM" if any(os.environ.get(v) for v in ("SLURM_MEM_PER_NODE","SLURM_MEM_PER_CPU")) else "psutil/proc"}]')
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={len(fitslist)})]')
    return workers


# -----------------------------------------------------------------------
# Per-file workers
# -----------------------------------------------------------------------

def process_stokes_I(fitsfile, spatial_mode, spatial_data):
    """Zero negative pixels inside the spatial region in a Stokes I model image."""
    mask = build_combined_mask(fitsfile, spatial_mode, spatial_data)
    img  = get_image(fitsfile)
    neg  = mask & (img < 0)
    n    = np.count_nonzero(neg)
    if n == 0:
        return f'I skipped (no negative pixels in region): {fitsfile}'
    img[neg] = 0.0
    flush_fits(img, fitsfile)
    return f'I zeroed {n} negative pixel(s) in region: {fitsfile}'


def process_stokes_I_all(fitsfile):
    """Zero ALL negative pixels in a Stokes I model image, regardless of region."""
    img = get_image(fitsfile)
    neg = img < 0
    n   = np.count_nonzero(neg)
    if n == 0:
        return f'I skipped (no negative pixels anywhere): {fitsfile}'
    img[neg] = 0.0
    flush_fits(img, fitsfile)
    return f'I zeroed {n} negative pixel(s) globally: {fitsfile}'


def process_stokes_V(fitsfile, spatial_mode, spatial_data):
    """Zero all non-zero pixels inside the spatial region in a Stokes V model image."""
    mask = build_combined_mask(fitsfile, spatial_mode, spatial_data)
    img  = get_image(fitsfile)
    n    = np.count_nonzero(img[mask])
    if n == 0:
        return f'V skipped (no non-zero pixels in region): {fitsfile}'
    img[mask] = 0.0
    flush_fits(img, fitsfile)
    return f'V zeroed {n} pixel(s) in region: {fitsfile}'


# -----------------------------------------------------------------------
# Parallel runner
# -----------------------------------------------------------------------

def run_parallel(fitslist, worker_fn, worker_args):
    max_workers = compute_max_workers(fitslist)
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(worker_fn, f, *worker_args): f for f in fitslist}
        for future in as_completed(futures):
            try:
                msg(future.result())
            except Exception as e:
                msg(f'ERROR processing {futures[future]}: {e}')


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Modify WSClean model images within a spatial region prior to selfcal predict.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Stokes modes:
  I   Zero negative pixels in *-I-*model.fits inside spatial input
  V   Zero all non-zero pixels in *-V-*model.fits inside spatial input
  VI  Both I and V

Spatial input (--spatial):
  .reg   DS9 region file (circle/ellipse, image or sky coords, auto-detected)
  .fits  Binary mask FITS — must match model image dimensions exactly
        """)

    parser.add_argument('--identifier',     required=True,
                        help='Glob prefix to find *model.fits files')
    parser.add_argument('--stokes',         required=True,
                        help='Stokes to modify: I, V, or VI')
    parser.add_argument('--spatial',        required=True,
                        help='Spatial input: DS9 .reg file or binary mask .fits file')
    parser.add_argument('--zero-all-neg-I', action='store_true', default=False,
                        help='Zero ALL negative Stokes I pixels globally (ignores --spatial)')

    args = parser.parse_args()

    identifier     = args.identifier
    stokes         = args.stokes.upper()
    spatial_path   = args.spatial
    zero_all_neg_I = args.zero_all_neg_I

    if stokes not in ('I', 'V', 'IV', 'VI'):
        parser.error(f'--stokes must be I, V, or VI — got "{stokes}"')
    if not os.path.exists(spatial_path):
        parser.error(f'--spatial file not found: {spatial_path}')

    try:
        spatial_mode, spatial_data = load_spatial(spatial_path)
    except ValueError as e:
        sys.exit(f'ERROR loading spatial input: {e}')

    msg(f'Stokes mode      : {stokes}')
    msg(f'Zero all neg I   : {zero_all_neg_I}')
    msg(f'Spatial input    : {spatial_path} ({spatial_mode})')
    if spatial_mode == 'reg':
        for i, region in enumerate(spatial_data):
            shape  = region.get('shape', 'circle')
            suffix = ' [larger semi-axis used]' if shape == 'ellipse' else ''
            if region['coord_type'] == 'sky':
                msg(f'  Region {i+1} ({shape}): RA={region["x"]:.6f} deg, '
                    f'Dec={region["y"]:.6f} deg, r={region["r"]:.2f} arcsec{suffix}')
            else:
                msg(f'  Region {i+1} ({shape}): x={region["x"]:.2f} px, '
                    f'y={region["y"]:.2f} px, r={region["r"]:.2f} px{suffix}')

    # ---- Stokes I -------------------------------------------------------
    if 'I' in stokes:
        pattern_I  = identifier + '*-I-*model.fits'
        fitslist_I = sorted(glob.glob(pattern_I))
        if not fitslist_I:
            msg(f'WARNING: No Stokes I model images found matching: {pattern_I}')
        else:
            msg(f'Found {len(fitslist_I)} Stokes I model image(s)')
            if zero_all_neg_I:
                msg('Stokes I mode: zeroing ALL negative pixels globally')
                run_parallel(fitslist_I, process_stokes_I_all, [])
            else:
                msg('Stokes I mode: zeroing negative pixels inside region only')
                run_parallel(fitslist_I, process_stokes_I, [spatial_mode, spatial_data])

    # ---- Stokes V -------------------------------------------------------
    if 'V' in stokes:
        pattern_V  = identifier + '*-V-*model.fits'
        fitslist_V = sorted(glob.glob(pattern_V))
        if not fitslist_V:
            msg(f'WARNING: No Stokes V model images found matching: {pattern_V}')
        else:
            msg(f'Found {len(fitslist_V)} Stokes V model image(s)')
            run_parallel(fitslist_V, process_stokes_V, [spatial_mode, spatial_data])

    msg('Done.')


if __name__ == '__main__':
    main()
