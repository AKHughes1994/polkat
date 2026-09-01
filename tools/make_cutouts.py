#!/usr/bin/env python3
"""
Script to create zoomed cutouts of WSCLEAN FITS images.
Iterates through files matching specified patterns and creates central cutouts.
"""

import os
import glob
import time
from astropy.io import fits
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


def ensure_list(item):
    """Convert string to list if needed, otherwise return as-is."""
    if isinstance(item, str):
        return [item]
    return item


def get_available_memory_bytes():
    """Return available memory in bytes, preferring SLURM allocation over node-wide figures."""
    slurm_mem_node = os.environ.get('SLURM_MEM_PER_NODE')
    if slurm_mem_node:
        try:
            allocated = int(slurm_mem_node) * 1024 ** 2
            msg(f'SLURM memory detected: SLURM_MEM_PER_NODE={slurm_mem_node} MB -> {allocated / 1024**3:.1f} GB')
            return allocated
        except ValueError:
            pass

    slurm_mem_cpu  = os.environ.get('SLURM_MEM_PER_CPU')
    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _alloc_cpus    = (int(_cpus_per_task) * int(_ntasks) if _cpus_per_task and _ntasks
                      else int(_cpus_per_task) if _cpus_per_task
                      else int(_ntasks) if _ntasks
                      else None)
    if slurm_mem_cpu and _alloc_cpus:
        try:
            allocated = int(slurm_mem_cpu) * _alloc_cpus * 1024 ** 2
            msg(f'SLURM memory detected: SLURM_MEM_PER_CPU={slurm_mem_cpu} MB x {_alloc_cpus} CPUs -> {allocated / 1024**3:.1f} GB')
            return allocated
        except ValueError:
            pass

    try:
        import psutil
        available = psutil.virtual_memory().available
        msg(f'No SLURM memory allocation detected — using node available memory: {available / 1024**3:.1f} GB')
        return available
    except ImportError:
        pass
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    available = int(line.split()[1]) * 1024
                    msg(f'No SLURM memory allocation detected — using /proc/meminfo: {available / 1024**3:.1f} GB')
                    return available
    except OSError:
        pass
    msg('WARNING: Could not determine available memory. Assuming 4 GB.')
    return 4 * 1024 ** 3


def estimate_memory_per_worker_bytes(fitsfile):
    """Conservative peak: 2 array copies (input + cutout) + 100 MB overhead."""
    with fits.open(fitsfile) as hdul:
        return 2 * hdul[0].data.nbytes + 100 * 1024 ** 2


def compute_max_workers(file_list):
    available  = get_available_memory_bytes()
    per_worker = estimate_memory_per_worker_bytes(file_list[0])
    mem_cap    = max(1, int(available // per_worker))

    import re as _re
    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _job_cpus      = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    _job_cpus_m    = _re.match(r'(\d+)', _job_cpus)

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

    workers = min(mem_cap, cpu_cap, len(file_list))
    mem_src = 'SLURM' if any(os.environ.get(v) for v in ('SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU')) else 'psutil/proc'
    msg(f'Memory available   : {available / 1024**3:.2f} GB  [{mem_src}]')
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={len(file_list)})]')
    return workers


def create_zoom_cutout(input_file, output_file, cutout_size=512):
    """
    Create a cutout of the central region of a FITS image.

    Parameters:
    -----------
    input_file : str
        Path to input FITS file
    output_file : str
        Path to output FITS file
    cutout_size : int
        Number of pixels for the square cutout (default: 512)
    """
    with fits.open(input_file) as hdul:
        header = hdul[0].header.copy()
        data   = hdul[0].data
        original_dtype = data.dtype

        shape = data.shape
        if len(shape) < 2:
            raise ValueError(f"Image must have at least 2 dimensions, got {len(shape)}")

        ny, nx = shape[-2], shape[-1]

        new_nx = min(cutout_size, nx)
        new_ny = min(cutout_size, ny)

        # Ensure even dimensions for clean centering
        if new_nx % 2 != 0:
            new_nx -= 1
        if new_ny % 2 != 0:
            new_ny -= 1

        start_x = (nx - new_nx) // 2
        start_y = (ny - new_ny) // 2

        if len(shape) == 4:
            cutout_data = data[:, :, start_y:start_y+new_ny, start_x:start_x+new_nx]
        elif len(shape) == 3:
            cutout_data = data[:, start_y:start_y+new_ny, start_x:start_x+new_nx]
        elif len(shape) == 2:
            cutout_data = data[start_y:start_y+new_ny, start_x:start_x+new_nx]
        else:
            raise ValueError(f"Unsupported number of dimensions: {len(shape)}")

        cutout_data = cutout_data.astype(original_dtype)

        if 'CRPIX1' in header and 'CRPIX2' in header:
            header['CRPIX1'] = header['CRPIX1'] - start_x
            header['CRPIX2'] = header['CRPIX2'] - start_y

        header['NAXIS1'] = new_nx
        header['NAXIS2'] = new_ny

        fits.writeto(output_file, cutout_data, header, overwrite=True)

    return output_file, cutout_data.shape


def process_one(input_file, identifier, full_dir, cutout_size, delete_originals):
    """Worker: create cutout for a single file and optionally delete the original."""
    base_name   = os.path.basename(input_file)
    output_name = base_name.replace(identifier, f'{identifier}_zoom')
    output_file = os.path.join(full_dir, output_name)

    output_file, shape = create_zoom_cutout(input_file, output_file, cutout_size)

    # Delete original unless it is an MFS image
    if delete_originals and not ('MFS-' in base_name and '-image' in base_name):
        os.remove(input_file)

    return f'Created {output_file} {shape}'


def main():

    # ================== MODIFIABLE VARIABLES ==================

    DIR              = 'IMAGES'
    identifiers      = ['datablind', 'datamask', 'diagnostic', 'pcalmask', 'uniform']
    suffixes         = ['image.fits', 'residual.fits', 'model.fits', 'dirty.fits','image.homogenized.fits']
    cutout_size      = 512
    delete_originals = True  # set False to keep original files after cutouts are made

    # ===========================================================

    identifiers = ensure_list(identifiers)
    suffixes    = ensure_list(suffixes)
    full_dir    = os.path.join(os.getcwd(), DIR)

    if not os.path.exists(full_dir):
        msg(f'WARNING: Directory {full_dir} does not exist!')
        return

    # Collect all files to process across all identifier/suffix combinations
    all_jobs = []  # list of (input_file, identifier)
    for identifier in identifiers:
        for suffix in suffixes:
            pattern       = os.path.join(full_dir, f'*{identifier}*{suffix}')
            matching      = [f for f in glob.glob(pattern) if '_zoom' not in os.path.basename(f)]
            all_jobs     += [(f, identifier) for f in matching]

    if not all_jobs:
        msg('No files found matching any pattern.')
        return

    msg(f'Found {len(all_jobs)} file(s) to process')

    max_workers = compute_max_workers([j[0] for j in all_jobs])

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(process_one, f, ident, full_dir, cutout_size, delete_originals): f
            for f, ident in all_jobs
        }
        for future in as_completed(futures):
            try:
                msg(future.result())
            except Exception as e:
                msg(f'ERROR processing {futures[future]}: {e}')

    msg('Processing complete!')


if __name__ == '__main__':
    main()
