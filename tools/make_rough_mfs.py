#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import os
import re
import sys
import glob
import time
import numpy as np
import bottleneck as bn
from astropy.io import fits
from concurrent.futures import ProcessPoolExecutor, as_completed

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


# Peak memory per worker, in units of one full channel stack (n_chan x ny x nx float32).
# The stack itself, plus the single working copy bn.nanmedian partitions in place.
STACK_ARRAY_EQUIVALENTS = 2
PROCESS_OVERHEAD_BYTES = 100 * 1024 ** 2  # 100 MB fixed per process

# Channels whose BMAJ is at or below this are treated as empty and excluded from the median
MIN_VALID_BMAJ = 1e-14


# Helper functions
def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush=True)

# Fits manipulation functions
def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
            image = np.array(input_hdu.data[:,:])
    elif len(input_hdu.data.shape) == 3:
            image = np.array(input_hdu.data[0,:,:])
    else:
            image = np.array(input_hdu.data[0,0,:,:])
    return image


# Resource detection
def get_available_memory_bytes():
    """Return available memory in bytes, preferring SLURM allocation over node-wide figures."""
    slurm_mem_node = os.environ.get('SLURM_MEM_PER_NODE')
    if slurm_mem_node:
        try:
            allocated = int(slurm_mem_node) * 1024 ** 2  # MB -> bytes
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
            allocated = int(slurm_mem_cpu) * _alloc_cpus * 1024 ** 2  # MB -> bytes
            msg(f'SLURM memory detected: SLURM_MEM_PER_CPU={slurm_mem_cpu} MB x {_alloc_cpus} CPUs -> {allocated / 1024**3:.1f} GB')
            return allocated
        except ValueError:
            pass

    if HAS_PSUTIL:
        available = psutil.virtual_memory().available
        msg(f'No SLURM memory allocation detected — using node available memory: {available / 1024**3:.1f} GB')
        return available
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    available = int(line.split()[1]) * 1024
                    msg(f'No SLURM memory allocation detected — using /proc/meminfo: {available / 1024**3:.1f} GB')
                    return available
    except OSError:
        pass
    msg('WARNING: Could not determine available memory (install psutil for accuracy). Assuming 4 GB.')
    return 4 * 1024 ** 3


def estimate_memory_per_worker_bytes(identifier, stoke):
    """
    Estimate peak memory for one (identifier, stokes) job: the full float32
    channel stack times STACK_ARRAY_EQUIVALENTS, plus a fixed process overhead.
    """
    images = get_channel_images(identifier, stoke)
    if not images:
        return PROCESS_OVERHEAD_BYTES
    header = fits.getheader(images[0])
    plane_bytes = header['NAXIS1'] * header['NAXIS2'] * 4
    return STACK_ARRAY_EQUIVALENTS * len(images) * plane_bytes + PROCESS_OVERHEAD_BYTES


def compute_max_workers(jobs):
    """
    Cap worker count by available memory (using the largest job's estimate),
    allocated CPUs, and the number of jobs. Returns at least 1.
    """
    available = get_available_memory_bytes()
    per_worker = max(estimate_memory_per_worker_bytes(ident, stoke) for ident, stoke in jobs)

    mem_cap = max(1, int(available // per_worker))
    # SLURM_CPUS_ON_NODE is excluded: it reports total node CPUs, not the allocation.
    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _job_cpus      = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    _job_cpus_m    = re.match(r'(\d+)', _job_cpus)

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
    n_jobs = len(jobs)

    workers = min(mem_cap, cpu_cap, n_jobs)

    msg(f'Memory available   : {available / 1024**3:.2f} GB')
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB  [largest channel stack]')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, jobs={n_jobs})]')

    return workers


# Naming and job construction
def get_identifiers(prefix):
    """
    Read in the input prefix and return the image name identifiers
    for each image (or set of images).

    Arguments:
        prefix: a string of the input prefix for the wsclean imaging call
    Returns:
        identifiers: an array containing the identifiers, where each identifier
        corresponds to a set of channelized images that are combined into a
        single MFS image
    """
    msg('Getting naming identifier "groups"')
    image_arr = sorted(glob.glob(f'{prefix}-t*'))
    suffix = np.unique([im.split(f'{prefix}-t')[-1].split('-')[0] for im in image_arr])
    return [f'{prefix}-t{s}' for s in suffix]


def get_stokes(identifier):
    """Return the Stokes suffixes present for an identifier, or [''] if none are labelled."""
    stokes = [s for s in ['-I', '-Q', '-U', '-V'] if glob.glob(f'{identifier}*{s}-*')]
    return stokes or ['']


def get_channel_images(identifier, stoke):
    return sorted(glob.glob(f'{identifier}-[!MFS]*{stoke}-image.fits'))


def get_mfs_name(identifier, stoke):
    return f'{identifier}-MFS{stoke}-image.homogenized.fits'


def build_jobs(prefix):
    """
    Return (identifier, stokes) pairs that need a pseudo-MFS image.
    Identifiers with a WSClean MFS PSF (shared across all Stokes) already have
    a real MFS and are skipped, as are Stokes whose pseudo-MFS output exists.
    """
    jobs = []
    for identifier in get_identifiers(prefix):
        if glob.glob(f'{identifier}-MFS-psf.fits'):
            msg(f'Skipped (WSClean MFS exists): {identifier}')
            continue
        for stoke in get_stokes(identifier):
            if os.path.exists(get_mfs_name(identifier, stoke)):
                msg(f'Skipped (pseudo-MFS exists): {get_mfs_name(identifier, stoke)}')
                continue
            jobs.append((identifier, stoke))
    return jobs


# Worker
def make_mfs(identifier, stoke):
    """
    Worker function: make a (pseudo-) MFS image for one identifier and Stokes
    parameter by taking the per-pixel median over channels with a valid beam.
    """
    valid = []
    for im in get_channel_images(identifier, stoke):
        header = fits.getheader(im)
        if header.get('BMAJ', 0) > MIN_VALID_BMAJ:
            valid.append((im, header))

    if not valid:
        return f'Skipped (no channels with a valid beam): {identifier}{stoke}'

    header = valid[0][1].copy()
    stack = np.empty((len(valid), header['NAXIS2'], header['NAXIS1']), dtype=np.float32)
    for k, (im, _) in enumerate(valid):
        stack[k] = get_image(im)

    # Adopt median values for each pixel and output MFS image
    data = bn.nanmedian(stack, axis=0)
    del stack
    for _ in range(4 - len(data.shape)): # Make (1, 1, N, N)-shaped
        data = data[None, :]

    header['CRVAL3'] = np.nanmean([h['CRVAL3'] for _, h in valid])
    header['BMAJ'] = np.nanmean([h['BMAJ'] for _, h in valid])
    header['BMIN'] = np.nanmean([h['BMIN'] for _, h in valid])
    header['BPA'] = np.nanmedian([h['BPA'] for _, h in valid])
    mfs_name = get_mfs_name(identifier, stoke)
    fits.PrimaryHDU(data=data, header=header).writeto(mfs_name, overwrite=True)
    return f'Done ({len(valid)} channels): {mfs_name}'


def main():

    # Read in prefix return error if missing it
    if len(sys.argv) != 2:
        msg('ERROR: Missing image prefix')
        sys.exit()
    prefix = sys.argv[-1]

    jobs = build_jobs(prefix)
    if not jobs:
        msg('No pseudo-MFS images to make')
        return

    msg(f'Making {len(jobs)} (pseudo-) MFS image(s); take these images with less than a grain of salt')
    max_workers = compute_max_workers(jobs)

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(make_mfs, identifier, stoke): (identifier, stoke)
            for identifier, stoke in jobs
        }
        for future in as_completed(futures):
            identifier, stoke = futures[future]
            try:
                msg(future.result())
            except Exception as e:
                msg(f'ERROR processing {identifier}{stoke}: {e}')

if __name__ == '__main__':
    main()
