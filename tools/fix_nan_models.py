import glob
import os
import sys
import time
import numpy
from astropy.io import fits
from concurrent.futures import ProcessPoolExecutor, as_completed


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
        image = numpy.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        image = numpy.array(input_hdu.data[0, :, :])
    else:
        image = numpy.array(input_hdu.data[0, 0, :, :])
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
    """Conservative peak: 2 array copies (read + zeroed) + 100 MB process overhead."""
    with fits.open(fitsfile) as hdul:
        return 2 * hdul[0].data.nbytes + 100 * 1024 ** 2


def compute_max_workers(fitslist):
    available  = get_available_memory_bytes()
    per_worker = estimate_memory_per_worker_bytes(fitslist[0])
    mem_cap    = max(1, int(available // per_worker))

    # --- SLURM-aware CPU cap ---
    # SLURM_CPUS_ON_NODE intentionally excluded (reflects total node, not allocation)
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

    workers = min(mem_cap, cpu_cap, len(fitslist))
    mem_src = 'SLURM' if any(os.environ.get(v) for v in ('SLURM_MEM_PER_NODE', 'SLURM_MEM_PER_CPU')) else 'psutil/proc'
    msg(f'Memory available   : {available / 1024**3:.2f} GB  [{mem_src}]')
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={len(fitslist)})]')
    return workers


def process_one(fitsfile):
    """Check and fix NaN values in a single model FITS file."""
    img    = get_image(fitsfile)
    maxval = numpy.max(img)
    if numpy.isnan(maxval):
        new_img = numpy.zeros((img.shape[0], img.shape[1]))
        flush_fits(new_img, fitsfile)
        return f'Zeroed NaN model: {fitsfile}'
    return f'OK (max={maxval:.6g}): {fitsfile}'


if __name__ == '__main__':

    if len(sys.argv) != 2:
        sys.exit('Usage: python fix_nan_models.py <image_prefix>')

    pattern  = sys.argv[1]
    fitslist = sorted(glob.glob(pattern + '*model.fits'))

    if not fitslist:
        sys.exit(f'No model FITS files found matching: {pattern}*model.fits')

    msg(f'Found {len(fitslist)} model file(s) to check')

    max_workers = compute_max_workers(fitslist)

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(process_one, f): f for f in fitslist}
        for future in as_completed(futures):
            try:
                msg(future.result())
            except Exception as e:
                msg(f'ERROR processing {futures[future]}: {e}')
