import glob, os, subprocess, sys, time
from astropy.io import fits
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import os.path as o

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


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


def estimate_memory_per_worker_bytes(image_Q):
    """
    Estimate peak memory consumed by one worker process for a given Q image.

    Per worker we load up to 3 input arrays (Q, U, V) and compute up to 2
    output arrays (Plin, Ptot) of the same shape and dtype.  astropy also
    holds an internal copy while the HDU is open, so we use a factor of 6
    array-equivalents as a conservative peak estimate, plus a fixed overhead
    for the Python process itself (~100 MB).
    """
    with fits.open(image_Q) as hdul:
        data = hdul[0].data
        array_bytes = data.nbytes  # single array footprint in bytes

    N_ARRAY_EQUIVALENTS = 6         # Q, U, V, Plin, Ptot + astropy internal copy
    PROCESS_OVERHEAD_BYTES = 100 * 1024 ** 2  # 100 MB fixed per process

    return N_ARRAY_EQUIVALENTS * array_bytes + PROCESS_OVERHEAD_BYTES


def get_available_memory_bytes():
    """Return available memory in bytes, preferring SLURM allocation over node-wide figures."""
    # Prefer SLURM allocation over node-wide figures so we don't over-subscribe
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
    # Fallback: parse /proc/meminfo on Linux
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


def compute_max_workers(image_list):
    """
    Cap worker count based on:
      - Available system memory vs. per-worker memory estimate
      - CPU count
    Uses the first image in the list as representative.
    Returns at least 1.
    """
    available = get_available_memory_bytes()
    per_worker = estimate_memory_per_worker_bytes(image_list[0])

    mem_cap = max(1, int(available // per_worker))
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
    n_images = len(image_list)

    workers = min(mem_cap, cpu_cap, n_images)

    msg(f'Memory available   : {available / 1024**3:.2f} GB  [{"SLURM" if any(os.environ.get(v) for v in ("SLURM_MEM_PER_NODE","SLURM_MEM_PER_CPU")) else "psutil/proc"}]')
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={n_images})]')

    return workers


def process_image(image_Q, directory, plin_only):
    """Worker function: compute Plin (and optionally Ptot) for one Q image."""
    image_U    = image_Q.replace('-Q-', '-U-')
    image_V    = image_Q.replace('-Q-', '-V-')
    image_Plin = image_Q.replace('-Q-', '-Plin-')
    image_Ptot = image_Q.replace('-Q-', '-Ptot-')

    plin_exists = os.path.exists(image_Plin)
    ptot_exists = os.path.exists(image_Ptot)

    if plin_only and plin_exists:
        return f'Skipped (Plin exists): {image_Q}'
    if not plin_only and plin_exists and ptot_exists:
        return f'Skipped (Plin+Ptot exist): {image_Q}'

    # --- Plin ---
    subprocess.run([f'cp {image_Q} {image_Plin}'], shell=True, check=True)
    flux_Q = get_image(image_Q)
    flux_U = get_image(image_U)
    flush_fits((flux_Q ** 2 + flux_U ** 2) ** 0.5, image_Plin)

    # --- Ptot (optional) ---
    if not plin_only:
        subprocess.run([f'cp {image_Q} {image_Ptot}'], shell=True, check=True)
        try:
            flux_V = get_image(image_V)
            flush_fits((flux_Q ** 2 + flux_U ** 2 + flux_V ** 2) ** 0.5, image_Ptot)
        except (FileNotFoundError, OSError) as e:
            return (f'Done Plin (Ptot skipped — Stokes V unavailable: {e}): '
                    f'{image_Q}')

    return f'Done: {image_Q}'


def main():

    # --- Argument parsing ---
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        sys.exit('ERROR: Usage: python make_pol_images.py <directory> [target_name] [plin_only]\n'
                 '  directory   - Directory to parse for images (e.g., IMAGES or INTERVALS)\n'
                 '  target_name - (Optional) Target name to filter images\n'
                 '  plin_only   - (Optional) If "true", only create Plin images (not Ptot)')

    directory   = sys.argv[1]
    target_name = sys.argv[2] if len(sys.argv) >= 3 else None
    plin_only   = len(sys.argv) == 4 and sys.argv[3].lower() == 'true'

    if not os.path.isdir(directory):
        sys.exit(f'ERROR: Directory "{directory}" does not exist')

    if target_name:
        msg(f'Searching for images matching target: {target_name}')
    else:
        msg(f'Searching for all images in: {directory}')

    pattern1 = f'{directory}/*-Q-*image.fits'
    pattern2 = f'{directory}/*-Q-*image.homogenized.fits'
    image_list = sorted(glob.glob(pattern1) + glob.glob(pattern2))

    if not image_list:
        sys.exit(f'ERROR: No Q images found in {directory}')

    msg(f'Found {len(image_list)} Q image(s) to process')

    # --- Worker count ---
    max_workers = compute_max_workers(image_list)

    # --- Parallel execution ---
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(process_image, q, directory, plin_only): q
            for q in image_list
        }
        for future in as_completed(futures):
            try:
                msg(future.result())
            except Exception as e:
                msg(f'ERROR processing {futures[future]}: {e}')


if __name__ == '__main__':
    main()
