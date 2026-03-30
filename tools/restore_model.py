#!/usr/bin/env python
#andrew.hughes@physics.ox.ac.uk

import glob
import logging
import numpy
import os
import subprocess
import random 
import scipy.signal
import shutil
import time
import string
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

from astropy.io import fits
from astropy.convolution import convolve,Gaussian2DKernel
from itertools import repeat
from multiprocessing import Pool

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)


def get_available_memory_bytes():
    # Prefer SLURM allocation over node-wide figures so we don't over-subscribe
    slurm_mem_node = os.environ.get('SLURM_MEM_PER_NODE')
    if slurm_mem_node:
        try:
            allocated = int(slurm_mem_node) * 1024 ** 2  # MB -> bytes
            msg(f'SLURM memory detected: SLURM_MEM_PER_NODE={slurm_mem_node} MB -> {allocated / 1024**3:.1f} GB')
            return allocated
        except ValueError:
            pass

    slurm_mem_cpu = os.environ.get('SLURM_MEM_PER_CPU')
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
    """Conservative peak estimate: 4 array copies + 100 MB process overhead."""
    with fits.open(fitsfile) as hdul:
        return 4 * hdul[0].data.nbytes + 100 * 1024 ** 2


def compute_max_workers(image_list):
    available  = get_available_memory_bytes()
    per_worker = estimate_memory_per_worker_bytes(image_list[0])
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
    workers    = min(mem_cap, cpu_cap, len(image_list))
    msg(f'Memory available   : {available / 1024**3:.2f} GB  [{"SLURM" if any(os.environ.get(v) for v in ("SLURM_MEM_PER_NODE","SLURM_MEM_PER_CPU")) else "psutil/proc"}]')
    msg(f'Est. memory/worker : {per_worker / 1024**3:.2f} GB')
    msg(f'Workers (mem cap)  : {mem_cap}')
    msg(f'Workers (CPU cap)  : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected) : {workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={len(image_list)})]')
    return workers

def get_image(fitsfile):
        input_hdu = fits.open(fitsfile)[0]
        if len(input_hdu.data.shape) == 2:
                image = numpy.array(input_hdu.data[:,:])
        elif len(input_hdu.data.shape) == 3:
                image = numpy.array(input_hdu.data[0,:,:])
        elif len(input_hdu.data.shape) == 4:
                image = numpy.array(input_hdu.data[0,0,:,:])
        else:
                image = numpy.array(input_hdu.data[0,0,0,:,:])
        return image


def flush_fits(newimage,fitsfile):
        f = fits.open(fitsfile,mode='update')
        input_hdu = f[0]
        if len(input_hdu.data.shape) == 2:
                input_hdu.data[:,:] = newimage
        elif len(input_hdu.data.shape) == 3:
                input_hdu.data[0,:,:] = newimage_
        elif len(input_hdu.data.shape) == 4:
                input_hdu.data[0,0,:,:] = newimage
        else:
                input_hdu.data[0,0,0,:,:] = newimage
        f.flush()


def deg2rad(xx):
    return numpy.pi*xx/180.0


def get_header(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    hdr = input_hdu.header
    bmaj = hdr.get('BMAJ')
    bmin = hdr.get('BMIN')
    bpa = hdr.get('BPA')
    pixscale = hdr.get('CDELT2')
    return bmaj, bmin, bpa, pixscale


def beam_header(fitsfile,bmaj,bmin,bpa):
        outhdu = fits.open(fitsfile,mode='update')
        outhdr = outhdu[0].header
        outhdr.set('BMAJ',bmaj,after='BUNIT')
        outhdr.set('BMIN',bmin,after='BMAJ')
        outhdr.set('BPA',bpa,after='BMIN')
        outhdr.remove('HISTORY')
        outhdu.flush()  


def restore_fits(modelsub_image, model_data, psf):
                
        restored_image = modelsub_image.replace('modelsub', 'restored')
        shutil.copyfile(modelsub_image, restored_image)

        # Get the fitted beam
        bmaj, bmin, bpa, pixscale = psf

        # Create restoring beam image
        xstd = bmaj / (2.3548 * pixscale)
        ystd = bmin / (2.3548 * pixscale)
        theta = deg2rad(bpa - 90.0)
        restoring = Gaussian2DKernel(x_stddev = xstd, y_stddev = ystd, theta=theta, x_size=cropsize, y_size=cropsize, mode='center')
        restoring_beam_image = restoring.array
        restoring_beam_image = restoring_beam_image / numpy.max(restoring_beam_image)

        # Convolve model with restoring beam
        modelconv_image = scipy.signal.fftconvolve(model_data, restoring_beam_image, mode='same')

        # Open model subtracted image and add convolved model
        modelsub_data = get_image(modelsub_image)
        restored_data = modelsub_data + modelconv_image

        # Flush restored FITS file and fix the header
        flush_fits(restored_data, restored_image)
        msg('Restoring complete: ' + restored_image.split(cfg.INTERVALS + '/')[-1] + '\n')

def process_one(modelsub_image, model_images, truncate_model, modelsub_shape, cropsize, target_prefix):
    """Worker: restore a single modelsub image. Returns a status string."""
    suffix    = '-'.join(modelsub_image.split(target_prefix + '-t')[-1].split('-')[1:]).replace('image.fits', 'model.fits')
    psf_image = modelsub_image[:]
    for substring in ['-I-image.fits', '-Q-image.fits', '-U-image.fits', '-V-image.fits']:
        psf_image = psf_image.replace(substring, '-image.fits')
    psf_image = psf_image.replace('-image.fits', '-psf.fits')
    psf = get_header(psf_image)

    model_image = [im for im in model_images if im.endswith(suffix)]
    if model_image:
        model_image = model_image[0]
    else:
        model_image = [im for im in model_images if im.endswith(suffix.replace('-model.fits', '-I-model.fits'))][0]

    bmaj = psf[0]
    if bmaj is None:
        bmaj = 0.0

    if bmaj > 1.0e-14:
        subprocess.run([f'cp {psf_image} {psf_image.replace("modelsub", "restored")}'], shell=True)

        if truncate_model:
            temp_model = model_image.replace('.fits', '.tmp.fits')
            subprocess.run([f'fitstool.py -f -z {modelsub_shape[0]} -o {temp_model} {model_image}'], shell=True)
            model_data = get_image(temp_model)
            subprocess.run([f'rm -rf {temp_model}'], shell=True)
        else:
            model_data = get_image(model_image)

        restore_fits(modelsub_image, model_data, psf)
        return f'Restored: {modelsub_image}'
    else:
        return f'Skipped (flagged channel): {modelsub_image}'


if __name__ == '__main__':  
    
        if len(sys.argv) != 4:
                print('Incorrect number of inputs should just be (in order) model prefix, target image prefix, stokes parameters (I or IQUV)')
                sys.exit()

        # Inputs
        mfs_only = '--mfs-only' in sys.argv
        args = [a for a in sys.argv[1:] if not a.startswith('--')]
        model_prefix = args[0]
        target_prefix = args[1]
        pol = args[2]
       
        # Convolution parameters 
        cropsize = 51 # Size of convolving kernel in pixels 

        # Get the array of modelsub images and restoring model images
        modelsub_images = sorted(glob.glob(f'{target_prefix}*image.fits'))
        if mfs_only:
            modelsub_images = [im for im in modelsub_images if '-MFS-' in im]
        model_images = sorted(glob.glob(f'{model_prefix}*-model.fits'))

        # Check to see if the model has different dimension than the snapshot imaging
        truncate_model = False
        model_shape    = get_image(model_images[0]).shape
        modelsub_shape = get_image(modelsub_images[0]).shape
    
        # Check that images are square
        if model_shape[0] != model_shape[1] or modelsub_shape[0] != modelsub_shape[1]: 
            sys.exit('ERROR: Both model and snapshot images must be square')

        if model_shape != modelsub_shape:
            truncate_model = True

        # Parallelise restoration across modelsub images
        max_workers = compute_max_workers(modelsub_images)

        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(process_one, im, model_images, truncate_model, modelsub_shape, cropsize, target_prefix): im
                for im in modelsub_images
            }
            for future in as_completed(futures):
                try:
                    msg(future.result())
                except Exception as e:
                    msg(f'ERROR processing {futures[future]}: {e}')
