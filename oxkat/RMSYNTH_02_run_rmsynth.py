import glob,os,datetime, subprocess, sys, json, shutil
import numpy as np
import os.path as o
import time
from functools import partial
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg
from oxkat import generate_jobs as gen

# Reinitialize print to flush in real-time for SLURM logging
print = partial(print, flush=True)

def mad_clip_rmsynth_file(fname, nsigma=5.0):
    '''
    Drop badly corrupted channels from an RM synthesis input file
    (freq, I, Q, U, dI, dQ, dU), with the same tests as the RMSYNTH_01
    spectrum plots:
      1. any of dI, dQ, dU more than nsigma robust sigmas (1.4826 x MAD) from
         its median, or non-finite
      2. Stokes I more than RMSYN_SPEC_INDEX_MAD_CLIP robust sigmas from the
         median, repeated until nothing more is clipped (Stokes I only: Q and
         U oscillate with Faraday rotation, so a clip on their values would
         remove real signal)

    The full spectrum is kept as <base>_unclipped.txt and the clip always
    starts from it, so re-running never clips twice. The copy is refreshed
    whenever fname is newer (i.e. RMSYNTH_01 has rewritten it).
    '''
    backup = fname[:-len('.txt')] + '_unclipped.txt'
    if not os.path.isfile(backup) or os.path.getmtime(fname) > os.path.getmtime(backup):
        shutil.copy(fname, backup)
    data = np.loadtxt(backup, ndmin=2)

    good = np.ones(len(data), dtype=bool)
    for rms in data[:, 4:7].T:
        med = np.nanmedian(rms)
        sig = 1.4826 * np.nanmedian(np.abs(rms - med))
        if sig > 0:
            good &= np.abs(rms - med) <= nsigma * sig
    n_rms = int(np.sum(~good))

    while np.sum(good) >= 3:
        I = data[good, 1]
        med = np.nanmedian(I)
        sig = 1.4826 * np.nanmedian(np.abs(I - med))
        if not sig > 0:
            break
        clip = np.abs(I - med) > cfg.RMSYN_SPEC_INDEX_MAD_CLIP * sig
        if not clip.any():
            break
        good[np.flatnonzero(good)[clip]] = False

    np.savetxt(fname, data[good])
    os.utime(backup)   # keep the copy newer than the clipped file
    print(gen.col('MAD clip')+f'{os.path.basename(fname)}: dropped {n_rms} (rms) + '
          f'{int(np.sum(~good)) - n_rms} (Stokes I) of {len(good)} channel(s)')

def main():

    # The usual rmsynth input files, plus the U/V-merged companions that
    # RMSYNTH_01 writes when there is no polarization angle calibrator. There
    # are none for Stokes I only extractions.
    fnames = sorted(glob.glob(cfg.RESULTS + '/*_rmsynth.txt') + glob.glob(cfg.RESULTS + '/*_rmsynth_UVfix.txt'))

    if len(fnames) == 0:
        print(gen.col('RM Synthesis')+'No *_rmsynth.txt or *_rmsynth_UVfix.txt files in '+cfg.RESULTS+' -- nothing to do')
        return

    # Build rmsynth1d flags from config
    rmsyn_flags  = f'-S -v -l {cfg.RMSYN_FARADAY_RANGE} -o {cfg.RMSYN_POLY_ORDER} -s {cfg.RMSYN_RMSF_SAMPLES}'
    if cfg.RMSYN_SUPER_RESOLUTION:
        rmsyn_flags += ' --super-resolution'

    # Build rmclean1d flags from config
    rmclean_flags = f'-S -v -c {cfg.RMCLEAN_CUTOFF} -w {cfg.RMCLEAN_WINDOW}'

    print(gen.col('RM Synthesis')+'Parameters used for rmsynth1d / rmclean1d (see oxkat/config.py to change):')
    print(gen.col('  Faraday range (-l)')   +f'+/- {cfg.RMSYN_FARADAY_RANGE} rad/m^2')
    print(gen.col('  Poly order (-o)')      +f'{cfg.RMSYN_POLY_ORDER}')
    print(gen.col('  RMSF samples (-s)')    +f'{cfg.RMSYN_RMSF_SAMPLES}')
    print(gen.col('  Super-resolution')     +f'{cfg.RMSYN_SUPER_RESOLUTION}')
    print(gen.col('  Clean cutoff (-c)')    +f'{cfg.RMCLEAN_CUTOFF}')
    print(gen.col('  Clean window (-w)')    +f'{cfg.RMCLEAN_WINDOW}')
    print(gen.col('rmsynth1d flags')+rmsyn_flags)
    print(gen.col('rmclean1d flags')+rmclean_flags)
    gen.print_spacer()

    for fname in fnames:
        # UVfix files are already MAD-clipped (including V) by RMSYNTH_01
        if not fname.endswith('_UVfix.txt'):
            mad_clip_rmsynth_file(fname)
        subprocess.run([f'rmsynth1d {fname} {rmsyn_flags}'],  shell=True)
        subprocess.run([f'rmclean1d {fname} {rmclean_flags}'], shell=True)

if __name__  == "__main__":
    main()
    
