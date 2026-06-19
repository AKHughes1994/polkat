import glob,os,datetime, subprocess, sys, json
import numpy as np
import os.path as o
import time
from functools import partial
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg
from oxkat import generate_jobs as gen

# Reinitialize print to flush in real-time for SLURM logging
print = partial(print, flush=True)

def main():

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

    # Iterate through all rmsynth.txt files in the RESULTS directory
    fnames = glob.glob(cfg.RESULTS + '/*_rmsynth.txt')

    for fname in fnames:
        subprocess.run([f'rmsynth1d {fname} {rmsyn_flags}'],  shell=True)
        subprocess.run([f'rmclean1d {fname} {rmclean_flags}'], shell=True)

if __name__  == "__main__":
    main()
    
