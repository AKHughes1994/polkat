import glob,os,datetime, subprocess, sys, json
import numpy as np
import os.path as o
import time
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

def main():

    # Build rmsynth1d flags from config
    rmsyn_flags  = f'-S -v -l {cfg.RMSYN_FARADAY_RANGE} -o {cfg.RMSYN_POLY_ORDER} -s {cfg.RMSYN_RMSF_SAMPLES}'
    if cfg.RMSYN_SUPER_RESOLUTION:
        rmsyn_flags += ' --super-resolution'

    # Build rmclean1d flags from config
    rmclean_flags = f'-S -v -c {cfg.RMCLEAN_CUTOFF} -w {cfg.RMCLEAN_WINDOW}'

    # Iterate through all rmsynth.txt files in the RESULTS directory
    fnames = glob.glob(cfg.RESULTS + '/*_rmsynth.txt')

    for fname in fnames:
        subprocess.run([f'rmsynth1d {fname} {rmsyn_flags}'],  shell=True)
        subprocess.run([f'rmclean1d {fname} {rmclean_flags}'], shell=True)

if __name__  == "__main__":
    main()
    
