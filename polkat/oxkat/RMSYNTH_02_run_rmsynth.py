import glob,os,datetime, subprocess, sys, json
import numpy as np
import os.path as o
import time
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

def main():

    # Iterate through all rmsynth.txt files in the RESULTS directory
    fnames = glob.glob(cfg.RESULTS + '/*_rmsynth.txt')

    for fname in fnames:
        syscall = f'rmsynth1d {fname} -S -v -o 4'
        subprocess.run([syscall], shell=True)

        syscall = f'rmclean1d {fname} -S -v'
        subprocess.run([syscall], shell=True)

    # Iterate through all rmsynth.txt files in the RESULTS directory
    fnames = glob.glob(cfg.RESULTS + '/*_rmsynth_sys.txt')

    for fname in fnames:
        syscall = f'rmsynth1d {fname} -S -v -o 4'
        subprocess.run([syscall], shell=True)

        syscall = f'rmclean1d {fname} -S -v'
        subprocess.run([syscall], shell=True)

if __name__  == "__main__":
    main()
    
