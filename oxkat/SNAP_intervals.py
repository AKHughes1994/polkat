#!/usr/bin/python
# andrew.hughes@physics.ox.ac.uk

import glob
import os
import pickle
import time
import sys
import os.path as o
import subprocess
import numpy as np

from pyrap.tables import table

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)


def get_total_ints(myms):

    tt = table(myms,readonly=True)

    scan_numbers = list(set(tt.getcol('SCAN_NUMBER')))
    exposure = round(np.mean(tt.getcol('EXPOSURE')),4)

    int_arr = []

    for scan in scan_numbers:
        subtab = tt.query(query='SCAN_NUMBER=='+str(scan))
        times = subtab.getcol('TIME')
        n_int = round((times[-1] - times[0]) / exposure, 0)
        int_arr.append(int(n_int))

    return int_arr

def main():
    
    intbin    = cfg.SNAP_INTBIN
    chanout = cfg.SNAP_CHANNELSOUT
    imsize   = cfg.SNAP_IMSIZE

    # Stokes choices
    pol = 'I'
    if cfg.SNAP_POL:
        pol = 'IQUV'

    # Deconvolution choice -- don't do this unless absolutely necessary
    niter = 0
    if cfg.SNAP_DECONV:
        niter = 100
              
    mask = False
    if cfg.SNAP_DECONVMASK:
        mask = cfg.SNAP_DECONVMASK  

    #img_size = cfg.SNAP_SIZE
    #snap_mask = cfg.SNAP_MASK_PATH
    
    if len(sys.argv) == 1:
        print('Please specify an MS file for interval imaging')
        sys.exit()
    else:
        myms = sys.argv[1]
    
    # Get number of integrations in each scan
    int_arr = get_total_ints(myms)
    nint_arr = []

    # Iterate through scans:
    scan_n = 0
    for ints in int_arr:

        # Starting interval        
        int0 = int(np.sum(int_arr[:scan_n]))
        
        # Determine if we want to truncate the end intervals or not
        if cfg.SNAP_INTEND:
            int1 = ints - ints % intbin + int0
            nint = int((int1 - int0) / intbin)
        else:
            int1 = int0 + ints
            nint = int(np.ceil((int1 - int0) / intbin))

        nint_arr.append(nint)

        # Temporary name for scan image
        image_prefix = cfg.INTERVALS+f'/img_{myms}_modelsub_scan{scan_n:04d}'

        # Generate imaging call and run
        imcall = gen.generate_syscall_wsclean(mslist = [myms],
                        imgname = image_prefix,
                        datacol = 'DATA',
                        chanout = chanout,
                        imsize = imsize,
                        niter = niter,
                        nomodel = True,
                        nodirty = True,
                        makepsf = True,
                        pol = pol,
                        intervalsout = nint,
                        interval0 = int0,
                        interval1 = int1,
                        field='0',
                        mask = mask)

        for syscall in imcall:
            subprocess.run([syscall], shell=True)

        # Fix channel naming (necessary when breaking up frequency images due to memory constraints)
        syscall = f'python3 {cfg.TOOLS}/fix_image_naming.py {chanout} {image_prefix}'
        subprocess.run([syscall], shell=True)

        # Fix 'scan' naming to obey the 't0...' wsclean naming convention -- again has to be a bettter way to do this but this'll work
        # Basic idea is to add the scan number to the t-label; e.g., if scan 0 has 100 integrations:
        # scan0001-t0000 becomes t100 (easier indexing imo and avoids scan boundaries without the need to split MS files)
        images= sorted(glob.glob(f'{image_prefix}*'))
        
        msg('Removing "scan" label in image names to standardize the snapshots')
        for image in images:
            image_prefix_fix = image_prefix.replace(f'_scan{scan_n:04d}','') # remove scan string
            image_fix = image.replace(image_prefix, image_prefix_fix) # replace it in image name
            suffix = image_fix.split(image_prefix_fix + '-t')[-1] # get suffix with improper t-label number
            suffix_fix = suffix.split('-') 
            suffix_fix[0] = '{:04d}'.format(int(suffix_fix[0]) + (nint_arr[scan_n - 1] if scan_n > 0 else 0)) # add past ints to t-label
            image_fix = image_fix.replace(suffix, '-'.join(suffix_fix)) # replace t-label in image
            syscall = f'mv {image} {image_fix}'
            subprocess.run([syscall], shell = True)

        scan_n += 1
        
if __name__ == "__main__":
    main()
