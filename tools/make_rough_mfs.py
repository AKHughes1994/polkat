#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import os
import sys
import glob
import subprocess
import time
import scipy
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from scipy.spatial import ConvexHull
import bottleneck as bn

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

# Helper functions
def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)

# Fits manipulation functions
def flush_fits(newimage,newheader, fitsfile):
    f = fits.open(fitsfile,mode='update')
    input_hdu = f[0]
    input_hdu.header = newheader
    if len(input_hdu.data.shape) == 2:
            input_hdu.data[:,:] = newimage
    elif len(input_hdu.data.shape) == 3:
            input_hdu.data[0,:,:] = newimage
    else:
            input_hdu.data[0,0,:,:] = newimage
    f.flush()

def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
            image = np.array(input_hdu.data[:,:])
    elif len(input_hdu.data.shape) == 3:
            image = np.array(input_hdu.data[0,:,:])
    else:
            image = np.array(input_hdu.data[0,0,:,:])
    return image


# Functions to fix names and actually run spatial homogenization
def make_mfs(prefix):
    """
    Combine all channels for each Stokes parameter into a single MFS image
    
    Arguments:
        prefix: a string of the input prefix (e.g., img_1697283077_sdp_l0_1024ch_SwiftJ1727.ms_pcalmask_zoom)
    """

    # Process all Stokes parameters
    stokes = ['-I', '-Q', '-U', '-V', '-Plin', '-Ptot']

    # Make a MFS image for each stokes parameter by combining ALL channels
    for stoke in stokes:
        # Find all channel images for this Stokes parameter, excluding MFS images
        images = sorted([im for im in glob.glob(f'{prefix}-*{stoke}-image.fits') if 'MFS' not in im])
        
        if images == []:
            msg(f'No images found for Stokes {stoke}, skipping')
            continue
            
        msg(f'Making MFS for Stokes{stoke} by combining {len(images)} channels')
        z = 0
        data = []
        freq  = []
        bmaj = []
        bmin = []
        bpa = []
        
        for k, im in enumerate(images[:]):
            if fits.getheader(im)['BMAJ'] > 1e-14:
                freq.append(fits.getheader(im)['CRVAL3'])
                bmaj.append(fits.getheader(im)['BMAJ'])
                bmin.append(fits.getheader(im)['BMIN'])
                bpa.append(fits.getheader(im)['BPA'])
                if z == 0:
                    header = fits.getheader(im)
                    z += 1
                data.append(get_image(im))

        # Skip if no valid data found
        if len(data) == 0:
            msg(f'No valid data for Stokes {stoke}, skipping')
            continue
            
        # Adopt median values for each pixel and output MFS image
        data = bn.nanmedian(data, axis = 0)
        for _ in range(4 - len(data.shape)): # Make (1, 1, N, N)-shaped
            data = data[None, :]

        header['CRVAL3'] = np.nanmean(freq)
        header['BMAJ'] = np.nanmean(bmaj)
        header['BMIN'] = np.nanmean(bmin)
        header['BPA'] = np.nanmedian(bpa)
        mfs_name = f'{prefix}-MFS{stoke}-image.fits'
        mfs_fits = fits.PrimaryHDU(data=data, header=header)
        mfs_fits.writeto(mfs_name, overwrite=True)
        msg(f'Wrote {mfs_name}')

def main():

    # Read in prefix return error if missing it
    if len(sys.argv) != 2:
        msg('ERROR: Missing image prefix')
        sys.exit()
    prefix = sys.argv[-1]

    # Make MFS image by combining all channels for each Stokes parameter
    make_mfs(prefix)
    
if __name__ in '__main__':
    main()
