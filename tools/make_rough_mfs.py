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
def get_identifiers(prefix):
    """
    Read in the input prefix and return the image name identifiers
    for each image (or set of images).

    Arguments:
        prefix: a string of the input prefix for the wsclean imaging call
    Returns:
        identifiers: an array containing the identifiers, where each identifier
        corresponds to a set of channelized images where we will calculate
        a single beam and homogenize the resolution
    """

    
    # Get identifiers should
    msg('Getting naming identifier "groups"')
    image_arr = sorted(glob.glob(f'{prefix}-t*'))
    suffix = np.unique([im.split(f'{prefix}-t')[-1].split('-')[0] for im in image_arr])
    return [f'{prefix}-t{s}' for s in suffix]

def make_mfs(identifier):
    """
    Simplified function to make (shitty MFS file if one does not exist)
    """

    # Check if MFS image exists, if not make a pseudo MFS image
    if glob.glob(f'{identifier}-MFS-psf.fits') == []:
          
        # Figure out what the Stokes parameters are included in the images
        stokes = []
        for stoke in ['-I', '-Q', '-U', '-V']:
            if glob.glob(f'{identifier}*{stoke}-*') != []:
                stokes.append(stoke)
            if stokes == []:        
                stokes = ['']

        # Make a MFS image for each stokes parameter
        for stoke in stokes:
            msg(f'Making a (pseudo-) Stokes{stoke} MFS image; take these images with less than a grain of salt')
            z = 0
            data = []
            freq  = []
            bmaj = []
            bmin = []
            bpa = []
            images = sorted(glob.glob(f'{identifier}-[!MFS]*{stoke}-image.fits'))
            for k, im in enumerate(images[:]):
                freq.append(fits.getheader(im)['CRVAL3'])
                bmaj.append(fits.getheader(im)['BMAJ'])
                bmin.append(fits.getheader(im)['BMIN'])
                bpa.append(fits.getheader(im)['BPA'])
                msg(im)
                if z == 0:
                    header = fits.getheader(im)
                    z += 1
                data.append(get_image(im))

            # Adopt median values for each pixel and output MFS image
            data = np.median(data, axis = 0)
            header['CRVAL3'] = np.nanmean(freq)
            header['BMAJ'] = np.nanmean(bmaj)
            header['BMIN'] = np.nanmean(bmin)
            header['BPA'] = np.nanmedian(bpa)
            mfs_name = f'{identifier}-MFS{stoke}-image.homogenized.fits'
            mfs_fits = fits.PrimaryHDU(data=data, header=header)
            mfs_fits.writeto(mfs_name, overwrite=True)

def main():

    # Read in prefix return error if missing it
    if len(sys.argv) != 2:
        msg('ERROR: Missing image prefix')
        sys.exit()
    prefix = sys.argv[-1]

    # Make MFS image
    image_identifiers = get_identifiers(prefix)
    for identifier in image_identifiers[:]:    
        make_mfs(identifier)
    
if __name__ in '__main__':
    main()
