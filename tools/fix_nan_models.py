#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
Fix NaN models and optionally zero inner regions of model FITS files.

This script processes wsclean model FITS files to:
1. Replace any NaN-valued models with zeros (prevents prediction issues)
2. Optionally zero out an inner box region of specified size (useful for removing central source)

Usage:
    python fix_nan_models.py <pattern> [inner_box_size]
    
Arguments:
    pattern: File pattern prefix for model files (e.g., 'img_myms_datamask')
    inner_box_size: (optional) Integer size of inner box to zero out (in pixels)

Example:
    python fix_nan_models.py img_myms_datamask
    python fix_nan_models.py img_myms_datamask 512
"""

import sys
import glob
import numpy
from astropy.io import fits


def get_image(fitsfile):
    """
    Extract 2D image data from a FITS file.
    
    Handles FITS files with 2, 3, or 4 dimensions by extracting the 2D spatial slice.
    
    Parameters:
        fitsfile (str): Path to the FITS file
        
    Returns:
        numpy.ndarray: 2D image array
    """
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
            image = numpy.array(input_hdu.data[:,:])
    elif len(input_hdu.data.shape) == 3:
            image = numpy.array(input_hdu.data[0,:,:])
    else:
            image = numpy.array(input_hdu.data[0,0,:,:])
    return image


def flush_fits(newimage, fitsfile):
    """
    Write modified 2D image data back to a FITS file.
    
    Updates the FITS file in place, handling 2, 3, or 4 dimensional data structures.
    
    Parameters:
        newimage (numpy.ndarray): 2D image array to write
        fitsfile (str): Path to the FITS file to update
    """
    f = fits.open(fitsfile, mode='update')
    input_hdu = f[0]
    if len(input_hdu.data.shape) == 2:
            input_hdu.data[:,:] = newimage
    elif len(input_hdu.data.shape) == 3:
            input_hdu.data[0,:,:] = newimage
    else:
            input_hdu.data[0,0,:,:] = newimage
    f.flush()


def zero_inner_box(image, box_size):
    """
    Zero out a centered square region in an image.
    
    Parameters:
        image (numpy.ndarray): 2D image array
        box_size (int): Size of the square box to zero (in pixels)
        
    Returns:
        numpy.ndarray: Modified image with zeroed inner box
    """
    ny, nx = image.shape
    center_y, center_x = ny // 2, nx // 2
    half_box = box_size // 2
    
    y_min = max(0, center_y - half_box)
    y_max = min(ny, center_y + half_box)
    x_min = max(0, center_x - half_box)
    x_max = min(nx, center_x + half_box)
    
    image[y_min:y_max, x_min:x_max] = 0.0
    return image


# Check command line arguments
if len(sys.argv) < 2 or len(sys.argv) > 3:
    print('ERROR: Invalid number of arguments', flush=True)
    print('Usage: python fix_nan_models.py <pattern> [inner_box_size]', flush=True)
    print('  pattern: File pattern prefix for model files', flush=True)
    print('  inner_box_size: (optional) Integer size of inner box to zero (in pixels)', flush=True)
    sys.exit(1)

# Parse command line arguments
pattern = sys.argv[1]
inner_box_size = None
if len(sys.argv) > 2:
    inner_box_size = int(sys.argv[2])
    print(f'Will zero inner {inner_box_size}x{inner_box_size} pixel box in all models', flush=True)

# Process all model FITS files matching the pattern
fitslist = sorted(glob.glob(pattern+'*model.fits'))
for fitsfile in fitslist:
        img = get_image(fitsfile)
        maxval = numpy.max(img)
        modified = False
        
        # Check for NaN models and replace with zeros
        if numpy.isnan(maxval):
                img = numpy.zeros((img.shape[0], img.shape[1]))
                print(fitsfile, maxval, 'zeroing NaN model', flush=True)
                modified = True
        
        # Zero inner box if requested
        if inner_box_size is not None:
                img = zero_inner_box(img, inner_box_size)
                print(fitsfile, maxval, f'zeroed inner {inner_box_size}x{inner_box_size} box', flush=True)
                modified = True
        
        # Write back to FITS file if modifications were made
        if modified:
                flush_fits(img, fitsfile)
        else:
                print(fitsfile, maxval, flush=True)