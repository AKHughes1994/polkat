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

        # Get the array of modelsub images and restoring model iamges
        modelsub_images = sorted(glob.glob(f'{target_prefix}*image.fits'))
        if mfs_only:
            modelsub_images = [im for im in modelsub_images if '-MFS-' in im]
        model_images = sorted(glob.glob(f'{model_prefix}*-model.fits'))

        # Check to see if the model has different dimension than the snapshot imaging
        truncate_model = False
        model_shape = get_image(model_images[0]).shape
        modelsub_shape = get_image(modelsub_images[0]).shape
    
        # Check that images are square
        if model_shape[0] != model_shape[1] or modelsub_shape[0] != modelsub_shape[1]: 
            sys.exit('ERROR: Both model and snapshot images must be square')

        if model_shape != modelsub_shape:
            truncate_model = True

        # Match the modelsub image with its static counterpart 
        for modelsub_image in modelsub_images[:]:

            # A bunch of string manipulation to get the relevant images
            suffix = '-'.join(modelsub_image.split(target_prefix + '-t')[-1].split('-')[1:]).replace('image.fits', 'model.fits') # this will extract frequency/stokes suffix for model matching
            psf_image = modelsub_image[:]
            for substring in ['-I-image.fits','-Q-image.fits','-U-image.fits','-V-image.fits']:
                psf_image = psf_image.replace(substring, '-image.fits')
            psf_image = psf_image.replace('-image.fits','-psf.fits')
            psf = get_header(psf_image)

            # Get the model image
            model_image = [im for im in model_images if im.endswith(suffix)]
            if model_image != []:
                model_image = model_image[0]
            else:
                # This is necessary if you did IQUV imaging for pcalmask, but only I imaging for snapshots
                model_image = [im for im in model_images if im.endswith(suffix.replace('-model.fits', '-I-model.fits'))][0]

            # Finally -- Restore images

            # Check channel is flagged (i.e., BEAM information is 0.0)
            bmaj = psf[0]
            if bmaj == None:
                bmaj = 0.0

            if bmaj > 1.0e-14:
                
                # Copy psf naming for restored image -- necessary for beam homogenizing routine
                subprocess.run([f'cp {psf_image} {psf_image.replace("modelsub", "restored")}'], shell = True)

                msg('Restoring image: ' + modelsub_image.split(cfg.INTERVALS + '/')[-1])
                msg('Restoring model: ' + model_image)
            
                # Get model data and truncate as necessary
                if truncate_model:
                    temp_model = model_image.replace('.fits', '.tmp.fits')
                    subprocess.run([f'fitstool.py -f -z {modelsub_shape[0]} -o {temp_model} {model_image}'], shell = True)
                    model_data = get_image(temp_model)
                    subprocess.run([f'rm -rf {temp_model}'], shell = True)
                else:
                    model_data = get_image(model_image)

                # Restore static model
                restore_fits(modelsub_image, model_data, psf)

            else:
                msg('Cannot Restore Image: ' + modelsub_image.split(cfg.INTERVALS + '/')[-1] + '; Channel Likely Flagged\n')            
                
