#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import glob
import logging
import numpy
import os
import random 
import scipy.signal
import shutil
import string
import sys

from astropy.io import fits
from astropy.convolution import convolve,Gaussian2DKernel
from itertools import repeat
from multiprocessing import Pool

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

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
    return bmaj,bmin,bpa,pixscale


def beam_header(fitsfile,bmaj,bmin,bpa):
        outhdu = fits.open(fitsfile,mode='update')
        outhdr = outhdu[0].header
        outhdr.set('BMAJ',bmaj,after='BUNIT')
        outhdr.set('BMIN',bmin,after='BMAJ')
        outhdr.set('BPA',bpa,after='BMIN')
        outhdr.remove('HISTORY')
        outhdu.flush()  


def convolve_fits(residual_fits,model_image,proc_id):

        proc_id = str(proc_id)

        # Set up FITS files
        #Accounts for all possible images and makes psf in any case
        psf_fits = residual_fits.replace('image','psf')
        psf_fits = residual_fits.replace('I-image','psf')
        psf_fits = residual_fits.replace('Q-image','psf')
        psf_fits = residual_fits.replace('U-image','psf')
        psf_fits = residual_fits.replace('V-image','psf')
        
        if not o.isdir(cfg.INTERVALS+'/restored'):
            os.mkdir(cfg.INTERVALS+'/restored')
        
        #restored_fits = residual_fits.replace('image','image-restored')
        restored_fits = cfg.INTERVALS+'/restored/'+residual_fits.split(cfg.INTERVALS + '/')[-1]
        shutil.copyfile(residual_fits,restored_fits)


        # Get the fitted beam
        bmaj,bmin,bpa,pixscale = get_header(psf_fits)

        logging.info('(File '+proc_id+') Residual    : '+residual_fits)
        logging.info('(File '+proc_id+') PSF         : '+psf_fits)
        logging.info('(File '+proc_id+') Restored    : '+restored_fits)
        logging.info('(File '+proc_id+') Fitted bmaj : '+str(bmaj*3600))
        logging.info('(File '+proc_id+') Fitted bmin : '+str(bmin*3600))
        logging.info('(File '+proc_id+') Fitted bpa  : '+str(bpa))

        # Create restoring beam image
        xstd = bmin/(2.3548*pixscale)
        ystd = bmaj/(2.3548*pixscale)
        theta = deg2rad(bpa)
        restoring = Gaussian2DKernel(x_stddev=xstd,y_stddev=ystd,theta=theta,x_size=cropsize,y_size=cropsize,mode='center')
        restoring_beam_image = restoring.array
        restoring_beam_image = restoring_beam_image / numpy.max(restoring_beam_image)

        # Convolve model with restoring beam
        model_conv_image = scipy.signal.fftconvolve(model_image, restoring_beam_image, mode='same')

        # Open residual image and add convolved model
        residual_image = get_image(residual_fits)
        restored_image = residual_image + model_conv_image

        # Flush restored FITS file and fix the header
        flush_fits(restored_image,restored_fits)
        beam_header(restored_fits,bmaj,bmin,bpa)


if __name__ == '__main__':  
    
        if len(sys.argv) < 5:
                print('Please specify model FITS image and a target field name and number of channels to restore')
                sys.exit()

        logfile = 'restore_model.log'
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s |  %(message)s', datefmt='%d/%m/%Y %H:%M:%S ')

        # Values of interest
        model_fits = sys.argv[1]
        targetname = sys.argv[2]
        num_channels = sys.argv[3]
        model_freq_splits = int(sys.argv[4])
        pol_list = sys.argv[5]
        
       
        
        print('Restoring model for polarisations: ' + str(pol_list))
        
        j = 1 # Number of parallel convolutions
        cropsize = 51 # Size of kernel thumbnail
        
        if int(num_channels) != 1:
            print('Overwriting input model image, looking for files createad during snapmask process')
            #SHOULD MAYBE PUT IN A CHECK HERE THAT NUM CHANNELS MATCHES THE NUMBER OF FOUDN MODEL IMAGES
            
            
            
           
            if model_freq_splits == 1:
                #Do the MFS images here
                for ij in pol_list:
                    
                    
                    # print(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits')
                    # print(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-000'+str(ei)+'-'+ij+'-model.fits')
                    try:
                        model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-MFS-'+ij+'-model.fits')[0]
                        fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-MFS-'+ij+'-image.fits')) # List of residuals
                        #fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits'))
                    except:
                        model_fits = []
                        fits_list = []
                        
                    #print(fits_list[0])
                    
                    #Checks in case -I- not in image name as only imaged I pol
                    if len(model_fits) == 0 and pol_list == 'I':
                        model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-MFS-'+'model.fits')[0]
                        fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-MFS'+'-image.fits')) # List of residuals
                    
                   
                    
                    ids = numpy.arange(0,len(fits_list))
            
                    # Get the image size from first image in list and create matched model image
                    img0 = get_image(fits_list[0])
                    nx,ny = numpy.shape(img0)
                    if nx != ny:
                            print('Only square images are supported at present')
                            sys.exit()
                    tmpmodel_fits = 'temp_model_'+''.join(random.choices(string.ascii_uppercase + string.digits, k=16))+'.fits'
                    os.system('fitstool.py -z '+str(nx)+' -o '+tmpmodel_fits+' '+model_fits)
                    model_image = get_image(tmpmodel_fits)
            
                    for i in range(0,len(fits_list)):
                        # print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                        # convolve_fits(fits_list[i],model_image,ids[i])
                        try: 
                            print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                            convolve_fits(fits_list[i],model_image,ids[i])
                        except:
                            print('Convolving Failed Skipping:' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                
                if num_channels != 1:
                    for ei in range(int(num_channels)):
                            
                        for ij in pol_list:
                            
                            
                            # print(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits')
                            # print(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-000'+str(ei)+'-'+ij+'-model.fits')
                            try:
                                model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-'+str(ei).zfill(4)+'-'+ij+'-model.fits')[0]
                                fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-'+str(ei).zfill(4)+'-'+ij+'-image.fits')) # List of residuals
                                #fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits'))
                            except:
                                model_fits = []
                                fits_list = []
                                
                            #print(fits_list[0])
                            
                            #Checks in case -I- not in image name as only imaged I pol
                            if len(model_fits) == 0 and pol_list == 'I':
                                model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-'+str(ei).zfill(4)+'-model.fits')[0]
                                fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-'+str(ei).zfill(4)+'-image.fits')) # List of residuals
                           
                            
                            ids = numpy.arange(0,len(fits_list))
                    
                            # Get the image size from first image in list and create matched model image
                            img0 = get_image(fits_list[0])
                            nx,ny = numpy.shape(img0)
                            if nx != ny:
                                    print('Only square images are supported at present')
                                    sys.exit()
                            tmpmodel_fits = 'temp_model_'+''.join(random.choices(string.ascii_uppercase + string.digits, k=16))+'.fits'
                            os.system('fitstool.py -z '+str(nx)+' -o '+tmpmodel_fits+' '+model_fits)
                            model_image = get_image(tmpmodel_fits)
                    
                            for i in range(0,len(fits_list)):
                                # print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                                # convolve_fits(fits_list[i],model_image,ids[i])
                                try: 
                                    print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                                    convolve_fits(fits_list[i],model_image,ids[i])
                                except:
                                    print('Convolving Failed Skipping:' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
             
            else:
                #Do the MFS images here, CREATING SUDO MFS IMAGE
                print('Frequency splits in model detected, creating sudo MFS image...these ARE NOT ACCURATE FLUX MEASUREMENTS')
                for ij in pol_list:
                    
                    mid_part = str(int(model_freq_splits/2))
                    # print(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits')
                    # print(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-000'+str(ei)+'-'+ij+'-model.fits')
                    try:
                        model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan'+'part'+mid_part+'-MFS-'+ij+'-model.fits')[0]
                        fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan-t*'+'-MFS-'+ij+'-image.fits')) # List of residuals
                        #fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits'))
                    except:
                        model_fits = []
                        fits_list = []
                        
                    #print(fits_list[0])
                    
                    #Checks in case -I- not in image name as only imaged I pol
                    if len(model_fits) == 0 and pol_list == 'I':
                        model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chanpart'+mid_part+'-MFS-'+'model.fits')[0]
                        fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan-t*'+'-MFS'+'-image.fits')) # List of residuals
                   
                    
                    ids = numpy.arange(0,len(fits_list))
            
                    # Get the image size from first image in list and create matched model image
                    img0 = get_image(fits_list[0])
                    nx,ny = numpy.shape(img0)
                    if nx != ny:
                            print('Only square images are supported at present')
                            sys.exit()
                    tmpmodel_fits = 'temp_model_'+''.join(random.choices(string.ascii_uppercase + string.digits, k=16))+'.fits'
                    os.system('fitstool.py -z '+str(nx)+' -o '+tmpmodel_fits+' '+model_fits)
                    model_image = get_image(tmpmodel_fits)
            
                    for i in range(0,len(fits_list)):
                        # print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                        # convolve_fits(fits_list[i],model_image,ids[i])
                        try: 
                            print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                            convolve_fits(fits_list[i],model_image,ids[i])
                        except:
                            print('Convolving Failed Skipping:' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                
                
                for nm in range(int(model_freq_splits)):
                    for ei in range(int(int(num_channels)/model_freq_splits)):
                        for ij in pol_list:
                            
                            real_chan_number = int(ei +(nm * (int(num_channels)/model_freq_splits)))
                            
                            #print('Looking for images like: '+cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-'+str(real_chan_number).zfill(4)+'-'+ij+'-image.fits')
                        
                            # print(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits')
                            # print(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-000'+str(ei)+'-'+ij+'-model.fits')
                            try:
                                model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chanpart'+str(nm)+'-'+str(ei).zfill(4)+'-'+ij+'-model.fits')[0]
                                fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-'+str(real_chan_number).zfill(4)+'-'+ij+'-image.fits')) # List of residuals
                                #fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+'-t*'+'-000'+str(ei)+'-'+ij+'-image.fits'))
                            except:
                                model_fits = []
                                fits_list = []
                                
                            print(model_fits)
                            
                            
                                
                            #print(fits_list[0])
                            
                            #Checks in case -I- not in image name as only imaged I pol
                            if len(model_fits) == 0 and pol_list == 'I':
                                model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chanpart'+str(nm)+'-'+str(ei).zfill(4)+'-model.fits')[0]
                                fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-'+str(real_chan_number).zfill(4)+'-image.fits')) # List of residuals
                           
                            
                            ids = numpy.arange(0,len(fits_list))
                    
                            # Get the image size from first image in list and create matched model image
                            img0 = get_image(fits_list[0])
                            nx,ny = numpy.shape(img0)
                            if nx != ny:
                                    print('Only square images are supported at present')
                                    sys.exit()
                            tmpmodel_fits = 'temp_model_'+''.join(random.choices(string.ascii_uppercase + string.digits, k=16))+'.fits'
                            os.system('fitstool.py -z '+str(nx)+' -o '+tmpmodel_fits+' '+model_fits)
                            model_image = get_image(tmpmodel_fits)
                    
                            for i in range(0,len(fits_list)):
                                # print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                                # convolve_fits(fits_list[i],model_image,ids[i])
                                try: 
                                    print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                                    convolve_fits(fits_list[i],model_image,ids[i])
                                except:
                                    print('Convolving Failed Skipping:' + fits_list[i].split(cfg.INTERVALS + '/')[-1])

    
#######################This is now 1 channel only, no need to use freqeuncy splits
            
        else:
            
            for ij in pol_list:
                
                try:
                    model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-'+ij+'-model.fits')[0]
                    fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-'+ij+'-image.fits')) # List of residuals
                except:
                    model_fits = []
                
                #Checks in case -I- not in image name as only imaged I pol
                if len(model_fits) == 0 and pol_list == 'I':
                    model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask'+str(num_channels)+'chan-'+'model.fits')[0]
                    fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*modelsub'+str(num_channels)+'chan'+'-t*'+'-image.fits')) # List of residuals
                    
                    # model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask*-model.fits')[0]
                    # print(cfg.INTERVALS + f'/*{targetname}*-image.fits')
                    # fits_list = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*-image.fits')) # List of residuals
                
                ids = numpy.arange(0,len(fits_list))
        
                # Get the image size from first image in list and create matched model image
                img0 = get_image(fits_list[0])
                nx,ny = numpy.shape(img0)
                if nx != ny:
                        print('Only square images are supported at present')
                        sys.exit()
                tmpmodel_fits = 'temp_model_'+''.join(random.choices(string.ascii_uppercase + string.digits, k=16))+'.fits'
                os.system('fitstool.py -z '+str(nx)+' -o '+tmpmodel_fits+' '+model_fits)
                model_image = get_image(tmpmodel_fits)
        
                for i in range(0,len(fits_list)):
                    try: 
                        print('Convolving: ' + fits_list[i].split(cfg.INTERVALS + '/')[-1])
                        convolve_fits(fits_list[i],model_image,ids[i])
                    except:
                        print('Convolving Failed Skipping:' + fits_list[i].split(cfg.INTERVALS + '/')[-1])

        # pool = Pool(processes=j)
        # pool.starmap(convolve_fits,zip(fits_list,repeat(model_image),ids))

