import glob, os, subprocess,sys, time
from astropy.io import fits
import numpy as np
import os.path as o
import sys

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush = True)


def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
            image = np.array(input_hdu.data[:,:])
    elif len(input_hdu.data.shape) == 3:
            image = np.array(input_hdu.data[0,:,:])
    else:
            image = np.array(input_hdu.data[0,0,:,:])
    return image


def flush_fits(newimage,fitsfile):
    f = fits.open(fitsfile,mode='update')
    input_hdu = f[0]
    if len(input_hdu.data.shape) == 2:
            input_hdu.data[:,:] = newimage
    elif len(input_hdu.data.shape) == 3:
            input_hdu.data[0,:,:] = newimage
    else:
            input_hdu.data[0,0,:,:] = newimage
    f.flush()

def main():

    # Read in directory
    if len(sys.argv) != 2:
        sys.exit('ERROR: Please include directory to parse for images (e.g., IMAGES or INTERVALS)')
    
    directory = sys.argv[-1]
    
    for image_Q in sorted(glob.glob(directory + '/*-Q-*image.fits') + glob.glob(directory + '/*-Q-*image.homogenized.fits')):

        # Get the other image names
        image_U = image_Q.replace('-Q-', '-U-')
        image_V = image_Q.replace('-Q-', '-V-')
        image_Plin = image_Q.replace('-Q-', '-Plin-')
        image_Ptot = image_Q.replace('-Q-', '-Ptot-')

        # Skip if both Plin and Ptot images already exist
        if os.path.exists(image_Plin) and os.path.exists(image_Ptot):
            msg(f"Skipping {image_Q}: Plin and Ptot images already exist.")
            continue

        # Initialize the P image by duplicating the Q image
        subprocess.run([f'cp {image_Q} {image_Plin}'], shell = True)
        subprocess.run([f'cp {image_Q} {image_Ptot}'], shell = True)
        msg(f'Making image: {image_Plin.split(directory + "/")[-1]}')

        # Run calculations and return P images
        flux_Q = get_image(image_Q)
        flux_U = get_image(image_U)
        flux_Plin = (flux_Q ** 2 + flux_U ** 2 ) ** (0.5)
        flush_fits(flux_Plin, image_Plin)

        # Incase Stokes V has been deleted 
        try:
            flux_V = get_image(image_V)
            flux_Ptot = (flux_Q ** 2 + flux_U ** 2 + flux_V ** 2) ** (0.5)
            flush_fits(flux_Ptot, image_Ptot)
  
        except:
            pass

if __name__ == "__main__":
    main()
