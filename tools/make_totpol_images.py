import glob, os, subprocess,sys, time
from astropy.io import fits
import numpy as np
import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)


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
    
    for image_Q in sorted(glob.glob(cfg.IMAGES + '/*-Q-*image.fits')):

        # Get the other image names
        image_U = image_Q.replace('-Q-', '-U-')
        image_V = image_Q.replace('-Q-', '-V-')
        image_P = image_Q.replace('-Q-', '-P-')

        # Initialize the P image by duplicating the Q image
        subprocess.run([f'cp {image_Q} {image_P}'], shell = True)
        msg(f'Making image: {image_P.split(cfg.IMAGES)[-1]}')

        # Run calculations and return P image
        flux_Q = get_image(image_Q)
        flux_U = get_image(image_U)
        flux_V = get_image(image_V)
        flux_P = (flux_Q ** 2 + flux_U ** 2 + flux_V ** 2) ** (0.5)
        flush_fits(flux_P, image_P)
  
if __name__ == "__main__":
    main()
