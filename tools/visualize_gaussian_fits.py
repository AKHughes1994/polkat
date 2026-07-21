import json
import glob
import subprocess
import numpy as np
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, Ellipse
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u


def make_image(image, coords, rms, ref = None):
    '''
    Make zoomed in PNG image near the source of interest
    and overplot the fitted positions with beam shapes
    Inputs:
        image  = string containing the name of the image
        coords = a single of an array of SkyCoord object coordinates
        rms    = rms noise 
    Outputs:
        None (Saves image to local directory)
    '''

    if ref is None:
        ref = coords

    # Load in the image
    hdu = fits.open(image)[0]
    
    # Get the beam parameters
    bmaj = hdu.header['BMAJ']
    bmin = hdu.header['BMIN']
    bpa  = hdu.header['BPA']
    pixel_size = abs(hdu.header['CDELT2'])

    # Correct the header file
    hdu.header['NAXIS'] = 2
    hdu.data = hdu.data[0,0,:,:]
    del hdu.header['NAXIS3']
    del hdu.header['NAXIS4']
    del hdu.header['CTYPE3']
    del hdu.header['CRPIX3']
    del hdu.header['CRVAL3']
    del hdu.header['CDELT3']
    del hdu.header['CUNIT3']
    del hdu.header['CTYPE4']
    del hdu.header['CRPIX4']
    del hdu.header['CRVAL4']
    del hdu.header['CDELT4']
    del hdu.header['CUNIT4']

    # Get the WCS header
    wcs = WCS(hdu.header)
    xpos, ypos = wcs.wcs_world2pix(ref.ra[0], ref.dec[0], 1)
  
    # Make the cutout
    cutout = Cutout2D(hdu.data * 1.0e3, position=(xpos,ypos), size=(80,80), wcs=wcs)
    
    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data
    
    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())
    wcs = WCS(hdu.header)

    # Set mpl Defaults
    majorw = 5.0; majorl = 20 #major tick width and length (pixels)

#print(mpl.rcParams.keys())

    mpl.rcParams['font.size'] = 20
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['figure.facecolor'] = 'white'
    mpl.rcParams['lines.color'] = 'black'
    mpl.rcParams['lines.markeredgecolor'] = 'black'
    mpl.rcParams["axes.edgecolor"] = "black"
    mpl.rcParams["axes.linewidth"] = 5.0
	
    # Initialize plot
    fig = plt.figure(figsize=(15,15))
    fig.set_facecolor('white')
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs) 

    # Make the image
    im = ax.imshow(hdu.data)
    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
    ax.tick_params(axis='both', which='major', direction='in', length=majorl, width=majorw, color='black')
    ax.contour(hdu.data, levels=[3.5 * rms,5 * rms, 10 * rms, 20 * rms, 40 * rms], colors='white', linewidths=2.5)
    #ax.grid('major', ls=':', lw=5.0, color='k', alpha=0.5)
    cbar = plt.colorbar(im, label="Pixel Flux Density (mJy / beam)", shrink=0.8)

    # Plot the best fit positions of the components
    for coord in coords:
        x, y = wcs.wcs_world2pix(coord.ra, coord.dec, 0)
        circ = Circle((x,y), 0.3, color='k')
        ax.add_patch(circ)
        ellipse = Ellipse((x,y),bmaj/pixel_size, bmin/pixel_size, angle=(bpa - 90.0), lw=5, ls='-', color = None, facecolor = None, edgecolor='red', fill = False, zorder=1000)
        ax.add_patch(ellipse)
    

    # Include the beam shpe in the image
    ellipse = Ellipse((5,5), bmaj/pixel_size, bmin/pixel_size, angle=(bpa - 90.0), \
                        lw=5, ls='-', color ='white', zorder=1000)
    ax.add_patch(ellipse)

    plotfile = image.split('.fits')[0].split('IMAGES/')[-1] + '_zoom.png'
    plt.savefig(plotfile)
    plt.clf()
    plt.close()

def main():

    # Load in the Arguments
    parser = ArgumentParser(description='TROLOLOLO')
    parser.add_argument('-i', '--image', dest='image',
                        help='Path to IQUV Images')
    parser.add_argument('-j', '--json-file', dest='jfile',
                        help="Path to .json dictionary fit information")

    args = parser.parse_args()
    image = str(args.image)
    jfile = str(args.jfile)

    # Get the I and P image paths
    image_I = glob.glob(image + '*I-image.fits')[0]
    image_P = glob.glob(image + '*P-image.fits')[0]

    # Open up the input .json dictionary and get the Stokes I and Stokes P positions to feed into visualization function
    with open(jfile, 'r') as f:
        json_dict = json.load(f)


    ra_I  = []
    dec_I = []
    ra_P  = []
    dec_P = []
    rms_I = []
    rms_P = []

    for key in [keys for keys in json_dict['MFS'].keys() if 'component' in keys]:
        ra_I.append(json_dict['MFS'][key]['I_RA_deg'])
        dec_I.append(json_dict['MFS'][key]['I_Dec_deg'])
        ra_P.append(json_dict['MFS'][key]['P_RA_deg'])
        dec_P.append(json_dict['MFS'][key]['P_Dec_deg'])
        rms_I.append(json_dict['MFS'][key]['I_rms_mJy'])  
        rms_P.append(json_dict['MFS'][key]['P_rms_mJy'])

    rms_I = np.mean(rms_I)
    rms_P = np.mean(rms_P)
    coords_I = SkyCoord(ra = ra_I * u.deg, dec = dec_I * u.deg)
    coords_P = SkyCoord(ra = ra_P * u.deg, dec = dec_P * u.deg)
           

    # Make (zoomed) image for Stokes I
    make_image(image_I, coords_I, rms_I, ref = coords_I) 
    make_image(image_P, coords_P, rms_P, ref = coords_I) 

if __name__ in "__main__":
    main()

