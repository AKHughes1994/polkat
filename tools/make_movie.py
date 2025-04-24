#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import logging
import os
import random
import numpy
import string
import sys

from astropy.io import fits
from astropy.time import Time
from multiprocessing import Pool
from PIL import Image, ImageDraw, ImageFont

# Path to the font used for overlaying text on images
fontPath = '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'
sans30 = ImageFont.truetype(fontPath, 30)

def generate_temp(k=16):
    """
    Generate a temporary filename for intermediate FITS files.

    Parameters:
    - k (int): Length of the random string appended to the filename.

    Returns:
    - str: Temporary FITS filename.
    """
    tmpfits = 'temp_' + ''.join(random.choices(string.ascii_uppercase + string.digits, k=k)) + '.fits'
    return tmpfits

def make_png(ff, i):
    """
    Process a FITS file to generate a PNG image with metadata overlays.

    Parameters:
    - ff (str): Path to the input FITS file.
    - i (int): Frame index for the current FITS file.

    Returns:
    - None
    """
    # Generate a temporary FITS file for processing
    tmpfits = generate_temp()

    logging.info(f' | File {i} | Input image {ff}')
    logging.info(f' | File {i} | Temp image  {tmpfits}')

    # Downsample the FITS file using mShrink
    os.system(f'mShrink {ff} {tmpfits} 2')

    # Extract metadata from the FITS header
    input_hdu = fits.open(ff)[0]
    hdr = input_hdu.header
    map_date = hdr.get('DATE-OBS')
    t_mjd = Time(map_date, format='isot', scale='utc').mjd
    tt = f'{map_date} | {t_mjd}'

    # Generate the output PNG filename
    pp = f'pic_{str(i).zfill(4)}.png'
    logging.info(f' | File {i} | PNG         {pp}')

    # Generate the PNG image using mViewer
    syscall = f'mViewer -ct 0 -gray {tmpfits} -0.0008 0.0018 -out {pp}'
    os.system(syscall)

    logging.info(f' | File {i} | Time        {tt}')

    # Open the generated PNG and overlay metadata
    img = Image.open(pp)
    xx, yy = img.size
    draw = ImageDraw.Draw(img)
    draw.text((0.03 * xx, 0.90 * yy), f'Frame : {str(i).zfill(len(str(nframes)))} / {nframes}', fill='white', font=sans30)
    draw.text((0.03 * xx, 0.93 * yy), f'Time  : {tt}', fill='white', font=sans30)
    draw.text((0.03 * xx, 0.96 * yy), f'Image : {ff}', fill='white', font=sans30)
    img.save(pp)

    # Remove the temporary FITS file
    os.system(f'rm {tmpfits}')
    logging.info(f' | File {i} | Done')

if __name__ == '__main__':
    """
    Main script execution:
    - Processes FITS files in the specified directory.
    - Generates PNG images with metadata overlays.
    - Compiles the PNG images into a video using FFmpeg.
    """
    # Get the input directory from the command-line arguments
    interval = sys.argv[-1]

    # Find unique prefixes for FITS files in the directory
    suffixes = glob.glob(f'{interval}/*_restored-*')
    suffixes = numpy.unique([x.split('_restored')[0] for x in suffixes])

    # Set up logging
    logfile = 'make_movie.log'
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s |  %(message)s', datefmt='%d/%m/%Y %H:%M:%S ')

    # Process each suffix
    for suffix in suffixes:

        for stoke in ['', '-I', '-Q', '-U', '-V', '-Ptot', '-Plin']:
            # Find FITS files for the current suffix and stoke
            fitslist = sorted(glob.glob(f'{suffix}_restored-t*-MFS{stoke}-image.fits'))
            if not fitslist:
                fitslist = sorted(glob.glob(f'{suffix}_restored-t*-MFS{stoke}-image.homogenized.fits'))

            if fitslist:
                ids = numpy.arange(0, len(fitslist))
                nframes = len(fitslist)
                j = 8  # Number of parallel processes

                # Process FITS files in parallel to generate PNGs
                pool = Pool(processes=j)
                pool.starmap(make_png, zip(fitslist, ids))

                # Compile the PNGs into a video using FFmpeg
                frame = '2340x2340'
                fps = 10
                opmovie = fitslist[0].split('-t')[0] + f'{stoke}.mp4'
                os.system(f'ffmpeg -y -r {fps} -f image2 -s {frame} -i pic_%04d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p {opmovie}')
