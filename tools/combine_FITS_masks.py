#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk

# Usage: combine_FITS_masks.py [--border N] <output.fits> <mask1.fits> <mask2.fits> [...]
# Writes the union (logical OR) of the input masks to output.fits, using the
# header of the first input mask. All inputs must share the same pixel grid.
# --border N zeroes the outermost N pixels along every image edge of the
# combined mask, so nothing near the edge (where the image is least reliable)
# is cleaned.


from astropy.io import fits
import numpy
import sys


def combine_masks(opfits, maskfits, border=0):
    with fits.open(maskfits[0]) as hdul:
        header = hdul[0].header
        combined = hdul[0].data != 0
    for mask in maskfits[1:]:
        data = fits.getdata(mask)
        if data.shape != combined.shape:
            sys.exit(f'Mask shape mismatch: {mask} {data.shape} vs {maskfits[0]} {combined.shape}')
        combined |= data != 0
    if border > 0:
        # The last two axes are the image (y, x); any leading axes are
        # degenerate Stokes/frequency axes
        combined[..., :border, :] = False
        combined[..., -border:, :] = False
        combined[..., :, :border] = False
        combined[..., :, -border:] = False
    fits.writeto(opfits, combined.astype(numpy.float32), header, overwrite=True)


def main():

    args = sys.argv[1:]
    border = 0
    if args[:1] == ['--border']:
        if len(args) < 2:
            sys.exit('--border needs a number of pixels')
        border = int(args[1])
        args = args[2:]
    if len(args) < 3:
        sys.exit('Usage: combine_FITS_masks.py [--border N] <output.fits> <mask1.fits> <mask2.fits> [...]')
    opfits = args[0]
    maskfits = args[1:]
    combine_masks(opfits, maskfits, border)


if __name__ == "__main__":


    main()
