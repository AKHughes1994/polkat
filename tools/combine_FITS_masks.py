#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk

# Usage: combine_FITS_masks.py <output.fits> <mask1.fits> <mask2.fits> [...]
# Writes the union (logical OR) of the input masks to output.fits, using the
# header of the first input mask. All inputs must share the same pixel grid.


from astropy.io import fits
import numpy
import sys


def combine_masks(opfits, maskfits):
    with fits.open(maskfits[0]) as hdul:
        header = hdul[0].header
        combined = hdul[0].data != 0
    for mask in maskfits[1:]:
        data = fits.getdata(mask)
        if data.shape != combined.shape:
            sys.exit(f'Mask shape mismatch: {mask} {data.shape} vs {maskfits[0]} {combined.shape}')
        combined |= data != 0
    fits.writeto(opfits, combined.astype(numpy.float32), header, overwrite=True)


def main():

    opfits = sys.argv[1]
    maskfits = sys.argv[2:]
    if len(maskfits) < 2:
        sys.exit('Usage: combine_FITS_masks.py <output.fits> <mask1.fits> <mask2.fits> [...]')
    combine_masks(opfits, maskfits)


if __name__ == "__main__":


    main()
