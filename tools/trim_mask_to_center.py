#!/usr/bin/env python3
"""
Trim a FITS mask file to the central NxN pixels.

This tool is useful when using the zoom step (smaller image size in WSCMA).
It takes a mask file created with the full image size and trims it to match
the smaller zoom image size, keeping only the central NxN pixels.

Usage:
    python trim_mask_to_center.py <input_mask.fits> <output_size> [output_mask.fits]

Arguments:
    input_mask.fits  : Input FITS mask file to trim
    output_size      : Size of the output image (NxN pixels)
    output_mask.fits : Optional output filename. If not provided, will use
                       input filename with '_trimmed' suffix

Example:
    python trim_mask_to_center.py my_mask-MFS-image.mask.fits 2560
    # Creates: my_mask-MFS-image.mask_trimmed.fits

    python trim_mask_to_center.py my_mask.fits 2560 my_mask_small.fits
    # Creates: my_mask_small.fits
"""

import sys
import os
import numpy as np
from astropy.io import fits


def trim_mask_to_center(input_file, output_size, output_file=None):
    """
    Trim a FITS mask file to the central NxN pixels.
    
    Parameters:
    -----------
    input_file : str
        Path to input FITS mask file
    output_size : int
        Size of the output square image (NxN pixels)
    output_file : str, optional
        Path to output FITS file. If None, creates a file with '_trimmed' suffix
    
    Returns:
    --------
    output_file : str
        Path to the created output file
    """
    
    # Read the input FITS file
    with fits.open(input_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
        
        # Get the shape of the input data
        if data.ndim == 2:
            input_shape = data.shape
        elif data.ndim == 3:
            # Handle 3D data (e.g., with frequency axis)
            input_shape = data.shape[-2:]
        elif data.ndim == 4:
            # Handle 4D data (e.g., with stokes and frequency axes)
            input_shape = data.shape[-2:]
        else:
            raise ValueError(f"Unexpected data dimensions: {data.ndim}. Expected 2, 3, or 4.")
        
        input_height, input_width = input_shape
        
        # Check that output size is not larger than input
        if output_size > min(input_height, input_width):
            raise ValueError(f"Output size {output_size} is larger than input size {input_shape}")
        
        # Calculate the center coordinates
        center_y = input_height // 2
        center_x = input_width // 2
        
        # Calculate the extraction boundaries
        half_size = output_size // 2
        y_start = center_y - half_size
        y_end = y_start + output_size
        x_start = center_x - half_size
        x_end = x_start + output_size
        
        # Extract the central region
        if data.ndim == 2:
            trimmed_data = data[y_start:y_end, x_start:x_end]
        elif data.ndim == 3:
            trimmed_data = data[:, y_start:y_end, x_start:x_end]
        elif data.ndim == 4:
            trimmed_data = data[:, :, y_start:y_end, x_start:x_end]
        
        # Update the header with new dimensions
        # Update NAXIS1 and NAXIS2 (or their equivalents)
        if 'NAXIS1' in header:
            header['NAXIS1'] = output_size
        if 'NAXIS2' in header:
            header['NAXIS2'] = output_size
            
        # Update CRPIX values to account for the shift in reference pixel
        if 'CRPIX1' in header:
            header['CRPIX1'] = header['CRPIX1'] - x_start
        if 'CRPIX2' in header:
            header['CRPIX2'] = header['CRPIX2'] - y_start
        
        # Add history information
        header.add_history(f'Trimmed to central {output_size}x{output_size} pixels')
        header.add_history(f'Original size: {input_width}x{input_height}')
        header.add_history(f'Trimmed using trim_mask_to_center.py')
        
        # Create output filename if not provided
        if output_file is None:
            base, ext = os.path.splitext(input_file)
            output_file = f"{base}_trimmed{ext}"
        
        # Write the trimmed data to a new FITS file
        hdu = fits.PrimaryHDU(data=trimmed_data, header=header)
        hdu.writeto(output_file, overwrite=True)
        
        print(f"Input file: {input_file}")
        print(f"  Input shape: {input_shape} (HxW)")
        print(f"  Output shape: {output_size}x{output_size}")
        print(f"  Extracted region: Y[{y_start}:{y_end}], X[{x_start}:{x_end}]")
        print(f"Output file: {output_file}")
        
        return output_file


def main():
    """Main function to parse command line arguments and run the trimming."""
    
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    try:
        output_size = int(sys.argv[2])
    except ValueError:
        print(f"Error: output_size must be an integer, got '{sys.argv[2]}'")
        sys.exit(1)
    
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Check that input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    try:
        trim_mask_to_center(input_file, output_size, output_file)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
