#!/usr/bin/env python3
"""
Script to create zoomed cutouts of WSCLEAN FITS images.
Iterates through files matching specified patterns and creates central cutouts.
"""

import os
import glob
from astropy.io import fits
from astropy import wcs
import numpy as np

def ensure_list(item):
    """Convert string to list if needed, otherwise return as-is."""
    if isinstance(item, str):
        return [item]
    return item

def create_zoom_cutout(input_file, output_file, cutout_size=1024):
    """
    Create a cutout of the central region of a FITS image.
    
    Parameters:
    -----------
    input_file : str
        Path to input FITS file
    output_file : str
        Path to output FITS file
    cutout_size : int
        Number of pixels for the square cutout (default: 1024)
    """
    
    # Open the FITS file
    with fits.open(input_file) as hdul:
        header = hdul[0].header.copy()
        data = hdul[0].data
        
        # Preserve original data type
        original_dtype = data.dtype
        
        # Get image dimensions
        # For WSCLEAN images, typical structure is [Stokes, Freq, Dec, RA]
        shape = data.shape
        
        # Find spatial axes (assuming last two dimensions are spatial)
        if len(shape) >= 2:
            ny, nx = shape[-2], shape[-1]
        else:
            raise ValueError(f"Image must have at least 2 dimensions, got {len(shape)}")
        
        # Calculate cutout region (central portion)
        new_nx = cutout_size
        new_ny = cutout_size
        
        # Ensure cutout doesn't exceed image dimensions
        if new_nx > nx:
            new_nx = nx
        if new_ny > ny:
            new_ny = ny
        
        # Ensure even dimensions for clean centering
        if new_nx % 2 != 0:
            new_nx -= 1
        if new_ny % 2 != 0:
            new_ny -= 1
        
        # Calculate starting indices for cutout
        start_x = (nx - new_nx) // 2
        start_y = (ny - new_ny) // 2
        
        # Create cutout - handle different numbers of dimensions
        if len(shape) == 4:  # [Stokes, Freq, Dec, RA]
            cutout_data = data[:, :, start_y:start_y+new_ny, start_x:start_x+new_nx]
        elif len(shape) == 3:  # [Freq, Dec, RA] or similar
            cutout_data = data[:, start_y:start_y+new_ny, start_x:start_x+new_nx]
        elif len(shape) == 2:  # [Dec, RA]
            cutout_data = data[start_y:start_y+new_ny, start_x:start_x+new_nx]
        else:
            raise ValueError(f"Unsupported number of dimensions: {len(shape)}")
        
        # Ensure data type is preserved
        cutout_data = cutout_data.astype(original_dtype)
        
        # Update WCS information in header
        if 'CRPIX1' in header and 'CRPIX2' in header:
            # Update reference pixel coordinates
            header['CRPIX1'] = header['CRPIX1'] - start_x
            header['CRPIX2'] = header['CRPIX2'] - start_y
            
        # Update NAXIS values
        if len(shape) >= 2:
            header['NAXIS1'] = new_nx
            header['NAXIS2'] = new_ny
        
        # Write the cutout to new file
        fits.writeto(output_file, cutout_data, header, overwrite=True)
        print(f"Created zoom cutout: {output_file}")
        print(f"  New image shape: {cutout_data.shape}")

def main():
    """Main function to process FITS images."""
    
    # ================== MODIFIABLE VARIABLES ==================
    
    # Directory containing FITS files (relative to current working directory)
    DIR = 'IMAGES'
    
    # Identifier(s) to match in filenames - can be string or list of strings
    identifiers = ['datablind', 'datamask', 'diagnostic', 'pcalmask', 'uniform']
    
    # Suffix(es) to match - can be string or list of strings  
    suffixes = ['image.fits', 'residual.fits', 'model.fits', 'image.homogenized.fits']
    
    # Cutout size in pixels (square cutout of this size)
    cutout_size = 256  # e.g., 256, 512, 1024
    
    # ===========================================================
    
    # Convert to lists if needed
    identifiers = ensure_list(identifiers)
    suffixes = ensure_list(suffixes)
    
    # Build full directory path
    full_dir = os.path.join(os.getcwd(), DIR)
    
    if not os.path.exists(full_dir):
        print(f"Warning: Directory {full_dir} does not exist!")
        return
    
    # Process each combination of identifier and suffix
    for identifier in identifiers[:]:
        for suffix in suffixes[:]:
            # Create search pattern: *identifier*suffix
            pattern = os.path.join(full_dir, f"*{identifier}*{suffix}")
            
            # Find matching files
            matching_files = glob.glob(pattern)
            
            if not matching_files:
                print(f"No files found matching pattern: {pattern}")
                continue
                
            print(f"\nProcessing {len(matching_files)} files for identifier '{identifier}' with suffix '{suffix}':")
            
            for input_file in matching_files[:]:
                try:
                    # Create output filename by replacing identifier with identifier_zoom
                    base_name = os.path.basename(input_file)
                    output_name = base_name.replace(identifier, f"{identifier}_zoom")
                    output_file = os.path.join(full_dir, output_name)
                    
                    # Create zoom cutout
                    create_zoom_cutout(input_file, output_file, cutout_size)
                    
                    # Delete original file unless it contains 'MFS-*-image' pattern
                    base_name = os.path.basename(input_file)
                    if not ('MFS-' in base_name and '-image' in base_name):
                        os.remove(input_file)
                    
                except Exception as e:
                    print(f"Error processing {input_file}: {str(e)}")
                    continue
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    main()
