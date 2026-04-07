import os
import sys
import glob
import subprocess
import time
import numpy as np

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt)


def fix_names(prefix, totchan):
    """
    Read in the input prefix and fix the names of the images;
    this is necessary if imaging the bandwidth has been broken up due to
    memory constraints.

    Infers the number of channels per part directly from the files found,
    rather than relying on the global WSC_MAX_CHANNELS config value.

    Arguments:
        prefix:  string, the input prefix for the wsclean imaging call
        totchan: int, total number of channels across all parts
    Returns:
        NOTHING
    """

    part_pattern = f'{prefix}_part*.fits'
    part_files   = glob.glob(part_pattern)

    if not part_files:
        msg(f'No "part" files found matching pattern: {part_pattern}')
        msg('Nothing to rename — exiting')
        return

    # -------------------------------------------------------------------------
    # Infer number of parts and channels per part from what is on disk
    # -------------------------------------------------------------------------
    part_numbers = set()
    for f in part_files:
        basename = os.path.basename(f)
        # extract the part index from '_partXXXX'
        part_str = basename.replace(os.path.basename(prefix) + '_part', '').split('-')[0]
        try:
            part_numbers.add(int(part_str))
        except ValueError:
            msg(f'WARNING: could not parse part index from filename: {basename}')

    num_parts = len(part_numbers)
    msg(f'Found {num_parts} image part(s) covering {totchan} total channels '
        f'(parts: {sorted(part_numbers)})')

    if totchan % num_parts != 0:
        msg(f'ERROR: totchan ({totchan}) is not evenly divisible by the number '
            f'of parts found ({num_parts}) — check inputs and rerun')
        sys.exit(1)

    maxchan = totchan // num_parts
    msg(f'Inferred {maxchan} channel(s) per part')

    # -------------------------------------------------------------------------
    # Remove MFS images produced during the parted run (they are invalid)
    # -------------------------------------------------------------------------
    mfs_files = glob.glob(f'{prefix}_part*-MFS-*')
    if mfs_files:
        msg(f'Removing {len(mfs_files)} invalid MFS file(s) from parted run')
        subprocess.run([f'rm -rf {prefix}_part*-MFS-*'], shell=True)
    else:
        msg('No MFS files to remove')

    # -------------------------------------------------------------------------
    # Re-glob non-MFS fits files after cleanup
    # -------------------------------------------------------------------------
    images = glob.glob(f'{prefix}_part*.fits')
    if not images:
        msg('ERROR: No part FITS files remain after MFS cleanup — nothing to rename')
        sys.exit(1)

    msg(f'Preparing to rename {len(images)} image file(s)')
    images = np.array([images, images])   # row 0 = old names, row 1 = new names

    # -------------------------------------------------------------------------
    # Build new filenames by offsetting channel indices per part
    # -------------------------------------------------------------------------
    rename_count = 0
    for k in range(images.shape[1]):
        image_new = images[1, k][:]

        # Identify which part this file belongs to
        factor = None
        leading_substring = None
        for part in range(num_parts):
            candidate = f'{prefix}_part{part:04d}'
            if candidate in image_new:
                leading_substring = candidate
                factor = part
                break

        if factor is None or leading_substring is None:
            msg(f'WARNING: could not identify part for file: {image_new} — skipping')
            continue

        # Handle time-interval naming (-tXXXX-YYYY-) vs plain channel naming
        if glob.glob(f'{leading_substring}-t*'):
            # e.g. prefix_part0001-t0000-0003-image.fits
            indexing_parts = image_new.split(f'{leading_substring}-t')[-1].split('-')[:2]
            indexing_old   = '-'.join(indexing_parts)
            indexing_parts[1] = '{:04d}'.format(int(indexing_parts[1]) + factor * maxchan)
            indexing_new   = '-'.join(indexing_parts)
        else:
            # e.g. prefix_part0001-0003-image.fits
            indexing_old = image_new.split(f'{leading_substring}-')[-1].split('-')[0]
            indexing_new = '{:04d}'.format(int(indexing_old) + factor * maxchan)

        image_new = (image_new
                     .replace(leading_substring, prefix)
                     .replace(indexing_old, indexing_new))
        images[1, k] = image_new
        rename_count += 1

    # -------------------------------------------------------------------------
    # Execute renames
    # -------------------------------------------------------------------------
    msg(f'Renaming {rename_count} file(s)...')
    for old, new in images.T:
        if old != new:
            subprocess.run([f'mv {old} {new}'], shell=True)

    msg('Image renaming complete')


def main():

    if len(sys.argv) != 3:
        msg('ERROR: Expected exactly 2 arguments: <totchan> <prefix>')
        msg('Usage: fix_image_naming.py <totchan> <prefix>')
        sys.exit(1)

    totchan = int(sys.argv[-2])
    prefix  = sys.argv[-1]

    msg(f'Running fix_image_naming for prefix: {prefix}')
    msg(f'Total channels expected: {totchan}')

    fix_names(prefix, totchan)


if __name__ in '__main__':
    main()
