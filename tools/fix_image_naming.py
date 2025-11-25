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
    print(stamp+txt)

def fix_names(prefix, totchan):
    """
    Read in the input prefix and fix the names of the images;
    this is necessary if imaging the bandwidth has been broken up due to
    memory constraints

    Arguments:
        prefix: a string of the input prefix for the wsclean imaging call
    Returns:
        NOTHING
    """

    maxchan = cfg.WSC_MAX_CHANNELS

    # Remove the 'parts' component to homogenize the naming conventions
    if maxchan < totchan and len(glob.glob(f'{prefix}_part*')) != 0: 
        msg('Found "part" naming; removing MFS files and standardizing naming convention')

        # Remove the bad MFS images:
        subprocess.run([f'rm -rf {prefix}_part*-MFS-*'], shell = True)

        images = glob.glob(f'{prefix}_part*.fits')
        images = np.array([images, images])  

        # There has to be a better way to do this -- AKH
        parts = int(totchan/maxchan)
        for k in range(images.shape[1]):
            image_new = images[1,k][:]
            for part in range(parts):
                if f'{prefix}_part{part:04d}' in image_new:
                    leading_substring = f'{prefix}_part{part:04d}'
                    factor = part

            # Check if there are intervals != 1 and thus -t000- strings in names
            if glob.glob(f'{leading_substring}-t*') != []:
                indexing_substring = image_new.split(f'{leading_substring}-t')[-1].split('-')[:2]
                indexing_substring_old = '-'.join(indexing_substring)
                indexing_substring[1] = '{:04d}'.format(int(indexing_substring[1]) + factor * maxchan)
                indexing_substring = '-'.join(indexing_substring)

            else:
                indexing_substring = image_new.split(f'{leading_substring}-')[-1].split('-')[0]
                indexing_substring_old = indexing_substring[:]
                indexing_substring = '{:04d}'.format(int(indexing_substring) + factor * maxchan)
                    
            image_new = image_new.replace(leading_substring, prefix).replace(indexing_substring_old, indexing_substring)
            images[1, k] = image_new

        for image_pair in images.T:
            syscall = f'mv {image_pair[0]} {image_pair[1]}'
            subprocess.run([syscall], shell=True)

    msg('Image names have been corrected')


def main():

    # Read in prefix and total number of channels return error if missing it
    if len(sys.argv) != 3:
        msg('ERROR: Missing image prefix, total number of channels, or both')
        sys.exit()
    prefix = sys.argv[-1]
    totchan = int(sys.argv[-2])

    fix_names(prefix, totchan)

if __name__ in '__main__':
    main()
