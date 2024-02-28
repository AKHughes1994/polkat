import sys
import glob
import os.path as o
import subprocess

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg

def main():

    if len(sys.argv) == 1:
        print('Please specify the name of the field to perform uvsubtraction')
        sys.exit()
    else:
        targetname = sys.argv[-1]

    # Look for all scan ms files associated with the target of interest
    mslist = glob.glob(f'*{targetname}*scan*.ms')

    for myms in mslist:
        syscall = f'python3 {cfg.TOOLS}/sum_MS_columns.py --src=MODEL_DATA --dest=DATA --subtract '+myms
        subprocess.run([syscall], shell=True)

if __name__ == "__main__":
    main()
