import glob
import os
import os.path as o
import sys
import subprocess
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg

def main():

    if len(sys.argv) == 1:
        print('Please specify the name of the target field for model restoration')
        sys.exit()
    else:
        targetname = sys.argv[1]
    
    # Get the full path for the restoring model
    model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask-MFS-model.fits')[0]

    syscall = 'python3 '+ cfg.TOOLS+f'/restore_model.py {model_fits} {targetname}'
    subprocess.run([syscall], shell=True)


if __name__ == "__main__":
    main()
