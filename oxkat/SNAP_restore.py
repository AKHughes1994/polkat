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
    
    num_channels = cfg.SNAP_CHANS
    pol = cfg.SNAP_POL
    model_freq_splits = cfg.SNAP_FREQ_SPLITS
    
    # Get the full path for the restoring model
    try:
        model_fits = glob.glob(cfg.IMAGES + f'/*{targetname}*snapmask-MFS-model.fits')[0]
    except:
        model_fits = 'default'
    
    if pol:
        print('doing pol')
        syscall = 'python3 '+ cfg.TOOLS+f'/restore_model.py {model_fits} {targetname} {num_channels} {model_freq_splits} IQUV'
    else:
        syscall = 'python3 '+ cfg.TOOLS+f'/restore_model.py {model_fits} {targetname} {num_channels} {model_freq_splits} I'
        
    subprocess.run([syscall], shell=True)


if __name__ == "__main__":
    main()
