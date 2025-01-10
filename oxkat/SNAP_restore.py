import glob
import os
import os.path as o
import sys
import subprocess
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg

def main():

    if len(sys.argv) != 3:
        print('Please specify the path to static model for restoration and target MS name')
        sys.exit()
    else:
        model_prefix = sys.argv[-2]
        myms = sys.argv[-1]

    # Define target_prefix
    target_prefix = cfg.INTERVALS + f'/img_{myms}_modelsub'

    # Configuration options    
    chanout = cfg.SNAP_CHANNELSOUT
    pol = 'I'
    if cfg.SNAP_POL:
        pol = 'IQUV'
    
    syscall='python3 '+ cfg.TOOLS+f'/restore_model.py {model_prefix} {target_prefix} {pol}'        
    subprocess.run([syscall], shell=True)


if __name__ == "__main__":
    main()
