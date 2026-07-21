import glob
import os
import os.path as o
import sys
import subprocess
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg

def main():

    if len(sys.argv) < 3:
        print('Please specify the path to static model for restoration, target MS name, and optionally target name')
        sys.exit()
    else:
        model_prefix = sys.argv[1]
        myms = sys.argv[2]
        
        # Get target name from command line or derive from MS name
        if len(sys.argv) >= 4:
            targetname = sys.argv[3]
        else:
            # Try to extract target name from MS filename
            targetname = myms.split('/')[-1].replace('_snapshot.ms', '').split('_', 1)[-1]
    
    # Create target-specific subdirectory under INTERVALS
    # (should already exist from SNAP_intervals.py, but check to be safe)
    filename_targetname = gen.scrub_target_name(targetname)
    intervals_subdir = cfg.INTERVALS + f'/{filename_targetname}'
    
    if not os.path.exists(intervals_subdir):
        print(f'ERROR: Intervals subdirectory does not exist: {intervals_subdir}')
        print('This should have been created by SNAP_intervals.py')
        sys.exit(1)
    
    # Define target_prefix in subdirectory
    target_prefix = intervals_subdir + f'/img_{myms}_modelsub'

    # Configuration options    
    chanout = cfg.SNAP_CHANNELSOUT
    pol = 'I'
    if cfg.SNAP_POL:
        pol = 'IQUV'
    
    syscall='python3 '+ cfg.TOOLS+f'/restore_model.py {model_prefix} {target_prefix} {pol}'        
    subprocess.run([syscall], shell=True)


if __name__ == "__main__":
    main()
