import glob
import shutil
import time
import datetime
import subprocess
import sys
import os
import os.path as o

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if len(sys.argv) == 1:
    print('Please specify a field to split')
    sys.exit()
else:
    targetname = sys.argv[-1]

# Load in the target ms files
target_ms_files = project_info['target_ms']

# If this is a target split from the self-calibrated SPLIT ms file
if any([targetname in target_ms for target_ms in target_ms_files]):
    myms = target_ms_files[target_names.index(targetname)]
    target_ms = myms.replace('.ms', '_snapshot.ms')
else:
    target_ms = myms.replace('.ms', f'_{targetname}_snapshot.ms')

if not o.isdir(myms):
    sys.exit('ERROR: The MS file you are trying to split from does not exist!')


# This is to check if you've already done snapshot imaging
dirs = glob.glob(target_ms.strip('.ms')+'*'+'.ms')
n = len(dirs)

if o.isdir(target_ms):
    os.rename(target_ms, myms.replace('.ms', f'_{targetname}_snapshot'+str(n)+'.ms'))

# Split out corrected_data column of the target of interest
mstransform(vis = myms, 
                        outputvis=target_ms, 
                        field=targetname,
                        datacolumn='corrected')

