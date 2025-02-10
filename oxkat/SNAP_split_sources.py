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

target_ms = myms.replace('.ms', f'_{targetname}_snapshot.ms')

dirs = glob.glob(target_ms.strip('.ms')+'*'+'.ms')
n = len(dirs)

if o.isdir(target_ms):
    os.rename(target_ms, myms.replace('.ms', f'_{targetname}_snapshot'+str(n)+'.ms'))

# Split out corrected_data column of the target of interest
mstransform(vis = myms, 
                        outputvis=target_ms, 
                        field=targetname,
                        datacolumn='corrected')
    

