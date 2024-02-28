import glob
import shutil
import time
import datetime
import subprocess
import sys


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if len(sys.argv) == 1:
    print('Please specify a field to split')
    sys.exit()
else:
    targetname = sys.argv[-1]

target_ms = myms.replace('.ms', f'_{targetname}.ms')

# Split out corrected_data column of the target of interest
mstransform(vis = myms, 
                        outputvis=target_ms, 
                        field=targetname,
                        datacolumn='corrected')
    

