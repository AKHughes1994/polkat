import glob
import shutil
import time
import datetime
import subprocess
import sys
import os
import os.path as o
import re


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if len(sys.argv) == 1:
    print('Please specify a field to split')
    sys.exit()
else:
    targetname = sys.argv[-1]

# Build list of all MS files to search through
all_ms_files = []

# Add target MS files
all_ms_files.extend(project_info['target_ms'])

# Add primary calibrator MS files
all_ms_files.extend(project_info['primary_ms'])

# Add secondary calibrator MS files
for pcal_ms_list in project_info['secondary_ms']:
    all_ms_files.extend(pcal_ms_list)

# Add polarization angle calibrator MS files if present
if 'polang_ms' in project_info and project_info['polang_ms']:
    all_ms_files.extend(project_info['polang_ms'])

# Search through all MS files to find the one containing this field
myms = None
for ms_file in all_ms_files:
    # Check if targetname matches the field in this MS
    # For targets: exact match with target name
    # For calibrators: match with fieldname_scanXX pattern
    if targetname in ms_file or f'_{targetname}_' in ms_file or ms_file.endswith(f'_{targetname}.ms'):
        myms = ms_file
        break

if myms is None:
    sys.exit(f'ERROR: Could not find MS file containing field {targetname}')

# Remove _scanXX suffix from fieldname if present (for calibrators)
fieldname_no_scan = re.sub(r'_scan\d+$', '', targetname)

target_ms = myms.replace('.ms', '_snapshot.ms')

if not o.isdir(myms):
    sys.exit('ERROR: The MS file you are trying to split from does not exist!')


# This is to check if you've already done snapshot imaging
# If target_ms exists, rename it to preserve the old snapshot before creating a new one
if o.isdir(target_ms):
    print(f'NOTE: {target_ms} already exists. Renaming it to preserve the old snapshot run.')
    dirs = glob.glob(target_ms.strip('.ms')+'*'+'.ms')
    n = len(dirs)
    os.rename(target_ms, target_ms.replace('.ms', f'_old{n:02d}.ms'))

# Split out corrected_data column of the target of interest
mstransform(vis = myms, 
                        outputvis=target_ms, 
                        field=fieldname_no_scan,
                        datacolumn='corrected')

