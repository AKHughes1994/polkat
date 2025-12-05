# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import subprocess
import sys


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

# Apply PRE_FIELDS filtering if specified
if PRE_FIELDS != '':
    target_names = user_targets
    pcal_names = user_pcals

# Split target fields (integrated over all scans)
for target in target_names:

    opms = ''

    for mm in target_ms:
        if target in mm:
            opms = mm

    if opms != '':

        mstransform(vis=myms,
            outputvis=opms,
            field=target,
            usewtspectrum=True,
            realmodelcol=True,
            datacolumn='corrected')

        flagmanager(vis=opms,
            mode='save',
            versionname='post-1GC')

    else:
        print('Target/MS mismatch in project info for '+target+', please check.')


# Split secondary calibrator fields (per scan)
for i, pcal in enumerate(pcal_names):

    # pcal_ms[i] is a list of MS files for each scan of this secondary
    for opms in pcal_ms[i]:
        
        # Extract scan number from MS filename (format: _scanXX.ms)
        scan_str = opms.split('_scan')[-1].replace('.ms', '')
        scan_num = str(int(scan_str))  # Remove zero-padding for CASA field selection
        
        mstransform(vis=myms,
            outputvis=opms,
            field=pcal,
            scan=scan_num,
            usewtspectrum=True,
            realmodelcol=True,
            datacolumn='corrected')

        flagmanager(vis=opms,
            mode='save',
            versionname='post-1GC')
