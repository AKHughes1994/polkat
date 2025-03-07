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

if PRE_FIELDS != '':
    target_names = user_targets

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
