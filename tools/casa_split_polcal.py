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

if pacal_name == '':
    sys.exit('No polarization calibration MS specified in project info; cannot proceed with polcal splitting.')


# Split the polcal MS into a separate MS for the specified field
try:
    opms = myms.replace('.ms', f'_{pacal_name}.ms')
    
    mstransform(vis=myms,
        outputvis=opms,
        field=pacal_name,
        usewtspectrum=True,
        realmodelcol=True,
        datacolumn='corrected')

    flagmanager(vis=opms,
        mode='save',
        versionname='post-1GC')

except Exception as e:
    sys.exit('Error during polcal MS transformation: ' + str(e))
