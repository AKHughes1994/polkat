# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
import numpy as np
import glob, json

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

myfields = PRE_FIELDS
myscans = PRE_SCANS
myoutputchans = int(PRE_NCHANS)
mytimebins = PRE_TIMEBIN

# Load in the master ms
opms = master_ms.replace('.ms','_'+str(myoutputchans)+'ch.ms')

flagmanager(opms, versionname='after_swap_and_zero', mode='restore')

# Remove old flag versions
for fname in ['1GC_flags', 'after_pcal_SwiftJ1727', 'applycal_1', 'autoflag_cals_data', 'basic', 'bpcal_residual_flags']:
    try:
        flagmanager(opms, versionname=fname, mode='delete')
    except:
        pass

# Get  names and field IDs for sources that are 
tb.open(opms+'/FIELD')
names = tb.getcol('NAME')
ids   = tb.getcol('SOURCE_ID')
tb.done()

# Append the working names and IDs to project info as mstranform will modify the field IDs if PRE_FIELDS != ''
with open('project_info.json','r') as j:
    project_info = json.load(j)

project_info['working_names'] = names.tolist()
project_info['working_ids'] = ids.tolist()


with open('project_info.json','w') as j:
    json.dump(project_info, j, indent=4, sort_keys = True)
