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


master_ms = glob.glob('*.ms')[0]
opms = master_ms.replace('.ms','_'+str(myoutputchans)+'ch.ms')

tb.open(master_ms+'/SPECTRAL_WINDOW')
nchan = tb.getcol('NUM_CHAN')[0]
tb.done()


mychanbin = int(nchan/myoutputchans)
if mychanbin <= 1:
	mychanave = False
else:
	mychanave = True

# Remove short scans that arise from metadata error from 2s integration observations
bad_scans = []
good_scans = []

tb.open(master_ms)
scans = np.unique(tb.getcol('SCAN_NUMBER'))
for scan in scans:
    subtab = tb.query(query='SCAN_NUMBER=='+str(scan)) # scan info
    scan_times = np.unique(subtab.getcol('TIME')) # scan integration times
    scan_dt = scan_times[-1] - scan_times[0] # total scan length (s)
    integration = scan_times[1] - scan_times[0] # integration length (s)
    if scan_dt < 10.0 and integration < 2.5:
        bad_scans.append(str(scan))
    else:
        good_scans.append(str(scan))
tb.close()

if myscans != '':
    myscans = myscans.split(',')
    myscans = ','.join([scan for scan in myscans if scan not in bad_scans])
    
else:
    myscans = ','.join(good_scans)

# Transform MS
mstransform(vis = master_ms,
	outputvis = opms,
	field = myfields,
	scan = myscans,
	datacolumn = 'data',
	chanaverage = mychanave,
	chanbin = mychanbin,
	# timeaverage = True,
	# timebin = '8s',
	realmodelcol = True,
	usewtspectrum = True)

# Save flags
flagmanager(vis = opms, mode = 'save', versionname = 'observatory')
clearcal(vis = opms, addmodel = True)

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


clearstat()
clearstat()
