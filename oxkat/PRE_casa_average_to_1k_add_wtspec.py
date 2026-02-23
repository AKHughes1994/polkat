# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
import numpy as np
import glob, json, re
import functools
print = functools.partial(print, flush=True)

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
    print('PRE: Manual scan selection enabled (PRE_SCANS != "")')
    print(f'PRE: Requested scans: {myscans}')
    myscans = myscans.split(',')
    removed = [scan for scan in myscans if scan in bad_scans]
    myscans = ','.join([scan for scan in myscans if scan not in bad_scans])
    if removed:
        print(f'PRE: Removed {len(removed)} bad/short scan(s) from selection: {removed}')
    else:
        print('PRE: No bad/short scans found in selection')
    print(f'PRE: Final scan selection: {myscans}')
    # Fix bug where CASA retains field names no longer present after scan selection:
    # derive myfields directly from the fields that appear in the selected scans
    tb.open(master_ms)
    scan_list = [int(s) for s in myscans.split(',') if s != '']
    field_ids = set()
    for scan in scan_list:
        subtab = tb.query(query='SCAN_NUMBER=='+str(scan))
        field_ids.update(np.unique(subtab.getcol('FIELD_ID')).tolist())
        subtab.close()
    tb.close()
    # Get the field names corresponding to these IDs
    tb.open(master_ms+'/FIELD')
    all_field_names = tb.getcol('NAME')
    tb.close()
    myfields = ','.join([all_field_names[fid] for fid in sorted(field_ids)])
    print(f'PRE: Derived fields from scan selection (overriding PRE_FIELDS): {myfields}')
    # Write derived fields back into config.py so all downstream scripts see them as PRE_FIELDS
    config_path = OXKAT+'/config.py'
    with open(config_path, 'r') as _f:
        _cfg_content = _f.read()
    _cfg_content = re.sub(r"^PRE_FIELDS\s*=\s*'[^']*'", f"PRE_FIELDS = '{myfields}'", _cfg_content, flags=re.MULTILINE)
    with open(config_path, 'w') as _f:
        _f.write(_cfg_content)
    print(f'PRE: Updated PRE_FIELDS in config.py to: {myfields}')

else:
    myscans = ','.join(good_scans)
    if bad_scans:
        print(f'PRE: Removed {len(bad_scans)} bad/short scan(s) from full scan list: {bad_scans}')
    print(f'PRE: No manual scan selection; using all {len(good_scans)} good scan(s)')

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
