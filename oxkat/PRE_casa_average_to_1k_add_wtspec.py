# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
import numpy as np
import glob, json

# Flush immediately, so tailing this stage's log shows target selection as it happens
import functools
print = functools.partial(print, flush=True)

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

myfields = PRE_FIELDS
myscans = PRE_SCANS
myoutputchans = int(PRE_NCHANS)
mytimebins = PRE_TIMEBIN


# ------- Target selection
#
# Targets are limited to the sources named in the first column of the XRB
# name-matching list. Any other target field is left out of the averaged MS,
# and out of the target lists in project_info that the rest of the pipeline
# reads, so nothing downstream looks for a field that is not there.
# Calibrators are never excluded here.

XRB_NAME_LIST = DATA+'/positions/XRB_name_matching.txt'

known_targets = set()
if os.path.isfile(XRB_NAME_LIST):
    with open(XRB_NAME_LIST) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                known_targets.add(line.split()[0].lower())
    print(f'\nTarget list: {len(known_targets)} name(s) read from {XRB_NAME_LIST}\n')
    print(f'\nName list contents: {sorted(known_targets)}\n')
else:
    print(f'\nWARNING: {XRB_NAME_LIST} not found -- keeping every target field\n')

# Names are compared without regard to case or surrounding whitespace
if known_targets:
    keep_target = [name.strip().lower() in known_targets for name in target_names]
else:
    keep_target = [True]*len(target_names)

dropped_names = [n for n, k in zip(target_names, keep_target) if not k]
dropped_ids = [i for i, k in zip(targets, keep_target) if not k]
kept_names = [n for n, k in zip(target_names, keep_target) if k]

n_total = len(target_names)
n_kept = len(kept_names)
frac_kept = (n_kept/n_total) if n_total else 0.0
print(f'\n{n_kept} out of {n_total} target field(s) in the MS are in the name '
      f'list ({frac_kept:.1%}): {kept_names}\n')

if dropped_names:
    print(f'\nExcluding {len(dropped_names)} target field(s) absent from the '
          f'name list: {dropped_names}\n')
    if not any(keep_target):
        print('\nWARNING: no target field matched the name list -- the averaged '
              'MS will hold calibrators only\n')

    # Build an explicit field selection. With PRE_FIELDS set, the excluded
    # targets are removed from what the user asked for; otherwise every field
    # but those targets is selected.
    if myfields != '':
        requested = [f.strip() for f in myfields.split(',') if f.strip()]
        myfields = ','.join([f for f in requested if f not in dropped_ids])
    else:
        cal_ids = [bpcal]+pcals
        if pacal not in ('', 'None'):
            cal_ids.append(pacal)
        kept_ids = [i for i, k in zip(targets, keep_target) if k]
        myfields = ','.join(list(dict.fromkeys(cal_ids+kept_ids)))

    print(f'\nmstransform field selection: {myfields}\n')


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

# Drop the excluded targets from the parallel target lists, so the field a
# later stage reads from project_info is one the averaged MS contains
if dropped_names:
    print(f'\nproject_info: trimming target lists to the {len(kept_names)} kept '
          f'field(s), dropping {dropped_names}\n')
    for key in ('target_names','target_ids','target_dirs','target_cal_map','target_ms'):
        values = project_info.get(key,[])
        if len(values) == len(keep_target):
            project_info[key] = [v for v,k in zip(values,keep_target) if k]
        else:
            print(f'\nWARNING: project_info["{key}"] has {len(values)} entries for '
                  f'{len(keep_target)} target(s) -- left untouched\n')

with open('project_info.json','w') as j:
    json.dump(project_info, j, indent=4, sort_keys = True)


clearstat()
clearstat()
