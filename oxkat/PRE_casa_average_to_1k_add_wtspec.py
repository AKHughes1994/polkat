# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import csv
import datetime
import functools
import glob
import json
import os
import re

import numpy as np

print = functools.partial(print, flush=True)

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

myfields      = PRE_FIELDS
myscans       = PRE_SCANS
myoutputchans = int(PRE_NCHANS)
mytimebins    = PRE_TIMEBIN


master_ms = glob.glob('*.ms')[0]
opms      = master_ms.replace('.ms', '_' + str(myoutputchans) + 'ch.ms')

tb.open(master_ms + '/SPECTRAL_WINDOW')
nchan = tb.getcol('NUM_CHAN')[0]
tb.done()

mychanbin = int(nchan / myoutputchans)
mychanave = mychanbin > 1

# -----------------------------------------------------------------------
# Remove short scans arising from metadata errors in 2s integration data
# -----------------------------------------------------------------------

bad_scans  = []
good_scans = []

tb.open(master_ms)
scans = np.unique(tb.getcol('SCAN_NUMBER'))
for scan in scans:
    subtab      = tb.query(query='SCAN_NUMBER==' + str(scan))
    scan_times  = np.unique(subtab.getcol('TIME'))
    scan_dt     = scan_times[-1] - scan_times[0]   # total scan duration (s)
    integration = scan_times[1]  - scan_times[0]   # integration length (s)
    if scan_dt < 10.0 and integration < 2.5:
        bad_scans.append(str(scan))
    else:
        good_scans.append(str(scan))
tb.close()

if myscans != '':
    print('PRE: Manual scan selection enabled (PRE_SCANS != "")')
    print(f'PRE: Requested scans: {myscans}')
    myscans = myscans.split(',')
    removed = [s for s in myscans if s in bad_scans]
    myscans = ','.join([s for s in myscans if s not in bad_scans])
    if removed:
        print(f'PRE: Removed {len(removed)} bad/short scan(s) from selection: {removed}')
    else:
        print('PRE: No bad/short scans found in selection')
    print(f'PRE: Final scan selection: {myscans}')

    # Fix CASA bug where field names linger after scan selection —
    # derive myfields from the fields actually present in the selected scans
    tb.open(master_ms)
    scan_list = [int(s) for s in myscans.split(',') if s != '']
    field_ids = set()
    for scan in scan_list:
        subtab = tb.query(query='SCAN_NUMBER==' + str(scan))
        field_ids.update(np.unique(subtab.getcol('FIELD_ID')).tolist())
        subtab.close()
    tb.close()

    tb.open(master_ms + '/FIELD')
    all_field_names = tb.getcol('NAME')
    tb.close()

    myfields = ','.join([all_field_names[fid] for fid in sorted(field_ids)])
    print(f'PRE: Derived fields from scan selection (overriding PRE_FIELDS): {myfields}')

    # Write derived fields back into config.py so all downstream scripts see them
    config_path = OXKAT + '/config.py'
    with open(config_path, 'r') as _f:
        _cfg_content = _f.read()
    _cfg_content = re.sub(
        r"^PRE_FIELDS\s*=\s*'[^']*'",
        f"PRE_FIELDS = '{myfields}'",
        _cfg_content,
        flags=re.MULTILINE)
    with open(config_path, 'w') as _f:
        _f.write(_cfg_content)
    print(f'PRE: Updated PRE_FIELDS in config.py to: {myfields}')

else:
    myscans = ','.join(good_scans)
    if bad_scans:
        print(f'PRE: Removed {len(bad_scans)} bad/short scan(s) from full scan list: {bad_scans}')
    print(f'PRE: No manual scan selection; using all {len(good_scans)} good scan(s)')


# -----------------------------------------------------------------------
# Transform MS
# -----------------------------------------------------------------------

mstransform(
    vis           = master_ms,
    outputvis     = opms,
    field         = myfields,
    scan          = myscans,
    datacolumn    = 'data',
    chanaverage   = mychanave,
    chanbin       = mychanbin,
    # timeaverage = True,
    # timebin     = '8s',
    realmodelcol  = True,
    usewtspectrum = True)

# Save observatory flags and initialise model column
flagmanager(vis=opms, mode='save', versionname='observatory')
clearcal(vis=opms, addmodel=True)

# Read field names and IDs from the output MS —
# mstransform renumbers field IDs when PRE_FIELDS != '', so we re-derive them here
tb.open(opms + '/FIELD')
names = tb.getcol('NAME')
ids   = tb.getcol('SOURCE_ID')
tb.done()

# Persist working names and IDs into project_info.json for downstream scripts
with open('project_info.json', 'r') as j:
    project_info = json.load(j)

project_info['working_names'] = names.tolist()
project_info['working_ids']   = ids.tolist()

with open('project_info.json', 'w') as j:
    json.dump(project_info, j, indent=4, sort_keys=True)

clearstat()
clearstat()


# -----------------------------------------------------------------------
# Write comprehensive per-field/scan timing CSV to RESULTS
# -----------------------------------------------------------------------

def mjd_to_utc(mjd):
    """Convert MJD (days) to a UTC ISO string."""
    unix = (mjd - 40587.0) * 86400.0
    return datetime.datetime.utcfromtimestamp(unix).strftime('%Y-%m-%dT%H:%M:%S')


if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)

# Load field name lookup from the output MS
tb.open(opms + '/FIELD')
field_names = tb.getcol('NAME')
tb.close()

# Read all field IDs, scan numbers, and timestamps from the output MS in one pass
tb.open(opms)
all_field_ids = tb.getcol('FIELD_ID')
all_scan_nums = tb.getcol('SCAN_NUMBER')
all_times     = tb.getcol('TIME')
tb.close()

# Iterate over every unique field+scan combination and extract timing
time_info = []
for field_id in np.unique(all_field_ids):
    for scan in np.unique(all_scan_nums[all_field_ids == field_id]):
        mask     = (all_field_ids == field_id) & (all_scan_nums == scan)
        t_times  = all_times[mask]
        t_start  = t_times.min() / 86400.0
        t_end    = t_times.max() / 86400.0
        duration = round((t_end - t_start) * 1440.0, 3)
        time_info.append((
            opms,
            field_names[field_id],
            int(field_id),
            int(scan),
            f'{t_start:.10f}',
            f'{t_end:.10f}',
            mjd_to_utc(t_start),
            mjd_to_utc(t_end),
            duration,
        ))

csv_path = os.path.join(RESULTS, opms.rstrip('/') + '_time_info.csv')

with open(csv_path, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['ms', 'field', 'field_id', 'scan', 'start_mjd', 'end_mjd',
                     'start_utc', 'end_utc', 'duration_min'])
    writer.writerows(time_info)

print(f'Timing information ({len(time_info)} field/scan entries) saved to: {csv_path}')
