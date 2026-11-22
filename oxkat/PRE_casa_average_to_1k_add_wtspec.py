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
    integration = float(np.median(subtab.getcol('EXPOSURE')))   # integration length (s)
    scan_dt     = len(scan_times) * integration                # total scan duration (s)
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


# -----------------------------------------------------------------------
# Weight check
# -----------------------------------------------------------------------
#
# Some MSs come out of mstransform with most of WEIGHT_SPECTRUM zero or
# negative, which silently down-weights or drops data in every later solve
# and image. The primary's weights are checked: a weight is good if it is
# positive, and only unflagged samples count, so a heavily flagged primary
# whose flagged samples carry zero weight does not trigger a rebuild.
#
#   WEIGHT_SPECTRUM good                  leave the weights as they are
#   WEIGHT_SPECTRUM bad, WEIGHT good      rebuild WEIGHT_SPECTRUM from WEIGHT
#   both bad                              rebuild both from channel width and
#                                         integration time (wtmode 'nyq'), and warn
#
# "Good" means at least PRE_MIN_GOOD_WEIGHT_FRAC (config.py) of the unflagged
# samples are positive. A rebuild applies to the whole MS, every field.

def good_weight_fractions(ms, field, chunk=2000):
    """
    Positive, zero and negative fractions of the unflagged WEIGHT and
    WEIGHT_SPECTRUM samples of one field, read in chunks of rows.

    A WEIGHT_SPECTRUM sample (corr, chan, row) is unflagged where FLAG is
    False; a WEIGHT sample (corr, row) where at least one of its channels is
    unflagged.

    Returns
    -------
    dict of {str: (float, float, float) or None}
        Column name to its (positive, zero, negative) fractions, or None if
        the field has no unflagged samples.
    """
    # Look the field up by name: mstransform may have renumbered the fields
    tb.open(ms + '/FIELD')
    fid = list(tb.getcol('NAME')).index(field)
    tb.close()

    # positive, zero, negative, total unflagged
    counts = {col: np.zeros(4) for col in ('WEIGHT', 'WEIGHT_SPECTRUM')}
    tb.open(ms)
    t = tb.query(f'FIELD_ID == {fid}')
    for start in range(0, t.nrows(), chunk):
        n = min(chunk, t.nrows() - start)
        flag = t.getcol('FLAG', start, n)                     # (corr, chan, row)
        ok = {'WEIGHT_SPECTRUM': ~flag,
              'WEIGHT': ~flag.all(axis=1)}                    # (corr, row)
        for col, c in counts.items():
            w = t.getcol(col, start, n)[ok[col]]
            c += [np.sum(w > 0), np.sum(w == 0), np.sum(w < 0), w.size]
    t.close()
    tb.close()

    return {col: (tuple(c[:3] / c[3]) if c[3] > 0 else None) for col, c in counts.items()}

def report_weights(ms, label):
    """
    Print the primary's unflagged weight fractions and return the good
    (positive) fraction of WEIGHT_SPECTRUM and of WEIGHT, or None for a
    column with no unflagged samples.
    """
    fractions = good_weight_fractions(ms, bpcal_name)
    print(f'\nWeights of the primary {bpcal_name} ({label}, unflagged samples only):')
    for col, frac in fractions.items():
        if frac is None:
            print(f'  {col:15s}  no unflagged samples')
        else:
            pos, zero, neg = frac
            print(f'  {col:15s}  positive {pos:.2%}  zero {zero:.2%}  negative {neg:.2%}')
    good = {col: (frac[0] if frac is not None else None) for col, frac in fractions.items()}
    return good['WEIGHT_SPECTRUM'], good['WEIGHT']

good_wtsp, good_wt = report_weights(opms, 'after mstransform')
if good_wtsp is None:
    # Nothing unflagged on the primary to judge by; later stages will fail on
    # it anyway, so the weights are left alone rather than rebuilt blind
    print(f'\n!!! WARNING: the primary {bpcal_name} has no unflagged data -- weights '
          f'not checked and left as they are !!!\n')
elif good_wtsp < PRE_MIN_GOOD_WEIGHT_FRAC:
    if good_wt is not None and good_wt >= PRE_MIN_GOOD_WEIGHT_FRAC:
        print(f'\nOnly {good_wtsp:.1%} of WEIGHT_SPECTRUM is positive but WEIGHT is good ({good_wt:.1%}) '
              f'-- reinitialising WEIGHT_SPECTRUM from WEIGHT\n')
        initweights(vis=opms, wtmode='weight', dowtsp=True)
    else:
        print(f'\n!!! WARNING: the weights in this MS look broken: only {good_wtsp:.1%} of WEIGHT_SPECTRUM and '
              f'{good_wt:.1%} of WEIGHT are positive (need {PRE_MIN_GOOD_WEIGHT_FRAC:.0%}). Reinitialising '
              f'both from channel width and integration time (wtmode nyq) -- check this MS !!!\n')
        initweights(vis=opms, wtmode='nyq', dowtsp=True)
    report_weights(opms, 'after initweights')
else:
    print(f'\nWeights look good: {good_wtsp:.1%} of WEIGHT_SPECTRUM and {good_wt:.1%} of WEIGHT are positive '
          f'(need {PRE_MIN_GOOD_WEIGHT_FRAC:.0%}) -- leaving them as they are\n')


# Read field names and IDs from the output MS. mstransform renumbers fields
# when PRE_FIELDS != '', and FIELD_ID in the main table is the row number of
# the FIELD table (not its SOURCE_ID, which can differ), so the IDs are
# re-derived here from the row order
tb.open(opms + '/FIELD')
names = tb.getcol('NAME')
ids   = np.arange(len(names))
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
all_exposures = tb.getcol('EXPOSURE')
tb.close()

# Iterate over every unique field+scan combination and extract timing
time_info = []
for field_id in np.unique(all_field_ids):
    for scan in np.unique(all_scan_nums[all_field_ids == field_id]):
        mask     = (all_field_ids == field_id) & (all_scan_nums == scan)
        t_times  = all_times[mask]
        t_start  = t_times.min() / 86400.0
        t_end    = t_times.max() / 86400.0
        # TIME is the integration midpoint, so the scan lasts one integration more than t_end - t_start
        duration = round(((t_end - t_start) * 86400.0 + np.median(all_exposures[mask])) / 60.0, 3)
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
