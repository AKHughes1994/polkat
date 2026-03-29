# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import sys

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

# -----------------------------------------------------------------------
# Read calibrator imaging settings from config
# -----------------------------------------------------------------------

cals_to_image = CALS_TO_IMAGE.strip()
combine_scan  = CAL_COMBINESCAN

if not cals_to_image:
    sys.exit('CALS_TO_IMAGE is empty in config — nothing to split.')

field_list = [f.strip() for f in cals_to_image.split(',') if f.strip()]
print(f'Fields to split: {field_list}')
print(f'Combine scans  : {combine_scan}')

# -----------------------------------------------------------------------
# Build a lookup of field name -> list of scan numbers from the working MS
# -----------------------------------------------------------------------

tb.open(myms)
all_fields = tb.getcol('FIELD_ID')
all_scans  = tb.getcol('SCAN_NUMBER')
tb.close()

tb.open(myms + '/FIELD')
field_names = tb.getcol('NAME').tolist()
tb.close()

# Map field name -> set of scan numbers present in the MS
field_scan_map = {}
for fid, scan in zip(all_fields, all_scans):
    name = field_names[fid]
    field_scan_map.setdefault(name, set()).add(int(scan))

# -----------------------------------------------------------------------
# Split
# -----------------------------------------------------------------------

for field in field_list:

    if field not in field_scan_map:
        print(f'WARNING: Field "{field}" not found in {myms} — skipping.')
        continue

    scans = sorted(field_scan_map[field])

    if combine_scan:
        # One MS per field, all scans combined
        scan_str = ','.join(str(s) for s in scans)
        opms = myms.replace('.ms', f'_{field}.ms')
        print(f'Splitting field {field} (all scans) -> {opms}')
        try:
            mstransform(
                vis           = myms,
                outputvis     = opms,
                field         = field,
                scan          = scan_str,
                usewtspectrum = True,
                realmodelcol  = True,
                datacolumn    = 'corrected')
            flagmanager(vis=opms, mode='save', versionname='post-1GC')
        except Exception as e:
            print(f'ERROR splitting {field}: {e}')

    else:
        # One MS per field per scan
        for scan in scans:
            scan_pad = f'{scan:03d}'
            opms = myms.replace('.ms', f'_{field}_scan{scan_pad}.ms')
            print(f'Splitting field {field} scan {scan} -> {opms}')
            try:
                mstransform(
                    vis           = myms,
                    outputvis     = opms,
                    field         = field,
                    scan          = str(scan),
                    usewtspectrum = True,
                    realmodelcol  = True,
                    datacolumn    = 'corrected')
                flagmanager(vis=opms, mode='save', versionname='post-1GC')
            except Exception as e:
                print(f'ERROR splitting {field} scan {scan}: {e}')

print('casa_split_cals.py complete.')
