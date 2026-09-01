import glob
import os
import os.path as o
import sys

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if len(sys.argv) < 3:
    print('Usage: SNAP_split_sources.py <targetname> <source_ms>')
    sys.exit()
else:
    targetname = sys.argv[-2]
    source_ms  = sys.argv[-1]

# The source MS priority (stage2 -> project target MS -> working MS) is
# resolved once, at setup time, by SNAP.py -- it's passed in directly here
# so a missing MS is caught before any jobs are generated, not at runtime.
if not o.isdir(source_ms):
    sys.exit(f'ERROR: source MS "{source_ms}" for target "{targetname}" does not exist.')

# -----------------------------------------------------------------------
# Determine the output snapshot MS name
# -----------------------------------------------------------------------

# Snapshot MS is always named relative to the working MS so that SNAP.py
# can locate it consistently regardless of which source MS was used
target_ms = myms.replace('.ms', f'_{targetname}_snapshot.ms')

# If a snapshot MS already exists from a previous run, archive it with a
# numeric suffix so the new split doesn't collide with it
if o.isdir(target_ms):
    existing = glob.glob(target_ms.rstrip('.ms') + '*.ms')
    n = len(existing)
    archived = myms.replace('.ms', f'_{targetname}_snapshot{n}.ms')
    print(f'SNAP_split_sources: Archiving existing snapshot MS to {archived}')
    os.rename(target_ms, archived)

# -----------------------------------------------------------------------
# Split the CORRECTED_DATA column for the requested field
# -----------------------------------------------------------------------

print(f'SNAP_split_sources: Splitting field "{targetname}" from {source_ms} -> {target_ms}')

mstransform(
    vis        = source_ms,
    outputvis  = target_ms,
    field      = targetname,
    datacolumn = 'corrected')
