import glob
import os
import os.path as o
import sys

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if len(sys.argv) == 1:
    print('Please specify a field to split')
    sys.exit()
else:
    targetname = sys.argv[-1]

# -----------------------------------------------------------------------
# Determine the best source MS to split from, in priority order:
#   1. stage2 MS  (<target>_stage2.ms) — produced by two-stage selfcal
#   2. Target MS  from project_info    — produced by 1GC_08 split
#   3. Working MS (working_ms)         — the averaged master MS, as a last resort
# -----------------------------------------------------------------------

target_ms_files = project_info['target_ms']

# Build the expected stage2 MS name: same convention as setup_2GC_twostage.py
if any(targetname in tms for tms in target_ms_files):
    project_target_ms = target_ms_files[target_names.index(targetname)]
    stage2_ms = project_target_ms.replace('.ms', '_stage2.ms')
else:
    project_target_ms = None
    stage2_ms = None

# Walk the priority list and pick the first MS that exists on disk
if stage2_ms and o.isdir(stage2_ms):
    source_ms = stage2_ms
    print(f'SNAP_split_sources: Using stage2 MS: {source_ms}')
elif project_target_ms and o.isdir(project_target_ms):
    source_ms = project_target_ms
    print(f'SNAP_split_sources: Using project target MS: {source_ms}')
elif o.isdir(myms):
    source_ms = myms
    print(f'SNAP_split_sources: Falling back to working MS: {source_ms}')
else:
    sys.exit(f'ERROR: No suitable MS found for target "{targetname}". '
             f'Checked: {stage2_ms}, {project_target_ms}, {myms}')

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
