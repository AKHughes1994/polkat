#!/usr/bin/env python
# Check whether a QuartiCal self-calibration solve succeeded. If it failed
# (as detected by absence of the expected log marker), re-run without
# parallactic angle correction, then apply parang via CASA applycal as a
# separate final step.
#
# Exception: with --parangmodel set (i.e. CAL_2GC_PARANGMODEL was True for
# this run), a failed solve almost always means a known CASA ephemeris quirk
# has snuck into QuartiCal's deterministic parang solver and made it choke.
# Every self-cal model image (WSDMA/WSCMI) already relied on that same
# solver, so the fallback below isn't used; this prints an error and exits
# non-zero instead.
#
# Usage:
#   python check_and_fix_parang_selfcal.py <log_outdir> <myms> <fallback_qc_command> [--parangmodel]
#
#   log_outdir          - Directory containing QuartiCal log output
#   myms                - Target MS path (used for split and applycal)
#   fallback_qc_command - Full QuartiCal command to run WITHOUT apply_p_jones_inv
#                         (passed as a single quoted string)
#   --parangmodel        - Pass this when CAL_2GC_PARANGMODEL was True for the
#                         run being checked (see exception above)

import glob
import os
import shutil
import subprocess
import sys
import time
import functools

print = functools.partial(print, flush=True)

SUCCESS_MARKER = 'Final post-solve chi-squared summary'


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt)


# -----------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------

parangmodel = '--parangmodel' in sys.argv
positional = [a for a in sys.argv[1:] if a != '--parangmodel']

if len(positional) != 3:
    sys.exit(
        'Usage: python check_and_fix_parang_selfcal.py '
        '<log_outdir> <myms> <fallback_qc_command> [--parangmodel]'
    )

log_outdir        = positional[0]
myms              = positional[1]
fallback_qc_cmd   = positional[2]

# Strip any singularity exec prefix upfront — this script runs inside the
# container already, so QuartiCal must be called directly
import re as _re
fallback_qc_cmd = _re.sub(r'^singularity\s+exec\s+\S+\s+', '', fallback_qc_cmd.strip())
msg(f'Fallback QC command: {fallback_qc_cmd}')


# -----------------------------------------------------------------------
# Check QuartiCal log for success marker
# -----------------------------------------------------------------------

# Prefer QuartiCal's native .log.qc files; fall back to plain .log/.txt
log_files = sorted(glob.glob(os.path.join(log_outdir, '*.log.qc')))
if not log_files:
    log_files = sorted(glob.glob(os.path.join(log_outdir, '*.log')) +
                       glob.glob(os.path.join(log_outdir, '*.txt')))

solve_succeeded = False

if not log_files:
    msg(f'WARNING: No log files found in {log_outdir}')
else:
    for logf in log_files:
        try:
            with open(logf) as f:
                content = f.read()
            if SUCCESS_MARKER in content:
                solve_succeeded = True
                msg(f'Solve succeeded (marker found in {os.path.basename(logf)})')
                break
        except OSError as e:
            msg(f'WARNING: Could not read {logf}: {e}')

if solve_succeeded:
    msg('QuartiCal parang solve verified OK — no intervention needed.')
    sys.exit(0)


# -----------------------------------------------------------------------
# Solve failed — re-run without apply_p_jones_inv
# -----------------------------------------------------------------------

msg('SUCCESS MARKER NOT FOUND — solve likely failed.')

if parangmodel:
    msg('ERROR: The deterministic parallactic angle solve failed — a known '
        'CASA ephemeris quirk sneaking into QuartiCal is the usual culprit. '
        'Set CAL_1GC_APPLYPARANG = True or CAL_2GC_PARANGMODEL = False in '
        'config.py, then re-run 2GC.')
    sys.exit(1)

msg('Re-running QuartiCal without parallactic angle correction...')

ret = subprocess.run(fallback_qc_cmd, shell=True)
if ret.returncode != 0:
    sys.exit(f'ERROR: Fallback QuartiCal command failed with return code {ret.returncode}')

msg('Fallback QuartiCal solve complete.')


# -----------------------------------------------------------------------
# Rename original MS, split CORRECTED_DATA into a fresh MS at original path
# -----------------------------------------------------------------------

pre_parang_ms = myms.rstrip('/').replace('.ms', '_pre_parang.ms')

msg(f'Renaming {myms} -> {pre_parang_ms}')
if os.path.exists(pre_parang_ms):
    shutil.rmtree(pre_parang_ms)
shutil.move(myms, pre_parang_ms)

msg(f'Splitting CORRECTED_DATA from {pre_parang_ms} into {myms}')

try:
    from casatasks import mstransform, applycal
    use_casatasks = True
except ImportError:
    use_casatasks = False
    msg('WARNING: casatasks not importable — falling back to subprocess CASA call')

if use_casatasks:
    mstransform(
        vis          = pre_parang_ms,
        outputvis    = myms,
        datacolumn   = 'corrected',
        usewtspectrum = True,
        realmodelcol  = True,
    )
else:
    casa_script = (
        f"mstransform(vis='{pre_parang_ms}', outputvis='{myms}', "
        f"datacolumn='corrected', usewtspectrum=True, realmodelcol=True)"
    )
    ret = subprocess.run(['casa', '--nologger', '--nogui', '-c', casa_script])
    if ret.returncode != 0:
        sys.exit(f'ERROR: mstransform failed with return code {ret.returncode}')


# -----------------------------------------------------------------------
# Apply parallactic angle correction via CASA applycal (no gaintables)
# -----------------------------------------------------------------------

msg(f'Applying parallactic angle correction to {myms} via applycal (parang=True)...')

if use_casatasks:
    applycal(
        vis      = myms,
        parang   = True,
        flagbackup = False,
    )
else:
    casa_script = (
        f"applycal(vis='{myms}', parang=True, flagbackup=False)"
    )
    ret = subprocess.run(['casa', '--nologger', '--nogui', '-c', casa_script])
    if ret.returncode != 0:
        sys.exit(f'ERROR: applycal failed with return code {ret.returncode}')

msg('Parallactic angle correction applied successfully.')
msg(f'Pipeline can continue — {myms} is ready.')
