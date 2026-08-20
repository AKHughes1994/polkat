#!/usr/bin/env python
# Apply a parallactic angle correction to an MS using CASA's applycal, with
# no gaintables involved (parang=True only). Writes/updates CORRECTED_DATA
# in place; DATA and the MS itself are left untouched.
#
# Usage:
#   python correct_parang_casa.py <myms>

import sys
import subprocess
import time
import functools

print = functools.partial(print, flush=True)


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt)


# -----------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------

if len(sys.argv) != 2:
    sys.exit('Usage: python correct_parang_casa.py <myms>')

myms = sys.argv[1]


# -----------------------------------------------------------------------
# Apply parallactic angle correction via CASA applycal (no gaintables)
# -----------------------------------------------------------------------

msg(f'Applying parallactic angle correction to {myms} via applycal (parang=True)...')

try:
    from casatasks import applycal
    use_casatasks = True
except ImportError:
    use_casatasks = False
    msg('WARNING: casatasks not importable — falling back to subprocess CASA call')

if use_casatasks:
    applycal(
        vis        = myms,
        parang     = True,
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
msg(f'{myms} CORRECTED_DATA now holds parang-corrected data.')
