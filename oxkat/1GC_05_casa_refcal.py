# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
1GC reference calibration: delays, bandpass, gains, leakage and cross-hand
phase from the calibrators, applied to every field.

Overview of the workflow, in the order it runs

    Models          primary model (crystalball or setjy), 1 Jy placeholders
                    for the secondaries and the polang cal
    Primary, stage 0
                    K0, Gp0, B0, G0, Df0 on the primary, applied to it and used
                    to flag its residuals (CORRECTED - MODEL)
    Primary, working tables
                    K, Gp, B, G, Df re-solved on the freshly flagged primary
    Secondaries + polang cal, stage 0
                    Gp0, G0 (and K0) solved against the primary's K, B, Df,
                    applied and used to flag each calibrator
    Secondaries + polang cal, working tables
                    Gp, G (and K) re-solved on the flagged data, then G is
                    bootstrapped onto the primary's flux scale (F table)
    Cross-hand      KCROSS and Xf on the polang cal (CASA solver, or the manual
                    solver in tools/manual_XF_solver.py)
    Applycal        every field, then target flagging

Which field each table is taken from:

    K   the field's own delays, or the primary's with CAL_1GC_DELAY_FROM_PRIMARY
    Gp  the field's own phases (targets: their paired secondary)
    B   always the primary
    G/F the field's own amplitudes (targets: their paired secondary), or the
        primary's with CAL_1GC_AMP_FROM_PRIMARY
    Df  always the primary
    Xf  always the polang cal
"""

import glob
import shutil
import time
import datetime
import subprocess
import sys
import numpy as np

# Flush immediately for better logging
import functools
print = functools.partial(print, flush=True)

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    targets = user_targets
    pcal_names = user_pcals
    target_cal_map = user_cal_map

# Targets paired with a CAL_1GC_BAD_SECONDARIES secondary take their gains from
# the primary: all of them, or with CAL_1GC_REPLACE_SECONDARY_AMP_ONLY just the
# amplitudes. target_cal_map sets the delay and phase field for each target,
# target_amp_cal_map the amplitude field. The secondary itself is still solved,
# flagged and calibrated with its own solutions. tools/manual_XF_solver.py, run
# from this script, re-reads the project info and applies the same mapping.
target_amp_cal_map = list(target_cal_map)
if CAL_1GC_REPLACE_SECONDARY_WITH_PRIMARY:
    replaced_what = 'amplitude gains' if CAL_1GC_REPLACE_SECONDARY_AMP_ONLY else 'all gains'
    for bad_pcal in CAL_1GC_BAD_SECONDARIES:
        if bad_pcal in target_cal_map:
            print(f'Replacing secondary {bad_pcal} with primary {bpcal_name} for target {replaced_what}')
    remapped_cal_map = [bpcal_name if c in CAL_1GC_BAD_SECONDARIES else c for c in target_cal_map]
    target_amp_cal_map = remapped_cal_map
    if not CAL_1GC_REPLACE_SECONDARY_AMP_ONLY:
        target_cal_map = remapped_cal_map

# Amplitude gains from the primary for every field (CAL_1GC_AMP_FROM_PRIMARY in
# config.py). Overrides target_amp_cal_map above; delays and phases unaffected.
if CAL_1GC_AMP_FROM_PRIMARY:
    print(f'Amplitudes: applying the primary ({bpcal_name}) amplitude gains to every field')
    target_amp_cal_map = [bpcal_name] * len(targets)

def amp_field(field):
    """Field whose amplitude (G) solutions are applied to `field`."""
    return bpcal_name if CAL_1GC_AMP_FROM_PRIMARY else field


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------- MANUAL FLAGS ---------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# Last-minute manual flags, applied before any solve. Each entry is one
# flagdata command in CASA list syntax (as in a flagdata inpfile), selecting
# whatever needs to go: antennas, fields, scans, times, channels, baselines.
# They are applied right after the 'autoflag_cals_data' flag version is
# restored below, so they survive into every table and every field, and they
# are kept in every flag version saved from there on. Leave empty for none.
#
# Examples:
#   "antenna='m029'"                                 # an antenna everywhere
#   "antenna='m038' field='J2011-0644' scan='12'"    # one antenna, one scan
#   "antenna='m010&m022'"                            # one baseline
#   "timerange='2026/06/02/22:10:00~22:25:00'"       # a time range, all fields
#   "spw='*:1419.8~1421.3MHz'"                       # a frequency range
MANUAL_FLAGS = []


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# --------------------------- FLAGGING PARAMETERS ---------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# The values below are the settled defaults from a long investigation into a
# break in the Stokes I spectrum above ~1.6 GHz in L-band. Simply raising
# BPCAL_CPARAM_RFLAG_FREQDEVSCALE to 7.0 turned out to be what fixes it.

DEBUG_PRINT_FLAGS = False

# Extend flags across correlations (flagdata rflag/tfcrop 'extendflags' param) --
# applied to calibrators only, not targets.
EXTEND_AUTO = False

# Target visibility rflag/tfcrop thresholds
TARGET_VIS_RFLAG_TIMEDEVSCALE = 5.0   # rflag threshold (sigma), time direction
TARGET_VIS_RFLAG_FREQDEVSCALE = 5.0   # rflag threshold (sigma), freq direction
TARGET_VIS_TFCROP_TIMECUTOFF  = 4.0   # tfcrop threshold (sigma), time direction
TARGET_VIS_TFCROP_FREQCUTOFF  = 3.0   # tfcrop threshold (sigma), freq direction
TARGET_VIS_EXTEND             = False  # whether the mode='extend' pass runs
TARGET_VIS_EXTENDPOLS         = False  # mode='extend' extendpols
TARGET_VIS_EXTEND_TIME        = 90.0  # mode='extend' growtime (% already flagged)
TARGET_VIS_EXTEND_FREQ        = 90.0  # mode='extend' growfreq (% already flagged)
TARGET_VIS_EXTEND_GROWAROUND  = False  # mode='extend' growaround
TARGET_VIS_EXTEND_FLAGNEARTIME = False  # mode='extend' flagneartime
TARGET_VIS_EXTEND_FLAGNEARFREQ = False  # mode='extend' flagnearfreq

# B0/B bandpass table (CPARAM) rflag/tfcrop thresholds
BPCAL_CPARAM_RFLAG_TIMEDEVSCALE = 5.0   # rflag threshold (sigma), time direction
BPCAL_CPARAM_RFLAG_FREQDEVSCALE = 7.0 if band == 'L' else 5.0   # rflag threshold (sigma), freq direction
BPCAL_CPARAM_TFCROP_TIMECUTOFF  = 4.0   # tfcrop threshold (sigma), time direction
BPCAL_CPARAM_TFCROP_FREQCUTOFF  = 3.0   # tfcrop threshold (sigma), freq direction
BPCAL_CPARAM_TFCROP_MAXNPIECES  = 7     # tfcrop maxnpieces, amplitude pass
BPCAL_CPARAM_EXTENDFLAGS        = False # rflag/tfcrop 'extendflags' on CPARAM

# Extra L-band-only CPARAM rflag over a restricted frequency range, run just
# before the main CPARAM rflag on both the B0 and B tables
BPCAL_CPARAM_LBAND_RFLAG              = False
BPCAL_CPARAM_LBAND_RFLAG_SPW          = '*:0.5GHz~1.6GHz'
BPCAL_CPARAM_LBAND_RFLAG_TIMEDEVSCALE = 5.0   # rflag threshold (sigma), time direction
BPCAL_CPARAM_LBAND_RFLAG_FREQDEVSCALE = 5.0   # rflag threshold (sigma), freq direction

# Extra CPARAM pass: tfcrop on bandpass phase only (ARG_Sol1,Sol2)
BPCAL_EXTRA_TFCROP = False
BPCAL_ARG_TFCROP_MAXNPIECES = 3     # tfcrop maxnpieces, phase-only (ARG) pass
BPCAL_ARG_TFCROP_TIMECUTOFF = 5.0   # tfcrop threshold (sigma), time direction
BPCAL_ARG_TFCROP_FREQCUTOFF = 5.0   # tfcrop threshold (sigma), freq direction

# CPARAM extend pass, propagating flags across polarisation. Independent of the
# phase-only tfcrop above; either pass can run without the other. growtime and
# growfreq at 100 confine the extension to the polarisation axis.
BPCAL_CPARAM_EXTEND      = False  # whether the mode='extend' pass runs
BPCAL_CPARAM_EXTENDPOLS  = True   # mode='extend' extendpols
BPCAL_CPARAM_EXTEND_TIME = 100.0  # mode='extend' growtime (% already flagged)
BPCAL_CPARAM_EXTEND_FREQ = 100.0  # mode='extend' growfreq (% already flagged)

# Calibrator (bpcal/pacal/secondary) visibility rflag/tfcrop thresholds
CAL_VIS_RFLAG_TIMEDEVSCALE = 5.0   # rflag threshold (sigma), time direction
CAL_VIS_RFLAG_FREQDEVSCALE = 5.0   # rflag threshold (sigma), freq direction
CAL_VIS_TFCROP_TIMECUTOFF  = 4.0   # tfcrop threshold (sigma), time direction
CAL_VIS_TFCROP_FREQCUTOFF  = 3.0   # tfcrop threshold (sigma), freq direction
CAL_VIS_TFCROP_MAXNPIECES  = 7     # tfcrop maxnpieces
CAL_VIS_EXTEND             = False  # whether the mode='extend' pass runs
CAL_VIS_EXTENDPOLS         = False  # mode='extend' extendpols
CAL_VIS_EXTEND_TIME        = 90.0  # mode='extend' growtime (% already flagged)
CAL_VIS_EXTEND_FREQ        = 90.0  # mode='extend' growfreq (% already flagged)
CAL_VIS_EXTEND_GROWAROUND  = False # mode='extend' growaround
CAL_VIS_EXTEND_FLAGNEARTIME = False # mode='extend' flagneartime
CAL_VIS_EXTEND_FLAGNEARFREQ = False # mode='extend' flagnearfreq

# Andrew's custom hacky flagger (tools/basedflagger.py) -- extra pass run on the
# primary (after its first-stage flag), the polang cal, each secondary, and
# each target. Off by default.
BASED_FLAGGER = False
BASED_FLAGGER_CORRELATION_PRODUCTS = 'XX,YY,XY,YX'
BASED_FLAGGER_ANTENNA_FLAG_CAP = 2
BASED_FLAGGER_OUTLIER_MODE = 'mixed'

# Test phase as well as amplitude. basedflagger runs the phase pass on the
# parallel hands only, whatever correlations are selected above.
BASED_FLAGGER_PHASE = True

# The --phase switch is appended to every basedflagger call from one place.
BASED_FLAGGER_PHASE_ARG = ' --phase' if BASED_FLAGGER_PHASE else ''

# DEBUGGING: shadems amp/phase vs baseline of the field basedflagger just ran
# on, once its flags are applied -- i.e. the CORRECTED_DATA it tested, before
# any later solve or applycal overwrites it.
DEBUG_BASED_SHADEMS = True

# Which fields to run BASED flagger on: comma-separated combination of
# 'primary', 'pacal', 'secondary', 'target'
BASED_FLAGGER_FIELDS = 'primary,pacal'
BASED_FLAGGER_FIELDS_LIST = [f.strip().lower() for f in BASED_FLAGGER_FIELDS.split(',') if f.strip()]

# Outlier mode for targets specifically, overriding BASED_FLAGGER_OUTLIER_MODE:
# targets are typically fainter fields, where a low outlier is more likely to
# be noise than a defect, so only high outliers are flagged by default
BASED_FLAGGER_TARGET_OUTLIER_MODE = 'high'

def based_flag_fields_arg(fields):
    """
    basedflagger --flag-fields argument: the fields a calibrator's flags are
    passed on to. Each flagged calibrator scan also flags those fields' scans
    between it and that calibrator's previous and next scans.
    """
    return f" --flag-fields '{','.join(fields)}'"

# Primary and polang cal flags go to every field; a secondary's go only to the
# targets paired with it. Targets pass their flags on to nothing.
BASED_FLAGGER_ALL_FIELDS = list(dict.fromkeys(
    [bpcal_name] + ([pacal_name] if pacal_name != '' else []) + pcal_names + targets))
BASED_FLAGGER_ALL_FIELDS_ARG = based_flag_fields_arg(BASED_FLAGGER_ALL_FIELDS)

def apply_based_flags(myms, inpfile, label, returncode=0):
    """
    Apply a basedflagger flag list, reporting the before/after flag fraction.

    A calibrator's flag commands select its own scans plus the scans of the
    other fields its flags are passed on to, so they remove more than the field they were
    derived from. basedflagger can only estimate its own field's share; this
    is where the amount actually removed is measured.

    basedflagger is an optional extra pass, so it never stops the calibration:
    a non-zero `returncode` (the flagger's exit status), a missing list, or a
    list with no commands in it all leave the flags as they are and say why.
    flagdata(mode='list') raises on an empty list, and basedflagger creates its
    list empty before it starts, so an empty list is the normal result for a
    field with nothing to flag as well as what a failed run leaves behind.
    """
    if returncode != 0:
        print(f'WARNING: BASED flagger [{label}]: exited with status {returncode} '
              f'-- no flags applied from {inpfile}')
        return

    if not os.path.isfile(inpfile):
        print(f'WARNING: BASED flagger [{label}]: {inpfile} not found '
              f'-- no flags applied')
        return

    with open(inpfile) as f:
        n_commands = sum(1 for line in f
                         if line.strip() and not line.lstrip().startswith('#'))
    if n_commands == 0:
        print(f'BASED flagger [{label}]: no baselines flagged ({inpfile} is empty)')
        return

    before = flagdata(vis=myms, mode='summary')
    flagdata(vis=myms, mode='list', inpfile=inpfile, flagbackup=False)
    after = flagdata(vis=myms, mode='summary')

    try:
        f0 = before['flagged'] / before['total']
        f1 = after['flagged'] / after['total']
        print(f'BASED flagger [{label}]: flagged {f0:.4%} -> {f1:.4%} '
              f'(+{f1 - f0:.4%}, {int(after["flagged"] - before["flagged"]):,} '
              f'visibilities) from {inpfile}')
    except (KeyError, TypeError, ZeroDivisionError):
        print(f'BASED flagger [{label}]: applied {inpfile} '
              f'(flag summary unavailable)')

def plot_post_based(myms, field, label):
    """
    Plot parallel-hand amp and phase vs baseline for `field` straight after
    its basedflagger flags are applied.
    """
    if not DEBUG_BASED_SHADEMS:
        return
    shadems_cmd = (f"shadems --dir {VISPLOTS} "
                   f"--xaxis BASELINE,BASELINE,BASELINE,BASELINE "
                   f"--yaxis CORRECTED_DATA:amp:XX,CORRECTED_DATA:amp:YY,"
                   f"CORRECTED_DATA:phase:XX,CORRECTED_DATA:phase:YY "
                   f"--colour-by ANTENNA1 --cnum 64 "
                   f"--png 'based_post{label}_{{ms}}_{{field}}_{{label}}.png' "
                   f"--field {field} {myms}")
    print(f'DEBUG: plotting amp/phase vs baseline after BASED flagger [{field}]: {shadems_cmd}')
    subprocess.run([shadems_cmd], shell=True)

def run_based_flagger(field, suffix, description,
                      outlier_mode=BASED_FLAGGER_OUTLIER_MODE, flag_fields_arg=''):
    """
    Run basedflagger on `field` and apply its flag list.

    `suffix` names the outputs (baseline_flags<suffix>.txt and the plots),
    `description` is what the log calls the field, and `flag_fields_arg` the
    --flag-fields argument for the fields its flags are passed on to.
    """
    based_cmd = (f"python_dask {TOOLS}/basedflagger.py "
                 f"{myms} {field} "
                 f"--correlation-products {BASED_FLAGGER_CORRELATION_PRODUCTS} "
                 f"--antenna-flag-cap {BASED_FLAGGER_ANTENNA_FLAG_CAP} --outlier-mode {outlier_mode} "
                 f"--save-dir {VISPLOTS}/ --output-suffix {suffix}"
                 + BASED_FLAGGER_PHASE_ARG + flag_fields_arg)
    print(f'Running BASED flagger on {description}: {based_cmd}')
    based_run = subprocess.run([based_cmd], shell=True)
    apply_based_flags(myms, f'baseline_flags{suffix}.txt', field, based_run.returncode)
    plot_post_based(myms, field, suffix)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------------------------- CALIBRATION PARAMETERS --------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

gapfill = CAL_1GC_FILLGAPS
myuvrange = CAL_1GC_UVRANGE
myspw = CAL_1GC_FREQRANGE

# Optional: uvrange for the primary's solves only (leave as '' to use myuvrange,
# CAL_1GC_UVRANGE in config.py, as every other calibrator does)
# primary_uvrange = '>1000m'
primary_uvrange = ''
primary_uvrange_use = primary_uvrange if primary_uvrange != '' else myuvrange

# Delay (K) handling for the polarisation angle calibrator and the secondaries
# is set by CAL_1GC_DELAY_FROM_PRIMARY in config.py. When on, no K is solved on
# those fields and the primary's delays are applied to them instead, so the K
# tables hold primary solutions only.

# interp for the K table wherever it is applied -- solves and applycal alike,
# including the cross-hand (KCROSS/Xf) path.
DELAY_INTERP_OWN = 'nearest'      # field applying its own delays
DELAY_INTERP_CARRIED = 'linear'   # delays carried over from another field

def delay_interp(applied_to, solutions_from):
    """
    interp for a K table holding `solutions_from`'s delays, applied to
    `applied_to`.

    A field applying its own delays has a solution at its own scans, so it
    snaps to the nearest. Delays carried over from another field have to span
    the gap in time between the two, so they are interpolated.
    """
    return DELAY_INTERP_OWN if solutions_from == applied_to else DELAY_INTERP_CARRIED

def delay_field(field):
    """
    Which field's K solutions to apply to `field`.

    With CAL_1GC_DELAY_FROM_PRIMARY the K tables carry no entry for the secondaries or
    the polang cal, so every field takes the primary's delays.
    """
    return bpcal_name if CAL_1GC_DELAY_FROM_PRIMARY else field


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------------------------------ HELPERS ------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

def tables(*terms):
    """
    The gaintable, gainfield and interp arguments of a CASA solve or applycal,
    from (table, field, interp) terms in apply order. Keeping each table next
    to the field and interp it is applied with stops the three lists drifting
    out of step.
    """
    return dict(gaintable=[t[0] for t in terms],
                gainfield=[t[1] for t in terms],
                interp=[t[2] for t in terms])

def debug_flags():
    """DEBUGGING: print the flag summary of the whole MS (DEBUG_PRINT_FLAGS)."""
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

def flag_bandpass(bptable):
    """Flag outliers in a bandpass table (CPARAM) before it is used."""

    # tfcrop on bandpass amplitude (CPARAM), all correlations
    flagdata(vis=bptable, mode='tfcrop', datacolumn='CPARAM',
        maxnpieces=BPCAL_CPARAM_TFCROP_MAXNPIECES,
        timecutoff=BPCAL_CPARAM_TFCROP_TIMECUTOFF, freqcutoff=BPCAL_CPARAM_TFCROP_FREQCUTOFF,
        extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

    # L-band only: rflag on bandpass amplitude (CPARAM) over a restricted frequency range
    if BPCAL_CPARAM_LBAND_RFLAG and band == 'L':
        flagdata(vis=bptable, mode='rflag', datacolumn='CPARAM', spw=BPCAL_CPARAM_LBAND_RFLAG_SPW,
            timedevscale=BPCAL_CPARAM_LBAND_RFLAG_TIMEDEVSCALE, freqdevscale=BPCAL_CPARAM_LBAND_RFLAG_FREQDEVSCALE,
            extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

    # rflag on bandpass amplitude (CPARAM), all correlations
    flagdata(vis=bptable, mode='rflag', datacolumn='CPARAM',
        timedevscale=BPCAL_CPARAM_RFLAG_TIMEDEVSCALE, freqdevscale=BPCAL_CPARAM_RFLAG_FREQDEVSCALE,
        extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

    # tfcrop on bandpass phase only (ARG_Sol1,Sol2)
    if BPCAL_EXTRA_TFCROP:
        flagdata(vis=bptable, mode='tfcrop', datacolumn='CPARAM', correlation='ARG_Sol1,Sol2',
            maxnpieces=BPCAL_ARG_TFCROP_MAXNPIECES,
            timecutoff=BPCAL_ARG_TFCROP_TIMECUTOFF, freqcutoff=BPCAL_ARG_TFCROP_FREQCUTOFF,
            extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

    # extend on the bandpass table (CPARAM)
    if BPCAL_CPARAM_EXTEND:
        flagdata(vis=bptable, mode='extend', extendpols=BPCAL_CPARAM_EXTENDPOLS,
            growtime=BPCAL_CPARAM_EXTEND_TIME, growfreq=BPCAL_CPARAM_EXTEND_FREQ,
            flagbackup=False)

def flag_calibrator(field, datacolumn):
    """
    rflag, tfcrop and (optionally) extend a calibrator's visibilities.
    `datacolumn` is 'residual' (CORRECTED - MODEL) for the primary, whose model
    is known, and 'corrected' for the 1 Jy placeholder calibrators.
    """
    flagdata(vis=myms,
        mode='rflag',
        datacolumn=datacolumn,
        field=field,
        timedevscale=CAL_VIS_RFLAG_TIMEDEVSCALE,
        freqdevscale=CAL_VIS_RFLAG_FREQDEVSCALE,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn=datacolumn,
        field=field,
        maxnpieces=CAL_VIS_TFCROP_MAXNPIECES,
        timecutoff=CAL_VIS_TFCROP_TIMECUTOFF,
        freqcutoff=CAL_VIS_TFCROP_FREQCUTOFF,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    # extend on this calibrator only (CAL_VIS_EXTEND)
    if CAL_VIS_EXTEND:
        flagdata(vis=myms,
            mode='extend',
            field=field,
            growtime=CAL_VIS_EXTEND_TIME,
            growfreq=CAL_VIS_EXTEND_FREQ,
            growaround=CAL_VIS_EXTEND_GROWAROUND,
            flagneartime=CAL_VIS_EXTEND_FLAGNEARTIME,
            flagnearfreq=CAL_VIS_EXTEND_FLAGNEARFREQ,
            extendpols=CAL_VIS_EXTENDPOLS,
            flagbackup=False)

    debug_flags()

def flag_target(target):
    """rflag, tfcrop, (optionally) extend and BASED-flag a calibrated target."""
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='corrected',
        field=target,
        timedevscale=TARGET_VIS_RFLAG_TIMEDEVSCALE,
        freqdevscale=TARGET_VIS_RFLAG_FREQDEVSCALE,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='corrected',
        field=target,
        timecutoff=TARGET_VIS_TFCROP_TIMECUTOFF,
        freqcutoff=TARGET_VIS_TFCROP_FREQCUTOFF,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    if TARGET_VIS_EXTEND:
        flagdata(vis=myms, mode='extend', growtime=TARGET_VIS_EXTEND_TIME, growfreq=TARGET_VIS_EXTEND_FREQ,
            growaround=TARGET_VIS_EXTEND_GROWAROUND, flagneartime=TARGET_VIS_EXTEND_FLAGNEARTIME,
            flagnearfreq=TARGET_VIS_EXTEND_FLAGNEARFREQ, extendpols=TARGET_VIS_EXTENDPOLS,
            field=target, flagbackup=False)

    if BASED_FLAGGER and 'target' in BASED_FLAGGER_FIELDS_LIST:
        run_based_flagger(target, f'_{target}', f'target {target}',
                          outlier_mode=BASED_FLAGGER_TARGET_OUTLIER_MODE)

def clip_leakage(dftable):
    """Flag leakage (Df) solutions with |D| > 0.1 (10%), far above MeerKAT's typical leakage."""
    flagdata(vis=dftable, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')

def move_aside(table, message):
    """
    Rename `table` to its _bad name so it is never applied, replacing any
    earlier _bad copy. `message` is logged with {} standing for the new name.
    """
    root, ext = os.path.splitext(table)
    bad_table = f'{root}_bad{ext}'
    print('  ' + message.format(bad_table))
    if os.path.isdir(bad_table):
        shutil.rmtree(bad_table)
    shutil.move(table, bad_table)

def save_flags(versionname):
    """Save the current flags as `versionname`, replacing any earlier copy."""
    flagmanager(vis=myms,
        mode='delete',
        versionname=versionname)

    flagmanager(vis=myms,
        mode='save',
        versionname=versionname)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------------------------------- WORKFLOW ----------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

if CAL_1GC_DELAY_FROM_PRIMARY:
    print(f'Delays: solving K on the primary ({bpcal_name}) only; secondaries '
          f'and polang cal take those solutions, applied with '
          f'interp={DELAY_INTERP_CARRIED}')


# ------- Gain table names

def caltable(suffix):
    return f'{GAINTABLES}/cal_1GC_{myms}.{suffix}'

# Stage 0 tables: solved first, used only to flag the calibrators
ktab0, bptab0, gptab0, gtab0, dftab0 = (caltable(s) for s in ('K0', 'B0', 'Gp0', 'G0', 'Df0'))

# Working tables: re-solved on the flagged data, applied to every field.
# ftab is gtab bootstrapped onto the primary's flux scale by fluxscale.
ktab, bptab, gptab, gtab, ftab, dftab = (caltable(s) for s in ('K', 'B', 'Gp', 'G', 'F', 'Df'))
dftab_original = dftab   # dftab is repointed below if Df is re-solved with KCROSS

# Cross-hand tables, and the Df re-solved with KCROSS applied
kcross, xftab, dftab6 = (caltable(s) for s in ('KCROSS', 'Xf', 'Df6'))

# Start from the flags 1GC_04 left (autoflagging of the calibrators' DATA), so a
# re-run of this script never stacks on flags from an earlier run
flagmanager(vis=myms,
        mode='restore',
        versionname='autoflag_cals_data')

# Last-minute manual flags (MANUAL_FLAGS above)
if MANUAL_FLAGS:
    print(f'Applying {len(MANUAL_FLAGS)} manual flag command(s):')
    for cmd in MANUAL_FLAGS:
        print(f'  {cmd}')
    flagdata(vis=myms, mode='list', inpfile=list(MANUAL_FLAGS), flagbackup=False)

# fluxscale fails if the F table already exists, so remove one left by an earlier run
if os.path.isdir(ftab):
    print(f"Removing: {ftab}")
    shutil.rmtree(ftab)


# ------- Calibrator models

# The primary's model goes into MODEL_DATA, which may not exist yet
tb.open(myms)
has_model_data = 'MODEL_DATA' in tb.colnames()
tb.close()
if not has_model_data:
    # Dummy setjy call to initialize non-existing model data column
    setjy(vis=myms,
        standard='manual',
        field=bpcal_name,
        fluxdensity=[1.0, 0, 0, 0],
        reffreq='1000MHz',
        usescratch=True)

# MeerKAT specific crystalball models from B. Hugo, for L and UHF band:
#   PKS 1934-638: https://archive-gw-1.kat.ac.za/public/repository/10.48479/hhhy-4r55/index.html
#   PKS 0408-65:  https://archive-gw-1.kat.ac.za/public/repository/10.48479/ez63-vx81/index.html
CRYSTALBALL_MODELS = {
    ('1934', 'L'):   'fitted.PKS1934.LBand.wsclean.cat.txt',
    ('1934', 'UHF'): 'fitted.PKS1934.UBand.wsclean.cat.txt',
    ('0408', 'L'):   'fitted.PKS0407.LBand.wsclean.cat.txt',
    ('0408', 'UHF'): 'fitted.PKS0407.UBand.wsclean.cat.txt',
}

if (primary_tag, band) in CRYSTALBALL_MODELS:
    syscall = (f"crystalball {myms} -f {bpcal_name} "
               f"-sm {DATA}/crystalball/{CRYSTALBALL_MODELS[(primary_tag, band)]}")
    subprocess.run([syscall],shell=True)

elif primary_tag == '1934':
    setjy(vis=myms,
        field=bpcal_name,
        standard='Stevens-Reynolds 2016',
        scalebychan=True,
        usescratch=True)

elif primary_tag == '0408':
    # Recommendation from SARAO, i.e., https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/1481408634/Flux+and+bandpass+calibration
    # See script in tools/SARAO_0408_model.py to see where CASA parameters come from
    setjy(vis=myms,
        field=bpcal_name,
        standard='manual',
        fluxdensity=[6.9862, 0.0, 0.0, 0.0],
        spix=[-1.2897, -0.2353, 0.0861],
        reffreq='2.7GHz',
        scalebychan=True,
        usescratch=True)

elif primary_tag == 'other':
    setjy(vis=myms,
        field=bpcal_name,
        standard='Perley-Butler 2013',
        scalebychan=True,
        usescratch=True)

# The secondaries (and, for the parallel-hand solves, the polang cal) get a
# 1 Jy point-source placeholder; fluxscale puts their gains on the primary's
# scale. The polang cal's real model (POLANG_MOD) is set before the cross-hand
# solves further down.
for pcal in pcal_names:
    if pcal != bpcal_name:
        setjy(vis =myms,
            field = pcal,
            standard = 'manual',
            fluxdensity = [1.0,0,0,0],
            reffreq = '1000MHz',
            usescratch = True)

if pacal_name != '' and pacal_name != bpcal_name and pacal_name not in pcal_names:
    setjy(vis=myms,
        field=pacal_name,
        standard='manual',
        fluxdensity=[1.0,0,0,0],
        reffreq='1000MHz',
        usescratch=True)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# --------------------------- PRIMARY (K, B, Df) ----------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

def solve_primary(k, gp, bp, g, df, df_append):
    """
    The full chain on the primary: K, then Gp (applying K), B (K, Gp),
    G (K, Gp, B) and Df (K, Gp, B, G). The bandpass is flagged before G is
    solved against it, and Df solutions above 0.1 are clipped.
    """

    # ------- K (primary; no prior calibration). Delays use the full band and
    # all baselines, so the uvrange and spw selections are left off
    gaincal(vis=myms,
        field=bpcal_name,
        caltable=k,
        # uvrange=primary_uvrange_use,
        # spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf')

    k_term = (k, bpcal_name, delay_interp(bpcal_name, bpcal_name))

    # ------- Gp (primary; apply K) -- Type G, phase-only
    gaincal(vis=myms,
        field=bpcal_name,
        uvrange=primary_uvrange_use,
        spw=myspw,
        caltable=gp,
        refant=str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        **tables(k_term))

    # ------- B (primary; apply K, Gp)
    bandpass(vis=myms,
        field=bpcal_name,
        uvrange=primary_uvrange_use,
        caltable=bp,
        refant=str(ref_ant),
        solint='inf',
        combine='',
        solnorm=False,
        minblperant=4,
        minsnr=3.0,
        bandtype='B',
        fillgaps=gapfill,
        **tables(k_term, (gp, bpcal_name, 'linear')))

    flag_bandpass(bp)
    debug_flags()

    # ------- G (primary; apply K, Gp, B) -- Type T, amp-only
    gaincal(vis=myms,
        field=bpcal_name,
        uvrange=primary_uvrange_use,
        spw=myspw,
        caltable=g,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        **tables(k_term, (gp, bpcal_name, 'linear'), (bp, bpcal_name, 'linear')))

    # ------- Df (primary; apply K, Gp, B, G)
    # The primary is unpolarised, so all of its cross-hand signal is leakage
    polcal(vis=myms,
        field=bpcal_name,
        uvrange=primary_uvrange_use,
        caltable=df,
        refant=str(ref_ant),
        solint='inf',
        poltype='Df',
        combine='scan',
        **tables(k_term, (gp, bpcal_name, 'linear'), (bp, bpcal_name, 'linear'),
                 (g, bpcal_name, 'linear')),
        **({'append': False} if df_append is False else {}))

    clip_leakage(df)
    debug_flags()


# --------------------------- Stage 0 ------------------------------------- #
# Solve on the primary as delivered, apply, and flag its residuals against the
# model. These tables are only used for this flagging.

solve_primary(ktab0, gptab0, bptab0, gtab0, dftab0, df_append=False)

# ------- Correct primary data with K0, Gp0, B0, G0, Df0
applycal(vis=myms,
    field=bpcal_name,
    parang=False,
    flagbackup=False,
    **tables((ktab0, bpcal_name, delay_interp(bpcal_name, bpcal_name)),
             (gptab0, bpcal_name, 'linear'),
             (bptab0, bpcal_name, 'linear'),
             (gtab0, bpcal_name, 'linear'),
             (dftab0, bpcal_name, 'linear')))

debug_flags()

# ------- Flag primary on CORRECTED_DATA - MODEL_DATA
flag_calibrator(bpcal_name, 'residual')
save_flags('bpcal_residual_flags')

if BASED_FLAGGER and 'primary' in BASED_FLAGGER_FIELDS_LIST:
    run_based_flagger(bpcal_name, '_bpcal', 'primary',
                      flag_fields_arg=BASED_FLAGGER_ALL_FIELDS_ARG)


# --------------------------- Working tables ------------------------------ #
# Re-solve the full chain on the flagged primary. K, B and Df from here are
# applied to every field.

solve_primary(ktab, gptab, bptab, gtab, dftab, df_append=None)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------- SECONDARIES + POLANG CAL (Gp, G, K) ----------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# Every calibrator other than the primary: the polang cal first, then the
# secondaries (a secondary that is also the primary or the polang cal is only
# solved once).
cal_fields = (([pacal_name] if pacal_name != '' else [])
              + [pcal for pcal in pcal_names if pcal != bpcal_name and pcal != pacal_name])

def solve_calibrator(field, gp, g, k):
    """
    Phase (Gp), then amplitude (G), then delay (K) solutions on a secondary or
    the polang cal, against the primary's working K, B and Df. All three are
    appended to tables already holding the primary's solutions. K is skipped
    with CAL_1GC_DELAY_FROM_PRIMARY.
    """
    from_primary = [(ktab, bpcal_name, delay_interp(field, bpcal_name)),
                    (bptab, bpcal_name, 'linear'),
                    (dftab, bpcal_name, 'linear')]

    # ------- Gp (apply K, B, Df from primary) -- Type G, phase-only
    gaincal(vis=myms,
        field=field,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gp,
        refant=str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        append=True,
        **tables(*from_primary))

    # ------- G (apply K, B, Df from primary, Gp from this field) -- Type T, amp-only
    gaincal(vis=myms,
        field=field,
        uvrange=myuvrange,
        spw=myspw,
        caltable=g,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        append=True,
        **tables(*from_primary, (gp, field, 'linear')))

    # ------- K (apply B, Df from primary, Gp, G from this field), full band and
    # all baselines as for the primary
    if not CAL_1GC_DELAY_FROM_PRIMARY:
        gaincal(vis=myms,
            field=field,
            caltable=k,
            # uvrange=myuvrange,
            # spw=myspw,
            refant=str(ref_ant),
            gaintype='K',
            solint='inf',
            append=True,
            **tables((bptab, bpcal_name, 'linear'),
                     (dftab, bpcal_name, 'linear'),
                     (gp, field, 'linear'),
                     (g, amp_field(field), 'linear')))


# --------------------------- Stage 0 ------------------------------------- #
# Solve, apply and flag each calibrator. These tables are only used for the
# flagging.

for field in cal_fields:
    solve_calibrator(field, gptab0, gtab0, ktab0)

for field in cal_fields:

    # -------- Applycal (B, Df from primary, K0, Gp0, G0 from this field) and Flag
    applycal(vis=myms,
        field=field,
        parang=False,
        flagbackup=False,
        **tables((ktab0, delay_field(field), delay_interp(field, delay_field(field))),
                 (bptab, bpcal_name, 'linear'),
                 (gptab0, field, 'linear'),
                 (gtab0, amp_field(field), 'linear'),
                 (dftab, bpcal_name, 'linear')))

    debug_flags()
    flag_calibrator(field, 'corrected')

    # The polang cal's flags are passed on to every field; a secondary's only
    # to the targets paired with it
    if field == pacal_name:
        if BASED_FLAGGER and 'pacal' in BASED_FLAGGER_FIELDS_LIST:
            run_based_flagger(field, '_pacal', 'pacal',
                              flag_fields_arg=BASED_FLAGGER_ALL_FIELDS_ARG)
    elif BASED_FLAGGER and 'secondary' in BASED_FLAGGER_FIELDS_LIST:
        paired_targets = [t for t, c in zip(targets, target_cal_map) if c == field]
        run_based_flagger(field, f'_{field}', f'secondary {field}',
                          flag_fields_arg=based_flag_fields_arg(paired_targets))


# --------------------------- Working tables ------------------------------ #
# Re-solve on the flagged calibrators. These are the solutions applied.

for field in cal_fields:
    solve_calibrator(field, gptab, gtab, ktab)

# ------- Flux scale: bootstrap every calibrator's G onto the primary's scale.
# With the primary as the only calibrator there is nothing to scale, so G is
# applied as it is.
if cal_fields:
    fluxscale(vis=myms,
        caltable=gtab,
        fluxtable=ftab,
        reference=bpcal_name,
        append=False,
        transfer='')
else:
    ftab = gtab


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------- CROSS-HAND (KCROSS, Xf; POLANG CAL) ----------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

def solve_kcross(scan=None):
    """
    KCROSS (cross-hand delay) on the polang cal, applying its K, Gp, G and the
    primary's B. `scan` restricts the solve to a comma-separated scan list.
    """
    gaincal(vis = myms,
        field = pacal_name,
        **({'scan': scan} if scan is not None else {}),
        spw = kcross_spw,
        caltable = kcross,
        refant = str(ref_ant),
        solint = 'inf',
        combine = combine_param,
        gaintype='KCROSS',
        parang = True,
        append = False,
        **tables((ktab, delay_field(pacal_name), delay_interp(pacal_name, delay_field(pacal_name))),
                 (gptab, pacal_name, 'linear'),
                 (bptab, bpcal_name, 'linear'),
                 (gtab, amp_field(pacal_name), 'linear')))

def solve_df_with_kcross(df):
    """
    Re-solve Df on the primary with KCROSS applied, into `df`.

    CASA enforces the apply order KCROSS -> Df -> Xf, and these tables do not
    commute, so Df must be re-solved whenever KCROSS is applied; otherwise the
    original (KCROSS-free) Df is incoherent with the applied KCROSS/Xf terms.
    """
    polcal(vis=myms,
        field=bpcal_name,
        uvrange=primary_uvrange_use,
        caltable=df,
        refant=str(ref_ant),
        solint='inf',
        poltype='Df',
        combine='scan',
        **tables((ktab, bpcal_name, delay_interp(bpcal_name, bpcal_name)),
                 (gptab, bpcal_name, 'linear'),
                 (bptab, bpcal_name, 'linear'),
                 (gtab, bpcal_name, 'linear'),
                 (kcross, pacal_name, 'linear')))

    clip_leakage(df)

def solve_xf(xf_terms, scan=None):
    """Xf (cross-hand phase) on the polang cal; `scan` restricts the solve to a scan list."""
    polcal(vis = myms,
        field = pacal_name,
        **({'scan': scan} if scan is not None else {}),
        uvrange = myuvrange,
        caltable = xftab,
        refant = str(ref_ant),
        solint = f'inf,{XF_CHANINT}ch',
        poltype='Xf',
        combine = combine_param,
        append = False,
        **tables(*xf_terms))


if pacal_name != '':

    # ------- Polang cal model: replaces the 1 Jy placeholder used for its
    # parallel-hand solves with the polarised model (POLANG_MOD in config.py)
    setjy(vis=myms,
        field=pacal_name,
        standard='manual',
        fluxdensity = POLANG_MOD,
        usescratch=True)

    # -------- Cross-hand terms. XF_MODE picks the solver: 'casa' (KCROSS + Xf
    # here), 'manual' (tools/manual_XF_solver.py), or 'auto' (CASA, falling back
    # to manual if no scan passes the phase-continuity check below)
    if XF_MODE not in ['casa', 'manual', 'auto']:
        print(f"Cross-hand phase mode ({XF_MODE}) not a valid option, defauling to 'auto'")
        XF_MODE = 'auto'

    manual_XF = False

    # Force manual mode for UHF band
    if band == 'UHF':
        print("UHF band detected, forcing manual XF mode")
        manual_XF = True
        XF_MODE = 'manual'

    if XF_MODE == 'casa' or XF_MODE == 'auto':

        print(f"Cross-hand phase mode is {XF_MODE}; solving CASA based solutions")

        combine_param = 'scan' if XF_AVG_SCAN else ''
        kcross_spw = '*:1.0~1.65GHz' if band == 'L' else ''

        # ------- KCROSS (polang cal; apply its K, Gp, G and the primary's B)
        # Always solve KCROSS (cheap computation), skip flag only controls application
        print(f"  Solving for KCROSS (cross-hand delay) {'[will NOT be applied]' if XF_SKIP_KCROSS else '[will be applied]'}")
        solve_kcross()

        if XF_SKIP_KCROSS:
            # Moved aside straight away, so it is never applied
            move_aside(kcross, 'Moving KCROSS to: {} (will not be used)')
        else:
            print(f"  Re-solving Df (leakage) with KCROSS included: {dftab6}")
            solve_df_with_kcross(dftab6)
            dftab = dftab6   # active Df is now the KCROSS-informed table

        # -------- Xf (polang cal; apply its K, Gp, G, the primary's B, Df, and KCROSS if used)
        print(f"  Solving for Xf (cross-hand phase) {'WITHOUT' if XF_SKIP_KCROSS else 'with'} KCROSS")
        xf_terms = [(ktab, delay_field(pacal_name), delay_interp(pacal_name, delay_field(pacal_name))),
                    (gptab, pacal_name, 'linear'),
                    (bptab, bpcal_name, 'linear'),
                    (gtab, amp_field(pacal_name), 'linear'),
                    (dftab, bpcal_name, 'linear')]
        if not XF_SKIP_KCROSS:
            xf_terms.append((kcross, pacal_name, 'linear'))
        solve_xf(xf_terms)

        # Cross-hand tables, fields and interp for the final applycal
        cross_table = [xftab] if XF_SKIP_KCROSS else [kcross, xftab]
        cross_field = [pacal_name] if XF_SKIP_KCROSS else [pacal_name, pacal_name]
        cross_interp = ['linear'] if XF_SKIP_KCROSS else ['linear', 'linear']

        # ------- Continuity check: a cross-hand phase that jumps by more than
        # XF_AUTO_ANG_JUMP between channels within a scan has gone wrong
        # (usually a 180 deg ambiguity), so that scan is not trusted.
        tb.open(xftab)
        gains = tb.getcol('CPARAM')
        flags = tb.getcol('FLAG')
        scans = tb.getcol('SCAN_NUMBER')
        tb.close()

        print(f"Gains shape: {gains.shape}")
        print(f"Unique scans: {np.unique(scans)}")

        unique_scans = np.unique(scans)
        manual_XF_per_scan = {}

        for scan in unique_scans:
            scan_mask = (scans == scan)

            # Median over time within the scan, unflagged channels only
            scan_gains = np.nanmedian(gains[0, :, scan_mask].T, axis=-1)
            scan_flags = np.nanmedian(flags[0, :, scan_mask].T, axis=-1).astype(bool)
            scan_gains = scan_gains[~scan_flags]

            if len(scan_gains) < 2:
                manual_XF_per_scan[scan] = True
                print(f"Scan {scan}: fewer than 2 unflagged channels, continuous = False")
                continue

            # Largest channel-to-channel phase jump, wrapped to +-180 deg
            phase_diffs = np.diff(np.angle(scan_gains))
            phase_diffs = np.arctan2(np.sin(phase_diffs), np.cos(phase_diffs))
            max_jump = np.max(np.abs(np.degrees(phase_diffs)))

            is_continuous = max_jump < XF_AUTO_ANG_JUMP
            manual_XF_per_scan[scan] = not is_continuous

            print(f"Scan {scan}: max jump = {max_jump:.2f}°, continuous = {is_continuous}")

        continuous_scans = [scan for scan in unique_scans if not manual_XF_per_scan[scan]]
        discontinuous_scans = [scan for scan in unique_scans if manual_XF_per_scan[scan]]

        print("\n" + "="*60)

        if len(continuous_scans) == len(unique_scans):
            # Case 1: All scans are continuous
            print(f"✓ All {len(unique_scans)} scan(s) are phase continuous")
            print("  Proceeding with standard XF calibration")

        elif len(continuous_scans) > 0:
            # Case 2: Some scans are continuous -- re-solve on those only
            print(f"✓ {len(continuous_scans)} scan(s) are phase continuous: {continuous_scans}")
            print(f"✗ {len(discontinuous_scans)} scan(s) failed continuity check: {discontinuous_scans}")
            good_scans = ','.join(map(str, continuous_scans))

            move_aside(xftab, 'Moving original Xf table to: {}')

            # Only handle KCROSS if it's being used (not already moved to _bad)
            if not XF_SKIP_KCROSS:
                move_aside(kcross, 'Moving original KCROSS table to: {}')

                print(f"  Remaking KCROSS table using only continuous scans: {continuous_scans}")
                solve_kcross(scan=good_scans)

                # dftab already points at dftab6 from the first KCROSS-informed solve
                # above -- overwrite it again with the continuous-scans-only KCROSS.
                print(f"  Re-solving Df (leakage) with KCROSS included, overwriting: {dftab}")
                solve_df_with_kcross(dftab)

            print(f"  Remaking Xf table using only continuous scans: {continuous_scans}")
            solve_xf(xf_terms, scan=good_scans)

        else:
            # Case 3: No scans are continuous -- hand over to the manual solver
            print(f"✗ None of the {len(unique_scans)} scan(s) are phase continuous")

            move_aside(xftab, 'Moving original Xf table to: {}')

            # Only move KCROSS if it's being used (not already moved to _bad)
            if not XF_SKIP_KCROSS:
                move_aside(kcross, 'Moving original KCROSS table to: {}')

                print(f"  Reverting Df to original table (KCROSS abandoned): {dftab_original}")
                dftab = dftab_original

            print(f"  Falling back to manual XF solver ({len(unique_scans)} scan(s))")
            manual_XF = True

        print("="*60)

    if XF_MODE == 'manual' or manual_XF:
        print(f"Cross-hand phase mode is {XF_MODE} (or auto detected a large phase jump? {manual_XF}); solving manual solutions")
        # Run in this namespace: the solver reads (and may reset) its tables here
        exec(open('tools/manual_XF_solver.py').read())

        # The manual solver writes a single Xf table, no KCROSS
        cross_table = [xftab]
        cross_field = [pacal_name]
        cross_interp = ['linear']

    # --- Optional: replace the active Xf table with an external one ---
    # When overriding, drop KCROSS from the gain tables since the external
    # table was not solved against it.
    if str(XF_OVERRIDE_TABLE) != '':
        print(f'  [XF] Overriding active Xf table with external: {XF_OVERRIDE_TABLE}')
        xftab        = XF_OVERRIDE_TABLE
        cross_table  = [xftab]
        cross_field  = ['']
        cross_interp = ['nearest,linear']
        dftab        = dftab_original
        print(f'  [XF] KCROSS removed from gain tables; cross_table = {cross_table}')

    # --- Optional: resolve the +/-180 deg degeneracy to a reference phase ---
    # Applied to whichever Xf table is active (solved or overridden above).
    if str(XF_REFERENCE_PHASE) != '':
        _ref_rad = np.deg2rad(float(XF_REFERENCE_PHASE))
        _n_flip  = 0
        tb.open(xftab, nomodify=False)
        _cparam = tb.getcol('CPARAM')
        _flags  = tb.getcol('FLAG')
        for _corr in range(_cparam.shape[0]):
            _ph      = np.angle(_cparam[_corr])
            _diff    = np.abs(np.angle(np.exp(1j * (_ph - _ref_rad))))
            _diff_fl = np.abs(np.angle(np.exp(1j * (_ph + np.pi - _ref_rad))))
            _flip    = (_diff_fl < _diff) & (~_flags[_corr])
            _cparam[_corr][_flip] *= -1
            _n_flip += int(np.sum(_flip))
        tb.putcol('CPARAM', _cparam)
        tb.close()
        print(f'  [XF] Phase degeneracy resolved to reference {XF_REFERENCE_PHASE} deg: '
              f'{_n_flip} channel(s) flipped.')

    # CAL_1GC_APPLYPARANG_CAL controls parang for the polang cal and the
    # secondaries, CAL_1GC_APPLYPARANG for the targets
    parang_cal = CAL_1GC_APPLYPARANG_CAL
    parang_target = CAL_1GC_APPLYPARANG

else:

    # --- No polang cal, so no Xf is solved. Optionally apply an external Xf
    # table to the secondaries and targets (never the primary). Without one,
    # nothing rotates the data into the sky frame, so parang is not applied.
    cross_table = []
    cross_field = []
    cross_interp = []
    parang_cal = False
    parang_target = False
    if str(XF_OVERRIDE_TABLE) != '':
        print(f'  [XF] No polarization angle calibrator, but XF_OVERRIDE_TABLE is set -- '
              f'applying external Xf table to secondaries and targets only: {XF_OVERRIDE_TABLE}')
        cross_table = [XF_OVERRIDE_TABLE]
        cross_field = ['']
        cross_interp = ['nearest,linear']
        parang_cal = CAL_1GC_APPLYPARANG_CAL
        parang_target = CAL_1GC_APPLYPARANG


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------- APPLYCAL (ALL FIELDS) ------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

def apply_all(field, phase_from, amp_from, parang):
    """
    Apply the working tables to `field`, taking phases from `phase_from` and
    amplitudes from `amp_from`; delays follow `phase_from` unless
    CAL_1GC_DELAY_FROM_PRIMARY. B and Df always come from the primary, the
    cross-hand terms from the polang cal.
    """
    applycal(vis=myms,
        field=field,
        parang=parang,
        flagbackup=False,
        gaintable=[ktab, gptab, bptab, ftab, dftab] + cross_table,
        gainfield=[delay_field(phase_from), phase_from, bpcal_name, amp_from, bpcal_name] + cross_field,
        interp=[delay_interp(field, delay_field(phase_from)), 'linear', 'linear', 'linear', 'linear'] + cross_interp)

# ------- Primary: its own tables, the original Df, no cross-hand terms
applycal(vis=myms,
    field=bpcal_name,
    parang=False,
    flagbackup=False,
    **tables((ktab, bpcal_name, delay_interp(bpcal_name, bpcal_name)),
             (gptab, bpcal_name, 'linear'),
             (bptab, bpcal_name, 'linear'),
             (ftab, bpcal_name, 'linear'),
             (dftab_original, bpcal_name, 'linear')))

# ------- Polang cal and secondaries: their own solutions
if pacal_name != '':
    apply_all(pacal_name, pacal_name, amp_field(pacal_name), parang_cal)

for pcal in pcal_names:
    if pcal == bpcal_name or pcal == pacal_name:
        continue   # already calibrated above
    apply_all(pcal, pcal, amp_field(pcal), parang_cal)

# ------- Targets: delays and phases from their paired secondary (target_cal_map),
# amplitudes per target_amp_cal_map, then target flagging
for i, target in enumerate(targets):
    apply_all(target, target_cal_map[i], target_amp_cal_map[i], parang_target)
    flag_target(target)

# ---- Aggressive flagging of CAL_1GC_BL_FREQS (CAL_1GC_AGGRESSIVE_FLAGS in config.py)
if CAL_1GC_AGGRESSIVE_FLAGS and CAL_1GC_BL_FREQS != []:
    flagdata(vis = myms,
        mode = 'manual',
        spw = ','.join(CAL_1GC_BL_FREQS),
        flagbackup=False)

save_flags('1GC_flags')

if pacal_name == '':
    sys.exit('Ending Early! No polarization angle calibrator')
