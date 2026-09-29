# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import glob
import json
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

def stamp():
    now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
    return now
    

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


gapfill = CAL_1GC_FILLGAPS
myuvrange = CAL_1GC_UVRANGE 
myspw = CAL_1GC_FREQRANGE

# Optional: Override uvrange for primary calibrator only (leave as '' to use myuvrange)
# primary_uvrange = '>1000m'
primary_uvrange = ''

# Use primary_uvrange if set, otherwise use myuvrange
if primary_uvrange != '':
    primary_uvrange_use = primary_uvrange
else:
    primary_uvrange_use = myuvrange

# --- Global cross-hand phase term
ALT_PACAL = ''         # Field name to use for Xf applycal when pacal_name is empty
ALT_TYPE  = 'combine'  # Which global_crosshand_phase table to apply: 'combine' or 'perscan'

# Validate
if ALT_TYPE not in ('combine', 'perscan'):
    raise ValueError(f"ALT_TYPE='{ALT_TYPE}' is not valid; must be 'combine' or 'perscan'.")

if ALT_PACAL != '':
    msmd.open(myms)
    _valid_fields = msmd.fieldnames()
    msmd.close()
    if ALT_PACAL not in _valid_fields:
        raise ValueError(
            f"ALT_PACAL='{ALT_PACAL}' not found in MS field list: {_valid_fields}")


# ------- Setup names


tt = stamp()

# Calibrator tables
ktab0 = GAINTABLES+'/cal_1GC_'+myms+'.K0'
bptab0 = GAINTABLES+'/cal_1GC_'+myms+'.B0'
gptab0 = GAINTABLES+'/cal_1GC_'+myms+'.Gp0'
gtab0 = GAINTABLES+'/cal_1GC_'+myms+'.G0'
dftab0  = GAINTABLES+'/cal_1GC_'+myms+'.Df0'

ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gptab = GAINTABLES+'/cal_1GC_'+myms+'.Gp'
gtab = GAINTABLES+'/cal_1GC_'+myms+'.G'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'
dftab_original = dftab

kcross  = GAINTABLES+'/cal_1GC_'+myms+'.KCROSS'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'
dftab6 = GAINTABLES+'/cal_1GC_'+myms+'.Df6'


# Restore the auto_cal flag version
flagmanager(vis=myms,
    mode='restore',
    versionname='autoflag_cals_data')

# Remove ftabs if they exist to prevent code from breaking
if os.path.isdir(ftab):
    print(f"Removing: {ftab}")
    shutil.rmtree(ftab)

# ------- Set BP calibrator models

# Check if MODEL_DATA column exists, if not initialize it
tb.open(myms)
if 'MODEL_DATA' not in tb.colnames():
    tb.close()
    # Dummy setjy call to initialize non-existing model data column
    setjy(vis=myms,
        standard='manual',
        field=bpcal_name,
        fluxdensity=[1.0, 0, 0, 0],
        reffreq='1000MHz',
        usescratch=True)
else:
    tb.close()

if primary_tag == '1934':
    
    # MeerKAT specific crystalball models for 1939 from B.Hugo: https://archive-gw-1.kat.ac.za/public/repository/10.48479/hhhy-4r55/index.htmlV
    if band == 'L':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS1934.LBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)

    elif band == 'UHF':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS1934.UBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)

    else:
        setjy(vis=myms,
            field=bpcal_name,
            standard='Stevens-Reynolds 2016',
            scalebychan=True,
            usescratch=True)

elif primary_tag == '0408':

    # MeerKAT specific crystalball models for 0408 from B.Hugo: https://archive-gw-1.kat.ac.za/public/repository/10.48479/ez63-vx81/index.html
    if band == 'L':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS0407.LBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)

    elif band == 'UHF':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS0407.UBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)
    
    else:
        # OXKAT Version that uses config.py input
        #bpcal_mod = CAL_1GC_0408_MODEL
        #setjy(vis=myms,
        #    field=bpcal_name,
        #    standard='manual',
        #    fluxdensity=bpcal_mod[0],
        #    spix=bpcal_mod[1],
        #    reffreq=bpcal_mod[2],
        #    scalebychan=True,
        #    usescratch=True)

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

for i in range(0,len(pcal_names)):
    pcal = pcal_names[i]
    if pcal != bpcal_name:
        setjy(vis =myms,
            field = pcal,
            standard = 'manual',
            fluxdensity = [1.0,0,0,0],
            reffreq = '1000MHz',
            usescratch = True)

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------- STAGE 0 ----------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #


# ------- K0 (primary)

gaincal(vis=myms,
    field=bpcal_name,
    caltable=ktab0,
    uvrange=primary_uvrange_use,
    # spw=myspw,
    refant=str(ref_ant),
    gaintype='K',
    solint='inf')


# ------- Gp0 (primary; apply K0) -- Type G, phase-only

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    spw=myspw,
    caltable=gptab0,
    refant=str(ref_ant),
    gaintype='G',
    solint='inf',
    calmode='p',
    minsnr=3,
    gainfield=[bpcal_name],
    interp=['nearest'],
    gaintable=[ktab0])


# ------- B0 (primary; apply K0, Gp0)

bandpass(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=bptab0,
    refant=str(ref_ant),
    solint='inf',
    combine='',
    solnorm=False,
    minblperant=4,
    minsnr=3.0,
    bandtype='B',
    fillgaps=gapfill,
    gainfield=[bpcal_name, bpcal_name],
    interp=['nearest', 'linear'],
    gaintable=[ktab0, gptab0])

# tfcrop on bandpass amplitude (CPARAM), all correlations
flagdata(vis=bptab0, mode='tfcrop', datacolumn='CPARAM',
    maxnpieces=BPCAL_CPARAM_TFCROP_MAXNPIECES,
    timecutoff=BPCAL_CPARAM_TFCROP_TIMECUTOFF, freqcutoff=BPCAL_CPARAM_TFCROP_FREQCUTOFF,
    extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# L-band only: rflag on bandpass amplitude (CPARAM) over a restricted frequency range
if BPCAL_CPARAM_LBAND_RFLAG and band == 'L':
    flagdata(vis=bptab0, mode='rflag', datacolumn='CPARAM', spw=BPCAL_CPARAM_LBAND_RFLAG_SPW,
        timedevscale=BPCAL_CPARAM_LBAND_RFLAG_TIMEDEVSCALE, freqdevscale=BPCAL_CPARAM_LBAND_RFLAG_FREQDEVSCALE,
        extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# rflag on bandpass amplitude (CPARAM), all correlations
flagdata(vis=bptab0, mode='rflag', datacolumn='CPARAM',
    timedevscale=BPCAL_CPARAM_RFLAG_TIMEDEVSCALE, freqdevscale=BPCAL_CPARAM_RFLAG_FREQDEVSCALE,
    extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# tfcrop on bandpass phase only (ARG_Sol1,Sol2)
if BPCAL_EXTRA_TFCROP:
    flagdata(vis=bptab0, mode='tfcrop', datacolumn='CPARAM', correlation='ARG_Sol1,Sol2',
        maxnpieces=BPCAL_ARG_TFCROP_MAXNPIECES,
        timecutoff=BPCAL_ARG_TFCROP_TIMECUTOFF, freqcutoff=BPCAL_ARG_TFCROP_FREQCUTOFF,
        extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# extend on the bandpass table (CPARAM)
if BPCAL_CPARAM_EXTEND:
    flagdata(vis=bptab0, mode='extend', extendpols=BPCAL_CPARAM_EXTENDPOLS,
        growtime=BPCAL_CPARAM_EXTEND_TIME, growfreq=BPCAL_CPARAM_EXTEND_FREQ,
        flagbackup=False)

# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

# ------- G0 (primary; apply K0, Gp0, B0) -- Type T, amp-only

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    spw=myspw,
    caltable=gtab0,
    refant=str(ref_ant),
    gaintype='T',
    solint='inf',
    calmode='a',
    minsnr=3,
    gainfield=[bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear'],
    gaintable=[ktab0, gptab0, bptab0])


# -------- Solve for Df0 (apply K0, Gp0, B0, G0)

polcal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=dftab0,
    refant=str(ref_ant),
    solint='inf',
    poltype='Df',
    combine='scan',
    gaintable=[ktab0, gptab0, bptab0, gtab0],
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear', 'linear'],
    append=False)

flagdata(vis=dftab0, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

# ------- Correct primary data with K0, Gp0, B0, G0, Df0

applycal(vis=myms,
    gaintable=[ktab0, gptab0, bptab0, gtab0, dftab0],
    field=bpcal_name,
    parang=False,
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear', 'linear', 'linear'],
    flagbackup=False)
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')


# ------- Flag primary on CORRECTED_DATA - MODEL_DATA

# rflag on CORRECTED_DATA - MODEL_DATA
flagdata(vis=myms,
    mode='rflag',
    datacolumn='residual',
    field=bpcal_name,
    timedevscale=CAL_VIS_RFLAG_TIMEDEVSCALE,
    freqdevscale=CAL_VIS_RFLAG_FREQDEVSCALE,
    extendflags=EXTEND_AUTO,
    flagbackup=False)

# tfcrop on CORRECTED_DATA - MODEL_DATA
flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='residual',
    field=bpcal_name,
    maxnpieces=CAL_VIS_TFCROP_MAXNPIECES,
    timecutoff=CAL_VIS_TFCROP_TIMECUTOFF,
    freqcutoff=CAL_VIS_TFCROP_FREQCUTOFF,
    extendflags=EXTEND_AUTO,
    flagbackup=False)

# extend on bpcal
if CAL_VIS_EXTEND:
    flagdata(vis=myms,
        mode='extend',
        growtime=CAL_VIS_EXTEND_TIME,
        growfreq=CAL_VIS_EXTEND_FREQ,
        growaround=CAL_VIS_EXTEND_GROWAROUND,
        flagneartime=CAL_VIS_EXTEND_FLAGNEARTIME,
        flagnearfreq=CAL_VIS_EXTEND_FLAGNEARFREQ,
        extendpols=CAL_VIS_EXTENDPOLS)

# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

flagmanager(vis=myms,
        mode='delete',
        versionname='bpcal_residual_flags')

flagmanager(vis=myms,
        mode='save',
        versionname='bpcal_residual_flags')


# ------- BASED flagger (primary)


if BASED_FLAGGER and 'primary' in BASED_FLAGGER_FIELDS_LIST:
    based_cmd = (f"python_dask {TOOLS}/basedflagger.py "
                 f"{myms} {bpcal_name} "
                 f"--correlation-products {BASED_FLAGGER_CORRELATION_PRODUCTS} "
                 f"--antenna-flag-cap {BASED_FLAGGER_ANTENNA_FLAG_CAP} --outlier-mode {BASED_FLAGGER_OUTLIER_MODE} "
                 f"--save-dir {VISPLOTS}/ --output-suffix _bpcal"
                 + BASED_FLAGGER_PHASE_ARG + BASED_FLAGGER_ALL_FIELDS_ARG)
    print(f'Running BASED flagger on primary: {based_cmd}')
    based_run = subprocess.run([based_cmd], shell=True)
    apply_based_flags(myms, 'baseline_flags_bpcal.txt', bpcal_name, based_run.returncode)
    plot_post_based(myms, bpcal_name, '_bpcal')


# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# --------------------------- Working Table (Primary)  --------------------------- #
# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #


# ------- K (primary; no prior calibration)

gaincal(vis=myms,
    field=bpcal_name,
    caltable=ktab,
    uvrange=primary_uvrange_use,
    # spw=myspw,
    refant=str(ref_ant),
    gaintype='K',
    solint='inf')


# -------- Gp (primary; apply K) -- Type G, phase-only

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    spw=myspw,
    caltable=gptab,
    refant=str(ref_ant),
    gaintype='G',
    solint='inf',
    calmode='p',
    minsnr=3,
    gaintable=[ktab],
    gainfield=[bpcal_name],
    interp=['nearest'])


# ------- B (primary; apply K, Gp)

bandpass(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=bptab,
    refant=str(ref_ant),
    solint='inf',
    combine='',
    solnorm=False,
    minblperant=4,
    minsnr=3.0,
    bandtype='B',
    fillgaps=gapfill,
    gaintable=[ktab, gptab],
    gainfield=[bpcal_name, bpcal_name],
    interp=['nearest', 'linear'])

# tfcrop on bandpass amplitude (CPARAM), all correlations
flagdata(vis=bptab, mode='tfcrop', datacolumn='CPARAM',
    maxnpieces=BPCAL_CPARAM_TFCROP_MAXNPIECES,
    timecutoff=BPCAL_CPARAM_TFCROP_TIMECUTOFF, freqcutoff=BPCAL_CPARAM_TFCROP_FREQCUTOFF,
    extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# L-band only: rflag on bandpass amplitude (CPARAM) over a restricted frequency range
if BPCAL_CPARAM_LBAND_RFLAG and band == 'L':
    flagdata(vis=bptab, mode='rflag', datacolumn='CPARAM', spw=BPCAL_CPARAM_LBAND_RFLAG_SPW,
        timedevscale=BPCAL_CPARAM_LBAND_RFLAG_TIMEDEVSCALE, freqdevscale=BPCAL_CPARAM_LBAND_RFLAG_FREQDEVSCALE,
        extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# rflag on bandpass amplitude (CPARAM), all correlations
flagdata(vis=bptab, mode='rflag', datacolumn='CPARAM',
    timedevscale=BPCAL_CPARAM_RFLAG_TIMEDEVSCALE, freqdevscale=BPCAL_CPARAM_RFLAG_FREQDEVSCALE,
    extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# tfcrop on bandpass phase only (ARG_Sol1,Sol2)
if BPCAL_EXTRA_TFCROP:
    flagdata(vis=bptab, mode='tfcrop', datacolumn='CPARAM', correlation='ARG_Sol1,Sol2',
        maxnpieces=BPCAL_ARG_TFCROP_MAXNPIECES,
        timecutoff=BPCAL_ARG_TFCROP_TIMECUTOFF, freqcutoff=BPCAL_ARG_TFCROP_FREQCUTOFF,
        extendflags=BPCAL_CPARAM_EXTENDFLAGS, flagbackup=False)

# extend on the bandpass table (CPARAM)
if BPCAL_CPARAM_EXTEND:
    flagdata(vis=bptab, mode='extend', extendpols=BPCAL_CPARAM_EXTENDPOLS,
        growtime=BPCAL_CPARAM_EXTEND_TIME, growfreq=BPCAL_CPARAM_EXTEND_FREQ,
        flagbackup=False)

# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

# -------- Ga (primary; apply K, Gp, B) -- Type T, amp-only

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    spw=myspw,
    caltable=gtab,
    refant=str(ref_ant),
    gaintype='T',
    solint='inf',
    calmode='a',
    minsnr=3,
    gaintable=[ktab, gptab, bptab],
    gainfield=[bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear'])


# -------- Solve for Df (primary; apply K, Gp, B, Ga)

polcal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=dftab,
    refant=str(ref_ant),
    solint='inf',
    poltype='Df',
    combine='scan',
    gaintable=[ktab, gptab, bptab, gtab],
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear', 'linear'])

flagdata(vis=dftab, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')


# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- Initial Table (Secondary + Pol. Cal.)  ---------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #

if pacal_name != '':

    # ------- Gp0 (polcal; apply B, Df, K from primary) -- Type G, phase-only

    gaincal(vis=myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gptab0,
        refant=str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['linear', 'linear', 'linear'],
        append=True)

    # ------- G0 (polcal; apply B, Df, K from primary, Gp0 from polcal) -- Type T, amp-only

    gaincal(vis=myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab0,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, bptab, dftab, gptab0],
        gainfield=[bpcal_name, bpcal_name, bpcal_name, pacal_name],
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

    # ------- K0 (polcal; apply B, Df from primary, Gp0, G0 from polcal)

    gaincal(vis=myms,
        field=pacal_name,
        caltable=ktab0,
        uvrange=myuvrange,
        # spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gptab0, gtab0],
        gainfield=[bpcal_name, bpcal_name, pacal_name, pacal_name],
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

# ----- Loop over secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    # ------- Gp0 (pcal; apply B, Df, K from primary) -- Type G, phase-only

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gptab0,
        refant=str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear'],
        append=True)

    # ------- G0 (pcal; apply B, Df, K from primary, Gp0 from pcal) -- Type T, amp-only

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab0,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, bptab, dftab, gptab0],
        gainfield=[bpcal_name, bpcal_name, bpcal_name, pcal],
        interp=['nearest', 'linear', 'linear', 'linear'],
        append=True)

    # ------- K0 (pcal; apply B, Df from primary, Gp0, G0 from pcal)

    gaincal(vis=myms,
        field=pcal,
        caltable=ktab0,
        uvrange=myuvrange,
        # spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gptab0, gtab0],
        gainfield=[bpcal_name, bpcal_name, pcal, pcal],
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

if pacal_name != '':
    # -------- Applycal (polcal; B, Df from primary, K0, Gp0, G0 from polcal) and Flag
    
    applycal(vis=myms,
        gaintable=[ktab0, bptab, gptab0, gtab0, dftab],
        field=pacal_name,
        parang=False,
        gainfield=[pacal_name, bpcal_name, pacal_name, pacal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear', 'linear', 'linear'],
        flagbackup=False)
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

    # rflag on pacal
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='corrected',
        field=pacal_name,
        timedevscale=CAL_VIS_RFLAG_TIMEDEVSCALE,
        freqdevscale=CAL_VIS_RFLAG_FREQDEVSCALE,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    # tfcrop on pacal
    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='corrected',
        field=pacal_name,
        maxnpieces=CAL_VIS_TFCROP_MAXNPIECES,
        timecutoff=CAL_VIS_TFCROP_TIMECUTOFF,
        freqcutoff=CAL_VIS_TFCROP_FREQCUTOFF,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    # extend on pacal
    if CAL_VIS_EXTEND:
        flagdata(vis=myms,
            mode='extend',
            growtime=CAL_VIS_EXTEND_TIME,
            growfreq=CAL_VIS_EXTEND_FREQ,
            growaround=CAL_VIS_EXTEND_GROWAROUND,
            flagneartime=CAL_VIS_EXTEND_FLAGNEARTIME,
            flagnearfreq=CAL_VIS_EXTEND_FLAGNEARFREQ,
            extendpols=CAL_VIS_EXTENDPOLS)

    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

    # ------- BASED flagger (pacal)

    if BASED_FLAGGER and 'pacal' in BASED_FLAGGER_FIELDS_LIST:
        based_cmd = (f"python_dask {TOOLS}/basedflagger.py "
                     f"{myms} {pacal_name} "
                     f"--correlation-products {BASED_FLAGGER_CORRELATION_PRODUCTS} "
                     f"--antenna-flag-cap {BASED_FLAGGER_ANTENNA_FLAG_CAP} --outlier-mode {BASED_FLAGGER_OUTLIER_MODE} "
                     f"--save-dir {VISPLOTS}/ --output-suffix _pacal"
                     + BASED_FLAGGER_PHASE_ARG + BASED_FLAGGER_ALL_FIELDS_ARG)
        print(f'Running BASED flagger on pacal: {based_cmd}')
        based_run = subprocess.run([based_cmd], shell=True)
        apply_based_flags(myms, 'baseline_flags_pacal.txt', pacal_name, based_run.returncode)
        plot_post_based(myms, pacal_name, '_pacal')


# ----- Loop over secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    # -------- Applycal (pcal; B, Df from primary, K0, Gp0, G0 from pcal) and Flag

    applycal(vis=myms,
        gaintable=[ktab0, bptab, gptab0, gtab0, dftab],
        field=pcal,
        parang=False,
        gainfield=[pcal, bpcal_name, pcal, pcal, bpcal_name],
        interp=['nearest', 'linear', 'linear', 'linear', 'linear'],
        flagbackup=False)
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

    # rflag on secondary
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='corrected',
        field=pcal,
        timedevscale=CAL_VIS_RFLAG_TIMEDEVSCALE,
        freqdevscale=CAL_VIS_RFLAG_FREQDEVSCALE,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    # tfcrop on secondary
    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='corrected',
        field=pcal,
        maxnpieces=CAL_VIS_TFCROP_MAXNPIECES,
        timecutoff=CAL_VIS_TFCROP_TIMECUTOFF,
        freqcutoff=CAL_VIS_TFCROP_FREQCUTOFF,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    # extend on secondary
    if CAL_VIS_EXTEND:
        flagdata(vis=myms,
            mode='extend',
            growtime=CAL_VIS_EXTEND_TIME,
            growfreq=CAL_VIS_EXTEND_FREQ,
            growaround=CAL_VIS_EXTEND_GROWAROUND,
            flagneartime=CAL_VIS_EXTEND_FLAGNEARTIME,
            flagnearfreq=CAL_VIS_EXTEND_FLAGNEARFREQ,
            extendpols=CAL_VIS_EXTENDPOLS)

    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

    # ------- BASED flagger (secondary)

    if BASED_FLAGGER and 'secondary' in BASED_FLAGGER_FIELDS_LIST:
        based_cmd = (f"python_dask {TOOLS}/basedflagger.py "
                     f"{myms} {pcal} "
                     f"--correlation-products {BASED_FLAGGER_CORRELATION_PRODUCTS} "
                     f"--antenna-flag-cap {BASED_FLAGGER_ANTENNA_FLAG_CAP} --outlier-mode {BASED_FLAGGER_OUTLIER_MODE} "
                     f"--save-dir {VISPLOTS}/ --output-suffix _{pcal}"
                     + BASED_FLAGGER_PHASE_ARG
                     + based_flag_fields_arg([t for t, c in zip(targets, target_cal_map) if c == pcal]))
        print(f'Running BASED flagger on secondary {pcal}: {based_cmd}')
        based_run = subprocess.run([based_cmd], shell=True)
        apply_based_flags(myms, f'baseline_flags_{pcal}.txt', pcal, based_run.returncode)
        plot_post_based(myms, pcal, f'_{pcal}')

# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- Working Table (Secondary + Pol. Cal.)  ------------------------ #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #

if pacal_name != '':   
 
    # ------- Gp (polcal; apply B, Df, K from primary) -- Type G, phase-only

    gaincal(vis=myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gptab,
        refant=str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear'],
        append=True)

    # ------- G (polcal; apply B, Df, K from primary, Gp from polcal) -- Type T, amp-only

    gaincal(vis=myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, bptab, dftab, gptab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name, pacal_name],
        interp=['nearest', 'linear', 'linear', 'linear'],
        append=True)

    # ------- K (polcal; apply B, Df from primary, Gp, G from polcal)
    
    gaincal(vis=myms,
        field=pacal_name,
        caltable=ktab,
        uvrange=myuvrange,
        # spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gptab, gtab],
        gainfield=[bpcal_name, bpcal_name, pacal_name, pacal_name],
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

            
for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    # ------- Gp (pcal; apply B, Df, K from primary) -- Type G, phase-only

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gptab,
        refant=str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear'],
        append=True)

    # ------- G (pcal; apply B, Df, K from primary, Gp from pcal) -- Type T, amp-only

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, bptab, dftab, gptab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name, pcal],
        interp=['nearest', 'linear', 'linear', 'linear'],
        append=True)

    # ------- K (pcal; apply B, Df from primary, Gp, G from pcal)

    gaincal(vis=myms,
        field=pcal,
        caltable=ktab,
        uvrange=myuvrange,
        # spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gptab, gtab],
        gainfield=[bpcal_name, bpcal_name, pcal, pcal],
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

# ----- Diagnostic: solve cross-hand delay (KCROSS) per secondary. This table is never
# added to cross_table/applycal -- it has no effect on the pipeline -- it's only saved
# for comparing the cross-hand delay recovered per secondary against the one solved on
# pacal_name. Parallel-hand terms only (K, Gp, B, Ga); no Df.
kcross_secondaries = GAINTABLES+'/cal_1GC_'+myms+'.KCROSS_secondaries'
_kcross_diag_first = True
for i in range(0, len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    gaincal(vis=myms,
        field=pcal,
        caltable=kcross_secondaries,
        uvrange=myuvrange,
        # spw=myspw,
        refant=str(ref_ant),
        solint='inf',
        combine='',
        gaintype='KCROSS',
        parang=True,
        gaintable=[ktab, gptab, bptab, gtab],
        gainfield=[pcal, pcal, bpcal_name, pcal],
        interp=['nearest', 'linear', 'linear', 'linear'],
        append=not _kcross_diag_first)
    _kcross_diag_first = False

# --- Apply fluxscaling to G but only if there are calibration fields other than the primary
if len([pcal for pcal in pcal_names if pcal != bpcal_name]) > 0 or pacal_name != '':
    
    # --- Apply fluxscaling to G
    fluxscale(vis=myms,
        caltable=gtab,
        fluxtable=ftab,
        reference=bpcal_name,
        append=False,
        transfer='')

# If there is no need to apply flux scaling, we can set ftab to gtab as the only calibrator is the primary    
else:
    ftab = gtab

# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- CROSS-HAND Tables (PA CAL)  ------------------------------------------------ #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #


if pacal_name != '':
 
    # ------- Set PA calibrator models

    setjy(vis=myms,
        field=pacal_name,
        standard='manual',
        fluxdensity = POLANG_MOD,
        usescratch=True)

    # -------- Solve for Cross-hand phase terms
    if XF_MODE not in ['casa', 'manual', 'auto']:
        print(f"Cross-hand phase mode ({XF_MODE}) not a valid option, defauling to 'auto'")
        XF_MODE = 'auto'

    manual_XF = False
    if XF_MODE == 'casa' or XF_MODE == 'auto':

        print(f"Cross-hand phase mode if {XF_MODE}; solving CASA based solutions")
        
        # ------- KCROSS (polcal; apply K, Gp, B, Ga from polcal/primary -- no Df: the
        # pre-KCROSS Df would bias this delay solve, so it's re-solved below instead)
        gaincal(vis=myms,
            field=pacal_name,
            caltable=kcross,
            uvrange=myuvrange,
            # spw=myspw,
            refant=str(ref_ant),
            solint='inf',
            gaintype='KCROSS',
            parang=True,
            gaintable=[ktab, gptab, bptab, gtab],
            gainfield=[pacal_name, pacal_name, bpcal_name, pacal_name],
            interp=['nearest', 'linear', 'linear', 'linear'],
            append=False)

        if not XF_SKIP_KCROSS:
            # ------- Re-solve Df (primary; apply K, Gp, B, Ga, KCROSS) into a separate table
            # CASA enforces the apply order KCROSS -> Df -> Xf, and these tables do not
            # commute, so Df must be re-solved now that KCROSS exists -- otherwise the
            # original (KCROSS-free) Df would be incoherent with the applied KCROSS/Xf terms.
            print(f"  Re-solving Df (leakage) with KCROSS included: {dftab6}")
            polcal(vis=myms,
                field=bpcal_name,
                uvrange=primary_uvrange_use,
                caltable=dftab6,
                refant=str(ref_ant),
                solint='inf',
                poltype='Df',
                combine='scan',
                gaintable=[ktab, gptab, bptab, gtab, kcross],
                gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name, pacal_name],
                interp=['nearest', 'linear', 'linear', 'linear', 'linear'])

            flagdata(vis=dftab6, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')

            # Active Df is now the KCROSS-informed table
            dftab = dftab6

        # -------- Xf (polcal; apply K, Gp, B, Ga, Df, KCROSS from polcal/primary)

        polcal(vis=myms,
            field=pacal_name,
            uvrange=myuvrange,
            caltable=xftab,
            refant=str(ref_ant),
            solint=f'inf,{XF_CHANINT}ch',
            poltype='Xf',
            combine='',
            gaintable=[ktab, gptab, bptab, gtab, dftab, kcross],
            gainfield=[pacal_name, pacal_name, bpcal_name, pacal_name, bpcal_name, pacal_name],
            interp=['nearest', 'linear', 'linear', 'linear', 'linear', 'nearest'],
            append=False)

        # Cross hand calibration tables
        cross_table = [kcross, xftab]
        cross_field = [pacal_name, pacal_name]
        cross_interp = ['nearest', 'linear']

        # Always check if the cross-hand phase is continuous per scan
        # Get the cross-hand phase            
        tb.open(xftab)
        gains = tb.getcol('CPARAM')
        flags = tb.getcol('FLAG')
        scans = tb.getcol('SCAN_NUMBER')
        tb.close()

        print(f"Gains shape: {gains.shape}")
        print(f"Unique scans: {np.unique(scans)}")

        # Process per scan
        unique_scans = np.unique(scans)
        manual_XF_per_scan = {}

        for scan in unique_scans:
            scan_mask = (scans == scan)
           
            # Get gains for this scan (median over time axis)
            scan_gains = np.nanmedian(gains[0, :, scan_mask].T, axis=-1)
            scan_flags = np.nanmedian(flags[0, :, scan_mask].T, axis=-1).astype(bool)
            scan_gains = scan_gains[~scan_flags]
            
            # Check continuity
            phases = np.angle(scan_gains)
            phase_diffs = np.diff(phases)
            phase_diffs = np.arctan2(np.sin(phase_diffs), np.cos(phase_diffs))
            max_jump = np.max(np.abs(np.degrees(phase_diffs)))
            
            is_continuous = max_jump < XF_AUTO_ANG_JUMP
            manual_XF_per_scan[scan] = not is_continuous
            
            print(f"Scan {scan}: max jump = {max_jump:.2f}°, continuous = {is_continuous}")

        # Determine which scans are continuous
        continuous_scans = [scan for scan in unique_scans if not manual_XF_per_scan[scan]]
        discontinuous_scans = [scan for scan in unique_scans if manual_XF_per_scan[scan]]

        print("\n" + "="*60)
        
        if len(continuous_scans) == len(unique_scans):
            # Case 1: All scans are continuous
            print(f"✓ All {len(unique_scans)} scan(s) are phase continuous")
            
            if XF_SKIP_KCROSS:
                print("  XF_SKIP_KCROSS = True: Skipping KCROSS, using XF only")
                # Remove the KCROSS table and remake XF without it
                old_kcross = kcross.replace('.KCROSS', '_notused.KCROSS')
                old_xftab = xftab.replace('.Xf', '_withKCROSS.Xf')
                print(f"  Moving KCROSS table to: {old_kcross}")
                if os.path.isdir(old_kcross):
                    shutil.rmtree(old_kcross)
                shutil.move(kcross, old_kcross)
                print(f"  Moving KCROSS-based Xf table to: {old_xftab}")
                if os.path.isdir(old_xftab):
                    shutil.rmtree(old_xftab)
                shutil.move(xftab, old_xftab)
                
                # Remake Xf table without KCROSS
                print(f"  Remaking Xf table without KCROSS")
                polcal(vis=myms,
                    field=pacal_name,
                    uvrange=myuvrange,
                    caltable=xftab,
                    refant=str(ref_ant),
                    solint=f'inf,{XF_CHANINT}ch',
                    poltype='Xf',
                    combine='',
                    gaintable=[ktab, gptab, bptab, gtab, dftab],
                    gainfield=[pacal_name, pacal_name, bpcal_name, pacal_name, bpcal_name],
                    interp=['nearest', 'linear', 'linear', 'linear', 'linear'],
                    append=False)
                
                # Update cross-hand calibration tables (no KCROSS)
                cross_table = [xftab]
                cross_field = [pacal_name]
                cross_interp = ['linear']
            else:
                print("  Proceeding with standard XF calibration (with KCROSS)")
            
        elif len(continuous_scans) > 0:
            # Case 2: Some scans are continuous, some are not
            print(f"✓ {len(continuous_scans)} scan(s) are phase continuous: {continuous_scans}")
            print(f"✗ {len(discontinuous_scans)} scan(s) failed continuity check: {discontinuous_scans}")
            
            if XF_SKIP_KCROSS:
                print("  XF_SKIP_KCROSS = True: Skipping KCROSS, using XF only")
                # Remove old tables
                old_kcross = kcross.replace('.KCROSS', '_notused.KCROSS')
                old_xftab = xftab.replace('.Xf', '_withKCROSS.Xf')
                print(f"  Moving KCROSS table to: {old_kcross}")
                if os.path.isdir(old_kcross):
                    shutil.rmtree(old_kcross)
                shutil.move(kcross, old_kcross)
                print(f"  Moving KCROSS-based Xf table to: {old_xftab}")
                if os.path.isdir(old_xftab):
                    shutil.rmtree(old_xftab)
                shutil.move(xftab, old_xftab)
                
                # Remake Xf table with only continuous scans, without KCROSS
                print(f"  Remaking Xf table using only continuous scans: {continuous_scans}")
                polcal(vis=myms,
                    field=pacal_name,
                    scan=','.join(map(str, continuous_scans)),
                    uvrange=myuvrange,
                    caltable=xftab,
                    refant=str(ref_ant),
                    solint=f'inf,{XF_CHANINT}ch',
                    poltype='Xf',
                    combine='',
                    gaintable=[ktab, gptab, bptab, gtab, dftab],
                    gainfield=[pacal_name, pacal_name, bpcal_name, pacal_name, bpcal_name],
                    interp=['nearest', 'linear', 'linear', 'linear', 'linear'],
                    append=False)
                
                # Update cross-hand calibration tables (no KCROSS)
                cross_table = [xftab]
                cross_field = [pacal_name]
                cross_interp = ['linear']
            else:
                # Move bad tables
                bad_kcross = kcross.replace('.KCROSS', '_bad.KCROSS')
                bad_xftab = xftab.replace('.Xf', '_bad.Xf')
                print(f"  Moving original KCROSS table to: {bad_kcross}")
                if os.path.isdir(bad_kcross):
                    shutil.rmtree(bad_kcross)
                shutil.move(kcross, bad_kcross)
                print(f"  Moving original Xf table to: {bad_xftab}")
                if os.path.isdir(bad_xftab):
                    shutil.rmtree(bad_xftab)
                shutil.move(xftab, bad_xftab)
                
                # Remake KCROSS with only continuous scans (no Df -- see note above)
                print(f"  Remaking KCROSS table using only continuous scans: {continuous_scans}")
                gaincal(vis = myms,
                    field = pacal_name,
                    scan = ','.join(map(str, continuous_scans)),
                    caltable = kcross,
                    uvrange=myuvrange,
                    # spw=myspw,
                    refant = str(ref_ant),
                    solint = 'inf',
                    gaintype='KCROSS',
                    parang = True,
                    gaintable=[ktab, gptab, bptab, gtab],
                    gainfield=[pacal_name, pacal_name, bpcal_name, pacal_name],
                    interp = ['nearest','linear','linear','linear'],
                    append = False)

                # dftab already points at dftab6 from the first KCROSS-informed solve
                # above -- overwrite it again with the continuous-scans-only KCROSS.
                print(f"  Re-solving Df (leakage) with KCROSS included, overwriting: {dftab}")
                polcal(vis=myms,
                    field=bpcal_name,
                    uvrange=primary_uvrange_use,
                    caltable=dftab,
                    refant=str(ref_ant),
                    solint='inf',
                    poltype='Df',
                    combine='scan',
                    gaintable=[ktab, gptab, bptab, gtab, kcross],
                    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name, pacal_name],
                    interp=['nearest', 'linear', 'linear', 'linear', 'linear'])

                flagdata(vis=dftab, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')

                # Remake Xf table with only continuous scans
                print(f"  Remaking Xf table using only continuous scans: {continuous_scans}")
                polcal(vis = myms,
                    field = pacal_name,
                    scan = ','.join(map(str, continuous_scans)),
                    uvrange = myuvrange,
                    caltable = xftab,
                    refant = str(ref_ant),
                    solint = f'inf,{XF_CHANINT}ch',
                    poltype='Xf',
                    combine = '',
                    gaintable=[ktab, gptab, bptab, gtab, dftab, kcross],
                    gainfield=[pacal_name, pacal_name, bpcal_name, pacal_name, bpcal_name, pacal_name],
                    interp = ['nearest','linear','linear','linear','linear','nearest'],
                    append = False)
            
        else:
            # Case 3: No scans are continuous
            print(f"✗ None of the {len(unique_scans)} scan(s) are phase continuous")
            
            # Move bad tables
            bad_kcross = kcross.replace('.KCROSS', '_bad.KCROSS')
            bad_xftab = xftab.replace('.Xf', '_bad.Xf')
            print(f"  Moving original KCROSS table to: {bad_kcross}")
            if os.path.isdir(bad_kcross):
                shutil.rmtree(bad_kcross)
            shutil.move(kcross, bad_kcross)
            print(f"  Moving original Xf table to: {bad_xftab}")
            if os.path.isdir(bad_xftab):
                shutil.rmtree(bad_xftab)
            shutil.move(xftab, bad_xftab)

            if not XF_SKIP_KCROSS:
                # KCROSS is abandoned -- revert to the original (KCROSS-free) Df table
                print(f"  Reverting Df to original table (KCROSS abandoned): {dftab_original}")
                dftab = dftab_original

            if len(unique_scans) == 1:
                # Single scan - can use manual solver
                print("  Single scan detected - falling back to manual XF solver")
                manual_XF = True
            else:
                # Multiple scans - cannot handle
                print("\nERROR: Cannot currently handle multiple cross-hand phase scans in a single MS.")
                print("       Cross-hand phase is stable over week timescales, so unless this is a")
                print("       week-long MS file, one scan should be sufficient.")
                print(f"       Please split out an MS file with only one XF scan and re-run.")
                sys.exit(1)

        print("="*60)

    if XF_MODE == 'manual' or manual_XF or band == 'UHF':
        print(f"Cross-hand phase mode is {XF_MODE} (or auto detected a large phase jump? {manual_XF}); solving manual solutions")
        exec(open('tools/manual_XF_solver.py').read())

        # Cross hand calibration tables
        cross_table = [xftab]
        cross_field = [pacal_name]
        cross_interp = ['linear']


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
# ------------------------- Global cross-hand phase  ---------------------------- #
# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #

print('\n=== Running global_crosshand_phase.py ===')
subprocess.run(
    ['casa', '--nologger', '--log2term', '--nogui', '-c',
     'tools/global_crosshand_phase.py'],
    check=True)
print('=== global_crosshand_phase.py complete ===\n')


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
# --------------------------- Applycal (All Fields)  ----------------------- #
# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #

# ------- BPCAL

applycal(vis = myms,
    gaintable = [ktab, gptab, bptab, ftab, dftab_original],
 #  sapplymode='calflagstrict',
    field = bpcal_name,
    #calwt = False,
    parang = False,
    gainfield = [bpcal_name, bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp = ['nearest','linear','linear','linear','linear'],
    flagbackup=False)
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')


# ----- If no polarization angle calibrator apply subset of tables and kill script

if pacal_name == '':

    # Determine whether to include a global Xf table from ALT_PACAL
    if ALT_PACAL != '':
        _alt_xftab = GAINTABLES + '/cal_1GC_' + myms + (
            '_combineScan.Xf' if ALT_TYPE == 'combine' else '_perScan.Xf')
        _alt_cross_table  = [_alt_xftab]
        _alt_cross_field  = [ALT_PACAL]
        _alt_cross_interp = ['linear']
        print(f'  Applying global Xf table: {_alt_xftab}  (ALT_PACAL={ALT_PACAL}, ALT_TYPE={ALT_TYPE})')
    else:
        _alt_cross_table  = []
        _alt_cross_field  = []
        _alt_cross_interp = []

    # ------- Secondaries

    for i in range(0, len(pcal_names)):

        pcal = pcal_names[i]

        # ------- Check if pcal is the same as bpcal_name or pacal_name
        if pcal == bpcal_name or pcal == pacal_name:
            # If so, skip to the next iteration as it is already in the working tables
            continue

        applycal(vis = myms,
            gaintable = [ktab, gptab, bptab, ftab, dftab] + _alt_cross_table,
            # applymode='calflagstrict',
            field = pcal,
            #calwt = False,
            parang = False,
            gainfield = [pcal, pcal, bpcal_name, pcal, bpcal_name] + _alt_cross_field,
            interp = ['nearest','linear','linear','linear','linear'] + _alt_cross_interp,
            flagbackup=False)
        # DEBUGGING: summarize flags
        if DEBUG_PRINT_FLAGS:
            print('DEBUG: PRINTING FLAGS')
            flagdata(myms, mode='summary')

    # ------- Targets
    for i in range(0, len(targets)):

        target = targets[i]
        related_pcal = target_cal_map[i]

        applycal(vis=myms,
                #applymode='calflagstrict',
                gaintable = [ktab, gptab, bptab, ftab, dftab] + _alt_cross_table,
                field=target,
                #calwt=False,
                parang=False,
                gainfield = ['', related_pcal, bpcal_name, related_pcal, bpcal_name] + _alt_cross_field,
                interp = ['nearest','linear','linear','linear','linear'] + _alt_cross_interp,
                flagbackup=False)
        # DEBUGGING: summarize flags
        if DEBUG_PRINT_FLAGS:
            print('DEBUG: PRINTING FLAGS')
            flagdata(myms, mode='summary')

        # Flag target
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

        # ------- BASED flagger (target)

        if BASED_FLAGGER and 'target' in BASED_FLAGGER_FIELDS_LIST:
            based_cmd = (f"python_dask {TOOLS}/basedflagger.py "
                         f"{myms} {target} "
                         f"--correlation-products {BASED_FLAGGER_CORRELATION_PRODUCTS} "
                         f"--antenna-flag-cap {BASED_FLAGGER_ANTENNA_FLAG_CAP} --outlier-mode {BASED_FLAGGER_TARGET_OUTLIER_MODE} "
                         f"--save-dir {VISPLOTS}/ --output-suffix _{target}" + BASED_FLAGGER_PHASE_ARG)
            print(f'Running BASED flagger on target {target}: {based_cmd}')
            based_run = subprocess.run([based_cmd], shell=True)
            apply_based_flags(myms, f'baseline_flags_{target}.txt', target, based_run.returncode)
            plot_post_based(myms, target, f'_{target}')

    # ---- Save flags

    flagmanager(vis=myms,
        mode='delete',
        versionname='1GC_flags')

    flagmanager(vis=myms,
        mode='save',
        versionname='1GC_flags')

    sys.exit('Ending Early! No polarization angle calibrator')

# -------- Full polarization 

# ------- PACAL

applycal(vis = myms,
 #       applymode='calflagstrict',
        field = pacal_name,
        #calwt = False,
        parang = False,
        gaintable = [ktab, gptab, bptab, ftab, dftab] + cross_table,
        gainfield = [pacal_name, pacal_name, bpcal_name, pacal_name, bpcal_name] + cross_field,
        interp = ['nearest','linear','linear','linear','linear'] + cross_interp,
        flagbackup=False)
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

# ------- Secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    applycal(vis = myms,
        # applymode='calflagstrict',
        field = pcal,
        #calwt = False,
        parang = False,
        gaintable = [ktab, gptab, bptab, ftab, dftab] + cross_table,
        gainfield = [pcal, pcal, bpcal_name, pcal, bpcal_name] + cross_field,
        interp = ['nearest','linear','linear','linear','linear'] + cross_interp,
        flagbackup=False)
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

# ------- Targets 
for i in range(0,len(targets)):

    target = targets[i]
    related_pcal = target_cal_map[i]

    applycal(vis=myms,
                #applymode='calflagstrict',
                field=target,
                #calwt=False,
                parang=False,
                gaintable = [ktab, gptab, bptab, ftab, dftab] + cross_table,
                gainfield = ['', related_pcal, bpcal_name, related_pcal, bpcal_name] + cross_field,
                interp = ['nearest','linear','linear','linear','linear'] + cross_interp,
                flagbackup=False)
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

    # Flag target
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

    # ------- BASED flagger (target)

    if BASED_FLAGGER and 'target' in BASED_FLAGGER_FIELDS_LIST:
        based_cmd = (f"python_dask {TOOLS}/basedflagger.py "
                     f"{myms} {target} "
                     f"--correlation-products {BASED_FLAGGER_CORRELATION_PRODUCTS} "
                     f"--antenna-flag-cap {BASED_FLAGGER_ANTENNA_FLAG_CAP} --outlier-mode {BASED_FLAGGER_TARGET_OUTLIER_MODE} "
                     f"--save-dir {VISPLOTS}/ --output-suffix _{target}" + BASED_FLAGGER_PHASE_ARG)
        print(f'Running BASED flagger on target {target}: {based_cmd}')
        based_run = subprocess.run([based_cmd], shell=True)
        apply_based_flags(myms, f'baseline_flags_{target}.txt', target, based_run.returncode)
        plot_post_based(myms, target, f'_{target}')

# ---- Apply aggressive flags if desired
if CAL_1GC_AGGRESSIVE_FLAGS and CAL_1GC_BL_FREQS != []:

    flagspw = ','.join(CAL_1GC_BL_FREQS)

    flagdata(vis = myms,
        mode = 'manual',
        spw = flagspw,
        flagbackup=False)
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

# ---- Save flags

flagmanager(vis=myms,
    mode='delete',
    versionname='1GC_flags')

flagmanager(vis=myms,
    mode='save',
    versionname='1GC_flags')


# ---- Flag summary

# The flag state once 1GC is done, written per run so the fraction removed can
# be followed from one reduction to the next. The per-field breakdown is
# printed as well, since that is where a flagging problem shows up first.

flag_summary = flagdata(vis=myms, mode='summary')

if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)

flag_summary_path = RESULTS+'/flagsummary_1GC_'+myms+'_'+stamp()+'.json'

def jsonable(obj):
    """Fall-back encoder for the numpy scalars CASA puts in the summary."""
    return obj.tolist() if hasattr(obj,'tolist') else str(obj)

with open(flag_summary_path,'w') as j:
    json.dump(flag_summary, j, indent=4, sort_keys=True, default=jsonable)

print(f'Flag summary written to: {flag_summary_path}')

try:
    flagged_pc = 100.0*flag_summary['flagged']/flag_summary['total']
    print(f'Flagged after 1GC: {flagged_pc:.4f}% of '
          f'{int(flag_summary["total"]):,} visibilities')
    for field_name, counts in sorted(flag_summary.get('field',{}).items()):
        if counts['total']:
            print(f'  {field_name}: '
                  f'{100.0*counts["flagged"]/counts["total"]:.4f}% '
                  f'({int(counts["flagged"]):,}/{int(counts["total"]):,})')
        else:
            print(f'  {field_name}: no visibilities')
except (KeyError, TypeError, ZeroDivisionError):
    print('Flag summary: breakdown unavailable')
