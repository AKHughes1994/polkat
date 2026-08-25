# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import functools
import json

print = functools.partial(print, flush=True)

# --- RFlag parameters (flagdata mode='rflag') ---
RFLAG_TIMEDEVSCALE = 6.0   # Flagging threshold (sigma) for time-direction deviation
RFLAG_FREQDEVSCALE = 6.0   # Flagging threshold (sigma) for freq-direction deviation

# --- TFCrop parameters (flagdata mode='tfcrop') ---
TFCROP_TIMECUTOFF = 5.0    # Flagging threshold (sigma) in the time direction
TFCROP_FREQCUTOFF = 4.0    # Flagging threshold (sigma) in the freq direction

# --- Extend flags across correlations (flagdata rflag/tfcrop 'extendflags' param) ---
EXTEND_AUTO = False

# --- Extend-mode parameters (flagdata mode='extend' growtime/growfreq, % already flagged) ---
EXTEND_TIME = 80.0
EXTEND_FREQ = 80.0

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    pcal_names = user_pcals

clearstat()
clearstat()

# Load the basic flags saved by 1GC_03_casa_basic_flags.py
flagmanager(vis = myms,
        mode = 'restore',
        versionname = 'basic')

# Flagging the primary calibrator data
flagdata(vis=myms,
    mode='rflag',
    datacolumn='data',
    field=bpcal_name,
    timedevscale=RFLAG_TIMEDEVSCALE,
    freqdevscale=RFLAG_FREQDEVSCALE,
    extendflags=EXTEND_AUTO,
    flagbackup=False)

flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='data',
    field=bpcal_name,
    timecutoff=TFCROP_TIMECUTOFF,
    freqcutoff=TFCROP_FREQCUTOFF,
    extendflags=EXTEND_AUTO,
    flagbackup=False)

flagdata(vis=myms,mode='extend',growtime=EXTEND_TIME,growfreq=EXTEND_FREQ,growaround=True,flagneartime=True,flagnearfreq=True,field=bpcal_name, flagbackup=False)

# Flagging the polarisation angle calibrator data, if it exists
if pacal_name != '' :
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='data',
        field=pacal_name,
        timedevscale=RFLAG_TIMEDEVSCALE,
        freqdevscale=RFLAG_FREQDEVSCALE,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='data',
        field=pacal_name,
        timecutoff=TFCROP_TIMECUTOFF,
        freqcutoff=TFCROP_FREQCUTOFF,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    flagdata(vis=myms,mode='extend',growtime=EXTEND_TIME,growfreq=EXTEND_FREQ,growaround=True,flagneartime=True,flagnearfreq=True,field=pacal_name, flagbackup=False)

# For edge cases where a phase calibrator is either the primary or the polarisation angle calibrator, we need to ensure that we don't double flag it.
pcal_names = [pcal for pcal in pcal_names if pcal != bpcal_name and pcal != pacal_name]

for pcal in pcal_names:
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='data',
        field=pcal,
        timedevscale=RFLAG_TIMEDEVSCALE,
        freqdevscale=RFLAG_FREQDEVSCALE,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='data',
        field=pcal,
        timecutoff=TFCROP_TIMECUTOFF,
        freqcutoff=TFCROP_FREQCUTOFF,
        extendflags=EXTEND_AUTO,
        flagbackup=False)

    flagdata(vis=myms,mode='extend',growtime=EXTEND_TIME,growfreq=EXTEND_FREQ,growaround=True,flagneartime=True,flagnearfreq=True,field=pcal, flagbackup=False)

# Saving the flag state
flagmanager(vis=myms,mode='save',versionname='autoflag_cals_data')


# Refant selection by post-flagging delay-solve S/N on the primary
if CAL_1GC_REFANT_SNR_SELECT:

    print(f'Solving per-scan delay S/N on the primary ({bpcal_name}) for refant ranking...')
    exec(open(f'{TOOLS}/antenna_delay_snr.py').read())

    top_names = [name for name, snr in refant_snr_ranking[:7]]
    if 'm060' in top_names and top_names[0] != 'm060':
        top_names.remove('m060')
        top_names.insert(0, 'm060')

    tb.open(myms+'/ANTENNA')
    ant_names = list(tb.getcol('NAME'))
    tb.close()

    new_ref_ant = ','.join(str(ant_names.index(name)) for name in top_names)
    print(f'Top {len(top_names)} antennas by average delay-solve S/N: {top_names}')
    print(f'Replacing ref_ant with: {new_ref_ant}')

    with open('project_info.json') as f:
        project_info = json.load(f)
    project_info['ref_ant'] = new_ref_ant
    with open('project_info.json', 'w') as f:
        json.dump(project_info, f, indent=4, sort_keys=True)


clearstat()
clearstat()
