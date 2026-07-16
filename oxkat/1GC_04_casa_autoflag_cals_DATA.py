# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

# --- RFlag parameters (flagdata mode='rflag') ---
RFLAG_TIMEDEVSCALE = 5.0   # Flagging threshold (sigma) for time-direction deviation
RFLAG_FREQDEVSCALE = 5.0   # Flagging threshold (sigma) for freq-direction deviation

# --- TFCrop parameters (flagdata mode='tfcrop') ---
TFCROP_TIMECUTOFF = 4.0    # Flagging threshold (sigma) in the time direction
TFCROP_FREQCUTOFF = 3.0    # Flagging threshold (sigma) in the freq direction

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
    flagbackup=False)

flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='data',
    field=bpcal_name,
    timecutoff=TFCROP_TIMECUTOFF,
    freqcutoff=TFCROP_FREQCUTOFF,
    flagbackup=False)

flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=bpcal_name, flagbackup=False)

# Flagging the polarisation angle calibrator data, if it exists
if pacal_name != '' :
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='data',
        field=pacal_name,
        timedevscale=RFLAG_TIMEDEVSCALE,
        freqdevscale=RFLAG_FREQDEVSCALE,
        flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='data',
        field=pacal_name,
        timecutoff=TFCROP_TIMECUTOFF,
        freqcutoff=TFCROP_FREQCUTOFF,
        flagbackup=False)

    flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=pacal_name, flagbackup=False)

# For edge cases where a phase calibrator is either the primary or the polarisation angle calibrator, we need to ensure that we don't double flag it.
pcal_names = [pcal for pcal in pcal_names if pcal != bpcal_name and pcal != pacal_name]

for pcal in pcal_names:
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='data',
        field=pcal,
        timedevscale=RFLAG_TIMEDEVSCALE,
        freqdevscale=RFLAG_FREQDEVSCALE,
        flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='data',
        field=pcal,
        timecutoff=TFCROP_TIMECUTOFF,
        freqcutoff=TFCROP_FREQCUTOFF,
        flagbackup=False)

    flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=pcal, flagbackup=False)

# Saving the flag state
flagmanager(vis=myms,mode='save',versionname='autoflag_cals_data')


clearstat()
clearstat()
