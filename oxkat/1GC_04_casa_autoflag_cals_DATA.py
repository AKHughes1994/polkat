# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    pcal_names = user_pcals

clearstat()
clearstat()

# Flagging the primary calibrator data
flagdata(vis=myms,mode='rflag',datacolumn='data',field=bpcal_name, flagbackup=False)
flagdata(vis=myms,mode='tfcrop',datacolumn='data',field=bpcal_name, flagbackup=False)
flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=bpcal_name, flagbackup=False)

# Flagging the polarisation angle calibrator data, if it exists
if pacal_name != '' :
    flagdata(vis=myms,mode='rflag',datacolumn='data',field=pacal_name, flagbackup=False)
    flagdata(vis=myms,mode='tfcrop',datacolumn='data',field=pacal_name, flagbackup=False)
    flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=pacal_name, flagbackup=False)

# For edge cases where a phase calibrator is either the primary or the polarisation angle calibrator, we need to ensure that we don't double flag it.
pcal_names = [pcal for pcal in pcal_names if pcal != bpcal_name and pcal != pacal_name]

for pcal in pcal_names:
    flagdata(vis=myms,mode='rflag',datacolumn='data',field=pcal, flagbackup=False)
    flagdata(vis=myms,mode='tfcrop',datacolumn='data',field=pcal, flagbackup=False)
    flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=pcal, flagbackup=False)

# Saving the flag state
flagmanager(vis=myms,mode='save',versionname='autoflag_cals_data')


clearstat()
clearstat()
