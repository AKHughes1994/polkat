# ian.heywood@physics.ox.ac.uk

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    pcals = user_pcals

clearstat()
clearstat()

flagdata(vis=myms,mode='rflag',datacolumn='data',field=bpcal_name)
flagdata(vis=myms,mode='tfcrop',datacolumn='data',field=bpcal_name)
flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=bpcal_name)


flagdata(vis=myms,mode='rflag',datacolumn='data',field=pacal_name)
flagdata(vis=myms,mode='tfcrop',datacolumn='data',field=pacal_name)
flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=pacal_name)


for pcal in pcals:
    flagdata(vis=myms,mode='rflag',datacolumn='data',field=pcal)
    flagdata(vis=myms,mode='tfcrop',datacolumn='data',field=pcal)
    flagdata(vis=myms,mode='extend',growtime=90.0,growfreq=90.0,growaround=True,flagneartime=True,flagnearfreq=True,field=pcal)


if SAVE_FLAGS:
   flagmanager(vis=myms,mode='save',versionname='autoflag_cals_data')


clearstat()
clearstat()
