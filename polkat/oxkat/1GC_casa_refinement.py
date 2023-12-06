# ian.heywood@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import sys


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    targets = user_targets
    pcals = user_pcals
    target_cal_map = user_cal_map

def stamp():
    now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
    return now
    

# ------- Parameters


gapfill = CAL_1GC_FILLGAPS
myuvrange = CAL_1GC_UVRANGE 
myspw = CAL_1GC_FREQRANGE
solint = sys.argv[-1]
field  = sys.argv[-2]


# ------- Setup names

tt = stamp()

ktab = GAINTABLES+'/cal_1GC_'+myms+'_{}'.format(field)+'.K'
gptab = GAINTABLES+'/cal_1GC_'+myms+'_{}'.format(field)+'.Gp'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'

# ------- Reload flags

if SAVE_FLAGS:
    flagmanager(vis=myms,mode='restore',versionname='bpcal_final_flags')

# ------- Flag CORRECTED_DATA - MODEL data column

flagdata(vis=myms,
    mode='rflag',
    datacolumn='residual',
    field=field, flagbackup=False)

flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='residual',
    field=field, flagbackup=False)

# -------- Solve for the pol cal gain calibrators
gaincal(vis=myms,
    field=field,
    caltable=ktab,
    refant = str(ref_ant),
    gaintype = 'K',
    solint = 'inf')

gaincal(vis=myms,
    field=field,
    uvrange=myuvrange,
    caltable=gptab,
    gaintype='G',
    solint=solint,
    calmode='p',
    gaintable=[ktab],
    gainfield=[field],
    interp = ['nearest'])

applycal(vis=myms,
    gaintable=[ktab,gptab,bptab,gatab],
#    applymode='calonly',
    field=field,
#    calwt=False,
    parang=False,
    gainfield=[field,field,bpcal_name, bpcal_name],
    interp = ['nearest','nearest','linear','linear'], flagbackup=False)
