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
solint = str(sys.argv[-1])
field  = str(sys.argv[-2])
pcal   = str(sys.argv[-3])

# ------- Setup names

tt = stamp()

#ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
ktab = GAINTABLES+'/cal_1GC_'+myms+'_{}'.format(field)+'.K'
sctab = GAINTABLES+'/cal_1GC_'+myms+'_{}'.format(field)+'.Gp'
ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
gptab = GAINTABLES+'/cal_1GC_'+myms+'.Gp'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'
tatab = GAINTABLES+'/cal_1GC_'+myms+'.Ta'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'
kcross = GAINTABLES+'/cal_1GC_'+myms+'.Kcross'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'

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
#gaincal(vis=myms,
#    field=field,
#    caltable=ktab,
#    parang=True,
#    refant = str(ref_ant),
#    gaintype = 'K',
#    solint = 'inf')

gaincal(vis=myms,
    field=field,
    uvrange=myuvrange,
    caltable=sctab,
    gaintype='G',
    solint='int,128MHz',
    calmode='p',
    parang=True,
    gaintable=[ktab,gptab.replace('.Gp', '_{}.Gp'.format(pcal)),bptab,gatab, tatab, dftab, kcross, xftab],
    gainfield=[pcal,pcal, bpcal_name, bpcal_name, pcal, bpcal_name, pacal_name, pacal_name],
    interp = ['linear','linear','linear','linear', 'linear','linear','linear', 'linear'])

applycal(vis=myms,
#            applymode='calonly',
            field=field,
#            calwt=False,
            parang=True,
            gaintable=[ktab,gptab.replace('.Gp', '_{}.Gp'.format(pcal)),bptab,gatab, tatab, dftab, kcross, xftab,sctab],
            gainfield=[pcal,pcal, bpcal_name, bpcal_name, pcal, bpcal_name, pacal_name, pacal_name,field],
            interp = ['linear','linear','linear','linear', 'linear','linear','linear', 'linear', 'nearest'])
