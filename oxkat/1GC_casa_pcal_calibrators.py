# ian.heywood@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import sys
from argparse import ArgumentParser


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    targets = user_targets
    pcals = user_pcals
    target_cal_map = user_cal_map

def stamp():
    now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
    return now
    

# ------ Read in arguments

parser = ArgumentParser(description='phase self-calibrate an arbitrary field (Stokes I only) [options] field')
parser.add_argument('field',
                        help="target field to calibrate")
parser.add_argument('-s', '--solint', dest="solint", default='inf',
                        help="String of solution interval [seconds/int/inf] to solve for phase solutions")
parser.add_argument('-i', '--initial-flag', dest='initial_flag', default='parr_hand_flags',
                        help='Initial set of flags to restore the MS file to before calibration')
parser.add_argument('-f', '--final-flag', dest='final_flag', default='pcal_flag',
                        help='Final set of flags to save to the MS file post-calibration')

args = parser.parse_args()
field = str(args.field)
solint = str(args.solint)
initial_flag = str(args.initial_flag)
final_flag = str(args.final_flag)

# ------- Parameters


gapfill = CAL_1GC_FILLGAPS
myuvrange = CAL_1GC_UVRANGE 
myspw = CAL_1GC_FREQRANGE

# ------- Setup names

tt = stamp()

ktab = GAINTABLES+'/cal_1GC_'+myms+f'_{field}.K'
gptab = GAINTABLES+'/cal_1GC_'+myms+f'_{field}.Gp'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'

# ------- Reload flags

flagmanager(vis=myms,mode='restore',versionname=initial_flag)

# ------- Flag CORRECTED_DATA - MODEL data column

flagdata(vis=myms,
    mode='rflag',
    datacolumn='residual',
    field=field, flagbackup=False, correlation='ABS_I')

flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='residual',
    field=field, flagbackup=False, correlation='ABS_I')

# -------- Solve for the K Solutions using self-cal model

gaincal(vis=myms,
    field=field,
    #uvrange=myuvrange,
    caltable=ktab,
    refant = str(ref_ant),
    gaintype = 'K',
    solint = 'inf',
    gaintable=[bptab,gatab],
    interp = ['linear','linear'],
    gainfield=[bpcal_name, bpcal_name])

# -------- Solve for the Gp Solutions using self-cal model

gaincal(vis=myms,
    field=field,
    uvrange=myuvrange,
    caltable=gptab,
    gaintype='G',
    refant = str(ref_ant),
    solint=solint,
    minsnr=5,
    calmode='p',
    gaintable=[bptab,gatab, ktab],
    interp = ['linear','linear','nearest'],
    gainfield=[bpcal_name, bpcal_name, field])


applycal(vis=myms,
    gaintable=[ktab,gptab,bptab,gatab],
    #applymode='calflagstrict',
    field=field,
    calwt=False,
    parang=False,
    gainfield=[field,field,bpcal_name, bpcal_name],
    interp = ['nearest','nearest','linear','linear'], flagbackup=False)

# ----- Save final flags for selfcal iteration

flagmanager(vis=myms,
        mode='delete',
        versionname=final_flag)

flagmanager(vis=myms,
        mode='save',
        versionname=final_flag)
