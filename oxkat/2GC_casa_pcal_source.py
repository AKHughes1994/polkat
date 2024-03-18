# ian.heywood@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import sys
from argparse import ArgumentParser


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

def stamp():
    now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
    return now
    

# ------ Read in arguments

parser = ArgumentParser(description='phase self-calibrate an arbitrary field (Stokes I only) [options] field')
parser.add_argument('field',
                        help="target field to calibrate")
parser.add_argument('-s', '--solint', dest="solint", default='inf',
                        help="String of solution interval [seconds/int/inf] to solve for phase solutions")
parser.add_argument('-f', '--final-flag', dest='final_flag', default='pcal_flag',
                        help='Final set of flags to save to the MS file post-calibration')

args = parser.parse_args()
field = str(args.field)
solint = str(args.solint)
final_flag = str(args.final_flag)

# ------ Target information

pcal_index = target_names.index(field)
related_pcal = target_cal_map[pcal_index]

# ------- Parameters


gapfill = CAL_1GC_FILLGAPS
myuvrange = CAL_1GC_UVRANGE 
myspw = CAL_1GC_FREQRANGE

# ------- Setup names

tt = stamp()

ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
gptab = GAINTABLES+'/cal_1GC_'+myms+f'_{field}.Gp'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'
kcross  = GAINTABLES+'/cal_1GC_'+myms+'.KCROSS'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'

if pacal_name != '':
    gaintables = [ktab, bptab, ftab, dftab, kcross, xftab]
    gainfields = [related_pcal, bpcal_name, related_pcal, bpcal_name, pacal_name, pacal_name]
    interps = ['linear', 'linear','linear', 'nearest','nearest','nearest']
else:
    gaintables = [ktab, bptab, ftab, dftab]
    gainfields = [related_pcal, bpcal_name, related_pcal, bpcal_name]
    interps = ['linear', 'linear','linear', 'nearest']

# -------- Solve for the Gp Solutions using self-cal model

gaincal(vis=myms,
    field=field,
   # uvrange=myuvrange,
    caltable=gptab,
    gaintype='G',
    refant = str(ref_ant),
    solint=solint,
    minsnr=3,
    calmode='p',
    gaintable = gaintables,
    gainfield = gainfields,
    interp = interps)

applycal(vis=myms,
    #applymode='calflagstrict',
    field=field,
#    calwt=False,
    parang=True,
    gaintable = gaintables + [gptab],
    gainfield = gainfields + [field],
    interp = interps + ['linear'])

# ----- Save final flags for selfcal iteration

flagmanager(vis=myms,
        mode='delete',
        versionname=final_flag)

flagmanager(vis=myms,
        mode='save',
        versionname=final_flag)
