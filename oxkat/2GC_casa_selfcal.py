# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import datetime
import sys
import time


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())


def stamp():
    now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
    return now


# Calibration mode: P = phase, K = delay, A = amp+phase
CAL_MODE = 'P'

# Rotate the sky-frame model into the feed frame during the solve (solutions stay in the feed frame)
SOLVE_PARANG = False

# Rotate the corrected data back to the sky frame when applying (use on the final solve only)
APPLY_PARANG = False

# Solutions below this S/N are flagged
MINSNR = 5.0

# Per-mode gaincal settings
MODE_SETTINGS = {
    'P': {'gaintype': 'G', 'calmode': 'p',  'solint': 'int', 'suffix': 'GP_2GC'},
    'K': {'gaintype': 'K', 'calmode': 'p',  'solint': '64s', 'suffix': 'K_2GC'},
    'A': {'gaintype': 'G', 'calmode': 'ap', 'solint': 'inf', 'suffix': 'GA_2GC'},
}


def str2bool(s):
    return s.strip().lower() in ('true', 't', 'yes', 'y', '1')


# Read ms=<name>, mode=P|K|A, solve_parang=True|False and parang=True|False from the command line, overriding the defaults above
for item in sys.argv:
    parts = item.split('=')
    if parts[0] == 'ms':
        myms = parts[1]
    if parts[0] == 'mode':
        CAL_MODE = parts[1].upper()
    if parts[0] == 'solve_parang':
        SOLVE_PARANG = str2bool(parts[1])
    if parts[0] == 'parang':
        APPLY_PARANG = str2bool(parts[1])


# Stop early if the requested mode is not P, K or A
if CAL_MODE not in MODE_SETTINGS:
    raise ValueError('mode must be one of P, K or A, got: '+CAL_MODE)

# Pick the gaincal settings for the chosen mode
settings = MODE_SETTINGS[CAL_MODE]

# Baselines used for the solve: the imaging UV cut (WSC_MINUVL/WSC_MAXUVL) if set, else all
myuvrange = ''
if WSC_MINUVL != '' and WSC_MAXUVL != '':
    myuvrange = str(WSC_MINUVL)+'~'+str(WSC_MAXUVL)+'lambda'
elif WSC_MINUVL != '':
    myuvrange = '>'+str(WSC_MINUVL)+'lambda'
elif WSC_MAXUVL != '':
    myuvrange = '<'+str(WSC_MAXUVL)+'lambda'


# Timestamped output table in GAINTABLES, suffix set by the mode
gtab = GAINTABLES+'/cal_'+myms+'_'+stamp()+'.'+settings['suffix']


# Summary of the parameters used for this run
print('2GC self-calibration summary', flush=True)
print('  MS:          '+myms, flush=True)
print('  Mode:        '+CAL_MODE, flush=True)
print('  Gaintype:    '+settings['gaintype'], flush=True)
print('  Calmode:     '+settings['calmode'], flush=True)
print('  Solint:      '+settings['solint'], flush=True)
print('  Min S/N:     '+str(MINSNR), flush=True)
print('  UV range:    '+(myuvrange if myuvrange != '' else 'all baselines'), flush=True)
print('  Solve parang:'+str(SOLVE_PARANG), flush=True)
print('  Apply parang:'+str(APPLY_PARANG), flush=True)
print('  Refant:      '+str(ref_ant), flush=True)
print('  Caltable:    '+gtab, flush=True)
print('  Interp:      linear', flush=True)


# Solve for the chosen calibration on field 0
gaincal(vis=myms,
    field='0',
    uvrange=myuvrange,                  # baselines to include
    caltable=gtab,                      # solutions are written here
    refant = str(ref_ant),              # reference antenna from project_info.json
    gaintype=settings['gaintype'],      # G (gain) or K (delay)
    solint=settings['solint'],          # solution interval for this mode
    solnorm=False,                      # keep absolute amplitudes
    minsnr=MINSNR,                      # flag solutions below this S/N
    calmode=settings['calmode'],        # p = phase only, ap = amp + phase
    parang=SOLVE_PARANG,                # rotate the model by the parallactic angle
    gaintable=[],                       # no prior calibration applied
    gainfield=[],
    interp=[],
    append=False)


# Apply the new solutions to the data, writing to CORRECTED_DATA
applycal(vis=myms,
    gaintable=[gtab],
    field='0',
    parang=APPLY_PARANG,                # derotate the corrected data to the sky frame
    gainfield='0',
    interp = ['linear'])


# Optional: recompute data weights after calibration
# statwt(vis=myms,
#     field='0')


