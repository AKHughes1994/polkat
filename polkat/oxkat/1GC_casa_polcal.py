# ian.heywood@physics.ox.ac.uk


import glob,os,datetime, subprocess
import shutil
import time

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
which_ms = CAL_1GC_WHICH_MODE

if SAVE_FLAGS:
    flagmanager(vis=myms,
        mode='save',
        versionname='before_polcal')

# ------- Setup names

ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gptab = GAINTABLES+'/cal_1GC_'+myms+'.Gp'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'
tatab = GAINTABLES+'/cal_1GC_'+myms+'.Ta'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'
kcross = GAINTABLES+'/cal_1GC_'+myms+'.Kcross'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'

# ---- Solve for K for PA Cal

gaincal(vis= myms,
        field = pacal_name,
        caltable = ktab,
        refant = str(ref_ant),
        gaintype = 'K',
        solint = 'inf',
        append = True)

# ---- Iterate through Phase cals

for i in range(0,len(pcals)):

    pcal = pcals[i]

    # --- Solve for K for secondaries


    gaincal(vis= myms,
        field = pcal,
        caltable = ktab,
        refant = str(ref_ant),
        gaintype = 'K',
        solint = 'inf',
        append = True)


# ---- Solve for Gp for PA Cal

gaincal(vis=myms,
    field=pacal_name,
    uvrange=myuvrange,
    caltable=gptab,
    gaintype='G',
    solint='int',
    calmode='p',
    gaintable=[ktab],
    gainfield=[pacal_name],
    interp = ['nearest'],
    append = True)


# ---- Iterate through phase cals

for i in range(0,len(pcals)):

    pcal = pcals[i]

    # --- Solve for Gp for secondaries

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        caltable=gptab,
        gaintype='G',
        solint='int',
        calmode='p',
        gaintable=[ktab],
        gainfield=[pcal],
        interp = ['nearest'],
        append = True)

    # --- Solve for a long-solution interval for secondaries which will be applied to the source field

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        caltable=gptab.replace('.Gp', '_{}.Gp'.format(pcal)),
        gaintype='G',
        solint='inf',
        calmode='p',
        gaintable=[ktab],
        gainfield=[pcal],
        interp = ['nearest'],
        append = False)

# ---- Solve for Ta for PA Cal

gaincal(vis=myms,
    field=pacal_name,
    uvrange=myuvrange,
    caltable=tatab,
    gaintype='T',
    solint='inf',
    calmode='a',
    gaintable=[ktab,gptab,bptab,gatab],
    gainfield=[pacal_name, pacal_name, bpcal_name,  bpcal_name],
    interp = ['nearest', 'nearest', 'linear', 'linear'],
    append = False)


# ---- Iterate through phase cals

for i in range(0,len(pcals)):

    pcal = pcals[i]

    # --- Solve for Ta for secondaries

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        caltable=tatab,
        gaintype='T',
        solint='inf',
        calmode='a',
        gaintable=[ktab,gptab,bptab,gatab],
        gainfield=[pcal, pcal, bpcal_name,  bpcal_name],
        interp = ['nearest', 'nearest', 'linear', 'linear'],
        append = True)

# ---- Solve for KCROSS

gaincal(vis = myms,
    field = pacal_name,
    uvrange = myuvrange,
    caltable = kcross,
    refant = str(ref_ant),
    solint = 'inf,100MHz',
    gaintype='KCROSS',
    parang = True,
    combine = '',
    gaintable=[ktab,gptab,bptab,gatab,tatab,dftab],
    gainfield=[pacal_name,pacal_name,bpcal_name,bpcal_name,pacal_name,bpcal_name],
    interp = ['nearest','nearest','linear','linear','nearest','linear'],
    append = False)

# ---- Solve for Xf

polcal(vis = myms,
    field = pacal_name,
    uvrange = myuvrange,
    caltable = xftab,
    refant = str(ref_ant),
    solint = 'inf,100MHz',
    poltype='Xf',
    combine = '',
    gaintable=[ktab,gptab,bptab,gatab,tatab,dftab,kcross],
    gainfield=[pacal_name,pacal_name,bpcal_name,bpcal_name,pacal_name,bpcal_name,pacal_name],
    interp = ['nearest','nearest','linear','linear','nearest','linear','nearest'],
    append = False)

#--------------------------------------------------------------------------------#

# ------- Apply final tables to primary

applycal(vis = myms,
#        applymode='calonly',
        field = bpcal_name,
        calwt = False,
        parang = True,
        gaintable = [ktab,gptab,bptab,gatab, dftab, kcross, xftab],
        gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name, bpcal_name, pacal_name, pacal_name],
        interp = ['nearest','nearest','nearest','nearest', 'nearest','linear', 'linear'])

# ------- Apply final tables to pol cal

applycal(vis = myms,
#        applymode='calonly',
        field = pacal_name,
        calwt = False,
        parang = True,
        gaintable = [ktab,gptab,bptab,gatab, tatab, dftab, kcross, xftab],
        gainfield = [pacal_name,pacal_name, bpcal_name, bpcal_name, pacal_name, bpcal_name, pacal_name, pacal_name],
        interp = ['nearest','nearest','linear','linear', 'nearest','linear','nearest', 'nearest'])

# ------- Apply final tables to secondaries


for i in range(0,len(pcals)):


    pcal = pcals[i]


    # --- Correct secondaries with K, Gp, B, Ga, Ta, Df, Kcross, Xf


    applycal(vis = myms,
#            applymode='calonly',
            field = pcal,
            calwt = False,
            parang = True,
            gaintable = [ktab,gptab,bptab,gatab, tatab, dftab, kcross, xftab],
            gainfield = [pcal,pcal, bpcal_name, bpcal_name, pcal, bpcal_name, pacal_name, pacal_name],
            interp = ['nearest','nearest','linear','linear', 'nearest','linear','linear', 'linear'])


    # ------- Apply final tables to targets


for i in range(0,len(targets)):


    target = targets[i]
    related_pcal = target_cal_map[i]


    # --- Correct secondaries with K, Gp, B, Ga, Ta, Df, Kcross, Xf


    applycal(vis=myms,
#            applymode='calonly',
            field=target,
            calwt=False,
            parang=True,
            gaintable = [ktab,gptab.replace('.Gp', '_{}.Gp'.format(related_pcal)),bptab,gatab, tatab, dftab, kcross, xftab],
            gainfield = [related_pcal,related_pcal, bpcal_name, bpcal_name, related_pcal, bpcal_name, pacal_name, pacal_name],
            interp = ['linear','linear','linear','linear', 'linear','linear','linear', 'linear'])

    flagdata(vis=myms, field=target, mode='tfcrop', flagbackup=False)
    flagdata(vis=myms, field=target, mode='rflag', flagbackup=False)


if SAVE_FLAGS:
    flagmanager(vis=myms,
        mode='save',
        versionname='after_polcal')





















