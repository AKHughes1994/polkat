# ian.heywood@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import subprocess


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


# ------- Setup names


tt = stamp()


ktab0 = GAINTABLES+'/cal_1GC_'+myms+'.K0'
bptab0 = GAINTABLES+'/cal_1GC_'+myms+'.B0'
gptab0 = GAINTABLES+'/cal_1GC_'+myms+'.Gp0'


ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gptab = GAINTABLES+'/cal_1GC_'+myms+'.Gp'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'


# ------- Set calibrator models



if primary_tag == '1934':
    setjy(vis=myms,
        field=bpcal_name,
        standard='Stevens-Reynolds 2016',
        scalebychan=True,
        usescratch=True)
    
    
elif primary_tag == '0408':
    bpcal_mod = CAL_1GC_0408_MODEL
    setjy(vis=myms,
        field=bpcal_name,
        standard='manual',
        fluxdensity=bpcal_mod[0],
        spix=bpcal_mod[1],
        reffreq=bpcal_mod[2],
        scalebychan=True,
        usescratch=True)


elif primary_tag == 'other':
    setjy(vis=myms,
        field=bpcal_name,
        standard='Perley-Butler 2013',
        scalebychan=True,
        usescratch=True)


#for i in range(0,len(pcals)):
#    pcal = pcals[i]
#    setjy(vis =myms,
#        field = pcal,
#        standard = 'manual',
#        fluxdensity = [1.0,0,0,0],
#        reffreq = '1000MHz',
#        usescratch = True)


# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------- STAGE 0 --------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #


# ------- K0 (primary)


gaincal(vis=myms,
    field=bpcal_name,
    #uvrange=myuvrange,
    #spw=myspw,
    caltable=ktab0,
    refant = str(ref_ant),
    gaintype = 'K',
    solint = 'inf',
    parang=False)


# ------- Gp0 (primary; apply K0)


gaincal(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    caltable=gptab0,
    gaintype='G',
    solint='inf',
    calmode='p',
    minsnr=5,
    gainfield=[bpcal_name],
    interp = ['nearest'],
    gaintable=[ktab0])


# ------- B0 (primary; apply K0, Gp0)


bandpass(vis=myms,
    field=bpcal_name, 
    uvrange=myuvrange,
    caltable=bptab0,
    refant = str(ref_ant),
    solint='inf',
    combine='',
    solnorm=False,
    minblperant=4,
    minsnr=3.0,
    bandtype='B',
    fillgaps=gapfill,
    parang=False,
    gainfield=[bpcal_name,bpcal_name],
    interp = ['nearest','nearest'],
    gaintable=[ktab0,gptab0])


flagdata(vis=bptab0,mode='tfcrop',datacolumn='CPARAM')
flagdata(vis=bptab0,mode='rflag',datacolumn='CPARAM')


# ------- Correct primary data with K0,B0,Gp0


applycal(vis=myms,
    gaintable=[ktab0,gptab0,bptab0],
#    applymode='calonly',
    field=bpcal_name,
#    calwt=False,
    parang=False,
    gainfield=[bpcal_name,bpcal_name,bpcal_name],
    interp = ['nearest','nearest','nearest'])


# ------- Flag primary on CORRECTED_DATA - MODEL_DATA


flagdata(vis=myms,
    mode='rflag',
    datacolumn='residual',
    field=bpcal_name)


flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='residual',
    field=bpcal_name)

if SAVE_FLAGS:
    flagmanager(vis=myms,
        mode='save',
        versionname='bpcal_residual_flags')


# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------- STAGE 1 --------------------------- #
# --------------------------------------------------------------- #
# --------------------------------------------------------------- #


# ------- K (primary; apply B0, Gp0)


gaincal(vis=myms,
    field=bpcal_name,
    caltable=ktab,
    refant = str(ref_ant),
    gaintype = 'K',
    solint = 'inf',
    parang=False,
    gaintable=[bptab0,gptab0],
    gainfield=[bpcal_name,bpcal_name],
    interp=['nearest','nearest'])


# ------- Gp (primary; apply K,B0)


gaincal(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    caltable=gptab,
    gaintype='G',
    solint='inf',
    calmode='p',
    minsnr=5,
    gainfield=[bpcal_name,bpcal_name],
    interp = ['nearest','nearest'],
    gaintable=[ktab,bptab0])


# ------- B (primary; apply K, Gp)


bandpass(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    caltable=bptab,
    refant = str(ref_ant),
    solint='inf',
    combine='',
    solnorm=False,
    minblperant=4,
    minsnr=3.0,
    bandtype='B',
    fillgaps=gapfill,
    parang=False,
    gainfield=[bpcal_name,bpcal_name],
    interp = ['nearest','nearest'],
    gaintable=[ktab,gptab])


flagdata(vis=bptab,mode='tfcrop',datacolumn='CPARAM')
flagdata(vis=bptab,mode='rflag',datacolumn='CPARAM')

# -------- Ga (primary; apply K, Gp, BP)

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    caltable=gatab,
    gaintype='G',
    solint='inf',
    calmode='a',
    minsnr=5,
    gainfield=[bpcal_name,bpcal_name,bpcal_name],
    interp = ['nearest','nearest','nearest'],
    gaintable=[ktab,gptab,bptab])


# ----- Fluxscale the amplitude table 

# fluxscale(vis=myms,
#   caltable = gatab,
#    fluxtable = ftab,
#    reference = bpcal_name,
#    append = False)

# ------- Apply parrallel hand tables to primary, flag, and plot visibilities

applycal(vis = myms,
        gaintable = [ktab, gptab, bptab, gatab],
#        applymode='calonly',
        field = bpcal_name,
#        calwt = False,
        parang = True,
        gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name],
        interp = ['nearest','nearest','nearest','nearest'])

flagdata(vis=myms,
    mode='rflag',
    datacolumn='residual',
    field=bpcal_name)


flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='residual',
    field=bpcal_name)

syscall = "shadems --xaxis FREQ,FREQ --yaxis CORRECTED_DATA:amp:YX,CORRECTED_DATA:amp:XY --colour-by ANTENNA1 --cnum 64 --suffix pre_leakage --dir {} --field {} {}".format(VISPLOTS, bpcal_name, myms)
subprocess.run([syscall],shell=True)

# -------- Df (primary; apply K, Gp, BP, Ga)

polcal(vis = myms,
    field = bpcal_name,
    uvrange = myuvrange,
    caltable = dftab,
    refant = str(ref_ant),
    solint = 'inf',
    poltype='Df',
    combine = 'scan',
    gaintable=[ktab, gptab, bptab, gatab],
    gainfield=[bpcal_name,bpcal_name,bpcal_name, bpcal_name],
    interp = ['nearest','nearest','nearest','nearest'],
    append = False)

flagdata(vis=dftab,mode='clip', clipminmax=[0,0.1], flagbackup=False, datacolumn='CPARAM')

# ------- Apply Leakage gain tables to primary, flag, and plot visibilities

applycal(vis = myms,
        gaintable = [ktab,gptab,bptab,gatab,dftab],
#        applymode='calonly',
        field = bpcal_name,
#        calwt = False,
        parang = True,
        gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name, bpcal_name],
        interp = ['nearest','nearest','nearest','nearest', 'nearest'])

# ------- Initial apply cal to pol cal

applycal(vis = myms,
        gaintable = [ktab,gptab,bptab,gatab,dftab],
#        applymode='calonly',
        field = pacal_name,
#        calwt = False,
        parang = True,
        gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name, bpcal_name],
        interp = ['linear','linear','linear','linear', 'linear'])

# ------- Initial apply to secondaries
for i in range(0,len(pcals)):

    pcal = pcals[i]

    applycal(vis = myms,
        gaintable = [ktab, gptab,bptab,gatab,dftab],
#        applymode='calonly',
        field = pcal,
#        calwt = False,
        parang = True,
        gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name, bpcal_name],
        interp = ['linear','linear','linear','linear', 'linear'])


syscall = "shadems --xaxis FREQ,FREQ --yaxis CORRECTED_DATA:amp:YX,CORRECTED_DATA:amp:XY --colour-by ANTENNA1 --cnum 64 --suffix post_leakage --dir {} --field {} {}".format(VISPLOTS, bpcal_name,myms)
subprocess.run([syscall],shell=True)

if SAVE_FLAGS:
    flagmanager(vis=myms,
        mode='save',
        versionname='bpcal_final_flags')


# ---- Flag the polarization calibrator
flagdata(vis=myms,
    mode='rflag',
    datacolumn='corrected',
    field=pacal_name, flagbackup=False)

flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='corrected',
    field=pacal_name, flagbackup=False)

# --- Flag the secondaries
for i in range(0,len(pcals)):

    pcal = pcals[i]

    flagdata(vis=myms,
        mode='rflag',
        datacolumn='corrected',
        field=pacal, flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='corrected',
        field=pcal, flagbackup=False)


if SAVE_FLAGS:
    flagmanager(vis=myms,
        mode='save',
        versionname='parr_hand_flags')
