# ian.heywood@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import subprocess
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


# ------- Setup names


tt = stamp()


# Initial Calibrators tables for the primary
ktab0 = GAINTABLES+'/cal_1GC_'+myms+'.K0'
bptab0 = GAINTABLES+'/cal_1GC_'+myms+'.B0'
gptab0 = GAINTABLES+'/cal_1GC_'+myms+'.Gp0'
gatab0 = GAINTABLES+'/cal_1GC_'+myms+'.Ga0'
ftab0 = GAINTABLES+'/cal_1GC_'+myms+'.F0'
dftab0  = GAINTABLES+'/cal_1GC_'+myms+'.Df0'

ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gptab = GAINTABLES+'/cal_1GC_'+myms+'.Gp'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'


kcross  = GAINTABLES+'/cal_1GC_'+myms+'.KCROSS'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'

# Restore the auto_cal flag version
flagmanager(vis=myms,
        mode='restore',
        versionname='autoflag_cals_data')


# ------- Set calibrator models

if primary_tag == '1934':
    #setjy(vis=myms,
    #   field=bpcal_name,
    #  standard='Stevens-Reynolds 2016',
    # scalebychan=True,
    #usescratch=True)

    syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS1934.LBand.wsclean.cat.txt"
    subprocess.run([syscall],shell=True)
    
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

for i in range(0,len(pcals)):
    pcal = pcals[i]
    setjy(vis =myms,
        field = pcal,
        standard = 'manual',
        fluxdensity = [1.0,0,0,0],
        reffreq = '1000MHz',
        usescratch = True)

if pacal_name != ''
    # Initialize the Polarization Calibrator model to have 0.0 Stokes V and Non-zero Stokes U
    setjy(vis=myms,
        field=pacal_name,
        standard='manual',
        fluxdensity = POLANG_MOD,
        usescratch=True)

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #
# --------------------------- STAGE 0 ----------------------- #
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
    solint = 'inf')


# ------- Gp0 (primary; apply K0)


gaincal(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    caltable=gptab0,
    refant = str(ref_ant),
    gaintype='G',
    solint='inf',
    calmode='p',
    minsnr=5,
    gainfield=[bpcal_name],
    interp = ['linear'],
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
    gainfield=[bpcal_name,bpcal_name],
    interp = ['linear','linear'],
    gaintable=[ktab0,gptab0])


flagdata(vis=bptab0,mode='tfcrop',datacolumn='CPARAM', flagbackup=False)
flagdata(vis=bptab0,mode='rflag',datacolumn='CPARAM', flagbackup=False)

# ------- Ga0 (primary; apply K0, Gp0, BP0) -- Type T


gaincal(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    spw = myspw,
    caltable=gatab0,
    refant = str(ref_ant),
    gaintype='T',
    solint='inf',
    calmode='a',
    minsnr=3,
    gainfield=[bpcal_name,bpcal_name, bpcal_name],
    interp = ['linear','linear', 'linear'],
    gaintable=[ktab0,gptab0, bptab0])


# -------- Solve for Df0 (apply K0, Gp0, Bp0, Ga0)

polcal(vis = myms,
    field = bpcal_name,
    uvrange = myuvrange,
    caltable = dftab0,
    refant = str(ref_ant),
    solint = 'inf',
    poltype='Df',
    combine = 'scan',
    gaintable=[ktab0,gptab0, bptab0,gatab0],
    gainfield=[bpcal_name,bpcal_name,bpcal_name, bpcal_name],
    interp = ['linear','linear','linear','linear'],
    append = False)

flagdata(vis=dftab0,mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')

# ------- Correct primary data with K0,B0,Gp0,gatab0,dftab0


applycal(vis=myms,
    gaintable=[ktab0,gptab0,bptab0, gatab0, dftab0],
    #applymode='calflagstrict',
    field=bpcal_name,
    #calwt=False,
    parang=True,
    gainfield=[bpcal_name,bpcal_name,bpcal_name, bpcal_name, bpcal_name],
    interp = ['linear','linear','linear', 'linear', 'nearest'], flagbackup=False)


# ------- Flag primary on CORRECTED_DATA - MODEL_DATA

flagdata(vis=myms,
    mode='rflag',
    datacolumn='residual',
    field=bpcal_name, 
    flagbackup=False) 

flagdata(vis=myms,
    mode='tfcrop',
    datacolumn='residual',
    field=bpcal_name,
    flagbackup=False) 

flagmanager(vis=myms,
        mode='delete',
        versionname='bpcal_residual_flags')

flagmanager(vis=myms,
        mode='save',
        versionname='bpcal_residual_flags')


# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #
# --------------------------- Working Table (Primary)  --------------------------- #
# ---------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------- #


# ------- K (primary; apply B0, Gp0, Ga0, Df0)


gaincal(vis=myms,
    field=bpcal_name,
    caltable=ktab,
    refant = str(ref_ant),
    gaintype = 'K',
    solint = 'inf',
    gaintable=[bptab0,gptab0, gatab0, dftab0],
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['linear','linear', 'linear', 'nearest'])


# ------- Gp (primary; apply K,B0, Ga0, Df0)


gaincal(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    caltable=gptab,
    gaintype='G',
    refant = str(ref_ant),
    solint='inf',
    calmode='p',
    minsnr=5,
    gaintable=[bptab0,ktab, gatab0, dftab0],
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['linear','linear', 'linear', 'nearest'])


# ------- B (primary; apply K, Gp, Ga0, Df0)


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
    gaintable=[ktab, gptab, gatab0, dftab0],
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['linear','linear', 'linear', 'nearest'])


flagdata(vis=bptab,mode='tfcrop',datacolumn='CPARAM' , flagbackup=False)
flagdata(vis=bptab,mode='rflag',datacolumn='CPARAM' , flagbackup=False)

# -------- Ga (primary; apply K, Gp, BP, Df0) -- Gaintype 'T'

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=myuvrange,
    spw = myspw,
    caltable=gatab,
    refant = str(ref_ant),
    gaintype='T',
    solint='inf',
    calmode='a',
    minsnr=3,
    gaintable=[ktab, gptab, bptab, dftab0],
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['linear','linear', 'linear', 'nearest'])


# -------- Solve for Df (primary; apply K, Gp, BP, Ga)

polcal(vis = myms,
    field = bpcal_name,
    uvrange = myuvrange,
    caltable = dftab,
    refant = str(ref_ant),
    solint = 'inf',
    poltype='Df',
    combine = 'scan',
    gaintable=[ktab, gptab, bptab, gatab],
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['linear','linear', 'linear', 'linear'])

flagdata(vis=dftab, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')


# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- Initial Table (Secondary + Pol. Cal.)  ---------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #

if pacal_name != '':

    # ------- Gp0 (polcal; apply Bp, Df, K (primary))


    gaincal(vis = myms,
        field=pacal_name,
        uvrange=myuvrange,
        # spw = myspw,
        caltable=gptab0,
        refant = str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'nearest'],
        append=True)


    # ------- Ga0 (polcal; apply Bp, Df, K (primary) Gp0 (polcal))


    gaincal(vis = myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw = myspw,
        caltable=gatab0,
        refant = str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, gptab0, bptab, dftab],
        gainfield=[bpcal_name, pacal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'nearest', 'linear', 'nearest'],
        append=True)

    # ------- K0 (polcal; apply Bp, Df (primary), Gp0, Ga0 (polcal))


    gaincal(vis= myms,
        field = pacal_name,
        #   uvrange = myuvrange,
        #   spw=myspw,
        caltable = ktab0,
        refant = str(ref_ant),
        gaintype = 'K',
        solint='inf',
        gaintable=[gptab0, gatab0, bptab, dftab],
        gainfield=[pacal_name,pacal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'nearest', 'linear', 'nearest'],
        append=True)

# ----- Loop over secondaries

for i in range(0,len(pcals)):

    pcal = pcals[i]
    
    # ------- Gp0 (pcal; apply Bp, Df, K (primary))


    gaincal(vis = myms,
        field=pcal,
        uvrange=myuvrange,
        # spw = myspw,
        caltable=gptab0,
        refant = str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'nearest'],
        append=True)


    # ------- Ga0 (pcal; apply Bp, Df, K (primary) Gp0 (pcal))

    gaincal(vis = myms,
        field=pcal,
        uvrange=myuvrange,
        spw = myspw,
        caltable=gatab0,
        refant = str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, gptab0, bptab, dftab],
        gainfield=[bpcal_name, pcal, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear', 'nearest'],
        append=True)

    # ------- K0 (pcal; apply Bp, Df (primary), Gp0, Ga0 (pcal))


    gaincal(vis= myms,
        field = pcal,
        #   uvrange = myuvrange,
        #   spw=myspw,
        caltable = ktab0,
        refant = str(ref_ant),
        gaintype = 'K',
        solint='inf',
        gaintable=[gptab0, gatab0, bptab, dftab],
        gainfield=[pcal,pcal, bpcal_name, bpcal_name],
        interp=['linear', 'linear', 'linear', 'nearest'],
        append=True)

# --- Apply fluxscaling to Ga0

fluxscale(vis=myms,
    caltable = gatab0,
    fluxtable = ftab0,
    reference = bpcal_name,
    append = False,
    transfer = '')

if pacal_name != '':
    # -------- Applycal (polcal; Bp, Df (primary), Ga0, K0, Gp0 (polcal)) and Flag
    
    applycal(vis=myms,
        gaintable=[ktab0,gptab0, ftab0, bptab, dftab],
        #applymode='calflagstrict',
        field=pacal_name,
        #calwt=False,
        parang=True,
        gainfield=[pacal_name, pacal_name, pacal_name, bpcal_name, bpcal_name],
        interp = ['nearest','nearest','nearest', 'linear', 'nearest'], 
        flagbackup=False)

    flagdata(vis=myms,
        mode='rflag',
        datacolumn='corrected',
        field=pacal_name, 
        flagbackup=False) 

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='corrected',
        field=pacal_name,
        flagbackup=False) 


# ----- Loop over secondaries

for i in range(0,len(pcals)):

    pcal = pcals[i]


    # -------- Applycal (pcal; Bp, Df (primary), Ga0, K0, Gp0 (polcal)) and Flag

    applycal(vis=myms,
        gaintable=[ktab0,gptab0, ftab0, bptab, dftab],
        #applymode='calflagstrict',
        field=pcal,
        #calwt=False,
        parang=True,
        gainfield=[pcal, pcal, pcal, bpcal_name, bpcal_name],
        interp = ['linear','linear','linear', 'linear', 'nearest'], 
        flagbackup=False)

    flagdata(vis=myms,
        mode='rflag',
        datacolumn='corrected',
        field=pcal,
        flagbackup=False) 

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='corrected',
        field=pcal,
        flagbackup=False) 

# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- Working Table (Secondary + Pol. Cal.)  ------------------------ #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #

if pacal_name != '':   
 
    # ------- Gp (polcal; apply Bp, Df, K (primary))


    gaincal(vis = myms,
        field=pacal_name,
        uvrange=myuvrange,
        # spw = myspw,
        caltable=gptab,
        refant = str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'nearest'],
        append=True)


    # ------- Ga (polcal; apply Bp, Df, K (primary) Gp (polcal))


    gaincal(vis = myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw = myspw,
        caltable=gatab,
        refant = str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, gptab, bptab, dftab],
        gainfield=[bpcal_name, pacal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'nearest', 'linear', 'nearest'],
        append=True)

    # ------- K (polcal; apply Bp, Df (primary), Gp, Ga (polcal))
    

    gaincal(vis= myms,
        field = pacal_name,
        #   uvrange = myuvrange,
        #   spw=myspw,
        caltable = ktab,
        refant = str(ref_ant),
        gaintype = 'K',
        solint='inf',
        gaintable=[gptab, gatab, bptab, dftab],
        gainfield=[pacal_name,pacal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'nearest', 'linear', 'nearest'],
        append=True)


    # ------- KCROSS (polcal; apply Bp, Df (primary), Ga, K, Gp (polcal))

    gaincal(vis = myms,
        field = pacal_name,
    #    uvrange = myuvrange,
        caltable = kcross,
        #   spw=myspw,
        refant = str(ref_ant),
        solint = 'inf',
        gaintype='KCROSS',
        parang = True,
        gaintable=[ktab,gptab,bptab,gatab,dftab],
        gainfield=[pacal_name, pacal_name,bpcal_name, pacal_name, bpcal_name],
        interp = ['nearest','nearest','linear','nearest','nearest'],
        append = False)

    # -------- Xf (polcal; apply Bp, Df (primary), Ga, Gp, K, KCROSS (polcal))

    polcal(vis = myms,
        field = pacal_name,
        uvrange = myuvrange,
        caltable = xftab,
        refant = str(ref_ant),
        solint = 'inf,64ch', # 16-channels-across in frequency
        poltype='Xf',
        combine = '',
        gaintable=[ktab,gptab,bptab,gatab,dftab, kcross],
        gainfield=[pacal_name,pacal_name,bpcal_name, pacal_name, bpcal_name, pacal_name],
        interp = ['nearest','nearest','linear','nearest','nearest', 'nearest'],
        append = False)

for i in range(0,len(pcals)):

    pcal = pcals[i]

    # ------- Gp0 (pcal; apply Bp, Df, K (primary))

    gaincal(vis = myms,
        field=pcal,
        uvrange=myuvrange,
        # spw = myspw,
        caltable=gptab,
        refant = str(ref_ant),
        gaintype='G',
        solint='inf',
        calmode='p',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'nearest'],
        append=True)


    # ------- Ga (pcal; apply Bp, Df, K (primary) Gp (pcal))

    gaincal(vis = myms,
        field=pcal,
        uvrange=myuvrange,
        spw = myspw,
        caltable=gatab,
        refant = str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='a',
        minsnr=3,
        gaintable=[ktab, gptab0, bptab, dftab],
        gainfield=[bpcal_name, pcal, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear', 'nearest'],
        append=True)

    # ------- K (pcal; apply Bp, Df (primary), Gp, Ga (pcal))


    gaincal(vis= myms,
        field = pcal,
        #   uvrange = myuvrange,
        #   spw=myspw,
        caltable = ktab,
        refant = str(ref_ant),
        gaintype = 'K',
        solint='inf',
        gaintable=[gptab, gatab, bptab, dftab],
        gainfield=[pcal, pcal, bpcal_name, bpcal_name],
        interp=['linear', 'linear', 'linear', 'nearest'],
        append=True)

# --- Apply fluxscaling to Ga
fluxscale(vis=myms,
    caltable = gatab,
    fluxtable = ftab,
    reference = bpcal_name,
    append = False,
    transfer = '')

# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
# --------------------------- Applycal (All Fields)  ----------------------- #
# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #

# ----- If no polarization angle calibrator apply subset of tables and kill script

if pacal_name == '':   

    # ------- BPCAL

    applycal(vis = myms,
        gaintable = [ktab,gptab,bptab,ftab,dftab],
 #       applymode='calflagstrict',
        field = bpcal_name,
        #calwt = False,
        parang = True,
        gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name, bpcal_name],
        interp = ['linear','linear','linear','linear','nearest'],
        flagbackup=False)

    # ------- Secondaries 

    for i in range(0,len(pcals)):

        pcal = pcals[i]
    
        applycal(vis = myms,
            gaintable = [ktab,gptab,bptab,ftab,dftab],
            # applymode='calflagstrict',
            field = pcal,
            #calwt = False,
            parang = True,
            gainfield = [pcal,pcal, bpcal_name, pcal, bpcal_name],
            interp = ['linear','linear','linear','linear','nearest'],
            flagbackup=False)

    # ------- Targets 
    for i in range(0,len(targets)):

        target = targets[i]
        related_pcal = target_cal_map[i]

        applycal(vis=myms,
                #applymode='calflagstrict',
                gaintable = [ktab,gptab,bptab,ftab,dftab],
                field=target,
                #calwt=False,
                parang=True,
                gainfield = [related_pcal, related_pcal, bpcal_name, related_pcal, bpcal_name],
                interp = ['linear','linear','linear','linear','nearest'],
                flagbackup=False)

        # Flag target
        flagdata(vis=myms,
            mode='rflag',
            datacolumn='corrected',
            field=target, flagbackup=False)

        flagdata(vis=myms,
            mode='tfcrop',
            datacolumn='corrected',
            field=target, flagbackup=False)

    # ---- Save flags

    flagmanager(vis=myms,
        mode='delete',
        versionname='1GC_flags')

    flagmanager(vis=myms,
        mode='save',
        versionname='1GC_flags')

    sys.exit('Ending Early! No polarization angle calibrator')

# -------- Full polarization 

# ------- BPCAL

applycal(vis = myms,
        gaintable = [ktab,gptab,bptab,ftab,dftab, kcross, xftab],
 #       applymode='calflagstrict',
        field = bpcal_name,
        #calwt = False,
        parang = True,
        gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name, bpcal_name, pacal_name, pacal_name],
        interp = ['linear','linear','linear','linear','nearest','nearest', 'nearest'],
        flagbackup=False)

# ------- PACAL

applycal(vis = myms,
        gaintable = [ktab,gptab,bptab,ftab,dftab, kcross, xftab],
 #       applymode='calflagstrict',
        field = pacal_name,
        #calwt = False,
        parang = True,
        gainfield = [pacal_name,pacal_name, bpcal_name, pacal_name, bpcal_name, pacal_name, pacal_name],
        interp = ['nearest','nearest','linear','nearest','nearest','nearest', 'nearest'],
        flagbackup=False)

# ------- Secondaries

for i in range(0,len(pcals)):

    pcal = pcals[i]

    applycal(vis = myms,
        gaintable = [ktab,gptab,bptab,ftab,dftab, kcross, xftab],
        # applymode='calflagstrict',
        field = pcal,
        #calwt = False,
        parang = True,
        gainfield = [pcal,pcal, bpcal_name, pcal, bpcal_name, pacal_name, pacal_name],
        interp = ['linear','linear','linear','linear','nearest','nearest', 'nearest'],
        flagbackup=False)

# ------- Targets 
for i in range(0,len(targets)):

    target = targets[i]
    related_pcal = target_cal_map[i]

    applycal(vis=myms,
                #applymode='calflagstrict',
                gaintable = [ktab,gptab,bptab,ftab,dftab, kcross, xftab],
                field=target,
                #calwt=False,
                parang=True,
                gainfield = [related_pcal, related_pcal, bpcal_name, related_pcal, bpcal_name, pacal_name, pacal_name],
                interp = ['linear','linear','linear','linear','nearest','nearest', 'nearest'],
                flagbackup=False)

    # Flag target
    flagdata(vis=myms,
        mode='rflag',
        datacolumn='corrected',
        field=target, flagbackup=False)

    flagdata(vis=myms,
        mode='tfcrop',
        datacolumn='corrected',
        field=target, flagbackup=False)

# ---- Save flags

flagmanager(vis=myms,
    mode='delete',
    versionname='1GC_flags')

flagmanager(vis=myms,
    mode='save',
    versionname='1GC_flags')
