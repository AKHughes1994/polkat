# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import subprocess
import sys
import numpy as np

# Flush immediately for better logging
import functools
print = functools.partial(print, flush=True)

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    targets = user_targets
    pcal_names = user_pcals
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

# Calibrator tables
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

# Remove ftabs if they exist to prevent code from breaking
if os.path.isdir(ftab0):
    print(f"Removing: {ftab0}")
    shutil.rmtree(ftab0)

if os.path.isdir(ftab):
    print(f"Removing: {ftab}")
    shutil.rmtree(ftab)

# ------- Set BP calibrator models

if primary_tag == '1934':
    
    # MeerKAT specific crystalball models for 1939 from B.Hugo: https://archive-gw-1.kat.ac.za/public/repository/10.48479/hhhy-4r55/index.htmlV
    if band == 'L':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS1934.LBand.wsclean.cat.txt"
        #subprocess.run([syscall],shell=True)
        pass

    elif band == 'UHF':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS1934.UBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)

    else:
        setjy(vis=myms,
            field=bpcal_name,
            standard='Stevens-Reynolds 2016',
            scalebychan=True,
            usescratch=True)

elif primary_tag == '0408':

    # MeerKAT specific crystalball models for 0408 from B.Hugo: https://archive-gw-1.kat.ac.za/public/repository/10.48479/ez63-vx81/index.html
    if band == 'L':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS0407.LBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)

    elif band == 'UHF':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS0407.UBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)
    
    else:
        # OXKAT Version that uses config.py input
        #bpcal_mod = CAL_1GC_0408_MODEL
        #setjy(vis=myms,
        #    field=bpcal_name,
        #    standard='manual',
        #    fluxdensity=bpcal_mod[0],
        #    spix=bpcal_mod[1],
        #    reffreq=bpcal_mod[2],
        #    scalebychan=True,
        #    usescratch=True)

        # Recommendation from SARAO, i.e., https://skaafrica.atlassian.net/wiki/spaces/ESDKB/pages/1481408634/Flux+and+bandpass+calibration
        # See script in tools/SARAO_0408_model.py to see where CASA parameters come from
        setjy(vis=myms,
            field=bpcal_name,
            standard='manual',
            fluxdensity=[6.9862, 0.0, 0.0, 0.0],
            spix=[-1.2897, -0.2353, 0.0861],
            reffreq='2.7GHz',
            scalebychan=True,
            usescratch=True)


elif primary_tag == 'other':
    setjy(vis=myms,
        field=bpcal_name,
        standard='Perley-Butler 2013',
        scalebychan=True,
        usescratch=True)

for i in range(0,len(pcal_names)):
    pcal = pcal_names[i]
    if pcal != bpcal_name:
        setjy(vis =myms,
            field = pcal,
            standard = 'manual',
            fluxdensity = [1.0,0,0,0],
            reffreq = '1000MHz',
            usescratch = True)

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
    interp = ['linear','linear','linear', 'linear', 'linear'], flagbackup=False)


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
    interp=['linear','linear', 'linear', 'linear'])


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
    interp=['linear','linear', 'linear', 'linear'])


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
    interp=['linear','linear', 'linear', 'linear'])


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
    interp=['linear','linear', 'linear', 'linear'])


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
        interp=['linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

# ----- Loop over secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

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
        interp=['linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

# --- Apply fluxscaling to Ga0 but only if there are calibration fields other than the primary
if len([pcal for pcal in pcal_names if pcal != bpcal_name]) > 0 or pacal_name != '':
    
    # --- Apply fluxscaling to Ga0
    fluxscale(vis=myms,
        caltable = gatab0,
        fluxtable = ftab0,
        reference = bpcal_name,
        append = False,
        transfer = '')

# If there is no need to apply flux scaling, we can set ftab0 to gatab0 as the only calibrator is the primary    
else:
    ftab0 = gatab0

if pacal_name != '':
    # -------- Applycal (polcal; Bp, Df (primary), Ga0, K0, Gp0 (polcal)) and Flag
    
    applycal(vis=myms,
        gaintable=[ktab0,gptab0, ftab0, bptab, dftab],
        #applymode='calflagstrict',
        field=pacal_name,
        #calwt=False,
        parang=True,
        gainfield=[pacal_name, pacal_name, pacal_name, bpcal_name, bpcal_name],
        interp = ['linear','linear','linear', 'linear', 'linear'], 
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

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    # -------- Applycal (pcal; Bp, Df (primary), Ga0, K0, Gp0 (polcal)) and Flag

    applycal(vis=myms,
        gaintable=[ktab0,gptab0, ftab0, bptab, dftab],
        #applymode='calflagstrict',
        field=pcal,
        #calwt=False,
        parang=True,
        gainfield=[pcal, pcal, pcal, bpcal_name, bpcal_name],
        interp = ['linear','linear','linear', 'linear', 'linear'], 
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
        interp=['linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

            
for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

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
        interp=['linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
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
        interp=['linear', 'linear', 'linear', 'linear'],
        append=True)

# --- Apply fluxscaling to Ga but only if there are calibration fields other than the primary
if len([pcal for pcal in pcal_names if pcal != bpcal_name]) > 0 or pacal_name != '':
    
    # --- Apply fluxscaling to Ga
    fluxscale(vis=myms,
        caltable = gatab,
        fluxtable = ftab,
        reference = bpcal_name,
        append = False,
        transfer = '')

# If there is no need to apply flux scaling, we can set ftab0 to gatab0 as the only calibrator is the primary    
else:
    ftab = gatab

# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- CROSS-HAND Tables (PA CAL)  ------------------------------------------------ #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #


if pacal_name != '':
 
    # ------- Set PA calibrator models

    setjy(vis=myms,
        field=pacal_name,
        standard='manual',
        fluxdensity = POLANG_MOD,
        usescratch=True)

    # -------- Solve for Cross-hand phase terms
    if XF_MODE not in ['casa', 'manual', 'auto']:
        print(f"Cross-hand phase mode ({XF_MODE}) not a valid option, defauling to 'auto'")
        XF_MODE = 'auto'

    manual_XF = False
    if XF_MODE == 'casa' or XF_MODE == 'auto':

        print(f"Cross-hand phase mode if {XF_MODE}; solving CASA based solutions")
        
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
            interp = ['linear','linear','linear','linear','linear'],
            append = False)

        # -------- Xf (polcal; apply Bp, Df (primary), Ga, Gp, K, KCROSS (polcal))

        polcal(vis = myms,
            field = pacal_name,
            uvrange = myuvrange,
            caltable = xftab,
            refant = str(ref_ant),
            solint = f'inf,{XF_CHANINT}ch',
            poltype='Xf',
            combine = '',
            gaintable=[ktab,gptab,bptab,gatab,dftab, kcross],
            gainfield=[pacal_name,pacal_name,bpcal_name, pacal_name, bpcal_name, pacal_name],
            interp = ['linear','linear','linear','linear','linear', 'linear'],
            append = False)

        # Cross hand calibration tables
        cross_table = [kcross, xftab]
        cross_field = [pacal_name, pacal_name]
        cross_interp = ['linear', 'linear']

        # Always check if the cross-hand phase is continuous per scan
        # Get the cross-hand phase            
        tb.open(xftab)
        gains = tb.getcol('CPARAM')
        flags = tb.getcol('FLAG')
        scans = tb.getcol('SCAN_NUMBER')
        tb.close()

        print(f"Gains shape: {gains.shape}")
        print(f"Unique scans: {np.unique(scans)}")

        # Process per scan
        unique_scans = np.unique(scans)
        manual_XF_per_scan = {}

        for scan in unique_scans:
            scan_mask = (scans == scan)
           
            # Get gains for this scan (median over time axis)
            scan_gains = np.nanmedian(gains[0, :, scan_mask].T, axis=-1)
            scan_flags = np.nanmedian(flags[0, :, scan_mask].T, axis=-1).astype(bool)
            scan_gains = scan_gains[~scan_flags]
            
            # Check continuity
            phases = np.angle(scan_gains)
            phase_diffs = np.diff(phases)
            phase_diffs = np.arctan2(np.sin(phase_diffs), np.cos(phase_diffs))
            max_jump = np.max(np.abs(np.degrees(phase_diffs)))
            
            is_continuous = max_jump < XF_AUTO_ANG_JUMP
            manual_XF_per_scan[scan] = not is_continuous
            
            print(f"Scan {scan}: max jump = {max_jump:.2f}°, continuous = {is_continuous}")

        # Determine which scans are continuous
        continuous_scans = [scan for scan in unique_scans if not manual_XF_per_scan[scan]]
        discontinuous_scans = [scan for scan in unique_scans if manual_XF_per_scan[scan]]

        print("\n" + "="*60)
        
        if len(continuous_scans) == len(unique_scans):
            # Case 1: All scans are continuous
            print(f"✓ All {len(unique_scans)} scan(s) are phase continuous")
            print("  Proceeding with standard XF calibration")
            
        elif len(continuous_scans) > 0:
            # Case 2: Some scans are continuous, some are not
            print(f"✓ {len(continuous_scans)} scan(s) are phase continuous: {continuous_scans}")
            print(f"✗ {len(discontinuous_scans)} scan(s) failed continuity check: {discontinuous_scans}")
            
            # Move bad tables
            bad_kcross = kcross.replace('.KCROSS', '_bad.KCROSS')
            bad_xftab = xftab.replace('.Xf', '_bad.Xf')
            print(f"  Moving original KCROSS table to: {bad_kcross}")
            shutil.move(kcross, bad_kcross)
            print(f"  Moving original Xf table to: {bad_xftab}")
            shutil.move(xftab, bad_xftab)
            
            # Remake KCROSS with only continuous scans
            print(f"  Remaking KCROSS table using only continuous scans: {continuous_scans}")
            gaincal(vis = myms,
                field = pacal_name,
                scan = ','.join(map(str, continuous_scans)),
                caltable = kcross,
                refant = str(ref_ant),
                solint = 'inf',
                gaintype='KCROSS',
                parang = True,
                gaintable=[ktab,gptab,bptab,gatab,dftab],
                gainfield=[pacal_name, pacal_name,bpcal_name, pacal_name, bpcal_name],
                interp = ['linear','linear','linear','linear','linear'],
                append = False)
            
            # Remake Xf table with only continuous scans
            print(f"  Remaking Xf table using only continuous scans: {continuous_scans}")
            polcal(vis = myms,
                field = pacal_name,
                scan = ','.join(map(str, continuous_scans)),
                uvrange = myuvrange,
                caltable = xftab,
                refant = str(ref_ant),
                solint = f'inf,{XF_CHANINT}ch',
                poltype='Xf',
                combine = '',
                gaintable=[ktab,gptab,bptab,gatab,dftab, kcross],
                gainfield=[pacal_name,pacal_name,bpcal_name, pacal_name, bpcal_name, pacal_name],
                interp = ['linear','linear','linear','linear','linear', 'linear'],
                append = False)
            
        else:
            # Case 3: No scans are continuous
            print(f"✗ None of the {len(unique_scans)} scan(s) are phase continuous")
            
            # Move bad tables
            bad_kcross = kcross.replace('.KCROSS', '_bad.KCROSS')
            bad_xftab = xftab.replace('.Xf', '_bad.Xf')
            print(f"  Moving original KCROSS table to: {bad_kcross}")
            shutil.move(kcross, bad_kcross)
            print(f"  Moving original Xf table to: {bad_xftab}")
            shutil.move(xftab, bad_xftab)
            
            if len(unique_scans) == 1:
                # Single scan - can use manual solver
                print("  Single scan detected - falling back to manual XF solver")
                manual_XF = True
            else:
                # Multiple scans - cannot handle
                print("\nERROR: Cannot currently handle multiple cross-hand phase scans in a single MS.")
                print("       Cross-hand phase is stable over week timescales, so unless this is a")
                print("       week-long MS file, one scan should be sufficient.")
                print(f"       Please split out an MS file with only one XF scan and re-run.")
                sys.exit(1)

        print("="*60)

    if XF_MODE == 'manual' or manual_XF or band == 'UHF':
        print(f"Cross-hand phase mode is {XF_MODE} (or auto detected a large phase jump? {manual_XF}); solving manual solutions")
        exec(open('tools/manual_XF_solver.py').read())

        # Cross hand calibration tables
        cross_table = [xftab]
        cross_field = [pacal_name]
        cross_interp = ['linear']     


# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #
# --------------------------- Applycal (All Fields)  ----------------------- #
# ------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------ #

# ------- BPCAL

applycal(vis = myms,
    gaintable = [ktab,gptab,bptab,ftab,dftab],
 #  sapplymode='calflagstrict',
    field = bpcal_name,
    #calwt = False,
    parang = True,
    gainfield = [bpcal_name,bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp = ['linear','linear','linear','linear','linear'],
    flagbackup=False)


# ----- If no polarization angle calibrator apply subset of tables and kill script

if pacal_name == '':   

    # ------- Secondaries 

    for i in range(0,len(pcal_names)):

        pcal = pcal_names[i]

        # ------- Check if pcal is the same as bpcal_name or pacal_name
        if pcal == bpcal_name or pcal == pacal_name:
            # If so, skip to the next iteration as it is already in the working tables
            continue
    
        applycal(vis = myms,
            gaintable = [ktab,gptab,bptab,ftab,dftab],
            # applymode='calflagstrict',
            field = pcal,
            #calwt = False,
            parang = False,
            gainfield = [pcal,pcal, bpcal_name, pcal, bpcal_name],
            interp = ['linear','linear','linear','linear','linear'],
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
                parang=False,
                gainfield = [related_pcal, related_pcal, bpcal_name, related_pcal, bpcal_name],
                interp = ['linear','linear','linear','linear','linear'],
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

# ------- PACAL

applycal(vis = myms,
 #       applymode='calflagstrict',
        field = pacal_name,
        #calwt = False,
        parang = True,
        gaintable = [ktab,gptab,bptab,ftab,dftab] + cross_table,
        gainfield = [pacal_name,pacal_name, bpcal_name, pacal_name, bpcal_name] + cross_field,
        interp = ['linear','linear','linear','linear','linear'] + cross_interp,
        flagbackup=False)

# ------- Secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    applycal(vis = myms,
        # applymode='calflagstrict',
        field = pcal,
        #calwt = False,
        parang = True,
        gaintable = [ktab,gptab,bptab,ftab,dftab] + cross_table,
        gainfield = [pcal,pcal, bpcal_name, pcal, bpcal_name] + cross_field,
        interp = ['linear','linear','linear','linear','linear'] + cross_interp,
        flagbackup=False)

# ------- Targets 
for i in range(0,len(targets)):

    target = targets[i]
    related_pcal = target_cal_map[i]

    applycal(vis=myms,
                #applymode='calflagstrict',
                field=target,
                #calwt=False,
                parang=True,
                gaintable = [ktab,gptab,bptab,ftab,dftab] + cross_table,
                gainfield = [related_pcal, related_pcal, bpcal_name, related_pcal, bpcal_name] + cross_field,
                interp = ['linear','linear','linear','linear','linear'] + cross_interp,
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

# ---- Apply aggressive flags if desired
if CAL_1GC_AGGRESSIVE_FLAGS and CAL_1GC_BL_FREQS != []:

    flagspw = ','.join(CAL_1GC_BL_FREQS)

    flagdata(vis = myms,
        mode = 'manual',
        spw = flagspw,
        flagbackup=False)

# ---- Save flags

flagmanager(vis=myms,
    mode='delete',
    versionname='1GC_flags')

flagmanager(vis=myms,
    mode='save',
    versionname='1GC_flags')
