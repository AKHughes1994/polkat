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

DEBUG_PRINT_FLAGS = False

gapfill = CAL_1GC_FILLGAPS
myuvrange = CAL_1GC_UVRANGE 
myspw = CAL_1GC_FREQRANGE

# Optional: Override uvrange for primary calibrator only (leave as '' to use myuvrange)
# primary_uvrange = '>1000m'
primary_uvrange = ''

# Use primary_uvrange if set, otherwise use myuvrange
if primary_uvrange != '':
    primary_uvrange_use = primary_uvrange
else:
    primary_uvrange_use = myuvrange


# ------- Setup names


tt = stamp()

# Calibrator tables
ktab0 = GAINTABLES+'/cal_1GC_'+myms+'.K0'
bptab0 = GAINTABLES+'/cal_1GC_'+myms+'.B0'
gtab0 = GAINTABLES+'/cal_1GC_'+myms+'.G0'
dftab0  = GAINTABLES+'/cal_1GC_'+myms+'.Df0'

ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gtab = GAINTABLES+'/cal_1GC_'+myms+'.G'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'

kcross  = GAINTABLES+'/cal_1GC_'+myms+'.KCROSS'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'

# Restore the auto_cal flag version
flagmanager(vis=myms,
        mode='restore',
        versionname='autoflag_cals_data')

# Remove ftabs if they exist to prevent code from breaking
if os.path.isdir(ftab):
    print(f"Removing: {ftab}")
    shutil.rmtree(ftab)

# ------- Set BP calibrator models

# Check if MODEL_DATA column exists, if not initialize it
tb.open(myms)
if 'MODEL_DATA' not in tb.colnames():
    tb.close()
    # Dummy setjy call to initialize non-existing model data column
    setjy(vis=myms,
        standard='manual',
        field=bpcal_name,
        fluxdensity=[1.0, 0, 0, 0],
        reffreq='1000MHz',
        usescratch=True)
else:
    tb.close()

if primary_tag == '1934':
    
    # MeerKAT specific crystalball models for 1939 from B.Hugo: https://archive-gw-1.kat.ac.za/public/repository/10.48479/hhhy-4r55/index.htmlV
    if band == 'L':
        syscall = f"crystalball {myms} -f {bpcal_name} -sm {DATA}/crystalball/fitted.PKS1934.LBand.wsclean.cat.txt"
        subprocess.run([syscall],shell=True)

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
    caltable=ktab0,
    uvrange=primary_uvrange_use,
    spw=myspw,
    refant=str(ref_ant),
    gaintype='K',
    solint='inf')


# ------- B0 (primary; apply K0)

bandpass(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=bptab0,
    refant=str(ref_ant),
    solint='inf',
    combine='',
    solnorm=False,
    minblperant=4,
    minsnr=3.0,
    bandtype='B',
    fillgaps=gapfill,
    gainfield=[bpcal_name],
    interp=['nearest'],
    gaintable=[ktab0])

flagdata(vis=bptab0, mode='tfcrop', datacolumn='CPARAM', flagbackup=False)
flagdata(vis=bptab0, mode='rflag', datacolumn='CPARAM', flagbackup=False)
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

# ------- G0 (primary; apply K0, B0) -- Type T, amp+phase

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    spw=myspw,
    caltable=gtab0,
    refant=str(ref_ant),
    gaintype='T',
    solint='inf',
    calmode='ap',
    minsnr=3,
    gainfield=[bpcal_name, bpcal_name],
    interp=['nearest', 'linear'],
    gaintable=[ktab0, bptab0])


# -------- Solve for Df0 (apply K0, B0, G0)

polcal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=dftab0,
    refant=str(ref_ant),
    solint='inf',
    poltype='Df',
    combine='scan',
    gaintable=[ktab0, bptab0, gtab0],
    gainfield=[bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear'],
    append=False)

flagdata(vis=dftab0, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

# ------- Correct primary data with K0, B0, G0, Df0

applycal(vis=myms,
    gaintable=[ktab0, bptab0, gtab0, dftab0],
    field=bpcal_name,
    parang=False,
    gainfield=[bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear', 'linear'],
    flagbackup=False)
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')


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
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

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


# ------- K (primary; no prior calibration)

gaincal(vis=myms,
    field=bpcal_name,
    caltable=ktab,
    uvrange=primary_uvrange_use,
    spw=myspw,
    refant=str(ref_ant),
    gaintype='K',
    solint='inf')


# ------- B (primary; apply K)

bandpass(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=bptab,
    refant=str(ref_ant),
    solint='inf',
    combine='',
    solnorm=False,
    minblperant=4,
    minsnr=3.0,
    bandtype='B',
    fillgaps=gapfill,
    gaintable=[ktab],
    gainfield=[bpcal_name],
    interp=['nearest'])

flagdata(vis=bptab, mode='tfcrop', datacolumn='CPARAM', flagbackup=False)
flagdata(vis=bptab, mode='rflag', datacolumn='CPARAM', flagbackup=False)
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')

# -------- G (primary; apply K, B) -- Type T, amp+phase

gaincal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    spw=myspw,
    caltable=gtab,
    refant=str(ref_ant),
    gaintype='T',
    solint='inf',
    calmode='ap',
    minsnr=3,
    gaintable=[ktab, bptab],
    gainfield=[bpcal_name, bpcal_name],
    interp=['nearest', 'linear'])


# -------- Solve for Df (primary; apply K, B, G)

polcal(vis=myms,
    field=bpcal_name,
    uvrange=primary_uvrange_use,
    caltable=dftab,
    refant=str(ref_ant),
    solint='inf',
    poltype='Df',
    combine='scan',
    gaintable=[ktab, bptab, gtab],
    gainfield=[bpcal_name, bpcal_name, bpcal_name],
    interp=['nearest', 'linear', 'linear'])

flagdata(vis=dftab, mode='clip', clipminmax=[0.0,0.1], flagbackup=False, datacolumn='CPARAM')
# DEBUGGING: summarize flags
if DEBUG_PRINT_FLAGS:
    print('DEBUG: PRINTING FLAGS')
    flagdata(myms, mode='summary')


# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- Initial Table (Secondary + Pol. Cal.)  ---------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #

if pacal_name != '':

    # ------- G0 (polcal; apply B, Df, K from primary) -- Type T, amp+phase

    gaincal(vis=myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab0,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='ap',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['linear', 'linear', 'linear'],
        append=True)

    # ------- K0 (polcal; apply B, Df from primary, G0 from polcal)

    gaincal(vis=myms,
        field=pacal_name,
        caltable=ktab0,
        uvrange=myuvrange,
        spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gtab0],
        gainfield=[bpcal_name, bpcal_name, pacal_name],
        interp=['linear', 'linear', 'linear'],
        append=True)

# ----- Loop over secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    # ------- G0 (pcal; apply B, Df, K from primary) -- Type T, amp+phase

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab0,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='ap',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear'],
        append=True)

    # ------- K0 (pcal; apply B, Df from primary, G0 from pcal)

    gaincal(vis=myms,
        field=pcal,
        caltable=ktab0,
        uvrange=myuvrange,
        spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gtab0],
        gainfield=[bpcal_name, bpcal_name, pcal],
        interp=['linear', 'linear', 'linear'],
        append=True)

if pacal_name != '':
    # -------- Applycal (polcal; B, Df from primary, K0, G0 from polcal) and Flag
    
    applycal(vis=myms,
        gaintable=[ktab0, bptab, gtab0, dftab],
        field=pacal_name,
        parang=False,
        gainfield=[pacal_name, bpcal_name, pacal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear', 'linear'],
        flagbackup=False)
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

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
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')


# ----- Loop over secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    # -------- Applycal (pcal; B, Df from primary, K0, G0 from pcal) and Flag

    applycal(vis=myms,
        gaintable=[ktab0, bptab, gtab0, dftab],
        field=pcal,
        parang=False,
        gainfield=[pcal, bpcal_name, pcal, bpcal_name],
        interp=['nearest', 'linear', 'linear', 'linear'],
        flagbackup=False)
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

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
    # DEBUGGING: summarize flags
    if DEBUG_PRINT_FLAGS:
        print('DEBUG: PRINTING FLAGS')
        flagdata(myms, mode='summary')

# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #
# --------------------------- Working Table (Secondary + Pol. Cal.)  ------------------------ #
# -------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------- #

if pacal_name != '':   
 
    # ------- G (polcal; apply B, Df, K from primary) -- Type T, amp+phase

    gaincal(vis=myms,
        field=pacal_name,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='ap',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear'],
        append=True)

    # ------- K (polcal; apply B, Df from primary, G from polcal)
    
    gaincal(vis=myms,
        field=pacal_name,
        caltable=ktab,
        uvrange=myuvrange,
        spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gtab],
        gainfield=[bpcal_name, bpcal_name, pacal_name],
        interp=['linear', 'linear', 'linear'],
        append=True)

            
for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    # ------- G (pcal; apply B, Df, K from primary) -- Type T, amp+phase

    gaincal(vis=myms,
        field=pcal,
        uvrange=myuvrange,
        spw=myspw,
        caltable=gtab,
        refant=str(ref_ant),
        gaintype='T',
        solint='inf',
        calmode='ap',
        minsnr=3,
        gaintable=[ktab, bptab, dftab],
        gainfield=[bpcal_name, bpcal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear'],
        append=True)

    # ------- K (pcal; apply B, Df from primary, G from pcal)

    gaincal(vis=myms,
        field=pcal,
        caltable=ktab,
        uvrange=myuvrange,
        spw=myspw,
        refant=str(ref_ant),
        gaintype='K',
        solint='inf',
        gaintable=[bptab, dftab, gtab],
        gainfield=[bpcal_name, bpcal_name, pcal],
        interp=['linear', 'linear', 'linear'],
        append=True)

# --- Apply fluxscaling to G but only if there are calibration fields other than the primary
if len([pcal for pcal in pcal_names if pcal != bpcal_name]) > 0 or pacal_name != '':
    
    # --- Apply fluxscaling to G
    fluxscale(vis=myms,
        caltable=gtab,
        fluxtable=ftab,
        reference=bpcal_name,
        append=False,
        transfer='')

# If there is no need to apply flux scaling, we can set ftab to gtab as the only calibrator is the primary    
else:
    ftab = gtab

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
    
    # Force manual mode for UHF band
    if band == 'UHF':
        print("UHF band detected, forcing manual XF mode")
        manual_XF = True
        XF_MODE = 'manual'
    
    if XF_MODE == 'casa' or XF_MODE == 'auto':

        print(f"Cross-hand phase mode is {XF_MODE}; solving CASA based solutions")
        
        # Setup combine parameter
        combine_param = 'scan' if XF_AVG_SCAN else ''
        
        # ------- KCROSS (polcal; apply B, Df from primary, G, K from polcal)
        # Always solve KCROSS (cheap computation), skip flag only controls application
        print(f"  Solving for KCROSS (cross-hand delay) {'[will NOT be applied]' if XF_SKIP_KCROSS else '[will be applied]'}")
        gaincal(vis = myms,
            field = pacal_name,
            caltable = kcross,
            refant = str(ref_ant),
            solint = 'inf',
            combine = combine_param,
            gaintype='KCROSS',
            parang = True,
            gaintable=[ktab, bptab, gtab, dftab],
            gainfield=[pacal_name, bpcal_name, pacal_name, bpcal_name],
            interp = ['linear','linear','linear','linear'],
            append = False)
        
        # If skipping KCROSS, immediately move it to bad name and never use it
        if XF_SKIP_KCROSS:
            bad_kcross = kcross.replace('.KCROSS', '_bad.KCROSS')
            print(f"  Moving KCROSS to: {bad_kcross} (will not be used)")
            if os.path.isdir(bad_kcross):
                shutil.rmtree(bad_kcross)
            shutil.move(kcross, bad_kcross)

        # -------- Xf (polcal; apply Bp, Df (primary), Ga, Gp, K, [KCROSS] (polcal))
        # Build gaintable list based on skip flagdani katseye twerking
        if XF_SKIP_KCROSS:
            print("  Solving for Xf (cross-hand phase) WITHOUT KCROSS")
        else:
            print("  Solving for Xf (cross-hand phase) with KCROSS")
        xf_gaintable = [ktab, bptab, gtab, dftab] + ([] if XF_SKIP_KCROSS else [kcross])
        xf_gainfield = [pacal_name, bpcal_name, pacal_name, bpcal_name] + ([] if XF_SKIP_KCROSS else [pacal_name])
        xf_interp = ['linear','linear','linear','linear'] + ([] if XF_SKIP_KCROSS else ['linear'])
        
        polcal(vis = myms,
            field = pacal_name,
            uvrange = myuvrange,
            caltable = xftab,
            refant = str(ref_ant),
            solint = f'inf,{XF_CHANINT}ch',
            poltype='Xf',
            combine = combine_param,
            gaintable=xf_gaintable,
            gainfield=xf_gainfield,
            interp=xf_interp,
            append = False)

        # Cross hand calibration tables for applycal
        cross_table = [xftab] if XF_SKIP_KCROSS else [kcross, xftab]
        cross_field = [pacal_name] if XF_SKIP_KCROSS else [pacal_name, pacal_name]
        cross_interp = ['linear'] if XF_SKIP_KCROSS else ['linear', 'linear']

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
            
            # Move bad Xf table
            bad_xftab = xftab.replace('.Xf', '_bad.Xf')
            print(f"  Moving original Xf table to: {bad_xftab}")
            shutil.move(xftab, bad_xftab)
            
            # Only handle KCROSS if it's being used (not already moved to _bad)
            if not XF_SKIP_KCROSS:
                # Move bad KCROSS table
                bad_kcross = kcross.replace('.KCROSS', '_bad.KCROSS')
                print(f"  Moving original KCROSS table to: {bad_kcross}")
                if os.path.isdir(bad_kcross):
                    shutil.rmtree(bad_kcross)
                shutil.move(kcross, bad_kcross)
                
                # Remake KCROSS with only continuous scans
                print(f"  Remaking KCROSS table using only continuous scans: {continuous_scans}")
                gaincal(vis = myms,
                    field = pacal_name,
                    scan = ','.join(map(str, continuous_scans)),
                    caltable = kcross,
                    refant = str(ref_ant),
                    solint = 'inf',
                    combine = combine_param,
                    gaintype='KCROSS',
                    parang = True,
                    gaintable=[ktab, bptab, gtab, dftab],
                    gainfield=[pacal_name, bpcal_name, pacal_name, bpcal_name],
                    interp = ['linear','linear','linear','linear'],
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
                combine = combine_param,
                gaintable=xf_gaintable,
                gainfield=xf_gainfield,
                interp=xf_interp,
                append = False)
            
        else:
            # Case 3: No scans are continuous
            print(f"✗ None of the {len(unique_scans)} scan(s) are phase continuous")
            
            # Move bad Xf table
            bad_xftab = xftab.replace('.Xf', '_bad.Xf')
            print(f"  Moving original Xf table to: {bad_xftab}")
            shutil.move(xftab, bad_xftab)
            
            # Only move KCROSS if it's being used (not already moved to _bad)
            if not XF_SKIP_KCROSS:
                bad_kcross = kcross.replace('.KCROSS', '_bad.KCROSS')
                print(f"  Moving original KCROSS table to: {bad_kcross}")
                if os.path.isdir(bad_kcross):
                    shutil.rmtree(bad_kcross)
                shutil.move(kcross, bad_kcross)
            
            # Fall back to manual solver (can handle both single and multiple scans)
            print(f"  Falling back to manual XF solver ({len(unique_scans)} scan(s))")
            manual_XF = True

        print("="*60)

    if XF_MODE == 'manual' or manual_XF:
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
    gaintable = [ktab, bptab, ftab, dftab],
    field = bpcal_name,
    parang = False,
    gainfield = [bpcal_name, bpcal_name, bpcal_name, bpcal_name],
    interp = ['nearest','linear','linear','linear'],
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
            gaintable = [ktab, bptab, ftab, dftab],
            field = pcal,
            parang = False,
            gainfield = [pcal, bpcal_name, pcal, bpcal_name],
            interp = ['nearest','linear','linear','linear'],
            flagbackup=False)

    # ------- Targets 
    for i in range(0,len(targets)):

        target = targets[i]
        related_pcal = target_cal_map[i]

        applycal(vis=myms,
                gaintable = [ktab, bptab, ftab, dftab],
                field=target,
                parang=False,
                gainfield = [related_pcal, bpcal_name, related_pcal, bpcal_name],
                interp = ['nearest','linear','linear','linear'],
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

# -------- Full polarization  # CAL_1GC_APPLYPARANG controls parang from here; above is XF-less path (no sky-frame rotation needed)

# ------- PACAL

applycal(vis = myms,
        field = pacal_name,
        parang = CAL_1GC_APPLYPARANG,
        gaintable = [ktab, bptab, ftab, dftab] + cross_table,
        gainfield = [pacal_name, bpcal_name, pacal_name, bpcal_name] + cross_field,
        interp = ['nearest','linear','linear','linear'] + cross_interp,
        flagbackup=False)

# ------- Secondaries

for i in range(0,len(pcal_names)):

    pcal = pcal_names[i]

    # ------- Check if pcal is the same as bpcal_name or pacal_name
    if pcal == bpcal_name or pcal == pacal_name:
        # If so, skip to the next iteration as it is already in the working tables
        continue

    applycal(vis = myms,
        field = pcal,
        parang = CAL_1GC_APPLYPARANG,
        gaintable = [ktab, bptab, ftab, dftab] + cross_table,
        gainfield = [pcal, bpcal_name, pcal, bpcal_name] + cross_field,
        interp = ['nearest','linear','linear','linear'] + cross_interp,
        flagbackup=False)

# ------- Targets 
for i in range(0,len(targets)):

    target = targets[i]
    related_pcal = target_cal_map[i]

    applycal(vis=myms,
                field=target,
                parang=CAL_1GC_APPLYPARANG,
                gaintable = [ktab, bptab, ftab, dftab] + cross_table,
                gainfield = [related_pcal, bpcal_name, related_pcal, bpcal_name] + cross_field,
                interp = ['nearest','linear','linear','linear'] + cross_interp,
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
