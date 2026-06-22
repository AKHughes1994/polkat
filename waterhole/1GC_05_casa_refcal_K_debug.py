# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

# --------------------------------------------------------------- #
# DEBUG SCRIPT: K0/G0/B0 ordering, primary only.
#
# Five versions of K0/G0/B0 ordering, each ending right after its B0
# solve. Every table solved below gets a ragavi-gains plot (subprocess
# call, same pattern as oxkat/1GC_06_plot_gaintables.py) into GAINPLOTS,
# with an informative title.
# --------------------------------------------------------------- #

import shutil
import subprocess
import sys

# Flush immediately for better logging
import functools
print = functools.partial(print, flush=True)

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

gapfill = CAL_1GC_FILLGAPS
myspw = CAL_1GC_FREQRANGE

# Optional: Override uvrange for primary calibrator only (leave as '' to use myuvrange)
# primary_uvrange = '>1000m'
primary_uvrange = ''

# Use primary_uvrange if set, otherwise use myuvrange
if primary_uvrange != '':
    primary_uvrange_use = primary_uvrange
else:
    primary_uvrange_use = CAL_1GC_UVRANGE

# Restore the auto_cal flag version
flagmanager(vis=myms,
        mode='restore',
        versionname='autoflag_cals_data')

# ------- Set BP calibrator model

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
        # subprocess.run([syscall],shell=True)
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

# --------------------------------------------------------------- #
# --------------------- DEBUG: K/G/B ordering ------------------- #
# --------------------------------------------------------------- #

os.makedirs(GAINPLOTS, exist_ok=True)

def plot_caltable(caltab, title):
    # Same pattern as oxkat/1GC_06_plot_gaintables.py: gaintype is derived
    # from the caltable's own extension, not passed in. ragavi-gains has no
    # --title flag, so `title` is only printed here (shows up in the log).
    gaintype = caltab.split('.')[-1]
    gaintype = ''.join([i for i in gaintype if not i.isdigit()])
    if gaintype == 'Ga' or gaintype == 'Gp':
        gaintype = 'G'
    htmlname = GAINPLOTS+'/'+caltab.split('/')[-1]+'.html'
    plotname = GAINPLOTS+'/'+caltab.split('/')[-1]+'.png'
    if not os.path.isfile(htmlname):
        syscall = 'ragavi-gains -g '+gaintype+' -t '+caltab+' --htmlname='+htmlname+' --plotname='+plotname
        print('DEBUG PLOT ['+title+']: '+syscall)
        subprocess.run([syscall],shell=True)
    else:
        print(htmlname+' exists, skipping')

# ---- Version 1: K0 -> G0 (apply K0) -> B0 (apply K0, G0)

print('')
print('='*70)
print('[V1] K0 -> G0 -> B0: starting')
print('='*70)

ktab0_dbg1 = GAINTABLES+'/cal_1GC_'+myms+'_dbgKGB.K0'
gtab0_dbg1 = GAINTABLES+'/cal_1GC_'+myms+'_dbgKGB.G0'
bptab0_dbg1 = GAINTABLES+'/cal_1GC_'+myms+'_dbgKGB.B0'

print('[V1] Solving K0 (no prior) -> '+ktab0_dbg1)
gaincal(vis=myms, field=bpcal_name, caltable=ktab0_dbg1,
    uvrange=primary_uvrange_use, spw=myspw, refant=str(ref_ant),
    gaintype='K', solint='inf')

print('[V1] Solving G0 (apply K0) -> '+gtab0_dbg1)
gaincal(vis=myms, field=bpcal_name, caltable=gtab0_dbg1,
    uvrange=primary_uvrange_use, spw=myspw, refant=str(ref_ant),
    gaintype='G', solint='inf', calmode='p', minsnr=3,
    gaintable=[ktab0_dbg1], gainfield=[bpcal_name], interp=['nearest'])

print('[V1] Solving B0 (apply K0, G0) -> '+bptab0_dbg1)
bandpass(vis=myms, field=bpcal_name, caltable=bptab0_dbg1,
    uvrange=primary_uvrange_use, refant=str(ref_ant), solint='inf',
    combine='', solnorm=False, minblperant=4, minsnr=3.0, bandtype='B',
    fillgaps=gapfill,
    gaintable=[ktab0_dbg1, gtab0_dbg1], gainfield=[bpcal_name, bpcal_name],
    interp=['nearest', 'linear'])

print('[V1] Plotting K0, G0, B0')
plot_caltable(ktab0_dbg1, 'V1 K0->G0->B0: K0 (no prior)')
plot_caltable(gtab0_dbg1, 'V1 K0->G0->B0: G0 (apply K0)')
plot_caltable(bptab0_dbg1, 'V1 K0->G0->B0: B0 (apply K0, G0)')
print('[V1] Done')

# ---- Version 2: G0 -> K0 (apply G0) -> B0 (apply G0, K0)

print('')
print('='*70)
print('[V2] G0 -> K0 -> B0: starting')
print('='*70)

gtab0_dbg2 = GAINTABLES+'/cal_1GC_'+myms+'_dbgGKB.G0'
ktab0_dbg2 = GAINTABLES+'/cal_1GC_'+myms+'_dbgGKB.K0'
bptab0_dbg2 = GAINTABLES+'/cal_1GC_'+myms+'_dbgGKB.B0'

print('[V2] Solving G0 (no prior) -> '+gtab0_dbg2)
gaincal(vis=myms, field=bpcal_name, caltable=gtab0_dbg2,
    uvrange=primary_uvrange_use, spw=myspw, refant=str(ref_ant),
    gaintype='G', solint='inf', calmode='p', minsnr=3)

print('[V2] Solving K0 (apply G0) -> '+ktab0_dbg2)
gaincal(vis=myms, field=bpcal_name, caltable=ktab0_dbg2,
    uvrange=primary_uvrange_use, spw=myspw, refant=str(ref_ant),
    gaintype='K', solint='inf',
    gaintable=[gtab0_dbg2], gainfield=[bpcal_name], interp=['nearest'])

print('[V2] Solving B0 (apply G0, K0) -> '+bptab0_dbg2)
bandpass(vis=myms, field=bpcal_name, caltable=bptab0_dbg2,
    uvrange=primary_uvrange_use, refant=str(ref_ant), solint='inf',
    combine='', solnorm=False, minblperant=4, minsnr=3.0, bandtype='B',
    fillgaps=gapfill,
    gaintable=[gtab0_dbg2, ktab0_dbg2], gainfield=[bpcal_name, bpcal_name],
    interp=['nearest', 'nearest'])

print('[V2] Plotting G0, K0, B0')
plot_caltable(gtab0_dbg2, 'V2 G0->K0->B0: G0 (no prior)')
plot_caltable(ktab0_dbg2, 'V2 G0->K0->B0: K0 (apply G0)')
plot_caltable(bptab0_dbg2, 'V2 G0->K0->B0: B0 (apply G0, K0)')
print('[V2] Done')

# ---- Version 3: K0 -> B0 (apply K0 only)

print('')
print('='*70)
print('[V3] K0 -> B0: starting')
print('='*70)

ktab0_dbg3 = GAINTABLES+'/cal_1GC_'+myms+'_dbgKB.K0'
bptab0_dbg3 = GAINTABLES+'/cal_1GC_'+myms+'_dbgKB.B0'

print('[V3] Solving K0 (no prior) -> '+ktab0_dbg3)
gaincal(vis=myms, field=bpcal_name, caltable=ktab0_dbg3,
    uvrange=primary_uvrange_use, spw=myspw, refant=str(ref_ant),
    gaintype='K', solint='inf')

print('[V3] Solving B0 (apply K0 only) -> '+bptab0_dbg3)
bandpass(vis=myms, field=bpcal_name, caltable=bptab0_dbg3,
    uvrange=primary_uvrange_use, refant=str(ref_ant), solint='inf',
    combine='', solnorm=False, minblperant=4, minsnr=3.0, bandtype='B',
    fillgaps=gapfill,
    gaintable=[ktab0_dbg3], gainfield=[bpcal_name], interp=['nearest'])

print('[V3] Plotting K0, B0')
plot_caltable(ktab0_dbg3, 'V3 K0->B0: K0 (no prior)')
plot_caltable(bptab0_dbg3, 'V3 K0->B0: B0 (apply K0 only)')
print('[V3] Done')

# ---- Version 4: G0 -> B0 (apply G0 only)

print('')
print('='*70)
print('[V4] G0 -> B0: starting')
print('='*70)

gtab0_dbg4 = GAINTABLES+'/cal_1GC_'+myms+'_dbgGB.G0'
bptab0_dbg4 = GAINTABLES+'/cal_1GC_'+myms+'_dbgGB.B0'

print('[V4] Solving G0 (no prior) -> '+gtab0_dbg4)
gaincal(vis=myms, field=bpcal_name, caltable=gtab0_dbg4,
    uvrange=primary_uvrange_use, spw=myspw, refant=str(ref_ant),
    gaintype='G', solint='inf', calmode='p', minsnr=3)

print('[V4] Solving B0 (apply G0 only) -> '+bptab0_dbg4)
bandpass(vis=myms, field=bpcal_name, caltable=bptab0_dbg4,
    uvrange=primary_uvrange_use, refant=str(ref_ant), solint='inf',
    combine='', solnorm=False, minblperant=4, minsnr=3.0, bandtype='B',
    fillgaps=gapfill,
    gaintable=[gtab0_dbg4], gainfield=[bpcal_name], interp=['nearest'])

print('[V4] Plotting G0, B0')
plot_caltable(gtab0_dbg4, 'V4 G0->B0: G0 (no prior)')
plot_caltable(bptab0_dbg4, 'V4 G0->B0: B0 (apply G0 only)')
print('[V4] Done')

# ---- Version 5: B0 alone (no prior gain tables)

print('')
print('='*70)
print('[V5] B0 only: starting')
print('='*70)

bptab0_dbg5 = GAINTABLES+'/cal_1GC_'+myms+'_dbgBonly.B0'

print('[V5] Solving B0 (no prior calibration) -> '+bptab0_dbg5)
bandpass(vis=myms, field=bpcal_name, caltable=bptab0_dbg5,
    uvrange=primary_uvrange_use, refant=str(ref_ant), solint='inf',
    combine='', solnorm=False, minblperant=4, minsnr=3.0, bandtype='B',
    fillgaps=gapfill)

print('[V5] Plotting B0')
plot_caltable(bptab0_dbg5, 'V5 B0 only: no prior calibration')
print('[V5] Done')

print('')
print('='*70)
print('DEBUG_KGB_ORDERING complete -- 5 K/G/B versions solved and plotted in GAINPLOTS, primary only.')
print('='*70)
