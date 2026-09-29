# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import subprocess
import sys
import os
import numpy as np


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

# Apply PRE_FIELDS filtering if specified
if PRE_FIELDS != '':
    target_names = user_targets
    pcal_names = user_pcals

def parang_range_deg(opms):
    """
    Parallactic angle range (min, max, delta, degrees) across every row of
    `opms`, computed directly from the MS: the field's phase centre and the
    array's own reference position (its first antenna, ITRF), via CASA's
    measures tool -- the same posangle-between-field-and-zenith method
    tools/correct_parang.py uses per-antenna, and the same convention
    waterhole/manual_XF_solver_continuity.compute_parallactic_angle uses for
    a single mid-epoch value, extended here to every dump time and unwrapped
    across the +-180 deg branch so a swing through it doesn't corrupt the
    range.
    """
    tb.open(opms + '/FIELD')
    ra_rad = tb.getcol('PHASE_DIR')[0][0][0]
    dec_rad = tb.getcol('PHASE_DIR')[1][0][0]
    tb.close()

    tb.open(opms + '/ANTENNA')
    ant0_pos = tb.getcol('POSITION')[:, 0]
    tb.close()

    tb.open(opms)
    times = np.unique(tb.getcol('TIME'))
    tb.close()

    me.doframe(me.position('itrf', *[qa.quantity(x, 'm') for x in ant0_pos]))
    field_dir = me.direction('J2000', qa.quantity(ra_rad, 'rad'), qa.quantity(dec_rad, 'rad'))
    zenith = me.direction('AZELGEO', '0deg', '90deg')

    chi_deg = np.empty(len(times))
    for i, t in enumerate(times):
        me.doframe(me.epoch('utc', qa.quantity(t, 's')))
        chi_deg[i] = np.degrees(me.posangle(field_dir, zenith).get_value('rad'))

    chi_unwrapped = np.degrees(np.unwrap(np.radians(chi_deg)))
    pa_min = float(np.min(chi_unwrapped))
    pa_max = float(np.max(chi_unwrapped))
    return pa_min, pa_max, pa_max - pa_min

# Initialize list to store timing information
time_info = []

# Ensure RESULTS directory exists
if not os.path.exists(RESULTS):
    os.makedirs(RESULTS)

if CAL_SKIP_CALS:
    print('CAL_SKIP_CALS is True -- skipping secondary calibrator splitting (primary and polarization-angle calibrator still split)')

# Split target fields (integrated over all scans)
for target in target_names:

    opms = ''

    for mm in target_ms:
        if target in mm:
            opms = mm

    if opms != '':

        if os.path.isdir(opms):
            print('Output MS already exists for '+target+', skipping mstransform: '+opms)
            # Still extract timing information
            tb.open(opms)
            times = tb.getcol('TIME')
            tb.close()
            t_start_mjd = times.min() / 86400.0
            t_end_mjd = times.max() / 86400.0
            pa_min, pa_max, pa_delta = parang_range_deg(opms)
            time_info.append((opms, target, 'all', t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))
        else:
            mstransform(vis=myms,
                outputvis=opms,
                field=target,
                usewtspectrum=True,
                realmodelcol=True,
                datacolumn='corrected')

            flagmanager(vis=opms,
                mode='save',
                versionname='post-1GC')

            # Extract timing information
            tb.open(opms)
            times = tb.getcol('TIME')
            tb.close()
            t_start_mjd = times.min() / 86400.0
            t_end_mjd = times.max() / 86400.0
            
            pa_min, pa_max, pa_delta = parang_range_deg(opms)
            time_info.append((opms, target, 'all', t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))

    else:
        print('Target/MS mismatch in project info for '+target+', please check.')


# Split primary calibrator fields (per scan)
for opms in primary_ms:
    
    # Extract scan number from MS filename (format: _scanXX.ms)
    scan_str = opms.split('_scan')[-1].replace('.ms', '')
    scan_num = str(int(scan_str))  # Remove zero-padding for CASA field selection
    
    if os.path.isdir(opms):
        print('Output MS already exists for '+bpcal_name+' scan '+scan_str+', skipping mstransform: '+opms)
        # Still extract timing information
        tb.open(opms)
        times = tb.getcol('TIME')
        tb.close()
        t_start_mjd = times.min() / 86400.0
        t_end_mjd = times.max() / 86400.0
        pa_min, pa_max, pa_delta = parang_range_deg(opms)
        time_info.append((opms, bpcal_name, scan_str, t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))
    else:
        mstransform(vis=myms,
            outputvis=opms,
            field=bpcal_name,
            scan=scan_num,
            usewtspectrum=True,
            realmodelcol=True,
            datacolumn='corrected')

        flagmanager(vis=opms,
            mode='save',
            versionname='post-1GC')

        # Extract timing information
        tb.open(opms)
        times = tb.getcol('TIME')
        tb.close()
        t_start_mjd = times.min() / 86400.0
        t_end_mjd = times.max() / 86400.0
        
        pa_min, pa_max, pa_delta = parang_range_deg(opms)
        time_info.append((opms, bpcal_name, scan_str, t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))


# Split secondary calibrator fields (per scan)
for i, pcal in enumerate(pcal_names):

    if CAL_SKIP_CALS:
        continue

    # pcal_ms[i] is a list of MS files for each scan of this secondary
    for opms in pcal_ms[i]:

        # Extract scan number from MS filename (format: _scanXX.ms)
        scan_str = opms.split('_scan')[-1].replace('.ms', '')
        scan_num = str(int(scan_str))  # Remove zero-padding for CASA field selection
        
        if os.path.isdir(opms):
            print('Output MS already exists for '+pcal+' scan '+scan_str+', skipping mstransform: '+opms)
            # Still extract timing information
            tb.open(opms)
            times = tb.getcol('TIME')
            tb.close()
            t_start_mjd = times.min() / 86400.0
            t_end_mjd = times.max() / 86400.0
            pa_min, pa_max, pa_delta = parang_range_deg(opms)
            time_info.append((opms, pcal, scan_str, t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))
        else:
            mstransform(vis=myms,
                outputvis=opms,
                field=pcal,
                scan=scan_num,
                usewtspectrum=True,
                realmodelcol=True,
                datacolumn='corrected')

            flagmanager(vis=opms,
                mode='save',
                versionname='post-1GC')

            # Extract timing information
            tb.open(opms)
            times = tb.getcol('TIME')
            tb.close()
            t_start_mjd = times.min() / 86400.0
            t_end_mjd = times.max() / 86400.0
            
            pa_min, pa_max, pa_delta = parang_range_deg(opms)
            time_info.append((opms, pcal, scan_str, t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))


# Split polarization angle calibrator fields (per scan)
if pacal_name != '':
    for opms in polang_ms:
        
        # Extract scan number from MS filename (format: _scanXX.ms)
        scan_str = opms.split('_scan')[-1].replace('.ms', '')
        scan_num = str(int(scan_str))  # Remove zero-padding for CASA field selection
        
        if os.path.isdir(opms):
            print('Output MS already exists for '+pacal_name+' scan '+scan_str+', skipping mstransform: '+opms)
            # Still extract timing information
            tb.open(opms)
            times = tb.getcol('TIME')
            tb.close()
            t_start_mjd = times.min() / 86400.0
            t_end_mjd = times.max() / 86400.0
            pa_min, pa_max, pa_delta = parang_range_deg(opms)
            time_info.append((opms, pacal_name, scan_str, t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))
        else:
            mstransform(vis=myms,
                outputvis=opms,
                field=pacal_name,
                scan=scan_num,
                usewtspectrum=True,
                realmodelcol=True,
                datacolumn='corrected')

            flagmanager(vis=opms,
                mode='save',
                versionname='post-1GC')

            # Extract timing information
            tb.open(opms)
            times = tb.getcol('TIME')
            tb.close()
            t_start_mjd = times.min() / 86400.0
            t_end_mjd = times.max() / 86400.0
            
            pa_min, pa_max, pa_delta = parang_range_deg(opms)
            time_info.append((opms, pacal_name, scan_str, t_start_mjd, t_end_mjd, pa_min, pa_max, pa_delta))


# Write timing information to file
output_file = os.path.join(RESULTS, myms.rstrip('/') + '_time_info.txt')
with open(output_file, 'w') as f:
    f.write('# MS timing information\n')
    f.write('# MS_NAME                                          FIELD_NAME            SCAN      START_MJD            END_MJD              PA_MIN_DEG    PA_MAX_DEG    PA_DELTA_DEG\n')
    for entry in time_info:
        ms_name, field_name, scan, t_start, t_end, pa_min, pa_max, pa_delta = entry
        f.write('%-50s %-21s %-9s %.10f %.10f %-13.4f %-13.4f %-13.4f\n' %
               (ms_name, field_name, scan, t_start, t_end, pa_min, pa_max, pa_delta))

print('Timing information saved to: ' + output_file)
