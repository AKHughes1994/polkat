# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import glob
import shutil
import time
import datetime
import subprocess
import sys
import os


exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

# Apply PRE_FIELDS filtering if specified
if PRE_FIELDS != '':
    target_names = user_targets
    pcal_names = user_pcals

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
            time_info.append((opms, target, 'all', t_start_mjd, t_end_mjd))
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
            
            time_info.append((opms, target, 'all', t_start_mjd, t_end_mjd))

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
        time_info.append((opms, bpcal_name, scan_str, t_start_mjd, t_end_mjd))
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
        
        time_info.append((opms, bpcal_name, scan_str, t_start_mjd, t_end_mjd))


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
            time_info.append((opms, pcal, scan_str, t_start_mjd, t_end_mjd))
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
            
            time_info.append((opms, pcal, scan_str, t_start_mjd, t_end_mjd))


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
            time_info.append((opms, pacal_name, scan_str, t_start_mjd, t_end_mjd))
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
            
            time_info.append((opms, pacal_name, scan_str, t_start_mjd, t_end_mjd))


# Write timing information to file
output_file = os.path.join(RESULTS, myms.rstrip('/') + '_time_info.txt')
with open(output_file, 'w') as f:
    f.write('# MS timing information\n')
    f.write('# MS_NAME                                          FIELD_NAME            SCAN      START_MJD            END_MJD\n')
    for entry in time_info:
        ms_name, field_name, scan, t_start, t_end = entry
        f.write('%-50s %-21s %-9s %.10f %.10f\n' % (ms_name, field_name, scan, t_start, t_end))

print('Timing information saved to: ' + output_file)
