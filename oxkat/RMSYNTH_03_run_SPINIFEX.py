#!/usr/bin/python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import logging, json, sys
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astropy.time import Time
from pyrap.tables import table
import os.path as o
import shutil
from functools import partial
from spinifex import get_rm

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

# Reinitialize print to flush in real-time for SLURM logging
print = partial(print, flush=True)

def casa_to_skycoord(pos0):
    """Convert CASA phase direction to astropy SkyCoord"""
    ra = pos0[0] * 180.0 / np.pi
    dec = pos0[1] * 180.0 / np.pi
    return SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

def skycoord_to_hmsdms(coord):
    """Convert SkyCoord to hmsdms string format"""
    return coord.to_string('hmsdms').replace('d', ':').replace('h', ':').replace('m', ':').replace('s', '')

def main():
    
    # Load project info
    with open('project_info.json') as f:
        project_info = json.load(f)

    myms = project_info['working_ms']
    
    # Define MeerKAT location
    lon = 21.44326 * u.deg
    lat = -30.712455 * u.deg
    meerkat = EarthLocation(lon=lon, lat=lat, height=1059.662443 * u.m)
    
    # Load FIELD information
    fldtab = table(myms+'/FIELD', ack=False)
    names = fldtab.getcol('NAME')
    ids = fldtab.getcol('SOURCE_ID')
    dirs = fldtab.getcol('PHASE_DIR')
    fldtab.done()
    
    # Parse scan information
    maintab = table(myms, ack=False)
    scan_info = {}  # Will store {field_id: {scan_num: [start_time, end_time]}}
    
    for field_id in ids:
        scan_info[field_id] = {}
    
    scans = np.unique(maintab.getcol('SCAN_NUMBER'))
    
    for sc in scans:
        subtab = maintab.query(query='SCAN_NUMBER=='+str(sc))
        st0 = subtab.getcol('TIME')[0]
        st1 = subtab.getcol('TIME')[-1]
        sfield = np.unique(subtab.getcol('FIELD_ID'))[0]
        
        if sc not in scan_info[sfield]:
            scan_info[sfield][sc] = []
        scan_info[sfield][sc] = [st0, st1]
        subtab.done()
    
    maintab.done()
    
    # Open output file in results directory
    output_file = o.join(cfg.RESULTS, f"{o.basename(myms)}_spinifex_rm.txt")
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Field_Name Position_hmsdms Scan_Number Time_ISOT Time_MJD RM_rad/m2 RM_err_rad/m2\n")
        
        # Process each field
        for field_id, field_name, phase_dir in zip(ids, names, dirs):
            print(f"Processing field: {field_name} (ID: {field_id})")
            
            # Convert position to SkyCoord
            source = casa_to_skycoord(phase_dir[0])
            pos_string = skycoord_to_hmsdms(source)
            
            # Process each scan for this field
            if field_id in scan_info:
                for scan_num in sorted(scan_info[field_id].keys()):
                    print(f"  Scan {scan_num}")
                    
                    # Get scan start and end times (in MJD seconds)
                    t_start_mjd_sec = scan_info[field_id][scan_num][0]
                    t_end_mjd_sec = scan_info[field_id][scan_num][1]
                    
                    # Convert to MJD days
                    t_start = Time(t_start_mjd_sec / 86400.0, format='mjd', scale='utc')
                    t_end = Time(t_end_mjd_sec / 86400.0, format='mjd', scale='utc')
                    
                    # Create 1-minute interval time array
                    duration = (t_end - t_start).to(u.min)
                    n_intervals = int(np.ceil(duration.value)) + 1
                    times = t_start + np.arange(n_intervals) * 1 * u.min
                    
                    # Run spinifex
                    try:
                        rm = get_rm.get_rm_from_skycoord(
                            loc=meerkat,
                            times=times,
                            output_directory=cfg.RESULTS + '/meerkat_gps_data',
                            source=source,
                            iono_model_name='ionex_iri'

                        )
                        
                        # Write results
                        for myrm, rm_err, tm in zip(rm.rm, rm.rm_error, rm.times):
                            f.write(f"{field_name} {pos_string} {scan_num} {tm.isot} {tm.mjd} {myrm} {rm_err}\n")
                    
                    except Exception as e:
                        print(f"    Warning: Failed to get RM for scan {scan_num}: {e}")
                        continue
    
    print(f"\nResults written to: {output_file}")
    
    # Remove intermediate MeerKAT GPS data
    gps_data_dir = o.join(cfg.RESULTS, 'meerkat_gps_data')
    if o.exists(gps_data_dir):
        shutil.rmtree(gps_data_dir)
        print(f"Removed intermediate GPS data: {gps_data_dir}")

if __name__ == "__main__":
    main()

