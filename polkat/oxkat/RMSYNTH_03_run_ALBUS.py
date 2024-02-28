#!/usr/bin/python
# ian.heywood@physics.ox.ac.uk

import logging, json, sys, subprocess
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from pyrap.tables import table
import os.path as o

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

def casa_to_albus_pos(pos0):
    ra = pos0[0] * 180.0 / np.pi
    dec = pos0[1] * 180.0 / np.pi

    coord = SkyCoord(ra = ra * u.deg, dec = dec * u.deg)
    return coord.to_string('hmsdms').replace('d',':').replace('h',':').replace('m',':').replace('s','').replace('+','plus').replace('-','minus').split(' ')

def main():
    
    # Load in the the prefields and project infor .json dictionaries
    with open('project_info.json') as f:
        project_info = json.load(f)

    # Specify myms
    myms  = project_info['working_ms']

    # Load in the  FIELD identifier Information
    fldtab = table(myms+'/FIELD', ack=False)
    names  = fldtab.getcol('NAME')
    ids         = fldtab.getcol('SOURCE_ID')
    dirs        = fldtab.getcol('PHASE_DIR')
    time        = fldtab.getcol('TIME')
    fldtab.done()

    # Parse the scan list
    maintab = table(myms, ack=False)
    scan_info = {}

    # Initialize array to store START and END times for the observations of each source
    for field_id in ids:
        scan_info[f'{field_id}'] = []
        
    scans = np.unique(maintab.getcol('SCAN_NUMBER'))
    for sc in scans:
            subtab = maintab.query(query='SCAN_NUMBER=='+str(sc))
            st0 = subtab.getcol('TIME')[0]
            st1 = subtab.getcol('TIME')[-1]
            sfield = np.unique(subtab.getcol('FIELD_ID'))[0]
            scan_info[f'{sfield}'].extend([st0, st1]) # Group by Field IDs
    
    subtab = maintab.query(query='SCAN_NUMBER==1')
    maintab.done()


    # Get the START and FINISH times of each source
    start_time   = []
    finish_time = []
    ra    = []
    dec = []
    for field_id in ids:
        start_time.append(np.amin(scan_info[f'{field_id}']) /86400.0)
        finish_time.append(np.amax(scan_info[f'{field_id}']) /86400.0)

    # Format the positions/times according to the Requirements from albus
    start_time   = Time(start_time, format='mjd', scale='utc').iso
    finish_time = Time(finish_time, format='mjd', scale='utc').iso
    start_times   = [start.replace('-','/') for start in start_time]
    finish_times   = [finish.replace('-','/') for finish in finish_time]

    for src_dir in dirs:
        ra.append(casa_to_albus_pos(src_dir[0])[0])
        dec.append(casa_to_albus_pos(src_dir[0])[1])

    # Iterate through sources and run ALBUS on the polarization calibrator or any target source
    for k in range(len(ra)):
        if names[k] in project_info['target_names'] or names[k] == project_info['polang_name']:
            syscall = f'python3 {cfg.TOOLS}/ALBUS_get_ionosphere.py -r "{ra[k]}" -d "{dec[k]}" -s "{start_times[k]}" -f "{finish_times[k]}" "{names[k]}"'
            subprocess.run([syscall], shell=True)


if __name__ == "__main__":
    main()

