#!/usr/bin/python
# andrew.hughes@physics.ox.ac.uk

import sys
import numpy
import logging
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
from pyrap.tables import table
from optparse import OptionParser


def main():

    parser = OptionParser(usage='%prog [options] msname')
    (options, args) = parser.parse_args()

    if len(args) != 1:
        print('Please specify a Measurement Set')
        sys.exit()
    else:
        myms = args[0].rstrip('/')

    # Extract MS name for log file
    msname = myms.split('/')[-1]
    logfile = 'elevations_' + msname + '.log'

    # Setup logging
    logging.basicConfig(filename=logfile, level=logging.DEBUG, 
                       format='%(asctime)s |  %(message)s', 
                       datefmt='%d/%m/%Y %H:%M:%S ')
    stream = logging.StreamHandler()
    stream.setLevel(logging.DEBUG)
    streamformat = logging.Formatter('%(asctime)s |  %(message)s', 
                                    datefmt='%d/%m/%Y %H:%M:%S ')
    stream.setFormatter(streamformat)
    mylogger = logging.getLogger(__name__)
    mylogger.setLevel(logging.DEBUG)
    mylogger.addHandler(stream)

    mylogger.info('')
    mylogger.info('--MS: ' + myms)
    mylogger.info('')

    # MEERKAT telescope location
    meerkat_lon = "21:26:35.736"
    meerkat_lat = "-30:42:44.838"
    meerkat_height = 1059.662443  # metres

    # Convert MeerKAT location to EarthLocation
    meerkat = EarthLocation.from_geodetic(
        lon=meerkat_lon,
        lat=meerkat_lat,
        height=meerkat_height * u.m
    )

    # ------- GETTING INFORMATION FROM MS -------

    # Get field information
    fldtab = table(myms + '/FIELD', ack=False)
    names = fldtab.getcol('NAME')
    ids = fldtab.getcol('SOURCE_ID')
    dirs = fldtab.getcol('PHASE_DIR')
    fldtab.done()

    # Get time information
    maintab = table(myms, ack=False)
    scans = numpy.unique(maintab.getcol('SCAN_NUMBER'))

    # Build list of scan information
    scan_info = []
    for sc in scans:
        subtab = maintab.query(query='SCAN_NUMBER==' + str(sc))
        times = subtab.getcol('TIME')
        sfield = numpy.unique(subtab.getcol('FIELD_ID'))[0]
        scan_info.append({
            'scan': sc,
            'field_id': sfield,
            'times': times
        })
        subtab.done()

    maintab.done()

    # ------- CALCULATE ELEVATIONS PER SCAN -------

    mylogger.info('---- ELEVATION/AZIMUTH RANGES (PER SCAN):')
    mylogger.info('')
    mylogger.info('     SCAN      FIELD_ID  NAME                  RA[deg]       DEC[deg]      RA[hms]           DEC[dms]          EL_MIN[deg]   EL_MAX[deg]   EL_MEAN[deg]  AZ_MIN[deg]   AZ_MAX[deg]   AZ_MEAN[deg]  MJD_MIN       MJD_MAX       MJD_CENTRAL')
    mylogger.info('')

    for scan_data in scan_info:
        sc = scan_data['scan']
        fld_id = scan_data['field_id']
        times = scan_data['times']
        
        # Get field name and position
        name = names[fld_id]
        ra_rad = float(dirs[fld_id][0][0])
        dec_rad = float(dirs[fld_id][0][1])
        
        # Convert to degrees
        ra_deg = numpy.degrees(ra_rad)
        dec_deg = numpy.degrees(dec_rad)
        
        # Create SkyCoord for the source
        source_coord = SkyCoord(ra=ra_rad * u.rad, dec=dec_rad * u.rad, frame='icrs')
        
        # Calculate elevation for this scan
        obs_times_mjd = numpy.array(times) / 86400.0  # Convert to MJD
        obs_times = Time(obs_times_mjd, format='mjd')
        
        # Get time statistics in MJD
        mjd_min = numpy.min(obs_times_mjd)
        mjd_max = numpy.max(obs_times_mjd)
        mjd_central = (mjd_min + mjd_max) / 2.0
        
        # HMS/DMS (colon-separated)
        ra_hms = source_coord.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
        dec_dms = source_coord.dec.to_string(unit=u.deg, sep=':', precision=2, alwayssign=True, pad=True)

        # Calculate Alt/Az at each time point
        altaz_frame = AltAz(obstime=obs_times, location=meerkat)
        source_altaz = source_coord.transform_to(altaz_frame)
        
        # Get elevation statistics
        elevations = source_altaz.alt.degree
        el_min = numpy.min(elevations)
        el_max = numpy.max(elevations)
        el_mean = numpy.mean(elevations)

        # Get azimuth statistics
        azimuths = source_altaz.az.degree
        az_min = numpy.min(azimuths)
        az_max = numpy.max(azimuths)
        az_mean = numpy.mean(azimuths)
        
        mylogger.info('     %-10s%-10s%-22s%-14.4f%-14.4f%-18s%-18s%-14.2f%-14.2f%-14.2f%-14.2f%-14.2f%-14.2f%-14.6f%-14.6f%-14.6f' % 
                 (sc, fld_id, name, ra_deg, dec_deg, ra_hms, dec_dms, el_min, el_max, el_mean, az_min, az_max, az_mean, mjd_min, mjd_max, mjd_central))

    mylogger.info('')
    mylogger.info('Elevation ranges written to: ' + logfile)
    mylogger.info('')


if __name__ == '__main__':
    main()
