#
# A basic python script that tracks a specified position on the 
# sky over the time range from START_TIME to END_TIME from
# a specific location on the Earth's surface.

# The output is a text file giving Slant Tec (STEC) and
# ionosphere rotation measure (RM) as a function of time

import os, subprocess, sys
import time
from argparse import ArgumentParser
import os.path as o

import MS_Iono_functions as iono 
import math

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg


def main():
    
    # Load in the Arguments
    parser = ArgumentParser(description='Run ALBUS to get the Ionospheric RM along any arbitrary LOS and at any arbitrary time')
    parser.add_argument('field',
                        help="Name of target field")
    parser.add_argument('-r', '--RA', dest="ra",
                        help="Right acension of source (DD:MM:SS.SSS)")
    parser.add_argument('-d', '--DEC', dest='dec',
                        help="Declination of source (HH:MM:SS.SSS)")
    parser.add_argument('-s', '--start-time', dest='t0',
                        help="Start time of source (YYYY/MM/DD hh:mm:ss.ss)")
    parser.add_argument('-f', '--final-time', dest='tf',
                        help="Declination of source (YYYY/MM/DD hh:mm:ss.ss)")

    args = parser.parse_args()
    OBJECT = f"meerkat_"+ str(args.field)
    RA = str(args.ra) 
    DEC = str(args.dec).replace('minus','-').replace('plus','+')
    START_TIME = str(args.t0)
    END_TIME = str(args.tf)

    # MEERKAT telescope location
    LONG="21:26:35.736"
    LAT="-30:42:44.838"
    HEIGHT=1059.662443

    # Default Solution parmaters from ALBUS examples -- Advise against touching this
    RED_TYPE = 'RI_G01'
    TIME_STEP = 300
    MAX_DIST = 750E3
    NUM_PROCESSORS = 8
    DO_SER = 0

    # Intialize a directory to store the (larger number) of GPS output files 
    GPS_DATA_DIR = cfg.RESULTS + '/meerkat_gps_data'
    if os.path.exists(GPS_DATA_DIR) is False:
        os.mkdir(GPS_DATA_DIR)

    
    # Process Ionospheric data 
    iono.process_ionosphere(time_step=TIME_STEP,
                                                    object=OBJECT,
                                                    Ra=RA,Dec=DEC,
                                                    Lat=LAT,
                                                    Long=LONG,
                                                    Height=HEIGHT,
                                                    start_time=START_TIME,
                                                    end_time=END_TIME,
                                                    max_dist=MAX_DIST,
                                                    processing_option=RED_TYPE,
                                                    do_serial=DO_SER,
                                                    num_processors=NUM_PROCESSORS, 
                                                    gps_data_directory=GPS_DATA_DIR)
    

if __name__ == "__main__":
    main()

