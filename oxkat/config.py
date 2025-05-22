#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import json
import os
import sys


# ------------------------------------------------------------------------
#
# Check for project_info file and get band
#

if os.path.isfile('project_info.json'):
    with open('project_info.json') as f:
        project_info = json.load(f)
    BAND = project_info['band']
else:
    BAND = 'not yet determined'


# ------------------------------------------------------------------------
#
# Paths for components and OUTPUTS
#

CWD = os.getcwd()
HOME = os.path.expanduser('~')

OXKAT = CWD+'/oxkat'
DATA = CWD+'/data'
TOOLS = CWD+'/tools'

GAINPLOTS = CWD+'/GAINPLOTS'
GAINTABLES = CWD+'/GAINTABLES'
IMAGES = CWD+'/IMAGES'
RESULTS = CWD+'/RESULTS'
LOGS = CWD+'/LOGS'
SCRIPTS = CWD+'/SCRIPTS'
VISPLOTS = CWD+'/VISPLOTS'
INTERVALS = CWD+'/INTERVALS'


# ------------------------------------------------------------------------
#
# Singularity settings
#

# Set to False to disable singularity entirely
USE_SINGULARITY = True

# If your data are symlinked and located in a path that singularity
# cannot see by default then set BIND to that path.
# If you wish to bind multiple paths then use a comma-separated list.
BIND = ''
BINDPATH = '$PWD,'+CWD+','+BIND

IDIA_CONTAINER_PATH = ['/idia/software/containers/',HOME+'/containers/']
CHPC_CONTAINER_PATH = [HOME+'/containers/']
HIPPO_CONTAINER_PATH = None
NODE_CONTAINER_PATH = [HOME+'/containers/', '/mnt/ephem/containers']


PYTHON3_PATTERN = 'polkat-0.1.2'
CASA_PATTERN = 'polkat-0.1.2'
QUARTICAL_PATTERN = 'polkat-0.1.2'
WSCLEAN_PATTERN = 'polkat-0.1.2'
SHADEMS_PATTERN = 'polkat-0.1.2'
ALBUS_PATTERN = 'polkat-albus'
TRICOLOUR_PATTERN = 'polkat-0.1.2'
MOVIE_PATTERN = 'polkat-0.1.2'


# ------------------------------------------------------------------------
#
# Slurm resource settings
#

SLURM_ACCOUNT = '' # e.g. b09-mightee-ag, b24-thunderkat-ag
SLURM_RESERVATION = '' # e.g. lsp-mightee

SLURM_NODELIST = '' # Specify node(s) to use
SLURM_EXCLUDE = '' # Specify node(s) to exclude


SLURM_DEFAULTS = {
	'TIME': '12:00:00',
	'PARTITION': 'Main',
	'NTASKS': '1',
	'NODES': '1',
	'CPUS': '8',
	'MEM': '64GB',
}

SLURM_TRICOLOUR = {
    'TIME': '06:00:00',
    'PARTITION': 'Main',
    'NTASKS': '1',
    'NODES': '1',
    'CPUS': '32',
    'MEM': '230GB'
}

SLURM_RM = {
	'TIME': '24:00:00',
	'PARTITION': 'Main',
	'NTASKS': '1',
	'NODES': '1',
	'CPUS': '8',
	'MEM': '64GB',
}

SLURM_WSCLEAN = {
    'TIME': '36:00:00',
    'PARTITION': 'Main',
    'NTASKS': '1',
    'NODES': '1',
    'CPUS': '32',
    'MEM': '230GB'
}

SLURM_PREDICT= {
    'TIME': '12:00:00',
    'PARTITION': 'Main',
    'NTASKS': '1',
    'NODES': '1',
    'CPUS': '32',
    'MEM': '230GB'
}

SLURM_EXTRALONG = {
    'TIME': '48:00:00',
    'PARTITION': 'Main',
    'NTASKS': '1',
    'NODES': '1',
    'CPUS': '32',
    'MEM': '230GB'
}

SLURM_HIGHMEM = {
    'TIME': '36:00:00',
    'PARTITION': 'HighMem',
    'NTASKS': '1',
    'NODES': '1',
    'CPUS': '32',
    'MEM': '480GB'
}

# ------------------------------------------------------------------------
#
# PBS resource settings
#

CHPC_ALLOCATION = 'ASTR1301'

PBS_DEFAULTS = {
	'PROGRAM': CHPC_ALLOCATION,
	'WALLTIME': '12:00:00',
	'QUEUE': 'serial',
	'NODES': '1',
	'PPN': '8',
	'MEM': '64gb'
}

PBS_TRICOLOUR = {
	'PROGRAM': CHPC_ALLOCATION,
	'WALLTIME': '06:00:00',
	'QUEUE': 'serial',
	'NODES': '1',
	'PPN': '24',
	'MEM': '120gb'
}

PBS_WSCLEAN = {
	'PROGRAM': CHPC_ALLOCATION,
	'WALLTIME': '12:00:00',
	'QUEUE': 'serial',
	'NODES': '1',
	'PPN': '24',
	'MEM': '120gb'
}

PBS_EXTRALONG = {
    'PROGRAM': CHPC_ALLOCATION,
    'WALLTIME': '48:00:00',
    'QUEUE': 'serial',
    'NODES': '1',
    'PPN': '24',
    'MEM': '120gb'
}

# ------------------------------------------------------------------------
#
# 1GC settings
#

# Scan intents, for automatic identification of cals/targets
CAL_1GC_TARGET_INTENT = 'TARGET'     # (partial) string to match for target intents
CAL_1GC_PRIMARY_INTENT = 'BANDPASS'  # (partial) string to match for primary intents
CAL_1GC_SECONDARY_INTENT = 'PHASE'   # (partial) string to match for secondary intents
CAL_1GC_DIAGNOSTICS = True          #  Choose if you want to make diagnostic plots of the Leakage + Phase cal
CAL_1GC_AGGRESSIVE_FLAGS = False     #  Choose if you want to aggresively flag the visibilities -- required for high-precision polarimetry (i.e., <1%)
CAL_1GC_RENAME_FEEDS = True # Turn true [Default] if you need switch feed naming; the SARAO Archive by default mislabels X as Y; and Y as X.  

# Pre-processing, operations applied when master MS is split to working MS
PRE_FIELDS = ''  # Comma-separated list of fields to select from raw MS
                                     # Names or IDs, do not mix, do not use spaces

# Polarization calibrator info --- Must be in PRE_FIELDS, if left blank will assume no polarization angle calibration
POLANG_NAME = 'J1331+3030'         # Specify the name of the field you want to use as a Polarization angle calibrator -- 3C286
POLANG_DIR  = '13:31:08.2881,+30.30.32.959' # CASA Format
POLANG_MOD  = [1.0, 0.0, 0.5, 0.0]

#POLANG_NAME = 'J0521+1638'         # Specify the name of the field you want to use as a Polarization angle calibrator -- 3C138
#POLANG_DIR  = '05:21:09.890000,+16.38.22.10000' # CASA Format
#POLANG_MOD  = [1.0, 0.3, -0.05, 0.0]

PRE_SCANS = ''                       # Comma-separated list of scans to select from raw MS
PRE_NCHANS = 1024                    # Integer number of channels for working MS
PRE_TIMEBIN = '8s'                   # Integration time for working MS

# Reference antennas
CAL_1GC_REF_ANT = 'auto'             # Comma-separated list to manually specify refant(s)
CAL_1GC_REF_POOL = ['m000','m001','m002','m003','m004','m006'] 
                                     # Pool to re-order for reference antenna list for 'auto'

# Field selection, IDs only at present. (Use tools/ms_info.py.)
CAL_1GC_PRIMARY = 'auto'             # Primary calibrator field ID
CAL_1GC_TARGETS = 'auto'             # Comma-separated target field IDs
CAL_1GC_SECONDARIES = 'auto'         # Comma-separated secondary IDs

                                     # - Lists of equal length in targets and secondaries maps cals to targets
                                     # - A single ID in uses same secondary for all targets
                                     # - A length mismatch reverts to auto, so double check!

# Sky model for primary calibrator --- EXPERIMENTAL (use 1GC_primary_models.py setup)
CAL_1GC_PRIMARY_MODEL = 'auto'       # setjy = use setjy component model only
                                     # auto = try to find a suitable model of the field sources in data/calmodels, defer to setjy if not found
                                     # or specify the location/of/wsclean-prefix for an arbitrary model cube


# GBK settings
CAL_1GC_DELAYCUT = 2.5               # [now defunct] Jy at central freq. Do not solve for K on secondaries weaker than this
CAL_1GC_FILLGAPS = 24                # Maximum channel gap over which to interpolate bandpass solutions

# Band specific options

if BAND == 'UHF':       

    CAL_1GC_FREQRANGE = '*:850~900MHz'        # Clean part of the band to use for generating UHF 1GC G-solutions
    CAL_1GC_UVRANGE = '>150m'               # Selection for baselines to include during 1GC B/G solving (K excluded)
    CAL_1GC_0408_MODEL = ([27.907,0.0,0.0,0.0],[-1.205],'850MHz')

    CAL_1GC_BAD_FREQS = ['*:540~570MHz',      # Lower band edge 
                        '*:1010~1150MHz']     # Upper band edge

    CAL_1GC_BL_FLAG_UVRANGE = '<600'        # Baseline range for which BL_FREQS are flagged
    CAL_1GC_BL_FREQS = []            

elif BAND == 'L':

    CAL_1GC_FREQRANGE = '*:1300~1400MHz'
    CAL_1GC_UVRANGE = '>150m'
    CAL_1GC_0408_MODEL = ([17.066,0.0,0.0,0.0],[-1.179],'1284MHz')

    CAL_1GC_BAD_FREQS = ['*:850~900MHz',      # Lower band edge
                        '*:1650~1800MHz',     # Upper bandpass edge
                        '*:1419.8~1421.3MHz'] # Galactic HI 

    CAL_1GC_BL_FLAG_UVRANGE = '<600'
    CAL_1GC_BL_FREQS = ['*:900MHz~915MHz',    # GSM and aviation
                        '*:925MHz~960MHz',                
                        '*:1080MHz~1095MHz',
                        '*:1565MHz~1585MHz',  # GPS
                        '*:1217MHz~1237MHz',
                        '*:1375MHz~1387MHz',
                        '*:1166MHz~1186MHz',
                        '*:1592MHz~1610MHz',  # GLONASS
                        '*:1242MHz~1249MHz',
                        '*:1191MHz~1217MHz',  # Galileo
                        '*:1260MHz~1300MHz',
                        '*:1453MHz~1490MHz',  # Afristar
                        '*:1616MHz~1626MHz',  # Iridium
                        '*:1526MHz~1554MHz',  # Inmarsat
                        '*:1600MHz']                 # Alkantpan
                                            # https://github.com/ska-sa/MeerKAT-Cookbook/blob/master/casa/L-band%20RFI%20frequency%20flagging.ipynb

elif BAND == 'S0':

    CAL_1GC_FREQRANGE = '*:2300~2400MHz'
    CAL_1GC_UVRANGE = '>150m'
    CAL_1GC_0408_MODEL = ([9.193,0.0,0.0,0.0],[-1.144],'2187MHz')   
    CAL_1GC_BAD_FREQS = ['*:1700~1800MHz',    # Lower band edge 
                        '*:2500~2650MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S1':

    CAL_1GC_FREQRANGE = ''
    CAL_1GC_UVRANGE = '>150m'
    CAL_1GC_0408_MODEL = ([8.244,0.0,0.0,0.0],[-1.138],'2406MHz')   
    CAL_1GC_BAD_FREQS = ['*:1967~2056MHz',    # Lower band edge 
                        '*:2756~2845MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S2':

    CAL_1GC_FREQRANGE = ''
    CAL_1GC_UVRANGE = '>150m'
    CAL_1GC_0408_MODEL = ([7.468,0.0,0.0,0.0],[-1.133],'2625MHz')   
    CAL_1GC_BAD_FREQS = ['*:2187~2275MHz',    # Lower band edge 
                        '*:2975~3063MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S3':

    CAL_1GC_FREQRANGE = ''
    CAL_1GC_UVRANGE = '>150m'
    CAL_1GC_0408_MODEL = ([6.822,0.0,0.0,0.0],[-1.128],'2483MHz')   
    CAL_1GC_BAD_FREQS = ['*:2405~2493MHz',    # Lower band edge 
                        '*:3194~3282MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S4':

    CAL_1GC_FREQRANGE = '*:2900~3000MHz'
    CAL_1GC_UVRANGE = '>150m'     
    CAL_1GC_0408_MODEL = ([6.423,0.0,0.0,0.0],[-1.124],'3000MHz')   
    CAL_1GC_BAD_FREQS = ['*:2600~2690MHz',    # Lower band edge 
                        '*:3420~3600MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600'
    CAL_1GC_BL_FREQS = []


# LINE modifiers
CAL_1GC_LINE_FILLGAPS = 48

# ------------------------------------------------------------------------
#
# 2GC settings
#


# CASA gaincal settings
CAL_2GC_UVRANGE = '>150m'            # Selection for baselines to include during G solving
CAL_2GC_PSOLINT = '32s'              # Solution interval for phase-only selfcal
CAL_2GC_APSOLINT = 'inf'             # Solution interval for amplitude and phase selfcal

# Quartical
CAL_2GC_YAML = DATA+'/quartical/2GC_phase.yaml'  # Frequency-dependent, phase-only self-calibration (diagonal terms: XX/YY only)
                                                 # NOTE: This does NOT de-rotate or re-apply parallactic angle corrections.
                                                 # WARNING: Since Stokes Q depends on both XX/YY and parallactic angle,
                                                 # omitting parallactic angle correction could, in principle, affect Q.
                                                 # However, in practice, no significant impact on Stokes Q has been observed in tests.

# These YAML files perform amplitude and phase self-calibration on the full 2x2 Jones matrix (XX, YY, XY, YX).
# Extensive testing indicates that using these options can significantly alter the measured polarized fluxes,
# likely due to strong instrumental polarization effects from the primary beam response.
# Without an accurate, full polarization beam model, these options SHOULD NOT BE USED.
# Proceed with extreme caution—you have been warned!
# CAL_2GC_YAML = DATA+'/quartical/2GC_fullpol.yaml'
# CAL_2GC_YAML = DATA+'/quartical/2GC_amppol.yaml'


# ------------------------------------------------------------------------
#
# wsclean and 2GC defaults
#
# General
WSC_MEM = 90
WSC_ABSMEM = -1 # in GB; mem is used if absmem is negative, calculated automatically for HPC, see absmem_helper
WSC_CONTINUE = False
# Outputs
WSC_MAKEPSF = False
WSC_NODIRTY = False
WSC_SOURCELIST = False
# Data selection
WSC_FIELD = 0
WSC_STARTCHAN = -1
WSC_ENDCHAN = -1
WSC_MINUVL = ''
WSC_MAXUVL = ''
WSC_EVEN = False
WSC_ODD = False
WSC_TUKEYTAPER = False
WSC_TAPERMASK = False
WSC_INTERVAL0 = None
WSC_INTERVAL1 = None
WSC_INTERVALSOUT = False
WSC_PARALLELREORDERING = 8
# Image dimensions
WSC_IMSIZE = 10240
WSC_CAL_IMSIZE = 2560
WSC_CELLSIZE = '1.1asec'
# Gridding / degridding
WSC_USEWGRIDDER = True
WSC_WGRIDDERACCURACY = 5e-5
WSC_BDA = False
WSC_BDAFACTOR = 10
WSC_NOMODEL = False
WSC_NWLAYERSFACTOR = 5
WSC_PADDING = 1.2
WSC_USEIDG = False # use image-domain gridder (not useable yet)
WSC_IDGMODE = 'CPU'
WSC_PREDICTCHANNELS = 64
WSC_PARALLELGRIDDING = 8
# Weighting
WSC_WEIGHT = 'briggs 0.0'
WSC_WEIGHT_CAL = 'uniform'
WSC_TAPERGAUSSIAN = ''
WSC_MFWEIGHT = False
# HIGH RES IMAGING
WSC_UNIFORM_IMAGE = True
WSC_WEIGHT_HIGHRES = 'uniform' # pick a more uniform weighting then WSC_WEIGHT -- uniform weight by default
# Deconvolution
WSC_PARALLELDECONVOLUTION = 2560
WSC_MULTISCALE = False
WSC_SCALES = '0,3,9'
WSC_MULTISCALE_BIAS = 0.7
WSC_CHANDECONV = False
WSC_NITER = 800000
WSC_GAIN = 0.15
WSC_MGAIN = 0.9
WSC_MASK_CHANNELSOUT = 8
WSC_IMAGE_CHANNELSOUT = 8
WSC_CAL_CHANNELSOUT = 16
WSC_MAX_CHANNELS = 16
WSC_FITSPECTRALPOL = 4
WSC_JOINCHANNELS = True
WSC_NONEGATIVE = False
WSC_STOPNEGATIVE = False
WSC_CIRCULARBEAM = False
WSC_POL = 'IQUV'
WSC_SPLITPOL = False # Image V/I and Q/U separately (necessary for High RM and MFS fitting)
WSC_JOINPOLARIZATIONS = True
WSC_SQUAREPOLARIZATIONS = False
# Masking
WSC_MASK = False
WSC_THRESHOLD = False
WSC_AUTOMASK = 3.0
WSC_AUTOTHRESHOLD = 1.0
WSC_LOCALRMS = False
# Determines if you want to Homogenize the resolution
WSC_HOMOGENIZEBEAM = False # Homogenize in freq
WSC_HOMOGENIZETIME = False # Homogenize in freq AND time
# Determine if you want to match the large resolvable angular scale by frequency
WSC_MATCHSCALES = WSC_HOMOGENIZEBEAM
if WSC_MATCHSCALES:
    speed_of_light = 299792458. # m / s
    min_baseline = 29. # m 
    # Calculate the largest resovable angular scale using top of band + 10 wavelengths as padding
    if BAND == 'UHF':
        minuvl = min_baseline / (speed_of_light / 1088.0e6) + 10. 
        WSC_MINUVL = '{}.0'.format(round(minuvl))
    if BAND == 'L':
        minuvl = min_baseline / (speed_of_light / 1712.0e6) + 10. 
        WSC_MINUVL = '{}.0'.format(round(minuvl))
    if BAND == 'S0':
        minuvl = min_baseline / (speed_of_light / 2625.0e6) + 10. 
        WSC_MINUVL = '{}.0'.format(round(minuvl))
    if BAND == 'S1':
        minuvl = min_baseline / (speed_of_light / 2843.0e6) + 10. 
        WSC_MINUVL = '{}.0'.format(round(minuvl))
    if BAND == 'S2':
        minuvl = min_baseline / (speed_of_light / 3062.0e6) + 10. 
        WSC_MINUVL = '{}.0'.format(round(minuvl))
    if BAND == 'S3':
        minuvl = min_baseline / (speed_of_light / 3281.0e6) + 10. 
        WSC_MINUVL = '{}.0'.format(round(minuvl))
    if BAND == 'S4':
        minuvl = min_baseline / (speed_of_light / 3500.0e6) + 10. 
        WSC_MINUVL = '{}.0'.format(round(minuvl))

# Band modifiers
if BAND == 'UHF':
    WSC_CELLSIZE = '1.7asec'
    WSC_BRIGGS = -0.5
    WSC_BDAFACTOR = 4
    WSC_NWLAYERSFACTOR = 5
if BAND == 'S0':
    WSC_CELLSIZE = '0.65asec'
if BAND == 'S1':
    WSC_CELLSIZE = '0.61asec'
if BAND == 'S2':
    WSC_CELLSIZE = '0.58asec'
if BAND == 'S3':
    WSC_CELLSIZE = '0.54asec'    
if BAND == 'S4':
    WSC_CELLSIZE = '0.5asec'


# ------------------------------------------------------------------------
#
# 3GC peeling settings
#

CAL_3GC_PEEL_NCHAN = 32
CAL_3GC_PEEL_MAXCHAN = WSC_MAX_CHANNELS
CAL_3GC_PEEL_POL = WSC_POL
CAL_3GC_PEEL_BRIGGS = 'briggs -0.6'
CAL_3GC_PEEL_REGION = ''  # Specify DS9 peeling region 
                          # Leave blank to search for <fieldname>*peel*.reg in the current path
CAL_3GC_PEEL_YAML = DATA+'/quartical/3GC_peel.yaml' # Diagonol only entries for dE (direction dependant peeling solutions) + frequency-dependant phase for full field
                                                    # Note that this doesn't de-rotate the parallactic angle, may cause issues with Stokes Q (but I haven't seen any yet)



CAL_3GC_FACET_REGION = '' # Specify DS9 region to define tessel centres
                          # Leave blank to search for <fieldname>*facet*.reg in the current path
                          # Regions specified here and above will apply to all fields, and so can
                          # be used to e.g. peel the same source from a compact mosaic rather than
                          # having to provide multiple copies of the same region on a per-field basis

# ------------------------------------------------------------------------
#
# MakeMask defaults
#


MAKEMASK_THRESH = 6.0
MAKEMASK_BOXSIZE = 500
MAKEMASK_SMALLBOX = 50
MAKEMASK_ISLANDSIZE = 30000
MAKEMASK_DILATION = 3

# ------------------------------------------------------------------------
#
# BREIZORRO defaults
#

BREIZORRO_THRESH = 6.0
BREIZORRO_BOXSIZE = 50
BREIZORRO_FILLHOLES = True
BREIZORRO_DILATION = 3

# ------------------------------------------------------------------------
#
# DDFacet defaults
#


# [Data]
DDF_DDID = 'D*'
DDF_FIELD = 'F0'
DDF_COLNAME = 'CORRECTED_DATA'
DDF_CHUNKHOURS = 0.5
DDF_DATASORT = 1
# [Predict]
DDF_PREDICTCOLNAME = '' # MODEL_DATA or leave empty to disable predict
DDF_INITDICOMODEL = ''
# [Output]
DDF_OUTPUTALSO = 'oenNS'
DDF_OUTPUTIMAGES = 'DdPMmRrIikz' # add 'A' to re-include spectral index map
DDF_OUTPUTCUBES = 'MmRi' # output intrinsic and apparent resid and model cubes
# [Image]
DDF_NPIX = 10125
DDF_CELL = 1.1
# [Facets]
DDF_DIAMMAX = 0.25
DDF_DIAMMIN = 0.05
DDF_NFACETS = 4 # crank this up (32?) to get better beam resolution if FITS beam is used
DDF_PSFOVERSIZE = 1.5
DDF_PADDING = 3.0 # padding needs increasing from default if NFacets is raised to prevent aliasing
# [Weight]
DDF_ROBUST = 0.0
# [Convolution Functions]
# DDF_NW = 100 # Increase for strong off-axis sources
# [Comp]
DDF_SPARSIFICATION = '0' # [100,30,10] grids every 100th visibility on major cycle 1, every 30th on cycle 2, etc.
# [Parallel]
DDF_NCPU = 8
# [Cache]
DDF_CACHERESET = 0
DDF_CACHEDIR = '.'
DDF_CACHEHMP = 1
# [Beam]
DDF_BEAM = '' # specify beam cube of the form: meerkat_pb_jones_cube_95channels_$(xy)_$(reim).fits
DDF_BEAMNBAND= 10
DDF_DTBEAMMIN = 1
DDF_FITSPARANGLEINCDEG = 0.5
DDF_BEAMCENTRENORM = True
DDF_FEEDSWAP = 1
DDF_BEAMSMOOTH = False
# [Freq]
DDF_NBAND = 8
DDF_NDEGRIDBAND = 8
# [DDESolutions]
DDF_DDSOLS = ''
DDF_DDMODEGRID = 'AP'
DDF_DDMODEDEGRID = 'AP'
# [Deconv]
DDF_GAIN = 0.15
DDF_FLUXTHRESHOLD = 3e-6
DDF_CYCLEFACTOR = 0
DDF_RMSFACTOR = 3.0	
DDF_DECONVMODE = 'hogbom'
DDF_SSD_DECONVPEAKFACTOR = 0.001
DDF_SSD_MAXMAJORITER = 3
DDF_SSD_MAXMINORITER = 120000
DDF_SSD_ENLARGEDATA = 0
DDF_HOGBOM_DECONVPEAKFACTOR = 0.1
DDF_HOGBOM_MAXMAJORITER = 5
DDF_HOGBOM_MAXMINORITER = 100000
DDF_HOGBOM_POLYFITORDER = 4
# [Mask]
DDF_MASK = 'auto' # 'auto' enables automasking 
                  # 'fits' uses the first *.mask.fits in the current folder
                  # otherwise pass a filename to use a specific FITS image
# [Misc]
DDF_MASKSIGMA = 4.5
DDF_CONSERVEMEMORY = 1


# Band modifiers
if BAND == 'UHF':
    DDF_CELL = 1.7
    DDF_ROBUST = -0.5
if BAND == 'S0':
    DDF_CELL = 0.65
if BAND == 'S1':
    DDF_CELL = 0.61
if BAND == 'S2':
    DDF_CELL = 0.58
if BAND == 'S3':
    DDF_CELL = 0.54
if BAND == 'S4':
    DDF_CELL = 0.5


# ------------------------------------------------------------------------
#
# killMS defaults
#


# [VisData]
KMS_TCHUNK = 0.5
KMS_INCOL = 'CORRECTED_DATA'
KMS_OUTCOL = 'MODEL_DATA'
# [Beam]
KMS_BEAM = ''
KMS_BEAMAT = 'Facet'
KMS_DTBEAMMIN = 1
KMS_CENTRENORM = 1
KMS_NCHANBEAMPERMS = 95
KMS_FITSPARANGLEINCDEG = 0.5
KMS_FITSFEEDSWAP = 1
# [ImageSkyModel]
KMS_DICOMODEL = ''
KMS_MAXFACETSIZE = 0.25
# [DataSelection]
KMS_UVMINMAX = '0.15,8500.0'
KMS_FIELDID = 0
KMS_DDID = 0
# [Actions]
KMS_NCPU = 16
KMS_DOBAR = 0
KMS_DEBUGPDB = 0
# [Solvers]
KMS_SOLVERTYPE = 'kafca'
KMS_DT = 5
KMS_NCHANSOLS = 8
# [KAFCA]
KMS_NITERKF = 9
KMS_COVQ = 0.05

# ------------------------------------------------------------------------
#
# Snapshot imaging defaults
#

SNAP_FIELDS = '' # Comma-separated list of field to run snapshot imaging on
SNAP_CHANNELSOUT = WSC_IMAGE_CHANNELSOUT # Integer number of channels to perform snap-shot imaging -- by default the same as the WSC channels out
SNAP_INTBIN = 1 # Integer number of intervals to image together during snapshot imaging (default = 1 is per intergration imaging)
SNAP_INTEND = True # This will throw away the last interval if 'incomplete' compared to bin
    # E.G., total ints = 17, int bin = 4, will only image integrations 1 to 16 (4 total images) discarding 17
SNAP_POL = True # Bool for whether or not to do full polarisation snapshot images
SNAP_IMSIZE = 2560 # Image size of snapshot images (not the model image)
SNAP_MODELIDENTIFIER = 'pcalmask' # identifier for image name following oxkat/polkat conventions
SNAP_MODELMASK = '' # Point to mask for initial model creation (or just don't delete the pcalmask model)
SNAP_DECONV = False # Deconvolve during the snapshot imaging process? -- rarely necessary
SNAP_DECONVMASK = '' # mask to use when doing snapshot imaging and deconvolving

# ------------------------------------------------------------------------
#
# PyBDSF defaults
#


PYBDSF_THRESH_PIX = 5.0
PYBDSF_THRESH_ISL = 3.0
PYBDSF_CATALOGTYPE = 'srl'
PYBDSF_CATALOGFORMAT = 'fits'


# ------------------------------------------------------------------------
#
# ClusterCat defaults
#


CLUSTERCAT_NDIR = 8
CLUSTERCAT_CENTRALRADIUS = 0.0
CLUSTERCAT_NGEN = 100
CLUSTERCAT_FLUXMIN = 0.000001
CLUSTERCAT_NCPU = 8


# ------------------------------------------------------------------------
#
# MeerKAT primary beam models
#


BEAM_L = HOME+'/Beams/meerkat_pb_jones_cube_95channels_$(xy)_$(reim).fits'


# ------------------------------------------------------------------------
#
# ALBUS GPS Data extraction method
#


#RED_TYPE = 'RI_G01'
RED_TYPE = 'RI_G03'
