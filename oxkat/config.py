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
BINDPATH = '$PWD,'+CWD+',/mnt,'+BIND

IDIA_CONTAINER_PATH = ['/idia/software/containers/',HOME+'/containers/', CWD, HOME]
CHPC_CONTAINER_PATH = [HOME+'/containers/']
HIPPO_CONTAINER_PATH = None
NODE_CONTAINER_PATH = [HOME+'/containers/', CWD, HOME]
GLAM_CONTAINER_PATH = ['/mnt/users/hughesa', '/mnt/extraspace/thunderkat_pol', CWD]


PYTHON3_PATTERN = 'polkat-0.2.5'
CASA_PATTERN = PYTHON3_PATTERN
QUARTICAL_PATTERN = PYTHON3_PATTERN
WSCLEAN_PATTERN = PYTHON3_PATTERN
SHADEMS_PATTERN = PYTHON3_PATTERN
ALBUS_PATTERN = 'polkat-albus'
SPINIFEX_PATTERN = PYTHON3_PATTERN
TRICOLOUR_PATTERN = PYTHON3_PATTERN
MOVIE_PATTERN = PYTHON3_PATTERN


# ------------------------------------------------------------------------
#
# Slurm resource settings
#

SLURM_ACCOUNT = '' # e.g. b09-mightee-ag, b24-thunderkat-ag
SLURM_RESERVATION = '' # e.g. lsp-mightee

SLURM_NODELIST = '' # Specify node(s) to use
SLURM_EXCLUDE = '' # Specify node(s) to exclude


SLURM_DEFAULTS = {
	'TIME': '16:00:00',
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
# GLAM resource settings (simplified SLURM-like for Glamdring cluster)
#

GLAM_SMALL = {
    'CPUS': '2',
    'MEM': '16GB',
}

GLAM_STANDARD = {
    'CPUS': '2',
    'MEM': '64GB',
}

GLAM_MEDIUM = {
    'CPUS': '8',
    'MEM': '120GB',
}

GLAM_WSC_SOURCE = {
    'CPUS': '12',
    'MEM': '180GB',
}

GLAM_WSC_CAL = {
    'CPUS': '8',
    'MEM': '64GB',
}

GLAM_CASA = {
    'CPUS': '8',
    'MEM': '64GB',
}

GLAM_LARGE = {
    'CPUS': '24',
    'MEM': '216GB',
}

GLAM_COREHEAVY = {
    'CPUS': '24',
    'MEM': '72GB',
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
PRE_SCANS = ''                       # Comma-separated list of scans to select from raw MS
PRE_NCHANS = 1024                    # Integer number of channels for working MS
PRE_TIMEBIN = '8s'                   # Integration time for working MS



# Polarization calibrator info --- Must be in PRE_FIELDS, if left blank will assume no polarization angle calibration
# https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol
POLANG_NAME = ''         # Specify the name of the field you want to use as a Polarization angle calibrator -- 3C286
POLANG_DIR  = '13:31:08.2881,+30.30.32.959' # CASA Format
POLANG_MOD  = [1.0, 0.0, 0.5, 0.0]
XF_TARGET_POLANG = 30.0  # Expected INTRINSIC (i.e. correcting RM effects) linear polarization angle in degrees
XF_TARGET_RM = 0.0 # Guess for the intrinsic RM; for 'manual' XF determination and RM trialing

#POLANG_NAME = 'J0521+1638'         # Specify the name of the field you want to use as a Polarization angle calibrator -- 3C138
#POLANG_DIR  = '05:21:09.890000,+16.38.22.10000' # CASA Format
#POLANG_MOD  = [1.0, 0.3, -0.05, 0.0]
#XF_TARGET_POLANG = -10.0  # Expected INTRINSIC (i.e. correcting RM effects) linear polarization angle in degrees
#XF_TARGET_RM = 0.0 # Guess for the intrinsic RM; for 'manual' XF determination and RM trialing

# PARAMETERS TO DECIDE XF SOLVE METHOD -- DEFAULTS ARE LIKELY ALL GOOD
XF_MODE = 'auto' # options are: auto (RECOMMENDED: determine based on band and/or if there is large phase discontinuties), casa or  manual
XF_AUTO_ANG_JUMP = 60.0 # angle in degrees where, if the CASA XF solver has adjacent solution intervals that have a discontinuity
                        # larger than this value, it will solve XF manually.
XF_SKIP_KCROSS = True # weather to skip the cross-hand delay and only solve for the phase
                      # my experimentation has found that some times this exclusion can give better solutions.

# XF table solution interval
XF_CHANINT = 16             # Channels per polcal solution interval; e.g. 1024 ch / 16 = 64 XF phase solutions across the band

# Calibration and outlier thresholds
XF_GLOBAL_SIGMA_CLIP = 10.0 # N-sigma threshold for global circular phase clip (gross outlier removal, wrap-safe)
XF_POLY_ORDER        = 4    # Polynomial order for cos/sin(phase) baseline fit; set to 0 to disable
XF_MAD_WINDOW        = None # Sliding circular MAD window size in channels; None = auto (ceil_odd(N_chan/10))
XF_MAD_SIGMA_CLIP    = 3.5  # N-sigma threshold for sliding circular MAD clip (local, wrap-safe)
XF_AVG_SCAN          = False # If True, flux-weighted circular average across scans after per-scan solve (increases SNR, loses scan-to-scan variation)
XF_MAD_ITERATIVE = True  # Iterative poly-baseline + sliding circular MAD clip until convergence (refits baseline each pass); False = single pass only
XF_MIN_SN = 5.0  # Minimum cross-hand S/N per channel; channels below this are flagged (0 = disabled)

# Smoothing control
XF_USE_SMOOTHING     = False # Apply Savitzky-Golay smoothing to XF solutions before writing the table
XF_SAVGOL_WINDOW     = None  # Smoothing window length in channels (odd integer); None = auto
XF_SAVGOL_POLYORDER  = 3     # Polynomial order for Savitzky-Golay filter (must be < XF_SAVGOL_WINDOW)

# Post-solve adjustments applied in 1GC_05_casa_refcal.py after XF solving
XF_REFERENCE_PHASE = ''  # Resolve the +/-180 deg degeneracy: each unflagged channel is flipped
                          # if the flipped value is closer to this reference (degrees).
                          # Applied to whichever Xf table is active (solved or override).
                          # Empty string = skip.

XF_OVERRIDE_TABLE  = ''  # If non-empty, replace the solved Xf table with this path before
                          # applycal.  Useful for importing a pre-computed table from another
                          # run (e.g. from global_crosshand_phase.py).
                          # Empty string = use the solved table.

# -----------------------------------------------------------------------
# DEPRECATED — parameters below are no longer used in the current
# manual XF solver (tools/manual_XF_solver.py). Retained for
# compatibility with older workflows or the CASA-based XF path.
# -----------------------------------------------------------------------
XF_POLY_SIGMA_CLIP   = 3.0   # Was: sigma threshold for iterative polynomial outlier clipping
XF_MAX_AVG_CHANNELS  = None  # Was: upper limit on channel averaging for XF solutions
XF_MIN_CROSS_FLUX    = 0.07  # Was: minimum total cross-hand flux (Jy) threshold for reliable solutions
XF_CLIP_WINDOW       = 16    # Was: window size (channels) for local scatter analysis
XF_EX                = True  # Was: enable gap-filling by extrapolation into flagged regions
XF_EX_FRAC           = 0.2   # Was: fraction of good bandwidth to use for extrapolation
XF_RM_ANG_MAX        = 20.0  # Was: max angular separation [deg] from XF_TARGET_POLANG for a valid RM trial solution
XF_DELTARM_TRIALS    = [-7,7] # Was: [min, max] RM offsets (rad/m²) defining the atmospheric RM trial grid

# Reference antennas
CAL_1GC_REF_ANT = 'auto'             # Comma-separated list to manually specify refant(s)
CAL_1GC_REF_POOL = ['m060','m061','m059','m058', 'm048', 'm062', 'm063'] 
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
CAL_1GC_FILLGAPS = 8                 # Maximum channel gap over which to interpolate bandpass solutions

# Parallactic angle correction during 1GC.
# If True, CASA applies the P-Jones correction as part of the 1GC solve.
# This is safe for phase-only and standard polarisation calibration.
# For high-fidelity amplitude self-calibration, set False: the feed-frame
# solve should remain in the antenna frame, with parallactic angle applied
# only as the final step to rotate into the sky frame.
CAL_1GC_APPLYPARANG = True

# Band specific options
CAL_1GC_BL_MODE = 'freq' # options are 'freq' to flag specific frequencies the 'uvrange' baselin
                         # 'aggressive' to the uvrange baseline on all frequencis,
                         # 'none' to not do any baseline specific flagging (initially)
                         # I've found that baseline flagging can cause artificiall spectral jumps, 
                         # so 'none' or 'aggressive' are the way to go (I think) whereas freq might be best for MFS (continuum) imaging.
                         # if you give the wrong option it will default to 'freq' and print a warning

if BAND == 'UHF':       

    CAL_1GC_FREQRANGE = '*:850~900MHz'        # Clean part of the band to use for generating UHF 1GC G-solutions
    CAL_1GC_UVRANGE = '>300m'               # Selection for baselines to include during 1GC B/G solving (K excluded)
    CAL_1GC_0408_MODEL = ([27.907,0.0,0.0,0.0],[-1.205],'850MHz') # Defunct, Model now hardcoded into 1GC_0

    CAL_1GC_BAD_FREQS = ['*:540~570MHz',      # Lower band edge 
                        '*:1010~1150MHz']     # Upper band edge

    CAL_1GC_BL_FLAG_UVRANGE = '<600m'        # Baseline range for which BL_FREQS are flagged
    CAL_1GC_BL_FREQS = []            

elif BAND == 'L':

    CAL_1GC_FREQRANGE = '*:1326~1357MHz'
    CAL_1GC_UVRANGE = '>300m'
    CAL_1GC_0408_MODEL = ([17.066,0.0,0.0,0.0],[-1.179],'1284MHz') # Defunct, Model now hardcoded into 1GC_05

    CAL_1GC_BAD_FREQS = ['*:850~900MHz',      # Lower band edge
                        '*:1650~1800MHz',     # Upper bandpass edge
                        '*:1419.8~1421.3MHz'] # Galactic HI 

    CAL_1GC_BL_FLAG_UVRANGE = '<600m'
    CAL_1GC_BL_FREQS = ['*:900MHz~915MHz',    # GSM and aviation
                        '*:925MHz~960MHz',
                        '*:935.40~960.05 MHz', # GSM (from katdal)
                        '*:1080MHz~1095MHz',
                        '*:1145.85~1300.90 MHz', # GPS/GLONASS (from katdal)
                        '*:1166MHz~1186MHz',
                        '*:1191MHz~1217MHz',  # Galileo
                        '*:1217MHz~1237MHz',
                        '*:1242MHz~1249MHz',
                        '*:1260MHz~1300MHz',
                        '*:1375MHz~1387MHz',
                        '*:1453MHz~1490MHz',  # Afristar
                        '*:1519.35~1608.30 MHz', # GPS/GLONASS (from katdal)
                        '*:1526MHz~1554MHz',  # Inmarsat
                        '*:1565MHz~1585MHz',  # GPS
                        '*:1592MHz~1610MHz',  # GLONASS
                        '*:1600MHz',                 # Alkantpan
                        '*:1616MHz~1626MHz']  # Iridium
                                            # https://github.com/ska-sa/MeerKAT-Cookbook/blob/master/casa/L-band%20RFI%20frequency%20flagging.ipynb

elif BAND == 'S0':

    CAL_1GC_FREQRANGE = '*:2300~2400MHz'
    CAL_1GC_UVRANGE = '>300m'
    CAL_1GC_0408_MODEL = ([9.193,0.0,0.0,0.0],[-1.144],'2187MHz') # Defunct, Model now hardcoded into 1GC_0
    CAL_1GC_BAD_FREQS = ['*:1700~1800MHz',    # Lower band edge 
                        '*:2500~2650MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600m'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S1':

    CAL_1GC_FREQRANGE = ''
    CAL_1GC_UVRANGE = '>300m'
    CAL_1GC_0408_MODEL = ([8.244,0.0,0.0,0.0],[-1.138],'2406MHz')    # Defunct, Model now hardcoded into 1GC_0
    CAL_1GC_BAD_FREQS = ['*:1967~2056MHz',    # Lower band edge 
                        '*:2756~2845MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600m'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S2':

    CAL_1GC_FREQRANGE = ''
    CAL_1GC_UVRANGE = '>300m'
    CAL_1GC_0408_MODEL = ([7.468,0.0,0.0,0.0],[-1.133],'2625MHz')    # Defunct, Model now hardcoded into 1GC_0
    CAL_1GC_BAD_FREQS = ['*:2187~2275MHz',    # Lower band edge 
                        '*:2975~3063MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600m'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S3':

    CAL_1GC_FREQRANGE = ''
    CAL_1GC_UVRANGE = '>300m'
    CAL_1GC_0408_MODEL = ([6.822,0.0,0.0,0.0],[-1.128],'2483MHz')   # Defunct, Model now hardcoded into 1GC_0
    CAL_1GC_BAD_FREQS = ['*:2405~2493MHz',    # Lower band edge 
                        '*:3194~3282MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600m'
    CAL_1GC_BL_FREQS = []

elif BAND == 'S4':

    CAL_1GC_FREQRANGE = '*:2900~3000MHz'
    CAL_1GC_UVRANGE = '>300m'     
    CAL_1GC_0408_MODEL = ([6.423,0.0,0.0,0.0],[-1.124],'3000MHz')   # Defunct, Model now hardcoded into 1GC_0
    CAL_1GC_BAD_FREQS = ['*:2600~2790MHz',    # Lower band edge 
                        '*:3375~3600MHz']     # Upper band edge
    CAL_1GC_BL_FLAG_UVRANGE = '<600m'
    CAL_1GC_BL_FREQS = []


# LINE modifiers
CAL_1GC_LINE_FILLGAPS = 48

# ------------------------------------------------------------------------
#
# 2GC settings — QuartiCal self-calibration
#

# Legacy CASA gaincal solution intervals (not used by QuartiCal; kept for reference)
CAL_2GC_UVRANGE = '>150m'            # UV range for baseline selection during gain solving
CAL_2GC_PSOLINT = '32s'              # Phase-only solution interval (DEFUNCT in QuartiCal version)
CAL_2GC_APSOLINT = 'inf'             # Amplitude+phase solution interval (DEFUNCT in QuartiCal version)

# QuartiCal YAML configurations
# See the individual .yaml files and the QuartiCal docs for full parameter details.
# Phase-only: solves diagonal Jones terms (XX, YY) for phase only.
CAL_2GC_YAML         = DATA+'/quartical/2GC_delay.yaml'

# Complex: solves the full 2x2 Jones matrix (XX, YY, XY, YX) for amplitude and phase.
# See the .yaml file for the exact solver configuration.
CAL_2GC_YAML_COMPLEX = DATA+'/quartical/2GC_phase.yaml'
# To use phase-only calibration for both stages (i.e. override the complex solve), uncomment:
# CAL_2GC_YAML_COMPLEX = CAL_2GC_YAML

# Skip primary beam correction in post-processing
SKIP_PB = False

# ------------------------------------------------------------------------
#
# 2GC settings — WSClean imaging
#

# --- Memory ---
WSC_MEM    = 90    # Memory limit as a percentage of total system RAM (-mem flag)
WSC_ABSMEM = -1    # Absolute memory limit in GB (-abs-mem flag); if negative, WSC_MEM % is used instead
                   # (calculated automatically for HPC configs via absmem_helper)

# --- Run control ---
WSC_CONTINUE = False  # Continue a previous deconvolution run (-continue flag)

# --- Output products ---
WSC_MAKEPSF    = False  # Force PSF image to be written even when not otherwise required
WSC_NODIRTY    = False  # Suppress writing dirty (unconvolved) images
WSC_SOURCELIST = False  # Save a component source list alongside the restored image

# --- Data selection ---
WSC_FIELD     = 0      # MS field ID to image
WSC_STARTCHAN = -1     # First channel to image (-channel-range); -1 = use all
WSC_ENDCHAN   = -1     # Last channel to image (-channel-range); -1 = use all
WSC_MINUVL    = ''     # Minimum UV distance in wavelengths (-minuv-l); '' = no cut
WSC_MAXUVL    = ''     # Maximum UV distance in wavelengths (-maxuv-l); '' = no cut
WSC_EVEN      = False  # Image even time-slots only (-even-timesteps)
WSC_ODD       = False  # Image odd time-slots only (-odd-timesteps)
WSC_INTERVAL0 = None   # First time interval to image (-interval start)
WSC_INTERVAL1 = None   # Last time interval to image (-interval end)
WSC_INTERVALSOUT = False  # Number of output time intervals (-intervals-out); False = single image

# --- Inner UV taper (suppress short-spacing artefacts) ---
WSC_TUKEYTAPER  = 3000.0  # Inner Tukey taper in wavelengths (-taper-inner-tukey);
                          # suppresses emission on scales larger than ~1/taper radians
WSC_TAPERMASK   = True    # If False, WSC_MINUVL/MAXUVL/TUKEYTAPER applied only to the
                          # final self-calibrated image; masking and self-cal use full UV
                          # range. If True, UV cuts also applied during blind masking and
                          # self-calibration.

# --- Parallelism ---
WSC_PARALLELREORDERING = 8  # Number of parallel threads for MS reordering (-parallel-reordering)
WSC_PARALLELGRIDDING   = 8  # Number of parallel gridding threads (-parallel-gridding)

# --- Image dimensions ---
WSC_IMSIZE     = 8500       # Image size in pixels for science targets (must be even)
WSC_CAL_IMSIZE = 2560        # Image size in pixels for calibrators
WSC_CELLSIZE   = '1.1asec'   # Pixel scale (arcsec/pixel); overridden per-band below

# --- Gridding --- DEPRECATED
WSC_USEWGRIDDER    = True   # Use the w-gridder algorithm (-gridder wgridder); recommended
WSC_WGRIDDERACCURACY = 5e-5 # w-gridder accuracy parameter (lower = more accurate, slower)
WSC_BDA            = False  # Enable baseline-dependent averaging during gridding
WSC_BDAFACTOR      = 10     # BDA smearing factor (-baseline-averaging); only used if WSC_BDA=True
WSC_NOMODEL        = False  # Skip writing model visibilities (-no-update-model-required)
WSC_NWLAYERSFACTOR = 5      # w-layers scaling factor (-nwlayers-factor); w-projection only
WSC_PADDING        = 1.2    # Image padding factor (-padding); w-projection only
WSC_USEIDG         = False  # Use image-domain gridder (experimental, not yet usable)
WSC_IDGMODE        = 'CPU'  # IDG compute mode: CPU or GPU
WSC_PREDICTCHANNELS = 64    # Number of channels used during model prediction

# --- Weighting ---
WSC_WEIGHT       = 'briggs 0.0'  # Visibility weighting scheme for science targets
WSC_WEIGHT_CAL   = 'uniform'     # Visibility weighting scheme for calibrators
WSC_TAPERGAUSSIAN = ''           # Gaussian UV taper (-taper-gaussian); '' = disabled
WSC_MFWEIGHT     = False          # Apply multi-frequency weighting (-mf-weighting / -no-mf-weighting)

# High-resolution uniform-weight image (made alongside the main Briggs image)
WSC_UNIFORM_IMAGE  = True       # Generate an additional uniform-weight image
WSC_WEIGHT_HIGHRES = 'uniform'  # Weighting scheme for the high-resolution image

# --- Deconvolution ---
WSC_PARALLELDECONVOLUTION = 2560  # Tile size in pixels for parallel deconvolution (-parallel-deconvolution)
WSC_MULTISCALE       = False      # Enable multi-scale CLEAN (-multiscale)
WSC_SCALES           = '0,3,9'   # Multi-scale CLEAN scale sizes in pixels
WSC_MULTISCALE_BIAS  = 0.7       # Multi-scale bias parameter; higher = prefer smaller scales
WSC_CHANDECONV       = 16      # Number of sub-band images used during deconvolution
                                  # (--deconvolution-channels); False = same as channels-out
WSC_NITER  = 800000  # Maximum number of CLEAN minor-cycle iterations (-niter)
WSC_GAIN   = 0.15    # Minor-cycle loop gain; fraction of peak subtracted per iteration (-gain)
WSC_MGAIN  = 0.9     # Major-cycle loop gain; fraction of residual cleaned per major cycle (-mgain)
WSC_NONEGATIVE  = False  # Forbid negative CLEAN components (-no-negative)
WSC_STOPNEGATIVE = False # Stop deconvolution when the first negative component is found (-stop-negative)
WSC_CIRCULARBEAM = False # Force a circular (not elliptical) restoring beam (-circular-beam)

# --- Channel strategy ---
WSC_BLIND_CHANNELSOUT = 8     # Blind (shallow, no mask) image — used only to generate initial mask
WSC_PCAL_CHANNELSOUT  = 128     # Final self-calibrated (pcalmask) image
WSC_DMASK_CHANNELSOUT = 128  # Data-masked image (stage 1 selfcal model);
                                               # set independently if a different channelisation is needed
WSC_CAL_CHANNELSOUT   = WSC_PCAL_CHANNELSOUT    # Calibrator images (secondaries / primary)

# Memory-safe channel chunking: if WSC_MAX_CHANNELS < channels-out, wsclean is called
# multiple times in channel-range blocks of WSC_MAX_CHANNELS, then the results are
# assembled into a full cube. Must divide evenly into both PRE_NCHANS and channels-out.
WSC_MAX_CHANNELS = 64

# Spectral polynomial fitting across output channels (-fit-spectral-pol N)
# Applied during deconvolution to constrain spectral index. 0 = disabled.
WSC_FITSPECTRALPOL = 4

# --- Channel joining and polarisation ---
# Join all output channels during peak finding to improve sensitivity (-join-channels)
WSC_JOINCHANNELS = True

# Use sum-of-squares (squared channel joining) for peak finding across channels
# (-squared-channel-joining). Only active when WSC_JOINCHANNELS = True.
# Appropriate for Q/U which can rapidly oscillate positive/negative across channels
# (e.g. due to Faraday rotation); the sum-of-squares avoids cancellation.
WSC_SQUARECHANS  = True

# Controls QU deconvolution scaling when WSC_SPLITPOL is True. IV always keeps
# whatever auto-mask/auto-threshold was set above, regardless of this setting.
# Always a plain float multiplier on both auto-mask and auto-threshold; sign
# is a mode switch, magnitude is the scale factor:
#   0        : disabled — no changes to QU
#   negative : |value| x auto-mask/auto-threshold, nothing else changes
#   positive : value x auto-mask/auto-threshold, and if WSC_SQUARECHANS is
#              also True, local-rms is additionally turned on for QU
# NOTE: with WSC_SQUARECHANS True, a purely multiplicative QU auto-mask/
# threshold (i.e. a negative scale) tends to diverge hard -- likely the
# sum-of-squares peak-finding's non-Gaussian noise statistics. That's why
# local-rms defaults on for positive scales instead of being opt-in.
WSC_QU_AUTOMASK_SCALE = 1.0

# Stokes parameters to image. Standard options: 'I', 'IQUV', 'QU', etc.
WSC_POL = 'IQUV'

# Image QU and IV in separate wsclean calls:
#   QU: no fit-spectral-pol — MFS polynomial fitting is inappropriate for sources
#       with large RMs where Q/U vary rapidly and non-smoothly across the band
#   IV: no squared-channel-joining — not appropriate for Stokes V
# Recommended for targets with large rotation measures (>~100 rad/m^2).
WSC_SPLITPOL = True

# Join polarisations during peak finding (-join-polarizations).
# Useful for faint sources where combining Stokes improves peak detection.
# Avoid for dynamic-range-limited sources or when amplitude self-calibrating,
# as it can fold Q/U artefacts into the model.
WSC_JOINPOLARIZATIONS = False

# --- Masking: blind image (stage 0 — mask generation only) ---
WSC_AUTOMASK_BLIND      = 5.0   # Auto-mask threshold (sigma) for blind image
WSC_AUTOTHRESHOLD_BLIND = 1.0   # Auto-threshold (sigma) at which CLEAN stops for blind image
WSC_THRESHOLD_BLIND     = False  # Absolute flux threshold for blind image; False = use auto-threshold
WSC_LOCALRMS_BLIND      = 0.5    # Local RMS map for masking in the blind image; local-rms strength

# --- Masking: first self-calibration pass (two-stage workflow only) ---
# Restricts cleaning to high-significance pixels before any amplitude self-calibration.
# Has no effect in the standard single-round 2GC workflow.
WSC_SHALLOWMASK          = 30.0  # Auto-mask threshold (sigma) for the first deconvolution pass;
                                  # calibrators use 2× this value
WSC_SHALLOWMASK_LOCALRMS = 0.5   # Local RMS map for the first deconvolution pass; local-rms strength

# --- Model modification before predict (two-stage self-calibration only) ---
# If True, runs mod_model_selfcal.py after fix_nan_models.py before each predict step
# in the two-stage self-calibration workflow. Has no effect in a single-round workflow.
MOD_MODEL_SELFCAL = False

# --- Masking: intermediate image (two-stage workflow only) ---
# An intermediate image is made from CORRECTED_DATA after stage 1 phase self-calibration.
# This model is used for amplitude self-calibration in stage 2.
# Has no effect in the standard single-round 2GC workflow.
WSC_INTER_LOCALRMS      = 0.33   # Local RMS map for intermediate image masking; local-rms strength
WSC_INTER_AUTOMASK      = 3.0    # Auto-mask threshold (sigma) for intermediate image
WSC_INTER_AUTOTHRESHOLD = 1.0    # Auto-threshold (sigma) for intermediate image and shallow clean

# --- Masking: final image ---
# Standard 2GC: applied to both self-calibration imaging passes.
# Two-stage 2GC: applied to the final imaging pass only.
WSC_MASK          = False    # FITS mask to use; False = no mask (blind deconvolution)
                              # Can also be a list of per-field mask paths
WSC_THRESHOLD     = False    # Absolute flux threshold (Jy) at which to stop cleaning; False = use auto
WSC_AUTOMASK      = 3.0      # Auto-mask threshold (sigma) for main deconvolution
WSC_AUTOTHRESHOLD = 1.0      # Auto-threshold (sigma) at which CLEAN stops
WSC_LOCALRMS      = False    # Local RMS map for adaptive masking during main deconvolution:
                              # False disables it, True enables it at the default strength (0.5),
                              # a float enables it at that strength (see generate_syscall_wsclean())

# --- Minimum baseline UV cut ---
# If True, sets WSC_MINUVL by converting WSC_BASELINE_CUTLENGTH to wavelengths at the band centre.
# Useful for excluding short baselines that contribute large-scale confusion.
# WSC_BASELINE_CUTLENGTH can be:
#   '<value>m'  — a physical length in metres; converted to wavelengths per band (frequency-dependent)
#   '<value>'   — already in wavelengths; used directly with no frequency scaling
WSC_BASELINE_CUT = True
WSC_BASELINE_CUTLENGTH = '750m'  # Default: exclude baselines shorter than 750 m

if WSC_BASELINE_CUT:
    _cutlen = str(WSC_BASELINE_CUTLENGTH).strip()
    if _cutlen.endswith('m'):
        _baseline_m = float(_cutlen[:-1])
        _speed_of_light = 299792458.  # m / s
        _band_mid = {
            'UHF':  816.0e6,
            'L':   1284.0e6,
            'S0':  2187.0e6,
            'S1':  2406.0e6,
            'S2':  2625.0e6,
            'S3':  2843.0e6,
            'S4':  3062.0e6,
        }
        if BAND in _band_mid:
            _wavelength = _speed_of_light / _band_mid[BAND]
            WSC_MINUVL = '{:.1f}'.format(_baseline_m / _wavelength)
    else:
        # Already in wavelengths — use directly, no frequency scaling
        WSC_MINUVL = '{:.1f}'.format(float(_cutlen))


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

# --- Post-processing: beam homogenisation ---
WSC_HOMOGENIZEBEAM = False  # Homogenise across frequency channels only
WSC_HOMOGENIZETIME = False  # Homogenise across both frequency channels and time intervals

# [EXPERIMENTAL] Match the minimum-UV cut to the physical minimum baseline length per band,
# ensuring the large-scale sensitivity is consistent across frequency when homogenising the beam.
# Override WSC_MINUVL manually below the auto-calculated block if needed.
WSC_MATCHSCALES = WSC_HOMOGENIZEBEAM
if WSC_MATCHSCALES:
    speed_of_light = 299792458. # m / s
    min_baseline = 29. # m
    # Largest resolvable angular scale: top of band + 10 wavelengths padding
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


# ------------------------------------------------------------------------
#
# Calibrator imaging settings (setup_2GC_cals.py)
#

# Comma-separated list of field names to image and self-calibrate as calibrators.
# By default uses the polarization angle calibrator defined above.
CALS_TO_IMAGE = POLANG_NAME

# If True, combine all scans per field into a single MS before imaging.
# If False, split per field per scan (one MS per field+scan).
CAL_COMBINESCAN = False

# Blind mask threshold for calibrator imaging — by default 3x the target blind mask.
# Calibrators are typically brighter and more compact so a higher threshold
# avoids cleaning noise into the model before self-calibration.
WSC_CALS_BLINDMASK = WSC_AUTOMASK_BLIND * 3

# Self-calibration file: Default is full 2x2 complex solutions (assuming this is for a pol cal)
CAL_2GC_FILE = DATA+'/quartical/2GC_complex.yaml'

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
# RM Synthesis pipeline defaults (RMSYNTH_01 extraction + RM-Tools)
#

# --- Extraction (RMSYNTH_01_extract_fluxes.py) ---

RMSYN_OVERWRITE           = True    # Re-run fitting even when the output JSON already exists.
                                    # Set False to skip to plotting only.

RMSYN_FORCE_FIX_STOKES_V  = True   # Always anchor the Stokes V position to Stokes I regardless
                                    # of S/N. Recommended: leakage dominates V for most sources.

RMSYN_SPEC_PLOT_SNR_THRESH = 10     # Minimum MFS Stokes I S/N for a per-channel spectrum PNG.

RMSYN_SPEC_INDEX_SNR_THRESH = 25    # Minimum MFS Stokes I S/N to attempt the spectral index fit.
                                    # Below this alpha, chi2, ndof, chi2_red are stored as None.

RMSYN_SPEC_INDEX_MAD_CLIP  = 10     # Iterative MAD outlier rejection threshold (in sigma) applied
                                    # to channels during the spectral index power-law fit.

RMSYN_MAX_I_DRIFT_PIX      = False  # Maximum pixel drift allowed for a Stokes I fitted position
                                    # before the fit is re-run with the position fixed.
                                    # False = no drift check; float = threshold in pixels.

# --- rmsynth1d (RM-Tools) ---

RMSYN_FARADAY_RANGE        = 2500   # Half-range of the Faraday depth domain: |phi| <= this (rad/m²) (-l)
RMSYN_POLY_ORDER           = 4      # Polynomial order for the Stokes I spectral model fit (-o)
RMSYN_RMSF_SAMPLES         = 10     # Number of samples across the FWHM of the RMSF (-s)
RMSYN_SUPER_RESOLUTION     = True   # Enable super-resolution deconvolution (--super-resolution)

# --- rmclean1d (RM-Tools) ---

RMCLEAN_CUTOFF             = -9     # CLEAN stopping threshold; negative = multiples of the noise (-c)
RMCLEAN_WINDOW             = -4     # CLEAN window half-width; negative = multiples of RMSF half-width (-w)

# ------------------------------------------------------------------------
#
# Snapshot imaging defaults
#

SNAP_FIELDS = '' # Comma-separated list of field to run snapshot imaging on
SNAP_CHANNELSOUT = WSC_PCAL_CHANNELSOUT # Integer number of channels to perform snap-shot imaging -- by default the same as the WSC channels out
SNAP_INTBIN = 1 # Integer number of intervals to image together during snapshot imaging (default = 1 is per intergration imaging)
SNAP_INTEND = True # This will throw away the last interval if 'incomplete' compared to bin
    # E.G., total ints = 17, int bin = 4, will only image integrations 1 to 16 (4 total images) discarding 17
SNAP_POL = True # Bool for whether or not to do full polarisation snapshot images
SNAP_IMSIZE = 2560 # Image size of snapshot images (not the model image)
SNAP_MODELIDENTIFIER = 'pcalmask' # identifier for image name following oxkat/polkat conventions
SNAP_MODELMASK = '' # Point to mask for initial model creation (or just don't delete the pcalmask model)
SNAP_DECONV = False # Deconvolve during the snapshot imaging process? -- rarely necessary
SNAP_DECONVMASK = '' # mask to use when doing snapshot imaging and deconvolving
SNAP_RESTORE_MFS_ONLY = False # Only restore the MFS model image from the snapshot imaging, rather than also outputting the per-channel snapshot images

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
