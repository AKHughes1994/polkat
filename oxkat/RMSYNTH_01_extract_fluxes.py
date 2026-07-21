# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import os
import datetime
import subprocess
import sys
import json
import shutil
import re
import numpy as np
import os.path as o
import time
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush=True)


def get_source_from_ms(myms, source_name):
    """
    Query MS to find source by name and return position and intent.
    For per-scan sources (e.g., "J1744-5144_scan2"), strips the scan suffix
    to match the base field name in the MS.
    """
    msg(f"Querying MS for source: {source_name}")
    
    # Strip scan suffix if present (e.g., "J1744-5144_scan2" -> "J1744-5144")
    base_field_name = source_name
    scan_number = None
    scan_match = re.search(r'(.+)_scan(\d+)$', source_name)
    if scan_match:
        base_field_name = scan_match.group(1)
        scan_number = scan_match.group(2)
        msg(f"  Detected per-scan source: base={base_field_name}, scan={scan_number}")
    
    # Query FIELD table for source
    tb.open(myms + '/FIELD')
    field_names = tb.getcol('NAME')
    field_dirs = tb.getcol('PHASE_DIR')
    tb.close()
    
    # Find matching source (case-insensitive, using base field name)
    source_idx = None
    for i, fname in enumerate(field_names):
        if fname.lower() == base_field_name.lower():
            source_idx = i
            break

    if source_idx is None:
        msg(f"ERROR: Source '{base_field_name}' not found in MS FIELD table")
        msg(f"  Available fields: {', '.join(field_names)}")
        return {'name': source_name, 'ra_deg': None, 'dec_deg': None, 'intent': None, 'found': False}
    
    # Get position (convert radians to degrees)
    # PHASE_DIR shape is (2, 1, num_fields) where [0]=RA, [1]=DEC
    ra_rad = field_dirs[0][0][source_idx]
    dec_rad = field_dirs[1][0][source_idx]
    ra_deg = np.degrees(ra_rad)
    dec_deg = np.degrees(dec_rad)
    
    # Get intent from STATE table via main table
    tb.open(myms)
    field_col = tb.getcol('FIELD_ID')
    state_col = tb.getcol('STATE_ID')
    tb.close()
    
    # Find scans with this field
    field_mask = field_col == source_idx
    field_states = np.unique(state_col[field_mask])
    
    # Get intent from STATE table
    tb.open(myms + '/STATE')
    obs_modes = tb.getcol('OBS_MODE')
    tb.close()
    
    intent = obs_modes[field_states[0]] if len(field_states) > 0 else 'UNKNOWN'
    
    msg(f"Found source: {field_names[source_idx]}")
    msg(f"  Position: RA={ra_deg:.6f} deg, Dec={dec_deg:.6f} deg")
    msg(f"  Intent: {intent}")
    
    return {
        'name': field_names[source_idx],
        'ra_deg': ra_deg,
        'dec_deg': dec_deg,
        'intent': intent,
        'found': True
    }


def parse_casa_position(position_str):
    """
    Parse CASA format position string to decimal degrees.
    CASA format: RA is colon-separated (HH:MM:SS.SS), Dec is period-separated (±DD.MM.SS.SS)
    Example: '17:02:49.40,-48.47.23.40'
    
    Returns
    -------
    ra_deg, dec_deg : float, float
        Position in decimal degrees
    """
    try:
        # Split RA and Dec by comma
        ra_str, dec_str = position_str.split(',')
        
        # Parse RA (colon-separated: HH:MM:SS.SS)
        ra_parts = ra_str.strip().split(':')
        if len(ra_parts) != 3:
            msg(f"WARNING: Expected 3 parts in RA '{ra_str}', got {len(ra_parts)}")
            return None, None
        ra_hours = float(ra_parts[0])
        ra_mins = float(ra_parts[1])
        ra_secs = float(ra_parts[2])
        ra_deg = (ra_hours + ra_mins/60.0 + ra_secs/3600.0) * 15.0  # Convert hours to degrees
        
        # Parse Dec (period-separated: ±DD.MM.SS.SS)
        dec_str = dec_str.strip()
        dec_sign = -1 if dec_str.startswith('-') else 1
        dec_str = dec_str.lstrip('+-')  # Remove sign
        
        # Handle period-separated format
        dec_parts = dec_str.split('.')
        if len(dec_parts) < 3:
            msg(f"WARNING: Expected at least 3 parts in Dec '{dec_str}', got {len(dec_parts)}")
            return None, None
        dec_deg_part = float(dec_parts[0])
        dec_arcmin = float(dec_parts[1])
        # Reconstruct arcseconds from remaining parts
        dec_arcsec = float(dec_parts[2] + ('.' + dec_parts[3] if len(dec_parts) > 3 else ''))
        
        dec_deg = dec_sign * (abs(dec_deg_part) + dec_arcmin/60.0 + dec_arcsec/3600.0)
        
        return ra_deg, dec_deg
        
    except Exception as e:
        msg(f"ERROR parsing CASA position '{position_str}': {e}")
        return None, None


def read_position_file(position_file, source_name):
    """
    Read source position from text file for target sources.
    Format: FIELD_NAME  CASA_POSITION
    where CASA_POSITION is RA:colon-separated,Dec:period-separated
    """
    if not os.path.exists(position_file):
        msg(f"WARNING: Position file {position_file} not found")
        return None, None
    
    msg(f"Reading position from file: {position_file}")
    
    with open(position_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[0].lower() == source_name.lower():
                # Parse CASA format position
                position_str = parts[1]
                ra_deg, dec_deg = parse_casa_position(position_str)
                if ra_deg is not None:
                    msg(f"Found position in file: RA={ra_deg:.6f} deg, Dec={dec_deg:.6f} deg")
                    return ra_deg, dec_deg
    
    msg(f"WARNING: Source {source_name} not found in position file")
    return None, None


def read_time_info(time_file, field_name, scan=None):
    """
    Read time information from RESULTS directory.
    Parses format:
    # MS_NAME  FIELD_NAME  SCAN  START_MJD  END_MJD
    
    Parameters
    ----------
    time_file : str
        Path to time info file
    field_name : str
        Field name to match
    scan : str, optional
        Scan number to match (for per-scan calibrators)
        
    Returns
    -------
    dict : Time information with start_mjd, end_mjd, middle_mjd, duration_hours
    """
    if not os.path.exists(time_file):
        msg(f"WARNING: Time file {time_file} not found")
        return None
    
    msg(f"Reading time info from: {time_file}")
    
    time_data = None
    with open(time_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split()
            if len(parts) >= 5:
                ms_file = parts[0]
                field_name_file = parts[1]
                scan_file = parts[2]
                start_mjd = float(parts[3])
                end_mjd = float(parts[4])
                
                # Match by field name and optionally by scan
                if field_name_file.lower() == field_name.lower():
                    if scan is None or scan == scan_file:
                        time_data = {
                            'ms_name': ms_file,
                            'field_name': field_name_file,
                            'scan': scan_file,
                            'start_mjd': start_mjd,
                            'end_mjd': end_mjd,
                            'middle_mjd': (start_mjd + end_mjd) / 2.0,
                            'duration_hours': (end_mjd - start_mjd) * 24.0
                        }
                        msg(f"Matched field: {field_name_file}, Scan: {scan_file}")
                        break
    
    if time_data:
        msg(f"  Start MJD: {time_data['start_mjd']:.10f}")
        msg(f"  End MJD: {time_data['end_mjd']:.10f}")
        msg(f"  Middle MJD: {time_data['middle_mjd']:.10f}")
        msg(f"  Duration: {time_data['duration_hours']:.4f} hours")
    else:
        msg(f"WARNING: No time info found for source {source_name}")
    
    return time_data


def get_imfit_values(fname, image, xpix, ypix):

    '''
    Run imfit on a specific image (single component only)
    Inputs: 
        fname = string containing estimate file name
        image = string containing path to image to be fit
        xpix  = pixel coordinate (RA)
        ypix  = pixel coordinate (Dec)
    Returns:
        imfit dictionary
    '''

    # Get the beam parameters
    bmaj        = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin        = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa         = imhead(image, mode='get', hdkey = 'BPA')['value']

    # Single component fitting region
    x = xpix
    y = ypix
    r = 3 * bmaj
    src_region = f'circle[[{x}pix,{y}pix],{r}arcsec]'
   
    return imfit(image, estimates = fname, region = src_region)
    
    

def calculate_P0(flux_P, rms_Q, rms_U, rms_V, Aq = 0.8):
    '''
    Calculate the de-biased total polarized flux (Ptot = sqrt(Q^2 + U^2 + V^2))
    Always using total polarization (including V).
    '''    

    # Get the noise ratio coeffs, and calculate noise, from Hales 2012. https://arxiv.org/abs/1205.5310
    if rms_Q >= rms_U:
        A = Aq
        B = 1. - Aq 
    else: 
        B = Aq
        A = 1. - Aq 

    rms_P = (A * rms_Q ** 2 + B * rms_U ** 2) ** 0.5 

    # De-bias if SNR >= 3, following Vaillancourt 2006. https://arxiv.org/abs/astro-ph/0603110
    if flux_P / rms_P >= 3:
        flux_P0 = (flux_P ** 2 - rms_P ** 2) ** (0.5)
    else:
        flux_P0 = flux_P

    return flux_P0, rms_P

        

def return_max(im, region):
    '''
    Return the value that has the higher absolute magnitude
    Necessary for fluxes that are non-positive definite
    '''

    ims = imstat(im, region=region)
    fmax = ims['max'][0]
    fmin = ims['min'][0]

    if fmax > abs(fmin):
        return fmax
    else:
        return fmin


def get_imstat_values(image, xpix, ypix, manual_rms_region = False):
    '''
    Take in an image and a position, 
    return the max, max pixel location(s), and rms
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']    

    # Define the regions of interest (rms is ~100 beam area)
    r_in  = 2.0 * bmaj
    r_out = np.sqrt(100 * 0.25 * bmaj * bmin + r_in ** 2)
    src_region = f'circle[[{xpix}pix,{ypix}pix],3.0pix]'
    rms_region = f'annulus[[{xpix}pix,{ypix}pix],[{r_in}arcsec,{r_out}arcsec]]'
    if manual_rms_region:
        rms_region = manual_rms_region

    # Values of interest -- Source
    ims = imstat(image, region = src_region)
    flux = return_max(image, src_region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]

    # Extract RMS
    rms = imstat(image, region = rms_region)['rms'][0]

    return [flux, xpix, ypix, rms]
    
    
def check_position(fname, image, xpix, ypix, snr_thresh = 5.0, P_image = False, manual_rms_region = False, force_fix_to_I = False):

    '''
    Code to check whether there is sufficient flux at a position to allow 
    imfit to fit for position, or if said position should be frozen (single component)
    Inputs:
        fname = string containing name of output estimate file
        xpix = RA pixel position
        ypix = Dec pixel position
        image = name of the image to fit
        snr_thresh = 5.0
        P_image = Check if it's a P-image or not
        force_fix_to_I = If True, always fix position to Stokes I regardless of S/N
    Outputs:
        Nothing, but makes an estimate file
    '''
    
    # Get beam parameters
    bmaj = imhead(image, mode='get', hdkey='bmaj')['value']
    bmin = imhead(image, mode='get', hdkey='bmin')['value']
    bpa  = imhead(image, mode='get', hdkey='bpa')['value']

    region = f'circle[[{xpix}pix,{ypix}pix],{2 * bmaj}arcsec]'

    # Use same base name as fname for check_pos file to avoid conflicts
    check_file = fname.replace('estimate_', 'check_pos_')
    f = open(check_file, 'w')
    f.write(f'0.0,{xpix},{ypix},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, xyabp')
    f.close()

    # Get flux at test position
    test_flux = abs(imfit(image, region = region, estimates=check_file)['results']['component0']['peak']['value'])

    # Get the rms from the imstat rms calculation
    rms = get_imstat_values(image, xpix, ypix, manual_rms_region = manual_rms_region)[3]
    
    # If its a P-image don't use the image plane noise as the check criteria as it is (very) non-gaussian
    if P_image is True:
        image_Q = image.replace('-Plin-', '-Q-').replace('-Ptot-', '-Q-')
        image_U = image.replace('-Plin-', '-U-').replace('-Ptot-', '-U-')
        ims_Q = get_imstat_values(image_Q, xpix, ypix, manual_rms_region = manual_rms_region)
        ims_U = get_imstat_values(image_U, xpix, ypix, manual_rms_region = manual_rms_region)
        rms = np.amax((ims_Q[3],ims_U[3]))
    
    # Determine if we should fit or fix position
    if force_fix_to_I:
        fix_var = 'xyabp'
        msg(f'Forcing position fix to Stokes I (non-target calibrator): {image}')
    elif test_flux > snr_thresh * rms:
        fix_var = 'abp'
        msg(f'Free fitting position with S/N = {test_flux/rms:.2f}: {image}')
    else:
        fix_var = 'xyabp'
        msg(f'Fixing position due to S/N = {test_flux/rms:.2f}: {image}')

    # Make the estimate file
    make_estimate(fname, image, xpix, ypix, fix_var, manual_rms_region = manual_rms_region)


def make_estimate(fname, image, xpix, ypix, fix_var, manual_rms_region = False):
    '''
    Take in imstat values and create an CASA imfit estimate file (single component)
    Inputs:
        fname  = string containing name of estimate file
        image  = string containing name of image to fit
        xpix   = Right ascension (pixel) estimate of source
        ypix   = Declination (pixel) estimate of source
        fix    = parameters to fix, default is assume a point source (abp) or fixing position (xyabp)
        manual_rms_region = option to specify the rms region
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']
    
    # Make estimate file
    f = open(fname, 'w')
    ims = get_imstat_values(image, xpix, ypix, manual_rms_region = manual_rms_region)
    f.write(f'{ims[0]},{xpix},{ypix},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix_var}\n')
    f.close()            
    
    return 0    


def check_fit_offset_and_refix(imfit_result, image, ref_ra_pix, ref_dec_pix, bmaj,
                               src_name, stokes_label, manual_rms_region=False):
    '''
    After an imfit call, check whether the fitted position has drifted more than
    0.5 BMAJ (rough Euclidean pixel distance) from the reference position.
    If it has, re-fit with the position fixed to the reference and report the event.

    Parameters
    ----------
    imfit_result    : dict  - result from get_imfit_values
    image           : str   - image path (used for pixel-scale look-up via CDELT2)
    ref_ra_pix      : float - reference RA pixel (usually from Stokes I)
    ref_dec_pix     : float - reference Dec pixel (usually from Stokes I)
    bmaj            : float - beam major axis in arcsec (used as the distance threshold)
    src_name        : str   - source name (used for the temporary estimate filename)
    stokes_label    : str   - label for informative messages (e.g. "Q", "Ptot", "Plin")
    manual_rms_region : str/bool - passed through to make_estimate if a re-fit is needed

    Returns
    -------
    imfit_result : dict  - possibly updated if a re-fit was performed
    was_refitted : bool  - True if position was out of tolerance and re-fit was done
    '''
    component = 'component0'
    fitted_ra_pix  = imfit_result['results'][component]['pixelcoords'][0]
    fitted_dec_pix = imfit_result['results'][component]['pixelcoords'][1]

    # Euclidean pixel distance from the reference position
    pixel_dist = np.sqrt((fitted_ra_pix - ref_ra_pix)**2 +
                         (fitted_dec_pix - ref_dec_pix)**2)

    # Convert BMAJ (arcsec) to pixels so the threshold is physically meaningful
    # imhead returns CDELT2 in radians, so convert rad -> deg -> arcsec
    cell_rad    = abs(imhead(image, mode='get', hdkey='CDELT2')['value'])
    cell_arcsec = np.degrees(cell_rad) * 3600.0
    bmaj_pixels = bmaj / cell_arcsec
    offset_arcsec = pixel_dist * cell_arcsec

    # Use 1/3 of the beam (BMAJ) as the refit threshold instead of 1/2
    threshold_fraction = 1.0 / 3.0

    if pixel_dist > threshold_fraction * bmaj_pixels:
        msg(f'  WARNING: Stokes {stokes_label} fitted position offset = {offset_arcsec:.2f}" '
            f'({pixel_dist:.1f} pix) exceeds 1/3 BMAJ ({threshold_fraction * bmaj:.2f}"). '
            f'Re-fitting with position fixed to reference.')
        fixed_est = f'estimate_fixed_{stokes_label}_{src_name}.txt'
        make_estimate(fixed_est, image, ref_ra_pix, ref_dec_pix, 'xyabp',
                      manual_rms_region=manual_rms_region)
        imfit_result = get_imfit_values(fixed_est, image, ref_ra_pix, ref_dec_pix)
        return imfit_result, True
    else:
        msg(f'  Position check OK: Stokes {stokes_label} offset = {offset_arcsec:.2f}" '
            f'({pixel_dist:.1f} pix) < 1/3 BMAJ ({threshold_fraction * bmaj:.2f}")')
        return imfit_result, False


def extract_polarization_properties(src_name,
    src_im_identifier,
    src_im_suffix, 
    src_ra, 
    src_dec, 
    src_ulims,
    manual_rms_region, 
    image_directory,
    image_identifier,
    time_info = None,
    force_fix_to_I = False,
    use_plin = False):

    '''
    Fit the Stokes IQUV cube for single component in an image. This assumes
    that the Q, U, and V images follow the WSCLEAN naming
    convention.

    input parameters:
        src_name       = name of source 
        src_im_identifier  = identifier for images that will have flux extraction
        src_im_suffix  = image suffix, included to differentiate between standard WSCLEAN image products (image.fits) and homogenized products (image.homogenized.fits)
        src_ra   = Estimated right ascension of source (degrees)
        src_dec  = Estimated declination of  source (degrees)
        manual_rms_region = option to specify the rms region, otherwise use an annulus centered on the source(s) with a ~100xPSF area
        time_info = dictionary with timing information from ms_time_info.txt
        force_fix_to_I = If True, always fix all positions to Stokes I (for non-target calibrators)
        use_plin = If True, use Plin (linear: sqrt(Q^2+U^2)) instead of Ptot (sqrt(Q^2+U^2+V^2)).
                   Should be set True when polang_name is present in project_info.

    Output parametres:
        output_dictionary = Dictionary containing all MFS and per-channel information for the IQUV fluxes
    ''' 

    # Initialize dictionary to contain the output parameters
    output_dictionary = {'name' : src_name}
    
    # Add time information if available
    if time_info:
        output_dictionary['timing'] = time_info
    
    output_dictionary['MFS'] = {}
    output_dictionary['CHAN'] = {}

    # Determine which polarized intensity image type to use.
    # Plin = sqrt(Q^2+U^2) is preferred when a polarization angle calibrator is present
    # because V is not used in angle calibration and Plin has lower noise.
    # Ptot = sqrt(Q^2+U^2+V^2) is the more general total polarized intensity.
    pol_image_type    = 'Plin' if use_plin else 'Ptot'
    pol_image_exclude = '-Ptot-' if use_plin else '-Plin-'  # Filter out the non-selected type
    if use_plin:
        msg(f'Polarization image type: Plin (linear only) -- polang_name is set, V excluded from P estimate')
    else:
        msg(f'Polarization image type: Ptot (total, including V)')
    output_dictionary['pol_image_type'] = pol_image_type

    # Check for MFS images: prioritize image.fits, fall back to image.homogenized.fits
    msg("Checking for MFS images...")
    mfs_standard = glob.glob(f'{src_im_identifier}MFS*image.fits')
    mfs_homogenized = glob.glob(f'{src_im_identifier}MFS*image.homogenized.fits')
    
    if mfs_standard:
        mfs_im_suffix = 'image.fits'
        msg(f"Using standard MFS images: {mfs_im_suffix}")
    elif mfs_homogenized:
        mfs_im_suffix = 'image.homogenized.fits'
        msg(f"Standard MFS images not found, using homogenized MFS images: {mfs_im_suffix}")
    else:
        msg(f"WARNING: No MFS images found with pattern: {src_im_identifier}MFS*")
        mfs_im_suffix = None
    
    # Check for channelised images: prioritize image.fits, fall back to image.homogenized.fits
    msg("Checking for channelised images...")
    chan_standard = glob.glob(f'{src_im_identifier}*-I-image.fits')
    chan_standard = [im for im in chan_standard if '-MFS-' not in im]
    chan_homogenized = glob.glob(f'{src_im_identifier}*-I-image.homogenized.fits')
    chan_homogenized = [im for im in chan_homogenized if '-MFS-' not in im]
    
    if chan_standard:
        src_im_suffix = 'image.fits'
        msg(f"Using standard channelised images: {src_im_suffix}")
    elif chan_homogenized:
        src_im_suffix = 'image.homogenized.fits'
        msg(f"Standard channelised images not found, using homogenized channelised images: {src_im_suffix}")
    else:
        msg(f"ERROR: No channelised images found")
        raise FileNotFoundError(
            f"No channelised images found with pattern: "
            f"{src_im_identifier}*-I-image[.homogenized].fits")
    
    # Get unique prefixes - use channelised images to determine prefix
    # Extract prefix by removing the Stokes and channel parts
    # Pattern: prefix-####-STOKES-suffix.fits where #### is the channel number
    prefix_arr = chan_standard if chan_standard else chan_homogenized
    # Split from the right to handle source names with hyphens (e.g., J2011-0644)
    # Remove: -####-I-image.fits or -####-I-image.homogenized.fits
    example_img = prefix_arr[0]
    # Find the channel number pattern: -####- where # is a digit
    match = re.search(r'(.+)-\d{4}-[IQUV]-', example_img)
    if match:
        prefix = match.group(1)
    else:
        # Fallback: remove last 3 hyphen-separated parts
        prefix = example_img.rsplit('-', 3)[0]
    
    # Extract the MFS image parameters if they exist
    if mfs_im_suffix:
        msg(f'Fitting MFS image(s) for prefix: {prefix}')
        
        # Collect MFS images; order after sort will be I, P(lin/tot), Q, U, V.
        # Exclude whichever polarized intensity type is NOT being used (controlled by use_plin).
        MFS_images = glob.glob(f'{prefix}-MFS-*-{mfs_im_suffix}')
        MFS_images = sorted([im for im in MFS_images if pol_image_exclude not in im])
        msg(f'MFS images selected (using {pol_image_type}, excluding {pol_image_exclude.strip("-")}): {[os.path.basename(im) for im in MFS_images]}')

        # Output image name
        for MFS_image in MFS_images:
            msg(f'Fitting MFS image name(s): {MFS_image}')

        # Get generalized properties from the first image (i.e., Stokes I header)
        MFS_freq_GHz = imhead(MFS_images[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
        MFS_bmaj = imhead(MFS_images[0], mode='get', hdkey = 'bmaj')['value']
        MFS_bmin = imhead(MFS_images[0], mode='get', hdkey = 'bmin')['value']
        MFS_bpa   = imhead(MFS_images[0], mode='get', hdkey = 'bpa')['value']

        # Convert the RA/DEC guess to pixel coordinates (single component)
        region = 'circle[[{}deg,{}deg],1.0pix]'.format(src_ra, src_dec)
        pixel_imstat = imstat(MFS_images[0], region = region)
        src_ra_pix = pixel_imstat['maxpos'][0]
        src_dec_pix = pixel_imstat['maxpos'][1]
            
        # First fit Stokes I -- use check_position to determine if position should be fixed
        check_position(f'estimate_I_{src_name}.txt', MFS_images[0], 
                      src_ra_pix, src_dec_pix, P_image=False, 
                      manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
        MFS_I_imfit = get_imfit_values(f'estimate_I_{src_name}.txt', 
                                       MFS_images[0], src_ra_pix, src_dec_pix)

        # Get component (single component only)
        component = 'component0'           
        MFS_I_ra_pix = MFS_I_imfit['results'][component]['pixelcoords'][0]
        MFS_I_dec_pix = MFS_I_imfit['results'][component]['pixelcoords'][1]
       
        # Extract values for single component
        flux_I = MFS_I_imfit['results'][component]['peak']['value'] * 1e3
        err_I = MFS_I_imfit['results'][component]['peak']['error'] * 1e3
        RA_I   = MFS_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
        DEC_I = MFS_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
        RA_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][0]
        DEC_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][1]
        rms_I = get_imstat_values(MFS_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
        
        msg('')
        msg(f"--- MFS Stokes I Fitted Position (Reference) ---")
        msg(f"RA:  {RA_I:.6f} deg  ({RA_pix_I:.2f} pix)")
        msg(f"Dec: {DEC_I:.6f} deg ({DEC_pix_I:.2f} pix)")
        msg(f"Flux: {flux_I:.3f} +/- {err_I:.3f} mJy (S/N = {flux_I/rms_I:.1f})")
        msg('')
       
        # Initialize single-value MFS dictionary (not lists)
        output_dictionary['MFS'] = {}
        output_dictionary['MFS']['upperlimit'] = src_ulims[0]           
        output_dictionary['MFS']['freq_GHz'] = MFS_freq_GHz
        output_dictionary['MFS']['bmaj_asec'] = MFS_bmaj
        output_dictionary['MFS']['bmin_asec'] = MFS_bmin
        output_dictionary['MFS']['bpa_deg'] = MFS_bpa
        output_dictionary['MFS']['I_flux_mJy'] = flux_I
        output_dictionary['MFS']['I_err_mJy'] = err_I
        output_dictionary['MFS']['I_rms_mJy'] = rms_I
        output_dictionary['MFS']['I_RA_deg'] = RA_I
        output_dictionary['MFS']['I_DEC_deg'] = DEC_I
        
        # Next fit Q, U, V, and P(lin/tot) -- which P image is used is set by pol_image_type.
        # Fit P first (image 1 after I) - initially referenced to I
        check_position(f'estimate_P_{src_name}.txt', MFS_images[1], 
                      MFS_I_ra_pix, MFS_I_dec_pix, P_image=True, 
                      manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
        MFS_P_imfit  = get_imfit_values(f'estimate_P_{src_name}.txt', 
                                        MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
        # Check fit has not wandered more than 0.5 BMAJ from I reference
        MFS_P_imfit, _ = check_fit_offset_and_refix(
            MFS_P_imfit, MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix,
            MFS_bmaj, src_name, pol_image_type, manual_rms_region=manual_rms_region)
        
        # Get P flux and RMS to check significance
        flux_P = MFS_P_imfit['results'][component]['peak']['value'] * 1e3
        err_P = MFS_P_imfit['results'][component]['peak']['error'] * 1e3
        rms_P = get_imstat_values(MFS_images[1], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
        
        # Check if P is >10 sigma detection (only relevant if not forcing I position)
        P_snr = flux_P / rms_P
        if not force_fix_to_I and P_snr > 10.0:
            msg(f'{pol_image_type} has high SNR ({P_snr:.1f}), referencing Q/U against {pol_image_type} position')
            # Get P position to use as reference for Q/U
            MFS_P_ra_pix = MFS_P_imfit['results'][component]['pixelcoords'][0]
            MFS_P_dec_pix = MFS_P_imfit['results'][component]['pixelcoords'][1]
            ref_ra_pix = MFS_P_ra_pix
            ref_dec_pix = MFS_P_dec_pix
        else:
            if force_fix_to_I:
                msg(f'Referencing Q/U against Stokes I position (forced for non-target calibrator)')
            else:
                msg(f'{pol_image_type} has low SNR ({P_snr:.1f}), referencing Q/U against Stokes I position')
            ref_ra_pix = MFS_I_ra_pix
            ref_dec_pix = MFS_I_dec_pix
        
        # Fit Q (image 2) - referenced to chosen position
        check_position(f'estimate_Q_{src_name}.txt', MFS_images[2], 
                      ref_ra_pix, ref_dec_pix, P_image=False, 
                      manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
        MFS_Q_imfit  = get_imfit_values(f'estimate_Q_{src_name}.txt', 
                                        MFS_images[2], ref_ra_pix, ref_dec_pix)
        # Check fit has not wandered more than 0.5 BMAJ from reference
        MFS_Q_imfit, _ = check_fit_offset_and_refix(
            MFS_Q_imfit, MFS_images[2], ref_ra_pix, ref_dec_pix,
            MFS_bmaj, src_name, 'Q', manual_rms_region=manual_rms_region)

        # Fit U (image 3) - referenced to chosen position
        check_position(f'estimate_U_{src_name}.txt', MFS_images[3], 
                      ref_ra_pix, ref_dec_pix, P_image=False, 
                      manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
        MFS_U_imfit  = get_imfit_values(f'estimate_U_{src_name}.txt', 
                                        MFS_images[3], ref_ra_pix, ref_dec_pix)
        # Check fit has not wandered more than 0.5 BMAJ from reference
        MFS_U_imfit, _ = check_fit_offset_and_refix(
            MFS_U_imfit, MFS_images[3], ref_ra_pix, ref_dec_pix,
            MFS_bmaj, src_name, 'U', manual_rms_region=manual_rms_region)

        # Fit V (image 4) - ALWAYS referenced to Stokes I regardless of P/Q/U reference choice.
        # V is not used in linear polarization analysis and is always checked against I.
        check_position(f'estimate_V_{src_name}.txt', MFS_images[4], 
                      MFS_I_ra_pix, MFS_I_dec_pix, P_image=False, 
                      manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
        MFS_V_imfit  = get_imfit_values(f'estimate_V_{src_name}.txt', 
                                        MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix)
        # Check fit has not wandered more than 0.5 BMAJ from Stokes I reference
        MFS_V_imfit, _ = check_fit_offset_and_refix(
            MFS_V_imfit, MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix,
            MFS_bmaj, src_name, 'V', manual_rms_region=manual_rms_region)

        # Extract remaining flux values and errors
        flux_Q = MFS_Q_imfit['results'][component]['peak']['value'] * 1e3
        flux_U = MFS_U_imfit['results'][component]['peak']['value'] * 1e3
        flux_V = MFS_V_imfit['results'][component]['peak']['value'] * 1e3
        
        err_P = MFS_P_imfit['results'][component]['peak']['error'] * 1e3
        err_Q = MFS_Q_imfit['results'][component]['peak']['error'] * 1e3
        err_U = MFS_U_imfit['results'][component]['peak']['error'] * 1e3
        err_V = MFS_V_imfit['results'][component]['peak']['error'] * 1e3

        # Get remaining RMS values (P already computed above)
        rms_Q = get_imstat_values(MFS_images[2], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
        rms_U = get_imstat_values(MFS_images[3], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
        rms_V = get_imstat_values(MFS_images[4], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
        
        # Extract and display fitted positions for verification
        RA_Q = MFS_Q_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
        DEC_Q = MFS_Q_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
        RA_U = MFS_U_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
        DEC_U = MFS_U_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
        RA_V = MFS_V_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
        DEC_V = MFS_V_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
        
        msg('')
        msg(f"--- MFS Fitted Positions Summary ---")
        msg(f"Stokes Q: RA={RA_Q:.6f} deg, Dec={DEC_Q:.6f} deg (S/N={flux_Q/rms_Q:.1f})")
        msg(f"Stokes U: RA={RA_U:.6f} deg, Dec={DEC_U:.6f} deg (S/N={flux_U/rms_U:.1f})")
        msg(f"Stokes V: RA={RA_V:.6f} deg, Dec={DEC_V:.6f} deg (S/N={flux_V/rms_V:.1f})")
        
        # Calculate offsets from Stokes I to verify position fixing
        # Apply cos(dec) correction to RA term for proper angular separation
        cosdec = np.cos(np.radians(DEC_I))
        offset_Q = np.sqrt(((RA_Q - RA_I) * cosdec)**2 + (DEC_Q - DEC_I)**2) * 3600  # arcsec
        offset_U = np.sqrt(((RA_U - RA_I) * cosdec)**2 + (DEC_U - DEC_I)**2) * 3600  # arcsec
        offset_V = np.sqrt(((RA_V - RA_I) * cosdec)**2 + (DEC_V - DEC_I)**2) * 3600  # arcsec
        
        msg('')
        msg(f"Position offsets from Stokes I:")
        msg(f"Q offset: {offset_Q:.4f} arcsec")
        msg(f"U offset: {offset_U:.4f} arcsec")
        msg(f"V offset: {offset_V:.4f} arcsec")
        if force_fix_to_I:
            msg(f"(Should be ~0 since positions are fixed to Stokes I)")
        msg('')
        msg(f"--- MFS Flux Densities ---")
        msg(f"  I:    {flux_I:+10.4f} +/- {err_I:.4f} mJy  (S/N = {flux_I/rms_I:.1f})")
        msg(f"  {pol_image_type+':':5s} {flux_P:+10.4f} +/- {err_P:.4f} mJy  (S/N = {flux_P/rms_P:.1f})")
        msg(f"  Q:    {flux_Q:+10.4f} +/- {err_Q:.4f} mJy  (S/N = {flux_Q/rms_Q:.1f})")
        msg(f"  U:    {flux_U:+10.4f} +/- {err_U:.4f} mJy  (S/N = {flux_U/rms_U:.1f})")
        msg(f"  V:    {flux_V:+10.4f} +/- {err_V:.4f} mJy  (S/N = {flux_V/rms_V:.1f})")
        msg('')

        # Store in output dictionary
        output_dictionary['MFS'][f'{pol_image_type}_flux_mJy'] = flux_P
        output_dictionary['MFS']['Q_flux_mJy'] = flux_Q
        output_dictionary['MFS']['U_flux_mJy'] = flux_U
        output_dictionary['MFS']['V_flux_mJy'] = flux_V
        
        output_dictionary['MFS'][f'{pol_image_type}_err_mJy'] = err_P
        output_dictionary['MFS']['Q_err_mJy'] = err_Q
        output_dictionary['MFS']['U_err_mJy'] = err_U
        output_dictionary['MFS']['V_err_mJy'] = err_V

        output_dictionary['MFS'][f'{pol_image_type}_rms_mJy'] = rms_P
        output_dictionary['MFS']['Q_rms_mJy'] = rms_Q
        output_dictionary['MFS']['U_rms_mJy'] = rms_U
        output_dictionary['MFS']['V_rms_mJy'] = rms_V
    else:
        msg('Skipping MFS image processing (no MFS images found)')

    # Now process channelized data (simplified for single component, 1D arrays)
    msg('')
    msg(f"{'='*80}")
    msg('Extracting channelized data')
    msg(f"{'='*80}")
    msg('')
    
    # Find channelized images (not MFS) - always use full IQUV
    images_I = sorted([im for im in glob.glob(f'{prefix}-*-I-{src_im_suffix}') if '-MFS-' not in im])
    images_Q = sorted([im for im in glob.glob(f'{prefix}-*-Q-{src_im_suffix}') if '-MFS-' not in im])
    images_U = sorted([im for im in glob.glob(f'{prefix}-*-U-{src_im_suffix}') if '-MFS-' not in im])
    images_V = sorted([im for im in glob.glob(f'{prefix}-*-V-{src_im_suffix}') if '-MFS-' not in im])
    images_P = sorted([im for im in glob.glob(f'{prefix}-*-{pol_image_type}-{src_im_suffix}') if '-MFS-' not in im])
    msg(f'Using {pol_image_type} for channelised polarized intensity ({len(images_P)} images found)')

    msg(f'Found {len(images_I)} channelized Stokes I images')
    
    # Initialize 1D arrays for channelized data
    output_dictionary['CHAN'] = {}
    output_dictionary['CHAN']['freq_GHz'] = []
    output_dictionary['CHAN']['bmaj_asec'] = []
    output_dictionary['CHAN']['bmin_asec'] = []
    output_dictionary['CHAN']['bpa_deg'] = []
    output_dictionary['CHAN']['I_flux_mJy'] = []
    output_dictionary['CHAN']['I_err_mJy'] = []
    output_dictionary['CHAN']['I_rms_mJy'] = []
    output_dictionary['CHAN']['I_RA_deg'] = []
    output_dictionary['CHAN']['I_DEC_deg'] = []
    output_dictionary['CHAN']['Q_flux_mJy'] = []
    output_dictionary['CHAN']['U_flux_mJy'] = []
    output_dictionary['CHAN']['V_flux_mJy'] = []
    output_dictionary['CHAN'][f'{pol_image_type}_flux_mJy'] = []
    output_dictionary['CHAN']['Q_err_mJy'] = []
    output_dictionary['CHAN']['U_err_mJy'] = []
    output_dictionary['CHAN']['V_err_mJy'] = []
    output_dictionary['CHAN'][f'{pol_image_type}_err_mJy'] = []
    output_dictionary['CHAN']['Q_rms_mJy'] = []
    output_dictionary['CHAN']['U_rms_mJy'] = []
    output_dictionary['CHAN']['V_rms_mJy'] = []
    output_dictionary['CHAN'][f'{pol_image_type}_rms_mJy'] = []

    # Loop through channels
    for ch_idx in range(len(images_I)):
        
        try:
            # Get frequency and beam info
            freq_GHz = imhead(images_I[ch_idx], mode='get', hdkey='CRVAL3')['value'] / 1.0e9
            bmaj = imhead(images_I[ch_idx], mode='get', hdkey='bmaj')['value']
            bmin = imhead(images_I[ch_idx], mode='get', hdkey='bmin')['value']
            bpa = imhead(images_I[ch_idx], mode='get', hdkey='bpa')['value']

            msg('')
            msg(f'--- Channel {ch_idx} | {freq_GHz:.4f} GHz ---')

            # Check if beam is valid (not flagged channel)
            if bmaj == 0.0 or bmin == 0.0:
                msg(f'Channel {ch_idx} appears to be flagged (zero beam), skipping')
                continue
            
            # Check beam scaling - skip if beam is >10x larger than expected from frequency scaling
            # Expected beam scales as freq_MFS / freq_channel
            expected_bmaj = MFS_bmaj * (MFS_freq_GHz / freq_GHz)
            expected_bmin = MFS_bmin * (MFS_freq_GHz / freq_GHz)
            
            if bmaj > 10.0 * expected_bmaj or bmin > 10.0 * expected_bmin:
                msg(f'Channel {ch_idx} has anomalous beam '
                    f'(bmaj={bmaj:.2f}, bmin={bmin:.2f} vs '
                    f'expected {expected_bmaj:.2f}, {expected_bmin:.2f}), skipping')
                continue
        
            # Fit Stokes I
            check_position(f'estimate_ch_I_{src_name}.txt', images_I[ch_idx], 
                          MFS_I_ra_pix, MFS_I_dec_pix, P_image=False, 
                          manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
            ch_I_imfit = get_imfit_values(f'estimate_ch_I_{src_name}.txt', 
                                          images_I[ch_idx], MFS_I_ra_pix, MFS_I_dec_pix)
            
            flux_I = ch_I_imfit['results'][component]['peak']['value'] * 1e3
            err_I = ch_I_imfit['results'][component]['peak']['error'] * 1e3
            RA_I = ch_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_I = ch_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            RA_pix_ch = ch_I_imfit['results'][component]['pixelcoords'][0]
            DEC_pix_ch = ch_I_imfit['results'][component]['pixelcoords'][1]
            rms_I = get_imstat_values(images_I[ch_idx], RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3
            
            # Append to arrays
            output_dictionary['CHAN']['freq_GHz'].append(freq_GHz)
            output_dictionary['CHAN']['bmaj_asec'].append(bmaj)
            output_dictionary['CHAN']['bmin_asec'].append(bmin)
            output_dictionary['CHAN']['bpa_deg'].append(bpa)
            output_dictionary['CHAN']['I_flux_mJy'].append(flux_I)
            output_dictionary['CHAN']['I_err_mJy'].append(err_I)
            output_dictionary['CHAN']['I_rms_mJy'].append(rms_I)
            output_dictionary['CHAN']['I_RA_deg'].append(RA_I)
            output_dictionary['CHAN']['I_DEC_deg'].append(DEC_I)
            
            # Fit Q, U, V, Ptot for full polarization
            if ch_idx < len(images_Q):
                
                # Fit P(lin/tot) first to check significance (pol_image_type set at function top)
                check_position(f'estimate_ch_P_{src_name}.txt', images_P[ch_idx], 
                              RA_pix_ch, DEC_pix_ch, P_image=True, 
                              manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
                ch_P_imfit = get_imfit_values(f'estimate_ch_P_{src_name}.txt', 
                                              images_P[ch_idx], RA_pix_ch, DEC_pix_ch)
                # Check fit has not wandered more than 0.5 BMAJ from I reference
                ch_P_imfit, _ = check_fit_offset_and_refix(
                    ch_P_imfit, images_P[ch_idx], RA_pix_ch, DEC_pix_ch,
                    bmaj, src_name, pol_image_type, manual_rms_region=manual_rms_region)
                flux_P = ch_P_imfit['results'][component]['peak']['value'] * 1e3
                err_P = ch_P_imfit['results'][component]['peak']['error'] * 1e3
                rms_P = get_imstat_values(images_P[ch_idx], RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3
                
                # Check if P is >10 sigma - if so, use P position for Q/U reference (only if not forcing I position)
                P_snr_ch = flux_P / rms_P
                if not force_fix_to_I and P_snr_ch > 10.0:
                    msg(f'  Ch {ch_idx}: {pol_image_type} SNR={P_snr_ch:.1f} > 10, using {pol_image_type} as Q/U reference')
                    ref_ra_pix_ch = ch_P_imfit['results'][component]['pixelcoords'][0]
                    ref_dec_pix_ch = ch_P_imfit['results'][component]['pixelcoords'][1]
                else:
                    # Use I position for Q/U
                    ref_ra_pix_ch = RA_pix_ch
                    ref_dec_pix_ch = DEC_pix_ch
                
                # Fit Q - referenced to chosen position
                check_position(f'estimate_ch_Q_{src_name}.txt', images_Q[ch_idx], 
                              ref_ra_pix_ch, ref_dec_pix_ch, P_image=False, 
                              manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
                ch_Q_imfit = get_imfit_values(f'estimate_ch_Q_{src_name}.txt', 
                                              images_Q[ch_idx], ref_ra_pix_ch, ref_dec_pix_ch)
                # Check fit has not wandered more than 0.5 BMAJ from reference
                ch_Q_imfit, _ = check_fit_offset_and_refix(
                    ch_Q_imfit, images_Q[ch_idx], ref_ra_pix_ch, ref_dec_pix_ch,
                    bmaj, src_name, 'Q', manual_rms_region=manual_rms_region)
                flux_Q = ch_Q_imfit['results'][component]['peak']['value'] * 1e3
                err_Q = ch_Q_imfit['results'][component]['peak']['error'] * 1e3
                rms_Q = get_imstat_values(images_Q[ch_idx], RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3
                
                # Fit U - referenced to chosen position
                check_position(f'estimate_ch_U_{src_name}.txt', images_U[ch_idx], 
                              ref_ra_pix_ch, ref_dec_pix_ch, P_image=False, 
                              manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
                ch_U_imfit = get_imfit_values(f'estimate_ch_U_{src_name}.txt', 
                                              images_U[ch_idx], ref_ra_pix_ch, ref_dec_pix_ch)
                # Check fit has not wandered more than 0.5 BMAJ from reference
                ch_U_imfit, _ = check_fit_offset_and_refix(
                    ch_U_imfit, images_U[ch_idx], ref_ra_pix_ch, ref_dec_pix_ch,
                    bmaj, src_name, 'U', manual_rms_region=manual_rms_region)
                flux_U = ch_U_imfit['results'][component]['peak']['value'] * 1e3
                err_U = ch_U_imfit['results'][component]['peak']['error'] * 1e3
                rms_U = get_imstat_values(images_U[ch_idx], RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3
                
                # Fit V - ALWAYS referenced to Stokes I (independent of P/Q/U reference logic)
                check_position(f'estimate_ch_V_{src_name}.txt', images_V[ch_idx], 
                              RA_pix_ch, DEC_pix_ch, P_image=False, 
                              manual_rms_region=manual_rms_region, force_fix_to_I=force_fix_to_I)
                ch_V_imfit = get_imfit_values(f'estimate_ch_V_{src_name}.txt', 
                                              images_V[ch_idx], RA_pix_ch, DEC_pix_ch)
                # Check fit has not wandered more than 0.5 BMAJ from Stokes I reference
                ch_V_imfit, _ = check_fit_offset_and_refix(
                    ch_V_imfit, images_V[ch_idx], RA_pix_ch, DEC_pix_ch,
                    bmaj, src_name, 'V', manual_rms_region=manual_rms_region)
                flux_V = ch_V_imfit['results'][component]['peak']['value'] * 1e3
                err_V = ch_V_imfit['results'][component]['peak']['error'] * 1e3
                rms_V = get_imstat_values(images_V[ch_idx], RA_pix_ch, DEC_pix_ch, manual_rms_region)[3] * 1e3

                msg(f"  Flux densities:")
                msg(f"    I:    {flux_I:+10.4f} +/- {err_I:.4f} mJy  (S/N = {flux_I/rms_I:.1f})")
                msg(f"    {pol_image_type+':':5s} {flux_P:+10.4f} +/- {err_P:.4f} mJy  (S/N = {flux_P/rms_P:.1f})")
                msg(f"    Q:    {flux_Q:+10.4f} +/- {err_Q:.4f} mJy  (S/N = {flux_Q/rms_Q:.1f})")
                msg(f"    U:    {flux_U:+10.4f} +/- {err_U:.4f} mJy  (S/N = {flux_U/rms_U:.1f})")
                msg(f"    V:    {flux_V:+10.4f} +/- {err_V:.4f} mJy  (S/N = {flux_V/rms_V:.1f})")

                # Append to arrays
                output_dictionary['CHAN']['Q_flux_mJy'].append(flux_Q)
                output_dictionary['CHAN']['U_flux_mJy'].append(flux_U)
                output_dictionary['CHAN']['V_flux_mJy'].append(flux_V)
                output_dictionary['CHAN'][f'{pol_image_type}_flux_mJy'].append(flux_P)
                output_dictionary['CHAN']['Q_err_mJy'].append(err_Q)
                output_dictionary['CHAN']['U_err_mJy'].append(err_U)
                output_dictionary['CHAN']['V_err_mJy'].append(err_V)
                output_dictionary['CHAN'][f'{pol_image_type}_err_mJy'].append(err_P)
                output_dictionary['CHAN']['Q_rms_mJy'].append(rms_Q)
                output_dictionary['CHAN']['U_rms_mJy'].append(rms_U)
                output_dictionary['CHAN']['V_rms_mJy'].append(rms_V)
                output_dictionary['CHAN'][f'{pol_image_type}_rms_mJy'].append(rms_P)
        
        except Exception as e:
            msg(f'Fitting failed for channel {ch_idx}: {e}')
            msg('Channel is likely flagged, skipping')
            continue    
    
    # Create timestamp prefix from time_info if available
    timestamp_prefix = ""
    if time_info and 'start_mjd' in time_info:
        # Convert MJD to datetime and format as YYYYMMDDTHH
        # MJD epoch is November 17, 1858 at 00:00:00 UTC
        mjd_epoch = datetime.datetime(1858, 11, 17, 0, 0, 0)
        dt = mjd_epoch + datetime.timedelta(days=time_info['start_mjd'])
        timestamp_prefix = dt.strftime('%Y%m%dT%H_')
        msg(f'Using timestamp prefix: {timestamp_prefix}')
    
    # Save output with timestamp prefix
    output_file = cfg.RESULTS + f'/{timestamp_prefix}{src_name}_polarization.json'
    with open(output_file, 'w') as j:
        json.dump(output_dictionary, j, indent=4)
    
    msg(f'Saved output to: {output_file}')
    
    # Also save text file for IQUV with timestamp prefix
    txt_file = cfg.RESULTS + f'/{timestamp_prefix}{src_name}_iquv.txt'
    with open(txt_file, 'w') as f:
        f.write(f"# Source: {src_name}\n")
        f.write(f"# Position: RA={RA_I:.6f} deg, Dec={DEC_I:.6f} deg\n")
        if time_info:
            f.write(f"# MS: {time_info.get('ms_name', 'N/A')}\n")
            f.write(f"# Scan: {time_info.get('scan', 'N/A')}\n")
            f.write(f"# Start MJD: {time_info.get('start_mjd', 0):.10f}\n")
            f.write(f"# End MJD: {time_info.get('end_mjd', 0):.10f}\n")
            f.write(f"# Middle MJD: {time_info.get('middle_mjd', 0):.10f}\n")
            f.write(f"# Duration: {time_info.get('duration_hours', 0):.4f} hours\n")
        f.write(f"# Columns: Channel Freq[GHz] I[mJy] Q[mJy] U[mJy] V[mJy] "
                f"{pol_image_type}[mJy] rms_I rms_Q rms_U rms_V rms_{pol_image_type}\n")
        f.write("#\n")
        for i in range(len(output_dictionary['CHAN']['freq_GHz'])):
            f.write(f"{i:4d} "
                   f"{output_dictionary['CHAN']['freq_GHz'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['I_flux_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['Q_flux_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['U_flux_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['V_flux_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN'][f'{pol_image_type}_flux_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['I_rms_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['Q_rms_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['U_rms_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN']['V_rms_mJy'][i]:10.4f} "
                   f"{output_dictionary['CHAN'][f'{pol_image_type}_rms_mJy'][i]:10.4f}\n")
    msg(f'Saved text file to: {txt_file}')
    
    # Create additional raw data file for polang if applicable
    # This file contains: Freq(Hz), I, Q, U, dI, dQ, dU in Jy with no header
    return output_dictionary, timestamp_prefix


def create_polang_raw_file(project_info, src_name, output_dict, timestamp_prefix):
    """
    Create raw data file for polang analysis if target matches polang_name.
    Format: Freq(Hz) I Q U dI dQ dU (all fluxes in Jy, no header)
    """
    # Check if polang_name exists and is not empty
    if 'polang_name' not in project_info or not project_info['polang_name']:
        return
    
    polang_name = project_info['polang_name']
    
    # Check if target includes polang_name
    if polang_name.lower() not in src_name.lower():
        msg(f'Source {src_name} does not include polang_name {polang_name}, skipping raw data file')
        return
    
    msg(f'Creating raw data file for polang analysis (polang_name={polang_name})')
    
    # Create raw data file with timestamp prefix
    raw_file = cfg.RESULTS + f'/{timestamp_prefix}{src_name}_rmsynth.txt'
    
    with open(raw_file, 'w') as f:
        # No header - just data
        for i in range(len(output_dict['CHAN']['freq_GHz'])):
            freq_hz = output_dict['CHAN']['freq_GHz'][i] * 1e9  # Convert GHz to Hz
            I_jy = output_dict['CHAN']['I_flux_mJy'][i] / 1e3  # Convert mJy to Jy
            Q_jy = output_dict['CHAN']['Q_flux_mJy'][i] / 1e3
            U_jy = output_dict['CHAN']['U_flux_mJy'][i] / 1e3
            dI_jy = output_dict['CHAN']['I_rms_mJy'][i] / 1e3
            dQ_jy = output_dict['CHAN']['Q_rms_mJy'][i] / 1e3
            dU_jy = output_dict['CHAN']['U_rms_mJy'][i] / 1e3
            
            f.write(f"{freq_hz:.6e} {I_jy:.6e} {Q_jy:.6e} {U_jy:.6e} "
                   f"{dI_jy:.6e} {dQ_jy:.6e} {dU_jy:.6e}\n")
    
    msg(f'Saved polang raw data file to: {raw_file}')


def plot_stokes_spectrum(output_dictionary, src_name, timestamp_prefix):
    """
    Plot Stokes I spectrum and save to file.
    """
    # Plot Stokes I spectrum
    if len(output_dictionary['CHAN']['freq_GHz']) > 0:
        plot_file = cfg.VISPLOTS + f'/{timestamp_prefix}{src_name}_I_spectrum.png'
        
        # Calculate 5th and 95th percentile for y-axis limits
        flux_array = np.array(output_dictionary['CHAN']['I_flux_mJy'])
        flux_p05 = np.percentile(flux_array, 5)
        flux_p95 = np.percentile(flux_array, 95)
        flux_range = flux_p95 - flux_p05
        padding = 0.25 * flux_range
        
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        ax.errorbar(output_dictionary['CHAN']['freq_GHz'], 
                    output_dictionary['CHAN']['I_flux_mJy'],
                    yerr=output_dictionary['CHAN']['I_err_mJy'],
                    fmt='o', markersize=4, capsize=3, label='Stokes I')
        
        # Add MFS flux as horizontal line if available
        if 'MFS' in output_dictionary and 'I_flux_mJy' in output_dictionary['MFS']:
            mfs_flux = output_dictionary['MFS']['I_flux_mJy']
            freq_min = min(output_dictionary['CHAN']['freq_GHz'])
            freq_max = max(output_dictionary['CHAN']['freq_GHz'])
            ax.axhline(y=mfs_flux, color='k', linestyle='--', linewidth=2, 
                      label=f'MFS: {mfs_flux:.2f} mJy')
        
        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Flux Density (mJy)')
        ax.set_title(f'{src_name} - Stokes I Spectrum')
        ax.set_ylim(flux_p05 - padding, flux_p95 + padding)
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        plt.savefig(plot_file)
        plt.close()
        msg(f'Saved Stokes I spectrum plot to: {plot_file}')


def main():
    
    # Parse command line arguments
    if len(sys.argv) < 2:
        msg("ERROR: Usage: RMSYNTH_01_extract_fluxes.py <source_name>")
        msg("  source_name: Name of source to process")
        msg("    - For targets: source_name")
        msg("    - For per-scan calibrators: cal_name_scanN (e.g., J1939-6342_scan3)")
        sys.exit(1)
    
    source_name = sys.argv[1]
    
    # Parse source_name to detect per-scan calibrators
    scan_match = re.search(r'(.+)_scan(\d+)$', source_name)
    if scan_match:
        # Per-scan calibrator
        field_name = scan_match.group(1)
        scan_num = scan_match.group(2)
        msg(f"Per-scan calibrator detected: field={field_name}, scan={scan_num}")
    else:
        # Regular target or calibrator
        field_name = source_name
        scan_num = None
        msg(f"Regular source: {field_name}")
    
    # Construct image directory from source name
    image_dir = f'IMAGES/{source_name}'
    msg(f"Starting IQUV extraction for source: {source_name}")
    msg(f"Image directory: {image_dir}")
    
    # Get MS path from project info
    with open('project_info.json') as f:
        project_info = json.load(f)
    myms = project_info['working_ms']
    
    # Query MS for source information using field_name (not full source_name)
    source_info = get_source_from_ms(myms, field_name)
    
    if not source_info['found']:
        msg("ERROR: Cannot proceed without source information")
        sys.exit(1)
    
    # ========== CALIBRATOR POSITION FIXING CONTROL ==========
    # Set to True to enable forced position fixing for non-target calibrators (except excluded sources)
    # Set to False to use standard S/N-based position fitting for all sources
    USE_CALIBRATOR_POSITION_FIXING = True
    # =========================================================
    
    # Determine position based on intent
    is_target = 'TARGET' in source_info['intent'].upper()
    
    # Define list of excluded calibrator substrings (calibrators that should NOT have forced I-position fixing)
    excluded_calibrator_substrings = ['J1331', 'J1733', '1424']
    
    # Determine if we should force all positions to be fixed at Stokes I
    # This applies to: non-targets that DON'T contain excluded substrings (only if flag is enabled)
    force_fix_to_I = False
    if USE_CALIBRATOR_POSITION_FIXING and not is_target:
        # Check if field name contains any excluded substring
        contains_excluded = any(substring in field_name for substring in excluded_calibrator_substrings)
        if not contains_excluded:
            force_fix_to_I = True
            msg('')
            msg(f"{'='*80}")
            msg(f"NON-TARGET CALIBRATOR DETECTED: {field_name}")
            msg(f"All positions will be FIXED to Stokes I position (no free fitting)")
            msg(f"This behavior is disabled for sources containing: {', '.join(excluded_calibrator_substrings)}")
            msg(f"{'='*80}")
            msg('')
        else:
            msg('')
            msg(f"Non-target calibrator '{field_name}' contains excluded substring")
            msg(f"Using standard S/N-based position fitting logic")
            msg('')
    elif not USE_CALIBRATOR_POSITION_FIXING:
        msg('')
        msg(f"USE_CALIBRATOR_POSITION_FIXING is disabled - using standard S/N-based fitting for all sources")
        msg('')
    
    if is_target:
        msg("Source is a TARGET - reading position from file")
        position_file = cfg.DATA + '/positions/XRB_pos_list.txt'
        ra_deg, dec_deg = read_position_file(position_file, field_name)
        if ra_deg is None:
            msg("WARNING: Position not found in file, using MS position")
            ra_deg = source_info['ra_deg']
            dec_deg = source_info['dec_deg']
        else:
            # Calculate and report position difference
            ra_ms = source_info['ra_deg']
            dec_ms = source_info['dec_deg']
            delta_ra_arcsec = (ra_deg - ra_ms) * 3600.0
            delta_dec_arcsec = (dec_deg - dec_ms) * 3600.0
            msg(f"Position difference (file - MS):")
            msg(f"  Delta RA:  {delta_ra_arcsec:+.2f} arcsec")
            msg(f"  Delta Dec: {delta_dec_arcsec:+.2f} arcsec")
            source_info['ra_deg'] = ra_deg
            source_info['dec_deg'] = dec_deg
    else:
        msg("Source is NOT a target - using MS position")
        ra_deg = source_info['ra_deg']
        dec_deg = source_info['dec_deg']
    
    # Read time information
    # Match by field name and scan number (if per-scan)
    # Time file is named: {working_ms}_time_info.txt in RESULTS directory
    time_file = cfg.RESULTS + '/' + myms + '_time_info.txt'
    time_info = read_time_info(time_file, field_name, scan=scan_num)
    
    # Construct image identifier - check for non-zoom first, then fall back to zoom
    src_im_identifier_nozoom = f'{image_dir}/*{field_name}*.ms_pcalmask-'
    src_im_identifier_zoom = f'{image_dir}/*{field_name}*.ms_pcalmask_zoom-'
    
    # Check if non-zoom images exist
    test_images = glob.glob(f'{src_im_identifier_nozoom}0*-I-image.fits')
    if test_images:
        src_im_identifier = src_im_identifier_nozoom
        msg(f'Using non-zoom images: {src_im_identifier}')
    else:
        src_im_identifier = src_im_identifier_zoom
        msg(f'Non-zoom images not found, using zoom images: {src_im_identifier}')
    
    src_im_suffix = 'image.fits'  # Will be checked and adjusted inside the function
    
    # Call extraction function with time_info
    src_ulims = [False]  # Single component, not an upper limit
    manual_rms_region = False

    # Fix for J2011+7205 (known issue with position)
    #if 'J2011' in field_name:
    #    manual_rms_region = 'circle[[20:11:19.2635096009,-6.45.26.5310001268],60.0arcsec]'
    #    msg('Applying manual RMS region for J2011+7205 due to known position issues')

    if 'J1727' in field_name:
        manual_rms_region = "circle[[17:27:34.8888647559,-16.13.19.8402342899],75arcsec]"

    # Log the RMS region choice once here before extraction begins
    if manual_rms_region:
        msg(f'RMS region: MANUAL -- {manual_rms_region}')
    else:
        msg('RMS region: default annulus (~100 beam areas) centred on source')

    # Use Plin instead of Ptot when a polarization angle calibrator is present,
    # as Plin (sqrt(Q^2+U^2)) is more suitable for angle calibration than Ptot (sqrt(Q^2+U^2+V^2)).
    use_plin = bool(project_info.get('polang_name', ''))
    if use_plin:
        msg(f'polang_name is set ({project_info["polang_name"]}): will use Plin images for polarized intensity')
    else:
        msg('polang_name is not set: will use Ptot images for polarized intensity')

    msg(f'Starting property extraction for identifier: {src_im_identifier}')
    output_dict, timestamp_prefix = extract_polarization_properties(
        source_name,  # Use full source_name (includes _scanN for per-scan)
        src_im_identifier,
        src_im_suffix, 
        ra_deg,  # Single value, not array
        dec_deg,  # Single value, not array
        src_ulims,
        manual_rms_region, 
        image_dir,
        '',  # image_identifier not needed
        time_info = time_info,
        force_fix_to_I = force_fix_to_I,
        use_plin = use_plin)
    
    # Create polang raw data file if applicable
    create_polang_raw_file(project_info, source_name, output_dict, timestamp_prefix)
    
    # Create plots
    plot_stokes_spectrum(output_dict, source_name, timestamp_prefix)
    
    # Clean up temporary files created during extraction
    msg('Cleaning up temporary check_pos and estimate files')
    for temp_file in glob.glob(f'check_pos_*_{source_name}.txt'):
        os.remove(temp_file)
    for temp_file in glob.glob(f'estimate_*_{source_name}.txt'):
        os.remove(temp_file)



if __name__  == "__main__":
    main()
