import glob,os,datetime, subprocess, sys, json, time
import shutil
import numpy as np
import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush=True)

def calculate_P0(flux_P, rms_Q, rms_U, Aq = 0.8):
    '''
    Calculate the de-biased linearly polarized flux
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
    Necessary for fluxes that are non-positive definate
    
    input parameters:
        fmax = Positive maximum
        fmin = Negative maximum
    '''

    ims = imstat(im, region=region)
    fmax = ims['max'][0]
    fmin = ims['min'][0]

    if fmax > abs(fmin):
        return fmax
    else:
        return fmin

def get_imstat_values(image, x, y, n_beams = 4.0):
    '''
    Take in an image and a region, 
    return the max, rms, and max pixel location
    Inputs:
        image = string containing image name
        pos      = position for regions in CASA format
        n_beams = radius (inner radius) in number of beams for the extraction (rms) region 
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']
    r_in  = n_beams * bmaj
    r_out = np.sqrt(500 * 0.25 * bmaj * bmin + r_in ** 2)

    # Define regions
    region = f'circle[[{x}pix,{y}pix],{r_in}arcsec]'
    rms_region = f'annulus[[{x}pix,{y}pix],[{r_in}arcsec,{r_out}arcsec]]'

    # Values of interest
    ims = imstat(image, region=region)
    flux = return_max(image, region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]
    rms = imstat(image, region=rms_region)['rms'][0]

    return [flux, xpix, ypix, rms]


def make_estimate(fname, image, x, y, fix_var='abp'):
    '''
    Take in an array of imstat values from f(get_imstat_values)
    and return an CASA imfit estimate file name fname
    Inputs:
        fname  = string containing name of estimate file
        image  = string containing name of image to fit
        src_ra    = Right acension of peak of source(s)
        src_dec = Right acension of dec of source(s)
        fix = paramters to fix, default is assume a point source (abp) other revelant example is fixing position (xyabp)
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']
    
    # Get flux estimate
    region = f'circle[[{x}pix, {y}pix], {bmaj}arcsec]'
    flux_guess = return_max(image, region)
    
    # Make estimate file
    f = open(fname, 'w')
    f.write(f'{flux_guess},{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix_var}\n')
    f.close()            
    
    return 0    


def get_imfit_values(fname, image, x, y, n_beams = 5.0):
    '''
    Take in an image and an estimate and  return the CASA  imfit dictionary
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']

    # Run imfit
    r  = n_beams * bmaj
    region = f'circle[[{x}pix,{y}pix],{r}arcsec]'    
    imf = imfit(image, estimates= fname, region=region)

    return imf



def check_position(fname, image, x, y, snr_thresh = 5.0):

    '''
    Code to check whether there is sufficient flux at a position to allow 
    imfit to fit for position, or if said position should be frozen
    Inputs:
        fname = string containing name of ouput estimate file
        xpix = RA pixel position(s)
        ypix = Dec pixel position(s)
        image = name of the image to fit
        snr_thresh = 5.0
        P_image = Check if it's a P-image or
    Outputs:
        Nothing, but makes an estimate file
    '''
    
    # Get beam parameters
    bmaj = imhead(image, mode='get', hdkey='bmaj')['value']
    bmin = imhead(image, mode='get', hdkey='bmin')['value']
    bpa  = imhead(image, mode='get', hdkey='bpa')['value']


    region = f'circle[[{x}pix,{y}pix],{2 * bmaj}arcsec]'

    f = open('check_pos.txt', 'w')
    f.write(f'0.0,{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, xyabp')
    f.close()

    # Get flux at test position
    test_flux = abs(imfit(image, region = region, estimates='check_pos.txt')['results']['component0']['peak']['value'])

    # If its a P-image don't use the image plane noise as the check criteria as it is (very) non-gaussian
    rms = get_imstat_values(image, x, y)[3]
    if test_flux > snr_thresh * rms:
        fix_var = 'abp'
    else:
        fix_var = 'xyabp'

    # Make the estimate file
    make_estimate(fname, image, x, y, fix_var)


def find_mfs_images(pattern):
    '''
    Find MFS images with fallback to homogenized images and diagnostic_zoom
    Inputs:
        pattern = glob pattern to search for (should end with *image.fits)
    Returns:
        List of image paths, preferring *image.fits over *image.homogenized.fits or diagnostic_zoom
    '''
    
    # First try to find standard image.fits files
    images = glob.glob(pattern)
    
    if len(images) == 0:
        # If no standard images found, try homogenized versions
        homogenized_pattern = pattern.replace('image.fits', 'image.homogenized.fits')
        images = glob.glob(homogenized_pattern)
        
        if len(images) > 0:
            msg(f'Standard images not found, using homogenized images: {len(images)} files found')
        else:
            # Try diagnostic_zoom if diagnostic doesn't exist
            zoom_pattern = pattern.replace('diagnostic', 'diagnostic_zoom')
            images = glob.glob(zoom_pattern)
            
            if len(images) > 0:
                msg(f'Diagnostic images not found, using diagnostic_zoom images: {len(images)} files found')
            else:
                # Try homogenized version of diagnostic_zoom
                zoom_homogenized_pattern = zoom_pattern.replace('image.fits', 'image.homogenized.fits')
                images = glob.glob(zoom_homogenized_pattern)
                
                if len(images) > 0:
                    msg(f'Using diagnostic_zoom homogenized images: {len(images)} files found')
                else:
                    msg(f'No images found for pattern: {pattern}, homogenized variant, diagnostic_zoom, or diagnostic_zoom homogenized')
    else:
        msg(f'Using standard images: {len(images)} files found')
    
    return images

        
def get_polcal_polarization(pacal_name, pacal_pos, bpcal_sys):
    '''
    Get the MFS parameters for the polarization calibrator and ouput a dictionary
    Inputs: 
        pacal_pos     = position (in CASA format) of calibrator (e.g., 13:31:08.2881,+30.30.32.959)
        pacal_name = name of calibrator 
    Outputs:
        MFS_dict       = Dictionary containing MFS information for the fill stokes  
    '''

    # Initialize the output dictionary
    output_dictionary = {}
    output_dictionary['MFS'] = {}
    output_dictionary['CHAN'] = {}

    # Get the images -- these will be in the order of I, P, Q, U, V
    MFS_images = find_mfs_images(cfg.IMAGES + f'/*{pacal_name}*diagnostic-MFS-*-image.fits')
    MFS_images = sorted([im for im in MFS_images if '-Ptot-' not in im])

    # Extract the MFS image parameters
    msg(f'Fitting MFS image(s) for Pol. Ang. Cal.: {pacal_name}')
            
    # Get generalized properties from the first image (i.e., Stokes I header)
    freq_GHz = imhead(MFS_images[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
    date_obs = imhead(MFS_images[0], mode='get', hdkey = 'DATE-OBS').replace('/','-',2).replace('/','T')
    bmaj = imhead(MFS_images[0], mode='get', hdkey = 'bmaj')['value']
    bmin = imhead(MFS_images[0], mode='get', hdkey = 'bmin')['value']
    bpa   = imhead(MFS_images[0], mode='get', hdkey = 'bpa')['value']
    
    # Convert the target position to pixel coordinates
    region = f'circle[[{pacal_pos}],1.0pix]'
    pixel_imstat = imstat(MFS_images[0], region = region)
    src_ra_pix = pixel_imstat['maxpos'][0]
    src_dec_pix = pixel_imstat['maxpos'][1]
            
    # First fit Stokes I -- also extract pixel coordinates
    make_estimate('estimate_I.txt', MFS_images[0], src_ra_pix, src_dec_pix,  'abp')
    MFS_I_imfit = get_imfit_values('estimate_I.txt', MFS_images[0], src_ra_pix, src_dec_pix)['results']['component0']
    MFS_I_ra_pix = MFS_I_imfit['pixelcoords'][0] 
    MFS_I_dec_pix = MFS_I_imfit['pixelcoords'][1]

    # Fit Stokes Q, U (and 'P')
    check_position('estimate_P.txt', MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
    MFS_P_imfit  = get_imfit_values('estimate_P.txt', MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)['results']['component0']
    MFS_P_ra_pix = MFS_P_imfit['pixelcoords'][0] 
    MFS_P_dec_pix = MFS_P_imfit['pixelcoords'][1] 
            
    check_position('estimate_Q.txt', MFS_images[2], MFS_P_ra_pix, MFS_P_dec_pix)
    MFS_Q_imfit  = get_imfit_values('estimate_Q.txt', MFS_images[2], MFS_P_ra_pix, MFS_P_dec_pix)['results']['component0']
    MFS_Q_ra_pix = MFS_Q_imfit['pixelcoords'][0] 
    MFS_Q_dec_pix = MFS_Q_imfit['pixelcoords'][1] 
   
    check_position('estimate_U.txt', MFS_images[3], MFS_P_ra_pix, MFS_P_dec_pix)
    MFS_U_imfit  = get_imfit_values('estimate_U.txt', MFS_images[3], MFS_P_ra_pix, MFS_P_dec_pix)['results']['component0']
    MFS_U_ra_pix = MFS_U_imfit['pixelcoords'][0] 
    MFS_U_dec_pix = MFS_U_imfit['pixelcoords'][1]  

    # Fit Stokes V using the same procedure as Q and U
    check_position('estimate_V.txt', MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix)
    MFS_V_imfit  = get_imfit_values('estimate_V.txt', MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix)['results']['component0']
    MFS_V_ra_pix = MFS_V_imfit['pixelcoords'][0] 
    MFS_V_dec_pix = MFS_V_imfit['pixelcoords'][1]

    # For ease of readability define the desired quantities as variables
    flux_I = MFS_I_imfit['peak']['value'] * 1e3
    flux_P = MFS_P_imfit['peak']['value'] * 1e3
    flux_Q = MFS_Q_imfit['peak']['value'] * 1e3
    flux_U = MFS_U_imfit['peak']['value'] * 1e3
    flux_V = MFS_V_imfit['peak']['value'] * 1e3
    # flux_V = get_imstat_values(MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix, n_beams = 1.0)[0] * 1e3  # Old max/min approach

    err_I = MFS_I_imfit['peak']['error'] * 1e3
    err_P = MFS_P_imfit['peak']['error'] * 1e3
    err_Q = MFS_Q_imfit['peak']['error'] * 1e3
    err_U = MFS_U_imfit['peak']['error'] * 1e3
    err_V = MFS_V_imfit['peak']['error'] * 1e3
    
    rms_I = get_imstat_values(MFS_images[0], MFS_I_ra_pix, MFS_I_dec_pix)[3] * 1e3
    rms_Q = get_imstat_values(MFS_images[2], MFS_Q_ra_pix, MFS_Q_dec_pix)[3] * 1e3
    rms_U = get_imstat_values(MFS_images[3], MFS_U_ra_pix, MFS_U_dec_pix)[3] * 1e3
    rms_V = get_imstat_values(MFS_images[4], MFS_V_ra_pix, MFS_V_dec_pix)[3] * 1e3
    flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U)

    RA_I = MFS_I_imfit['shape']['direction']['m0']['value'] * 180 / np.pi
    RA_P = MFS_P_imfit['shape']['direction']['m0']['value'] * 180 / np.pi
    DEC_I = MFS_I_imfit['shape']['direction']['m1']['value'] * 180 / np.pi
    DEC_P = MFS_P_imfit['shape']['direction']['m1']['value'] * 180 / np.pi

    # Calculate extra parameters
    LP_frac     = flux_P0 / flux_I * 100.0
    LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2)

    LP_EVPA     = 0.5 * np.arctan2(flux_U, flux_Q) * 180.0 / np.pi 
    LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi

    # Put in MFS values
    output_dictionary['date_obs'] = date_obs

    output_dictionary['MFS']['freq_GHz'] = freq_GHz
    output_dictionary['MFS']['bmaj_asec'] = bmaj
    output_dictionary['MFS']['bmin_asec'] = bmin
    output_dictionary['MFS']['bpa_deg'] = bpa

    output_dictionary['MFS']['I_flux_mJy'] = flux_I
    output_dictionary['MFS']['P_flux_mJy'] = flux_P
    output_dictionary['MFS']['P0_flux_mJy'] = flux_P0
    output_dictionary['MFS']['Q_flux_mJy'] = flux_Q
    output_dictionary['MFS']['U_flux_mJy'] = flux_U
    output_dictionary['MFS']['V_flux_mJy'] = flux_V

    output_dictionary['MFS']['I_rms_mJy'] = rms_I
    output_dictionary['MFS']['P_rms_mJy'] = rms_P
    output_dictionary['MFS']['Q_rms_mJy'] = rms_Q
    output_dictionary['MFS']['U_rms_mJy'] = rms_U
    output_dictionary['MFS']['V_rms_mJy'] = rms_V

    output_dictionary['MFS']['I_err_mJy'] = err_I
    output_dictionary['MFS']['P_err_mJy'] = err_P
    output_dictionary['MFS']['Q_err_mJy'] = err_Q
    output_dictionary['MFS']['U_err_mJy'] = err_U
    output_dictionary['MFS']['V_err_mJy'] = err_V

    output_dictionary['MFS']['P_RA_deg'] = RA_P
    output_dictionary['MFS']['P_DEC_deg'] = DEC_P
    output_dictionary['MFS']['I_RA_deg'] = RA_I
    output_dictionary['MFS']['I_DEC_deg'] = DEC_I

    output_dictionary['MFS']['LP_frac'] = LP_frac
    output_dictionary['MFS']['LP_frac_err'] = LP_frac_err
    output_dictionary['MFS']['LP_EVPA'] = LP_EVPA
    output_dictionary['MFS']['LP_EVPA_err'] = LP_EVPA_err

    # Calculate systematic and also append 
    pacal_sys = flux_V / flux_I
    
    output_dictionary['date_obs'] = date_obs
    output_dictionary['BPCAL_RESIDUAL_POLFRAC'] = bpcal_sys
    output_dictionary['PACAL_RESIDUAL_VFRAC'] = pacal_sys

    # Extract channelized data
    CHAN_images = sorted(glob.glob(cfg.IMAGES + f'/*{pacal_name}*diagnostic-[!MFS]*-image.fits'))
    if len(CHAN_images) == 0:
        msg(f'Diagnostic channel images not found, trying diagnostic_zoom for {pacal_name}')
        CHAN_images = sorted(glob.glob(cfg.IMAGES + f'/*{pacal_name}*diagnostic_zoom-[!MFS]*-image.fits'))
    CHAN_images = sorted([im for im in CHAN_images if im.split('-')[-2] in ['I', 'Q', 'U', 'V', 'Plin']])
    CHAN_images_arr = np.array(CHAN_images).reshape(int(len(CHAN_images) / 5), 5) 

    for k, CHAN_images in enumerate(CHAN_images_arr[:]):

        msg(f'Fitting Channel for Pol. Ang. Cal. {pacal_name}: {CHAN_images[0]}')

        try:

            # Get generalized properties from the first image (i.e., Stokes I header)
            freq_GHz = imhead(CHAN_images[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
            date_obs = imhead(CHAN_images[0], mode='get', hdkey = 'DATE-OBS').replace('/','-',2).replace('/','T')
            bmaj = imhead(CHAN_images[0], mode='get', hdkey = 'bmaj')['value']
            bmin = imhead(CHAN_images[0], mode='get', hdkey = 'bmin')['value']
            bpa   = imhead(CHAN_images[0], mode='get', hdkey = 'bpa')['value']

            # Check if the beam is significantly elongated with respect to the expectation from MFS image (Channel is probably very flagged and should be omitted)
            beam_scaled = output_dictionary['MFS']['freq_GHz'] / freq_GHz   * output_dictionary['MFS']['bmaj_asec']
            if bmaj > 3.0 * beam_scaled :
                msg(f'Skipping Channel: BMAJ is different from expectation from MFS image (likely highly flagged; predicted = {beam_scaled:.2f}arcsec; observed = {bmaj:.2f}arcsec)')
                continue 
                        
            # First fit Stokes I -- also extract pixel coordinates
            check_position('estimate_I.txt', CHAN_images[0], MFS_I_ra_pix, MFS_I_dec_pix)
            CHAN_I_imfit = get_imfit_values('estimate_I.txt', CHAN_images[0], MFS_I_ra_pix, MFS_I_dec_pix)['results']['component0']
            CHAN_I_ra_pix = CHAN_I_imfit['pixelcoords'][0] 
            CHAN_I_dec_pix = CHAN_I_imfit['pixelcoords'][1]

            # Fit Stokes Q, U (and 'P')
            check_position('estimate_P.txt', CHAN_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
            CHAN_P_imfit  = get_imfit_values('estimate_P.txt', CHAN_images[1], MFS_I_ra_pix, MFS_I_dec_pix)['results']['component0']
            CHAN_P_ra_pix = CHAN_P_imfit['pixelcoords'][0] 
            CHAN_P_dec_pix = CHAN_P_imfit['pixelcoords'][1] 
                    
            check_position('estimate_Q.txt', CHAN_images[2], MFS_P_ra_pix, MFS_P_dec_pix)
            CHAN_Q_imfit  = get_imfit_values('estimate_Q.txt', CHAN_images[2], MFS_P_ra_pix, MFS_P_dec_pix)['results']['component0']
            CHAN_Q_ra_pix = CHAN_Q_imfit['pixelcoords'][0] 
            CHAN_Q_dec_pix = CHAN_Q_imfit['pixelcoords'][1] 
    
            check_position('estimate_U.txt', CHAN_images[3], MFS_P_ra_pix, MFS_P_dec_pix)
            CHAN_U_imfit  = get_imfit_values('estimate_U.txt', CHAN_images[3], MFS_P_ra_pix, MFS_P_dec_pix)['results']['component0']
            CHAN_U_ra_pix = CHAN_U_imfit['pixelcoords'][0] 
            CHAN_U_dec_pix = CHAN_U_imfit['pixelcoords'][1]  

            # Fit Stokes V using the same procedure as Q and U
            check_position('estimate_V.txt', CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix)
            CHAN_V_imfit  = get_imfit_values('estimate_V.txt', CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix)['results']['component0']
            CHAN_V_ra_pix = CHAN_V_imfit['pixelcoords'][0] 
            CHAN_V_dec_pix = CHAN_V_imfit['pixelcoords'][1]

            # For ease of readability define the desired quantities as variables
            flux_I = CHAN_I_imfit['peak']['value'] * 1e3
            flux_P = CHAN_P_imfit['peak']['value'] * 1e3
            flux_Q = CHAN_Q_imfit['peak']['value'] * 1e3
            flux_U = CHAN_U_imfit['peak']['value'] * 1e3
            flux_V = CHAN_V_imfit['peak']['value'] * 1e3
            # flux_V = get_imstat_values(CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix, n_beams = 1.0)[0] * 1e3  # Old max/min approach
    
            err_I = CHAN_I_imfit['peak']['error'] * 1e3
            err_P = CHAN_P_imfit['peak']['error'] * 1e3
            err_Q = CHAN_Q_imfit['peak']['error'] * 1e3
            err_U = CHAN_U_imfit['peak']['error'] * 1e3
            err_V = CHAN_V_imfit['peak']['error'] * 1e3
    
            rms_I = get_imstat_values(CHAN_images[0], CHAN_I_ra_pix, CHAN_I_dec_pix)[3] * 1e3
            rms_Q = get_imstat_values(CHAN_images[2], CHAN_Q_ra_pix, CHAN_Q_dec_pix)[3] * 1e3
            rms_U = get_imstat_values(CHAN_images[3], CHAN_U_ra_pix, CHAN_U_dec_pix)[3] * 1e3
            rms_V = get_imstat_values(CHAN_images[4], CHAN_V_ra_pix, CHAN_V_dec_pix)[3] * 1e3
            flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U)
        
            # Calculate extra parameters
            LP_frac     = flux_P0 / flux_I * 100.0
            LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2)

            LP_EVPA     = 0.5 * np.arctan2(flux_U, flux_Q) * 180.0 / np.pi 
            LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi

            # Try to append flux values to a list
            if 'freq_GHz' in output_dictionary['CHAN']:   
                output_dictionary['CHAN']['freq_GHz'].append(freq_GHz)
                output_dictionary['CHAN']['bmaj_asec'].append(bmaj)
                output_dictionary['CHAN']['bmin_asec'].append(bmin)
                output_dictionary['CHAN']['bpa_deg'].append(bpa)

                output_dictionary['CHAN']['I_flux_mJy'].append(flux_I)
                output_dictionary['CHAN']['P_flux_mJy'].append(flux_P)
                output_dictionary['CHAN']['P0_flux_mJy'].append(flux_P0)
                output_dictionary['CHAN']['Q_flux_mJy'].append(flux_Q)
                output_dictionary['CHAN']['U_flux_mJy'].append(flux_U)
                output_dictionary['CHAN']['V_flux_mJy'].append(flux_V)

                output_dictionary['CHAN']['I_rms_mJy'].append(rms_I)
                output_dictionary['CHAN']['P_rms_mJy'].append(rms_P)
                output_dictionary['CHAN']['Q_rms_mJy'].append(rms_Q)
                output_dictionary['CHAN']['U_rms_mJy'].append(rms_U)
                output_dictionary['CHAN']['V_rms_mJy'].append(rms_V)

                output_dictionary['CHAN']['I_err_mJy'].append(err_I)
                output_dictionary['CHAN']['P_err_mJy'].append(err_P)
                output_dictionary['CHAN']['Q_err_mJy'].append(err_Q)
                output_dictionary['CHAN']['U_err_mJy'].append(err_U)
                output_dictionary['CHAN']['V_err_mJy'].append(err_V)

                output_dictionary['CHAN']['LP_frac'].append(LP_frac)
                output_dictionary['CHAN']['LP_frac_err'].append(LP_frac_err)
                output_dictionary['CHAN']['LP_EVPA'].append(LP_EVPA)
                output_dictionary['CHAN']['LP_EVPA_err'].append(LP_EVPA_err)
                  
            # If it fails initialize the arrays       
            else: 
                msg('Initializing CHAN dictionary')
                output_dictionary['CHAN']['freq_GHz'] = [freq_GHz]
                output_dictionary['CHAN']['bmaj_asec'] = [bmaj]
                output_dictionary['CHAN']['bmin_asec'] = [bmin]
                output_dictionary['CHAN']['bpa_deg'] = [bpa]

                output_dictionary['CHAN']['I_flux_mJy'] = [flux_I]
                output_dictionary['CHAN']['P_flux_mJy'] = [flux_P]
                output_dictionary['CHAN']['P0_flux_mJy'] = [flux_P0]
                output_dictionary['CHAN']['Q_flux_mJy'] = [flux_Q]
                output_dictionary['CHAN']['U_flux_mJy'] = [flux_U]
                output_dictionary['CHAN']['V_flux_mJy'] = [flux_V]

                output_dictionary['CHAN']['I_rms_mJy'] = [rms_I]
                output_dictionary['CHAN']['P_rms_mJy'] = [rms_P]
                output_dictionary['CHAN']['Q_rms_mJy'] = [rms_Q]
                output_dictionary['CHAN']['U_rms_mJy'] = [rms_U]
                output_dictionary['CHAN']['V_rms_mJy'] = [rms_V]

                output_dictionary['CHAN']['I_err_mJy'] = [err_I]
                output_dictionary['CHAN']['P_err_mJy'] = [err_P]
                output_dictionary['CHAN']['Q_err_mJy'] = [err_Q]
                output_dictionary['CHAN']['U_err_mJy'] = [err_U]
                output_dictionary['CHAN']['V_err_mJy'] = [err_V]

                output_dictionary['CHAN']['LP_frac'] = [LP_frac]
                output_dictionary['CHAN']['LP_frac_err'] = [LP_frac_err]
                output_dictionary['CHAN']['LP_EVPA'] = [LP_EVPA]
                output_dictionary['CHAN']['LP_EVPA_err'] = [LP_EVPA_err]


        except:
            msg(f'Fitting Channel {k} failed')

    # Make the RM Synthesis files
    rmsynth_arr = np.array([np.array(output_dictionary['CHAN']['freq_GHz']) * 1e9, 
                                                    np.array(output_dictionary['CHAN']['I_flux_mJy']) / 1e3, 
                                                    np.array(output_dictionary['CHAN']['Q_flux_mJy']) / 1e3, 
                                                    np.array(output_dictionary['CHAN']['U_flux_mJy']) / 1e3,
                                                    np.array(output_dictionary['CHAN']['I_rms_mJy']) / 1e3,
                                                    np.array(output_dictionary['CHAN']['Q_rms_mJy']) / 1e3,
                                                    np.array(output_dictionary['CHAN']['U_rms_mJy']) / 1e3])

    # Save RM Synthesis file
    prefix = MFS_images[0].split('-MFS')[0].replace(cfg.IMAGES, cfg.RESULTS)
    np.savetxt(f'{prefix}_rmsynth.txt', rmsynth_arr.T)

    # Save output dictionary
    with open(f'{prefix}_polarization.json','w') as j:
        json.dump(output_dictionary, j, indent = 4)

    return pacal_sys


def get_primary_systematic(bpcal_name, bpcal_pos):
    '''
    Estimate the systematic error of the quality of the polarization leakage calibration
    using an the (unpolarized) primary calibrator (e.g., J1939-6342)

    Since the source is polarized any residual Q,U,V flux is likely due to unmodeled leakage near the phase center
    causing some Stokes I signal to appear polarized

    We can estimate the fractional systematic effect through the simple relation: 
    (Q + U + V) / (I) where, in the ideal case Q = U = V = 0.0!
    '''
    
    # Get the pixel coordinates
    MFS_I_images = find_mfs_images(cfg.IMAGES + f'/*{bpcal_name}*diagnostic-MFS-I-image.fits')
    MFS_I_image  = MFS_I_images[0]
    region = f'circle[[{bpcal_pos}],1.0pix]'
    pixel_imstat = imstat(MFS_I_image, region = region)
    x, y = pixel_imstat['maxpos'][0], pixel_imstat['maxpos'][1]

    # Fit Stokes I (source is very bright)
    estimate_I   = make_estimate('estimate_I.txt',  MFS_I_image, x, y, 'abp')
    MFS_I_imfit  = get_imfit_values('estimate_I.txt',  MFS_I_image,  x, y)
    flux_I = MFS_I_imfit['results']['component0']['peak']['value'] 

    # Get peak pixel Value from the total polarization image, theoretically should be zero, and thus, will quantify the systematic leakage
    MFS_P_images = find_mfs_images(cfg.IMAGES + f'/*{bpcal_name}*diagnostic-MFS-Ptot-image.fits')
    MFS_P_image = MFS_P_images[0] # Total polarization image 
    flux_P = get_imstat_values(MFS_P_image, x, y, n_beams = 1.0)[0]

    # Caculate systematic and return it
    bpcal_sys = flux_P / flux_I

    return bpcal_sys



def main():

    # Load in the project info dictionary
    with open('project_info.json') as f:
        project_info = json.load(f)

    # Get the coordinates for the primary from the project information file
    bpcal_name = project_info['primary_name']

    if bpcal_name == 'J1939-6342':
        bpcal_pos = '19:39:25.0264,-63.42.45.624'

    else: #J0408-6545
        bpcal_pos = '04:08:20.3782,-65.45.09.080'

    # Check if they made Stokes images of the calibrators
    bpcal_check_images = find_mfs_images(cfg.IMAGES + f'/*{bpcal_name}*diagnostic-MFS-I-image.fits')
    if len(bpcal_check_images) == 0:
        msg('ERROR: You do not have Stokes I,Q,U,V images of your calibrator(s)')
        sys.exit()

    # Load in the polarization angle calibrator name
    pacal_name  = cfg.POLANG_NAME
    pacal_pos   = cfg.POLANG_DIR

    # Check if the observations has a polarization angle calibrator or not
    pol_flag = False
    if pacal_name != '':
        pol_flag = True

    # Initialize dictionary to contain systematic terms
    # BPCAL_RESIDUAL_POLFRAC is the residual fraction (not %!) of Q^2 + U^2 + V^2 (which should be zero)
    # PACAL_RESIDUAL_VFRAC is the residual stokes V fraction (i.e., V / I, which should be zero)  
    systematics = {'BPCAL_RESIDUAL_POLFRAC': 0.0, 'PACAL_RESIDUAL_VFRAC': 0.0}

    # Primary Calibrator Systematic
    msg(f'Fitting Primary Systematics for BP Cal.: {bpcal_name}')
    bpcal_sys = get_primary_systematic(bpcal_name, bpcal_pos)
    systematics['BPCAL_RESIDUAL_POLFRAC'] = bpcal_sys

    # Check if pol. ang. systematics are necessary
    if pol_flag:
        
        # Polarization Calibrator Systematic
        msg(f'Fitting Polcal Systematics')

        fname = glob.glob(f'{cfg.RESULTS}/*{pacal_name}*_polarization.json')

        if fname == []:
            msg(f'Fitting Primary Systematics for PA Cal.: {pacal_name}')
            systematics['PACAL_RESIDUAL_VFRAC'] = get_polcal_polarization(pacal_name, pacal_pos, bpcal_sys)

        else:
            msg(f'Dictionary for {pacal_name} already exists; Opening; please delete if you want to re-calculate systematics')
            with open(fname[0], 'r') as j:
                pacal_dict = json.load(j)
            systematics['PACAL_RESIDUAL_VFRAC'] = pacal_dict['PACAL_RESIDUAL_VFRAC']

    # Append systematics to the polarization.json files
    for f in glob.glob(cfg.RESULTS + '/*polarization.json'):

        with open(f, 'r') as j:
            output_dictionary = json.load(j)

        output_dictionary['BPCAL_RESIDUAL_POLFRAC'] = systematics['BPCAL_RESIDUAL_POLFRAC'] 
        output_dictionary['PACAL_RESIDUAL_VFRAC'] = systematics['PACAL_RESIDUAL_VFRAC']

        with open(f, 'w') as j:
            json.dump(output_dictionary, j, indent = 4)


if __name__ == "__main__":
    main()


























