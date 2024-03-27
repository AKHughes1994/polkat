import glob,os,datetime, subprocess, sys, json, time
import shutil
import numpy as np
import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)

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

    # De-bias if SNR >= 4, following Vaillancourt 2006. https://arxiv.org/abs/astro-ph/0603110
    if flux_P / rms_P >= 4:
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

def get_IQUVP_names(im_I):
    '''
    Take in an image prefix and return the separated IQUVP images
    '''
    
    im_Q = im_I.replace('-I-', '-Q-')
    im_U = im_I.replace('-I-', '-U-')
    im_V = im_I.replace('-I-', '-V-')
    im_P = im_I.replace('-I-', '-Plin-') # This well measure linear polarization precisely (U^2 + Q^2) ** 0.5, with pol. ang.

    return [im_I, im_Q, im_U, im_V, im_P]


def get_imstat_values(image, pos, n_beams = 4.0):
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
    region = f'circle[[{pos}],{r_in}arcsec]'
    rms_region = f'annulus[[{pos}],[{r_in}arcsec,{r_out}arcsec]]'

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

def get_imfit_values(fname, image, x, y, n_beams = 4.0):
    '''
    Take in an image and an estimate and  return the CASA  imfit dictionary
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']

    # Fit region
    r  = n_beams * bmaj
    region = f'circle[[{x}pix,{y}pix],{r}arcsec]'
    
    imf = imfit(image, estimates= fname, region=region)

    return imf

def fit_channel(i_chan_image, pix_I, pix_Q, pix_U, pix_P, pos):
    '''
    Extract the IQUP parameters from the working channelized image
    Inputs:
        i_chan_image = working image
        xpix/ypix  = for the MFS image of each Stokes Parameter
    Outputs:
        Array containing the relevant flux/Polarization necessary for RM Synthesis
    '''

    # Get all of the image names:
    IQUVP_names = get_IQUVP_names(i_chan_image)

    # Get frequency
    freq_Hz  =  imhead(IQUVP_names[0], mode='get', hdkey = 'CRVAL3')['value']

    # Fit Stokes I
    check_position('estimate_I.txt', IQUVP_names[0], pix_I[0], pix_I[1], pos, snr_thresh=5.0)
    imf_I  = get_imfit_values('estimate_I.txt', IQUVP_names[0], pix_I[0], pix_I[1])

    # Fit Stokes Q
    check_position('estimate_Q.txt', IQUVP_names[1], pix_Q[0], pix_Q[1], pos, snr_thresh=5.0)
    imf_Q  = get_imfit_values('estimate_Q.txt', IQUVP_names[1], pix_Q[0], pix_Q[1])

    # Fit Stokes U
    check_position('estimate_U.txt', IQUVP_names[2], pix_U[0], pix_U[1], pos, snr_thresh=5.0)
    imf_U  = get_imfit_values('estimate_U.txt', IQUVP_names[2], pix_U[0], pix_U[1])

    # Fit Linear Intensity
    check_position('estimate_P.txt', IQUVP_names[4], pix_P[0], pix_P[1], pos, snr_thresh=5.0)
    imf_P  = get_imfit_values('estimate_P.txt', IQUVP_names[4], pix_P[0], pix_P[1])

    # Seperate fluxes
    flux_I = imf_I['results']['component0']['peak']['value'] 
    flux_Q = imf_Q['results']['component0']['peak']['value'] 
    flux_U = imf_U['results']['component0']['peak']['value']
    flux_P = imf_P['results']['component0']['peak']['value']  

    # RMS values
    rms_I =  get_imstat_values(IQUVP_names[0], pos)[3] 
    rms_Q =  get_imstat_values(IQUVP_names[1], pos)[3]
    rms_U =  get_imstat_values(IQUVP_names[2], pos)[3]
    flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U, Aq = 0.8)

    # Calculate the additional Linear Polarization parameters
    LP_frac     = flux_P0 / flux_I * 100.0
    LP_frac_err = LP_frac * ((rms_P / flux_P0) ** 2 + (rms_I / flux_I) ** 2) ** 0.5 

    return [freq_Hz, flux_I, flux_Q, flux_U, rms_I, rms_Q, rms_U, flux_P0, rms_P, LP_frac, LP_frac_err]   
        

def check_position(fname, image, x,y, pos, snr_thresh=5.0):
    '''
    Code to check whether there is sufficient flux at a position to allow 
    imfit to fit for position, or if said position should be frozen
    Inputs:
        fname = string containing name of ouput estimate file
        xpix = RA pixel position(s)
        ypix = Dec pixel position(s)
        image = name of the image to fit
        snr_thresh = 5.0
    Outputs:
        Nothing, but makes an estimate file
    '''

    # If a component is weak (i.e., < 5-sigma fix the position to the the reference otherwise fit for the position wonder
    ims = get_imstat_values(image, pos, n_beams=0.1)
    if abs(ims[0])/ims[3] > 5.0:
        fix_var = 'abp'
    else:
        fix_var = 'xyabp'

    # Make the estimate file
    make_estimate(fname, image, x, y, fix_var)

def calculate_sys_err(flux, rms, systematics, stokes = 'I'):
    '''
    Code to calculate a systematic error using the the residual polarized signals
    in the calibrators:
        For Bandpass Calibrator (BP) this is ANY polarization
        For Polartization Angle (PA) Calibrator this is ONLY Stokes V
    Inputs:
        flux = array containing [I,Q,U,V,P] fluxes
        rms  = array containing [I,Q,U,V,P] rms
        systematics = fractional systematic error from BP, PA (in that order)
        Stokes = Stokes parameter that will be solved (systematic is Stokes Depedant)
    Outputs:
        systematic error (single number or array)
    '''

    bpcal_sys, pacal_sys = systematics

    if type(flux) is list:
        flux = np.array(flux)
        rms  = np.array(rms)

    # Break apart the input arrays for readability
    I, Q, U, V, P      = flux
    dI, dQ, dU, dV, dP = rms
    
    if stokes == 'I':
        sys_err_sq = dI ** 2 +  (bpcal_sys * P) ** 2

    elif stokes == 'Q':
        sys_err_sq = dQ ** 2 +  (bpcal_sys * P) ** 2

    elif stokes == 'U':
        sys_err_sq = dU ** 2 +  (bpcal_sys * I) ** 2 + (pacal_sys * P) ** 2

    elif stokes == 'V':
        sys_err_sq = dV ** 2 +  (bpcal_sys * I) ** 2 + (pacal_sys * P) ** 2

    else:
        raise NameError("Sorry, must specify I, Q, U, or V")          
    
    return (sys_err_sq ** (0.5)).tolist()


def get_polcal_polarization(pacal_name, pacal_pos, bpcal_sys):
    '''
    Get the MFS parameters for the polarization calibrator and ouput a dictionary
    Inputs: 
        pacal_pos     = position (in CASA format) of calibrator (e.g., 13:31:08.2881,+30.30.32.959)
        pacal_name = name of calibrator 
    Outputs:
        MFS_dict       = Dictionary containing MFS information for the fill stokes  
    '''

    # Stokes I image
    image_I = glob.glob(cfg.IMAGES + f'/*{pacal_name}*postXf-MFS-I-image.fits')[0]

    # Now get the Full IQUV polarization properties for the polarization calibrator
    image_I, image_Q, image_U, image_V, image_P = get_IQUVP_names(image_I)

    # Get the basic properties from the image header
    freqMFS_GHz  =  imhead(image_I, mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
    date_obs  =  imhead(image_I, mode='get', hdkey = 'DATE-OBS').replace('/','-',2).replace('/','T')

    # Define the output file prefix
    file_prefix = image_I.split('-MFS')[0].replace('IMAGES', 'RESULTS')

    # Get imstat to get estimate of Stokes I position
    ims_I = get_imstat_values(image_I, pacal_pos)

    # Make estimate files for I, Q, U, and P
    estimate_I   = make_estimate('estimate_I.txt', image_I, ims_I[1], ims_I[2], 'abp')
    estimate_Q   = make_estimate('estimate_Q.txt', image_Q, ims_I[1], ims_I[2], 'abp')
    estimate_U   = make_estimate('estimate_U.txt', image_U, ims_I[1], ims_I[2], 'abp')
    estimate_P   = make_estimate('estimate_P.txt', image_P, ims_I[1], ims_I[2], 'abp')

    # Imfit estimate files for I, Q, U, and P
    imf_I = get_imfit_values('estimate_I.txt', image_I, ims_I[1], ims_I[2])
    imf_Q = get_imfit_values('estimate_Q.txt', image_Q, ims_I[1], ims_I[2])
    imf_U = get_imfit_values('estimate_U.txt', image_U, ims_I[1], ims_I[2])
    imf_P = get_imfit_values('estimate_P.txt', image_P, ims_I[1], ims_I[2])

    # Get imstat regions for I,Q,U,P,V
    ims_I = get_imstat_values(image_I, pacal_pos)
    ims_Q = get_imstat_values(image_Q, pacal_pos)
    ims_U = get_imstat_values(image_U, pacal_pos)
    ims_P = get_imstat_values(image_P, pacal_pos)
    ims_V = get_imstat_values(image_V, pacal_pos, n_beams = 1.0)

    # Get I, Q, U, P -- MFS pixel positions
    pix_I = imf_I['results']['component0']['pixelcoords']
    pix_Q = imf_Q['results']['component0']['pixelcoords']
    pix_U = imf_U['results']['component0']['pixelcoords']
    pix_P = imf_P['results']['component0']['pixelcoords']

    #  Separate out the fluxes and the RMS values
    flux_I  = imf_I['results']['component0']['peak']['value'] * 1e3
    rms_I  = ims_I[3] * 1e3

    flux_Q = imf_Q['results']['component0']['peak']['value'] * 1e3
    rms_Q = ims_Q[3] * 1e3

    flux_U = imf_U['results']['component0']['peak']['value'] * 1e3
    rms_U = ims_U[3] * 1e3

    flux_P = imf_P['results']['component0']['peak']['value'] * 1e3
    flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U, Aq = 0.8)

    flux_V = ims_V[0] * 1e3
    rms_V = ims_V[3] * 1e3

    # Since 3C286 is meant to be unpolarized in Stokes V, any residual Stokes V is from improper cross-hand phase calculate fractional systematic
    pacal_sys = abs(flux_V) / (flux_P ** 2 + flux_V **2) ** (0.5)
    systematics = [bpcal_sys, pacal_sys]
    flux = [flux_I, flux_Q, flux_U, flux_V, flux_P]
    rms  = [rms_I, rms_Q, rms_U, rms_V, rms_P]

    # Calculate error that include the systematic effects from insuffiicent leakage/cross-hand phase calibration
    sys_I = calculate_sys_err(flux, rms, systematics, stokes = 'I')
    sys_Q = calculate_sys_err(flux, rms, systematics, stokes = 'Q')
    sys_U = calculate_sys_err(flux, rms, systematics, stokes = 'U')
    sys_V = calculate_sys_err(flux, rms, systematics, stokes = 'V')
    sys_P = 0.5 * (sys_Q + sys_U)

    # Calculate other parameters
    LP_frac     = flux_P0 / flux_I * 100.0
    LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2)
    LP_frac_sys = LP_frac * np.sqrt( (sys_I / flux_I) ** 2 + (sys_P / flux_P0) ** 2)

    LP_EVPA     = 0.5 * np.arctan2(flux_U, flux_Q) * 180.0 / np.pi 
    LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi
    LP_EVPA_sys = LP_frac * np.sqrt( (sys_I / flux_I) ** 2 + (sys_P / flux_P0) ** 2)


    # Make MFS dictionary
    MFS_dict = {}
    MFS_dict['beam'] = {}
    MFS_dict['beam']['bmaj_asec'] = imf_I['results']['component0']['beam']['beamarcsec']['major']['value']
    MFS_dict['beam']['bmin_asec'] = imf_I['results']['component0']['beam']['beamarcsec']['minor']['value']
    MFS_dict['beam']['bpa_deg']   = imf_I['results']['component0']['beam']['beamarcsec']['positionangle']['value']
    MFS_dict['beam']['bpa_deg']   = imf_I['results']['component0']['beam']['beamarcsec']['positionangle']['value']
        
    # Fluxes
    MFS_dict['I_flux_mJy'] = flux_I
    MFS_dict['Q_flux_mJy'] = flux_Q
    MFS_dict['U_flux_mJy'] = flux_U
    MFS_dict['V_flux_mJy'] = flux_V
    MFS_dict['P_flux_mJy'] = flux_P
    MFS_dict['P0_flux_mJy'] = flux_P0

    # RMS
    MFS_dict['I_rms_mJy'] = rms_I
    MFS_dict['Q_rms_mJy'] = rms_Q
    MFS_dict['U_rms_mJy'] = rms_U
    MFS_dict['V_rms_mJy'] = rms_V
    MFS_dict['P_rms_mJy'] = rms_P

    # Systematic error
    MFS_dict['I_sys_mJy'] = sys_I
    MFS_dict['Q_sys_mJy'] = sys_Q
    MFS_dict['U_sys_mJy'] = sys_U
    MFS_dict['V_sys_mJy'] = sys_V
    MFS_dict['P_sys_mJy'] = sys_P

    # Position
    MFS_dict['I_RA_deg']   = imf_I['results']['component0']['shape']['direction']['m0']['value'] * 180 / np.pi
    MFS_dict['I_Dec_deg']  = imf_I['results']['component0']['shape']['direction']['m1']['value'] * 180 / np.pi
    MFS_dict['P_RA_deg']   = imf_P['results']['component0']['shape']['direction']['m0']['value'] * 180 / np.pi
    MFS_dict['P_Dec_deg']  = imf_P['results']['component0']['shape']['direction']['m1']['value'] * 180 / np.pi 

    # Polarization Parameters
    MFS_dict['LP_frac'] = LP_frac
    MFS_dict['LP_frac_err'] = LP_frac_err
    MFS_dict['LP_frac_sys'] = LP_frac_sys
    MFS_dict['LP_EVPA'] = LP_EVPA
    MFS_dict['LP_EVPA_err'] = LP_EVPA_err
    MFS_dict['LP_EVPA_sys'] = LP_EVPA_sys

    # Global Information
    MFS_dict['freqMFS_GHz'] = freqMFS_GHz
    MFS_dict['date_isot'] = date_obs

    # Initialize Channel dictionary
    chan_dict = {'freq_GHz' : [], 'I_flux_mJy': [], 'Q_flux_mJy': [], 'U_flux_mJy': [], 'P0_flux_mJy': [], 'LP_frac': [],  
                                         'I_rms_mJy': [], 'Q_rms_mJy': [], 'U_rms_mJy': [], 'P_rms_mJy': [], 'LP_frac_err': [],
                                         'I_sys_mJy': [], 'Q_sys_mJy': [], 'U_sys_mJy': [], 'P_sys_mJy': [], 'LP_frac_sys':[]}

    # Iterate through channelized images extracting fluxes to make RM Synthesis files
    rmsynth_arr = []
    for i_chan_image in sorted(glob.glob(cfg.IMAGES + f'/*{pacal_name}*postXf-[!MFS]*-I-image.fits')):
        try:                
            msg(f'Fitting Image: {i_chan_image.split("IMAGES/")[-1]}')
            # Measured Errors + values
            chan_data = fit_channel(i_chan_image, pix_I, pix_Q, pix_U, pix_P, pacal_pos)
            chan_dict['freq_GHz'].append(chan_data[0] / 1e9)
            chan_dict['I_flux_mJy'].append(chan_data[1] * 1e3)
            chan_dict['Q_flux_mJy'].append(chan_data[2] * 1e3)
            chan_dict['U_flux_mJy'].append(chan_data[3] * 1e3)
            chan_dict['I_rms_mJy'].append(chan_data[4] * 1e3)
            chan_dict['Q_rms_mJy'].append(chan_data[5] * 1e3)
            chan_dict['U_rms_mJy'].append(chan_data[6] * 1e3)

            chan_dict['P0_flux_mJy'].append(chan_data[7] * 1e3)
            chan_dict['P_rms_mJy'].append(chan_data[8] * 1e3)
            chan_dict['LP_frac'].append(chan_data[9])
            chan_dict['LP_frac_err'].append(chan_data[10])


            # Systematic Errors
            flux = [chan_data[1], chan_data[2], chan_data[3], chan_data[1], chan_data[7]] # Don't need stokes V -- dummy variable just make it I
            rms  = [chan_data[4], chan_data[5], chan_data[6], chan_data[1], chan_data[1]] # Don't need dV or dP -- dummy variables
            chan_sys_I  = calculate_sys_err(flux, rms, systematics, stokes = 'I')
            chan_sys_Q  = calculate_sys_err(flux, rms, systematics, stokes = 'Q')
            chan_sys_U  = calculate_sys_err(flux, rms, systematics, stokes = 'U')
            chan_sys_P  = (chan_sys_Q + chan_sys_U) * 0.5
            chan_sys_LP = chan_data[9] * (chan_sys_I ** 2 / chan_data[1] ** 2 + chan_sys_P ** 2 / chan_data[7] ** 2) ** 0.5

            chan_dict['I_sys_mJy'].append(chan_sys_I * 1e3)
            chan_dict['Q_sys_mJy'].append(chan_sys_Q * 1e3)
            chan_dict['U_sys_mJy'].append(chan_sys_U * 1e3)     
            chan_dict['P_sys_mJy'].append(chan_sys_P * 1e3)
            chan_dict['LP_frac_sys'].append(chan_sys_LP)
            
            rmsynth_arr.append(chan_data[:7])

        except:
            msg('Fitting Failed: Channel is likely flagged')
 
    # Save the RMSynth file
    np.savetxt(f'/{file_prefix}_rmsynth.txt', np.array(rmsynth_arr))


    # Make over arching dictionary
    pacal_dict = {}
    pacal_dict['date_isot'] = date_obs
    pacal_dict['Leakage_Fractional_Systematic']   = bpcal_sys
    pacal_dict['CrossHand_Fractional_Systematic'] = pacal_sys
    pacal_dict['MFS'] = MFS_dict
    pacal_dict['CHAN'] = chan_dict

    # Save dictionary
    with open(f'/{file_prefix}_polarization.json', 'w') as jfile:
        jfile.write(json.dumps(pacal_dict, indent=4, sort_keys=True))

    # Save RM files
    return pacal_sys

def update_rmsynth_file(fname, systematics):
    '''
    Load in RMsynth text files and update them, 
    such that we apply the systematic errors from the Leakage and Cross-hand
    calibration errors
    Input:
        fname    = Input string containing the path to the rmsynth file
        systematics = Array containing the bandpass and polarization angle calibrator fractional systematcis (in that order)
    Output:
        None -- But saves a new file
    '''

    # Load in the data
    data_arr = np.genfromtxt(fname)
    
    # Solve for the individual fluxes
    I, Q, U    = data_arr[:,1], data_arr[:,2], data_arr[:,3]
    dI, dQ, dU = data_arr[:,4], data_arr[:,5], data_arr[:,6]

    # Solve for polarized flux
    P = (Q ** 2 + U **2) ** (0.5)
    
    flux = [I, Q, U, I, P] # Don't need stokes V -- dummy variable just make it I
    rms = [dI, dQ, dU, dI, dI] # Don't need dV or dP -- dummy variables

    data_arr[:,4] = calculate_sys_err(flux, rms, systematics, stokes = 'I')
    data_arr[:,5] = calculate_sys_err(flux, rms, systematics, stokes = 'Q')
    data_arr[:,6] = calculate_sys_err(flux, rms, systematics, stokes = 'U')

    np.savetxt(fname.replace('.txt', '_sys.txt'), np.array(data_arr))

    return 0

def calc_sys_P(dQ, dU, Aq=0.8):
    if dQ >= dU:
        A = Aq
        B = 1. - Aq 
    else: 
        B = Aq
        A = 1. - Aq 

    dP = (A * rms_Q ** 2 + B * rms_U ** 2) ** 0.5     

    return dP


def update_src_pol_dict(fname, systematics, Aq = 0.8):

    with open(fname, 'r') as jfile:
        src_dict = json.load(jfile)

    
    bpcal_sys = systematics[0]
    pacal_sys = systematics[1] 

    src_dict['Leakage_Fractional_Systematic'] = bpcal_sys
    src_dict['CrossHand_Fractional_Systematic'] = pacal_sys

    #Append Systematic errors
    comps = [key for key in src_dict['MFS'].keys() if 'component' in key]
    for comp in comps:
        for subdict in ['MFS', 'CHAN']:
            flux = [src_dict[subdict][comp]['I_flux_mJy'], src_dict[subdict][comp]['Q_flux_mJy'], src_dict[subdict][comp]['U_flux_mJy'], src_dict[subdict][comp]['V_flux_mJy'], src_dict[subdict][comp]['P_flux_mJy']]
            rms = [src_dict[subdict][comp]['I_rms_mJy'], src_dict[subdict][comp]['Q_rms_mJy'], src_dict[subdict][comp]['U_rms_mJy'], src_dict[subdict][comp]['V_rms_mJy'], src_dict[subdict][comp]['P_rms_mJy']]

            # Append systematics to the dictionary
            src_dict[subdict][comp]['I_sys_mJy'] = calculate_sys_err(flux, rms, systematics, stokes = 'I')
            src_dict[subdict][comp]['Q_sys_mJy'] = calculate_sys_err(flux, rms, systematics, stokes = 'Q')
            src_dict[subdict][comp]['U_sys_mJy'] = calculate_sys_err(flux, rms, systematics, stokes = 'U')
            src_dict[subdict][comp]['V_sys_mJy'] = calculate_sys_err(flux, rms, systematics, stokes = 'V')

            # Define for variables for cleanliness of calculations
            flux_I = src_dict[subdict][comp]['I_flux_mJy'] 
            flux_Q = src_dict[subdict][comp]['Q_flux_mJy'] 
            flux_U = src_dict[subdict][comp]['U_flux_mJy'] 
            flux_V = src_dict[subdict][comp]['V_flux_mJy'] 
            flux_P = src_dict[subdict][comp]['P0_flux_mJy'] 
            LP_frac = src_dict[subdict][comp]['LP_frac'] 

            sys_I = src_dict[subdict][comp]['I_sys_mJy'] 
            sys_Q = src_dict[subdict][comp]['Q_sys_mJy'] 
            sys_U = src_dict[subdict][comp]['U_sys_mJy'] 
            sys_V = src_dict[subdict][comp]['V_sys_mJy'] 

            # Calculate pol. properties with systematic corrections
            
            # For CHAN apply list comprehension
            if type(sys_Q) is list:
                sys_P = [np.amax([dQ,dU,dV]) for dQ,dU,dV in zip(sys_Q, sys_U, sys_V)]
                LP_EVPA_err = [None] * len(sys_P)
                if pacal_sys != 0: # case with pol. cal.
                    sys_P       = [calc_sys_P(dQ, dU) for dQ,dU in zip(sys_Q, sys_U)]
                    LP_EVPA_err = [0.5 * np.sqrt(U ** 2 * dQ **2  + Q ** 2 * dU ** 2) / (U ** 2  + Q ** 2) * 180.0 / np.pi for U,dU,Q,dQ in zip(flux_U, sys_U, flux_Q, sys_Q)]
                LP_frac_sys = [LP * (dI ** 2 / I ** 2 + dP ** 2 / P ** 2) ** 0.5 for LP,I,dI,P,dP in zip(LP_frac, flux_I, sys_I, flux_P, sys_P)]
                

            # For MFS its just single numbers
            else:
                sys_P = np.amax([sys_Q, sys_U, sys_V])
                LP_EVPA_err = None
                if pacal_sys != 0:
                    sys_P = calc_sys_P(sys_Q, sys_U)
                    LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * sys_Q **2  + flux_Q ** 2 * sys_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi
                LP_frac_sys = LP_frac * np.sqrt( (sys_I / flux_I) ** 2 + (sys_P / flux_P) ** 2 )


            src_dict[subdict][comp]['P_sys_mJy'] = sys_P
            src_dict[subdict][comp]['LP_EVPA_sys'] = LP_EVPA_err
            src_dict[subdict][comp]['LP_frac_sys'] = LP_frac_sys

    # Resave the dictionary
    with open(fname, 'w') as jfile:
        jfile.write(json.dumps(src_dict, indent=4, sort_keys=True))

def get_primary_systematic(bpcal_name, bpcal_pos):
    '''
    Estimate the systematic error of the quality of the polarization leakage calibration
    using an the (unpolarized) primary calibrator (e.g., J1939-6342)

    Since the source is polarized any residual Q,U,V flux is likely due to unmodeled leakage near the phase center
    causing some Stokes I signal to appear polarized

    We can estimate the fractional systematic effect through the simple relation: 
    (Q + U + V) / (I) where, in the ideal case Q = U = V = 0.0!
    '''

    # Get imstat parameters from Stokes I image
    image_I = glob.glob(cfg.IMAGES + f'/*{bpcal_name}*postXf-MFS-I-image.fits')[0]
    ims_I       = get_imstat_values(image_I, bpcal_pos)
    
    # Fit Stokes I (source is very bright)
    estimate_I   = make_estimate('estimate_I.txt',  image_I,  ims_I[1], ims_I[2], 'abp')
    imf_I   = get_imfit_values('estimate_I.txt',  image_I,  ims_I[1], ims_I[2])
    flux_I = imf_I['results']['component0']['peak']['value'] 

    # Get peak pixel Value from the total polarization image, theoretically should be zero, and thus, will quantify the systematic leakage
    image_P = glob.glob(cfg.IMAGES + f'/*{bpcal_name}*postXf-MFS-Ptot-image.fits')[0] # Total polarization image 
    flux_P = get_imstat_values(image_P, bpcal_pos, n_beams=1.0)[0]

    # Caculate systematic and return it
    bpcal_sys = flux_P / (flux_I)

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

    # Load in the polarization angle calibrator name
    pacal_name  = cfg.POLANG_NAME
    pacal_pos   = cfg.POLANG_DIR

    # Check if the observations has a polarization angle calibrator or not
    pol_flag = False
    if pacal_name != '':
        pol_flag = True

    # Initialize array to contain systematic terms
    systematics = []

    # Primary Calibrator Systematic
    msg(f'Fitting Primary Systematics')
    bpcal_sys = get_primary_systematic(bpcal_name, bpcal_pos)
    systematics.append(bpcal_sys)

    if pol_flag:
        
        # Polarization Calibrator Systematic
        msg(f'Fitting Polcal MFS image')
        pacal_sys = get_polcal_polarization(pacal_name, pacal_pos, bpcal_sys)
        systematics.append(pacal_sys)


        # Update all rmsynth.txt files to include systematic errors
        for f in glob.glob(cfg.RESULTS + '/*_rmsynth.txt'):
            update_rmsynth_file(f, systematics)
    else:
        pacal_sys = 0 # No polarization calibrator
        systematics.append(pacal_sys)

    # Update the source dictionary with the systematic corrections
    rmsynth_info = np.genfromtxt(cfg.DATA + '/rmsynth/rmsynth_info.txt', skip_header = 2, dtype=str)

    # If this is a 1-D array convert to two 2-D
    rmsynth_info = np.atleast_2d(rmsynth_info)

    # Break the columns into the relevant properties
    src_names = rmsynth_info[:,0]
    src_im_identifiers = rmsynth_info[:,1]

    # Add systematic errors to src dictionary
    for src_name, src_im_identifier in zip(src_names, src_im_identifiers):
        fname = glob.glob(cfg.RESULTS + f'/*{src_name}*{src_im_identifier}*_polarization.json')[0]
        update_src_pol_dict(fname, systematics)


if __name__ == "__main__":
    main()


























