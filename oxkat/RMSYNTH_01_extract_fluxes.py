# VERY MESSY SCRIPT TO CALCULATE THE FLUXES OF POLARIZED COMPONENTS
# BUT IT WORKS

import glob,os,datetime, subprocess, sys, json
import shutil
import numpy as np
import os.path as o
import time
from scipy.spatial.distance import cdist
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)

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

def get_imfit_values(fname, image, xpix, ypix):
    '''
    Run imfit on a specific image
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
    pixel_asec  = imhead(image, mode='get', hdkey='cdelt2')['value'] * 3600 * 180.0 / np.pi # pixel size in asecs


    # Get the numnber of components
    n_comp = len(xpix)
   
    # For single components
    if n_comp == 1: 
        x = xpix[0]
        y = ypix[0]
        r = 3 * bmaj
        src_region = f'circle[[{x}pix,{y}pix],{r}arcsec]'
    
    # For multi-components (this won't work if there are components that are VERY far from eachother)
    else:
        # Define an array of coorindates [[x1,y1], [x2,y2], etc.])
        point_array =  np.array([xpix, ypix]).T

        # Get maximum distance between points
        max_dist = np.amax(cdist(point_array,point_array)) * pixel_asec
    
        # Take either twice the maximum distance or 3 times the bmaj axis as the bounding region radius
        x = xpix[0]
        y = ypix[0]
        r = np.amax((3 * bmaj, 2.0 * max_dist))
        src_region = f'circle[[{x}pix,{y}pix],{r}arcsec]'

    return imfit(image, estimates = fname, region = src_region)
        

def get_imstat_values(image, xpix, ypix):
    '''
    Take in an image and a position, 
    return the max, max pixel location(s), and rms
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']    

    # Define the regions of interest (rms is ~500 beam area)
    r_in  = 3.0 * bmaj
    r_out = np.sqrt(500 * 0.25 * bmaj * bmin + r_in ** 2)
    src_region = f'circle[[{xpix}pix,{ypix}pix],3.0pix]'
    rms_region = f'annulus[[{xpix}pix,{ypix}pix],[{r_in}arcsec,{r_out}arcsec]]'

    # Values of interest -- Source
    ims = imstat(image, region = src_region)
    flux = return_max(image, src_region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]

    # Extract RMS
    rms = imstat(image, region = rms_region)['rms'][0]

    return [flux, xpix, ypix, rms]


def make_estimate(fname, image, xpix, ypix, fix_var):
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
    
    # Make estimate file
    f = open(fname, 'w')
    for x, y, fix in zip(xpix, ypix, fix_var):
        ims = get_imstat_values(image, x, y)
        f.write(f'{ims[0]},{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix}\n')
    f.close()            
    
    return 0    

def get_IQUVP_names(im_I):
    '''
    Take in an image prefix and return the separated IQUVP images
    '''
    
    im_Q = im_I.replace('-I-', '-Q-')
    im_U = im_I.replace('-I-', '-U-')
    im_V = im_I.replace('-I-', '-V-')
    im_P = im_I.replace('-I-', '-P-')

    return [im_I, im_Q, im_U, im_V, im_P]

def check_position(fname, image, xpix, ypix, snr_thresh=5.0, P_image = False):
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

    fix_var = []

    # If a component is weak (i.e., < snr-thresh * sigma fix the position to the the reference otherwise fit for the position)
    for x, y in zip(xpix, ypix):
        ims   = get_imstat_values(image, x, y)

        # If its a P-image don't use the image plane noise as the check criteria as it is (very) non-gaussian
        if P_image:
            image_Q = image.replace('-P-', '-Q-')
            image_U = image.replace('-P-', '-U-')
            ims_Q = get_imstat_values(image_Q, x, y)
            ims_U = get_imstat_values(image_U, x, y)
            ims[3] = (ims_Q[3] + ims_U[3]) * 0.5

        if abs(ims[0])/ims[3] > snr_thresh:
            fix_var.append('abp')
        else:
            fix_var.append('xyabp')

    # Make the estimate file
    make_estimate(fname, image, xpix, ypix, fix_var)

def initialize_MFS_dict(i_image, imf_I, imf_Q, imf_U, imf_V, imf_P):
    '''
    Function that will take in the imfit dictionaries for all 4 (+ Lin. Pol.) Stokes Parameters, while
    returning a trimmed dictionary containing all of the parameters of interest
    '''

    IQUVP_names = get_IQUVP_names(i_image)

    # Initialize Dictionary
    MFS_dict = {}
    MFS_dict['beam'] = {}
    MFS_dict['beam']['bmaj_asec'] = imf_I['results']['component0']['beam']['beamarcsec']['major']['value']
    MFS_dict['beam']['bmin_asec'] = imf_I['results']['component0']['beam']['beamarcsec']['minor']['value']
    MFS_dict['beam']['bpa_deg']   = imf_I['results']['component0']['beam']['beamarcsec']['positionangle']['value']
   
    # Get the component keys
    comps = [key for key in imf_I['results'].keys() if 'component' in key]
    for comp in comps:
        MFS_dict[comp] = {}

        # Seperate out the fluxes and RMS for readability
        flux_I = imf_I['results'][comp]['peak']['value'] * 1e3
        flux_Q = imf_Q['results'][comp]['peak']['value'] * 1e3
        flux_U = imf_U['results'][comp]['peak']['value'] * 1e3
        flux_V = imf_V['results'][comp]['peak']['value'] * 1e3
        flux_P = imf_P['results'][comp]['peak']['value'] * 1e3

        rms_I = get_imstat_values(IQUVP_names[0], imf_I['results'][comp]['pixelcoords'][0], imf_I['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_Q = get_imstat_values(IQUVP_names[1], imf_Q['results'][comp]['pixelcoords'][0], imf_Q['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_U = get_imstat_values(IQUVP_names[2], imf_U['results'][comp]['pixelcoords'][0], imf_U['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_V = get_imstat_values(IQUVP_names[3], imf_V['results'][comp]['pixelcoords'][0], imf_V['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_P = 0.5 * (rms_Q + rms_U)        

        # Append fluxes and RMS to dictionary
        MFS_dict[comp]['I_flux_mJy'] = flux_I
        MFS_dict[comp]['Q_flux_mJy'] = flux_Q
        MFS_dict[comp]['U_flux_mJy'] = flux_U
        MFS_dict[comp]['V_flux_mJy'] = flux_V
        MFS_dict[comp]['P_flux_mJy'] = flux_P

        MFS_dict[comp]['I_rms_mJy'] = rms_I
        MFS_dict[comp]['Q_rms_mJy'] = rms_Q
        MFS_dict[comp]['U_rms_mJy'] = rms_U
        MFS_dict[comp]['V_rms_mJy'] = rms_V
        MFS_dict[comp]['P_rms_mJy'] = 0.5 * (rms_Q + rms_U)

        # Calculate the other Polarisation parameters -- Following George et al. 2012: arXiv:1106.5362
        if flux_P > 4.0 * rms_P:
            flux_P0 = np.sqrt(flux_P ** 2  - 2.3 * rms_P ** 2)
        else:
            flux_P0 = flux_P

        LP_frac     = flux_P0 / flux_I * 100.0
        LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2 )
        LP_EVPA     = np.arctan2(flux_U, flux_Q) * 180.0 / np.pi * 0.5
        LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi
    
        # Append additional polarisation parameters
        MFS_dict[comp]['P0_flux_mJy'] = flux_P0
        MFS_dict[comp]['LP_frac']     = LP_frac
        MFS_dict[comp]['LP_frac_err'] = LP_frac_err
        MFS_dict[comp]['LP_EVPA']     = LP_EVPA
        MFS_dict[comp]['LP_EVPA_err'] = LP_EVPA_err
    
        # Append the positions (Use Stokes I and Lin. Pol.)
        MFS_dict[comp]['I_RA_deg']  = imf_I['results'][comp]['shape']['direction']['m0']['value'] * 180 / np.pi
        MFS_dict[comp]['I_Dec_deg']  = imf_I['results'][comp]['shape']['direction']['m1']['value'] * 180 / np.pi
        MFS_dict[comp]['P_RA_deg'] = imf_P['results'][comp]['shape']['direction']['m0']['value'] * 180 / np.pi
        MFS_dict[comp]['P_Dec_deg'] = imf_P['results'][comp]['shape']['direction']['m1']['value'] * 180 / np.pi 

    return MFS_dict

def fit_channel(i_image, xpix_I, ypix_I, xpix_Q, ypix_Q, xpix_U, ypix_U, xpix_P, ypix_P, xpix_V, ypix_V):
    '''
    Extract the IQUP parameters from the working channelized image
    Inputs:
        i_chan_image = working image
        xpix/ypix  = for the MFS image of each Stokes Parameter
    Outputs:
        Array containing the relevant flux/Polarization parameters
    '''

    # Get all of the image names:
    IQUVP_names = get_IQUVP_names(i_image)

    # Fit Stokes I
    check_position('estimate_I.txt', IQUVP_names[0], xpix_I, ypix_I)
    imf_I  = get_imfit_values('estimate_I.txt', IQUVP_names[0], xpix_I, ypix_I)

    # Fit Stokes Q
    check_position('estimate_Q.txt', IQUVP_names[1], xpix_Q, ypix_Q)
    imf_Q  = get_imfit_values('estimate_Q.txt', IQUVP_names[1], xpix_Q, ypix_Q)

    # Fit Stokes U
    check_position('estimate_U.txt', IQUVP_names[2], xpix_U, ypix_U)
    imf_U  = get_imfit_values('estimate_U.txt', IQUVP_names[2], xpix_U, ypix_U)

    # Fit Stokes V
    check_position('estimate_V.txt', IQUVP_names[3], xpix_V, ypix_V)
    imf_V  = get_imfit_values('estimate_V.txt', IQUVP_names[3], xpix_V, ypix_V)

    # Fit Stokes P
    check_position('estimate_P.txt', IQUVP_names[4], xpix_P, ypix_P, P_image=True)
    imf_P  = get_imfit_values('estimate_P.txt', IQUVP_names[4], xpix_P, ypix_P)
       
    # Get the number of components
    comps = [key for key in imf_I['results'].keys() if 'component' in key]
    n_comps = len(comps)

    # Get the frequencies
    freq = np.ones(n_comps) * imhead(IQUVP_names[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9

    # Get the fluxes
    flux_I, flux_Q, flux_U, flux_V, flux_P, flux_P0 = np.ones((6, n_comps))
    rms_I,  rms_Q,  rms_U,  rms_V, rms_P, = np.ones((5, n_comps))

    # Get the component keys
    comps = [key for key in imf_I['results'].keys() if 'component' in key]
    for k, comp in enumerate(comps):

        # Fluxes
        flux_I[k] = imf_I['results'][comp]['peak']['value'] * 1e3
        flux_Q[k] = imf_Q['results'][comp]['peak']['value'] * 1e3
        flux_U[k] = imf_U['results'][comp]['peak']['value'] * 1e3
        flux_V[k] = imf_V['results'][comp]['peak']['value'] * 1e3
        flux_P[k] = imf_P['results'][comp]['peak']['value'] * 1e3

        # RMS
        rms_I[k] = get_imstat_values(IQUVP_names[0], imf_I['results'][comp]['pixelcoords'][0], imf_I['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_Q[k] = get_imstat_values(IQUVP_names[1], imf_Q['results'][comp]['pixelcoords'][0], imf_Q['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_U[k] = get_imstat_values(IQUVP_names[2], imf_U['results'][comp]['pixelcoords'][0], imf_U['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_V[k] = get_imstat_values(IQUVP_names[3], imf_V['results'][comp]['pixelcoords'][0], imf_V['results'][comp]['pixelcoords'][1])[3] * 1e3
        rms_P[k] = rms_Q[k] * 0.5 + rms_U[k] * 0.5

        # Bias corrected flux
        if flux_P[k] > 4.0 * rms_P[k]:
            flux_P0[k] = np.sqrt(flux_P[k] ** 2  - 2.3 * rms_P[k] ** 2)
        else:
            flux_P0[k] = flux_P[k]

    # Angle + Fraction calculations
    LP_frac     = flux_P0 / flux_I * 100.0
    LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2 )
    LP_EVPA     = np.arctan2(flux_U, flux_Q) * 180.0 / np.pi * 0.5
    LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi
    
    return [freq, flux_I, flux_Q, flux_U, flux_V, flux_P, flux_P0, LP_frac, LP_EVPA, rms_I, rms_Q, rms_U, rms_V, rms_P, LP_frac_err, LP_EVPA_err]


def make_rmsynth_files(src_header, chan_dict):
    '''
    Conver the channel dictionary to a text file with the following columns
    freq (Hz), I flux, Q flux, U flux, I rms, Q rms, U rms
    '''

    comps = [key for key in chan_dict.keys() if 'component' in key]
    for comp in comps:
        freq = np.array(chan_dict[comp]['freq_GHz']) * 1e9
        I   =   np.array(chan_dict[comp]['I_flux_mJy']) / 1e3
        Q  =  np.array(chan_dict[comp]['Q_flux_mJy']) / 1e3
        U  =  np.array(chan_dict[comp]['U_flux_mJy']) / 1e3
        dI  = np.array(chan_dict[comp]['I_rms_mJy']) / 1e3
        dQ = np.array(chan_dict[comp]['Q_rms_mJy']) / 1e3
        dU = np.array(chan_dict[comp]['U_rms_mJy']) / 1e3

        np.savetxt(cfg.RESULTS + f'/{src_header}_{comp}_rmsynth.txt', np.array([freq, I, Q, U, dI, dQ, dU]).T)    


def extract_polarization_properties(src_name,  src_im_identifier, src_ra, src_dec):

    '''
    Fit the Stokes IQUV cube for all components in an image. This assumes
    that the Q, U, and V images follow the WSCLEAN naming
    convention. This should work for slightly extended emission where
    you have two components that are yet to be separated by > 1 beam

    input parameters:
        src_name       = name of source 
        src_im_prefix = identifier for IQUV images that will have there fluxes extracted
        src_ra    = Estimated right acension of source(s) pixel units
        src_dec = Estimated declination of  source(s) pixel units

    Output parametres:
        flux_dict = Compehensive Dictionary containing all MFS and per-channel information for the IQUV fluxes
        rmsynth_arr = Array containing the necessary information to run RM synthesis on each image
    ''' 

    # Initilaize flux dictionary
    flux_dict = {'name': src_name}

    # Get the image prefixes with and without the full path
    src_im_prefix = glob.glob(cfg.IMAGES + f'/*{src_name}*{src_im_identifier}*')[0].split(f'{src_im_identifier}')[0] + src_im_identifier
    src_header = src_im_prefix.split(cfg.IMAGES + '/')[-1]
    msg(f'Fitting Image MFS images')
    
    # Begin by getting the properties of the MFS images
    IQUVP_names = get_IQUVP_names(glob.glob(f'{src_im_prefix}-MFS-I-*image.fits')[0])
    freqMFS_GHz = imhead(IQUVP_names[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
    date_obs  =  imhead(IQUVP_names[0], mode='get', hdkey = 'DATE-OBS').replace('/','-',2).replace('/','T')

    # The src_ra/decs can be a comma separed list for mutli_component fitting
    src_ra    = [ra for ra in src_ra.split(',')]
    src_dec = [dec for dec in src_dec.split(',')]

    # Convert from CASA hmsdms to pixel coordinates
    for k in range(len(src_ra)):
        region = f'circle[[{src_ra[k]},{src_dec[k]}],1.0pix]'
        ims = imstat(IQUVP_names[0], region = region)
        src_ra[k]  = ims['maxpos'][0]
        src_dec[k] = ims['maxpos'][1]

    # Make estimate files and fit for Stokes I
    make_estimate('estimate_I.txt', IQUVP_names[0], src_ra, src_dec,  ['abp'] * len(src_ra))
    imf_I  = get_imfit_values('estimate_I.txt', IQUVP_names[0], src_ra, src_dec)
    xpix_I = [imf_I['results'][key]['pixelcoords'][0] for key in imf_I['results'].keys() if 'component' in key]
    ypix_I = [imf_I['results'][key]['pixelcoords'][1] for key in imf_I['results'].keys() if 'component' in key]

    # Make estimate and fit for Lin. Pol. Intensity -- checking against position of stokes I components
    check_position('estimate_P.txt', IQUVP_names[4], xpix_I, ypix_I, P_image=True)
    imf_P  = get_imfit_values('estimate_P.txt', IQUVP_names[4], xpix_I, ypix_I)
    xpix_P = [imf_P['results'][key]['pixelcoords'][0] for key in imf_P['results'].keys() if 'component' in key]
    ypix_P = [imf_P['results'][key]['pixelcoords'][1] for key in imf_P['results'].keys() if 'component' in key]

    # Make estimate and fit for Q/U -- checking against position of Lin. Pol. components
    check_position('estimate_Q.txt', IQUVP_names[1], xpix_P, ypix_P)
    imf_Q  = get_imfit_values('estimate_Q.txt', IQUVP_names[1], xpix_P, ypix_P)
    xpix_Q = [imf_Q['results'][key]['pixelcoords'][0] for key in imf_Q['results'].keys() if 'component' in key]
    ypix_Q = [imf_Q['results'][key]['pixelcoords'][1] for key in imf_Q['results'].keys() if 'component' in key]
   
    check_position('estimate_U.txt', IQUVP_names[2], xpix_P, ypix_P)
    imf_U  = get_imfit_values('estimate_U.txt', IQUVP_names[2], xpix_P, ypix_P)
    xpix_U = [imf_U['results'][key]['pixelcoords'][0] for key in imf_U['results'].keys() if 'component' in key]
    ypix_U = [imf_U['results'][key]['pixelcoords'][1] for key in imf_U['results'].keys() if 'component' in key]

    # Make estimate and fit Stokes V  -- checking against position of stokes I components
    check_position('estimate_V.txt', IQUVP_names[3], xpix_I, ypix_I)
    imf_V  = get_imfit_values('estimate_V.txt', IQUVP_names[3], xpix_I, ypix_I)
    xpix_V = [imf_V['results'][key]['pixelcoords'][0] for key in imf_V['results'].keys() if 'component' in key]
    ypix_V = [imf_V['results'][key]['pixelcoords'][1] for key in imf_V['results'].keys() if 'component' in key]

    # Initialize the MFS parameters
    MFS_dict = initialize_MFS_dict(IQUVP_names[0], imf_I, imf_Q, imf_U, imf_V, imf_P)

    # Append the MFS imge parameters
    MFS_dict['freqMFS_GHz'] = freqMFS_GHz
    MFS_dict['date_isot']   = date_obs

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Now iterate through channelized images to get the frequecny evolution  #
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    
    # Initialize dictionary
    chan_dict = {}
    comps = [key for key in imf_I['results'].keys() if 'component' in key]
    for comp in comps:
        chan_dict[comp] = {'freq_GHz' : [], 'I_flux_mJy': [], 'Q_flux_mJy': [], 'U_flux_mJy': [], 'V_flux_mJy': [], 'P_flux_mJy': [], 'P0_flux_mJy': [], 'LP_frac': [], 'LP_EVPA' : [],  
                                            'I_rms_mJy': [], 'Q_rms_mJy': [], 'U_rms_mJy': [],  'V_rms_mJy': [], 'P_rms_mJy': [], 'LP_frac_err': [], 'LP_EVPA_err':[]}
    
    # Iterate through the channelized images    
    i_chan_images = sorted(glob.glob(f'{src_im_prefix}-00*-I-image.fits'))
    for i_chan_image in i_chan_images[:]:
        try: 
            msg(f'Fitting Image: {i_chan_image}')
            chan_data = fit_channel(i_chan_image, xpix_I, ypix_I, xpix_Q, ypix_Q, xpix_U, ypix_U, xpix_P, ypix_P, xpix_V, ypix_V)

            for k, comp in enumerate(comps):
                chan_dict[comp]['freq_GHz'].append(chan_data[0][k])
                chan_dict[comp]['I_flux_mJy'].append(chan_data[1][k])
                chan_dict[comp]['Q_flux_mJy'].append(chan_data[2][k])
                chan_dict[comp]['U_flux_mJy'].append(chan_data[3][k])
                chan_dict[comp]['V_flux_mJy'].append(chan_data[4][k])
                chan_dict[comp]['P_flux_mJy'].append(chan_data[5][k])
                chan_dict[comp]['P0_flux_mJy'].append(chan_data[6][k])
                chan_dict[comp]['LP_frac'].append(chan_data[7][k])
                chan_dict[comp]['LP_EVPA'].append(chan_data[8][k])
                chan_dict[comp]['I_rms_mJy'].append(chan_data[9][k])
                chan_dict[comp]['Q_rms_mJy'].append(chan_data[10][k])
                chan_dict[comp]['U_rms_mJy'].append(chan_data[11][k])
                chan_dict[comp]['V_rms_mJy'].append(chan_data[12][k])
                chan_dict[comp]['P_rms_mJy'].append(chan_data[13][k])
                chan_dict[comp]['LP_frac_err'].append(chan_data[14][k])
                chan_dict[comp]['LP_EVPA_err'].append(chan_data[15][k])
        except TypeError:
            msg('Fitting Failed: Channel is likely flagged')

    # Attach the sub_directories to the large directory + save as JSON
    flux_dict['MFS']  = MFS_dict
    flux_dict['CHAN'] = chan_dict

    with open(cfg.RESULTS + f'/{src_header}_polarization.json', 'w') as jfile:
        jfile.write(json.dumps(flux_dict, indent=4, sort_keys=True))

    # Make RM synthesis text files
    make_rmsynth_files(src_header, chan_dict)

    return 0

def main():
    
    # Load in the rmsynthesis data
    rmsynth_info = np.genfromtxt(cfg.DATA + '/rmsynth/rmsynth_info.txt', skip_header = 2, dtype=str)

    # If this is a 1-D array convert to two 2-D
    rmsynth_info = np.atleast_2d(rmsynth_info)

    # Break the columns into the relevant properties
    src_names, src_im_identifiers, src_ras, src_decs= rmsynth_info[:,0], rmsynth_info[:,1], rmsynth_info[:,2], rmsynth_info[:,3]

    for src_name, src_im_identifier, src_ra, src_dec in zip(src_names, src_im_identifiers, src_ras, src_decs):
        print('\n')
        msg(f'Starting Analysis of {src_name}')
        extract_polarization_properties(src_name, src_im_identifier, src_ra, src_dec)

if __name__  == "__main__":
    main()
