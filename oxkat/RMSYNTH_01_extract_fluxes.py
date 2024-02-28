# VERY MESSY SCRIPT TO CALCULATE THE FLUXES OF POLARIZED COMPONENTS

import glob,os,datetime, subprocess, sys, json
import shutil
import numpy as np
import os.path as o
import time
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

def get_imstat_values(image, region):
    '''
    Take in an image and a region, 
    return the max, rms, and max pixel location
    '''

    # Values of interest
    ims = imstat(image, region=region)
    flux = return_max(image, region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']

    # Get rms region around annulus centered on source (area will be ~500 beams)
    r_in  = 3.0 * bmaj
    r_out = np.sqrt(500 * 0.25 * bmaj * bmin + r_in ** 2)
    rms_region = f'annulus[[{xpix}pix,{ypix}pix],[{r_in}arcsec,{r_out}arcsec]]'
    rms = imstat(image, region=rms_region)['rms'][0]

    return [flux, xpix, ypix, rms]

def get_imfit_values(image, fname, src_ra, src_dec):
    '''
    Take in an image and an estimate and  return the CASA  imfit dictionary
    '''

    # Get the beam parameters
    bmaj = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey = 'BPA')['value']

    # Fit region
    r  = 5.0 * bmaj
    src_ra = np.mean(src_ra)
    src_dec = np.mean(src_dec)    

    region = f'circle[[{src_ra}pix,{src_dec}pix],{r}arcsec]'
    
    imf = imfit(image, estimates= fname, region=region)

    return imf


def make_estimate(fname, image, src_ra, src_dec, fix_var):
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
    for ra, dec, fix in zip(src_ra, src_dec, fix_var):
        region = f'circle[[{ra}pix, {dec}pix], 1.5pix]'
        ims = get_imstat_values(image, region)
        f.write(f'{ims[0]},{ra},{dec},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix}\n')
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

def check_position(fname, xpix, ypix, image, snr_thresh=5.0):
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

    # If a component is weak (i.e., < snr-thresh*sigma fix the position to the the reference otherwise fit for the position wonder)
    for x, y in zip(xpix, ypix):
        region = f'circle[[{x}pix, {y}pix], 1.5pix]'
        ims = get_imstat_values(image, region = region)
        if abs(ims[0])/ims[3] > snr_thresh:
            fix_var.append('abp')
        else:
            fix_var.append('xyabp')

    # Make the estimate file
    make_estimate(fname, image, xpix, ypix, fix_var)

def initialize_MFS_dict(imf_I, imf_Q, imf_U, imf_V, imf_P):
    '''
    Function that will take in the imfit dictionaries for all 4 (+ Lin. Pol.) Stokes Parameters, while
    returning a trimmed dictionary containing all of the parameters of interest
    '''

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
        
        # Raw fluxes
        MFS_dict[comp]['I_flux_mJy'] = imf_I['results'][comp]['peak']['value'] * 1e3
        MFS_dict[comp]['Q_flux_mJy'] = imf_Q['results'][comp]['peak']['value'] * 1e3
        MFS_dict[comp]['U_flux_mJy'] = imf_U['results'][comp]['peak']['value'] * 1e3
        MFS_dict[comp]['V_flux_mJy'] = imf_V['results'][comp]['peak']['value'] * 1e3
        MFS_dict[comp]['P_flux_mJy'] = imf_P['results'][comp]['peak']['value'] * 1e3

        # Position (Use Stokes I and Lin. Pol.)
        MFS_dict[comp]['I_RA_deg']  = imf_I['results'][comp]['shape']['direction']['m0']['value'] * 180 / np.pi
        MFS_dict[comp]['I_Dec_deg']  = imf_I['results'][comp]['shape']['direction']['m1']['value'] * 180 / np.pi
        MFS_dict[comp]['P_RA_deg'] = imf_P['results'][comp]['shape']['direction']['m0']['value'] * 180 / np.pi
        MFS_dict[comp]['P_Dec_deg'] = imf_P['results'][comp]['shape']['direction']['m1']['value'] * 180 / np.pi 

    return MFS_dict

def fit_channel(i_chan_image, xpix_I, ypix_I, xpix_Q, ypix_Q, xpix_U, ypix_U, xpix_P, ypix_P):
    '''
    Extract the IQUP parameters from the working channelized image
    Inputs:
        i_chan_image = working image
        xpix/ypix  = for the MFS image of each Stokes Parameter
    Outputs:
        Array containing the relevant flux/Polarization parameters
    '''

    # Get all of the image names:
    IQUVP_names = get_IQUVP_names(i_chan_image)

    # Fit Stokes I
    check_position('estimate_I.txt', xpix_I, ypix_I, IQUVP_names[0], snr_thresh=5.0)
    imf_I  = get_imfit_values(IQUVP_names[0], 'estimate_I.txt', xpix_I, ypix_I)

    # Fit Stokes Q
    check_position('estimate_Q.txt', xpix_Q, ypix_Q, IQUVP_names[1], snr_thresh=5.0)
    imf_Q  = get_imfit_values(IQUVP_names[1], 'estimate_Q.txt', xpix_Q, ypix_Q)

    # Fit Stokes U
    check_position('estimate_U.txt', xpix_U, ypix_U, IQUVP_names[2], snr_thresh=5.0)
    imf_U  = get_imfit_values(IQUVP_names[2], 'estimate_U.txt', xpix_U, ypix_U)

    # Fit Stokes V
    make_estimate('estimate_V.txt', IQUVP_names[3], xpix_I, ypix_I,  ['xyabp'] * len(xpix_I))
    imf_V  = get_imfit_values(IQUVP_names[3], 'estimate_V.txt', xpix_I, ypix_I)

    # Fit Stokes P
    check_position('estimate_P.txt', xpix_P, ypix_P, IQUVP_names[4], snr_thresh=4.0)
    imf_P  = get_imfit_values(IQUVP_names[4], 'estimate_P.txt', xpix_P, ypix_P)
       
    # Get the number of components
    comps = [key for key in imf_I['results'].keys() if 'component' in key]
    n_comps = len(comps)

    # Get the frequencies
    freq = np.ones(n_comps) * imhead(IQUVP_names[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9

    # Get the rms values
    region = f'circle[[{np.mean(xpix_I)}pix,{np.mean(ypix_I)}pix],1.5arcsec]'

    rms_I   =  get_imstat_values(IQUVP_names[0], region)[3] * 1e3 * np.ones(n_comps)
    rms_Q =  get_imstat_values(IQUVP_names[1], region)[3] * 1e3 * np.ones(n_comps)
    rms_U =  get_imstat_values(IQUVP_names[2], region)[3] * 1e3 * np.ones(n_comps)
    rms_V =  get_imstat_values(IQUVP_names[3], region)[3] * 1e3 * np.ones(n_comps)
    rms_P = 0.5* (rms_Q + rms_U) * np.ones(n_comps)

    # Get the fluxes
    flux_I, flux_Q, flux_U, flux_V, flux_P, flux_P0 = np.ones((6, n_comps))

    for k, comp in enumerate(comps):
        flux_I[k] = imf_I['results'][comp]['peak']['value'] * 1e3
        flux_Q[k] = imf_Q['results'][comp]['peak']['value'] * 1e3
        flux_U[k] = imf_U['results'][comp]['peak']['value'] * 1e3
        flux_V[k] = imf_V['results'][comp]['peak']['value'] * 1e3
        flux_P[k] = imf_P['results'][comp]['peak']['value'] * 1e3

        if flux_P[k] > 4.0 * rms_P[k]:
            flux_P0[k] = np.sqrt(flux_P[k] ** 2  - 2.3 * rms_P[k] ** 2)
        else:
            flux_P0[k] = flux_P[k] 

    # Calculate other parameters
    LP             = flux_P0 / flux_I * 100.0
    LP_err      =  LP * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2 )
    EVPA        = np.arctan2(flux_U, flux_Q) * 180.0 / np.pi * 0.5
    EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi
    
    return [freq, flux_I, flux_Q, flux_U, flux_V, flux_P, flux_P0, LP, EVPA, rms_I, rms_Q, rms_U, rms_V, rms_P, LP_err, EVPA_err]


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
    src_ra    = [float(ra) for ra in src_ra.split(',')]
    src_dec = [float(dec) for dec in src_dec.split(',')]

    # Make estimate files and fit for Stokes I
    make_estimate('estimate_I.txt', IQUVP_names[0], src_ra, src_dec,  ['abp'] * len(src_ra))
    imf_I  = get_imfit_values(IQUVP_names[0], 'estimate_I.txt', src_ra, src_dec)
    xpix_I = [imf_I['results'][key]['pixelcoords'][0] for key in imf_I['results'].keys() if 'component' in key]
    ypix_I = [imf_I['results'][key]['pixelcoords'][1] for key in imf_I['results'].keys() if 'component' in key]

    # Make estimate and fit for Lin. Pol. Intensity -- checking against position of stokes I components
    check_position('estimate_P.txt', xpix_I, ypix_I, IQUVP_names[4], snr_thresh=4.0) # lower threshold due to elevated P-noise
    imf_P  = get_imfit_values(IQUVP_names[4], 'estimate_P.txt', src_ra, src_dec)
    xpix_P = [imf_P['results'][key]['pixelcoords'][0] for key in imf_P['results'].keys() if 'component' in key]
    ypix_P = [imf_P['results'][key]['pixelcoords'][1] for key in imf_P['results'].keys() if 'component' in key]

    # Make estimate and fit for Q/U -- checking against position of Lin. Pol. components
    check_position('estimate_Q.txt', xpix_P, ypix_P, IQUVP_names[1], snr_thresh=5.0)
    imf_Q  = get_imfit_values(IQUVP_names[1], 'estimate_Q.txt', src_ra, src_dec)
    xpix_Q = [imf_Q['results'][key]['pixelcoords'][0] for key in imf_Q['results'].keys() if 'component' in key]
    ypix_Q = [imf_Q['results'][key]['pixelcoords'][1] for key in imf_Q['results'].keys() if 'component' in key]
   
    check_position('estimate_U.txt', xpix_P, ypix_P, IQUVP_names[2], snr_thresh=5.0)
    imf_U  = get_imfit_values(IQUVP_names[2], 'estimate_U.txt', src_ra, src_dec)
    xpix_U = [imf_U['results'][key]['pixelcoords'][0] for key in imf_U['results'].keys() if 'component' in key]
    ypix_U = [imf_U['results'][key]['pixelcoords'][1] for key in imf_U['results'].keys() if 'component' in key]

    # Make estimate for V, but fix always fix at the Stokes I position (this is just to estimate the Stokes V level)
    make_estimate('estimate_V.txt', IQUVP_names[3], xpix_I, ypix_I,  ['xyabp'] * len(src_ra))
    imf_V  = get_imfit_values(IQUVP_names[3], 'estimate_V.txt', src_ra, src_dec)

    # Initialize the MFS parameters
    MFS_dict = initialize_MFS_dict(imf_I, imf_Q, imf_U, imf_V, imf_P)

    # Get image plane RMS noise and add to dictionary
    region = f'circle[[{np.mean(src_ra)}pix,{np.mean(src_dec)}pix],1.5arcsec]'
    rms_I     =  get_imstat_values(IQUVP_names[0], region)[3] * 1e3
    rms_Q   =  get_imstat_values(IQUVP_names[1], region)[3] * 1e3
    rms_U   =  get_imstat_values(IQUVP_names[2], region)[3] * 1e3
    rms_V   =  get_imstat_values(IQUVP_names[3], region)[3] * 1e3
    rms_P   = 0.5 * (rms_Q + rms_U)

    MFS_dict['freqMFS_GHz'] = freqMFS_GHz
    MFS_dict['date_isot'] = date_obs

    # Iterate through components
    comps = [key for key in MFS_dict.keys() if 'component' in key]
    for comp in comps:
        
        # Seperate out the fluxes for readability
        MFS_dict[comp]['I_rms_mJy']   = rms_I 
        MFS_dict[comp]['Q_rms_mJy'] = rms_Q
        MFS_dict[comp]['U_rms_mJy'] = rms_U 
        MFS_dict[comp]['V_rms_mJy'] = rms_V 
        MFS_dict[comp]['P_rms_mJy'] = (rms_Q + rms_U) * 0.5     

        flux_I = MFS_dict[comp]['I_flux_mJy']
        flux_Q = MFS_dict[comp]['Q_flux_mJy']
        flux_U = MFS_dict[comp]['U_flux_mJy']
        flux_V = MFS_dict[comp]['V_flux_mJy']
        flux_P = MFS_dict[comp]['P_flux_mJy']    

        # Check if S/N of P is >=4, if so apply bias correction according to George et al. 2012: 1106.5362 and calculate Lin. Pol. Frac + Angle
        if flux_P >= 4.0 *  rms_P :
            flux_P0 = np.sqrt(flux_P ** 2 - 2.3 * rms_P ** 2)

        else: 
            flux_P0 = flux_P

        LP_frac        = flux_P0 / flux_I * 100.0
        LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2)

        LP_EVPA         = 0.5 * np.arctan2(flux_U, flux_Q) * 180.0 / np.pi 
        LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi

        # Append to dictionary 
        MFS_dict[comp]['P0_flux_mJy'] =  flux_P0

        MFS_dict[comp]['LP_frac']        = LP_frac
        MFS_dict[comp]['LP_frac_err'] = LP_frac_err

        MFS_dict[comp]['LP_EVPA']        = LP_EVPA
        MFS_dict[comp]['LP_EVPA_err'] = LP_EVPA_err
    

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Now iterate through channelized images to get the frequecny evolution   #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    # Initialize dictionary
    chan_dict = {}
    for comp in comps:
        chan_dict[comp] = {'freq_GHz' : [], 'I_flux_mJy': [], 'Q_flux_mJy': [], 'U_flux_mJy': [], 'V_flux_mJy': [], 'P_flux_mJy': [], 'P0_flux_mJy': [], 'LP_frac': [], 'LP_EVPA' : [],  
                                             'I_rms_mJy': [], 'Q_rms_mJy': [], 'U_rms_mJy': [],  'V_rms_mJy': [], 'P_rms_mJy': [], 'LP_frac_err': [], 'LP_EVPA_err':[]}
    
    # Iterate through the channelized images    
    i_chan_images = sorted(glob.glob(f'{src_im_prefix}-00*-I-image.fits'))
    for i_chan_image in i_chan_images[:]:
        try: 
            msg(f'Fitting Image: {i_chan_image}')
            chan_data = fit_channel(i_chan_image, xpix_I, ypix_I, xpix_Q, ypix_Q, xpix_U, ypix_U, xpix_P, ypix_P)

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
        except:
            msg('Fitting Failed: Channel is likely flagged')

    # Attach the sub_directories to the large directory + save as JSON
    flux_dict ['MFS'] = MFS_dict
    flux_dict ['CHAN'] = chan_dict

    with open(cfg.RESULTS + f'/{src_header}_polarization.json', 'w') as jfile:
        jfile.write(json.dumps(flux_dict, indent=4, sort_keys=True))

    # Make RM synthesis text files
    make_rmsynth_files(src_header, chan_dict)

    return 0

def main():
    
    # Load in the rmsynthesis data
    rmsynth_info = np.genfromtxt(cfg.DATA + '/rmsynth/rmsynth_info.txt', skip_header = 3, dtype=str)

    # If this is a 1-D array convert to two 2-D
    if len(rmsynth_info.shape) == 1:
        rmsynth_info = rmsynth_info[np.newaxis, :]

    # Break the columns into the relevant properties
    src_names, src_im_identifiers, src_ras, src_decs= rmsynth_info[:,0], rmsynth_info[:,1], rmsynth_info[:,2], rmsynth_info[:,3]

    for src_name, src_im_identifier, src_ra, src_dec in zip(src_names, src_im_identifiers, src_ras, src_decs):
        print('\n')
        msg(f'Starting Analysis of {src_name}')
        extract_polarization_properties(src_name, src_im_identifier, src_ra, src_dec)
   

if __name__  == "__main__":
    main()
