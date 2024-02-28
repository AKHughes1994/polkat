import glob,os,datetime, subprocess, sys, json, time
import shutil
import numpy as np
import os.path as o
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

def get_IQUVP_names(im_I):
    '''
    Take in an image prefix and return the separated IQUVP images
    '''
    
    im_Q = im_I.replace('-I-', '-Q-')
    im_U = im_I.replace('-I-', '-U-')
    im_V = im_I.replace('-I-', '-V-')
    im_P = im_I.replace('-I-', '-P-')

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

def fit_channel(i_chan_image, pix_I, pix_Q, pix_U, pos):
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

    # Seperate fluxes
    flux_I   = imf_I['results']['component0']['peak']['value'] 
    flux_Q = imf_Q['results']['component0']['peak']['value'] 
    flux_U = imf_U['results']['component0']['peak']['value'] 

    # RMS values
    rms_I   =  get_imstat_values(IQUVP_names[0], pos)[3] 
    rms_Q =  get_imstat_values(IQUVP_names[1], pos)[3]
    rms_U =  get_imstat_values(IQUVP_names[2], pos)[3]

    return [freq_Hz, flux_I, flux_Q, flux_U, rms_I, rms_Q, rms_U]   
        

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

def calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes_I = False):

    if type(flux) is list:
        flux = np.array(flux)
        rms = np.array(rms)
    
    sys_err_sq = rms ** 2 + (flux * bpcal_sys) ** 2
    if not stokes_I:
        sys_err_sq += (flux * pacal_sys) ** 2
    
    return (sys_err_sq ** (0.5)).tolist()


def get_polcal_polarization(pacal_name, pacal_pos, bpcal_systematic):
    '''
    Get the MFS parameters for the polarization calibrator and ouput a dictionary
    Inputs: 
        pacal_pos     = position (in CASA format) of calibrator (e.g., 13:31:08.2881,+30.30.32.959)
        pacal_name = name of calibrator 
    Outputs:
        MFS_dict       = Dictionary containing MFS information for the fill stokes  
    '''

    # Get imstat parameters from Stokes I image
    image_I = glob.glob(cfg.IMAGES + f'/*{pacal_name}_postXf-MFS-I-image.fits')[0]
    freqMFS_GHz  =  imhead(image_I, mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
    date_obs  =  imhead(image_I, mode='get', hdkey = 'DATE-OBS').replace('/','-',2).replace('/','T')
    ims_I       = get_imstat_values(image_I, pacal_pos)
    file_prefix = image_I.split('-MFS')[0].replace('IMAGES', 'RESULTS')

    # Now get the Full IQUV polarization properties for the polarization calibrator
    image_I, image_Q, image_U, image_V, image_P = get_IQUVP_names(image_I)

    # Make estimate files for I, Q, U, and P
    estimate_I   = make_estimate('estimate_I.txt',  image_I,  ims_I[1], ims_I[2], 'abp')
    estimate_Q = make_estimate('estimate_Q.txt', image_Q, ims_I[1], ims_I[2], 'abp')
    estimate_U = make_estimate('estimate_U.txt', image_U, ims_I[1], ims_I[2], 'abp')
    estimate_P  = make_estimate('estimate_P.txt', image_P, ims_I[1], ims_I[2], 'xyabp')

    # Imfit estimate files for I, Q, U, and P
    imf_I   = get_imfit_values('estimate_I.txt',  image_I,  ims_I[1], ims_I[2])
    imf_Q = get_imfit_values('estimate_Q.txt', image_Q, ims_I[1], ims_I[2])
    imf_U = get_imfit_values('estimate_U.txt', image_U, ims_I[1], ims_I[2])
    imf_P = get_imfit_values('estimate_P.txt', image_P, ims_I[1], ims_I[2])

    # Get imstat regions for I,Q,U,P,V
    ims_I   = get_imstat_values(image_I,  pacal_pos)
    ims_Q = get_imstat_values(image_Q, pacal_pos)
    ims_U = get_imstat_values(image_U, pacal_pos)
    ims_P = get_imstat_values(image_P, pacal_pos)
    ims_V = get_imstat_values(image_V, pacal_pos, n_beams = 1.0)

    # Get I, Q, U, P -- MFS pixel positions
    pix_I = imf_I['results']['component0']['pixelcoords']
    pix_Q = imf_I['results']['component0']['pixelcoords']
    pix_U = imf_I['results']['component0']['pixelcoords']

    #  Separate out the fluxes and the RMS values
    flux_I  = imf_I['results']['component0']['peak']['value'] * 1e3
    rms_I  = ims_I[3] * 1e3

    flux_Q = imf_Q['results']['component0']['peak']['value'] * 1e3
    rms_Q = ims_Q[3] * 1e3

    flux_U = imf_U['results']['component0']['peak']['value'] * 1e3
    rms_U = ims_U[3] * 1e3

    flux_P = imf_P['results']['component0']['peak']['value'] * 1e3
    rms_P = (rms_Q + rms_U) * 0.5

    flux_V = ims_V[0] * 1e3
    rms_V = ims_V[3] * 1e3

    # Since 3C286 is meant to be unpolarized in Stokes V, any residual Stokes V is from improper cross-hand phase calculate fractional systematic
    pacal_systematic = abs(flux_V) / (flux_P ** 2 + flux_V **2) ** (0.5)

    # Calculate error that include the systematic effects from insuffiicent leakage/cross-hand phase calibration
    sys_I = calculate_sys_err(flux_I, rms_I, bpcal_systematic, pacal_systematic, stokes_I = True)
    sys_Q = calculate_sys_err(flux_Q, rms_Q, bpcal_systematic, pacal_systematic, stokes_I = False)
    sys_U = calculate_sys_err(flux_U, rms_U, bpcal_systematic, pacal_systematic, stokes_I = False)
    sys_V = calculate_sys_err(flux_V, rms_V, bpcal_systematic, pacal_systematic, stokes_I = False)
    sys_P = 0.5 * (sys_Q + sys_U)

    # Calculate other parameters
    flux_P0 = np.sqrt(flux_P ** 2 - 2.3 * rms_P ** 2)
   
    LP_frac        = flux_P0 / flux_I * 100.0
    LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2)

    LP_EVPA         = 0.5 * np.arctan2(flux_U, flux_Q) * 180.0 / np.pi 
    LP_EVPA_err = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi


    # Make MFS dictionary
    MFS_dict = {}
    MFS_dict['beam'] = {}
    MFS_dict['beam']['bmaj_asec'] = imf_I['results']['component0']['beam']['beamarcsec']['major']['value']
    MFS_dict['beam']['bmin_asec'] = imf_I['results']['component0']['beam']['beamarcsec']['minor']['value']
    MFS_dict['beam']['bpa_deg']   = imf_I['results']['component0']['beam']['beamarcsec']['positionangle']['value']
    MFS_dict['beam']['bpa_deg']   = imf_I['results']['component0']['beam']['beamarcsec']['positionangle']['value']
        
    # Image information
    MFS_dict['freqMFS_GHz'] = freqMFS_GHz
    MFS_dict['date_isot'] = date_obs

   #  Fluxes
    MFS_dict['I_flux_mJy']   = flux_I
    MFS_dict['Q_flux_mJy'] = flux_Q
    MFS_dict['U_flux_mJy'] = flux_U
    MFS_dict['V_flux_mJy'] = flux_V
    MFS_dict['P_flux_mJy'] = flux_P
    MFS_dict['P0_flux_mJy'] = flux_P0

   #  RMS
    MFS_dict['I_rms_mJy']   = rms_I
    MFS_dict['Q_rms_mJy'] = rms_Q
    MFS_dict['U_rms_mJy'] = rms_U
    MFS_dict['V_rms_mJy'] = rms_V
    MFS_dict['P_rms_mJy'] = rms_P

   #  Systematic error
    MFS_dict['I_sys_mJy']   = sys_I
    MFS_dict['Q_sys_mJy'] = sys_Q
    MFS_dict['U_sys_mJy'] = sys_U
    MFS_dict['V_sys_mJy'] = sys_V
    MFS_dict['P_sys_mJy'] = sys_P

    # Position (Use Stokes I and Lin. Pol.)
    MFS_dict['I_RA_deg']  = imf_I['results']['component0']['shape']['direction']['m0']['value'] * 180 / np.pi
    MFS_dict['I_Dec_deg']  = imf_I['results']['component0']['shape']['direction']['m1']['value'] * 180 / np.pi
    MFS_dict['P_RA_deg'] = imf_P['results']['component0']['shape']['direction']['m0']['value'] * 180 / np.pi
    MFS_dict['P_Dec_deg'] = imf_P['results']['component0']['shape']['direction']['m1']['value'] * 180 / np.pi 

    # Polarization Parameters
    MFS_dict['LP_frac'] = LP_frac
    MFS_dict['LP_frac_err'] = LP_frac_err
    MFS_dict['LP_EVPA'] = LP_EVPA
    MFS_dict['LP_EVPA_err'] = LP_EVPA_err
    MFS_dict['Polarization_Calibrator_Fractional_Systematic'] = pacal_systematic
    MFS_dict['Primary_Calibrator_Fractional_Systematic'] = bpcal_systematic

    # Iterate through channelized images extracting fluxes to make RM Synthesis files
    rmsynth_arr = []
    for i_chan_image in sorted(glob.glob(cfg.IMAGES + f'/*{pacal_name}_postXf-00*-I-image.fits')):
        try:                
            msg(f'Fitting Image: {i_chan_image}')
            rmsynth_arr.append(fit_channel(i_chan_image, pix_I, pix_Q, pix_U, pacal_pos))

        except:
            msg('Fitting Failed: Channel is likely flagged')
 
    # Save the RMSynth file
    np.savetxt(f'/{file_prefix}_rmsynth.txt', np.array(rmsynth_arr))

    # Make a new array with the systematic noise corrections
    #rmsynth_sys_arr = np.array(rmsynth_arr)
    #rmsynth_sys_arr[:,4] = calculate_sys_err(rmsynth_sys_arr[:,1] , rmsynth_sys_arr[:,4] , bpcal_systematic, pacal_systematic, stokes_I = True)
    #rmsynth_sys_arr[:,5] = calculate_sys_err(rmsynth_sys_arr[:,2] , rmsynth_sys_arr[:,5] , bpcal_systematic, pacal_systematic, stokes_I = False)
    #rmsynth_sys_arr[:,6] = calculate_sys_err(rmsynth_sys_arr[:,3] , rmsynth_sys_arr[:,6] , bpcal_systematic, pacal_systematic, stokes_I = False)

    # Save dictionary
    with open(f'/{file_prefix}_polarization.json', 'w') as jfile:
        jfile.write(json.dumps(MFS_dict, indent=4, sort_keys=True))

    # Save RM files
    return MFS_dict

def update_rmsynth_file(fname, MFS_dict):

    bpcal_systematic = MFS_dict['Primary_Calibrator_Fractional_Systematic'] 
    pacal_systematic = MFS_dict['Polarization_Calibrator_Fractional_Systematic'] 
    
    data_arr = np.genfromtxt(fname)
    data_arr[:,4] = calculate_sys_err(data_arr[:,1] , data_arr[:,4] , bpcal_systematic, pacal_systematic, stokes_I = True)
    data_arr[:,5] = calculate_sys_err(data_arr[:,2] , data_arr[:,5] , bpcal_systematic, pacal_systematic, stokes_I = False)
    data_arr[:,6] = calculate_sys_err(data_arr[:,3] , data_arr[:,6] , bpcal_systematic, pacal_systematic, stokes_I = False)

    np.savetxt(fname.replace('.txt', '_sys.txt'), np.array(data_arr))

def update_src_pol_dict(fname, MFS_dict):

    with open(fname, 'r') as jfile:
        src_dict = json.load(jfile)

    bpcal_systematic = MFS_dict['Primary_Calibrator_Fractional_Systematic'] 
    pacal_systematic = MFS_dict['Polarization_Calibrator_Fractional_Systematic'] 

    src_dict['Polarization_Calibrator_Fractional_Systematic'] = pacal_systematic
    src_dict['Primary_Calibrator_Fractional_Systematic'] = bpcal_systematic

    #Append Systematic errors
    comps = [key for key in src_dict['MFS'].keys() if 'component' in key]
    for comp in comps:
        for subdict in ['MFS', 'CHAN']:
            src_dict[subdict][comp]['I_sys_mJy'] = calculate_sys_err(src_dict[subdict][comp]['I_flux_mJy'], src_dict[subdict][comp]['I_rms_mJy'], bpcal_systematic, pacal_systematic, stokes_I = True)
            src_dict[subdict][comp]['Q_sys_mJy'] = calculate_sys_err(src_dict[subdict][comp]['Q_flux_mJy'], src_dict[subdict][comp]['Q_rms_mJy'], bpcal_systematic, pacal_systematic, stokes_I = False)
            src_dict[subdict][comp]['U_sys_mJy'] = calculate_sys_err(src_dict[subdict][comp]['U_flux_mJy'], src_dict[subdict][comp]['U_rms_mJy'], bpcal_systematic, pacal_systematic, stokes_I = False)
            src_dict[subdict][comp]['V_sys_mJy'] = calculate_sys_err(src_dict[subdict][comp]['V_flux_mJy'], src_dict[subdict][comp]['V_rms_mJy'], bpcal_systematic, pacal_systematic, stokes_I = False)
            sys_Q = src_dict[subdict][comp]['Q_sys_mJy'] 
            sys_U = src_dict[subdict][comp]['U_sys_mJy'] 
            if type(sys_Q) is list:
                sys_P = (0.5 * np.array(sys_Q) + np.array(sys_U)).tolist()
            else:
                sys_P = 0.5 * (sys_Q + sys_U)
            src_dict[subdict][comp]['P_sys_mJy'] = sys_P

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
    image_I = glob.glob(cfg.IMAGES + f'/*{bpcal_name}_postXf-MFS-I-image.fits')[0]
    image_I, image_Q, image_U, image_V, image_P = get_IQUVP_names(image_I)
    ims_I       = get_imstat_values(image_I, bpcal_pos)
    
    # Fit Stokes I (source is very bright)
    estimate_I   = make_estimate('estimate_I.txt',  image_I,  ims_I[1], ims_I[2], 'abp')
    imf_I   = get_imfit_values('estimate_I.txt',  image_I,  ims_I[1], ims_I[2])
    flux_I = imf_I['results']['component0']['peak']['value'] 

    # Get peak pixel of QUV image
    flux_Q = get_imstat_values(image_Q, bpcal_pos, n_beams=1.0)[0]
    flux_U = get_imstat_values(image_U, bpcal_pos, n_beams=1.0)[0]
    flux_V = get_imstat_values(image_V, bpcal_pos, n_beams=1.0)[0]

    # Caculate systematic and return it
    bpcal_systematic = (flux_Q ** 2 + flux_U ** 2 + flux_V ** 2) ** 0.5 / (flux_I)

    return bpcal_systematic

def main():

    # Modify these parameters to match your calibrators
    pacal_name  =  'J1331+3030'
    pacal_pos      = '13:31:08.2881,+30.30.32.959'    

    bpcal_name = 'J1939-6342'
    bpcal_pos     =  '19:39:25.0264,-63.42.45.624'

    # Primary Calibrator
    bpcal_systematic = get_primary_systematic(bpcal_name, bpcal_pos)

    # Polarization Calibrator
    pacal_MFS_dict = get_polcal_polarization(pacal_name, pacal_pos, bpcal_systematic)

    # Update all rmsynth.txt files to include systematic errors
    for f in glob.glob(cfg.RESULTS + '/*_rmsynth.txt'):
        update_rmsynth_file(f, pacal_MFS_dict)

    # Update the source dictionary with the sysematic corrections
    rmsynth_info = np.genfromtxt(cfg.DATA + '/rmsynth/rmsynth_info.txt', skip_header = 3, dtype=str)

    # If this is a 1-D array convert to two 2-D
    if len(rmsynth_info.shape) == 1:
        rmsynth_info = rmsynth_info[np.newaxis, :]

    # Break the columns into the relevant properties
    src_names = rmsynth_info[:,0]
    src_im_identifiers = rmsynth_info[:,1]

    # Add systematic errors to src dictionary
    for src_name, src_im_identifier in zip(src_names, src_im_identifiers):
        fname = glob.glob(cfg.RESULTS + f'/*{src_name}*{src_im_identifier}*_polarization.json')[0]
        update_src_pol_dict(fname, pacal_MFS_dict)


if __name__ == "__main__":
    main()


























