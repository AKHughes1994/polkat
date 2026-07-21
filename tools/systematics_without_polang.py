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

def calculate_sys_err(flux, rms, bpcal_sys, stokes = 'I'):
    '''
    Code to calculate a systematic error using the the residual polarized signals
    in the calibrators:
        For Bandpass Calibrator (BP) this is ANY polarization
        For Polartization Angle (PA) Calibrator this is ONLY Stokes V
    Inputs:
        flux = array containing [I,Q,U,V,P] fluxes
        rms  = array containing [I,Q,U,V,P] rms
        bpcal_sys = fractional systematic error from BP
        pacal_sys = fractional systematic error from PA
        Stokes = Stokes parameter that will be solved (systematic is Stokes Depedant)
    Outputs:
        systematic error (single number or array)
    '''

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
        sys_err_sq = dU ** 2 +  (bpcal_sys * I) ** 2

    elif stokes == 'V':
        sys_err_sq = dV ** 2 +  (bpcal_sys * I) ** 2

    else:
        raise NameError("Sorry, must specify I, Q, U, or V")          
    
    return (sys_err_sq ** (0.5)).tolist()

s
def update_rmsynth_file(fname, bpcal_sys):
    '''
    Load in RMsynth text files and update them, 
    such that we apply the systematic errors from the Leakage and Cross-hand
    calibration errors
    Input:
        fname    = Input string containing the path to the rmsynth file
        MFS_dict = Dictionary containing the fractional systematic terms
    Output:
        None -- But saves a new file
    '''


    # Load in the data
    data_arr = np.genfromtxt(fname)
    
    # Solve for the individual fluxes
    I, Q, U    = data_arr[:,1], data_arr[:,2], data_arr[:,3]
    dI, dQ, dU = data_arr[:,4], data_arr[:,5], data_arr[:,6]

    # Solve for P flux
    P = (Q ** 2 + U ** 2) ** (0.5)
    
    flux = [I, Q, U, I, P] # Don't need stokes V -- dummy variable just make it I
    rms = [dI, dQ, dU, dI, dI] # Don't need dV or dP -- dummy variables

    data_arr[:,4] = calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes = 'I')
    data_arr[:,5] = calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes = 'Q')
    data_arr[:,6] = calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes = 'U')

    np.savetxt(fname.replace('.txt', '_sys.txt'), np.array(data_arr))

    return 0

def update_src_pol_dict(fname, bpcal_sys):

    with open(fname, 'r') as jfile:
        src_dict = json.load(jfile)

    src_dict['Leakage_Fractional_Systematic'] = bpcal_sys

    #Append Systematic errors
    comps = [key for key in src_dict['MFS'].keys() if 'component' in key]
    for comp in comps:
        for subdict in ['MFS', 'CHAN']:
            flux = [src_dict[subdict][comp]['I_flux_mJy'], src_dict[subdict][comp]['Q_flux_mJy'], src_dict[subdict][comp]['U_flux_mJy'], src_dict[subdict][comp]['V_flux_mJy'], src_dict[subdict][comp]['P_flux_mJy']]
            rms = [src_dict[subdict][comp]['I_rms_mJy'], src_dict[subdict][comp]['Q_rms_mJy'], src_dict[subdict][comp]['U_rms_mJy'], src_dict[subdict][comp]['V_rms_mJy'], src_dict[subdict][comp]['P_rms_mJy']]

            src_dict[subdict][comp]['I_sys_mJy'] = calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes = 'I')
            src_dict[subdict][comp]['Q_sys_mJy'] = calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes = 'Q')
            src_dict[subdict][comp]['U_sys_mJy'] = calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes = 'U')
            src_dict[subdict][comp]['V_sys_mJy'] = calculate_sys_err(flux, rms, bpcal_sys, pacal_sys, stokes = 'V')
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
    image_I = glob.glob(cfg.IMAGES + f'/*{bpcal_name}*postXf-MFS-I-image.fits')[0]
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
    bpcal_sys = (flux_Q ** 2 + flux_U ** 2 + flux_V ** 2) ** 0.5 / (flux_I)

    return bpcal_sys

def main():

    # Modify these parameters to match your calibrators
    bpcal_name = 'J1939-6342'
    bpcal_pos     =  '19:39:25.0264,-63.42.45.624'

    # Primary Calibrator
    msg(f'Fitting Primary Systematics')
    bpcal_sys = get_primary_systematic(bpcal_name, bpcal_pos)

    # Update all rmsynth.txt files to include systematic errors
    for f in glob.glob(cfg.RESULTS + '/*_rmsynth.txt'):
        update_rmsynth_file(f, bpcal_sys)

    # Update the source dictionary with the sysematic corrections
    rmsynth_info = np.genfromtxt(cfg.DATA + '/rmsynth/rmsynth_info.txt', skip_header = 2, dtype=str)

    # If this is a 1-D array convert to two 2-D
    rmsynth_info = np.atleast_2d(rmsynth_info)

    # Break the columns into the relevant properties
    src_names = rmsynth_info[:,0]
    src_im_identifiers = rmsynth_info[:,1]

    # Add systematic errors to src dictionary
    for src_name, src_im_identifier in zip(src_names, src_im_identifiers):
        fname = glob.glob(cfg.RESULTS + f'/*{src_name}*{src_im_identifier}*_polarization.json')[0]
        update_src_pol_dict(fname, bpcal_sys)


if __name__ == "__main__":
    main()


























