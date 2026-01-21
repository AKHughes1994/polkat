# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import os
import datetime
import subprocess
import sys
import json
import shutil
import numpy as np
import os.path as o
import time
from scipy.spatial.distance import cdist
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush=True)


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
    
    # For multi-components (this won't work if there are components that are VERY far from eachother largely because it will just time out)
    else:
        # Define an array of coorindates [[x1,y1], [x2,y2], etc.])
        point_array =  np.array([xpix, ypix]).T

        # Get maximum distance between points
        max_dist = np.amax(cdist(point_array,point_array)) * pixel_asec
    
        # Take either twice the maximum distance or 10 times the bmaj axis as the bounding region radius
        x = xpix[0]
        y = ypix[0]
        r = np.amax((5 * bmaj, 2.0 * max_dist))
        src_region = f'circle[[{x}pix,{y}pix],{r}arcsec]'
   
    return imfit(image, estimates = fname, region = src_region)
    
    

def calculate_P0(flux_P, rms_Q, rms_U, rms_V, pol_flag, Aq = 0.8):
    '''
    Calculate the de-biased linearly polarized flux
    '''    

    # If polarization angle calibrator incluided
    if pol_flag:

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

    # If there is no polarization angle calibrator
    else:
        
        # This doesn't have studies that I can find (Last Update: Mar 27, 2024) -- Going conservative until I do this properly
        rms_P = np.amax([rms_Q, rms_U, rms_V]) # adopt maximum
        
        # Always de-bias - from my own Mote Carlo experiments it seems like the bias correction becomes a factor of 2 for P^2 = Q^2 + U^2 + V^2  
        P0 =  (flux_P ** 2 - 2.0 * rms_P ** 2) ** (0.5)

    return flux_P0, rms_P

        

def return_max(im, region):
    '''
    Return the value that has the higher absolute magnitude
    Necessary for fluxes that are non-positive definate
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
    r_in  = 3.0 * bmaj
    r_out = np.sqrt(100 * 0.25 * bmaj * bmin + r_in ** 2)
    src_region = f'circle[[{xpix}pix,{ypix}pix],3.0pix]'
    rms_region = f'annulus[[{xpix}pix,{ypix}pix],[{r_in}arcsec,{r_out}arcsec]]'
    if manual_rms_region:
        rms_region = manual_rms_region
    #else:
    #    msg('Using default annular RMS region')

    # Values of interest -- Source
    ims = imstat(image, region = src_region)
    flux = return_max(image, src_region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]

    # Extract RMS
    rms = imstat(image, region = rms_region)['rms'][0]

    return [flux, xpix, ypix, rms]
    
    
def check_position(fname, image, xpix, ypix, snr_thresh = 5.0, P_image = False, fix_additional_comps = False, manual_rms_region = False):

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

    fix_var = []

    # If a component is weak (i.e., < snr-thresh * sigma fix the position to the the reference otherwise fit for the position)
    k = 0
    
    # Get beam parameters
    bmaj = imhead(image, mode='get', hdkey='bmaj')['value']
    bmin = imhead(image, mode='get', hdkey='bmin')['value']
    bpa  = imhead(image, mode='get', hdkey='bpa')['value']

    for x, y in zip(xpix, ypix):
        region = f'circle[[{x}pix,{y}pix],{2 * bmaj}arcsec]'

        f = open('check_pos.txt', 'w')
        f.write(f'0.0,{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, xyabp')
        f.close()

        # Get flux at test position
        test_flux = abs(imfit(image, region = region, estimates='check_pos.txt')['results']['component0']['peak']['value'])

        # If its a P-image don't use the image plane noise as the check criteria as it is (very) non-gaussian
        rms = get_imstat_values(image, x, y, manual_rms_region = manual_rms_region)[3]
        
        if P_image is True:
            # msg(f'Check image: {image}')
            image_Q = image.replace('-Plin-', '-Q-').replace('-Ptot-', '-Q-')
            image_U = image.replace('-Plin-', '-U-').replace('-Ptot-', '-U-')
            ims_Q = get_imstat_values(image_Q, x, y, manual_rms_region = manual_rms_region)
            ims_U = get_imstat_values(image_U, x, y, manual_rms_region = manual_rms_region)
            rms = np.amax((ims_Q[3],ims_U[3]))
        if fix_additional_comps is True and k > 0:
            fix_var.append('xyabp')
            msg(f'Fixing position of component {k} because manual flag set: {image}')
        elif fix_additional_comps is True and k == 0 and test_flux > snr_thresh * rms:
            fix_var.append('abp')
            msg(f'Free fitting position of component {k} with S/N = {test_flux/rms:.2f}: {image}')
        elif fix_additional_comps is False and test_flux > snr_thresh * rms:
            fix_var.append('abp')
            msg(f'Free fitting position of component {k} with S/N = {test_flux/rms:.2f}: {image}')
        else:
            fix_var.append('xyabp')
            msg(f'Fixing position of component {k} due to S/N = {test_flux/rms:.2f}: {image}')
        k+=1

    # Make the estimate file
    make_estimate(fname, image, xpix, ypix, fix_var)


def make_estimate(fname, image, xpix, ypix, fix_var):
    '''
    Take in an array of imstat values from f(get_imstat_values)
    and return an CASA imfit estimate file name fname
    Inputs:
        fname  = string containing name of estimate file
        image  = string containing name of image to fit
        xpix    = Right acension (pixel) estimate of source(s)
        ypix  = Declination (pixel) estimate of source(s)
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


def extract_polarization_properties(src_name,
    src_im_identifier,
    src_im_suffix, 
    src_ra, 
    src_dec, 
    src_ulims,
    pol_flag, 
    manual_rms_region, 
    image_directory,
    image_identifier,
    fix_additional_comps = False):

    '''
    Fit the Stokes IQUV cube for all components in an image. This assumes
    that the Q, U, and V images follow the WSCLEAN naming
    convention. This should work for slightly extended emission where
    you have two components that are yet to be separated by > 1 beam

    input parameters:
        src_name       = name of source 
        src_im_prefix  = identifier for images that will have fluxe extraction
        src_im_suffix  = image suffix, included to differentiate between standard WSCLEAN image products (image.fits) and homogenized produces (image.homogenized.fits)
        src_ra   = Estimated right acension of source(s) pixel units
        src_dec  = Estimated declination of  source(s) pixel units
        pol_flag = determines whether or not to solve for total or linear/circular polarization (depending on if cross-hand calibrator was included)
        manual_rms_region = option to specify the rms region, otherwise use an annuls centerd on the source(s) with a ~100xPSF area
        fix_addition_comps[default=True] = option to fix the Q,U,V,P position to Stokes I (primarily for faint partially unresolved ejecta)

    Output parametres:
        flux_dict = Dictionary containing all MFS and per-channel information for the IQUV fluxes
        rmsynth_arr = Array containing the necessary information to run RM synthesis on each image
    ''' 

    # Initilaize dictionary to contain the output parameters
    output_dictionary = {'name' : src_name}
    output_dictionary['MFS'] = {}
    output_dictionary['CHAN'] = {}

    # Check if the MFS image exists for the suffix, if not, use homgenized one
    mfs_im_suffix = src_im_suffix[:]
    if glob.glob(f'{src_im_identifier}*-MFS*{src_im_suffix}') == []:
        msg(f'WARNING: {src_im_suffix} does not exist in MFS image due to channel splitting; Using: image.homogenized.fits')
        mfs_im_suffix = 'image.homogenized.fits' # if this doesn't work you don't have the images required for this analysis
    
    # Determine if there are multiple stokes parameters or if its strictly Stokes I (suffixes will be, e.g., MFS-image)
    if glob.glob(f'{src_im_identifier}*-MFS-I-{mfs_im_suffix}') == []:
        only_intensity = True
    else:
        only_intensity = False
    
    # Get unique prefixes, this will now iterate in time (if applicable) if not it will just return arrays of length 1:
    prefix_arr = glob.glob(f'{src_im_identifier}*-MFS*{mfs_im_suffix}')
    prefix_arr = sorted(list(set([x.split('-MFS')[0] for x in prefix_arr])))
    
    for k, prefix in enumerate(prefix_arr[:]):
    
        # Extract the MFS image parameters
        msg(f'Fitting MFS image(s) for prefix {k}: {prefix}')
        
        # This is to get the image into an array
        if only_intensity:
            MFS_images = [f'{prefix}-MFS-{mfs_im_suffix}']
        else:
            # The sorted function will make order the images as I, P, Q, U, V
            MFS_images = glob.glob(f'{prefix}-MFS-*-{mfs_im_suffix}')
            if pol_flag:
                MFS_images = sorted([im for im in MFS_images if '-Ptot-' not in im])
            else:
                MFS_images = sorted([im for im in MFS_images if '-Plin-' not in im])

        # Output image name
        for MFS_image in MFS_images:
            msg(f'Fitting MFS image name(s): {MFS_image}')
    
        # Get generalized properties from the first image (i.e., Stokes I header)
        freq_GHz = imhead(MFS_images[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
        date_obs = imhead(MFS_images[0], mode='get', hdkey = 'DATE-OBS').replace('/','-',2).replace('/','T')
        bmaj = imhead(MFS_images[0], mode='get', hdkey = 'bmaj')['value']
        bmin = imhead(MFS_images[0], mode='get', hdkey = 'bmin')['value']
        bpa   = imhead(MFS_images[0], mode='get', hdkey = 'bpa')['value']
    
        # Convert the RA/DEC guesses to pixel coordinates
        src_ra_pix = []
        src_dec_pix = []
        for z in range(src_ra.size):
            region = 'circle[[{}deg,{}deg],1.0pix]'.format(src_ra[z], src_dec[z])
            pixel_imstat = imstat(MFS_images[0], region = region)
            src_ra_pix.append(pixel_imstat['maxpos'][0])
            src_dec_pix.append(pixel_imstat['maxpos'][1])
            
        # First fit Stokes I -- also extract pixel coordinates
        fix = []
        for z, boolean in enumerate(src_ulims):        
            msg(f'For source {src_name} component {z} has upper limit = {boolean}')
            if boolean:
                fix.append('xyabp')
            else:
                fix.append('abp')
    
        make_estimate('estimate_I.txt', MFS_images[0], src_ra_pix, src_dec_pix,  fix)
        MFS_I_imfit = get_imfit_values('estimate_I.txt', MFS_images[0], src_ra_pix, src_dec_pix)

        # Get number of components
        components = [key for key in MFS_I_imfit['results'].keys() if 'component' in key]              
        MFS_I_ra_pix = [MFS_I_imfit['results'][key]['pixelcoords'][0] for key in components]
        MFS_I_dec_pix = [MFS_I_imfit['results'][key]['pixelcoords'][1] for key in components]
       
        # Initialize arrays if they don't exist, else append values to existsing arrays
        for z, component in enumerate(components):
        
            # For ease of readability define the desired quantities as variables
            flux_I = MFS_I_imfit['results'][component]['peak']['value'] * 1e3
            err_I = MFS_I_imfit['results'][component]['peak']['error'] * 1e3
            RA_I   = MFS_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
            DEC_I = MFS_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
            RA_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][0]
            DEC_pix_I = MFS_I_imfit['results'][component]['pixelcoords'][1]
            rms_I = get_imstat_values(MFS_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
       
            if component not in output_dictionary['MFS']:

                output_dictionary['MFS'][component] = {}

                output_dictionary['MFS'][component]['upperlimit'] = [src_ulims[z]]           
                output_dictionary['MFS'][component]['freq_GHz'] = [freq_GHz]
                output_dictionary['MFS'][component]['date_isot'] = [date_obs]
                output_dictionary['MFS'][component]['bmaj_asec'] = [bmaj]
                output_dictionary['MFS'][component]['bmin_asec'] = [bmin]
                output_dictionary['MFS'][component]['bpa_deg'] = [bpa]
                output_dictionary['MFS'][component]['I_flux_mJy'] = [flux_I]
                output_dictionary['MFS'][component]['I_err_mJy'] = [err_I]
                output_dictionary['MFS'][component]['I_rms_mJy'] = [rms_I]
                output_dictionary['MFS'][component]['I_RA_deg'] = [RA_I]
                output_dictionary['MFS'][component]['I_DEC_deg'] = [DEC_I]
                         
            else:
                output_dictionary['MFS'][component]['freq_GHz'].append(freq_GHz)
                output_dictionary['MFS'][component]['date_isot'].append(date_obs)
                output_dictionary['MFS'][component]['bmaj_asec'].append(bmaj)
                output_dictionary['MFS'][component]['bmin_asec'].append(bmin)
                output_dictionary['MFS'][component]['bpa_deg'].append(bpa)
                output_dictionary['MFS'][component]['I_flux_mJy'].append(flux_I)
                output_dictionary['MFS'][component]['I_err_mJy'].append(err_I)
                output_dictionary['MFS'][component]['I_rms_mJy'].append(rms_I)
                output_dictionary['MFS'][component]['I_RA_deg'].append(RA_I)
                output_dictionary['MFS'][component]['I_DEC_deg'].append(DEC_I)
            
        
        # Next fit Polarization intensity and check if source is bright enough to allow positions to vary freely
        if not only_intensity:
            
            # Run a first pass
            try:
                check_position('estimate_P.txt', MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix, P_image = True, fix_additional_comps = False, manual_rms_region = manual_rms_region)
                MFS_P_imfit  = get_imfit_values('estimate_P.txt', MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
                MFS_P_ra_pix = [MFS_P_imfit['results'][key]['pixelcoords'][0] for key in components]
                MFS_P_dec_pix = [MFS_P_imfit['results'][key]['pixelcoords'][1] for key in components]   
           
                # If any of the P components are very far from Stokes I refit fixing the additional comps 
                dist_check = False
                dist_max = 0.0
                for z in range(len(MFS_P_ra_pix)):
                    dist = ((MFS_I_ra_pix[z] - MFS_P_ra_pix[z]) ** 2 + (MFS_I_dec_pix[z] - MFS_P_dec_pix[z]) ** 2) ** 0.5
                    if dist > 3.0: # should not be >3 pixels
                        dist_check = True
                        if dist > dist_max:
                            dist_max = dist
                            
            except Exception as e:
                msg(f'Free fitting failed, assuming multi-component fit to complex re-trying fixing the secondary component positions to Stokes I: {e}')
                dist_max = -1.0
                dist_check = True

            if dist_check:
                msg(f'Significant offset between Stokes I and P {dist_max:.1f} pix, re-fitting while fixing multi-component')
                check_position('estimate_P.txt', MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix, P_image = True, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
                MFS_P_imfit  = get_imfit_values('estimate_P.txt', MFS_images[1], MFS_I_ra_pix, MFS_I_dec_pix)
                MFS_P_ra_pix = [MFS_P_imfit['results'][key]['pixelcoords'][0] for key in components]
                MFS_P_dec_pix = [MFS_P_imfit['results'][key]['pixelcoords'][1] for key in components]   

            # Make estimate and fit for Q/U -- checking against position of Lin. Pol. components
            check_position('estimate_Q.txt', MFS_images[2], MFS_P_ra_pix, MFS_P_dec_pix, P_image = False, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
            MFS_Q_imfit  = get_imfit_values('estimate_Q.txt', MFS_images[2], MFS_P_ra_pix, MFS_P_dec_pix)
            MFS_Q_ra_pix = [MFS_Q_imfit['results'][key]['pixelcoords'][0] for key in components]
            MFS_Q_dec_pix = [MFS_Q_imfit['results'][key]['pixelcoords'][1] for key in components]   
   
            check_position('estimate_U.txt', MFS_images[3], MFS_P_ra_pix, MFS_P_dec_pix, P_image = False, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
            MFS_U_imfit  = get_imfit_values('estimate_U.txt', MFS_images[3], MFS_P_ra_pix, MFS_P_dec_pix)
            MFS_U_ra_pix = [MFS_U_imfit['results'][key]['pixelcoords'][0] for key in components]
            MFS_U_dec_pix = [MFS_U_imfit['results'][key]['pixelcoords'][1] for key in components]   

            # Make estimate and fit Stokes V  -- checking against position of stokes I components
            check_position('estimate_V.txt', MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix, P_image = False, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
            MFS_V_imfit  = get_imfit_values('estimate_V.txt', MFS_images[4], MFS_I_ra_pix, MFS_I_dec_pix)
            MFS_V_ra_pix = [MFS_V_imfit['results'][key]['pixelcoords'][0] for key in components]
            MFS_V_dec_pix = [MFS_V_imfit['results'][key]['pixelcoords'][1] for key in components]   
  
            # Once again initialize arrays that don't exist
            for component in components[:]:
            
                # For ease and readability separate out the parameters for calculations and rms extraction
                flux_I = MFS_I_imfit['results'][component]['peak']['value'] * 1e3
                flux_P = MFS_P_imfit['results'][component]['peak']['value'] * 1e3
                flux_Q = MFS_Q_imfit['results'][component]['peak']['value'] * 1e3
                flux_U = MFS_U_imfit['results'][component]['peak']['value'] * 1e3
                flux_V = MFS_V_imfit['results'][component]['peak']['value'] * 1e3
                
                err_P = MFS_P_imfit['results'][component]['peak']['error'] * 1e3
                err_Q = MFS_Q_imfit['results'][component]['peak']['error'] * 1e3
                err_U = MFS_U_imfit['results'][component]['peak']['error'] * 1e3
                err_V = MFS_V_imfit['results'][component]['peak']['error'] * 1e3

                RA_P   = MFS_P_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
                DEC_P = MFS_P_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
                
                rms_I = get_imstat_values(MFS_images[1], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                rms_Q = get_imstat_values(MFS_images[2], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                rms_U = get_imstat_values(MFS_images[3], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                rms_V = get_imstat_values(MFS_images[4], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                

                # Calculate additional linear polarisation parameters
                flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U, rms_V, pol_flag, Aq = 0.8)
                
                # Calculate the other Polarisation parameters
                LP_frac     = flux_P0 / flux_I * 100.0
                LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2 )

                if pol_flag:
                    LP_EVPA     = np.arctan2(flux_U, flux_Q) * 180.0 / np.pi * 0.5
                    LP_EVPA_err = np.sqrt(flux_U ** 2 * rms_Q ** 2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi * 0.5
   
                else:
                    LP_EVPA = None
                    LP_EVPA_err = None  

                # Append values to dictionary -- Initialize if it doesn't exist
                if 'P_flux_mJy' not in output_dictionary['MFS'][component]:

                    output_dictionary['MFS'][component]['Q_flux_mJy'] = [flux_Q]            
                    output_dictionary['MFS'][component]['U_flux_mJy'] = [flux_U]                
                    output_dictionary['MFS'][component]['V_flux_mJy'] = [flux_V]                
                    output_dictionary['MFS'][component]['P_flux_mJy'] = [flux_P]           
                    output_dictionary['MFS'][component]['P0_flux_mJy'] = [flux_P0]    

                    output_dictionary['MFS'][component]['Q_err_mJy'] = [err_Q]            
                    output_dictionary['MFS'][component]['U_err_mJy'] = [err_U]                
                    output_dictionary['MFS'][component]['V_err_mJy'] = [err_V]                
                    output_dictionary['MFS'][component]['P_err_mJy'] = [err_P]           

                    output_dictionary['MFS'][component]['Q_rms_mJy'] = [rms_Q]            
                    output_dictionary['MFS'][component]['U_rms_mJy'] = [rms_U]                
                    output_dictionary['MFS'][component]['V_rms_mJy'] = [rms_V]                
                    output_dictionary['MFS'][component]['P_rms_mJy'] = [rms_P]        

                    output_dictionary['MFS'][component]['P_RA_deg'] = [RA_P]
                    output_dictionary['MFS'][component]['P_DEC_deg'] = [DEC_P]    

                    output_dictionary['MFS'][component]['LP_frac'] = [LP_frac]
                    output_dictionary['MFS'][component]['LP_frac_err'] = [LP_frac_err]
                    output_dictionary['MFS'][component]['LP_EVPA'] = [LP_EVPA]
                    output_dictionary['MFS'][component]['LP_EVPA_err'] = [LP_EVPA_err]        
            
                else:

                    output_dictionary['MFS'][component]['Q_flux_mJy'].append(flux_Q)         
                    output_dictionary['MFS'][component]['U_flux_mJy'].append(flux_U)               
                    output_dictionary['MFS'][component]['V_flux_mJy'].append(flux_V)               
                    output_dictionary['MFS'][component]['P_flux_mJy'].append(flux_P)           
                    output_dictionary['MFS'][component]['P0_flux_mJy'].append(flux_P0)    

                    output_dictionary['MFS'][component]['Q_err_mJy'].append(err_Q)            
                    output_dictionary['MFS'][component]['U_err_mJy'].append(err_U)               
                    output_dictionary['MFS'][component]['V_err_mJy'].append(err_V)                
                    output_dictionary['MFS'][component]['P_err_mJy'].append(err_P)           

                    output_dictionary['MFS'][component]['Q_rms_mJy'].append(rms_Q)            
                    output_dictionary['MFS'][component]['U_rms_mJy'].append(rms_U)                
                    output_dictionary['MFS'][component]['V_rms_mJy'].append(rms_V)                
                    output_dictionary['MFS'][component]['P_rms_mJy'].append(rms_P)        

                    output_dictionary['MFS'][component]['P_RA_deg'].append(RA_P)
                    output_dictionary['MFS'][component]['P_DEC_deg'].append(DEC_P)    

                    output_dictionary['MFS'][component]['LP_frac'].append(LP_frac)
                    output_dictionary['MFS'][component]['LP_frac_err'].append(LP_frac_err)
                    output_dictionary['MFS'][component]['LP_EVPA'].append(LP_EVPA)
                    output_dictionary['MFS'][component]['LP_EVPA_err'].append(LP_EVPA_err)     
                

        ###########################
        # Extract the CHAN image parameters  #
        ###########################

        msg(f'Fitting CHAN image(s) for prefix {k} with {src_im_suffix}: {prefix}')

        # Glob the images
        CHAN_images = sorted(glob.glob(f'{prefix}-[!MFS]*-{src_im_suffix}'))

        # Reshape to so that each component is a set of Stokes parameters
        if only_intensity:
            CHAN_images_arr = np.array(CHAN_images).reshape(len(CHAN_images), 1)

        else:
            if pol_flag:
               CHAN_images = sorted([im for im in CHAN_images if '-Ptot-' not in im])
            else:
               CHAN_images = sorted([im for im in CHAN_images if '-Plin-' not in im])

            CHAN_images_arr = np.array(CHAN_images).reshape(int(len(CHAN_images) / 5), 5) # reshape to group in frequency for each set of Stokes paramete

        # Iterate through the frequency channels, fit, and append to output dictionary
        for CHAN_images in CHAN_images_arr[:]:
         
            try:

                msg(f'Fitting image: {CHAN_images[0]}')

                # First fit Stokes I amd get the fix coordinates
                bmaj = imhead(CHAN_images[0], mode='get', hdkey = 'bmaj')['value']
                bmin = imhead(CHAN_images[0], mode='get', hdkey = 'bmin')['value']
                bpa   = imhead(CHAN_images[0], mode='get', hdkey = 'bpa')['value']           
                freq_GHz = imhead(CHAN_images[0], mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9     

                # Check if the beam is significantly elongated with respect to the expectation from MFS image (Channel is probably very flagged and should be omitted)
                beam_scaled = output_dictionary['MFS']['component0']['freq_GHz'][k] / freq_GHz   * output_dictionary['MFS']['component0']['bmaj_asec'][k]
                if bmaj > 3.0 * beam_scaled:
                    msg('Skipping Channel: BMAJ is very different from expectation based on MFS image (likely high flagged)')
                    continue 

                check_position('estimate_I.txt', CHAN_images[0], MFS_I_ra_pix, MFS_I_dec_pix, manual_rms_region = manual_rms_region)
                CHAN_I_imfit = get_imfit_values('estimate_I.txt', CHAN_images[0], MFS_I_ra_pix, MFS_I_dec_pix)
                CHAN_I_ra_pix = [CHAN_I_imfit['results'][key]['pixelcoords'][0] for key in components]
                CHAN_I_dec_pix = [CHAN_I_imfit['results'][key]['pixelcoords'][1] for key in components]

                # Initialize arrays if they don't exist, else append values to existsing arrays
                for component in components:
        
                    # For ease of readability define the desired quantities as variables
                    flux_I = CHAN_I_imfit['results'][component]['peak']['value'] * 1e3
                    err_I = CHAN_I_imfit['results'][component]['peak']['error'] * 1e3
                    RA_I   = CHAN_I_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
                    DEC_I = CHAN_I_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
                    RA_pix_I = CHAN_I_imfit['results'][component]['pixelcoords'][0]
                    DEC_pix_I = CHAN_I_imfit['results'][component]['pixelcoords'][1]
                    rms_I = get_imstat_values(CHAN_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3
        
                    # Append values to the ouput dictionary
                    if component not in output_dictionary['CHAN']:

                        output_dictionary['CHAN'][component] = {}
                
                        output_dictionary['CHAN'][component]['freq_GHz'] = [[freq_GHz]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['bmaj_asec'] = [[bmaj]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['bmin_asec'] = [[bmin]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['bpa_deg'] = [[bpa]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['I_flux_mJy'] = [[flux_I]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['I_err_mJy'] = [[err_I]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['I_rms_mJy'] = [[rms_I]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['I_RA_deg'] = [[RA_I]] + (len(prefix_arr) - 1) * [[]]
                        output_dictionary['CHAN'][component]['I_DEC_deg'] = [[DEC_I]] + (len(prefix_arr) - 1) * [[]]
                         
                    else:
                        output_dictionary['CHAN'][component]['freq_GHz'][k].append(freq_GHz)
                        output_dictionary['CHAN'][component]['bmaj_asec'][k].append(bmaj)
                        output_dictionary['CHAN'][component]['bmin_asec'][k].append(bmin)
                        output_dictionary['CHAN'][component]['bpa_deg'][k].append(bpa)
                        output_dictionary['CHAN'][component]['I_flux_mJy'][k].append(flux_I)
                        output_dictionary['CHAN'][component]['I_err_mJy'][k].append(err_I)
                        output_dictionary['CHAN'][component]['I_rms_mJy'][k].append(rms_I)
                        output_dictionary['CHAN'][component]['I_RA_deg'][k].append(RA_I)
                        output_dictionary['CHAN'][component]['I_DEC_deg'][k].append(DEC_I)


            
                # Next fit Polarization intensity and check if source is bright enough to allow positions to vary freely
                if not only_intensity:

                    check_position('estimate_P.txt', CHAN_images[1], MFS_P_ra_pix, MFS_P_dec_pix, P_image = True, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
                    CHAN_P_imfit  = get_imfit_values('estimate_P.txt', CHAN_images[1], CHAN_I_ra_pix, CHAN_I_dec_pix)
                    CHAN_P_ra_pix = [CHAN_P_imfit['results'][key]['pixelcoords'][0] for key in components]
                    CHAN_P_dec_pix = [CHAN_P_imfit['results'][key]['pixelcoords'][1] for key in components]   
                
                    # Make estimate and fit for Q/U -- checking against position of Lin. Pol. components
                    check_position('estimate_Q.txt', CHAN_images[2], CHAN_P_ra_pix, CHAN_P_dec_pix, P_image = False, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
                    CHAN_Q_imfit  = get_imfit_values('estimate_Q.txt', CHAN_images[2], CHAN_P_ra_pix, CHAN_P_dec_pix)
                    CHAN_Q_ra_pix = [CHAN_Q_imfit['results'][key]['pixelcoords'][0] for key in components]
                    CHAN_Q_dec_pix = [CHAN_Q_imfit['results'][key]['pixelcoords'][1] for key in components]   
   
                    check_position('estimate_U.txt', CHAN_images[3], CHAN_P_ra_pix, CHAN_P_dec_pix, P_image = False, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
                    CHAN_U_imfit  = get_imfit_values('estimate_U.txt', CHAN_images[3], CHAN_P_ra_pix, CHAN_P_dec_pix)
                    CHAN_U_ra_pix = [CHAN_U_imfit['results'][key]['pixelcoords'][0] for key in components]
                    CHAN_U_dec_pix = [CHAN_U_imfit['results'][key]['pixelcoords'][1] for key in components]   

                    # Make estimate and fit Stokes V  -- checking against position of stokes I components
                    check_position('estimate_V.txt', CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix, P_image = False, fix_additional_comps = fix_additional_comps, manual_rms_region = manual_rms_region)
                    CHAN_V_imfit  = get_imfit_values('estimate_V.txt', CHAN_images[4], CHAN_I_ra_pix, CHAN_I_dec_pix)
                    CHAN_V_ra_pix = [CHAN_V_imfit['results'][key]['pixelcoords'][0] for key in components]
                    CHAN_V_dec_pix = [CHAN_V_imfit['results'][key]['pixelcoords'][1] for key in components]   
  
                    # Once again initialize arrays that don't exist
                    for component in components:

                        RA_pix_I = CHAN_I_imfit['results'][component]['pixelcoords'][0]
                        DEC_pix_I = CHAN_I_imfit['results'][component]['pixelcoords'][1]
            
                        # For ease and readability separate out the parameters for calculations and rms extraction
                        flux_I = CHAN_I_imfit['results'][component]['peak']['value'] * 1e3
                        flux_P = CHAN_P_imfit['results'][component]['peak']['value'] * 1e3
                        flux_Q = CHAN_Q_imfit['results'][component]['peak']['value'] * 1e3
                        flux_U = CHAN_U_imfit['results'][component]['peak']['value'] * 1e3
                        flux_V = CHAN_V_imfit['results'][component]['peak']['value'] * 1e3
                
                        err_P = CHAN_P_imfit['results'][component]['peak']['error'] * 1e3
                        err_Q = CHAN_Q_imfit['results'][component]['peak']['error'] * 1e3
                        err_U = CHAN_U_imfit['results'][component]['peak']['error'] * 1e3
                        err_V = CHAN_V_imfit['results'][component]['peak']['error'] * 1e3
    
                        RA_P   = CHAN_P_imfit['results'][component]['shape']['direction']['m0']['value'] * 180 / np.pi
                        DEC_P = CHAN_P_imfit['results'][component]['shape']['direction']['m1']['value'] * 180 / np.pi
                    
                        rms_I = get_imstat_values(CHAN_images[0], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                        rms_Q = get_imstat_values(CHAN_images[2], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                        rms_U = get_imstat_values(CHAN_images[3], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
                        rms_V = get_imstat_values(CHAN_images[4], RA_pix_I, DEC_pix_I, manual_rms_region)[3] * 1e3                
    
                        # Calculate additional linear polarisation parameters
                        flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U, rms_V, pol_flag, Aq = 0.8)
                    
                        # Calculate the other Polarisation parameters
                        LP_frac     = flux_P0 / flux_I * 100.0
                        LP_frac_err = LP_frac * np.sqrt( (rms_I / flux_I) ** 2 + (rms_P / flux_P0) ** 2 )

                        if pol_flag:
                            LP_EVPA     = np.arctan2(flux_U, flux_Q) * 180.0 / np.pi * 0.5
                            LP_EVPA_err = np.sqrt(flux_U ** 2 * rms_Q ** 2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi * 0.5
                        else:
                            LP_EVPA = None
                            LP_EVPA_err = None  

                        # Append values to dictionary -- Initialize if it doesn't exist
                        if 'P_flux_mJy' not in output_dictionary['CHAN'][component]:
    
                            output_dictionary['CHAN'][component]['Q_flux_mJy'] = [[flux_Q]] + (len(prefix_arr) - 1) * [[]]      
                            output_dictionary['CHAN'][component]['U_flux_mJy'] = [[flux_U]] + (len(prefix_arr) - 1) * [[]]              
                            output_dictionary['CHAN'][component]['V_flux_mJy'] = [[flux_V]] + (len(prefix_arr) - 1) * [[]]              
                            output_dictionary['CHAN'][component]['P_flux_mJy'] = [[flux_P]] + (len(prefix_arr) - 1) * [[]]         
                            output_dictionary['CHAN'][component]['P0_flux_mJy'] = [[flux_P0]] + (len(prefix_arr) - 1) * [[]]   
    
                            output_dictionary['CHAN'][component]['Q_err_mJy'] = [[err_Q]] + (len(prefix_arr) - 1) * [[]]            
                            output_dictionary['CHAN'][component]['U_err_mJy'] = [[err_U]]  + (len(prefix_arr) - 1) * [[]]              
                            output_dictionary['CHAN'][component]['V_err_mJy'] = [[err_V]]  + (len(prefix_arr) - 1) * [[]]              
                            output_dictionary['CHAN'][component]['P_err_mJy'] = [[err_P]]  + (len(prefix_arr) - 1) * [[]]         

                            output_dictionary['CHAN'][component]['Q_rms_mJy'] = [[rms_Q]] + (len(prefix_arr) - 1) * [[]]            
                            output_dictionary['CHAN'][component]['U_rms_mJy'] = [[rms_U]]  + (len(prefix_arr) - 1) * [[]]              
                            output_dictionary['CHAN'][component]['V_rms_mJy'] = [[rms_V]]  + (len(prefix_arr) - 1) * [[]]              
                            output_dictionary['CHAN'][component]['P_rms_mJy'] = [[rms_P]]  + (len(prefix_arr) - 1) * [[]]      

                            output_dictionary['CHAN'][component]['P_RA_deg'] = [[RA_P]] + (len(prefix_arr) - 1) * [[]]
                            output_dictionary['CHAN'][component]['P_DEC_deg'] = [[DEC_P]] + (len(prefix_arr) - 1) * [[]]    

                            output_dictionary['CHAN'][component]['LP_frac'] = [[LP_frac]] + (len(prefix_arr) - 1) * [[]]
                            output_dictionary['CHAN'][component]['LP_frac_err'] = [[LP_frac_err]] + (len(prefix_arr) - 1) * [[]]
                            output_dictionary['CHAN'][component]['LP_EVPA'] = [[LP_EVPA]] + (len(prefix_arr) - 1) * [[]]
                            output_dictionary['CHAN'][component]['LP_EVPA_err'] = [[LP_EVPA_err]] + (len(prefix_arr) - 1) * [[]]      
                
                        else:

                            output_dictionary['CHAN'][component]['Q_flux_mJy'][k].append(flux_Q)         
                            output_dictionary['CHAN'][component]['U_flux_mJy'][k].append(flux_U)               
                            output_dictionary['CHAN'][component]['V_flux_mJy'][k].append(flux_V)               
                            output_dictionary['CHAN'][component]['P_flux_mJy'][k].append(flux_P)           
                            output_dictionary['CHAN'][component]['P0_flux_mJy'][k].append(flux_P0)    
    
                            output_dictionary['CHAN'][component]['Q_err_mJy'][k].append(err_Q)            
                            output_dictionary['CHAN'][component]['U_err_mJy'][k].append(err_U)               
                            output_dictionary['CHAN'][component]['V_err_mJy'][k].append(err_V)                
                            output_dictionary['CHAN'][component]['P_err_mJy'][k].append(err_P)           

                            output_dictionary['CHAN'][component]['Q_rms_mJy'][k].append(rms_Q)            
                            output_dictionary['CHAN'][component]['U_rms_mJy'][k].append(rms_U)                
                            output_dictionary['CHAN'][component]['V_rms_mJy'][k].append(rms_V)                
                            output_dictionary['CHAN'][component]['P_rms_mJy'][k].append(rms_P)        

                            output_dictionary['CHAN'][component]['P_RA_deg'][k].append(RA_P)
                            output_dictionary['CHAN'][component]['P_DEC_deg'][k].append(DEC_P)    
        
                            output_dictionary['CHAN'][component]['LP_frac'][k].append(LP_frac)
                            output_dictionary['CHAN'][component]['LP_frac_err'][k].append(LP_frac_err)
                            output_dictionary['CHAN'][component]['LP_EVPA'][k].append(LP_EVPA)
                            output_dictionary['CHAN'][component]['LP_EVPA_err'][k].append(LP_EVPA_err)     

            except:
                msg('Fitting Failed: Channel is likely flagged')


        # Write RM Synthesis files (freq, I, Q, U, dI, dQ, dU)
        if not only_intensity:
            for component in components:
                rmsynth_arr = np.array([np.array(output_dictionary['CHAN'][component]['freq_GHz'][k]) * 1e9, 
                                                    np.array(output_dictionary['CHAN'][component]['I_flux_mJy'][k]) / 1e3, 
                                                    np.array(output_dictionary['CHAN'][component]['Q_flux_mJy'][k]) / 1e3, 
                                                    np.array(output_dictionary['CHAN'][component]['U_flux_mJy'][k]) / 1e3,
                                                    np.array(output_dictionary['CHAN'][component]['I_rms_mJy'][k]) / 1e3,
                                                    np.array(output_dictionary['CHAN'][component]['Q_rms_mJy'][k]) / 1e3,
                                                    np.array(output_dictionary['CHAN'][component]['U_rms_mJy'][k]) / 1e3])

                np.savetxt('{}_{}_rmsynth.txt'.format(prefix, component).replace(image_directory, 'RESULTS'), rmsynth_arr.T)

    # Save the full dictionary with all times, etc.
    if len(prefix.split(f'{image_identifier}-t')) > 1:
        file_name = prefix.split(f'{image_identifier}-t')[0] + image_identifier
    else:
        file_name = prefix
    with open('{}_polarization.json'.format(file_name).replace(image_directory, 'RESULTS'), 'w') as j:
        json.dump(output_dictionary, j, indent = 4)

    return 0

def main():
    
    # Load in the rmsynthesis data
    with open(cfg.DATA + '/rmsynth/rmsynth_info.json', 'r') as j:
        rmsynth_info = json.load(j)

    # Check to see if there is a Polarization angle calibrator -- pol_flag = Trye means that you do have
    pol_flag=False
    if cfg.POLANG_NAME != '':
        pol_flag = True
        
    # Iterate through sources as specified in rmsynth_info.json
    for k in range(len(rmsynth_info['image_directory']))[:]:
   
        if not os.path.exists(rmsynth_info['image_directory'][k]):
            msg('Skipping {} as the directory does not exist'.format(rmsynth_info['image_directory'][k]))
            continue

        # Construct image identifier based on input options
        if rmsynth_info["image_timing"][k]:
            src_im_identifier = cfg.CWD +'/{}/*{}*.ms_{}-t'.format(rmsynth_info["image_directory"][k], 
                                                                                            rmsynth_info["source_name"][k], 
                                                                                            rmsynth_info["image_identifier"][k])
        else:
            src_im_identifier = cfg.CWD +'/{}/*{}*.ms_{}'.format(rmsynth_info["image_directory"][k], 
                                                                                            rmsynth_info["source_name"][k], 
                                                                                           rmsynth_info["image_identifier"][k])                                                                                           
        # Separate out source positions
        src_ra   = []
        src_dec = []
        for pos in rmsynth_info['source_pos'][k]:
        
            # If the positions are given in CASA h:m:s d.m.s format covert to degrees
            if ':' in pos:
                pos_i = pos.split(',')
                pos_i[0], pos_i[1] = pos_i[0].split(':'), pos_i[1].split('.')
                
                if len(pos_i[1]) > 3:
                    pos_i[1][2] = pos_i[1][2] + '.' + pos_i[1][3]
                    pos_i[1] = pos_i[1][:3]
                
                sign = 1.0
                if '-' in pos_i[1][0]:
                    sign = -1.0
                
                pos_i[0] = float(pos_i[0][0]) * 15 + float(pos_i[0][1]) / 4. + float(pos_i[0][2]) / 240.
                pos_i[1] = (abs(float(pos_i[1][0])) + float(pos_i[1][1]) / 60.0 + float(pos_i[1][2]) / 3600.0) * sign

            # Else just remove 'deg' string (if there) and convert to floats
            else:
                pos_i.replace('deg', '').split(',')
            
            src_ra.append(pos_i[0])
            src_dec.append(pos_i[1])      
            
        src_ra, src_deg = np.array(src_ra), np.array(src_dec)  

        ###########################
        # Add the position prediction from here #
        ###########################
        
        msg(f'Starting property extraction for identifier: {src_im_identifier}')
        extract_polarization_properties(rmsynth_info["source_name"][k],
            src_im_identifier,
            rmsynth_info['image_suffix'][k], 
            src_ra, 
            src_dec, 
            rmsynth_info['source_ulim'][k],
            pol_flag, 
            rmsynth_info['rms_region'][k],
            rmsynth_info['image_directory'][k],
            rmsynth_info["image_identifier"][k],
            fix_additional_comps = True)



if __name__  == "__main__":
    main()
