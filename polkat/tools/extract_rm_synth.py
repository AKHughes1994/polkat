import glob,sys, time, argparse
import numpy as np
from astropy.time import Time

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

def fit_cube(image_I, region, output=False):

    '''
    Fit the Stokes IQUV cube for an image. This assumes
    that the Q, U, and V images follow the WSCLEAN naming
    convention

    input parameters:
        images_I = Path the Stokes I image
        region   = Fiting region (should be cenetered around the source)
        output (optional) = Output the Polarization fluxes, fraction, and angle
    ''' 

    # Get the names from the Q, U and V image names from the input I image
    image_Q = image_I.replace('-I-', '-Q-')
    image_U = image_I.replace('-I-', '-U-')
    image_V = image_I.replace('-I-', '-V-')
    image_P = image_I.replace('-I-', '-P-') # "Stokes P" is the linear polarization intensity

    # Load in the beam parameters -- these should be ~constant across stokes parameters... if not... we have some problems  
    bmaj = imhead(image_I, mode='get', hdkey = 'BMAJ')['value']
    bmin = imhead(image_I, mode='get', hdkey = 'BMIN')['value']
    bpa  = imhead(image_I, mode='get', hdkey = 'BPA')['value']
    freq = imhead(image_I, mode='get', hdkey = 'CRVAL3')['value']

    # Get inital guesses for the positions and fluxex from CASA task imstat (no fitting)
    imstat_I = imstat(image_I, region=region)
    imstat_P = imstat(image_P, region=region)
    flux_I_guess, x_guess, y_guess   = imstat_I['max'][0], imstat_I['maxpos'][0], imstat_I['maxpos'][1]
    flux_P_guess, Px_guess, Py_guess = imstat_P['max'][0], imstat_P['maxpos'][0], imstat_P['maxpos'][1]
    flux_Q_guess  = return_max(image_Q, region = 'circle[[{}pix, {}pix],{}arcsec]'.format(Px_guess, Py_guess, 0.25 * bmaj))
    flux_U_guess  = return_max(image_U, region = 'circle[[{}pix, {}pix],{}arcsec]'.format(Px_guess, Py_guess, 0.25 * bmaj))
    flux_V_guess  = return_max(image_V, region = 'circle[[{}pix, {}pix],{}arcsec]'.format(x_guess, y_guess, 0.25 * bmaj))

    # Define the rms extraction region -- by default has an area of ~100 beams centered around peak pixel of source
    # Should be ~20uJy/beam for a 15-min L-band observation
    # May vary for bright sources, if it varies, tweak the region
    r_in  = 2.0 * bmaj
    r_out = np.sqrt(100 * 0.25 * bmaj * bmin + r_in ** 2)
    rms_region = 'annulus[[%spix,%spix],[%sarcsec,%sarcsec]]' %(x_guess,y_guess,r_in,r_out)

    # Extract the rms for each image
    rms_I = imstat(image_I, region=rms_region)['rms'][0]
    rms_Q = imstat(image_Q, region=rms_region)['rms'][0]
    rms_U = imstat(image_U, region=rms_region)['rms'][0]
    rms_V = imstat(image_V, region=rms_region)['rms'][0]
    rms_P = 0.5 * (rms_Q + rms_U) # This is the standard RMS assumption for Stokes P (could change to a "max" statement if conservative)

    # Make the estimate files
    f = open('estimate_I.txt', 'w')
    f.write('{},{},{},{}arcsec,{}arcsec,{}deg, abp'.format(flux_I_guess, x_guess, y_guess, bmaj, bmin, bpa))
    f.close()    
    
    f = open('estimate_P.txt', 'w')
    f.write('{},{},{},{}arcsec,{}arcsec,{}deg, abp'.format(flux_P_guess, Px_guess, Py_guess, bmaj, bmin, bpa))
    f.close()
 
    f = open('estimate_Q.txt', 'w')
    f.write('{},{},{},{}arcsec,{}arcsec,{}deg, xyabp'.format(flux_Q_guess, Px_guess, Py_guess, bmaj, bmin, bpa))
    f.close()

    f = open('estimate_U.txt', 'w')
    f.write('{},{},{},{}arcsec,{}arcsec,{}deg, xyabp'.format(flux_U_guess, Px_guess, Py_guess, bmaj, bmin, bpa))
    f.close()

    f = open('estimate_V.txt', 'w')
    f.write('{},{},{},{}arcsec,{}arcsec,{}deg, xyabp'.format(flux_V_guess, x_guess, y_guess, bmaj, bmin, bpa))
    f.close()

    # Fit each image
    fit_I = imfit(image_I, region = region, estimates='estimate_I.txt')
    fit_P = imfit(image_P, region = region, estimates='estimate_P.txt')  
    fit_Q = imfit(image_Q, region = region, estimates='estimate_Q.txt')  
    fit_U = imfit(image_U, region = region, estimates='estimate_U.txt')  
    fit_V = imfit(image_V, region = region, estimates='estimate_V.txt')  
    
    # Get the fluxes
    flux_I = fit_I['results']['component0']['peak']['value'] 
    flux_P = fit_P['results']['component0']['peak']['value']
    flux_Q = fit_Q['results']['component0']['peak']['value']
    flux_U = fit_U['results']['component0']['peak']['value']
    flux_V = fit_V['results']['component0']['peak']['value']
    
    # Calculate the Pol. Fracs and angles
    lin_pol    = np.sqrt(flux_P ** 2 - 2.3 * rms_P **2) / flux_I * 100. # This correction is standard (adopted from George, Stil, and Keller, 2012)
    lin_pol_err    = lin_pol * np.sqrt((rms_I/flux_I) ** 2 + (rms_P/flux_P) ** 2)
 
    EVPA       = 0.5 * np.arctan(flux_U / flux_Q) * 180.0 / np.pi
    EVPA_err       = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi

    circ_pol    = flux_V / flux_I * 100.
    circ_pol_err    = circ_pol * np.sqrt((rms_I/flux_I) ** 2 + (rms_V/flux_V) ** 2)

    # If ouput is desired print the ouput: 
    if output:
        msg('The Stokes I Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_I, rms_I))
        msg('The Stokes Q Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_Q, rms_Q))
        msg('The Stokes U Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_U, rms_U))
        msg('The Stokes V Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_V, rms_V))
        msg('The Lin. Pol. Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_P, rms_P))
        msg('Corresponding to a Lin. Pol. Frac of {:.5f} +/- {:.5f}% and an EVPA of {:.3f} +/- {:.3f} deg'.format(lin_pol, lin_pol_err, EVPA, EVPA_err))
        msg('Corresponding to a Circ. Pol. Frac of {:.5f} +/- {:.5f}%\n'.format(circ_pol, circ_pol_err))
    

    return [freq,flux_I,flux_Q,flux_U,rms_I,rms_Q,rms_U]
    
 

def main():

    # Load in the arguments 
    parser = argparse.ArgumentParser(description='Input parameters for fitting Stokes I, Q, and U for RM Synthesis fitting')    
    parser.add_argument('--im', '-i', help='model image prefix')
    parser.add_argument('--prefix', '-p', help='prefix for ouput file')
    parser.add_argument('--region', '-r',  help='approximate coorindates (in CASA format, e.g., "circle[[00:00:00.0,00.00.00.0],15pix]") for the source being fit')
    args = parser.parse_args()
    im   = str(args.im)
    prefix   = str(args.prefix)
    region = str(args.region)

    # This array will contain the full RM Synthesis data txt file
    data = []

    # First fit the MFS image to get an estimate of the polarization properties
    image_I = glob.glob(im + '*MFS-I-image.fits')[0]
    print('\n')
    msg('MFS image parameters:')
    fit_cube(image_I, region, output=True)

    # Also get the date from the MFS header file
    time = imhead(image_I, mode='get', hdkey='DATE-OBS')
    time_isot = time.replace('/','-',2).replace('/','T')
    time_mjd = Time(time_isot, format='isot').mjd

    # Now iterate through channelized images
    for image_I in sorted(glob.glob(im + '*-00*-I-image.fits')):

        # The Try statment is so that flagged channel fitting will fail without crashing the script
        try:    
            msg('Fitting channelized image: ' +  image_I)
            data_i = fit_cube(image_I, region, output=False)
            data.append(data_i)
 
        except:
            msg('Channelized image fitting failed (likely flagged): ' + image_I)

    np.savetxt('{}_{}.txt'.format(prefix, time_mjd), data)

if __name__ in "__main__":
    main()



