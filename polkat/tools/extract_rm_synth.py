import glob,sys
import numpy as np
from astropy.time import Time

im  = sys.argv[-1]
pos = sys.argv[-2]

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

print('\n\n')
data = []

# MFS image polarization fraction
image_I = glob.glob(im + '*MFS-I-image.fits')[0]
image_P = image_I.replace('-I-','-P-')    
image_Q = image_I.replace('-I-','-Q-')    
    
# MFS Images
bmaj = imhead(image_I, mode='get', hdkey = 'BMAJ')['value']
bmin = imhead(image_I, mode='get', hdkey = 'BMIN')['value']
bpa  = imhead(image_I, mode='get', hdkey = 'BPA')['value']

region = 'circle[[{}],{}arcsec]'.format(pos,1.5 * bmaj)

# Get initial estimate 
imstat_I = imstat(image_I, region=region)
imstat_P = imstat(image_P, region=region)
flux_I_guess, x_guess, y_guess   = imstat_I['max'][0], imstat_I['maxpos'][0], imstat_I['maxpos'][1]
flux_P_guess, Px_guess, Py_guess = imstat_P['max'][0], imstat_P['maxpos'][0], imstat_P['maxpos'][1]

# Make the rms region
r_in  = 2.0 * bmaj
r_out = np.sqrt(100 * 0.25 * bmaj * bmin + r_in ** 2)
rms_region = 'annulus[[%spix,%spix],[%sarcsec,%sarcsec]]' %(x_guess,y_guess,r_in,r_out)

# Get the RMS
rms_I = imstat(image_I, region=rms_region)['rms'][0]
rms_P = imstat(image_Q, region=rms_region)['rms'][0]

# Fit the fluxes
fit_I = imfit(image_I, region = region, estimates='estimate_I.txt')
fit_P = imfit(image_P, region = region, estimates='estimate_P.txt')  

# Get the fluxes
flux_I = fit_I['results']['component0']['peak']['value'] 
flux_P = fit_P['results']['component0']['peak']['value']

# Calculate the Pol. Fracs
lin_pol        = np.sqrt(flux_P ** 2 - 2.3 * rms_P **2) / flux_I * 100.       # This correction is standard (adopted from George, Stil, and Keller, 2012) 
lin_pol_err    = lin_pol * np.sqrt((rms_I/flux_I) ** 2 + (rms_P/flux_P) ** 2)

# Print Relevant Values
print('For frequency channel corresponding to Stokes I image:', image_I)
print('The Stokes I Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_I, rms_I))
print('The Lin. Pol. Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_P, rms_P))
print('Corresponding to a Lin. Pol. Frac of {:.5f} +/- {:.5f}%'.format(lin_pol, lin_pol_err))
print('\n'

# Iterate through all of the images
for image_I in sorted(glob.glob(im + '*-00*-I-image.fits')):
    image_Q = image_I.replace('-I-', '-Q-')
    image_U = image_I.replace('-I-', '-U-')
    image_V = image_I.replace('-I-', '-V-')
    image_P = image_I.replace('-I-', '-P-')

    try:    
        # Load in the global parameters from the header file
        bmaj = imhead(image_I, mode='get', hdkey = 'BMAJ')['value']
        bmin = imhead(image_I, mode='get', hdkey = 'BMIN')['value']
        bpa  = imhead(image_I, mode='get', hdkey = 'BPA')['value']
        freq = imhead(image_I, mode='get', hdkey='CRVAL3')['value']
        time = imhead(image_I, mode='get', hdkey='DATE-OBS')
        time_isot = time.replace('/','-',2).replace('/','T')
        time_mjd = Time(time_isot, format='isot').mjd
        
        # Get the extents along RA and Dec
        dec_extent = 2.0 * np.sqrt((bmaj * 0.5) ** 2 * np.cos(bpa) ** 2  + (bmin * 0.5) ** 2 * np.sin(bpa) ** 2)
        ra_extent  = 2.0 * np.sqrt((bmaj * 0.5) ** 2 * np.sin(bpa) ** 2  + (bmin * 0.5) ** 2 * np.cos(bpa) ** 2)
       
    
        # Get an initial postion/flux estimates
        imstat_I = imstat(image_I, region=region)
        imstat_P = imstat(image_P, region=region)
        flux_I_guess, x_guess, y_guess   = imstat_I['max'][0], imstat_I['maxpos'][0], imstat_I['maxpos'][1]
        flux_P_guess, Px_guess, Py_guess = imstat_P['max'][0], imstat_P['maxpos'][0], imstat_P['maxpos'][1]
        flux_Q_guess  = return_max(image_Q, region = 'circle[[{}pix, {}pix],{}arcsec]'.format(x_guess, y_guess, 0.25 * bmaj))
        flux_U_guess  = return_max(image_U, region = 'circle[[{}pix, {}pix],{}arcsec]'.format(x_guess, y_guess, 0.25 * bmaj))
        flux_V_guess  = return_max(image_V, region = 'circle[[{}pix, {}pix],{}arcsec]'.format(x_guess, y_guess, 0.25 * bmaj))
       
        # Make the rms region
        r_in  = 2.0 * bmaj
        r_out = np.sqrt(100 * 0.25 * bmaj * bmin + r_in ** 2)
        rms_region = 'annulus[[%spix,%spix],[%sarcsec,%sarcsec]]' %(x_guess,y_guess,r_in,r_out)
        
        # Extract the rms for each image
        rms_I = imstat(image_I, region=rms_region)['rms'][0]
        rms_Q = imstat(image_Q, region=rms_region)['rms'][0]
        rms_U = imstat(image_U, region=rms_region)['rms'][0]
        rms_V = imstat(image_V, region=rms_region)['rms'][0]
        rms_P = 0.5 * (rms_Q + rms_U)
        
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
    
        # Calculate the Pol. Fracs
        lin_pol        = np.sqrt(flux_P ** 2 - 2.3 * rms_P **2) / flux_I * 100.
        lin_pol_err    = lin_pol * np.sqrt((rms_I/flux_I) ** 2 + (rms_P/flux_P) ** 2)
 
        EVPA           = 0.5 * np.arctan(flux_U / flux_Q) * 180.0 / np.pi
        EVPA_err       = 0.5 * np.sqrt(flux_U ** 2 * rms_Q **2  + flux_Q ** 2 * rms_U ** 2) / (flux_U ** 2  + flux_Q ** 2) * 180.0 / np.pi

        circ_pol        = flux_V / flux_I * 100.
        circ_pol_err    = circ_pol * np.sqrt((rms_I/flux_I) ** 2 + (rms_V/flux_V) ** 2)

        data.append([freq,flux_I,flux_Q,flux_U,rms_I,rms_Q,rms_U])

        # Print Relevant Values
        print('For frequency channel corresponding to Stokes I image:', image_I)
        print('The Stokes I Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_I, rms_I))
        print('The Stokes Q Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_Q, rms_Q))
        print('The Stokes U Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_U, rms_U))
        print('The Stokes V Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_V, rms_V))
        print('The Lin. Pol. Flux density is {:.5f} +/ {:.5f} Jy/beam'.format(flux_P, rms_P))
        print('Corresponding to a Lin. Pol. Frac of {:.5f} +/- {:.5f}% and an EVPA of {:.3f} +/- {:.3f} deg'.format(lin_pol, lin_pol_err, EVPA, EVPA_err))
        print('Corresponding to a Circ. Pol. Frac of {:.5f} +/- {:.5f}%\n'.format(circ_pol, circ_pol_err))
    except:
        print('For frequency channel corresponding to Stokes I image:', image_I)
        print('FIT FAILED: Very likely the channel is flagged')

np.savetxt('{}_rmsynth_data.txt'.format(im), data)
