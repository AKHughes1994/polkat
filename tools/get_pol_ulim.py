import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.special import gammaincinv, gammainc, iv
from scipy.integrate import trapz

fracpols = []
iPeaks = []


def main():

    # Read in the Stokes I image 
    if len(sys.argv) < 2:
        print('Please specify model FITS image and a target field name')
        sys.exit() 
    im_I = sys.argv[-1]
    im_P = im_I.replace('-I-','-P-')
    im_Q = im_I.replace('-I-','-Q-')
    im_U = im_I.replace('-I-','-U-')

    # Source Coorindates
    src_coord = '17:00:58.46,-46.11.08.6' # CASA format

    # Read in beam information from images
    bmaj_I = imhead(im_I, mode='get',hdkey='bmaj')['value']
    bmin_I = imhead(im_I, mode='get',hdkey='bmin')['value']    
    bpa_I  = imhead(im_I, mode='get',hdkey='bpa')['value']        

    bmaj_P = imhead(im_P, mode='get',hdkey='bmaj')['value']
    bmin_P = imhead(im_P, mode='get',hdkey='bmin')['value']    
    bpa_P  = imhead(im_P, mode='get',hdkey='bpa')['value']    

    src_region = f'circle[[{src_coord}], {3 * bmaj_I}arcsec]'    

    # Get prelim stokes I fit
    iFit = imfit(im_I, region = src_region)
    iPeak = iFit['results']['component0']['peak']['value']
    x = iFit['results']['component0']['pixelcoords'][0]
    y = iFit['results']['component0']['pixelcoords'][1] 

    # Make Stokes I estimate file
    myfile = open('estimate_I.txt','w')
    myfile.write(f'{iPeak},{x},{y},{bmaj_I}arcsec,{bmin_I}arcsec,{bpa_I}deg, abp')
    myfile.close()

    # Get stokes I fit
    iFit = imfit(im_I, region = src_region, estimates='estimate_I.txt')
    iPeak = iFit['results']['component0']['peak']['value']
    x = iFit['results']['component0']['pixelcoords'][0]
    y = iFit['results']['component0']['pixelcoords'][1] 

    # Make estimate file for Stokes P (forced aperature)
    myfile = open('estimate_P.txt','w')
    myfile.write(f'0.0,{x},{y},{bmaj_P}arcsec,{bmin_P}arcsec,{bpa_P}deg, xyabp')
    myfile.close()

    # Extract Linear Polarization flux density using a force aperature
    pFit = imfit(im_P, region = src_region, estimates ='estimate_P.txt')
    pPeak = pFit['results']['component0']['peak']['value']

    # Extract RMS
    r = np.sqrt(100 * bmaj_I ** 2) 

    di = imstat(im_I.replace('image.fits','residual.fits'),region = f'circle[[{src_coord}],{r}arcsec]')['rms'][0]
    du = imstat(im_U.replace('image.fits','residual.fits'),region = f'circle[[{src_coord}],{r}arcsec]')['rms'][0]
    dq = imstat(im_Q.replace('image.fits','residual.fits'),region = f'circle[[{src_coord}],{r}arcsec]')['rms'][0]
    
    # Solve for error in polarization intensity note that dq ~ du for all observations
    dp = 0.5 * du + 0.5 * dq

    # Solve for 99.73% CI upper limit following Vaillancourt (2006)
    def F(p,p0):
           return p * np.exp(-(p**2 + p0**2)/2) * iv(0,p*p0)

    def N(p):
            return p * np.sqrt(np.pi/2.) * np.exp(-p**2/4) * iv(0,p**2/4)

    for ulim in np.linspace(0.2,10,500):
            p0_arr = np.linspace(0,ulim,100000)
            frac = trapz(F(pPeak/dp,p0_arr)/N(pPeak/dp),x = p0_arr) * 100.0
            if frac >= 99.73:
                p0_ulim = ulim * dp
                print(f'Fraction: {frac:.2f}')
                break

    # Print Upper Limit in polarization fraction
    print(f'Stokes I: {iPeak * 1e3} +/- {di * 1e3} mJy')
    print(f'Stokes P (3-sigma Upper limit): {p0_ulim * 1e3} +/- {dp * 1e3} mJy')
    print(f'Polarization Fraction (3-sigma Upper limit): {p0_ulim/iPeak * 100}%')

if __name__ == "__main__":
    main()
