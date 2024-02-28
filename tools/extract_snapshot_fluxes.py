import numpy as np
import glob
import sys
import json
import os
import os.path as o
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from argparse import ArgumentParser

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)

def fit_image(image, pos):

    # Get the image parameters
    bmaj  = imhead(image, mode='get', hdkey = 'BMAJ')['value']
    bmin  = imhead(image, mode='get', hdkey = 'BMIN')['value']
    bpa    = imhead(image, mode='get', hdkey = 'BPA')['value']
    freq_GHz  = imhead(image, mode='get', hdkey = 'CRVAL3')['value'] / 1.0e9
    date_obs  =  imhead(image, mode='get', hdkey = 'DATE-OBS').replace('/','-',2).replace('/','T')

    # Get and initial imstat fit for the peak position and flux as well as the rms
    r_in  = 3.0 * bmaj
    r_out = np.sqrt(500 * 0.25 * bmaj * bmin + r_in ** 2)
    
    # Define regions
    region         = f'circle[[{pos}], {r_in}arcsec]'
    rms_region = f'annulus[[{pos}], [{r_in}arcsec, {r_out}arcsec]]'
   
    # Get rms
    rms = imstat(image, region = rms_region)['rms'][0] * 1e3

    # Run imstat
    ims = imstat(image, region=region)

    # Get the intial guesses
    flux_guess = ims['max'][0]
    x_guess = ims['maxpos'][0]
    y_guess = ims['maxpos'][1]

    # Make an estimate file
    f = open('estimate.txt', 'w')
    f.write(f'{flux_guess},{x_guess},{y_guess},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, abp')
    f.close()
    
    # Fit source using imfit
    imf = imfit(image, region=region, estimates='estimate.txt')

    # Seperate the imfit values
    flux = imf['results']['component0']['peak']['value'] * 1e3
    ra    = imf['results']['component0']['shape']['direction']['m0']['value'] * 180 / np.pi
    dec = imf['results']['component0']['shape']['direction']['m1']['value'] * 180 / np.pi

    return [flux, rms, ra, dec, bmaj, bmin, bpa, freq_GHz, date_obs]



def main():

    # Read in the arguments
    parser = ArgumentParser(description='CASA script to iteratively solve for snapshot light cruves of a (very) bright source, this script assumes the standard polkat directory structure + naming convention')
    parser.add_argument('targetname',
                        help="target field to calibrate")
    parser.add_argument('-c', '--coordinates', dest="coordinates",
                        help="string containing the (comma-seperated) RA and dec estimate in a CASA compatible format")
    parser.add_argument('-p', '--plot', dest="plot", action='store_true', default=False,
                        help="boolean flag to determine whether the routine will output diagnostic plots")

    args = parser.parse_args()
    targetname = str(args.targetname)
    plot = args.plot
    coordinates = str(args.coordinates)

    # Read in the INTERVAL images from the source of interest
    restored_images = sorted(glob.glob(cfg.INTERVALS + f'/*{targetname}*image-restored.fits'))

    # Initialize storage directoyr:
    fit_dict = {'flux_mjy':[], 'time_isot':[], 'freq_GHz':[], 'rms_mjy':[], 'ra_deg':[], 'dec_deg':[], 'bmaj_arcsec':[], 'bmin_arcsec':[], 'bpa_deg':[]}

    # Iterate through images getting fluxes
    for restored_image in restored_images[:]:
        try:
            msg('Fitting image: ' + restored_image.split(cfg.INTERVALS + '/')[-1])
            data_arr = fit_image(restored_image, coordinates)
            fit_dict['flux_mjy'].append(data_arr[0])
            fit_dict['rms_mjy'].append(data_arr[1])
            fit_dict['ra_deg'].append(data_arr[2])
            fit_dict['dec_deg'].append(data_arr[3])
            fit_dict['bmaj_arcsec'].append(data_arr[4])
            fit_dict['bmin_arcsec'].append(data_arr[5])
            fit_dict['bpa_deg'].append(data_arr[6])
            fit_dict['freq_GHz'].append(data_arr[7])
            fit_dict['time_isot'].append(data_arr[8])
        except:
            msg('Fitting Failed: Skipping')


    # Sort dictionary by time
    index = np.argsort(fit_dict['time_isot'])
    for key in fit_dict.keys():
         fit_dict[key] = (np.array(fit_dict[key])[index]).tolist()

    # Save the dictionary
    json.dump(fit_dict, open(cfg.RESULTS+ f'/{targetname}_snapshot.json', "w"), indent=4)


    if plot:
    
        # Caclulate the weighted average + error on the average of the flux
        flux_avg        = np.average(fit_dict['flux_mjy'], weights = np.array(fit_dict['rms_mjy']) ** (-2))
        flux_avg_err = np.sum(np.array(fit_dict['rms_mjy']) ** (-2)) ** (-1)

        print(flux_avg, flux_avg_err)

        lw = 2 #2.5  # border width (pixels)
        majorw = 2.5; majorl = 12.5 #major tick width and length (pixels)
        minorw = 2.5; minorl = 6.25 #minor tick width and length (pixels)

        mpl.rcParams['font.size'] = 20
        mpl.rcParams['font.family'] = 'serif'
        mpl.rcParams['figure.facecolor'] = 'white'

        indexes = np.arange(len(fit_dict['flux_mjy']))

        fig, ax = plt.subplots(2, figsize=(12,24), sharex=True, gridspec_kw={'hspace': 0.03})

        ax[0].set_title('Source: ' + targetname)

        # Flux
        ax[0].errorbar(indexes, fit_dict['flux_mjy'], fit_dict['rms_mjy'] , fmt='o--', c='k') 
        ax[0].axhline(flux_avg, c ='k', ls=':')
        ax[0].axhspan(flux_avg - flux_avg_err, flux_avg + flux_avg_err, color='k', alpha=0.25)
        ax[0].set_ylabel('Flux Density\n(mJy/beam)')

        # Positional offset
        ra_offset = (np.array(fit_dict['ra_deg']) - np.mean(fit_dict['ra_deg'])) * 3600.0
        dec_offset = (np.array(fit_dict['dec_deg']) - np.mean(fit_dict['dec_deg'])) * 3600.0
        ax[1].errorbar(indexes, ra_offset, yerr = np.array(fit_dict['bmaj_arcsec']) * 0.1, c='r', fmt='o--', label='RA')
        ax[1].errorbar(indexes, dec_offset, yerr = np.array(fit_dict['bmaj_arcsec']) * 0.1, c='b', fmt='o--', label='Dec')
        ax[1].axhline(0.0, c ='k', ls=':')
        ax[1].set_ylabel('Positional Offset\n(arcsec)')
        ax[1].set_xlabel('Image # (Arbitrary)')

        # Make the plot look pretty
        for ax_i in ax:
            ax_i.minorticks_on()
            ax_i.tick_params(axis='both', which='major', direction='in',length=majorl,width=majorw,top=True,right=False)
            ax_i.tick_params(axis='both', which='minor', direction='in',length=minorl,width=minorw,top=True,right=False)

            for axis in ['top','bottom','left','right']: 
                ax_i.spines[axis].set_linewidth(lw)       #set border thickness

        # Save figure
        plt.savefig(f'{targetname}_snapshot.png')
        plt.close()

if __name__ == "__main__":
    main()

        
        
        
        

        
        
