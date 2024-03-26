--- 

### What is this?

This modified version of the semi-automated routine [oxkat](https://github.com/IanHeywood/oxkat) includes polarisation calibration and imaging. This assumes you are familiar with the oxkat workflow:


---
##### Change of workflow

* `INFO` -- 

* There are a lot of similarities between the two containers; however, `polkat-0.0.1.sif` now uses [`casa`](https://casa.nrao.edu/) version 6.x., and I have the intention of (eventually) converting all of the [`cubical`](https://github.com/ratt-ru/CubiCal) scripts to use the updates [`quartical`](https://github.com/ratt-ru/CubiCal) package. 


---
##### Major changes with respect to oxkat

* The ms file is averaged to 1024 channels in the INFO step (oxkat does this in 1GC). This involves getting the numeric indexes of the (sub-)fields used in each calibration run. INFO also identifies the polarization calibrator automatically—currently, MeerKAT defaults to specifying the scan intent as "UNKNOWN" for the polarization calibrator.

* By default, any `wsclean` call uses its auto-threshold and auto-masking routine (alongside manual masking). The auto-routines use Gaussian kernels when looking for clean components and thus will often find single-pixel clean components for point sources. Anecdotally, I find that this results in better self-calibration, as standard thresholding will clean every pixel within the manual mask (but the difference is often non-specific). This should also help for messy fields where the RMS is well above the systematic threshold of the instrument due to, for example, diffuse emission!

* 1GC is the most significant change to the calibration routine. The 1GC calibration now solves for complete polarization solutions following the below workflow:
  1. Solve for K (Parr-hand Delay), Gp (Gain-phase), Ga (Gain-amplitude), Bp (Bandpass), and Df (Leakage) solutions using the (unpolarized) primary (this is typically J1939)
  2. Solve for K, Gp, Ga, Bp, and Df for the polarization angle and secondary calibrators -- **NOTE** all Gain amplitude solutions adopt a gain type of 'T', which means that there is one solution for both antenna-feed polarizations. Long story short, this preserves the polarization information for polarized calibrators, as CASA, by default, will assume an unpolarized calibrator for gain type 'G.' 
  3. Solve for KCROSS (Cross-hand delay) and Xf (cross-hand phase) using the polarization angle calibrator. Solving for Xf for a linear-feed system is (pseudo)-model-independent. The only requirement is that the model for the polarization angle calibrator be initialized such that the flux has U > Q and a V = 0.0. By default, this is done with `setjy` with the parameter `fluxdensity=[1.0,0.0,0.5,0.0].`
  6. At the end of 1GC, by default, you will image the (unpolarized) primary and the (known polarization) polarization angle calibrator. This can be turned off in `config.py` (setting `CAL_1GC_DIAGNOSTICS = False`). I would suggest against turning this off, as these act as built-in check sources to probe the accuracy of our calibration. Any excess polarization in the primary can be used as a systematic error (probably caused by residual leakage), and the polarization calibrator can have its angle/fraction(s) checked against their known values as a further systematic. I'll talk about this more later on.
 
* FLAG and 2GC have been combined into a single step. The major differences:
  1. "datamask" and "pcalmask" images are full stokes IQUV.
  2. The routine no longer splits out the calibrated source(s) into separate ms file(s). The reason is that, for linear feeds, the parallactic angle corrections due not commute. More simply, What this means is that if you perform full polarization calibration (with parallactic angle corrections) and then split out the data for self-calibrating, you will "over-correct" the parallactic angle. Further gain solutions will be incorrect as you are now solving for the values against the rotated visibilities. 
  3. As a result of 2. the routine is now self-calibration using `casa` rather than `cubical.` I am investigating whether we can use Quartical for full-polarization self-calibration, but this is a large work in progress. 
 
* 3GC has yet to work; don't try it. 
 
### VERY IMPORTANT: CHECK YOUR VISIBILITES + GAIN TABLES BEFORE MOVING ON FROM 1GC to 2GC; POLARIZATION CALIBRATION CAN DO WEIRD THINGS! WHAT YOU WANT TO ENSURE IS THAT THE STOKES V FOR THE POLARIZATION CALIBRATOR IS ~0.IN VISIBILITY SPACE, THIS CORRESPONDS TO THE IMAGINARY COMPONENT OF XY AND YX 

---
##### The New Setups

There are two new setups:

###### SNAP

This routine added the Heywood snapshot imaging routine to the general workflow. After 2GC you can run:
   ```
   $ python setups/SNAP.py idia
   $ ./submit_snap_job.sh
   ```
By default, the routine will make snapshot images of all Targets in the MS file. However, you can modify the `SNAP_FIELDS` variable inside of `config.py` if you want to specify the sources for snapshot imaging (e.g., want to snapshot image a calibrator as a check source)

###### RMSYNTH

This semi-automated routine should fit/extract the properties of an arbitrary number of Gaussian components for your IQUV imaging cube(s).
   ```
   $ python setups/RMSYNTH.py idia
   $ ./submit_rmsynth_job.sh
   ```
The only thing that needs to be modified is the file `data/rmsynth/rmsynth_info.txt` below is an example file:

```
# Text file containing information to feed into the RMSYNTH_01_extract_fluxes.py routine columns are:
# Field name, image identifier (for IQUV cube, e.g., "pcalmask"), ra, dec (pixels) seperated by spaces (can include multiple RA/DEC for one image)
Field1 pcalmask 17:27:43.3346781657,17:27:43.3378907659 -16.12.19.5120108691,-16.12.26.3558690263
Field2 datamask 18:00:00 -17:00:00
```

Ignoring the preamble, the code will look for the four columns to get the necessary information for fitting:
  1. The first column is the field name (as seen in the ms file)
  2. The second column is the image identifier; if you want to fit the self-called images, use 'pcalmask'
  3. The third/fourth column is the RA/Dec guess(es) in the standard CASA format (note the period separators for declination). These can be single-coordinates or comma-separated lists for multi-component fitting (e.g., if you have core + jet ejecta)

For this example, you will fit the Field1 self-calibrated image with a two-component fit, and the Field2 masked image (no self-cal) with a one-component fit.

The RMSYNTH step will also measure the polarization systematics using the calibrators, measure RM/pol angle using RM Synthesis, and get estimates of the ionospheric RM contributions using ALBUS! I will go into detail about the fitting approaches soon!

---
##### Example run

## Quick start

1. Navigate to a working area / scratch space:

   ```
   $ cd /scratch/users/emperor-zerg/
   ```

2. Clone the root contents of this repo into it:

   ```
   $ git clone -b master https://github.com/AKHughes1994/polkat.git .
   ```

3. Make a symlink to your MeerKAT Measurement Set (or place it in the working folder, it will not be modified at all):

   ```
   $ ln -s /idia/raw/xkat/SCI-20230907-RF-01/1708829174/1708829174_sdp_l0.ms .
   ```

4. The first step is to run a script that gathers some required information about the observation:

   ```
   $ python setups/0_GET_INFO.py idia
   $ ./submit_info_job.sh
   ```

5. Once this is complete, you can generate and submit the jobs required for the reference calibration (1GC):

   ```
   $ python setups/1GC.py idia
   $ ./submit_1GC_jobs.sh
   ```

6. Once this is complete, you can generate and submit the jobs required for the target imaging and self-calibration (2GC):

   ```
   $ python setups/2GC.py idia
   $ ./submit_2GC_jobs.sh
   ```

6. If something goes wrong, you can kill the running and queued jobs on a cluster with, e.g.,:

   ```
   $ source SCRIPTS/kill_1GC_jobs.sh
   ```

7. Once all the jobs have been completed, you can examine the products and move on to the setup for the next steps, namely RMSYNTH and SNAP for data analysis.

Please see the [setups README](setups/README.md) for more details about the general workflow. Most settings can be tuned via the [`config.py`](oxkat/config.py) file.

---
##### To-Do List

1. Figure out if there is a significant difference between the `-mf-weighting` and `-no-mf-weighting` channelized images
1. Investigate full self-calibration polarisation routines.
3. Add peeling capabilities

