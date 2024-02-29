--- 

### What is this?

This is a modified version of the semi-automated routine [oxkat](https://github.com/IanHeywood/oxkat) designed to make polarization calibration as hands-off as possible. I direct the reader to the original oxkat documentation for a more detailed pipeline description. Here, I will highlight some key differences.


---
##### Overarching Changes

* This uses a mix of containers, namely `oxkat-0.41.sif`, `polkat-0.0.2.sif`, `polkat-albus.sif`. These are all contained on the idia cluster (`/idia/software/containers`), and, as a result, if you want to use polkat on other clusters (or locally) you have to download the container and put it in the appropriate directories.

* There are a lot of similarities between the two containers; however, `polkat-0.0.1.sif` now uses [`casa`](https://casa.nrao.edu/) version 6.x., and I have the intention of (eventually) converting all of the [`cubical`](https://github.com/ratt-ru/CubiCal) scripts to use the updates [`quartical`](https://github.com/ratt-ru/CubiCal) package. 


---
##### Major changes with respect to oxkat

* The ms file is averaged to 1024 channels in the INFO step (oxkat does this in 1GC). This involves getting the numeric indexes of the (sub-)fields used in each calibration run. INFO also identifies the polarization calibrator automatically—currently, MeerKAT defaults to specifying the scan intent as "UNKNOWN" for the polarization calibrator.

* By default, `clean` uses its auto-threshold and auto-masking routine (alongside manual masking). The auto-routines use Gaussian kernels when looking for clean components and thus will often find single-pixel clean components for point sources. Anecdotally, I find that this results in better self-calibration, as standard thresholding will clean every pixel within the manual mask (but the difference is often non-specific). This should also help for messy fields where the RMS is well above the systematic threshold of the instrument due to, for example, diffuse emission!

* 1GC is the most significant change to the calibration routine. The 1GC calibration now solves for complete polarization solutions following the below workflow:
  1. Solve for K (Parr-hand Delay), Gp (Gain-phase), Ga (Gain-amplitude), Bp (Bandpass), and Df (Leakage) solutions using the (unpolarized) primary (this is typically J1939)
  2. Solve for K, Gp, Ga, Bp, and Df for the polarization angle and secondary calibrators -- **NOTE** all Gain amplitude solutions adopt a gain type of 'T', which means that there is one solution for both antenna-feed polarizations. Long story short, this preserves the polarization information for polarized calibrators, as CASA, by default, will assume an unpolarized calibrator for gain type 'G.' 
  3. Solve for KCROSS (Cross-hand delay) and Xf (cross-hand phase) using the polarization angle calibrator. Solving for Xf for a linear-feed system is (pseudo)-model-independent. The only requirement is that the model for the polarization angle calibrator be initialized such that the flux has U > Q and a V = 0.0. By default, this is done with `setjy` with the parameter `fluxdensity=[1.0,0.0,0.5,0.0].`
  6. At the end of 1GC, by default, you will image the (unpolarized) primary and the (known polarization) polarization angle calibrator. This can be turned off in `config.py` (setting `CAL_1GC_DIAGNOSTICS = False`). I would suggest against turning this off, as these act as built-in check sources to probe the accuracy of our calibration. Any excess polarization in the primary can be used as a systematic error (probably caused by residual leakage), and the polarization calibrator can have its angle/fraction(s) checked against their known values as a further systematic. I'll talk about this more later on.
 
* FLAG and 2GC have been combined into a single step. The major differences:
  1. "datamask" and "pcalmask" images are full stokes IQUV.
  2. The routine no longer splits out the calibrated source(s) into separate ms file(s). The reason is that, for linear feeds, the parallactic angle corrections due not commute. More simply, What this means is that if you perform full polarization calibration (with parallactic angle corrections) and then split out the data for self-calibrating, you will "over-correct" the parallactic angle. Further gain solutions will be incorrect as you are now solving for the values against the rotated visibilities. 
  3. As a result of 2. the routine is now self-calibration using `casa` rather than `cubical.` I am investigating whether we can use Quartical for full-polarization self-calibration, but this is a large work in progress. 

* There are two new setups, namely RMSYNTH and SNAP.
  1. SNAP added the Heywood snapshot imaging routine into the standard workflow (only Stokes-I currently)
  2. RMSYNTH (pseudo-)automatically extracts full IQUV images from the source of interest
 
* 3GC has yet to work; don't try it. 
 
[AD MORE INFO ON SNAP AND RMSYNTH]


### VERY IMPORTANT: CHECK YOUR VISIBILITES + GAIN TABLES BEFORE MOVING ON FROM 1GC to 2GC; POLARIZATION CALIBRATION CAN DO WEIRD THINGS! WHAT YOU WANT TO ENSURE IS THAT THE STOKES V FOR THE POLARIZATION CALIBRATOR IS ~0.IN VISIBILITY SPACE, THIS CORRESPONDS TO THE IMAGINARY COMPONENT OF XY AND YX 

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

1. Investigate full self-calibration polarisation routines.
2. Add peeling capabilities

