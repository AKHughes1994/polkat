--- 

### What is this?

This modified version of the semi-automated routine [oxkat](https://github.com/IanHeywood/oxkat) includes polarisation calibration and imaging. This assumes you are familiar with the oxkat workflow


---
##### Change of workflow when compared to Oxkat

* INFO — The ms file is averaged in this step; creates a `json` dictionary called "pre_fields" that contains index-to-name mapping
* 1GC  — Calibration now included leakage and cross-hand phase calibration
* 2GC  — FLAG and 2GC have been combined into a single setup; all imaging/masking/self-cal is now in 2GC

---
##### Changes to `config.py` — new parameters

* `POLANG_NAME = 'J1331+3030'` — This is the source's name to be used as a polarization angle calibrator. The default is J1331+3030 (3C286). If blank, the routine will skip the polarization angle calibration. 
* `POLANG_DIR  = '13:31:08.2881,+30.30.32.959'` — Coordinates of polarization angle calibrator. The default is J1331+3030 (3C286) coorindates
* `PRE_FIELDS = ''` — If this is not blank (''), must contain `POLANG_NAME`
* `POLANG_MOD = [1.0, 0.0, 0.5, 0.0]` — an array to initialize a (quasi-)arbitrary polarization model. The default works for 3C286.
* `UNIFORM_IMAGE = True` — Flag to make an additional uniform weighted image. Useful if you want high angular resolution for ejecta and high sensitivity for polarisation. Images are made by default.
* `SNAP_FIELDS = ''` — Name of field(s) to perform Heywood snapshot imaging routine.
* `RED_TYPE = RI_G03` — Type of GPS fitting when determining ionospheric RM using [ALBUS](https://github.com/twillis449/ALBUS_ionosphere/blob/master/python_scripts/MS_Iono_functions.py).
*  `CAL_1GC_DIAGNOSTICS = True` — Flag to determine if you will image/fit leakage and polarisation angle calibrator to quantify systematics. Just leave this turned on. 


---
##### Standard Workflow

If you are working on IDIA/ILIFU, here would be a standard workflow from start to finish:

###### INFO

Get ms info, and average channels (to 1024 by default)

   ```
   $ python setups/INFO.py idia
   $ ./submit_info_job.sh
   ```

###### 1GC

After INFO is complete, you can perform reference 1GC calibration using calibrator fields.

   ```
   $ python setups/1GC.py idia
   $ ./submit_1GC_job.sh
   ```

please inspect the visibility/gain solutions as polarization can be finicky

###### 2GC

After 1GC is complete, perform final flagging, imaging, and direction-independent phase self-calibration

   ```
   $ python setups/2GC.py idia
   $ ./submit_2GC_job.sh
   ```

At this point, you should have all the imaging data products. Here, if you want to run snapshot imaging on a field and have set the `SNAP_FIELDS` parameter in `config.py` you can run, 

   ```
   $ python setups/SNAP.py idia
   $ ./submit_snap_job.sh
   ```

###### RMSYNTH

This is a new setup that will fit the full Stokes I, Q, U, and V images, before running it open the images with the 'Plin' or 'Ptot' identifiers and see if there is a polarization detection.  If there is you can run:

   ```
   $ python setups/RMSYNTH.py idia
   $ ./submit_rmsynth_job.sh
   ```

What this will do is: 

   * Fit the source with an arbitrary number of Gaussians using the `casa` task `imfit`
   * If `CAL_1GC_DIAGNOSTICS = True` (which it should be) it will quantify the systematic calibration effects through image plane analysis of the calibrators
   * Run RM Synthesis on every source/component extracting polarisation angles/rotation measures
   * Run ALBUS to get ionospheric RM for post-processing corrections. This randomly fails sometimes you can (i) run it again; (ii) Change `RED_TYPE = RI_G03` to `RED_TYPE = RI_G01`, the former will fail if only one GPS station is operational

The only thing that needs to be modified is the file `data/rmsynth/rmsynth_info.txt`. Below is an example file:

```
# Text file containing information to feed into the RMSYNTH_01_extract_fluxes.py routine columns are:
# Field name, image identifier (for IQUV cube, e.g., "pcalmask"), ra, dec (pixels) separated by spaces (can include multiple RA/DEC for one image)
Field1 pcalmask 17:27:43.3346781657,17:27:43.3378907659 -16.12.19.5120108691,-16.12.26.3558690263
Field2 datamask 18:00:00 -17:00:00
```
Ignoring the preamble, the code will look for the four columns to get the necessary information for fitting:
  1. The first column is the field name (as seen in the ms file)
  2. The second column is the image identifier; if you want to fit the self-called images, use 'pcalmask'
  3. The third/fourth column is the RA/Dec guess(es) in the standard CASA format (note the period separators for declination). These can be single-coordinates or comma-separated lists for multi-component fitting (e.g. if you have core + jet ejecta)

For this example, you will fit the Field1 self-calibrated image with a two-component fit, and the Field2 masked image (no self-cal) with a one-component fit.

---

##### ThunderKAT (No Polarization Angle Calibration)

polkat will work without a polarization angle calibrator. As long as the primary if unpolarized leakage calibration will be sufficient to measure the total polarization, you will not be able to distinguish between circular and linearly polarized emission. However, given that XRBs are synchrotron sources, the VAST majority will be circularly unpolarized, and thus, the total polarization is a good measure of linear polarization in the absence of a polarization calibrator. 

Leave `POLANG_NAME = ''`, RMSYNTH will skip the RM Synthesis steps, and only quantify the systematics using the leakage calibrator

---
##### To-Do List

1. Investigate systematic offset of 3C286 properties vs. expectation
2. Investigate full self-calibration polarisation routines
3. Add peeling capabilities

