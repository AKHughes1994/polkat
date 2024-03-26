--- 

### What is this?

This modified version of the semi-automated routine [oxkat](https://github.com/IanHeywood/oxkat) includes polarisation calibration and imaging. This assumes you are familiar with the oxkat workflow:


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

