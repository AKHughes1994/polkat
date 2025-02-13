--- 

### What is this?

This modified version of the MeerKAT semi-automated data processing routine [oxkat](https://github.com/IanHeywood/oxkat) that was tweaked to include full polarisation calibration and for Stokes I, Q, U, V imaging. We assume you are familiar with the oxkat workflow and its file system and that data processing options are contained in and changed by editing `oxkat/config.py`. This guide will walk you through a standard use case and highlight some changes to the `config.py` file and new options. 


---
##### Before We Start 

This routine has been designed primarily for use on the ILIFU clusters operated by The Inter-university Institute for Data Intensive Astronomy (IDIA); however, you can run it locally if you have the comprising software. The necessary software has been combined into containers using [apptainer](https://apptainer.org/) (previously known as singularity). The containers can be found in the `/software/containers` directory on ILIFU and follow the naming convention of `polkat-[version].sif`. If you do not have access to ILIFU but still want to use polkat and thus require the container, please email 'hughesakh [at] gmail [dot] com', and a download link can be made available.

---
#### Standard Workflow

Henceforth, we will assume a Linux-based operating system (I use Ubuntu). 

You can begin by initializing a working directory, moving or symlinking your CASA measurement file (i.e., ms-file), and cloning this repository. For example:

```
mkdir working_directory
cd working_directory
git clone -b polkat_casa https://github.com/AKHughes1994/polkat.git
ln -s /idia/raw/point/to/your/file.ms .
```

The above example makes a symbolic link following the directory structure of ILIFU; if you are running locally instead, you would `mv` your file into `working_directory/`. Furthermore, in this example, the `git` call will make a directory called `polkat/` inside `working_directory/`. You will need to move the contents of `polkat/' one-level up, e.g.,:

```
mv polkat/* .
```

You run commands in `working_dir/` for the remainder of the data processing. It is always a good idea to check your ms-file to know your calibrator/target field names, observing band, frequency resolution, etc. This can done through the following command:

```
singularity exec /point/to/container/polkat-[version].sif python3 tools/ms_info.py [ms-file.ms]
```

Once you know your observations, you may begin processing. 

---

##### INFO

The first step, 'INFO,' can be run with the following commands: 

```
python3 setups/INFO.py idia
./submit_info_job.sh
```

If you are running polkat locally, replace `idia` with `node` and make sure the variable `NODE_CONTAINER_PATH` in `config.py` includes wherever you have the `polkat-[version].sif`

This step goes through your ms-file and extracts the info for the targets/calibrators, storing them in `project_info.json`. Furthermore, one change that has been made is INFO now splits out your desired fields and averages your ms-file down to (by default) 1024 frequency channels (oxkat did this at the start of the next step, 1GC). Moreover, the current version flaggs any short scans (< 10 seconds) if the dump time of your observations is 2 seconds as there is a known META data bug that will create short scans and mislabel the pointing direction (for most use cases, this removal will be desired unless you have target scans that are actually < 10 seconds). 

There are some `config.py` variables that will be commonly customised to fit your specific observations:

* `POLANG_NAME = 'J1331+3030'` — This is the name of the polarization angle calibrator as seen in the ms-file. The default is J1331+3030 (3C286). If left blank, the routine assumes no polarization angle calibration. 
* `POLANG_DIR  = '13:31:08.2881,+30.30.32.959'` — Coordinates of polarization angle calibrator. The default is J1331+3030 (3C286) coordindates
* `PRE_FIELDS = ''` — This is where you specify the fields of interest. For multi-target ms-files, it will run polkat on each target. If you are only interested in one target, you'll want to specify the names of the target, the target phase calibrator, primary, and (if applicable) the polarization angle calibrator

The two polarisation angle calibrators provided by SARAO/MeerKAT are J1331+3030 (3C286) and J0521+1638 (3C138), the default assumes 3C286 but in the config file. The correct parameters for 3C138 are also included. There are ways to use other calibrators if you are being created, but for if you do use a non-standard calibrator you likely an expert user and, as a result, are on your own ;). 

At the end of INFO, you should have your working ms-file that will (most likely) follow the naming covention `[ms-file]_1024ch.ms`. You are now ready to proceed.

---

##### 1GC

The second step, '1GC', performs reference calibration (i.e., it uses the calibrator fields to calibrate your target) using [casa](https://casa.nrao.edu/). Like 'INFO', it is run with: 

```
python3 setups/1GC.py idia
./submit_1GC_job.sh
```

At the time of writing this (Feb 13, 2025), MeerKAT ms-files that were made using the SARAO archive (i.e., pressing this button ![image](https://github.com/user-attachments/assets/05d49a3a-b4cf-42af-9f6c-9db278b647fe)) will have mislabelled the there X/Y feeds. This mislabeling will result in WRONG polarisation properties if not corrected (see EVLA Memo 219). polkat corrects this mislabelling in its first two steps. CAUTION: Other (less common) ways to pull data from the archive may also correct for this, and you don't want to double correct; polkat assumes you use the button, but you should make sure by knowing who or what you got your data from! A good way to know you didn't correct this effect is if the polarization properties of the diagnostic images (more on that later) strongly disagree with the archival values see [here](https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol). For most use cases, polkat should handle all this for you! 

Here are the most import `config.py` variables that will be commonly customised to fit specific observations:

* `CAL_1GC_DIAGNOSTICS = True` — This flag (on by default) will image your calibrators in Stokes I, Q, U, and V. The primary (also the leakage calibrator, typically J1939) should be unpolarized (Q = U = V = 0) and the polarization angle calibrator should have properties consistent with the catalogue values. Please leave this on; the extra time is worth knowing whether the calibration worked adequately.
* `CAL_1GC_AGGRESSIVE_FLAGS = False` — This flag (off by default) will more aggressively flag the baselines. By default, some RFI is only flagged on shorter baselines as they are more significantly affected. Turning this on flags that RFI on all baselines is necessary for high-precision polarimetry. The default flagging will result in spurious polarisation at the ~0.3% level, whereas turning it on the spurious signal is often <0.1% (this comes at the cost of some signal-to-noise).
* `POLANG_MOD  = [1.0, 0.0, 0.5, 0.0]` — This list contains a quasi-arbitrary initialization model for the polarisation angle calibrator (default is 3C286); it is fed into the casa command `setjy`. The `config.py` file also contains the model tested for 3C138.

The last parameter, `POLANG_MOD', is why I recommend leaving diagnostic imaging turned on. It turns out that for linear feed instruments, the casa cross-hand phase solver does not need the correct model; all it requires is to get the 'quadrant' of the angle approximately correct. This model independency is discussed in detail in EVLA Memo 219, and (I think) it needs an initialization because multiple solutions exist for the (model-independent) cross-hand phase. As a result of this pseudo-model independence, the polarisation angle calibrator also acts as a polarisation check source! All my testing converges on the correct solution despite the input model being junk, allowing us to get a feel for the systematic errors by comparing the measured values to the expected ones!

If 1GC is successful, please take a look at the visibility/gain solutions; as polarization can be finicky, you can now move on to 2GC (self-calibration and target imaging). 

---

##### 2GC

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

