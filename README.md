If you make use of this software please cite: 

```
@software{2025ascl.soft02026H,
       author = {{Hughes}, Andrew K. and {Cowie}, Fraser J. and {Heywood}, Ian and {Hugo}, Ben},
        title = "{polkat: Semi-automate full polarization of MeerKAT observations}",
 howpublished = {Astrophysics Source Code Library, record ascl:2502.026},
         year = 2025,
        month = feb,
          eid = {ascl:2502.026},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2025ascl.soft02026H},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```

### Caution

Currently, this routine has been extensively tested for MeerKATs L- and S-bands. We are looking into UHF band, and while the current version seems to get the correct polarisation degrees, the polarization angle < 800 MHz is providing problems. As a result, you should be careful (and check your data) for UHF band. We will update this documentation as we figure out what is going on.

--- 

### What is this?

This modified version of the MeerKAT semi-automated data processing routine [oxkat](https://github.com/IanHeywood/oxkat) that was tweaked to include full polarisation calibration and for Stokes I, Q, U, V imaging. We assume you are familiar with the oxkat workflow and its file system and that data processing options are contained in and changed by editing `oxkat/config.py`. This guide will walk you through a standard use case and highlight some changes to the `config.py` file and new options. 


---
##### Before We Start 

This routine has been designed primarily for use on the ILIFU clusters operated by The Inter-university Institute for Data Intensive Astronomy (IDIA); however, you can run it locally if you have the comprising software. The necessary software has been combined into containers using [apptainer](https://apptainer.org/) (previously known as singularity). The containers can be found in the `/software/containers` directory on ILIFU and follow the naming convention of `polkat-[version].sif`. If you do not have access to ILIFU but still want to use polkat and thus require the container, you can pull them from dockerhub:

```
# The main container -- CURRENTLY DOESN'T WORK
singularity pull polkat-0.1.2.sif docker://hughesakh/polkat:0.1.2

# ALBUS container for ionospheric corrections
singularity pull polkat-albus.sif docker://hughesakh/polkat_albus:latest
```

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

After 1GC is complete, final flagging, imaging (with WSCLEAN), and direction-independent phase self-calibration are performed. The image products will include channelised images, as well as a single MFS image that combines the data from the different frequency channels into a single Multi-frequency synthesis image (maximizing sensitivity, but may cause bandwidth depolarisation)

   ```
   $ python setups/2GC.py idia
   $ ./submit_2GC_job.sh
   ```

There are too many bells and whistles to list all of them (feel free to discuss them if there are any questions). Some key ones are:

* `WSC_WEIGHT  = briggs 0.0` — The weighting for the imaging; by default it uses a [Briggs](https://casa.nrao.edu/Documents/Briggs-PhD.pdf) robustness of zero. This is the most natural (i.e., sensitivity maximizing) weighting before the MeerKAT synthesized beam becomes extremely non-Gaussian (and thus hard to deconvolve via CLEAN-based imagers)
* `WSC_UNIFORM_IMAGE = True` — This flag will make a high-angular resolution image in addition to the standard images using the weighting set by the variable `WSC_WEIGHT_HIGHRES`. It is helpful if you have a weakly polarised point source where you want natural weighting to maximize the sensitivity for pol. Detections, by high angular resolution for Stokes I astrometry.
* `WSC_POL = 'IQUV'` — Specify which Stokes parameters you want images; just leave this alone. 
* `WSC_IMAGE_CHANNELSOUT = 8` — This sets the frequency channels to be images, which may need to be increased if you are doing polarisation imaging of a source with a considerable rotation measure.
* `WSC_MAX_CHANNELS = 16` — This is new to polkat; it says a maximum amount of frequency channels to be imaged at a given time as you tend to run out of memory if you go above 16 channels with IQUV imaging. If `WSC_IMAGE_CHANNELSOUT` > `WSC_MAX_CHANNELS`, the channelized imaging is broken into multiple steps; the downside is that you no longer get an MFS image stacked in visibility space. CAUTION: polkat will make an MFS image by homogenizing the beam in frequency and stacking in the image plane, but this should be taken as qualitative; DO NOT REPORT THE FLUXES FROM THE MFS IMAGE IN THIS CASE.
* `WSC_HOMOGENIZEBEAM = False`  — This is new to polkat; it produces a new set of images with the suffix `.homogenized.fits`. It uses the Welzl algorithm to solve for the [minimum enclosing eclipse](https://github.com/dorshaviv/lowner-john-ellipse) and then enforces that as the synthesized beam shape, thereby homogenizing angular resolution in frequency (necessary to make, for example, spectral index maps from the channelized images). It also homogenizes the residual images using [pypher](https://pypher.readthedocs.io/en/latest/). CAUTION: This should be seen as an experimental as the residual homogenization is a pseudo-averaging that drastically reduces the RMS in the upper-frequency band; more thought needs to go into this to understand the correct way to interpret the signal-to-noise of detections post-homogenization. 

Also, you can control the solution interval for self-calibration using `CAL_2GC_PSOLINT=32s`; by default, this solves in 32s intervals (i.e. four 8s dump times following the oxkat default). Decrease this number if you have a good model of a bright target field.

Please don't hesitate to reach out if you'd like to discuss any of the other parameters. Most `config.py` imaging parameters are found after the comment `# wsclean and 2GC defaults` and are just links to the various [wsclean](https://wsclean.readthedocs.io/en/latest/) options.

For most people, the final image(s) you will be working with will have suffixes like `pcalmask-MFS-I-image.fits` or `pcalmask-[CHAN_NUMBER]-I-image.fits`; where in this case 'I' identifies the image as a Stokes I (or total intensity image). One step of 2GC will also produce linear polarisation intensity images (`Plin`, $\sqrt{Q^2+U^2}$) and total polarisation intensity images (`Ptot`; $\sqrt{Q^2+U^2+V^2}$) images.

---

##### SNAP

This new step will perform snapshot imaging following the Heywood-ian approach. It will (efficiently) make short-timescale images by subtracting a time-averaged model and then looking for image plane variability in the post-subtracted images. It can be run with. 

   ```
   $ python setups/SNAP.py idia
   $ ./submit_snap_job.sh
   ```

The relevant parameters in `config.py` can be found after the comment `Snapshot imaging defaults' below are the key variables:
   
   * `SNAP_FIELDS = ''` — Names of the fields to perform snapshot imaging; the default will perform the snapshot imaging on every field (so you'll probably want to specify your target)
   * `SNAP_INTBIN = 1` — Integer number of dump times per snapshot image; default of 1 will make one snapshot image per dump time; 4 would make combine 4 dump times in each snapshot image
   * `SNAP_INTEND = True`—This will throw away the edge bins if they don't equal `SNAP_INTBIN`; e.g., for observation with 11 total dump times and 'SNAP_INTBIN = 2', `SNAP_INTEND = True` will result in 5 snapshot images being made using dump times 1 to 10, but it will throw away the last dump time (as it can't be grouped into 2) to ensure equispaced (in time) snapshot images. This works across scan boundaries.
   * `SNAP_Pol = True` — Default makes IQUV snapshots; turn it off to make only Stokes I snapshots.

The rest of the options should be reasonably self-explanatory! Reach out if anything is confusing.

---

##### RMSYNTH

This new step seeks to automate going from images to fluxes for some arbitrary number of point sources; please verify that the fluxes make sense;

What it does is:

   1. Fit the target MFS/channelized/snapshot images with some user-specified number of Gaussians using the casa task `imfit`
   2. If `CAL_1GC_DIAGNOSTICS = True` (which it should be), it will quantify the systematic calibration effects through image plane analysis of the calibrators
   3. Run RM Synthesis on every Gaussian component from every target, extracting polarisation angles/rotation measures
   4. Run ALBUS to get ionospheric RM for post-processing corrections. This randomly fails sometimes. You can (i) rerun it; (ii) Change `RED_TYPE = RI_G03` to `RED_TYPE = RI_G01` in `config.py`; the former will fail if only one GPS station is operational

The only thing that needs to be modified is the file `data/rmsynth/rmsynth_info.json`. Below is the default file contents that I left as an example:

```
{
    "image_directory": ["IMAGES", "INTERVALS"],
    "image_identifier": ["pcalmask", "restored"],
    "image_suffix": ["image.fits", "image.fits"],
    "image_timing": [false, true],
    "source_name": ["SwiftJ1727", "SwiftJ1727"],
    "source_ulim":[[false], [false]],
    "rms_region": [false, false],
    "source_pos": [["17:27:43.307,-16.12.17.619"], ["17:27:43.307,-16.12.17.619"]],
    "pos_coeff":[[[]], [[]]],
    "time0":[[], []]
}
```

It will extract the properties for N sets of images where N is set by the length of "image_directory" (and all the other parameters) 

Here is a list of what each entry does:

   * "image_directory": Specify the directory where the images you are going to fit are being held; IMAGES will contain 2GC images; INTERVALS will contain snapshot images
   * "image_identifier": Specify the identifier used to locate the images; pcalmask points to the 2GC images products and restored points to the snapshot images
   * "image_suffix": Specify the type of image, "image.fits" is the standard wsclean output; the other option would be "image.homogenized.fits" if you want the homogenized resolution images
   * "image_timing": Specify if the images have been broken up in time (e.g., snapshot images), this will change the wsclean output names to include a '-t[SNAPSHOT_NUMBER]' in the file name.
   * "source_name": Self-explanatory
   * "source_ulim": Determines whether of not to fix the position of the component to the value given by "source_pos"
   * "rms_region": Option to specify a manual RMS region (in standard CASA region format), if not, it will be calculated within the software as an annulus
   * "source_pos": List of component positions in standard CASA format. In this example, we are only fitting one component per image set. If you want to fit multiple components, just add another entry, e.g., "source_pos": [["17:27:43.307,-16.12.17.619", "17:27:44.307,-16.12.20.619"], ["17:27:43.307,-16.12.17.619", "17:27:44.307,-16.12.20.619"]],

You don't need to worry about the bottom two; I haven't implemented them yet. The default should be good for most targetted point source observations; all you need to do is change "source_name" and "source_pos". If you didn't do snapshot imaging, remove the second entry from each of the lists!

---

##### Note: No Polarization Angle Calibrator

polkat will work without a polarization angle calibrator. As long as the primary can work as an unpolarized leakage calibrator, it will be sufficient to measure the total polarization; you cannot distinguish between circular and linearly polarized emission. However, given that a lot of radio is synchrotron sources, the VAST majority will be circularly unpolarized, and thus, the total polarization is a good measure of linear polarization without a polarization angle calibrator. 

Leave `POLANG_NAME = ''` for these cases

---
## THE COMPRISING PACKAGES

polkat is fundamentally a compilation of actual software; these packages are the actual lifeblood and deserve all the credit!

| Package | Stage | Purpose | Reference |
| --- | --- | --- | --- | 
| [`astropy`](https://www.astropy.org/) | 1GC, 2GC, 3GC | Coordinates, time standards, FITS file manipulation | [Astropy Collaboration, 2013](https://ui.adsabs.harvard.edu/abs/2013A%26A...558A..33A/abstract), [Astropy Collaboration, 2018](https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A/abstract)|
| [`CASA`](https://casa.nrao.edu/) | 1GC, 2GC | Averaging, splitting, cross calibration, DI self-calibration, flagging | [The CASA Team, et al.](https://ui.adsabs.harvard.edu/abs/2022PASP..134k4501C/abstract)|
| [`Quartical`](https://github.com/ratt-ru/CubiCal) | 2GC, 3GC | DI / DD self-calibration | [Kenyon et al., 2018](https://ui.adsabs.harvard.edu/abs/2024arXiv241210072K/abstract)|
| [`DDFacet`](https://github.com/saopicc/DDFacet) | 3GC | Imaging with direction-dependent corrections | [Tasse et al., 2018](https://ui.adsabs.harvard.edu/abs/2018A%26A...611A..87T/abstract) | 
| [`killMS`](https://github.com/saopicc/killMS) | 3GC | DD self-calibration| [Tasse, 2014](https://ui.adsabs.harvard.edu/abs/2014A%26A...566A.127T/abstract); [Smirnov & Tasse, 2014](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.2668S/abstract) |
| [`owlcat`](https://github.com/ska-sa/owlcat/) | 2GC, 3GC | FITS file manipulation | - |
| [`ragavi`](https://github.com/ratt-ru/ragavi/) | 1GC, 2GC | Plotting gain solutions | - |
| [`shadeMS`](https://github.com/ratt-ru/shadeMS/) | 1GC | Plotting visibilities | [Smirnov et al., 2022](https://ui.adsabs.harvard.edu/abs/2022ASPC..532..385S/abstract) |
| [`Singularity and/or Apptainer`](https://apptainer.org/) | 1GC, 2GC, 3GC | Containerisation | [Kurtzer, Sochat & Bauer, 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459) |
| [`tricolour`](https://github.com/ska-sa/tricolour) | 2GC | Flagging | [Hugo et al., 2022](https://ui.adsabs.harvard.edu/abs/2022arXiv220609179H/abstract) |
| [`wsclean`](https://gitlab.com/aroffringa/wsclean) | 2GC, 3GC | Imaging, model prediction | [Offringa et al., 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.444..606O/abstract)|
| [`pypher`](https://pypher.readthedocs.io/en/latest/) | 2GC, 3GC | Residual Homogenization | [Boucaud et al., 2016](https://ui.adsabs.harvard.edu/abs/2016A%26A...596A..63B/abstract)|
| [`RMTools`](https://github.com/CIRADA-Tools/RM-Tools?tab=readme-ov-file) | RMSynth | Rotation Measure Synthesis | Van Eck et al.|
| [`ALBUS`](https://github.com/twillis449/ALBUS_ionosphere) | RMSynth | Ionospheric RM Estimation | Willis et al.|

---

## TO DO LIST

1. Further improve documentation and include a PDF with some examples
2. Make a new version called polkat_QC that uses Quartical for self-calibration following a similar routine as oxkat
3. Add 3GC, specifically after 2 as it will be easier to adapt



