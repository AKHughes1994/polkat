
If you make use of this software, please cite:

```bibtex
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

---


### Latest modification (April, 2026)

#### Update: UHF Band and Manual Cross-Hand Solver

Full UHF band support is now available. The cross-hand phase (XF) calibration has been extended with a **manual solver** that handles cases where CASA’s XF solver introduces large discontinuities, discussed more under INFO.

This routine has been tested on standard polarization calibrators **3C286** and **3C138** for **MeerKAT L-, S-, and UHF-bands**, and now recovers stable and correct polarization angles across the full frequency range.  


---

### What is this?

This repository contains a modified version of the MeerKAT semi-automated data processing routine [oxkat](https://github.com/IanHeywood/oxkat), enhanced to support full polarization calibration and Stokes I, Q, U, V imaging. It is assumed that you are already familiar with the oxkat workflow and file system, and that you configure data processing options by editing `oxkat/config.py`. This guide walks you through a standard use case, highlighting changes to `config.py` and introducing new options.

This branch, `polkat_QC_selfcal`, is designed to closely mimic the main `oxkat` branch, with the primary change being the upgrade from Cubical to Quartical for self-calibration. Unlike the main `polkat` branch (which uses CASA for self-calibration and does not split out the target field), this version is a near one-to-one adaptation, but with additional features and the ability to handle full polarization observations. An added benefit is that this branch includes the capability to perform the 3GC peeling step.

**IMPORTANT: For time series or multi-epoch campaigns**, it is strongly recommended that you download a local branch and use it consistently throughout your project. This repository is actively maintained, and improvements to flux density calibration or other processing steps may introduce global offsets or systematic differences between epochs processed with different versions. To ensure consistency across your dataset, fix your version at the start of your campaign.

---

#### Before We Start

This routine is designed primarily for use on the ILIFU clusters operated by The Inter-university Institute for Data Intensive Astronomy (IDIA), but you can run it locally if you have the required software. The necessary software is bundled into containers using [apptainer](https://apptainer.org/) (formerly singularity). Containers are available in `/software/containers` on ILIFU and follow the naming convention `polkat-[version].sif`. If you do not have ILIFU access but want to use polkat, you can get container access via pulling from dockerhub:

```
# The main container 
singularity pull polkat-0.2.4.sif docker://hughesakh/polkat:0.2.4

# ALBUS container for ionospheric corrections (no longer needed as the main has SPINIFEX)
singularity pull polkat-albus.sif docker://hughesakh/polkat_albus:latest
```

---

#### Standard Workflow

Assume a Linux-based OS (e.g., Ubuntu).

1. **Initialize a working directory and prepare your data:**

   ```bash
   mkdir working_directory
   cd working_directory
   git clone -b polkat_QC_selfcal https://github.com/AKHughes1994/polkat.git
   ln -s /idia/raw/point/to/your/file.ms .
   ```

   If running locally, use `mv` to move your file into `working_directory/`. The `git` command creates a `polkat/` directory; move its contents up one level:

   ```bash
   mv polkat/* .
   ```

   All commands should be run from `working_dir/` for the remainder of the process.

2. **Check your ms-file** to know your calibrator/target field names, observing band, frequency resolution, etc. For example:

   ```bash
   singularity exec /point/to/container/polkat-[version].sif python3 tools/ms_info.py [ms-file.ms]
   ```

---

## INFO

The first step, `INFO`, is run with:

```bash
python3 setups/INFO.py idia
./submit_info_job.sh
```

If running locally, replace `idia` with `node` and ensure `NODE_CONTAINER_PATH` in `config.py` points to your `polkat-[version].sif`.

This step extracts info for targets/calibrators from your ms-file, storing it in `project_info.json`. INFO also splits out your desired fields and averages your ms-file down to (by default) 1024 frequency channels (oxkat did this at the start of 1GC). 


**Commonly customized `config.py` variables:**

```python
POLANG_NAME = 'J1331+3030'  # Name of the polarization angle calibrator as it appears in the ms-file
                             # (default: J1331+3030/3C286). Leave blank if you do not wish to 
                             # perform polarization angle calibration.
POLANG_DIR = '13:31:08.2881,+30.30.32.959'  # Coordinates of the polarization angle calibrator 
                                             # (default: 3C286)
PRE_FIELDS = ''  # Specify fields of interest. For multi-target ms-files, polkat will process 
                 # each target. To process only a specific target, provide the names of the 
                 # target, phase calibrator, primary calibrator, and (if applicable) polarization 
                 # angle calibrator. If you use PRE_FIELDS, be sure to include the polarization 
                 # calibrator in the string.

# Cross-hand phase (XF) calibration parameters (must be set per calibrator; defaults for 3C286 and 3C138 provided)
XF_TARGET_POLANG = 30.0  # Expected intrinsic (RM-corrected) linear polarization angle [deg]
XF_TARGET_RM = 0.0       # Initial guess for the intrinsic rotation measure [rad/m^2]
XF_MODE = 'auto'         # Options: 'auto', 'casa', or 'manual'
                         # 'auto' (RECOMMENDED): uses CASA solver by default, switches to manual 
                         # if large phase discontinuities are detected (almost certainly due to 
                         # low cross-hand flux in some channels)
XF_AUTO_ANG_JUMP = 60.0  # Threshold [deg]; if adjacent CASA XF solutions differ by more than this,
                         # a manual re-solve is triggered
```

**Cross-hand phase (XF) calibration parameters** have been introduced to [`oxkat/config.py`](oxkat/config.py) and must be set **per calibrator** (defaults provided for 3C286 and 3C138).

The manual solver is implemented in [`tools/manual_XF_solver.py`](tools/manual_XF_solver.py).  
For most use cases, the default settings will suffice for UHF-through-S band --- especially if you're not deeply familiar with cross-hand calibration.
If you develop or know of a more robust cross-hand phase solver, please get in touch!

- ### Note: No Polarization Angle Calibrator?

If no polang calibrator is available, make sure `POLANG_NAME = ''`. The primary (e.g., J1939) must be unpolarized, so you can still solve for leakage terms.

SARAO/MeerKAT provides two polarization angle calibrators: J1331+3030 (3C286, default) and J0521+1638 (3C138). Using a non-standard calibrator strategies is possible, even without a standard calibrator but will be added later and is intended for expert users.

**Metadata bug workaround:**

A block of code has been added to `oxkat/PRE_casa_average_to_1k_add_wtspec.py` to address a known metadata issue with 2-second integration observations:

```python
# Remove short scans that arise from metadata error from 2s integration observations
bad_scans = []
good_scans = []

tb.open(master_ms)
scans = np.unique(tb.getcol('SCAN_NUMBER'))
for scan in scans:
    subtab = tb.query(query='SCAN_NUMBER=='+str(scan))  # scan info
    scan_times = np.unique(subtab.getcol('TIME'))  # scan integration times
    scan_dt = scan_times[-1] - scan_times[0]  # total scan length (s)
    integration = scan_times[1] - scan_times[0]  # integration length (s)
    if scan_dt < 10.0 and integration < 2.5:
        bad_scans.append(str(scan))
    else:
        good_scans.append(str(scan))
tb.close()

if myscans != '':
    myscans = myscans.split(',')
    myscans = ','.join([scan for scan in myscans if scan not in bad_scans])
else:
    myscans = ','.join(good_scans)
```

This code flags any short scans (less than 10 seconds) when the dump time is 2 seconds, due to a known metadata bug that can mislabel the pointing direction. For most users, this removal is desirable unless you have genuine target scans shorter than 10 seconds. If this behavior is not desired, you can manually comment out this code block—it will not affect the rest of the workflow.

At the end of the INFO step, you should have a working ms-file, typically named `[ms-file]_1024ch.ms`. You are now ready to proceed to the next stage.

---

## 1GC

The second step, `1GC`, performs reference calibration (using calibrator fields to calibrate your target) with [casa](https://casa.nrao.edu/):

```bash
python3 setups/1GC.py idia
./submit_1GC_job.sh
```

**Important:** As of LAST UPDATE, MeerKAT ms-files made using the SARAO archive (i.e., using the download button) have mislabelled X/Y feeds. This results in incorrect polarization properties if not corrected (see EVLA Memo 219). polkat corrects this mislabelling in its first two steps. **CAUTION:** Other (less common) archive download methods may already correct this; do not double-correct. polkat assumes you used the button. If polarization properties of diagnostic images strongly disagree with archival values (see [here](https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol)), you may have missed or double-applied this correction.

**Commonly customized `config.py` variables:**

```python
CAL_1GC_DIAGNOSTICS = True  # (default: on) Images calibrators in Stokes I, Q, U, V. The primary 
                             # (typically J1939) should be unpolarized; the polarization angle 
                             # calibrator should match catalogue values. Leave this on for 
                             # calibration checks.
CAL_1GC_AGGRESSIVE_FLAGS = False  # (default: off) More aggressive baseline flagging. Default flags 
                                   # RFI on short baselines; turning this on flags RFI on all 
                                   # baselines, can remove RFI inducing low-fractional polarisation 
                                   # structure.
POLANG_MOD = [1.0, 0.0, 0.5, 0.0]  # Initialization model for the polarization angle calibrator 
                                    # (default: 3C286); used by setjy in CASA. The config file also 
                                    # includes a model for 3C138. This model is used by the default 
                                    # CASA cross-hand phase (XF) solver.
```

**Note:** For linear feed instruments, the CASA cross-hand phase solver only needs the angle quadrant approximately correct (see EVLA Memo 219). The polarization angle calibrator also acts as a polarization check source. All testing converges on the correct solution despite the input model, allowing you to estimate systematic errors by comparing measured and expected values.

If 1GC is successful, check the visibility/gain solutions. Polarization can be finicky. You can now move on to 2GC (self-calibration and target imaging), and the pipeline will split out individual MS files for your target field(s). 

---

## 2GC

After completing 1GC, move on to the 2GC step. This stage performs final flagging, imaging (using [wsclean](https://wsclean.readthedocs.io/en/latest/)), and direction-independent phase self-calibration. The 2GC process produces both channelized images and a single MFS (multi-frequency synthesis) image. The MFS image maximizes sensitivity, but may be affected by bandwidth depolarization.

**To run 2GC:**

```bash
python3 setups/2GC.py idia
./submit_2GC_job.sh
```

For fields with clear artefacts around bright sources (possibly due to calibration errors or residual RFI), an alternative two-stage self-calibration approach is available via `python3 waterhole/setup_2GC_twostage.py idia`. This performs iterative self-calibration first on pixels meeting the `WSC_SHALLOWMASK` threshold before proceeding to the full imaging workflow.

Most imaging parameters for 2GC are found in `config.py` after the comment `# wsclean and 2GC defaults`, and correspond directly to [wsclean](https://wsclean.readthedocs.io/en/latest/) options. **Please consult the wsclean documentation for detailed explanations of all available parameters.** The variables listed below are the ones most commonly adjusted for typical use cases (but note that the defaults are generally suitable for standard MeerKAT observations).:

**Commonly customized `config.py` variables:**

```python
# Data selection
WSC_MINUVL = ''  # Minimum uv-distance in wavelengths (leave empty for no restriction)
WSC_MAXUVL = ''  # Maximum uv-distance in wavelengths (leave empty for no restriction)
                 # Useful for matching angular scales across frequency
WSC_TUKEYTAPER = False  # Apply Tukey taper to reduce edge effects in uv-coverage
WSC_TAPERMASK = False   # Apply taper/uvcuts to blind and mask images (not just final images)

# Weighting
WSC_WEIGHT = 'briggs 0.0'  # Imaging weighting scheme (default: Briggs robustness 0.0)
                            # Maximizes sensitivity before MeerKAT beam becomes non-Gaussian
WSC_UNIFORM_IMAGE = True   # Also generate high-resolution image using WSC_WEIGHT_HIGHRES
                           # Useful for tracking proper motion

# Polarization and spectral settings
WSC_POL = 'IQUV'  # Stokes parameters to image. Pipeline designed for either Stokes I or 
                  # full IQUV. Atypical subsets (e.g., QUV, UV) may work but other pipeline 
                  # steps (e.g., RMSYNTH) may fail.
WSC_SPLITPOL = False  # Image V/I and Q/U separately (necessary for high RM and MFS fitting)
WSC_JOINPOLARIZATIONS = True  # Join polarizations during deconvolution
WSC_SQUAREPOLARIZATIONS = False  # Use squared polarization during deconvolution

# Frequency channels
WSC_BLIND_CHANNELSOUT = 8  # Number of channels for initial blind imaging
WSC_PCAL_CHANNELSOUT = 8   # Number of channels for final post-calibration imaging
                           # Increase for sources with high rotation measure
WSC_DMASK_CHANNELSOUT = WSC_PCAL_CHANNELSOUT  # Channels for datamask image (in case you want 
                                               # different channelization than pcalmask)
WSC_CAL_CHANNELSOUT = 16   # Number of channels for calibrator imaging
WSC_MAX_CHANNELS = 16      # Maximum channels imaged at once (memory limit). If 
                           # WSC_PCAL_CHANNELSOUT > WSC_MAX_CHANNELS, imaging is performed 
                           # in steps and MFS is created by stacking in image plane 
                           # (qualitative only; do not report fluxes from this MFS image)
WSC_FITSPECTRALPOL = 4     # Polynomial order for spectral polarization fitting
WSC_JOINCHANNELS = True    # Join channels during deconvolution

# Deconvolution constraints
WSC_NONEGATIVE = False     # Enforce non-negative constraint during deconvolution
WSC_STOPNEGATIVE = False   # Stop deconvolution when negative components appear
WSC_CIRCULARBEAM = False   # Force circular restoring beam

# Masking for blind imaging (initial image to create mask)
WSC_AUTOMASK_BLIND = 10.0        # Auto-masking sigma threshold for blind image
WSC_AUTOTHRESHOLD_BLIND = 2.0    # Auto-threshold sigma for blind cleaning
WSC_THRESHOLD_BLIND = False      # Manual threshold for blind image (False = auto)
WSC_LOCALRMS_BLIND = False       # Use local RMS for blind image masking (False = global RMS)

# Masking for science imaging (final images)
WSC_MASK = False            # Use existing mask file (False = auto-generate from blind image)
WSC_THRESHOLD = False       # Manual cleaning threshold (False = auto)
WSC_SHALLOWMASK = 50.0      # Shallow mask sigma threshold for initial mask creation
WSC_AUTOMASK = 4.0          # Auto-masking sigma threshold for science imaging
WSC_AUTOTHRESHOLD = 1.0     # Auto-threshold sigma for final cleaning
WSC_LOCALRMS = False        # Use local RMS for masking (False = global RMS)
```


**Multi-frequency weighting (IMPORTANT):**

This pipeline supports spectrally resolved functionality (e.g., Rotation Measure synthesis). The `WSC_MFWEIGHT` parameter controls whether multi-frequency weighting is applied during imaging. It is set to `True` by default, optimized for a single combined MFS image with the lowest RMS noise. **If you plan to use spectral functionality, set this to `False` to avoid introducing artificial spectral structure** (see [wsclean MF weighting documentation](https://wsclean.readthedocs.io/en/latest/mf_weighting.html)).

**Self-calibration control:**

Self-calibration solutions are controlled via Quartical YAML files in `data/quartical`. The default configuration (`2GC_phase.yaml`) performs phase-only self-calibration with a frequency-dependent slope. See the [Quartical documentation](https://quartical.readthedocs.io/en/latest/) for details. The key parameter is `time_interval='1'`, which sets the solution interval to 1 integration time by default.

**Note on parallactic angle:** The QC routine does not de-rotate or re-apply parallactic angle corrections during self-calibration, as the DATA column already has parallactic angle rotation applied after 1GC. While there was initial concern this could reduce calibration fidelity (due to feed-frame vs. sky-frame issues), in practice this effect appears minor—much smaller than measurement errors or other systematics. If you encounter a case where this makes a significant difference, please report it.

**Output:**  
Final images will have suffixes such as `pcalmask-MFS-I-image.fits` or `pcalmask-[CHAN_NUMBER]-I-image.fits` (`I` for Stokes I). In addition to Stokes images, 2GC also produces linear polarization intensity (`Plin`, $\sqrt{Q^2+U^2}$) and total polarization intensity (`Ptot`, $\sqrt{Q^2+U^2+V^2}$) images.

---

## 3GC (Peeling)

The 3GC (Third Generation Calibration) "peeling" step is an advanced calibration and imaging process designed to further improve image fidelity by correcting for direction-dependent effects (DDEs) around bright sources in your field. Peeling is especially useful for fields with strong off-axis sources or complex extended emission, where standard direction-independent calibration is insufficient.

**How does 3GC_peel work?**

- Peeling directions are defined using DS9 region files (one region per source/direction). You can specify these manually by setting the `CAL_3GC_PEEL_REGION` variable in `config.py` to a comma-separated list of region file paths. Alternatively, you can place region files in your working directory using the naming convention `[SOURCE_NAME]_peel[PEEL_NUMBER].reg` (e.g., `CygX1_peel1.reg`, `CygX1_peel2.reg`, etc.).
- The pipeline supports an arbitrary number of peel directions, but note that each direction adds a new data column to your MS file, which can rapidly increase its size.
- Peeling involves amplitude and phase self-calibration in each direction and therefore may have weird interactions with instrumental, off-axis leakage (which MeerKAT has a lot of). In my tests, peeling generally preserves source fidelity (including polarization fluxes), but if you notice significant changes in target flux densities, especially for polarization, please report these cases.

**What does 3GC_peel do?**

- Identifies and calibrates bright sources that may cause significant DDEs.
- Iteratively subtracts these sources from the data (the "peeling" process).
- Applies direction-dependent calibration solutions to improve the dynamic range and fidelity of the final images.
- Produces both residual and "peeled" images for scientific analysis.

**To run 3GC_peel:**

```bash
python3 setups/3GC_peel.py idia
./submit_3GC_peel_job.sh
```

**Output:**

- Peeled measurement sets and images for each direction/source.
- Residual images with the peeled sources removed.
- Diagnostic plots and logs for each peeling iteration.

**Notes and Best Practices:**

- Peeling is computationally intensive. Ensure you have sufficient resources, especially for large fields or many peel directions.
- For most science cases, peeling is only necessary if you see strong artefacts around bright sources after 2GC.
- All 3GC_peel parameters can be found and adjusted in `config.py` under the section labeled `# 3GC peeling defaults`.

For more details on the theory and practice of peeling and direction-dependent calibration, see [arxiv:1101.1765](https://arxiv.org/abs/1101.1765) and the [oxkat documentation](https://github.com/IanHeywood/oxkat).

---

## SNAP

This step performs snapshot imaging following the Heywood-ian approach, efficiently producing short-timescale images by subtracting a time-averaged model and searching for image-plane variability. This is particularly useful for detecting transient or variable sources.

**To run SNAP:**

```bash
python3 setups/SNAP.py idia
./submit_snap_job.sh
```

**Commonly customized `config.py` variables** (see `# Snapshot imaging defaults` section):

```python
SNAP_FIELDS = ''  # Fields to perform snapshot imaging on (comma-separated). 
                  # Leave empty to process all target fields; specify field names 
                  # to process only specific targets.
SNAP_CHANNELSOUT = WSC_PCAL_CHANNELSOUT  # Number of frequency channels for snapshot images 
                                          # (default: same as WSC_PCAL_CHANNELSOUT)
SNAP_INTBIN = 1   # Number of integration times to bin together per snapshot image 
                  # (default: 1 = per-integration imaging)
SNAP_INTEND = True  # If True, discard incomplete bins at the end of the observation 
                    # to ensure all snapshot images have equal time spacing
                    # Example: 17 integrations with INTBIN=4 → images 1-16, discard 17
SNAP_POL = True   # If True, produce full IQUV snapshot images; if False, Stokes I only
SNAP_IMSIZE = 2560  # Image size in pixels for snapshot images (not the model image)
SNAP_MODELIDENTIFIER = 'pcalmask'  # Identifier for the model image used for subtraction 
                                    # (follows oxkat/polkat naming conventions)
SNAP_MODELMASK = ''  # Path to mask for initial model creation (or leave empty to use 
                     # existing pcalmask model from 2GC)
SNAP_DECONV = False  # If True, perform deconvolution during snapshot imaging 
                     # (rarely necessary for typical variability studies)
SNAP_DECONVMASK = ''  # Mask to use if SNAP_DECONV=True
```

**Optional post-processing:**

Once snapshot imaging completes, you can generate movies to visually inspect for variables:

```bash
python3 waterhole/setup_movie.py idia
```

This creates MP4 files in your `INTERVALS/` directories, providing a convenient way to identify variable sources in your snapshot images.

---

## RMSYNTH

The RMSYNTH step automates the process of extracting fluxes and polarization properties from images for arbitrary point sources. **Please verify that the resulting fluxes and measurements are physically reasonable for your science case.**

**What this step does:**

1. Fits target MFS, channelized, and/or snapshot images with user-specified Gaussians using the CASA task `imfit`.
2. If `CAL_1GC_DIAGNOSTICS = True`, quantifies systematic calibration effects via image-plane analysis of calibrators.
3. Runs RM Synthesis on every Gaussian component from every target, extracting polarization angles and rotation measures.
4. Runs ALBUS to estimate the ionospheric RM for post-processing corrections. (Note: ALBUS may fail randomly; if so, rerun or change `RED_TYPE = RI_G03` to `RED_TYPE = RI_G01` in `config.py`.)

You may need to modify `data/rmsynth/rmsynth_info.json` to match your dataset. Example configuration:

```json
{
    "image_directory": ["IMAGES", "INTERVALS"],
    "image_identifier": ["pcalmask", "restored"],
    "image_suffix": ["image.fits", "image.fits"],
    "image_timing": [false, true],
    "source_name": ["SwiftJ1727", "SwiftJ1727"],
    "source_ulim": [[false], [false]],
    "rms_region": [false, false],
    "source_pos": [["17:27:43.307,-16.12.17.619"], ["17:27:43.307,-16.12.17.619"]],
    "pos_coeff": [[[]], [[]]],
    "time0": [[], []]
}
```

**Parameter explanations:**

- `"image_directory"`: Directory containing images to fit (`IMAGES` for 2GC, `INTERVALS` for snapshots).
- `"image_identifier"`: Identifier for images (`pcalmask` for 2GC, `restored` for snapshots).
- `"image_suffix"`: Image file type (`image.fits` for standard, `image.homogenized.fits` for homogenized images).
- `"image_timing"`: Set to `true` if images are time-split (e.g., snapshots).
- `"source_name"`: Names of sources to fit.
- `"source_ulim"`: Whether to fix the component position to `"source_pos"`.
- `"rms_region"`: Manual RMS region (CASA format); otherwise, an annulus is used.
- `"source_pos"`: List of component positions (CASA format). Add more entries for multiple components.

You can ignore `"pos_coeff"` and `"time0"` for most use cases. The default configuration should work for typical targeted point source observations—just update `"source_name"` and `"source_pos"` as needed. If you did not perform snapshot imaging, remove the second entry from each list.

**Warning:**  
While this step streamlines the process, it is important to understand the underlying tools and concepts—especially Faraday Rotation and RM Synthesis. For more background, see [Brentjens & de Bruyn 2005](https://arxiv.org/abs/astro-ph/0507349).

---

## THE COMPRISING PACKAGES

polkat is fundamentally a compilation of existing software; these packages are the lifeblood and deserve all the credit!

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
| [`RMTools`](https://github.com/CIRADA-Tools/RM-Tools?tab=readme-ov-file) | RMSynth | Rotation Measure Synthesis | [Van Eck et al.](https://ui.adsabs.harvard.edu/abs/2020ascl.soft05003P/abstract)|
| [`ALBUS`](https://github.com/twillis449/ALBUS_ionosphere) | RMSynth | Ionospheric RM Estimation | Willis et al.|


---

## TO DO LIST

1. Further improve documentation and include a PDF with some examples.
2. Make a Quartical only branch


