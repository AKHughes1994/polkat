### Caution

This routine has been extensively tested for MeerKAT L- and S-bands. UHF band support is experimental: while the current version recovers correct polarization degrees, the polarization angle below 800 MHz is problematic. Please check your data carefully for UHF band reductions. We will update this documentation as we resolve these issues.

---

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

### What is this?

This repository contains a modified version of the MeerKAT semi-automated data processing routine [oxkat](https://github.com/IanHeywood/oxkat), enhanced to support full polarization calibration and Stokes I, Q, U, V imaging. It is assumed that you are already familiar with the oxkat workflow and file system, and that you configure data processing options by editing `oxkat/config.py`. This guide walks you through a standard use case, highlighting changes to `config.py` and introducing new options.

This branch, `polkat_QC_selfcal`, is designed to closely mimic the main `oxkat` branch, with the primary change being the upgrade from Cubical to Quartical for self-calibration. Unlike the main `polkat` branch (which uses CASA for self-calibration and does not split out the target field), this version is a near one-to-one adaptation, but with additional features and the ability to handle full polarization observations. An added benefit is that this branch includes the capability to perform the 3GC peeling step.

---

#### Before We Start

This routine is designed primarily for use on the ILIFU clusters operated by The Inter-university Institute for Data Intensive Astronomy (IDIA), but you can run it locally if you have the required software. The necessary software is bundled into containers using [apptainer](https://apptainer.org/) (formerly singularity). Containers are available in `/software/containers` on ILIFU and follow the naming convention `polkat-[version].sif`. If you do not have ILIFU access but want to use polkat, email 'hughesakh [at] gmail [dot] com' for a download link.

---

#### Standard Workflow

Assume a Linux-based OS (e.g., Ubuntu).

1. **Initialize a working directory and prepare your data:**

   ```bash
   mkdir working_directory
   cd working_directory
   git clone -b polkat_casa https://github.com/AKHughes1994/polkat.git
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

- `POLANG_NAME = 'J1331+3030'` — Name of the polarization angle calibrator as it appears in the ms-file (default: J1331+3030/3C286). Leave blank if you do not wish to perform polarization angle calibration.
- `POLANG_DIR  = '13:31:08.2881,+30.30.32.959'` — Coordinates of the polarization angle calibrator (default: 3C286).
- `PRE_FIELDS = ''` — Specify fields of interest. For multi-target ms-files, polkat will process each target. To process only a specific target, provide the names of the target, phase calibrator, primary calibrator, and (if applicable) polarization angle calibrator. **If you use `PRE_FIELDS`, be sure to include the polarization calibrator in the string.**

SARAO/MeerKAT provides two polarization angle calibrators: J1331+3030 (3C286, default) and J0521+1638 (3C138). Parameters for 3C138 are also included in the config file. Using a non-standard calibrator is possible, but intended for expert users.

A block of code has been added to `oxkat/PRE_casa_average_to_1k_add_wtspec.py` to address a known metadata issue:

```
# Remove short scans that arise from metadata error from 2s integration observations
bad_scans = []
good_scans = []

tb.open(master_ms)
scans = np.unique(tb.getcol('SCAN_NUMBER'))
for scan in scans:
    subtab = tb.query(query='SCAN_NUMBER=='+str(scan)) # scan info
    scan_times = np.unique(subtab.getcol('TIME')) # scan integration times
    scan_dt = scan_times[-1] - scan_times[0] # total scan length (s)
    integration = scan_times[1] - scan_times[0] # integration length (s)
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

**Important:** As of May 22, 2025, MeerKAT ms-files made using the SARAO archive (i.e., using the download button) have mislabelled X/Y feeds. This results in incorrect polarization properties if not corrected (see EVLA Memo 219). polkat corrects this mislabelling in its first two steps. **CAUTION:** Other (less common) archive download methods may already correct this; do not double-correct. polkat assumes you used the button. If polarization properties of diagnostic images strongly disagree with archival values (see [here](https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol)), you may have missed or double-applied this correction.

**Commonly customized `config.py` variables:**

- `CAL_1GC_DIAGNOSTICS = True` — (default: on) Images calibrators in Stokes I, Q, U, V. The primary (typically J1939) should be unpolarized; the polarization angle calibrator should match catalogue values. Leave this on for calibration checks.
- `CAL_1GC_AGGRESSIVE_FLAGS = False` — (default: off) More aggressive baseline flagging. Default flags RFI on short baselines; turning this on flags RFI on all baselines, improving polarimetry at the cost of some S/N.
- `POLANG_MOD  = [1.0, 0.0, 0.5, 0.0]` — Initialization model for the polarization angle calibrator (default: 3C286); used by `setjy` in casa. The config file also includes a model for 3C138.

**Note:** For linear feed instruments, the casa cross-hand phase solver only needs the angle quadrant approximately correct (see EVLA Memo 219). The polarization angle calibrator also acts as a polarization check source. All testing converges on the correct solution despite the input model, allowing you to estimate systematic errors by comparing measured and expected values.

If 1GC is successful, check the visibility/gain solutions. Polarization can be finicky. You can now move on to 2GC (self-calibration and target imaging), and the pipeline will split out inidvidual MS files for you target field(s). 

---

## 2GC

After completing 1GC, move on to the 2GC step. This stage performs final flagging, imaging (using WSCLEAN), and direction-independent phase self-calibration. The 2GC process produces both channelized images and a single MFS (multi-frequency synthesis) image. The MFS image maximizes sensitivity, but may be affected by bandwidth depolarization.

**To run 2GC:**

```bash
python setups/2GC.py idia
./submit_2GC_job.sh
```

**Key configurable options for 2GC include:**

- `WSC_WEIGHT = briggs 0.0` — Imaging weighting scheme (default: Briggs robustness 0.0). This maximizes sensitivity before the MeerKAT synthesized beam becomes non-Gaussian.
- `WSC_UNIFORM_IMAGE = True` — Also generates a high-angular-resolution image using `WSC_WEIGHT_HIGHRES`. Useful for tracking proper motion.
- `WSC_POL = 'IQUV'` — Stokes parameters to image. The pipeline is designed for either Stokes I or full IQUV imaging. Atypical subsets (e.g., QUV, UV) may work, but other parts of the pipeline (such as RMSYNTH) may fail.
- `WSC_IMAGE_CHANNELSOUT = 8` — Number of frequency channels to image. Increase this value for sources with high rotation measure.
- `WSC_MAX_CHANNELS = 16` — Maximum number of frequency channels imaged at once (memory limit). If `WSC_IMAGE_CHANNELSOUT > WSC_MAX_CHANNELS`, imaging is performed in steps. In this case, the MFS image is created by stacking in the image plane (qualitative only; **do not report fluxes from this MFS image**).
- `WSC_HOMOGENIZEBEAM = False` — Produces `.homogenized.fits` images with a frequency-homogenized beam using the Welzl algorithm. Also homogenizes residuals using [pypher](https://pypher.readthedocs.io/en/latest/).  
  **Experimental:** Residual homogenization can reduce RMS in the upper-frequency band, but interpret S/N with caution. In most cases, it is preferable to specify the maximum uvdist with `WSC_MAXUVL = ''` (automation for this may be added in the future).

**Self-calibration control:**  
You can adjust self-calibration solutions via the Quartical YAML files in `data/quartical`. The default file is phase-only self-calibration, which includes a frequency slope (e.g., `2GC_phase.yaml`).  
Quartical parameter documentation can be found [here](https://quartical.readthedocs.io/en/latest/). The key parameter is `time_interval='4'`, which sets the number of integration times per solution interval. By default, this is 4 (matching `oxkat`), but you may decrease it to 1 for well-behaved fields with bright point sources.

**Warning:**  
The QC routine does not account for parallactic angle rotation/de-rotation during self-calibration, since the DATA column has the parallactic angle rotation applied after 1GC. While there was initial concern this could reduce the fidelity of self-calibrated data (due to feed-frame vs. sky-frame calibration issues), in practice this effect appears to be minor—much smaller than measurement errors or other systematics. If you ever find a case where this makes a significant difference, please let me know.

Most imaging parameters for 2GC are found in `config.py` after the comment `# wsclean and 2GC defaults`, and correspond directly to [wsclean](https://wsclean.readthedocs.io/en/latest/) options.

**Output:**  
Final images will have suffixes such as `pcalmask-MFS-I-image.fits` or `pcalmask-[CHAN_NUMBER]-I-image.fits` (`I` for Stokes I). In addition to Stokes images, 2GC also produces linear polarization intensity (`Plin`, $\sqrt{Q^2+U^2}$) and total polarization intensity (`Ptot`, $\sqrt{Q^2+U^2+V^2}$) images.

---

## 3GC (Peeling)

The 3GC (Third Generation Calibration) "peeling" step is an advanced calibration and imaging process designed to further improve image fidelity by correcting for direction-dependent effects (DDEs) around bright sources in your field. Peeling is especially useful for fields with strong off-axis sources or complex extended emission, where standard direction-independent calibration is insufficient.

**How does 3GC_peel work?**

- Peeling directions are defined using DS9 region files (one region per source/direction). You can specify these manually by setting the `CAL_3GC_PEEL_REGION` variable in `config.py` to a comma-separated list of region file paths. Alternatively, you can place region files in your working directory using the naming convention `[SOURCE_NAME]_peel[PEEL_NUMBER].reg` (e.g., `CygX1_peel1.reg`, `CygX1_peel2.reg`, etc.).
- The pipeline supports an arbitrary number of peel directions, but note that each direction adds a new data column to your MS file, which can rapidly increase its size.
- Peeling involves amplitude and phase self-calibration in each direction and therefore may have weird interactions with instrumental, off-axis Leakage (which MeerKAT has a lot of). In my tests, peeling generally preserves source fidelity (including polarization fluxes), but if you notice significant changes in target flux densities—especially for polarization—please report these cases.

**What does 3GC_peel do?**

- Identifies and calibrates bright sources that may cause significant DDEs.
- Iteratively subtracts these sources from the data (the "peeling" process).
- Applies direction-dependent calibration solutions to improve the dynamic range and fidelity of the final images.
- Produces both residual and "peeled" images for scientific analysis.

**To run 3GC_peel:**

```bash
python setups/3GC_peel.py idia
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

This step performs snapshot imaging following the Heywood-ian approach, efficiently making short-timescale images by subtracting a time-averaged model and searching for image-plane variability. Run with:

```bash
python setups/SNAP.py idia
./submit_snap_job.sh
```

Key `config.py` variables (see `Snapshot imaging defaults`):

- `SNAP_FIELDS = ''` — Fields to perform snapshot imaging on (default: all fields; specify your target).
- `SNAP_INTBIN = 1` — Number of dump times per snapshot image (default: 1).
- `SNAP_INTEND = True` — If true, discards edge bins that don't fit `SNAP_INTBIN` for equispaced snapshots.
- `SNAP_Pol = True` — Default: IQUV snapshots; set to False for Stokes I only.

Other options are self-explanatory. Reach out if anything is unclear.

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

#### Note: No Polarization Angle Calibrator

polkat works without a polarization angle calibrator. As long as the primary can serve as an unpolarized leakage calibrator, you can measure total polarization, but cannot distinguish between circular and linear polarization. For most synchrotron sources (typically circularly unpolarized), total polarization is a good proxy for linear polarization. Leave `POLANG_NAME = ''` in these cases.

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
| [`RMTools`](https://github.com/CIRADA-Tools/RM-Tools?tab=readme-ov-file) | RMSynth | Rotation Measure Synthesis | Van Eck et al.|
| [`ALBUS`](https://github.com/twillis449/ALBUS_ionosphere) | RMSynth | Ionospheric RM Estimation | Willis et al.|


---

## TO DO LIST

1. Further improve documentation and include a PDF with some examples.
2. Make a Quartical only branch


