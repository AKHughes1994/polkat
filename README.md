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

This is a modified version of the MeerKAT semi-automated data processing routine [oxkat](https://github.com/IanHeywood/oxkat), enhanced to include full polarization calibration and Stokes I, Q, U, V imaging. We assume you are familiar with the oxkat workflow and its file system, and that data processing options are set by editing `oxkat/config.py`. This guide walks you through a standard use case, highlighting changes to `config.py` and new options.

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

#### INFO

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

#### 1GC

The second step, `1GC`, performs reference calibration (using calibrator fields to calibrate your target) with [casa](https://casa.nrao.edu/):

```bash
python3 setups/1GC.py idia
./submit_1GC_job.sh
```

**Important:** As of Feb 13, 2025, MeerKAT ms-files made using the SARAO archive (i.e., using the download button) have mislabelled X/Y feeds. This results in incorrect polarization properties if not corrected (see EVLA Memo 219). polkat corrects this mislabelling in its first two steps. **CAUTION:** Other (less common) archive download methods may already correct this; do not double-correct. polkat assumes you used the button. If polarization properties of diagnostic images strongly disagree with archival values (see [here](https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol)), you may have missed or double-applied this correction.

**Commonly customized `config.py` variables:**

- `CAL_1GC_DIAGNOSTICS = True` — (default: on) Images calibrators in Stokes I, Q, U, V. The primary (typically J1939) should be unpolarized; the polarization angle calibrator should match catalogue values. Leave this on for calibration checks.
- `CAL_1GC_AGGRESSIVE_FLAGS = False` — (default: off) More aggressive baseline flagging. Default flags RFI on short baselines; turning this on flags RFI on all baselines, improving polarimetry at the cost of some S/N.
- `POLANG_MOD  = [1.0, 0.0, 0.5, 0.0]` — Initialization model for the polarization angle calibrator (default: 3C286); used by `setjy` in casa. The config file also includes a model for 3C138.

**Note:** For linear feed instruments, the casa cross-hand phase solver only needs the angle quadrant approximately correct (see EVLA Memo 219). The polarization angle calibrator also acts as a polarization check source. All testing converges on the correct solution despite the input model, allowing you to estimate systematic errors by comparing measured and expected values.

If 1GC is successful, check the visibility/gain solutions. Polarization can be finicky. You can now move on to 2GC (self-calibration and target imaging).

---

#### 2GC

After 1GC, final flagging, imaging (with WSCLEAN), and direction-independent phase self-calibration are performed. Image products include channelized images and a single MFS image (maximizing sensitivity, but may cause bandwidth depolarization):

```bash
python setups/2GC.py idia
./submit_2GC_job.sh
```

There are many options; here are some key ones:

- `WSC_WEIGHT  = briggs 0.0` — Imaging weighting (default: Briggs robustness 0.0). Maximizes sensitivity before the MeerKAT synthesized beam becomes non-Gaussian.
- `WSC_UNIFORM_IMAGE = True` — Also makes a high-angular resolution image using `WSC_WEIGHT_HIGHRES`. Useful for weakly polarized point sources.
- `WSC_POL = 'IQUV'` — Stokes parameters to image (leave as is).
- `WSC_IMAGE_CHANNELSOUT = 8` — Number of frequency channels to image. Increase for sources with high rotation measure.
- `WSC_MAX_CHANNELS = 16` — Maximum frequency channels imaged at once (memory limit). If `WSC_IMAGE_CHANNELSOUT > WSC_MAX_CHANNELS`, imaging is split into steps. In this case, the MFS image is made by stacking in the image plane (qualitative only; **do not report fluxes from this MFS image**).
- `WSC_HOMOGENIZEBEAM = False` — Produces `.homogenized.fits` images with frequency-homogenized beam using the Welzl algorithm. Also homogenizes residuals using [pypher](https://pypher.readthedocs.io/en/latest/). **Experimental:** Residual homogenization reduces RMS in the upper-frequency band; interpret S/N with caution.

You can control the self-calibration solution interval with `CAL_2GC_PSOLINT=32s` (default: 32s). Decrease for bright targets.

Most `config.py` imaging parameters are after the comment `# wsclean and 2GC defaults` and link to [wsclean](https://wsclean.readthedocs.io/en/latest/) options.

Final images will have suffixes like `pcalmask-MFS-I-image.fits` or `pcalmask-[CHAN_NUMBER]-I-image.fits` (`I` for Stokes I). 2GC also produces linear polarization intensity (`Plin`, $\sqrt{Q^2+U^2}$) and total polarization intensity (`Ptot`, $\sqrt{Q^2+U^2+V^2}$) images.

---

#### SNAP

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

#### RMSYNTH

This step automates going from images to fluxes for arbitrary point sources. Please verify the fluxes make sense.

**What it does:**

1. Fits target MFS/channelized/snapshot images with user-specified Gaussians using the casa task `imfit`.
2. If `CAL_1GC_DIAGNOSTICS = True`, quantifies systematic calibration effects via image-plane analysis of calibrators.
3. Runs RM Synthesis on every Gaussian component from every target, extracting polarization angles/rotation measures.
4. Runs ALBUS to get ionospheric RM for post-processing corrections. ALBUS may fail randomly; rerun or change `RED_TYPE = RI_G03` to `RED_TYPE = RI_G01` in `config.py`.

Modify `data/rmsynth/rmsynth_info.json` as needed. Example:

```json
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

- `"image_directory"`: Directory for images to fit (`IMAGES` for 2GC, `INTERVALS` for snapshots).
- `"image_identifier"`: Identifier for images (`pcalmask` for 2GC, `restored` for snapshots).
- `"image_suffix"`: Image type (`image.fits` standard; `image.homogenized.fits` for homogenized images).
- `"image_timing"`: True if images are time-split (e.g., snapshots).
- `"source_name"`: Source names.
- `"source_ulim"`: Whether to fix component position to `"source_pos"`.
- `"rms_region"`: Manual RMS region (CASA format); otherwise, calculated as an annulus.
- `"source_pos"`: List of component positions (CASA format). Add more entries for multiple components.

You can ignore `"pos_coeff"` and `"time0"` for now. The default should work for most targeted point source observations; just change `"source_name"` and `"source_pos"`. If you didn't do snapshot imaging, remove the second entry from each list.

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


