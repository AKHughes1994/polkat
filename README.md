### Caution

This routine has been extensively tested for MeerKAT L- and S-bands. UHF band support is currently experimental. While current versions recover the correct polarization degrees, the polarization angle below 800 MHz remains problematic. Please check your data carefully when reducing UHF band observations. This documentation will be updated as issues are resolved.

---

If you make use of this software, please cite:

```bibtex
@software{2025ascl.soft02026H,
  author       = {Hughes, Andrew K. and Cowie, Fraser J. and Heywood, Ian and Hugo, Ben},
  title        = {polkat: Semi-automate full polarization of MeerKAT observations},
  howpublished = {Astrophysics Source Code Library, record ascl:2502.026},
  year         = {2025},
  month        = feb,
  eid          = {ascl:2502.026},
  adsurl       = {https://ui.adsabs.harvard.edu/abs/2025ascl.soft02026H},
  adsnote      = {Provided by the SAO/NASA Astrophysics Data System}
```

---

### What is this?

This repository contains a modified version of the MeerKAT semi-automated data processing routine [oxkat](https://github.com/IanHeywood/oxkat), extended to enable full polarization calibration and imaging of Stokes I, Q, U, V. This version uses CASA for calibration and assumes familiarity with the oxkat workflow, including configuration via `oxkat/config.py`.

This branch is distinct from `polkat_QC_selfcal`, which uses Quartical for self-calibration and includes support for 3GC (peeling). The current version documented here does not include 3GC functionality.

---

### Before We Start

polkat was designed for ILIFU clusters (hosted by IDIA), but it can also be run locally if the required software is available. All dependencies are bundled using [apptainer](https://apptainer.org/) (previously known as Singularity). Containers are available in `/software/containers` on ILIFU with filenames like `polkat-[version].sif`.

To run outside ILIFU, containers can be pulled from DockerHub:

```bash
singularity pull polkat-0.1.2.sif docker://hughesakh/polkat:0.1.2
singularity pull polkat-albus.sif docker://hughesakh/polkat_albus:fixed
```

---

### Standard Workflow

Assuming a Linux-based OS:

```bash
mkdir working_directory
cd working_directory
git clone -b polkat_casa https://github.com/AKHughes1994/polkat.git
ln -s /idia/raw/point/to/your/file.ms .
mv polkat/* .
```

Use `singularity exec` to run helper tools and ensure `config.py` is updated for your system. Look at your ms file before you begin:

```bash
singularity exec /path/to/polkat-[version].sif python3 tools/ms_info.py yourfile.ms
```

---

### INFO

Run INFO with:

```bash
python3 setups/INFO.py idia
./submit_info_job.sh
```

Replace `idia` with `node` if running locally and ensure `NODE_CONTAINER_PATH` is correctly set in `config.py`.

This step extracts observation metadata and writes `project_info.json`. It also splits fields and averages the ms-file to 1024 channels. Scans <10s (from 2s integration data) are flagged to address known metadata issues.

Key `config.py` options:

- `POLANG_NAME`: Default `'J1331+3030'` (3C286); blank disables polarization angle calibration.
- `POLANG_DIR`: Coordinates of the polang calibrator.
- `PRE_FIELDS`: List of fields (target, calibrator, etc.).



### Note: No Polarization Angle Calibrator?

If no polang calibrator is available, leave `POLANG_NAME = ''`. The primary (e.g., J1939) must be unpolarized. You can recover total polarization (but not separate linear from circular).

---

### 1GC

The second step, `1GC`, performs reference calibration (using calibrator fields to calibrate your target) with [casa](https://casa.nrao.edu/):

```bash
python3 setups/1GC.py idia
./submit_1GC_job.sh
```

**Important:** As of May 22, 2025, MeerKAT ms-files made using the SARAO archive (i.e., using the download button) have mislabelled X/Y feeds. This results in incorrect polarization properties if not corrected (see EVLA Memo 219). polkat corrects this mislabelling in its first two steps. **CAUTION:** Other (less common) archive download methods may already correct this; do not double-correct. polkat assumes you used the button. If polarization properties of diagnostic images strongly disagree with archival values (see [here](https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/pol)), you may have missed or double-applied this correction.

**Commonly customized **``** variables:**

- `CAL_1GC_DIAGNOSTICS = True` — (default: on) Images calibrators in Stokes I, Q, U, V. The primary (typically J1939) should be unpolarized; the polarization angle calibrator should match catalogue values. Leave this on for calibration checks.
- `CAL_1GC_AGGRESSIVE_FLAGS = False` — (default: off) More aggressive baseline flagging. Default flags RFI on short baselines; turning this on flags RFI on all baselines, improving polarimetry at the cost of some S/N.
- `POLANG_MOD  = [1.0, 0.0, 0.5, 0.0]` — Initialization model for the polarization angle calibrator (default: 3C286); used by `setjy` in casa. The config file also includes a model for 3C138.

**Note:** For linear feed instruments, the casa cross-hand phase solver only needs the angle quadrant approximately correct (see EVLA Memo 219). The polarization angle calibrator also acts as a polarization check source. All testing converges on the correct solution despite the input model, allowing you to estimate systematic errors by comparing measured and expected values.

If 1GC is successful, check the visibility/gain solutions. Polarization can be finicky. You can now move on to 2GC (self-calibration and target imaging).

---

### 2GC

After completing 1GC, move on to the 2GC step. This stage performs final flagging, imaging (using WSCLEAN), and direction-independent phase self-calibration. The 2GC process produces both channelized images and a single MFS (multi-frequency synthesis) image. The MFS image maximizes sensitivity, but may be affected by bandwidth depolarization.

**To run 2GC:**

```bash
python3 setups/2GC.py idia
./submit_2GC_job.sh
```

**Key configurable options for 2GC include:**

- `WSC_WEIGHT = briggs 0.0` — Imaging weighting scheme (default: Briggs robustness 0.0). This maximizes sensitivity before the MeerKAT synthesized beam becomes non-Gaussian.
- `WSC_UNIFORM_IMAGE = True` — Also generates a high-angular-resolution image using `WSC_WEIGHT_HIGHRES`. Useful for tracking proper motion.
- `WSC_POL = 'IQUV'` — Stokes parameters to image. The pipeline is designed for either Stokes I or full IQUV imaging.
- `WSC_IMAGE_CHANNELSOUT = 8` — Number of frequency channels to image. Increase this value for sources with high rotation measure.
- `WSC_MAX_CHANNELS = 16` — Maximum number of frequency channels imaged at once (memory limit). If `WSC_IMAGE_CHANNELSOUT > WSC_MAX_CHANNELS`, imaging is performed in steps. In this case, the MFS image is created by stacking in the image plane (qualitative only; **do not report fluxes from this MFS image**).
- `WSC_HOMOGENIZEBEAM = False` — Produces `.homogenized.fits` images with a frequency-homogenized beam using the Welzl algorithm. Also homogenizes residuals using [pypher](https://pypher.readthedocs.io/en/latest/).\
  **Experimental:** Residual homogenization can reduce RMS in the upper-frequency band, but interpret S/N with caution.

**Self-calibration control:**\
Unlike the `polkat_QC_selfcal` branch, this version uses CASA for self-calibration. The solution interval is set using `CAL_2GC_PSOLINT`, which controls the number of seconds over which phase solutions are solved. The default is `32s`, matching the oxkat default. Reduce this value to improve calibration for high-S/N datasets.

Most imaging parameters for 2GC are found in `config.py` after the comment `# wsclean and 2GC defaults`, and correspond directly to [wsclean](https://wsclean.readthedocs.io/en/latest/) options.

**Output:**\
Final images will have suffixes such as `pcalmask-MFS-I-image.fits` or `pcalmask-[CHAN_NUMBER]-I-image.fits` (`I` for Stokes I). In addition to Stokes images, 2GC also produces linear polarization intensity (`Plin`, sqrt(Q^2 + U^2)) and total polarization intensity (`Ptot`, sqrt(Q^2 + U^2 + V^2)) images.

---

### SNAP

This step performs snapshot imaging following the Heywood-ian approach, efficiently making short-timescale images by subtracting a time-averaged model and searching for image-plane variability. Run with:

```bash
python3 setups/SNAP.py idia
./submit_snap_job.sh
```

Key `config.py` variables (see `Snapshot imaging defaults`):

- `SNAP_FIELDS = ''` — Fields to perform snapshot imaging on (default: all fields; specify your target).
- `SNAP_INTBIN = 1` — Number of dump times per snapshot image (default: 1).
- `SNAP_INTEND = True` — If true, discards edge bins that don't fit `SNAP_INTBIN` for equispaced snapshots.
- `SNAP_Pol = True` — Default: IQUV snapshots; set to False for Stokes I only.

Once snapshot finishing running there is an optional setup:

```bash
python3 waterhole/setup_movie.py idia
```

This script will make mp4s in your `INTERVAL` directories so you can visually look for variables in your snapshot images.

---

### RMSYNTH

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

You can ignore `"pos_coeff"` and `"time0"` as these don't do anything currently. The default configuration should work for typical targeted point source observations—just update `"source_name"` and `"source_pos"` as needed. If you did not perform snapshot imaging, remove the second entry from each list.

---

### THE COMPRISING PACKAGES

polkat depends on many existing tools; the heavy lifting is done by them:

| Package                                                   | Stage    | Purpose                 | Reference                                                                                                                                                    |
| --------------------------------------------------------- | -------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| [`astropy`](https://www.astropy.org/)                     | 1GC, 2GC | FITS, coordinates       | [2013A&A...558A..33A](https://ui.adsabs.harvard.edu/abs/2013A%26A...558A..33A), [2018AJ....156..123A](https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A) |
| [`CASA`](https://casa.nrao.edu/)                          | 1GC, 2GC | Calibration, flagging   | [2022PASP..134k4501C](https://ui.adsabs.harvard.edu/abs/2022PASP..134k4501C)                                                                                 |
| [`Quartical`](https://github.com/ratt-ru/CubiCal)         | 2GC      | Self-cal                | [2024arXiv241210072K](https://ui.adsabs.harvard.edu/abs/2024arXiv241210072K)                                                                                 |
| [`owlcat`](https://github.com/ska-sa/owlcat/)             | 2GC      | FITS manipulation       | -                                                                                                                                                            |
| [`ragavi`](https://github.com/ratt-ru/ragavi/)            | 1GC, 2GC | Gain plotting           | -                                                                                                                                                            |
| [`shadeMS`](https://github.com/ratt-ru/shadeMS/)          | 1GC      | Visibility plotting     | [2022ASPC..532..385S](https://ui.adsabs.harvard.edu/abs/2022ASPC..532..385S)                                                                                 |
| [`tricolour`](https://github.com/ska-sa/tricolour)        | 2GC      | Flagging                | [2022arXiv220609179H](https://ui.adsabs.harvard.edu/abs/2022arXiv220609179H)                                                                                 |
| [`wsclean`](https://gitlab.com/aroffringa/wsclean)        | 2GC      | Imaging                 | [2014MNRAS.444..606O](https://ui.adsabs.harvard.edu/abs/2014MNRAS.444..606O)                                                                                 |
| [`pypher`](https://pypher.readthedocs.io/en/latest/)      | 2GC      | Residual homogenization | [2016A&A...596A..63B](https://ui.adsabs.harvard.edu/abs/2016A%26A...596A..63B)                                                                               |
| [`RMTools`](https://github.com/CIRADA-Tools/RM-Tools)     | RMSYNTH  | RM synthesis            | [2020ascl.soft05003P](https://ui.adsabs.harvard.edu/abs/2020ascl.soft05003P)                                                                                 |
| [`ALBUS`](https://github.com/twillis449/ALBUS_ionosphere) | RMSYNTH  | Ionosphere RM           | Willis et al.                                                                                                                                                |

---

### TODO

1. Expand snapshot imaging options.
2. Assess beam homogenization post-processing.

