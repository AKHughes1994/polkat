> **NOTE:** Following a large merge of changes from `dev`, `main` (this branch) is now the default branch going forward. `polkat_QC_selfcal` is no longer the default — it is preserved as-is for posterity.

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


### Latest modification (July, 2026)

#### Update: Image-Plane UV Cuts for Short Baselines (X-KAT Tuning)

The default `config.py` is now tuned for X-KAT observations, where preserving spectral fidelity across the band is the priority. Short baselines are now excluded via an **image-plane** UV cut (`WSC_BASELINE_CUT = True`, `WSC_BASELINE_CUTLENGTH = '750m'` by default, `WSC_TAPERMASK = True`) rather than by flagging the data. This addresses the same short-baseline artefacts that motivate `CAL_1GC_BL_MODE`, but without discarding any visibilities from the MS — the cut is only applied at imaging/deconvolution time, so you can always re-image the same MS with the full baseline range afterward (e.g. by setting `WSC_BASELINE_CUT = False`) to search for extended structure.

# **IF YOU WANT THE FULL ARRAY / EXTENDED STRUCTURE, MAKE SURE TO SET `WSC_BASELINE_CUT = False`, OR RUN A SEPARATE IMAGING PASS WITHOUT THE CUT AFTER 2GC / YOUR DESIRED SELF-CALIBRATION — OTHERWISE SHORT BASELINES ARE SILENTLY EXCLUDED BY DEFAULT.**

---

#### Update: Default Channels-Out Increased to 64; Joint-Polarization Cleaning Removed

`WSC_PCAL_CHANNELSOUT` (and the calibrator / data-masked channel counts, which follow it by default) now default to **64** channels, up from 8, and `WSC_JOINPOLARIZATIONS` now defaults to `False`. Sixty-four channels is about the highest channelisation a standard IDIA node can image without either running into memory errors or having to split the channelisation into blocks (`WSC_MAX_CHANNELS`), which loses the MFS (continuum) image. Please read through `config.py` — particularly the 2GC section — to see the full set of changes, and raise a discussion if anything is unclear on why a particular default was chosen.

---

#### Update: New `extra/` Directory for Advanced Imaging Routines

A new [`extra/`](extra/) directory has been added for advanced/experimental imaging scripts that sit outside the standard pipeline stages. The two-stage self-calibration script has moved here, from `waterhole/setup_2GC_twostage.py` to [`extra/setup_2GC_twostage.py`](extra/setup_2GC_twostage.py) — update any submission scripts that reference the old path. `extra/` will be populated with further routines over time as they mature enough for general use, and is intended to serve advanced users who want to go beyond the standard workflow.

---

#### Update: New `extra/setup_extended_imaging.py` for Full-UV-Range Extended Imaging

[`extra/setup_extended_imaging.py`](extra/setup_extended_imaging.py) jointly images an arbitrary set of measurement sets (given as a glob pattern or a `.txt` list) at the **full UV range** — it forces mf-weighting on and forces the inner Tukey taper / `minuv-l` / `maxuv-l` cuts off, overriding `config.py` regardless of the short-baseline defaults described in the 2GC section below. It performs no calibration: blind image → automask → final masked image (or final image only, if an existing FITS mask is supplied via `-m`). This is intended for recovering extended/diffuse structure that the default image-plane uv-cut and taper are designed to suppress during standard 2GC point-source imaging.

---

#### Update: Bring Your Own Xf Table (`XF_OVERRIDE_TABLE`)

For observations without a dedicated cross-hand calibrator, `XF_OVERRIDE_TABLE` lets you supply a pre-solved cross-hand phase (Xf) table to use instead of solving one from the current dataset. Cross-hand phase appears to be very stable over timescales of months, provided the same reference antenna is used, so a table solved on a different observation with the same refant can often be reused directly. A helper function that builds an override table from a set of input tables, or from a text file of phase/frequency pairings, is planned. This is intended for **advanced users** — open a discussion on the repository if you'd like to see how this might work for your use case.

---

#### Update: Ionospheric RM — Migrated from ALBUS to Spinifex

Ionospheric RM estimation has been migrated from ALBUS to [Spinifex](https://spinifex.readthedocs.io/en/latest/), which is now bundled in the main container. The separate `polkat-albus.sif` container is no longer required.

---

#### Update: New 1GC Config Parameters (`CAL_1GC_BL_MODE`, `CAL_1GC_APPLYPARANG`)

**`CAL_1GC_BL_MODE`** — a new flagging mode during the deterministic flagging step of 1GC, replacing the older `CAL_1GC_AGGRESSIVE_FLAGS` boolean. The available options are:

- `’freq’` *(default)* — flags specific RFI-affected frequency bands only on short baselines (< 600 m). Recommended for most observations.
- `’none’` — disables all baseline-length-dependent flagging entirely.
- `’aggressive’` — flags *all* data on the affected short baselines, regardless of frequency. This was found to remove discontinuities in the Stokes I spectrum arising from dynamic-range-limited data. If using `’aggressive’`, also consider disabling multi-frequency weighting by setting `WSC_MFWEIGHT = False` in `config.py`.

**`CAL_1GC_APPLYPARANG`** (`True` by default) — when `True`, the parallactic angle correction is applied to `CORRECTED_DATA` before splitting to the target MS, placing the data in the sky frame and meaning self-calibration corrections are also solved in the sky frame. When `False`, the correction is instead applied after the final self-calibration step, keeping the solve in the feed frame — useful for amplitude self-calibration of the off-diagonal Jones terms, but with caveats for longer tracks. See the 1GC section below, under "Parallactic angle and self-calibration frame", for the full details and the sky-model workaround needed for longer tracks.

---

#### Update: Two-Stage Self-Calibration and wsclean 3.5

An alternative 2GC script, [`extra/setup_2GC_twostage.py`](extra/setup_2GC_twostage.py), performs two rounds of self-calibration: a first pass on high-significance pixels only (controlled by `WSC_SHALLOWMASK`), followed by amplitude self-calibration on an intermediate image, then a final imaging pass. wsclean has been upgraded to **version 3.5**, adding `--local-rms-strength`; every `localrms`-style parameter (e.g. `WSC_INTER_LOCALRMS`, `WSC_SHALLOWMASK_LOCALRMS`) is now a single knob for both on/off and strength -- `False` disables it, `True` enables it at the default strength (0.5), a float enables it at that strength. Intermediate masking parameters (`WSC_INTER_AUTOMASK`, `WSC_INTER_AUTOTHRESHOLD`, `WSC_INTER_LOCALRMS`) apply only to this workflow. See the 2GC section below for full parameter listings.

**Amplitude self-calibration** YAML files have been added to `data/quartical/`. [`2GC_complex.yaml`](data/quartical/2GC_complex.yaml) performs direction-independent amplitude self-calibration; [`2GC_complex_2dir.yaml`](data/quartical/2GC_complex_2dir.yaml) extends this to two directions and is used to peel a problematic in-field source. Place a DS9 region file named `DIR.reg` in your working directory to activate the peeling. This is intended for **advanced users**; inspect your solutions carefully.

---

#### Update: New Config Parameter (`MOD_MODEL_SELFCAL`)

`MOD_MODEL_SELFCAL` (`False` by default) zeros Stokes V in the sky model before each predict step in the two-stage self-calibration workflow. This is designed for polarization observations of dynamic-range-limited fields (e.g., X-ray binaries) where the model could otherwise accumulate a spurious non-zero Stokes V component. Has no effect in the standard single-round 2GC workflow.

---

### What is this?

This repository contains a modified version of the MeerKAT semi-automated data processing routine [oxkat](https://github.com/IanHeywood/oxkat), enhanced to support full polarization calibration and Stokes I, Q, U, V imaging. It is assumed that you are already familiar with the oxkat workflow and file system, and that you configure data processing options by editing `oxkat/config.py`. This guide walks you through a standard use case, highlighting changes to `config.py` and introducing new options.

**IMPORTANT: For time series or multi-epoch campaigns**, it is strongly recommended that you download a local branch and use it consistently throughout your project. This repository is actively maintained, and improvements to flux density calibration or other processing steps may introduce global offsets or systematic differences between epochs processed with different versions. To ensure consistency across your dataset, fix your version at the start of your campaign.

---

#### Before We Start

This routine is designed primarily for use on the ILIFU clusters operated by The Inter-university Institute for Data Intensive Astronomy (IDIA), but you can run it locally if you have the required software. The necessary software is bundled into containers using [apptainer](https://apptainer.org/) (formerly singularity). Containers are available in `/software/containers` on ILIFU and follow the naming convention `polkat-[version].sif`. If you do not have ILIFU access but want to use polkat, you can get container access via pulling from dockerhub:

**Note:** the image bundles CASA, wsclean, QuartiCal, Spinifex, and tricolour, so it is several GB, and `singularity pull` doesn't just download it — it also has to unpack and re-squash every layer into a single SIF file, which is CPU/I/O-bound rather than network-bound. Do this from an interactive/compute-node session with local scratch space rather than a login node with a networked home directory, or it can take a very long time.

```
# The main container 
singularity pull polkat-0.2.5.sif docker://hughesakh/polkat:0.2.5

# If you need to point to particular storage areas for build (i.e., can't use default home/tmp)
mkdir -p .singularity_tmp .singularity_cache
SINGULARITY_TMPDIR="$PWD/.singularity_tmp" \
SINGULARITY_CACHEDIR="$PWD/.singularity_cache" \
singularity pull "$PWD/polkat-0.2.5.sif" docker://hughesakh/polkat:0.2.5

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
   git clone https://github.com/AKHughes1994/polkat.git
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
CAL_1GC_APPLYPARANG = True  # Apply parallactic angle correction after 1GC (sky-frame self-cal).
                             # Set False for feed-frame amplitude self-calibration; recommended 
                             # only for short snapshots (≲1 hour). See 2GC notes for full details.
CAL_1GC_BL_MODE = 'freq'    # Baseline-dependent flagging mode during deterministic flags.
                             # 'freq'       (default) — flags RFI-affected bands on baselines < 600 m only.
                             # 'none'                 — disables all baseline-dependent flagging.
                             # 'aggressive'           — flags all data on affected short baselines.
                             #   The 'aggressive' mode removes Stokes I spectral discontinuities
                             #   in dynamic-range-limited data; recommended for in-band spectral
                             #   work alongside WSC_MFWEIGHT = False.
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

By default, 2GC applies an image-plane UV cut and taper to suppress short-baseline artefacts (`WSC_BASELINE_CUT = True`, `WSC_TAPERMASK = True` — see "Image-Plane UV Cuts for Short Baselines" above) — X-KAT observations are mostly focused on point sources, so this is the appropriate default. If you instead need to recover extended/diffuse structure, re-image with `WSC_BASELINE_CUT = False`, or use [`extra/setup_extended_imaging.py`](extra/setup_extended_imaging.py) to jointly image across the full UV range without re-running calibration.

**To run 2GC:**

```bash
python3 setups/2GC.py idia
./submit_2GC_job.sh
```

For fields with clear artefacts around bright sources (possibly due to calibration errors or residual RFI), an alternative two-stage self-calibration approach is available:

```bash
python3 extra/setup_2GC_twostage.py idia
./submit_2GC_twostage_job.sh
```

This performs two rounds of self-calibration: a first pass restricted to high-significance pixels (controlled by `WSC_SHALLOWMASK`), followed by amplitude self-calibration on an intermediate image, before proceeding to the final imaging pass. wsclean has been upgraded to **version 3.5**, which introduces the `--local-rms-strength` option. Every `localrms`-style parameter is now a single knob for both on/off and strength: `False` disables it, `True` enables it at the default strength (0.5), and a float enables it at that strength directly. The relevant `config.py` parameters are:

```python
# Intermediate image (two-stage self-calibration only)
WSC_INTER_LOCALRMS      = 0.33   # Local RMS map for intermediate image masking; local-rms strength
WSC_INTER_AUTOMASK      = 3.0    # Auto-mask threshold (sigma) for intermediate image
WSC_INTER_AUTOTHRESHOLD = 1.0    # Auto-threshold (sigma) for intermediate image

WSC_SHALLOWMASK          = 30.0  # Initial auto-mask threshold (sigma) for the first deconvolution pass;
                                  # in the two-stage workflow, calibrators use 2× this value
WSC_SHALLOWMASK_LOCALRMS = 0.5   # Local RMS map for the first deconvolution pass; local-rms strength
```

**Amplitude self-calibration** YAML files are available in `data/quartical/`. [`2GC_complex.yaml`](data/quartical/2GC_complex.yaml) performs direction-independent amplitude self-calibration; [`2GC_complex_2dir.yaml`](data/quartical/2GC_complex_2dir.yaml) extends this to two directions to peel a problematic in-field source. To use the peeling configuration, create a DS9 region file named `DIR.reg` in your working directory specifying the source. This is intended for **advanced users**.

A new `config.py` parameter, `MOD_MODEL_SELFCAL`, controls whether Stokes V is zeroed in the sky model before each predict step in the two-stage workflow. This is designed for polarization observations of dynamic-range-limited fields (e.g., X-ray binaries), where the model could otherwise accumulate a spurious non-zero Stokes V component:

```python
MOD_MODEL_SELFCAL = False  # If True, zeros Stokes V in the self-calibration sky model before 
                            # each predict step (two-stage workflow only). Recommended for  
                            # polarization observations of dynamic-range-limited sources.
```

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

# Masking: blind image (mask generation only)
WSC_AUTOMASK_BLIND      = 5.0   # Auto-mask threshold (sigma) for blind image
WSC_AUTOTHRESHOLD_BLIND = 1.0   # Auto-threshold (sigma) at which CLEAN stops for blind image
WSC_THRESHOLD_BLIND     = False  # Absolute flux threshold for blind image; False = use auto-threshold
WSC_LOCALRMS_BLIND      = 0.5    # Local RMS map for masking in the blind image; local-rms strength

# Masking: final image
# Standard 2GC: applied to both self-calibration imaging passes.
# Two-stage 2GC: applied to the final imaging pass only.
# (For two-stage shallow/intermediate masking params, see WSC_SHALLOWMASK / WSC_INTER_* above.)
WSC_MASK          = False   # FITS mask to use; False = no mask (blind deconvolution)
WSC_THRESHOLD     = False   # Absolute flux threshold (Jy) to stop cleaning; False = use auto
WSC_AUTOMASK      = 2.0     # Auto-mask threshold (sigma) for main deconvolution
WSC_AUTOTHRESHOLD = 0.5     # Auto-threshold (sigma) at which CLEAN stops
WSC_LOCALRMS      = False   # Local RMS map for adaptive masking: False disables, True enables at
                             # the default strength (0.5), a float enables at that strength
```


**Multi-frequency weighting (IMPORTANT):**

This pipeline supports spectrally resolved functionality (e.g., Rotation Measure synthesis). The `WSC_MFWEIGHT` parameter controls whether multi-frequency weighting is applied during imaging. It is set to `True` by default, optimized for a single combined MFS image with the lowest RMS noise. **If you plan to use spectral functionality, set this to `False` to avoid introducing artificial spectral structure** (see [wsclean MF weighting documentation](https://wsclean.readthedocs.io/en/latest/mf_weighting.html)).

**Self-calibration control:**

Self-calibration solutions are controlled via Quartical YAML files in `data/quartical`. The default configuration (`2GC_phase.yaml`) performs phase-only self-calibration with a frequency-dependent slope. See the [Quartical documentation](https://quartical.readthedocs.io/en/latest/) for details. The key parameter is `time_interval='1'`, which sets the solution interval to 1 integration time by default.

**Parallactic angle and self-calibration frame (`CAL_1GC_APPLYPARANG`):**

```python
CAL_1GC_APPLYPARANG = True  # When True, the parallactic angle correction is applied to 
                             # CORRECTED_DATA before splitting to the target MS (sky-frame 
                             # self-calibration). When False, the correction is applied after 
                             # the final self-calibration step (feed-frame self-calibration).
```

When `True` (default), `CORRECTED_DATA` is in the sky frame before self-calibration, meaning the diagonal (XX/YY) and off-diagonal (XY/YX) Jones terms may be mixed by parallactic angle rotation. For phase-only self-calibration this has negligible effect. If you are performing **amplitude self-calibration — especially of the off-diagonal Jones terms** — consider setting this to `False` to keep the solve in the feed frame, with the parallactic angle correction applied after the final self-calibration step instead (`output.apply_p_jones_inv=true` in the Quartical YAML, applied automatically by `setups/2GC.py` whenever `CAL_1GC_APPLYPARANG = False`). This is recommended only for **short-track observations (≲1 hour)**, where parallactic angle rotation over the track is usually small enough that Q and U are not significantly mixed — but verify this for your own field rather than treating it as a blanket rule, since a source transiting near zenith can accumulate parallactic angle quickly even in under an hour, depending on declination. For longer tracks where phase-only calibration is sufficient, leave this `True`.

If you have a longer track and need feed-frame amplitude self-calibration anyway, you cannot just image/time-average the raw feed-frame data to get a model — the source's polarization angle is fixed in the sky frame but appears to rotate over the track in the feed frame, so averaging feed-frame data over a long track smears the polarized signal, the same way bandwidth-smearing washes out signal across frequency. The model has to come from the sky frame instead:

1. Parang-correct a copy of the data and split it out separately; image/self-calibrate that copy to build a clean sky-frame model (safe to time-average, since polarization angle is static in the sky frame for a non-variable source).
2. Predict that sky-frame model onto the dataset you'll actually run the feed-frame amplitude self-calibration on (make sure this dataset has *not* already had a parallactic angle correction applied — it needs to genuinely still be in the feed frame, or the forward/inverse P-Jones steps below will be solving against the wrong frame).
3. Calibrate with **both** `input_model.apply_p_jones=true` (forward rotation, applied to the model per-timestep at solve time, projecting the static sky-frame model into the instantaneous feed frame to match the uncorrected visibilities) and `output.apply_p_jones_inv=true` (inverse rotation, applied once at the end to bring `CORRECTED_DATA` back into the sky frame). This per-timestep rotation is deterministic and doesn't smear anything — smearing only happens if the model itself was built by averaging in the wrong frame.

Watch for **double-correction** too: if `CAL_1GC_APPLYPARANG=True` (data already de-rotated into the sky frame at 1GC) but `output.apply_p_jones_inv` is still left `True`, the inverse rotation gets applied a second time to data that's already sky-frame.

Finally, if the source is intrinsically variable over the track (in flux, polarization, or both — not just apparently rotating due to parallactic angle, which the forward/inverse P-Jones steps above already handle exactly), a single static model — even a correctly-built sky-frame one — is not enough; you'll need a time-resolved model as well. Between the double-correction and smearing failure modes and the extra model-building work, this whole path warrants care and should only be attempted if you're confident in your approach. Feel free to raise a discussion on the repository if this applies to your use case.

The corresponding Quartical YAML controls for the parallactic angle are:

```yaml
input_model:
  apply_p_jones: false      # Apply P-Jones rotation to the model before solving.
                             # Set True when the data is in the feed frame (CAL_1GC_APPLYPARANG = False)
                             # so that the model is rotated to match.
output:
  apply_p_jones_inv: true   # Apply the inverse P-Jones rotation when writing CORRECTED_DATA.
                             # Set True to rotate the output back into the sky frame.
```

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
4. Runs Spinifex to estimate the ionospheric RM for post-processing corrections.

You may need to modify `data/rmsynth/rmsynth_info.json` to match your dataset. Example configuration:

```json
{
    "image_directory": ["IMAGES", "INTERVALS"],
    "image_identifier": ["pcalmask", "restored"],
    "image_suffix": ["image.fits", "image.fits"],
    "image_timing": [false, true],
    "source_name": ["SwiftJ1727", "SwiftJ1727"],
    "label": ["", ""],
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
- `"label"`: Optional output name override, one per entry. Leave as `""` to use `"source_name"` as-is (the default, unchanged behaviour). If you're fitting multiple components/entries against the same field/`"source_name"` — e.g. separate jet ejecta or repeated fits with different `"source_pos"` — give each entry its own `"label"` so their outputs (JSON, plots) get distinct filenames instead of overwriting the previous entry's results. `"source_name"` itself still drives image discovery and MS timing lookups either way, so `"label"` never affects which input files are found.
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
| [`Spinifex`](https://spinifex.readthedocs.io/en/latest/) | RMSynth | Ionospheric RM Estimation | - |


---

## TO DO LIST

1. A dedicated processing guide for the X-KAT/ThunderKAT collaboration is currently in preparation and will be appended to this repository when complete.
2. Further improve documentation and include a PDF with worked examples.
3. Make a Quartical-only branch.


