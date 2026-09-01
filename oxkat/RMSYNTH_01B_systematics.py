# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import os
import re
import sys
import json
import time
import numpy as np
import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from oxkat import config as cfg


# =============================================================================
# Global fitting options
# =============================================================================

# Set True to always re-run the full imfit loop even if a polarization JSON
# already exists for the pol-angle calibrator.  Set False to skip fitting and
# load the existing JSON directly (useful for re-running main() after editing
# the downstream appending logic without redoing slow imfit calls).
OVERWRITE = True

# =============================================================================


def msg(txt):
    """Print a timestamped log message to stdout."""
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


# =============================================================================
# Polarization utilities
# =============================================================================

def calculate_P0(flux_P, rms_Q, rms_U, Aq=0.8):
    """
    Calculate the de-biased linearly polarized flux density.

    Noise combination follows Hales et al. (2012), doi:10.48550/arXiv.1205.5310.
    De-biasing is applied when S/N >= 3, following Vaillancourt (2006),
    doi:10.48550/arXiv.astro-ph/0603110.

    Parameters
    ----------
    flux_P : float  -- measured P = sqrt(Q^2 + U^2), any flux unit
    rms_Q  : float  -- RMS noise in Stokes Q, same unit
    rms_U  : float  -- RMS noise in Stokes U, same unit
    Aq     : float  -- noise-ratio coefficient for the dominant Stokes (default 0.8)

    Returns
    -------
    flux_P0 : float  -- de-biased linear polarization flux density
    rms_P   : float  -- combined RMS noise for linear polarization
    """
    if rms_Q >= rms_U:
        A, B = Aq, 1.0 - Aq
    else:
        A, B = 1.0 - Aq, Aq

    rms_P   = (A * rms_Q**2 + B * rms_U**2)**0.5
    flux_P0 = (flux_P**2 - rms_P**2)**0.5 if flux_P / rms_P >= 3.0 else flux_P

    return flux_P0, rms_P


def return_max(im, region):
    """
    Return the signed peak value with the larger absolute magnitude from a region.

    Necessary for Stokes parameters that are not positive-definite (Q, U, V).

    Parameters
    ----------
    im     : str  -- CASA image path
    region : str  -- CASA region string

    Returns
    -------
    float : signed peak (positive or negative) with the largest |value|
    """
    ims  = imstat(im, region=region)
    fmax = ims['max'][0]
    fmin = ims['min'][0]
    return fmax if fmax > abs(fmin) else fmin


# =============================================================================
# CASA imfit wrappers
# =============================================================================

def get_imstat_values(image, x, y, n_beams=4.0):
    """
    Extract peak flux, peak pixel location, and local RMS from a CASA image.

    The signal region is a circle of radius n_beams * BMAJ centred on (x, y).
    The noise is measured in a surrounding annulus enclosing ~500 beam areas.

    Parameters
    ----------
    image   : str   -- CASA image path
    x       : float -- RA pixel coordinate
    y       : float -- Dec pixel coordinate
    n_beams : float -- signal region radius in units of BMAJ (default 4)

    Returns
    -------
    [flux, xpix, ypix, rms] : list of float
        flux  -- signed peak flux density in image units
        xpix  -- RA pixel of the peak
        ypix  -- Dec pixel of the peak
        rms   -- RMS noise in the surrounding annulus
    """
    bmaj = imhead(image, mode='get', hdkey='BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey='BMIN')['value']

    r_in       = n_beams * bmaj
    r_out      = np.sqrt(500 * 0.25 * bmaj * bmin + r_in**2)
    region     = f'circle[[{x}pix,{y}pix],{r_in}arcsec]'
    rms_region = f'annulus[[{x}pix,{y}pix],[{r_in}arcsec,{r_out}arcsec]]'

    ims  = imstat(image, region=region)
    flux = return_max(image, region)
    xpix = ims['maxpos'][0]
    ypix = ims['maxpos'][1]
    rms  = imstat(image, region=rms_region)['rms'][0]

    return [flux, xpix, ypix, rms]


def make_estimate(fname, image, x, y, fix_var='abp'):
    """
    Write a CASA imfit estimate file for a single point-source component.

    The amplitude seed is the signed peak within one BMAJ of (x, y).
    The Gaussian shape is initialised to the restoring beam.

    Parameters
    ----------
    fname   : str   -- output estimate file path (read by imfit)
    image   : str   -- CASA image from which to read beam and peak flux
    x       : float -- RA pixel coordinate
    y       : float -- Dec pixel coordinate
    fix_var : str   -- CASA constraint string:
                       'abp'   -- fit amplitude, beam axes, PA; position free
                       'xyabp' -- additionally freeze the sky position

    Returns
    -------
    0 : int  (success indicator)
    """
    bmaj = imhead(image, mode='get', hdkey='BMAJ')['value']
    bmin = imhead(image, mode='get', hdkey='BMIN')['value']
    bpa  = imhead(image, mode='get', hdkey='BPA')['value']

    region     = f'circle[[{x}pix,{y}pix],{bmaj}arcsec]'
    flux_guess = return_max(image, region)

    with open(fname, 'w') as f:
        f.write(f'{flux_guess},{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, {fix_var}\n')

    return 0


def get_imfit_values(fname, image, x, y, n_beams=5.0):
    """
    Run CASA imfit on a single-component source in a circular fit region.

    Parameters
    ----------
    fname   : str   -- estimate file produced by make_estimate()
    image   : str   -- CASA image path to fit
    x       : float -- RA pixel coordinate (centre of the fit region)
    y       : float -- Dec pixel coordinate (centre of the fit region)
    n_beams : float -- fit region radius in units of BMAJ (default 5)

    Returns
    -------
    dict : full CASA imfit result dictionary
    """
    bmaj   = imhead(image, mode='get', hdkey='BMAJ')['value']
    r      = n_beams * bmaj
    region = f'circle[[{x}pix,{y}pix],{r}arcsec]'
    return imfit(image, estimates=fname, region=region)


def check_position(fname, image, x, y, snr_thresh=5.0):
    """
    Decide whether the source position should be freed or frozen, then write
    the imfit estimate file.

    If the absolute peak flux at (x, y) exceeds snr_thresh * local_rms, the
    position is freed ('abp'); otherwise it is fixed ('xyabp').

    Parameters
    ----------
    fname      : str   -- output estimate file path (passed to imfit)
    image      : str   -- CASA image in which to test the S/N
    x          : float -- RA pixel coordinate
    y          : float -- Dec pixel coordinate
    snr_thresh : float -- minimum S/N ratio for a free position fit (default 5)
    """
    bmaj = imhead(image, mode='get', hdkey='bmaj')['value']
    bmin = imhead(image, mode='get', hdkey='bmin')['value']
    bpa  = imhead(image, mode='get', hdkey='bpa')['value']

    check_file = f'check_pos_{os.path.basename(fname)}.txt'
    with open(check_file, 'w') as cf:
        cf.write(f'0.0,{x},{y},{bmaj}arcsec,{bmin}arcsec,{bpa}deg, xyabp')

    region    = f'circle[[{x}pix,{y}pix],{2 * bmaj}arcsec]'
    test_flux = abs(imfit(image, region=region,
                          estimates=check_file)['results']['component0']['peak']['value'])
    rms = get_imstat_values(image, x, y)[3]
    snr = test_flux / rms if rms > 0 else 0.0

    if snr >= snr_thresh:
        fix_var = 'abp'
        msg(f'  check_position: FREE  | S/N = {snr:.1f} (>= {snr_thresh}) '
            f'| {os.path.basename(image)}')
    else:
        fix_var = 'xyabp'
        msg(f'  check_position: FIXED | S/N = {snr:.1f} (<  {snr_thresh}) '
            f'| {os.path.basename(image)}')

    make_estimate(fname, image, x, y, fix_var)


# =============================================================================
# Image discovery
# =============================================================================

def find_cal_images(cal_name, stokes='I'):
    """
    Locate calibrator MFS diagnostic FITS images, prioritising per-scan images.

    Search order:
      1. Images containing 'scan' in the filename  -->  preferred
      2. All other matching images                  -->  fallback only

    If both types exist, scan images are used exclusively and the non-scan
    images are explicitly logged and ignored.

    Parameters
    ----------
    cal_name : str  -- calibrator field name as it appears in image filenames
    stokes   : str  -- Stokes label to search for (e.g. 'I', 'Q', 'V')

    Returns
    -------
    prefixes      : list of str  -- sorted image prefixes (everything before '-MFS-'),
                                    one per scan or image set found
    used_fallback : bool         -- True if the non-scan fallback path was taken
    """
    scan_pattern    = cfg.IMAGES + f'/*{cal_name}*scan*diagnostic-MFS-{stokes}-image.fits'
    nonscan_pattern = cfg.IMAGES + f'/*{cal_name}*diagnostic-MFS-{stokes}-image.fits'

    scan_images    = sorted(glob.glob(scan_pattern))
    nonscan_images = sorted(glob.glob(nonscan_pattern))
    used_fallback  = False

    if scan_images:
        non_scan_only = [im for im in nonscan_images if im not in scan_images]

        msg(f'Found {len(scan_images)} scan image(s) for {cal_name} Stokes {stokes} '
            f'-- scan images take priority')
        if non_scan_only:
            msg(f'  Ignoring {len(non_scan_only)} non-scan prefix image(s):')
            for im in non_scan_only:
                msg(f'    Ignored: {os.path.basename(im)}')

        prefixes = sorted(list(set([im.split('-MFS-')[0] for im in scan_images])))

    elif nonscan_images:
        msg(f'WARNING: No scan images found for {cal_name} Stokes {stokes} '
            f'-- falling back to {len(nonscan_images)} non-scan image(s)')
        prefixes      = sorted(list(set([im.split('-MFS-')[0] for im in nonscan_images])))
        used_fallback = True

    else:
        msg(f'No images found for {cal_name} Stokes {stokes} '
            f'(tried scan and non-scan patterns)')
        prefixes = []

    return prefixes, used_fallback


def cal_json_path(cal_name, first_prefix):
    """
    Build the per-calibrator output JSON path.

    Strips any scan index token (e.g. '_scan001') from the prefix so the JSON
    is named per-calibrator rather than per-scan.  For example:

        .../RESULTS/1727_J1939-6342_scan001_polarization.json  (WRONG)
        .../RESULTS/1727_J1939-6342_polarization.json          (correct)

    Parameters
    ----------
    cal_name     : str  -- calibrator name (used only for the log message)
    first_prefix : str  -- the first scan prefix string from find_cal_images()

    Returns
    -------
    str : absolute path to the output JSON file inside cfg.RESULTS
    """
    # Strip everything from -MFS- onward, then remove any _scanNNN token
    base = first_prefix.split('-MFS-')[0]
    base = re.sub(r'_scan\d+', '', base)
    base = base.replace(cfg.IMAGES, cfg.RESULTS)
    return f'{base}_polarization.json'


# =============================================================================
# Plotting helpers
# =============================================================================

def _mad_ylim(ax, flux, rms, floor_lo=None):
    """
    Set y-axis limits using MAD-based 5-sigma outlier rejection.

    The limits are driven by the actual data range (min(flux - rms),
    max(flux + rms)) of inlier points, with 15% padding.  Symmetry around
    zero is NOT enforced -- if the data run from -2 to 20000, the axis
    will be roughly -2 to 20000, not -20000 to 20000.

    Parameters
    ----------
    ax       : matplotlib Axes
    flux     : array-like of float  -- central values
    rms      : array-like of float  -- per-point uncertainties
    floor_lo : float or None        -- if given, clamp the lower limit to this value
                                       (e.g. floor_lo=0.0 for Stokes I)
    """
    flux = np.asarray(flux, dtype=float)
    rms  = np.asarray(rms,  dtype=float)

    med   = np.nanmedian(flux)
    mad   = np.nanmedian(np.abs(flux - med))
    sigma = 1.4826 * mad
    inlier = (np.abs(flux - med) <= 5.0 * sigma) if sigma > 0 else np.ones(len(flux), dtype=bool)

    lo = np.nanmin(flux[inlier] - rms[inlier])
    hi = np.nanmax(flux[inlier] + rms[inlier])

    if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
        span = hi - lo
        pad  = 0.15 * span if span > 0 else 0.5
        lo   = max(floor_lo, lo - pad) if floor_lo is not None else lo - pad
        ax.set_ylim(lo, hi + pad)


def plot_cal_spectrum(output_dictionary, cal_name, scan_key):
    """
    Plot per-channel Stokes I, Q, U, V spectra for one calibrator scan.

    Four panels stacked vertically, sharing the frequency axis:
      Panel 1 (top): Stokes I -- floor at zero; MFS dashed reference;
                     spectral index power-law fit when S/N > 30 and >= 3 valid channels.
      Panels 2-4:   Stokes Q, U, V -- non-positive-definite; data-driven y-range
                     via MAD-based outlier rejection; no enforced symmetry; MFS dashed.

    Parameters
    ----------
    output_dictionary : dict  -- polarization dictionary from get_polcal_polarization()
    cal_name          : str   -- calibrator name for title and filename
    scan_key          : str   -- scan key, e.g. 'scan000'
    """
    plot_dir = os.path.join('RESULTS', 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)

    chan = output_dictionary['CHAN'].get(scan_key, {})
    mfs  = output_dictionary['MFS'].get(scan_key, {})

    if not chan or not chan.get('freq_GHz'):
        msg(f'  plot_cal_spectrum: no channel data for {cal_name} {scan_key} -- skipping')
        return

    freq_arr = np.array(chan['freq_GHz'])
    I_flux   = np.array(chan['I_flux_mJy'])
    I_rms    = np.array(chan['I_rms_mJy'])
    Q_flux   = np.array(chan['Q_flux_mJy'])
    Q_rms    = np.array(chan['Q_rms_mJy'])
    U_flux   = np.array(chan['U_flux_mJy'])
    U_rms    = np.array(chan['U_rms_mJy'])
    V_flux   = np.array(chan['V_flux_mJy'])
    V_rms    = np.array(chan['V_rms_mJy'])

    mfs_I    = mfs.get('I_flux_mJy', 0.0)
    mfs_Q    = mfs.get('Q_flux_mJy', 0.0)
    mfs_U    = mfs.get('U_flux_mJy', 0.0)
    mfs_V    = mfs.get('V_flux_mJy', 0.0)
    mfs_freq = mfs.get('freq_GHz',   float(np.nanmean(freq_arr)))
    mfs_rms  = mfs.get('I_rms_mJy',  1.0)
    mfs_snr  = mfs_I / mfs_rms if mfs_rms > 0 else 0.0

    # ------------------------------------------------------------------
    # Spectral index fit on Stokes I (log-log, S/N-weighted)
    # ------------------------------------------------------------------
    alpha = alpha_err = fit_nu = fit_S = None
    valid = (I_flux > 0) & (I_rms > 0) & (I_flux / I_rms >= 3.0)

    if mfs_snr > 30.0 and np.sum(valid) >= 3:
        nu_v   = freq_arr[valid]
        S_v    = I_flux[valid]
        w_v    = S_v / I_rms[valid]
        log_nu = np.log(nu_v / mfs_freq)
        log_S  = np.log(S_v)
        try:
            coeffs, cov = np.polyfit(log_nu, log_S, 1, w=w_v, cov=True)
            alpha     = coeffs[0]
            alpha_err = np.sqrt(cov[0, 0])
            nu_plot   = np.linspace(nu_v.min(), nu_v.max(), 200)
            fit_S     = np.exp(coeffs[1]) * (nu_plot / mfs_freq) ** alpha
            fit_nu    = nu_plot
            msg(f'  Spectral index ({cal_name} {scan_key}): '
                f'alpha = {alpha:.4f} +/- {alpha_err:.4f}  '
                f'(MFS S/N = {mfs_snr:.1f}, N_chan = {int(np.sum(valid))})')
        except Exception as e:
            msg(f'  WARNING: Spectral index fit failed for {cal_name} {scan_key}: {e}')
    else:
        reason = (f'MFS S/N = {mfs_snr:.1f} (< 30)' if mfs_snr <= 30.0
                  else f'only {int(np.sum(valid))} valid channels (need >= 3)')
        msg(f'  Spectral index fit skipped ({cal_name} {scan_key}): {reason}')

    # ------------------------------------------------------------------
    # Figure: 4 panels stacked vertically, shared x-axis
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    # --- Panel 0: Stokes I ---
    ax = axes[0]
    ax.errorbar(freq_arr, I_flux, yerr=I_rms,
                fmt='o', markersize=4, capsize=3, color='royalblue',
                label='Stokes I', zorder=3)
    ax.axhline(y=mfs_I, color='k', linestyle='--', linewidth=1.5,
               label=f'MFS: {mfs_I:.4f} mJy  (S/N = {mfs_snr:.1f})', zorder=2)
    if fit_nu is not None:
        ax.plot(fit_nu, fit_S, color='tomato', linewidth=1.8,
                label=rf'$\alpha$ = {alpha:.2f} $\pm$ {alpha_err:.2f}', zorder=4)
    _mad_ylim(ax, I_flux, I_rms, floor_lo=0.0)
    ax.set_ylabel('I  (mJy/beam)', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.legend(fontsize=9, loc='best')

    # --- Panels 1-3: Stokes Q, U, V (data-driven range, no forced symmetry) ---
    for ax, flux, rms, label, mfs_val, colour in zip(
            axes[1:],
            [Q_flux,       U_flux,        V_flux],
            [Q_rms,        U_rms,         V_rms],
            ['Q',          'U',           'V'],
            [mfs_Q,        mfs_U,         mfs_V],
            ['darkorange', 'forestgreen', 'mediumpurple']):

        ax.errorbar(freq_arr, flux, yerr=rms,
                    fmt='o', markersize=4, capsize=3, color=colour,
                    label=f'Stokes {label}', zorder=3)
        ax.axhline(y=mfs_val, color='k', linestyle='--', linewidth=1.5,
                   label=f'MFS: {mfs_val:.4f} mJy', zorder=2)
        ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8, zorder=1)
        _mad_ylim(ax, flux, rms)   # no floor -- Q/U/V can be negative
        ax.set_ylabel(f'{label}  (mJy/beam)', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle=':')
        ax.legend(fontsize=9, loc='best')

    axes[-1].set_xlabel('Frequency (GHz)', fontsize=11)
    fig.suptitle(f'{cal_name}  |  {scan_key}  |  IQUV spectrum', fontsize=12, y=1.01)
    plt.tight_layout()

    plot_fname = os.path.join(plot_dir, f'{cal_name}_{scan_key}_IQUV_spectrum.png')
    plt.savefig(plot_fname, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  IQUV spectrum plot saved: {plot_fname}')


def plot_cal_lightcurve(output_dictionary, cal_name):
    """
    Plot per-scan MFS Stokes I, Q, U, V (and optional derived quantities) for a
    calibrator across all available scans.

    Layout (rows shown only when data exist):
      Row 0 : Stokes I  (mJy/beam; floor at zero)
      Row 1 : Stokes Q  (mJy/beam; data-driven, no forced symmetry)
      Row 2 : Stokes U  (mJy/beam; data-driven)
      Row 3 : Stokes V  (mJy/beam; data-driven)
      Row 4 : LP fraction (%) -- only if 'LP_frac' key present in MFS scan dict
      Row 5 : V/I systematic  -- only if 'VFRAC' key present in MFS scan dict

    All y-axes use MAD-based 5-sigma outlier rejection.  Symmetry around zero
    is NOT enforced (the data range drives the limits in all panels except I).

    Parameters
    ----------
    output_dictionary : dict  -- polarization dictionary; MFS sub-dict keyed by scan_key
    cal_name          : str   -- calibrator name for plot title and filename
    """
    plot_dir = os.path.join('RESULTS', 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)

    mfs_dict  = output_dictionary.get('MFS', {})
    scan_keys = sorted(mfs_dict.keys())

    if not scan_keys:
        msg(f'  plot_cal_lightcurve: no MFS scan data for {cal_name} -- skipping')
        return

    scan_idx = np.arange(len(scan_keys))

    def _get(key, default=np.nan):
        return np.array([mfs_dict[sk].get(key, default) for sk in scan_keys], dtype=float)

    I_flux   = _get('I_flux_mJy');  I_rms  = _get('I_rms_mJy')
    Q_flux   = _get('Q_flux_mJy');  Q_rms  = _get('Q_rms_mJy')
    U_flux   = _get('U_flux_mJy');  U_rms  = _get('U_rms_mJy')
    V_flux   = _get('V_flux_mJy');  V_rms  = _get('V_rms_mJy')

    has_quv   = not np.all(np.isnan(Q_flux))
    has_lp    = 'LP_frac' in mfs_dict[scan_keys[0]]
    has_vfrac = 'VFRAC'   in mfs_dict[scan_keys[0]]

    n_rows = 4 + int(has_lp) + int(has_vfrac)

    fig, axes = plt.subplots(n_rows, 1, figsize=(10, 2.8 * n_rows),
                             sharex=True, squeeze=False)
    axes = [ax[0] for ax in axes]   # flatten from (n,1) to list

    ax_idx = 0

    # Row 0: Stokes I
    ax = axes[ax_idx]; ax_idx += 1
    ax.errorbar(scan_idx, I_flux, yerr=I_rms,
                fmt='o', markersize=5, capsize=3, color='royalblue')
    ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8)
    _mad_ylim(ax, I_flux, I_rms, floor_lo=0.0)
    ax.set_ylabel('I  (mJy/beam)', fontsize=10)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.set_title(f'{cal_name}  |  per-scan MFS IQUV', fontsize=12)

    # Rows 1-3: Q, U, V (no floor, data-driven range)
    for flux, rms, label, colour in [
            (Q_flux, Q_rms, 'Q', 'darkorange'),
            (U_flux, U_rms, 'U', 'forestgreen'),
            (V_flux, V_rms, 'V', 'mediumpurple')]:
        ax = axes[ax_idx]; ax_idx += 1
        ax.errorbar(scan_idx, flux, yerr=rms,
                    fmt='o', markersize=5, capsize=3, color=colour)
        ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8)
        _mad_ylim(ax, flux, rms)
        ax.set_ylabel(f'{label}  (mJy/beam)', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle=':')

    # Row 4 (optional): LP fraction
    if has_lp:
        LP_frac = _get('LP_frac');  LP_err = _get('LP_frac_err')
        ax = axes[ax_idx]; ax_idx += 1
        ax.errorbar(scan_idx, LP_frac, yerr=LP_err,
                    fmt='s', markersize=5, capsize=3, color='teal')
        ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8)
        _mad_ylim(ax, LP_frac, LP_err, floor_lo=0.0)
        ax.set_ylabel('LP fraction (%)', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle=':')

    # Row 5 (optional): V/I systematic
    if has_vfrac:
        vfrac = _get('VFRAC')
        ax = axes[ax_idx]; ax_idx += 1
        ax.plot(scan_idx, vfrac, 'D', markersize=5, color='crimson', label='V/I')
        ax.axhline(y=0, color='grey', linestyle=':', linewidth=0.8)
        _mad_ylim(ax, vfrac, np.ones_like(vfrac) * 1e-6)
        ax.set_ylabel('V/I  (systematic)', fontsize=10)
        ax.grid(True, alpha=0.3, linestyle=':')
        ax.legend(fontsize=9)

    axes[-1].set_xlabel('Scan index', fontsize=11)
    axes[-1].set_xticks(scan_idx)
    axes[-1].set_xticklabels(scan_keys, rotation=30, ha='right', fontsize=8)

    plt.tight_layout()

    plot_fname = os.path.join(plot_dir, f'{cal_name}_MFS_lightcurve.png')
    plt.savefig(plot_fname, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  MFS light curve plot saved: {plot_fname}')


# =============================================================================
# Systematic estimation functions
# =============================================================================

def get_primary_systematic(bpcal_name, bpcal_pos):
    """
    Estimate instrumental polarization leakage systematics from the bandpass
    calibrator (J1939-6342 or J0408-6545).

    The bandpass calibrator is assumed to be intrinsically unpolarised.  Any
    residual Stokes Q, U, or V signal at the source position therefore reflects
    unmodelled leakage from Stokes I.  The systematic fractions are:

        BPCAL_RESIDUAL_QFRAC = Q / I
        BPCAL_RESIDUAL_UFRAC = U / I
        BPCAL_RESIDUAL_VFRAC = V / I

    In the ideal case all three are exactly zero.

    Fitting strategy
    ----------------
    Stokes I is fitted with a free position ('abp') -- the source is bright
    and we want the best astrometric reference for the leakage fits.

    Stokes Q, U, and V positions are ALWAYS frozen to the Stokes I pixel
    position ('xyabp').  The source is intrinsically unpolarised, so any
    signal is leakage from I; fitting the position freely in this
    noise-dominated regime produces spurious offsets that contaminate the
    systematic estimate.  check_position() is NOT called for Q, U, or V.

    When multiple scan images are found, systematics are computed per scan
    and averaged, giving a more robust estimate than a single combined image.
    A per-scan MFS light curve plot is produced after fitting.

    If per-channel Stokes I images exist for a scan (same naming convention
    as the pol-angle calibrator's channel images), they are fitted too --
    Q, U, and V again with position frozen to Stokes I -- and a per-scan
    IQUV spectrum plot is written to RESULTS/fitting_plots, exactly as for
    the pol-angle calibrator.  Scans without channel images are still
    processed for their MFS systematic; only the channel spectrum plot is
    skipped for those.

    Parameters
    ----------
    bpcal_name : str  -- calibrator field name as it appears in image filenames
    bpcal_pos  : str  -- CASA-format sky position, e.g. '19:39:25.0264,-63.42.45.624'

    Returns
    -------
    dict with keys:
        BPCAL_RESIDUAL_QFRAC : float  -- mean Q/I leakage fraction
        BPCAL_RESIDUAL_UFRAC : float  -- mean U/I leakage fraction
        BPCAL_RESIDUAL_VFRAC : float  -- mean V/I leakage fraction
        _mfs_dict            : dict   -- per-scan MFS data (internal; for plotting)
        _chan_dict           : dict   -- per-scan channel data (internal; for plotting)
    """
    msg('')
    msg(f'{"="*70}')
    msg(f'  BP calibrator systematic estimation: {bpcal_name}')
    msg(f'{"="*70}')

    prefixes, _ = find_cal_images(bpcal_name, stokes='I')

    if not prefixes:
        msg(f'WARNING: No Stokes I images found for {bpcal_name} -- '
            f'skipping BP systematic estimation, returning zeros')
        return {'BPCAL_RESIDUAL_QFRAC': 0.0,
                'BPCAL_RESIDUAL_UFRAC': 0.0,
                'BPCAL_RESIDUAL_VFRAC': 0.0,
                '_mfs_dict': {}, '_chan_dict': {}}

    msg(f'Looping over {len(prefixes)} scan prefix(es) for {bpcal_name}')

    Q_fracs, U_fracs, V_fracs = [], [], []
    mfs_dict  = {}
    chan_dict = {}

    for k, prefix in enumerate(prefixes):

        scan_key = f'scan{k:03d}'
        msg('')
        msg(f'  --- Scan {k+1}/{len(prefixes)}: {os.path.basename(prefix)} ---')

        im_I = f'{prefix}-MFS-I-image.fits'
        im_Q = f'{prefix}-MFS-Q-image.fits'
        im_U = f'{prefix}-MFS-U-image.fits'
        im_V = f'{prefix}-MFS-V-image.fits'

        missing = [s for s, im in [('I', im_I), ('Q', im_Q), ('U', im_U), ('V', im_V)]
                   if not os.path.exists(im)]
        if missing:
            msg(f'  WARNING: Missing Stokes {missing} image(s) for this scan -- skipping')
            continue

        # Convert sky position to pixel coordinates in this scan's Stokes I image
        region       = f'circle[[{bpcal_pos}],1.0pix]'
        pixel_imstat = imstat(im_I, region=region)
        x_seed       = pixel_imstat['maxpos'][0]
        y_seed       = pixel_imstat['maxpos'][1]

        # ------------------------------------------------------------------
        # Stokes I -- free position fit ('abp')
        # ------------------------------------------------------------------
        make_estimate('estimate_I_bpcal.txt', im_I, x_seed, y_seed, 'abp')
        I_imfit = get_imfit_values('estimate_I_bpcal.txt', im_I, x_seed, y_seed)

        flux_I = I_imfit['results']['component0']['peak']['value'] * 1e3
        err_I  = I_imfit['results']['component0']['peak']['error'] * 1e3
        I_xpix = I_imfit['results']['component0']['pixelcoords'][0]
        I_ypix = I_imfit['results']['component0']['pixelcoords'][1]
        rms_I  = get_imstat_values(im_I, I_xpix, I_ypix)[3] * 1e3

        msg(f'  Stokes I : flux = {flux_I:.4f} mJy  err = {err_I:.4f} mJy  '
            f'rms = {rms_I:.4f} mJy  S/N = {flux_I / rms_I:.1f}')
        msg(f'  Position : x = {I_xpix:.1f}, y = {I_ypix:.1f} pix  '
            f'[astrometric reference for Q, U, V]')

        freq_GHz = imhead(im_I, mode='get', hdkey='CRVAL3')['value'] / 1.0e9
        bmaj     = imhead(im_I, mode='get', hdkey='bmaj')['value']

        scan_mfs = {
            'I_flux_mJy': flux_I, 'I_rms_mJy': rms_I, 'I_err_mJy': err_I,
            'freq_GHz': freq_GHz, 'bmaj_asec': bmaj,
        }

        # ------------------------------------------------------------------
        # Stokes Q, U, V -- position ALWAYS frozen to Stokes I ('xyabp')
        #
        # J1939/J0408 is intrinsically unpolarised; any Q/U/V is pure leakage
        # from I.  Free-fitting the position in this noise-dominated regime
        # introduces spurious offsets that bias the systematic fractions.
        # check_position() is deliberately bypassed here.
        # ------------------------------------------------------------------
        for stokes, im in [('Q', im_Q), ('U', im_U), ('V', im_V)]:

            est_fname = f'estimate_{stokes}_bpcal.txt'
            make_estimate(est_fname, im, I_xpix, I_ypix, 'xyabp')
            pol_imfit = get_imfit_values(est_fname, im, I_xpix, I_ypix)

            flux = pol_imfit['results']['component0']['peak']['value'] * 1e3
            err  = pol_imfit['results']['component0']['peak']['error'] * 1e3
            rms  = get_imstat_values(im, I_xpix, I_ypix)[3] * 1e3
            frac = flux / flux_I

            msg(f'  Stokes {stokes} : flux = {flux:.4f} mJy  err = {err:.4f} mJy  '
                f'rms = {rms:.4f} mJy  S/N = {abs(flux) / rms:.1f}  '
                f'{stokes}/I = {frac:.6f}  [position FIXED to Stokes I]')

            scan_mfs[f'{stokes}_flux_mJy'] = flux
            scan_mfs[f'{stokes}_rms_mJy']  = rms
            scan_mfs[f'{stokes}_err_mJy']  = err

            if stokes == 'Q':
                Q_fracs.append(frac)
            elif stokes == 'U':
                U_fracs.append(frac)
            elif stokes == 'V':
                V_fracs.append(frac)

        mfs_dict[scan_key] = scan_mfs

        # ------------------------------------------------------------------
        # Per-channel Stokes I images (optional) -- fit and plot spectra
        #
        # Q, U, V positions are ALWAYS frozen to the Stokes I pixel position,
        # for the same reason as the MFS fit above: J1939/J0408 is
        # intrinsically unpolarised, so check_position() is bypassed for them.
        # ------------------------------------------------------------------
        CHAN_I_images = sorted(glob.glob(f'{prefix}-[!M]*-I-image.fits'))
        if not CHAN_I_images:
            prefix_zoom   = prefix.replace('diagnostic', 'diagnostic_zoom')
            CHAN_I_images = sorted(glob.glob(f'{prefix_zoom}-[!M]*-I-image.fits'))
            if CHAN_I_images:
                msg(f'  Channel I images not found at prefix; '
                    f'using diagnostic_zoom ({len(CHAN_I_images)} found)')

        scan_chan = {
            'freq_GHz'  : [], 'I_flux_mJy': [], 'Q_flux_mJy': [], 'U_flux_mJy': [],
            'V_flux_mJy': [], 'I_rms_mJy' : [], 'Q_rms_mJy' : [], 'U_rms_mJy' : [],
            'V_rms_mJy' : [],
        }

        if not CHAN_I_images:
            msg(f'  No per-channel Stokes I images found for this scan -- '
                f'skipping channel fits')
        else:
            msg(f'  Found {len(CHAN_I_images)} per-channel I image(s) -- looping')

            for n, ch_I in enumerate(CHAN_I_images):

                ch_Q = ch_I.replace('-I-image.fits', '-Q-image.fits')
                ch_U = ch_I.replace('-I-image.fits', '-U-image.fits')
                ch_V = ch_I.replace('-I-image.fits', '-V-image.fits')

                missing_chan = [s for s, im in [('Q', ch_Q), ('U', ch_U), ('V', ch_V)]
                                if not os.path.exists(im)]
                if missing_chan:
                    msg(f'    Channel {n+1}: missing Stokes {missing_chan} -- skipping')
                    continue

                try:
                    c_freq_GHz = imhead(ch_I, mode='get', hdkey='CRVAL3')['value'] / 1.0e9
                    c_bmaj     = imhead(ch_I, mode='get', hdkey='bmaj')['value']

                    # Skip heavily flagged channels (beam >> MFS beam scaled by freq ratio)
                    beam_scaled = freq_GHz / c_freq_GHz * bmaj
                    if c_bmaj > 3.0 * beam_scaled:
                        msg(f'    Channel {n+1}: BMAJ = {c_bmaj:.2f}" >> 3 x scaled MFS '
                            f'beam = {3.0 * beam_scaled:.2f}" -- likely heavily flagged, '
                            f'skipping')
                        continue

                    # Stokes I -- check_position
                    check_position('estimate_I_bpcal_chan.txt', ch_I, I_xpix, I_ypix)
                    cI_imfit = get_imfit_values('estimate_I_bpcal_chan.txt', ch_I,
                                                 I_xpix, I_ypix)
                    c_flux_I = cI_imfit['results']['component0']['peak']['value'] * 1e3
                    cI_xpix  = cI_imfit['results']['component0']['pixelcoords'][0]
                    cI_ypix  = cI_imfit['results']['component0']['pixelcoords'][1]
                    c_rms_I  = get_imstat_values(ch_I, cI_xpix, cI_ypix)[3] * 1e3

                    # Stokes Q, U, V -- position FROZEN to Stokes I ('xyabp')
                    c_flux = {'I': c_flux_I}
                    c_rms  = {'I': c_rms_I}
                    for stokes, im in [('Q', ch_Q), ('U', ch_U), ('V', ch_V)]:
                        est_fname = f'estimate_{stokes}_bpcal_chan.txt'
                        make_estimate(est_fname, im, I_xpix, I_ypix, 'xyabp')
                        c_imfit = get_imfit_values(est_fname, im, I_xpix, I_ypix)
                        c_flux[stokes] = c_imfit['results']['component0']['peak']['value'] * 1e3
                        c_rms[stokes]  = get_imstat_values(im, I_xpix, I_ypix)[3] * 1e3

                    msg(f'    Channel {n+1}/{len(CHAN_I_images)}: freq = {c_freq_GHz:.4f} GHz  '
                        f'I={c_flux["I"]:.3f}  Q={c_flux["Q"]:.3f}  '
                        f'U={c_flux["U"]:.3f}  V={c_flux["V"]:.3f}  [mJy]')

                    scan_chan['freq_GHz'  ].append(c_freq_GHz)
                    scan_chan['I_flux_mJy'].append(c_flux['I'])
                    scan_chan['Q_flux_mJy'].append(c_flux['Q'])
                    scan_chan['U_flux_mJy'].append(c_flux['U'])
                    scan_chan['V_flux_mJy'].append(c_flux['V'])
                    scan_chan['I_rms_mJy' ].append(c_rms['I'])
                    scan_chan['Q_rms_mJy' ].append(c_rms['Q'])
                    scan_chan['U_rms_mJy' ].append(c_rms['U'])
                    scan_chan['V_rms_mJy' ].append(c_rms['V'])

                except Exception as e:
                    msg(f'    Channel {n+1} fitting failed: {e}')

        chan_dict[scan_key] = scan_chan

        # Per-scan IQUV spectrum plot (MFS + channel data collected above)
        plot_cal_spectrum({'MFS': mfs_dict, 'CHAN': chan_dict}, bpcal_name, scan_key)

    # ------------------------------------------------------------------
    # Average systematics across all successfully processed scans
    # ------------------------------------------------------------------
    n_good = len(Q_fracs)

    if n_good == 0:
        msg('ERROR: No valid scans successfully processed -- returning zero systematics')
        return {'BPCAL_RESIDUAL_QFRAC': 0.0,
                'BPCAL_RESIDUAL_UFRAC': 0.0,
                'BPCAL_RESIDUAL_VFRAC': 0.0,
                '_mfs_dict': mfs_dict, '_chan_dict': chan_dict}

    Q_sys = float(np.mean(Q_fracs))
    U_sys = float(np.mean(U_fracs))
    V_sys = float(np.mean(V_fracs))

    msg('')
    msg(f'  BP cal systematic fractions averaged over {n_good} scan(s):')
    msg(f'    Per-scan Q/I : {[f"{v:.6f}" for v in Q_fracs]}')
    msg(f'    Per-scan U/I : {[f"{v:.6f}" for v in U_fracs]}')
    msg(f'    Per-scan V/I : {[f"{v:.6f}" for v in V_fracs]}')
    msg(f'    BPCAL_RESIDUAL_QFRAC = {Q_sys:.6f}')
    msg(f'    BPCAL_RESIDUAL_UFRAC = {U_sys:.6f}')
    msg(f'    BPCAL_RESIDUAL_VFRAC = {V_sys:.6f}')

    # MFS light curve plot across scans
    if mfs_dict:
        plot_cal_lightcurve({'MFS': mfs_dict}, bpcal_name)

    return {'BPCAL_RESIDUAL_QFRAC': Q_sys,
            'BPCAL_RESIDUAL_UFRAC': U_sys,
            'BPCAL_RESIDUAL_VFRAC': V_sys,
            '_mfs_dict': mfs_dict, '_chan_dict': chan_dict}


def get_polcal_polarization(pacal_name, pacal_pos, bpcal_sys):
    """
    Extract full-Stokes MFS and per-channel flux densities from the pol-angle
    calibrator (e.g. J1331+3030), and estimate the Stokes V leakage systematic.

    For each scan image set found:
      - Stokes I is fitted with a free position ('abp') to establish the
        astrometric reference.
      - Stokes Q and U positions are checked against local S/N via
        check_position() and either freed or fixed accordingly.
      - Stokes V is always fitted with its position frozen to Stokes I
        ('xyabp').  V is intrinsically weak; free-fitting the position in this
        noise-dominated regime biases the V/I systematic estimate.

    Ptot images are not used anywhere in this function.

    The leakage systematic returned is:
        PACAL_RESIDUAL_VFRAC = mean(V / I)  averaged across scans

    A per-calibrator JSON (scan indices stripped from the filename), per-scan
    IQUV spectrum plots, and a cross-scan MFS light curve are all produced.

    Parameters
    ----------
    pacal_name : str  -- pol-angle calibrator field name as in image filenames
    pacal_pos  : str  -- CASA-format sky position
    bpcal_sys  : dict -- bpcal systematic fractions from get_primary_systematic(),
                         stored in the output JSON for traceability

    Returns
    -------
    pacal_sys : float  -- PACAL_RESIDUAL_VFRAC (mean V/I across all good scans)
    """
    msg('')
    msg(f'{"="*70}')
    msg(f'  Pol-angle calibrator systematic & flux extraction: {pacal_name}')
    msg(f'{"="*70}')

    prefixes, _ = find_cal_images(pacal_name, stokes='I')

    if not prefixes:
        msg(f'WARNING: No Stokes I images found for {pacal_name} -- '
            f'skipping pol-cal extraction, returning zero systematic')
        return 0.0

    msg(f'Looping over {len(prefixes)} scan prefix(es) for {pacal_name}')

    output_dictionary = {
        'name'                  : pacal_name,
        'MFS'                   : {},
        'CHAN'                  : {},
        'BPCAL_RESIDUAL_QFRAC'  : bpcal_sys.get('BPCAL_RESIDUAL_QFRAC', 0.0),
        'BPCAL_RESIDUAL_UFRAC'  : bpcal_sys.get('BPCAL_RESIDUAL_UFRAC', 0.0),
        'BPCAL_RESIDUAL_VFRAC'  : bpcal_sys.get('BPCAL_RESIDUAL_VFRAC', 0.0),
        'PACAL_RESIDUAL_VFRAC'  : 0.0,
    }

    V_fracs = []

    for k, prefix in enumerate(prefixes):

        scan_key = f'scan{k:03d}'
        msg('')
        msg(f'  --- Scan {k+1}/{len(prefixes)}: {os.path.basename(prefix)} ---')

        im_I = f'{prefix}-MFS-I-image.fits'
        im_Q = f'{prefix}-MFS-Q-image.fits'
        im_U = f'{prefix}-MFS-U-image.fits'
        im_V = f'{prefix}-MFS-V-image.fits'

        missing = [s for s, im in [('I', im_I), ('Q', im_Q), ('U', im_U), ('V', im_V)]
                   if not os.path.exists(im)]
        if missing:
            msg(f'  WARNING: Missing Stokes {missing} image(s) for this scan -- skipping')
            continue

        # Observation metadata
        freq_GHz = imhead(im_I, mode='get', hdkey='CRVAL3')['value'] / 1.0e9
        date_obs = (imhead(im_I, mode='get', hdkey='DATE-OBS')
                    .replace('/', '-', 2).replace('/', 'T'))
        bmaj     = imhead(im_I, mode='get', hdkey='bmaj')['value']
        bmin     = imhead(im_I, mode='get', hdkey='bmin')['value']
        bpa      = imhead(im_I, mode='get', hdkey='bpa')['value']

        msg(f'  Metadata: date={date_obs}  freq={freq_GHz:.4f} GHz  '
            f'beam={bmaj:.2f}"x{bmin:.2f}" PA={bpa:.1f} deg')

        region       = f'circle[[{pacal_pos}],1.0pix]'
        pixel_imstat = imstat(im_I, region=region)
        x_seed       = pixel_imstat['maxpos'][0]
        y_seed       = pixel_imstat['maxpos'][1]

        # ------------------------------------------------------------------
        # MFS Stokes I -- free position fit
        # ------------------------------------------------------------------
        make_estimate('estimate_I_pacal.txt', im_I, x_seed, y_seed, 'abp')
        I_imfit = get_imfit_values('estimate_I_pacal.txt', im_I, x_seed, y_seed)

        flux_I = I_imfit['results']['component0']['peak']['value'] * 1e3
        err_I  = I_imfit['results']['component0']['peak']['error'] * 1e3
        I_xpix = I_imfit['results']['component0']['pixelcoords'][0]
        I_ypix = I_imfit['results']['component0']['pixelcoords'][1]
        RA_I   = (I_imfit['results']['component0']['shape']['direction']['m0']['value']
                  * 180.0 / np.pi)
        DEC_I  = (I_imfit['results']['component0']['shape']['direction']['m1']['value']
                  * 180.0 / np.pi)
        rms_I  = get_imstat_values(im_I, I_xpix, I_ypix)[3] * 1e3

        msg(f'  Stokes I : flux = {flux_I:.4f} mJy  err = {err_I:.4f} mJy  '
            f'rms = {rms_I:.4f} mJy  S/N = {flux_I / rms_I:.1f}')
        msg(f'  Position : RA = {RA_I:.6f} deg  Dec = {DEC_I:.6f} deg  '
            f'(x={I_xpix:.1f}, y={I_ypix:.1f} pix)  [astrometric reference]')

        # ------------------------------------------------------------------
        # MFS Stokes Q -- check_position (S/N-dependent free vs fixed)
        # ------------------------------------------------------------------
        check_position('estimate_Q_pacal.txt', im_Q, I_xpix, I_ypix)
        Q_imfit = get_imfit_values('estimate_Q_pacal.txt', im_Q, I_xpix, I_ypix)
        flux_Q  = Q_imfit['results']['component0']['peak']['value'] * 1e3
        err_Q   = Q_imfit['results']['component0']['peak']['error'] * 1e3
        Q_xpix  = Q_imfit['results']['component0']['pixelcoords'][0]
        Q_ypix  = Q_imfit['results']['component0']['pixelcoords'][1]
        RA_Q    = (Q_imfit['results']['component0']['shape']['direction']['m0']['value']
                   * 180.0 / np.pi)
        DEC_Q   = (Q_imfit['results']['component0']['shape']['direction']['m1']['value']
                   * 180.0 / np.pi)
        rms_Q   = get_imstat_values(im_Q, Q_xpix, Q_ypix)[3] * 1e3

        msg(f'  Stokes Q : flux = {flux_Q:.4f} mJy  err = {err_Q:.4f} mJy  '
            f'rms = {rms_Q:.4f} mJy  S/N = {abs(flux_Q) / rms_Q:.1f}')

        # ------------------------------------------------------------------
        # MFS Stokes U -- check_position (S/N-dependent free vs fixed)
        # ------------------------------------------------------------------
        check_position('estimate_U_pacal.txt', im_U, I_xpix, I_ypix)
        U_imfit = get_imfit_values('estimate_U_pacal.txt', im_U, I_xpix, I_ypix)
        flux_U  = U_imfit['results']['component0']['peak']['value'] * 1e3
        err_U   = U_imfit['results']['component0']['peak']['error'] * 1e3
        U_xpix  = U_imfit['results']['component0']['pixelcoords'][0]
        U_ypix  = U_imfit['results']['component0']['pixelcoords'][1]
        rms_U   = get_imstat_values(im_U, U_xpix, U_ypix)[3] * 1e3

        msg(f'  Stokes U : flux = {flux_U:.4f} mJy  err = {err_U:.4f} mJy  '
            f'rms = {rms_U:.4f} mJy  S/N = {abs(flux_U) / rms_U:.1f}')

        # ------------------------------------------------------------------
        # MFS Stokes V -- position FROZEN to Stokes I ('xyabp')
        # ------------------------------------------------------------------
        make_estimate('estimate_V_pacal.txt', im_V, I_xpix, I_ypix, 'xyabp')
        V_imfit = get_imfit_values('estimate_V_pacal.txt', im_V, I_xpix, I_ypix)
        flux_V  = V_imfit['results']['component0']['peak']['value'] * 1e3
        err_V   = V_imfit['results']['component0']['peak']['error'] * 1e3
        rms_V   = get_imstat_values(im_V, I_xpix, I_ypix)[3] * 1e3

        msg(f'  Stokes V : flux = {flux_V:.4f} mJy  err = {err_V:.4f} mJy  '
            f'rms = {rms_V:.4f} mJy  S/N = {abs(flux_V) / rms_V:.1f}  '
            f'[position FIXED to Stokes I]')

        # Derived polarization quantities
        flux_P         = np.sqrt(flux_Q**2 + flux_U**2)
        flux_P0, rms_P = calculate_P0(flux_P, rms_Q, rms_U)
        err_P          = ((np.sqrt((flux_Q * err_Q)**2 + (flux_U * err_U)**2) / flux_P)
                          if flux_P > 0 else 0.0)

        LP_frac        = flux_P0 / flux_I * 100.0
        LP_frac_err    = LP_frac * np.sqrt((rms_I  / flux_I)**2 +
                                            (rms_P  / max(flux_P0, 1e-10))**2)
        LP_EVPA        = 0.5 * np.arctan2(flux_U, flux_Q) * 180.0 / np.pi
        LP_EVPA_err    = (0.5 * np.sqrt(flux_U**2 * rms_Q**2 + flux_Q**2 * rms_U**2) /
                          max(flux_U**2 + flux_Q**2, 1e-30) * 180.0 / np.pi)

        V_frac = flux_V / flux_I
        V_fracs.append(V_frac)

        msg(f'  V/I this scan = {V_frac:.6f}  '
            f'LP frac = {LP_frac:.4f}%  EVPA = {LP_EVPA:.4f} deg')

        # ------------------------------------------------------------------
        # Store MFS results for this scan
        # ------------------------------------------------------------------
        output_dictionary['MFS'][scan_key] = {
            'date_obs'    : date_obs,
            'freq_GHz'    : freq_GHz,
            'bmaj_asec'   : bmaj,
            'bmin_asec'   : bmin,
            'bpa_deg'     : bpa,
            'I_flux_mJy'  : flux_I,
            'P_flux_mJy'  : flux_P,
            'P0_flux_mJy' : flux_P0,
            'Q_flux_mJy'  : flux_Q,
            'U_flux_mJy'  : flux_U,
            'V_flux_mJy'  : flux_V,
            'I_rms_mJy'   : rms_I,
            'P_rms_mJy'   : rms_P,
            'Q_rms_mJy'   : rms_Q,
            'U_rms_mJy'   : rms_U,
            'V_rms_mJy'   : rms_V,
            'I_err_mJy'   : err_I,
            'P_err_mJy'   : err_P,
            'Q_err_mJy'   : err_Q,
            'U_err_mJy'   : err_U,
            'V_err_mJy'   : err_V,
            'I_RA_deg'    : RA_I,
            'I_DEC_deg'   : DEC_I,
            'Q_RA_deg'    : RA_Q,
            'Q_DEC_deg'   : DEC_Q,
            'LP_frac'     : LP_frac,
            'LP_frac_err' : LP_frac_err,
            'LP_EVPA'     : LP_EVPA,
            'LP_EVPA_err' : LP_EVPA_err,
            'VFRAC'       : V_frac,
        }

        # ------------------------------------------------------------------
        # Per-channel fits for this scan
        # ------------------------------------------------------------------
        CHAN_I_images = sorted(glob.glob(f'{prefix}-[!M]*-I-image.fits'))

        if not CHAN_I_images:
            prefix_zoom   = prefix.replace('diagnostic', 'diagnostic_zoom')
            CHAN_I_images = sorted(glob.glob(f'{prefix_zoom}-[!M]*-I-image.fits'))
            if CHAN_I_images:
                msg(f'  Channel I images not found at prefix; '
                    f'using diagnostic_zoom ({len(CHAN_I_images)} found)')
                prefix = prefix_zoom
            else:
                msg(f'  WARNING: No per-channel images found for scan {k+1} -- '
                    f'skipping channel fits')
                output_dictionary['CHAN'][scan_key] = {}
                # Plot MFS spectrum (no channels) then continue
                plot_cal_spectrum(output_dictionary, pacal_name, scan_key)
                continue

        msg(f'  Found {len(CHAN_I_images)} per-channel I image(s) -- looping')

        output_dictionary['CHAN'][scan_key] = {
            'freq_GHz'    : [], 'bmaj_asec'  : [], 'bmin_asec'   : [], 'bpa_deg'     : [],
            'I_flux_mJy'  : [], 'Q_flux_mJy' : [], 'U_flux_mJy'  : [], 'V_flux_mJy'  : [],
            'I_rms_mJy'   : [], 'Q_rms_mJy'  : [], 'U_rms_mJy'   : [], 'V_rms_mJy'   : [],
            'I_err_mJy'   : [], 'Q_err_mJy'  : [], 'U_err_mJy'   : [], 'V_err_mJy'   : [],
            'LP_frac'     : [], 'LP_frac_err': [], 'LP_EVPA'      : [], 'LP_EVPA_err' : [],
        }

        for n, chan_I in enumerate(CHAN_I_images):

            chan_Q = chan_I.replace('-I-image.fits', '-Q-image.fits')
            chan_U = chan_I.replace('-I-image.fits', '-U-image.fits')
            chan_V = chan_I.replace('-I-image.fits', '-V-image.fits')

            missing_chan = [s for s, im in [('Q', chan_Q), ('U', chan_U), ('V', chan_V)]
                            if not os.path.exists(im)]
            if missing_chan:
                msg(f'    Channel {n+1}: missing Stokes {missing_chan} -- skipping')
                continue

            try:
                msg('')
                msg(f'    Channel {n+1}/{len(CHAN_I_images)}: {os.path.basename(chan_I)}')

                c_freq_GHz = imhead(chan_I, mode='get', hdkey='CRVAL3')['value'] / 1.0e9
                c_bmaj     = imhead(chan_I, mode='get', hdkey='bmaj')['value']
                c_bmin     = imhead(chan_I, mode='get', hdkey='bmin')['value']
                c_bpa      = imhead(chan_I, mode='get', hdkey='bpa')['value']

                # Skip heavily flagged channels (beam >> MFS beam scaled by freq ratio)
                beam_scaled = (output_dictionary['MFS'][scan_key]['freq_GHz'] /
                               c_freq_GHz * output_dictionary['MFS'][scan_key]['bmaj_asec'])
                if c_bmaj > 3.0 * beam_scaled:
                    msg(f'    Skipping: BMAJ = {c_bmaj:.2f}" >> 3 x scaled MFS beam '
                        f'= {3.0 * beam_scaled:.2f}" -- likely heavily flagged channel')
                    continue

                # Stokes I -- check_position
                check_position('estimate_I_chan.txt', chan_I, I_xpix, I_ypix)
                cI_imfit  = get_imfit_values('estimate_I_chan.txt', chan_I, I_xpix, I_ypix)
                c_flux_I  = cI_imfit['results']['component0']['peak']['value'] * 1e3
                c_err_I   = cI_imfit['results']['component0']['peak']['error'] * 1e3
                cI_xpix   = cI_imfit['results']['component0']['pixelcoords'][0]
                cI_ypix   = cI_imfit['results']['component0']['pixelcoords'][1]
                c_rms_I   = get_imstat_values(chan_I, cI_xpix, cI_ypix)[3] * 1e3

                # Stokes Q -- check_position
                check_position('estimate_Q_chan.txt', chan_Q, I_xpix, I_ypix)
                cQ_imfit  = get_imfit_values('estimate_Q_chan.txt', chan_Q, I_xpix, I_ypix)
                c_flux_Q  = cQ_imfit['results']['component0']['peak']['value'] * 1e3
                c_err_Q   = cQ_imfit['results']['component0']['peak']['error'] * 1e3
                cQ_xpix   = cQ_imfit['results']['component0']['pixelcoords'][0]
                cQ_ypix   = cQ_imfit['results']['component0']['pixelcoords'][1]
                c_rms_Q   = get_imstat_values(chan_Q, cQ_xpix, cQ_ypix)[3] * 1e3

                # Stokes U -- check_position
                check_position('estimate_U_chan.txt', chan_U, I_xpix, I_ypix)
                cU_imfit  = get_imfit_values('estimate_U_chan.txt', chan_U, I_xpix, I_ypix)
                c_flux_U  = cU_imfit['results']['component0']['peak']['value'] * 1e3
                c_err_U   = cU_imfit['results']['component0']['peak']['error'] * 1e3
                cU_xpix   = cU_imfit['results']['component0']['pixelcoords'][0]
                cU_ypix   = cU_imfit['results']['component0']['pixelcoords'][1]
                c_rms_U   = get_imstat_values(chan_U, cU_xpix, cU_ypix)[3] * 1e3

                # Stokes V -- FROZEN to Stokes I position ('xyabp')
                make_estimate('estimate_V_chan.txt', chan_V, I_xpix, I_ypix, 'xyabp')
                cV_imfit  = get_imfit_values('estimate_V_chan.txt', chan_V, I_xpix, I_ypix)
                c_flux_V  = cV_imfit['results']['component0']['peak']['value'] * 1e3
                c_err_V   = cV_imfit['results']['component0']['peak']['error'] * 1e3
                c_rms_V   = get_imstat_values(chan_V, I_xpix, I_ypix)[3] * 1e3

                # Derived polarization quantities
                c_flux_P           = np.sqrt(c_flux_Q**2 + c_flux_U**2)
                c_flux_P0, c_rms_P = calculate_P0(c_flux_P, c_rms_Q, c_rms_U)
                c_LP_frac          = c_flux_P0 / c_flux_I * 100.0
                c_LP_frac_err      = c_LP_frac * np.sqrt(
                    (c_rms_I / c_flux_I)**2 +
                    (c_rms_P / max(c_flux_P0, 1e-10))**2)
                c_LP_EVPA          = 0.5 * np.arctan2(c_flux_U, c_flux_Q) * 180.0 / np.pi
                c_LP_EVPA_err      = (0.5 * np.sqrt(c_flux_U**2 * c_rms_Q**2 +
                                                     c_flux_Q**2 * c_rms_U**2) /
                                      max(c_flux_U**2 + c_flux_Q**2, 1e-30) * 180.0 / np.pi)

                msg(f'      I={c_flux_I:.3f}  Q={c_flux_Q:.3f}  '
                    f'U={c_flux_U:.3f}  V={c_flux_V:.3f}  [mJy]  '
                    f'LP={c_LP_frac:.4f}%  EVPA={c_LP_EVPA:.4f} deg')

                cd = output_dictionary['CHAN'][scan_key]
                cd['freq_GHz'   ].append(c_freq_GHz)
                cd['bmaj_asec'  ].append(c_bmaj)
                cd['bmin_asec'  ].append(c_bmin)
                cd['bpa_deg'    ].append(c_bpa)
                cd['I_flux_mJy' ].append(c_flux_I)
                cd['Q_flux_mJy' ].append(c_flux_Q)
                cd['U_flux_mJy' ].append(c_flux_U)
                cd['V_flux_mJy' ].append(c_flux_V)
                cd['I_rms_mJy'  ].append(c_rms_I)
                cd['Q_rms_mJy'  ].append(c_rms_Q)
                cd['U_rms_mJy'  ].append(c_rms_U)
                cd['V_rms_mJy'  ].append(c_rms_V)
                cd['I_err_mJy'  ].append(c_err_I)
                cd['Q_err_mJy'  ].append(c_err_Q)
                cd['U_err_mJy'  ].append(c_err_U)
                cd['V_err_mJy'  ].append(c_err_V)
                cd['LP_frac'    ].append(c_LP_frac)
                cd['LP_frac_err'].append(c_LP_frac_err)
                cd['LP_EVPA'    ].append(c_LP_EVPA)
                cd['LP_EVPA_err'].append(c_LP_EVPA_err)

            except Exception as e:
                msg(f'    Channel {n+1} fitting failed: {e}')

        # Per-scan IQUV spectrum plot (after channel data are collected for this scan)
        plot_cal_spectrum(output_dictionary, pacal_name, scan_key)

    # ------------------------------------------------------------------
    # Average V/I systematic across all successfully processed scans
    # ------------------------------------------------------------------
    if not V_fracs:
        msg('ERROR: No valid scan V/I measurements obtained -- returning zero systematic')
        pacal_sys = 0.0
    else:
        pacal_sys = float(np.mean(V_fracs))
        msg('')
        msg(f'  Pol-angle cal V/I fractions per scan: {[f"{v:.6f}" for v in V_fracs]}')
        msg(f'  PACAL_RESIDUAL_VFRAC (mean across {len(V_fracs)} scan(s)) = {pacal_sys:.6f}')

    output_dictionary['PACAL_RESIDUAL_VFRAC'] = pacal_sys

    # Cross-scan MFS light curve plot
    if output_dictionary['MFS']:
        plot_cal_lightcurve(output_dictionary, pacal_name)

    # ------------------------------------------------------------------
    # RM synthesis input file (all scan channel data concatenated)
    # ------------------------------------------------------------------
    all_freq, all_I, all_Q, all_U = [], [], [], []
    all_rmsI, all_rmsQ, all_rmsU = [], [], []

    for sk in output_dictionary['CHAN']:
        cd = output_dictionary['CHAN'][sk]
        if cd.get('freq_GHz'):
            all_freq += cd['freq_GHz']
            all_I    += cd['I_flux_mJy']
            all_Q    += cd['Q_flux_mJy']
            all_U    += cd['U_flux_mJy']
            all_rmsI += cd['I_rms_mJy']
            all_rmsQ += cd['Q_rms_mJy']
            all_rmsU += cd['U_rms_mJy']

    json_out  = cal_json_path(pacal_name, prefixes[0])
    out_prefix = json_out.replace('_polarization.json', '')
    os.makedirs(os.path.dirname(json_out), exist_ok=True)

    if all_freq:
        rmsynth_arr = np.array([
            np.array(all_freq)  * 1e9,
            np.array(all_I)     / 1e3,
            np.array(all_Q)     / 1e3,
            np.array(all_U)     / 1e3,
            np.array(all_rmsI)  / 1e3,
            np.array(all_rmsQ)  / 1e3,
            np.array(all_rmsU)  / 1e3,
        ])
        rmsynth_path = f'{out_prefix}_rmsynth.txt'
        np.savetxt(rmsynth_path, rmsynth_arr.T)
        msg(f'  RM synthesis input file saved: {rmsynth_path}')
    else:
        msg('  WARNING: No channel data accumulated -- rmsynth file not written')

    # Per-calibrator JSON (scan indices stripped from filename)
    with open(json_out, 'w') as j:
        json.dump(output_dictionary, j, indent=4)
    msg(f'  Polarization JSON saved: {json_out}')

    return pacal_sys


# =============================================================================
# Entry point
# =============================================================================

def main():
    """
    Orchestrate calibrator systematic estimation and propagate results to all
    target polarization JSON files in cfg.RESULTS.

    Workflow
    --------
    1. Load project_info.json to identify the bandpass calibrator.

    2. Attempt Step 1 (bpcal): get_primary_systematic() estimates Q/I, U/I, V/I
       leakage fractions for J1939-6342 (or J0408-6545) from per-scan images.
       Step 1 is SKIPPED (with a warning, not an exit) if:
         - bpcal_name is not in the known-positions dict; or
         - no bpcal images are found in cfg.IMAGES.

    3. Attempt Step 2 (pacal): get_polcal_polarization() measures V/I for the
       pol-angle calibrator and extracts the full MFS + per-channel flux dictionary.
       Step 2 is SKIPPED if:
         - cfg.POLANG_NAME is empty; or
         - no pacal images are found in cfg.IMAGES.

    4. If BOTH steps were skipped there is nothing to do and main() returns
       after logging a clear diagnostic message.

    5. Append the four systematic terms to every target *_polarization.json in
       cfg.RESULTS.  Calibrator JSONs are excluded.  If no target JSONs exist,
       a msg() is printed and the append step is skipped cleanly.

    Terms appended to each target JSON:
        BPCAL_RESIDUAL_QFRAC  -- Q/I leakage from bandpass calibrator (0 if skipped)
        BPCAL_RESIDUAL_UFRAC  -- U/I leakage from bandpass calibrator (0 if skipped)
        BPCAL_RESIDUAL_VFRAC  -- V/I leakage from bandpass calibrator (0 if skipped)
        PACAL_RESIDUAL_VFRAC  -- V/I residual from pol-angle calibrator (0 if skipped)
    """
    msg('')
    msg(f'{"="*70}')
    msg('  RMSYNTH_01B_systematics.py -- calibrator systematic estimation')
    msg(f'{"="*70}')

    with open('project_info.json') as f:
        project_info = json.load(f)

    bpcal_name = project_info['primary_name']

    bpcal_positions = {
        'J1939-6342': '19:39:25.0264,-63.42.45.624',
        'J0408-6545': '04:08:20.3782,-65.45.09.080',
    }

    # PKS B1934-638 is the B1950 name for J1939-6342 (same source).  Some
    # project_info.json files record the calibrator this way even though the
    # image/JSON filenames on disk still use the field's real name -- so we
    # resolve the POSITION via this alias without renaming bpcal_name itself,
    # keeping find_cal_images() and every filename-based lookup below working
    # against whatever name actually appears on disk.
    if '1934' in bpcal_name and bpcal_name not in bpcal_positions:
        msg(f'  Bandpass calibrator "{bpcal_name}" recognised as J1939-6342 '
            f'(B1934-638 / J1939-6342 are the same source) -- using its position')
        bpcal_positions[bpcal_name] = bpcal_positions['J1939-6342']

    systematics = {
        'BPCAL_RESIDUAL_QFRAC': 0.0,
        'BPCAL_RESIDUAL_UFRAC': 0.0,
        'BPCAL_RESIDUAL_VFRAC': 0.0,
        'PACAL_RESIDUAL_VFRAC': 0.0,
    }

    ran_bpcal = False
    ran_pacal = False

    # ------------------------------------------------------------------
    # Step 1: BP calibrator systematics
    # ------------------------------------------------------------------
    msg('')
    msg(f'Step 1: BP calibrator systematics ({bpcal_name})')

    if bpcal_name not in bpcal_positions:
        msg(f'WARNING: Unknown bandpass calibrator "{bpcal_name}" -- '
            f'add its position to bpcal_positions in main().  Skipping Step 1.')
    else:
        bpcal_pos      = bpcal_positions[bpcal_name]
        bpcal_prefixes, _ = find_cal_images(bpcal_name, stokes='I')

        if not bpcal_prefixes:
            msg(f'WARNING: No images found for {bpcal_name} -- skipping Step 1')
        else:
            bpcal_result = get_primary_systematic(bpcal_name, bpcal_pos)
            systematics['BPCAL_RESIDUAL_QFRAC'] = bpcal_result['BPCAL_RESIDUAL_QFRAC']
            systematics['BPCAL_RESIDUAL_UFRAC'] = bpcal_result['BPCAL_RESIDUAL_UFRAC']
            systematics['BPCAL_RESIDUAL_VFRAC'] = bpcal_result['BPCAL_RESIDUAL_VFRAC']
            ran_bpcal = True

    # ------------------------------------------------------------------
    # Step 2: Pol-angle calibrator systematics and flux extraction
    # ------------------------------------------------------------------
    pacal_name = cfg.POLANG_NAME
    pacal_pos  = cfg.POLANG_DIR
    pol_flag   = pacal_name != ''

    if not pol_flag:
        msg('')
        msg('Step 2: No pol-angle calibrator configured '
            '(cfg.POLANG_NAME is empty) -- skipping')
    else:
        msg('')
        msg(f'Step 2: Pol-angle calibrator systematics & extraction ({pacal_name})')

        pacal_prefixes, _ = find_cal_images(pacal_name, stokes='I')

        if not pacal_prefixes:
            msg(f'WARNING: No images found for {pacal_name} -- skipping Step 2')
        else:
            pacal_json_matches = glob.glob(
                f'{cfg.RESULTS}/*{pacal_name}*_polarization.json')

            if pacal_json_matches and not OVERWRITE:
                msg(f'  OVERWRITE=False and JSON exists -- loading: '
                    f'{os.path.basename(pacal_json_matches[0])}')
                with open(pacal_json_matches[0], 'r') as j:
                    pacal_dict = json.load(j)
                systematics['PACAL_RESIDUAL_VFRAC'] = pacal_dict.get(
                    'PACAL_RESIDUAL_VFRAC', 0.0)
                msg(f'  Loaded PACAL_RESIDUAL_VFRAC = '
                    f'{systematics["PACAL_RESIDUAL_VFRAC"]:.6f}')
            else:
                bpcal_sys_for_pacal = {
                    k: systematics[k] for k in
                    ['BPCAL_RESIDUAL_QFRAC', 'BPCAL_RESIDUAL_UFRAC', 'BPCAL_RESIDUAL_VFRAC']
                }
                pacal_sys = get_polcal_polarization(
                    pacal_name, pacal_pos, bpcal_sys_for_pacal)
                systematics['PACAL_RESIDUAL_VFRAC'] = pacal_sys

            ran_pacal = True

    # ------------------------------------------------------------------
    # Bail out if neither step ran
    # ------------------------------------------------------------------
    if not ran_bpcal and not ran_pacal:
        msg('')
        msg('WARNING: Both Step 1 (bpcal) and Step 2 (pacal) were skipped -- '
            'no calibrator images found.  Have you imaged the calibrators '
            'with full Stokes IQUV?  Exiting without modifying any JSONs.')
        return

    # ------------------------------------------------------------------
    # Step 3: Propagate systematics to all target polarization JSON files
    # ------------------------------------------------------------------
    msg('')
    msg('Step 3: Appending systematic terms to target polarization JSON files')

    all_jsons = glob.glob(cfg.RESULTS + '/*_polarization.json')

    # Exclude calibrator JSONs -- they already embed their own systematics
    exclude_names = [bpcal_name] + ([pacal_name] if pol_flag else [])
    target_jsons  = [f for f in all_jsons
                     if not any(ex in os.path.basename(f) for ex in exclude_names)]

    if not target_jsons:
        msg('  No target polarization JSON files found in '
            f'{cfg.RESULTS} -- skipping the append step.  '
            'Run RMSYNTH_01_extract_fluxes.py first to generate target JSONs.')
    else:
        n_updated = 0
        for fpath in target_jsons:
            fname = os.path.basename(fpath)
            with open(fpath, 'r') as j:
                out_dict = json.load(j)

            out_dict['BPCAL_RESIDUAL_QFRAC'] = systematics['BPCAL_RESIDUAL_QFRAC']
            out_dict['BPCAL_RESIDUAL_UFRAC'] = systematics['BPCAL_RESIDUAL_UFRAC']
            out_dict['BPCAL_RESIDUAL_VFRAC'] = systematics['BPCAL_RESIDUAL_VFRAC']
            out_dict['PACAL_RESIDUAL_VFRAC'] = systematics['PACAL_RESIDUAL_VFRAC']

            with open(fpath, 'w') as j:
                json.dump(out_dict, j, indent=4)

            msg(f'  Updated: {fname}')
            n_updated += 1

        msg(f'Step 3 complete: {n_updated} target JSON file(s) updated')

    # ------------------------------------------------------------------
    # Final summary
    # ------------------------------------------------------------------
    msg('')
    msg(f'{"="*70}')
    msg('  Final systematic values:')
    msg(f'    BPCAL_RESIDUAL_QFRAC = {systematics["BPCAL_RESIDUAL_QFRAC"]:.6f}'
        + ('' if ran_bpcal else '  [not measured -- bpcal skipped]'))
    msg(f'    BPCAL_RESIDUAL_UFRAC = {systematics["BPCAL_RESIDUAL_UFRAC"]:.6f}'
        + ('' if ran_bpcal else '  [not measured -- bpcal skipped]'))
    msg(f'    BPCAL_RESIDUAL_VFRAC = {systematics["BPCAL_RESIDUAL_VFRAC"]:.6f}'
        + ('' if ran_bpcal else '  [not measured -- bpcal skipped]'))
    msg(f'    PACAL_RESIDUAL_VFRAC = {systematics["PACAL_RESIDUAL_VFRAC"]:.6f}'
        + ('' if ran_pacal else '  [not measured -- pacal skipped]'))
    msg(f'{"="*70}')
    msg('RMSYNTH_01B_systematics.py complete.')


if __name__ == '__main__':
    main()
