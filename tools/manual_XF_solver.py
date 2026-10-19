# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
Manual cross-hand phase (XF) solver for the polarization angle calibrator.

Run by 1GC_05_casa_refcal.py when XF_MODE = 'manual', or when 'auto' finds a
phase jump in the CASA solution. Writes GAINTABLES/cal_1GC_<ms>.Xf and
diagnostic plots to GAINPLOTS/manualXF/.

CASA's polcal solves the XF phase only modulo 180 deg. For each channel, the
calibrator's Q, U, V are corrected for the XF phase (with the polcal sign kept
and flipped), the parallactic angle and Faraday rotation, and the two EVPAs
are compared with the calibrator's known polarization angle. Flipping the sign
flips the cross-hand signal, so the two EVPAs differ in proportion to it: where
the cross-hand signal is a small part of the polarized flux they cannot be
told apart. The signs
of all channels are chosen together (Viterbi): channels where the two EVPAs
differ by more than the confidence threshold decide by the EVPA, and the rest
follow continuity of the XF phase from them.

FARADAY ROTATION
    The ionospheric RM comes from Spinifex, corrected for its known bias
    (Perley & Butler 2026) and held at or below zero:
        RM_iono = min(RM_spinifex - SPINIFEX_RM_BIAS, 0)

    The RM is searched over RM_intrinsic + RM_iono +/- 2 x Spinifex error,
    with the upper edge capped at RM_intrinsic: first on a coarse grid, then
    per channel.

    3C286 / 3C138      EVPA model from Perley et al. (2026); RM_intrinsic = 0,
                       as the model already includes the source's own rotation.
    Other calibrators  Constant EVPA of XF_TARGET_POLANG; RM_intrinsic =
                       XF_TARGET_RM.

WORKFLOW
    1.  applycal K, B, Gp, G, Df to the calibrator (no XF)
    2.  Average IQUV per scan (cached in stokes_perscan.npz)
    3.  Parallactic angle per scan
    4.  Ionospheric RM per scan (Spinifex)
    5.  polcal Xf, XF_CHANINT channels per solution
    6.  Cross-hand S/N and imaginary fraction per channel; flag channels below
        XF_MIN_SN or above XF_MAX_IMAG_FRAC
    7.  RM and the sign of every channel fitted together (EVPA where it can
        tell the signs apart, continuity elsewhere)
    8.  Local-trend clip: flag channels off a straight line through their
        neighbours by more than XF_TREND_CLIP_NSIGMA x max(MAD, 1/S/N)
    9.  Optional Savitzky-Golay smoothing (XF_USE_SMOOTHING)
    10. Write the XF table; optional cross-scan average (XF_AVG_SCAN)
    11. Diagnostics; optional applycal check on the MS (XF_APPLY_TO_MS)

SETTINGS (oxkat/config.py)
    XF_TARGET_POLANG, XF_TARGET_RM   EVPA and RM for calibrators other than
                                     3C286 / 3C138
    XF_CHANINT                       Channels per polcal solution
    XF_MIN_SN                        Cross-hand S/N floor (0 = off)
    XF_MAX_IMAG_FRAC                 Cross-hand |Im| / |Re| ceiling (0 = off)
    XF_SIGN_SYS_FLOOR, XF_SIGN_NSIGMA
                                     Confidence threshold of the sign choice:
                                     max(floor deg, nsigma x EVPA noise)
    XF_SIGN_JUMP_COST                Cost of a pi jump between neighbouring
                                     channels, in units of (margin/threshold)^2
    XF_TREND_CLIP_NSIGMA, XF_TREND_CLIP_WINDOW
                                     Local-trend clip (XF_TREND_CLIP_NSIGMA = 0
                                     turns it off; window None = auto)
    XF_USE_SMOOTHING, XF_SAVGOL_WINDOW, XF_SAVGOL_POLYORDER
                                     Savitzky-Golay smoothing
    XF_AVG_SCAN                      One flux-weighted solution for all scans

OUTPUTS (GAINPLOTS/manualXF/)
    Plots are numbered 01-11 in the order they are made (see OUTPUT_PLOTS);
    00_README_plots.txt describes each one. Data files:
    stokes_perscan.npz                  Cached per-scan IQUV
    final_corrected_stokes_spectra.npz  IQUV after XF correction: feed frame
                                        (*_corrected_scans) and sky frame
                                        (Q_sky_scans, U_sky_scans), with
                                        chi_deg_scans, rm_scans (median search
                                        RM per scan) and model_evpa_deg
"""

import functools
import os
import shutil
import subprocess

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

print = functools.partial(print, flush=True)

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    targets = user_targets
    pcal_names = user_pcals
    target_cal_map = user_cal_map


# ============================================================================
# SETTINGS
# ============================================================================

# Spinifex ionospheric RM bias in rad/m^2 (Perley & Butler 2026); subtracted
# from every Spinifex value
SPINIFEX_RM_BIAS = -0.3

# Re-apply all tables plus the new XF table to the MS and plot IQUV from
# CORRECTED_DATA as a check on the analytic correction
XF_APPLY_TO_MS = False

C_LIGHT = 2.998e8   # m/s

# Plot style: larger text, inward major and minor ticks on all four sides,
# thicker axes. Applied only inside the plot functions (see styled), so it
# does not change plots made by the calling script.
PLOT_STYLE = {
    'font.size': 14, 'axes.labelsize': 16, 'axes.titlesize': 15,
    'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 11,
    'axes.linewidth': 1.8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.minor.visible': True, 'ytick.minor.visible': True,
    'xtick.major.size': 7, 'ytick.major.size': 7, 'xtick.major.width': 1.6, 'ytick.major.width': 1.6,
    'xtick.minor.size': 4, 'ytick.minor.size': 4, 'xtick.minor.width': 1.2, 'ytick.minor.width': 1.2,
}


# Output plots (GAINPLOTS/manualXF/), numbered in the order they are made.
# The numbers are fixed, so a plot that is not made leaves a gap.
# key: (file name, description written to 00_README_plots.txt)
OUTPUT_PLOTS = {
    'evpa_model':     ('01_calibrator_evpa_model.png',
                       'Model EVPA of the polarization angle calibrator against frequency '
                       '(3C286 / 3C138 only, Perley et al. 2026), with XF_TARGET_POLANG for reference.'),
    'stokes_pre':     ('02_stokes_IQUV_before_xf.png',
                       'Per-scan I, Q, U, V spectra of the calibrator with K, B, Gp, G, Df applied '
                       'and no XF correction (feed frame).'),
    'sign_choice':    ('03_xf_sign_choice_kept_vs_rejected.png',
                       'Raw polcal XF phase; EVPA residual from the model for the kept sign (on the '
                       'y-scale of the rejected sign, should sit near 0; inset: zoom) and for the '
                       'rejected sign; and the sign margin |res flip - res keep| per channel against '
                       'the confidence threshold (filled: decided by the EVPA, open: decided by '
                       'continuity, red x: the joint choice went against the channel\'s own better sign).'),
    'resid_signed':   ('04_evpa_residuals_after_sign_choice.png',
                       'De-rotated EVPA (XF, parallactic angle and RM removed) against the model, its '
                       'residual against frequency, and a histogram of residuals, before flagging.'),
    'resid_flagged':  ('05_evpa_residuals_after_flagging.png',
                       'As 04, for the channels that survive the S/N cut and the local-trend clip.'),
    'rm_histogram':   ('06_rm_search_per_channel_histogram.png',
                       'Per-channel RM from the fine RM search, with the coarse best RM (dashed) '
                       'and the median (dotted), per scan.'),
    'xf_phase':       ('07_xf_phase_raw_to_final.png',
                       'XF phase at each step: raw polcal; after the sign choice (every channel, '
                       'filled where the polcal sign was kept, open where it was flipped); after '
                       'flagging (marking the S/N, imaginary-fraction and local-trend clip flags); '
                       'the written solution; '
                       'cross-hand S/N per channel against XF_MIN_SN.'),
    'resid_opposite': ('08_evpa_residuals_final_plus_pi.png',
                       'As 04, for the written solution shifted by pi (the other polcal solution). '
                       'Residuals differ from 05 in proportion to the cross-hand signal, so the two '
                       'look alike where it is a small part of the polarized flux.'),
    'stokes_feed':    ('09_stokes_IQUV_after_xf_feed_frame.png',
                       'Per-scan I, Q, U, V after the analytic XF correction, feed frame '
                       '(no parallactic angle correction).'),
    'stokes_sky':     ('10_stokes_IQUV_after_xf_sky_frame.png',
                       'As 09, with Q and U rotated to the sky frame by the parallactic angle; '
                       'the scans should agree with each other.'),
    'psi0_v':         ('11_psi0_and_V_vs_lambda2.png',
                       'Final data against lambda^2: intrinsic angle psi0 = 0.5 arctan2(U_sky, Q_sky) '
                       '- RM lambda^2 (sky frame, median search RM of each scan) with the model EVPA '
                       '(dashed, should match), and Stokes V (should be near 0).'),
    'stokes_casa':    ('12_stokes_IQUV_after_xf_casa_applycal.png',
                       'Per-scan I, Q, U, V read back from CORRECTED_DATA after applycal with the '
                       'new XF table and parang=True (sky frame); should match 10. '
                       'Only with XF_APPLY_TO_MS.'),
}
PLOT = {key: name for key, (name, _) in OUTPUT_PLOTS.items()}


def write_plot_readme(output_dir):
    """Write 00_README_plots.txt describing every output plot and data file."""
    path = os.path.join(output_dir, '00_README_plots.txt')
    with open(path, 'w') as f:
        f.write('Manual XF solver outputs (manual_XF_solver_v2.py)\n')
        f.write('Plots are numbered in the order they are made; a missing number is a plot\n'
                'that was not needed for this run.\n\n')
        for name, description in OUTPUT_PLOTS.values():
            exists = '' if os.path.exists(os.path.join(output_dir, name)) else '  [not made in this run]'
            f.write(f'{name}{exists}\n    {description}\n\n')
        f.write('stokes_perscan.npz\n    Cached per-scan averaged IQUV (complex), reused when the '
                'channels and scans match.\n\n')
        f.write('xf_sign_choice.npz\n    XF gains before (raw_gains) and after (signed_gains) the '
                'sign choice, per scan and XF channel, with its inputs: raw_flags, per-channel '
                'and coarse RM, Q/U/V on the XF grid, chi_deg_scans, model_evpa_deg, iono_rm(_err), '
                'rm_intrinsic, the cross-hand S/N (crosshand_snr, NaN where already flagged), '
                'S/N flags (snr_flagged) and noise per XF channel (crosshand_noise), the cross-hand '
                '|Im| / |Re| (imag_frac) and its flags (imag_flagged), '
                'and the sign evidence: '
                'res_keep / res_flip (|EVPA - model| with the polcal sign kept / flipped, at the '
                'scan RM), threshold, confident and usable.\n\n')
        f.write('final_corrected_stokes_spectra.npz\n    I, Q, U, V after the analytic XF correction: '
                'feed frame (*_corrected_scans) and sky frame (Q_sky_scans, U_sky_scans), with '
                'chi_deg_scans, rm_scans (median search RM per scan) and model_evpa_deg.\n')
    print(f"Saved: {path}")


def styled(plot_function):
    """Run plot_function with PLOT_STYLE applied."""
    @functools.wraps(plot_function)
    def wrapper(*args, **kwargs):
        with matplotlib.rc_context(PLOT_STYLE):
            return plot_function(*args, **kwargs)
    return wrapper


# ============================================================================
# DATA LOADING
# ============================================================================

def read_chan_freq_ghz(table):
    """Channel frequencies (GHz) of an MS or caltable."""
    tb.open(table + '/SPECTRAL_WINDOW')
    chan_freq = tb.getcol('CHAN_FREQ')
    tb.close()
    return (chan_freq / 1.0e9).flatten()


def channel_average_ms(vis, field, max_chan=1024):
    """
    Write a channel-averaged copy of `field` (CORRECTED_DATA -> DATA) when vis
    has more than max_chan channels, and return its path. Returns None when no
    averaging is needed.
    """
    nchan = len(read_chan_freq_ghz(vis))
    if nchan <= max_chan:
        return None

    factor = int(np.ceil(nchan / max_chan))
    tmpms = vis.rstrip('/') + '_tmp_chanavg.ms'
    print(f"{nchan} channels: averaging by {factor} into {tmpms}")
    if os.path.isdir(tmpms):
        shutil.rmtree(tmpms)
    mstransform(vis=vis, outputvis=tmpms, field=field,
                chanaverage=True, chanbin=factor,
                datacolumn='corrected', keepflags=True)
    return tmpms


def build_iquv_flag(raw_flag, nchan, nrow):
    """
    Flags of shape (4, nchan, nrow) for I, Q, U, V. I and Q are flagged where
    either parallel hand is; U and V where either cross hand is.
    """
    if raw_flag.ndim == 3 and raw_flag.shape[0] == 4:
        f_iq = raw_flag[0] | raw_flag[3]
        f_uv = raw_flag[1] | raw_flag[2]
        return np.stack([f_iq, f_iq, f_uv, f_uv], axis=0)
    if raw_flag.ndim == 3 and raw_flag.shape[0] == 1:
        return np.broadcast_to(raw_flag, (4, nchan, nrow))
    if raw_flag.shape == (nchan, nrow):
        return np.broadcast_to(raw_flag[None], (4, nchan, nrow))
    raise ValueError(f"Unrecognised FLAG shape: {raw_flag.shape}")


def build_iquv_weights(d, nchan, nrow):
    """
    Weights of shape (4, nchan, nrow) for I, Q, U, V from WEIGHT_SPECTRUM, or
    WEIGHT when the MS has no WEIGHT_SPECTRUM. I and Q take the mean of the
    parallel-hand weights; U and V the mean of the cross-hand weights.
    """
    if 'weight_spectrum' in d:
        W = np.asarray(d['weight_spectrum'])
        if W.shape == (nchan, nrow):
            return np.broadcast_to(W[None], (4, nchan, nrow))
        if W.ndim == 3 and W.shape[0] == 4:
            w_iq = 0.5 * (W[0] + W[3])
            w_uv = 0.5 * (W[1] + W[2])
            return np.stack([w_iq, w_iq, w_uv, w_uv], axis=0)
        if W.ndim == 3 and W.shape[0] == 1:
            return np.broadcast_to(W, (4, nchan, nrow))
        raise ValueError(f"Unrecognised WEIGHT_SPECTRUM shape: {W.shape}")

    W = np.asarray(d['weight'])
    if W.shape == (nrow,):
        W = W[None]
    if W.ndim == 2 and W.shape[-1] == nrow:
        if W.shape[0] == 4:
            w_iq = 0.5 * (W[0] + W[3])
            w_uv = 0.5 * (W[1] + W[2])
            per_stokes = np.stack([w_iq, w_iq, w_uv, w_uv], axis=0)
        elif W.shape[0] == 1:
            per_stokes = np.repeat(W, 4, axis=0)
        else:
            raise ValueError(f"Unrecognised WEIGHT shape: {W.shape}")
        return np.broadcast_to(per_stokes[:, None, :], (4, nchan, nrow))
    raise ValueError(f"Unrecognised WEIGHT shape: {W.shape}")


def compute_weighted_averages(data, weights, flags):
    """Weighted average over rows of (4, nchan, nrow) data, ignoring flagged samples."""
    W_eff = np.where(flags, 0.0, weights)
    num = np.sum(W_eff * data, axis=-1, dtype=np.complex128)
    den = np.sum(W_eff, axis=-1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return num / den


def read_scan_iquv(vis, field, scan, datacol):
    """
    Average one scan's IQUV visibilities over baselines and time.

    Returns (vis_avg, time, frac_unflagged): vis_avg is complex (4, nchan),
    time is the scan's median time (MJD seconds), frac_unflagged the fraction
    of samples with non-zero weight. Returns None if the scan cannot be read.
    """
    ms.open(vis)
    ms.selectinit(reset=True)
    if not ms.msselect({'field': field, 'scan': str(scan)}):
        ms.close()
        print(f"  WARNING: could not select field {field}, scan {scan}")
        return None
    ms.selectpolarization(['I', 'Q', 'U', 'V'])
    d = ms.getdata([datacol, 'flag', 'weight_spectrum', 'weight', 'time'])
    ms.close()

    if datacol not in d or 'flag' not in d:
        print(f"  WARNING: scan {scan} has no {datacol.upper()} or FLAG")
        return None
    data = d[datacol]
    ncorr, nchan, nrow = data.shape
    if ncorr != 4:
        print(f"  WARNING: scan {scan} has {ncorr} correlations (expected 4)")
        return None

    flags = build_iquv_flag(np.asarray(d['flag']), nchan, nrow)
    weights = np.where(flags, 0.0, build_iquv_weights(d, nchan, nrow))
    vis_avg = compute_weighted_averages(data, weights, flags)
    return vis_avg, float(np.median(d['time'])), float(np.mean(weights > 0))


def load_stokes_per_scan(vis, field, scan_numbers, freq_ghz, cache_path):
    """
    Per-scan averaged IQUV, complex (4, n_scans, nchan), and each scan's median
    time. Reuses cache_path when its channel count and scans match; otherwise
    reads the MS and rewrites the cache.
    """
    if os.path.isfile(cache_path):
        try:
            cache = np.load(cache_path)
            if (cache['freq_ghz'].shape[0] == len(freq_ghz) and
                    np.array_equal(cache['scan_numbers'], scan_numbers)):
                print(f"Using cached IQUV: {cache_path}")
                return cache['vis_avg'], cache['scan_times']
            print(f"Cached IQUV does not match this MS; re-reading: {cache_path}")
        except Exception as e:
            print(f"Could not read cached IQUV ({e}); re-reading the MS")

    tb.open(vis)
    datacol = 'corrected_data' if 'CORRECTED_DATA' in tb.colnames() else 'data'
    tb.close()
    print(f"Averaging IQUV per scan from {datacol.upper()}")

    vis_avg = np.full((4, len(scan_numbers), len(freq_ghz)), np.nan, dtype=complex)
    scan_times = np.zeros(len(scan_numbers))
    for i, scan in enumerate(scan_numbers):
        result = read_scan_iquv(vis, field, scan, datacol)
        if result is None:
            continue
        vis_avg[:, i, :], scan_times[i], frac_good = result
        print(f"  Scan {scan}: {100 * frac_good:.1f}% unflagged, "
              f"median time {scan_times[i]:.1f} s (MJD)")

    np.savez(cache_path, freq_ghz=freq_ghz, scan_numbers=scan_numbers,
             scan_times=scan_times, vis_avg=vis_avg)
    print(f"Saved: {cache_path}")
    return vis_avg, scan_times


def report_iquv_quality(vis_avg, scan_numbers):
    """
    Print per-scan statistics of Re(IQUV) and the fraction of flux in the
    imaginary part, which should be well under 10% for a calibrated point
    source at the phase centre.
    """
    for i, scan in enumerate(scan_numbers):
        print(f"  Scan {scan}:")
        for s, avg in zip('IQUV', vis_avg[:, i, :]):
            re = np.real(avg)
            re = re[np.isfinite(re)]
            with np.errstate(invalid='ignore', divide='ignore'):
                im_frac = 100 * np.nanmedian(np.abs(np.imag(avg)) / np.abs(avg))
            if re.size == 0:
                print(f"    {s}: no finite values")
                continue
            med = np.median(re)
            mad = np.median(np.abs(re - med))
            status = 'WARNING: high' if im_frac > 10.0 else 'ok'
            print(f"    {s}: median {med:+.4g} Jy, MAD {mad:.4g} Jy, "
                  f"imaginary {im_frac:.1f}% ({status})")


# ============================================================================
# CALIBRATOR EVPA MODELS
# ============================================================================

def compute_3c286_evpa(freq_ghz):
    """
    Intrinsic EVPA (deg) of 3C286 (J1331+3030), Perley et al. (2026, ApJS 283,
    82), Sec. 8.1, Eqs. 9-10, with x = log10(nu / GHz):

        0.5 <= nu <= 1.0 GHz :  26.0 + 57.0 x + 615 x^2 + 3790 x^3
        1.0 <  nu <= 50  GHz :  26.1 + 17.1 x - 16.1 x^2 + 5.75 x^3

    Channels outside 0.5-50 GHz are NaN.
    """
    freq_ghz = np.asarray(freq_ghz, dtype=float)
    x = np.log10(freq_ghz)
    evpa = np.full_like(freq_ghz, np.nan)

    low = (freq_ghz >= 0.5) & (freq_ghz <= 1.0)
    high = (freq_ghz > 1.0) & (freq_ghz <= 50.0)
    evpa[low] = 26.0 + 57.0 * x[low] + 615.0 * x[low]**2 + 3790.0 * x[low]**3
    evpa[high] = 26.1 + 17.1 * x[high] - 16.1 * x[high]**2 + 5.75 * x[high]**3

    n_out = int(np.sum(~(low | high)))
    if n_out:
        print(f"WARNING: {n_out} channel(s) outside the 0.5-50 GHz range of the "
              f"3C286 model -> NaN")
    return evpa


def compute_3c138_evpa(freq_ghz):
    """
    Intrinsic EVPA (deg) of 3C138 (J0521+1638), Perley et al. (2026, ApJS 283,
    82), Sec. 8.2, Eqs. 11-12, with x = log10(nu / GHz):

        0.5 <= nu <= 1.0 GHz : -21.9 + 71.1 x - 1435 x^2 - 11110 x^3 - 23190 x^4
        1.0 <  nu <= 4.0 GHz : -22.0 + 95.4 x -  340 x^2 +   534 x^3 -   308 x^4

    Channels outside 0.5-4 GHz are NaN: above 4 GHz the EVPA oscillates and
    3C138 is not a reliable angle calibrator.
    """
    freq_ghz = np.asarray(freq_ghz, dtype=float)
    x = np.log10(freq_ghz)
    evpa = np.full_like(freq_ghz, np.nan)

    low = (freq_ghz >= 0.5) & (freq_ghz <= 1.0)
    high = (freq_ghz > 1.0) & (freq_ghz <= 4.0)
    evpa[low] = (-21.9 + 71.1 * x[low] - 1435.0 * x[low]**2
                 - 11110.0 * x[low]**3 - 23190.0 * x[low]**4)
    evpa[high] = (-22.0 + 95.4 * x[high] - 340.0 * x[high]**2
                  + 534.0 * x[high]**3 - 308.0 * x[high]**4)

    n_out = int(np.sum(~(low | high)))
    if n_out:
        print(f"WARNING: {n_out} channel(s) outside the 0.5-4 GHz range of the "
              f"3C138 model -> NaN")
    return evpa


def calibrator_model_name(field_name):
    """'3C286' or '3C138' when field_name is one of the modelled calibrators, else None."""
    name = field_name.lower()
    if '3c286' in name or 'j1331' in name:
        return '3C286'
    if '3c138' in name or 'j0521' in name:
        return '3C138'
    return None


def calibrator_evpa(model_name, freq_ghz):
    """EVPA (deg) per channel: the model for 3C286 / 3C138, else XF_TARGET_POLANG."""
    if model_name == '3C286':
        return compute_3c286_evpa(freq_ghz)
    if model_name == '3C138':
        return compute_3c138_evpa(freq_ghz)
    return np.full(len(freq_ghz), float(XF_TARGET_POLANG))


@styled
def plot_evpa_model(model_name, freq_ghz, evpa, output_dir):
    """Plot the calibrator EVPA model against frequency."""
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    ax.plot(freq_ghz, evpa, 'b-', linewidth=2, label=f'{model_name} EVPA model')
    ax.axhline(XF_TARGET_POLANG, color='red', linestyle='--', linewidth=1.5,
               alpha=0.7, label=f'XF_TARGET_POLANG = {XF_TARGET_POLANG} deg')
    ax.set_xlabel('Frequency [GHz]')
    ax.set_ylabel('EVPA [deg]')
    ax.set_title(f'{model_name} EVPA model (Perley et al. 2026)')
    _legend(ax)
    ax.grid(True, alpha=0.3)
    _save(fig, output_dir, PLOT['evpa_model'])


# ============================================================================
# POLARIZATION GEOMETRY
# ============================================================================

def _quantity_deg(q):
    """CASA quantity in degrees, as a float."""
    return float(np.atleast_1d(qa.convert(q, 'deg')['value'])[0])


def compute_parallactic_angle(vis, field_name, time_mjd_s):
    """
    Parallactic angle of field_name at time_mjd_s (MJD seconds): the position
    angle of the zenith seen from the source, in deg, wrapped to [-180, 180).
    Returns (chi_deg, az_deg, el_deg).
    """
    msmd.open(vis)
    phase_dir = msmd.phasecenter(int(msmd.fieldsforname(field_name)[0]))
    msmd.close()

    tb.open(vis + '/OBSERVATION')
    observatory = str(tb.getcol('TELESCOPE_NAME')[0])
    tb.close()

    me.doframe(me.observatory(observatory))
    me.doframe(me.epoch('utc', qa.quantity(time_mjd_s, 's')))
    zenith = me.direction('AZELGEO', '0deg', '90deg')
    chi_deg = _quantity_deg(me.posangle(phase_dir, zenith))
    chi_deg = ((chi_deg + 180.0) % 360.0) - 180.0

    azel = me.measure(phase_dir, 'azelgeo')
    return chi_deg, _quantity_deg(azel['m0']), _quantity_deg(azel['m1'])


def correct_crosshand_phase(u, v, rho):
    """Rotate U, V by the cross-hand phase rho (rad)."""
    cos_rho, sin_rho = np.cos(rho), np.sin(rho)
    return u * cos_rho - v * sin_rho, u * sin_rho + v * cos_rho


def correct_parallactic_angle(q, u, chi_deg):
    """Rotate Q, U from the feed frame to the sky frame (by -2 chi)."""
    ang = np.radians(-2 * chi_deg)
    return q * np.cos(ang) + u * np.sin(ang), -q * np.sin(ang) + u * np.cos(ang)


def calculate_derotated_angle(q, u, rm, lambda_sq):
    """
    EVPA (deg, wrapped to [-90, 90)) after removing Faraday rotation RM * lambda^2.
    rm may be a scalar or one value per channel.
    """
    ang = 2 * rm * lambda_sq
    q_derot = q * np.cos(ang) + u * np.sin(ang)
    u_derot = -q * np.sin(ang) + u * np.cos(ang)
    evpa = np.degrees(0.5 * np.arctan2(u_derot, q_derot))
    return ((evpa + 90.0) % 180.0) - 90.0


def derotated_evpa(gains, Q, U, V, chi_deg, rm, freq_ghz):
    """
    Sky-frame EVPA (deg) per channel after removing the XF phase of `gains`,
    the parallactic angle and Faraday rotation `rm` (scalar or per channel).
    The XF phase is removed as CASA applies it, with conj(gain).
    """
    U_xh, _ = correct_crosshand_phase(U, V, np.angle(np.conj(gains)))
    Q_sky, U_sky = correct_parallactic_angle(Q, U_xh, chi_deg)
    lambda_sq = (C_LIGHT / (freq_ghz * 1e9)) ** 2
    return calculate_derotated_angle(Q_sky, U_sky, rm, lambda_sq)


def evpa_residual(evpa, target):
    """EVPA minus target, wrapped to [-90, 90) deg."""
    return (evpa - target + 90) % 180 - 90


# ============================================================================
# IONOSPHERE
# ============================================================================

def get_spinifex_rm_per_scan(field_name, scan_numbers):
    """
    Ionospheric RM and its error per scan (rad/m^2) for field_name.

    Runs RMSYNTH_03_run_SPINIFEX.py unless RESULTS/<ms>_spinifex_rm.txt already
    exists. Per scan, RM_spinifex is the median over the scan's time steps,
    corrected for SPINIFEX_RM_BIAS and capped at 0:
        RM_iono = min(RM_spinifex - SPINIFEX_RM_BIAS, 0)
    The error combines the median Spinifex model error with the RM drift
    (standard deviation over the scan's time steps) in quadrature.
    A scan with no Spinifex value gets RM 0 and error 0.
    """
    n_scans = len(scan_numbers)
    iono_rm = np.zeros(n_scans)
    iono_rm_err = np.zeros(n_scans)

    script = os.path.join(OXKAT, 'RMSYNTH_03_run_SPINIFEX.py')
    output_file = os.path.join(RESULTS, f"{os.path.basename(myms)}_spinifex_rm.txt")

    if os.path.isfile(output_file):
        print(f"Using existing Spinifex output: {output_file}")
    else:
        print(f"Running Spinifex: {script}")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        result = subprocess.run(['python-spinifex', script], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"WARNING: Spinifex failed (return code {result.returncode}); "
                  f"ionospheric RM set to 0 for every scan")
            print(f"  stderr: {result.stderr[:500]}")
            return iono_rm, iono_rm_err

    # Columns: field, RA, Dec, scan, time (ISOT), time (MJD), RM, RM error
    rows = {int(scan): [] for scan in scan_numbers}
    try:
        with open(output_file) as f:
            for line in f:
                parts = line.split()
                if len(parts) < 8 or parts[0].startswith('#') or parts[0] != field_name:
                    continue
                scan = int(parts[3])
                if scan in rows:
                    rows[scan].append((float(parts[6]), float(parts[7])))
    except (OSError, ValueError) as e:
        print(f"WARNING: could not read Spinifex output ({e}); "
              f"ionospheric RM set to 0 for every scan")
        return iono_rm, iono_rm_err

    print(f"Ionospheric RM (bias {SPINIFEX_RM_BIAS:+.2f} rad/m^2 removed, capped at 0):")
    for i, scan in enumerate(scan_numbers):
        if not rows[int(scan)]:
            print(f"  Scan {scan}: no Spinifex value for {field_name}; RM set to 0")
            continue
        rm_vals, err_vals = np.array(rows[int(scan)]).T
        rm_raw = np.median(rm_vals)
        model_err = np.median(err_vals)
        drift = np.std(rm_vals)
        rm_corr = rm_raw - SPINIFEX_RM_BIAS
        iono_rm[i] = min(rm_corr, 0.0)
        iono_rm_err[i] = np.hypot(model_err, drift)
        capped = '  (capped at 0)' if rm_corr > 0 else ''
        print(f"  Scan {scan}: Spinifex {rm_raw:+.3f} -> corrected {rm_corr:+.3f} "
              f"-> RM_iono {iono_rm[i]:+.3f} +/- {iono_rm_err[i]:.3f} rad/m^2 "
              f"(model {model_err:.3f}, drift {drift:.3f}){capped}")
    return iono_rm, iono_rm_err


# ============================================================================
# 180 DEG AMBIGUITY
# ============================================================================

def _valid(*arrays):
    """Channels where every array is finite."""
    return np.logical_and.reduce([np.isfinite(a) for a in arrays])


def rm_search_window(iono_rm, iono_rm_err, rm_intrinsic):
    """
    RM search range (lo, hi) in rad/m^2: rm_intrinsic + RM_iono +/- 2 x error,
    with the ionospheric part capped at 0.
    """
    lo = iono_rm - 2 * iono_rm_err
    hi = min(iono_rm + 2 * iono_rm_err, 0.0)
    return rm_intrinsic + lo, rm_intrinsic + hi


def _rm_grid(lo, hi, step):
    """Evenly spaced trial RMs from lo to hi inclusive, spacing at most step."""
    n = int(np.ceil((hi - lo) / step - 1e-9)) + 1
    return np.linspace(lo, hi, max(n, 1))


def sign_evidence(gains, Q, U, V, chi_deg, rm, freq_ghz, target):
    """
    |EVPA - model| (deg) per channel with the polcal sign kept (res_keep) and
    flipped (res_flip), at one RM for the whole scan. The two differ in
    proportion to the cross-hand signal, so they are close where it is a small
    part of the polarized flux. NaN where invalid.
    """
    res_keep = np.full(len(gains), np.nan)
    res_flip = np.full(len(gains), np.nan)
    ok = _valid(gains, Q, U, V, target)
    args = (Q[ok], U[ok], V[ok], chi_deg, rm, freq_ghz[ok])
    res_keep[ok] = np.abs(evpa_residual(derotated_evpa(gains[ok], *args), target[ok]))
    res_flip[ok] = np.abs(evpa_residual(derotated_evpa(-gains[ok], *args), target[ok]))
    return res_keep, res_flip


def evpa_noise_deg(Q, U, noise):
    """
    EVPA noise (deg) per XF channel: 0.5 x noise / P, with P = sqrt(Q^2 + U^2)
    and noise the cross-hand noise of the XF channel.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.degrees(0.5 * noise / np.hypot(Q, U))


def resolve_signs_joint(gains, usable, res_keep, res_flip, threshold, jump_cost):
    """
    Choose the sign of every usable channel of one scan together, as the
    lowest-cost path through the channels in frequency order (Viterbi), with
    two states per channel: polcal sign kept or flipped.

    Data cost: the worse sign costs (|res_flip - res_keep| / threshold)^2, the
    better sign 0. Channels whose margin is below the threshold (the EVPA
    cannot tell the signs apart) therefore cost almost nothing either way.
    Jump cost: a pi jump between neighbouring channels costs jump_cost (less
    when their phase step is near 90 deg and so ambiguous itself, and spread
    over any flagged gap). Confident channels set the sign; ambiguous ones
    follow continuity from them.

    Returns (flip, confident): bool arrays over all channels (False where not
    usable). confident marks channels whose margin exceeds the threshold.
    """
    n = len(gains)
    flip = np.zeros(n, dtype=bool)
    margin = res_flip - res_keep                           # > 0: keeping is better
    confident = usable & (np.abs(margin) > threshold)
    idx = np.where(usable)[0]
    if idx.size == 0:
        return flip, confident

    data = np.zeros((idx.size, 2))
    data[:, 0] = (np.clip(-margin[idx], 0, None) / threshold[idx]) ** 2   # keep
    data[:, 1] = (np.clip(margin[idx], 0, None) / threshold[idx]) ** 2    # flip

    phase = np.angle(gains[idx])
    cost = data[0].copy()
    back = np.zeros((idx.size, 2), dtype=int)
    for k in range(1, idx.size):
        step = phase[k] - phase[k - 1]
        gap = idx[k] - idx[k - 1]
        # Same state: a jump if the step is near +/-pi; different states: if near 0
        same = jump_cost * (1.0 - np.cos(step)) / 2.0 / gap
        diff = jump_cost * (1.0 + np.cos(step)) / 2.0 / gap
        trans = np.array([[same, diff], [diff, same]])     # [from, to]
        total = cost[:, None] + trans
        back[k] = np.argmin(total, axis=0)
        cost = total[back[k], [0, 1]] + data[k]

    state = int(np.argmin(cost))
    for k in range(idx.size - 1, -1, -1):
        flip[idx[k]] = state == 1
        state = back[k, state]
    return flip, confident


def rm_search(gains, usable, Q, U, V, chi_deg, freq_ghz, target,
              iono_rm, iono_rm_err, rm_intrinsic, threshold, jump_cost):
    """
    RM and XF signs of one scan, fitted together over rm_search_window(iono_rm,
    iono_rm_err, rm_intrinsic).

    Coarse grid, step max(0.05, error / 3): at each trial RM the signs of all
    usable channels are chosen jointly (resolve_signs_joint), and the trial
    whose chosen signs give the smallest median |EVPA - model| is the scan's
    RM; its signs are the ones applied. Fine grid, step 0.04 within +/-0.25 of
    it (and inside the window): each channel's own best RM with its sign held
    fixed, used for the EVPA diagnostics.

    Returns (rm_chan, coarse_rm, flip, confident, res_keep, res_flip), the
    last four at coarse_rm; rm_chan is NaN outside `usable`.
    """
    rm_chan = np.full(len(gains), np.nan)
    lo, hi = rm_search_window(iono_rm, iono_rm_err, rm_intrinsic)

    def signs_at(rm):
        res_keep, res_flip = sign_evidence(gains, Q, U, V, chi_deg, rm, freq_ghz, target)
        flip, confident = resolve_signs_joint(gains, usable, res_keep, res_flip, threshold, jump_cost)
        return flip, confident, res_keep, res_flip

    if np.sum(usable) < 5:
        coarse_rm = 0.5 * (lo + hi)
        print(f"    only {int(np.sum(usable))} usable channels; RM set to the window centre")
        rm_chan[usable] = coarse_rm
        return (rm_chan, coarse_rm) + signs_at(coarse_rm)

    coarse = _rm_grid(lo, hi, max(0.05, iono_rm_err / 3.0))
    scores = []
    for rm in coarse:
        flip, _, res_keep, res_flip = signs_at(rm)
        scores.append(np.median(np.where(flip, res_flip, res_keep)[usable]))
    coarse_rm = float(coarse[int(np.argmin(scores))])
    flip, confident, res_keep, res_flip = signs_at(coarse_rm)
    print(f"    Coarse: {len(coarse)} trials over [{lo:+.3f}, {hi:+.3f}] rad/m^2; best "
          f"{coarse_rm:+.3f} (median |EVPA - model| {min(scores):.2f} deg with the joint signs)")

    s = np.where(flip, -gains, gains)[usable]
    args = (Q[usable], U[usable], V[usable], chi_deg)
    f, t = freq_ghz[usable], target[usable]
    fine = _rm_grid(max(lo, coarse_rm - 0.25), min(hi, coarse_rm + 0.25), 0.04)
    best_dev = np.full(s.size, np.inf)
    best_rm = np.full(s.size, coarse_rm)
    for rm in fine:
        dev = np.abs(evpa_residual(derotated_evpa(s, *args, rm, f), t))
        better = dev < best_dev
        best_dev[better], best_rm[better] = dev[better], rm

    rm_chan[usable] = best_rm
    print(f"    Fine: {len(fine)} trials; per-channel RM {best_rm.min():+.3f} to "
          f"{best_rm.max():+.3f}, median {np.median(best_rm):+.3f} rad/m^2; "
          f"median |EVPA - model| {np.median(best_dev):.2f} deg")
    return rm_chan, coarse_rm, flip, confident, res_keep, res_flip


def contiguous_ranges(freq_ghz, mask):
    """Frequency ranges 'lo-hi GHz' of the runs of True in mask (channels in order)."""
    ranges, start = [], None
    for k in range(len(mask) + 1):
        on = k < len(mask) and mask[k]
        if on and start is None:
            start = k
        elif not on and start is not None:
            ranges.append(f'{freq_ghz[start]:.3f}-{freq_ghz[k - 1]:.3f} GHz')
            start = None
    return ranges


def print_sign_table(freq_ghz, raw, flip, usable, res_keep, res_flip, threshold, confident, target):
    """
    Per usable channel: frequency, raw polcal XF phase, sign applied, the
    resulting XF phase, model EVPA, EVPA residual with the polcal sign kept
    and flipped, the confidence threshold, and what decided the sign (EVPA,
    or continuity where the margin is below the threshold).
    """
    print(f"    {'GHz':>8} {'raw XF':>8} {'sign':>5} {'new XF':>8} {'model':>8} {'res keep':>9} "
          f"{'res flip':>9} {'thresh':>7} {'decided by':>11}")
    for c in np.where(usable)[0]:
        decided = 'EVPA' if confident[c] else 'continuity'
        new = -raw[c] if flip[c] else raw[c]
        print(f"    {freq_ghz[c]:>8.4f} {np.degrees(np.angle(raw[c])):>8.2f} "
              f"{'-1' if flip[c] else '+1':>5} {np.degrees(np.angle(new)):>8.2f} "
              f"{target[c]:>8.2f} {res_keep[c]:>9.2f} "
              f"{res_flip[c]:>9.2f} {threshold[c]:>7.2f} {decided:>11}")


# ============================================================================
# FLAGGING
# ============================================================================

def bin_to_grid(values, n_out):
    """
    Mean of the finite values in each of n_out equal blocks of channels, and
    the number of finite values in each block.
    """
    edges = np.round(np.linspace(0, len(values), n_out + 1)).astype(int)
    out = np.full(n_out, np.nan)
    count = np.zeros(n_out, dtype=int)
    for i in range(n_out):
        block = values[edges[i]:edges[i + 1]]
        block = block[np.isfinite(block)]
        count[i] = block.size
        if block.size:
            out[i] = block.mean()
    return out, count


def estimate_crosshand_noise(U_scans, V_scans):
    """
    Noise (Jy) on the native-resolution cross-hand flux sqrt(U^2 + V^2),
    pooled over scans. Uses the MAD of first differences in windows of
    max(4, nchan / 100) channels, so the Faraday oscillation of the spectrum
    is not counted as noise. Returns 0 when it cannot be estimated.
    """
    nchan = U_scans.shape[1]
    win = max(4, nchan // 100)
    mads = []
    for xh in np.sqrt(U_scans**2 + V_scans**2):
        for start in range(0, nchan, win):
            seg = xh[start:start + win]
            seg = seg[np.isfinite(seg)]
            if seg.size >= 4:
                d = np.diff(seg)
                mads.append(np.median(np.abs(d - np.median(d))) / np.sqrt(2.0))
    return 1.4826 * float(np.median(mads)) if mads else 0.0


def snr_flag(U_scans, V_scans, flags, scan_numbers, n_xf, min_sn):
    """
    Cross-hand S/N per XF channel: sqrt(U^2 + V^2) binned onto the XF grid,
    over the noise of that bin (native-channel noise / sqrt(channels in the
    bin), as polcal solves on the same bins); NaN where already flagged. When
    min_sn > 0, flags the channels below it: there the polcal phase is set by
    noise. Updates `flags` in place and returns (snr, newly_flagged,
    noise_xf), dicts over scans; noise_xf is the noise (Jy) per XF channel,
    NaN where it could not be estimated.
    """
    snr = {scan: np.full(n_xf, np.nan) for scan in scan_numbers}
    new = {scan: np.zeros(n_xf, dtype=bool) for scan in scan_numbers}
    noise_xf = {scan: np.full(n_xf, np.nan) for scan in scan_numbers}
    noise = estimate_crosshand_noise(U_scans, V_scans)
    if noise <= 0:
        print("  WARNING: could not estimate the cross-hand noise; no S/N available")
        return snr, new, noise_xf
    print(f"  Cross-hand noise: {noise:.4e} Jy per native channel")
    if min_sn <= 0:
        print("  S/N flagging disabled (XF_MIN_SN <= 0)")

    for i, scan in enumerate(scan_numbers):
        xh, count = bin_to_grid(np.sqrt(U_scans[i]**2 + V_scans[i]**2), n_xf)
        with np.errstate(divide='ignore'):
            noise_xf[scan] = np.where(count > 0, noise / np.sqrt(count), np.nan)
        snr[scan] = np.where(flags[scan], np.nan, xh / noise_xf[scan])
        if min_sn > 0:
            new[scan] = np.isfinite(snr[scan]) & (snr[scan] < min_sn)
            flags[scan] |= new[scan]
        finite = snr[scan][np.isfinite(snr[scan])]
        rng = f"{finite.min():.1f} to {finite.max():.1f}" if finite.size else "none"
        print(f"  Scan {scan}: {int(new[scan].sum())} channels below S/N {min_sn:.1f} "
              f"(S/N {rng})")
    return snr, new, noise_xf


def imag_fraction_flag(vis_avg, flags, scan_numbers, n_xf, max_frac):
    """
    Imaginary fraction of the cross-hand per XF channel, |Im(U, V)| / |Re(U, V)|,
    with U and V binned onto the XF grid (NaN where there is no data). A point
    source at the phase centre has real Stokes visibilities and the XF phase
    only mixes Re(U) and Re(V), so the imaginary part measures what else is in
    the cross-hands (off-axis sources, residual phase errors). Where it is a
    large fraction of the real part, the polcal phase follows that
    contamination rather than the calibrator. When max_frac > 0, flags the
    unflagged channels above it. Updates `flags` in place and returns (frac,
    newly_flagged), dicts over scans.
    """
    frac = {scan: np.full(n_xf, np.nan) for scan in scan_numbers}
    new = {scan: np.zeros(n_xf, dtype=bool) for scan in scan_numbers}
    if max_frac <= 0:
        print("  Imaginary-fraction flagging disabled (XF_MAX_IMAG_FRAC <= 0)")
    for i, scan in enumerate(scan_numbers):
        re_u, im_u, re_v, im_v = (bin_to_grid(part, n_xf)[0] for part in
                                  (vis_avg[2, i].real, vis_avg[2, i].imag,
                                   vis_avg[3, i].real, vis_avg[3, i].imag))
        with np.errstate(divide='ignore', invalid='ignore'):
            frac[scan] = np.hypot(im_u, im_v) / np.hypot(re_u, re_v)
        if max_frac > 0:
            new[scan] = ~flags[scan] & np.isfinite(frac[scan]) & (frac[scan] > max_frac)
            flags[scan] |= new[scan]
        finite = frac[scan][np.isfinite(frac[scan])]
        rng = f"{finite.min():.2f} to {finite.max():.2f}" if finite.size else "none"
        print(f"  Scan {scan}: {int(new[scan].sum())} channels above imaginary fraction "
              f"{max_frac:.2f} (fraction {rng})")
    return frac, new


def print_crosshand_table(freq_ghz, raw_flags, snr, snr_fl, frac, imag_fl):
    """
    Per XF channel of one scan: frequency, cross-hand S/N, imaginary fraction
    |Im| / |Re|, and the flag applied (polcal / no data, S/N, imaginary
    fraction, or none).
    """
    print(f"    {'ch':>3} {'GHz':>8} {'S/N':>8} {'|Im|/|Re|':>10} {'flag':>10}")
    for c in range(len(freq_ghz)):
        flag = ('polcal' if raw_flags[c] else 'S/N' if snr_fl[c]
                else 'imag frac' if imag_fl[c] else '-')
        print(f"    {c:>3} {freq_ghz[c]:>8.4f} {snr[c]:>8.1f} {frac[c]:>10.3f} {flag:>10}")


def local_trend_clip(gains, flags, freq_ghz, snr, window, nsigma, max_iter=10):
    """
    Flag XF channels whose phase departs from the local trend of their
    neighbours. For each unflagged channel, a straight line is fitted
    (weighted by S/N) to the phases of its `window` nearest unflagged
    neighbours, so a phase slope is never read as a deviation, also at the
    band edges. The channel is flagged when its offset from the line exceeds

        nsigma x max(1.4826 x MAD of the neighbours about the line, 1 / S/N)

    The 1/S/N floor is the expected phase noise (rad), so at high S/N a
    sub-degree wiggle is not flagged just because the neighbours are smooth.
    Repeats until nothing new is flagged. Returns the updated flags.
    """
    flags = flags.copy()
    phase_noise = np.where(np.isfinite(snr) & (snr > 0), 1.0 / snr, np.inf)
    weight = np.where(np.isfinite(snr), snr, 0.0)

    for iteration in range(1, max_iter + 1):
        good_idx = np.where(~flags & np.isfinite(gains))[0]
        new = np.zeros(len(gains), dtype=bool)
        for ch in good_idx:
            others = good_idx[good_idx != ch]
            if others.size < 3:
                continue
            nearest = np.sort(others[np.argsort(np.abs(others - ch))[:window]])
            # Phases relative to the channel, unwrapped along frequency
            rel = np.unwrap(np.angle(gains[nearest] * np.conj(gains[ch])))
            w = weight[nearest]
            if not np.any(w > 0):
                continue
            slope, offset = np.polyfit(freq_ghz[nearest], rel, 1, w=w)
            resid = rel - (slope * freq_ghz[nearest] + offset)
            scale = max(1.4826 * np.median(np.abs(resid - np.median(resid))), phase_noise[ch])
            # The channel's own relative phase is 0; its offset is minus the line there
            new[ch] = abs(np.angle(np.exp(-1j * (slope * freq_ghz[ch] + offset)))) > nsigma * scale
        flags |= new
        print(f"    Pass {iteration}: {int(new.sum())} flagged, "
              f"{int(np.sum(~flags & np.isfinite(gains)))} usable channels remain")
        if not new.any():
            break
    return flags


# ============================================================================
# XF TABLE
# ============================================================================

def load_xf_gains(xftab, scan_numbers):
    """
    Median XF gain over antennas for each scan. Returns (gains, flags), dicts
    over scans; a channel is flagged when polcal flagged it on every antenna.
    """
    tb.open(xftab)
    cparam = tb.getcol('CPARAM')[0]      # (n_chan, n_rows)
    flag = tb.getcol('FLAG')[0]
    scans_tab = tb.getcol('SCAN_NUMBER')
    tb.close()

    gains, flags = {}, {}
    for scan in scan_numbers:
        rows = scans_tab == scan
        if not rows.any():
            print(f"  WARNING: scan {scan} is not in {xftab}")
            gains[scan] = np.full(cparam.shape[0], np.nan, dtype=complex)
            flags[scan] = np.ones(cparam.shape[0], dtype=bool)
            continue
        g = np.where(flag[:, rows], np.nan, cparam[:, rows])
        with np.errstate(all='ignore'):
            med = np.nanmedian(g.real, axis=1) + 1j * np.nanmedian(g.imag, axis=1)
        flags[scan] = np.all(flag[:, rows], axis=1) | ~np.isfinite(med)
        gains[scan] = np.where(flags[scan], np.nan, med)
        print(f"  Scan {scan}: {int(rows.sum())} antennas, "
              f"{int(np.sum(~flags[scan]))}/{len(med)} channels unflagged")
    return gains, flags


def smooth_gains(gains, flags, scan_numbers):
    """Savitzky-Golay smooth the real and imaginary parts of each scan's unflagged gains, in place."""
    for scan in scan_numbers:
        good = ~flags[scan] & np.isfinite(gains[scan])
        n_good = int(good.sum())
        if n_good < 5:
            continue
        window = XF_SAVGOL_WINDOW if XF_SAVGOL_WINDOW is not None else min(5, n_good)
        window = min(window | 1, n_good if n_good % 2 else n_good - 1)
        order = min(XF_SAVGOL_POLYORDER, window - 1)
        g = gains[scan][good]
        smoothed = (savgol_filter(g.real, window, order)
                    + 1j * savgol_filter(g.imag, window, order))
        gains[scan][good] = smoothed / np.abs(smoothed)
        print(f"  Scan {scan}: {n_good} channels (window {window}, order {order})")


def write_xf_gains(xftab, gains, flags, scan_numbers, scan_times):
    """
    Write each scan's gains (unit amplitude) to every antenna of that scan in
    xftab; channels flagged here become 1+0j and FLAG=True. Each scan's TIME is
    set to its median time. Returns the written gains per scan with flagged
    channels as NaN, for plotting.
    """
    tb.open(xftab, nomodify=False)
    cparam = tb.getcol('CPARAM')
    flag = tb.getcol('FLAG')
    time = tb.getcol('TIME')
    scans_tab = tb.getcol('SCAN_NUMBER')

    written = {}
    for i, scan in enumerate(scan_numbers):
        rows = np.where(scans_tab == scan)[0]
        if rows.size == 0:
            print(f"  WARNING: scan {scan} is not in {xftab}; not written")
            continue
        amp = np.abs(gains[scan])
        f = flags[scan] | ~np.isfinite(gains[scan]) | (amp == 0)
        g = np.where(f, np.nan, gains[scan] / np.where(f, 1.0, amp))
        cparam[0][:, rows] = np.where(f[:, None], 1.0 + 0.0j, g[:, None])
        flag[0][:, rows] |= f[:, None]
        time[rows] = scan_times[i]
        written[scan] = np.where(f | flag[0][:, rows[0]], np.nan, g)
        print(f"  Scan {scan}: {int(np.sum(~f))}/{len(f)} channels written")

    tb.putcol('CPARAM', cparam)
    tb.putcol('FLAG', flag)
    tb.putcol('TIME', time)
    tb.flush()
    tb.close()
    return written


def average_across_scans(xftab, gains, flags, U_xf, V_xf, scan_numbers):
    """
    Overwrite every row of xftab with one solution: the circular mean of the
    per-scan gains, weighted by cross-hand flux sqrt(U^2 + V^2). A channel is
    flagged only when no scan has an unflagged value. Returns that solution
    with flagged channels as NaN.
    """
    phasor_sum = np.zeros(U_xf.shape[1], dtype=complex)
    for i, scan in enumerate(scan_numbers):
        w = np.sqrt(U_xf[i]**2 + V_xf[i]**2)
        good = ~flags[scan] & np.isfinite(gains[scan]) & np.isfinite(w)
        phasor_sum[good] += w[good] * np.exp(1j * np.angle(gains[scan][good]))

    flagged = np.abs(phasor_sum) <= 1e-10
    combined = np.where(flagged, 1.0 + 0.0j, phasor_sum / np.where(flagged, 1.0, np.abs(phasor_sum)))
    print(f"  Combined solution: {int(np.sum(~flagged))}/{len(combined)} channels unflagged")

    tb.open(xftab, nomodify=False)
    cparam = tb.getcol('CPARAM')
    flag = tb.getcol('FLAG')
    cparam[0] = combined[:, None]
    flag[0] = flagged[:, None]
    tb.putcol('CPARAM', cparam)
    tb.putcol('FLAG', flag)
    tb.flush()
    tb.close()
    return np.where(flagged, np.nan, combined)


# ============================================================================
# PLOTS
# ============================================================================

def _scan_colors(n_scans):
    return ['C0'] if n_scans == 1 else plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))


def _legend(ax):
    """Legend outside the panel, to its right, so it never covers data, insets or text boxes."""
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)


def _save(fig, output_dir, filename):
    path = os.path.join(output_dir, filename)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


@styled
def plot_stokes_spectra(freq_ghz, I, Q, U, V, scan_numbers, output_dir, filename, title):
    """I, Q, U, V against frequency, one colour per scan."""
    colors = _scan_colors(len(scan_numbers))
    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True, constrained_layout=True)
    for i, scan in enumerate(scan_numbers):
        label = f'Scan {scan}' if len(scan_numbers) > 1 else None
        for ax, data in zip(axes, (I, Q, U, V)):
            ax.plot(freq_ghz, data[i], '.-', color=colors[i % len(colors)], alpha=0.7,
                    markersize=6, linewidth=0.8, label=label)
    for ax, s in zip(axes, 'IQUV'):
        ax.axhline(0.0, color='gray', linestyle='--', alpha=0.7)
        ax.set_ylabel(f'Stokes {s} [Jy]')
        ax.set_title(f'{title} - {s}')
        ax.grid(True, alpha=0.3)
    if len(scan_numbers) > 1:
        _legend(axes[0])
    axes[3].set_xlabel('Frequency [GHz]')
    _save(fig, output_dir, filename)


@styled
def plot_sign_choice(freq_ghz, raw, signed, Q, U, V, chi, rm_chan, target,
                     scan_numbers, output_dir, evidence):
    """
    Sign choice. 1: raw polcal XF phase. 2: EVPA residual of the kept sign
    (should sit near 0), on the y-scale of the rejected sign, with an inset
    zooming in on the kept residuals. 3: residual of the rejected sign.
    4: the sign margin |res_flip - res_keep| against the confidence threshold;
    filled = decided by the EVPA, open = decided by continuity, red x = the
    joint choice went against the channel's own better sign.
    evidence: {scan: (res_keep, res_flip, threshold, confident, usable)}.
    """
    colors = _scan_colors(len(scan_numbers))
    fig, axes = plt.subplots(4, 1, figsize=(14, 18), sharex=True, constrained_layout=True)
    kept_f, kept_r, kept_c = [], [], []

    for i, scan in enumerate(scan_numbers):
        color = colors[i % len(colors)]
        label = f'Scan {scan}' if len(scan_numbers) > 1 else None
        valid = _valid(signed[scan], Q[i], U[i], V[i], target, rm_chan[scan])
        if valid.sum() < 5:
            continue
        f, s = freq_ghz[valid], signed[scan][valid]
        args = (Q[i][valid], U[i][valid], V[i][valid], chi[i], rm_chan[scan][valid], f)
        res_kept = evpa_residual(derotated_evpa(s, *args), target[valid])
        res_rej = evpa_residual(derotated_evpa(-s, *args), target[valid])

        axes[0].plot(f, np.degrees(np.angle(raw[scan][valid])), '.', color=color,
                     alpha=0.85, markersize=6, label=label)
        axes[1].plot(f, res_kept, '.', color=color, alpha=0.85, markersize=6, label=label)
        axes[2].plot(f, res_rej, '.', color=color, alpha=0.85, markersize=6, label=label)
        kept_f.extend(f)
        kept_r.extend(res_kept)
        kept_c.extend([color] * len(f))

        # Sign margin against the threshold
        res_keep, res_flip, threshold, confident, usable = evidence[scan]
        margin = np.maximum(np.abs(res_flip - res_keep), 1e-3)   # log axis
        flipped =np.real(signed[scan] * np.conj(raw[scan])) < 0
        against = usable & (flipped != (res_flip < res_keep))
        sel_conf, sel_amb = usable & confident, usable & ~confident
        axes[3].plot(freq_ghz[sel_conf], margin[sel_conf], 'o', color=color, markersize=6,
                     label=f'{label + ": " if label else ""}decided by EVPA')
        axes[3].plot(freq_ghz[sel_amb], margin[sel_amb], 'o', color=color, markersize=6,
                     markerfacecolor='none', label=f'{label + ": " if label else ""}decided by continuity')
        axes[3].plot(freq_ghz[usable], threshold[usable], '-', color=color, linewidth=1.2, alpha=0.8,
                     label=f'{label + ": " if label else ""}confidence threshold')
        if against.any():
            axes[3].plot(freq_ghz[against], margin[against], 'x', color='red', markersize=9,
                         markeredgewidth=2, label='Against own better sign' if i == 0 else None)

    axes[0].set_ylabel('XF phase [deg]')
    axes[0].set_title('Raw polcal XF phase')
    axes[1].set_title('EVPA residual, kept sign (y-scale of the rejected sign; inset: zoom)')
    axes[2].set_title('EVPA residual, rejected sign')
    axes[3].set_title('Sign margin |res flip - res keep| and confidence threshold')
    axes[3].set_ylabel('Margin [deg]')
    axes[3].set_yscale('log')
    axes[3].set_xlabel('Frequency [GHz]')
    for ax in axes[1:3]:
        ax.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.set_ylabel('EVPA - model [deg]')
    axes[1].set_ylim(axes[2].get_ylim())
    for ax in axes:
        ax.grid(True, alpha=0.3)
    if len(scan_numbers) > 1:
        _legend(axes[0])
    _legend(axes[3])
    if kept_r:
        # The kept residuals sit near 0, so the inset goes in the half of the
        # panel away from 0
        ylo, yhi = axes[1].get_ylim()
        zero_frac = (0.0 - ylo) / (yhi - ylo)
        inset = axes[1].inset_axes([0.62, 0.06 if zero_frac > 0.5 else 0.54, 0.36, 0.40])
        inset.scatter(kept_f, kept_r, c=kept_c, alpha=0.6, s=8)
        inset.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.6)
        inset.text(0.02, 0.04, 'Kept sign, zoom', transform=inset.transAxes,
                   va='bottom', fontsize=10, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        inset.tick_params(labelsize=9)
        inset.grid(True, alpha=0.2)
    _save(fig, output_dir, PLOT['sign_choice'])


@styled
def plot_evpa_residuals(freq_ghz, gains, flags, Q, U, V, chi, rm_chan, target,
                        scan_numbers, output_dir, filename, stage):
    """
    De-rotated EVPA against the model, its residual against frequency, and a
    histogram of residuals, using the unflagged channels of `gains`.
    """
    colors = _scan_colors(len(scan_numbers))
    fig, axes = plt.subplots(3, 1, figsize=(14, 14), constrained_layout=True)
    all_res = []

    for i, scan in enumerate(scan_numbers):
        color = colors[i % len(colors)]
        label = f'Scan {scan}' if len(scan_numbers) > 1 else None
        valid = ~flags[scan] & _valid(gains[scan], Q[i], U[i], V[i], target, rm_chan[scan])
        if valid.sum() < 5:
            print(f"  WARNING: scan {scan} has fewer than 5 valid channels to plot")
            continue
        f = freq_ghz[valid]
        evpa = derotated_evpa(gains[scan][valid], Q[i][valid], U[i][valid], V[i][valid],
                              chi[i], rm_chan[scan][valid], f)
        res = evpa_residual(evpa, target[valid])
        axes[0].plot(f, evpa, '.', color=color, alpha=0.85, markersize=6, label=label)
        axes[1].plot(f, res, '.', color=color, alpha=0.85, markersize=6, label=label)
        axes[2].hist(res, bins=50, alpha=0.6, color=color, edgecolor='black', linewidth=0.5,
                     label=f'Scan {scan} (median RM {np.median(rm_chan[scan][valid]):+.2f} rad/m^2)')
        all_res.extend(res)

    axes[0].plot(freq_ghz, target, 'k--', linewidth=1.5, alpha=0.8, label='Model EVPA')
    axes[0].set_ylabel('De-rotated EVPA [deg]')
    axes[0].set_title(f'{stage}: EVPA after XF, parallactic angle and Faraday correction')
    _legend(axes[0])
    axes[1].axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.8)
    axes[1].set_ylabel('EVPA - model [deg]')
    axes[1].set_xlabel('Frequency [GHz]')
    axes[1].set_title(f'{stage}: EVPA residual (should be near 0)')
    axes[2].set_xlabel('EVPA - model [deg]')
    axes[2].set_ylabel('Channels')
    axes[2].set_title(f'{stage}: EVPA residuals (should be centred on 0)')

    if all_res:
        all_res = np.array(all_res)
        mean, med, std = all_res.mean(), np.median(all_res), all_res.std()
        axes[2].axvline(0, color='black', linestyle='--', linewidth=2, zorder=10)
        # The median is the reference: robust to the few outlier channels
        axes[2].axvline(med, color='green', linestyle=':', linewidth=2,
                        label=f'Median = {med:+.2f} deg')
        # Statistics as a legend entry, outside the panel
        axes[2].plot([], [], ' ', label=f'Mean = {mean:+.2f} deg\nStd = {std:.2f} deg\nN = {len(all_res)}')
        print(f"  EVPA residual: median {med:+.2f}, mean {mean:+.2f}, std {std:.2f} deg")
    _legend(axes[2])
    for ax in axes:
        ax.grid(True, alpha=0.3)
    _save(fig, output_dir, filename)


@styled
def plot_rm_histogram(rm_chan, coarse_rm, flags, scan_numbers, output_dir):
    """Per-channel RM from the fine search, with the coarse RM (dashed) and median (dotted)."""
    colors = _scan_colors(len(scan_numbers))
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    for i, scan in enumerate(scan_numbers):
        color = colors[i % len(colors)]
        vals = rm_chan[scan][~flags[scan] & np.isfinite(rm_chan[scan])]
        if vals.size == 0:
            continue
        prefix = f'Scan {scan} ' if len(scan_numbers) > 1 else ''
        ax.hist(vals, bins=40, alpha=0.6, color=color, edgecolor='black', linewidth=0.4,
                label=f'{prefix}per-channel RM (N={vals.size})')
        ax.axvline(coarse_rm[scan], color=color, linestyle='--', linewidth=2,
                   label=f'{prefix}coarse RM = {coarse_rm[scan]:+.3f}')
        ax.axvline(np.median(vals), color=color, linestyle=':', linewidth=1.5,
                   label=f'{prefix}median = {np.median(vals):+.3f} +/- {np.std(vals):.3f}')
    ax.set_xlabel('Per-channel RM [rad/m^2]')
    ax.set_ylabel('Channels')
    ax.set_title('Per-channel RM from the fine search')
    _legend(ax)
    ax.grid(True, alpha=0.3)
    _save(fig, output_dir, PLOT['rm_histogram'])


@styled
def plot_psi0_and_v(freq_ghz, Q_sky, U_sky, V, rm_scans, model_evpa, scan_numbers, output_dir):
    """
    Final data against lambda^2: the intrinsic angle psi0 (sky frame, each
    scan's RM removed) with the model EVPA, and Stokes V.
    """
    lambda_sq = (C_LIGHT / (freq_ghz * 1e9)) ** 2
    colors = _scan_colors(len(scan_numbers))
    fig, (ax_psi, ax_v) = plt.subplots(2, 1, figsize=(12, 9), sharex=True, constrained_layout=True)
    shown = [model_evpa]
    for i, scan in enumerate(scan_numbers):
        color = colors[i % len(colors)]
        psi0 = calculate_derotated_angle(Q_sky[i], U_sky[i], rm_scans[i], lambda_sq)
        shown.append(psi0)
        ax_psi.plot(lambda_sq, psi0, '.', color=color, markersize=6,
                    label=f'Scan {scan} (RM {rm_scans[i]:+.2f} rad/m^2)')
        ax_v.plot(lambda_sq, V[i] * 1e3, '.', color=color, markersize=6, label=f'Scan {scan}')
    ax_psi.plot(lambda_sq, model_evpa, 'k--', linewidth=1.5, label='Model EVPA')
    # Zoom on the bulk of the data and the model: 1st-99th percentile plus 25% padding
    vals = np.concatenate(shown)
    lo, hi = np.nanpercentile(vals[np.isfinite(vals)], [1, 99])
    pad = max(0.25 * (hi - lo), 1.0)
    ax_psi.set_ylim(max(lo - pad, -95), min(hi + pad, 95))
    ax_psi.set_ylabel(r'$\psi_0$ [deg]')
    ax_psi.set_title(r'Intrinsic angle $\psi_0$ (sky frame, RM removed) and Stokes V')
    ax_v.axhline(0, color='grey', linestyle=':', linewidth=1.0)
    ax_v.set_ylabel('Stokes V [mJy]')
    ax_v.set_xlabel(r'$\lambda^2$ [m$^2$]')
    for ax in (ax_psi, ax_v):
        ax.grid(True, alpha=0.3)
        _legend(ax)
    _save(fig, output_dir, PLOT['psi0_v'])


@styled
def plot_xf_phase_diagnostic(freq_ghz, raw, signed, cleaned, flags, clip_fl,
                             snr_fl, imag_fl, final, snr, scan_numbers, output_dir):
    """
    XF phase through each step: raw polcal; after the sign choice (filled =
    polcal sign kept, open = flipped); after flagging (x local-trend clip,
    square S/N, triangle imaginary fraction); the written solution; and the
    cross-hand S/N per channel against XF_MIN_SN.
    """
    colors = _scan_colors(len(scan_numbers))
    fig, axes = plt.subplots(5, 1, figsize=(14, 20), sharex=True, constrained_layout=True)
    markers = (('x', 'red', 'Trend clip'), ('s', 'purple', 'S/N flag'),
               ('^', 'darkorange', 'Imaginary-fraction flag'))
    p4 = []

    def deg(g):
        return np.degrees(np.angle(g))

    for i, scan in enumerate(scan_numbers):
        color = colors[i % len(colors)]
        label = f'Scan {scan}' if len(scan_numbers) > 1 else None
        first = i == 0

        ok = np.isfinite(raw[scan])
        axes[0].plot(freq_ghz[ok], deg(raw[scan][ok]), '.', color=color,
                     alpha=0.85, markersize=6, label=label)

        # Every channel's phase after the sign choice, marked by the sign applied
        ok = np.isfinite(signed[scan])
        flipped = ok & (np.real(signed[scan] * np.conj(raw[scan])) < 0)
        prefix = f'Scan {scan}: ' if len(scan_numbers) > 1 else ''
        axes[1].plot(freq_ghz[ok & ~flipped], deg(signed[scan][ok & ~flipped]), 'o', color=color,
                     markersize=5, label=f'{prefix}polcal sign kept')
        axes[1].plot(freq_ghz[flipped], deg(signed[scan][flipped]), 'o', color=color,
                     markersize=5, markerfacecolor='none', label=f'{prefix}flipped')

        good = ~flags[scan] & np.isfinite(cleaned[scan])
        axes[2].plot(freq_ghz[good], deg(cleaned[scan][good]), '.', color=color,
                     alpha=0.85, markersize=6, label=label)

        for fl, (m, c, name) in zip((clip_fl[scan], snr_fl[scan], imag_fl[scan]), markers):
            sel = fl & ok
            if sel.any():
                axes[2].plot(freq_ghz[sel], deg(cleaned[scan][sel]), m, color=c, alpha=0.8,
                             markersize=6, label=f'{name} ({int(sel.sum())} ch)' if first else None)

        ok = np.isfinite(final[scan])
        axes[3].plot(freq_ghz[ok], deg(final[scan][ok]), '-o', color=color, alpha=0.8,
                     markersize=6, linewidth=1, label=label)
        p4.extend(deg(final[scan][ok]))

        ok = np.isfinite(snr[scan])
        axes[4].plot(freq_ghz[ok], snr[scan][ok], '-o', color=color, alpha=0.8,
                     markersize=6, linewidth=1, label=label)
        if (snr_fl[scan] & ok).any():
            axes[4].plot(freq_ghz[snr_fl[scan] & ok], snr[scan][snr_fl[scan] & ok], 's',
                         color='purple', alpha=0.8, markersize=6,
                         label='S/N flagged' if first else None)

    if XF_MIN_SN > 0:
        axes[4].axhline(XF_MIN_SN, color='red', linestyle='--', linewidth=1.5,
                        label=f'XF_MIN_SN = {XF_MIN_SN:.1f}')

    titles = ['Raw polcal XF phase',
              'After sign choice (filled = polcal sign kept, open = flipped)',
              'After flagging (x = trend clip, square = S/N, triangle = imaginary fraction)',
              'Written XF solution',
              'Cross-hand S/N per channel']
    for ax, title in zip(axes, titles):
        ax.set_title(title)
        ax.set_ylabel('XF phase [deg]')
        ax.grid(True, alpha=0.3)
    axes[4].set_ylabel('S/N')
    for ax in axes[:4]:
        ax.set_ylim(-190, 190)
    if p4:
        lo, hi = np.nanmin(p4), np.nanmax(p4)
        margin = max(0.1 * (hi - lo), 2.0)
        axes[3].set_ylim(lo - margin, hi + margin)
    axes[4].set_xlabel('Frequency [GHz]')
    axes[4].set_xlim(freq_ghz.min(), freq_ghz.max())
    axes[4].set_ylim(bottom=0)
    if len(scan_numbers) > 1:
        for ax in (axes[0], axes[3]):
            _legend(ax)
    _legend(axes[1])
    _legend(axes[2])
    _legend(axes[4])
    _save(fig, output_dir, PLOT['xf_phase'])


def caltable_has_field(caltable, vis, field_name):
    """Whether caltable holds a solution for field_name (matched by name via the MS FIELD table)."""
    tb.open(vis + '/FIELD')
    field_names = list(tb.getcol('NAME'))
    tb.close()
    tb.open(caltable)
    field_ids = set(int(i) for i in tb.getcol('FIELD_ID'))
    tb.close()
    return field_name in {field_names[i] for i in field_ids if 0 <= i < len(field_names)}


# ============================================================================
# 1. APPLY CALIBRATION (NO XF)
# ============================================================================

xfdir = GAINPLOTS + '/manualXF'
os.makedirs(xfdir, exist_ok=True)

ktab = GAINTABLES + '/cal_1GC_' + myms + '.K'
bptab = GAINTABLES + '/cal_1GC_' + myms + '.B'
gptab = GAINTABLES + '/cal_1GC_' + myms + '.Gp'
gtab = GAINTABLES + '/cal_1GC_' + myms + '.G'
dftab = GAINTABLES + '/cal_1GC_' + myms + '.Df'
xftab = GAINTABLES + '/cal_1GC_' + myms + '.Xf'

# The K table holds the calibrator's own delays unless 1GC_05 ran with
# DELAY_FROM_PRIMARY, in which case the primary's are interpolated across.
if caltable_has_field(ktab, myms, pacal_name):
    k_field, k_interp = pacal_name, 'nearest'
else:
    k_field, k_interp = bpcal_name, 'linear'
print(f"Delays from {k_field} (interp={k_interp})")

pre_xf_tables = dict(gaintable=[ktab, bptab, gptab, gtab, dftab],
                     gainfield=[k_field, bpcal_name, pacal_name, pacal_name, bpcal_name],
                     interp=[k_interp, 'linear', 'linear', 'linear', 'linear'])

applycal(vis=myms, field=pacal_name, parang=False, flagbackup=False, **pre_xf_tables)

# The 3C286 / 3C138 models are intrinsic EVPAs, so only the ionosphere rotates them
model_name = calibrator_model_name(pacal_name)
rm_intrinsic = 0.0 if model_name else float(XF_TARGET_RM)

print("\n=== Manual XF solve ===")
print(f"Calibrator: {pacal_name}")
if model_name:
    print(f"EVPA: {model_name} model (Perley et al. 2026)")
else:
    print(f"EVPA: XF_TARGET_POLANG = {XF_TARGET_POLANG} deg")
print(f"RM search: {rm_intrinsic:+.2f} rad/m^2 + RM_iono +/- 2 x Spinifex error")
print(f"polcal solution interval: {XF_CHANINT} channels; "
      f"average across scans: {XF_AVG_SCAN}")
print(f"Sign choice: EVPA decides where the signs differ by > max({XF_SIGN_SYS_FLOOR} deg, "
      f"{XF_SIGN_NSIGMA} x EVPA noise), else continuity (jump cost {XF_SIGN_JUMP_COST})")
print(f"Flagging: S/N < {XF_MIN_SN}; |Im| / |Re| > {XF_MAX_IMAG_FRAC}; local-trend clip {XF_TREND_CLIP_NSIGMA} sigma "
      f"(window {XF_TREND_CLIP_WINDOW or 'auto'})")
print(f"Smoothing: {XF_USE_SMOOTHING}")


# ============================================================================
# 2. PER-SCAN IQUV
# ============================================================================

print("\n=== Per-scan IQUV ===")
msmd.open(myms)
try:
    scan_numbers = msmd.scansforfield(msmd.fieldsforname(pacal_name)[0])
finally:
    msmd.close()
n_scans = len(scan_numbers)
print(f"{n_scans} scan(s) of {pacal_name}: {scan_numbers}")

tmpms = channel_average_ms(myms, pacal_name)
vis_to_load = tmpms or myms
freq_ghz = read_chan_freq_ghz(vis_to_load)
print(f"{len(freq_ghz)} channels, {freq_ghz.min():.3f} - {freq_ghz.max():.3f} GHz")

vis_avg, scan_times = load_stokes_per_scan(vis_to_load, pacal_name, scan_numbers, freq_ghz,
                                           os.path.join(xfdir, 'stokes_perscan.npz'))
if tmpms:
    shutil.rmtree(tmpms)

# A calibrated point source at the phase centre has real Stokes visibilities
I_scans, Q_scans, U_scans, V_scans = np.real(vis_avg)
report_iquv_quality(vis_avg, scan_numbers)

if model_name:
    plot_evpa_model(model_name, freq_ghz, calibrator_evpa(model_name, freq_ghz), xfdir)

if not os.path.isfile(os.path.join(xfdir, PLOT['stokes_pre'])):
    plot_stokes_spectra(freq_ghz, I_scans, Q_scans, U_scans, V_scans, scan_numbers, xfdir,
                        PLOT['stokes_pre'], 'Pre-XF Stokes')


# ============================================================================
# 3. PARALLACTIC ANGLE
# ============================================================================

print("\n=== Parallactic angle ===")
chi_deg_scans = np.zeros(n_scans)
for i, scan in enumerate(scan_numbers):
    chi_deg_scans[i], az, el = compute_parallactic_angle(myms, pacal_name, scan_times[i])
    print(f"  Scan {scan}: chi {chi_deg_scans[i]:+7.2f} deg (az {az:6.1f}, el {el:5.1f})")


# ============================================================================
# 4. IONOSPHERIC RM
# ============================================================================

print("\n=== Ionospheric RM (Spinifex) ===")
iono_rm_per_scan, iono_rm_err_per_scan = get_spinifex_rm_per_scan(pacal_name, scan_numbers)


# ============================================================================
# 5. POLCAL XF
# ============================================================================

print(f"\n=== polcal Xf ({XF_CHANINT} channels per solution) ===")
if os.path.isdir(xftab):
    shutil.rmtree(xftab)
polcal(vis=myms, field=pacal_name, caltable=xftab, refant=str(ref_ant),
       solint=f'inf,{XF_CHANINT}ch', poltype='Xf', combine='', append=False,
       **pre_xf_tables)

freq_xf_ghz = read_chan_freq_ghz(xftab)
n_xf_chan = len(freq_xf_ghz)
raw_gains, raw_flags = load_xf_gains(xftab, scan_numbers)

# IQUV and the model EVPA on the XF channel grid
Q_xf, U_xf, V_xf = (np.array([np.interp(freq_xf_ghz, freq_ghz[np.isfinite(s)], s[np.isfinite(s)])
                              if np.sum(np.isfinite(s)) > 1 else np.full(n_xf_chan, np.nan)
                              for s in stokes])
                    for stokes in (Q_scans, U_scans, V_scans))
target_xf = calibrator_evpa(model_name, freq_xf_ghz)

# Channels without IQUV or a model EVPA cannot be sign-checked
for i, scan in enumerate(scan_numbers):
    unusable = ~_valid(Q_xf[i], U_xf[i], V_xf[i], target_xf)
    n_new = int(np.sum(unusable & ~raw_flags[scan]))
    if n_new:
        print(f"  Scan {scan}: flagging {n_new} channels with no IQUV or model EVPA")
    raw_flags[scan] |= unusable
    raw_gains[scan][unusable] = np.nan


# ============================================================================
# 6. CROSS-HAND S/N AND CONTAMINATION
# ============================================================================
# Flags the channels where the polcal phase is set by noise (S/N) or by
# contamination of the cross-hands (imaginary fraction), and gives the EVPA
# noise that sets how far apart the two signs must be to be told apart

print(f"\n=== Cross-hand S/N (XF_MIN_SN = {XF_MIN_SN}) ===")
flags = {scan: raw_flags[scan].copy() for scan in scan_numbers}
snr_per_scan, snr_flagged, xh_noise = snr_flag(U_scans, V_scans, flags, scan_numbers,
                                               n_xf_chan, XF_MIN_SN)
for scan in scan_numbers:
    print(f"  Scan {scan}: noise per XF channel {np.nanmedian(xh_noise[scan]):.4e} Jy (median)")

print(f"\n=== Cross-hand imaginary fraction (XF_MAX_IMAG_FRAC = {XF_MAX_IMAG_FRAC}) ===")
imag_frac_per_scan, imag_flagged = imag_fraction_flag(vis_avg, flags, scan_numbers,
                                                      n_xf_chan, XF_MAX_IMAG_FRAC)
for scan in scan_numbers:
    print(f"\n  Scan {scan}:")
    print_crosshand_table(freq_xf_ghz, raw_flags[scan], snr_per_scan[scan], snr_flagged[scan],
                          imag_frac_per_scan[scan], imag_flagged[scan])


# ============================================================================
# 7. RESOLVE THE 180 DEG AMBIGUITY
# ============================================================================
# Each channel's EVPA with the polcal sign kept and flipped is compared with
# the model. Where they differ by more than the threshold, max(XF_SIGN_SYS_FLOOR,
# XF_SIGN_NSIGMA x EVPA noise), the EVPA decides; elsewhere (weak cross-hand
# signal, so the sign barely moves the EVPA) the sign follows continuity. The
# RM is fitted together with the signs, so a sign pattern cannot absorb an
# RM error.

print("\n=== RM and sign of the XF solution ===")
signed_gains, rm_chan_per_scan, coarse_rm_per_scan, evidence = {}, {}, {}, {}
for i, scan in enumerate(scan_numbers):
    print(f"\n  Scan {scan} (chi {chi_deg_scans[i]:+.2f} deg, "
          f"RM_iono {iono_rm_per_scan[i]:+.3f} +/- {iono_rm_err_per_scan[i]:.3f} rad/m^2)")
    usable = ~flags[scan] & _valid(raw_gains[scan], Q_xf[i], U_xf[i], V_xf[i], target_xf)
    threshold = np.fmax(XF_SIGN_SYS_FLOOR,
                        XF_SIGN_NSIGMA * evpa_noise_deg(Q_xf[i], U_xf[i], xh_noise[scan]))
    (rm_chan_per_scan[scan], coarse_rm_per_scan[scan],
     flip, confident, res_keep, res_flip) = rm_search(
        raw_gains[scan], usable, Q_xf[i], U_xf[i], V_xf[i], chi_deg_scans[i], freq_xf_ghz,
        target_xf, iono_rm_per_scan[i], iono_rm_err_per_scan[i], rm_intrinsic,
        threshold, XF_SIGN_JUMP_COST)
    # Flagged channels take their own better sign, so they can still be plotted
    flip = np.where(usable, flip, res_flip < res_keep)
    signed_gains[scan] = np.where(flip, -raw_gains[scan], raw_gains[scan])
    evidence[scan] = (res_keep, res_flip, threshold, confident, usable)

    n_use, n_conf = int(usable.sum()), int(confident.sum())
    n_over = int(np.sum(confident & (flip != (res_flip < res_keep))))
    print(f"    {n_use} usable channels: {n_conf} decided by the EVPA, {n_use - n_conf} by "
          f"continuity; {int(np.sum(flip & usable))} flipped; {n_over} confident channel(s) "
          f"overruled by continuity")
    ambiguous = contiguous_ranges(freq_xf_ghz[usable], ~confident[usable])
    if n_use and (n_use - n_conf) / n_use >= 0.25:
        print(f"    WARNING: the EVPA cannot tell the two signs apart over "
              f"{100 * (n_use - n_conf) / n_use:.0f}% of the band ({', '.join(ambiguous)}). "
              f"The cross-hand signal there is a small part of the polarized flux, so the "
              f"XF sign follows continuity only.")
    elif ambiguous:
        print(f"    Decided by continuity: {', '.join(ambiguous)}")
    print_sign_table(freq_xf_ghz, raw_gains[scan], flip, usable, res_keep, res_flip,
                     threshold, confident, target_xf)


def save_sign_choice():
    """
    Save the XF phases before and after the sign choice, with every input to
    it, to xf_sign_choice.npz (arrays are scans x XF channels), so the sign
    resolution can be rerun and tested offline.
    """
    def stack(per_scan):
        return np.array([per_scan[scan] for scan in scan_numbers])
    path = os.path.join(xfdir, 'xf_sign_choice.npz')
    ev = {name: np.array([evidence[scan][k] for scan in scan_numbers])
          for k, name in enumerate(('res_keep', 'res_flip', 'threshold', 'confident', 'usable'))}
    np.savez(path, freq_xf_ghz=freq_xf_ghz, scan_numbers=scan_numbers,
             raw_gains=stack(raw_gains), raw_flags=stack(raw_flags),
             signed_gains=stack(signed_gains), rm_chan=stack(rm_chan_per_scan),
             coarse_rm=stack(coarse_rm_per_scan), Q_xf=Q_xf, U_xf=U_xf, V_xf=V_xf,
             chi_deg_scans=chi_deg_scans, model_evpa_deg=target_xf,
             iono_rm=iono_rm_per_scan, iono_rm_err=iono_rm_err_per_scan,
             rm_intrinsic=rm_intrinsic, crosshand_snr=stack(snr_per_scan),
             snr_flagged=stack(snr_flagged), crosshand_noise=stack(xh_noise),
             imag_frac=stack(imag_frac_per_scan), imag_flagged=stack(imag_flagged), **ev)
    print(f"Saved: {path}")


save_sign_choice()
geometry = (Q_xf, U_xf, V_xf, chi_deg_scans, rm_chan_per_scan, target_xf, scan_numbers, xfdir)
plot_sign_choice(freq_xf_ghz, raw_gains, signed_gains, *geometry, evidence)
plot_evpa_residuals(freq_xf_ghz, signed_gains, raw_flags, *geometry,
                    PLOT['resid_signed'], 'After sign choice')


# ============================================================================
# 8. LOCAL-TREND CLIP
# ============================================================================

gains = {scan: signed_gains[scan].copy() for scan in scan_numbers}
clip_window = int(XF_TREND_CLIP_WINDOW) if XF_TREND_CLIP_WINDOW else max(5, int(round(n_xf_chan / 10.0)))
clip_flagged = {scan: np.zeros(n_xf_chan, dtype=bool) for scan in scan_numbers}
print(f"\n=== Local-trend clip (window {clip_window} channels, {XF_TREND_CLIP_NSIGMA} sigma) ===")
if XF_TREND_CLIP_NSIGMA > 0:
    for scan in scan_numbers:
        print(f"  Scan {scan}:")
        before = flags[scan].copy()
        flags[scan] = local_trend_clip(gains[scan], flags[scan], freq_xf_ghz, snr_per_scan[scan],
                                       clip_window, XF_TREND_CLIP_NSIGMA)
        clip_flagged[scan] = flags[scan] & ~before
else:
    print("  Disabled (XF_TREND_CLIP_NSIGMA <= 0)")

cleaned_gains = {scan: gains[scan].copy() for scan in scan_numbers}
plot_evpa_residuals(freq_xf_ghz, gains, flags, *geometry,
                    PLOT['resid_flagged'], 'After flagging')
plot_rm_histogram(rm_chan_per_scan, coarse_rm_per_scan, flags, scan_numbers, xfdir)


# ============================================================================
# 9. SMOOTHING
# ============================================================================

if XF_USE_SMOOTHING:
    print("\n=== Savitzky-Golay smoothing ===")
    smooth_gains(gains, flags, scan_numbers)


# ============================================================================
# 10. WRITE THE XF TABLE
# ============================================================================

print(f"\n=== Writing {xftab} ===")
final_gains = write_xf_gains(xftab, gains, flags, scan_numbers, scan_times)

if XF_AVG_SCAN and n_scans > 1:
    print("\n=== Flux-weighted average across scans ===")
    combined = average_across_scans(xftab, gains, flags, U_xf, V_xf, scan_numbers)
    final_gains = {scan: combined for scan in scan_numbers}


# ============================================================================
# 11. DIAGNOSTICS
# ============================================================================

print("\n=== Diagnostics ===")
plot_xf_phase_diagnostic(freq_xf_ghz, raw_gains, signed_gains, cleaned_gains, flags,
                         clip_flagged, snr_flagged, imag_flagged, final_gains, snr_per_scan,
                         scan_numbers, xfdir)

# The written solution shifted by pi in every channel (the other polcal
# solution), for comparison with plot 05
plot_evpa_residuals(freq_xf_ghz, {scan: -final_gains[scan] for scan in scan_numbers},
                    {scan: ~np.isfinite(final_gains[scan]) for scan in scan_numbers},
                    *geometry, PLOT['resid_opposite'],
                    'Final + pi (for comparison)')

# Analytic XF correction of the per-scan IQUV; CASA applies conj(gain)
U_corr = np.full_like(U_scans, np.nan)
V_corr = np.full_like(V_scans, np.nan)
for i, scan in enumerate(scan_numbers):
    ok = np.isfinite(final_gains[scan])
    if ok.sum() < 2:
        continue
    g = final_gains[scan][ok]
    g_interp = (np.interp(freq_ghz, freq_xf_ghz[ok], g.real)
                + 1j * np.interp(freq_ghz, freq_xf_ghz[ok], g.imag))
    U_corr[i], V_corr[i] = correct_crosshand_phase(U_scans[i], V_scans[i],
                                                   np.angle(np.conj(g_interp)))
    with np.errstate(invalid='ignore', divide='ignore'):
        p_frac = 100 * np.nanmedian(np.hypot(Q_scans[i], U_corr[i]) / I_scans[i])
        v_frac = 100 * np.nanmedian(np.abs(V_corr[i]) / I_scans[i])
    print(f"  Scan {scan}: P/I {p_frac:.3f}%, |V|/I {v_frac:.3f}%")

plot_stokes_spectra(freq_ghz, I_scans, Q_scans, U_corr, V_corr, scan_numbers, xfdir,
                    PLOT['stokes_feed'], 'Post-XF Stokes (analytic, feed frame)')
# Saved in the feed frame (XF phase removed only) and in the sky frame (Q, U
# also rotated by -2 chi), with what is needed to compare the sky-frame EVPA,
# 0.5 arctan2(U_sky, Q_sky) - RM lambda^2, with the model
Q_sky = np.full_like(Q_scans, np.nan)
U_sky = np.full_like(U_scans, np.nan)
for i in range(n_scans):
    Q_sky[i], U_sky[i] = correct_parallactic_angle(Q_scans[i], U_corr[i], chi_deg_scans[i])
rm_scans = np.array([np.nanmedian(rm_chan_per_scan[scan]) if np.any(np.isfinite(rm_chan_per_scan[scan]))
                     else np.nan for scan in scan_numbers])
plot_stokes_spectra(freq_ghz, I_scans, Q_sky, U_sky, V_corr, scan_numbers, xfdir,
                    PLOT['stokes_sky'], 'Post-XF Stokes (analytic, sky frame)')
np.savez(os.path.join(xfdir, 'final_corrected_stokes_spectra.npz'),
         freq_ghz=freq_ghz, scan_numbers=scan_numbers, I_corrected_scans=I_scans,
         Q_corrected_scans=Q_scans, U_corrected_scans=U_corr, V_corrected_scans=V_corr,
         chi_deg_scans=chi_deg_scans, Q_sky_scans=Q_sky, U_sky_scans=U_sky,
         rm_scans=rm_scans, model_evpa_deg=calibrator_evpa(model_name, freq_ghz))
print(f"Saved: {xfdir}/final_corrected_stokes_spectra.npz")
plot_psi0_and_v(freq_ghz, Q_sky, U_sky, V_corr, rm_scans, calibrator_evpa(model_name, freq_ghz),
                scan_numbers, xfdir)

if XF_APPLY_TO_MS:
    print("\n=== Applying the XF table to the MS ===")
    # parang=True: CORRECTED_DATA is rotated to the sky frame, for comparison with plot 10
    applycal(vis=myms, field=pacal_name, parang=True, flagbackup=False,
             gaintable=pre_xf_tables['gaintable'] + [xftab],
             gainfield=pre_xf_tables['gainfield'] + [pacal_name],
             interp=pre_xf_tables['interp'] + ['linear'])
    freq_ms_ghz = read_chan_freq_ghz(myms)
    iquv_ms = np.full((4, n_scans, len(freq_ms_ghz)), np.nan)
    for i, scan in enumerate(scan_numbers):
        result = read_scan_iquv(myms, pacal_name, scan, 'corrected_data')
        if result is not None:
            iquv_ms[:, i, :] = np.real(result[0])
    plot_stokes_spectra(freq_ms_ghz, *iquv_ms, scan_numbers, xfdir,
                        PLOT['stokes_casa'], 'Post-XF Stokes (CORRECTED_DATA, sky frame)')

write_plot_readme(xfdir)
print("\n=== Manual XF solve complete ===")
