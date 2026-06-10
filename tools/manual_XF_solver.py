# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
Manual Cross-Hand Phase (XF) Polarization Calibration Solver

Solves for the instrumental cross-hand phase (XF) using a known polarization
angle calibrator (e.g. 3C286). Operates directly at the working channelisation
(XF_CHANINT channels) — CASA's polcal solver maximises SNR within each solution
interval. No full-resolution solve or manual frequency averaging is performed.

SCIENTIFIC BACKGROUND:
The XF term represents the phase difference between the X and Y feed receptors.
It rotates observed Stokes U and V in the cross-hand correlations, and must be
removed before accurate polarization measurements can be made.

The ±π degeneracy in CASA's polcal Xf solution is resolved by comparing the
de-rotated EVPA (after applying the trial gain, parallactic angle correction,
and Faraday de-rotation) to the known source EVPA model. The ionospheric RM is
estimated per scan using Spinifex and folded into the RM trial grid.

WORKFLOW:
1.  applycal on myms with [K, B, Gp, G, Df] tables (no XF yet)
2.  Load IQUV visibilities per scan via weighted average of CORRECTED_DATA
3.  Compute parallactic angles at each scan's mid-time
4.  Run Spinifex to estimate ionospheric RM per scan
5.  polcal at XF_CHANINT ch resolution → xftab
6.  Load xftab gains; flag channels with no valid IQUV
7.  Interpolate IQUV onto xftab frequency grid for ±π resolution
8.  Coarse RM grid search (XF_TARGET_RM + iono_RM ± 2.0 rad/m²) to find
    global best RM per scan
9.  Fine RM grid (±0.25 rad/m²) per channel independently — each channel
    adopts the (RM, sign) pair minimising |EVPA − target|
10. Population sign check: flip channels inconsistent with band-wide median
11. 80% sanity check: force-flip scan if <80% of channels agree with target
12. Global XF_GLOBAL_SIGMA_CLIP σ circular phase clip (wrap-safe)
13. Polynomial baseline + sliding circular MAD clip (wrap-safe):
      - Fit cos/sin(phase) vs normalised channel index (poly order XF_POLY_ORDER)
      - Subtract baseline phasor; flag channels where baseline-subtracted
        deviation > XF_MAD_SIGMA_CLIP σ × 1.4826 × local circular MAD
        (window = XF_MAD_WINDOW channels)
14. Optional Savitzky-Golay smoothing
15. Write corrected gains back to xftab (FLAG column preserved from polcal)
16. Interpolate xftab gains onto full MS channel grid for IQUV diagnostic
    correction and Stokes before/after plots

DIAGNOSTIC OUTPUTS (all in GAINPLOTS/manualXF/):
    stokes_perscan.npz              Cached IQUV visibilities
    3c286_evpa_model.png            EVPA model vs frequency
    stokes_spectra_preXF.png        Pre-XF IQUV spectra
    xf_phase_stage1_raw.png         Raw polcal XF phases + EVPA residuals
    xf_phase_stage2_post_pi.png     Post ±π resolution residuals
    xf_phase_stage3_post_flagging.png  Post-clipping residuals
    xf_phase_diagnostic.png         4-panel: raw → ±π → clipped → final
    xf_per_channel_rm_histogram.png Per-channel adopted RM distribution
    stokes_spectra_postXF_analytic.png  Post-XF IQUV (analytic correction)
"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import shutil
import datetime
import subprocess

# Smoothing library
from scipy.signal import savgol_filter

# Configure immediate output flushing for better logging
import functools
print = functools.partial(print, flush=True)

# Load oxkat configuration and project information
exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

# Handle user-defined fields if provided
if PRE_FIELDS != '':
    targets = user_targets
    pcal_names = user_pcals
    target_cal_map = user_cal_map

# ============================
# USER GLOBALS
# ============================

# Whether to run a final applycal with the new xftab and extract IQUV from the
# MS for the post-XF Stokes diagnostic plot. False = use analytic correction only.
XF_APPLY_TO_MS = False


# ============================
# HELPER FUNCTIONS
# ============================


def build_iquv_flag(raw_flag, nchan, nrow):
    """
    Convert raw FLAG array to IQUV-compatible flag array.

    Inputs:
        raw_flag : array - Raw FLAG array from MS (various shapes supported)
        nchan : int - Number of frequency channels
        nrow : int - Number of rows (baselines × times)

    Outputs:
        flag_iquv : array - Flag array with shape (4, nchan, nrow) for I,Q,U,V

    Maps correlation flags to Stokes parameter flags:
        - I,Q flagged if XX or YY flagged
        - U,V flagged if XY or YX flagged
    """
    if raw_flag.ndim == 3:
        c = raw_flag.shape[0]
        print(f"FLAG has {c} correlation planes with shape {raw_flag.shape}.")
        if c == 4:
            print("Converting 4-corr flags (XX,XY,YX,YY) → I,Q,U,V flags.")
            f_xx = raw_flag[0, :, :]
            f_xy = raw_flag[1, :, :]
            f_yx = raw_flag[2, :, :]
            f_yy = raw_flag[3, :, :]
            f_iq = f_xx | f_yy
            f_uv = f_xy | f_yx
            return np.stack([f_iq, f_iq, f_uv, f_uv], axis=0)
        elif c == 1:
            print("Single-plane flags; broadcasting to I,Q,U,V.")
            return np.broadcast_to(raw_flag, (4, nchan, nrow))
        else:
            raise ValueError(f"Unexpected number of correlation planes: {c}. Expected 1 or 4.")
    elif raw_flag.ndim == 2 and raw_flag.shape == (nchan, nrow):
        print("Per-channel flags (nchan,nrow); broadcasting to I,Q,U,V.")
        return np.broadcast_to(raw_flag[None, :, :], (4, nchan, nrow))
    else:
        raise ValueError(f"Unrecognised FLAG shape: {raw_flag.shape}.")


def build_iquv_weights(weight_data, have_wspec, nchan, nrow):
    """
    Convert raw WEIGHT or WEIGHT_SPECTRUM to IQUV-compatible weight array.

    Inputs:
        weight_data : dict - Dictionary containing 'weight_spectrum' or 'weight'
        have_wspec : bool - True if WEIGHT_SPECTRUM available
        nchan : int - Number of frequency channels
        nrow : int - Number of rows

    Outputs:
        W : array - Weight array with shape (4, nchan, nrow) for I,Q,U,V
    """
    if have_wspec:
        W = np.asarray(weight_data['weight_spectrum'])
        print(f"Using WEIGHT_SPECTRUM with original shape {W.shape}")
        if W.ndim == 2 and W.shape == (nchan, nrow):
            print("Single weight per channel → broadcasting to (4, nchan, nrow).")
            W = np.broadcast_to(W[np.newaxis, :, :], (4, nchan, nrow))
        elif W.ndim == 3:
            c = W.shape[0]
            if c == 4:
                print("Converting 4-corr WEIGHT_SPECTRUM → I,Q,U,V weights.")
                w_xx = W[0, :, :]
                w_xy = W[1, :, :]
                w_yx = W[2, :, :]
                w_yy = W[3, :, :]
                w_iq = 0.5 * (w_xx + w_yy)
                w_uv = 0.5 * (w_xy + w_yx)
                W = np.stack([w_iq, w_iq, w_uv, w_uv], axis=0)
            elif c == 1:
                print("Single-pol WEIGHT_SPECTRUM → broadcasting to (4, nchan, nrow).")
                W = np.broadcast_to(W, (4, nchan, nrow))
            else:
                raise ValueError(f"Unexpected correlations in WEIGHT_SPECTRUM: {c}.")
        else:
            raise ValueError(f"Unexpected WEIGHT_SPECTRUM dimensions: {W.shape}.")
    else:
        Wraw = np.asarray(weight_data['weight'])
        print(f"Using WEIGHT with original shape {Wraw.shape}")
        if Wraw.ndim == 1 and Wraw.shape[0] == nrow:
            print("Single weight per row → broadcasting to (4, nchan, nrow).")
            W = np.ones((4, nchan, nrow), dtype=float) * Wraw[None, None, :]
        elif Wraw.ndim == 2 and Wraw.shape[-1] == nrow:
            c = Wraw.shape[0]
            if c == 4:
                print("Converting 4-corr WEIGHT → I,Q,U,V weights.")
                w_iq = 0.5 * (Wraw[0, :] + Wraw[3, :])
                w_uv = 0.5 * (Wraw[1, :] + Wraw[2, :])
                W = np.zeros((4, nchan, nrow), dtype=float)
                W[0, :, :] = w_iq[None, :]
                W[1, :, :] = w_iq[None, :]
                W[2, :, :] = w_uv[None, :]
                W[3, :, :] = w_uv[None, :]
            elif c == 1:
                print("Single correlation WEIGHT → broadcasting to (4, nchan, nrow).")
                W = np.ones((4, nchan, nrow), dtype=float) * Wraw[0, :][None, None, :]
            else:
                raise ValueError(f"Unexpected correlations in WEIGHT: {c}.")
        else:
            raise ValueError(f"Unexpected WEIGHT shape: {Wraw.shape}.")
    return W


def compute_weighted_averages(data, weights, flags):
    """
    Compute weighted vector averages per channel, accounting for flags.

    Inputs:
        data : array - Complex visibility data (4, nchan, nrow)
        weights : array - Weight array (4, nchan, nrow)
        flags : array - Flag array (4, nchan, nrow)

    Outputs:
        vis_avg : array - Complex averaged visibilities (4, nchan)
    """
    W_eff = np.where(flags, 0.0, weights)
    den = np.sum(W_eff, axis=-1)
    num = np.sum(W_eff * data, axis=-1, dtype=np.complex128)
    with np.errstate(invalid='ignore', divide='ignore'):
        vis_avg = num / den
    return vis_avg


def compute_3c286_evpa(freq_ghz):
    """
    Compute frequency-dependent EVPA for 3C286 (J1331+3030).

    Reference: Perley & Butler (2013) / Hugo & Perley (2024)

    EVPA(ν) [deg] =
        32.64 - 85.37λ²                              for ν ∈ [1.7, 12] GHz
        29.53 + λ²(4005.88(log₁₀ν)³ - 39.38)        for ν < 1.7 GHz
    """
    c = 2.99792458e8
    freq_hz = freq_ghz * 1e9
    lambda_m = c / freq_hz
    lambda_sq = lambda_m ** 2
    evpa_deg = np.zeros_like(freq_ghz)
    high_freq_mask = freq_ghz >= 1.7
    evpa_deg[high_freq_mask] = 32.64 - 85.37 * lambda_sq[high_freq_mask]
    low_freq_mask = freq_ghz < 1.7
    if np.any(low_freq_mask):
        log_nu_cubed = np.log10(freq_ghz[low_freq_mask]) ** 3
        evpa_deg[low_freq_mask] = 29.53 + lambda_sq[low_freq_mask] * (4005.88 * log_nu_cubed - 39.38)
    return evpa_deg


def quick_stats(name, arr):
    """Print quick statistics (median, MAD, std) for a numerical array."""
    finite = np.isfinite(arr)
    if not np.any(finite):
        print(f'{name}: no finite values')
        return
    med = np.nanmedian(arr[finite])
    mad = np.nanmedian(np.abs(arr[finite] - med))
    std = np.nanstd(arr[finite])
    print(f'{name}: median={med:.6g}, MAD={mad:.6g}, std={std:.6g}')


def calc_im_fraction(complex_avg):
    """Calculate fraction of flux density in imaginary component (%)."""
    real_part = np.abs(np.real(complex_avg))
    imag_part = np.abs(np.imag(complex_avg))
    total_amp = np.sqrt(real_part ** 2 + imag_part ** 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        fraction = imag_part / total_amp
    return np.nanmedian(fraction) * 100


def q_to_rad_scalar(q):
    """Convert CASA quantity to radians, returning a Python float."""
    return float(np.atleast_1d(qa.getvalue(qa.convert(q, 'rad')))[0])


def rad_to_deg_scalar(xrad):
    """Convert radians to degrees using CASA, returning a Python float."""
    return float(np.atleast_1d(qa.getvalue(qa.convert({'value': xrad, 'unit': 'rad'}, 'deg')))[0])


def compute_parallactic_angle(vis, field_name, time_mjd, field_id=None):
    """
    Compute parallactic angle at specified time using CASA AZ/EL method.

    χ = atan2(-sin(A), tan(φ)cos(e) - cos(A)sin(e))
    where A=azimuth, e=elevation, φ=site latitude.

    Returns chi_deg (float) and diagnostics (dict).
    """
    msmd.open(vis)
    try:
        fids = msmd.fieldsforname(field_name)
        field_id = int(fids[0]) if len(fids) else None
    except Exception:
        field_id = None

    if field_id is None:
        tb.open(vis)
        field_ids_all = tb.getcol('FIELD_ID')
        tb.close()
        vals, cnts = np.unique(field_ids_all, return_counts=True)
        field_id = int(vals[np.argmax(cnts)])
        print(f'WARNING: field "{field_name}" not found; using modal FIELD_ID {field_id}')

    phase_dir = msmd.phasecenter(field_id)
    msmd.close()

    tb.open(vis + '/OBSERVATION')
    tel_names = tb.getcol('TELESCOPE_NAME')
    tb.close()
    obsname = str(tel_names[0]) if tel_names.size > 0 else 'UNKNOWN'

    pos_meas = me.observatory(obsname)
    pos_wgs = me.measure(pos_meas, 'wgs84')
    lat_rad = q_to_rad_scalar(pos_wgs['m1'])

    me.doframe(pos_meas)
    me.doframe(me.epoch('utc', qa.quantity(time_mjd, 's')))
    azel = me.measure(phase_dir, 'azel')

    az_rad = q_to_rad_scalar(azel['m0'])
    el_rad = q_to_rad_scalar(azel['m1'])

    num = -np.sin(az_rad)
    den = np.tan(lat_rad) * np.cos(el_rad) - np.cos(az_rad) * np.sin(el_rad)
    chi_rad = np.arctan2(num, den)
    chi_deg = rad_to_deg_scalar(chi_rad)
    chi_deg = ((chi_deg + 180.0) % 360.0) - 180.0

    diagnostics = {
        'az_deg': rad_to_deg_scalar(az_rad),
        'el_deg': rad_to_deg_scalar(el_rad),
        'lat_deg': rad_to_deg_scalar(lat_rad),
        'time': time_mjd,
        'observatory': obsname
    }
    return chi_deg, diagnostics


def correct_crosshand_phase(u_prime, v_prime, rho):
    """
    Apply cross-hand phase correction.

    U_corr = U'cos(ρ) - V'sin(ρ)
    V_corr = U'sin(ρ) + V'cos(ρ)
    """
    cos_rho = np.cos(rho)
    sin_rho = np.sin(rho)
    u_corrected = u_prime * cos_rho - v_prime * sin_rho
    v_corrected = u_prime * sin_rho + v_prime * cos_rho
    return u_corrected, v_corrected


def correct_parallactic_angle(q, u, parang_deg):
    """
    Apply parallactic angle correction (feed → sky frame).

    Rotates Q,U from feed frame to sky frame by -2χ.
    """
    parang_rad = np.radians(-2 * parang_deg)
    q_corrected = q * np.cos(parang_rad) + u * np.sin(parang_rad)
    u_corrected = -q * np.sin(parang_rad) + u * np.cos(parang_rad)
    return q_corrected, u_corrected


def calculate_derotated_angle(q, u, rm, lambda_sq_val):
    """
    Calculate de-rotated polarization angle after RM correction.

    Applies Faraday de-rotation: θ = 2·RM·λ², then returns
    0.5·arctan2(U', Q') in degrees, wrapped to (-90, 90].
    """
    rot_angle = 2 * rm * lambda_sq_val
    cos_rot = np.cos(rot_angle)
    sin_rot = np.sin(rot_angle)
    q_derot = q * cos_rot + u * sin_rot
    u_derot = -q * sin_rot + u * cos_rot
    pol_angle_rad = 0.5 * np.arctan2(u_derot, q_derot)
    pol_angle_deg = np.degrees(pol_angle_rad)
    pol_angle_deg = ((pol_angle_deg + 90.0) % 180.0) - 90.0
    return pol_angle_deg


def plot_stokes_spectra(freq, I_scans, Q_scans, U_scans, V_scans, scan_numbers,
                        output_dir, filename='stokes_spectra.png',
                        title='Stokes Parameters', zoom_percentile=None,
                        zoom_stokes=None):
    """
    Generate per-channel spectra for all Stokes parameters with multiple scans.

    Creates 4-panel plot (I, Q, U, V) vs frequency. Each scan in a different
    colour. Optional percentile zoom on selected Stokes panels.
    """
    n_scans = I_scans.shape[0]
    colors = ['C0'] if n_scans == 1 else plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))

    fig, ax = plt.subplots(4, 1, figsize=(12, 14), sharex=True, constrained_layout=True)

    for scan_idx in range(n_scans):
        color = colors[scan_idx % len(colors)]
        label = f'Scan {scan_numbers[scan_idx]}' if n_scans > 1 else None
        kw = dict(marker='.', linestyle='-', alpha=0.7, color=color,
                  label=label, markersize=6, linewidth=0.8)
        ax[0].plot(freq, I_scans[scan_idx], **kw)
        ax[1].plot(freq, Q_scans[scan_idx], **kw)
        ax[2].plot(freq, U_scans[scan_idx], **kw)
        ax[3].plot(freq, V_scans[scan_idx], **kw)

    if zoom_percentile is not None:
        lower_p = (100 - zoom_percentile) / 2
        upper_p = 100 - lower_p
        stokes_names = ['I', 'Q', 'U', 'V']
        for idx, (data, sname) in enumerate(
                zip([I_scans, Q_scans, U_scans, V_scans], stokes_names)):
            if zoom_stokes is None or sname in zoom_stokes:
                valid_data = data[np.isfinite(data)]
                if len(valid_data) > 0:
                    y_min = np.percentile(valid_data, lower_p)
                    y_max = np.percentile(valid_data, upper_p)
                    y_range = y_max - y_min
                    ax[idx].set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

    for a, ylabel, t in zip(ax,
                             ['Stokes I [Jy]', 'Stokes Q [Jy]', 'Stokes U [Jy]', 'Stokes V [Jy]'],
                             ['I', 'Q', 'U', 'V']):
        a.axhline(0.0, color='gray', linestyle='--', alpha=0.7)
        a.set_ylabel(ylabel)
        a.set_title(f'{title} - {t}')
        a.grid(True, alpha=0.3)

    if n_scans > 1:
        ax[0].legend(fontsize=8, ncol=min(3, n_scans))

    ax[3].set_xlabel('Frequency [GHz]')
    plt.savefig(os.path.join(output_dir, filename), dpi=150)
    plt.close(fig)


def get_spinifex_rm_per_scan(pacal_name, scan_numbers):
    """
    Run Spinifex to estimate ionospheric RM and parse the output per scan.

    If the output file already exists it is parsed directly; otherwise
    Spinifex is launched via subprocess before parsing.

    Inputs:
        pacal_name  : str   - Name of the polarization angle calibrator field
        scan_numbers: array - MS scan numbers

    Outputs:
        iono_rm     : array - Median ionospheric RM per scan (rad/m²)
        iono_rm_err : array - Median RM uncertainty per scan (rad/m²)
    """
    spinifex_script = os.path.join(OXKAT, 'RMSYNTH_03_run_SPINIFEX.py')
    output_file = os.path.join(RESULTS,
                               f"{os.path.basename(myms)}_spinifex_rm.txt")

    iono_rm = np.zeros(len(scan_numbers), dtype=float)
    iono_rm_err = np.zeros(len(scan_numbers), dtype=float)

    # Run Spinifex if output not already present
    if not os.path.isfile(output_file):
        print(f"  Spinifex output not found — running: {spinifex_script}")
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        result = subprocess.run(
            ['python-spinifex', spinifex_script],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            print(f"  WARNING: Spinifex subprocess failed (returncode={result.returncode})")
            print(f"  stderr: {result.stderr[:500]}")
            print("  Ionospheric RM will be treated as 0 rad/m² for all scans.")
            return iono_rm, iono_rm_err
        print("  Spinifex completed successfully.")
    else:
        print(f"  Using existing Spinifex output: {output_file}")

    # Parse output file
    # Columns: Field_Name Position_hmsdms Scan_Number Time_ISOT Time_MJD RM RM_err
    scan_rm_data = {int(sc): [] for sc in scan_numbers}
    scan_rm_err_data = {int(sc): [] for sc in scan_numbers}
    n_field_rows = 0

    try:
        with open(output_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) < 8:
                    continue
                field = parts[0]
                if field != pacal_name:
                    continue
                n_field_rows += 1
                scan_num = int(parts[3])
                rm_val = float(parts[6])
                rm_err_val = float(parts[7])
                if scan_num in scan_rm_data:
                    scan_rm_data[scan_num].append(rm_val)
                    scan_rm_err_data[scan_num].append(rm_err_val)

        if n_field_rows == 0:
            print(f"  WARNING: Field '{pacal_name}' not found in Spinifex output.")
            print("  Ionospheric RM will be treated as 0 rad/m² for all scans.")
            return iono_rm, iono_rm_err

        for scan_idx, scan in enumerate(scan_numbers):
            sc = int(scan)
            if scan_rm_data[sc]:
                iono_rm[scan_idx] = np.median(scan_rm_data[sc])
                iono_rm_err[scan_idx] = np.median(scan_rm_err_data[sc])
                print(f"  Scan {scan}: ionospheric RM = "
                      f"{iono_rm[scan_idx]:+.3f} ± {iono_rm_err[scan_idx]:.3f} rad/m²")
            else:
                print(f"  Scan {scan}: not found in Spinifex output — using 0 rad/m²")

    except Exception as e:
        print(f"  WARNING: Failed to parse Spinifex output: {e}")
        print("  Ionospheric RM will be treated as 0 rad/m² for all scans.")

    return iono_rm, iono_rm_err


def resolve_pi_degeneracy_scan(median_gains, Q_scan, U_scan, V_scan,
                                chi_deg, freq_ghz, target_polang_array,
                                iono_rm, iono_rm_err=0.0, rm_range=2.0):
    """
    Resolve the ±π degeneracy in CASA polcal Xf gains for a single scan.

    For each trial RM in the grid:
        - Per channel, pick the sign (+1 or -1) that minimises |EVPA − target|
        - Compute robust stats across all valid channels: MAD residual and
          number of channels preferring each sign
    The trial RM with the lowest MAD residual is adopted globally.
    Per-channel signs are then re-evaluated at that single adopted RM.

    Grid step = iono_rm_err / 3 (minimum 0.05, fallback 0.2 if err is zero).

    Inputs:
        median_gains    : complex array (nchan,)
        Q_scan, U_scan, V_scan : float arrays (nchan,)
        chi_deg         : float — parallactic angle (degrees)
        freq_ghz        : float array (nchan,)
        target_polang_array : float array (nchan,)
        iono_rm         : float — Spinifex ionospheric RM (rad/m²)
        iono_rm_err     : float — Spinifex RM uncertainty (rad/m²)
        rm_range        : float — half-width of search grid (rad/m²)

    Outputs:
        signed_gains : complex array (nchan,) — per-channel sign applied
        best_rm      : float — adopted RM
        grid_centre  : float — XF_TARGET_RM + iono_rm
        best_delta_rm: float — grid offset from centre
        med_dev      : float — MAD residual at best RM
    """
    c = 2.998e8

    grid_centre = XF_TARGET_RM + iono_rm
    rm_step = max(0.05, iono_rm_err / 3.0) if iono_rm_err > 0 else 0.2
    rm_trials = np.arange(grid_centre - rm_range,
                          grid_centre + rm_range + rm_step / 2,
                          rm_step)

    signed_gains = np.full(len(median_gains), np.nan, dtype=complex)

    valid = (np.isfinite(median_gains) &
             np.isfinite(Q_scan) & np.isfinite(U_scan) & np.isfinite(V_scan))

    if np.sum(valid) < 5:
        print(f"    WARNING: Only {np.sum(valid)} valid channels — skipping grid search")
        signed_gains[valid] = median_gains[valid]
        return signed_gains, grid_centre, grid_centre, 0.0, np.inf

    freq_hz   = freq_ghz[valid] * 1e9
    lambda_sq = (c / freq_hz) ** 2
    Q_v  = Q_scan[valid]
    U_v  = U_scan[valid]
    V_v  = V_scan[valid]
    g_v  = median_gains[valid]
    target_v = target_polang_array[valid]
    n_valid  = int(np.sum(valid))

    # Pre-compute corrected Q_sky, U_sky for both signs
    sign_results = {}
    for sign in [+1, -1]:
        rho = np.angle(sign * g_v)
        U_xh, V_xh = correct_crosshand_phase(U_v, V_v, rho)
        Q_sky, U_sky = correct_parallactic_angle(Q_v, U_xh, chi_deg)
        sign_results[sign] = (Q_sky, U_sky)

    # Grid search
    print(f"\n    Grid search: {len(rm_trials)} trials, "
          f"step={rm_step:.3f} rad/m², "
          f"range=[{rm_trials[0]:.3f}, {rm_trials[-1]:.3f}] rad/m²")
    print(f"    {'Trial RM':>12}  {'ΔRM':>8}  {'MAD resid(°)':>14}  "
          f"{'N(+sign)':>10}  {'N(-sign)':>10}")
    print(f"    {'-'*60}")

    best_mad   = np.inf
    best_rm    = grid_centre

    for rm in rm_trials:
        delta_rm = rm - grid_centre
        dev_pos = np.abs(calculate_derotated_angle(
            sign_results[+1][0], sign_results[+1][1], rm, lambda_sq) - target_v)
        dev_neg = np.abs(calculate_derotated_angle(
            sign_results[-1][0], sign_results[-1][1], rm, lambda_sq) - target_v)
        dev_pos = np.where(dev_pos > 90, 180 - dev_pos, dev_pos)
        dev_neg = np.where(dev_neg > 90, 180 - dev_neg, dev_neg)

        # Per channel: pick best sign
        best_dev_per_ch = np.minimum(dev_pos, dev_neg)
        n_pos = int(np.sum(dev_pos <= dev_neg))
        n_neg = n_valid - n_pos

        # Robust stat: MAD of per-channel minimum residuals
        med  = float(np.median(best_dev_per_ch))
        mad  = float(np.median(np.abs(best_dev_per_ch - med)))

        marker = ' <<<' if mad < best_mad else ''
        print(f"    {rm:>12.4f}  {delta_rm:>+8.4f}  {mad:>14.4f}  "
              f"{n_pos:>10d}  {n_neg:>10d}{marker}")

        if mad < best_mad:
            best_mad = mad
            best_rm  = rm

    best_delta_rm = best_rm - grid_centre

    # ── Fine grid: per-channel RM and sign ─────────────────────────────────
    # Each channel independently finds the (RM, sign) pair within ±0.2 of the
    # coarse best_rm that minimises |EVPA − target| for that channel.
    fine_step   = 0.04
    fine_range  = 0.25
    fine_trials = np.arange(best_rm - fine_range,
                            best_rm + fine_range + fine_step / 2,
                            fine_step)

    print(f"\n    Fine grid (per-channel): {len(fine_trials)} trials, "
          f"step={fine_step:.3f} rad/m², "
          f"range=[{fine_trials[0]:.4f}, {fine_trials[-1]:.4f}] rad/m²")

    # Pre-compute corrected skies for both signs (same as coarse, reuse sign_results)
    per_ch_rm   = np.full(n_valid, best_rm)
    per_ch_sign = np.ones(n_valid, dtype=int)
    per_ch_dev  = np.full(n_valid, np.inf)

    for rm in fine_trials:
        for sign in [+1, -1]:
            angles = calculate_derotated_angle(
                sign_results[sign][0], sign_results[sign][1], rm, lambda_sq)
            dev = np.abs(angles - target_v)
            dev = np.where(dev > 90, 180 - dev, dev)
            improved = dev < per_ch_dev
            per_ch_rm[improved]   = rm
            per_ch_sign[improved] = sign
            per_ch_dev[improved]  = dev[improved]

    coarse_rm     = best_rm
    best_rm       = float(np.median(per_ch_rm))
    best_delta_rm = best_rm - grid_centre
    med_dev       = float(np.median(per_ch_dev))

    n_pos = int(np.sum(per_ch_sign == +1))
    n_neg = n_valid - n_pos
    print(f"    Per-channel RM range: [{per_ch_rm.min():.4f}, {per_ch_rm.max():.4f}] rad/m²")
    print(f"    Median per-channel RM: {best_rm:.4f} rad/m²  (MAD dev = {med_dev:.4f}°)")
    print(f"    Per-channel sign: {n_pos} channels +1, {n_neg} channels −1")

    # Compare fine-grid per-channel sign to what the coarse RM alone would give
    dev_pos_coarse = np.abs(calculate_derotated_angle(
        sign_results[+1][0], sign_results[+1][1], coarse_rm, lambda_sq) - target_v)
    dev_neg_coarse = np.abs(calculate_derotated_angle(
        sign_results[-1][0], sign_results[-1][1], coarse_rm, lambda_sq) - target_v)
    dev_pos_coarse = np.where(dev_pos_coarse > 90, 180 - dev_pos_coarse, dev_pos_coarse)
    dev_neg_coarse = np.where(dev_neg_coarse > 90, 180 - dev_neg_coarse, dev_neg_coarse)
    coarse_sign_per_ch = np.where(dev_pos_coarse <= dev_neg_coarse, +1, -1)

    n_swapped   = int(np.sum(per_ch_sign != coarse_sign_per_ch))
    n_unchanged = n_valid - n_swapped
    print(f"\n    Fine grid sign changes vs coarse RM sign:")
    print(f"      Unchanged : {n_unchanged:5d} / {n_valid} channels "
          f"({100.*n_unchanged/n_valid:.1f}%)")
    print(f"      Swapped   : {n_swapped:5d} / {n_valid} channels "
          f"({100.*n_swapped/n_valid:.1f}%)")

    print(f"\n    ── Adopted RM breakdown ──────────────────────────────")
    print(f"    Target RM (config):    XF_TARGET_RM       = {XF_TARGET_RM:+.4f} rad/m²")
    print(f"    Ionospheric RM:        iono_rm            = {iono_rm:+.4f} rad/m²")
    print(f"    Grid offset (ΔRM):     median(fine)−centre= {best_delta_rm:+.4f} rad/m²")
    print(f"    ──────────────────────────────────────────────────────")
    print(f"    Adopted RM (median of per-channel fine):   {best_rm:+.4f} rad/m²"
          f"  (MAD dev = {med_dev:.4f}°)")

    # Apply per-channel sign from fine grid and map RM back to full channel array
    valid_idx = np.where(valid)[0]
    for i, ch in enumerate(valid_idx):
        signed_gains[ch] = per_ch_sign[i] * median_gains[ch]

    per_ch_rm_full = np.full(len(median_gains), np.nan)
    for i, ch in enumerate(valid_idx):
        per_ch_rm_full[ch] = per_ch_rm[i]

    return signed_gains, best_rm, grid_centre, best_delta_rm, med_dev, coarse_rm, per_ch_rm_full


def global_sigma_clip_gains(per_scan_gains, per_scan_flags, scan_numbers, sigma=10.0):
    """
    Global σ-clip on XF phase, using circular statistics to handle the
    -180°/+180° wrap correctly.

    For each channel, the gain lies on the unit circle. We compute the
    circular mean phase across all unflagged channels, then measure each
    channel's angular deviation from that mean (wrapped to (-180, 180]).
    Channels where |deviation| > sigma × 1.4826 × MAD(deviations) are flagged.

    Inputs:
        per_scan_gains : dict {scan: complex array (nchan,)}
        per_scan_flags : dict {scan: bool array (nchan,)}
        scan_numbers   : array
        sigma          : float — clipping threshold (default 10σ)

    Outputs:
        per_scan_flags       : dict — updated flag arrays
        global_flagged_per_scan : dict {scan: bool array} — channels newly
                                  flagged by this step (not pre-existing flags)
    """
    # Pool all unflagged phases
    all_phases_deg = []
    for scan in scan_numbers:
        g = per_scan_gains[scan]
        f = per_scan_flags[scan]
        good = ~f & np.isfinite(g)
        all_phases_deg.extend(np.degrees(np.angle(g[good])))

    all_phases_deg = np.array(all_phases_deg)

    global_flagged_per_scan = {scan: np.zeros(len(per_scan_flags[scan]), dtype=bool)
                               for scan in scan_numbers}

    if len(all_phases_deg) < 10:
        print("  WARNING: Too few unflagged gains for global clip — skipping")
        return per_scan_flags, global_flagged_per_scan

    # Circular mean via unit-vector averaging
    mean_complex = np.mean(np.exp(1j * np.radians(all_phases_deg)))
    circ_mean_deg = np.degrees(np.angle(mean_complex))

    # Angular deviations, wrapped to (-180, 180]
    deviations = all_phases_deg - circ_mean_deg
    deviations = (deviations + 180) % 360 - 180

    mad_dev = np.median(np.abs(deviations - np.median(deviations)))
    sig_dev = 1.4826 * mad_dev
    threshold = sigma * sig_dev

    print(f"  Circular mean phase: {circ_mean_deg:+.3f}°, "
          f"MAD={mad_dev:.3f}°, threshold={threshold:.3f}° ({sigma:.0f}σ)")

    total_new = 0
    for scan in scan_numbers:
        g = per_scan_gains[scan]
        f = per_scan_flags[scan].copy()
        good = ~f & np.isfinite(g)

        phases = np.degrees(np.angle(g))
        dev = (phases - circ_mean_deg + 180) % 360 - 180
        clip = np.abs(dev) > threshold

        new_flags = good & clip
        n_new = int(np.sum(new_flags))
        total_new += n_new
        per_scan_flags[scan] = f | new_flags
        global_flagged_per_scan[scan] = new_flags

        if n_new > 0:
            print(f"    Scan {scan}: {n_new} channels flagged by {sigma:.0f}σ global clip")
        else:
            print(f"    Scan {scan}: 0 channels flagged by global clip")

    print(f"  Total newly flagged by global clip: {total_new}")
    return per_scan_flags, global_flagged_per_scan


def _bin_qu_native_to_xf(qu_native, n_xf):
    """Bin a native-resolution QU array to n_xf channels by nanmean."""
    nchan    = len(qu_native)
    bin_size = nchan / n_xf
    qu_xf    = np.full(n_xf, np.nan)
    for i in range(n_xf):
        lo    = int(round(i * bin_size))
        hi    = int(round((i + 1) * bin_size))
        valid = qu_native[lo:hi]
        valid = valid[np.isfinite(valid)]
        if len(valid):
            qu_xf[i] = np.nanmean(valid)
    return qu_xf


def poly_baseline_mad_clip(gains, flags, n_chan, poly_order, window, sigma,
                           iterative=True, max_iter=20, weights=None):
    """
    Polynomial-baseline sliding circular MAD clip on XF gains.

    Completely wrap-safe — all operations in the phasor/complex domain:

    Step 1 — Polynomial baseline fit:
        Fit poly_order polynomial to cos(phase) and sin(phase) vs normalised
        channel index on [−1, +1]. This is wrap-safe because cos/sin are
        continuous everywhere. Subtract the baseline phasor from each channel's
        gain to give baseline-subtracted residual phasors.
        If weights are provided (log1p-compressed cross-hand flux per channel),
        the polyfit is weighted so high-S/N channels anchor the baseline.

    Step 2 — Sliding circular MAD:
        For each unflagged channel, gather the nearest unflagged neighbours
        (by index). Compute their circular median (minimises sum of |angle|
        differences). Compute circular MAD = median(|angle(sig × conj(med))|).
        Flag the channel if |angle(sig_ch × conj(med))| > sigma × 1.4826 × circ_MAD.

    If iterative=True, steps 1 and 2 repeat until no new channels are flagged
    or max_iter is reached. Each iteration refits the polynomial on surviving
    channels, giving a cleaner baseline for the next MAD pass.

    Inputs:
        gains      : complex array (nchan,) — unit-amplitude gains
        flags      : bool array (nchan,)    — True = already flagged
        n_chan     : int
        poly_order : int — polynomial degree (use XF_POLY_ORDER)
        window     : int — sliding window size in channels
        sigma      : float — MAD threshold (use XF_MAD_SIGMA_CLIP)
        iterative  : bool — refit and re-clip until convergence (use XF_MAD_ITERATIVE)
        max_iter   : int — maximum iterations
        weights    : float array (nchan,) or None — per-channel polyfit weights,
                     e.g. log1p(sqrt(Q²+U²) + 1e-8) on xftab frequency grid.
                     None = uniform weighting.

    Outputs:
        flags     : bool array — updated
        p_cos     : polynomial coefficients for cos(phase) baseline (or None)
        p_sin     : polynomial coefficients for sin(phase) baseline (or None)
        n_flagged : int — total channels newly flagged across all iterations
    """
    flags     = flags.copy()
    chan_idx  = np.arange(n_chan, dtype=float)
    chan_norm = 2.0 * chan_idx / max(n_chan - 1, 1) - 1.0
    n_flagged_total = 0
    p_cos = None
    p_sin = None

    # log1p-compressed weights: sqrt(Q²+U²) is positive-definite so log1p is
    # safe with a 1e-8 floor. No normalisation needed — polyfit only cares
    # about relative weights, and log1p naturally compresses dynamic range
    # without setting an artificial floor.
    if weights is not None:
        w_norm = np.log1p(np.maximum(weights, 1e-8))
        w_norm = np.where(np.isfinite(w_norm), w_norm, 0.0)
    else:
        w_norm = np.ones(n_chan)

    for iteration in range(max_iter):
        # ── Step 1: polynomial baseline ─────────────────────────────────────
        good   = ~flags & np.isfinite(gains)
        n_good = int(np.sum(good))

        if n_good > poly_order + 1:
            ph       = np.angle(gains[good])
            c_norm   = chan_norm[good]
            w_good   = w_norm[good]
            try:
                p_cos = np.polyfit(c_norm, np.cos(ph), poly_order, w=w_good)
                p_sin = np.polyfit(c_norm, np.sin(ph), poly_order, w=w_good)
            except (np.linalg.LinAlgError, ValueError):
                p_cos = None
                p_sin = None

        if p_cos is not None:
            baseline       = np.polyval(p_cos, chan_norm) + 1j * np.polyval(p_sin, chan_norm)
            baseline_phase = np.angle(baseline)
            signal         = gains * np.exp(-1j * baseline_phase)
        else:
            signal = gains.copy()

        # ── Step 2: sliding circular MAD ────────────────────────────────────
        good_idx = np.where(~flags & np.isfinite(signal))[0]
        n_valid  = len(good_idx)

        chan_dev_sigma = np.full(n_chan, np.nan)
        chan_off_deg   = np.full(n_chan, np.nan)
        chan_mad_deg   = np.full(n_chan, np.nan)
        chan_med_deg   = np.full(n_chan, np.nan)
        chan_win_n     = np.zeros(n_chan, dtype=int)
        new_flag       = np.zeros(n_chan, dtype=bool)

        if n_valid >= 2:
            for pos, ch in enumerate(good_idx):
                others = good_idx[good_idx != ch]
                if len(others) < 2:
                    continue
                dists   = np.abs(others - ch)
                win_idx = others[np.argsort(dists)[:window]]
                win_sig = signal[win_idx]
                chan_win_n[ch] = len(win_sig)

                win_angles  = np.angle(win_sig)
                best_sum    = np.inf
                circ_med_ph = win_angles[0]
                for _cand in win_angles:
                    _s = float(np.sum(np.abs(np.angle(np.exp(1j * (win_angles - _cand))))))
                    if _s < best_sum:
                        best_sum    = _s
                        circ_med_ph = _cand
                circ_med = np.exp(1j * circ_med_ph)

                circ_mad = float(np.median(np.abs(np.angle(win_sig * np.conj(circ_med)))))
                if circ_mad < 1e-10:
                    continue
                circ_stdev = 1.4826 * circ_mad

                off_ch    = np.angle(signal[ch] * np.conj(circ_med))
                dev_sigma = abs(off_ch) / circ_stdev

                chan_dev_sigma[ch] = dev_sigma
                chan_off_deg[ch]   = np.degrees(off_ch)
                chan_mad_deg[ch]   = np.degrees(circ_mad)
                chan_med_deg[ch]   = np.degrees(circ_med_ph)

                if dev_sigma > sigma:
                    new_flag[ch] = True

        n_new = int(np.sum(new_flag))
        n_flagged_total += n_new

        iter_label = f"iteration {iteration + 1}" if iterative else "single pass"
        print(f"\n    [{iter_label}] {n_new} channels flagged "
              f"({int(np.sum(~flags)) - n_new}/{n_chan} remaining after this pass)")

        # ── Logging ─────────────────────────────────────────────────────────
        has_no_data = ~(~flags & np.isfinite(gains))
        print(f"    {'ch':>5}  {'N':>4}  {'median(deg)':>11}  {'MAD(deg)':>9}  "
              f"{'offset(deg)':>11}  {'dev(sigma)':>10}  {'thr':>6}  {'flag':>6}")
        print(f"    {'-'*72}")
        for ch in range(n_chan):
            if has_no_data[ch]:
                print(f"    {ch:>5}  {'---':>4}  {'pre-flagged':>11}  "
                      f"{'---':>9}  {'---':>11}  {'---':>10}  {sigma:>6.2f}  {'(pre)':>6}")
            elif np.isnan(chan_dev_sigma[ch]):
                print(f"    {ch:>5}  {chan_win_n[ch]:>4}  {'(no window)':>11}  "
                      f"{'---':>9}  {'---':>11}  {'---':>10}  {sigma:>6.2f}  {'---':>6}")
            else:
                flag_str = 'CLIP' if new_flag[ch] else 'ok'
                print(f"    {ch:>5}  {chan_win_n[ch]:>4}  "
                      f"{chan_med_deg[ch]:>+11.3f}  "
                      f"{chan_mad_deg[ch]:>9.3f}  "
                      f"{chan_off_deg[ch]:>+11.3f}  "
                      f"{chan_dev_sigma[ch]:>10.3f}  "
                      f"{sigma:>6.2f}  "
                      f"{flag_str:>6}")

        flags = flags | new_flag

        if n_new == 0 or not iterative:
            break

    if iterative:
        print(f"\n    Converged after {iteration + 1} iteration(s), "
              f"{n_flagged_total} total channels flagged")

    return flags, p_cos, p_sin, n_flagged_total


def plot_xf_phase_diagnostic(freq_ghz, chan_freq_xf_ghz,
                              pre_pi_gains_per_scan, resolved_gains_per_scan,
                              flagged_gains_per_scan, poly_coeffs_per_scan,
                              final_gains_per_scan, snr_per_scan,
                              scan_numbers, output_dir):
    """
    Five-panel diagnostic plot of XF gains through the full processing pipeline.

    Panel 1: Raw CASA polcal output (before ±π resolution)
    Panel 2: After ±π resolution — flagged channels highlighted
    Panel 3: After MAD clipping — surviving channels, flagged markers, poly baseline
    Panel 4: Final averaged XF solution (xftab)
    Panel 5: Per-channel cross-hand SNR with XF_MIN_SN threshold

    Inputs:
        freq_ghz               : float array (nchan,)
        chan_freq_xf_ghz       : float array (n_xf,)
        pre_pi_gains_per_scan  : dict {scan: complex (nchan,)}
        resolved_gains_per_scan: dict {scan: complex (nchan,)}
        flagged_gains_per_scan : dict {scan: 5-tuple (gains, flags, global_fl, poly_fl, snr_fl)}
        poly_coeffs_per_scan   : dict {scan: (p_cos, p_sin) or (None, None)}
        final_gains_per_scan   : dict {scan: complex (n_xf,)}
        snr_per_scan           : dict {scan: float array (n_xf,)}
        scan_numbers           : array
        output_dir             : str
    """
    n_scans = len(scan_numbers)
    colors = ['C0'] if n_scans == 1 else plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))

    fig, axes = plt.subplots(5, 1, figsize=(14, 20), sharex=True,
                             constrained_layout=True)

    p3_unflagged_phases = []
    p4_phases = []

    for scan_idx, scan in enumerate(scan_numbers):
        color = colors[scan_idx % len(colors)]
        label = f'Scan {scan}' if n_scans > 1 else None

        # --- Panel 1: raw CASA polcal ---
        pre_g = pre_pi_gains_per_scan[scan]
        valid_pre = np.isfinite(pre_g)
        axes[0].plot(freq_ghz[valid_pre], np.degrees(np.angle(pre_g[valid_pre])),
                     '.', color=color, alpha=0.85, markersize=6, label=label)

        # --- Panel 2: after ±π resolution, flagged channels highlighted ---
        res_g = resolved_gains_per_scan[scan]
        valid_res = np.isfinite(res_g)
        gains_fl, flags_fl, global_fl, poly_fl, snr_fl = flagged_gains_per_scan[scan]
        will_be_flagged = (global_fl | poly_fl | snr_fl) & valid_res
        stays_good = valid_res & ~will_be_flagged

        axes[1].plot(freq_ghz[stays_good], np.degrees(np.angle(res_g[stays_good])),
                     '.', color=color, alpha=0.85, markersize=6, label=label)
        if np.any(will_be_flagged & global_fl):
            axes[1].plot(freq_ghz[will_be_flagged & global_fl],
                         np.degrees(np.angle(res_g[will_be_flagged & global_fl])),
                         'D', color='darkorange', alpha=0.8, markersize=6,
                         label='Will be global-clipped' if scan_idx == 0 else None)
        if np.any(will_be_flagged & poly_fl):
            axes[1].plot(freq_ghz[will_be_flagged & poly_fl],
                         np.degrees(np.angle(res_g[will_be_flagged & poly_fl])),
                         'x', color='red', alpha=0.8, markersize=6, linewidth=0.8,
                         label='Will be MAD-clipped' if scan_idx == 0 else None)
        if np.any(will_be_flagged & snr_fl):
            axes[1].plot(freq_ghz[will_be_flagged & snr_fl],
                         np.degrees(np.angle(res_g[will_be_flagged & snr_fl])),
                         's', color='purple', alpha=0.8, markersize=6,
                         label='Will be SNR-flagged' if scan_idx == 0 else None)

        # --- Panel 3: after clipping ---
        good = ~flags_fl & np.isfinite(gains_fl)
        good_phases = np.degrees(np.angle(gains_fl[good]))
        axes[2].plot(freq_ghz[good], good_phases,
                     '.', color=color, alpha=0.85, markersize=6, label=label)
        p3_unflagged_phases.extend(good_phases)

        global_bad = global_fl & np.isfinite(gains_fl)
        if np.any(global_bad):
            axes[2].plot(freq_ghz[global_bad], np.degrees(np.angle(gains_fl[global_bad])),
                         'D', color='darkorange', alpha=0.8, markersize=6,
                         label=f'Global clip ({int(np.sum(global_bad))} ch)'
                               if scan_idx == 0 else None)
        poly_bad = poly_fl & np.isfinite(gains_fl)
        if np.any(poly_bad):
            axes[2].plot(freq_ghz[poly_bad], np.degrees(np.angle(gains_fl[poly_bad])),
                         'x', color='red', alpha=0.8, markersize=6,
                         label=f'MAD clip ({int(np.sum(poly_bad))} ch)'
                               if scan_idx == 0 else None)
        snr_bad = snr_fl & np.isfinite(gains_fl)
        if np.any(snr_bad):
            axes[2].plot(freq_ghz[snr_bad], np.degrees(np.angle(gains_fl[snr_bad])),
                         's', color='purple', alpha=0.8, markersize=6,
                         label=f'SNR flag ({int(np.sum(snr_bad))} ch)'
                               if scan_idx == 0 else None)

        # Polynomial baseline overlay
        p_cos, p_sin = poly_coeffs_per_scan[scan]
        if p_cos is not None and p_sin is not None:
            n_plt    = len(freq_ghz)
            c_norm_plt = 2.0 * np.arange(n_plt) / max(n_plt - 1, 1) - 1.0
            baseline_complex = (np.polyval(p_cos, c_norm_plt) +
                                1j * np.polyval(p_sin, c_norm_plt))
            axes[2].plot(freq_ghz, np.degrees(np.angle(baseline_complex)),
                         '-', color='black', linewidth=2.0, zorder=100,
                         label='Poly baseline' if scan_idx == 0 else None)

        # --- Panel 4: final averaged solution ---
        final_g = final_gains_per_scan[scan]
        valid_final = np.isfinite(final_g)
        final_phases = np.degrees(np.angle(final_g[valid_final]))
        axes[3].plot(chan_freq_xf_ghz[valid_final], final_phases,
                     '-o', color=color, alpha=0.8, markersize=6,
                     linewidth=1, label=label)
        p4_phases.extend(final_phases)

        # --- Panel 5: per-channel cross-hand SNR ---
        snr = snr_per_scan[scan]
        valid_snr = np.isfinite(snr)
        axes[4].plot(chan_freq_xf_ghz[valid_snr], snr[valid_snr],
                     '-o', color=color, alpha=0.8, markersize=6,
                     linewidth=1, label=label)
        # Mark SNR-flagged channels
        if np.any(snr_fl & valid_snr):
            axes[4].plot(chan_freq_xf_ghz[snr_fl & valid_snr],
                         snr[snr_fl & valid_snr],
                         's', color='purple', alpha=0.8, markersize=6,
                         label=f'SNR flagged' if scan_idx == 0 else None)

    # SNR threshold line
    if XF_MIN_SN > 0:
        axes[4].axhline(XF_MIN_SN, color='red', linestyle='--', linewidth=1.5,
                        label=f'Threshold (SNR={XF_MIN_SN:.1f})')

    titles = [
        'Raw CASA polcal Xf (before ±π resolution)',
        'After ±π resolution (◆ = global, × = MAD, ■ = SNR flagged)',
        'After clipping (◆ = global, × = MAD, ■ = SNR, black = poly baseline)',
        'Final averaged XF solution (xftab)',
        f'Cross-hand SNR per channel (threshold = {XF_MIN_SN:.1f})',
    ]
    ylabels = ['XF Phase [deg]', 'XF Phase [deg]', 'XF Phase [deg]',
               'XF Phase [deg]', 'SNR']
    for ax, title, ylabel in zip(axes, titles, ylabels):
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

    for ax in axes[:4]:
        ax.set_ylim(-190, 190)

    # Panel 3 zoom
    if len(p3_unflagged_phases) > 0:
        p3_min = np.nanmin(p3_unflagged_phases)
        p3_max = np.nanmax(p3_unflagged_phases)
        p3_pad = (p3_max - p3_min) * 0.10
        axes[2].set_ylim(p3_min - p3_pad, p3_max + p3_pad)

    # Panel 4 zoom
    if len(p4_phases) > 0:
        p4_min = np.nanmin(p4_phases)
        p4_max = np.nanmax(p4_phases)
        p4_pad = max((p4_max - p4_min) * 0.10, 2.0)
        axes[3].set_ylim(p4_min - p4_pad, p4_max + p4_pad)

    axes[4].set_xlabel('Frequency [GHz]')
    axes[4].set_xlim(freq_ghz.min(), freq_ghz.max())
    axes[4].set_ylim(bottom=0)

    # Legends
    if n_scans > 1:
        for ax in (axes[0], axes[1], axes[3], axes[4]):
            ax.legend(fontsize=8, ncol=min(3, n_scans))
    axes[2].legend(fontsize=8, ncol=min(3, n_scans + 2))
    axes[4].legend(fontsize=8)

    plt.savefig(os.path.join(output_dir, 'xf_phase_diagnostic.png'), dpi=150)
    plt.close(fig)


def plot_derotated_angle_residuals(freq_ghz, Q_scans, U_scans, V_scans,
                                   resolved_gains_per_scan, raw_flags_per_scan,
                                   chi_deg_scans, best_rm_per_scan,
                                   target_polang_array, scan_numbers, output_dir,
                                   filename='derotated_angle_residuals.png',
                                   flag_override=None, per_ch_rm_per_scan=None):
    """
    Three-panel diagnostic plot of de-rotated polarization angle after ±π resolution.

    Panel 1: De-rotated EVPA vs frequency, overlaid with target.
    Panel 2: Residual deviation (EVPA − target) vs frequency, wrapped to (−90, 90].
    Panel 3: Histogram of EVPA residuals centred on zero.

    flag_override      : dict {scan: bool array} — post-clipping flags for stage 3.
    per_ch_rm_per_scan : dict {scan: float array (nchan,)} — per-channel fine grid RM.
                         If supplied, Faraday de-rotation uses the per-channel value
                         rather than the single best_rm. Falls back to best_rm_per_scan
                         for channels where per_ch_rm is NaN.
    """
    c = 2.998e8
    n_scans = len(scan_numbers)
    colors = ['C0'] if n_scans == 1 else plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))

    fig, axes = plt.subplots(3, 1, figsize=(14, 14), constrained_layout=True)

    all_evpa = []
    all_scan_labels = []

    for scan_idx, scan in enumerate(scan_numbers):
        color = colors[scan_idx % len(colors)]
        label = f'Scan {scan}' if n_scans > 1 else None

        g = resolved_gains_per_scan[scan]
        flags = flag_override[scan] if flag_override is not None else raw_flags_per_scan[scan]
        valid = (~flags & np.isfinite(g) &
                 np.isfinite(Q_scans[scan_idx]) & np.isfinite(U_scans[scan_idx]) &
                 np.isfinite(V_scans[scan_idx]))

        if np.sum(valid) < 5:
            print(f"  WARNING: Scan {scan} has fewer than 5 valid channels for residuals plot")
            continue

        freq_hz = freq_ghz[valid] * 1e9
        lambda_sq = (c / freq_hz) ** 2

        Q_v = Q_scans[scan_idx][valid]
        U_v = U_scans[scan_idx][valid]
        V_v = V_scans[scan_idx][valid]
        g_v = g[valid]

        # Apply XF correction using gain phase directly
        rho = np.angle(g_v)
        U_xh, V_xh = correct_crosshand_phase(U_v, V_v, rho)

        # Parallactic angle correction (feed → sky)
        Q_sky, U_sky = correct_parallactic_angle(Q_v, U_xh, chi_deg_scans[scan_idx])

        # Faraday de-rotation: per-channel RM if available, else global best_rm
        if per_ch_rm_per_scan is not None:
            rm_ch = per_ch_rm_per_scan[scan][valid]
            fallback = ~np.isfinite(rm_ch)
            rm_ch[fallback] = best_rm_per_scan[scan]
            # Compute per-channel EVPA using per-channel RM
            evpa_deg = np.array([
                calculate_derotated_angle(
                    Q_sky[i:i+1], U_sky[i:i+1], rm_ch[i], lambda_sq[i:i+1]
                )[0]
                for i in range(len(rm_ch))
            ])
        else:
            evpa_deg = calculate_derotated_angle(Q_sky, U_sky, best_rm_per_scan[scan], lambda_sq)

        target_v = target_polang_array[valid]
        deviation = (evpa_deg - target_v + 90) % 180 - 90

        axes[0].plot(freq_ghz[valid], evpa_deg, '.', color=color,
                     alpha=0.85, markersize=6, label=label)
        axes[1].plot(freq_ghz[valid], deviation, '.', color=color,
                     alpha=0.85, markersize=6, label=label)

        all_evpa.extend(evpa_deg)
        all_scan_labels.extend([scan] * len(evpa_deg))

    axes[0].plot(freq_ghz, target_polang_array, 'k--', linewidth=1.5,
                 alpha=0.8, label='Target EVPA')
    axes[1].axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.8)

    axes[0].set_ylabel('De-rotated EVPA [deg]')
    axes[0].set_title('De-rotated polarization angle after XF + parang + Faraday correction')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8, ncol=min(3, n_scans + 1))

    axes[1].set_ylabel('EVPA − Target [deg]')
    axes[1].set_xlabel('Frequency [GHz]')
    axes[1].set_title('Residual deviation from target EVPA (should be ≈ 0)')
    axes[1].grid(True, alpha=0.3)

    # Panel 3: histogram of residuals (EVPA − target), centred on zero
    all_evpa = np.array(all_evpa)
    all_scan_labels = np.array(all_scan_labels)

    # Recompute residuals for the histogram
    all_residuals = []
    all_residual_labels = []
    for scan_idx, scan in enumerate(scan_numbers):
        g = resolved_gains_per_scan[scan]
        flags = flag_override[scan] if flag_override is not None else raw_flags_per_scan[scan]
        valid = (~flags & np.isfinite(g) &
                 np.isfinite(Q_scans[scan_idx]) & np.isfinite(U_scans[scan_idx]) &
                 np.isfinite(V_scans[scan_idx]))
        if np.sum(valid) < 5:
            continue
        freq_hz = freq_ghz[valid] * 1e9
        lambda_sq = (c / freq_hz) ** 2
        Q_v = Q_scans[scan_idx][valid]
        U_v = U_scans[scan_idx][valid]
        V_v = V_scans[scan_idx][valid]
        g_v = g[valid]
        rho = np.angle(g_v)
        U_xh, V_xh = correct_crosshand_phase(U_v, V_v, rho)
        Q_sky, U_sky = correct_parallactic_angle(Q_v, U_xh, chi_deg_scans[scan_idx])
        if per_ch_rm_per_scan is not None:
            rm_ch = per_ch_rm_per_scan[scan][valid].copy()
            rm_ch[~np.isfinite(rm_ch)] = best_rm_per_scan[scan]
            evpa = np.array([
                calculate_derotated_angle(
                    Q_sky[i:i+1], U_sky[i:i+1], rm_ch[i], lambda_sq[i:i+1]
                )[0]
                for i in range(len(rm_ch))
            ])
        else:
            evpa = calculate_derotated_angle(Q_sky, U_sky, best_rm_per_scan[scan], lambda_sq)
        residual = (evpa - target_polang_array[valid] + 90) % 180 - 90
        all_residuals.extend(residual)
        all_residual_labels.extend([scan] * len(residual))

    all_residuals = np.array(all_residuals)
    all_residual_labels = np.array(all_residual_labels)

    if len(all_residuals) > 0:
        unique_scans = np.unique(all_residual_labels)
        for i, scan in enumerate(unique_scans):
            color = colors[i % len(colors)]
            mask = all_residual_labels == scan
            rm = best_rm_per_scan[scan]
            axes[2].hist(all_residuals[mask], bins=50, alpha=0.6,
                         label=f'Scan {scan} (RM={rm:.2f} rad/m²)',
                         color=color, edgecolor='black', linewidth=0.5)

        mean_res  = float(np.mean(all_residuals))
        median_res = float(np.median(all_residuals))
        std_res   = float(np.std(all_residuals))

        axes[2].axvline(0, color='black', linestyle='--', linewidth=2,
                        label='Zero (perfect correction)', zorder=10)
        axes[2].axvline(mean_res, color='green', linestyle=':', linewidth=2,
                        label=f'Mean = {mean_res:+.2f}°', alpha=0.8)

        textstr = (f'Mean:   {mean_res:+.2f}°\n'
                   f'Median: {median_res:+.2f}°\n'
                   f'Std:    {std_res:.2f}°\n'
                   f'N: {len(all_residuals)}')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        axes[2].text(0.02, 0.98, textstr, transform=axes[2].transAxes, fontsize=10,
                     verticalalignment='top', bbox=props)

        print(f"  EVPA residual statistics:")
        print(f"    Mean={mean_res:+.2f}°, Median={median_res:+.2f}°, Std={std_res:.2f}°")

    axes[2].set_xlabel('EVPA − Target [deg]')
    axes[2].set_ylabel('Count')
    axes[2].set_title('Histogram of EVPA residuals (should be centred on 0)')
    axes[2].legend(fontsize=8)
    axes[2].grid(True, alpha=0.3)

    plt.savefig(os.path.join(output_dir, filename), dpi=150)
    plt.close(fig)


def plot_raw_xf_phase_panels(freq_ghz, raw_gains_per_scan, raw_flags_per_scan,
                              pre_pi_gains_per_scan, best_sign_per_scan,
                              Q_scans, U_scans, V_scans, chi_deg_scans,
                              best_rm_per_scan, target_polang_array,
                              scan_numbers, output_dir):
    """
    Stage 1 three-panel diagnostic.

    Panel 1: Raw XF phase vs frequency (pre-resolution)
    Panel 2: Accepted per-channel residuals (EVPA_accepted − target) vs frequency
    Panel 3: Rejected per-channel residuals vs frequency, y-range locked to panel 2.
             Inset (top-right) shows full y-range so divergence is visible even
             when main panel is clipped. Title notes ideal case: blank main panel.
    """
    c = 2.998e8
    n_scans = len(scan_numbers)
    colors = ['C0'] if n_scans == 1 else plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))

    fig, axes = plt.subplots(3, 1, figsize=(14, 14), sharex=True,
                             constrained_layout=True)

    rej_all_freq = []
    rej_all_dev  = []
    rej_all_col  = []

    for scan_idx, scan in enumerate(scan_numbers):
        color = colors[scan_idx % len(colors)]
        label = f'Scan {scan}' if n_scans > 1 else None

        g_pre  = pre_pi_gains_per_scan[scan]   # raw, before per-channel sign
        g_acc  = raw_gains_per_scan[scan]       # after per-channel sign
        flags  = raw_flags_per_scan[scan]
        valid  = (~flags & np.isfinite(g_pre) &
                  np.isfinite(Q_scans[scan_idx]) & np.isfinite(U_scans[scan_idx]) &
                  np.isfinite(V_scans[scan_idx]))

        if np.sum(valid) < 5:
            continue

        # Panel 1: raw phase
        axes[0].plot(freq_ghz[valid], np.degrees(np.angle(g_pre[valid])),
                     '.', color=color, alpha=0.85, markersize=6, label=label)

        freq_hz   = freq_ghz[valid] * 1e9
        lambda_sq = (c / freq_hz) ** 2
        Q_v = Q_scans[scan_idx][valid]
        U_v = U_scans[scan_idx][valid]
        V_v = V_scans[scan_idx][valid]
        target_v = target_polang_array[valid]
        rm = best_rm_per_scan[scan]

        # Accepted: use per-channel signed gains
        g_acc_v = g_acc[valid]
        rho_acc = np.angle(g_acc_v)
        U_xh, V_xh = correct_crosshand_phase(U_v, V_v, rho_acc)
        Q_sky, U_sky = correct_parallactic_angle(Q_v, U_xh, chi_deg_scans[scan_idx])
        evpa_acc = calculate_derotated_angle(Q_sky, U_sky, rm, lambda_sq)
        dev_acc  = (evpa_acc - target_v + 90) % 180 - 90
        axes[1].plot(freq_ghz[valid], dev_acc, '.', color=color,
                     alpha=0.85, markersize=6, label=label)

        # Rejected: negate per-channel signed gains
        rho_rej = np.angle(-g_acc_v)
        U_xh, V_xh = correct_crosshand_phase(U_v, V_v, rho_rej)
        Q_sky, U_sky = correct_parallactic_angle(Q_v, U_xh, chi_deg_scans[scan_idx])
        evpa_rej = calculate_derotated_angle(Q_sky, U_sky, rm, lambda_sq)
        dev_rej  = (evpa_rej - target_v + 90) % 180 - 90
        axes[2].plot(freq_ghz[valid], dev_rej, '.', color=color,
                     alpha=0.85, markersize=6, label=label)

        rej_all_freq.extend(freq_ghz[valid])
        rej_all_dev.extend(dev_rej)
        rej_all_col.extend([color] * int(np.sum(valid)))

    axes[0].set_ylabel('XF Phase [deg]')
    axes[0].set_title('Stage 1: Raw CASA polcal Xf phase (before ±π resolution)')
    axes[0].grid(True, alpha=0.3)
    if n_scans > 1:
        axes[0].legend(fontsize=8, ncol=min(3, n_scans))

    axes[1].axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    axes[1].set_ylabel('EVPA − Target [deg]')
    axes[1].set_title('Stage 1: Accepted per-channel residuals (should cluster near 0)')
    axes[1].grid(True, alpha=0.3)

    # Lock rejected y-range to accepted; add inset showing full range
    acc_ylim = axes[1].get_ylim()
    axes[2].axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    axes[2].set_ylim(acc_ylim)
    axes[2].set_ylabel('EVPA − Target [deg]')
    axes[2].set_xlabel('Frequency [GHz]')
    axes[2].set_title(
        'Stage 1: Rejected per-channel residuals (same y-scale; ideally blank — '
        'divergence clipped to accepted range, see inset for full range)')
    axes[2].grid(True, alpha=0.3)

    # Inset: full y-range of rejected residuals (top-right of panel 3)
    if len(rej_all_dev) > 0:
        ax_inset = axes[2].inset_axes([0.62, 0.55, 0.36, 0.42])
        ax_inset.scatter(rej_all_freq, rej_all_dev,
                         c=rej_all_col, alpha=0.3, s=6)
        ax_inset.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.6)
        ax_inset.set_title('Full range', fontsize=7)
        ax_inset.tick_params(labelsize=6)
        ax_inset.grid(True, alpha=0.2)

    plt.savefig(os.path.join(output_dir, 'xf_phase_stage1_raw.png'), dpi=150)
    plt.close(fig)


def print_per_channel_resolution(freq_ghz, signed_gains, pre_pi_gains, flags,
                                  Q_scan, U_scan, V_scan, chi_deg,
                                  best_rm, target_polang_array):
    """
    Print per-channel resolution table after per-channel ±π decision is made.

    For every valid channel prints: frequency, expected EVPA, measured EVPA
    for the accepted (signed) and rejected (negated) solutions, and the delta.
    """
    c = 2.998e8
    valid = (~flags & np.isfinite(signed_gains) &
             np.isfinite(Q_scan) & np.isfinite(U_scan) & np.isfinite(V_scan))

    if np.sum(valid) < 1:
        print("    (no valid channels)")
        return

    freq_hz   = freq_ghz[valid] * 1e9
    lambda_sq = (c / freq_hz) ** 2
    Q_v  = Q_scan[valid]
    U_v  = U_scan[valid]
    V_v  = V_scan[valid]
    g_acc = signed_gains[valid]
    target_v = target_polang_array[valid]
    freq_v   = freq_ghz[valid]

    results = {}
    for label, g_use in [('acc', g_acc), ('rej', -g_acc)]:
        rho = np.angle(g_use)
        U_xh, V_xh = correct_crosshand_phase(U_v, V_v, rho)
        Q_sky, U_sky = correct_parallactic_angle(Q_v, U_xh, chi_deg)
        evpa = calculate_derotated_angle(Q_sky, U_sky, best_rm, lambda_sq)
        delta = (evpa - target_v + 90) % 180 - 90
        results[label] = (evpa, delta)

    evpa_acc, delta_acc = results['acc']
    evpa_rej, delta_rej = results['rej']

    # Per-channel XF phase (from signed gains) and sign (+1 or -1)
    g_pre_v = pre_pi_gains[valid]
    xf_phase_deg = np.degrees(np.angle(g_acc))
    # Sign = +1 if accepted gain matches pre-pi gain, -1 if negated
    sign_per_ch = np.where(np.real(g_acc * np.conj(g_pre_v)) >= 0, '+1', '-1')

    print(f"\n    {'Chan_freq(GHz)':<16} {'XF_Phase(°)':<13} {'Sign':<6} {'Target(°)':<12} "
          f"{'Acc_EVPA(°)':<14} {'Acc_Δ(°)':<12} "
          f"{'Rej_EVPA(°)':<14} {'Rej_Δ(°)':<10}")
    print(f"    {'-'*97}")
    for i in range(len(freq_v)):
        print(f"    {freq_v[i]:<16.4f} {xf_phase_deg[i]:<13.3f} {sign_per_ch[i]:<6} "
              f"{target_v[i]:<12.3f} "
              f"{evpa_acc[i]:<14.3f} {delta_acc[i]:<12.3f} "
              f"{evpa_rej[i]:<14.3f} {delta_rej[i]:<10.3f}")


def plot_per_channel_rm_histogram(freq_ghz, per_ch_rm_per_scan, coarse_rm_per_scan,
                                   raw_flags_per_scan, scan_numbers, output_dir):
    """
    Histogram of per-channel adopted RMs from the fine grid search.

    Shows the distribution of per-channel best RMs across all valid channels,
    with a vertical dashed line at the coarse grid RM for reference, and a
    dotted line at the per-channel median.
    """
    n_scans = len(scan_numbers)
    colors  = ['C0'] if n_scans == 1 else plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))

    fig, ax = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)

    for scan_idx, scan in enumerate(scan_numbers):
        color    = colors[scan_idx % len(colors)]
        flags    = raw_flags_per_scan[scan]
        per_ch   = per_ch_rm_per_scan[scan]
        valid    = ~flags & np.isfinite(per_ch)
        if not np.any(valid):
            continue

        rm_vals   = per_ch[valid]
        coarse_rm = coarse_rm_per_scan[scan]
        med_rm    = float(np.median(rm_vals))
        std_rm    = float(np.std(rm_vals))

        ax.hist(rm_vals, bins=40, alpha=0.6, color=color,
                edgecolor='black', linewidth=0.4,
                label=f'Scan {scan} (N={int(np.sum(valid))})'
                      if n_scans > 1 else f'Per-channel RM (N={int(np.sum(valid))})')
        ax.axvline(coarse_rm, color=color, linestyle='--', linewidth=2,
                   alpha=0.9,
                   label=f'Scan {scan} coarse RM = {coarse_rm:.4f} rad/m²'
                         if n_scans > 1 else f'Coarse RM = {coarse_rm:.4f} rad/m²')
        ax.axvline(med_rm, color=color, linestyle=':', linewidth=1.5, alpha=0.8,
                   label=f'Scan {scan} median = {med_rm:.4f} ± {std_rm:.4f} rad/m²'
                         if n_scans > 1 else f'Median = {med_rm:.4f} ± {std_rm:.4f} rad/m²')

    ax.set_xlabel('Per-channel adopted RM [rad/m²]')
    ax.set_ylabel('Number of channels')
    ax.set_title('Per-channel adopted RM distribution (fine grid)\n'
                 'Dashed = coarse grid RM, Dotted = per-channel median')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.savefig(os.path.join(output_dir, 'xf_per_channel_rm_histogram.png'), dpi=150)
    plt.close(fig)


# ============================
# Main Script
# ============================
xfdir = GAINPLOTS + '/manualXF'
os.makedirs(xfdir, exist_ok=True)

# Define calibration tables (matching 1GC_05_casa_refcal.py)
ktab = GAINTABLES + '/cal_1GC_' + myms + '.K'
bptab = GAINTABLES + '/cal_1GC_' + myms + '.B'
gptab = GAINTABLES + '/cal_1GC_' + myms + '.Gp'
gtab = GAINTABLES + '/cal_1GC_' + myms + '.G'
dftab = GAINTABLES + '/cal_1GC_' + myms + '.Df'

xftab = GAINTABLES + '/cal_1GC_' + myms + '.Xf'

# Apply all calibration except XF (which we are solving for) to write
# CORRECTED_DATA into myms. IQUV will be read from CORRECTED_DATA directly.
applycal(vis=myms,
         field=pacal_name,
         parang=False,
         gaintable=[ktab, bptab, gptab, gtab, dftab],
         gainfield=[pacal_name, bpcal_name, pacal_name, pacal_name, bpcal_name],
         interp=['linear', 'linear', 'linear', 'linear', 'linear'],
         flagbackup=False)

# ============================
# Display Analysis Parameters
# ============================
print("\n=== Cross-Hand Phase (XF) Calibration Parameters ===")
print(f"Polarization calibrator: {pacal_name}")
print(f"Known source polarization angle: {XF_TARGET_POLANG}°")
print(f"Known source rotation measure: {XF_TARGET_RM} rad/m²")
print(f"RM coarse grid: XF_TARGET_RM + iono_RM ± 2.0 rad/m², step = iono_err/3")
print(f"RM fine grid: coarse_best ± 0.25 rad/m², step = 0.04 rad/m² (per channel)")
print(f"Channel averaging interval (polcal solint): {XF_CHANINT} channels")
print(f"Scan-averaged solution: {XF_AVG_SCAN}")
print(f"Global circular clip: {XF_GLOBAL_SIGMA_CLIP}σ")
print(f"Polynomial baseline order: {XF_POLY_ORDER}")
print(f"Sliding MAD window: {XF_MAD_WINDOW} channels (None = auto ceil(N/10))")
print(f"MAD sigma clip: {XF_MAD_SIGMA_CLIP}σ")
print(f"MAD iterative: {XF_MAD_ITERATIVE}")
print(f"SNR threshold: {XF_MIN_SN} (0 = disabled)")
print(f"Apply smoothing: {XF_USE_SMOOTHING}")
if XF_USE_SMOOTHING:
    print(f"  Savitzky-Golay window: {XF_SAVGOL_WINDOW} (None = auto)")
    print(f"  Savitzky-Golay polynomial order: {XF_SAVGOL_POLYORDER}")
print("=" * 60)

# ============================
# Extract Channel Frequencies
# ============================
print("\nExtracting channel frequencies from MS...")
tb.open(myms + '/SPECTRAL_WINDOW')
chan_freq = tb.getcol('CHAN_FREQ')
tb.close()
freq_ghz = (chan_freq / 1.0e9).flatten()
nchan = len(freq_ghz)
freq_min = freq_ghz.min()
freq_max = freq_ghz.max()
print(f"Spectral window: {nchan} channels")
print(f"Frequency range: {freq_min:.3f} - {freq_max:.3f} GHz")

# ============================
# Channel Averaging (if input > 1024 channels)
# ============================
# Build a channel-averaged temporary MS so downstream code always works on
# ≤1024 channels.  The temp MS writes the averaged CORRECTED_DATA into its
# DATA column; it is deleted after the per-scan visibility loading is done.
vis_to_load = myms
tmpms = None

if nchan > 1024:
    avg_factor = int(np.ceil(nchan / 1024))
    nchan_avg = int(np.ceil(nchan / avg_factor))
    print(f"\nInput has {nchan} channels (> 1024) — channel-averaging by factor "
          f"{avg_factor} to {nchan_avg} channels before extracting visibilities.")
    tmpms = myms.rstrip('/') + '_tmp_chanavg.ms'
    if os.path.isdir(tmpms):
        shutil.rmtree(tmpms)
    mstransform(vis=myms, outputvis=tmpms,
                field=pacal_name,
                chanaverage=True, chanbin=avg_factor,
                datacolumn='corrected',
                keepflags=True)
    tb.open(tmpms + '/SPECTRAL_WINDOW')
    chan_freq = tb.getcol('CHAN_FREQ')
    tb.close()
    freq_ghz = (chan_freq / 1.0e9).flatten()
    nchan = len(freq_ghz)
    freq_min = freq_ghz.min()
    freq_max = freq_ghz.max()
    vis_to_load = tmpms
    print(f"Averaged spectral window: {nchan} channels, "
          f"{freq_min:.3f} — {freq_max:.3f} GHz")

# ============================
# Retrieve Scan Information
# ============================
print("\nRetrieving scan information from MS...")
msmd.open(myms)
try:
    scan_numbers = msmd.scansforfield(msmd.fieldsforname(pacal_name)[0])
except Exception:
    msmd.close()
    sys.exit(f"ERROR: Could not retrieve scans for field {pacal_name}")

print(f"Found {len(scan_numbers)} scan(s) for {pacal_name}: {scan_numbers}")
n_scans = len(scan_numbers)
msmd.close()

# ============================
# Frequency-Dependent EVPA Model
# ============================
if '3c286' in pacal_name.lower() or 'j1331' in pacal_name.lower():
    print("\n" + "=" * 70)
    print("DETECTED 3C286 (J1331+3030) — APPLYING FREQUENCY-DEPENDENT EVPA MODEL")
    print("=" * 70)
    target_polang_array = compute_3c286_evpa(freq_ghz)
    print(f"EVPA range: {target_polang_array.min():.2f}° — {target_polang_array.max():.2f}°")
    print(f"EVPA at band centre ({np.median(freq_ghz):.2f} GHz): "
          f"{np.median(target_polang_array):.2f}°")

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(freq_ghz, target_polang_array, 'b-', linewidth=2, label='3C286 EVPA model')
    ax.axhline(XF_TARGET_POLANG, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
               label=f'Config value: {XF_TARGET_POLANG}°')
    ax.set_xlabel('Frequency [GHz]')
    ax.set_ylabel('EVPA [deg]')
    ax.set_title('3C286 Frequency-Dependent EVPA Model (Perley & Butler 2013)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(xfdir, '3c286_evpa_model.png'), dpi=150)
    plt.close(fig)
    print(f"Saved: {xfdir}/3c286_evpa_model.png\n" + "=" * 70)
else:
    print(f"\nUsing constant EVPA from config: {XF_TARGET_POLANG}° (not 3C286/J1331)")
    target_polang_array = np.full_like(freq_ghz, XF_TARGET_POLANG)

# ============================
# Load and Average Visibility Data Per Scan
# ============================
# If stokes_perscan.npz already exists from a previous run and is consistent
# with the current MS (same channel count and scan numbers), load from it
# directly and skip the expensive per-baseline weighted averaging step.

stokes_cache_path = os.path.join(xfdir, 'stokes_perscan.npz')
stokes_cached = False

if os.path.isfile(stokes_cache_path):
    print(f"\nFound cached Stokes data: {stokes_cache_path}")
    try:
        cache = np.load(stokes_cache_path)
        cached_nchan = int(cache['freq_ghz'].shape[0])
        cached_scans = cache['scan_numbers']
        if cached_nchan == nchan and np.array_equal(cached_scans, scan_numbers):
            scan_times    = cache['scan_times']
            vis_avg       = cache['vis_avg']          # (4, n_scans, nchan) complex
            stokes_cached = True
            print(f"  Cache consistent ({n_scans} scan(s) × {nchan} channels) — "
                  f"skipping MS visibility loading")
        else:
            print(f"  WARNING: Cache mismatch — cached {cached_nchan} ch / "
                  f"{len(cached_scans)} scan(s) vs current {nchan} ch / "
                  f"{n_scans} scan(s). Reloading from MS.")
    except Exception as e:
        print(f"  WARNING: Failed to load cache ({e}) — reloading from MS")

if not stokes_cached:
    print(f"\nAveraging visibilities per scan (collapse baseline/time axes)...")

    # Check once whether CORRECTED_DATA exists in this MS (split/averaged MSes
    # have DATA only; the original calibrated MS has CORRECTED_DATA).
    tb.open(vis_to_load)
    _has_corrected = 'CORRECTED_DATA' in tb.colnames()
    tb.close()
    _data_col = 'corrected_data' if _has_corrected else 'data'
    print(f"  Using data column: {_data_col.upper()}")

    vis_avg = np.zeros((4, n_scans, nchan), dtype=complex)
    scan_times = np.zeros(n_scans, dtype=float)

    total_rows = 0
    total_unflagged = 0

    for scan_idx, scan in enumerate(scan_numbers):
        ms.open(vis_to_load)
        ms.selectinit(reset=True)
        ok = ms.msselect({'field': pacal_name, 'scan': str(scan)})
        if not ok:
            ms.close()
            print(f"WARNING: Could not select scan {scan}, skipping")
            vis_avg[:, scan_idx, :] = np.nan
            continue

        ms.selectpolarization(['I', 'Q', 'U', 'V'])
        d = ms.getdata([_data_col, 'flag', 'weight_spectrum', 'weight', 'time'])
        ms.close()

        # Normalise to 'data' key regardless of which column was read
        if _data_col != 'data' and _data_col in d:
            d['data'] = d.pop(_data_col)

        missing = [k for k in ('data', 'flag') if k not in d]
        if missing:
            print(f"WARNING: Scan {scan} missing {missing}, skipping")
            vis_avg[:, scan_idx, :] = np.nan
            del d
            continue

        data_scan = d['data']
        raw_flag_scan = np.asarray(d['flag'])
        time_scan = d['time']
        ncorr, nchan_scan, nrow_scan = data_scan.shape

        if ncorr != 4:
            print(f"WARNING: Scan {scan} has {ncorr} correlations (expected 4), skipping")
            vis_avg[:, scan_idx, :] = np.nan
            del d, data_scan, raw_flag_scan, time_scan
            continue

        have_wspec = 'weight_spectrum' in d
        flag_iquv_scan = build_iquv_flag(raw_flag_scan, nchan_scan, nrow_scan)
        W_scan = build_iquv_weights(d, have_wspec, nchan_scan, nrow_scan)
        del d

        vis_avg_scan = compute_weighted_averages(data_scan, W_scan, flag_iquv_scan)
        vis_avg[:, scan_idx, :] = vis_avg_scan
        scan_times[scan_idx] = np.median(time_scan)

        W_eff_scan = np.where(flag_iquv_scan, 0.0, W_scan)
        unflagged_scan = int(np.sum(W_eff_scan > 0))
        total_scan = int(W_eff_scan.size)
        total_rows += nrow_scan
        total_unflagged += unflagged_scan

        print(f"  Scan {scan:3d}: {nrow_scan:5d} rows, "
              f"{unflagged_scan:8d}/{total_scan:8d} unflagged "
              f"({100.0 * unflagged_scan / total_scan:5.1f}%), "
              f"time={scan_times[scan_idx]:.1f} MJD")

        del data_scan, raw_flag_scan, time_scan, flag_iquv_scan, W_scan, W_eff_scan

    print(f"\nScan averaging complete. Total rows: {total_rows:,d}, "
          f"unflagged data: {total_unflagged:,d}")

    np.savez(stokes_cache_path,
             freq_ghz=freq_ghz, scan_numbers=scan_numbers, scan_times=scan_times,
             vis_avg=vis_avg)
    print(f"Saved: {stokes_cache_path}")

    if tmpms is not None and os.path.isdir(tmpms):
        print(f"\nRemoving temporary channel-averaged MS: {tmpms}")
        shutil.rmtree(tmpms)
        tmpms = None

# ============================
# Extract Stokes and Run Quality Checks (always, even on cached runs)
# ============================
# For a calibrated point source at the phase centre, Stokes flux is a real
# physical quantity. The real part of the averaged visibility carries the
# signal; the imaginary part reflects calibration residuals and should be
# small (<10%). We use np.real() throughout for computations (EVPA comparison,
# V→0 check) because EVPA and Stokes V are defined in terms of real flux
# densities. The complex vis_avg is retained in the cache solely so the
# imaginary fraction diagnostic can be re-run at any time.
I_scans = np.real(vis_avg[0])   # (n_scans, nchan)
Q_scans = np.real(vis_avg[1])
U_scans = np.real(vis_avg[2])
V_scans = np.real(vis_avg[3])

print("\n=== Data Quality Assessment ===")
for scan_idx, scan in enumerate(scan_numbers):
    print(f"\n--- Scan {scan} ---")
    for stokes_label, avg in zip('IQUV', vis_avg[:, scan_idx, :]):
        quick_stats(f'Scan {scan} Re({stokes_label})', np.real(avg))
    im_fracs = {s: calc_im_fraction(vis_avg[i, scan_idx, :])
                for i, s in enumerate('IQUV')}
    print(f"  Imaginary flux fractions (expect <<10% for point source):")
    for s, frac in im_fracs.items():
        status = "WARNING: HIGH" if frac > 10.0 else "OK"
        print(f"    {s}: {frac:.1f}% {status}")

# Pre-correction Stokes spectra (regenerate if not present, otherwise skip)
preXF_plot = os.path.join(xfdir, 'stokes_spectra_preXF.png')
if not os.path.isfile(preXF_plot):
    print("\nGenerating pre-correction Stokes spectra plot...")
    plot_stokes_spectra(freq_ghz, I_scans, Q_scans, U_scans, V_scans,
                        scan_numbers, xfdir,
                        filename='stokes_spectra_preXF.png',
                        title='Pre-XF Stokes Parameters')
    print(f"Saved: {preXF_plot}")
else:
    print(f"\nPre-XF Stokes plot already exists: {preXF_plot}")

# ============================
# Compute Parallactic Angles Per Scan
# ============================
print("\n=== Computing Parallactic Angles Per Scan ===")
chi_deg_scans = np.zeros(n_scans)

for scan_idx, scan in enumerate(scan_numbers):
    chi_deg, diagnostics = compute_parallactic_angle(myms, pacal_name,
                                                     scan_times[scan_idx])
    chi_deg_scans[scan_idx] = chi_deg
    print(f"  Scan {scan}: χ = {chi_deg:+7.3f}° "
          f"(AZ={diagnostics['az_deg']:6.1f}°, EL={diagnostics['el_deg']:5.1f}°)")

print(f"\nParallactic angle range: "
      f"{chi_deg_scans.min():+.3f}° to {chi_deg_scans.max():+.3f}°")

# ============================
# Spinifex Ionospheric RM
# ============================
print("\n=== Estimating Ionospheric RM via Spinifex ===")
iono_rm_per_scan, iono_rm_err_per_scan = get_spinifex_rm_per_scan(
    pacal_name, scan_numbers)

print(f"\nIonospheric RM summary:")
for scan_idx, scan in enumerate(scan_numbers):
    rm_centre = XF_TARGET_RM + iono_rm_per_scan[scan_idx]
    print(f"  Scan {scan}: iono_RM={iono_rm_per_scan[scan_idx]:+.3f} ± "
          f"{iono_rm_err_per_scan[scan_idx]:.3f} rad/m²  "
          f"→ RM grid centre = {rm_centre:.3f} rad/m²  "
          f"(range [{rm_centre - 2.0:.2f}, {rm_centre + 2.0:.2f}])")

# Weighted mean across scans (weight by 1/err² where err > 0)
with np.errstate(divide='ignore', invalid='ignore'):
    weights = np.where(iono_rm_err_per_scan > 0, 1.0 / iono_rm_err_per_scan**2, 0.0)
if np.sum(weights) > 0:
    mean_iono_rm = np.sum(weights * iono_rm_per_scan) / np.sum(weights)
    mean_iono_rm_err = 1.0 / np.sqrt(np.sum(weights))
else:
    mean_iono_rm = float(np.mean(iono_rm_per_scan))
    mean_iono_rm_err = float(np.std(iono_rm_per_scan))
print(f"\n  Weighted mean ionospheric RM: {mean_iono_rm:+.3f} ± {mean_iono_rm_err:.3f} rad/m²"
      f"  (XF_TARGET_RM + iono = {XF_TARGET_RM + mean_iono_rm:.3f} rad/m²)")

# ============================
# CASA polcal: Working-Channelisation XF Table
# ============================
print(f"\n=== Solving XF Table ({XF_CHANINT}ch intervals) ===")
print(f"combine='' (always per-scan; XF_AVG_SCAN={XF_AVG_SCAN} controls manual post-averaging)")

if os.path.isdir(xftab):
    shutil.rmtree(xftab)

polcal(vis=myms,
       field=pacal_name,
       caltable=xftab,
       refant=str(ref_ant),
       solint=f'inf,{XF_CHANINT}ch',
       poltype='Xf',
       combine='',
       gaintable=[ktab, bptab, gptab, gtab, dftab],
       gainfield=[pacal_name, bpcal_name, pacal_name, pacal_name, bpcal_name],
       interp=['linear', 'linear', 'linear', 'linear', 'linear'],
       append=False)

print(f"Created XF table: {xftab}")

# Read xftab frequency grid
tb.open(xftab + '/SPECTRAL_WINDOW')
chan_freq_xf = tb.getcol('CHAN_FREQ').flatten()
tb.close()
chan_freq_xf_ghz = chan_freq_xf / 1e9
n_xf_chan = len(chan_freq_xf_ghz)
print(f"xftab frequency grid: {n_xf_chan} channels "
      f"({chan_freq_xf_ghz.min():.3f} — {chan_freq_xf_ghz.max():.3f} GHz)")

# ============================
# Load XF Gains from Working-Channelisation Table
# ============================
print("\n=== Loading XF Gains ===")

tb.open(xftab)
gains_raw = tb.getcol('CPARAM')        # (1, n_xf_chan, n_rows)
flags_raw = tb.getcol('FLAG')          # (1, n_xf_chan, n_rows)
scans_tab = tb.getcol('SCAN_NUMBER')   # (n_rows,)
times_tab = tb.getcol('TIME')          # (n_rows,)
tb.close()

print(f"Table shape: {gains_raw.shape} (pol, chan, rows)")
print(f"Unique scans in table: {np.unique(scans_tab)}")

raw_gains_per_scan = {}   # {scan: complex (n_xf_chan,)}
raw_flags_per_scan = {}   # {scan: bool (n_xf_chan,)}

# IQUV interpolated onto xftab channel grid for ±π resolution
Q_xf_scans = np.zeros((n_scans, n_xf_chan))
U_xf_scans = np.zeros((n_scans, n_xf_chan))
V_xf_scans = np.zeros((n_scans, n_xf_chan))

for scan_idx, scan in enumerate(scan_numbers):
    scan_mask = scans_tab == scan
    if not np.any(scan_mask):
        print(f"WARNING: Scan {scan} not found in xftab")
        raw_gains_per_scan[scan] = np.full(n_xf_chan, np.nan, dtype=complex)
        raw_flags_per_scan[scan] = np.ones(n_xf_chan, dtype=bool)
        Q_xf_scans[scan_idx] = np.nan
        U_xf_scans[scan_idx] = np.nan
        V_xf_scans[scan_idx] = np.nan
        continue

    ant_indices = np.where(scan_mask)[0]
    g_scan = gains_raw[0, :, :][:, ant_indices]   # (n_xf_chan, n_ant)
    f_scan = flags_raw[0, :, :][:, ant_indices]   # (n_xf_chan, n_ant)

    g_masked = np.where(f_scan, np.nan + 0j, g_scan)
    with np.errstate(all='ignore'):
        med_real = np.nanmedian(np.real(g_masked), axis=-1)
        med_imag = np.nanmedian(np.imag(g_masked), axis=-1)
    median_g = med_real + 1j * med_imag

    all_flagged = np.all(f_scan, axis=-1)
    good = ~all_flagged & np.isfinite(median_g)
    n_good = int(np.sum(good))
    n_polcal_flagged = int(np.sum(all_flagged))
    print(f"  Scan {scan}: {int(np.sum(scan_mask))} antennas, "
          f"{n_good}/{n_xf_chan} good channels "
          f"({n_polcal_flagged} flagged by polcal)")

    median_g[~good] = np.nan
    raw_gains_per_scan[scan] = median_g
    raw_flags_per_scan[scan] = ~good

    # Interpolate IQUV from MS channels onto xftab channels
    for arr_ms, arr_xf in [(Q_scans[scan_idx], Q_xf_scans),
                            (U_scans[scan_idx], U_xf_scans),
                            (V_scans[scan_idx], V_xf_scans)]:
        finite = np.isfinite(arr_ms)
        if np.sum(finite) > 1:
            arr_xf[scan_idx] = np.interp(chan_freq_xf_ghz, freq_ghz[finite], arr_ms[finite])
        else:
            arr_xf[scan_idx] = np.nan

    # Flag xftab channels where IQUV is NaN after interpolation
    iquv_bad = (~np.isfinite(Q_xf_scans[scan_idx]) |
                ~np.isfinite(U_xf_scans[scan_idx]) |
                ~np.isfinite(V_xf_scans[scan_idx]))
    if np.any(iquv_bad):
        n_iquv = int(np.sum(iquv_bad & ~raw_flags_per_scan[scan]))
        if n_iquv > 0:
            print(f"  Scan {scan}: flagging {n_iquv} xftab channels with no valid IQUV")
        raw_flags_per_scan[scan] = raw_flags_per_scan[scan] | iquv_bad
        raw_gains_per_scan[scan][iquv_bad] = np.nan

# Also build a target_polang array on the xftab frequency grid
if '3c286' in pacal_name.lower() or 'j1331' in pacal_name.lower():
    target_polang_xf = compute_3c286_evpa(chan_freq_xf_ghz)
else:
    target_polang_xf = np.full(n_xf_chan, XF_TARGET_POLANG)

# ============================
# ±π Degeneracy Resolution Per Scan
# ============================
print("\n=== Resolving ±π Degeneracy Per Scan ===")

# Snapshot the raw CASA polcal output before any sign flip is applied.
# This is used in panel 1 of the diagnostic plot to show the ±π ambiguity.
pre_pi_gains_per_scan = {scan: raw_gains_per_scan[scan].copy() for scan in scan_numbers}
best_rm_per_scan    = {}
best_sign_per_scan  = {}
coarse_rm_per_scan  = {}
per_ch_rm_per_scan  = {}

for scan_idx, scan in enumerate(scan_numbers):
    print(f"\n  Scan {scan} (χ={chi_deg_scans[scan_idx]:+.2f}°, "
          f"iono_RM={iono_rm_per_scan[scan_idx]:+.3f} rad/m²)")

    median_g = raw_gains_per_scan[scan]

    signed_gains, best_rm, grid_centre_rm, delta_rm, med_dev, coarse_rm, per_ch_rm = \
        resolve_pi_degeneracy_scan(
            median_g,
            Q_xf_scans[scan_idx], U_xf_scans[scan_idx], V_xf_scans[scan_idx],
            chi_deg_scans[scan_idx], chan_freq_xf_ghz, target_polang_xf,
            iono_rm_per_scan[scan_idx],
            iono_rm_err=iono_rm_err_per_scan[scan_idx])

    print(f"    Best RM={best_rm:.4f} rad/m² (Δ={delta_rm:+.4f}), MAD deviation={med_dev:.3f}°")

    print_per_channel_resolution(
        chan_freq_xf_ghz, signed_gains, pre_pi_gains_per_scan[scan],
        raw_flags_per_scan[scan],
        Q_xf_scans[scan_idx], U_xf_scans[scan_idx], V_xf_scans[scan_idx],
        chi_deg_scans[scan_idx], best_rm, target_polang_xf)

    raw_gains_per_scan[scan] = signed_gains
    best_rm_per_scan[scan]   = best_rm
    best_sign_per_scan[scan] = None
    coarse_rm_per_scan[scan] = coarse_rm
    per_ch_rm_per_scan[scan] = per_ch_rm

# Rename for clarity in downstream code and diagnostic plots
resolved_gains_per_scan = raw_gains_per_scan

# Stage 1: raw phase + accepted/rejected residuals
print("\nGenerating Stage 1 raw XF phase plot...")
plot_raw_xf_phase_panels(chan_freq_xf_ghz, raw_gains_per_scan, raw_flags_per_scan,
                          pre_pi_gains_per_scan, best_sign_per_scan,
                          Q_xf_scans, U_xf_scans, V_xf_scans, chi_deg_scans,
                          best_rm_per_scan, target_polang_xf,
                          scan_numbers, xfdir)
print(f"Saved: {xfdir}/xf_phase_stage1_raw.png")

# Stage 2: post ±π resolution residuals (pre-flagging)
print("\nGenerating Stage 2 derotated angle residuals plot (post ±π resolution)...")
plot_derotated_angle_residuals(chan_freq_xf_ghz, Q_xf_scans, U_xf_scans, V_xf_scans,
                               resolved_gains_per_scan, raw_flags_per_scan,
                               chi_deg_scans, best_rm_per_scan,
                               target_polang_xf, scan_numbers, xfdir,
                               filename='xf_phase_stage2_post_pi.png',
                               per_ch_rm_per_scan=per_ch_rm_per_scan)
print(f"Saved: {xfdir}/xf_phase_stage2_post_pi.png")

# ============================
# Global Sign Sanity Check (80% rule)
# ============================
# For each scan, compute the fraction of valid channels whose de-rotated EVPA
# lies within 45° of the target. If fewer than 80% pass, the ±π resolution
# chose the wrong sign — force-flip that scan's gains now, before any
# clipping or table writing occurs.
print("\n=== Global Sign Sanity Check (80% rule) ===")
_c = 2.998e8
_any_flipped = False

for scan_idx, scan in enumerate(scan_numbers):
    g = resolved_gains_per_scan[scan]
    flags = raw_flags_per_scan[scan]
    valid = (~flags & np.isfinite(g) &
             np.isfinite(Q_xf_scans[scan_idx]) & np.isfinite(U_xf_scans[scan_idx]) &
             np.isfinite(V_xf_scans[scan_idx]))

    if np.sum(valid) < 5:
        print(f"  Scan {scan}: insufficient channels for sanity check — skipping")
        continue

    freq_hz = chan_freq_xf_ghz[valid] * 1e9
    lambda_sq = (_c / freq_hz) ** 2
    Q_v = Q_xf_scans[scan_idx][valid]
    U_v = U_xf_scans[scan_idx][valid]
    V_v = V_xf_scans[scan_idx][valid]
    g_v = g[valid]

    rho = np.angle(g_v)
    U_xh, V_xh = correct_crosshand_phase(U_v, V_v, rho)
    Q_sky, U_sky = correct_parallactic_angle(Q_v, U_xh, chi_deg_scans[scan_idx])
    evpa_deg = calculate_derotated_angle(Q_sky, U_sky, best_rm_per_scan[scan], lambda_sq)

    deviation = (evpa_deg - target_polang_xf[valid] + 90) % 180 - 90
    frac_good = float(np.sum(np.abs(deviation) < 45)) / float(np.sum(valid))

    if frac_good < 0.80:
        print(f"  Scan {scan}: only {frac_good*100:.1f}% of channels within 45° of target "
              f"— FORCE-FLIPPING SIGN")
        resolved_gains_per_scan[scan] = -resolved_gains_per_scan[scan]
        _any_flipped = True
    else:
        print(f"  Scan {scan}: {frac_good*100:.1f}% of channels within 45° of target — OK")

if _any_flipped:
    print("\nRe-generating Stage 2 plot after force-flip(s)...")
    plot_derotated_angle_residuals(chan_freq_xf_ghz, Q_xf_scans, U_xf_scans, V_xf_scans,
                                   resolved_gains_per_scan, raw_flags_per_scan,
                                   chi_deg_scans, best_rm_per_scan,
                                   target_polang_xf, scan_numbers, xfdir,
                                   filename='xf_phase_stage2_post_pi.png',
                                   per_ch_rm_per_scan=per_ch_rm_per_scan)
    print(f"Saved (updated): {xfdir}/xf_phase_stage2_post_pi.png")

# ============================
# Per-Channel Population Sign Check
# ============================
# Post sign-resolution but pre-clipping. For each channel, compute its
# circular distance from the band-wide median XF phase. If flipping the sign
# reduces that distance by more than 30° (well below π/2), it's almost
# certainly a misresolved ±π channel rather than noise — adopt the flip.
# Channels already flagged (NaN) are skipped.
print("\n=== Per-Channel Population Sign Check (threshold=30°) ===")
_flip_threshold = 30.0
_total_ch_flipped = 0

for scan_idx, scan in enumerate(scan_numbers):
    g = resolved_gains_per_scan[scan].copy()
    flags = raw_flags_per_scan[scan]
    valid = ~flags & np.isfinite(g)

    if np.sum(valid) < 10:
        print(f"  Scan {scan}: too few valid channels — skipping")
        continue

    # Circular median phase of the whole band
    phases = np.degrees(np.angle(g[valid]))
    circ_mean = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(phases)))))

    n_flipped = 0
    valid_idx = np.where(valid)[0]

    for ch in valid_idx:
        phase_ch = np.degrees(np.angle(g[ch]))
        phase_fl = np.degrees(np.angle(-g[ch]))

        # Circular distance to population median, wrapped to [0, 180]
        dist_now  = abs(((phase_ch - circ_mean + 180) % 360) - 180)
        dist_flip = abs(((phase_fl - circ_mean + 180) % 360) - 180)

        if dist_flip < dist_now - _flip_threshold:
            g[ch] = -g[ch]
            n_flipped += 1

    if n_flipped > 0:
        resolved_gains_per_scan[scan] = g
        print(f"  Scan {scan}: {n_flipped} channels flipped "
              f"(reduced distance to population by >{_flip_threshold}°)")
    else:
        print(f"  Scan {scan}: 0 channels flipped — all consistent with population")

    _total_ch_flipped += n_flipped

print(f"  Total channels flipped by population check: {_total_ch_flipped}")

# ============================
# Global 10σ Circular Phase Clip
# ============================
# Work on copies so we can show resolved vs post-clipping in the diagnostic plot
working_gains = {scan: resolved_gains_per_scan[scan].copy() for scan in scan_numbers}
working_flags = {scan: raw_flags_per_scan[scan].copy() for scan in scan_numbers}

print(f"\n=== Global {XF_GLOBAL_SIGMA_CLIP}σ Circular Phase Clip ===")
working_flags, global_flagged_per_scan = global_sigma_clip_gains(
    working_gains, working_flags, scan_numbers, sigma=float(XF_GLOBAL_SIGMA_CLIP))

# Snapshot flags after global clip so SNR and MAD flags can be identified separately
post_global_flags = {scan: working_flags[scan].copy() for scan in scan_numbers}

# ============================
# Cross-Hand SNR Flagging
# ============================
print(f"\n=== Cross-Hand SNR Flagging (threshold = {XF_MIN_SN}) ===")

snr_flagged_per_scan = {scan: np.zeros(n_xf_chan, dtype=bool) for scan in scan_numbers}
snr_per_scan         = {scan: np.full(n_xf_chan, np.nan) for scan in scan_numbers}

if XF_MIN_SN > 0:
    # Estimate noise at native resolution — more data points and finer windows
    # than operating on the XF-averaged grid.
    win_size_nat = max(3, nchan // 100)
    all_mads = []
    for scan_idx, scan in enumerate(scan_numbers):
        qu_nat = np.sqrt(U_scans[scan_idx]**2 + V_scans[scan_idx]**2)
        for start in range(0, nchan, win_size_nat):
            seg   = qu_nat[start:start + win_size_nat]
            valid = seg[np.isfinite(seg)]
            if len(valid) >= 3:
                all_mads.append(float(np.median(np.abs(valid - np.median(valid)))))
    if all_mads:
        noise = float(np.median(all_mads)) * 1.4826
        print(f"  Estimated noise: {noise:.4e} Jy "
              f"(pooled {len(all_mads)} windows of {win_size_nat} native ch, "
              f"{n_scans} scan(s))")
    else:
        noise = 0.0
        print("  WARNING: could not estimate noise — SNR flagging skipped")

    if noise > 0:
        for scan_idx, scan in enumerate(scan_numbers):
            qu_nat = np.sqrt(U_scans[scan_idx]**2 + V_scans[scan_idx]**2)
            qu_xf  = _bin_qu_native_to_xf(qu_nat, n_xf_chan)
            snr    = np.where(np.isfinite(qu_xf), qu_xf / noise, np.nan)
            # NaN any channel already flagged for any reason at this point
            already_flagged = raw_flags_per_scan[scan] | working_flags[scan]
            snr = np.where(already_flagged, np.nan, snr)
            snr_per_scan[scan] = snr
            low_snr = np.isfinite(snr) & (snr < XF_MIN_SN) & ~working_flags[scan]
            n_new = int(np.sum(low_snr))
            working_flags[scan] = working_flags[scan] | low_snr
            snr_flagged_per_scan[scan] = low_snr

            _snr_fin = snr[np.isfinite(snr)]
            _snr_rng = (f"[{_snr_fin.min():.2f}, {_snr_fin.max():.2f}]"
                        if len(_snr_fin) else "[all NaN]")
            print(f"  Scan {scan}: {n_new} channels flagged (SNR < {XF_MIN_SN:.1f}), "
                  f"SNR range {_snr_rng}")
            print(f"    {'ch':>5}  {'Freq(GHz)':>10}  {'SNR':>10}  {'flag':>8}")
            print(f"    {'-'*40}")
            for ch in range(n_xf_chan):
                if np.isnan(snr[ch]):
                    print(f"    {ch:>5}  {chan_freq_xf_ghz[ch]:>10.4f}  {'---':>10}  {'(pre)':>8}")
                else:
                    flag_str = 'FLAGGED' if low_snr[ch] else 'ok'
                    print(f"    {ch:>5}  {chan_freq_xf_ghz[ch]:>10.4f}  {snr[ch]:>10.2f}  {flag_str:>8}")
else:
    print("  XF_MIN_SN <= 0 — SNR flagging disabled")

# Snapshot flags after SNR clip so MAD flags can be identified separately
post_snr_flags = {scan: working_flags[scan].copy() for scan in scan_numbers}

# ============================
# Polynomial Baseline + Sliding Circular MAD Clip
# ============================
poly_order = int(XF_POLY_ORDER) if XF_POLY_ORDER and int(XF_POLY_ORDER) > 0 else 3
import math
_auto_window = math.ceil(n_xf_chan / 5)
if _auto_window % 2 == 0:
    _auto_window += 1
mad_window = int(XF_MAD_WINDOW) if XF_MAD_WINDOW else _auto_window
mad_sigma  = float(XF_MAD_SIGMA_CLIP) if XF_MAD_SIGMA_CLIP else 5.0
print(f"  MAD window: {mad_window} channels "
      f"({'user-specified' if XF_MAD_WINDOW else f'auto: ceil_odd(ceil({n_xf_chan}/5))'})")

iter_str = f"iterative, max {20} passes" if XF_MAD_ITERATIVE else "single pass"
print(f"\n=== Polynomial Baseline + Sliding Circular MAD Clip "
      f"(poly_order={poly_order}, window={mad_window}, {mad_sigma}σ, {iter_str}) ===")

poly_coeffs_per_scan  = {scan: (None, None) for scan in scan_numbers}
poly_flagged_per_scan = {scan: np.zeros(n_xf_chan, dtype=bool) for scan in scan_numbers}

for scan_idx, scan in enumerate(scan_numbers):
    n_before = int(np.sum(~working_flags[scan]))
    print(f"\n  Scan {scan} ({n_before}/{n_xf_chan} unflagged channels):")
    U_w = U_xf_scans[scan_idx]
    V_w = V_xf_scans[scan_idx]
    poly_weights = np.log1p(np.sqrt(U_w**2 + V_w**2) + 1e-8)
    poly_weights = np.where(np.isfinite(poly_weights) & ~working_flags[scan],
                            poly_weights, 0.0)
    working_flags[scan], p_cos, p_sin, n_clipped = poly_baseline_mad_clip(
        working_gains[scan], working_flags[scan],
        n_xf_chan, poly_order, mad_window, mad_sigma,
        iterative=bool(XF_MAD_ITERATIVE), weights=poly_weights)
    poly_coeffs_per_scan[scan]  = (p_cos, p_sin)
    # Only channels newly flagged after the SNR step are counted as MAD-flagged
    poly_flagged_per_scan[scan] = working_flags[scan] & ~post_snr_flags[scan]
    n_after = int(np.sum(~working_flags[scan]))
    print(f"  Scan {scan}: {n_clipped} channels flagged "
          f"({n_after}/{n_xf_chan} remaining)")

# Store flagged gains and flag sets for diagnostic plot (5-tuple: global + poly + snr)
flagged_gains_per_scan = {
    scan: (working_gains[scan].copy(),
           working_flags[scan].copy(),
           global_flagged_per_scan[scan].copy(),
           poly_flagged_per_scan[scan].copy(),
           snr_flagged_per_scan[scan].copy())
    for scan in scan_numbers}

# Stage 3: post flagging residuals (surviving channels only)
print("\nGenerating Stage 3 derotated angle residuals plot (post flagging)...")
plot_derotated_angle_residuals(chan_freq_xf_ghz, Q_xf_scans, U_xf_scans, V_xf_scans,
                               resolved_gains_per_scan, raw_flags_per_scan,
                               chi_deg_scans, best_rm_per_scan,
                               target_polang_xf, scan_numbers, xfdir,
                               filename='xf_phase_stage3_post_flagging.png',
                               flag_override=working_flags,
                               per_ch_rm_per_scan=per_ch_rm_per_scan)
print(f"Saved: {xfdir}/xf_phase_stage3_post_flagging.png")

# Per-channel RM histogram
print("\nGenerating per-channel RM histogram...")
plot_per_channel_rm_histogram(chan_freq_xf_ghz, per_ch_rm_per_scan, coarse_rm_per_scan,
                               raw_flags_per_scan, scan_numbers, xfdir)
print(f"Saved: {xfdir}/xf_per_channel_rm_histogram.png")

# ============================
# Optional Savitzky-Golay Smoothing
# ============================
if XF_USE_SMOOTHING:
    print("\n=== Applying Savitzky-Golay Smoothing ===")
    for scan in scan_numbers:
        good = ~working_flags[scan] & np.isfinite(working_gains[scan])
        n_good = int(np.sum(good))
        if n_good < 5:
            continue
        window_length = XF_SAVGOL_WINDOW if XF_SAVGOL_WINDOW is not None else min(5, n_good)
        if window_length % 2 == 0:
            window_length += 1
        window_length = min(window_length, n_good)
        if window_length % 2 == 0:
            window_length -= 1
        polyorder = min(XF_SAVGOL_POLYORDER, window_length - 1)
        g_good = working_gains[scan][good]
        re_smooth = savgol_filter(np.real(g_good), window_length, polyorder)
        im_smooth = savgol_filter(np.imag(g_good), window_length, polyorder)
        smoothed = re_smooth + 1j * im_smooth
        smoothed = smoothed / np.abs(smoothed)
        working_gains[scan][good] = smoothed
        print(f"  Scan {scan}: smoothed {n_good} channels "
              f"(window={window_length}, poly={polyorder})")

# ============================
# Write Solutions to xftab
# ============================
print("\n=== Writing Solutions to xftab ===")

tb.open(xftab, nomodify=False)
gains_out = tb.getcol('CPARAM')     # (1, n_xf_chan, n_rows)
flags_out = tb.getcol('FLAG')       # (1, n_xf_chan, n_rows)
time_out  = tb.getcol('TIME')       # (n_rows,)
scans_out = tb.getcol('SCAN_NUMBER')
tb.close()

print(f"xftab dimensions: {gains_out.shape}")

final_gains_per_scan = {}

for scan_idx, scan in enumerate(scan_numbers):
    scan_mask_out = scans_out == scan
    if not np.any(scan_mask_out):
        print(f"WARNING: Scan {scan} not found in xftab rows — skipping")
        continue

    ant_idx_out = np.where(scan_mask_out)[0]
    g_final = working_gains[scan].copy()
    f_final = working_flags[scan].copy()

    # Enforce unit amplitude on all unflagged channels
    valid_g = ~f_final & np.isfinite(g_final)
    if np.any(valid_g):
        g_final[valid_g] = g_final[valid_g] / np.abs(g_final[valid_g])

    # NaN flagged channels for diagnostic plot
    g_diag = g_final.copy()
    g_diag[f_final] = np.nan

    # Also NaN channels flagged in polcal table itself
    polcal_flagged = flags_out[0, :, :][:, ant_idx_out[0]].astype(bool)
    g_diag[polcal_flagged] = np.nan

    final_gains_per_scan[scan] = g_diag

    gains_out[0, :, :][:, ant_idx_out] = np.where(
        f_final[:, np.newaxis], 1.0 + 0.0j, g_final[:, np.newaxis])
    flags_out[0, :, :][:, ant_idx_out] = np.where(
        f_final[:, np.newaxis], True, flags_out[0, :, :][:, ant_idx_out])

    time_out[scan_mask_out] = scan_times[scan_idx]

    n_good = int(np.sum(~f_final))
    n_flagged = int(np.sum(f_final))
    print(f"  Scan {scan}: wrote {n_good}/{n_xf_chan} valid channels "
          f"({n_flagged} flagged)")

tb.open(xftab, nomodify=False)
tb.putcol('CPARAM', gains_out)
tb.putcol('FLAG', flags_out)
tb.putcol('TIME', time_out)
tb.flush()
tb.close()

print(f"\nSUCCESS: XF table populated: {xftab}")

# ============================
# Cross-Scan Flux-Weighted Circular Average
# ============================
# When XF_AVG_SCAN=True, compute a flux-weighted circular phasor average across
# all per-scan solutions and overwrite every row with the combined solution.
# Weight per (scan, channel) = sqrt(Q²+U²) on the xftab frequency grid.
if XF_AVG_SCAN and n_scans > 1:
    print("\n=== Cross-Scan Flux-Weighted Circular Average (XF_AVG_SCAN=True) ===")

    # Build per-scan QU weights on the xftab frequency grid
    scan_weights_xf = {}
    for scan_idx, scan in enumerate(scan_numbers):
        U_xf = U_xf_scans[scan_idx]
        V_xf = V_xf_scans[scan_idx]
        w = np.where(np.isfinite(U_xf) & np.isfinite(V_xf),
                     np.sqrt(U_xf**2 + V_xf**2), 0.0)
        w = np.where(w > 0, w, 0.0)
        scan_weights_xf[scan] = w
        print(f"  Scan {scan}: median cross-hand weight = {np.nanmedian(w[w>0]):.4f} Jy"
              if np.any(w > 0) else f"  Scan {scan}: no valid cross-hand flux")

    # Accumulate weighted phasors across scans
    wsum_re = np.zeros(n_xf_chan)
    wsum_im = np.zeros(n_xf_chan)
    wsum    = np.zeros(n_xf_chan)
    flag_all = np.ones(n_xf_chan, dtype=bool)

    for scan_idx, scan in enumerate(scan_numbers):
        g = working_gains[scan]
        f = working_flags[scan]
        w = scan_weights_xf[scan]
        good = ~f & np.isfinite(g)
        ph   = np.angle(g)
        eff_w = np.where(good, w, 0.0)
        wsum_re  += eff_w * np.cos(ph)
        wsum_im  += eff_w * np.sin(ph)
        wsum     += eff_w
        flag_all &= f   # channel flagged only if flagged in ALL scans

    has_data = wsum > 0.0
    amp      = np.sqrt(wsum_re**2 + wsum_im**2)
    valid    = has_data & (amp > 1e-10)
    safe_amp = np.where(valid, amp, 1.0)
    avg_re   = np.where(valid, wsum_re / safe_amp, 0.0)
    avg_im   = np.where(valid, wsum_im / safe_amp, 0.0)
    avg_gains_combined = (avg_re + 1j * avg_im).astype(complex)
    avg_flag_combined  = ~has_data | flag_all

    n_avg_good = int(np.sum(~avg_flag_combined))
    print(f"  Combined solution: {n_avg_good}/{n_xf_chan} valid channels")

    # Overwrite every row in the table with the combined solution
    tb.open(xftab, nomodify=False)
    gains_out2 = tb.getcol('CPARAM')
    flags_out2 = tb.getcol('FLAG')
    n_rows_tab = gains_out2.shape[2]
    for row in range(n_rows_tab):
        gains_out2[0, :, row] = np.where(
            avg_flag_combined, 1.0 + 0.0j, avg_gains_combined)
        flags_out2[0, :, row] = avg_flag_combined
    tb.putcol('CPARAM', gains_out2)
    tb.putcol('FLAG',   flags_out2)
    tb.flush()
    tb.close()

    # Update final_gains_per_scan for diagnostic plots
    avg_diag = avg_gains_combined.copy()
    avg_diag[avg_flag_combined] = np.nan
    for scan in scan_numbers:
        final_gains_per_scan[scan] = avg_diag

    print(f"  All {n_rows_tab} table rows overwritten with combined solution.")

# ============================
# Diagnostic Plots
# ============================
print("\n=== Generating Diagnostic Plots ===")

# XF phase diagnostic (raw → flagged → final)
print("Generating XF phase diagnostic plot...")
plot_xf_phase_diagnostic(chan_freq_xf_ghz, chan_freq_xf_ghz,
                          pre_pi_gains_per_scan,
                          resolved_gains_per_scan,
                          flagged_gains_per_scan,
                          poly_coeffs_per_scan,
                          final_gains_per_scan,
                          snr_per_scan,
                          scan_numbers, xfdir)
print(f"Saved: {xfdir}/xf_phase_diagnostic.png")

# Post-correction Stokes spectra and V-zeroing check
print("\nApplying solutions to IQUV for post-correction diagnostic...")

I_corrected_scans = np.zeros((n_scans, nchan))
Q_corrected_scans = np.zeros((n_scans, nchan))
U_corrected_scans = np.zeros((n_scans, nchan))
V_corrected_scans = np.zeros((n_scans, nchan))

for scan_idx, scan in enumerate(scan_numbers):
    g_xf = working_gains[scan].copy()
    f_xf = working_flags[scan].copy()
    valid_xf = ~f_xf & np.isfinite(g_xf)

    if np.sum(valid_xf) < 2:
        I_corrected_scans[scan_idx] = np.nan
        Q_corrected_scans[scan_idx] = np.nan
        U_corrected_scans[scan_idx] = np.nan
        V_corrected_scans[scan_idx] = np.nan
        continue

    # Interpolate xftab gains onto MS full-resolution frequency grid
    freq_src  = chan_freq_xf_ghz[valid_xf]
    gains_src = g_xf[valid_xf]
    gains_src = gains_src / np.abs(gains_src)  # ensure unit amplitude

    real_interp = np.interp(freq_ghz, freq_src, np.real(gains_src))
    imag_interp = np.interp(freq_ghz, freq_src, np.imag(gains_src))
    complex_interp = real_interp + 1j * imag_interp
    complex_interp = complex_interp / np.abs(complex_interp)

    # CASA applies conj(gain) to data
    rho_interp = np.angle(np.conj(complex_interp))

    U_corr, V_corr = correct_crosshand_phase(
        U_scans[scan_idx], V_scans[scan_idx], rho_interp)

    I_corrected_scans[scan_idx] = I_scans[scan_idx]
    Q_corrected_scans[scan_idx] = Q_scans[scan_idx]
    U_corrected_scans[scan_idx] = U_corr
    V_corrected_scans[scan_idx] = V_corr

    P_corr = np.sqrt(Q_scans[scan_idx] ** 2 + U_corr ** 2)
    P_frac = np.nanmedian(P_corr / I_scans[scan_idx]) * 100
    V_frac = np.nanmedian(np.abs(V_corr) / I_scans[scan_idx]) * 100
    print(f"  Scan {scan}: P/I = {P_frac:.3f}%, |V|/I = {V_frac:.3f}%")

# Post-correction check: optionally apply xftab to myms and extract full IQUV
I_ms_scans = None
Q_ms_scans = None
U_ms_scans = None
V_ms_scans = None

if XF_APPLY_TO_MS:
    print("\n=== Applying XF table to myms and extracting MS IQUV ===")

    applycal(vis=myms,
             field=pacal_name,
             parang=False,
             gaintable=[ktab, bptab, gptab, gtab, dftab, xftab],
             gainfield=[pacal_name, bpcal_name, pacal_name, pacal_name,
                        bpcal_name, pacal_name],
             interp=['linear', 'linear', 'linear', 'linear', 'linear', 'linear'],
             flagbackup=False)

    print("  Extracting IQUV from CORRECTED_DATA in myms...")
    I_ms_scans = np.full((n_scans, nchan), np.nan)
    Q_ms_scans = np.full((n_scans, nchan), np.nan)
    U_ms_scans = np.full((n_scans, nchan), np.nan)
    V_ms_scans = np.full((n_scans, nchan), np.nan)

    for scan_idx, scan in enumerate(scan_numbers):
        ms.open(myms)
        ms.selectinit(reset=True)
        ok = ms.msselect({'field': pacal_name, 'scan': str(scan)})
        if not ok:
            ms.close()
            print(f"  WARNING: Could not select field={pacal_name}, scan={scan}")
            continue

        ms.selectpolarization(['I', 'Q', 'U', 'V'])
        d = ms.getdata(['corrected_data', 'flag', 'weight_spectrum', 'weight'])
        ms.close()

        if 'corrected_data' not in d or 'flag' not in d:
            print(f"  WARNING: Scan {scan} missing corrected_data")
            continue

        d['data'] = d.pop('corrected_data')
        data_scan = d['data']
        raw_flag_scan = np.asarray(d['flag'])
        ncorr, nchan_scan, nrow_scan = data_scan.shape

        if ncorr != 4:
            print(f"  WARNING: Scan {scan} has {ncorr} correlations, skipping")
            continue

        have_wspec = 'weight_spectrum' in d
        flag_iquv_scan = build_iquv_flag(raw_flag_scan, nchan_scan, nrow_scan)
        W_scan = build_iquv_weights(d, have_wspec, nchan_scan, nrow_scan)
        del d

        vis_avg_scan = compute_weighted_averages(data_scan, W_scan, flag_iquv_scan)
        I_ms_scans[scan_idx] = np.real(vis_avg_scan[0])
        Q_ms_scans[scan_idx] = np.real(vis_avg_scan[1])
        U_ms_scans[scan_idx] = np.real(vis_avg_scan[2])
        V_ms_scans[scan_idx] = np.real(vis_avg_scan[3])

        with np.errstate(all='ignore'):
            V_med = np.nanmedian(np.abs(V_ms_scans[scan_idx]))
        print(f"  Scan {scan}: median |V_MS| = {V_med:.4f} Jy")

        del data_scan, raw_flag_scan, flag_iquv_scan, W_scan

# Post-XF Stokes spectra — analytic correction
plot_stokes_spectra(freq_ghz, I_corrected_scans, Q_corrected_scans,
                    U_corrected_scans, V_corrected_scans,
                    scan_numbers, xfdir,
                    filename='stokes_spectra_postXF_analytic.png',
                    title='Post-XF Stokes Parameters (Analytic correction)')
print(f"Saved: {xfdir}/stokes_spectra_postXF_analytic.png")

# Post-XF Stokes spectra — CASA MS CORRECTED_DATA (if available)
if I_ms_scans is not None:
    plot_stokes_spectra(freq_ghz, I_ms_scans, Q_ms_scans,
                        U_ms_scans, V_ms_scans,
                        scan_numbers, xfdir,
                        filename='stokes_spectra_postXF_casa.png',
                        title='Post-XF Stokes Parameters (CASA CORRECTED_DATA)')
    print(f"Saved: {xfdir}/stokes_spectra_postXF_casa.png")

# Save corrected spectra
np.savez(os.path.join(xfdir, 'final_corrected_stokes_spectra.npz'),
         freq_ghz=freq_ghz,
         scan_numbers=scan_numbers,
         I_corrected_scans=I_corrected_scans,
         Q_corrected_scans=Q_corrected_scans,
         U_corrected_scans=U_corrected_scans,
         V_corrected_scans=V_corrected_scans)
print(f"Saved: {xfdir}/final_corrected_stokes_spectra.npz")

print("\n=== Script Complete ===")
