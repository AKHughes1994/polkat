# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
Manual Cross-Hand Phase (XF) Polarization Calibration Solver - Per-Scan Version

This script performs sophisticated polarization calibration analysis on radio 
interferometric data to solve for instrumental cross-hand phase corrections
on a per-scan basis, properly accounting for time-variable parallactic angles.

SCIENTIFIC BACKGROUND:
The cross-hand phase (XF) represents the phase offset between the XY and YX 
correlations in a dual-polarization feed system. This instrumental effect 
rotates the observed Stokes U and V parameters in the UV plane, which must 
be corrected before accurate polarization measurements can be made.

For a linearly polarized source, the cross-hand phase ρ can be determined from:
    ρ = arctan2(-V_obs, U_obs)

However, this measurement has a ±π degeneracy. This script resolves the 
degeneracy by comparing Faraday-derotated polarization angles to a known 
source polarization angle across a range of rotation measures (RM).

WORKFLOW SUMMARY:
1. Load and average visibility data per scan (collapse baseline/time axes)
2. Compute parallactic angles at each scan's mid-time
3. Calculate raw cross-hand phase: ρ = arctan2(-V, U)
4. Three-stage RM trial analysis to find global RM and resolve ±π per channel
5. Flag outliers based on cross-hand flux and local phase scatter
6. Apply frequency averaging if needed (complex-plane vector averaging)
7. Generate CASA XF calibration table (scan-averaged or per-scan)
8. Apply corrections and generate diagnostic plots

DIAGNOSTIC OUTPUTS:
- Pre-XF Stokes I,Q,U,V spectra (multi-scan)
- Cross-hand flux and phase before correction
- Delta_RM analysis (trial distribution and angle convergence)
- Final rho values with outlier flagging
- Post-XF corrected Stokes spectra (zoomed to show structure)
"""

# Standard library imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys
import os
import glob
import shutil
import time
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

def stamp():
    """Generate timestamp string for logging."""
    now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
    return now


# ============================
# FUNCTIONS
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
            f_iq = f_xx | f_yy  # I,Q from XX|YY
            f_uv = f_xy | f_yx  # U,V from XY|YX
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
        raise ValueError(f"Unrecognised FLAG shape: {raw_flag.shape}. Expected (4,nchan,nrow), (1,nchan,nrow), or (nchan,nrow).")


def build_iquv_weights(weight_data, have_wspec, nchan, nrow):
    """
    Convert raw WEIGHT or WEIGHT_SPECTRUM to IQUV-compatible weight array.
    
    Inputs:
        weight_data : dict - Dictionary containing 'weight_spectrum' or 'weight' from MS
        have_wspec : bool - True if WEIGHT_SPECTRUM available, False if only WEIGHT
        nchan : int - Number of frequency channels
        nrow : int - Number of rows
    
    Outputs:
        W : array - Weight array with shape (4, nchan, nrow) for I,Q,U,V
    
    Maps correlation weights to Stokes parameter weights:
        - I,Q weights = average of XX and YY weights
        - U,V weights = average of XY and YX weights
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
                print("Converting 4-corr WEIGHT_SPECTRUM (XX,XY,YX,YY) → I,Q,U,V weights.")
                w_xx = W[0, :, :]
                w_xy = W[1, :, :]
                w_yx = W[2, :, :]
                w_yy = W[3, :, :]
                w_iq = 0.5*(w_xx + w_yy)
                w_uv = 0.5*(w_xy + w_yx)
                W = np.stack([w_iq, w_iq, w_uv, w_uv], axis=0)
            elif c == 1:
                print("Single-pol WEIGHT_SPECTRUM → broadcasting to (4, nchan, nrow).")
                W = np.broadcast_to(W, (4, nchan, nrow))
            else:
                raise ValueError(f"Unexpected number of correlation planes in WEIGHT_SPECTRUM: {c}. Expected 1 or 4.")
        else:
            raise ValueError(f"Unexpected WEIGHT_SPECTRUM dimensions: {W.shape}. Expected (nchan,nrow) or (corr,nchan,nrow).")
    else:
        Wraw = np.asarray(weight_data['weight'])
        print(f"Using WEIGHT with original shape {Wraw.shape}")

        if Wraw.ndim == 1 and Wraw.shape[0] == nrow:
            print("Single weight per row → broadcasting to (4, nchan, nrow).")
            wrow = Wraw
            W = np.ones((4, nchan, nrow), dtype=float) * wrow[None, None, :]
        elif Wraw.ndim == 2 and Wraw.shape[-1] == nrow:
            c = Wraw.shape[0]
            if c == 4:
                print("Converting 4-corr WEIGHT (XX,XY,YX,YY) → I,Q,U,V weights.")
                w_iq = 0.5*(Wraw[0, :] + Wraw[3, :])
                w_uv = 0.5*(Wraw[1, :] + Wraw[2, :])
                W = np.zeros((4, nchan, nrow), dtype=float)
                W[0, :, :] = w_iq[None, :]  # I
                W[1, :, :] = w_iq[None, :]  # Q
                W[2, :, :] = w_uv[None, :]  # U
                W[3, :, :] = w_uv[None, :]  # V
            elif c == 1:
                print("Single correlation WEIGHT → broadcasting to (4, nchan, nrow).")
                W = np.ones((4, nchan, nrow), dtype=float) * Wraw[0, :][None, None, :]
            else:
                raise ValueError(f"Unexpected number of correlations in WEIGHT: {c}. Expected 1 or 4.")
        else:
            raise ValueError(f"Unexpected WEIGHT shape: {Wraw.shape}. Expected (nrow,) or (corr,nrow).")
    
    return W


def compute_weighted_averages(data, weights, flags):
    """
    Compute weighted vector averages per channel, accounting for flags.
    
    Inputs:
        data : array - Complex visibility data with shape (4, nchan, nrow)
        weights : array - Weight array with shape (4, nchan, nrow)
        flags : array - Flag array with shape (4, nchan, nrow)
    
    Outputs:
        vis_avg : array - Complex averaged visibilities with shape (4, nchan)
        sigma_proxy : array - Uncertainty proxy (1/sqrt(sum(weights))) with shape (4, nchan)
    
    Applies flags by zeroing weights, then computes weighted mean per channel.
    Returns both averaged visibilities and uncertainty estimates.
    """
    W_eff = np.where(flags, 0.0, weights)
    
    den = np.sum(W_eff, axis=-1)  # (4, nchan)
    num = np.sum(W_eff * data, axis=-1, dtype=np.complex128)  # (4, nchan)
    
    with np.errstate(invalid='ignore', divide='ignore'):
        vis_avg = num / den
        sigma_proxy = 1.0 / np.sqrt(den)
    
    return vis_avg, sigma_proxy


def align_stokes_arrays(freq_ghz, I_flux, Q_flux, U_flux, V_flux):
    """
    Align Stokes parameter arrays to common minimum length.
    
    Inputs:
        freq_ghz : array - Frequency array in GHz
        I_flux, Q_flux, U_flux, V_flux : arrays - Stokes parameter flux arrays
    
    Outputs:
        f, I, Q, U, V : arrays - Aligned arrays truncated to common length
    
    Handles cases where frequency and Stokes arrays may have different lengths
    due to processing steps. Truncates all to minimum common length.
    """
    m = min(freq_ghz.size, I_flux.size, Q_flux.size, U_flux.size, V_flux.size)
    if m != freq_ghz.size:
        print(f"Note: Aligning arrays to {m} entries (was {freq_ghz.size}).")
    return freq_ghz[:m], I_flux[:m], Q_flux[:m], U_flux[:m], V_flux[:m]


def quick_stats(name, arr):
    """
    Print quick statistics for an array.
    
    Inputs:
        name : str - Name of the array for display
        arr : array - Numerical array to analyze
    
    Outputs:
        None (prints to console)
    
    Displays median, MAD (Median Absolute Deviation), and standard deviation
    for finite values in the array.
    """
    finite = np.isfinite(arr)
    if not np.any(finite):
        print(f'{name}: no finite values')
        return
    med = np.nanmedian(arr[finite])
    mad = np.nanmedian(np.abs(arr[finite] - med))
    std = np.nanstd(arr[finite])
    print(f'{name}: median={med:.6g}, MAD={mad:.6g}, std={std:.6g}')


def calc_im_fraction(complex_avg):
    """
    Calculate fraction of flux density in imaginary component.
    
    Inputs:
        complex_avg : array - Complex visibility array
    
    Outputs:
        im_frac : float - Median percentage of total amplitude in imaginary part
    
    For a well-calibrated point source, imaginary flux should be negligible (<10%).
    High imaginary fractions indicate calibration problems.
    """
    real_part = np.abs(np.real(complex_avg))
    imag_part = np.abs(np.imag(complex_avg))
    total_amp = np.sqrt(real_part**2 + imag_part**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        fraction = imag_part / total_amp
    return np.nanmedian(fraction) * 100


def q_to_rad_scalar(q):
    """
    Convert CASA quantity to radians, returning a Python float.
    
    Inputs:
        q : CASA quantity - Angle quantity from CASA
    
    Outputs:
        rad : float - Angle in radians as Python scalar
    
    Avoids NumPy 1.25+ warnings by ensuring scalar return type.
    """
    return float(np.atleast_1d(qa.getvalue(qa.convert(q, 'rad')))[0])


def rad_to_deg_scalar(xrad):
    """
    Convert radians to degrees using CASA, returning a Python float.
    
    Inputs:
        xrad : float - Angle in radians
    
    Outputs:
        deg : float - Angle in degrees as Python scalar
    
    Uses CASA's quantity conversion for consistency with CASA conventions.
    """
    return float(np.atleast_1d(qa.getvalue(qa.convert({'value': xrad, 'unit': 'rad'}, 'deg')))[0])


def compute_parallactic_angle(temp_ms, field_name, time_mjd, field_id=None):
    """
    Compute parallactic angle at specified time using CASA AZ/EL method.
    
    Inputs:
        temp_ms : str - Path to measurement set
        field_name : str - Name of the field/source
        time_mjd : float - Time in MJD seconds
        field_id : int or None - Field ID (auto-determined if None)
    
    Outputs:
        chi_deg : float - Parallactic angle at specified time in degrees
        diagnostics : dict - Dictionary with AZ, EL, LAT, time, and observatory name
    
    Computes parallactic angle using: χ = atan2(-sin(A), tan(φ)cos(e) - cos(A)sin(e))
    where A=azimuth, e=elevation, φ=site latitude (CASA convention: A=0 at N, +E).
    """
    # Resolve field ID
    msmd.open(temp_ms)
    try:
        fids = msmd.fieldsforname(field_name)
        field_id = int(fids[0]) if len(fids) else None
    except Exception:
        field_id = None

    if field_id is None:
        tb.open(temp_ms)
        field_ids_all = tb.getcol('FIELD_ID')
        tb.close()
        vals, cnts = np.unique(field_ids_all, return_counts=True)
        field_id = int(vals[np.argmax(cnts)])
        print(f'WARNING: field "{field_name}" not found; using modal FIELD_ID {field_id}')
    
    phase_dir = msmd.phasecenter(field_id)
    msmd.close()

    tb.open(temp_ms + '/OBSERVATION')
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

    # Compute parallactic angle: χ = atan2(-sin(A), tan(φ)cos(e) - cos(A)sin(e))
    num = -np.sin(az_rad)
    den = np.tan(lat_rad)*np.cos(el_rad) - np.cos(az_rad)*np.sin(el_rad)
    chi_rad = np.arctan2(num, den)
    chi_deg = rad_to_deg_scalar(chi_rad)
    chi_deg = ((chi_deg + 180.0) % 360.0) - 180.0  # Wrap to (-180, 180]

    az_deg = rad_to_deg_scalar(az_rad)
    el_deg = rad_to_deg_scalar(el_rad)
    lt_deg = rad_to_deg_scalar(lat_rad)
    
    diagnostics = {
        'az_deg': az_deg, 
        'el_deg': el_deg, 
        'lat_deg': lt_deg,
        'time': time_mjd,
        'observatory': obsname
    }
    
    return chi_deg, diagnostics


def correct_crosshand_phase(u_prime, v_prime, rho):
    """
    Apply cross-hand phase correction using inverse rotation.
    
    Inputs:
        u_prime : float/array - Uncorrected Stokes U
        v_prime : float/array - Uncorrected Stokes V
        rho : float/array - Cross-hand phase in radians
    
    Outputs:
        u_corrected : float/array - Corrected Stokes U
        v_corrected : float/array - Corrected Stokes V
    
    Applies rotation matrix to remove instrumental cross-hand phase:
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
    
    Inputs:
        q : float/array - Stokes Q in feed frame
        u : float/array - Stokes U in feed frame
        parang_deg : float - Parallactic angle in degrees
    
    Outputs:
        q_corrected : float/array - Stokes Q in sky frame
        u_corrected : float/array - Stokes U in sky frame
    
    Rotates Q,U from feed frame to sky frame using inverse rotation by -2χ:
        Q_sky = Q_feed·cos(-2χ) + U_feed·sin(-2χ)
        U_sky = -Q_feed·sin(-2χ) + U_feed·cos(-2χ)
    """
    parang_rad = np.radians(-2 * parang_deg)
    
    q_corrected = q * np.cos(parang_rad) + u * np.sin(parang_rad)
    u_corrected = -q * np.sin(parang_rad) + u * np.cos(parang_rad)
    
    return q_corrected, u_corrected


def calculate_derotated_angle(q, u, rm, lambda_sq_val):
    """
    Calculate de-rotated polarization angle after RM correction.
    
    Inputs:
        q : float/array - Stokes Q
        u : float/array - Stokes U
        rm : float - Rotation Measure in rad/m²
        lambda_sq_val : float - Wavelength squared in m²
    
    Outputs:
        pol_angle_deg : float/array - Polarization angle in degrees, range (-90, 90]
    
    Applies Faraday rotation correction and computes intrinsic polarization angle:
        Rotation: θ = 2·RM·λ²
        Q' = Q·cos(θ) + U·sin(θ), U' = -Q·sin(θ) + U·cos(θ)
        Angle = 0.5·arctan2(U', Q')
    """
    rot_angle = 2 * rm * lambda_sq_val
    cos_rot = np.cos(rot_angle)
    sin_rot = np.sin(rot_angle)
    
    q_derot = q * cos_rot + u * sin_rot
    u_derot = -q * sin_rot + u * cos_rot
    
    pol_angle_rad = 0.5 * np.arctan2(u_derot, q_derot)
    pol_angle_deg = np.degrees(pol_angle_rad)
    
    # Wrap to (-90, 90] for polarization angle
    pol_angle_deg = ((pol_angle_deg + 90.0) % 180.0) - 90.0
    
    return pol_angle_deg


def angle_difference(angle1, angle2):
    """
    Calculate smallest angular difference between two angles.
    
    Inputs:
        angle1 : float - First angle in degrees
        angle2 : float - Second angle in degrees
    
    Outputs:
        abs_diff : float - Absolute angular difference in degrees
    
    Handles 2π wrapping to find the smallest angular separation,
    accounting for angles wrapping around ±180°.
    """
    diff = angle1 - angle2
    diff = ((diff + 180.0) % 360.0) - 180.0
    return abs(diff)


def plot_stokes_spectra(freq, I_scans, Q_scans, U_scans, V_scans, scan_numbers, output_dir, filename='stokes_spectra.png', title='Stokes Parameters', zoom_percentile=None, zoom_stokes=None):
    """
    Generate per-channel spectra for all Stokes parameters with multiple scans.
    
    Inputs:
        freq : array - Frequency array in GHz (nchan,)
        I_scans, Q_scans, U_scans, V_scans : arrays - Stokes flux densities (n_scans, nchan)
        scan_numbers : array - Scan number labels (n_scans,)
        output_dir : str - Directory path for saving plot
        filename : str - Output filename
        title : str - Plot title
        zoom_percentile : float or None - If provided, zoom to inner percentile range (e.g., 90 = 5th to 95th percentile)
        zoom_stokes : list or None - List of Stokes parameters to zoom (e.g., ['V'] to zoom only Stokes V). If None, zoom all.
    
    Outputs:
        None (saves plot to disk)
    
    Creates 4-panel plot showing flux density vs frequency for I, Q, U, V.
    Each scan is plotted with a different color. For single scan (n_scans=1), 
    plots a single trace.
    """
    n_scans = I_scans.shape[0]
    
    if n_scans == 1:
        colors = ['C0']
    else:
        colors = plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))
    
    fig, ax = plt.subplots(4, 1, figsize=(12, 14), sharex=True, constrained_layout=True)
    
    for scan_idx in range(n_scans):
        color = colors[scan_idx % len(colors)]
        label = f'Scan {scan_numbers[scan_idx]}' if n_scans > 1 else None
        
        ax[0].plot(freq, I_scans[scan_idx], marker='.', linestyle='-', alpha=0.6, 
                   color=color, label=label, markersize=3, linewidth=0.5)
        ax[1].plot(freq, Q_scans[scan_idx], marker='.', linestyle='-', alpha=0.6, 
                   color=color, label=label, markersize=3, linewidth=0.5)
        ax[2].plot(freq, U_scans[scan_idx], marker='.', linestyle='-', alpha=0.6, 
                   color=color, label=label, markersize=3, linewidth=0.5)
        ax[3].plot(freq, V_scans[scan_idx], marker='.', linestyle='-', alpha=0.6, 
                   color=color, label=label, markersize=3, linewidth=0.5)
    
    if zoom_percentile is not None:
        lower_p = (100 - zoom_percentile) / 2
        upper_p = 100 - lower_p
        
        stokes_names = ['I', 'Q', 'U', 'V']
        for idx, (data, stokes_name) in enumerate(zip([I_scans, Q_scans, U_scans, V_scans], stokes_names)):
            # Only zoom if zoom_stokes is None (zoom all) or if this Stokes parameter is in the list
            if zoom_stokes is None or stokes_name in zoom_stokes:
                valid_data = data[np.isfinite(data)]
                if len(valid_data) > 0:
                    y_min = np.percentile(valid_data, lower_p)
                    y_max = np.percentile(valid_data, upper_p)
                    # Add 10% padding
                    y_range = y_max - y_min
                    ax[idx].set_ylim(y_min - 0.1*y_range, y_max + 0.1*y_range)
    
    ax[0].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[0].set_ylabel('Stokes I [Jy]')
    ax[0].set_title(f'{title} - I')
    ax[0].grid(True, alpha=0.3)
    if n_scans > 1:
        ax[0].legend(fontsize=8, ncol=min(3, n_scans))
    
    ax[1].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[1].set_ylabel('Stokes Q [Jy]')
    ax[1].set_title(f'{title} - Q')
    ax[1].grid(True, alpha=0.3)
    
    ax[2].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[2].set_ylabel('Stokes U [Jy]')
    ax[2].set_title(f'{title} - U')
    ax[2].grid(True, alpha=0.3)
    
    ax[3].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[3].set_ylabel('Stokes V [Jy]')
    ax[3].set_xlabel('Frequency [GHz]')
    ax[3].set_title(f'{title} - V')
    ax[3].grid(True, alpha=0.3)
    
    plt.savefig(os.path.join(output_dir, filename), dpi=150)
    plt.close(fig)


def plot_polang_deviation(freq_ghz, Q_scans, U_scans, V_scans, rho_rad_scans, chi_deg_scans, 
                          selected_solution_scans, scan_numbers, global_rm, target_polang, output_dir):
    """
    Plot deviation of de-rotated polarization angle from target (one plot per scan).
    
    Shows BOTH ±π options for each channel - which allows visualization of which
    option (positive or negative ρ) deviates less from the target. The option with
    smaller deviation (closer to zero) is the one selected in the calibration.
    
    Inputs:
        freq_ghz : array - Frequency in GHz
        Q_scans, U_scans, V_scans : arrays - Stokes parameters per scan
        rho_rad_scans : array - Raw cross-hand phase per scan
        chi_deg_scans : array - Parallactic angles per scan
        selected_solution_scans : array - Solution labels ('positive'=ρ≥0, 'negative'=ρ<0)
        scan_numbers : array - Scan numbers
        global_rm : float - Global rotation measure
        target_polang : float - Target polarization angle
        output_dir : str - Output directory
    """
    c = 2.99792458e8  # Speed of light
    
    for scan_idx, scan in enumerate(scan_numbers):
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        
        Q_scan = Q_scans[scan_idx]
        U_scan = U_scans[scan_idx]
        V_scan = V_scans[scan_idx]
        rho_rad_scan = rho_rad_scans[scan_idx]
        sol_labels = selected_solution_scans[scan_idx]
        
        # Only channels with solutions (not empty, not outlier)
        has_solution = (sol_labels != '') & (sol_labels != 'outlier')
        solution_idx = np.where(has_solution)[0]
        
        if len(solution_idx) == 0:
            plt.close(fig)
            continue
        
        freq_hz = freq_ghz[solution_idx] * 1e9
        lambda_sq = (c / freq_hz)**2
        
        Q_channels = Q_scan[solution_idx]
        U_channels = U_scan[solution_idx]
        V_channels = V_scan[solution_idx]
        rho_raw = rho_rad_scan[solution_idx]
        sol_names = sol_labels[solution_idx]
        
        # Compute both ±π options
        rho_option1 = ((rho_raw + np.pi) % (2*np.pi)) - np.pi
        rho_option2 = ((rho_raw + np.pi + np.pi) % (2*np.pi)) - np.pi
        
        # Label by sign: determine which option is positive and which is negative
        rho_positive = np.where(rho_option1 >= 0, rho_option1, rho_option2)
        rho_negative = np.where(rho_option1 < 0, rho_option1, rho_option2)
        
        # Apply corrections for both sign options
        U_xh_positive, V_xh_positive = correct_crosshand_phase(U_channels, V_channels, rho_positive)
        U_xh_negative, V_xh_negative = correct_crosshand_phase(U_channels, V_channels, rho_negative)
        
        Q_final_positive, U_final_positive = correct_parallactic_angle(Q_channels, U_xh_positive, chi_deg_scans[scan_idx])
        Q_final_negative, U_final_negative = correct_parallactic_angle(Q_channels, U_xh_negative, chi_deg_scans[scan_idx])
        
        # De-rotate at global RM for both sign options
        angles_positive = calculate_derotated_angle(Q_final_positive, U_final_positive, global_rm, lambda_sq)
        angles_negative = calculate_derotated_angle(Q_final_negative, U_final_negative, global_rm, lambda_sq)
        
        # Compute deviations for BOTH options (not just selected)
        deviation_positive = angles_positive - target_polang
        deviation_positive = np.where(deviation_positive > 90, deviation_positive - 180, deviation_positive)
        deviation_positive = np.where(deviation_positive < -90, deviation_positive + 180, deviation_positive)
        
        deviation_negative = angles_negative - target_polang
        deviation_negative = np.where(deviation_negative > 90, deviation_negative - 180, deviation_negative)
        deviation_negative = np.where(deviation_negative < -90, deviation_negative + 180, deviation_negative)
        
        # Convert rho to degrees for colormap
        rho_positive_deg = np.degrees(rho_positive)
        rho_negative_deg = np.degrees(rho_negative)
        
        # Plot BOTH options for ALL channels (not just selected)
        # Positive option (circles)
        scatter1 = ax.scatter(freq_ghz[solution_idx], deviation_positive, 
                             c=rho_positive_deg, marker='o', s=50, 
                             cmap='viridis', alpha=0.8, edgecolors='black', linewidths=0.5,
                             vmin=-180, vmax=180)
        
        # Negative option (squares)
        scatter2 = ax.scatter(freq_ghz[solution_idx], deviation_negative, 
                             c=rho_negative_deg, marker='s', s=50,
                             cmap='viridis', alpha=0.8, edgecolors='black', linewidths=0.5,
                             vmin=-180, vmax=180)
        
        # Add colorbar (showing full ρ range from -180 to 180)
        cbar = plt.colorbar(scatter2, ax=ax)
        cbar.set_label('ρ [deg]', rotation=270, labelpad=20)
        
        # Add reference line at zero
        ax.axhline(0, color='black', linestyle='--', alpha=0.5, linewidth=1)
        
        # Custom legend for markers
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
                   markeredgecolor='black', markersize=8, label='Positive ρ option (ρ ≥ 0)'),
            Line2D([0], [0], marker='s', color='w', markerfacecolor='gray',
                   markeredgecolor='black', markersize=8, label='Negative ρ option (ρ < 0)')
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
        
        ax.set_xlabel('Frequency [GHz]')
        ax.set_ylabel('Deviation [deg]')
        ax.set_title(f'Scan {scan}: Both ±π Options (RM={global_rm:.2f} rad/m²)')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'polang_deviation_scan{scan}.png'), dpi=150)
        plt.close(fig)


def plot_crosshand_phase(freq, U_scans, V_scans, rho_deg_scans, scan_numbers, min_flux, output_dir):
    """
    Generate cross-hand phase diagnostic plot for multiple scans.
    
    Inputs:
        freq : array - Frequency array in GHz (nchan,)
        U_scans : array - Stokes U flux (n_scans, nchan)
        V_scans : array - Stokes V flux (n_scans, nchan)
        rho_deg_scans : array - Cross-hand phase in degrees (n_scans, nchan)
        scan_numbers : array - Scan number labels (n_scans,)
        min_flux : float - Minimum flux threshold for reference line
        output_dir : str - Directory path for saving plot
    
    Outputs:
        None (saves plot to disk)
    
    Top panel shows cross-hand flux amplitude, bottom panel shows phase.
    Each scan plotted in different color. For single scan, plots single trace.
    """
    n_scans = U_scans.shape[0]
    
    if n_scans == 1:
        colors = ['red']
    else:
        colors = plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))
    
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    for scan_idx in range(n_scans):
        color = colors[scan_idx % len(colors)]
        label = f'Scan {scan_numbers[scan_idx]}' if n_scans > 1 else None
        
        U = U_scans[scan_idx]
        V = V_scans[scan_idx]
        rho_deg = rho_deg_scans[scan_idx]
        
        valid = np.isfinite(U) & np.isfinite(V) & np.isfinite(freq) & np.isfinite(rho_deg)
        
        # Top: cross-hand flux
        crosshand_flux = np.sqrt(U**2 + V**2)
        ax_top.scatter(freq[valid], crosshand_flux[valid], s=8, c=[color], alpha=0.6, label=label)
        
        # Bottom: cross-hand phase
        ax_bottom.scatter(freq[valid], rho_deg[valid], s=8, c=[color], alpha=0.6, label=label)
    
    ax_top.axhline(min_flux, color='gray', linestyle='--', alpha=0.7, 
                   label=f'{min_flux} Jy reference')
    ax_top.set_ylabel('Cross-hand flux [Jy] (√(U² + V²))')
    ax_top.set_title('Per-channel cross-hand flux (before correction)')
    ax_top.grid(True, alpha=0.3)
    ax_top.legend(fontsize=8, ncol=min(3, n_scans+1))
    
    ax_bottom.set_xlabel('Frequency [GHz]')
    ax_bottom.set_ylabel('Cross-hand phase ρ [deg]')
    ax_bottom.set_title('Per-channel cross-hand phase (arctan2(-V, U))')
    ax_bottom.grid(True, alpha=0.3)
    if n_scans > 1:
        ax_bottom.legend(fontsize=8, ncol=min(3, n_scans))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'crosshand_phase.png'), dpi=150)
    plt.close(fig)


def plot_rm_sweep(rm_trials, angles_per_rm_positive, angles_per_rm_negative, target_angle, output_dir, title_suffix=''):
    """
    Plot polarization angle vs RM for both sign options (positive and negative ρ).
    
    Inputs:
        rm_trials : array - RM values trialed
        angles_per_rm_positive : array - Angles for positive ρ at each RM
        angles_per_rm_negative : array - Angles for negative ρ at each RM
        target_angle : float - Target polarization angle
        output_dir : str - Output directory
        title_suffix : str - Additional text for title
    
    Outputs:
        None (saves plot to disk)
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    ax.plot(rm_trials, angles_per_rm_positive, 'o-', label='Positive ρ', alpha=0.7)
    ax.plot(rm_trials, angles_per_rm_negative, 's-', label='Negative ρ', alpha=0.7)
    ax.axhline(target_angle, color='red', linestyle='--', label=f'Target = {target_angle}°')
    ax.set_xlabel('RM [rad/m²]')
    ax.set_ylabel('Polarization Angle [deg]')
    ax.set_title(f'RM Sweep {title_suffix}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'rm_sweep{title_suffix.replace(" ", "_")}.png'), dpi=150)
    plt.close(fig)


def plot_delta_rm_diagnostics(all_delta_rms, all_min_deviations, all_best_angles, filtered_delta_rms, 
                                filtered_deviations, filtered_angles, global_delta_rm, scan_labels, 
                                filtered_scan_labels, target_angle, output_dir):
    """
    Generate 3-panel diagnostic plot for delta_RM analysis.
    
    Inputs:
        all_delta_rms : array - All delta_RM values from Stage 1
        all_min_deviations : array - Corresponding minimum deviations
        all_best_angles : array - Actual de-rotated angles at best delta_RM
        filtered_delta_rms : array - Good delta_RM values after filtering
        filtered_deviations : array - Corresponding deviations
        filtered_angles : array - Corresponding actual angles
        global_delta_rm : float - Selected global delta_RM
        scan_labels : array - Scan labels for all points
        filtered_scan_labels : array - Scan labels for filtered points
        target_angle : float - Target polarization angle
        output_dir : str - Output directory
    
    Outputs:
        None (saves plot to disk)
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12))
    
    # Panel 1: All delta_RM trials
    unique_scans = np.unique(scan_labels)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_scans)))
    
    for i, scan in enumerate(unique_scans):
        mask = scan_labels == scan
        ax1.scatter(all_delta_rms[mask], all_min_deviations[mask], 
                   c=[colors[i]], label=f'Scan {scan}', alpha=0.5, s=10)
    
    ax1.axvline(global_delta_rm, color='red', linestyle='--', linewidth=2, 
               label=f'Global δRM = {global_delta_rm:.2f}')
    ax1.set_xlabel('Delta RM [rad/m²]')
    ax1.set_ylabel('Min Deviation from Target [deg]')
    ax1.set_title('Stage 1: All Delta_RM Trials')
    ax1.legend(fontsize=8, ncol=2)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Filtered (good) delta_RM trials
    for i, scan in enumerate(unique_scans):
        if len(filtered_scan_labels) > 0:
            mask = filtered_scan_labels == scan
            if np.any(mask):
                ax2.scatter(filtered_delta_rms[mask], filtered_deviations[mask], 
                           c=[colors[i]], label=f'Scan {scan}', alpha=0.6, s=15)
    
    ax2.axvline(global_delta_rm, color='red', linestyle='--', linewidth=2, 
               label=f'Global δRM = {global_delta_rm:.2f}')
    ax2.set_xlabel('Delta RM [rad/m²]')
    ax2.set_ylabel('Min Deviation from Target [deg]')
    ax2.set_title('Stage 2: Filtered (Good Quality) Delta_RM Trials')
    ax2.legend(fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Scatter plot of actual de-rotated angles for good values (colored by scan)
    for i, scan in enumerate(unique_scans):
        if len(filtered_scan_labels) > 0:
            mask = filtered_scan_labels == scan
            if np.any(mask):
                ax3.scatter(filtered_delta_rms[mask], filtered_angles[mask], 
                           c=[colors[i]], label=f'Scan {scan}', alpha=0.6, s=15)
    
    ax3.axvline(global_delta_rm, color='red', linestyle='--', linewidth=2, 
               label=f'Median δRM = {global_delta_rm:.2f}')
    ax3.axhline(target_angle, color='green', linestyle=':', linewidth=2, alpha=0.7,
               label=f'Target = {target_angle}°')
    ax3.set_xlabel('Delta RM [rad/m²]')
    ax3.set_ylabel('De-rotated Polarization Angle [deg]')
    ax3.set_title('Stage 2: De-rotated Angles at Good Delta_RM Values')
    ax3.legend(fontsize=8, ncol=2)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'delta_rm_diagnostics.png'), dpi=150)
    plt.close(fig)


def plot_final_rho_values(freq_ghz, U_scans, V_scans, rho_raw_scans, rho_final_scans, 
                          freq_averaged_data_scans, selected_solution_scans, scan_numbers, output_dir, min_flux):
    """
    Generate 3-panel plot of final rho values with flagging information.
    
    Inputs:
        freq_ghz : array - Frequency in GHz
        U_scans : array - Stokes U per scan (n_scans, nchan)
        V_scans : array - Stokes V per scan (n_scans, nchan)
        rho_raw_scans : array - Raw rho values per scan in radians (n_scans, nchan)
        rho_final_scans : array - Final selected rho values in radians (n_scans, nchan)
        freq_averaged_data_scans : list - Frequency averaged data per scan (or None)
        selected_solution_scans : array - Solution labels per scan (n_scans, nchan)
        scan_numbers : array - Scan number labels
        output_dir : str - Output directory
        min_flux : float - Minimum flux threshold for reference
    
    Outputs:
        None (saves plot to disk)
    """
    n_scans = len(scan_numbers)
    colors = plt.cm.tab10(np.linspace(0, 1, min(n_scans, 10)))
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12), sharex=True)
    
    # Panel 1: Cross-hand flux with flagging
    for scan_idx in range(n_scans):
        color = colors[scan_idx % len(colors)]
        label = f'Scan {scan_numbers[scan_idx]}'
        
        crosshand_flux = np.sqrt(U_scans[scan_idx]**2 + V_scans[scan_idx]**2)
        sol_labels = selected_solution_scans[scan_idx]
        
        # Good data
        good_mask = (sol_labels != 'outlier') & (sol_labels != '') & np.isfinite(crosshand_flux)
        if np.any(good_mask):
            ax1.plot(freq_ghz[good_mask], crosshand_flux[good_mask], '.', color=color, alpha=0.6, 
                    label=label, markersize=4)
        
        # Flagged data
        flagged_mask = (sol_labels == 'outlier') & np.isfinite(crosshand_flux)
        if np.any(flagged_mask):
            ax1.plot(freq_ghz[flagged_mask], crosshand_flux[flagged_mask], 'x', color='red', 
                    alpha=0.6, markersize=4)
    
    # Add threshold line
    ax1.axhline(min_flux, color='gray', linestyle='--', alpha=0.7, linewidth=1.5,
                label=f'Threshold ({min_flux:.3f} Jy)')
    
    ax1.set_ylabel('Cross-hand Flux [Jy]')
    ax1.set_title('Cross-hand Flux √(U² + V²)')
    ax1.legend(fontsize=8, ncol=min(3, n_scans))
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Raw cross-hand phase with flagging
    for scan_idx in range(n_scans):
        color = colors[scan_idx % len(colors)]
        rho_deg = np.degrees(rho_raw_scans[scan_idx])
        sol_labels = selected_solution_scans[scan_idx]
        
        # Good data
        good_mask = (sol_labels != 'outlier') & (sol_labels != '') & np.isfinite(rho_deg)
        if np.any(good_mask):
            ax2.plot(freq_ghz[good_mask], rho_deg[good_mask], '.', color=color, alpha=0.6, markersize=4)
        
        # Flagged data
        flagged_mask = (sol_labels == 'outlier') & np.isfinite(rho_deg)
        if np.any(flagged_mask):
            ax2.plot(freq_ghz[flagged_mask], rho_deg[flagged_mask], 'x', color='red', 
                    alpha=0.6, markersize=4)
    
    ax2.set_ylabel('Raw ρ [deg]')
    ax2.set_title('Raw Cross-hand Phase: ρ = arctan2(-V, U)')
    ax2.grid(True, alpha=0.3)
    
    # Panel 3: Final cross-hand phase with flagging
    # Check if scan-averaged mode (all freq_averaged_data are same object)
    is_scan_averaged = False
    if len(freq_averaged_data_scans) > 1 and freq_averaged_data_scans[0] is not None:
        is_scan_averaged = all(freq_averaged_data_scans[i] is freq_averaged_data_scans[0] 
                               for i in range(1, len(freq_averaged_data_scans)))
    
    # Plot per-channel data for each scan
    for scan_idx in range(n_scans):
        color = colors[scan_idx % len(colors)]
        rho_final_deg = np.degrees(rho_final_scans[scan_idx])
        sol_labels = selected_solution_scans[scan_idx]
        
        # Unflagged (good) points
        good_mask = (sol_labels != 'outlier') & (sol_labels != '') & np.isfinite(rho_final_deg)
        if np.any(good_mask):
            ax3.scatter(freq_ghz[good_mask], rho_final_deg[good_mask], 
                       c=[color], s=20, alpha=0.7, marker='o', edgecolors='black', linewidths=0.5)
        
        # Flagged (outlier) points
        outlier_mask = (sol_labels == 'outlier') & np.isfinite(rho_final_deg)
        if np.any(outlier_mask):
            ax3.scatter(freq_ghz[outlier_mask], rho_final_deg[outlier_mask], 
                       c='red', s=30, alpha=0.8, marker='x', linewidths=1.5)
    
    # Plot frequency-averaged solution
    if is_scan_averaged and freq_averaged_data_scans[0] is not None:
        # SCAN-AVERAGED MODE: Plot once in BLACK
        avg_data = freq_averaged_data_scans[0]
        avg_freq = avg_data['freq']
        avg_rho_deg = avg_data['rho_deg']
        valid_avg = np.isfinite(avg_rho_deg)
        
        # Plot complete line
        ax3.plot(avg_freq[valid_avg], avg_rho_deg[valid_avg], 
                linestyle='-', color='black', linewidth=2, alpha=0.9)
        
        if 'extrapolated' in avg_data:
            extrap_bins = avg_data['extrapolated']
            
            # Filled diamonds for good bins
            good_bins = valid_avg & ~extrap_bins
            if np.any(good_bins):
                ax3.scatter(avg_freq[good_bins], avg_rho_deg[good_bins], 
                           marker='D', s=64, color='black', edgecolors='black', 
                           linewidths=0.5, zorder=5)
            
            # Open diamonds for extrapolated bins
            extrap_bins_plot = valid_avg & extrap_bins
            if np.any(extrap_bins_plot):
                ax3.scatter(avg_freq[extrap_bins_plot], avg_rho_deg[extrap_bins_plot],
                           marker='D', s=64, facecolors='none', edgecolors='black',
                           linewidths=2, zorder=5)
        else:
            # No extrapolation info
            ax3.scatter(avg_freq[valid_avg], avg_rho_deg[valid_avg],
                       marker='D', s=64, color='black', edgecolors='black',
                       linewidths=0.5, zorder=5)
    else:
        # PER-SCAN MODE: Plot each scan in its color
        for scan_idx in range(n_scans):
            if freq_averaged_data_scans[scan_idx] is not None:
                color = colors[scan_idx % len(colors)]
                avg_data = freq_averaged_data_scans[scan_idx]
                avg_freq = avg_data['freq']
                avg_rho_deg = avg_data['rho_deg']
                valid_avg = np.isfinite(avg_rho_deg)
                
                # Plot complete line
                ax3.plot(avg_freq[valid_avg], avg_rho_deg[valid_avg], 
                        linestyle='-', color=color, linewidth=2, alpha=0.9)
                
                if 'extrapolated' in avg_data:
                    extrap_bins = avg_data['extrapolated']
                    
                    # Filled diamonds for good bins
                    good_bins = valid_avg & ~extrap_bins
                    if np.any(good_bins):
                        ax3.scatter(avg_freq[good_bins], avg_rho_deg[good_bins], 
                                   marker='D', s=64, color=color, edgecolors='black', 
                                   linewidths=0.5, zorder=5)
                    
                    # Open diamonds for extrapolated bins
                    extrap_bins_plot = valid_avg & extrap_bins
                    if np.any(extrap_bins_plot):
                        ax3.scatter(avg_freq[extrap_bins_plot], avg_rho_deg[extrap_bins_plot],
                                   marker='D', s=64, facecolors='none', edgecolors=color,
                                   linewidths=2, zorder=5)
                else:
                    # No extrapolation info
                    ax3.scatter(avg_freq[valid_avg], avg_rho_deg[valid_avg],
                               marker='D', s=64, color=color, edgecolors='black',
                               linewidths=0.5, zorder=5)
    
    # Add legend elements
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
               markeredgecolor='black', markersize=8, label='Good'),
        Line2D([0], [0], marker='x', color='w', markerfacecolor='red', 
               markeredgecolor='red', markersize=8, label='Flagged', linewidth=1.5),
        Line2D([0], [0], marker='D', color='gray', markerfacecolor='gray',
               markeredgecolor='black', markersize=8, linestyle='-', linewidth=2, 
               label='Freq. Averaged'),
        Line2D([0], [0], marker='D', color='w', markerfacecolor='none',
               markeredgecolor='gray', markersize=8, markeredgewidth=2, 
               label='Extrapolated')
    ]
    ax3.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    ax3.set_xlabel('Frequency [GHz]')
    ax3.set_ylabel('Final ρ [deg]')
    ax3.set_title('Final Cross-hand Phase (±π resolved, outliers flagged)')
    ax3.grid(True, alpha=0.3)
    
    # Adjust panel 3 (final rho) to zoom on unflagged data range
    all_final_rho_deg = []
    for scan_idx in range(n_scans):
        rho_final_deg = np.degrees(rho_final_scans[scan_idx])
        sol_labels = selected_solution_scans[scan_idx]
        good_mask = (sol_labels != 'outlier') & (sol_labels != '') & np.isfinite(rho_final_deg)
        if np.any(good_mask):
            all_final_rho_deg.extend(rho_final_deg[good_mask])
    
    if len(all_final_rho_deg) > 0:
        rho_min = np.min(all_final_rho_deg)
        rho_max = np.max(all_final_rho_deg)
        rho_range = rho_max - rho_min
        margin = rho_range * 0.1 if rho_range > 0 else 10
        ax3.set_ylim(rho_min - margin, rho_max + margin)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'final_rho_values.png'), dpi=150)
    plt.close(fig)


# ============================
# MAIN SCRIPT
# ============================

# ============================
# Configuration and Setup
# ============================
temp_ms = 'pol_ang_temp.ms'
xfdir = GAINPLOTS + '/manualXF'
os.makedirs(xfdir, exist_ok=True)

# ====================================
# ====================================
# We apply all calibration solutions up to and including the D-term (leakage)
# corrections, but NOT the cross-hand phase (XF), since that's what we're solving for.
# This gives us calibrated I,Q,U,V visibilities with only the XF correction missing.

# Define calibration tables (primary calibrator solutions)
ktab0 = GAINTABLES+'/cal_1GC_'+myms+'.K0'
bptab0 = GAINTABLES+'/cal_1GC_'+myms+'.B0'
gptab0 = GAINTABLES+'/cal_1GC_'+myms+'.Gp0'
gatab0 = GAINTABLES+'/cal_1GC_'+myms+'.Ga0'
ftab0 = GAINTABLES+'/cal_1GC_'+myms+'.F0'
dftab0  = GAINTABLES+'/cal_1GC_'+myms+'.Df0'

ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gptab = GAINTABLES+'/cal_1GC_'+myms+'.Gp'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'

kcross  = GAINTABLES+'/cal_1GC_'+myms+'.KCROSS'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'

applycal(vis = myms,
        field = pacal_name,
        parang = False,
        gainfield = [pacal_name,pacal_name, bpcal_name, pacal_name, bpcal_name],
        gaintable = [ktab,gptab,bptab,ftab,dftab],
        interp = ['linear','linear','linear','linear','linear'],
        flagbackup = False)

# Split out polarization calibrator to temporary MS for analysis
mstransform(vis = myms,
        outputvis = temp_ms,
        field = pacal_name,
        usewtspectrum = False,
        datacolumn='corrected')

# ===============================
# Display Analysis Parameters
# ===============================
print("\n=== Cross-Hand Phase (XF) Calibration Parameters ===")
print(f"Polarization calibrator: {pacal_name}")
print(f"Known source polarization angle: {XF_TARGET_POLANG}°")
print(f"Known source rotation measure: {XF_TARGET_RM} rad/m²")
print(f"RM trial search range: {XF_DELTARM_TRIALS} rad/m² (delta from target)")
print(f"Channel averaging interval: {XF_CHANINT} channels")
print(f"Minimum cross-hand flux for valid measurement: {XF_MIN_CROSS_FLUX} Jy")
print(f"Outlier detection sigma threshold: {XF_SIGMA_CLIP}")
print(f"Outlier detection window size: {XF_CLIP_WINDOW} channels")
print(f"Apply smoothing during table creation: {XF_USE_SMOOTHING}")
if XF_USE_SMOOTHING:
    print(f"  Savitzky-Golay filter window: {XF_SAVGOL_WINDOW} (None = auto-calculate)")
    print(f"  Savitzky-Golay polynomial order: {XF_SAVGOL_POLYORDER}")
print(f"Gap filling / extrapolation enabled: {XF_EX}")
if XF_EX:
    print(f"  Extrapolation bandwidth fraction: {XF_EX_FRAC} ({XF_EX_FRAC*100:.0f}% of good bandwidth)")
print(f"Generate scan-averaged solution: {XF_AVG_SCAN}")
print("="*60)

# ============================
# ============================
print("\nExtracting channel frequencies from MS...")
tb.open(temp_ms + '/SPECTRAL_WINDOW')
chan_freq = tb.getcol('CHAN_FREQ')
tb.close()
freq_ghz = (chan_freq / 1.0e9).flatten()
print(f"Spectral window: {freq_ghz.shape[0]} channels")
print(f"Frequency range: {freq_ghz.min():.3f} - {freq_ghz.max():.3f} GHz")

# If not specified, calculate from channel interval and total channels
if XF_MAX_AVG_CHANNELS is None:
    XF_MAX_AVG_CHANNELS = int(len(chan_freq) / XF_CHANINT)
    print(f"Frequency averaging: {XF_MAX_AVG_CHANNELS} channels (auto-calculated from N_CHAN/XF_CHANINT)")
elif not isinstance(XF_MAX_AVG_CHANNELS, int):
    XF_MAX_AVG_CHANNELS = int(len(chan_freq) / XF_CHANINT)
    print(f"Frequency averaging: {XF_MAX_AVG_CHANNELS} channels (defaulted to N_CHAN/XF_CHANINT)")
else:
    print(f"Frequency averaging: {XF_MAX_AVG_CHANNELS} channels (user-specified)")

# ============================
# Load and Average Visibility Data Per Scan
# ============================
# We load data one scan at a time for memory efficiency, averaging over
# baselines and time within each scan to produce per-channel Stokes I,Q,U,V.
# This preserves the scan structure needed for parallactic angle corrections.

print("\nRetrieving scan information from MS...")

msmd.open(temp_ms)
try:
    scan_numbers = msmd.scansforfield(msmd.fieldsforname(pacal_name)[0])
except:
    msmd.close()
    sys.exit(f"ERROR: Could not retrieve scans for field {pacal_name}")

print(f"Found {len(scan_numbers)} scan(s) for {pacal_name}: {scan_numbers}")

nchan = len(freq_ghz) 
print(f"Processing {nchan} frequency channels")

msmd.close()

# Initialize arrays to store scan-averaged visibilities and weights
# Shape: (N_stokes=4, N_scans, N_channels)
n_scans = len(scan_numbers)
vis_avg = np.zeros((4, n_scans, nchan), dtype=complex)
sigma_proxy = np.zeros((4, n_scans, nchan), dtype=float)
scan_times = np.zeros(n_scans, dtype=float)

# Process each scan: average over baselines and time to get per-channel Stokes
print(f"\nAveraging visibilities per scan (collapse baseline/time axes)...")
total_rows = 0
total_unflagged = 0

for scan_idx, scan in enumerate(scan_numbers):
    # Select this scan from the MS
    ms.open(temp_ms)
    ms.selectinit(reset=True)
    ok = ms.msselect({'field': pacal_name, 'scan': str(scan)})
    if not ok:
        ms.close()
        print(f"WARNING: Could not select scan {scan}, skipping")
        vis_avg[:, scan_idx, :] = np.nan
        sigma_proxy[:, scan_idx, :] = np.nan
        continue
    
    # Request Stokes I,Q,U,V data
    ms.selectpolarization(['I','Q','U','V'])
    
    d = ms.getdata(['data', 'flag', 'weight_spectrum', 'weight', 'time'])
    ms.close()
    
    # Verify required data is present
    missing = [k for k in ('data', 'flag') if k not in d]
    if missing:
        print(f"WARNING: Scan {scan} missing {missing}, skipping")
        vis_avg[:, scan_idx, :] = np.nan
        sigma_proxy[:, scan_idx, :] = np.nan
        del d
        continue
    
    data_scan = d['data']
    raw_flag_scan = np.asarray(d['flag'])
    time_scan = d['time']
    ncorr, nchan_scan, nrow_scan = data_scan.shape
    
    # Verify we have 4 Stokes parameters
    if ncorr != 4:
        print(f"WARNING: Scan {scan} has {ncorr} correlations (expected 4), skipping")
        vis_avg[:, scan_idx, :] = np.nan
        sigma_proxy[:, scan_idx, :] = np.nan
        del d, data_scan, raw_flag_scan, time_scan
        continue
    
    have_wspec = 'weight_spectrum' in d
    flag_iquv_scan = build_iquv_flag(raw_flag_scan, nchan_scan, nrow_scan)
    W_scan = build_iquv_weights(d, have_wspec, nchan_scan, nrow_scan)
    
    del d  # Free memory
    
    # Compute weighted average across baselines and time
    # This collapses (4, nchan, nrow) → (4, nchan)
    vis_avg_scan, sigma_scan = compute_weighted_averages(data_scan, W_scan, flag_iquv_scan)
    
    vis_avg[:, scan_idx, :] = vis_avg_scan
    sigma_proxy[:, scan_idx, :] = sigma_scan
    scan_times[scan_idx] = np.median(time_scan)  # Use mid-time of scan
    
    # Track flagging statistics
    W_eff_scan = np.where(flag_iquv_scan, 0.0, W_scan)
    unflagged_scan = int(np.sum(W_eff_scan > 0))
    total_scan = int(W_eff_scan.size)
    total_rows += nrow_scan
    total_unflagged += unflagged_scan
    
    print(f"  Scan {scan:3d}: {nrow_scan:5d} rows, {unflagged_scan:8d}/{total_scan:8d} unflagged ({100.0*unflagged_scan/total_scan:5.1f}%), time={scan_times[scan_idx]:.1f} MJD")
    
    del data_scan, raw_flag_scan, time_scan, flag_iquv_scan, W_scan, W_eff_scan, vis_avg_scan, sigma_scan

print(f"\nScan averaging complete:")
print(f"  Total rows processed: {total_rows:,d}")
print(f"  Total unflagged data points: {total_unflagged:,d}")
print(f"  Output array shape: {vis_avg.shape} (Stokes, scans, channels)")

# For a properly calibrated point source, Stokes should be real-valued
I_scans = np.real(vis_avg[0])  # (n_scans, nchan)
Q_scans = np.real(vis_avg[1])
U_scans = np.real(vis_avg[2])
V_scans = np.real(vis_avg[3])

I_sigma_scans = sigma_proxy[0]
Q_sigma_scans = sigma_proxy[1]
U_sigma_scans = sigma_proxy[2]
V_sigma_scans = sigma_proxy[3]

print(f"\nExtracted per-scan Stokes parameters: {I_scans.shape} (scans, channels)")

# ============================
# Data Quality Checks
# ============================
# Verify that visibilities are well-calibrated by checking:
# 1. Real vs imaginary components (point source should be predominantly real)
# 2. Statistical properties (median, scatter)

print("\n=== Data Quality Assessment ===")

for scan_idx, scan in enumerate(scan_numbers):
    print(f"\n--- Scan {scan} ---")
    
    I_avg_scan = vis_avg[0, scan_idx, :]
    Q_avg_scan = vis_avg[1, scan_idx, :]
    U_avg_scan = vis_avg[2, scan_idx, :]
    V_avg_scan = vis_avg[3, scan_idx, :]
    
    # Display statistics for real and imaginary components
    quick_stats(f'Scan {scan} Re(I_avg)', np.real(I_avg_scan))
    quick_stats(f'Scan {scan} Im(I_avg)', np.imag(I_avg_scan))
    quick_stats(f'Scan {scan} Re(Q_avg)', np.real(Q_avg_scan))
    quick_stats(f'Scan {scan} Im(Q_avg)', np.imag(Q_avg_scan))
    quick_stats(f'Scan {scan} Re(U_avg)', np.real(U_avg_scan))
    quick_stats(f'Scan {scan} Im(U_avg)', np.imag(U_avg_scan))
    quick_stats(f'Scan {scan} Re(V_avg)', np.real(V_avg_scan))
    quick_stats(f'Scan {scan} Im(V_avg)', np.imag(V_avg_scan))
    
    # Calculate fraction of flux in imaginary component
    # For a well-calibrated point source, this should be << 10%
    im_fracs = {
        'I': calc_im_fraction(I_avg_scan),
        'Q': calc_im_fraction(Q_avg_scan),
        'U': calc_im_fraction(U_avg_scan),
        'V': calc_im_fraction(V_avg_scan)
    }
    
    print(f'Scan {scan} imaginary flux fractions (expect <<10% for point source):')
    for stokes, frac in im_fracs.items():
        status = "WARNING: HIGH" if frac > 10.0 else "OK"
        print(f'  {stokes}: {frac:.1f}% {status}')
    
    # Warn if any Stokes parameter has high imaginary component
    high_im_stokes = [f'{s}({f:.1f}%)' for s, f in im_fracs.items() if f > 10.0]
    if high_im_stokes:
        print(f'WARNING: Scan {scan} has high imaginary flux in {", ".join(high_im_stokes)}')

# ============================
# Generate Pre-Correction Diagnostic Plots
# ============================
print("\nGenerating pre-correction Stokes spectra plots...")

plot_stokes_spectra(freq_ghz, I_scans, Q_scans, U_scans, V_scans, 
                               scan_numbers, xfdir, 
                               filename='stokes_spectra_preXF.png',
                               title='Pre-XF Stokes Parameters')

print(f"Saved: {xfdir}/stokes_spectra_preXF.png")

# Save per-scan data arrays for later analysis
np.savez(os.path.join(xfdir, 'stokes_perscan.npz'),
         freq_ghz=freq_ghz,
         scan_numbers=scan_numbers,
         scan_times=scan_times,
         I_scans=I_scans,
         Q_scans=Q_scans,
         U_scans=U_scans,
         V_scans=V_scans,
         I_sigma_scans=I_sigma_scans,
         Q_sigma_scans=Q_sigma_scans,
         U_sigma_scans=U_sigma_scans,
         V_sigma_scans=V_sigma_scans)

print(f"Saved: {xfdir}/stokes_perscan.npz")

# ============================
# Measure Raw Cross-Hand Phase
# ============================
# The cross-hand phase ρ is the phase offset between XY and YX correlations.
# For a linearly polarized source, it can be measured from Stokes U and V:
#     ρ = arctan2(-V_obs, U_obs)
# This measurement has a ±π ambiguity that must be resolved using RM analysis.

print("\n=== Measuring Raw Cross-Hand Phase: ρ = arctan2(-V, U) ===")

# Calculate cross-hand phase per scan
rho_rad_scans = np.full((n_scans, nchan), np.nan, dtype=float)  # (n_scans, nchan)
rho_deg_scans = np.full((n_scans, nchan), np.nan, dtype=float)

for scan_idx, scan in enumerate(scan_numbers):
    U_scan = U_scans[scan_idx]
    V_scan = V_scans[scan_idx]
    
    valid = np.isfinite(U_scan) & np.isfinite(V_scan)
    num_valid = int(np.sum(valid))
    
    rho_rad = np.full_like(U_scan, np.nan, dtype=float)
    rho_rad[valid] = np.arctan2(-V_scan[valid], U_scan[valid])
    rho_deg = np.degrees(rho_rad)
    rho_deg_wrapped = ((rho_deg + 180.0) % 360.0) - 180.0
    
    rho_rad_scans[scan_idx] = rho_rad
    rho_deg_scans[scan_idx] = rho_deg_wrapped
    
    if num_valid > 0:
        med_rho = np.nanmedian(rho_deg_wrapped)
        std_rho = np.nanstd(rho_deg_wrapped)
        print(f"  Scan {scan}: {num_valid}/{nchan} valid, median={med_rho:.3f} deg, std={std_rho:.3f} deg")

# Plot multi-scan cross-hand phase
plot_crosshand_phase(freq_ghz, U_scans, V_scans, rho_deg_scans, 
                                scan_numbers, XF_MIN_CROSS_FLUX, xfdir)

print(f"Saved: {xfdir}/crosshand_phase.png")

# Save cross-hand phase measurements
np.savez(os.path.join(xfdir, 'crosshand_phase_perscan.npz'),
         freq_ghz=freq_ghz,
         scan_numbers=scan_numbers,
         rho_rad_scans=rho_rad_scans,
         rho_deg_scans=rho_deg_scans)

print(f"Saved: {xfdir}/crosshand_phase_perscan.npz")

# ============================
# Compute Parallactic Angles
# ============================
# The parallactic angle χ describes the rotation of the antenna feed frame
# relative to the sky frame. It varies with time as the source moves across
# the sky. We compute χ at each scan's mid-time for use in ±π disambiguation.
#
# Parallactic angle is given by: χ = atan2(-sin(A), tan(φ)cos(e) - cos(A)sin(e))
# where A = azimuth, e = elevation, φ = observatory latitude

print("\n=== Computing Parallactic Angles Per Scan ===")

chi_deg_scans = np.zeros(n_scans)
parang_diagnostics_scans = []

for scan_idx, scan in enumerate(scan_numbers):
    chi_deg, diagnostics = compute_parallactic_angle(temp_ms, pacal_name, scan_times[scan_idx])
    chi_deg_scans[scan_idx] = chi_deg
    parang_diagnostics_scans.append(diagnostics)
    
    print(f"  Scan {scan}: χ = {chi_deg:+7.3f}° (AZ={diagnostics['az_deg']:6.1f}°, EL={diagnostics['el_deg']:5.1f}°, time={scan_times[scan_idx]:.1f} MJD)")

print(f"\nParallactic angle range: {chi_deg_scans.min():+.3f}° to {chi_deg_scans.max():+.3f}°")
print(f"Parallactic angle span: {chi_deg_scans.max() - chi_deg_scans.min():.3f}°")
print(f"Observatory: {parang_diagnostics_scans[0]['observatory']}")
print(f"Latitude: {parang_diagnostics_scans[0]['lat_deg']:.3f}°")

# Print time differences between consecutive scans
if n_scans > 1:
    print(f"\nTime differences between consecutive scans:")
    for scan_idx in range(1, n_scans):
        time_diff_seconds = scan_times[scan_idx] - scan_times[scan_idx-1]
        time_diff_minutes = time_diff_seconds / 60.0
        print(f"  Scan {scan_numbers[scan_idx-1]} → Scan {scan_numbers[scan_idx]}: {time_diff_minutes:.2f} minutes")
    
    total_time_minutes = (scan_times[-1] - scan_times[0]) / 60.0
    print(f"  Total time span: {total_time_minutes:.2f} minutes")
else:
    print("\nSingle scan observation")

# ============================
# Rotation Measure Trial Analysis (±π Degeneracy Resolution)
# ============================
# The cross-hand phase measurement ρ = arctan2(-V, U) has a ±π ambiguity.
# We resolve this by comparing de-rotated polarization angles to the known
# source angle across a range of rotation measures.
#
# PHYSICS:
# Faraday rotation causes polarization angle to rotate with wavelength:
#     Δχ = RM × λ²
# where RM is the rotation measure in rad/m².
#
# THREE-STAGE APPROACH:
# Stage 1: For each scan/channel, trial both ±π options across a range of RM
#          values, tracking which delta_RM gives the best match to target angle
# Stage 2: Pool all results, filter outliers, find global median delta_RM
# Stage 3: Use global RM to select the correct ±π option per scan/channel
#
# This approach is robust to:
# - Time-variable parallactic angles (per-scan chi values)
# - Uncertain source RM (searches around target RM)
# - Channel-to-channel variations

print("\n=== ±π Degeneracy Resolution via RM Trial Analysis ===")
print("Stage 1: Trial delta_RM for each scan/channel/option")
print("Stage 2: Pool results across scans, find global median delta_RM")
print("Stage 3: Use global RM to select ±π solution per scan/channel")

# Initialize per-scan solution arrays
selected_rho_per_channel_scans = np.full((n_scans, nchan), np.nan)  # (n_scans, nchan)
selected_solution_per_channel_scans = np.full((n_scans, nchan), '', dtype='U8')

delta_rm_trials = np.linspace(XF_DELTARM_TRIALS[0], XF_DELTARM_TRIALS[-1], 30)
rm_trials = delta_rm_trials + XF_TARGET_RM
c = 2.998e8  # Speed of light (m/s)

print(f"\nTrialing {len(rm_trials)} RM values from {rm_trials[0]:.1f} to {rm_trials[-1]:.1f} rad/m²")
print(f"Delta RM search range: {delta_rm_trials[0]:.2f} to {delta_rm_trials[-1]:.2f} rad/m²")

# Initialize global collection arrays for Stage 2 pooling
all_delta_rms = []
all_min_deviations = []
all_best_angles = []  # Actual de-rotated angles at best delta_RM
all_scan_labels = []  # Track which scan each result came from

# ============================
# STAGE 1: Collect best delta_RM for each scan/channel/option
# ============================
print("\n--- Stage 1: Collecting best delta_RM per scan/channel/option ---")

for scan_idx, scan in enumerate(scan_numbers):
    print(f"\nProcessing Scan {scan} (χ = {chi_deg_scans[scan_idx]:.3f}°)")
    
    Q_scan = Q_scans[scan_idx]
    U_scan = U_scans[scan_idx]
    V_scan = V_scans[scan_idx]
    rho_rad_scan = rho_rad_scans[scan_idx]
    
    # Find unflagged channels for this scan
    unflagged_mask = np.isfinite(Q_scan) & np.isfinite(U_scan) & np.isfinite(V_scan) & np.isfinite(rho_rad_scan)
    unflagged_idx = np.where(unflagged_mask)[0]
    
    if len(unflagged_idx) == 0:
        print(f"  WARNING: No unflagged channels for scan {scan}")
        continue
    
    print(f"  Found {len(unflagged_idx)} unflagged channels")
    
    freq_hz = freq_ghz[unflagged_idx] * 1e9
    lambda_sq = (c / freq_hz)**2
    
    Q_channels = Q_scan[unflagged_idx]
    U_channels = U_scan[unflagged_idx]
    V_channels = V_scan[unflagged_idx]
    
    rho_raw = rho_rad_scan[unflagged_idx]
    rho_option1 = ((rho_raw + np.pi) % (2*np.pi)) - np.pi
    rho_option2 = ((rho_raw + np.pi + np.pi) % (2*np.pi)) - np.pi
    
    # Label by sign: determine which option is positive and which is negative
    rho_positive = np.where(rho_option1 >= 0, rho_option1, rho_option2)
    rho_negative = np.where(rho_option1 < 0, rho_option1, rho_option2)
    
    U_xh_positive, V_xh_positive = correct_crosshand_phase(U_channels, V_channels, rho_positive)
    U_xh_negative, V_xh_negative = correct_crosshand_phase(U_channels, V_channels, rho_negative)
    
    Q_final_positive, U_final_positive = correct_parallactic_angle(Q_channels, U_xh_positive, chi_deg_scans[scan_idx])
    Q_final_negative, U_final_negative = correct_parallactic_angle(Q_channels, U_xh_negative, chi_deg_scans[scan_idx])
    
    # Trial RMs for this scan
    n_valid = len(unflagged_idx)
    best_delta_rm_positive = np.zeros(n_valid)
    best_delta_rm_negative = np.zeros(n_valid)
    min_dev_positive = np.full(n_valid, np.inf)
    min_dev_negative = np.full(n_valid, np.inf)
    best_angle_positive = np.zeros(n_valid)
    best_angle_negative = np.zeros(n_valid)
    
    for rm, delta_rm in zip(rm_trials, delta_rm_trials):
        angles_positive = calculate_derotated_angle(Q_final_positive, U_final_positive, rm, lambda_sq)
        angles_negative = calculate_derotated_angle(Q_final_negative, U_final_negative, rm, lambda_sq)
        
        dev_positive = np.abs(angles_positive - XF_TARGET_POLANG)
        dev_positive = np.where(dev_positive > 90, 180 - dev_positive, dev_positive)
        
        dev_negative = np.abs(angles_negative - XF_TARGET_POLANG)
        dev_negative = np.where(dev_negative > 90, 180 - dev_negative, dev_negative)
        
        better_mask_positive = dev_positive < min_dev_positive
        min_dev_positive[better_mask_positive] = dev_positive[better_mask_positive]
        best_delta_rm_positive[better_mask_positive] = delta_rm
        best_angle_positive[better_mask_positive] = angles_positive[better_mask_positive]
        
        better_mask_negative = dev_negative < min_dev_negative
        min_dev_negative[better_mask_negative] = dev_negative[better_mask_negative]
        best_delta_rm_negative[better_mask_negative] = delta_rm
        best_angle_negative[better_mask_negative] = angles_negative[better_mask_negative]
    
    # Add to global collection with scan labels
    all_delta_rms.extend(best_delta_rm_positive)
    all_delta_rms.extend(best_delta_rm_negative)
    all_min_deviations.extend(min_dev_positive)
    all_min_deviations.extend(min_dev_negative)
    all_best_angles.extend(best_angle_positive)
    all_best_angles.extend(best_angle_negative)
    all_scan_labels.extend([scan] * (2 * n_valid))
    
    print(f"  Collected {2*n_valid} delta_RM values (2 options × {n_valid} channels)")

# Convert to arrays
all_delta_rms = np.array(all_delta_rms)
all_min_deviations = np.array(all_min_deviations)
all_best_angles = np.array(all_best_angles)
all_scan_labels = np.array(all_scan_labels)

print(f"\n--- Stage 1 Complete ---")
print(f"Total collected: {len(all_delta_rms)} delta_RM values from {n_scans} scans")
print(f"Delta_RM range: {all_delta_rms.min():.2f} to {all_delta_rms.max():.2f} rad/m²")
print(f"Min deviation range: {all_min_deviations.min():.2f}° to {all_min_deviations.max():.2f}°")

# ============================
# STAGE 2: Quality Filtering and Global Mode Selection
# ============================
# We now have delta_RM estimates from all scans/channels/options.
# common value) which represents the global RM offset from our target.
#
# FILTERING APPROACH:
# - Remove solutions with deviation > median + 3×MAD (outlier rejection)
# - Cap at 75th percentile of deviations (remove worst quartile)
# - Cap at absolute threshold (20°)

print("\n--- Stage 2: Quality filtering and global mode selection ---")

median_dev = np.median(all_min_deviations)
mad_dev = np.median(np.abs(all_min_deviations - median_dev))
threshold_dev = median_dev + 3.0 * 1.4826 * mad_dev

percentile_threshold = np.percentile(all_min_deviations, 75)
absolute_threshold = 20.0

quality_threshold = min(threshold_dev, percentile_threshold, absolute_threshold)

print(f"Deviation statistics:")
print(f"  Median deviation: {median_dev:.2f}°")
print(f"  MAD-based threshold: {threshold_dev:.2f}°")
print(f"  75th percentile: {percentile_threshold:.2f}°")
print(f"  Absolute threshold: {absolute_threshold:.2f}°")
print(f"  Using quality threshold: {quality_threshold:.2f}° (most conservative)")

good_solutions = all_min_deviations < quality_threshold
filtered_delta_rms = all_delta_rms[good_solutions]
filtered_deviations = all_min_deviations[good_solutions]

n_rejected = np.sum(~good_solutions)
n_accepted = np.sum(good_solutions)

print(f"\nQuality filtering results:")
print(f"  Accepted: {n_accepted}/{len(all_delta_rms)} solutions ({n_accepted/len(all_delta_rms)*100:.1f}%)")
print(f"  Rejected: {n_rejected}/{len(all_delta_rms)} solutions ({n_rejected/len(all_delta_rms)*100:.1f}%)")

if n_accepted < 10:
    print(f"WARNING: Very few good solutions ({n_accepted}). Using all solutions instead.")
    filtered_delta_rms = all_delta_rms
    filtered_deviations = all_min_deviations
    quality_threshold = np.inf

global_delta_rm = np.median(filtered_delta_rms)
global_rm = XF_TARGET_RM + global_delta_rm

print(f"\nGlobal median selection:")
print(f"  Median delta_RM: {global_delta_rm:+.2f} rad/m²")
print(f"  Global RM: {global_rm:.2f} rad/m²")

# Additional delta_RM statistics
print(f"\nDetailed delta_RM statistics:")
print(f"  All solutions: mean={np.mean(all_delta_rms):+.3f}, std={np.std(all_delta_rms):.3f}, median={np.median(all_delta_rms):+.3f} rad/m²")
print(f"  Filtered (good) solutions: mean={np.mean(filtered_delta_rms):+.3f}, std={np.std(filtered_delta_rms):.3f}, median={np.median(filtered_delta_rms):+.3f} rad/m²")
print(f"  Delta_RM range (all): [{all_delta_rms.min():+.3f}, {all_delta_rms.max():+.3f}] rad/m²")
print(f"  Delta_RM range (filtered): [{filtered_delta_rms.min():+.3f}, {filtered_delta_rms.max():+.3f}] rad/m²")

# Generate delta_RM diagnostic plot
print("\nGenerating delta_RM diagnostic plot...")

good_solutions = all_min_deviations < quality_threshold
filtered_angles = all_best_angles[good_solutions]

# Match scan labels to filtered data
filtered_scan_labels = all_scan_labels[good_solutions]

plot_delta_rm_diagnostics(all_delta_rms, all_min_deviations, all_best_angles, 
                          filtered_delta_rms, filtered_deviations, filtered_angles,
                          global_delta_rm, all_scan_labels, filtered_scan_labels,
                          XF_TARGET_POLANG, xfdir)
print(f"Saved: {xfdir}/delta_rm_diagnostics.png")

# ============================
# Generate Polarization Angle Histogram at Median Delta_RM
# ============================
# Calculate the de-rotated polarization angles for each scan using the global_rm
# and create a histogram showing the distribution colored by scan
print("\nGenerating polarization angle histogram at median delta_RM...")

polang_at_median_rm = []
polang_scan_labels = []

for scan_idx, scan in enumerate(scan_numbers):
    Q_scan = Q_scans[scan_idx]
    U_scan = U_scans[scan_idx]
    V_scan = V_scans[scan_idx]
    rho_rad_scan = rho_rad_scans[scan_idx]
    
    unflagged_mask = np.isfinite(Q_scan) & np.isfinite(U_scan) & np.isfinite(V_scan) & np.isfinite(rho_rad_scan)
    unflagged_idx = np.where(unflagged_mask)[0]
    
    if len(unflagged_idx) == 0:
        continue
    
    freq_hz = freq_ghz[unflagged_idx] * 1e9
    lambda_sq = (c / freq_hz)**2
    
    Q_channels = Q_scan[unflagged_idx]
    U_channels = U_scan[unflagged_idx]
    V_channels = V_scan[unflagged_idx]
    
    rho_raw = rho_rad_scan[unflagged_idx]
    rho_option1 = ((rho_raw + np.pi) % (2*np.pi)) - np.pi
    rho_option2 = ((rho_raw + np.pi + np.pi) % (2*np.pi)) - np.pi
    
    # Label by sign: determine which option is positive and which is negative
    rho_positive = np.where(rho_option1 >= 0, rho_option1, rho_option2)
    rho_negative = np.where(rho_option1 < 0, rho_option1, rho_option2)
    
    U_xh_positive, V_xh_positive = correct_crosshand_phase(U_channels, V_channels, rho_positive)
    U_xh_negative, V_xh_negative = correct_crosshand_phase(U_channels, V_channels, rho_negative)
    
    Q_final_positive, U_final_positive = correct_parallactic_angle(Q_channels, U_xh_positive, chi_deg_scans[scan_idx])
    Q_final_negative, U_final_negative = correct_parallactic_angle(Q_channels, U_xh_negative, chi_deg_scans[scan_idx])
    
    # Calculate angles at global RM for both sign options
    angles_positive = calculate_derotated_angle(Q_final_positive, U_final_positive, global_rm, lambda_sq)
    angles_negative = calculate_derotated_angle(Q_final_negative, U_final_negative, global_rm, lambda_sq)
    
    # Select the option closer to target
    dev_positive = np.abs(angles_positive - XF_TARGET_POLANG)
    dev_positive = np.where(dev_positive > 90, 180 - dev_positive, dev_positive)
    
    dev_negative = np.abs(angles_negative - XF_TARGET_POLANG)
    dev_negative = np.where(dev_negative > 90, 180 - dev_negative, dev_negative)
    
    select_positive = dev_positive <= dev_negative
    selected_angles = np.where(select_positive, angles_positive, angles_negative)
    
    polang_at_median_rm.extend(selected_angles)
    polang_scan_labels.extend([scan] * len(selected_angles))

polang_at_median_rm = np.array(polang_at_median_rm)
polang_scan_labels = np.array(polang_scan_labels)

fig, ax = plt.subplots(1, 1, figsize=(10, 6))

unique_scans = np.unique(polang_scan_labels)
colors = plt.cm.tab10(np.linspace(0, 1, len(unique_scans)))

# Plot histogram for each scan
for i, scan in enumerate(unique_scans):
    mask = polang_scan_labels == scan
    ax.hist(polang_at_median_rm[mask], bins=50, alpha=0.6, label=f'Scan {scan}', 
            color=colors[i], edgecolor='black', linewidth=0.5)

# Add vertical line for target angle
ax.axvline(XF_TARGET_POLANG, color='red', linestyle='--', linewidth=2, 
          label=f'Target = {XF_TARGET_POLANG}°', zorder=10)

# Calculate and display statistics
mean_angle = np.mean(polang_at_median_rm)
median_angle = np.median(polang_at_median_rm)
std_angle = np.std(polang_at_median_rm)

ax.axvline(mean_angle, color='green', linestyle=':', linewidth=2, 
          label=f'Mean = {mean_angle:.2f}°', alpha=0.7)

ax.set_xlabel('Polarization Angle [deg]')
ax.set_ylabel('Count')
ax.set_title(f'Polarization Angle Distribution at Median RM = {global_rm:.2f} rad/m²\n' + 
            f'(δRM = {global_delta_rm:+.2f} rad/m²)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Add text box with statistics
textstr = f'Mean: {mean_angle:.2f}°\nMedian: {median_angle:.2f}°\nStd: {std_angle:.2f}°\nN: {len(polang_at_median_rm)}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

plt.tight_layout()
plt.savefig(os.path.join(xfdir, 'polang_histogram_at_median_rm.png'), dpi=150)
plt.close(fig)

print(f"Saved: {xfdir}/polang_histogram_at_median_rm.png")
print(f"  Polarization angle statistics at median RM:")
print(f"    Mean: {mean_angle:.2f}°, Median: {median_angle:.2f}°, Std: {std_angle:.2f}°")
print(f"    Deviation from target: {mean_angle - XF_TARGET_POLANG:+.2f}°")

# ============================
# STAGE 3: Use global RM to select ±π solution per scan/channel
# ============================
print("\n--- Stage 3: Selecting ±π solutions using global RM ---")

for scan_idx, scan in enumerate(scan_numbers):
    print(f"\nScan {scan} (χ = {chi_deg_scans[scan_idx]:.3f}°)")
    
    Q_scan = Q_scans[scan_idx]
    U_scan = U_scans[scan_idx]
    V_scan = V_scans[scan_idx]
    rho_rad_scan = rho_rad_scans[scan_idx]
    
    unflagged_mask = np.isfinite(Q_scan) & np.isfinite(U_scan) & np.isfinite(V_scan) & np.isfinite(rho_rad_scan)
    unflagged_idx = np.where(unflagged_mask)[0]
    
    if len(unflagged_idx) == 0:
        continue
    
    freq_hz = freq_ghz[unflagged_idx] * 1e9
    lambda_sq = (c / freq_hz)**2
    
    Q_channels = Q_scan[unflagged_idx]
    U_channels = U_scan[unflagged_idx]
    V_channels = V_scan[unflagged_idx]
    
    rho_raw = rho_rad_scan[unflagged_idx]
    rho_option1 = ((rho_raw + np.pi) % (2*np.pi)) - np.pi
    rho_option2 = ((rho_raw + np.pi + np.pi) % (2*np.pi)) - np.pi
    
    # Label by sign: determine which option is positive and which is negative
    rho_positive = np.where(rho_option1 >= 0, rho_option1, rho_option2)
    rho_negative = np.where(rho_option1 < 0, rho_option1, rho_option2)
    
    U_xh_positive, V_xh_positive = correct_crosshand_phase(U_channels, V_channels, rho_positive)
    U_xh_negative, V_xh_negative = correct_crosshand_phase(U_channels, V_channels, rho_negative)
    
    Q_final_positive, U_final_positive = correct_parallactic_angle(Q_channels, U_xh_positive, chi_deg_scans[scan_idx])
    Q_final_negative, U_final_negative = correct_parallactic_angle(Q_channels, U_xh_negative, chi_deg_scans[scan_idx])
    
    # Calculate angles at global RM for both sign options
    angles_positive_fixed = calculate_derotated_angle(Q_final_positive, U_final_positive, global_rm, lambda_sq)
    angles_negative_fixed = calculate_derotated_angle(Q_final_negative, U_final_negative, global_rm, lambda_sq)
    
    dev_positive_fixed = np.abs(angles_positive_fixed - XF_TARGET_POLANG)
    dev_positive_fixed = np.where(dev_positive_fixed > 90, 180 - dev_positive_fixed, dev_positive_fixed)
    
    dev_negative_fixed = np.abs(angles_negative_fixed - XF_TARGET_POLANG)
    dev_negative_fixed = np.where(dev_negative_fixed > 90, 180 - dev_negative_fixed, dev_negative_fixed)
    
    # Select best option per channel based on deviation from target
    select_positive = dev_positive_fixed <= dev_negative_fixed
    
    selected_rho = np.where(select_positive, rho_positive, rho_negative)
    solution_names = np.where(select_positive, 'positive', 'negative')
    
    for i, ch in enumerate(unflagged_idx):
        selected_rho_per_channel_scans[scan_idx, ch] = selected_rho[i]
        selected_solution_per_channel_scans[scan_idx, ch] = solution_names[i]
    
    n_positive = np.sum(selected_rho >= 0)
    n_negative = np.sum(selected_rho < 0)
    
    print(f"  Selected solutions: {n_positive} positive, {n_negative} negative")

print("\n=== Per-Channel RM Trial Analysis Complete ===")

# ============================
# Plot Polarization Angle Deviation
# ============================
print("\nGenerating polarization angle deviation diagnostic plots (per scan)...")
plot_polang_deviation(freq_ghz, Q_scans, U_scans, V_scans, rho_rad_scans, chi_deg_scans,
                     selected_solution_per_channel_scans, scan_numbers, global_rm, XF_TARGET_POLANG, xfdir)
print(f"Saved: {xfdir}/polang_deviation_scan*.png")

# ============================
# Two-Stage Outlier Detection
# ============================
# Flag channels with unreliable cross-hand phase measurements using:
# 
# STAGE 1: Flux-based flagging
#   - Flag channels where cross-hand flux √(U² + V²) < threshold
#   - These have insufficient SNR for reliable phase measurement
#
# STAGE 2: Local scatter analysis
#   - For each channel, compute local phase scatter in frequency
#   - Flag channels with anomalously high local scatter
#   - Uses sliding window to detect localized RFI or calibration issues

print("\n=== Two-Stage Outlier Detection ===")

for scan_idx, scan in enumerate(scan_numbers):
    print(f"\n--- Scan {scan} Outlier Detection ---")
    
    U_scan = U_scans[scan_idx]
    V_scan = V_scans[scan_idx]
    rho_per_ch = selected_rho_per_channel_scans[scan_idx]
    sol_per_ch = selected_solution_per_channel_scans[scan_idx]
    
    has_solution = ~np.isnan(rho_per_ch) & (sol_per_ch != '')
    n_with_solution = np.sum(has_solution)
    
    if n_with_solution == 0:
        print(f"  No solutions found for scan {scan}")
        continue
    
    # STAGE 1: Low flux flagging
    print(f"  Stage 1: Low flux flagging (< {XF_MIN_CROSS_FLUX} Jy)")
    low_flux_mask = np.sqrt(U_scan**2 + V_scan**2) < XF_MIN_CROSS_FLUX
    low_flux_channels = np.where(low_flux_mask)[0]
    n_low_flux = len(low_flux_channels)
    print(f"    Flagged {n_low_flux} channels")
    
    for ch in low_flux_channels:
        selected_solution_per_channel_scans[scan_idx, ch] = 'outlier'
    
    # STAGE 2: Local scatter flagging
    print(f"  Stage 2: Local scatter flagging")
    
    # Recalculate has_solution after low flux flagging
    has_solution = ~np.isnan(rho_per_ch) & (sol_per_ch != 'outlier') & (sol_per_ch != '')
    n_with_solution = np.sum(has_solution)
    
    if n_with_solution > 10:
        freq_for_outliers = freq_ghz[has_solution]
        rho_deg_for_outliers = np.degrees(rho_per_ch[has_solution])
        channels_for_outliers = np.where(has_solution)[0]
        
        freq_sort_idx = np.argsort(freq_for_outliers)
        sorted_rho_deg = rho_deg_for_outliers[freq_sort_idx]
        
        local_scatter_outliers = np.zeros(len(rho_deg_for_outliers), dtype=bool)
        local_stdevs = np.full(len(sorted_rho_deg), np.nan)
        
        for i in range(len(sorted_rho_deg)):
            start_idx = max(0, i - XF_CLIP_WINDOW // 2)
            end_idx = min(len(sorted_rho_deg), i + XF_CLIP_WINDOW // 2 + 1)
            
            if end_idx - start_idx >= 5:
                window_data = sorted_rho_deg[start_idx:end_idx]
                window_complex = np.exp(1j * np.radians(window_data))
                window_centroid = np.mean(window_complex)
                window_centroid = window_centroid / np.abs(window_centroid)
                angular_deviations = np.abs(window_complex - window_centroid)
                local_stdevs[i] = np.std(angular_deviations)
        
        valid_stdevs = local_stdevs[~np.isnan(local_stdevs)]
        
        if len(valid_stdevs) > 5:
            stdev_mean = np.mean(valid_stdevs)
            stdev_std = np.std(valid_stdevs)
            stdev_median = np.median(valid_stdevs)
            stdev_mad = np.median(np.abs(valid_stdevs - stdev_median))
            
            sigma_threshold = stdev_mean + XF_SIGMA_CLIP * stdev_std
            mad_threshold = stdev_median + XF_SIGMA_CLIP * 1.4826 * stdev_mad
            threshold = max(sigma_threshold, mad_threshold)
            
            print(f"    Local stdev threshold: {threshold:.4f}")
            
            n_flagged = 0
            for i in range(len(sorted_rho_deg)):
                if not np.isnan(local_stdevs[i]) and local_stdevs[i] > threshold:
                    original_idx = freq_sort_idx[i]
                    local_scatter_outliers[original_idx] = True
                    n_flagged += 1
            
            print(f"    Flagged {n_flagged} channels ({n_flagged/len(valid_stdevs)*100:.1f}%)")
            
            outlier_channels = channels_for_outliers[local_scatter_outliers]
            for ch in outlier_channels:
                if selected_solution_per_channel_scans[scan_idx, ch] != 'outlier':
                    selected_solution_per_channel_scans[scan_idx, ch] = 'outlier'
        else:
            print(f"    Insufficient data for threshold calculation")
    else:
        print(f"    Insufficient channels ({n_with_solution}) for local scatter analysis")

print("\n=== Two-Stage Outlier Detection Complete ===")

# ============================
# Frequency Averaging with Linear Extrapolation
# ============================

print("\n=== Frequency Averaging ===")

freq_averaged_data_scans = []

# Define bins across full frequency range
total_channels = len(freq_ghz)
freq_min = freq_ghz.min()
freq_max = freq_ghz.max()

if XF_AVG_SCAN:
    # ============================
    # SCAN-AVERAGED MODE: Pool all scans into single binned solution
    # ============================
    print("Scan-averaged mode: Pooling data from all scans")
    
    # Collect all good solutions from all scans
    all_good_freqs = []
    all_good_rhos = []
    
    for scan_idx, scan in enumerate(scan_numbers):
        rho_per_ch = selected_rho_per_channel_scans[scan_idx]
        sol_per_ch = selected_solution_per_channel_scans[scan_idx]
        
        has_solution = ~np.isnan(rho_per_ch) & (sol_per_ch != 'outlier') & (sol_per_ch != '')
        
        if np.sum(has_solution) > 0:
            all_good_freqs.append(freq_ghz[has_solution])
            all_good_rhos.append(rho_per_ch[has_solution])
    
    if len(all_good_freqs) == 0:
        print("ERROR: No good solutions found across any scan!")
        sys.exit(1)
    
    # Concatenate all scans
    all_good_freqs = np.concatenate(all_good_freqs)
    all_good_rhos = np.concatenate(all_good_rhos)
    
    n_scans = len(scan_numbers)
    max_data_points = total_channels * n_scans
    print(f"  Total good data points: {len(all_good_freqs)} / {max_data_points}")
    
    # Create evenly-spaced bins
    n_freq_bins = XF_MAX_AVG_CHANNELS
    freq_bin_edges = np.linspace(freq_min, freq_max, n_freq_bins + 1)
    freq_bin_centers = 0.5 * (freq_bin_edges[:-1] + freq_bin_edges[1:])
    
    averaged_rho_rad = np.full(n_freq_bins, np.nan)
    bin_channel_counts = np.zeros(n_freq_bins, dtype=int)
    
    # Fill bins with pooled data
    for i in range(n_freq_bins):
        in_bin = (all_good_freqs >= freq_bin_edges[i]) & (all_good_freqs < freq_bin_edges[i+1])
        if i == n_freq_bins - 1:
            in_bin = (all_good_freqs >= freq_bin_edges[i]) & (all_good_freqs <= freq_bin_edges[i+1])
        
        if np.sum(in_bin) > 0:
            bin_rhos = all_good_rhos[in_bin]
            complex_vectors = np.exp(1j * bin_rhos)
            mean_complex = np.mean(complex_vectors)
            averaged_rho_rad[i] = np.angle(mean_complex)
            bin_channel_counts[i] = np.sum(in_bin)
    
    # Sparsity check: expected = (nchan × nscan) / nbin × 0.25
    n_scans = len(scan_numbers)
    expected_per_bin = (total_channels * n_scans) / n_freq_bins
    min_data_points = max(1, int(0.25 * expected_per_bin))
    
    # Report minimum sparsity
    filled_bins_mask = bin_channel_counts > 0
    if np.any(filled_bins_mask):
        min_bin_count = np.min(bin_channel_counts[filled_bins_mask])
        print(f"  Min bin occupancy: {min_bin_count} data points (expected: {expected_per_bin:.1f}, threshold: {min_data_points})")
    
    sparse_bins = (bin_channel_counts > 0) & (bin_channel_counts < min_data_points)
    if np.any(sparse_bins):
        averaged_rho_rad[sparse_bins] = np.nan
        n_sparse = np.sum(sparse_bins)
        print(f"  Flagged {n_sparse} sparse bins (<25% filled)")
    
    good_bins = np.isfinite(averaged_rho_rad)
    n_filled_bins = np.sum(good_bins)
    
    print(f"  Filled {n_filled_bins} / {n_freq_bins} bins")
    
    # Extrapolation
    extrap_bin_mask = np.zeros(n_freq_bins, dtype=bool)
    
    if XF_EX and n_filled_bins >= 2:
        good_bin_indices = np.where(good_bins)[0]
        first_good = good_bin_indices[0]
        last_good = good_bin_indices[-1]
        
        n_fit_bins = max(2, int(np.ceil(XF_EX_FRAC * n_filled_bins)))
        
        if first_good > 0:
            fit_indices = good_bin_indices[:n_fit_bins]
            freq_fit = freq_bin_centers[fit_indices]
            complex_fit = np.exp(1j * averaged_rho_rad[fit_indices])
            
            p_real = np.polyfit(freq_fit, np.real(complex_fit), 1)
            p_imag = np.polyfit(freq_fit, np.imag(complex_fit), 1)
            
            for i in range(first_good):
                freq_extrap = freq_bin_centers[i]
                real_extrap = np.polyval(p_real, freq_extrap)
                imag_extrap = np.polyval(p_imag, freq_extrap)
                complex_extrap = (real_extrap + 1j * imag_extrap) / np.abs(real_extrap + 1j * imag_extrap)
                averaged_rho_rad[i] = np.angle(complex_extrap)
                extrap_bin_mask[i] = True
            
            print(f"  Extrapolated {first_good} bins at lower edge")
        
        if last_good < n_freq_bins - 1:
            fit_indices = good_bin_indices[-n_fit_bins:]
            freq_fit = freq_bin_centers[fit_indices]
            complex_fit = np.exp(1j * averaged_rho_rad[fit_indices])
            
            p_real = np.polyfit(freq_fit, np.real(complex_fit), 1)
            p_imag = np.polyfit(freq_fit, np.imag(complex_fit), 1)
            
            for i in range(last_good + 1, n_freq_bins):
                freq_extrap = freq_bin_centers[i]
                real_extrap = np.polyval(p_real, freq_extrap)
                imag_extrap = np.polyval(p_imag, freq_extrap)
                complex_extrap = (real_extrap + 1j * imag_extrap) / np.abs(real_extrap + 1j * imag_extrap)
                averaged_rho_rad[i] = np.angle(complex_extrap)
                extrap_bin_mask[i] = True
            
            n_extrap_upper = n_freq_bins - 1 - last_good
            print(f"  Extrapolated {n_extrap_upper} bins at upper edge")
    
    valid_avg = np.isfinite(averaged_rho_rad)
    freq_for_interp = freq_bin_centers[valid_avg]
    rho_for_interp = averaged_rho_rad[valid_avg]
    extrap_for_plot = extrap_bin_mask[valid_avg]
    
    # Store as single solution replicated for all scans
    freq_averaged_pooled = {
        'freq': freq_for_interp,
        'rho_rad': rho_for_interp,
        'rho_deg': np.degrees(rho_for_interp),
        'n_bins': len(freq_for_interp),
        'n_good': len(freq_for_interp),
        'extrapolated': extrap_for_plot
    }
    
    for _ in range(n_scans):
        freq_averaged_data_scans.append(freq_averaged_pooled)

else:
    # ============================
    # PER-SCAN MODE: Do frequency averaging with extrapolation per scan
    # ============================
    for scan_idx, scan in enumerate(scan_numbers):
        rho_per_ch = selected_rho_per_channel_scans[scan_idx]
        sol_per_ch = selected_solution_per_channel_scans[scan_idx]
    
        has_solution = ~np.isnan(rho_per_ch) & (sol_per_ch != 'outlier') & (sol_per_ch != '')
        n_with_solution = np.sum(has_solution)
    
        if n_with_solution == 0:
            print(f"  Scan {scan}: No valid solutions, skipping")
            freq_averaged_data_scans.append(None)
            continue
    
        if n_with_solution > XF_MAX_AVG_CHANNELS:
            print(f"  Scan {scan}: Averaging {total_channels} total channels → {XF_MAX_AVG_CHANNELS} bins ({n_with_solution} / {total_channels} good)")
        
            freq_for_avg = freq_ghz[has_solution]
            rho_for_avg = rho_per_ch[has_solution]
        
            # Create evenly-spaced bins across full frequency range
            freq_bin_edges = np.linspace(freq_min, freq_max, XF_MAX_AVG_CHANNELS + 1)
            freq_bin_centers = 0.5 * (freq_bin_edges[:-1] + freq_bin_edges[1:])
        
            freq_averaged_rho_rad = np.full(XF_MAX_AVG_CHANNELS, np.nan)
            bin_channel_counts = np.zeros(XF_MAX_AVG_CHANNELS, dtype=int)
        
            for i in range(XF_MAX_AVG_CHANNELS):
                bin_low = freq_bin_edges[i]
                bin_high = freq_bin_edges[i + 1]
            
                in_bin = (freq_for_avg >= bin_low) & (freq_for_avg < bin_high)
                if i == XF_MAX_AVG_CHANNELS - 1:
                    in_bin = (freq_for_avg >= bin_low) & (freq_for_avg <= bin_high)
            
                if np.any(in_bin):
                    bin_rho_rad = rho_for_avg[in_bin]
                    complex_vectors = np.exp(1j * bin_rho_rad)
                    mean_complex = np.mean(complex_vectors)
                    averaged_rho_rad = np.angle(mean_complex)
                
                    freq_averaged_rho_rad[i] = averaged_rho_rad
                    bin_channel_counts[i] = np.sum(in_bin)
        
            # Sparsity check: expected = (nchan × 1) / nbin × 0.25
            expected_per_bin = total_channels / XF_MAX_AVG_CHANNELS
            min_channels = max(1, int(0.25 * expected_per_bin))
            
            # Report minimum sparsity
            filled_bins_mask = bin_channel_counts > 0
            if np.any(filled_bins_mask):
                min_bin_count = np.min(bin_channel_counts[filled_bins_mask])
                print(f"    Min bin occupancy: {min_bin_count} channels (expected: {expected_per_bin:.1f}, threshold: {min_channels})")
            
            sparse_bins = (bin_channel_counts > 0) & (bin_channel_counts < min_channels)
            if np.any(sparse_bins):
                freq_averaged_rho_rad[sparse_bins] = np.nan
                n_sparse = np.sum(sparse_bins)
                print(f"    Flagged {n_sparse} sparse bins (<25% filled)")
        
            good_bins = np.isfinite(freq_averaged_rho_rad)
            n_filled_bins = np.sum(good_bins)
        
            print(f"    Filled bins: {n_filled_bins} / {XF_MAX_AVG_CHANNELS}")
        
            extrap_bin_mask = np.zeros(XF_MAX_AVG_CHANNELS, dtype=bool)
        
            if XF_EX and n_filled_bins >= 2:
                good_bin_indices = np.where(good_bins)[0]
                first_good = good_bin_indices[0]
                last_good = good_bin_indices[-1]
            
                n_fit_bins = max(2, int(np.ceil(XF_EX_FRAC * n_filled_bins)))
            
                if first_good > 0:
                    fit_indices = good_bin_indices[:n_fit_bins]
                    freq_fit = freq_bin_centers[fit_indices]
                    complex_fit = np.exp(1j * freq_averaged_rho_rad[fit_indices])
                
                    p_real = np.polyfit(freq_fit, np.real(complex_fit), 1)
                    p_imag = np.polyfit(freq_fit, np.imag(complex_fit), 1)
                
                    for i in range(first_good):
                        freq_extrap = freq_bin_centers[i]
                        real_extrap = np.polyval(p_real, freq_extrap)
                        imag_extrap = np.polyval(p_imag, freq_extrap)
                        complex_extrap = (real_extrap + 1j * imag_extrap) / np.abs(real_extrap + 1j * imag_extrap)
                        freq_averaged_rho_rad[i] = np.angle(complex_extrap)
                        extrap_bin_mask[i] = True
                
                    if first_good > 0:
                        print(f"    Extrapolated {first_good} bins at lower edge")
            
                if last_good < XF_MAX_AVG_CHANNELS - 1:
                    fit_indices = good_bin_indices[-n_fit_bins:]
                    freq_fit = freq_bin_centers[fit_indices]
                    complex_fit = np.exp(1j * freq_averaged_rho_rad[fit_indices])
                
                    p_real = np.polyfit(freq_fit, np.real(complex_fit), 1)
                    p_imag = np.polyfit(freq_fit, np.imag(complex_fit), 1)
                
                    for i in range(last_good + 1, XF_MAX_AVG_CHANNELS):
                        freq_extrap = freq_bin_centers[i]
                        real_extrap = np.polyval(p_real, freq_extrap)
                        imag_extrap = np.polyval(p_imag, freq_extrap)
                        complex_extrap = (real_extrap + 1j * imag_extrap) / np.abs(real_extrap + 1j * imag_extrap)
                        freq_averaged_rho_rad[i] = np.angle(complex_extrap)
                        extrap_bin_mask[i] = True
                
                    n_extrap_upper = XF_MAX_AVG_CHANNELS - 1 - last_good
                    if n_extrap_upper > 0:
                        print(f"    Extrapolated {n_extrap_upper} bins at upper edge")
        
            freq_averaged_rho_deg = np.degrees(freq_averaged_rho_rad)
        
            freq_averaged_data_scans.append({
                'freq': freq_bin_centers,
                'rho_deg': freq_averaged_rho_deg,
                'rho_rad': freq_averaged_rho_rad,
                'n_bins': XF_MAX_AVG_CHANNELS,
                'n_good': n_filled_bins,
                'extrapolated': extrap_bin_mask,
                'bin_edges': freq_bin_edges
            })
        else:
            print(f"  Scan {scan}: No averaging needed ({n_with_solution} ≤ {XF_MAX_AVG_CHANNELS})")
            freq_averaged_data_scans.append(None)

print("\n=== Frequency Averaging Complete ===")

# ============================
# Generate CASA XF Calibration Table
# ============================

print("\n=== Generating CASA XF Calibration Table ===")

combine_param = 'scan' if XF_AVG_SCAN else ''
print(f"Creating XF table with combine='{combine_param}'")

polcal(vis=temp_ms,
       field=pacal_name,
       caltable=xftab,
       refant=str(ref_ant),
       solint=f'inf,{XF_CHANINT}ch',
       poltype='Xf',
       combine=combine_param,
       append=False)

print(f"Created XF table structure: {xftab}")

# Read XF table frequency grid
tb.open(xftab + '/SPECTRAL_WINDOW')
chan_freq_xf = tb.getcol('CHAN_FREQ').flatten()
tb.close()

chan_freq_xf_ghz = chan_freq_xf / 1e9

print(f"XF table frequency grid: {len(chan_freq_xf_ghz)} channels")
print(f"  Range: {chan_freq_xf_ghz.min():.3f} - {chan_freq_xf_ghz.max():.3f} GHz")

# Read table structure
tb.open(xftab, nomodify=False)
gains = tb.getcol('CPARAM')
flag_original = tb.getcol('FLAG')
time_col = tb.getcol('TIME')
tb.close()

print(f"  Table dimensions: {gains.shape}")

if XF_AVG_SCAN:
    # ============================
    # Scan-Averaged XF Table
    # ============================
    
    print("\n=== Scan-Averaged XF Table ===")
    
    # Use pre-computed frequency-averaged solution
    scan_sol_data = freq_averaged_data_scans[0]  # All scans point to same solution
    freq_for_interp = scan_sol_data['freq']
    rho_for_interp = scan_sol_data['rho_rad']
    
    print(f"  Using {len(freq_for_interp)} frequency bins for interpolation")
    
    # Smooth if requested
    if XF_USE_SMOOTHING and len(freq_for_interp) >= 5:
        complex_for_smooth = np.exp(1j * rho_for_interp)
        
        window_length = XF_SAVGOL_WINDOW if XF_SAVGOL_WINDOW is not None else min(5, len(freq_for_interp))
        if window_length % 2 == 0:
            window_length += 1
        window_length = min(window_length, len(freq_for_interp))
        if window_length % 2 == 0:
            window_length -= 1
        
        polyorder = min(XF_SAVGOL_POLYORDER, window_length - 1)
        
        real_smooth = savgol_filter(np.real(complex_for_smooth), window_length, polyorder)
        imag_smooth = savgol_filter(np.imag(complex_for_smooth), window_length, polyorder)
        complex_smooth = real_smooth + 1j * imag_smooth
        complex_smooth = complex_smooth / np.abs(complex_smooth)
        print(f"  Applied smoothing (window={window_length}, poly={polyorder})")
    else:
        complex_smooth = np.exp(1j * rho_for_interp)
    
    # Interpolate onto XF table frequencies
    real_interp = np.interp(chan_freq_xf_ghz, freq_for_interp, np.real(complex_smooth))
    imag_interp = np.interp(chan_freq_xf_ghz, freq_for_interp, np.imag(complex_smooth))
    
    complex_interp = real_interp + 1j * imag_interp
    complex_interp = complex_interp / np.abs(complex_interp)
    
    valid_range = (chan_freq_xf_ghz >= freq_for_interp.min()) & (chan_freq_xf_ghz <= freq_for_interp.max())
    
    final_rho_for_corrections = complex_interp
    complex_interp_casa = np.conj(complex_interp)
    
    # Write to table (single time stamp = mean of all scans)
    tb.open(xftab, nomodify=False)
    gains = tb.getcol('CPARAM')
    flags = tb.getcol('FLAG')
    
    n_pol, n_freq, n_antennas = gains.shape
    
    print(f"  XF table structure: {gains.shape} (pol, freq, ant)")
    
    if n_freq != len(chan_freq_xf_ghz):
        print(f"WARNING: Expected {len(chan_freq_xf_ghz)} frequency channels, got {n_freq}")
    
    # Broadcast solution across all antennas: (freq,) -> (1, freq, ant)
    gains[0, :, :] = complex_interp_casa[:, np.newaxis]
    
    # Flag out-of-range channels for all antennas
    flags[0, ~valid_range, :] = True
    
    # Set flagged gains to unity
    gains[0, ~valid_range, :] = 1.0 + 0.0j
    
    # Update time stamp to mean of all scans
    mean_time = np.mean(scan_times)
    time_col[:] = mean_time
    
    tb.putcol('CPARAM', gains)
    tb.putcol('FLAG', flags)
    tb.putcol('TIME', time_col)
    tb.flush()
    tb.close()
    
    print(f"  Wrote scan-averaged XF solutions (time={mean_time:.1f} MJD)")
    print(f"  Flagged {np.sum(~valid_range)} / {len(chan_freq_xf_ghz)} out-of-range channels")
    
    print("\nGenerating final rho values diagnostic plot...")
    
    # freq_averaged_data_scans already contains the pre-computed solution
    plot_final_rho_values(freq_ghz, U_scans, V_scans, rho_rad_scans, selected_rho_per_channel_scans,
                         freq_averaged_data_scans, selected_solution_per_channel_scans, scan_numbers, xfdir, XF_MIN_CROSS_FLUX)
    print(f"Saved: {xfdir}/final_rho_values.png")

else:
    # ============================
    # Per-Scan XF Table
    # ============================
    
    print("\n=== Per-Scan XF Table ===")
    
    tb.open(xftab, nomodify=False)
    gains = tb.getcol('CPARAM')
    flags = tb.getcol('FLAG')
    time_col = tb.getcol('TIME')
    
    n_pol, n_freq, n_ant_total = gains.shape
    
    print(f"  Table dimensions: {gains.shape} (pol, freq, ant×scans)")
    
    n_antennas_per_scan = n_ant_total // n_scans
    
    if n_antennas_per_scan * n_scans != n_ant_total:
        print(f"WARNING: Antenna dimension {n_ant_total} not evenly divisible by {n_scans} scans")
        n_antennas_per_scan = n_ant_total // n_scans
    
    print(f"  Detected {n_antennas_per_scan} antennas per scan, {n_scans} scans")
    print(f"  Scan solutions will be written to antenna indices: ", end='')
    for i in range(n_scans):
        print(f"scan{i}=[{i*n_antennas_per_scan}:{(i+1)*n_antennas_per_scan}]", end=' ')
    print()
    
    per_scan_rho_for_corrections = []
    
    for scan_idx in range(n_scans):
        scan = scan_numbers[scan_idx]
        
        print(f"  Processing scan {scan} (index {scan_idx})...")
        
        scan_sol_data = freq_averaged_data_scans[scan_idx]
        
        if scan_sol_data is not None:
            model_freq_ghz = scan_sol_data['freq']
            model_rho_rad = scan_sol_data['rho_rad']
        else:
            model_freq_ghz = freq_ghz
            model_rho_rad = selected_rho_per_channel_scans[scan_idx]
        
        good_mask = np.isfinite(model_rho_rad)
        
        if np.sum(good_mask) < 2:
            print(f"    Scan {scan}: Insufficient good data, flagging all")
            ant_start = scan_idx * n_antennas_per_scan
            ant_end = (scan_idx + 1) * n_antennas_per_scan
            flags[:, :, ant_start:ant_end] = True
            gains[:, :, ant_start:ant_end] = 1.0 + 0.0j
            per_scan_rho_for_corrections.append(None)
            continue
        
        freq_good = model_freq_ghz[good_mask]
        rho_good = model_rho_rad[good_mask]
        
        # Convert to complex unit vectors
        complex_good = np.exp(1j * rho_good)
        
        # Smooth if requested
        if XF_USE_SMOOTHING and len(freq_good) >= 5:
            window_length = XF_SAVGOL_WINDOW if XF_SAVGOL_WINDOW is not None else min(5, len(freq_good))
            if window_length % 2 == 0:
                window_length += 1
            window_length = min(window_length, len(freq_good))
            if window_length % 2 == 0:
                window_length -= 1
            
            polyorder = min(XF_SAVGOL_POLYORDER, window_length - 1)
            
            real_smooth = savgol_filter(np.real(complex_good), window_length, polyorder)
            imag_smooth = savgol_filter(np.imag(complex_good), window_length, polyorder)
            complex_smooth = real_smooth + 1j * imag_smooth
            complex_smooth = complex_smooth / np.abs(complex_smooth)
        else:
            complex_smooth = complex_good
        
        # Interpolate onto XF table frequencies
        real_interp = np.interp(chan_freq_xf_ghz, freq_good, np.real(complex_smooth))
        imag_interp = np.interp(chan_freq_xf_ghz, freq_good, np.imag(complex_smooth))
        
        complex_interp = real_interp + 1j * imag_interp
        complex_interp = complex_interp / np.abs(complex_interp)
        
        valid_range = (chan_freq_xf_ghz >= freq_good.min()) & (chan_freq_xf_ghz <= freq_good.max())
        
        per_scan_rho_for_corrections.append(complex_interp)
        
        complex_interp_casa = np.conj(complex_interp)
        
        # Write to gains for this scan's antennas
        ant_start = scan_idx * n_antennas_per_scan
        ant_end = (scan_idx + 1) * n_antennas_per_scan
        
        # Broadcast across all antennas for this scan: (1, n_freq, n_ant_per_scan)
        gains_scan = np.tile(complex_interp_casa.reshape(1, -1, 1), (1, 1, n_antennas_per_scan))
        gains[:, :, ant_start:ant_end] = gains_scan
        
        # Flag out-of-range channels for this scan
        flags[0, ~valid_range, ant_start:ant_end] = True
        
        # Set flagged gains to unity
        gains[:, :, ant_start:ant_end][flags[:, :, ant_start:ant_end]] = 1.0 + 0.0j
        
        print(f"    Scan {scan}: Wrote solutions to antenna indices {ant_start}-{ant_end-1}")
        print(f"    Flagged {np.sum(~valid_range)} / {len(chan_freq_xf_ghz)} out-of-range channels")
    
    tb.putcol('CPARAM', gains)
    tb.putcol('FLAG', flags)
    
    for scan_idx in range(n_scans):
        ant_start = scan_idx * n_antennas_per_scan
        ant_end = (scan_idx + 1) * n_antennas_per_scan
        time_col[ant_start:ant_end] = scan_times[scan_idx]
    
    tb.putcol('TIME', time_col)
    tb.flush()
    tb.close()
    
    print(f"  Wrote per-scan XF solutions for {n_scans} scans")
    
    print("\nGenerating final rho values diagnostic plot...")
    
    plot_final_rho_values(freq_ghz, U_scans, V_scans, rho_rad_scans, selected_rho_per_channel_scans,
                         freq_averaged_data_scans, selected_solution_per_channel_scans, scan_numbers, xfdir, XF_MIN_CROSS_FLUX)
    print(f"Saved: {xfdir}/final_rho_values.png")

print(f"SUCCESS: XF table populated: {xftab}")

# ============================
# ============================
# for diagnostic plots.

print("\n=== Applying Cross-Hand Phase Corrections for Diagnostic Plots ===")

if XF_AVG_SCAN:
    print("Scan-averaged solution: applying same rho to each scan independently")
    
    # Interpolate scan-averaged rho solution onto original frequency grid
    rho_interp_rad = np.angle(final_rho_for_corrections)
    rho_interp = np.interp(freq_ghz, chan_freq_xf_ghz, rho_interp_rad)
    
    I_corrected_scans = np.zeros((n_scans, nchan))
    Q_corrected_scans = np.zeros((n_scans, nchan))
    U_corrected_scans = np.zeros((n_scans, nchan))
    V_corrected_scans = np.zeros((n_scans, nchan))
    
    for scan_idx in range(n_scans):
        U_corrected, V_corrected = correct_crosshand_phase(U_scans[scan_idx], V_scans[scan_idx], rho_interp)
        
        I_corrected_scans[scan_idx] = I_scans[scan_idx]
        Q_corrected_scans[scan_idx] = Q_scans[scan_idx]
        U_corrected_scans[scan_idx] = U_corrected
        V_corrected_scans[scan_idx] = V_corrected
    
    print(f"  Applied scan-averaged correction to {n_scans} scan(s)")
    
    # Calculate polarization statistics per scan
    for scan_idx, scan in enumerate(scan_numbers):
        P_corr = np.sqrt(Q_corrected_scans[scan_idx]**2 + U_corrected_scans[scan_idx]**2)
        P_frac = np.nanmedian(P_corr / I_corrected_scans[scan_idx]) * 100
        V_frac = np.nanmedian(np.abs(V_corrected_scans[scan_idx]) / I_corrected_scans[scan_idx]) * 100
        print(f"  Scan {scan}: P/I = {P_frac:.2f}%, |V|/I = {V_frac:.2f}%")
    
    # Generate corrected Stokes spectra plot (zoomed on Stokes V only)
    plot_stokes_spectra(freq_ghz, I_corrected_scans, Q_corrected_scans, 
                                   U_corrected_scans, V_corrected_scans, 
                                   scan_numbers, xfdir, 
                                   filename='stokes_spectra_postXF.png',
                                   title='Post-XF Stokes Parameters (Scan-Averaged Rho)',
                                   zoom_percentile=90,
                                   zoom_stokes=['V'])
    
    # Save corrected data
    np.savez(os.path.join(xfdir, 'final_corrected_stokes_spectra.npz'),
             freq_ghz=freq_ghz,
             scan_numbers=scan_numbers,
             I_corrected_scans=I_corrected_scans,
             Q_corrected_scans=Q_corrected_scans,
             U_corrected_scans=U_corrected_scans,
             V_corrected_scans=V_corrected_scans,
             rho_applied=rho_interp)
    
    print(f"  Saved: {xfdir}/final_corrected_stokes_spectra.npz")

else:
    # Per-scan mode: apply each scan's individual rho solution
    print("Per-scan solutions: applying individual rho per scan")
    
    I_corrected_scans = np.zeros((n_scans, nchan))
    Q_corrected_scans = np.zeros((n_scans, nchan))
    U_corrected_scans = np.zeros((n_scans, nchan))
    V_corrected_scans = np.zeros((n_scans, nchan))
    
    for scan_idx in range(n_scans):
        if per_scan_rho_for_corrections[scan_idx] is None:
            I_corrected_scans[scan_idx] = np.nan
            Q_corrected_scans[scan_idx] = np.nan
            U_corrected_scans[scan_idx] = np.nan
            V_corrected_scans[scan_idx] = np.nan
            continue
        
        # Interpolate this scan's solution onto original frequency grid
        rho_interp_rad = np.angle(per_scan_rho_for_corrections[scan_idx])
        rho_interp = np.interp(freq_ghz, chan_freq_xf_ghz, rho_interp_rad)
        
        U_corr, V_corr = correct_crosshand_phase(U_scans[scan_idx], V_scans[scan_idx], rho_interp)
        
        I_corrected_scans[scan_idx] = I_scans[scan_idx]
        Q_corrected_scans[scan_idx] = Q_scans[scan_idx]
        U_corrected_scans[scan_idx] = U_corr
        V_corrected_scans[scan_idx] = V_corr
    
    print(f"  Applied per-scan corrections to {n_scans} scans")
    
    # Calculate per-scan corrected polarization statistics
    for scan_idx, scan in enumerate(scan_numbers):
        if np.all(np.isnan(U_corrected_scans[scan_idx])):
            continue
        P_corr = np.sqrt(Q_corrected_scans[scan_idx]**2 + U_corrected_scans[scan_idx]**2)
        P_frac = np.nanmedian(P_corr / I_corrected_scans[scan_idx]) * 100
        V_frac = np.nanmedian(np.abs(V_corrected_scans[scan_idx]) / I_corrected_scans[scan_idx]) * 100
        print(f"  Scan {scan}: P = {P_frac:.3f}%, |V| = {V_frac:.3f}%")
    
    # Plot final corrected spectra (zoomed on Stokes V only)
    plot_stokes_spectra(freq_ghz, I_corrected_scans, Q_corrected_scans, 
                                   U_corrected_scans, V_corrected_scans, 
                                   scan_numbers, xfdir, 
                                   filename='stokes_spectra_postXF.png',
                                   title='Post-XF Stokes Parameters (Per-Scan)',
                                   zoom_percentile=90,
                                   zoom_stokes=['V'])
    
    # Save corrected spectra
    np.savez(os.path.join(xfdir, 'final_corrected_stokes_spectra.npz'),
             freq_ghz=freq_ghz,
             scan_numbers=scan_numbers,
             I_corrected_scans=I_corrected_scans,
             Q_corrected_scans=Q_corrected_scans,
             U_corrected_scans=U_corrected_scans,
             V_corrected_scans=V_corrected_scans)
    
    print(f"  Saved: {xfdir}/final_corrected_stokes_spectra.npz")

print(f"Saved: {xfdir}/stokes_spectra_postXF.png")

# Cleanup
shutil.rmtree(temp_ms)
print("\n=== Script Complete ===")
