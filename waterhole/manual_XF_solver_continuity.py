# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
Manual Cross-Hand Phase (XF) Polarization Calibration Solver

This script performs sophisticated polarization calibration analysis on radio 
interferometric data to solve for instrumental cross-hand phase corrections.
It handles the ±π degeneracy, applies parallactic angle corrections, and 
performs full-bandwidth pseudo-continuity analysis with outlier detection.
"""

# Default imports
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import shutil
import time
import datetime
import subprocess
import sys

# Imports for smoothing
from scipy.signal import savgol_filter

# Flush immediately for better logging
import functools
print = functools.partial(print, flush=True)

# Load in config
exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

if PRE_FIELDS != '':
    targets = user_targets
    pcal_names = user_pcals
    target_cal_map = user_cal_map

def stamp():
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


def compute_parallactic_angle(temp_ms, field_name, field_id=None):
    """
    Compute parallactic angle at mid-time using CASA AZ/EL method.
    
    Inputs:
        temp_ms : str - Path to measurement set
        field_name : str - Name of the field/source
        field_id : int or None - Field ID (auto-determined if None)
    
    Outputs:
        chi_deg : float - Parallactic angle at mid-observation time in degrees
        diagnostics : dict - Dictionary with AZ, EL, LAT, mid_time, and observatory name
    
    Computes parallactic angle using: χ = atan2(-sin(A), tan(φ)cos(e) - cos(A)sin(e))
    where A=azimuth, e=elevation, φ=site latitude (CASA convention: A=0 at N, +E).
    """
    print("\n--- Computing parallactic angle (CASA AZ/EL method) ---")
    
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
    else:
        print(f'Field "{field_name}" → FIELD_ID {field_id}')

    # Get phase center and timestamps
    phase_dir = msmd.phasecenter(field_id)
    t_arr = np.array(msmd.timesforfield(field_id), float)
    msmd.close()

    if t_arr.size == 0:
        tb.open(temp_ms)
        time_all = tb.getcol('TIME')
        field_all = tb.getcol('FIELD_ID')
        tb.close()
        t_arr = time_all[field_all == field_id]
    if t_arr.size == 0:
        raise RuntimeError("No timestamps found for selected field.")

    mid_time = float(np.median(t_arr))
    print(f"Mid time: {mid_time:.3f} s (MJD) | N times: {t_arr.size}")

    # Get observatory information
    tb.open(temp_ms + '/OBSERVATION')
    tel_names = tb.getcol('TELESCOPE_NAME')
    tb.close()
    obsname = str(tel_names[0]) if tel_names.size > 0 else 'UNKNOWN'
    print(f"Observatory: {obsname}")

    # Get site geodetic latitude
    pos_meas = me.observatory(obsname)
    pos_wgs = me.measure(pos_meas, 'wgs84')
    lat_rad = q_to_rad_scalar(pos_wgs['m1'])

    # Convert source to AZ/EL at mid-time
    me.doframe(pos_meas)
    me.doframe(me.epoch('utc', qa.quantity(mid_time, 's')))
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

    print(f"AZ={az_deg:.3f}°, EL={el_deg:.3f}°, LAT={lt_deg:.3f}°")
    print(f"Parallactic angle at mid-time: {chi_deg:.3f}°")
    print("--- End parallactic angle computation ---\n")
    
    diagnostics = {
        'az_deg': az_deg, 
        'el_deg': el_deg, 
        'lat_deg': lt_deg,
        'mid_time': mid_time,
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


def plot_stokes_histograms(I_avg, Q_avg, U_avg, V_avg, output_dir):
    """
    Generate histograms of real and imaginary Stokes parameters.
    
    Inputs:
        I_avg, Q_avg, U_avg, V_avg : arrays - Complex averaged Stokes parameters
        output_dir : str - Directory path for saving plot
    
    Outputs:
        None (saves plot to disk)
    
    Creates 2×4 histogram grid showing distributions of real and imaginary
    components for all Stokes parameters.
    """
    fig, ax = plt.subplots(2, 4, figsize=(15, 7), constrained_layout=True)
    
    ax[0,0].hist(np.real(I_avg), bins=40); ax[0,0].set_title('Re(I_avg) per channel')
    ax[0,1].hist(np.real(Q_avg), bins=40); ax[0,1].set_title('Re(Q_avg) per channel')
    ax[0,2].hist(np.real(U_avg), bins=40); ax[0,2].set_title('Re(U_avg) per channel')
    ax[0,3].hist(np.real(V_avg), bins=40); ax[0,3].set_title('Re(V_avg) per channel')
    
    ax[1,0].hist(np.imag(I_avg), bins=40); ax[1,0].set_title('Im(I_avg) per channel')
    ax[1,1].hist(np.imag(Q_avg), bins=40); ax[1,1].set_title('Im(Q_avg) per channel')
    ax[1,2].hist(np.imag(U_avg), bins=40); ax[1,2].set_title('Im(U_avg) per channel')
    ax[1,3].hist(np.imag(V_avg), bins=40); ax[1,3].set_title('Im(V_avg) per channel')
    
    for a in ax.ravel():
        a.set_xlabel('Value')
    
    plt.savefig(os.path.join(output_dir, 'iquv_preXF_avg_hist.png'))
    plt.close(fig)


def plot_stokes_spectra(freq, I_flux, Q_flux, U_flux, V_flux, output_dir, filename, title):
    """
    Generate per-channel spectra for all Stokes parameters.
    
    Inputs:
        freq : array - Frequency array in GHz
        I_flux, Q_flux, U_flux, V_flux : arrays - Stokes flux densities in Jy
        output_dir : str - Directory path for saving plot
        filename : str - Output filename (e.g., 'final_corrected_stokes.png')
        title : str - Main title for the plot
    
    Outputs:
        None (saves plot to disk)
    
    Creates 4-panel plot showing flux density vs frequency for I, Q, U, V.
    """
    fig, ax = plt.subplots(4, 1, figsize=(10, 12), sharex=True, constrained_layout=True)
    
    # Main title for the entire figure
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    ax[0].plot(freq, I_flux, marker='.', linestyle='none', alpha=0.7)
    ax[0].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[0].set_ylabel('Stokes I [Jy]')
    ax[0].grid(True, alpha=0.3)
    
    ax[1].plot(freq, Q_flux, marker='.', linestyle='none', alpha=0.7)
    ax[1].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[1].set_ylabel('Stokes Q [Jy]')
    ax[1].grid(True, alpha=0.3)
    
    ax[2].plot(freq, U_flux, marker='.', linestyle='none', alpha=0.7)
    ax[2].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[2].set_ylabel('Stokes U [Jy]')
    ax[2].grid(True, alpha=0.3)
    
    ax[3].plot(freq, V_flux, marker='.', linestyle='none', alpha=0.7)
    ax[3].axhline(0.0, color='gray', linestyle='--', alpha=0.7)
    ax[3].set_ylabel('Stokes V [Jy]')
    ax[3].set_xlabel('Frequency [GHz]')
    ax[3].grid(True, alpha=0.3)
    
    plt.savefig(os.path.join(output_dir, filename), dpi=150)
    plt.close(fig)

def plot_crosshand_phase_basic(freq, U, V, rho_deg, min_flux, output_dir):
    """
    Generate basic two-panel cross-hand phase diagnostic plot.
    
    Inputs:
        freq : array - Frequency array in GHz
        U : array - Stokes U flux in Jy
        V : array - Stokes V flux in Jy
        rho_deg : array - Cross-hand phase in degrees
        min_flux : float - Minimum flux threshold for reference line
        output_dir : str - Directory path for saving plot
    
    Outputs:
        None (saves plot to disk)
    
    Top panel shows cross-hand flux amplitude, bottom panel shows phase.
    """
    valid = np.isfinite(U) & np.isfinite(V) & np.isfinite(freq) & np.isfinite(rho_deg)
    
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Top: cross-hand flux
    crosshand_flux = np.sqrt(U**2 + V**2)
    ax_top.scatter(freq[valid], crosshand_flux[valid], s=8, c='red', alpha=0.7)
    ax_top.axhline(min_flux, color='gray', linestyle='--', alpha=0.7, 
                   label=f'{min_flux} Jy reference')
    ax_top.set_ylabel('Cross-hand flux [Jy] (√(U² + V²))')
    ax_top.set_title('Per-channel cross-hand flux (before correction)')
    ax_top.grid(True, alpha=0.3)
    ax_top.legend()
    
    # Bottom: cross-hand phase
    ax_bottom.scatter(freq[valid], rho_deg[valid], s=8, c='blue', alpha=0.7)
    ax_bottom.set_xlabel('Frequency [GHz]')
    ax_bottom.set_ylabel('Cross-hand phase ρ [deg]')
    ax_bottom.set_title('Per-channel cross-hand phase (arctan2(-V, U))')
    ax_bottom.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'crosshand_phase.png'))
    plt.close(fig)


# ============================
# MAIN SCRIPT
# ============================

# ============================
# Configuration Parameters
# ============================
temp_ms = 'pol_ang_temp.ms'
xfdir = GAINPLOTS + '/manualXF'
os.makedirs(xfdir, exist_ok=True)

# ====================================
# SPLIT OUT POLARIZATION CALIBRATOR MS
# ====================================

# Calibration tables
ktab = GAINTABLES+'/cal_1GC_'+myms+'.K'
bptab = GAINTABLES+'/cal_1GC_'+myms+'.B'
gptab = GAINTABLES+'/cal_1GC_'+myms+'.Gp'
gatab = GAINTABLES+'/cal_1GC_'+myms+'.Ga'
ftab = GAINTABLES+'/cal_1GC_'+myms+'.F'
dftab  = GAINTABLES+'/cal_1GC_'+myms+'.Df'


kcross  = GAINTABLES+'/cal_1GC_'+myms+'.KCROSS'
xftab  = GAINTABLES+'/cal_1GC_'+myms+'.Xf'

# ------- Apply calibration solutions up to Df (leakage) then split out temporary files

applycal(vis = myms,
 #       applymode='calflagstrict',
        field = pacal_name,
        #calwt = False,
        parang = False,
        gainfield = [pacal_name,pacal_name, bpcal_name, pacal_name, bpcal_name],
        gaintable = [ktab,gptab,bptab,ftab,dftab],
        interp = ['linear','linear','linear','linear','linear'],
        flagbackup = False)

mstransform(vis = myms,
        outputvis = temp_ms,
        field = pacal_name,
        usewtspectrum = False,
        datacolumn='corrected')

# ===============================
# Print INPUTS
# ===============================

print(f"Field: {pacal_name}")
print(f"Target lin. pol. angle: {XF_TARGET_POLANG} deg")
print(f"Channel averaging interval: {XF_CHANINT} channels")
print(f"Min cross-hand flux threshold: {XF_MIN_CROSS_FLUX} Jy")
print(f"Local scatter sigma threshold: {XF_SIGMA_CLIP}")
print(f"Local scatter window size: {XF_CLIP_WINDOW} channels")
print(f"Smooth during XF table creation?: {XF_USE_SMOOTHING}")
if XF_USE_SMOOTHING:
    print(f"Savgol filter window size = {XF_SAVGOL_WINDOW} (None - Calculate within code)")
    print(f"Savgol polynomial order = {XF_SAVGOL_POLYORDER}")

# ============================
# Load Frequency Information
# ============================

print("\nReading channel frequencies...")
tb.open(temp_ms + '/SPECTRAL_WINDOW')
chan_freq = tb.getcol('CHAN_FREQ')
tb.close()
freq_ghz = (chan_freq / 1.0e9).flatten()
print(f"Frequency array shape: {freq_ghz.shape}")

if XF_MAX_AVG_CHANNELS is None:
    XF_MAX_AVG_CHANNELS = int(len(chan_freq) / XF_CHANINT)
    print(f"No channel averaging specified, adopting N_CHAN / XF_CHAN_INT = {XF_MAX_AVG_CHANNELS}")

elif not isinstance(XF_MAX_AVG_CHANNELS, int):
    XF_MAX_AVG_CHANNELS = int(len(chan_freq) / XF_CHANINT)
    print(f"XF_MAX_AVG_CHANNELS is neither None nor an integer, defaulting to N_CHAN / XF_CHAN_INT = {XF_MAX_AVG_CHANNELS}")

else:
    print(f"XF_MAX_AVG_CHANNELS is specifiec as: {XF_MAX_AVG_CHANNELS}")  

# ============================
# Load Visibility Data
# ============================
print("\nOpening MS and selecting field + pseudo-Stokes I,Q,U,V...")
ms.open(temp_ms)
ms.selectinit(reset=True)
ok = ms.msselect({'field': pacal_name})
print(f'msselect(field="{pacal_name}") ok? {ok}')
ms.selectpolarization(['I','Q','U','V'])

d = ms.getdata(['data', 'flag', 'weight_spectrum', 'weight'])

missing = [k for k in ('data', 'flag') if k not in d]
if missing:
    ms.close()
    sys.exit(f"ERROR: ms.getdata() did not return required key(s): {missing}")

have_wspec = 'weight_spectrum' in d
if have_wspec:
    print("Weights: using WEIGHT_SPECTRUM from ms.getdata().")
elif 'weight' in d:
    print("Weights: WEIGHT_SPECTRUM not present; using WEIGHT from ms.getdata().")
else:
    ms.close()
    sys.exit("ERROR: Neither WEIGHT_SPECTRUM nor WEIGHT present in ms.getdata() result.")

ms.close()
print("MS closed.")

data = d['data']
raw_flag = np.asarray(d['flag'])
ncorr, nchan, nrow = data.shape
print(f"DATA shape ([I,Q,U,V],chan,row): {data.shape}")

if ncorr != 4:
    sys.exit(f"ERROR: Expected 4 correlations (I,Q,U,V), got {ncorr}.")

# ============================
# Process Flags and Weights
# ============================
print("\n--- Flag handling ---")
flag_iquv = build_iquv_flag(raw_flag, nchan, nrow)
print(f"Final FLAG_IQUV shape: {flag_iquv.shape}")
print("--- End flag handling ---\n")

print("--- Weight handling ---")
W = build_iquv_weights(d, have_wspec, nchan, nrow)
W_eff = np.where(flag_iquv, 0.0, W)

nonzero = int(np.sum(W_eff > 0))
total = int(W_eff.size)
print(f"Unflagged (non-zero) weights: {nonzero:,d} / {total:,d} ({100.0*nonzero/total:.1f}%)")
print(f"Weight stats: min={np.nanmin(W_eff):.6g}, median={np.nanmedian(W_eff):.6g}, max={np.nanmax(W_eff):.6g}")
print("--- End weight handling ---\n")

# ============================
# Compute Weighted Averages
# ============================
print("Computing weighted vector averages per channel...")
vis_avg, sigma_proxy = compute_weighted_averages(data, W, flag_iquv)

I_avg = vis_avg[0]
Q_avg = vis_avg[1]
U_avg = vis_avg[2]
V_avg = vis_avg[3]

I_flux = np.real(I_avg)
Q_flux = np.real(Q_avg)
U_flux = np.real(U_avg)
V_flux = np.real(V_avg)

I_sigma = sigma_proxy[0]
Q_sigma = sigma_proxy[1]
U_sigma = sigma_proxy[2]
V_sigma = sigma_proxy[3]

print("Averaging complete.")
print(f"Median I_flux={np.nanmedian(I_flux):.6g} Jy, Median Q_flux={np.nanmedian(Q_flux):.6g} Jy")
print(f"Median U_flux={np.nanmedian(U_flux):.6g} Jy, Median V_flux={np.nanmedian(V_flux):.6g} Jy")

# ============================
# Consistency Checks
# ============================
print("\n=== Consistency checks (per-channel) ===")
quick_stats('Re(I_avg)', np.real(I_avg))
quick_stats('Im(I_avg)', np.imag(I_avg))
quick_stats('Re(Q_avg)', np.real(Q_avg))
quick_stats('Im(Q_avg)', np.imag(Q_avg))
quick_stats('Re(U_avg)', np.real(U_avg))
quick_stats('Im(U_avg)', np.imag(U_avg))
quick_stats('Re(V_avg)', np.real(V_avg))
quick_stats('Im(V_avg)', np.imag(V_avg))

im_fracs = {
    'I': calc_im_fraction(I_avg),
    'Q': calc_im_fraction(Q_avg),
    'U': calc_im_fraction(U_avg),
    'V': calc_im_fraction(V_avg)
}

print(f'\nImaginary flux fractions (should be <<10% for good calibration):')
for stokes, frac in im_fracs.items():
    status = "⚠️ HIGH" if frac > 10.0 else "✓ OK"
    print(f'  {stokes}: {frac:.1f}% {status}')

high_im_stokes = [f'{s}({f:.1f}%)' for s, f in im_fracs.items() if f > 10.0]
if high_im_stokes:
    print(f'WARNING: High imaginary amplitude in {", ".join(high_im_stokes)} - check prior calibration!')

# ============================
# Generate Diagnostic Plots
# ============================
print("\nGenerating diagnostic plots...")

plot_stokes_histograms(I_avg, Q_avg, U_avg, V_avg, xfdir)

f, I_plot, Q_plot, U_plot, V_plot = align_stokes_arrays(freq_ghz, I_flux, Q_flux, U_flux, V_flux)
plot_stokes_spectra(f, I_plot, Q_plot, U_plot, V_plot, xfdir, 'iquv_pre_XF.png', title='Pre-XF Stokes Parameters')

np.savez(os.path.join(xfdir, 'iquv_vector_avgs_per_channel.npz'),
         freq_ghz=f,
         I_avg=I_avg[:len(f)], Q_avg=Q_avg[:len(f)], 
         U_avg=U_avg[:len(f)], V_avg=V_avg[:len(f)],
         I_flux=I_plot, Q_flux=Q_plot, U_flux=U_plot, V_flux=V_plot,
         I_sigma=I_sigma[:len(f)], Q_sigma=Q_sigma[:len(f)], 
         U_sigma=U_sigma[:len(f)], V_sigma=V_sigma[:len(f)])

# ============================
# Cross-Hand Phase Solve
# ============================
print("\n=== Solving cross-hand phase: ρ = arctan2(-V, U) ===")

# Use aligned arrays
f_use = f
I_use = I_plot
Q_use = Q_plot
U_use = U_plot
V_use = V_plot

valid = np.isfinite(U_use) & np.isfinite(V_use) & np.isfinite(f_use)
num_valid = int(np.sum(valid))
print(f"Valid channels for cross-hand phase: {num_valid} / {f_use.size}")

rho_rad = np.full_like(U_use, np.nan, dtype=float)
rho_rad[valid] = np.arctan2(-V_use[valid], U_use[valid])
rho_deg = np.degrees(rho_rad)
rho_deg_wrapped = ((rho_deg + 180.0) % 360.0) - 180.0

if num_valid > 0:
    med_rho = np.nanmedian(rho_deg_wrapped)
    std_rho = np.nanstd(rho_deg_wrapped)
    print(f"Cross-hand phase: median={med_rho:.3f} deg, std={std_rho:.3f} deg")

plot_crosshand_phase_basic(f_use, U_use, V_use, rho_deg_wrapped, XF_MIN_CROSS_FLUX, xfdir)

np.savez(os.path.join(xfdir, 'crosshand_phase_per_channel.npz'),
         freq_ghz=f_use,
         rho_deg=rho_deg_wrapped,
         rho_rad=rho_rad)

print(f"Saved: {xfdir}/crosshand_phase.png, {xfdir}/crosshand_phase_per_channel.npz")

# ============================
# Parallactic Angle
# ============================
chi_deg, parang_diagnostics = compute_parallactic_angle(temp_ms, pacal_name)

# ============================
# Top 5% Bandwidth Analysis
# ============================
print("\n=== Top 5% bandwidth analysis ===")

unflagged_mask = np.isfinite(Q_use) & np.isfinite(U_use) & np.isfinite(V_use) & np.isfinite(rho_rad)
unflagged_idx = np.where(unflagged_mask)[0]

if len(unflagged_idx) == 0:
    sys.exit("ERROR: No unflagged channels found!")

print(f"Found {len(unflagged_idx)} unflagged channels out of {len(f_use)} total channels")

# Select top 5% by frequency
n_top5 = max(1, int(0.05 * len(unflagged_idx)))
unflagged_freqs = f_use[unflagged_idx]
freq_sorted_idx = np.argsort(unflagged_freqs)
top5_relative_idx = freq_sorted_idx[-n_top5:]
top5_idx = unflagged_idx[top5_relative_idx]

print(f"Using top {n_top5} unflagged channels ({100*n_top5/len(unflagged_idx):.1f}% of unflagged bandwidth)")
print(f"Frequency range: {f_use[top5_idx].min():.3f} - {f_use[top5_idx].max():.3f} GHz")

# Calculate averages for top 5%
I_top5 = np.mean(I_use[top5_idx])
Q_top5 = np.mean(Q_use[top5_idx])
U_top5 = np.mean(U_use[top5_idx])
V_top5 = np.mean(V_use[top5_idx])

print(f"\nRaw averages (top 5% unflagged bandwidth):")
print(f"  I = {I_top5:.6f} Jy")
print(f"  Q = {Q_top5:.6f} Jy, U = {U_top5:.6f} Jy, V = {V_top5:.6f} Jy")
print(f"  Q/I = {Q_top5/I_top5*100:.3f}%, U/I = {U_top5/I_top5*100:.3f}%, V/I = {V_top5/I_top5*100:.3f}%")

# Calculate cross-hand phase for top 5%
rho_top5_rad = np.mean(rho_rad[top5_idx])
rho_top5_deg = np.degrees(rho_top5_rad)
print(f"\nCross-hand phase (top 5% average): {rho_top5_deg:.3f} deg ({rho_top5_rad:.6f} rad)")

# Define ±π solutions
rho_pos = rho_top5_rad
rho_neg = rho_top5_rad + np.pi

# Wrap to [-π, π]
rho_pos = ((rho_pos + np.pi) % (2*np.pi)) - np.pi
rho_neg = ((rho_neg + np.pi) % (2*np.pi)) - np.pi

# Ensure opposite signs
if rho_pos * rho_neg > 0:
    if rho_pos > 0:
        rho_neg = rho_pos - np.pi
    else:
        rho_neg = rho_pos + np.pi
    rho_neg = ((rho_neg + np.pi) % (2*np.pi)) - np.pi

if rho_pos < 0 and rho_neg > 0:
    rho_pos, rho_neg = rho_neg, rho_pos

print(f"Positive solution: {np.degrees(rho_pos):.3f} deg (rho = {rho_pos:.6f} rad)")
print(f"Negative solution: {np.degrees(rho_neg):.3f} deg (rho = {rho_neg:.6f} rad)")

assert -np.pi <= rho_pos <= np.pi, f"rho_pos ({rho_pos:.6f}) not in [-π, π]"
assert -np.pi <= rho_neg <= np.pi, f"rho_neg ({rho_neg:.6f}) not in [-π, π]"
assert rho_pos >= 0, f"rho_pos ({rho_pos:.6f}) should be non-negative"
assert rho_neg <= 0, f"rho_neg ({rho_neg:.6f}) should be non-positive"
print("✓ Both solutions are within [-π, π] with opposite signs")

# Apply corrections for both solutions
U_xh_pos, V_xh_pos = correct_crosshand_phase(U_top5, V_top5, rho_pos)
U_xh_neg, V_xh_neg = correct_crosshand_phase(U_top5, V_top5, rho_neg)

print(f"\nAfter cross-hand phase correction:")
print(f"  Positive: U = {U_xh_pos:.6f} Jy, V = {V_xh_pos:.6f} Jy (V/I = {V_xh_pos/I_top5*100:.3f}%)")
print(f"  Negative: U = {U_xh_neg:.6f} Jy, V = {V_xh_neg:.6f} Jy (V/I = {V_xh_neg/I_top5*100:.3f}%)")

# Apply parallactic angle correction
Q_final_pos, U_final_pos = correct_parallactic_angle(Q_top5, U_xh_pos, chi_deg)
Q_final_neg, U_final_neg = correct_parallactic_angle(Q_top5, U_xh_neg, chi_deg)

print(f"\nAfter parallactic angle correction (χ = {chi_deg:.3f}°):")
print(f"  Positive: Q = {Q_final_pos:.6f} Jy, U = {U_final_pos:.6f} Jy, V = {V_xh_pos:.6f} Jy")
print(f"  Negative: Q = {Q_final_neg:.6f} Jy, U = {U_final_neg:.6f} Jy, V = {V_xh_neg:.6f} Jy")

P_pos = np.sqrt(Q_final_pos**2 + U_final_pos**2) / I_top5 * 100
P_neg = np.sqrt(Q_final_neg**2 + U_final_neg**2) / I_top5 * 100

print(f"\nTotal linear polarization fraction:")
print(f"  Positive solution: P = {P_pos:.3f}%")
print(f"  Negative solution: P = {P_neg:.3f}%")

# ============================
# RM Trial Analysis
# ============================
print("\n=== RM Trial Analysis ===")
print(f"Target angle: {XF_TARGET_POLANG:.1f} degrees")
print("Trialing 21 RM values from -7 to +7 rad/m² to select best cross-hand phase solution")

rm_trials = np.linspace(-7.0, 7.0, 21)
c = 2.998e8  # Speed of light in m/s

freq_top5_hz = np.mean(f_use[top5_idx]) * 1e9
lambda_sq = (c / freq_top5_hz)**2
print(f"Average frequency (top 5%): {freq_top5_hz/1e9:.3f} GHz")
print(f"Average wavelength squared: {lambda_sq:.6e} m²")

results_pos = []
results_neg = []

print(f"\n{'RM':>8s} {'Pos_Angle':>10s} {'Neg_Angle':>10s} {'Pos_Diff':>10s} {'Neg_Diff':>10s}")
print("-" * 50)

for rm in rm_trials:
    angle_pos = calculate_derotated_angle(Q_final_pos, U_final_pos, rm, lambda_sq)
    angle_neg = calculate_derotated_angle(Q_final_neg, U_final_neg, rm, lambda_sq)
    
    diff_pos = abs(angle_pos - XF_TARGET_POLANG)
    diff_neg = abs(angle_neg - XF_TARGET_POLANG)
    
    if diff_pos > 90:
        diff_pos = 180 - diff_pos
    if diff_neg > 90:
        diff_neg = 180 - diff_neg
    
    results_pos.append({'rm': rm, 'angle': angle_pos, 'diff': diff_pos})
    results_neg.append({'rm': rm, 'angle': angle_neg, 'diff': diff_neg})
    
    print(f"{rm:8.2f} {angle_pos:10.3f} {angle_neg:10.3f} {diff_pos:10.3f} {diff_neg:10.3f}")

total_diff_pos = sum(r['diff'] for r in results_pos)
total_diff_neg = sum(r['diff'] for r in results_neg)

print(f"\nTotal differences across all RM trials:")
print(f"  Positive solution: {total_diff_pos:.3f}° (sum of {len(results_pos)} trials)")
print(f"  Negative solution: {total_diff_neg:.3f}° (sum of {len(results_neg)} trials)")

# Select best solution
if total_diff_pos < total_diff_neg:
    best_overall = 'positive'
    best_rho = rho_pos
    best_Q = Q_final_pos
    best_U = U_final_pos
    best_V = V_xh_pos
    selected_total_diff = total_diff_pos
else:
    best_overall = 'negative'
    best_rho = rho_neg
    best_Q = Q_final_neg
    best_U = U_final_neg
    best_V = V_xh_neg
    selected_total_diff = total_diff_neg

print(f"\n=== SOLUTION SELECTION RESULTS ===")
print(f"Winner: {best_overall} solution")
print(f"Selected cross-hand phase: {np.degrees(best_rho):.3f}° ({best_rho:.6f} rad)")
print(f"Final: Q/I = {best_Q/I_top5*100:.3f}%, U/I = {best_U/I_top5*100:.3f}%, "
      f"V/I = {best_V/I_top5*100:.3f}%, P = {np.sqrt(best_Q**2+best_U**2)/I_top5*100:.3f}%")

# ============================
# Full-Bandwidth Pseudo-Continuity Analysis
# ============================
print("\n=== Pseudo-Continuity Analysis Across Full Bandwidth ===")
print("Determining cross-hand phase solutions in ~5% frequency chunks with continuity constraint")

n_channels = len(f_use)
selected_rho_per_channel = np.full(n_channels, np.nan)
selected_solution_per_channel = np.full(n_channels, '', dtype='U8')

valid_channels = unflagged_idx

# Define frequency chunks
target_channels_per_chunk = max(1, int(0.05 * len(valid_channels)))
n_chunks = max(5, min(40, len(valid_channels) // target_channels_per_chunk))

unflagged_freqs = f_use[valid_channels]
freq_min, freq_max = unflagged_freqs.min(), unflagged_freqs.max()
freq_edges = np.linspace(freq_min, freq_max, n_chunks + 1)

print(f"Dividing {len(valid_channels)} unflagged channels into {n_chunks} frequency chunks")
print(f"Frequency range: {freq_min:.3f} - {freq_max:.3f} GHz")

previous_rho = best_rho
previous_solution = best_overall
print(f"Starting with top-of-band solution: {previous_solution} (rho = {np.degrees(previous_rho):.3f}°)")

chunk_results = []

# Process from high to low frequency
for i in range(n_chunks-1, -1, -1):
    freq_low = freq_edges[i]
    freq_high = freq_edges[i+1]
    
    chunk_mask = (f_use >= freq_low) & (f_use < freq_high) & unflagged_mask
    chunk_channels = np.where(chunk_mask)[0]
    
    if len(chunk_channels) == 0:
        print(f"Chunk {i+1:2d}: {freq_low:.3f}-{freq_high:.3f} GHz - No valid channels, skipping")
        continue
    
    # Calculate chunk averages
    Q_chunk = np.mean(Q_use[chunk_channels])
    U_chunk = np.mean(U_use[chunk_channels])
    V_chunk = np.mean(V_use[chunk_channels])
    
    rho_chunk_raw = np.arctan2(-V_chunk, U_chunk)
    
    # Define ±π solutions
    rho_chunk_option1 = ((rho_chunk_raw + np.pi) % (2*np.pi)) - np.pi
    rho_chunk_option2 = ((rho_chunk_raw + np.pi + np.pi) % (2*np.pi)) - np.pi
    
    # Calculate proximity to previous solution
    prev_rho_deg = np.degrees(previous_rho)
    proximity_option1 = angle_difference(np.degrees(rho_chunk_option1), prev_rho_deg)
    proximity_option2 = angle_difference(np.degrees(rho_chunk_option2), prev_rho_deg)
    
    # Select based on proximity
    if proximity_option1 <= proximity_option2:
        use_pi_correction = False
        selected_rho_chunk = rho_chunk_option1
        selected_proximity = proximity_option1
    else:
        use_pi_correction = True
        selected_rho_chunk = rho_chunk_option2
        selected_proximity = proximity_option2
    
    # Apply to individual channels
    for ch in chunk_channels:
        rho_ch_raw = rho_rad[ch]
        if use_pi_correction:
            rho_ch_corrected = rho_ch_raw + np.pi
        else:
            rho_ch_corrected = rho_ch_raw
        rho_ch_corrected = ((rho_ch_corrected + np.pi) % (2*np.pi)) - np.pi
        selected_rho_per_channel[ch] = rho_ch_corrected
    
    solution_name = 'no_pi' if not use_pi_correction else 'plus_pi'
    
    chunk_info = {
        'chunk_id': i+1,
        'freq_range': (freq_low, freq_high),
        'n_channels': len(chunk_channels),
        'use_pi_correction': use_pi_correction,
        'chunk_avg_rho_deg': np.degrees(selected_rho_chunk),
        'proximity_deg': selected_proximity
    }
    chunk_results.append(chunk_info)
    
    for ch in chunk_channels:
        selected_solution_per_channel[ch] = solution_name
    
    pi_status = "+π" if use_pi_correction else "raw"
    print(f"Chunk {i+1:2d}: {freq_low:.3f}-{freq_high:.3f} GHz ({len(chunk_channels):3d} ch) - "
          f"{pi_status:4s} (avg_rho={np.degrees(selected_rho_chunk):7.3f}°, "
          f"prev={prev_rho_deg:7.3f}°, diff={selected_proximity:5.1f}°)")
    
    previous_rho = selected_rho_chunk
    previous_solution = solution_name

chunk_results.reverse()
print(f"\nPseudo-continuity analysis complete: Processed {len(chunk_results)} chunks")

# ============================
# Two-Stage Outlier Detection
# ============================
print("\n=== Two-Stage Outlier Detection ===")

has_solution = ~np.isnan(selected_rho_per_channel)
n_with_solution = np.sum(has_solution)

# STAGE 1: Low flux flagging
print(f"\nStage 1: Low flux flagging (< {XF_MIN_CROSS_FLUX} Jy)")
low_flux_mask = np.sqrt(U_flux**2 + V_flux**2) < XF_MIN_CROSS_FLUX
low_flux_channels = np.where(low_flux_mask)[0]
n_low_flux = len(low_flux_channels)
print(f"  Found {n_low_flux} channels with cross-hand flux < {XF_MIN_CROSS_FLUX} Jy")

for ch in low_flux_channels:
    selected_solution_per_channel[ch] = 'outlier'

# STAGE 2: Local scatter flagging
print(f"\nStage 2: Local scatter flagging")

if n_with_solution > 10:
    freq_for_outliers = f_use[has_solution]
    rho_deg_for_outliers = np.degrees(selected_rho_per_channel[has_solution])
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
        
        print(f"  Local stdev threshold: {threshold:.4f}")
        
        n_flagged = 0
        for i in range(len(sorted_rho_deg)):
            if not np.isnan(local_stdevs[i]) and local_stdevs[i] > threshold:
                original_idx = freq_sort_idx[i]
                local_scatter_outliers[original_idx] = True
                n_flagged += 1
        
        print(f"  Flagged {n_flagged} channels ({n_flagged/len(valid_stdevs)*100:.1f}%)")
    else:
        n_flagged = 0
        print(f"  Insufficient data for threshold calculation")
    
    n_stage1_outliers = np.sum(local_scatter_outliers)
    
    outlier_channels = channels_for_outliers[local_scatter_outliers]
    for ch in outlier_channels:
        if selected_solution_per_channel[ch] != 'outlier':
            selected_solution_per_channel[ch] = 'outlier'
    
    n_total_outliers = n_low_flux + n_stage1_outliers
    print(f"\nTotal flagged: {n_total_outliers} / {n_with_solution} channels")
else:
    print(f"Insufficient data for outlier detection ({n_with_solution} channels)")
    n_total_outliers = n_low_flux

print("\n=== Two-Stage Outlier Detection Complete ===")

# ============================
# Frequency Averaging
# ============================
print("\n=== Frequency Averaging ===")

freq_averaged_data = None

if n_with_solution > XF_MAX_AVG_CHANNELS:
    print(f"Averaging {n_with_solution} channels → {XF_MAX_AVG_CHANNELS} bins")
    
    freq_for_avg = f_use[has_solution]
    rho_for_avg = selected_rho_per_channel[has_solution]
    solution_for_avg = selected_solution_per_channel[has_solution]
    
    freq_min, freq_max = freq_for_avg.min(), freq_for_avg.max()
    freq_bin_edges = np.linspace(freq_min, freq_max, XF_MAX_AVG_CHANNELS + 1)
    freq_bin_centers = 0.5 * (freq_bin_edges[:-1] + freq_bin_edges[1:])
    
    freq_averaged_rho_rad = np.full(XF_MAX_AVG_CHANNELS, np.nan)
    freq_averaged_solution = np.full(XF_MAX_AVG_CHANNELS, '', dtype='U8')
    
    n_good_bins = 0
    n_flagged_bins = 0
    
    for i in range(XF_MAX_AVG_CHANNELS):
        bin_low = freq_bin_edges[i]
        bin_high = freq_bin_edges[i + 1]
        
        in_bin = (freq_for_avg >= bin_low) & (freq_for_avg < bin_high)
        if i == XF_MAX_AVG_CHANNELS - 1:
            in_bin = (freq_for_avg >= bin_low) & (freq_for_avg <= bin_high)
        
        bin_channels = np.where(in_bin)[0]
        
        if len(bin_channels) == 0:
            continue
        
        bin_rho_rad = rho_for_avg[bin_channels]
        bin_solutions = solution_for_avg[bin_channels]
        
        good_channels = bin_solutions != 'outlier'
        
        if np.any(good_channels):
            good_rho_rad = bin_rho_rad[good_channels]
            complex_vectors = np.exp(1j * good_rho_rad)
            mean_complex = np.mean(complex_vectors)
            averaged_rho_rad = np.angle(mean_complex)
            
            freq_averaged_rho_rad[i] = averaged_rho_rad
            
            unique_solutions, counts = np.unique(bin_solutions[good_channels], return_counts=True)
            predominant_solution = unique_solutions[np.argmax(counts)]
            freq_averaged_solution[i] = predominant_solution
            
            n_good_bins += 1
        else:
            freq_averaged_rho_rad[i] = np.nan
            freq_averaged_solution[i] = 'outlier'
            n_flagged_bins += 1
    
    freq_averaged_rho_deg = np.degrees(freq_averaged_rho_rad)
    
    print(f"  Good bins: {n_good_bins} / {XF_MAX_AVG_CHANNELS}")
    print(f"  Flagged bins: {n_flagged_bins} / {XF_MAX_AVG_CHANNELS}")
    
    freq_averaged_data = {
        'freq': freq_bin_centers,
        'rho_deg': freq_averaged_rho_deg,
        'solution': freq_averaged_solution,
        'n_bins': XF_MAX_AVG_CHANNELS,
        'n_good': n_good_bins,
        'n_flagged': n_flagged_bins
    }
else:
    print(f"Skipping frequency averaging ({n_with_solution} ≤ {XF_MAX_AVG_CHANNELS} channels)")

print("\n=== Frequency Averaging Complete ===")

# ============================
# Three-Panel Diagnostic Plot
# ============================
print("\n=== Generating three-panel diagnostic plot ===")

if n_with_solution > 0:
    freq_plot = f_use[has_solution]
    rho_deg_plot = np.degrees(selected_rho_per_channel[has_solution])
    solution_plot = selected_solution_per_channel[has_solution]
    
    fig, (ax_top, ax_middle, ax_bottom) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)
    
    # Top: Cross-hand flux
    crosshand_flux = np.sqrt(U_flux[has_solution]**2 + V_flux[has_solution]**2)
    raw_mask = solution_plot == 'no_pi'
    plus_pi_mask = solution_plot == 'plus_pi'
    outlier_mask = solution_plot == 'outlier'
    good_mask = ~outlier_mask
    
    if np.any(good_mask):
        ax_top.scatter(freq_plot[good_mask], crosshand_flux[good_mask], 
                      c='blue', alpha=0.6, label='Cross-hand flux')
    if np.any(outlier_mask):
        ax_top.scatter(freq_plot[outlier_mask], crosshand_flux[outlier_mask], 
                      c='grey', alpha=0.3, label=f'Outliers ({np.sum(outlier_mask)})')
    
    ax_top.axhline(y=XF_MIN_CROSS_FLUX, color='red', linestyle='--', alpha=0.7, 
                   label=f'{XF_MIN_CROSS_FLUX} Jy reference')
    ax_top.set_ylabel('Cross-Hand Flux (Jy)')
    ax_top.set_title('Cross-Hand Flux vs Frequency')
    ax_top.grid(True, alpha=0.3)
    ax_top.legend()
    
    # Function to plot phase data
    def plot_phase_panel(ax, title_suffix=""):
        if np.any(raw_mask):
            ax.scatter(freq_plot[raw_mask], rho_deg_plot[raw_mask], 
                      c='blue', s=15, alpha=0.8, label='Raw (no π correction)', zorder=3)
        if np.any(plus_pi_mask):
            ax.scatter(freq_plot[plus_pi_mask], rho_deg_plot[plus_pi_mask], 
                      c='red', s=15, alpha=0.8, label='+π correction', zorder=3)
        if np.any(outlier_mask):
            ax.scatter(freq_plot[outlier_mask], rho_deg_plot[outlier_mask], 
                      c='black', s=25, alpha=0.9, marker='x', linewidths=2,
                      label=f'Outliers ({np.sum(outlier_mask)})', zorder=4)
        
        for cr in chunk_results[1:]:
            ax.axvline(cr['freq_range'][0], color='gray', alpha=0.3, 
                      linestyle='--', linewidth=0.5, zorder=1)
        
        if freq_averaged_data is not None:
            avg_freq = freq_averaged_data['freq']
            avg_rho = freq_averaged_data['rho_deg']
            valid_avg_mask = ~np.isnan(avg_rho)
            
            if np.any(valid_avg_mask):
                ax.scatter(avg_freq[valid_avg_mask], avg_rho[valid_avg_mask],
                          c='purple', s=60, alpha=0.9, marker='s', edgecolors='black', linewidths=1,
                          label=f'Averaged ({freq_averaged_data["n_good"]} bins)', zorder=5)
        
        ax.set_ylabel('Cross-hand phase ρ [deg]')
        ax.grid(True, alpha=0.3, zorder=0)
        ax.set_title(f'Cross-hand phase {title_suffix}')
        ax.legend(fontsize=9)
    
    # Middle: Full range
    plot_phase_panel(ax_middle, "(Full Range)")
    y_min_full, y_max_full = rho_deg_plot.min(), rho_deg_plot.max()
    y_margin_full = 0.05 * (y_max_full - y_min_full)
    ax_middle.set_ylim(y_min_full - y_margin_full, y_max_full + y_margin_full)
    
    # Bottom: Zoomed to good data
    plot_phase_panel(ax_bottom, "(Zoomed to Good Data)")
    
    if np.any(good_mask):
        good_rho = rho_deg_plot[good_mask]
        y_min_good, y_max_good = good_rho.min(), good_rho.max()
        
        if freq_averaged_data is not None:
            avg_good_data = freq_averaged_data['rho_deg'][~np.isnan(freq_averaged_data['rho_deg'])]
            if len(avg_good_data) > 0:
                y_min_good = min(y_min_good, avg_good_data.min())
                y_max_good = max(y_max_good, avg_good_data.max())
        
        y_margin_good = max(0.5 * (y_max_good - y_min_good), 1.0)
        ax_bottom.set_ylim(y_min_good - y_margin_good, y_max_good + y_margin_good)
    else:
        ax_bottom.set_ylim(y_min_full - y_margin_full, y_max_full + y_margin_full)
    
    ax_bottom.set_xlabel('Frequency [GHz]')
    
    plt.tight_layout()
    plt.savefig(os.path.join(xfdir, 'crosshand_phase_fullband_flagged.png'), dpi=150)
    plt.close(fig)
    
    print(f"Saved: {xfdir}/crosshand_phase_fullband_flagged.png")
    
    # Save results
    save_dict = {
        'freq_ghz': freq_plot,
        'rho_deg': rho_deg_plot,
        'rho_rad': selected_rho_per_channel[has_solution],
        'solution_type': solution_plot,
        'outlier_flags': outlier_mask,
        'chunk_info': [cr for cr in chunk_results]
    }
    
    if freq_averaged_data is not None:
        save_dict['freq_averaged_ghz'] = freq_averaged_data['freq']
        save_dict['freq_averaged_rho_deg'] = freq_averaged_data['rho_deg']
        save_dict['freq_averaged_solution'] = freq_averaged_data['solution']
    
    np.savez(os.path.join(xfdir, 'crosshand_phase_fullband_flagged.npz'), **save_dict)
    print(f"Saved: {xfdir}/crosshand_phase_fullband_flagged.npz")
else:
    print("ERROR: No channels have cross-hand phase solutions!")

# ============================
# Generate CASA XF Table with Computed Solutions
# ============================
print("\n=== Generating CASA XF calibration table ===")

# Create initial XF table structure with polcal
polcal(vis=temp_ms,
       field=pacal_name,
       caltable=xftab,
       refant=str(ref_ant),
       solint=f'inf,{XF_CHANINT}ch',
       poltype='Xf',
       combine='',
       append=False)

print(f"Created XF table structure: {xftab}")

# Read XF table frequencies
tb.open(xftab + '/SPECTRAL_WINDOW')
chan_freq_xf = tb.getcol('CHAN_FREQ').flatten()  # Hz
tb.close()

chan_freq_xf_ghz = chan_freq_xf / 1e9  # Convert to GHz

print(f"XF table: {len(chan_freq_xf_ghz)} frequency channels")
print(f"  Frequency range: {chan_freq_xf_ghz.min():.3f} - {chan_freq_xf_ghz.max():.3f} GHz")

# Decide which data to use for interpolation
use_averaged = True  # Prefer averaged data if available
if use_averaged and freq_averaged_data is not None:
    # Use frequency-averaged data (preferred for smoother solutions)
    model_freq_ghz = freq_averaged_data['freq']
    model_rho_rad = np.deg2rad(freq_averaged_data['rho_deg'])
    model_solution = freq_averaged_data['solution']
    
    # Filter to good (non-outlier) data
    good_mask = (model_solution != 'outlier') & np.isfinite(model_rho_rad)
    
    print(f"Using frequency-averaged data: {np.sum(good_mask)}/{len(model_freq_ghz)} good bins")
else:
    # Use per-channel data
    model_freq_ghz = f_use
    model_rho_rad = selected_rho_per_channel
    model_solution = selected_solution_per_channel
    
    # Filter to good (non-outlier) data
    good_mask = (model_solution != 'outlier') & np.isfinite(model_rho_rad)
    
    print(f"Using per-channel data: {np.sum(good_mask)}/{len(model_freq_ghz)} good channels")

if np.sum(good_mask) < 2:
    print("ERROR: Insufficient good data points for interpolation (need at least 2)")
    print("XF table created but not populated with computed solutions")
else:
    # Extract good data for interpolation
    freq_good = model_freq_ghz[good_mask]
    rho_good = model_rho_rad[good_mask]
    
    freq_min, freq_max = freq_good.min(), freq_good.max()
    print(f"  Model frequency range: {freq_min:.3f} - {freq_max:.3f} GHz")
    
    # Convert to complex representation (handles 2π wrapping)
    complex_good = np.exp(1j * rho_good)
    real_good = np.real(complex_good)
    imag_good = np.imag(complex_good)
    
    # ===== Apply smoothing if enabled =====
    if XF_USE_SMOOTHING and len(freq_good) >= 5:
        # Auto-calculate window length if not specified
        if XF_SAVGOL_WINDOW is None:
            window_length = min(11, len(freq_good) if len(freq_good) % 2 == 1 else len(freq_good) - 1)
            if window_length < 5:
                window_length = 5 if len(freq_good) >= 5 else len(freq_good)
                if window_length % 2 == 0:
                    window_length -= 1
        else:
            window_length = XF_SAVGOL_WINDOW
            # Ensure window is odd and valid
            if window_length % 2 == 0:
                window_length += 1
            window_length = min(window_length, len(freq_good))
        
        polyorder = min(XF_SAVGOL_POLYORDER, window_length - 1)
        
        # Apply Savitzky-Golay filter to real and imaginary parts
        real_smooth = savgol_filter(real_good, window_length, polyorder)
        imag_smooth = savgol_filter(imag_good, window_length, polyorder)
        
        # Reconstruct smoothed complex gains
        complex_smooth = real_smooth + 1j * imag_smooth
        complex_smooth = complex_smooth / np.abs(complex_smooth)  # Normalize to unit magnitude
        
        print(f"  Applied Savitzky-Golay smoothing (window={window_length}, poly={polyorder})")
    else:
        # No smoothing - use original data
        complex_smooth = complex_good
        if XF_USE_SMOOTHING:
            print(f"  Smoothing disabled: insufficient points ({len(freq_good)} < 5)")
        else:
            print(f"  No smoothing applied (XF_USE_SMOOTHING=False)")
    
    # ===== Interpolate onto original data frequencies for final correction =====
    # Interpolate the final cross-hand phases back onto the original f_use grid
    real_for_data = np.real(complex_smooth)
    imag_for_data = np.imag(complex_smooth)
    
    real_interp_data = np.interp(f_use, freq_good, real_for_data)
    imag_interp_data = np.interp(f_use, freq_good, imag_for_data)
    
    complex_interp_data = real_interp_data + 1j * imag_interp_data
    complex_interp_data = complex_interp_data / np.abs(complex_interp_data)
    
    # Extract interpolated phases in radians
    rho_final_interp = np.angle(complex_interp_data)
    
    # Apply final corrections to Q, U, V
    print("\n=== Applying final corrections to Q, U, V ===")
    
    # Apply cross-hand phase correction
    U_corrected, V_corrected = correct_crosshand_phase(U_use, V_use, rho_final_interp)
    
    # Apply parallactic angle correction
    Q_corrected, U_corrected = correct_parallactic_angle(Q_use, U_corrected, chi_deg)
    
    print(f"  Applied cross-hand phase correction (interpolated from model)")
    print(f"  Applied parallactic angle correction (χ = {chi_deg:.3f}°)")
    
    # Calculate corrected polarization statistics
    P_corrected = np.sqrt(Q_corrected**2 + U_corrected**2)
    P_frac_corrected = np.nanmedian(P_corrected / I_use) * 100
    V_frac_corrected = np.nanmedian(np.abs(V_corrected) / I_use) * 100
    
    print(f"  Corrected polarization: P = {P_frac_corrected:.3f}%, |V| = {V_frac_corrected:.3f}%")

    # Plot final corrected I, Q, U, V spectra
    plot_stokes_spectra(f, I_use, Q_corrected, U_corrected, V_corrected, xfdir, 'iquv_post_XF.png', title='Post-XF Stokes Parameters')
        
    # Save corrected spectra to file
    np.savez(os.path.join(xfdir, 'final_corrected_stokes_spectra.npz'),
             freq_ghz=f_use,
             I=I_use,
             Q_corrected=Q_corrected,
             U_corrected=U_corrected,
             V_corrected=V_corrected,
             P_corrected=P_corrected,
             rho_applied=rho_final_interp)
    
    print(f"  Saved: {xfdir}/final_corrected_stokes_spectra.npz")
    
    # ===== Interpolate onto XF table frequencies =====
    # Extract phases for XF table interpolation
    real_for_interp = np.real(complex_smooth)
    imag_for_interp = np.imag(complex_smooth)
    
    # Interpolate real and imaginary parts separately
    real_interp = np.interp(chan_freq_xf_ghz, freq_good, real_for_interp)
    imag_interp = np.interp(chan_freq_xf_ghz, freq_good, imag_for_interp)
    
    # Reconstruct complex gains
    complex_interp = real_interp + 1j * imag_interp
    
    # Normalize to unit magnitude (phase-only correction)
    complex_interp = complex_interp / np.abs(complex_interp)
    
    # NEGATE THE PHASES because CASA applies them with wrong sign convention
    complex_interp_casa = np.conj(complex_interp)  # Complex conjugate negates the phase
    
    print(f"  Negated phases for CASA convention (complex conjugate applied)")
    
    # Determine which channels are outside valid frequency range (extrapolated)
    valid_range = (chan_freq_xf_ghz >= freq_min) & (chan_freq_xf_ghz <= freq_max)
    n_good = np.sum(valid_range)
    n_flagged_outofrange = np.sum(~valid_range)
    
    print(f"  XF table channels: {n_good} interpolated, {n_flagged_outofrange} out of range")
    
    if n_flagged_outofrange > 0:
        flagged_freqs = chan_freq_xf_ghz[~valid_range]
        print(f"  Out-of-range frequency spans: {flagged_freqs.min():.3f} - {flagged_freqs.max():.3f} GHz")
    
    # Open XF table and read existing data
    tb.open(xftab, nomodify=False)
    gains = tb.getcol('CPARAM')        # Shape: (1, n_freq_solutions, n_antennas)
    flag_original = tb.getcol('FLAG')  # Shape: (1, n_freq_solutions, n_antennas)
    
    print(f"  XF CPARAM shape: {gains.shape}")
    print(f"  Original flags: {np.sum(flag_original)} / {flag_original.size} flagged")
    
    # Tile interpolated gains to all antennas (cross-hand phase is antenna-independent)
    n_antennas = gains.shape[2]
    gains_new = np.tile(complex_interp_casa.reshape(1, -1, 1), (1, 1, n_antennas))
    
    # Create new flag array, preserving original flags
    flag_new = flag_original.copy()
    
    # Add flags for out-of-range channels
    if n_flagged_outofrange > 0:
        flag_new[0, ~valid_range, :] = True
    
    # Set gains to unity (1+0j) for ALL flagged channels (CASA standard)
    # Flagged channels get identity gains (no correction applied)
    gains_new[flag_new] = 1.0 + 0.0j
    
    print(f"  New flags: {np.sum(flag_new)} / {flag_new.size} flagged")
    print(f"  Added {np.sum(flag_new) - np.sum(flag_original)} new flags for out-of-range channels")
    print(f"  Set {np.sum(flag_new)} gain values to unity (all flagged channels)")
    
    # Write back to table
    tb.putcol('CPARAM', gains_new)
    tb.putcol('FLAG', flag_new)
    tb.flush()
    tb.close()
    
    print(f"✓ XF table populated with negated interpolated cross-hand phases: {xftab}")

shutil.rmtree(temp_ms)
print("\n=== Script Complete ===")
