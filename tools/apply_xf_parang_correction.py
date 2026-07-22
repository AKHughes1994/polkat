# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

"""
Post-hoc Xf (cross-hand phase) + parallactic angle correction for per-source
IQUV spectrum files, with Monte Carlo error propagation.

Iterates through RESULTS for RMSYNTH per-source '*_iquv.txt' files, applies:
  1. Cross-hand phase correction:  (U,V) rotated by +rho, rho = angle(Xf gain)
     [identical convention to manual_XF_solver.correct_crosshand_phase]
  2. Parallactic angle correction: (Q,U) rotated feed -> sky by -2*chi
     [identical convention to manual_XF_solver.correct_parallactic_angle]
     chi computed at the epoch middle MJD via the CASA AZ/EL method using the
     source position from the file header (no MS required).

Xf solutions are drawn from the perScan Xf table produced by
global_crosshand_phase.py, restricted to PA-calibrator fields whose names
contain any of XF_FIELD_IDENTIFIERS (case-insensitive substring match).

ERROR PROPAGATION (Monte Carlo, default 1000 iterations):
  Per channel, per iteration:
    - draw one (field, scan) Xf solution, categorically weighted by the
      cross-hand flux S/N at that channel (weights from the xf_fluxspec_*
      files written by global_crosshand_phase.py; noise via pooled
      first-difference window MADs — immune to the RM oscillation gradient).
      With a single available solution the sampling is degenerate (no-op),
    - jitter the drawn phase by that solution's own statistical uncertainty
      sigma_phi ~ 1/SNR rad (ensemble circular-std fallback when no S/N info
      exists; floored at 0.5 deg, capped at 30 deg) — this propagates each
      solution's error and smooths the K-lumpy solution mixture,
    - perturb I, Q, U, V by their per-channel Gaussian errors,
    - apply the Xf rotation, then the parallactic angle rotation,
    - accumulate corrected I, Q, U, V and P = sqrt(Q^2 + U^2).
  Adopt the per-channel MEDIAN as the value and the 16th/84th percentiles as
  asymmetric errors.  The standard rms_* columns carry the LARGER of the two
  one-sided errors (conservative); the full asymmetric pairs are appended as
  extra columns.  An asymmetry summary is printed per Stokes.

OUTPUTS (per input file):
  RESULTS/<prefix>_<src>_<ident>_iquv_XFcorr.txt      corrected spectrum
  RESULTS/fitting_plots/<prefix>_<src>_IQUV_spectrum_XFcorr.png
      plot matching RMSYNTH_01's plot_stokes_spectrum styling exactly, with
      translucent P = sqrt(Q^2+U^2) overlaid on the Q and U panels.
      (MFS reference lines and the spectral-index fit are omitted: neither
       the MFS fluxes nor the fit configuration live in the txt files.)

Files already ending in the correction suffix are skipped, as are I-only
files (no Q/U/V columns).

RUNTIME: plain Python with python-casacore (pyrap) — no CASA installation or
casatools instance is required.  Table access uses casacore.tables and the
parallactic angle uses casacore.measures (legacy pyrap namespaces are
supported as a fallback).
"""

import argparse
import glob
import os
import re
import sys
import datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from oxkat import config as cfg

# --- casacore (pyrap) backend: plain-python, no CASA installation needed ---
# Prefers the modern python-casacore namespace, falls back to legacy pyrap.
try:
    from casacore.tables import table as casatable
    from casacore.measures import measures as _measures_ctor
    from casacore import quanta as dq
except ImportError:                                    # legacy pyrap installs
    from pyrap.tables import table as casatable
    from pyrap.measures import measures as _measures_ctor
    from pyrap import quanta as dq

dm = _measures_ctor()
# python-casacore names the frame setter do_frame; old pyrap used doframe
_do_frame = getattr(dm, 'do_frame', None) or getattr(dm, 'doframe')


def _qval(q, unit):
    """Value of a casacore quantity (dict) in the requested unit."""
    return dq.quantity(q).get_value(unit)


# =============================================================================
# USER GLOBALS
# =============================================================================

# PA-calibrator field identifiers: an Xf table field contributes solutions if
# its name contains ANY of these substrings (case-insensitive).
XF_FIELD_IDENTIFIERS = ['1424', '1733', '1331', '3c286']

N_MC_DEFAULT   = 1000
SEED_DEFAULT   = 42
CORR_SUFFIX    = '_XFcorr'

# A data channel can use a given Xf solution only if its nearest valid
# (unflagged) Xf table channel is within MAX_GAP_FACTOR x the median Xf
# channel spacing — prevents interpolating a solution across wide flagged gaps.
MAX_GAP_FACTOR = 2.0

# Warn if the parallactic angle swings by more than this over the epoch
# (single mid-time chi is applied to the whole file).
CHI_SWING_WARN_DEG = 5.0


def msg(txt=''):
    ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'{ts} | {txt}', flush=True)


# =============================================================================
# Correction conventions — verbatim from manual_XF_solver.py
# =============================================================================

def correct_crosshand_phase(u_prime, v_prime, rho):
    """U_corr = U'cos(rho) - V'sin(rho);  V_corr = U'sin(rho) + V'cos(rho)."""
    cos_rho = np.cos(rho)
    sin_rho = np.sin(rho)
    u_corrected = u_prime * cos_rho - v_prime * sin_rho
    v_corrected = u_prime * sin_rho + v_prime * cos_rho
    return u_corrected, v_corrected


def correct_parallactic_angle(q, u, parang_deg):
    """Rotate Q,U from feed frame to sky frame by -2*chi."""
    parang_rad = np.radians(-2 * parang_deg)
    q_corrected = q * np.cos(parang_rad) + u * np.sin(parang_rad)
    u_corrected = -q * np.sin(parang_rad) + u * np.cos(parang_rad)
    return q_corrected, u_corrected


def compute_parang_from_header(ra_deg, dec_deg, mjd_days,
                               observatory='MeerKAT'):
    """Parallactic angle (deg) at mjd_days for (ra, dec) using the CASA
    AZ/EL method — same formula as manual_XF_solver.compute_parallactic_angle,
    but driven by the file-header position instead of an MS, and evaluated
    with casacore measures (no CASA instance required)."""
    ra_deg = float(ra_deg) % 360.0
    pos_meas = dm.observatory(observatory)
    pos_wgs  = dm.measure(pos_meas, 'WGS84')
    lat_rad  = _qval(pos_wgs['m1'], 'rad')

    _do_frame(pos_meas)
    _do_frame(dm.epoch('utc', dq.quantity(mjd_days, 'd')))
    dirm = dm.direction('J2000',
                        dq.quantity(ra_deg, 'deg'),
                        dq.quantity(dec_deg, 'deg'))
    azel = dm.measure(dirm, 'AZEL')

    az_rad = _qval(azel['m0'], 'rad')
    el_rad = _qval(azel['m1'], 'rad')

    num = -np.sin(az_rad)
    den = np.tan(lat_rad) * np.cos(el_rad) - np.cos(az_rad) * np.sin(el_rad)
    chi_deg = np.degrees(np.arctan2(num, den))
    return ((chi_deg + 180.0) % 360.0) - 180.0


# =============================================================================
# Input parsing
# =============================================================================

def parse_iquv_file(fpath):
    """Parse an RMSYNTH '*_iquv.txt' file.

    Returns (header_dict, data_dict) or (None, None) if the file lacks
    polarisation columns.  Data arrays are in mJy, freq in GHz.
    """
    hdr = {'source': None, 'ra_deg': None, 'dec_deg': None, 'ms': None,
           'start_mjd': None, 'end_mjd': None, 'mid_mjd': None,
           'header_lines': []}
    rows = []
    ncol_expected = 12  # Channel + 11 numeric columns

    with open(fpath, 'r') as fh:
        for line in fh:
            raw = line.rstrip('\n')
            s = raw.strip()
            if s.startswith('#'):
                hdr['header_lines'].append(raw)
                m = re.match(r'#\s*Source:\s*(\S+)', s)
                if m:
                    hdr['source'] = m.group(1)
                m = re.match(r'#\s*Position:\s*RA=([-\d.]+)\s*deg,\s*Dec=([-\d.]+)\s*deg', s)
                if m:
                    hdr['ra_deg']  = float(m.group(1))
                    hdr['dec_deg'] = float(m.group(2))
                m = re.match(r'#\s*MS:\s*(\S+)', s)
                if m:
                    hdr['ms'] = m.group(1)
                m = re.match(r'#\s*Start MJD:\s*([\d.]+)', s)
                if m:
                    hdr['start_mjd'] = float(m.group(1))
                m = re.match(r'#\s*End MJD:\s*([\d.]+)', s)
                if m:
                    hdr['end_mjd'] = float(m.group(1))
                m = re.match(r'#\s*Middle MJD:\s*([\d.]+)', s)
                if m:
                    hdr['mid_mjd'] = float(m.group(1))
                continue
            if not s:
                continue
            parts = s.split()
            if len(parts) < ncol_expected:
                return None, None   # I-only or unknown format — skip file
            try:
                rows.append([float(p) for p in parts[:ncol_expected]])
            except ValueError:
                continue

    if not rows:
        return None, None

    arr = np.array(rows, dtype=float)
    data = {
        'chan':     arr[:, 0].astype(int),
        'freq_ghz': arr[:, 1],
        'I':  arr[:, 2],  'Q':  arr[:, 3],  'U':  arr[:, 4],  'V': arr[:, 5],
        'P':  arr[:, 6],
        'rms_I': arr[:, 7], 'rms_Q': arr[:, 8], 'rms_U': arr[:, 9],
        'rms_V': arr[:, 10], 'rms_P': arr[:, 11],
    }
    return hdr, data


# =============================================================================
# Xf solutions + cross-hand S/N weights
# =============================================================================

def load_xf_solutions(xf_table_path, identifiers):
    """Read the perScan Xf table and return per-(field, scan) solutions for
    fields matching the identifier list.

    Returns
    -------
    solutions : dict {(field_name, scan): (phasor complex (nchan,), valid bool (nchan,))}
        Phasor = circular mean of corr-0 gains across unflagged antenna rows.
    freq_ghz_xf : ndarray — Xf table channel frequencies (GHz).
    """
    idents = [i.lower() for i in identifiers]

    t = casatable(os.path.join(xf_table_path, 'FIELD'), ack=False)
    field_names = list(t.getcol('NAME'))
    t.close()

    t = casatable(os.path.join(xf_table_path, 'SPECTRAL_WINDOW'), ack=False)
    freq_ghz_xf = np.asarray(t.getcol('CHAN_FREQ')).flatten() / 1.0e9
    t.close()

    t = casatable(xf_table_path, ack=False)
    # casacore returns (nrows, nchan, ncorr); transpose to the CASA-tool
    # convention (ncorr, nchan, nrows) so the indexing below matches the
    # sibling scripts.
    cparam = np.asarray(t.getcol('CPARAM')).T
    flags  = np.asarray(t.getcol('FLAG')).T
    fids   = np.asarray(t.getcol('FIELD_ID'))
    scans  = (np.asarray(t.getcol('SCAN_NUMBER'))
              if 'SCAN_NUMBER' in t.colnames() else None)
    t.close()

    n_rows = cparam.shape[2]
    matched_fields = {}
    for row in range(n_rows):
        fid = int(fids[row])
        fn  = field_names[fid] if fid < len(field_names) else None
        if fn is None or not any(i in fn.lower() for i in idents):
            continue
        sn  = int(scans[row]) if scans is not None else 0
        matched_fields.setdefault((fn, sn), []).append(row)

    if not matched_fields:
        raise RuntimeError(
            f'No Xf table fields match identifiers {identifiers}. '
            f'Table fields: {field_names}')

    solutions = {}
    for (fn, sn), rows in sorted(matched_fields.items()):
        sum_ph = np.zeros(cparam.shape[1], dtype=complex)
        counts = np.zeros(cparam.shape[1], dtype=int)
        for row in rows:
            good = ~flags[0, :, row]
            ph   = np.exp(1j * np.angle(cparam[0, :, row]))
            sum_ph[good] += ph[good]
            counts[good] += 1
        valid  = counts > 0
        phasor = np.where(valid, sum_ph / np.maximum(counts, 1), np.nan + 0j)
        # renormalise circular mean to unit amplitude where defined
        amp = np.abs(phasor)
        phasor = np.where(valid & (amp > 0), phasor / np.where(amp > 0, amp, 1.0),
                          np.nan + 0j)
        solutions[(fn, sn)] = (phasor, valid)
        msg(f'  Xf solution [{fn} scan {sn}]: '
            f'{int(valid.sum())}/{len(valid)} valid channels '
            f'({len(rows)} antenna rows pooled)')

    return solutions, freq_ghz_xf


def _parse_fluxspec_file(fpath):
    """Parse an xf_fluxspec_* file, locating the cross-hand flux column by
    HEADER NAME (robust to both the legacy 4/5-column and current 7-column
    IQUV formats).  Returns (freq_ghz, crosshand_flux_jy) or (None, None)."""
    col_names = None
    rows = []
    with open(fpath, 'r') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('#'):
                if 'freq_GHz' in s:
                    col_names = s.lstrip('#').split()
                continue
            if not s:
                continue
            parts = s.split()
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    if col_names is None or not rows:
        return None, None
    try:
        i_freq = col_names.index('freq_GHz')
        i_flux = col_names.index('crosshand_flux_Jy')
    except ValueError:
        return None, None
    arr = np.array([r for r in rows if len(r) >= len(col_names)], dtype=float)
    if arr.size == 0:
        return None, None
    return arr[:, i_freq], arr[:, i_flux]


def load_crosshand_weights(results_dir, solutions):
    """For each (field, scan) solution, find the matching xf_fluxspec_* file
    and return {key: (freq_ghz, snr)} where snr = crosshand flux / noise.

    Noise per field: pooled first-difference window MADs across that field's
    scans (same estimator as global_crosshand_phase._flag_sn_threshold).
    Solutions without a fluxspec file get None (uniform weight downstream).
    """
    # Locate candidate files per (field, scan)
    spec_map = {}
    for (fn, sn) in solutions:
        safe_fn = fn.replace('/', '_').replace(' ', '_')
        pattern = os.path.join(
            results_dir, f'xf_fluxspec_*_{safe_fn}_IQUV_scan{sn:04d}.txt')
        hits = sorted(glob.glob(pattern))
        if not hits:
            msg(f'  WARNING: no fluxspec file for [{fn} scan {sn}] '
                f'(pattern {os.path.basename(pattern)}) — uniform weight')
            spec_map[(fn, sn)] = None
            continue
        if len(hits) > 1:
            msg(f'  WARNING: {len(hits)} fluxspec files match [{fn} scan {sn}]; '
                f'using {os.path.basename(hits[-1])}')
        freq, flux = _parse_fluxspec_file(hits[-1])
        spec_map[(fn, sn)] = None if freq is None else (freq, flux)
        if freq is None:
            msg(f'  WARNING: could not parse {os.path.basename(hits[-1])} — '
                f'uniform weight')

    # Per-field noise via diff-MAD windows
    field_specs = {}
    for (fn, sn), v in spec_map.items():
        if v is not None:
            field_specs.setdefault(fn, []).append(v)
    field_noise = {}
    for fn, specs in field_specs.items():
        all_mads = []
        for freq, flux in specs:
            n = len(flux)
            win = max(4, n // 100)
            for start in range(0, n, win):
                seg   = flux[start:start + win]
                valid = seg[np.isfinite(seg)]
                if len(valid) >= 4:
                    d = np.diff(valid)
                    all_mads.append(np.median(np.abs(d - np.median(d))) / np.sqrt(2.0))
        if all_mads:
            noise = float(np.median(all_mads)) * 1.4826
            if noise > 0:
                field_noise[fn] = noise
                msg(f'  Field "{fn}": cross-hand noise = {noise:.4e} Jy '
                    f'({len(all_mads)} diff-MAD windows, {len(specs)} scan(s))')

    weights = {}
    for key, v in spec_map.items():
        fn = key[0]
        if v is None or fn not in field_noise:
            weights[key] = None
        else:
            freq, flux = v
            snr = np.where(np.isfinite(flux), flux / field_noise[fn], 0.0)
            weights[key] = (freq, np.maximum(snr, 0.0))
    return weights


def interp_solution(phasor, valid, freq_xf, freq_data, max_gap_factor):
    """Interpolate one Xf solution's phase onto the data frequency grid.

    Interpolation is done in the phasor (re/im) domain over valid channels
    only.  Data channels farther than max_gap_factor x median Xf channel
    spacing from the nearest valid Xf channel are marked unavailable.
    Returns (phase_rad (n_data,), available bool (n_data,)).
    """
    v_freq = freq_xf[valid]
    if len(v_freq) < 2:
        return np.zeros(len(freq_data)), np.zeros(len(freq_data), dtype=bool)
    v_ph = phasor[valid]
    re_i = np.interp(freq_data, v_freq, np.real(v_ph))
    im_i = np.interp(freq_data, v_freq, np.imag(v_ph))
    phase = np.arctan2(im_i, re_i)

    spacing = np.median(np.diff(np.sort(freq_xf)))
    max_gap = max_gap_factor * spacing
    nearest = np.min(np.abs(freq_data[:, None] - v_freq[None, :]), axis=1)
    inside  = (freq_data >= v_freq.min() - max_gap) & \
              (freq_data <= v_freq.max() + max_gap)
    return phase, (nearest <= max_gap) & inside


# =============================================================================
# Monte Carlo
# =============================================================================

# Per-solution phase-uncertainty bounds (radians): the statistical floor stops
# sigma collapsing to zero on very bright channels, the cap stops it exploding
# where the cross-hand S/N is marginal (those channels are dominated by the
# ensemble scatter anyway).
PHASE_SIGMA_FLOOR_RAD = np.deg2rad(0.5)
PHASE_SIGMA_CAP_RAD   = np.deg2rad(30.0)


def resolve_xf_sign(data, phase_mat, avail_mat, weight_mat, forced=None):
    """Determine the sign with which the table Xf phase must be applied.

    Whether the tabulated phase corresponds to the corruption (+rho, undo
    with -rho) or to the correction (+rho applied directly) is a convention
    that can differ between solvers and data paths.  The data settle it:
    corrupted linear polarisation leaks into Stokes V, and the correct sign
    collapses |V| toward the noise while the wrong sign inflates it.

    Uses the weight-averaged solution phase per channel, applies both signs,
    and scores each by the robust median of |V_corrected| / rms_V over valid
    channels.  Returns (sign, diagnostics dict).  forced = +1/-1 skips the
    choice but still reports both scores.
    """
    w = np.where(avail_mat, np.maximum(weight_mat, 1e-6), 0.0)
    wsum = w.sum(axis=0)
    ok = wsum > 0
    mean_phasor = np.zeros(phase_mat.shape[1], dtype=complex)
    mean_phasor[ok] = (w[:, ok] * np.exp(1j * phase_mat[:, ok])).sum(axis=0) / wsum[ok]
    mean_phase = np.angle(mean_phasor)

    valid = ok & np.isfinite(data['U']) & np.isfinite(data['V']) & \
            np.isfinite(data['rms_V']) & (data['rms_V'] > 0)

    scores = {}
    for s in (+1, -1):
        _, V_c = correct_crosshand_phase(data['U'], data['V'], s * mean_phase)
        scores[s] = float(np.median(np.abs(V_c[valid]) / data['rms_V'][valid]))
    score0 = float(np.median(np.abs(data['V'][valid]) / data['rms_V'][valid]))

    sign = forced if forced in (+1, -1) else (+1 if scores[+1] <= scores[-1] else -1)
    diag = {'uncorrected': score0, 'plus': scores[+1], 'minus': scores[-1],
            'forced': forced is not None, 'n_valid': int(valid.sum())}
    return sign, diag


def build_phase_sigma(phase_mat, avail_mat, weight_mat, has_snr):
    """Per-(solution, channel) phase uncertainty (rad) for the MC jitter.

    Each solved cross-hand phase carries its own statistical error
    sigma_phi ~ 1/SNR rad (phase error of a phasor at amplitude S/N).  Adding
    this jitter to each drawn phase (i) propagates the per-solution
    uncertainty the categorical sampling alone cannot see, and (ii) smooths
    the K-lumpy ensemble mixture so the output percentiles are not quantised
    to the K discrete solutions.

    Solutions without fluxspec S/N information (has_snr False) fall back to
    the per-channel circular standard deviation across the available
    solutions (the ensemble jitter itself), or the floor if only one
    solution exists.
    """
    K, n_chan = phase_mat.shape

    # Ensemble circular std per channel (fallback + single-solution case)
    ph_masked = np.where(avail_mat, np.exp(1j * phase_mat), 0.0)
    counts    = avail_mat.sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'):
        R = np.abs(ph_masked.sum(axis=0)) / np.maximum(counts, 1)
    ens_sigma = np.sqrt(np.maximum(-2.0 * np.log(np.clip(R, 1e-6, 1.0)), 0.0))
    ens_sigma = np.where(counts >= 2, ens_sigma, PHASE_SIGMA_FLOOR_RAD)

    sigma_mat = np.empty((K, n_chan))
    for k in range(K):
        if has_snr[k]:
            with np.errstate(divide='ignore'):
                sigma_mat[k] = 1.0 / np.maximum(weight_mat[k], 1e-3)
        else:
            sigma_mat[k] = ens_sigma
    return np.clip(sigma_mat, PHASE_SIGMA_FLOOR_RAD, PHASE_SIGMA_CAP_RAD)


def run_monte_carlo(data, phase_mat, avail_mat, weight_mat, sigma_mat,
                    chi_deg, n_mc, seed):
    """Vectorised MC over (n_mc, n_chan).

    phase_mat, avail_mat, weight_mat, sigma_mat : (K, n_chan) — per-solution
        phase, availability, cross-hand S/N weight, and phase uncertainty
        (rad) on the data frequency grid.

    Each iteration draws a solution index per channel (categorical, S/N
    weighted), then jitters the drawn phase by its own sigma — so the output
    errors carry both the between-solution scatter and the per-solution
    statistical uncertainty.

    Returns dict of {stokes: (median, err_lo, err_hi)} for I, Q, U, V, P,
    plus 'n_sol' (solutions available per channel) and 'dead' mask.
    """
    rng    = np.random.default_rng(seed)
    n_chan = len(data['freq_ghz'])
    K      = phase_mat.shape[0]

    w = np.where(avail_mat, weight_mat, 0.0)
    # If every available solution has zero weight at a channel (e.g. missing
    # fluxspec), fall back to uniform weights over the available solutions.
    col_wsum = w.sum(axis=0)
    uniform  = (col_wsum <= 0) & (avail_mat.sum(axis=0) > 0)
    w[:, uniform] = avail_mat[:, uniform].astype(float)
    col_wsum = w.sum(axis=0)

    dead = col_wsum <= 0          # no usable Xf solution at this channel
    n_sol = avail_mat.sum(axis=0)

    # Categorical sampling per channel via inverse-CDF
    w_safe = np.where(dead[None, :], 1.0, w)
    cdf = np.cumsum(w_safe / w_safe.sum(axis=0, keepdims=True), axis=0)
    r   = rng.random((n_mc, n_chan))
    idx = (r[:, None, :] > cdf[None, :, :]).sum(axis=1)     # (n_mc, n_chan)
    phase_draws = np.take_along_axis(
        np.broadcast_to(phase_mat[None, :, :], (n_mc, K, n_chan)),
        idx[:, None, :], axis=1)[:, 0, :]
    sigma_draws = np.take_along_axis(
        np.broadcast_to(sigma_mat[None, :, :], (n_mc, K, n_chan)),
        idx[:, None, :], axis=1)[:, 0, :]
    # Jitter each drawn phase by that solution's own statistical uncertainty
    phase_draws = phase_draws + rng.standard_normal((n_mc, n_chan)) * sigma_draws

    def draw(key):
        return data[key][None, :] + \
               rng.standard_normal((n_mc, n_chan)) * data['rms_' + key][None, :]

    I_d = draw('I')
    Q_d = draw('Q')
    U_d = draw('U')
    V_d = draw('V')

    U_x, V_x = correct_crosshand_phase(U_d, V_d, phase_draws)
    Q_c, U_c = correct_parallactic_angle(Q_d, U_x, chi_deg)
    V_c = V_x
    P_c = np.hypot(Q_c, U_c)

    out = {'n_sol': n_sol, 'dead': dead, 'chi_deg': chi_deg}
    for name, mat in (('I', I_d), ('Q', Q_c), ('U', U_c),
                      ('V', V_c), ('P', P_c)):
        p16, p50, p84 = np.percentile(mat, [15.865, 50.0, 84.135], axis=0)
        med    = np.where(dead, np.nan, p50)
        err_lo = np.where(dead, np.nan, p50 - p16)
        err_hi = np.where(dead, np.nan, p84 - p50)
        out[name] = (med, err_lo, err_hi)

    # Circular statistics of the sampled + jittered Xf phase itself — the
    # actual per-channel distribution the corrections were drawn from
    # (mixture over solutions + per-solution 1/SNR jitter).  Percentiles are
    # taken on wrapped deviations about the circular mean so they are safe
    # against the ±180 deg branch cut.
    ref = np.angle(np.mean(np.exp(1j * phase_draws), axis=0))
    dev = np.angle(np.exp(1j * (phase_draws - ref[None, :])))
    d16, d50, d84 = np.percentile(dev, [15.865, 50.0, 84.135], axis=0)
    ph_med = np.degrees(np.angle(np.exp(1j * (ref + d50))))
    out['phase'] = (np.where(dead, np.nan, ph_med),
                    np.where(dead, np.nan, np.degrees(d50 - d16)),
                    np.where(dead, np.nan, np.degrees(d84 - d50)))
    return out


def print_asymmetry_summary(mc, warn_frac=0.20):
    msg('  Asymmetry summary  [(hi - lo) / mean(hi, lo)]:')
    for name in ('I', 'Q', 'U', 'V', 'P'):
        med, lo, hi = mc[name]
        good = np.isfinite(lo) & np.isfinite(hi) & ((lo + hi) > 0)
        if not good.any():
            msg(f'    {name}: no valid channels')
            continue
        asym = (hi[good] - lo[good]) / (0.5 * (hi[good] + lo[good]))
        n_warn = int(np.sum(np.abs(asym) > warn_frac))
        msg(f'    {name}:  median = {np.median(asym):+.3f}   '
            f'max|.| = {np.max(np.abs(asym)):.3f}   '
            f'channels >{warn_frac*100:.0f}%: {n_warn}/{int(good.sum())}')


# =============================================================================
# Output: corrected txt + plot
# =============================================================================

def write_corrected_file(out_path, hdr, data, mc, xf_table, sol_keys,
                         n_mc, seed):
    med = {k: mc[k][0] for k in ('I', 'Q', 'U', 'V', 'P')}
    lo  = {k: mc[k][1] for k in ('I', 'Q', 'U', 'V', 'P')}
    hi  = {k: mc[k][2] for k in ('I', 'Q', 'U', 'V', 'P')}
    # Standard rms columns: LARGER of the two one-sided 68% errors (conservative)
    rms = {k: np.maximum(lo[k], hi[k]) for k in med}

    with open(out_path, 'w') as fh:
        for line in hdr['header_lines']:
            if line.strip().startswith('# Columns:'):
                break
            fh.write(line + '\n')
        fh.write('#\n')
        fh.write('# --- Xf + parallactic angle correction (this file) ---\n')
        fh.write(f'# Generated:        {datetime.datetime.now().isoformat()}\n')
        fh.write(f'# Xf table:         {xf_table}\n')
        fh.write('# Solutions used:   ' +
                 ', '.join(f'{fn}:scan{sn}' for fn, sn in sol_keys) + '\n')
        fh.write(f'# Parallactic chi:  {mc["chi_deg"]:+.4f} deg (epoch mid-time)\n')
        fh.write(f'# Xf sign applied:  {mc.get("xf_sign", +1):+d} '
                 '(auto-resolved by minimising median |V|/rms unless forced)\n')
        fh.write(f'# Monte Carlo:      {n_mc} iterations, seed {seed}\n')
        fh.write('# Phase jitter:     each drawn Xf phase perturbed by its own\n')
        fh.write('#                   sigma_phi ~ 1/SNR rad (ensemble circular\n')
        fh.write('#                   std fallback; floor 0.5 deg, cap 30 deg)\n')
        fh.write('# Values:           per-channel MC median; rms_* = larger of\n')
        fh.write('#                   the asymmetric 68% errors; *_lo/*_hi give\n')
        fh.write('#                   the (16th, 84th) percentile offsets.\n')
        fh.write('# Columns: Channel Freq[GHz] I[mJy] Q[mJy] U[mJy] V[mJy] '
                 'Plin[mJy] rms_I rms_Q rms_U rms_V rms_Plin '
                 'errI_lo errI_hi errQ_lo errQ_hi errU_lo errU_hi '
                 'errV_lo errV_hi errPlin_lo errPlin_hi\n')
        fh.write('#\n')
        for i in range(len(data['freq_ghz'])):
            vals = [f'{int(data["chan"][i]):4d}', f'{data["freq_ghz"][i]:10.4f}']
            vals += [f'{med[k][i]:12.4f}' for k in ('I', 'Q', 'U', 'V', 'P')]
            vals += [f'{rms[k][i]:10.4f}' for k in ('I', 'Q', 'U', 'V', 'P')]
            for k in ('I', 'Q', 'U', 'V', 'P'):
                vals += [f'{lo[k][i]:10.4f}', f'{hi[k][i]:10.4f}']
            fh.write('  '.join(vals) + '\n')
    msg(f'  Saved corrected spectrum: {out_path}')


def _mad_ylim(arr, err=None, k=10.0, pad=0.10):
    """Robust MAD-based y-axis limits — verbatim from RMSYNTH_01."""
    arr = np.asarray(arr, dtype=float)
    err = np.zeros_like(arr) if err is None else np.asarray(err, dtype=float)
    finite = np.isfinite(arr) & np.isfinite(err)
    if not finite.any():
        return (-1.0, 1.0)
    arr = arr[finite]
    err = err[finite]

    med_f   = np.median(arr)
    mad_f   = np.median(np.abs(arr - med_f))
    sigma_f = mad_f * 1.4826
    if sigma_f == 0:
        sigma_f = np.std(arr) if len(arr) > 1 else 1.0
    flux_mask = (arr >= med_f - k * sigma_f) & (arr <= med_f + k * sigma_f)

    med_e   = np.median(err)
    mad_e   = np.median(np.abs(err - med_e))
    sigma_e = mad_e * 1.4826
    if sigma_e == 0:
        sigma_e = np.std(err) if len(err) > 1 else 1.0
    err_mask = err <= med_e + k * sigma_e

    mask = flux_mask & err_mask
    if not mask.any():
        mask = np.ones(len(arr), dtype=bool)

    lo = np.min(arr[mask])
    hi = np.max(arr[mask])
    if lo == hi:
        lo -= 1.0
        hi += 1.0
    margin = pad * (hi - lo)
    return (lo - margin, hi + margin)


def plot_corrected_spectrum(plot_path, src_name, data, mc):
    """4x1 IQUV plot matching RMSYNTH_01.plot_stokes_spectrum styling, with
    asymmetric MC errorbars and translucent P on the Q and U panels."""
    freq = data['freq_ghz']
    P_med = mc['P'][0]

    stokes = [('I', 'tab:blue'), ('Q', 'tab:orange'),
              ('U', 'tab:green'), ('V', 'tab:red')]

    fig, axes = plt.subplots(4, 1, figsize=(12, 3.5 * 4),
                             sharex=True, squeeze=False)
    axes = axes[:, 0]

    for ax, (label, colour) in zip(axes, stokes):
        med, lo, hi = mc[label]

        ax.errorbar(freq, med, yerr=[lo, hi],
                    fmt='o', markersize=3, capsize=2, color=colour,
                    elinewidth=0.8, linewidth=0, label=f'Stokes {label}')

        if label in ('Q', 'U', 'V'):
            ax.axhline(0, color='grey', linestyle=':', linewidth=0.8)

        # Translucent linear polarised intensity on Q and U panels
        if label in ('Q', 'U'):
            ax.plot(freq, P_med, color='purple', alpha=0.30,
                    linewidth=1.2, label=r'$P=\sqrt{Q^2+U^2}$')

        sym_err = np.where(np.isfinite(lo) & np.isfinite(hi),
                           np.maximum(lo, hi), np.nan)
        ylim = _mad_ylim(med, err=sym_err)
        if label in ('Q', 'U'):
            p_hi = _mad_ylim(P_med)[1]
            ylim = (ylim[0], max(ylim[1], p_hi))
        ax.set_ylim(ylim)

        ax.set_ylabel(f'S (mJy) [Stokes {label}]', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.tick_params(direction='in', which='both', top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

    axes[-1].set_xlabel('Frequency (GHz)', fontsize=11)
    axes[0].set_title(f'{src_name}  —  Stokes IQUV Spectra (Xf + parang corrected)',
                      fontsize=12)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  Saved plot: {plot_path}')


def plot_xf_solutions(plot_path, src_name, freq_data, phase_mat, avail_mat,
                      sol_keys, mc):
    """Diagnostic plot of the Xf terms used for the correction: each
    contributing (field, scan) solution's phase, overlaid with the median and
    68% CI of the sampled + jittered distribution actually applied — makes
    the adopted Xf jitter (between-solution scatter + per-solution 1/SNR
    smoothing) directly visible."""
    ph_med, ph_lo, ph_hi = mc['phase']

    fig, ax = plt.subplots(1, 1, figsize=(12, 4.5))

    cmap = plt.cm.tab10
    for k, (fn, sn) in enumerate(sol_keys):
        ph = np.where(avail_mat[k], np.degrees(phase_mat[k]), np.nan)
        ax.plot(freq_data, ph, 'o', markersize=2.5, alpha=0.55,
                color=cmap(k % 10), label=f'{fn}: scan {sn}')

    ax.fill_between(freq_data, ph_med - ph_lo, ph_med + ph_hi,
                    color='k', alpha=0.18, linewidth=0,
                    label='MC 68% CI (mixture + jitter)')
    ax.plot(freq_data, ph_med, '-', color='k', linewidth=1.4,
            label='MC median')

    ax.set_xlabel('Frequency (GHz)', fontsize=11)
    ax.set_ylabel('Xf phase (deg)', fontsize=10)
    ax.set_title(f'{src_name}  —  Xf solutions and adopted MC distribution',
                 fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc='upper right')
    ax.tick_params(direction='in', which='both', top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    msg(f'  Saved Xf solution/jitter plot: {plot_path}')


# =============================================================================
# Per-file driver
# =============================================================================

def process_file(fpath, solutions, weights, freq_ghz_xf, xf_table,
                 n_mc, seed, max_gap_factor, make_xf_plot=False,
                 xf_sign_forced=None):
    base = os.path.basename(fpath)
    msg()
    msg(f'=== {base} ===')

    hdr, data = parse_iquv_file(fpath)
    if hdr is None:
        msg('  SKIP: not a full-Stokes iquv file (I-only or unknown format).')
        return False

    missing = [k for k in ('ra_deg', 'dec_deg', 'mid_mjd') if hdr[k] is None]
    if missing:
        msg(f'  SKIP: header missing {missing} — cannot compute parang.')
        return False

    src = hdr['source'] or base.split('_')[1]
    freq_data = data['freq_ghz']
    n_chan = len(freq_data)
    msg(f'  Source {src}: {n_chan} channels, '
        f'{freq_data.min():.4f}-{freq_data.max():.4f} GHz')

    # Parallactic angle at start/mid/end — apply mid, warn on large swing
    chi_mid = compute_parang_from_header(hdr['ra_deg'], hdr['dec_deg'],
                                         hdr['mid_mjd'])
    if hdr['start_mjd'] is not None and hdr['end_mjd'] is not None:
        chi_s = compute_parang_from_header(hdr['ra_deg'], hdr['dec_deg'],
                                           hdr['start_mjd'])
        chi_e = compute_parang_from_header(hdr['ra_deg'], hdr['dec_deg'],
                                           hdr['end_mjd'])
        swing = abs(((chi_e - chi_s + 180.0) % 360.0) - 180.0)
        msg(f'  Parallactic angle: chi = {chi_mid:+.3f} deg at mid-time '
            f'(start {chi_s:+.3f}, end {chi_e:+.3f}, swing {swing:.3f} deg)')
        if swing > CHI_SWING_WARN_DEG:
            msg(f'  WARNING: chi swings {swing:.2f} deg over the epoch but a '
                f'single mid-time value is applied — consider per-scan files.')
    else:
        msg(f'  Parallactic angle: chi = {chi_mid:+.3f} deg at mid-time')

    # Interpolate each solution + weight onto the data grid
    keys = sorted(solutions.keys())
    K = len(keys)
    phase_mat  = np.zeros((K, n_chan))
    avail_mat  = np.zeros((K, n_chan), dtype=bool)
    weight_mat = np.zeros((K, n_chan))
    has_snr    = np.zeros(K, dtype=bool)
    for k, key in enumerate(keys):
        phasor, valid = solutions[key]
        phase_mat[k], avail_mat[k] = interp_solution(
            phasor, valid, freq_ghz_xf, freq_data, max_gap_factor)
        wv = weights.get(key)
        if wv is None:
            weight_mat[k] = 1.0
        else:
            w_freq, w_snr = wv
            weight_mat[k] = np.interp(freq_data, w_freq, w_snr,
                                      left=0.0, right=0.0)
            has_snr[k] = True

    # Resolve the sign with which the table phase must be applied — the
    # correct sign collapses the leaked linear polarisation in Stokes V.
    xf_sign, sd = resolve_xf_sign(data, phase_mat, avail_mat, weight_mat,
                                  forced=xf_sign_forced)
    msg(f'  Xf sign resolution (median |V|/rms over {sd["n_valid"]} channels): '
        f'uncorrected = {sd["uncorrected"]:.2f}, '
        f'+rho = {sd["plus"]:.2f}, -rho = {sd["minus"]:.2f}  '
        f'-> applying sign {xf_sign:+d}'
        f'{" (forced by --xf-sign)" if sd["forced"] else ""}')
    if not sd['forced'] and min(sd['plus'], sd['minus']) > 0.8 * sd['uncorrected']:
        msg('  WARNING: neither sign reduces |V|/rms by >20% — the Xf table '
            'may not match this data, or the source has genuine circular '
            'polarisation at the leakage level.  Inspect before trusting.')
    phase_mat = xf_sign * phase_mat

    sigma_mat = build_phase_sigma(phase_mat, avail_mat, weight_mat, has_snr)

    mc = run_monte_carlo(data, phase_mat, avail_mat, weight_mat, sigma_mat,
                         chi_mid, n_mc, seed)
    mc['xf_sign'] = xf_sign

    # Confirm the correction on V explicitly (this is the whole point)
    _vok = np.isfinite(mc['V'][0]) & (data['rms_V'] > 0)
    v_pre  = np.median(np.abs(data['V'][_vok]) / data['rms_V'][_vok])
    v_post = np.median(np.abs(mc['V'][0][_vok]) / data['rms_V'][_vok])
    msg(f'  Stokes V check: median |V|/rms  {v_pre:.2f} (pre) -> '
        f'{v_post:.2f} (post-correction)')
    n_dead = int(mc['dead'].sum())
    msg(f'  MC complete: {n_mc} iterations, {K} candidate solution(s); '
        f'median available per channel: {int(np.median(mc["n_sol"]))}; '
        f'{n_dead} channel(s) with no usable solution (NaN in output).')
    print_asymmetry_summary(mc)

    # Outputs
    out_txt = fpath[:-len('_iquv.txt')] + f'_iquv{CORR_SUFFIX}.txt'
    write_corrected_file(out_txt, hdr, data, mc, xf_table, keys, n_mc, seed)

    plot_dir = os.path.join(os.path.dirname(fpath), 'fitting_plots')
    os.makedirs(plot_dir, exist_ok=True)
    prefix = base.split('_')[0]
    plot_path = os.path.join(
        plot_dir, f'{prefix}_{src}_IQUV_spectrum{CORR_SUFFIX}.png')
    plot_corrected_spectrum(plot_path, src, data, mc)

    # One-off diagnostic: Xf solutions + adopted MC phase distribution,
    # produced for the first file corrected in this run.
    if make_xf_plot:
        xf_plot_path = os.path.join(
            plot_dir, f'{prefix}_{src}_XFsolutions{CORR_SUFFIX}.png')
        plot_xf_solutions(xf_plot_path, src, freq_data,
                          phase_mat, avail_mat, keys, mc)

    return True


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Apply Xf + parallactic angle corrections to RMSYNTH '
                    'per-source iquv files with MC error propagation.')
    parser.add_argument('--xf-table', type=str, default='',
                        help='Path to the perScan Xf table.  Default: '
                             'auto-discover a single *_perScan.Xf in GAINTABLES.')
    parser.add_argument('--results', type=str, default=cfg.RESULTS,
                        help='Directory to scan for *_iquv.txt files.')
    parser.add_argument('--identifiers', nargs='*',
                        default=XF_FIELD_IDENTIFIERS,
                        help='PA-cal field identifier substrings '
                             '(case-insensitive).')
    parser.add_argument('--n-mc', type=int, default=N_MC_DEFAULT)
    parser.add_argument('--seed', type=int, default=SEED_DEFAULT)
    parser.add_argument('--max-gap-factor', type=float, default=MAX_GAP_FACTOR)
    parser.add_argument('--xf-sign', type=str, default='auto',
                        choices=['auto', '+1', '-1', '1'],
                        help='Sign with which the table Xf phase is applied. '
                             'auto (default): choose per file by minimising '
                             'median |V|/rms; +1/-1: force.')
    args, _ = parser.parse_known_args()

    xf_table = args.xf_table
    if not xf_table:
        hits = sorted(glob.glob(os.path.join(cfg.GAINTABLES, '*_perScan.Xf')))
        if len(hits) != 1:
            raise RuntimeError(
                f'Found {len(hits)} *_perScan.Xf tables in {cfg.GAINTABLES}; '
                f'specify one with --xf-table.  Hits: {hits}')
        xf_table = hits[0]
    msg(f'Xf table: {xf_table}')
    msg(f'Field identifiers: {args.identifiers}')

    solutions, freq_ghz_xf = load_xf_solutions(xf_table, args.identifiers)
    weights = load_crosshand_weights(args.results, solutions)

    pattern = os.path.join(args.results, '*_iquv.txt')
    files = sorted(f for f in glob.glob(pattern)
                   if not f.endswith(f'_iquv{CORR_SUFFIX}.txt'))
    if not files:
        msg(f'No *_iquv.txt files found in {args.results}.')
        return
    msg(f'Found {len(files)} iquv file(s) to correct.')

    xf_plot_pending = True
    for fpath in files:
        try:
            _forced = None if args.xf_sign == 'auto' else int(args.xf_sign)
            ok = process_file(fpath, solutions, weights, freq_ghz_xf, xf_table,
                              args.n_mc, args.seed, args.max_gap_factor,
                              make_xf_plot=xf_plot_pending,
                              xf_sign_forced=_forced)
            if ok:
                xf_plot_pending = False
        except Exception as e:
            msg(f'  ERROR processing {os.path.basename(fpath)}: {e}')

    msg()
    msg('Xf + parang correction complete.')


if __name__ == '__main__':
    main()
