# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
#
# Standalone Xf (global cross-hand phase) solver.
# Assumes K, B, Gp, G, and Df tables have already been produced.
# This script ONLY solves for Xf (no KCROSS).  It always produces:
#   [PREFIX]_combineScan.Xf      -- solved with combine='scan', solint='inf'
#   [PREFIX]_perScan.Xf          -- solved per scan (only if any PA cal field
#                                   has more than one scan; do_per_scan=True)
#   [PREFIX]_combineScanField.Xf -- circular phasor average of combineScan
#                                   across all PA cal fields (only when
#                                   do_multi_field=True, i.e. >1 PA cal field)
#
# Fields are ALWAYS handled separately through all flagging and statistics.
# The only cross-field step is the final combineScanField circular average.
#
# Post-processing pipeline per table:
#   1. +/-180 deg phase degeneracy resolution against a reference phase.
#   --- flagging steps (skipped for phasor-averaged tables; mad_only for
#       combineScanField) ---
#   2. Gross MAD pre-clip (per-field, XF_GROSS_CLIP sigma).
#   3. Per-field flag: channels where cross-hand fraction < XF_MIN_CROSS_FRAC.
#      Fraction is computed per (field, scan); the per-field MAD flag mask is
#      broadcast to all rows of that field.
#   4. Per-field, per-scan S/N threshold: field noise = median(all per-window
#      MADs from all scans of that field) * 1.4826; channels where any scan
#      has S/N < XF_MIN_SN are flagged.
#   5. Fine MAD clip (per-field polynomial baseline + sliding circular MAD).
#      Also applied to combineScanField after the cross-field average.
#   ---
#   6. Diagnostic PNG plot (phases + cross-hand fraction + S/N panels).
#      Y-axis limits are saved from the first plot produced and reused in
#      all subsequent plots for direct comparison.
#   7. Optional manual SPW:channel flagging (XF_FLAG_FREQS).
#
# Per-field-scan flux spectrum text files are saved to RESULTS.

import sys
import os
import datetime
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def msg(text=''):
    """Print a timestamped message to stdout, flushed immediately."""
    ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'{ts} | {text}', flush=True)

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())


# -----------------------------------------------------------------------
#  User-defined parameters (argparse)
#
#  myuvrange and myspw are drawn from config.py (CAL_1GC_UVRANGE /
#  CAL_1GC_FREQRANGE) and are not exposed here.
#  POLANG_MOD, XF_GLOBAL_SIGMA_CLIP (gross clip) and XF_MAD_SIGMA_CLIP
#  (local MAD clip) defaults are drawn from config.py.
# -----------------------------------------------------------------------

import argparse

_parser = argparse.ArgumentParser(
    description='Standalone Xf (global cross-hand phase) solver.')

_parser.add_argument(
    '--pacal-name',
    type=str,
    default='',
    help='PA calibrator name(s), comma-separated.  Empty string (default) '
         'triggers auto-discovery from MS scan intents.')
_parser.add_argument(
    '--polang-mod',
    nargs=4,
    type=float,
    default=POLANG_MOD,
    metavar=('I', 'Q', 'U', 'V'),
    help='IQUV Stokes model for the PA calibrator (Jy).')
_parser.add_argument(
    '--xf-reference-phase',
    type=float,
    default=0.0,
    metavar='DEG',
    help='Reference phase (deg) for resolving the +/-180 deg degeneracy.')
_parser.add_argument(
    '--xf-flag-freqs',
    nargs='*',
    default=[],
    metavar='SPW:CHAN',
    help='SPW:channel range(s) to flag in the Xf table after solving '
         '(e.g. "0:0~10").  Leave empty to skip.')
_parser.add_argument(
    '--xf-min-cross-frac',
    type=float,
    default=0.0015,
    metavar='FRAC',
    help='Min cross-hand fraction sqrt(U^2+V^2)/I; channels below this '
         'are flagged in the Xf table.')
_parser.add_argument(
    '--xf-gross-clip',
    type=float,
    default=XF_GLOBAL_SIGMA_CLIP,
    metavar='SIGMA',
    help='Sigma threshold for the single-pass gross-outlier pre-clip (run before '
         'fractional and S/N cuts; catches wildly wrong channels).')
_parser.add_argument(
    '--xf-poly-order',
    type=int,
    default=4,
    metavar='ORDER',
    help='Polynomial order for the global phase baseline fit.')
_parser.add_argument(
    '--xf-mad-window',
    type=int,
    default=0,
    metavar='NCHAN',
    help='Sliding channel window width for local MAD clipping.  0 (default) '
         'auto-sets to the nearest odd integer of nchan_xf / 10.')
_parser.add_argument(
    '--xf-mad-sigma',
    type=float,
    default=XF_MAD_SIGMA_CLIP,
    metavar='SIGMA',
    help='Sigma threshold for local circular MAD channel clipping.')
_parser.add_argument(
    '--xf-field-mad-sigma',
    type=float,
    default=3.0,
    metavar='SIGMA',
    help='Sigma threshold for the final fine MAD clip applied to the '
         'combineScanField table after cross-field averaging (default 3.0; '
         'tighter than the per-field clip because the averaged solution is '
         'smoother).')
_parser.add_argument(
    '--xf-use-poly',
    dest='xf_use_poly',
    action='store_true',
    default=True,
    help='Fit a complex-plane polynomial baseline before MAD clip; MAD then '
         'operates on residuals rather than raw phases (default True).')
_parser.add_argument(
    '--xf-no-poly',
    dest='xf_use_poly',
    action='store_false',
    help='Disable polynomial baseline; MAD operates directly on raw phases.')
_parser.add_argument(
    '--xf-verbose',
    dest='xf_verbose',
    action='store_true',
    default=False,
    help='Print per-scan MAD clipping table to the log.')
_parser.add_argument(
    '--xf-no-verbose',
    dest='xf_verbose',
    action='store_false',
    help='Suppress per-scan MAD clipping table output.')
_parser.add_argument(
    '--xf-min-sn',
    type=float,
    default=10.0,
    metavar='SN',
    help='Minimum cross-hand power S/N; channels below this are flagged in '
         'the Xf table.  0 disables the cut.')

_args, _ = _parser.parse_known_args()

pacal_name         = _args.pacal_name
POLANG_MOD         = _args.polang_mod
XF_REFERENCE_PHASE = _args.xf_reference_phase
XF_FLAG_FREQS      = _args.xf_flag_freqs
XF_FLAG_CHANS      = []   # set to '0.92~0.94GHz' or ['0.92~0.94GHz'] to override XF_FLAG_FREQS
XF_MIN_CROSS_FRAC  = _args.xf_min_cross_frac
XF_GROSS_CLIP      = _args.xf_gross_clip
XF_POLY_ORDER      = _args.xf_poly_order
XF_MAD_WINDOW      = _args.xf_mad_window
XF_MAD_SIGMA       = _args.xf_mad_sigma
XF_FIELD_MAD_SIGMA = _args.xf_field_mad_sigma
XF_USE_POLY        = _args.xf_use_poly
XF_VERBOSE         = _args.xf_verbose
XF_MIN_SN          = _args.xf_min_sn


# -----------------------------------------------------------------------
#  Derived parameters
# -----------------------------------------------------------------------

myuvrange = CAL_1GC_UVRANGE
myspw     = CAL_1GC_FREQRANGE


# -----------------------------------------------------------------------
#  Calibration table paths  (must already exist)
# -----------------------------------------------------------------------

ktab  = GAINTABLES + '/cal_1GC_' + myms + '.K'
bptab = GAINTABLES + '/cal_1GC_' + myms + '.B'
gptab = GAINTABLES + '/cal_1GC_' + myms + '.Gp'
gtab  = GAINTABLES + '/cal_1GC_' + myms + '.G'
ftab  = GAINTABLES + '/cal_1GC_' + myms + '.F'
dftab = GAINTABLES + '/cal_1GC_' + myms + '.Df'
# Xf output tables are defined in theain flow after scan-count determination


# -----------------------------------------------------------------------
#  Output directories
# -----------------------------------------------------------------------

# Diagnostic plots go into a dedicated subdirectory under GAINPLOTS
xf_plot_dir = GAINPLOTS + '/XF'
os.makedirs(xf_plot_dir, exist_ok=True)

# Per-source/scan flux spectrum files go into RESULTS
os.makedirs(RESULTS, exist_ok=True)


# -----------------------------------------------------------------------
#  Resolve PA calibrator field name(s)
#
#  If pacal_name is empty, auto-discover fields from the MS by inspecting
#  scan intents.  Fields whose intent(s) match any of the following
#  criteria (case-insensitive) are included:
#    - intent == 'UNKNOWN'         (exact)
#    - intent == 'TARGET'          (exact)
#    - intent starts with 'CALIBRATE_AMP'  (partial)
#    - intent starts with 'CALIBRATE_PHA'  (partial)
#  The result is a comma-separated string of matching field names.
# -----------------------------------------------------------------------

if pacal_name == '':
    msg('pacal_name is empty — auto-discovering PA calibrator fields from MS intents.')

    _matched_fields = []
    msmd.open(myms)
    try:
        _all_names = msmd.fieldnames()
        for _fid in range(len(_all_names)):
            _fname   = _all_names[_fid]
            _intents = msmd.intentsforfield(_fid)
            for _intent in _intents:
                _il = _intent.lower()
                if (_il == 'unknown' or
                        _il.startswith('calibrate_amp') or
                        _il.startswith('calibrate_pha')):
                    _matched_fields.append(_fname)
                    break   # one matching intent is enough per field
    finally:
        msmd.close()

    pacal_name = ','.join(_matched_fields)

    if pacal_name:
        msg(f'  Auto-discovered PA calibrator field(s): {pacal_name}')
    else:
        msg('WARNING: No fields matched the PA calibrator intent criteria. '
            'Xf solving will be skipped.')


# -----------------------------------------------------------------------
#  Helper functions
# -----------------------------------------------------------------------

# MJD epoch: MJD 0 = 1858-11-17 00:00:00 UTC.
# CASA TIME columns are stored in MJD seconds.
_MJD_EPOCH = datetime.datetime(1858, 11, 17, 0, 0, 0)


def mjd_sec_to_isot(mjd_sec):
    """Convert a CASA MJD-seconds value to an ISO 8601 string (ms precision).

    Parameters
    ----------
    mjd_sec : float
        Time in seconds since MJD epoch (as stored in CASA TIME columns).

    Returns
    -------
    str
        ISO 8601 string, e.g. '2023-03-15T12:34:56.789'.
    """
    dt = _MJD_EPOCH + datetime.timedelta(seconds=float(mjd_sec))
    # Format to millisecond precision
    return dt.strftime('%Y-%m-%dT%H:%M:%S.') + '{:03d}'.format(dt.microsecond // 1000)


def rad_to_hms(rad):
    """Convert right ascension in radians to an HH:MM:SS.SSS string.

    Parameters
    ----------
    rad : float
        RA in radians.  Wrapped to [0, 2*pi) before conversion.

    Returns
    -------
    str
    """
    hours = np.degrees(rad) / 15.0
    hours = hours % 24.0          # Ensure [0, 24)
    h = int(hours)
    m_float = (hours - h) * 60.0
    m = int(m_float)
    s = (m_float - m) * 60.0
    return '{:02d}:{:02d}:{:06.3f}'.format(h, m, s)


def rad_to_dms(rad):
    """Convert declination in radians to a +/-DD:MM:SS.SS string.

    Parameters
    ----------
    rad : float
        Dec in radians.

    Returns
    -------
    str
    """
    deg_total = np.degrees(rad)
    sign = '+' if deg_total >= 0.0 else '-'
    deg_total = abs(deg_total)
    d = int(deg_total)
    m_float = (deg_total - d) * 60.0
    m = int(m_float)
    s = (m_float - m) * 60.0
    return '{}{:02d}:{:02d}:{:05.2f}'.format(sign, d, m, s)


def time_average_iquv(field_name, scan_num):
    """Read visibilities for a single field-scan and average over time.

    Uses CASA's ms tool to select the specified field and scan, requests
    Stokes IQUV (CASA converts from native correlations on the fly), then
    forms a weighted mean over all rows (baselines x integrations) to
    produce per-channel Stokes spectra.

    Falls back to reading raw correlations (XX,XY,YX,YY) and computing
    Stokes manually if the Stokes selection fails.

    Parameters
    ----------
    field_name : str
        Field name as it appears in the MS FIELD table.
    scan_num : int or str
        Scan number to select.

    Returns
    -------
    I_chan, Q_chan, U_chan, V_chan : ndarray or None
        Per-channel real flux densities (Jy), length = nchan_native.
        None for all four if selection fails.
    weight_uv : ndarray or None
        Per-channel effective weight sum for the cross-hand (U/V)
        correlations after flagging.  Shape (nchan_native,).  Used as
        the binning weight so that flagged channels are naturally
        down-weighted when averaging to the Xf table resolution.
        None if selection fails.
    t_start, t_end : float or None
        Scan start and end times in MJD seconds.
        None if selection fails.
    """
    # ---- First attempt: request Stokes IQUV directly ----
    ms.open(myms)
    ms.selectinit(reset=True)
    ok = ms.msselect({'field': field_name, 'scan': str(scan_num)})
    if not ok:
        ms.close()
        msg(f'    WARNING: ms.msselect failed for field={field_name} scan={scan_num}')
        return None, None, None, None, None, None

    stokes_mode = False
    try:
        ms.selectpolarization(['I', 'Q', 'U', 'V'])
        d = ms.getdata(['corrected_data', 'flag', 'weight', 'time'])
        stokes_mode = True
    except Exception as e:
        msg(f'    NOTE: Stokes selectpolarization failed ({e}); '
              f'falling back to raw correlations.')
        d = ms.getdata(['corrected_data', 'flag', 'weight', 'time'])
    ms.close()

    if 'corrected_data' not in d or 'flag' not in d:
        msg(f'    WARNING: getdata returned no corrected_data for field={field_name} scan={scan_num}')
        return None, None, None, None, None, None

    data_arr  = d['corrected_data']
    raw_flag  = np.asarray(d['flag'])        # bool, shape varies
    time_arr  = d['time']                    # (nrow,)
    ncorr_read, nchan_read, nrow = data_arr.shape

    t_start = float(time_arr.min())
    t_end   = float(time_arr.max())

    # ---- Build per-correlation weights ----
    # Use the scalar WEIGHT column (ncorr, nrow) and broadcast to channels.
    Wraw = np.asarray(d.get('weight', np.ones((ncorr_read, nrow), dtype=float)))
    if Wraw.ndim == 1 and Wraw.shape[0] == nrow:
        # Single weight per row
        W = np.ones((ncorr_read, nchan_read, nrow), dtype=float)
        for c in range(ncorr_read):
            W[c] = Wraw[np.newaxis, :]
    elif Wraw.ndim == 2 and Wraw.shape[-1] == nrow:
        # (ncorr, nrow)
        W = np.zeros((ncorr_read, nchan_read, nrow), dtype=float)
        for c in range(min(ncorr_read, Wraw.shape[0])):
            W[c] = Wraw[c, np.newaxis, :]
    else:
        W = np.ones((ncorr_read, nchan_read, nrow), dtype=float)

    # ---- Expand flags to (ncorr, nchan, nrow) ----
    if raw_flag.ndim == 3 and raw_flag.shape == (ncorr_read, nchan_read, nrow):
        flags = raw_flag.astype(bool)
    elif raw_flag.ndim == 3 and raw_flag.shape[0] == 1:
        flags = np.broadcast_to(
            raw_flag.astype(bool), (ncorr_read, nchan_read, nrow)).copy()
    elif raw_flag.ndim == 2 and raw_flag.shape == (nchan_read, nrow):
        flags = np.broadcast_to(
            raw_flag[np.newaxis].astype(bool), (ncorr_read, nchan_read, nrow)).copy()
    else:
        msg(f'    WARNING: Unrecognised FLAG shape {raw_flag.shape}; assuming no flags.')
        flags = np.zeros((ncorr_read, nchan_read, nrow), dtype=bool)

    # ---- Weighted time-average: (ncorr, nchan, nrow) -> (ncorr, nchan) ----
    W_eff = np.where(flags, 0.0, W)
    den   = np.sum(W_eff, axis=-1)                               # (ncorr, nchan)
    num   = np.sum(W_eff * data_arr, axis=-1, dtype=np.complex128)  # (ncorr, nchan)

    with np.errstate(invalid='ignore', divide='ignore'):
        vis_avg = np.where(den > 0, num / den, np.nan + 0j)      # (ncorr, nchan)

    # ---- Extract cross-hand effective weight per channel ----
    # den[corr, chan] = sum of unflagged weights for that correlation plane.
    # For the cross-hand flux (sqrt(U^2+V^2)) we want the weight of the
    # U/V-sensitive correlations:
    #   Stokes mode  -> planes 2 (U) and 3 (V)
    #   Raw XY/YX    -> planes 1 (XY) and 2 (YX)
    #   2-corr       -> no cross-hand available; set to zero
    if stokes_mode and ncorr_read >= 4:
        weight_uv = 0.5 * (den[2] + den[3])   # mean of U, V weight sums
    elif (not stokes_mode) and ncorr_read >= 4:
        weight_uv = 0.5 * (den[1] + den[2])   # mean of XY, YX weight sums
    else:
        weight_uv = np.zeros(nchan_read, dtype=float)

    del d, data_arr, raw_flag, time_arr, W, W_eff, den, num

    # ---- Convert to Stokes IQUV ----
    if stokes_mode:
        # CASA already returned IQUV order
        I_chan = np.real(vis_avg[0])
        Q_chan = np.real(vis_avg[1])
        U_chan = np.real(vis_avg[2])
        V_chan = np.real(vis_avg[3])
    else:
        # Convert from linear feed correlations (XX,XY,YX,YY)
        if ncorr_read >= 4:
            XX = vis_avg[0]; XY = vis_avg[1]
            YX = vis_avg[2]; YY = vis_avg[3]
            I_chan = np.real((XX + YY) / 2.0)
            Q_chan = np.real((XX - YY) / 2.0)
            U_chan = np.real((XY + YX) / 2.0)
            # V from imaginary part of (XY - YX): sign follows IAU convention
            V_chan = np.imag((YX - XY) / 2.0)
        elif ncorr_read == 2:
            # Parallel hands only
            I_chan = np.real((vis_avg[0] + vis_avg[1]) / 2.0)
            Q_chan = np.real((vis_avg[0] - vis_avg[1]) / 2.0)
            U_chan = np.full(nchan_read, np.nan)
            V_chan = np.full(nchan_read, np.nan)
        else:
            I_chan = np.real(vis_avg[0])
            Q_chan = np.full(nchan_read, np.nan)
            U_chan = np.full(nchan_read, np.nan)
            V_chan = np.full(nchan_read, np.nan)

    del vis_avg
    return I_chan, Q_chan, U_chan, V_chan, weight_uv, t_start, t_end


def bin_channels_weighted(arr_native, weights_native, n_out):
    """Bin a native-resolution 1-D array to n_out channels using weights.

    Each output channel is the weighted mean of the input channels that
    fall within its proportional bin.  Channels whose weight is zero or
    whose value is NaN are excluded; if all inputs in a bin are excluded
    the output bin is NaN with zero weight.

    Using the effective weight sum (from the time-averaging step) as the
    bin weight naturally propagates flagging: channels where most
    baselines were flagged contribute little to the bin average, and
    completely flagged channels (weight=0 / value=NaN) contribute nothing.

    Parameters
    ----------
    arr_native : ndarray, shape (nchan_native,)
        Values to bin (e.g. cross-hand flux per native channel).
    weights_native : ndarray, shape (nchan_native,)
        Per-channel effective weight sums from the time-averaging step.
        Must be >= 0; zero means fully flagged.
    n_out : int
        Number of output channels.

    Returns
    -------
    result : ndarray, shape (n_out,)
        Weighted-mean binned values.  NaN where no valid data in bin.
    result_weights : ndarray, shape (n_out,)
        Sum of weights in each output bin (useful for downstream
        weighting or diagnostics).
    """
    n_in = len(arr_native)

    if n_in == n_out:
        # No binning required; zero-weight channels become NaN
        result = arr_native.copy().astype(float)
        result_w = weights_native.copy().astype(float)
        bad = (result_w <= 0.0) | ~np.isfinite(result)
        result[bad] = np.nan
        result_w[bad] = 0.0
        return result, result_w

    result   = np.full(n_out, np.nan, dtype=float)
    result_w = np.zeros(n_out, dtype=float)

    for i in range(n_out):
        i_start = int(round(i       * n_in / float(n_out)))
        i_end   = int(round((i + 1) * n_in / float(n_out)))
        i_end   = max(i_start + 1, i_end)

        chunk_vals = arr_native[i_start:i_end]
        chunk_wts  = weights_native[i_start:i_end]

        # Valid: finite value AND positive weight
        valid = np.isfinite(chunk_vals) & (chunk_wts > 0.0)
        if not np.any(valid):
            continue

        w_sum = float(np.sum(chunk_wts[valid]))
        result[i]   = float(np.sum(chunk_wts[valid] * chunk_vals[valid])) / w_sum
        result_w[i] = w_sum

    return result, result_w


def flux_textfile_path(field_name, scan_num, outdir):
    """Return the expected path for a flux spectrum text file.

    Uses the same naming convention as save_flux_textfile so the two
    functions always agree on the filename.
    """
    safe_name = field_name.replace('/', '_').replace(' ', '_')
    safe_ms   = myms.replace('/', '_').replace(' ', '_').rstrip('_')
    fname = 'xf_fluxspec_{}_{}_{}_scan{:04d}.txt'.format(
        safe_ms, safe_name, 'IQUV', int(scan_num))
    return outdir.rstrip('/') + '/' + fname


def load_flux_textfile(fpath, nchan_xf):
    """Load cross-hand flux and fraction from an existing flux spectrum text file.

    Reads the four-column text file written by save_flux_textfile and
    returns the cross-hand flux and fraction columns binned to the Xf
    table channel resolution.  Uses equal weights since original
    per-channel weights are not stored; NaN entries are excluded.

    Parameters
    ----------
    fpath : str
        Full path to the text file.
    nchan_xf : int
        Number of channels in the Xf calibration table.

    Returns
    -------
    cross_flux_xf : ndarray, shape (nchan_xf,)
        Cross-hand flux binned to the Xf table resolution.
    cross_frac_xf : ndarray, shape (nchan_xf,)
        Cross-hand fraction (flux / I) binned to the Xf table resolution.
    """
    xflux = []
    xfrac = []
    with open(fpath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                xflux.append(float(parts[2]))
                xfrac.append(float(parts[3]) if len(parts) >= 4 else float('nan'))
            except ValueError:
                continue

    cross_flux_native = np.array(xflux, dtype=float)
    cross_frac_native = np.array(xfrac, dtype=float)
    weights_native    = np.where(np.isfinite(cross_flux_native), 1.0, 0.0)

    cross_flux_xf, _ = bin_channels_weighted(cross_flux_native, weights_native, nchan_xf)
    cross_frac_xf, _ = bin_channels_weighted(cross_frac_native, weights_native, nchan_xf)

    return cross_flux_xf, cross_frac_xf, cross_flux_native


def save_flux_textfile(field_name, ra_rad, dec_rad, scan_num,
                       t_start_mjd, t_end_mjd,
                       freq_ghz, I_chan, cross_flux_chan,
                       outdir):
    """Save per-channel total and cross-hand flux spectrum to a text file.

    One file is written per field-scan combination to `outdir`.  The
    header records the field name, J2000 direction, and scan time range
    in both MJD seconds and ISO 8601 (ISOT) format.

    Parameters
    ----------
    field_name : str
    ra_rad, dec_rad : float
        Field J2000 direction in radians (from FIELD subtable).
    scan_num : int
    t_start_mjd, t_end_mjd : float
        Scan start/end in MJD seconds.
    freq_ghz : ndarray
        Channel frequencies in GHz.
    I_chan : ndarray
        Time-averaged Stokes I per channel (Jy).
    cross_flux_chan : ndarray
        Time-averaged sqrt(U^2 + V^2) per channel (Jy).
    outdir : str
        Directory in which to write the file.
    """
    fpath = flux_textfile_path(field_name, scan_num, outdir)

    ra_str  = rad_to_hms(ra_rad)
    dec_str = rad_to_dms(dec_rad)

    t_start_isot = mjd_sec_to_isot(t_start_mjd)
    t_end_isot   = mjd_sec_to_isot(t_end_mjd)

    # Cross-hand fraction = sqrt(U^2+V^2) / I per channel
    with np.errstate(invalid='ignore', divide='ignore'):
        cross_frac_chan = np.where(
            np.isfinite(I_chan) & (I_chan > 0.0),
            cross_flux_chan / I_chan,
            np.nan)

    with open(fpath, 'w') as fh:
        fh.write('# Xf calibration flux spectrum\n')
        fh.write('# MS:                {}\n'.format(myms))
        fh.write('# Field:              {}\n'.format(field_name))
        fh.write('# Direction (J2000):  {} {}\n'.format(ra_str, dec_str))
        fh.write('# Scan:               {}\n'.format(scan_num))
        fh.write('# Time start (MJD s): {:.3f}\n'.format(t_start_mjd))
        fh.write('# Time start (ISOT):  {}\n'.format(t_start_isot))
        fh.write('# Time end   (MJD s): {:.3f}\n'.format(t_end_mjd))
        fh.write('# Time end   (ISOT):  {}\n'.format(t_end_isot))
        fh.write('# XF_MIN_CROSS_FRAC threshold: {:.6f} ({:.4f}%)\n'.format(
            XF_MIN_CROSS_FRAC, XF_MIN_CROSS_FRAC * 100.0))
        fh.write('#\n')
        fh.write('# {:>14s}  {:>14s}  {:>20s}  {:>20s}\n'.format(
            'freq_GHz', 'I_Jy', 'crosshand_flux_Jy', 'crosshand_frac'))

        for ch in range(len(freq_ghz)):
            i_val = I_chan[ch]          if np.isfinite(I_chan[ch])          else float('nan')
            x_val = cross_flux_chan[ch] if np.isfinite(cross_flux_chan[ch]) else float('nan')
            f_val = cross_frac_chan[ch] if np.isfinite(cross_frac_chan[ch]) else float('nan')
            fh.write('  {:>16.8f}  {:>14.6f}  {:>20.6f}  {:>20.8f}\n'.format(
                freq_ghz[ch], i_val, x_val, f_val))

    msg(f'    Saved flux spectrum: {fpath}')


# -----------------------------------------------------------------------
#  Processing functions
#
#  These are called once per Xf table (combineScan and, if applicable,
#  perScan).  They rely on module-level variables (myms, tb, msmd,
#  all_field_names, field_info, XF_* thresholds) that are defined
#  before any call is made.
# -----------------------------------------------------------------------

def _build_field_mean_frac(scan_xflux_map):
    """Mean cross-hand fraction across scans per field (for combineScan flagging)."""
    field_fracs = {}
    for (fn, sn), v in scan_xflux_map.items():
        field_fracs.setdefault(fn, []).append(v['frac'])
    result = {}
    for fn, fracs in field_fracs.items():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            result[fn] = np.nanmean(np.vstack(fracs), axis=0)
    return result


def _resolve_degeneracy(tab_path):
    """Flip each unflagged Xf channel by +/-180 deg to minimise distance to
    XF_REFERENCE_PHASE.  Modifies CPARAM in-place."""
    ref_rad = np.deg2rad(XF_REFERENCE_PHASE)
    tb.open(tab_path, nomodify=False)
    cparam = tb.getcol('CPARAM')
    flags  = tb.getcol('FLAG')
    n_flipped = 0
    for corr in range(cparam.shape[0]):
        phases    = np.angle(cparam[corr])
        diff      = np.abs(np.angle(np.exp(1j * (phases - ref_rad))))
        diff_flip = np.abs(np.angle(np.exp(1j * (phases + np.pi - ref_rad))))
        flip_mask = (diff_flip < diff) & (~flags[corr])
        cparam[corr][flip_mask] *= -1
        n_flipped += int(np.sum(flip_mask))
    tb.putcol('CPARAM', cparam)
    tb.close()
    msg(f'  Phase degeneracy resolved: reference = {XF_REFERENCE_PHASE:.1f} deg, '
        f'{n_flipped} channel(s) flipped.')


def _flag_flux_threshold(tab_path, field_frac_map):
    """Flag Xf table channels where cross-hand fraction < XF_MIN_CROSS_FRAC.

    Always operates per-field: all rows (scans × antennas) belonging to the
    same field receive the same flag mask, derived from the field-mean
    cross-hand fraction across scans.

    Parameters
    ----------
    tab_path : str
    field_frac_map : dict
        {field_name: frac_array} — per-field mean cross-hand fraction.
    """
    tb.open(tab_path, nomodify=False)
    cparam        = tb.getcol('CPARAM')
    tab_flag      = tb.getcol('FLAG')
    tab_field_ids = tb.getcol('FIELD_ID')

    n_rows_tab    = cparam.shape[2]
    n_corr_tab    = cparam.shape[0]
    n_chan_tab     = cparam.shape[1]
    n_total_chrows = n_rows_tab * n_corr_tab * n_chan_tab

    n_flagged_before = int(np.sum(tab_flag))
    msg(f'  Before flagging: {n_flagged_before} / {n_total_chrows} channel-rows already flagged')

    fs_stats = {}

    for row_idx in range(n_rows_tab):
        fid = int(tab_field_ids[row_idx])
        fn  = all_field_names[fid] if fid < len(all_field_names) else None

        if fn is None or fn not in field_info or fn not in field_frac_map:
            continue

        cross_frac_xf  = field_frac_map[fn]
        low_flux_mask  = np.isfinite(cross_frac_xf) & (cross_frac_xf < XF_MIN_CROSS_FRAC)
        low_flux_mask |= ~np.isfinite(cross_frac_xf)

        if fn not in fs_stats:
            n_before = int(np.sum(tab_flag[0, :, row_idx]))
            fs_stats[fn] = {'n_before': n_before, 'n_after': n_before}

        for corr in range(n_corr_tab):
            tab_flag[corr, low_flux_mask, row_idx] = True

        fs_stats[fn]['n_after'] = int(np.sum(tab_flag[0, :, row_idx]))

    for fn in sorted(fs_stats.keys()):
        st = fs_stats[fn]
        msg(f'  Field "{fn}":  {st["n_before"]} / {n_chan_tab} flagged  →  '
            f'{st["n_after"]} / {n_chan_tab} flagged  '
            f'(+{st["n_after"] - st["n_before"]})')

    n_flagged_after = int(np.sum(tab_flag))
    msg(f'  Total newly flagged: +{n_flagged_after - n_flagged_before} channel-rows  '
        f'(cross-hand fraction < {XF_MIN_CROSS_FRAC*100:.4f}% of Stokes I)')

    tb.putcol('FLAG', tab_flag)
    tb.close()


def _flag_mad_clip(tab_path, verbose=True, sigma=None, single_pass=False):
    """Sliding-window circular MAD clip on Xf phases.

    All operations are performed in the complex (phasor) domain so the
    result is completely insensitive to 2π phase wrapping.

    By default (single_pass=False) rows are grouped by field; all rows of a
    field are pooled for the polynomial fit and circular-mean MAD signal, and
    the resulting per-channel flag mask is broadcast to every row of that field.

    When single_pass=True all rows are treated as a single group (used for
    combineScanField, which carries the same cross-field average in every row).

    Returns
    -------
    poly_coeffs_by_field : dict
        {group_label: (coeffs_re, coeffs_im)} for groups where a polynomial
        was successfully fitted.  Empty dict when XF_USE_POLY is False.
    """
    _sigma   = sigma if sigma is not None else XF_MAD_SIGMA
    mode_tag = 'single-pass' if single_pass else 'per-field'
    mode_str = f'poly-residual [{mode_tag}]' if XF_USE_POLY else f'direct [{mode_tag}]'
    msg()
    msg(f'=== Flagging Xf table: circular MAD clip '
        f'(window={XF_MAD_WINDOW} ch, >{_sigma} sigma, mode={mode_str}) ===')

    tb.open(tab_path, nomodify=False)
    cparam        = tb.getcol('CPARAM')
    tab_flag      = tb.getcol('FLAG')
    tab_field_ids = tb.getcol('FIELD_ID')

    n_rows_tab    = cparam.shape[2]
    n_corr_tab    = cparam.shape[0]
    n_chan_tab     = cparam.shape[1]
    n_total_chrows = n_rows_tab * n_corr_tab * n_chan_tab

    chan_idx  = np.arange(n_chan_tab, dtype=float)
    chan_norm = 2.0 * chan_idx / max(n_chan_tab - 1, 1) - 1.0

    n_flagged_before = int(np.sum(tab_flag))
    msg(f'  Before flagging: {n_flagged_before} / {n_total_chrows} '
        f'channel-rows already flagged')

    # Group rows: all into one group (single_pass) or by field
    if single_pass:
        field_rows = {'all': list(range(n_rows_tab))}
    else:
        field_rows = {}
        for row_idx in range(n_rows_tab):
            fid = int(tab_field_ids[row_idx])
            fn  = all_field_names[fid] if fid < len(all_field_names) else 'unknown'
            if fn not in field_rows:
                field_rows[fn] = []
            field_rows[fn].append(row_idx)

    n_newly_total        = 0
    poly_coeffs_by_field = {}

    for fn in sorted(field_rows.keys()):
        rows_in_field = field_rows[fn]

        # Per-field polynomial fit from all rows/corrs of this field
        cr, ci = None, None
        if XF_USE_POLY:
            pool_c, pool_re, pool_im = [], [], []
            for row_idx in rows_in_field:
                for corr in range(n_corr_tab):
                    good = ~tab_flag[corr, :, row_idx]
                    if int(np.sum(good)) < 2:
                        continue
                    ph = np.angle(cparam[corr, :, row_idx])
                    pool_c.append(chan_norm[good])
                    pool_re.append(np.cos(ph[good]))
                    pool_im.append(np.sin(ph[good]))
            if pool_c:
                all_c  = np.concatenate(pool_c)
                all_re = np.concatenate(pool_re)
                all_im = np.concatenate(pool_im)
                if len(all_c) > XF_POLY_ORDER + 1:
                    cr = np.polyfit(all_c, all_re, XF_POLY_ORDER)
                    ci = np.polyfit(all_c, all_im, XF_POLY_ORDER)
                    poly_coeffs_by_field[fn] = (cr, ci)
                    msg(f'  [{fn}] Fitted degree-{XF_POLY_ORDER} polynomial '
                        f'({len(all_c)} samples, {len(rows_in_field)} rows pooled).')
                else:
                    msg(f'  [{fn}] WARNING: insufficient phasors for polynomial fit; '
                        'falling back to direct circular MAD.')

        # Circular mean phasor (corr=0) across all rows, baseline-subtracted
        poly_re_v = np.polyval(cr, chan_norm) if cr is not None else None
        poly_im_v = np.polyval(ci, chan_norm) if ci is not None else None

        sum_re = np.zeros(n_chan_tab)
        sum_im = np.zeros(n_chan_tab)
        counts = np.zeros(n_chan_tab, dtype=int)
        for row_idx in rows_in_field:
            flags_r = tab_flag[0, :, row_idx]
            good    = ~flags_r
            ph      = np.angle(cparam[0, :, row_idx])
            phasors = np.exp(1j * ph)
            if poly_re_v is not None:
                res = phasors * np.exp(-1j * np.arctan2(poly_im_v, poly_re_v))
            else:
                res = phasors
            sum_re[good] += np.real(res[good])
            sum_im[good] += np.imag(res[good])
            counts[good] += 1

        has_data = counts > 0
        safe_cnt = np.where(has_data, counts.astype(float), 1.0)
        signal   = np.where(has_data,
                            sum_re / safe_cnt + 1j * sum_im / safe_cnt,
                            0.0 + 0.0j)
        good_idx = np.where(has_data)[0]

        chan_dev_sigma = np.full(n_chan_tab, np.nan)
        chan_med_deg   = np.full(n_chan_tab, np.nan)
        chan_stdev_deg = np.full(n_chan_tab, np.nan)
        chan_off_deg   = np.full(n_chan_tab, np.nan)
        chan_win_n     = np.zeros(n_chan_tab, dtype=int)
        new_flag       = np.zeros(n_chan_tab, dtype=bool)

        if len(good_idx) >= 2:
            for ch in good_idx:
                others = good_idx[good_idx != ch]
                if len(others) < 2:
                    continue
                dists      = np.abs(others - ch)
                win_idx    = others[np.argsort(dists)[:XF_MAD_WINDOW]]
                win_sig    = signal[win_idx]
                chan_win_n[ch] = len(win_sig)

                win_angles  = np.angle(win_sig)
                best_sum    = np.inf
                circ_med_ph = win_angles[0]
                for _cand in win_angles:
                    _s = np.sum(np.abs(np.angle(np.exp(1j * (win_angles - _cand)))))
                    if _s < best_sum:
                        best_sum    = _s
                        circ_med_ph = _cand
                circ_med = np.exp(1j * circ_med_ph)

                circ_mad   = np.median(np.abs(np.angle(win_sig * np.conj(circ_med))))
                if circ_mad < 1e-10:
                    continue
                circ_stdev = 1.4826 * circ_mad

                off_ch    = np.angle(signal[ch] * np.conj(circ_med))
                dev_sigma = abs(off_ch) / circ_stdev

                chan_dev_sigma[ch] = dev_sigma
                chan_med_deg[ch]   = np.rad2deg(circ_med_ph)
                chan_stdev_deg[ch] = np.rad2deg(circ_stdev)
                chan_off_deg[ch]   = np.rad2deg(off_ch)

                if dev_sigma > _sigma:
                    new_flag[ch] = True

        n_newly = int(np.sum(new_flag))
        n_newly_total += n_newly

        # Broadcast flag mask to every row and corr in this field
        for row_idx in rows_in_field:
            for corr in range(n_corr_tab):
                tab_flag[corr, :, row_idx] |= new_flag

        if verbose:
            flags_row = ~has_data
            msg(f'  [{fn}  all scans pooled]  threshold={_sigma:.2f}σ')
            msg(f'  {"ch":>5}  {"N":>4}  {"median(deg)":>11}  '
                f'{"stdev(deg)":>10}  {"offset(deg)":>11}  '
                f'{"dev(σ)":>7}  {"thr(σ)":>7}  {"flag":>6}')
            for _c in range(n_chan_tab):
                if flags_row[_c]:
                    msg(f'  {_c:>5}  {"---":>4}  {"pre-flagged":>11}  '
                        f'{"---":>10}  {"---":>11}  {"---":>7}  '
                        f'{_sigma:>7.2f}  {"(pre)":>6}')
                elif np.isnan(chan_dev_sigma[_c]):
                    msg(f'  {_c:>5}  {chan_win_n[_c]:>4}  {"(no window)":>11}  '
                        f'{"---":>10}  {"---":>11}  {"---":>7}  '
                        f'{_sigma:>7.2f}  {"---":>6}')
                else:
                    flag_str = 'CLIP' if new_flag[_c] else 'ok'
                    msg(f'  {_c:>5}  {chan_win_n[_c]:>4}  '
                        f'{chan_med_deg[_c]:>+11.3f}  '
                        f'{chan_stdev_deg[_c]:>10.3f}  '
                        f'{chan_off_deg[_c]:>+11.3f}  '
                        f'{chan_dev_sigma[_c]:>7.3f}  '
                        f'{_sigma:>7.2f}  '
                        f'{flag_str:>6}')
        msg(f'  [{fn}  all scans pooled]  '
            f'newly flagged: {n_newly} / {n_chan_tab} channels')

    n_flagged_after = int(np.sum(tab_flag))
    msg(f'  Total newly flagged: +{n_newly_total} channels  '
        f'({n_flagged_after - n_flagged_before} channel-rows)')

    tb.putcol('FLAG', tab_flag)
    tb.close()

    return poly_coeffs_by_field


def _make_diagnostic_plot(tab_path, flag_polcal_snap, freq_ghz_xf,
                           poly_coeffs_by_field, frac_lookup,
                           plot_path, combine=False,
                           combine_xf_tab=None,
                           native_qu_lookup=None, freq_ghz_native=None,
                           phase_ylim=None, sn_ylim=None,
                           field_combine=False):
    """Diagnostic PDF: three-panel plot.

    Top panel    — Xf phases after flagging, coloured by (field, scan);
                   our-flagged channels shown as grey x markers;
                   per-field polynomial baselines overlaid when available.
    Middle panel — per-channel fractional cross-hand flux (sqrt(U²+V²)/I)
                   used for the flux-threshold flagging step, with the
                   XF_MIN_CROSS_FRAC threshold marked as a dashed line.
    Bottom panel — S/N of native-resolution sqrt(Q²+U²), where the global
                   noise is estimated from the median windowed MAD across
                   the band (log y-scale).
    """
    msg()
    msg(f'=== Producing diagnostic phase plot ===')
    msg(f'  Output: {plot_path}')

    tb.open(tab_path)
    cparam_plot = tb.getcol('CPARAM')
    flag_final  = tb.getcol('FLAG')
    fid_plot    = tb.getcol('FIELD_ID')
    scan_col_p  = tb.getcol('SCAN_NUMBER') if 'SCAN_NUMBER' in tb.colnames() else None
    tb.close()

    n_rows_p = cparam_plot.shape[2]
    n_corr_p = cparam_plot.shape[0]
    freq_x   = freq_ghz_xf
    chan_x   = np.arange(len(freq_x), dtype=float)
    chan_n   = 2.0 * chan_x / max(len(freq_x) - 1, 1) - 1.0

    our_flags = flag_final & ~flag_polcal_snap

    fs_keys_ordered = []
    fs_keys_set     = set()
    key_rep_row     = {}   # key -> first row index (corr=0 representative)
    # row_phase_key maps each row to its phase-panel key:
    # field_combine=True collapses all rows to a single 'combineScanField' key
    # for the phase panel only; fraction/S/N panels use per-field keys as normal.
    row_phase_key = {}
    _FIELD_COMBINE_KEY = ('combineScanField', -1)
    for row_idx in range(n_rows_p):
        fid = int(fid_plot[row_idx])
        fn  = all_field_names[fid] if fid < len(all_field_names) else 'unknown'
        sn  = int(scan_col_p[row_idx]) if scan_col_p is not None else -1
        key = _FIELD_COMBINE_KEY if field_combine else (fn, sn)
        row_phase_key[row_idx] = key
        if key not in fs_keys_set:
            fs_keys_ordered.append(key)
            fs_keys_set.add(key)
            key_rep_row[key] = row_idx

    n_keys     = len(fs_keys_ordered)
    cmap       = plt.get_cmap('tab20' if n_keys <= 20 else 'hsv')
    key_colour = {k: cmap(i / max(n_keys - 1, 1)) for i, k in enumerate(fs_keys_ordered)}

    # For fraction/S/N panels always use per-field keys (not collapsed)
    fs_perfield_keys_ordered = []
    fs_perfield_keys_set     = set()
    perfield_key_rep_row     = {}
    for row_idx in range(n_rows_p):
        fid = int(fid_plot[row_idx])
        fn  = all_field_names[fid] if fid < len(all_field_names) else 'unknown'
        sn  = int(scan_col_p[row_idx]) if scan_col_p is not None else -1
        pfkey = (fn, sn)
        if pfkey not in fs_perfield_keys_set:
            fs_perfield_keys_ordered.append(pfkey)
            fs_perfield_keys_set.add(pfkey)
            perfield_key_rep_row[pfkey] = row_idx

    fs_good_traces  = {k: [] for k in fs_keys_ordered}
    fs_ourflag_freq = {k: [] for k in fs_keys_ordered}
    fs_ourflag_ph   = {k: [] for k in fs_keys_ordered}

    for row_idx in range(n_rows_p):
        key = row_phase_key[row_idx]
        for corr in range(n_corr_p):
            phases     = np.rad2deg(np.angle(cparam_plot[corr, :, row_idx])).astype(float)
            pre_flag   = flag_polcal_snap[corr, :, row_idx]
            our_flag   = our_flags[corr, :, row_idx]
            good_trace = phases.copy()
            good_trace[pre_flag | our_flag] = np.nan
            fs_good_traces[key].append(good_trace)
            ours_only = our_flag & ~pre_flag
            fs_ourflag_freq[key].extend(freq_x[ours_only].tolist())
            fs_ourflag_ph[key].extend(phases[ours_only].tolist())

    # poly_coeffs_by_field: {field_name: (coeffs_re, coeffs_im)}
    # One dashed overlay per field is drawn after the scatter plots.


    BORDER_LW = 2.0
    TICK_MAJ  = 8
    TICK_MIN  = 4
    TICK_LW   = 1.5

    # Determine whether we can build the S/N panel
    _have_sn_data = (native_qu_lookup is not None) and (freq_ghz_native is not None)

    if _have_sn_data:
        fig, (ax, ax_frac, ax_sn) = plt.subplots(
            3, 1, figsize=(12, 11),
            gridspec_kw={'height_ratios': [3, 1.5, 1.5], 'hspace': 0.25})
    else:
        fig, (ax, ax_frac) = plt.subplots(
            2, 1, figsize=(12, 9),
            gridspec_kw={'height_ratios': [3, 1.5], 'hspace': 0.25})
        ax_sn = None

    # ---- Top panel: Xf phases ----
    plotted_line         = set()
    all_flagged_freq     = []
    all_flagged_ph       = []
    all_unflagged_phases = []

    for key in fs_keys_ordered:
        fn, sn = key
        colour = key_colour[key]
        if field_combine:
            label = 'combineScanField Xf'
        elif combine:
            label = f'{fn}: combineScan'
        else:
            label = f'{fn}: scan {sn}'
        for trace in fs_good_traces[key]:
            valid = np.isfinite(trace)
            if not np.any(valid):
                continue
            all_unflagged_phases.extend(trace[valid].tolist())
            lbl = label if key not in plotted_line else '_nolegend_'
            ax.scatter(freq_x[valid], trace[valid],
                       s=18, color=colour, alpha=0.7,
                       edgecolors='black', linewidths=0.4,
                       rasterized=True, label=lbl)
            plotted_line.add(key)
        all_flagged_freq.extend(fs_ourflag_freq[key])
        all_flagged_ph.extend(fs_ourflag_ph[key])

    poly_drawn = set()
    for _pkey in fs_keys_ordered:
        _pfn, _ = _pkey
        # field_combine: poly is stored under 'all'; label as combineScanField
        _poly_key = 'all' if field_combine else _pfn
        _poly_lbl = (f'combineScanField poly-{XF_POLY_ORDER} baseline' if field_combine
                     else f'{_pfn}: poly-{XF_POLY_ORDER} baseline')
        if _poly_key in poly_drawn or _poly_key not in poly_coeffs_by_field:
            continue
        _cr, _ci = poly_coeffs_by_field[_poly_key]
        _pre = np.polyval(_cr, chan_n)
        _pim = np.polyval(_ci, chan_n)
        _pph = np.rad2deg(np.arctan2(_pim, _pre))
        ax.plot(freq_x, _pph, color=key_colour[_pkey], lw=1.5, ls='--',
                zorder=10, label=_poly_lbl)
        poly_drawn.add(_poly_key)

    if all_flagged_freq:
        ax.scatter(all_flagged_freq, all_flagged_ph,
                   marker='x', s=18, linewidths=0.8,
                   color='#999999', alpha=0.6, zorder=3, label='Flagged')

    if all_unflagged_phases:
        y_lo  = min(all_unflagged_phases)
        y_hi  = max(all_unflagged_phases)
        y_pad = 0.5 * (y_hi - y_lo) if y_hi > y_lo else 10.0
        ax.set_ylim(y_lo - y_pad, y_hi + y_pad)
    if phase_ylim is not None:
        ax.set_ylim(phase_ylim)

    # ---- combineScan overlay on perScan plot ----
    # Computed directly from the per-scan data already in memory: phasor-average
    # cos(ph) and sin(ph) across all rows per channel, using flag_final to exclude
    # flagged channels. Only channels where at least one row is unflagged are shown.
    if combine_xf_tab is not None and not combine:
        n_chan_ov = cparam_plot.shape[1]
        avg_re_ov = np.zeros(n_chan_ov)
        avg_im_ov = np.zeros(n_chan_ov)
        counts_ov = np.zeros(n_chan_ov, dtype=int)
        for _ri in range(cparam_plot.shape[2]):
            _good_ch = ~flag_final[0, :, _ri]
            _ph      = np.angle(cparam_plot[0, :, _ri])
            avg_re_ov[_good_ch] += np.cos(_ph[_good_ch])
            avg_im_ov[_good_ch] += np.sin(_ph[_good_ch])
            counts_ov[_good_ch] += 1
        has_ov   = counts_ov > 0
        ov_phase = np.where(has_ov,
                            np.rad2deg(np.arctan2(avg_im_ov, avg_re_ov)),
                            np.nan)
        if np.any(has_ov):
            ax.scatter(freq_x[has_ov], ov_phase[has_ov],
                       s=35, marker='s', color='#444444',
                       alpha=0.45, edgecolors='none',
                       zorder=25, rasterized=True, label='combineScan solution')

    ax.set_ylabel('Xf phase (deg)', fontsize=13, labelpad=8)
    ax.set_title(f'Xf cross-hand phases after flagging  —  {myms}',
                 fontsize=13, pad=10)
    ax.set_xlim(freq_x[0], freq_x[-1])
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(3))
    for spine in ax.spines.values():
        spine.set_linewidth(BORDER_LW)
    ax.tick_params(which='major', direction='in',
                   length=TICK_MAJ, width=TICK_LW,
                   top=True, bottom=True, left=True, right=True,
                   labelsize=11)
    ax.tick_params(which='minor', direction='in',
                   length=TICK_MIN, width=TICK_LW * 0.75,
                   top=True, bottom=True, left=True, right=True)
    legend = ax.legend(
        loc='upper left', bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0, fontsize=9, frameon=True, framealpha=0.95,
        edgecolor='#666666', handlelength=2.2, handleheight=1.2,
        borderpad=0.8, labelspacing=0.5)
    legend.get_frame().set_linewidth(1.2)

    # ---- Middle panel: fractional cross-hand flux ----
    # Always per-field — use perfield keys (ignores phase-panel collapsing)
    _pf_n_keys   = len(fs_perfield_keys_ordered)
    _pf_cmap     = plt.get_cmap('tab20' if _pf_n_keys <= 20 else 'hsv')
    _pf_colour   = {k: _pf_cmap(i / max(_pf_n_keys - 1, 1))
                    for i, k in enumerate(fs_perfield_keys_ordered)}
    frac_plotted     = set()
    frac_flagged_freq = []
    frac_flagged_val  = []
    for key in fs_perfield_keys_ordered:
        fn, sn   = key
        colour   = _pf_colour[key]
        label    = f'{fn}: combineScan' if combine else f'{fn}: scan {sn}'
        if combine or field_combine:
            frac_arr = frac_lookup.get(fn, None)
        else:
            entry    = frac_lookup.get((fn, sn), None)
            frac_arr = entry['frac'] if isinstance(entry, dict) else entry
        if frac_arr is None:
            continue
        frac_arr = np.asarray(frac_arr, dtype=float)
        if len(frac_arr) != len(freq_x):
            continue

        # Per-channel flag masks (corr=0, representative row)
        rep = perfield_key_rep_row.get(key, 0)
        our_f = our_flags[0, :, rep]
        pre_f = flag_polcal_snap[0, :, rep]
        good  = ~our_f & ~pre_f
        bad   = our_f & ~pre_f   # we flagged it

        lbl = label if key not in frac_plotted else '_nolegend_'
        if np.any(good):
            ax_frac.scatter(freq_x[good], frac_arr[good] * 100.0,
                            s=18, color=colour, alpha=0.75,
                            edgecolors='black', linewidths=0.4,
                            rasterized=True, label=lbl)
            frac_plotted.add(key)
        if np.any(bad):
            frac_flagged_freq.extend(freq_x[bad].tolist())
            frac_flagged_val.extend((frac_arr[bad] * 100.0).tolist())

    if frac_flagged_freq:
        ax_frac.scatter(frac_flagged_freq, frac_flagged_val,
                        marker='x', s=18, linewidths=0.8,
                        color='#999999', alpha=0.6, zorder=3, label='Flagged')

    ax_frac.axhline(XF_MIN_CROSS_FRAC * 100.0,
                    color='red', lw=1.5, ls='--',
                    label=f'Threshold ({XF_MIN_CROSS_FRAC*100:.2f}%)')
    if ax_sn is None:
        ax_frac.set_xlabel('Frequency (GHz)', fontsize=12, labelpad=8)
    ax_frac.set_ylabel('Cross-hand fraction (%)', fontsize=12, labelpad=8)
    ax_frac.set_title('Fractional cross-hand flux  sqrt(U²+V²)/I',
                      fontsize=12, pad=6)
    ax_frac.set_xlim(freq_x[0], freq_x[-1])
    ax_frac.set_ylim(bottom=0)
    ax_frac.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax_frac.yaxis.set_minor_locator(ticker.AutoMinorLocator(3))
    for spine in ax_frac.spines.values():
        spine.set_linewidth(BORDER_LW)
    ax_frac.tick_params(which='major', direction='in',
                        length=TICK_MAJ, width=TICK_LW,
                        top=True, bottom=True, left=True, right=True,
                        labelsize=10)
    ax_frac.tick_params(which='minor', direction='in',
                        length=TICK_MIN, width=TICK_LW * 0.75,
                        top=True, bottom=True, left=True, right=True)
    ax_frac.legend(
        loc='upper left', bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0, fontsize=9, frameon=True, framealpha=0.95,
        edgecolor='#666666', handlelength=2.2, borderpad=0.8,
        labelspacing=0.5).get_frame().set_linewidth(1.2)

    # Suppress x-tick labels on panels that aren't the bottom-most
    if ax_sn is not None:
        ax.tick_params(labelbottom=False)
        ax_frac.tick_params(labelbottom=False)
    else:
        ax.tick_params(labelbottom=False)

    # ---- Bottom panel: sqrt(Q²+U²) S/N binned to Xf grid ----
    if ax_sn is not None and _have_sn_data:
        freq_nat   = np.asarray(freq_ghz_native, dtype=float)
        n_native_p = len(freq_nat)
        win_size   = max(3, n_native_p // 100)
        nchan_xf_p = len(freq_x)

        # Build field-pooled noise using per-field keys (always per-field)
        field_noise_map = {}
        all_field_keys = {}
        for key in fs_perfield_keys_ordered:
            fn, _ = key
            all_field_keys.setdefault(fn, []).append(key)
        for fn, fkeys in all_field_keys.items():
            all_mads = []
            for key in fkeys:
                qu_k = (native_qu_lookup.get(fn, None) if (combine or field_combine)
                        else native_qu_lookup.get(key, None))
                if qu_k is None:
                    continue
                qu_k = np.asarray(qu_k, dtype=float)
                if len(qu_k) != n_native_p:
                    continue
                for start in range(0, n_native_p, win_size):
                    window = qu_k[start:start + win_size]
                    valid  = window[np.isfinite(window)]
                    if len(valid) < 3:
                        continue
                    all_mads.append(np.nanmedian(np.abs(valid - np.nanmedian(valid))))
            if all_mads:
                field_noise_map[fn] = np.nanmedian(all_mads) * 1.4826

        sn_plotted      = set()
        sn_flagged_freq = []
        sn_flagged_val  = []

        for key in fs_perfield_keys_ordered:
            fn, sn   = key
            colour   = _pf_colour[key]
            label    = f'{fn}: combineScan' if combine else f'{fn}: scan {sn}'

            qu_arr = (native_qu_lookup.get(fn, None) if (combine or field_combine)
                      else native_qu_lookup.get((fn, sn), None))
            if qu_arr is None:
                continue
            qu_arr = np.asarray(qu_arr, dtype=float)
            if len(qu_arr) != n_native_p:
                continue

            field_noise = field_noise_map.get(fn, 0.0)
            if field_noise <= 0.0:
                continue

            # Bin to Xf grid (same as frac)
            weights_nat = np.where(np.isfinite(qu_arr), 1.0, 0.0)
            qu_xf, _    = bin_channels_weighted(qu_arr, weights_nat, nchan_xf_p)
            sn_xf       = qu_xf / field_noise

            # Per-channel flag masks
            rep   = perfield_key_rep_row.get(key, 0)
            our_f = our_flags[0, :, rep]
            pre_f = flag_polcal_snap[0, :, rep]
            good  = ~our_f & ~pre_f
            bad   = our_f & ~pre_f

            good_sn = good & np.isfinite(sn_xf) & (sn_xf > 0.0)
            bad_sn  = bad  & np.isfinite(sn_xf)

            lbl = label if key not in sn_plotted else '_nolegend_'
            if np.any(good_sn):
                ax_sn.scatter(freq_x[good_sn], sn_xf[good_sn],
                              s=18, color=colour, alpha=0.75,
                              edgecolors='black', linewidths=0.4,
                              rasterized=True, label=lbl)
                sn_plotted.add(key)
            if np.any(bad_sn):
                sn_flagged_freq.extend(freq_x[bad_sn].tolist())
                sn_flagged_val.extend(sn_xf[bad_sn].tolist())

        if sn_flagged_freq:
            ax_sn.scatter(sn_flagged_freq, sn_flagged_val,
                          marker='x', s=18, linewidths=0.8,
                          color='#999999', alpha=0.6, zorder=3, label='Flagged')

        ax_sn.set_yscale('log')
        ax_sn.set_ylim(bottom=1.0)
        if sn_ylim is not None:
            ax_sn.set_ylim(sn_ylim)
        if XF_MIN_SN > 0.0:
            ax_sn.axhline(XF_MIN_SN, color='red', lw=1.5, ls='--',
                          label=f'S/N threshold ({XF_MIN_SN:.1f})')
        ax_sn.set_xlabel('Frequency (GHz)', fontsize=12, labelpad=8)
        ax_sn.set_ylabel('S/N  [sqrt(Q²+U²) / noise]', fontsize=12, labelpad=8)
        ax_sn.set_title('Cross-hand power S/N  (Xf grid)',
                        fontsize=12, pad=6)
        ax_sn.set_xlim(freq_x[0], freq_x[-1])
        ax_sn.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        for spine in ax_sn.spines.values():
            spine.set_linewidth(BORDER_LW)
        ax_sn.tick_params(which='major', direction='in',
                          length=TICK_MAJ, width=TICK_LW,
                          top=True, bottom=True, left=True, right=True,
                          labelsize=10)
        ax_sn.tick_params(which='minor', direction='in',
                          length=TICK_MIN, width=TICK_LW * 0.75,
                          top=True, bottom=True, left=True, right=True)
        ax_sn.legend(
            loc='upper left', bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0, fontsize=9, frameon=True, framealpha=0.95,
            edgecolor='#666666', handlelength=2.2, borderpad=0.8,
            labelspacing=0.5).get_frame().set_linewidth(1.2)
    elif ax_sn is not None:
        ax_sn.set_visible(False)

    plt.savefig(plot_path, bbox_inches='tight', dpi=200)
    phase_ylim_out = ax.get_ylim()
    sn_ylim_out    = ax_sn.get_ylim() if ax_sn is not None else None
    plt.close(fig)
    msg(f'  Saved: {plot_path}')
    return phase_ylim_out, sn_ylim_out


def _flag_sn_threshold(tab_path, scan_qu_map, freq_ghz_native):
    """Flag Xf table channels where the binned cross-hand S/N < XF_MIN_SN.

    Noise is estimated per field by pooling ALL per-window MADs from every
    window in every scan of that field, then taking the median × 1.4826.
    S/N is then computed per (field, scan): binned scan QU / field noise.
    Channels below XF_MIN_SN in any scan are flagged.  Only active when
    XF_MIN_SN > 0.
    """
    if XF_MIN_SN <= 0.0:
        return

    msg()
    msg(f'=== Flagging Xf table: cross-hand S/N < {XF_MIN_SN:.2f} ===')

    tb.open(tab_path, nomodify=False)
    tab_flag      = tb.getcol('FLAG')
    tab_field_ids = tb.getcol('FIELD_ID')
    tab_scan_s    = tb.getcol('SCAN_NUMBER') if 'SCAN_NUMBER' in tb.colnames() else None
    n_rows_tab    = tab_flag.shape[2]
    n_corr_tab    = tab_flag.shape[0]
    n_chan_tab     = tab_flag.shape[1]

    win_size = max(3, nchan_native // 100)

    # Step 1: per-field noise — pool every window MAD from every scan
    field_noise_map = {}
    field_scans = {}
    for (fn, sn), qu_arr in scan_qu_map.items():
        field_scans.setdefault(fn, []).append(sn)
    for fn, scans in field_scans.items():
        all_mads = []
        for sn in scans:
            qu_arr = np.asarray(scan_qu_map[(fn, sn)], dtype=float)
            if len(qu_arr) != nchan_native:
                continue
            for start in range(0, nchan_native, win_size):
                window = qu_arr[start:start + win_size]
                valid  = window[np.isfinite(window)]
                if len(valid) >= 3:
                    all_mads.append(np.nanmedian(np.abs(valid - np.nanmedian(valid))))
        if all_mads:
            noise = np.nanmedian(all_mads) * 1.4826
            if noise > 0.0:
                field_noise_map[fn] = noise
                msg(f'  Field "{fn}": noise = {noise:.4e}  '
                    f'(pooled {len(all_mads)} windows from {len(scans)} scan(s))')

    # Step 2: group table rows by (field, scan) when possible, else by field
    key_rows = {}
    for row_idx in range(n_rows_tab):
        fid = int(tab_field_ids[row_idx])
        fn  = all_field_names[fid] if fid < len(all_field_names) else None
        if fn is None or fn not in field_info:
            continue
        if tab_scan_s is not None:
            key = (fn, int(tab_scan_s[row_idx]))
        else:
            key = (fn, None)
        key_rows.setdefault(key, []).append(row_idx)

    # Step 3: flag per key
    n_new_total = 0
    for key, rows in key_rows.items():
        fn, sn = key
        if fn not in field_noise_map:
            continue
        noise = field_noise_map[fn]

        if sn is not None and (fn, sn) in scan_qu_map:
            qu_arr = np.asarray(scan_qu_map[(fn, sn)], dtype=float)
        else:
            # Fall back to field mean
            arrs = [np.asarray(scan_qu_map[(fn, s)], dtype=float)
                    for s in field_scans.get(fn, [])
                    if (fn, s) in scan_qu_map]
            if not arrs:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                qu_arr = np.nanmean(np.array(arrs), axis=0)

        if len(qu_arr) != nchan_native:
            continue

        weights_nat = np.where(np.isfinite(qu_arr), 1.0, 0.0)
        qu_xf, _    = bin_channels_weighted(qu_arr, weights_nat, n_chan_tab)
        sn_xf       = qu_xf / noise

        new_flag = np.isfinite(sn_xf) & (sn_xf < XF_MIN_SN)
        n_new    = int(np.sum(new_flag))
        if n_new == 0:
            continue

        for row_idx in rows:
            for corr in range(n_corr_tab):
                tab_flag[corr, :, row_idx] |= new_flag
        n_new_total += n_new
        lbl = f'Field "{fn}" scan {sn}' if sn is not None else f'Field "{fn}"'
        msg(f'  {lbl}: {n_new} / {n_chan_tab} channels below S/N={XF_MIN_SN:.2f}')

    msg(f'  Total newly flagged by S/N cut: {n_new_total} channels')
    tb.putcol('FLAG', tab_flag)
    tb.close()


def _avg_perscan_to_combine(perscan_tab, combine_tab):
    """Overwrite combineScan CPARAM with a cross-hand flux weighted phasor average
    of per-scan solutions, renormalised to unit amplitude per channel.

    Weight for each scan at each channel = binned cross-hand flux (sqrt(Q²+U²))
    for that (field, scan) from scan_native_qu_map.  Higher-S/N scans contribute
    more.  Channels flagged in ALL contributing scans remain flagged.
    """
    msg()
    msg('=== Building combineScan solution (cross-hand flux weighted per-scan average) ===')

    tb.open(perscan_tab)
    per_cparam = tb.getcol('CPARAM')
    per_flag   = tb.getcol('FLAG')
    per_fid    = tb.getcol('FIELD_ID')
    per_ant    = tb.getcol('ANTENNA1')
    per_scan   = tb.getcol('SCAN_NUMBER') if 'SCAN_NUMBER' in tb.colnames() else None
    tb.close()

    tb.open(combine_tab, nomodify=False)
    comb_cparam = tb.getcol('CPARAM')
    comb_flag   = tb.getcol('FLAG')
    comb_fid    = tb.getcol('FIELD_ID')
    comb_ant    = tb.getcol('ANTENNA1')
    n_corr      = comb_cparam.shape[0]
    n_chan       = comb_cparam.shape[1]
    n_comb_rows = comb_cparam.shape[2]

    # Pre-bin per-scan QU weights to the Xf channel grid
    scan_w_xf = {}
    for (fn, sn), qu_nat in scan_native_qu_map.items():
        qu_nat = np.asarray(qu_nat, dtype=float)
        w_nat  = np.where(np.isfinite(qu_nat), 1.0, 0.0)
        qu_xf, _ = bin_channels_weighted(qu_nat, w_nat, n_chan)
        w = np.where(np.isfinite(qu_xf) & (qu_xf > 0.0), qu_xf, 0.0)
        scan_w_xf[(fn, sn)] = w

    # Log per-field scan weights (median weight as proxy for relative S/N)
    for fn in sorted(set(k[0] for k in scan_w_xf)):
        parts = []
        for (f2, sn), w in sorted(scan_w_xf.items()):
            if f2 != fn:
                continue
            parts.append(f'scan {sn}: {np.nanmedian(w):.4f}')
        msg(f'  [{fn}] cross-hand flux weights — ' + ',  '.join(parts))

    for comb_row in range(n_comb_rows):
        fid = int(comb_fid[comb_row])
        ant = int(comb_ant[comb_row])
        fn  = all_field_names[fid] if fid < len(all_field_names) else None
        per_rows_f = [r for r in range(per_cparam.shape[2])
                      if int(per_fid[r]) == fid and int(per_ant[r]) == ant]
        if not per_rows_f:
            continue

        for corr in range(n_corr):
            wsum_re = np.zeros(n_chan)
            wsum_im = np.zeros(n_chan)
            wsum    = np.zeros(n_chan)

            for pr in per_rows_f:
                sn   = int(per_scan[pr]) if per_scan is not None else None
                w    = scan_w_xf.get((fn, sn), np.ones(n_chan)) if fn else np.ones(n_chan)
                good = ~per_flag[corr, :, pr]
                ph   = np.angle(per_cparam[corr, :, pr])
                eff_w    = np.where(good, w, 0.0)
                wsum_re += eff_w * np.cos(ph)
                wsum_im += eff_w * np.sin(ph)
                wsum    += eff_w

            has_data = wsum > 0.0
            amp      = np.sqrt(wsum_re**2 + wsum_im**2)
            valid    = has_data & (amp > 1e-10)
            safe_amp = np.where(valid, amp, 1.0)
            norm_re  = np.where(valid, wsum_re / safe_amp, 0.0)
            norm_im  = np.where(valid, wsum_im / safe_amp, 0.0)

            comb_cparam[corr, :, comb_row] = (norm_re + 1j * norm_im).astype(complex)
            comb_flag[corr, :, comb_row]   = ~has_data

    tb.putcol('CPARAM', comb_cparam)
    tb.putcol('FLAG',   comb_flag)
    tb.close()
    msg('  combineScan CPARAM replaced with cross-hand flux weighted per-scan average.')


def _avg_combine_fields_to_global(combine_tab, field_tab, field_qu_weights):
    """Build combineScanField: cross-hand flux weighted circular average across fields.

    For each antenna and channel, the phasor from each field row is weighted by
    the mean cross-hand flux (sqrt(Q²+U²)) at that channel for that field.
    This gives brighter channels more influence on the combined solution.
    Channels flagged in every contributing field row remain flagged.

    field_qu_weights: {field_name: 1-D array of length n_chan} — per-channel
        cross-hand flux at Xf resolution, one entry per PA cal field.
    """
    import shutil
    msg()
    msg('=== Building combineScanField solution (cross-hand flux weighted average) ===')

    if os.path.exists(field_tab):
        shutil.rmtree(field_tab)
    shutil.copytree(combine_tab, field_tab)
    msg(f'  Copied {combine_tab} → {field_tab}')

    tb.open(field_tab, nomodify=False)
    cparam   = tb.getcol('CPARAM')
    flag     = tb.getcol('FLAG')
    ant      = tb.getcol('ANTENNA1')
    fid_col  = tb.getcol('FIELD_ID')
    n_corr   = cparam.shape[0]
    n_chan    = cparam.shape[1]
    n_rows    = cparam.shape[2]
    tb.close()

    # Build field_id → per-channel weight array (fallback: uniform weight 1)
    fid_to_weight = {}
    for fn, fid_tuple in field_info.items():
        fid = fid_tuple[0]
        w   = field_qu_weights.get(fn)
        if w is not None:
            w = np.asarray(w, dtype=float)
            w = np.where(np.isfinite(w) & (w > 0.0), w, 0.0)
        else:
            w = np.ones(n_chan, dtype=float)
        fid_to_weight[fid] = w
        msg(f'  [{fn}] median cross-hand flux weight: {np.nanmedian(w):.4f}')

    unique_ants = sorted(set(int(ant[r]) for r in range(n_rows)))

    for a in unique_ants:
        ant_rows = [r for r in range(n_rows) if int(ant[r]) == a]
        for corr in range(n_corr):
            wsum_re  = np.zeros(n_chan)
            wsum_im  = np.zeros(n_chan)
            wsum     = np.zeros(n_chan)
            for r in ant_rows:
                fid = int(fid_col[r])
                w   = fid_to_weight.get(fid, np.ones(n_chan))
                good = ~flag[corr, :, r]
                ph   = np.angle(cparam[corr, :, r])
                eff_w            = np.where(good, w, 0.0)
                wsum_re         += eff_w * np.cos(ph)
                wsum_im         += eff_w * np.sin(ph)
                wsum            += eff_w
            has_data = wsum > 0.0
            amp      = np.sqrt(wsum_re**2 + wsum_im**2)
            valid    = has_data & (amp > 1e-10)
            safe_amp = np.where(valid, amp, 1.0)
            norm_re  = np.where(valid, wsum_re / safe_amp, 0.0)
            norm_im  = np.where(valid, wsum_im / safe_amp, 0.0)
            for r in ant_rows:
                cparam[corr, :, r] = (norm_re + 1j * norm_im).astype(complex)
                flag[corr, :, r]   = ~has_data

    tb.open(field_tab, nomodify=False)
    tb.putcol('CPARAM', cparam)
    tb.putcol('FLAG',   flag)
    tb.close()
    msg('  combineScanField CPARAM set to cross-hand flux weighted cross-field average.')


def _freq_range_to_spw(freq_str, freq_ghz_xf):
    """Convert a frequency range string to a CASA SPW:chan selection for a cal table.

    Accepts either an already-resolved 'SPW:chan' string (passed through as-is)
    or a frequency range like '0.92~0.94GHz' / '920~940MHz' / '9.2e8~9.4e8Hz'.
    Returns None if the range falls entirely outside the table's frequency axis.
    """
    if ':' in freq_str:
        return freq_str  # already channel-based, pass through

    _unit_scales = {'ghz': 1.0, 'mhz': 1e-3, 'khz': 1e-6, 'hz': 1e-9}
    _s = freq_str.strip()
    _scale = 1.0
    for _unit, _sc in _unit_scales.items():
        if _s.lower().endswith(_unit):
            _scale = _sc
            _s = _s[:-len(_unit)]
            break
    try:
        _lo_str, _hi_str = _s.split('~')
        _lo_ghz = float(_lo_str) * _scale
        _hi_ghz = float(_hi_str) * _scale
    except (ValueError, TypeError):
        msg(f'  WARNING: Could not parse XF_FLAG_FREQS entry "{freq_str}" — skipping.')
        return None

    _chans = np.where((freq_ghz_xf >= _lo_ghz) & (freq_ghz_xf <= _hi_ghz))[0]
    if len(_chans) == 0:
        msg(f'  WARNING: Frequency range "{freq_str}" ({_lo_ghz:.4f}~{_hi_ghz:.4f} GHz) '
            f'has no matching channels in this table — skipping.')
        return None
    return f'0:{_chans[0]}~{_chans[-1]}'


def _flag_freq_ranges(tab_path, freq_ghz_xf):
    """Apply manual SPW:channel flagging to an Xf table (XF_FLAG_FREQS).
    Set XF_FLAG_CHANS to override with direct channel indices ('0:10~25').
    """
    if XF_FLAG_CHANS:
        _entries = XF_FLAG_CHANS if isinstance(XF_FLAG_CHANS, list) else [XF_FLAG_CHANS]
        _spws = [_freq_range_to_spw(f, freq_ghz_xf) for f in _entries]
        _spws = [s for s in _spws if s is not None]
        if _spws:
            flagdata(vis=tab_path, mode='manual', spw=','.join(_spws), flagbackup=False)
            msg(f'  XF_FLAG_CHANS override applied: {_spws}')
        return
    if not XF_FLAG_FREQS:
        msg('  XF_FLAG_FREQS is empty — no table-level frequency flagging applied.')
        return

    spw_strs = [_freq_range_to_spw(f, freq_ghz_xf) for f in XF_FLAG_FREQS]
    spw_strs = [s for s in spw_strs if s is not None]
    if not spw_strs:
        msg('  No valid frequency ranges to flag after conversion — skipping.')
        return

    flagspw = ','.join(spw_strs)
    flagdata(vis=tab_path, mode='manual', spw=flagspw, flagbackup=False)
    msg(f'  Flagged SPW ranges in Xf table: {flagspw}  (from {XF_FLAG_FREQS})')


def _process_xftab(tab_path, label, frac_lookup, use_field_mean,
                   field_frac_map, scan_qu_map,
                   flag_polcal_snap, freq_ghz_xf, plot_path,
                   combine_xf_tab=None,
                   native_qu_lookup=None, freq_ghz_native=None,
                   skip_flagging=False,
                   phase_ylim=None, sn_ylim=None,
                   mad_only=False,
                   fine_sigma=None):
    """Full post-polcal processing pipeline for one Xf table:
      1. Resolve +/-180 deg phase degeneracy
      2. (skip if skip_flagging/mad_only) Gross MAD pre-clip (per-field)
      3. (skip if skip_flagging/mad_only) Flag channels below cross-hand flux threshold
      4. (skip if skip_flagging/mad_only) S/N threshold flagging (per-field)
      5. Fine MAD clip (per-field)  [skipped only if skip_flagging=True]
      6. Diagnostic plot
      7. Manual frequency flagging

    scan_qu_map: {(field_name, scan_num): qu_array} — per-scan native-res QU arrays
    field_frac_map: {field_name: frac_array}

    skip_flagging=True: only degeneracy + plot + freq flags (phasor-averaged tables)
    mad_only=True: degeneracy + fine MAD only (no gross/frac/S/N) + plot + freq flags
                   used for combineScanField (solved by polcal combine='scan,field')
    Returns (phase_ylim, sn_ylim) from the plot.
    """
    msg()
    msg(f'--- Degeneracy resolution [{label}] ---')
    _resolve_degeneracy(tab_path)

    if not skip_flagging and not mad_only:
        # 1. Gross outlier pass — catches wildly wrong channels before other cuts
        _flag_mad_clip(tab_path, verbose=False, sigma=XF_GROSS_CLIP)

        # 2. Fractional cross-hand flux threshold (per-field)
        msg()
        msg(f'=== Flagging Xf table [{label}]: cross-hand fraction < '
            f'{XF_MIN_CROSS_FRAC*100:.4f}% of Stokes I ===')
        _flag_flux_threshold(tab_path, field_frac_map)

        # 3. S/N threshold (per-field, per-scan)
        if XF_MIN_SN > 0.0:
            _flag_sn_threshold(tab_path, scan_qu_map, freq_ghz_native)

    if not skip_flagging:
        # Fine MAD clip — single pass for combineScanField (all rows carry the
        # same cross-field average), per-field for all other tables.
        _fine_sigma = fine_sigma if fine_sigma is not None else XF_MAD_SIGMA
        poly_coeffs = _flag_mad_clip(tab_path, verbose=XF_VERBOSE,
                                     sigma=_fine_sigma, single_pass=mad_only)
    else:
        msg(f'  [skip_flagging=True — solution is phasor average; flagging skipped]')
        poly_coeffs = {}

    phase_ylim_out, sn_ylim_out = _make_diagnostic_plot(
        tab_path, flag_polcal_snap, freq_ghz_xf,
        poly_coeffs, frac_lookup,
        plot_path, combine=use_field_mean,
        combine_xf_tab=combine_xf_tab,
        native_qu_lookup=native_qu_lookup,
        freq_ghz_native=freq_ghz_native,
        phase_ylim=phase_ylim, sn_ylim=sn_ylim,
        field_combine=mad_only)

    _flag_freq_ranges(tab_path, freq_ghz_xf)
    return phase_ylim_out, sn_ylim_out


# -----------------------------------------------------------------------
#  Set model for PA calibrator
# -----------------------------------------------------------------------

_tmp_fields = [f.strip() for f in pacal_name.split(',') if f.strip()]

for _pf in _tmp_fields:
    setjy(vis=myms,
        field=_pf,
        standard='manual',
        fluxdensity=POLANG_MOD,
        usescratch=True)


# -----------------------------------------------------------------------
#  Determine if any pacal field has more than one scan
# -----------------------------------------------------------------------

_field_scan_counts = {}
msmd.open(myms)
try:
    _all_names_tmp = msmd.fieldnames()
    for _fn in _tmp_fields:
        if _fn in _all_names_tmp:
            _fid_tmp = _all_names_tmp.index(_fn)
            _field_scan_counts[_fn] = len(list(msmd.scansforfield(_fid_tmp)))
finally:
    msmd.close()

do_per_scan    = any(c > 1 for c in _field_scan_counts.values())
do_multi_field = len(_tmp_fields) > 1
msg(f'Scan counts per PA cal field: {_field_scan_counts}')
msg(f'do_per_scan    = {do_per_scan}')
msg(f'do_multi_field = {do_multi_field}')

xf_prefix      = GAINTABLES + '/cal_1GC_' + myms
xftab_combine  = xf_prefix + '_combineScan.Xf'
xftab_per      = xf_prefix + '_perScan.Xf'
xftab_field    = xf_prefix + '_combineScanField.Xf'

# Tables that will be produced (conditionally) in this run
_expected_tabs = [xftab_combine]
if do_per_scan:
    _expected_tabs.append(xftab_per)
if do_multi_field:
    _expected_tabs.append(xftab_field)

_tabs_exist = all(os.path.exists(t) for t in _expected_tabs)

msg()
msg('=== Xf table existence check ===')
for _t in _expected_tabs:
    msg(f'  {"[OK]" if os.path.exists(_t) else "[MISSING]"}  {_t}')

if _tabs_exist:
    msg()
    msg('  All Xf tables already exist — skipping polcal solves.')
    msg('  To force a complete reprocess, run:')
    msg()
    msg('    rm -rf ' + ' '.join(_expected_tabs))
    msg()
else:
    msg()
    msg('  One or more tables missing — running polcal solves.')


# -----------------------------------------------------------------------
#  Solve for Xf  (no KCROSS)
#
#  XF_AVG_PERSCAN=True (default):  solve perScan first, then build the
#    combineScan table by averaging per-scan phasors.  PerScan must be
#    available (do_per_scan must be True); if not, fall back to CASA
#    combine='scan' with a warning.
#  XF_AVG_PERSCAN=False: CASA polcal combine='scan' for combineScan,
#    then perScan if applicable (original behaviour).
#  combineScanField (do_multi_field): solved directly with combine='scan,field'
#    — CASA pools all scans of all PA cal fields into a single solution.
# -----------------------------------------------------------------------

_polcal_kwargs = dict(
    vis      = myms,
    field    = pacal_name,
    uvrange  = myuvrange,
    refant   = str(ref_ant),
    solint   = f'inf,{XF_CHANINT}ch',
    poltype  = 'Xf',
    gaintable = [ktab, bptab, gptab, gtab, dftab],
    gainfield = [pacal_name, bpcal_name, pacal_name, pacal_name, bpcal_name],
    interp    = ['nearest', 'linear', 'linear', 'linear', 'linear'],
    append    = False,
)

if not _tabs_exist:
    if do_per_scan:
        msg()
        msg('=== Solving Xf (perScan) ===')
        polcal(caltable=xftab_per, combine='', **_polcal_kwargs)

    # combineScan template: solved with combine='scan' to get the correct table
    # structure/dimensions.  CPARAM is overwritten after per-scan flagging by
    # _avg_perscan_to_combine (cross-hand flux weighted average of per-scan solutions).
    # If only one scan exists, the per-scan solution IS the combine solution.
    msg()
    msg('=== Solving Xf (combineScan template) ===')
    polcal(caltable=xftab_combine, combine='scan', **_polcal_kwargs)

    # combineScanField template is built from xftab_combine via copytree inside
    # _avg_combine_fields_to_global — no separate polcal solve needed.

# Read FLAG snapshots from whichever tables now exist
if do_per_scan:
    tb.open(xftab_per)
    flag_polcal_per = tb.getcol('FLAG').copy()
    tb.close()

tb.open(xftab_combine)
flag_polcal_combine = tb.getcol('FLAG').copy()
tb.close()


# -----------------------------------------------------------------------
#  Read spectral axis information (from combineScan table)
# -----------------------------------------------------------------------

msg()
msg('Reading spectral window information...')

tb.open(myms + '/SPECTRAL_WINDOW')
chan_freq_native = tb.getcol('CHAN_FREQ').flatten()   # Hz
tb.close()
freq_ghz_native = chan_freq_native / 1.0e9
nchan_native    = len(freq_ghz_native)
msg(f'  MS native channels   : {nchan_native}  '
      f'({freq_ghz_native.min():.4f} - {freq_ghz_native.max():.4f} GHz)')

tb.open(xftab_combine + '/SPECTRAL_WINDOW')
chan_freq_xf = tb.getcol('CHAN_FREQ').flatten()       # Hz
tb.close()
freq_ghz_xf = chan_freq_xf / 1.0e9
nchan_xf    = len(freq_ghz_xf)
msg(f'  Xf table channels    : {nchan_xf}  '
      f'({freq_ghz_xf.min():.4f} - {freq_ghz_xf.max():.4f} GHz)')

# Auto-set MAD window if not specified (nearest odd integer of nchan_xf / 10)
if XF_MAD_WINDOW <= 0:
    _auto_w = max(3, int(round(nchan_xf / 10.0)))
    if _auto_w % 2 == 0:
        _auto_w += 1
    XF_MAD_WINDOW = _auto_w
    msg(f'  XF_MAD_WINDOW auto-set to {XF_MAD_WINDOW} '
        f'(nearest odd to {nchan_xf}/10 = {nchan_xf/10:.1f})')


# -----------------------------------------------------------------------
#  Read FIELD subtable  (names + J2000 directions)
# -----------------------------------------------------------------------

tb.open(myms + '/FIELD')
all_field_names = list(tb.getcol('NAME'))
ref_dirs        = tb.getcol('REFERENCE_DIR')   # shape (2, 1, nfields) radians
tb.close()

pacal_fields = [f.strip() for f in pacal_name.split(',')]
field_info   = {}

for fn in pacal_fields:
    if fn in all_field_names:
        fid     = all_field_names.index(fn)
        ra_rad  = float(ref_dirs[0, 0, fid])
        dec_rad = float(ref_dirs[1, 0, fid])
        field_info[fn] = (fid, ra_rad, dec_rad)
        msg(f'  Field "{fn}": FIELD_ID={fid}  '
              f'RA={rad_to_hms(ra_rad)}  Dec={rad_to_dms(dec_rad)}')
    else:
        msg(f'  WARNING: Field "{fn}" not found in MS FIELD table; skipping.')

do_multi_field = len(field_info) > 1
msg(f'do_multi_field = {do_multi_field}')


# -----------------------------------------------------------------------
#  Check whether flux spectrum text files already exist for every
#  field-scan combination.  If they all do, skip the applycal and
#  visibility extraction and load the cross-hand flux directly from
#  the cached files — saving significant time on re-runs.
# -----------------------------------------------------------------------

msg()
msg('=== Checking for cached flux spectrum files ===')

all_field_scans = []
for fn in list(field_info.keys()):
    fid = field_info[fn][0]
    msmd.open(myms)
    try:
        scan_nums_fn = list(msmd.scansforfield(fid))
    except Exception as e:
        msmd.close()
        msg(f'  WARNING: Could not retrieve scans for field "{fn}": {e}')
        continue
    msmd.close()
    for sn in scan_nums_fn:
        all_field_scans.append((fn, sn))

missing_files = []
for fn, sn in all_field_scans:
    fpath = flux_textfile_path(fn, sn, RESULTS)
    if not os.path.exists(fpath):
        missing_files.append((fn, sn, fpath))

scan_xflux_map    = {}
scan_native_qu_map = {}

if not missing_files:
    msg(f'  All {len(all_field_scans)} flux spectrum file(s) found — '
        f'loading from cache, skipping applycal and extraction.')
    msg()
    for fn, sn in all_field_scans:
        fpath = flux_textfile_path(fn, sn, RESULTS)
        msg(f'  Loading: {fpath}')
        cross_flux_xf, cross_frac_xf, cross_flux_native_c = load_flux_textfile(fpath, nchan_xf)
        scan_xflux_map[(fn, sn)]     = {'flux': cross_flux_xf, 'frac': cross_frac_xf}
        scan_native_qu_map[(fn, sn)] = cross_flux_native_c   # sqrt(U²+V²) as proxy
    msg()
    msg(f'  Loaded cross-hand flux for {len(scan_xflux_map)} field-scan(s).')

else:
    msg(f'  {len(missing_files)} / {len(all_field_scans)} flux spectrum file(s) missing '
        f'— running applycal and extraction.')
    for fn, sn, fp in missing_files:
        msg(f'    Missing: {fp}')

    msg()
    msg('=== Applying calibration tables (K, B, Gp, F, Df) for flux extraction ===')

    applycal(vis=myms,
        field=pacal_name,
        parang=False,
        gaintable=[ktab, bptab, gptab, ftab, dftab],
        gainfield=[pacal_name, bpcal_name, pacal_name, pacal_name, bpcal_name],
        interp=['nearest', 'linear', 'linear', 'linear', 'linear'],
        flagbackup=False)

    msg('  applycal complete — reading from CORRECTED_DATA column.')

    msg()
    msg('=== Extracting per-field-scan cross-hand flux spectra ===')

    for field_name, (fid, ra_rad, dec_rad) in field_info.items():

        scan_nums = [sn for fn, sn in all_field_scans if fn == field_name]

        msg()
        msg(f'  Field "{field_name}" — {len(scan_nums)} scan(s): {scan_nums}')

        for scan_num in scan_nums:

            msg(f'    Scan {scan_num}: time-averaging visibilities...')

            I_chan, Q_chan, U_chan, V_chan, weight_uv, t_start, t_end = \
                time_average_iquv(field_name, scan_num)

            if I_chan is None:
                msg(f'    SKIP: no data returned for field={field_name} scan={scan_num}')
                continue

            cross_flux_native = np.sqrt(
                np.where(np.isfinite(U_chan), U_chan, 0.0)**2 +
                np.where(np.isfinite(V_chan), V_chan, 0.0)**2)
            uv_all_nan = (~np.isfinite(U_chan)) & (~np.isfinite(V_chan))
            cross_flux_native[uv_all_nan] = np.nan

            with np.errstate(invalid='ignore', divide='ignore'):
                cross_frac_native = np.where(
                    np.isfinite(I_chan) & (I_chan > 0.0),
                    cross_flux_native / I_chan,
                    np.nan)

            n_below = int(np.sum(
                np.isfinite(cross_frac_native) &
                (cross_frac_native < XF_MIN_CROSS_FRAC)))
            n_valid = int(np.sum(np.isfinite(cross_frac_native)))
            msg(f'    Valid channels: {n_valid}  |  '
                f'Below {XF_MIN_CROSS_FRAC*100:.4f}% threshold: {n_below} / {n_valid}')

            save_flux_textfile(
                field_name, ra_rad, dec_rad, scan_num,
                t_start, t_end,
                freq_ghz_native, I_chan, cross_flux_native,
                RESULTS)

            cross_flux_xf, _ = bin_channels_weighted(cross_flux_native, weight_uv, nchan_xf)
            cross_frac_xf, _ = bin_channels_weighted(cross_frac_native, weight_uv, nchan_xf)
            scan_xflux_map[(field_name, scan_num)] = {'flux': cross_flux_xf,
                                                      'frac': cross_frac_xf}

            # Native-resolution sqrt(Q²+U²) for S/N diagnostic panel
            cross_qu_native = np.sqrt(
                np.where(np.isfinite(Q_chan), Q_chan, 0.0)**2 +
                np.where(np.isfinite(U_chan), U_chan, 0.0)**2)
            qu_all_nan = (~np.isfinite(Q_chan)) & (~np.isfinite(U_chan))
            cross_qu_native[qu_all_nan] = np.nan
            scan_native_qu_map[(field_name, scan_num)] = cross_qu_native


# -----------------------------------------------------------------------
#  Build combine frac map (mean across scans per field)
# -----------------------------------------------------------------------

combine_frac_map = _build_field_mean_frac(scan_xflux_map)

# -----------------------------------------------------------------------
#  Build combine native-QU map (nanmean across scans per field)
# -----------------------------------------------------------------------

combine_native_qu_map = {}
for fn in field_info:
    scans_qu = [scan_native_qu_map[(fn, sn)]
                for (f2, sn) in scan_native_qu_map if f2 == fn]
    if scans_qu:
        stacked = np.array(scans_qu)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            combine_native_qu_map[fn] = np.nanmean(stacked, axis=0)

safe_ms_stem = myms.replace('/', '_').replace(' ', '_').rstrip('_')


# -----------------------------------------------------------------------
#  Process Xf tables.
#  perScan always runs first (full flagging pipeline).  combineScan is then
#  built as a cross-hand flux weighted phasor average of the flagged per-scan
#  solutions (mad_only — fine MAD + plot only).  If only one scan exists,
#  combineScan runs the full flagging pipeline directly.
#  combineScanField (do_multi_field): cross-hand flux weighted average across
#  fields written into the polcal template; fine MAD single-pass + plot.
#  Y-axis limits are captured from the first plot and reused in all
#  subsequent plots for direct comparison.
# -----------------------------------------------------------------------

def _run_process_combine(skip_flagging=False, phase_ylim=None, sn_ylim=None):
    msg()
    msg('=' * 60)
    msg('Processing combineScan Xf table')
    msg('=' * 60)
    if skip_flagging:
        msg('  (Cross-hand flux weighted phasor average of flagged perScan solutions)')
    else:
        msg('  (Single scan — full flagging pipeline)')
    return _process_xftab(
        tab_path          = xftab_combine,
        label             = 'combineScan',
        frac_lookup       = combine_frac_map,
        use_field_mean    = True,
        field_frac_map    = combine_frac_map,
        scan_qu_map       = scan_native_qu_map,
        flag_polcal_snap  = flag_polcal_combine,
        freq_ghz_xf       = freq_ghz_xf,
        plot_path         = xf_plot_dir + '/xf_phases_{}_combineScan.png'.format(safe_ms_stem),
        native_qu_lookup  = combine_native_qu_map,
        freq_ghz_native   = freq_ghz_native,
        skip_flagging     = skip_flagging,
        phase_ylim        = phase_ylim,
        sn_ylim           = sn_ylim,
    )

def _run_process_per():
    msg()
    msg('=' * 60)
    msg('Processing perScan Xf table')
    msg('=' * 60)
    return _process_xftab(
        tab_path          = xftab_per,
        label             = 'perScan',
        frac_lookup       = scan_xflux_map,      # per-scan for plot colours
        use_field_mean    = False,                # plot: label by scan
        field_frac_map    = combine_frac_map,     # per-field for flagging
        scan_qu_map       = scan_native_qu_map,
        flag_polcal_snap  = flag_polcal_per,
        freq_ghz_xf       = freq_ghz_xf,
        plot_path         = xf_plot_dir + '/xf_phases_{}_perScan.png'.format(safe_ms_stem),
        combine_xf_tab    = xftab_combine,
        native_qu_lookup  = scan_native_qu_map,  # per-scan for plot S/N panel
        freq_ghz_native   = freq_ghz_native,
    )

def _run_process_field_combine(phase_ylim=None, sn_ylim=None):
    msg()
    msg('=' * 60)
    msg('Processing combineScanField Xf table')
    msg('=' * 60)
    tb.open(xftab_field)
    flag_polcal_field = tb.getcol('FLAG').copy()
    tb.close()
    return _process_xftab(
        tab_path          = xftab_field,
        label             = 'combineScanField',
        frac_lookup       = combine_frac_map,
        use_field_mean    = True,
        field_frac_map    = combine_frac_map,
        scan_qu_map       = scan_native_qu_map,
        flag_polcal_snap  = flag_polcal_field,
        freq_ghz_xf       = freq_ghz_xf,
        plot_path         = xf_plot_dir + '/xf_phases_{}_combineScanField.png'.format(safe_ms_stem),
        native_qu_lookup  = combine_native_qu_map,
        freq_ghz_native   = freq_ghz_native,
        mad_only          = True,
        fine_sigma        = XF_FIELD_MAD_SIGMA,
        phase_ylim        = phase_ylim,
        sn_ylim           = sn_ylim,
    )

_phase_ylim = None
_sn_ylim    = None

if do_per_scan:
    _phase_ylim, _sn_ylim = _run_process_per()
    # Build combineScan from cross-hand flux weighted average of flagged per-scan solutions
    _avg_perscan_to_combine(xftab_per, xftab_combine)
    tb.open(xftab_combine)
    flag_polcal_combine = tb.getcol('FLAG').copy()
    tb.close()
    _run_process_combine(skip_flagging=True,
                         phase_ylim=_phase_ylim, sn_ylim=_sn_ylim)
else:
    # Single scan — combineScan IS the per-scan solution; full flagging pipeline runs
    _phase_ylim, _sn_ylim = _run_process_combine()

if do_multi_field:
    # Build per-channel cross-hand flux weights (field-mean QU binned to Xf grid)
    _nchan_xf   = len(freq_ghz_xf)
    _field_qu_w = {}
    for _fn, _qu_nat in combine_native_qu_map.items():
        _qu_nat = np.asarray(_qu_nat, dtype=float)
        _w_nat  = np.where(np.isfinite(_qu_nat), 1.0, 0.0)
        _qu_xf, _ = bin_channels_weighted(_qu_nat, _w_nat, _nchan_xf)
        _field_qu_w[_fn] = np.where(np.isfinite(_qu_xf) & (_qu_xf > 0.0),
                                     _qu_xf, 0.0)
    # Overwrite polcal template with cross-hand flux weighted cross-field average
    _avg_combine_fields_to_global(xftab_combine, xftab_field, _field_qu_w)
    _run_process_field_combine(phase_ylim=_phase_ylim, sn_ylim=_sn_ylim)


msg()
msg('Xf solver complete.')

