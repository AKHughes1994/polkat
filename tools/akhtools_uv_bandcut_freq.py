#!/usr/bin/env python3
# andrew.hughes@physics.ox.ac.uk

"""
Suggest circular uv cuts (WSClean -minuv-l and -maxuv-l) that give the same uv
coverage at every frequency across the band.

The uv coordinates of a baseline in wavelengths scale with frequency, so the
bottom of the band has the smallest outer extent and the top of the band the
largest central hole. Cutting every channel to the ring both ends share makes
the PSF, and the range of angular scales sampled, the same at every frequency.
This matters for spectropolarimetry, where frequency-dependent uv coverage
shows up as spurious structure in Q, U and the RM spectrum.

Method. The unflagged visibilities in the bottom and top EDGE_FRAC of the
unflagged channels (every channel with at least one unflagged visibility;
fully flagged channels, at the band edges or inside the band, are skipped) are
taken, with their Hermitian mirrors (-u, -v), and

    outer cut   the semi-minor axis of the smallest ellipse enclosing the
                bottom-of-band uv points (Khachiyan's algorithm, centred on the
                origin). The longest unflagged bottom-of-band baseline is
                reported as an alternative
    inner cut   by default the band's PolKAT default minimum baseline
                (INNER_CUT_DEFAULT_M: 500 m at UHF and L, 100 m at S-band,
                where RFI is much lighter), converted to wavelengths at the
                band centre as oxkat/config.py does for WSC_MINUVL. Also
                reported: the calculated cut, the shortest unflagged
                baseline at the top of the band (the largest circle about
                the origin with no top-of-band data in it). INNER_CUT picks
                which is used; MINUVL_OVERRIDE replaces both

For each end of the band it reports how much of the unflagged data each cut
removes, on top of what is already flagged.

A visibility counts as unflagged only when none of its correlations is
flagged. Autocorrelations (zero uv length) are left out throughout.

Outputs, in directories created under the current working directory

    PLOTS/akhtools_uv_bandcut_<ms>_<field>.png
        2 x 2 panels. Columns are the bottom and top of the band, the top row
        the full uv plane and the bottom row a zoom on the central hole.
        grey, log scale     unflagged uv density (with mirrors)
        green, shaded       the kept uv range, between the two cuts
        green, dashed       outer cut (-maxuv-l)
        grey, dashed        alternative outer cut, longest bottom baseline
        orange, dashed      calculated inner cut (top-of-band shortest baseline)
        purple, dashed      band default inner cut (INNER_CUT_DEFAULT_M)
        black, dotted       MINUVL_OVERRIDE, if set
        the green shading stops at whichever inner cut is used
        blue, solid         enclosing ellipse, bottom of band
        red, solid          enclosing ellipse, top of band
    LOGS/akhtools_uv_bandcut_<field>_<timestamp>.log
        everything printed to the screen

All inputs are globals in the CONFIGURATION block below. Edit them and run
    python akhtools_uv_bandcut_freq.py
"""

import warnings
import os

warnings.filterwarnings('ignore')
os.environ['PYTHONWARNINGS'] = 'ignore'

import logging
import re
from datetime import datetime

import numpy as np
import dask
from daskms import xds_from_ms, xds_from_table
from scipy.spatial import ConvexHull, QhullError

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# =====================================================================
# CONFIGURATION (edit here)
# =====================================================================

# Input
MS_FILE = '/path/to/target.ms'
FIELD = None                      # field ID (int) or name; None uses the first field

# Band edges
EDGE_FRAC = 0.1                   # fraction of the unflagged channels taken at each end

# Inner cut. INNER_CUT picks which is used: 'default', the band's PolKAT
# default minimum baseline below (as oxkat/config.py WSC_BASELINE_CUTLENGTH),
# or 'calculated', the shortest unflagged baseline at the top of the band.
# Both are always reported and drawn.
INNER_CUT = 'default'
INNER_CUT_DEFAULT_M = {'UHF': 500.0, 'L': 500.0, 'S': 100.0}   # metres, per band

# Inner cut in wavelengths that replaces both, or None
MINUVL_OVERRIDE = None

# Plot
BINS = 800                        # uv density grid size per axis

# Output
SAVE_DIR = 'PLOTS/'
LOG_DIR = 'LOGS/'
PRINTING = True                   # False logs to file only

C_LIGHT = 299792458.0             # m/s


# =====================================================================
# LOGGING
# =====================================================================

def log(message=""):
    """Timestamped line to the screen (if PRINTING) and always to the log file."""
    logging.info(message)
    if PRINTING:
        if message == "":
            print(flush=True)
        else:
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def scrub_name(name):
    return re.sub(r'[^A-Za-z0-9_.+-]', '_', str(name))


def meerkat_band(freqs):
    """MeerKAT receiver band from the channel frequencies (Hz): UHF, L or S."""
    centre = 0.5 * (np.min(freqs) + np.max(freqs))
    if centre < 1.05e9:
        return 'UHF'
    if centre < 1.75e9:
        return 'L'
    return 'S'


# =====================================================================
# SMALLEST ENCLOSING ELLIPSE (Khachiyan, centred on the origin)
# =====================================================================
# The uv points always include their Hermitian mirrors (-u, -v), so the set is
# symmetric about the origin and so is its smallest enclosing ellipse. That
# ellipse, x^T A x <= 1, is found with Khachiyan's algorithm, which converges
# to the exact minimum-area solution with no recursion. (tools/uv_band_cut.py
# used Welzl's algorithm, as in tools/homogenize_beams.py; on some uv
# coverages, e.g. SS433's, it returns a degenerate ellipse with a semi-major
# axis of ~1e12 lambda, and elsewhere one that is not the smallest, which
# moves the outer cut.)
# An ellipse is (centre, semi-major, semi-minor, rotation of the major axis).

def center_form_to_geometric(F, c):
    """(centre, semi_major, semi_minor, rotation) from x^T F x <= 1 about c, or None."""
    try:
        w, V = np.linalg.eigh(F)
        if np.any(w <= 0):
            return None
        if V[1, 0] < 0:
            V[:, 0] = -V[:, 0]
        rotation_angle = np.arctan2(V[1, 0], V[0, 0])
        semi_major = 1.0 / np.sqrt(np.min(w))
        semi_minor = 1.0 / np.sqrt(np.max(w))
        return c, semi_major, semi_minor, rotation_angle
    except (np.linalg.LinAlgError, ValueError, ZeroDivisionError):
        return None


def enclosing_ellipse(points, tol=1e-10, max_iter=50000):
    """
    Minimum-area ellipse centred on the origin enclosing `points` (n, 2), or
    None if the points do not span two dimensions.

    Khachiyan's algorithm: weights u on the points define X = sum u_i x_i x_i^T;
    each step moves weight onto the point furthest outside the current ellipse
    x^T (X^-1 / 2) x <= 1 until the weights stop changing. The final ellipse is
    scaled up by the residual, if any, so every point is strictly inside.
    """
    P = np.asarray(points, dtype=float).T
    n = P.shape[1]
    if n < 2 or np.linalg.matrix_rank(P) < 2:
        return None
    u = np.full(n, 1.0 / n)
    for _ in range(max_iter):
        X = (P * u) @ P.T
        M = np.einsum('ij,ij->j', P, np.linalg.solve(X, P))   # x_i^T X^-1 x_i
        j = np.argmax(M)
        step = (M[j] - 2.0) / (2.0 * (M[j] - 1.0))
        new_u = (1.0 - step) * u
        new_u[j] += step
        converged = np.abs(new_u - u).max() < tol
        u = new_u
        if converged:
            break
    A = np.linalg.inv((P * u) @ P.T) / 2.0
    A /= max(1.0, np.einsum('ij,ij->j', P, A @ P).max())
    return center_form_to_geometric(A, np.zeros(2))


def ellipse_outline(ellipse, num_pts=400):
    """(num_pts, 2) points around an ellipse, for plotting."""
    c, a, b, t = ellipse
    theta = np.linspace(0, 2 * np.pi, num_pts)
    major = np.array([np.cos(t), np.sin(t)])
    minor = np.array([-np.sin(t), np.cos(t)])
    return c + np.outer(a * np.cos(theta), major) + np.outer(b * np.sin(theta), minor)


# =====================================================================
# MS READING
# =====================================================================
# dask-ms gives lazy dask arrays (ds.COLUMN.data); nothing is read until
# dask.compute(). Only .data and dask.compute() are used, never xarray methods
# such as .values, so this runs with or without xarray installed.

def _decode(values):
    return [v.decode() if isinstance(v, (bytes, np.bytes_)) else str(v) for v in values]


def resolve_field(ms, field_arg):
    """(field_id, field_name). With FIELD = None: the only field, or the first with a warning."""
    t = xds_from_table(f'{ms}::FIELD', columns=['NAME'])[0]
    names = _decode(t.NAME.data.compute())
    listing = ', '.join(f'{i}={n}' for i, n in enumerate(names))

    if field_arg is None:
        if len(names) > 1:
            log(f'WARNING: FIELD is None and the MS has {len(names)} fields '
                f'({listing}); using the first, {names[0]}')
        return 0, names[0]
    if str(field_arg).isdigit():
        fid = int(field_arg)
        if not 0 <= fid < len(names):
            raise ValueError(f'Field ID {fid} out of range. Available fields: {listing}')
        return fid, names[fid]
    for i, name in enumerate(names):
        if name.strip().lower() == str(field_arg).strip().lower():
            return i, name
    raise ValueError(f"Field '{field_arg}' not found. Available fields: {listing}")


def load_field(ms, field_id):
    """
    UVW (row, 3) in metres, per-channel 'any correlation flagged' (row, chan),
    and the channel frequencies (Hz) of one field's spectral window.
    """
    # One dataset per spectral window (DATA_DESC_ID). A grouping column is not
    # a data column: dask-ms exposes it on the dataset as a plain integer
    ds_list = xds_from_ms(ms, columns=['UVW', 'FLAG'],
                          taql_where=f'FIELD_ID=={field_id}',
                          group_cols=['DATA_DESC_ID'], chunks={'row': 100_000})
    if not ds_list:
        raise RuntimeError(f'No rows for FIELD_ID {field_id}')
    if len(ds_list) > 1:
        log(f'WARNING: {len(ds_list)} spectral windows; using the first only')
    ds = ds_list[0]

    # Reduce FLAG over correlations lazily, so only (row, chan) booleans are
    # held in memory rather than the full (row, chan, corr) array
    uvw, flagged = dask.compute(ds.UVW.data, ds.FLAG.data.any(axis=2))
    ddid = int(ds.DATA_DESC_ID)

    dd = xds_from_table(f'{ms}::DATA_DESCRIPTION', columns=['SPECTRAL_WINDOW_ID'])[0]
    spw_id = int(dd.SPECTRAL_WINDOW_ID.data.compute()[ddid])
    spw = xds_from_table(f'{ms}::SPECTRAL_WINDOW', columns=['CHAN_FREQ'])[0]
    freqs = np.asarray(spw.CHAN_FREQ.data.compute())[spw_id]
    return np.asarray(uvw), np.asarray(flagged), freqs


# =====================================================================
# UV GEOMETRY
# =====================================================================

def block_points(uvw, flagged, freqs, chans):
    """Per channel in `chans`, its unflagged (u, v) in wavelengths."""
    for ch in chans:
        ok = ~flagged[:, ch]
        yield ch, uvw[ok, :2] * (freqs[ch] / C_LIGHT)


def hull_vertices(points):
    """Convex-hull vertices of `points` (all of them if there are too few for a hull)."""
    if len(points) < 3:
        return points
    try:
        return points[ConvexHull(points).vertices]
    except QhullError:
        return points


def block_ellipse(uvw, flagged, freqs, chans):
    """
    Smallest ellipse enclosing a channel block's unflagged uv points and their
    mirrors. Only convex-hull vertices can touch the enclosing ellipse, so
    each channel is reduced to its hull first, then the block's combined hull.
    """
    verts = [hull_vertices(uv) for _, uv in block_points(uvw, flagged, freqs, chans) if len(uv)]
    if not verts:
        return None
    verts = np.vstack(verts)
    verts = hull_vertices(np.vstack([verts, -verts]))
    log(f'    {len(verts)} hull vertices')
    return enclosing_ellipse(verts)


def baseline_lengths(uvw, flagged, freqs, chans):
    """Unflagged uv distances (wavelengths) of every channel in `chans`, one array per channel."""
    return [np.hypot(uv[:, 0], uv[:, 1]) for _, uv in block_points(uvw, flagged, freqs, chans)]


# =====================================================================
# PLOT
# =====================================================================

def plot_coverage(uvw, flagged, freqs, bottom, top, ell_bottom, ell_top,
                  max_cut, max_cut_uvmax, min_cut, inner_cuts, title, path):
    """
    Unflagged uv density at each end of the band with the cuts and ellipses
    (see preamble). `min_cut` is the inner cut used; `inner_cuts` maps each
    candidate inner cut's name to (lambda, colour, linestyle).
    """
    # Full-plane view sized by the longest unflagged baseline at the top of the
    # band (the largest uv distance at either end), widened to show the
    # ellipses if they reach further, but never by more than 2x, so a bad fit
    # cannot blow the axes up
    uvmax_top = max(r.max() for r in baseline_lengths(uvw, flagged, freqs, top) if len(r))
    reach = [e[1] for e in (ell_bottom, ell_top) if e is not None and e[1] <= 2 * uvmax_top]
    extent = 1.05 * max([uvmax_top] + reach) / 1e3
    zoom = 4.0 * max(v[0] for v in inner_cuts.values()) / 1e3
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    circle = np.linspace(0, 2 * np.pi, 400)

    for row, lim in ((0, extent), (1, zoom)):
        edges = np.linspace(-lim, lim, BINS + 1)
        for col, (label, chans) in enumerate((('Bottom', bottom), ('Top', top))):
            ax = axes[row, col]

            # uv density in klambda, each visibility counted with its mirror
            hist = np.zeros((BINS, BINS))
            for _, uv in block_points(uvw, flagged, freqs, chans):
                uv = uv / 1e3
                hist += np.histogram2d(uv[:, 0], uv[:, 1], bins=[edges, edges])[0]
                hist += np.histogram2d(-uv[:, 0], -uv[:, 1], bins=[edges, edges])[0]
            hist[hist == 0] = np.nan
            ax.imshow(hist.T, origin='lower', extent=[-lim, lim, -lim, lim],
                      cmap='Greys', norm=LogNorm(), interpolation='nearest')

            # Kept uv range: the ring between the cuts (inner circle traced backwards)
            ring_u = np.concatenate([max_cut * np.cos(circle), min_cut * np.cos(circle[::-1])]) / 1e3
            ring_v = np.concatenate([max_cut * np.sin(circle), min_cut * np.sin(circle[::-1])]) / 1e3
            ax.fill(ring_u, ring_v, color='tab:green', alpha=0.15, lw=0, label='kept uv range')

            for ell, colour, name in ((ell_bottom, 'tab:blue', 'bottom ellipse'),
                                      (ell_top, 'tab:red', 'top ellipse')):
                if ell is not None:
                    outline = ellipse_outline(ell) / 1e3
                    ax.plot(outline[:, 0], outline[:, 1], color=colour, lw=1.5, label=name)

            cuts = [(max_cut, 'tab:green', '--', f'maxuv {max_cut/1e3:.2f} k$\\lambda$'),
                    (max_cut_uvmax, 'grey', '--', f'bottom uvmax {max_cut_uvmax/1e3:.2f} k$\\lambda$')]
            for name, (radius, colour, ls) in inner_cuts.items():
                used = ' (used)' if radius == min_cut else ''
                cuts.append((radius, colour, ls, f'minuv {name} {radius/1e3:.3f} k$\\lambda${used}'))
            for radius, colour, ls, name in cuts:
                ax.plot(radius / 1e3 * np.cos(circle), radius / 1e3 * np.sin(circle),
                        color=colour, lw=1.5, ls=ls, label=name)

            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_title(f'{label} {len(chans)} channel(s): '
                         f'{freqs[chans[0]]/1e9:.3f}-{freqs[chans[-1]]/1e9:.3f} GHz'
                         + (' (inner zoom)' if row == 1 else ''))
            ax.set_xlabel('u (k$\\lambda$)')
            ax.set_ylabel('v (k$\\lambda$)')
            ax.set_aspect('equal')
            ax.legend(loc='upper right')

    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=150)
    plt.close(fig)


# =====================================================================
# MAIN
# =====================================================================

def main():
    ms = MS_FILE.rstrip('/')
    field_id, field_name = resolve_field(ms, FIELD)
    tag = scrub_name(field_name)

    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(SAVE_DIR, exist_ok=True)
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(LOG_DIR, f'akhtools_uv_bandcut_{tag}_{stamp}.log')
    logging.basicConfig(filename=log_file, level=logging.INFO, format='[%(asctime)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S', filemode='w', force=True)

    log("=" * 70)
    log("UV BAND CUT")
    log("=" * 70)
    log(f"  MS: {ms}")
    log(f"  Field: {field_name} (ID {field_id})")
    log(f"  Edge fraction: {EDGE_FRAC}")
    log(f"  Inner cut: {INNER_CUT}"
        + (f"; MINUVL_OVERRIDE = {MINUVL_OVERRIDE} lambda" if MINUVL_OVERRIDE is not None else ''))
    log(f"  Log: {log_file}")
    log("")

    # ------- Read UVW and flags; drop autocorrelations (zero uv length)
    uvw, flagged, freqs = load_field(ms, field_id)
    cross = np.any(uvw[:, :2] != 0, axis=1)
    if not cross.all():
        log(f'Dropping {np.count_nonzero(~cross):,} autocorrelation row(s)')
        uvw, flagged = uvw[cross], flagged[cross]
    nrow, nchan = flagged.shape
    log(f'{nrow:,} rows, {nchan} channels, {freqs[0]/1e9:.4f}-{freqs[-1]/1e9:.4f} GHz')

    # ------- Edge blocks: EDGE_FRAC of the unflagged channels at each end of the
    # band. Fully flagged channels (band edges, RFI) are skipped, so each block
    # holds n_edge channels with data rather than n_edge of the channel range
    good = np.where((~flagged).any(axis=0))[0]
    if len(good) == 0:
        raise RuntimeError('Every channel is fully flagged for this field')
    c0, c1 = good[0], good[-1]
    n_edge = max(1, int(round(EDGE_FRAC * len(good))))
    bottom = good[:n_edge]
    top = good[-n_edge:]
    log(f'Unflagged channels: {len(good)} of {nchan}, between channels {c0} and {c1} '
        f'({freqs[c0]/1e9:.4f}-{freqs[c1]/1e9:.4f} GHz); {n_edge} at each end')
    log(f'  Bottom: channels {bottom[0]}-{bottom[-1]} ({freqs[bottom[0]]/1e9:.4f}-{freqs[bottom[-1]]/1e9:.4f} GHz)')
    log(f'  Top:    channels {top[0]}-{top[-1]} ({freqs[top[0]]/1e9:.4f}-{freqs[top[-1]]/1e9:.4f} GHz)')

    # ------- Outer cut: the bottom-of-band enclosing ellipse's semi-minor axis,
    # the largest circle about the centre that the bottom of the band fills
    log('Fitting enclosing ellipses...')
    log('  Bottom of band:')
    ell_bottom = block_ellipse(uvw, flagged, freqs, bottom)
    log('  Top of band:')
    ell_top = block_ellipse(uvw, flagged, freqs, top)
    if ell_bottom is None:
        raise RuntimeError('Could not fit an ellipse to the bottom of the band')

    max_cut = ell_bottom[2]
    # An origin-centred enclosing ellipse cannot have a semi-minor axis longer
    # than the longest baseline it encloses; flag it if it ever does
    # Alternative outer cut: the longest unflagged bottom-of-band baseline,
    # which keeps every baseline the bottom of the band has
    max_cut_uvmax = max(r.max() for r in baseline_lengths(uvw, flagged, freqs, bottom) if len(r))
    if max_cut > 1.001 * max_cut_uvmax:
        log(f'WARNING: bottom ellipse semi-minor ({max_cut:.0f} lambda) exceeds the longest '
            f'bottom baseline ({max_cut_uvmax:.0f} lambda); check the fit')
    for label, ell in (('Bottom', ell_bottom), ('Top', ell_top)):
        if ell is None:
            log(f'  {label}: ellipse fit failed')
            continue
        c, a, b, t = ell
        log(f'  {label}: semi-major {a/1e3:.2f} klambda, semi-minor {b/1e3:.2f} klambda, '
            f'PA {np.degrees(t):.1f} deg, centre ({c[0]:.1f}, {c[1]:.1f}) lambda')

    # ------- Inner cut candidates
    # calculated: the shortest unflagged top-of-band baseline, the largest hole
    # about the origin the top of the band leaves
    min_cut_calc = min((r.min() for r in baseline_lengths(uvw, flagged, freqs, top) if len(r)), default=None)
    if min_cut_calc is None:
        raise RuntimeError('No unflagged data at the top of the band')
    log(f'  Top: shortest unflagged baseline {min_cut_calc:.0f} lambda')

    # default: the band's minimum baseline in metres, in wavelengths at the
    # band centre (frequency-dependent, as oxkat/config.py does for WSC_MINUVL)
    band = meerkat_band(freqs)
    centre = 0.5 * (freqs.min() + freqs.max())
    default_m = INNER_CUT_DEFAULT_M[band]
    min_cut_default = default_m * centre / C_LIGHT
    log(f'  {band}-band default: {default_m:.0f} m at {centre/1e9:.4f} GHz = {min_cut_default:.0f} lambda')

    inner_cuts = {'calculated': (min_cut_calc, 'tab:orange', '--'),
                  'default': (min_cut_default, 'tab:purple', '--')}
    if MINUVL_OVERRIDE is not None:
        inner_cuts['override'] = (float(MINUVL_OVERRIDE), 'black', ':')
        inner_used = 'override'
    elif INNER_CUT in ('default', 'calculated'):
        inner_used = INNER_CUT
    else:
        raise ValueError(f"INNER_CUT must be 'default' or 'calculated', not {INNER_CUT!r}")
    min_cut = inner_cuts[inner_used][0]
    log(f'  Inner cut used: {inner_used} ({min_cut:.0f} lambda)')

    # ------- What each cut removes at each end of the band, on top of the
    # existing flags (fractions of the unflagged data)
    log('')
    log('=' * 70)
    for label, chans in (('Bottom', bottom), ('Top', top)):
        n_total = nrow * len(chans)
        n_unflagged = n_inner = n_outer = n_outer_uvmax = 0
        n_inner_each = {name: 0 for name in inner_cuts}
        for r in baseline_lengths(uvw, flagged, freqs, chans):
            n_unflagged += len(r)
            n_inner += int(np.count_nonzero(r < min_cut))
            for name, (radius, _, _) in inner_cuts.items():
                n_inner_each[name] += int(np.count_nonzero(r < radius))
            n_outer += int(np.count_nonzero(r > max_cut))
            n_outer_uvmax += int(np.count_nonzero(r > max_cut_uvmax))
        frac_flagged = 1.0 - n_unflagged / n_total
        frac_inner = n_inner / n_unflagged if n_unflagged else 0.0
        frac_outer = n_outer / n_unflagged if n_unflagged else 0.0
        frac_outer_uvmax = n_outer_uvmax / n_unflagged if n_unflagged else 0.0
        frac_both = frac_inner + frac_outer
        log(f'{label} {len(chans)} channel(s), {n_total:,} visibilities:')
        log(f'  Already flagged:                      {frac_flagged:.2%}')
        log(f'  Unflagged inside -minuv-l:            {n_inner:,} of {n_unflagged:,} ({frac_inner:.2%})')
        log(f'  Unflagged beyond -maxuv-l:            {n_outer:,} of {n_unflagged:,} ({frac_outer:.2%})')
        log(f'  Flagged + both cuts (total lost):     {frac_flagged + (1 - frac_flagged) * frac_both:.2%}')
        log(f'  Unflagged beyond bottom uvmax cut:    {n_outer_uvmax:,} of {n_unflagged:,} ({frac_outer_uvmax:.2%})')
        for name in inner_cuts:
            frac = n_inner_each[name] / n_unflagged if n_unflagged else 0.0
            log(f'  Unflagged inside {name + " inner cut:":22s}{n_inner_each[name]:,} of {n_unflagged:,} ({frac:.2%})')
        log('')
    log(f'Outer cut (bottom-of-band semi-minor axis): {max_cut:.0f} lambda ({max_cut/1e3:.2f} klambda)')
    log(f'Outer cut at the bottom-of-band uvmax:      {max_cut_uvmax:.0f} lambda ({max_cut_uvmax/1e3:.2f} klambda)')
    log(f'Inner cut (top-of-band shortest baseline):  {min_cut_calc:.0f} lambda ({min_cut_calc/1e3:.3f} klambda)')
    log(f'Inner cut ({band}-band default, {default_m:.0f} m):'.ljust(44)
        + f' {min_cut_default:.0f} lambda ({min_cut_default/1e3:.3f} klambda)')
    if MINUVL_OVERRIDE is not None:
        log('Inner cut (MINUVL_OVERRIDE):'.ljust(44) + f' {min_cut:.0f} lambda ({min_cut/1e3:.3f} klambda)')
    log(f'Inner cut used: {inner_used}')
    log(f'  -minuv-l {min_cut:.0f} -maxuv-l {max_cut:.0f}')
    log('=' * 70)

    # ------- Diagnostic plot
    png = os.path.join(SAVE_DIR, f'akhtools_uv_bandcut_{scrub_name(os.path.basename(ms))}_{tag}.png')
    plot_coverage(uvw, flagged, freqs, bottom, top, ell_bottom, ell_top,
                  max_cut, max_cut_uvmax, min_cut, inner_cuts,
                  f'{os.path.basename(ms)} -- {field_name}: unflagged uv coverage', png)
    log(f'Plot: {png}')
    log('Done')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        log(f"ERROR: {exc}")
        raise
