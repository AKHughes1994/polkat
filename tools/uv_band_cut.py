#!/usr/bin/env python3
# andrew.hughes@physics.ox.ac.uk

"""
Suggest circular -maxuv-l and -minuv-l cuts that match the uv coverage across
the band, at both the smallest and the largest angular scales.

The uv coverage in wavelengths scales with frequency, so the bottom of the band
has the smallest outer extent and the top of the band the largest inner hole.
This tool takes the unflagged visibilities in the bottom and top EDGE_FRAC of
the usable band (the channel range left once the fully flagged edge channels
are dropped) and:
  - maxuv: finds the smallest ellipse enclosing each end (Welzl algorithm, as
    in tools/homogenize_beams.py) and takes the semi-minor axis of the
    bottom-of-band ellipse as the outer cut;
  - minuv: takes the shortest unflagged baseline at the top of the band, the
    largest circle around the origin with no top-of-band data in it, as the
    inner cut, unless WSC_MINUVL is set in oxkat/config.py, in which case that
    is used instead (the calculated value is still printed).
It reports how much of the unflagged data at each end of the band each cut
removes, and plots the uv coverage of both ends with the ellipses and cuts.

A visibility counts as unflagged only when none of its correlations is
flagged. Autocorrelations (zero baseline length) are left out throughout.

Usage:
    python3 uv_band_cut.py <ms> [--field NAME_OR_ID] [--edge-frac 0.1] [--outdir DIR] [--bins 800]
Any option not given takes its value from the globals below.
"""

import argparse
import os
import sys
import threading

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.spatial import ConvexHull, QhullError
from daskms import xds_from_ms, xds_from_table

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

C_LIGHT = 299792458.0

FIELD = None      # Field name or ID; None uses the first field (warning if the MS has several)
EDGE_FRAC = 0.1   # Fraction of the usable band taken at each end
OUTDIR = '.'      # Directory for the diagnostic plot
BINS = 800        # uv plot grid size per axis


def msg(txt):
    print(txt, flush=True)


# ============================================================================
# WELZL ALGORITHM (copied from tools/homogenize_beams.py)
# ============================================================================

def is_singular(A, tolerance=1e-12):
    try:
        _, s, _ = np.linalg.svd(A)
        return np.min(s) < tolerance * np.max(s)
    except np.linalg.LinAlgError:
        return True


def center_form_to_geometric(F, c):
    """(center, semi_major, semi_minor, rotation of the major axis) or None."""
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


def ellipse_from_boundary5(S):
    if S.shape != (5, 2):
        return None
    try:
        x, y = S[:, 0], S[:, 1]
        A = np.column_stack((x**2, y**2, 2 * x * y, x, y))
        if is_singular(A):
            return None
        sol = np.linalg.solve(A, -np.ones(S.shape[0]))
        center_matrix = -2 * np.array([[sol[0], sol[2]], [sol[2], sol[1]]])
        if is_singular(center_matrix):
            return None
        c = np.linalg.solve(center_matrix, sol[3:5])
        A_ext = np.vstack([
            np.hstack([np.eye(3), -np.array([[sol[0], sol[2], sol[1]]]).T]),
            np.array([c[0]**2, 2 * c[0] * c[1], c[1]**2, -1])
        ])
        s = np.linalg.solve(A_ext, np.array([0, 0, 0, 1]))
        F = np.array([[s[0], s[1]], [s[1], s[2]]])
        return center_form_to_geometric(F, c)
    except (np.linalg.LinAlgError, ValueError):
        return None


def ellipse_from_boundary4(S):
    if S.shape != (4, 2):
        return None
    try:
        Sc = S - np.mean(S, axis=0)
        angles = np.arctan2(Sc[:, 1], Sc[:, 0])
        S = S[np.argsort(-angles), :]

        A = np.column_stack([S[2, :] - S[0, :], S[1, :] - S[3, :]])
        if is_singular(A):
            return None
        b = S[1, :] - S[0, :]
        s = np.linalg.solve(A, b)
        diag_intersect = S[0, :] + s[0] * (S[2, :] - S[0, :])
        S = S - diag_intersect

        AC = S[2, :] - S[0, :]
        theta = np.arctan2(AC[1], AC[0])
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        rot_mat = np.array([[cos_theta, sin_theta], [-sin_theta, cos_theta]])
        S = (rot_mat @ S.T).T

        denominator = S[3, 1] - S[1, 1]
        if abs(denominator) < 1e-12:
            return None
        m = (S[1, 0] - S[3, 0]) / denominator
        shear_mat = np.array([[1.0, m], [0.0, 1.0]])
        S = (shear_mat @ S.T).T

        b = np.linalg.norm(S, axis=1)
        if np.any(b == 0):
            return None
        d = (b[1] * b[3]) / (b[2] * b[0])
        if d <= 0:
            return None
        d_quarter = d ** 0.25
        stretch_mat = np.diag([d_quarter, 1.0 / d_quarter])
        S = (stretch_mat @ S.T).T

        a = np.linalg.norm(S, axis=1)
        coeff = np.zeros(4)
        coeff[0] = -4 * a[1]**2 * a[2] * a[0]
        coeff[1] = -4 * a[1] * (a[2] - a[0]) * (a[1]**2 - a[2] * a[0])
        coeff[2] = (3 * a[1]**2 * (a[1]**2 + a[2]**2) -
                    8 * a[1]**2 * a[2] * a[0] +
                    3 * (a[1]**2 + a[2]**2) * a[0]**2)
        coeff[3] = coeff[1] / 2.0
        roots = np.roots(coeff)
        valid_roots = roots[(-1 < np.real(roots)) & (np.real(roots) < 1) & (np.abs(np.imag(roots)) < 1e-10)]
        if len(valid_roots) == 0:
            return None
        theta = np.arcsin(np.real(valid_roots[0]))

        cos_theta_sqrt = np.cos(theta) ** 0.5
        if cos_theta_sqrt == 0:
            return None
        D_mat = np.array([[1.0 / cos_theta_sqrt, np.sin(theta) / cos_theta_sqrt],
                          [0.0, cos_theta_sqrt]])
        S = (D_mat @ S.T).T

        boundary = S[:-1, :]
        A = np.vstack([-2 * boundary.T, np.ones(boundary.shape[0])]).T
        if is_singular(A):
            return None
        b = -np.sum(boundary**2, axis=1)
        s = np.linalg.solve(A, b)
        circle_c = s[:2]
        circle_r_sq = np.sum(circle_c**2) - s[2]
        if circle_r_sq <= 0:
            return None
        circle_r = np.sqrt(circle_r_sq)

        T_mat = D_mat @ stretch_mat @ shear_mat @ rot_mat
        ellipse_c = np.linalg.solve(T_mat, circle_c) + diag_intersect
        ellipse_F = T_mat.T @ T_mat / (circle_r**2)
        return center_form_to_geometric(ellipse_F, ellipse_c)
    except (np.linalg.LinAlgError, ValueError, ZeroDivisionError):
        return None


def ellipse_from_boundary3(S):
    if S.shape != (3, 2):
        return None
    try:
        c = np.mean(S, axis=0)
        Sc = S - c
        if is_singular(Sc):
            return None
        StS = Sc.T @ Sc
        if is_singular(StS):
            return None
        F = 1.5 * np.linalg.inv(StS)
        return center_form_to_geometric(F, c)
    except (np.linalg.LinAlgError, ValueError):
        return None


def welzl(interior, boundary=None):
    if boundary is None:
        boundary = np.zeros((0, 2))

    if interior.shape[0] == 0 or boundary.shape[0] == 5:
        if boundary.shape[0] <= 2:
            return None
        elif boundary.shape[0] == 3:
            return ellipse_from_boundary3(boundary)
        elif boundary.shape[0] == 4:
            return ellipse_from_boundary4(boundary)
        else:
            return ellipse_from_boundary5(boundary)

    i = np.random.randint(interior.shape[0])
    p = interior[i, :]
    interior_wo_p = np.delete(interior, i, 0)

    ellipse = welzl(interior_wo_p, boundary)
    if is_in_ellipse(p, ellipse):
        return ellipse
    new_boundary = np.vstack([boundary, p]) if boundary.size > 0 else p.reshape(1, -1)
    return welzl(interior_wo_p, new_boundary)


def is_in_ellipse(point, ellipse):
    if ellipse is None or point is None:
        return False
    try:
        c, a, b, t = ellipse
        v = point - c
        cos_t, sin_t = np.cos(t), np.sin(t)
        rot_mat = np.array([[cos_t, sin_t], [-sin_t, cos_t]])
        F = rot_mat.T @ np.diag([1.0/(a**2), 1.0/(b**2)]) @ rot_mat
        return v.T @ F @ v <= 1.0 + 1e-10
    except (ValueError, TypeError, IndexError):
        return False


def enclosing_ellipse(points):
    """
    Smallest ellipse enclosing `points`. Welzl recurses once per point, so it
    runs on a thread with a large stack and a raised recursion limit.
    """
    np.random.seed(0)
    result = {}

    def run():
        result['ellipse'] = welzl(points)

    sys.setrecursionlimit(max(sys.getrecursionlimit(), 10 * len(points) + 1000))
    threading.stack_size(512 * 1024 * 1024)
    worker = threading.Thread(target=run)
    worker.start()
    worker.join()
    return result.get('ellipse')


def ellipse_outline(ellipse, num_pts=400):
    c, a, b, t = ellipse
    theta = np.linspace(0, 2 * np.pi, num_pts)
    major = np.array([np.cos(t), np.sin(t)])
    minor = np.array([-np.sin(t), np.cos(t)])
    return c + np.outer(a * np.cos(theta), major) + np.outer(b * np.sin(theta), minor)


# ============================================================================
# MS READING
# ============================================================================

def resolve_field(ms, field_arg):
    """(field_id, field_name). With no FIELD: the only field, or the first with a warning."""
    names = xds_from_table(f'{ms}::FIELD', columns=['NAME'])[0].NAME.values
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in names]

    if field_arg is None:
        if len(names) > 1:
            msg(f'WARNING: no FIELD given and the MS has {len(names)} fields '
                f'({", ".join(names)}) -- using the first, {names[0]}')
        return 0, names[0]
    if str(field_arg).isdigit():
        return int(field_arg), names[int(field_arg)]
    if field_arg not in names:
        sys.exit(f'ERROR: field {field_arg} not in the MS (fields: {", ".join(names)})')
    return names.index(field_arg), field_arg


def load_field(ms, field_id):
    """UVW (rows x 3) and per-channel 'any correlation flagged' (rows x chans) for one field."""
    ds_list = xds_from_ms(ms, columns=['UVW', 'FLAG'], taql_where=f'FIELD_ID=={field_id}',
                          group_cols=['DATA_DESC_ID'], chunks={'row': 100_000})
    if not ds_list:
        sys.exit(f'ERROR: no rows for FIELD_ID {field_id}')
    if len(ds_list) > 1:
        msg(f'WARNING: {len(ds_list)} spectral windows -- using the first only')
    ds = ds_list[0]

    ddid = ds.attrs['DATA_DESC_ID']
    spw_id = int(xds_from_table(f'{ms}::DATA_DESCRIPTION', columns=['SPECTRAL_WINDOW_ID'])[0]
                 .SPECTRAL_WINDOW_ID.values[ddid])
    freqs = xds_from_table(f'{ms}::SPECTRAL_WINDOW', columns=['CHAN_FREQ'])[0].CHAN_FREQ.values[spw_id]

    uvw = ds.UVW.data.compute()
    flagged = ds.FLAG.data.any(axis=2).compute()
    return uvw, flagged, freqs


# ============================================================================
# MAIN
# ============================================================================

def block_points(uvw, flagged, freqs, chans):
    """Per channel in `chans`, the unflagged (u, v) in wavelengths."""
    for ch in chans:
        ok = ~flagged[:, ch]
        yield ch, uvw[ok, :2] * (freqs[ch] / C_LIGHT)


def hull_vertices(points):
    if len(points) < 3:
        return points
    try:
        return points[ConvexHull(points).vertices]
    except QhullError:
        return points


def block_ellipse(uvw, flagged, freqs, chans):
    """Enclosing ellipse of a channel block's unflagged uv points and their mirrors."""
    verts = [hull_vertices(uv) for _, uv in block_points(uvw, flagged, freqs, chans) if len(uv)]
    if not verts:
        return None, 0
    verts = np.vstack(verts)
    verts = hull_vertices(np.vstack([verts, -verts]))
    msg(f'    {len(verts)} hull vertices')
    return enclosing_ellipse(verts), len(verts)


def main():
    parser = argparse.ArgumentParser(description='Suggest circular -minuv-l and -maxuv-l cuts matching uv coverage across the band')
    parser.add_argument('ms', help='Measurement set')
    parser.add_argument('--field', default=FIELD, help=f'Field name or ID (default: {FIELD}; None uses the first field)')
    parser.add_argument('--edge-frac', type=float, default=EDGE_FRAC,
                        help=f'Fraction of the usable band taken at each end (default: {EDGE_FRAC})')
    parser.add_argument('--outdir', default=OUTDIR, help=f'Directory for the diagnostic plot (default: {OUTDIR})')
    parser.add_argument('--bins', type=int, default=BINS, help=f'uv plot grid size per axis (default: {BINS})')
    args = parser.parse_args()

    ms = args.ms.rstrip('/')
    field_id, field_name = resolve_field(ms, args.field)
    msg(f'MS: {ms}')
    msg(f'Field: {field_name} (ID {field_id})')

    uvw, flagged, freqs = load_field(ms, field_id)
    cross = np.any(uvw[:, :2] != 0, axis=1)
    if not cross.all():
        msg(f'Dropping {np.count_nonzero(~cross):,} autocorrelation row(s)')
        uvw, flagged = uvw[cross], flagged[cross]
    nrow, nchan = flagged.shape
    msg(f'{nrow:,} rows, {nchan} channels, {freqs[0]/1e9:.4f}-{freqs[-1]/1e9:.4f} GHz')

    # Usable band: drop channels with no unflagged data at either end
    good = np.where((~flagged).any(axis=0))[0]
    if len(good) == 0:
        sys.exit('ERROR: every channel is fully flagged for this field')
    c0, c1 = good[0], good[-1]
    n_usable = c1 - c0 + 1
    n_edge = max(1, int(round(args.edge_frac * n_usable)))
    bottom = np.arange(c0, c0 + n_edge)
    top = np.arange(c1 - n_edge + 1, c1 + 1)
    msg(f'Usable band: channels {c0}-{c1} ({freqs[c0]/1e9:.4f}-{freqs[c1]/1e9:.4f} GHz); '
        f'{n_edge} channel(s) at each end')
    msg(f'  Bottom: channels {bottom[0]}-{bottom[-1]} ({freqs[bottom[0]]/1e9:.4f}-{freqs[bottom[-1]]/1e9:.4f} GHz)')
    msg(f'  Top:    channels {top[0]}-{top[-1]} ({freqs[top[0]]/1e9:.4f}-{freqs[top[-1]]/1e9:.4f} GHz)')

    msg('Fitting enclosing ellipses...')
    msg('  Bottom of band:')
    ell_bottom, _ = block_ellipse(uvw, flagged, freqs, bottom)
    msg('  Top of band:')
    ell_top, _ = block_ellipse(uvw, flagged, freqs, top)
    if ell_bottom is None:
        sys.exit('ERROR: could not fit an ellipse to the bottom of the band')

    max_cut = ell_bottom[2]
    # Alternative outer cut: the longest unflagged baseline at the bottom of the band,
    # which matches the maximum baseline length across the band
    max_cut_uvmax = max(np.hypot(uv[:, 0], uv[:, 1]).max()
                        for _, uv in block_points(uvw, flagged, freqs, bottom) if len(uv))
    for label, ell in (('Bottom', ell_bottom), ('Top', ell_top)):
        if ell is None:
            msg(f'  {label}: ellipse fit failed')
            continue
        c, a, b, t = ell
        msg(f'  {label}: semi-major {a/1e3:.2f} klambda, semi-minor {b/1e3:.2f} klambda, '
            f'PA {np.degrees(t):.1f} deg, centre ({c[0]:.1f}, {c[1]:.1f}) lambda')

    # Inner cut: the shortest unflagged baseline at the top of the band
    min_cut = min((np.hypot(uv[:, 0], uv[:, 1]).min()
                   for _, uv in block_points(uvw, flagged, freqs, top) if len(uv)), default=None)
    if min_cut is None:
        sys.exit('ERROR: no unflagged data at the top of the band')
    msg(f'  Top: shortest unflagged baseline {min_cut:.0f} lambda')
    min_cut_calc = min_cut
    if str(cfg.WSC_MINUVL).strip() != '':
        min_cut = float(cfg.WSC_MINUVL)
        msg(f'  Using WSC_MINUVL from config ({min_cut:.0f} lambda) as the inner cut instead')

    # What each cut removes from each end of the band, on top of the existing flags
    msg('')
    msg('=' * 70)
    for label, chans in (('Bottom', bottom), ('Top', top)):
        n_total = nrow * len(chans)
        n_unflagged = n_inner = n_outer = n_outer_uvmax = 0
        for _, uv in block_points(uvw, flagged, freqs, chans):
            r = np.hypot(uv[:, 0], uv[:, 1])
            n_unflagged += len(r)
            n_inner += int(np.count_nonzero(r < min_cut))
            n_outer += int(np.count_nonzero(r > max_cut))
            n_outer_uvmax += int(np.count_nonzero(r > max_cut_uvmax))
        frac_flagged = 1.0 - n_unflagged / n_total
        frac_inner = n_inner / n_unflagged if n_unflagged else 0.0
        frac_outer = n_outer / n_unflagged if n_unflagged else 0.0
        frac_both = frac_inner + frac_outer
        msg(f'{label} {len(chans)} channel(s), {n_total:,} visibilities:')
        msg(f'  Already flagged:                      {frac_flagged:.2%}')
        msg(f'  Unflagged inside -minuv-l:            {n_inner:,} of {n_unflagged:,} ({frac_inner:.2%})')
        msg(f'  Unflagged beyond -maxuv-l:            {n_outer:,} of {n_unflagged:,} ({frac_outer:.2%})')
        msg(f'  Flagged + both cuts (total lost):     {frac_flagged + (1 - frac_flagged) * frac_both:.2%}')
        frac_outer_uvmax = n_outer_uvmax / n_unflagged if n_unflagged else 0.0
        msg(f'  Unflagged beyond bottom uvmax cut:    {n_outer_uvmax:,} of {n_unflagged:,} ({frac_outer_uvmax:.2%})')
        msg('')
    msg(f'Outer cut (bottom-of-band semi-minor axis): {max_cut:.0f} lambda ({max_cut/1e3:.2f} klambda)')
    msg(f'Outer cut at the bottom-of-band uvmax:      {max_cut_uvmax:.0f} lambda ({max_cut_uvmax/1e3:.2f} klambda)')
    msg(f'Inner cut (top-of-band shortest baseline):  {min_cut_calc:.0f} lambda ({min_cut_calc/1e3:.3f} klambda)')
    if min_cut != min_cut_calc:
        msg(f'Inner cut used (WSC_MINUVL from config):    {min_cut:.0f} lambda ({min_cut/1e3:.3f} klambda)')
    msg(f'  -minuv-l {min_cut:.0f} -maxuv-l {max_cut:.0f}')
    msg('=' * 70)

    # Diagnostic plot: unflagged uv coverage (with mirrors) at each end of the band,
    # the full plane on the top row and the inner hole zoomed in on the bottom row
    extent = max(e[1] for e in (ell_bottom, ell_top) if e is not None) * 1.05 / 1e3
    zoom = 4.0 * min_cut / 1e3
    fig, axes = plt.subplots(2, 2, figsize=(16, 16))
    circle = np.linspace(0, 2 * np.pi, 400)
    for row, lim in ((0, extent), (1, zoom)):
        edges = np.linspace(-lim, lim, args.bins + 1)
        for col, (label, chans) in enumerate((('Bottom', bottom), ('Top', top))):
            ax = axes[row, col]
            hist = np.zeros((args.bins, args.bins))
            for _, uv in block_points(uvw, flagged, freqs, chans):
                uv = uv / 1e3
                hist += np.histogram2d(uv[:, 0], uv[:, 1], bins=[edges, edges])[0]
                hist += np.histogram2d(-uv[:, 0], -uv[:, 1], bins=[edges, edges])[0]
            hist[hist == 0] = np.nan
            ax.imshow(hist.T, origin='lower', extent=[-lim, lim, -lim, lim],
                      cmap='Greys', norm=LogNorm(), interpolation='nearest')
            # Kept uv range: the ring between the inner and outer cuts (inner circle traced backwards)
            ring_u = np.concatenate([max_cut * np.cos(circle), min_cut * np.cos(circle[::-1])]) / 1e3
            ring_v = np.concatenate([max_cut * np.sin(circle), min_cut * np.sin(circle[::-1])]) / 1e3
            ax.fill(ring_u, ring_v, color='tab:green', alpha=0.15, lw=0, label='kept uv range')
            for ell, colour, name in ((ell_bottom, 'tab:blue', 'bottom ellipse'), (ell_top, 'tab:red', 'top ellipse')):
                if ell is not None:
                    outline = ellipse_outline(ell) / 1e3
                    ax.plot(outline[:, 0], outline[:, 1], color=colour, lw=1.5, label=name)
            ax.plot(max_cut / 1e3 * np.cos(circle), max_cut / 1e3 * np.sin(circle), color='tab:green',
                    lw=1.5, ls='--', label=f'maxuv {max_cut/1e3:.2f} k$\\lambda$')
            ax.plot(max_cut_uvmax / 1e3 * np.cos(circle), max_cut_uvmax / 1e3 * np.sin(circle), color='grey',
                    lw=1.5, ls='--', label=f'bottom uvmax {max_cut_uvmax/1e3:.2f} k$\\lambda$')
            ax.plot(min_cut / 1e3 * np.cos(circle), min_cut / 1e3 * np.sin(circle), color='tab:orange',
                    lw=1.5, ls='--', label=f'minuv {min_cut/1e3:.3f} k$\\lambda$'
                    + (' (config)' if min_cut != min_cut_calc else ''))
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_title(f'{label} {len(chans)} channel(s): '
                         f'{freqs[chans[0]]/1e9:.3f}-{freqs[chans[-1]]/1e9:.3f} GHz'
                         + (' (inner zoom)' if row == 1 else ''))
            ax.set_xlabel('u (k$\\lambda$)')
            ax.set_ylabel('v (k$\\lambda$)')
            ax.set_aspect('equal')
            ax.legend(loc='upper right')
    fig.suptitle(f'{os.path.basename(ms)} -- {field_name}: unflagged uv coverage')
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(args.outdir, exist_ok=True)
    png = os.path.join(args.outdir, f'uv_band_cut_{os.path.basename(ms)}_{field_name}.png')
    fig.savefig(png, dpi=150)
    plt.close(fig)
    msg(f'Plot: {png}')


if __name__ == '__main__':
    main()
