"""
fix_delay_table.py

Check a K (delay) table after gaincal and, if CASA flagged most of one
field's solutions, replace them with delays from the delay spectra of the
baselines to the leading reference antenna.

For the field just solved:

1. The gaincal's prior tables (same gaintable, gainfield and interp) are
   applied to CORRECTED_DATA, so the data match what gaincal solved on.
   Without prior tables DATA is used.
2. For each scan, antenna and parallel hand, data / MODEL_DATA on the
   antenna's baseline to the refant is vector-averaged over the scan, with
   every row oriented (antenna, refant). Its delay spectrum
   (antenna_delay_snr.py) gives the delay and S/N. With CASA's convention
   V_ij ~ g_i conj(g_j), g_i = exp(2 pi i nu tau_i), the measured delay is
   tau_antenna - tau_refant, which is the K solution.
3. If CAL_1GC_KFIX_DIAGNOSTIC, these delays are plotted against CASA's in
   GAINPLOTS. On a good CASA solve they agree; a sign error shows as
   mirrored delays.
4. If at least CAL_1GC_KFIX_FLAG_FRAC of the field's solutions are flagged,
   the table is copied to <ktab>_bad_<field> and the field's rows are
   rewritten: the new delay where its S/N is at least CAL_1GC_KFIX_MINSNR,
   flagged otherwise.

Load with exec() inside CASA (needs `tb`, `applycal`, TOOLS, GAINPLOTS and
the CAL_1GC_KFIX_* settings), then call fix_delay_table() after each K
gaincal.
"""

import os
import re
import shutil

exec(open(f"{TOOLS}/antenna_delay_snr.py").read())


def uvrange_mask(uvw, uvrange):
    """
    Rows inside a CASA uvrange given in metres.

    Parameters
    ----------
    uvw : numpy.ndarray, shape (nrow, 3)
        UVW (m).
    uvrange : str
        '', '>X m', '<X m' or 'X~Y m'.

    Returns
    -------
    numpy.ndarray of bool, shape (nrow,)
    """
    if not uvrange:
        return np.ones(len(uvw), bool)
    # Only plain metre ranges are understood; anything else stops the run
    # rather than silently solving on different baselines from gaincal
    m = re.fullmatch(r"\s*([<>])?\s*([\d.]+)\s*(?:~\s*([\d.]+))?\s*m\s*", uvrange)
    if m is None or (m[1] is None) == (m[3] is None):
        raise ValueError(f"uvrange {uvrange!r} is not supported (metres only: '>X', '<X' or 'X~Y')")
    uv = np.hypot(uvw[:, 0], uvw[:, 1])
    if m[3] is not None:
        return (uv >= float(m[2])) & (uv <= float(m[3]))
    return uv > float(m[2]) if m[1] == ">" else uv < float(m[2])


def refant_spectra(ms, ant, ref, where, datacol, corr_idx, nchan, uvrange):
    """
    Scan-averaged data / MODEL_DATA spectrum of baseline (ant, ref), per correlation.

    Rows stored as (ref, ant) are conjugated, so every spectrum is oriented
    (ant, ref).

    Parameters
    ----------
    ms : str
    ant, ref : int
        Antenna and reference antenna indices.
    where : str
        TaQL selection (field and scan).
    datacol : str
        'DATA' or 'CORRECTED_DATA'.
    corr_idx : list of int
        Parallel-hand correlation indices.
    nchan : int
    uvrange : str
        See `uvrange_mask`.

    Returns
    -------
    spectra : numpy.ndarray of complex, shape (ncorr, nchan)
    frac : numpy.ndarray, shape (ncorr,)
        Unflagged fraction per correlation.
    """
    pair = f"((ANTENNA1=={ant} && ANTENNA2=={ref}) || (ANTENNA1=={ref} && ANTENNA2=={ant}))"
    cols = query_cols(ms, f"{pair} && {where}",
                      [datacol, "MODEL_DATA", "FLAG", "FLAG_ROW", "ANTENNA1", "UVW"])
    # Same baselines as gaincal: drop rows outside its uvrange
    if cols[0].size:
        keep = uvrange_mask(cols[5], uvrange)
        cols = [c[keep] for c in cols]
    data, model, flag, flag_row, ant1, _ = cols
    if not data.size:
        return np.zeros((len(corr_idx), nchan), complex), np.zeros(len(corr_idx))

    # Divide by the model as gaincal does, leaving only the antenna-based
    # phase; samples with no model are treated as flagged
    vis = np.divide(data, model, out=np.zeros_like(data), where=model != 0)
    # A (ref, ant) row is the conjugate of (ant, ref): flip it so the delay
    # always comes out as tau_antenna - tau_refant
    vis[ant1 == ref] = np.conj(vis[ant1 == ref])
    good = ~(flag | flag_row[:, None, None] | (model == 0))[:, :, corr_idx]
    vis = np.where(good, vis[:, :, corr_idx], 0.0)

    # Vector average over the scan, per channel and correlation
    n = good.sum(axis=0)
    spectra = np.divide(vis.sum(axis=0), n, out=np.zeros(n.shape, complex), where=n > 0)
    return spectra.T, good.mean(axis=(0, 1))


def refant_delays(ms, fid, scans, ants, ref, datacol, uvrange):
    """
    Delay and S/N of every antenna against the refant, per scan and parallel hand.

    Parameters
    ----------
    ms : str
    fid : int
        Field ID.
    scans, ants : list of int
        Scans and antennas to solve.
    ref : int
        Reference antenna (delay 0).
    datacol, uvrange : str
        See `refant_spectra`.

    Returns
    -------
    dict of {(scan, antenna): (delay_ns, snr)}
        numpy.ndarray pairs of shape (ncorr,). Unusable baselines have NaN
        delay and S/N 0; the refant has delay 0 and infinite S/N.
    """
    corr_idx = parallel_hands(ms)[0]
    chan_freq = read_col(f"{ms}/SPECTRAL_WINDOW", "CHAN_FREQ")[0]
    df = abs(float(np.median(np.diff(chan_freq))))
    fwhm, null = peak_widths(len(chan_freq))
    noise_half_width = null + int(round(NOISE_ANNULUS_FWHM * fwhm))

    # One solution per scan, antenna and polarization, like gaincal with
    # solint='inf' and no combine
    out = {}
    for scan in scans:
        where = f"FIELD_ID=={fid} && SCAN_NUMBER=={scan}"
        for ant in ants:
            # The refant defines zero delay
            if ant == ref:
                out[scan, ant] = (np.zeros(len(corr_idx)), np.full(len(corr_idx), np.inf))
                continue
            spectra, frac = refant_spectra(ms, ant, ref, where, datacol, corr_idx,
                                           len(chan_freq), uvrange)
            delay, snr = np.full(len(corr_idx), np.nan), np.zeros(len(corr_idx))
            for k, (spectrum, f) in enumerate(zip(spectra, frac)):
                # Too little unflagged data on the baseline to the refant:
                # left at S/N 0, so it ends up flagged
                if f >= MIN_UNFLAGGED_FRAC:
                    r = delay_and_snr(*delay_spectrum(spectrum, df), null, noise_half_width)
                    delay[k], snr[k] = r["delay_ns"], r["snr"]
            out[scan, ant] = (delay, snr)
    return out


def plot_delay_comparison(out_png, names, ants, casa_delay, casa_flag, delay, good, title):
    """
    Re-solved delays against CASA's, one column per parallel hand.

    Top: both sets of delays per antenna (CASA's unflagged solutions only).
    Bottom: re-solved minus CASA where both are good.

    Parameters
    ----------
    out_png : str
    names : list of str
        Antenna names.
    ants : numpy.ndarray of int, shape (nrow,)
        Antenna of each solution (one row per scan and antenna).
    casa_delay, delay : numpy.ndarray, shape (nrow, npol)
        CASA's and the re-solved delays (ns).
    casa_flag, good : numpy.ndarray of bool, shape (nrow, npol)
        CASA's flags, and which re-solved delays pass the S/N cut.
    title : str
    """
    order = sorted(set(ants), key=lambda a: names[a])
    x = np.array([order.index(a) for a in ants])
    npol = delay.shape[1]
    fig, axes = plt.subplots(2, npol, figsize=(max(8.0, min(40.0, 0.2 * len(order) * npol)), 7),
                             sharex=True, squeeze=False)

    for p in range(npol):
        casa_ok = ~casa_flag[:, p]
        both = casa_ok & good[:, p]
        top, bottom = axes[0, p], axes[1, p]
        top.plot(x[casa_ok], casa_delay[casa_ok, p], "o", mfc="none", mec="0.3", ms=6,
                 label=f"CASA ({casa_ok.sum()} unflagged)")
        top.plot(x[good[:, p]], delay[good[:, p], p], ".", color="crimson", ms=5,
                 label=f"re-solved ({good[:, p].sum()} good)")
        top.set_title(f"Sol{p + 1}")
        top.set_ylabel("Delay (ns)")
        top.legend(fontsize=7)

        diff = delay[both, p] - casa_delay[both, p]
        bottom.axhline(0.0, color="0.6", lw=1)
        bottom.plot(x[both], diff, ".", color="k", ms=5)
        if both.any():
            bottom.set_title(f"median |diff| = {np.median(np.abs(diff)):.3f} ns ({both.sum()} solutions)",
                             fontsize=8)
        bottom.set_ylabel("Re-solved - CASA (ns)")
        bottom.set_xticks(range(len(order)))
        bottom.set_xticklabels([names[a] for a in order], rotation=90, fontsize=6)
        for ax in (top, bottom):
            ax.grid(color="0.9", lw=0.6)

    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def fix_delay_table(ms, ktab, field, refant, gaintable=(), gainfield=(), interp=(), uvrange=""):
    """
    Check one field's solutions in a K table, and re-solve and rewrite them if they failed.

    Parameters
    ----------
    ms : str
    ktab : str
        K table just written by gaincal.
    field : str
        Field name solved.
    refant : str or int
        gaincal's refant; only the first entry (index or name) is used.
    gaintable, gainfield, interp : sequence of str
        gaincal's prior tables and how they were applied.
    uvrange : str
        gaincal's uvrange (metres).
    """
    # Leading refant of gaincal's list, given as an index or a name
    names = list(read_col(f"{ms}/ANTENNA", "NAME"))
    ref = str(refant).split(",")[0].strip()
    ref = int(ref) if ref.isdigit() else names.index(ref)
    fid = field_id(ms, field)

    # The table is appended across fields: only this field's rows are
    # checked and, if needed, rewritten
    rows = np.flatnonzero(read_col(ktab, "FIELD_ID") == fid)
    if not rows.size:
        print(f"fix_delay_table: no solutions for {field} in {ktab}")
        return
    ants = read_col(ktab, "ANTENNA1")[rows]
    scans = read_col(ktab, "SCAN_NUMBER")[rows]
    casa_delay = read_col(ktab, "FPARAM")[rows, 0, :]
    casa_flag = read_col(ktab, "FLAG")[rows, 0, :]
    casa_ref = sorted({names[a] for a in read_col(ktab, "ANTENNA2")[rows]})

    # Fraction of CASA's solutions (all scans, antennas, polarizations) flagged
    flagged = float(casa_flag.mean())
    rewrite = flagged >= CAL_1GC_KFIX_FLAG_FRAC
    print(f"fix_delay_table: {os.path.basename(ktab)} {field}: {100 * flagged:.0f}% of solutions flagged"
          + (f" (>= {100 * CAL_1GC_KFIX_FLAG_FRAC:.0f}%), re-solving against {names[ref]}" if rewrite else ""))
    if not (rewrite or CAL_1GC_KFIX_DIAGNOSTIC):
        return

    # Calibrate with the same prior tables gaincal applied, so the delays are
    # measured on the same data. applycal also applies those tables' flags
    # to this field in the MS, as the later applycals with them do.
    datacol = "DATA"
    if gaintable:
        applycal(vis=ms, field=field, gaintable=list(gaintable), gainfield=list(gainfield),
                 interp=list(interp), parang=False, calwt=False, flagbackup=False)
        datacol = "CORRECTED_DATA"

    # Re-solve, in the same row order as the table
    sol = refant_delays(ms, fid, sorted(set(scans)), sorted(set(ants)), ref, datacol, uvrange)
    delay = np.array([sol[s, a][0] for s, a in zip(scans, ants)])
    snr = np.array([sol[s, a][1] for s, a in zip(scans, ants)])
    good = snr >= CAL_1GC_KFIX_MINSNR

    # Re-solved vs CASA, made for every check so good CASA solves validate the method
    out_png = f"{GAINPLOTS}/kfix_{os.path.basename(ktab)}_{field}.png"
    plot_delay_comparison(out_png, names, ants, casa_delay, casa_flag, delay, good,
                          f"{os.path.basename(ktab)} {field}: re-solved (refant {names[ref]}) "
                          f"vs CASA (refant {', '.join(casa_ref)}), {100 * flagged:.0f}% CASA flagged")
    print(f"fix_delay_table: {good.sum()} of {good.size} re-solved delays pass S/N >= "
          f"{CAL_1GC_KFIX_MINSNR}; comparison plot {out_png}")
    if not rewrite:
        return

    # Keep CASA's table before touching it
    backup = f"{ktab}_bad_{field}"
    shutil.rmtree(backup, ignore_errors=True)
    shutil.copytree(ktab, backup)

    # Replace every solution of this field, so they all share one refant:
    # the re-solved delay (ns) where it passes the S/N cut, flagged otherwise.
    # The refant keeps CASA's SNR value; ANTENNA2 records the refant used.
    tb.open(ktab, nomodify=False)
    fparam, flag, snr_col, ant2 = (tb.getcol(c).T for c in ("FPARAM", "FLAG", "SNR", "ANTENNA2"))
    fparam[rows, 0, :] = np.where(good, delay, 0.0)
    flag[rows, 0, :] = ~good
    not_ref = rows[ants != ref]
    snr_col[not_ref, 0, :] = snr[ants != ref]
    ant2[rows] = ref
    for col, values in (("FPARAM", fparam), ("FLAG", flag), ("SNR", snr_col), ("ANTENNA2", ant2)):
        tb.putcol(col, values.T)
    tb.close()
    print(f"fix_delay_table: rewrote {field} in {ktab}; original kept as {backup}")
