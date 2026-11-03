"""
antenna_delay_snr.py

Per-antenna delay-solve S/N ranking for reference antenna selection, on the
primary calibrator.

For each baseline (a1 < a2) and scan, the unflagged parallel-hand
correlations (XX/YY or RR/LL) are vector-averaged over time into one
complex spectrum

    V_ij(nu) ~ A exp(+2 pi i nu tau_ij)

which is Hanning-tapered, zero-padded and FFT'd (numpy sign convention, so
a positive phase slope gives a positive delay). The delay is the
power-weighted centroid of the central peak, and the S/N is the peak
amplitude over the rms of a noise annulus around it. Each antenna's S/N is
the mean over every connected baseline in every scan.

Run with exec() inside CASA, which provides `tb`, `myms`, `bpcal_name` and
`GAINPLOTS`. Sets `refant_snr_ranking`, a list of (antenna name, S/N)
sorted best first.

Assumes a single spectral window with uniform channel spacing (MeerKAT).
"""

import functools
import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

print = functools.partial(print, flush=True)

ZERO_PAD_FACTOR = 8       # FFT oversampling factor
NOISE_ANNULUS_FWHM = 100  # noise annulus extent on each side of the peak, in FWHM
ZOOM_FWHM = 5             # half-width of the zoom inset in the FFT plots, in FWHM
MIN_UNFLAGGED_FRAC = 0.1  # skip baselines with less unflagged data than this
TOP_N = 7                 # antennas printed at each end of the ranking

PARALLEL_HANDS = {5: "RR", 8: "LL", 9: "XX", 12: "YY"}  # casacore Stokes codes


# --------------------------------------------------------------------------
# MS access
# --------------------------------------------------------------------------

def read_col(table, col):
    """
    Read one column of a table.

    Parameters
    ----------
    table : str
        Table path.
    col : str
        Column name.

    Returns
    -------
    numpy.ndarray
        Row axis first, as python-casacore returns it.
    """
    tb.open(table)
    values = tb.getcol(col).T
    tb.close()
    return values


def query_cols(table, where, cols):
    """
    Read columns of the rows of a table matching a TaQL condition.

    Parameters
    ----------
    table : str
        Table path.
    where : str
        TaQL condition.
    cols : list of str
        Column names.

    Returns
    -------
    list of numpy.ndarray
        One array per column, row axis first (empty when no row matches).
    """
    tb.open(table)
    sub = tb.query(where)
    values = [sub.getcol(c).T if sub.nrows() else np.array([]) for c in cols]
    sub.close()
    tb.close()
    return values


def field_id(ms, name):
    """
    ID of the field called `name`.

    Raises
    ------
    RuntimeError
        If the name is missing or appears more than once.
    """
    ids = [i for i, n in enumerate(read_col(f"{ms}/FIELD", "NAME")) if n == name]
    if len(ids) != 1:
        raise RuntimeError(f"Field '{name}' found {len(ids)} times in {ms}/FIELD")
    return ids[0]


def parallel_hands(ms):
    """
    Indices and labels of the parallel-hand correlations.

    Returns
    -------
    idx : list of int
    labels : list of str
    """
    corr_type = read_col(f"{ms}/POLARIZATION", "CORR_TYPE")[0]
    idx = [i for i, c in enumerate(corr_type) if c in PARALLEL_HANDS]
    if not idx:
        raise RuntimeError(f"No parallel-hand correlations in CORR_TYPE={corr_type}")
    return idx, [PARALLEL_HANDS[corr_type[i]] for i in idx]


def baseline_spectrum(ms, a1, a2, corr_idx, where):
    """
    Time- and correlation-averaged spectrum of baseline (a1, a2), a1 < a2.

    Rows stored as (a2, a1) are conjugated before averaging.

    Parameters
    ----------
    ms : str
        Measurement set.
    a1, a2 : int
        Antenna indices, a1 < a2.
    corr_idx : list of int
        Correlation indices to average.
    where : str
        Extra TaQL selection (field and scan).

    Returns
    -------
    spectrum : numpy.ndarray of complex, shape (nchan,), or None
        None when the baseline has no rows.
    unflagged_frac : float
        Fraction of selected samples that are unflagged.
    """
    pair = f"((ANTENNA1=={a1} && ANTENNA2=={a2}) || (ANTENNA1=={a2} && ANTENNA2=={a1}))"
    data, flag, flag_row, ant1 = query_cols(ms, f"{pair} && {where}",
                                            ["DATA", "FLAG", "FLAG_ROW", "ANTENNA1"])
    if not data.size:
        return None, 0.0

    data[ant1 == a2] = np.conj(data[ant1 == a2])
    good = ~(flag[:, :, corr_idx] | flag_row[:, None, None])
    data = np.where(good, data[:, :, corr_idx], 0.0)

    n = good.sum(axis=(0, 2))
    spectrum = np.divide(data.sum(axis=(0, 2)), n, out=np.zeros(n.shape, complex), where=n > 0)
    return spectrum, float(good.mean())


# --------------------------------------------------------------------------
# Delay and S/N
# --------------------------------------------------------------------------

def delay_spectrum(spectrum, df):
    """
    Oversampled delay spectrum of a channel spectrum.

    Parameters
    ----------
    spectrum : numpy.ndarray of complex, shape (nchan,)
    df : float
        Channel width (Hz).

    Returns
    -------
    delay_ns : numpy.ndarray
        Delay axis (ns), increasing.
    amp : numpy.ndarray
        |FFT| of the Hanning-tapered, zero-padded spectrum.
    """
    npad = len(spectrum) * ZERO_PAD_FACTOR
    amp = np.abs(np.fft.fftshift(np.fft.fft(spectrum * np.hanning(len(spectrum)), npad)))
    delay_ns = np.fft.fftshift(np.fft.fftfreq(npad, d=df)) * 1e9
    return delay_ns, amp


def peak_widths(nchan):
    """
    FWHM and first-null half-width of the noise-free delay peak.

    Measured on the delay spectrum of a flat (zero-delay) spectrum, so they
    are exact for this taper and padding.

    Parameters
    ----------
    nchan : int

    Returns
    -------
    fwhm : float
        Full width at half maximum (oversampled bins).
    null : int
        Bins from the peak to its first null.
    """
    amp = delay_spectrum(np.ones(nchan, complex), 1.0)[1]
    side = amp[int(np.argmax(amp)):]
    null = int(np.argmax(np.diff(side) >= 0))
    i = int(np.argmax(side < 0.5 * side[0]))
    fwhm = 2.0 * (i - 1 + (side[i - 1] - 0.5 * side[0]) / (side[i - 1] - side[i]))
    return fwhm, null


def delay_and_snr(delay_ns, amp, null, noise_half_width):
    """
    Delay and S/N of the strongest peak in a delay spectrum.

    Parameters
    ----------
    delay_ns, amp : numpy.ndarray
        Output of `delay_spectrum`.
    null : int
        Peak half-width out to its first null (bins).
    noise_half_width : int
        Half-width of the noise annulus (bins), measured from the peak.

    Returns
    -------
    dict
        delay_ns (power-weighted centroid of the peak), snr (peak over
        annulus rms), noise_rms, peak (index), window and noise_window
        (index ranges of the peak and the annulus).
    """
    n = len(amp)
    peak = int(np.argmax(amp))
    lo, hi = max(0, peak - null), min(n, peak + null + 1)
    w = amp[lo:hi] ** 2
    delay = float(np.sum(delay_ns[lo:hi] * w) / np.sum(w))

    nlo, nhi = max(0, peak - noise_half_width), min(n, peak + noise_half_width + 1)
    noise = np.zeros(n, bool)
    noise[nlo:nhi] = True
    noise[lo:hi] = False
    if noise.sum() < 5:  # peak near an edge: use everything outside the peak
        noise[:], noise[lo:hi] = True, False
        nlo, nhi = 0, n
    rms = float(np.std(amp[noise]))

    return {"delay_ns": delay, "snr": float(amp[peak]) / rms if rms > 0 else np.inf,
            "noise_rms": rms, "peak": peak, "window": (lo, hi), "noise_window": (nlo, nhi)}


def solve_scan(ms, fid, scan, corr_idx, df, null, noise_half_width):
    """
    Delay and S/N of every baseline of one field and scan.

    Parameters
    ----------
    ms : str
    fid, scan : int
        Field ID and scan number.
    corr_idx : list of int
        Parallel-hand correlation indices.
    df : float
        Channel width (Hz).
    null, noise_half_width : int
        See `delay_and_snr`.

    Returns
    -------
    list of dict
        One `delay_and_snr` result per usable baseline, with a1, a2 and the
        delay spectrum (delay_axis, amp) added.
    """
    where = f"FIELD_ID=={fid} && SCAN_NUMBER=={scan}"
    a1s, a2s = query_cols(ms, where, ["ANTENNA1", "ANTENNA2"])
    pairs = sorted({(min(i, j), max(i, j)) for i, j in zip(a1s, a2s) if i != j})

    results = []
    for a1, a2 in pairs:
        spectrum, frac = baseline_spectrum(ms, a1, a2, corr_idx, where)
        if spectrum is None or frac < MIN_UNFLAGGED_FRAC:
            continue
        delay_axis, amp = delay_spectrum(spectrum, df)
        results.append({**delay_and_snr(delay_axis, amp, null, noise_half_width),
                        "a1": a1, "a2": a2, "delay_axis": delay_axis, "amp": amp})
    return results


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------

def median_ci(values):
    """Median and its 68% interval (1.4826 x MAD) as (median, lo, hi)."""
    med = float(np.median(values))
    sig = 1.4826 * float(np.median(np.abs(np.asarray(values) - med)))
    return med, med - sig, med + sig


def antenna_axes(n):
    """Figure and axes sized for n antenna labels, at most 120 of them shown."""
    fig, ax = plt.subplots(figsize=(min(60.0, max(8.0, 0.35 * n)), 5.0))
    ticks = list(range(0, n, max(1, -(-n // 120))))
    return fig, ax, ticks


def plot_snr(snrs_by_scan, names, out_png, title):
    """
    Median S/N +/- 68% interval per antenna, one series per scan.

    Parameters
    ----------
    snrs_by_scan : dict of {label: dict of {antenna: list of float}}
    names : list of str
        Antenna names.
    out_png : pathlib.Path
    title : str
    """
    ants = sorted({a for s in snrs_by_scan.values() for a in s}, key=lambda a: names[a])
    fig, ax, ticks = antenna_axes(len(ants))
    x = np.arange(len(ants))
    cmap = plt.get_cmap("tab10" if len(snrs_by_scan) <= 10 else "tab20")
    offsets = np.linspace(-0.3, 0.3, len(snrs_by_scan)) if len(snrs_by_scan) > 1 else [0.0]

    for k, (label, snrs) in enumerate(snrs_by_scan.items()):
        stats = np.array([median_ci(snrs[a]) if snrs.get(a) else (np.nan,) * 3 for a in ants])
        ax.errorbar(x + offsets[k], stats[:, 0],
                    yerr=[stats[:, 0] - stats[:, 1], stats[:, 2] - stats[:, 0]],
                    fmt="o", ms=4, color=cmap(k % cmap.N), capsize=2, alpha=0.85, label=label)

    ax.grid(axis="y", color="0.85", lw=0.6, zorder=0)
    ax.set_xticks(ticks)
    ax.set_xticklabels([names[ants[i]] for i in ticks], rotation=90, fontsize=8)
    ax.set_ylabel("Delay-solve S/N (median, 68% interval)")
    ax.set_title(title)
    if len(snrs_by_scan) > 1:
        ax.legend(fontsize=7, ncol=min(len(snrs_by_scan), 6))
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_best_fft(best, names, out_dir, zoom):
    """
    Delay spectrum of each antenna's highest-S/N baseline.

    Parameters
    ----------
    best : dict of {antenna: dict}
        `solve_scan` result per antenna, with its baseline partner added.
    names : list of str
        Antenna names.
    out_dir : pathlib.Path
    zoom : int
        Half-width of the zoom inset (bins).
    """
    out_dir.mkdir(exist_ok=True)
    for ant, r in best.items():
        d, amp = r["delay_axis"], r["amp"]
        (lo, hi), (nlo, nhi) = r["window"], r["noise_window"]
        zlo, zhi = max(0, r["peak"] - zoom), min(len(d), r["peak"] + zoom + 1)

        fig, ax = plt.subplots(figsize=(7, 4.5))
        axins = ax.inset_axes([0.04, 0.54, 0.42, 0.42])
        for a in (ax, axins):
            a.plot(d, amp, color="0.4", lw=1)
            a.axvspan(d[nlo], d[nhi - 1], color="steelblue", alpha=0.12)
            a.axvspan(d[lo], d[hi - 1], color="orange", alpha=0.3)
            a.axvline(r["delay_ns"], color="crimson", lw=1.5)
            a.axhline(r["noise_rms"], color="steelblue", ls="--", lw=1)
        axins.set_xlim(d[zlo], d[zhi - 1])
        axins.tick_params(labelsize=6)
        ax.indicate_inset_zoom(axins, edgecolor="0.3")

        partner = names[r["partner"]]
        ax.set_xlabel("Delay (ns)")
        ax.set_ylabel("Amplitude")
        ax.set_title(f"{names[ant]}-{partner}: delay {r['delay_ns']:.4f} ns, S/N {r['snr']:.1f}")
        fig.tight_layout()
        fig.savefig(out_dir / f"fft_{names[ant]}_best_{partner}.png", dpi=150)
        plt.close(fig)


# --------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------

def rank_antennas(ms, field, plot_dir):
    """
    Rank antennas by mean delay-solve S/N on one field, and plot it.

    Parameters
    ----------
    ms : str
        Measurement set.
    field : str
        Field name (normally the primary calibrator).
    plot_dir : pathlib.Path
        Directory for refant_snr_ranking.json and the plots.

    Returns
    -------
    list of (str, float)
        Antenna name and mean S/N over every connected baseline in every
        scan, best first.
    """
    plot_dir.mkdir(parents=True, exist_ok=True)

    names = list(read_col(f"{ms}/ANTENNA", "NAME"))
    fid = field_id(ms, field)
    scans = sorted(set(query_cols(ms, f"FIELD_ID=={fid}", ["SCAN_NUMBER"])[0].tolist()))
    corr_idx, corr_labels = parallel_hands(ms)
    chan_freq = read_col(f"{ms}/SPECTRAL_WINDOW", "CHAN_FREQ")[0]
    df = abs(float(np.median(np.diff(chan_freq))))

    fwhm, null = peak_widths(len(chan_freq))
    noise_half_width = null + int(round(NOISE_ANNULUS_FWHM * fwhm))
    print(f"Delay S/N on {field} (field {fid}), scans {scans}, correlations {corr_labels}, "
          f"{len(chan_freq)} channels; peak FWHM {fwhm:.2f} bins, first null {null} bins")

    snrs_by_scan, best = {}, {}
    for scan in scans:
        snrs = {}
        for r in solve_scan(ms, fid, scan, corr_idx, df, null, noise_half_width):
            for ant, partner in ((r["a1"], r["a2"]), (r["a2"], r["a1"])):
                snrs.setdefault(ant, []).append(r["snr"])
                if ant not in best or r["snr"] > best[ant]["snr"]:
                    best[ant] = {**r, "partner": partner}
        print(f"  scan {scan}: {sum(map(len, snrs.values())) // 2} usable baseline(s)")
        if snrs:
            snrs_by_scan[f"scan {scan}"] = snrs

    if not snrs_by_scan:
        raise RuntimeError(f"No usable baselines on {field} in any scan")

    combined = {}
    for snrs in snrs_by_scan.values():
        for ant, values in snrs.items():
            combined.setdefault(ant, []).extend(values)

    ranking = sorted(((names[a], float(np.mean(v))) for a, v in combined.items()),
                     key=lambda r: r[1], reverse=True)

    n = min(TOP_N, len(ranking))
    print(f"Best {n}:  " + ", ".join(f"{a} {s:.1f}" for a, s in ranking[:n]))
    print(f"Worst {n}: " + ", ".join(f"{a} {s:.1f}" for a, s in ranking[::-1][:n]))

    with open(plot_dir / "refant_snr_ranking.json", "w") as f:
        json.dump([{"antenna": a, "avg_snr": s} for a, s in ranking], f, indent=2)

    plot_snr(snrs_by_scan, names, plot_dir / "refant_snr_per_scan.png",
             f"Per-scan delay-solve S/N by antenna ({field})")
    plot_snr({"all scans": combined}, names, plot_dir / "refant_snr_combined.png",
             f"Combined (all-scan) delay-solve S/N by antenna ({field})")
    plot_best_fft(best, names, plot_dir / "refant_fft_per_antenna", int(round(ZOOM_FWHM * fwhm)))
    print(f"Ranking and plots written to {plot_dir}")
    return ranking


refant_snr_ranking = rank_antennas(myms, bpcal_name, Path(GAINPLOTS))
