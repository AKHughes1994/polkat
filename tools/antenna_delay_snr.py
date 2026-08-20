#!/usr/bin/env python3
"""
antenna_delay_snr.py
=====================
Per-antenna delay-solve S/N ranking for reference-antenna selection, based
on a primary/delay-calibrator field.

Method
------
For each baseline (a1, a2), the parallel-hand correlations (XX/YY, or
RR/LL if the array is circular) are vector-averaged over time and over
all unflagged rows within one scan to build one complex spectrum per
baseline:

    V_ij(f) ~= A * exp(-2*pi*j*f*tau_ij) + noise

That spectrum is Hanning-tapered, zero-padded, and FFT'd into an
oversampled delay spectrum. Delay is the power-weighted centroid of
the central peak (sub-bin precision). The peak's FWHM and first-null
half-width are measured analytically (`analytic_peak_widths`) from the
window function's own noise-free peak shape. S/N is the peak amplitude
over the RMS of a local noise annulus that starts just outside the
peak (its first null) and extends NOISE_ANNULUS_FWHM FWHM further on
each side:

    SNR = peak_amplitude / local_annulus_rms

Baselines are grouped by antenna. For each antenna, the single baseline
with the highest SNR (among all baselines touching that antenna) is
taken as that antenna's representative measurement for that solve.

Reference-antenna ranking
--------------------------
This solve is repeated **independently per scan** on the given field
(e.g. the primary), giving each antenna one representative SNR value
per scan. Those per-scan values are then pooled across all scans
("combined") and averaged to give a single per-antenna S/N used to
rank candidate reference antennas -- an antenna whose delay solution
is only strong in one scan won't outrank one that's consistently good
throughout the track.

Assumptions
-----------
- Single spectral window (SPECTRAL_WINDOW row 0 / uniform channel
  spacing); like MeerKAT. Multi-SPW MSs are not handled.
"""

from __future__ import annotations

import functools
import json
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

print = functools.partial(print, flush=True)

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

ZERO_PAD_FACTOR = 8       # oversampling factor applied before the FFT
NOISE_ANNULUS_FWHM = 100  # noise annulus extent on EACH side of the excluded peak, in FWHM
ZOOM_FWHM = 5              # per-antenna FFT plot zoom half-width, in FWHM
MIN_UNFLAGGED_FRAC = 0.1  # skip a baseline if less than this fraction of samples are unflagged
TOP_N = 7                 # number of best/worst antennas to print/rank

# Debug/smoke-test switch: when True, stop after the first baseline that
# clears MIN_UNFLAGGED_FRAC and run the full pipeline (delay/SNR estimate +
# all plots) on just that one baseline as a quick example/sanity check.
SINGLE_BASELINE_TEST = False

# casacore Stokes enum codes for the parallel-hand correlations
PARALLEL_HAND_NAMES = {
    5: "RR", 8: "LL",    # circular
    9: "XX", 12: "YY",   # linear
}

# --------------------------------------------------------------------------
# MS access helpers
# --------------------------------------------------------------------------

def get_antenna_names(ms_path: Path) -> list[str]:
    tb.open(str(ms_path / "ANTENNA"))
    names = list(tb.getcol("NAME"))
    tb.close()
    return names


def get_field_id(ms_path: Path, field_name: str) -> int:
    tb.open(str(ms_path / "FIELD"))
    names = list(tb.getcol("NAME"))
    tb.close()
    matches = [i for i, n in enumerate(names) if n == field_name]
    if not matches:
        raise RuntimeError(f"Field '{field_name}' not found in {ms_path}::FIELD (have {names})")
    if len(matches) > 1:
        raise RuntimeError(f"Field '{field_name}' is ambiguous in {ms_path}::FIELD (IDs {matches})")
    return matches[0]


def get_scans_for_field(ms_path: Path, field_id: int) -> list[int]:
    tb.open(str(ms_path))
    sub = tb.query(f"FIELD_ID=={field_id}")
    scans = sorted({int(s) for s in sub.getcol("SCAN_NUMBER")})
    sub.close()
    tb.close()
    return scans


def get_parallel_hand_indices(ms_path: Path) -> tuple[list[int], list[str]]:
    tb.open(str(ms_path / "POLARIZATION"))
    corr_type = tb.getcol("CORR_TYPE").T[0]  # tb returns (ncorr, nrow); flip to (nrow, ncorr)
    tb.close()
    idx = [i for i, c in enumerate(corr_type) if c in PARALLEL_HAND_NAMES]
    if not idx:
        raise RuntimeError(
            f"No parallel-hand (XX/YY or RR/LL) correlations found in CORR_TYPE={corr_type}"
        )
    labels = [PARALLEL_HAND_NAMES[corr_type[i]] for i in idx]
    return idx, labels


def get_spw_chan_freq(ms_path: Path, spw_id: int = 0) -> np.ndarray:
    tb.open(str(ms_path / "SPECTRAL_WINDOW"))
    chan_freq = tb.getcol("CHAN_FREQ").T[spw_id]  # Hz; tb returns (nchan, nrow)
    tb.close()
    return np.asarray(chan_freq, dtype=float)


def get_present_baselines(
    ms_path: Path, field_id: int | None = None, scan_number: int | None = None
) -> list[tuple[int, int]]:
    tb.open(str(ms_path))
    conds = []
    if field_id is not None:
        conds.append(f"FIELD_ID=={field_id}")
    if scan_number is not None:
        conds.append(f"SCAN_NUMBER=={scan_number}")
    if conds:
        sub = tb.query(" && ".join(conds))
        a1 = sub.getcol("ANTENNA1")
        a2 = sub.getcol("ANTENNA2")
        sub.close()
    else:
        a1 = tb.getcol("ANTENNA1")
        a2 = tb.getcol("ANTENNA2")
    tb.close()
    pairs = {(int(i), int(j)) for i, j in zip(a1, a2) if i != j}
    pairs = {(min(i, j), max(i, j)) for i, j in pairs}
    return sorted(pairs)


def load_baseline_spectrum(
    ms_path: Path, a1: int, a2: int, corr_indices: list[int],
    field_id: int | None = None, scan_number: int | None = None,
) -> tuple[np.ndarray | None, float]:
    """Vector-average DATA over all rows in the given selection for baseline
    (a1, a2), using only unflagged samples in the given (parallel-hand)
    correlation indices.

    `a1` must be < `a2` -- callers always pass the canonical (lower, higher)
    antenna-index ordering (see `get_present_baselines`) so each physical
    baseline is only ever solved once. The TaQL query still matches both
    (ANTENNA1, ANTENNA2) storage orders so rows aren't silently dropped if
    the MS happens to store a pair the other way round.

    Returns (spectrum[nchan] complex, unflagged_fraction).
    """
    assert a1 < a2, f"expected canonical (lower, higher) antenna order, got ({a1}, {a2})"
    tb.open(str(ms_path))
    conds = [f"((ANTENNA1=={a1} && ANTENNA2=={a2}) || (ANTENNA1=={a2} && ANTENNA2=={a1}))"]
    if field_id is not None:
        conds.append(f"FIELD_ID=={field_id}")
    if scan_number is not None:
        conds.append(f"SCAN_NUMBER=={scan_number}")
    sub = tb.query(" && ".join(conds))
    if sub.nrows() == 0:
        sub.close()
        tb.close()
        return None, 0.0
    data = sub.getcol("DATA").T          # tb returns (ncorr, nchan, nrow) -> (nrow, nchan, ncorr)
    flag = sub.getcol("FLAG").T          # same axis flip
    flag_row = sub.getcol("FLAG_ROW")    # (nrow,)
    row_ant1 = sub.getcol("ANTENNA1")    # (nrow,)
    sub.close()
    tb.close()

    # Rows stored as (a2, a1) instead of the canonical (a1, a2) carry the
    # conjugate visibility (V_ji = conj(V_ij)) -- conjugate them back before
    # averaging in with the rest, so they don't corrupt the phase/delay.
    reversed_rows = row_ant1 == a2
    if reversed_rows.any():
        data[reversed_rows] = np.conj(data[reversed_rows])

    data = data[:, :, corr_indices]
    flag = flag[:, :, corr_indices] | flag_row[:, None, None]

    good = ~flag
    unflagged_frac = float(good.mean()) if good.size else 0.0

    data = np.where(good, data, 0.0)
    num = data.sum(axis=(0, 2))
    den = good.sum(axis=(0, 2))
    spectrum = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return spectrum, unflagged_frac


# --------------------------------------------------------------------------
# Delay / SNR estimation
# --------------------------------------------------------------------------

def compute_delay_spectrum(
    spectrum: np.ndarray, chan_freq: np.ndarray, zero_pad_factor: int = ZERO_PAD_FACTOR
) -> tuple[np.ndarray, np.ndarray]:
    """FFT a (Hanning-tapered, zero-padded) channel spectrum into an
    oversampled delay spectrum. Returns (delay_ns, amplitude)."""
    nchan = len(chan_freq)
    df = abs(float(np.median(np.diff(chan_freq))))  # Hz
    npad = int(nchan * zero_pad_factor)

    windowed = spectrum * np.hanning(nchan)
    padded = np.zeros(npad, dtype=complex)
    padded[:nchan] = windowed

    delay_spec = np.fft.fftshift(np.fft.fft(padded))
    delay_axis_s = np.fft.fftshift(np.fft.fftfreq(npad, d=df))
    return delay_axis_s * 1e9, np.abs(delay_spec)


def analytic_peak_widths(nchan: int, zero_pad_factor: int = ZERO_PAD_FACTOR) -> tuple[float, int]:
    """Measure the FWHM and first-null half-width (both in oversampled
    bins) of the delay-domain peak this pipeline's Hanning-tapered,
    zero-padded FFT produces for a perfect, noise-free zero-delay point
    source -- the exact peak shape, not an estimate from noisy data.
    Used to size the peak-exclusion window and noise annulus to the
    peak's real width instead of a fixed bin count."""
    window = np.hanning(nchan)
    npad = int(nchan * zero_pad_factor)
    padded = np.zeros(npad, dtype=complex)
    padded[:nchan] = window
    amp = np.abs(np.fft.fftshift(np.fft.fft(padded)))
    peak_idx = int(np.argmax(amp))
    peak_amp = amp[peak_idx]

    def _interp_crossing(threshold: float) -> float:
        """Bins from the peak (sub-bin, linearly interpolated) to where
        amp first drops below `threshold` walking outward from the peak."""
        i = peak_idx
        while i + 1 < len(amp) and amp[i + 1] >= threshold:
            i += 1
        if i + 1 >= len(amp):
            return float(i - peak_idx)
        a0, a1 = amp[i], amp[i + 1]
        frac = (a0 - threshold) / (a0 - a1) if a0 != a1 else 0.0
        return (i - peak_idx) + frac

    fwhm_bins = 2.0 * _interp_crossing(0.5 * peak_amp)

    i = peak_idx
    while i + 1 < len(amp) and amp[i + 1] < amp[i]:
        i += 1
    peak_half_width = int(np.ceil(i - peak_idx))

    return fwhm_bins, peak_half_width


def estimate_delay_and_snr(
    delay_ns: np.ndarray, amp: np.ndarray,
    half_width: int, noise_half_width: int,
) -> dict:
    n = len(amp)
    peak_idx = int(np.argmax(amp))
    lo = max(0, peak_idx - half_width)
    hi = min(n, peak_idx + half_width + 1)

    window_delay = delay_ns[lo:hi]
    window_amp = amp[lo:hi]
    weights = window_amp ** 2  # power-weighted centroid
    centroid_delay = (
        float(np.sum(window_delay * weights) / np.sum(weights))
        if weights.sum() > 0
        else float(delay_ns[peak_idx])
    )

    # Local noise annulus, excluding the full central peak out to its first null.
    noise_lo = max(0, peak_idx - noise_half_width)
    noise_hi = min(n, peak_idx + noise_half_width + 1)
    noise_mask = np.zeros(n, dtype=bool)
    noise_mask[noise_lo:noise_hi] = True
    noise_mask[lo:hi] = False

    if noise_mask.sum() < 5:
        # annulus too small (peak near an edge, or a short spectrum) --
        # fall back to everything outside the central peak
        noise_mask = np.ones(n, dtype=bool)
        noise_mask[lo:hi] = False
        noise_lo, noise_hi = 0, n
    noise_rms = float(np.std(amp[noise_mask])) if noise_mask.sum() >= 5 else float(np.std(amp))

    peak_amp = float(amp[peak_idx])
    snr = peak_amp / noise_rms if noise_rms > 0 else float("inf")

    return {
        "peak_idx": peak_idx,
        "peak_delay_ns": float(delay_ns[peak_idx]),
        "centroid_delay_ns": centroid_delay,
        "peak_amp": peak_amp,
        "noise_rms": noise_rms,
        "snr": snr,
        "window": (lo, hi),
        "noise_window": (noise_lo, noise_hi),
    }


def solve_baselines(
    ms_path: Path, pairs: list[tuple[int, int]], corr_indices: list[int],
    chan_freq: np.ndarray, antenna_names: list[str],
    half_width: int, noise_half_width: int,
    field_id: int | None = None, scan_number: int | None = None,
    min_unflagged_frac: float = MIN_UNFLAGGED_FRAC, single_baseline_test: bool = False,
    verbose: bool = True,
) -> list[dict]:
    """Solve delay/SNR for every baseline in `pairs`, within the given
    (field, scan) selection. Returns a list of per-baseline result dicts."""
    results: list[dict] = []
    for a1, a2 in pairs:
        spectrum, unflagged_frac = load_baseline_spectrum(
            ms_path, a1, a2, corr_indices, field_id=field_id, scan_number=scan_number)
        name1, name2 = antenna_names[a1], antenna_names[a2]
        if spectrum is None or unflagged_frac < min_unflagged_frac:
            if verbose:
                print(f"    skip {name1}-{name2}: unflagged fraction {unflagged_frac:.2f} "
                      f"< {min_unflagged_frac}")
            continue

        delay_ns, amp = compute_delay_spectrum(spectrum, chan_freq)
        res = estimate_delay_and_snr(delay_ns, amp, half_width, noise_half_width)
        res.update({
            "a1": a1, "a2": a2, "name1": name1, "name2": name2,
            "delay_ns_axis": delay_ns, "amp": amp,
            "unflagged_frac": unflagged_frac,
        })
        results.append(res)
        if verbose:
            print(f"    {name1}-{name2}: delay={res['centroid_delay_ns']:+.4f} ns, "
                  f"SNR={res['snr']:.1f}")

        if single_baseline_test:
            break

    return results


def best_per_antenna_from_baselines(
    baseline_results: list[dict],
) -> tuple[dict[int, dict], dict[int, list[float]]]:
    """Group baseline results by antenna, keeping the highest-SNR connected
    baseline per antenna. Returns (best_per_antenna, antenna_snrs) where
    antenna_snrs holds *all* connected baselines' SNRs per antenna."""
    best_per_antenna: dict[int, dict] = {}
    antenna_snrs: dict[int, list[float]] = {}
    for res in baseline_results:
        for ant, other, sign in ((res["a1"], res["a2"], +1.0), (res["a2"], res["a1"], -1.0)):
            antenna_snrs.setdefault(ant, []).append(res["snr"])
            cur = best_per_antenna.get(ant)
            if cur is None or res["snr"] > cur["snr"]:
                best_per_antenna[ant] = {
                    **res,
                    "signed_delay_ns": sign * res["centroid_delay_ns"],
                    "partner": other,
                }
    return best_per_antenna, antenna_snrs


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------

def _labeled_figure(
    n_items: int, per_item: float, min_width: float = 8.0, max_width: float = 60.0,
    height: float = 5.0, max_ticks: int = 120,
):
    """Create a fig/ax sized for n_items x-tick labels.

    The width grows with n_items but is capped at `max_width` inches so the
    rendered image can never approach matplotlib/Agg's 2**16 px-per-side
    limit (large antenna arrays can easily produce 1000+ baselines). If
    there'd be more labels than `max_ticks` can show legibly at that width,
    only every Nth one is drawn. Returns (fig, ax, tick_indices).
    """
    width = min(max_width, max(min_width, per_item * n_items))
    fig, ax = plt.subplots(figsize=(width, height))
    step = max(1, -(-n_items // max_ticks))  # ceil(n_items / max_ticks)
    tick_idx = list(range(0, n_items, step))
    return fig, ax, tick_idx


def median_mad_ci(values: list[float]) -> tuple[float, float, float]:
    """Median and 68% CI (1.4826x median absolute deviation). Returns
    (median, ci_lo, ci_hi)."""
    arr = np.asarray(values, dtype=float)
    med = float(np.median(arr))
    if arr.size < 2:
        return med, med, med
    sigma = 1.4826 * float(np.median(np.abs(arr - med)))
    return med, med - sigma, med + sigma


def plot_snr_per_scan(
    snr_by_scan: dict[int, dict[int, list[float]]], antenna_names: list[str], plot_dir: Path
) -> None:
    """SNR-vs-antenna, one series per scan: median +/- 68% CI over all
    baselines connected to that antenna in that scan."""
    ants = sorted({a for per_ant in snr_by_scan.values() for a in per_ant},
                  key=lambda i: antenna_names[i])
    labels = [antenna_names[i] for i in ants]
    scans = sorted(snr_by_scan.keys())

    fig, ax, tick_idx = _labeled_figure(len(labels), per_item=0.35)
    x = np.arange(len(labels))
    cmap = plt.get_cmap("tab10" if len(scans) <= 10 else "tab20")
    offsets = np.linspace(-0.3, 0.3, len(scans)) if len(scans) > 1 else [0.0]
    for si, scan in enumerate(scans):
        meds, err_lo, err_hi = [], [], []
        for a in ants:
            vals = snr_by_scan[scan].get(a)
            if not vals:
                meds.append(np.nan); err_lo.append(0.0); err_hi.append(0.0)
                continue
            med, lo, hi = median_mad_ci(vals)
            meds.append(med)
            err_lo.append(med - lo)
            err_hi.append(hi - med)
        ax.errorbar(x + offsets[si], meds, yerr=[err_lo, err_hi], fmt="o", ms=4,
                    color=cmap(si % cmap.N), ecolor=cmap(si % cmap.N),
                    elinewidth=1, capsize=2, alpha=0.85, zorder=3, label=f"scan {scan}")
    ax.grid(axis="y", color="0.85", lw=0.6, zorder=0)
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([labels[i] for i in tick_idx], rotation=90, fontsize=8)
    ax.set_ylabel("Delay-solve SNR (median ± 68% CI)")
    ax.set_title("Per-scan delay-solve SNR by antenna (primary calibrator)")
    ax.legend(fontsize=7, ncol=min(len(scans), 6))
    fig.tight_layout()
    fig.savefig(plot_dir / "refant_snr_per_scan.png", dpi=150)
    plt.close(fig)


def plot_antenna_snr_combined(
    antenna_snrs: dict[int, list[float]], antenna_names: list[str], plot_dir: Path
) -> None:
    """Median SNR +/- 68% CI over all connected-baseline SNRs, all scans
    pooled -- this is the population the ranking is drawn from."""
    ants = sorted(antenna_snrs.keys(), key=lambda i: antenna_names[i])
    labels = [antenna_names[i] for i in ants]

    medians, err_lo, err_hi = [], [], []
    for i in ants:
        med, lo, hi = median_mad_ci(antenna_snrs[i])
        medians.append(med)
        err_lo.append(med - lo)
        err_hi.append(hi - med)

    fig, ax, tick_idx = _labeled_figure(len(labels), per_item=0.35)
    x = np.arange(len(labels))
    ax.errorbar(x, medians, yerr=[err_lo, err_hi], fmt="none", ecolor="0.4",
                elinewidth=1, capsize=3, zorder=2)
    ax.scatter(x, medians, c="#2ca02c", s=70, zorder=3, edgecolor="k", linewidth=0.5)
    ax.grid(axis="y", color="0.85", lw=0.6, zorder=0)
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([labels[i] for i in tick_idx], rotation=90, fontsize=8)
    ax.set_ylabel("Median SNR, all scans combined (68% CI)")
    ax.set_title("Combined (all-scan) delay-solve SNR by antenna")
    fig.tight_layout()
    fig.savefig(plot_dir / "refant_snr_combined.png", dpi=150)
    plt.close(fig)


def plot_best_baseline_fft_per_antenna(
    best_per_antenna: dict[int, dict], antenna_names: list[str], plot_dir: Path,
    zoom_half_width: int,
) -> None:
    out_dir = plot_dir / "refant_fft_per_antenna"
    out_dir.mkdir(exist_ok=True)

    for ant_idx, res in best_per_antenna.items():
        if "delay_ns_axis" not in res:
            continue  # no baseline data at all for this antenna

        delay_ns = res["delay_ns_axis"]
        amp = res["amp"]
        lo, hi = res["window"]
        noise_lo, noise_hi = res["noise_window"]
        ant_name = antenna_names[ant_idx]
        partner_name = antenna_names[res["partner"]]

        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(delay_ns, amp, color="0.4", lw=1, label="|FFT|")
        ax.axvspan(delay_ns[noise_lo], delay_ns[noise_hi - 1], color="steelblue", alpha=0.12,
                    label="local noise annulus")
        ax.axvspan(delay_ns[lo], delay_ns[hi - 1], color="orange", alpha=0.3,
                    label="central peak")
        ax.axvline(res["centroid_delay_ns"], color="crimson", lw=1.5,
                    label=f"centroid = {res['centroid_delay_ns']:.4f} ns")
        ax.axhline(res["noise_rms"], color="steelblue", ls="--", lw=1,
                    label=f"local noise RMS (SNR={res['snr']:.1f})")
        ax.set_xlabel("Delay (ns)")
        ax.set_ylabel("Amplitude")
        ax.set_title(f"{ant_name}: best baseline -> {partner_name} (SNR={res['snr']:.1f})")
        ax.legend(fontsize=7)

        # Zoomed-in inset (top-left): peak +/- zoom_half_width (in FWHM),
        # so the sub-bin centroid location is visible even when it sits
        # well inside a single pixel of the full-width plot.
        peak_idx = res["peak_idx"]
        zlo = max(0, peak_idx - zoom_half_width)
        zhi = min(len(delay_ns), peak_idx + zoom_half_width + 1)

        axins = ax.inset_axes([0.04, 0.54, 0.42, 0.42])
        axins.plot(delay_ns[zlo:zhi], amp[zlo:zhi], color="0.4", lw=1.2, marker=".", ms=3)
        axins.axvspan(delay_ns[noise_lo], delay_ns[noise_hi - 1], color="steelblue", alpha=0.12)
        axins.axvspan(delay_ns[lo], delay_ns[hi - 1], color="orange", alpha=0.3)
        axins.axvline(res["centroid_delay_ns"], color="crimson", lw=1.5)
        axins.axhline(res["noise_rms"], color="steelblue", ls="--", lw=1)
        axins.set_xlim(delay_ns[zlo], delay_ns[zhi - 1])
        axins.set_title("zoom near peak", fontsize=7, pad=2)
        axins.tick_params(labelsize=6)
        ax.indicate_inset_zoom(axins, edgecolor="0.3")
        fig.tight_layout()
        fig.savefig(out_dir / f"fft_{ant_name}_best_{partner_name}.png", dpi=150)
        plt.close(fig)


# --------------------------------------------------------------------------
# Ranking / reporting
# --------------------------------------------------------------------------

def print_best_worst(
    avg_snr_by_ant: dict[int, float], antenna_names: list[str], top_n: int
) -> list[tuple[str, float]]:
    """Print the top-N and bottom-N antennas by average SNR. Returns the
    full ranking as (name, avg_snr) sorted best-first."""
    ranked = sorted(
        ((antenna_names[a], snr) for a, snr in avg_snr_by_ant.items()),
        key=lambda x: x[1], reverse=True,
    )
    n = min(top_n, len(ranked))
    print(f"\nAntenna delay-solve S/N ranking ({len(ranked)} antennas, "
          f"average over all connected baselines, all scans, of the primary):")
    print(f"  BEST {n}:")
    for name, snr in ranked[:n]:
        print(f"    {name}: {snr:.1f}")
    print(f"  WORST {n}:")
    for name, snr in ranked[-n:][::-1]:
        print(f"    {name}: {snr:.1f}")
    return ranked


def write_ranking_json(ranked: list[tuple[str, float]], plot_dir: Path) -> Path:
    out_path = plot_dir / "refant_snr_ranking.json"
    with open(out_path, "w") as f:
        json.dump([{"antenna": name, "avg_snr": snr} for name, snr in ranked], f, indent=2)
    return out_path


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

t_start = time.perf_counter()

ms_path = Path(myms)
plot_dir = Path(GAINPLOTS)
plot_dir.mkdir(parents=True, exist_ok=True)

antenna_names = get_antenna_names(ms_path)

field_id = get_field_id(ms_path, bpcal_name)
scans = get_scans_for_field(ms_path, field_id)
print(f"Field '{bpcal_name}' (ID {field_id}): {len(scans)} scan(s): {scans}")

corr_indices, corr_labels = get_parallel_hand_indices(ms_path)
print(f"Using parallel-hand correlations: {corr_labels}")

chan_freq = get_spw_chan_freq(ms_path)
print(f"{len(chan_freq)} channels, {chan_freq[0] / 1e9:.4f}-{chan_freq[-1] / 1e9:.4f} GHz")

fwhm_bins, peak_half_width = analytic_peak_widths(len(chan_freq), ZERO_PAD_FACTOR)
noise_half_width = peak_half_width + int(round(NOISE_ANNULUS_FWHM * fwhm_bins))
zoom_half_width = int(round(ZOOM_FWHM * fwhm_bins))
print(f"Analytic peak shape: FWHM = {fwhm_bins:.2f} bins, central peak half-width "
      f"(first null) = {peak_half_width} bins, noise annulus half-width = {noise_half_width} bins")

snr_by_scan: dict[int, dict[int, list[float]]] = {}
combined_snrs: dict[int, list[float]] = {}
combined_best_per_antenna: dict[int, dict] = {}  # for the FFT diagnostic plots

for scan in scans:
    t_scan = time.perf_counter()
    pairs = get_present_baselines(ms_path, field_id=field_id, scan_number=scan)
    print(f"\nScan {scan}: {len(pairs)} baseline(s)")
    if pairs:
        baseline_results = solve_baselines(
            ms_path, pairs, corr_indices, chan_freq, antenna_names,
            peak_half_width, noise_half_width,
            field_id=field_id, scan_number=scan,
            single_baseline_test=SINGLE_BASELINE_TEST)
        if not baseline_results:
            print(f"  no usable baselines in scan {scan}")
        else:
            best_this_scan, all_this_scan = best_per_antenna_from_baselines(baseline_results)
            snr_by_scan[scan] = all_this_scan
            for a, snr_list in all_this_scan.items():
                combined_snrs.setdefault(a, []).extend(snr_list)
            for a, r in best_this_scan.items():
                cur = combined_best_per_antenna.get(a)
                if cur is None or r["snr"] > cur["snr"]:
                    combined_best_per_antenna[a] = r

    print(f"  scan {scan} done in {time.perf_counter() - t_scan:.1f}s "
          f"(total elapsed {time.perf_counter() - t_start:.1f}s)")

    if SINGLE_BASELINE_TEST:
        break

if not combined_snrs:
    raise RuntimeError("No usable baselines found in any scan (check flags / correlation selection).")

avg_snr_by_ant = {a: float(np.mean(vals)) for a, vals in combined_snrs.items()}

refant_snr_ranking = print_best_worst(avg_snr_by_ant, antenna_names, TOP_N)
ranking_path = write_ranking_json(refant_snr_ranking, plot_dir)
print(f"\nRanking written to {ranking_path}")

t_plot = time.perf_counter()
plot_snr_per_scan(snr_by_scan, antenna_names, plot_dir)
plot_antenna_snr_combined(combined_snrs, antenna_names, plot_dir)
plot_best_baseline_fft_per_antenna(combined_best_per_antenna, antenna_names, plot_dir, zoom_half_width)
print(f"Plots written to {plot_dir} ({time.perf_counter() - t_plot:.1f}s)")

print(f"Total runtime: {time.perf_counter() - t_start:.1f}s")
