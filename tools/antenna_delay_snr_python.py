#!/usr/bin/env python3
"""
antenna_delay_snr.py
=====================
Per-antenna delay-solve S/N ranking for reference-antenna selection, based
on a primary/delay-calibrator field.

Method
------
For each baseline (a1, a2), the parallel-hand correlations (XX/YY, or
RR/LL if the array is circular) are vector-averaged over all unflagged
data in the selection (both parallel-hand pols) to build one complex
spectrum per baseline:

    V_ij(f) ~= A * exp(-2*pi*j*f*tau_ij) + noise

That spectrum is Hanning-tapered and zero-padded before being FFT'd
across frequency, which oversamples the resulting delay axis well
beyond the native 1/bandwidth resolution -- this is what lets very
small (near-zero) delays be resolved, not just large ones.

The delay is then estimated as the power-weighted **centroid** of
|FFT| within a small window around its peak (sub-bin precision),
rather than just the discrete peak bin.

A per-baseline signal-to-noise ratio is estimated from a *local* noise
annulus around the peak (+/- NOISE_HALF_WIDTH bins, excluding the
centroid window itself) rather than the whole spectrum, so far-off
structure -- bandpass ripple, RFI channels, edge/aliasing effects --
can't drag the RMS around:

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
  spacing). Multi-SPW MSs are not handled.
- DATA column (not CORRECTED_DATA).

Requires: python-casacore (pyrap.tables / casacore.tables), numpy,
matplotlib.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

try:
    from pyrap.tables import table
except ImportError:  # fall back to the modern package name
    from casacore.tables import table

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

ZERO_PAD_FACTOR = 8       # oversampling factor applied before the FFT
CENTROID_HALF_WIDTH = 4   # +/- bins (in the oversampled spectrum) for the centroid itself
NOISE_HALF_WIDTH_MULTIPLE = 200  # NOISE_HALF_WIDTH = this many multiples of CENTROID_HALF_WIDTH
NOISE_HALF_WIDTH = NOISE_HALF_WIDTH_MULTIPLE * CENTROID_HALF_WIDTH  # +/- bins for the local noise annulus (default 800)
ZOOM_HALF_WIDTH_MULTIPLE = 5    # per-antenna FFT plot zooms to peak +/- this many x CENTROID_HALF_WIDTH
MIN_UNFLAGGED_FRAC = 0.1  # skip a baseline if less than this fraction of samples are unflagged
TOP_N_DEFAULT = 7         # default number of best/worst antennas to print/rank

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
    with table(str(ms_path / "ANTENNA"), ack=False) as ant:
        return list(ant.getcol("NAME"))


def get_field_id(ms_path: Path, field_name: str) -> int:
    with table(str(ms_path / "FIELD"), ack=False) as fld:
        names = list(fld.getcol("NAME"))
    matches = [i for i, n in enumerate(names) if n == field_name]
    if not matches:
        raise RuntimeError(f"Field '{field_name}' not found in {ms_path}::FIELD (have {names})")
    if len(matches) > 1:
        raise RuntimeError(f"Field '{field_name}' is ambiguous in {ms_path}::FIELD (IDs {matches})")
    return matches[0]


def get_scans_for_field(ms_path: Path, field_id: int) -> list[int]:
    with table(str(ms_path), ack=False) as ms:
        sub = ms.query(f"FIELD_ID=={field_id}")
        scans = sorted({int(s) for s in sub.getcol("SCAN_NUMBER")})
        sub.close()
    return scans


def get_parallel_hand_indices(ms_path: Path) -> tuple[list[int], list[str]]:
    with table(str(ms_path / "POLARIZATION"), ack=False) as pol:
        corr_type = pol.getcol("CORR_TYPE")[0]
    idx = [i for i, c in enumerate(corr_type) if c in PARALLEL_HAND_NAMES]
    if not idx:
        raise RuntimeError(
            f"No parallel-hand (XX/YY or RR/LL) correlations found in CORR_TYPE={corr_type}"
        )
    labels = [PARALLEL_HAND_NAMES[corr_type[i]] for i in idx]
    return idx, labels


def get_spw_chan_freq(ms_path: Path, spw_id: int = 0) -> np.ndarray:
    with table(str(ms_path / "SPECTRAL_WINDOW"), ack=False) as spw:
        chan_freq = spw.getcol("CHAN_FREQ")[spw_id]  # Hz
    return np.asarray(chan_freq, dtype=float)


def get_present_baselines(
    ms_path: Path, field_id: int | None = None, scan_number: int | None = None
) -> list[tuple[int, int]]:
    with table(str(ms_path), ack=False) as ms:
        sub = ms
        conds = []
        if field_id is not None:
            conds.append(f"FIELD_ID=={field_id}")
        if scan_number is not None:
            conds.append(f"SCAN_NUMBER=={scan_number}")
        if conds:
            sub = ms.query(" && ".join(conds))
        a1 = sub.getcol("ANTENNA1")
        a2 = sub.getcol("ANTENNA2")
        if sub is not ms:
            sub.close()
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
    with table(str(ms_path), ack=False) as ms:
        conds = [f"((ANTENNA1=={a1} && ANTENNA2=={a2}) || (ANTENNA1=={a2} && ANTENNA2=={a1}))"]
        if field_id is not None:
            conds.append(f"FIELD_ID=={field_id}")
        if scan_number is not None:
            conds.append(f"SCAN_NUMBER=={scan_number}")
        sub = ms.query(" && ".join(conds))
        if sub.nrows() == 0:
            sub.close()
            return None, 0.0
        data = sub.getcol("DATA")            # (nrow, nchan, ncorr)
        flag = sub.getcol("FLAG")            # (nrow, nchan, ncorr)
        flag_row = sub.getcol("FLAG_ROW")    # (nrow,)
        row_ant1 = sub.getcol("ANTENNA1")    # (nrow,)
        sub.close()

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


def estimate_delay_and_snr(
    delay_ns: np.ndarray, amp: np.ndarray,
    half_width: int = CENTROID_HALF_WIDTH, noise_half_width: int = NOISE_HALF_WIDTH,
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

    # Noise is measured in a *local* annulus around the peak (excluding the
    # centroid window itself), not across the whole (heavily oversampled)
    # spectrum -- that keeps far-off structure (bandpass ripple, RFI
    # channels, edge/aliasing effects) from dragging the RMS around and
    # gives a noise floor that's actually representative of the region
    # right next to the peak.
    noise_lo = max(0, peak_idx - noise_half_width)
    noise_hi = min(n, peak_idx + noise_half_width + 1)
    noise_mask = np.zeros(n, dtype=bool)
    noise_mask[noise_lo:noise_hi] = True
    noise_mask[lo:hi] = False

    if noise_mask.sum() < 5:
        # annulus too small (peak near an edge, or a short spectrum) --
        # fall back to everything outside the centroid window
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
        res = estimate_delay_and_snr(delay_ns, amp)
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


def median_ci(
    values: list[float], ci: float = 0.95, n_boot: int = 2000, seed: int = 0
) -> tuple[float, float, float]:
    """Bootstrap the median of `values` and a percentile confidence interval.
    Returns (median, ci_lo, ci_hi). With fewer than 2 samples the CI just
    collapses onto the (single) value -- there's nothing to resample."""
    arr = np.asarray(values, dtype=float)
    med = float(np.median(arr))
    if arr.size < 2:
        return med, med, med
    rng = np.random.default_rng(seed)
    boot_idx = rng.integers(0, arr.size, size=(n_boot, arr.size))
    boot_medians = np.median(arr[boot_idx], axis=1)
    lo, hi = np.percentile(boot_medians, [100 * (1 - ci) / 2, 100 * (1 + ci) / 2])
    return med, float(lo), float(hi)


def plot_snr_per_scan(
    snr_by_scan: dict[int, dict[int, float]], antenna_names: list[str], plot_dir: Path
) -> None:
    """SNR-vs-antenna, one series per scan -- shows whether an antenna's
    delay S/N is consistent scan-to-scan or just a one-off."""
    ants = sorted({a for per_ant in snr_by_scan.values() for a in per_ant},
                  key=lambda i: antenna_names[i])
    labels = [antenna_names[i] for i in ants]
    scans = sorted(snr_by_scan.keys())

    fig, ax, tick_idx = _labeled_figure(len(labels), per_item=0.35)
    x = np.arange(len(labels))
    cmap = plt.get_cmap("tab10" if len(scans) <= 10 else "tab20")
    for si, scan in enumerate(scans):
        y = [snr_by_scan[scan].get(a, np.nan) for a in ants]
        ax.scatter(x, y, s=45, color=cmap(si % cmap.N), edgecolor="k", linewidth=0.4,
                   zorder=3, label=f"scan {scan}")
    ax.grid(axis="y", color="0.85", lw=0.6, zorder=0)
    ax.set_xticks(tick_idx)
    ax.set_xticklabels([labels[i] for i in tick_idx], rotation=90, fontsize=8)
    ax.set_ylabel("Delay-solve SNR")
    ax.set_title("Per-scan delay-solve SNR by antenna (primary calibrator)")
    ax.legend(fontsize=7, ncol=min(len(scans), 6))
    fig.tight_layout()
    fig.savefig(plot_dir / "refant_snr_per_scan.png", dpi=150)
    plt.close(fig)


def plot_antenna_snr_combined(
    antenna_snrs: dict[int, list[float]], antenna_names: list[str], plot_dir: Path
) -> None:
    """Median SNR (with a bootstrap CI) pooling every scan's per-antenna
    measurement together ("combine the two populations") -- this is the
    population the average-SNR reference-antenna ranking is drawn from."""
    ants = sorted(antenna_snrs.keys(), key=lambda i: antenna_names[i])
    labels = [antenna_names[i] for i in ants]

    medians, err_lo, err_hi = [], [], []
    for i in ants:
        med, lo, hi = median_ci(antenna_snrs[i])
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
    ax.set_ylabel("Median SNR, all scans combined (95% bootstrap CI)")
    ax.set_title("Combined (all-scan) delay-solve SNR by antenna")
    fig.tight_layout()
    fig.savefig(plot_dir / "refant_snr_combined.png", dpi=150)
    plt.close(fig)


def plot_best_baseline_fft_per_antenna(
    best_per_antenna: dict[int, dict], antenna_names: list[str], plot_dir: Path
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
                    label="centroid window")
        ax.axvline(res["centroid_delay_ns"], color="crimson", lw=1.5,
                    label=f"centroid = {res['centroid_delay_ns']:.4f} ns")
        ax.axhline(res["noise_rms"], color="steelblue", ls="--", lw=1,
                    label=f"local noise RMS (SNR={res['snr']:.1f})")
        ax.set_xlabel("Delay (ns)")
        ax.set_ylabel("Amplitude")
        ax.set_title(f"{ant_name}: best baseline -> {partner_name} (SNR={res['snr']:.1f})")
        ax.legend(fontsize=7)

        # Zoomed-in inset (top-left): peak +/- ZOOM_HALF_WIDTH_MULTIPLE x
        # CENTROID_HALF_WIDTH, so the sub-bin centroid location is visible
        # even when it sits well inside a single pixel of the full-width
        # plot. Deliberately sized off the centroid window, not the (much
        # wider) noise annulus -- the noise annulus may just clip in at the
        # edges of the zoom, or not appear at all, and that's fine.
        peak_idx = res["peak_idx"]
        zoom_half_width = ZOOM_HALF_WIDTH_MULTIPLE * CENTROID_HALF_WIDTH
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
          f"average over all scans of the primary):")
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

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("myms", type=Path, help="Path to the MS")
    parser.add_argument("field", type=str, help="Field name to solve on (e.g. the primary/bpcal)")
    parser.add_argument("--outdir", type=Path, required=True,
                         help="Directory for ranking JSON + diagnostic plots (e.g. GAINPLOTS)")
    parser.add_argument("--top-n", type=int, default=TOP_N_DEFAULT,
                         help=f"Number of best/worst antennas to print (default {TOP_N_DEFAULT})")
    args = parser.parse_args()

    t_start = time.perf_counter()

    ms_path: Path = args.myms
    plot_dir: Path = args.outdir
    plot_dir.mkdir(parents=True, exist_ok=True)

    antenna_names = get_antenna_names(ms_path)

    field_id = get_field_id(ms_path, args.field)
    scans = get_scans_for_field(ms_path, field_id)
    print(f"Field '{args.field}' (ID {field_id}): {len(scans)} scan(s): {scans}")

    corr_indices, corr_labels = get_parallel_hand_indices(ms_path)
    print(f"Using parallel-hand correlations: {corr_labels}")

    chan_freq = get_spw_chan_freq(ms_path)
    print(f"{len(chan_freq)} channels, {chan_freq[0] / 1e9:.4f}-{chan_freq[-1] / 1e9:.4f} GHz")

    # --- Solve independently per scan, keep each antenna's best-baseline
    # SNR per scan, and pool ("combine") those per-scan values across scans.
    snr_by_scan: dict[int, dict[int, float]] = {}
    combined_snrs: dict[int, list[float]] = {}
    combined_best_per_antenna: dict[int, dict] = {}  # for the FFT diagnostic plots

    for scan in scans:
        t_scan = time.perf_counter()
        pairs = get_present_baselines(ms_path, field_id=field_id, scan_number=scan)
        print(f"\nScan {scan}: {len(pairs)} baseline(s)")
        if not pairs:
            continue

        baseline_results = solve_baselines(
            ms_path, pairs, corr_indices, chan_freq, antenna_names,
            field_id=field_id, scan_number=scan,
            single_baseline_test=SINGLE_BASELINE_TEST)
        if not baseline_results:
            print(f"  no usable baselines in scan {scan}")
            continue

        best_this_scan, _ = best_per_antenna_from_baselines(baseline_results)
        snr_by_scan[scan] = {a: r["snr"] for a, r in best_this_scan.items()}
        for a, r in best_this_scan.items():
            combined_snrs.setdefault(a, []).append(r["snr"])
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

    ranked = print_best_worst(avg_snr_by_ant, antenna_names, args.top_n)
    ranking_path = write_ranking_json(ranked, plot_dir)
    print(f"\nRanking written to {ranking_path}")

    t_plot = time.perf_counter()
    plot_snr_per_scan(snr_by_scan, antenna_names, plot_dir)
    plot_antenna_snr_combined(combined_snrs, antenna_names, plot_dir)
    plot_best_baseline_fft_per_antenna(combined_best_per_antenna, antenna_names, plot_dir)
    print(f"Plots written to {plot_dir} ({time.perf_counter() - t_plot:.1f}s)")

    print(f"Total runtime: {time.perf_counter() - t_start:.1f}s")


if __name__ == "__main__":
    main()
