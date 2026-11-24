"""Benchmark of the pile-up methods in ``pileup.py`` (step 5).

Three tests, each comparing every method's estimate with **the same method's estimate of
the same pulse alone**, so differences in scale between methods (TM reads a sample, HYPER
extrapolates an exponential) cancel and only the pile-up error remains:

A. **Template pairs in real noise.** Two template pulses (random sub-sample phases) at a
   controlled separation and amplitude ratio are injected into real noise, *found by the
   real detector*, clustered as in the pipeline, and resolved with the detected positions.
B. **Real-pulse pairs.** Two real isolated pulses (with their real shape variations and
   noise) are added with a controlled separation and resolved with their known positions.
C. **Real pile-up events.** The ⁶⁰Co pile-up clusters stored in an event file are resolved
   by each method. The ⁶⁰Co peak positions and widths must match those the same method
   gives for isolated pulses.

Run:  python -m gamma_spectrometry evaluate-pileup STORE TEMPLATE EVENTS OUT_DIR
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import plotting as P
from .baseline import design_from_template, filter_template
from .detect import Detector, decay_constant, noise_autocovariance
from .detection_eval import noise_streams
from .events import FLAG_EDGE, FLAG_LARGE, FLAG_WARMUP, EventFile
from .pileup import METHODS, PileupContext, resolve
from .preprocess import PreprocessConfig, Preprocessor
from .spectrum import fit_peaks
from .store import CaptureStore
from .template import PulseTemplate

SEPARATIONS = (3, 4, 5, 6, 8, 10, 15, 20, 30, 45)
PAIRS_V = {"equal 300/300 mV": (0.3, 0.3), "big then small 300/30 mV": (0.3, 0.03),
           "small then big 30/300 mV": (0.03, 0.3)}
LABELS = {"tm": "template matching (paper)", "hyper": "HYPER", "ls": "least squares, integer times",
          "ls_refine": "least squares + timing", "gls_refine": "whitened LS (GLS) + timing",
          "wls_refine": "noise + mismatch weighted LS + timing"}
PRE, POST = 20, 80          # snippet: from first pulse - PRE to last pulse + POST (as in the pipeline)


def _robust(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict(n=0, bias=np.nan, spread=np.nan)
    med = float(np.median(x))
    return dict(n=int(x.size), bias=med, spread=float(1.4826 * np.median(np.abs(x - med))))


def make_context(template_path, noise=None):
    T = PulseTemplate.load(template_path)
    pz = design_from_template(T)
    Tp = filter_template(T, pz.b, pz.a)
    acov = None
    if noise is not None:
        acov = noise_autocovariance(noise, np.ones(noise.size, bool), 4096)
    return T, pz, Tp, PileupContext(Tp, decay_constant(Tp), acov)


# ----------------------------------------------------------------------------------------------
# Test A: template pairs in real noise, detected by the real detector
# ----------------------------------------------------------------------------------------------

def test_template_pairs(Tp, ctx, noise, rng, n_max=800):
    det = Detector(Tp)
    out = {}
    t_single = np.arange(PRE + POST + 1, dtype=float)
    for pname, (a1, a2) in PAIRS_V.items():
        out[pname] = {}
        for d in SEPARATIONS:
            spacing = 600
            starts = np.arange(400, noise.size - 600, spacing)[:n_max]
            c1 = starts + rng.uniform(0, 1, starts.size)
            c2 = c1 + d
            y = noise.copy()
            k = np.arange(-10, 300)
            for x1, x2 in zip(c1, c2):
                for c, a in ((x1, a1), (x2, a2)):
                    idx = int(np.floor(c)) + k
                    idx = idx[idx < y.size]
                    y[idx] += a * Tp(idx - c)
            found = det.detect(y).index
            errs = {m: ([], []) for m in METHODS}
            n_pair = 0
            for x1, x2 in zip(c1, c2):
                near = found[(found >= x1 - 40) & (found <= x2 + 60)]
                if near.size != 2 or abs(near[0] - x1) > 2 or abs(near[1] - x2) > 2:
                    continue
                n_pair += 1
                s0 = int(near[0]) - PRE
                seg = y[s0: int(near[1]) + POST + 1]
                pos = near - s0
                for m in METHODS:
                    est = resolve(m, seg, pos, ctx)
                    for j, (c, a, p) in enumerate(((x1, a1, pos[0]), (x2, a2, pos[1]))):
                        ref = resolve(m, a * Tp(t_single - (PRE + (c - np.round(c)))), np.array([PRE]), ctx)[0]
                        errs[m][j].append(est[j] / ref - 1)
            out[pname][str(d)] = dict(
                detected_as_pair=n_pair / starts.size,
                methods={m: dict(first=_robust(e[0]), second=_robust(e[1])) for m, e in errs.items()})
    return out


# ----------------------------------------------------------------------------------------------
# Real isolated pulses (for tests B and C)
# ----------------------------------------------------------------------------------------------

def real_isolated_pulses(store_path, Tp, pz, captures, a_lo=0.03, a_hi=0.7):
    """Windows y[p-PRE : p+POST+60] of isolated real pulses (clean signal) and their peak index PRE."""
    wins = []
    with CaptureStore(store_path) as s:
        pre = Preprocessor(1.0 / s.dt_s, PreprocessConfig(), pz)
        det = Detector(Tp)
        for i in captures:
            pre.reset()
            y = pre.process(s.volts(i))
            d = det.detect(y)
            idx, amp = d.index, d.filter_amplitude
            g_prev = np.diff(np.r_[-10**9, idx])
            g_next = np.diff(np.r_[idx, 10**12])
            ok = (g_prev > 300) & (g_next > 300) & (idx > 1500) & (idx < y.size - 400) & (amp > a_lo) & (amp < a_hi)
            for p in idx[ok]:
                wins.append(y[p - PRE: p + POST + 61])
    return np.array(wins)


# ----------------------------------------------------------------------------------------------
# Test B: pairs of real pulses
# ----------------------------------------------------------------------------------------------

def test_real_pairs(W, ctx, rng, n_per_d=600, n_pool=4000):
    L = W.shape[1]
    singles = {m: np.array([resolve(m, w, np.array([PRE]), ctx)[0] for w in W]) for m in METHODS}
    out = {}
    for d in SEPARATIONS:
        i1 = rng.integers(0, min(len(W), n_pool), n_per_d)
        i2 = rng.integers(0, min(len(W), n_pool), n_per_d)
        errs = {m: ([], []) for m in METHODS}
        for a, b in zip(i1, i2):
            if a == b:
                continue
            y = np.zeros(L + d)
            y[:L] += W[a]
            y[d:] += W[b]
            y = y[: PRE + d + POST + 1]
            pos = np.array([PRE, PRE + d])
            for m in METHODS:
                est = resolve(m, y, pos, ctx)
                errs[m][0].append(est[0] / singles[m][a] - 1)
                errs[m][1].append(est[1] / singles[m][b] - 1)
        out[str(d)] = {m: dict(first=_robust(e[0]), second=_robust(e[1])) for m, e in errs.items()}
    return out, singles


# ----------------------------------------------------------------------------------------------
# Test C: the real pile-up events
# ----------------------------------------------------------------------------------------------

def test_real_events(events_path, singles, ctx, co_guess=(0.316, 0.357)):
    ef = EventFile(events_path)
    ev = ef.events
    bad = FLAG_WARMUP | FLAG_EDGE | FLAG_LARGE
    amps = {m: [] for m in METHODS}
    gaps = []
    size = []
    clusters = np.unique(ev["cluster"][(ev["cluster_size"] >= 2) & (ev["cluster_size"] <= 4)])
    by_cluster = {}
    for i in np.flatnonzero(np.isin(ev["cluster"], clusters)):
        by_cluster.setdefault(int(ev["cluster"][i]), []).append(i)
    for cl, rows in by_cluster.items():
        e = ev[rows]
        if np.any(e["flags"] & bad):
            continue
        try:
            _, start, y = ef.snippet(cl)
        except KeyError:
            continue
        pos = e["index"] - start
        for m in METHODS:
            amps[m].extend(resolve(m, y, pos, ctx))
        g = np.minimum(np.where(e["gap_prev"] > 0, e["gap_prev"], 10**6), np.where(e["gap_next"] > 0, e["gap_next"], 10**6))
        gaps.extend(g)
        size.extend([e.size] * e.size)
    gaps = np.array(gaps)
    res = {}
    for m in METHODS:
        a = np.array(amps[m])
        iso = singles[m]
        scale = 1.0
        row = {}
        try:
            ref = fit_peaks(iso, co_guess, co_guess[0] * 0.85, co_guess[1] * 1.12, 0.002)
            scale = float(np.mean(np.asarray(co_guess) / ref.centers))
            row["isolated"] = dict(centers=ref.centers.tolist(), fwhm_pct=(100 * ref.resolution).tolist(), n=int(iso.size))
        except Exception as exc:  # noqa: BLE001
            row["isolated"] = dict(error=str(exc))
        for label, sel in (("all pile-up", np.ones(a.size, bool)), ("gap < 10", gaps < 10),
                           ("gap 10-30", (gaps >= 10) & (gaps < 30)), ("gap 30-60", (gaps >= 30) & (gaps < 60))):
            try:
                g = np.asarray(co_guess) / scale
                pf = fit_peaks(a[sel], g, g[0] * 0.85, g[1] * 1.12, 0.002)
                row[label] = dict(n=int(sel.sum()), centers=pf.centers.tolist(),
                                  shift_pct=(100 * (pf.centers / ref.centers - 1)).tolist(),
                                  fwhm_pct=(100 * pf.resolution).tolist())
            except Exception as exc:  # noqa: BLE001
                row[label] = dict(n=int(sel.sum()), error=str(exc))
        res[m] = row
    return res, {m: np.array(v) for m, v in amps.items()}, gaps


# ----------------------------------------------------------------------------------------------
# Driver and plots
# ----------------------------------------------------------------------------------------------

def evaluate(store_path, template_path, events_path, out_dir, captures=range(0, 100, 5), seed=0,
             verbose=True) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    T, pz, Tp, _ = make_context(template_path)
    real, _, _ = noise_streams(store_path, Tp, pz, list(captures), seed)
    _, _, _, ctx = make_context(template_path, real)
    if verbose:
        print("A: template pairs in real noise ...", flush=True)
    A = test_template_pairs(Tp, ctx, real, rng)
    if verbose:
        print("B: pairs of real pulses ...", flush=True)
    # all captures for the isolated-pulse references (their peak fits need the statistics)
    with CaptureStore(store_path) as st:
        all_caps = list(range(len(st)))
    W = real_isolated_pulses(store_path, Tp, pz, all_caps)
    B, singles = test_real_pairs(W, ctx, rng)
    if verbose:
        print(f"   {W.shape[0]} real isolated pulses; C: real pile-up events ...", flush=True)
    C, amps, gaps = test_real_events(events_path, singles, ctx)
    summary = dict(template_pairs=A, real_pairs=B, real_events=C, n_real_isolated=int(W.shape[0]),
                   speed_us_per_cluster=_speed(W, ctx))
    (out / "pileup_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _plots(out, A, B, C, singles, amps, gaps)
    return summary


def _speed(W, ctx, n=300):
    """Microseconds per two-pulse cluster for each method (single core)."""
    import time
    L = W.shape[1]
    ys = []
    for k in range(n):
        y = np.zeros(L + 8)
        y[:L] += W[k % len(W)]
        y[8:] += W[(k + 1) % len(W)]
        ys.append(y[: PRE + 8 + POST + 1])
    out = {}
    for m in METHODS:
        t0 = time.perf_counter()
        for y in ys:
            resolve(m, y, np.array([PRE, PRE + 8]), ctx)
        out[m] = 1e6 * (time.perf_counter() - t0) / n
    return out


def _plots(out, A, B, C, singles, amps, gaps):
    ds = [int(d) for d in SEPARATIONS]
    # A: bias and spread vs separation, per pair type (second pulse is the hard one)
    fig, ax = P.figure(2, 3, width=15, height_per_row=3.3)
    for j, pname in enumerate(PAIRS_V):
        for k, m in enumerate(METHODS):
            for row, key in ((0, "bias"), (1, "spread")):
                vals = [100 * A[pname][str(d)]["methods"][m]["second"][key] for d in ds]
                ax[row, j].plot(ds, vals, marker="o", markersize=3.5, linewidth=1.2, color=P.SERIES[k], label=LABELS[m])
        ax[0, j].set_title(f"A, second pulse bias: {pname}")
        ax[1, j].set_title(f"A, second pulse spread: {pname}")
        ax[0, j].axhline(0, color=P.AXIS, linewidth=0.8)
        ax[0, j].set_ylim(-25, 25)
        ax[1, j].set_yscale("log")
        for r in (0, 1):
            ax[r, j].set_xlabel("separation (samples)")
            ax[r, j].set_ylabel("% of true amplitude")
    ax[0, 0].legend()
    P.save(fig, out / "pileup_template_pairs.png")

    # B: real pairs
    fig, ax = P.figure(2, 2, height_per_row=3.3)
    for k, m in enumerate(METHODS):
        for j, which in enumerate(("first", "second")):
            ax[0, j].plot(ds, [100 * B[str(d)][m][which]["bias"] for d in ds], marker="o", markersize=3.5,
                          linewidth=1.2, color=P.SERIES[k], label=LABELS[m])
            ax[1, j].plot(ds, [100 * B[str(d)][m][which]["spread"] for d in ds], marker="o", markersize=3.5,
                          linewidth=1.2, color=P.SERIES[k], label=LABELS[m])
            ax[0, j].set_title(f"B, real pulse pairs: {which} pulse bias")
            ax[1, j].set_title(f"B, real pulse pairs: {which} pulse spread")
    for a in ax.ravel():
        a.set_xlabel("separation (samples)")
        a.set_ylabel("%")
    ax[0, 0].axhline(0, color=P.AXIS, linewidth=0.8)
    ax[0, 1].axhline(0, color=P.AXIS, linewidth=0.8)
    for a in ax[0]:
        a.set_ylim(-15, 15)
    for a in ax[1]:
        a.set_yscale("log")
    ax[0, 0].legend()
    P.save(fig, out / "pileup_real_pairs.png")

    # C: spectra of resolved pile-up events per method, scaled to the method's isolated peaks
    fig, ax = P.figure(1, 1, height_per_row=4.0)
    edges = np.linspace(0.2, 0.45, 220)
    for k, m in enumerate(METHODS):
        try:
            c_iso = np.array(C[m]["isolated"]["centers"])
            scale = float(np.mean(np.array([0.3157, 0.3571]) / c_iso))
        except Exception:  # noqa: BLE001
            scale = 1.0
        h, _ = np.histogram(amps[m] * scale, edges)
        ax[0, 0].stairs(h, edges, color=P.SERIES[k], linewidth=1.1, label=LABELS[m])
    ax[0, 0].set_title("C: real ⁶⁰Co pile-up events resolved by each method (scaled to its isolated-pulse peaks)")
    ax[0, 0].set_xlabel("amplitude (V, scaled)")
    ax[0, 0].set_ylabel("counts per bin")
    ax[0, 0].legend()
    P.save(fig, out / "pileup_real_events.png")
