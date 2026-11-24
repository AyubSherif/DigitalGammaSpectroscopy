"""Evaluate detection methods with synthetic pulses injected into real noise.

Noise sources (both from the pre-processed signal of real captures):

* **real noise**: pulse-free stretches of the clean signal, joined together. It can still
  contain real pulses too small for the mask, so it is used only for injection tests.
* **pure noise**: the same stream with randomised Fourier phases. Same power spectrum,
  no pulses. Used for the false-trigger rate.

Tests, for each method:

* output noise sigma and false triggers per second (pure noise);
* efficiency vs amplitude, and extra detections per injected pulse (artifacts);
* pair resolution: two pulses at separation d (equal, big-then-small, small-then-big).

Run:  python -m gammaspec evaluate-detection STORE TEMPLATE OUT_DIR
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import plotting as P
from .baseline import design_from_template, filter_template
from .detect import DetectConfig, Detector
from .noise import pulse_free_windows
from .preprocess import PreprocessConfig, Preprocessor
from .signal import pulse_free_mask
from .store import CaptureStore
from .template import PulseTemplate

AMPLITUDES_V = (0.001, 0.0015, 0.002, 0.003, 0.004, 0.006, 0.008, 0.012, 0.02, 0.05, 0.3, 1.0, 2.0)
SEPARATIONS = (2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 45)
PAIRS_V = {"equal 300/300 mV": (0.3, 0.3), "big then small 300/30 mV": (0.3, 0.03),
           "small then big 30/300 mV": (0.03, 0.3)}


def default_methods() -> dict[str, DetectConfig]:
    base = DetectConfig()
    return {
        "optimal + cancellation (default)": base,
        "optimal, prominence only": replace(base, cancel=False),
        "matched": replace(base, method="matched", cancel=False),
        "sharpened": replace(base, method="sharpened", cancel=False, threshold=6.0),
        "edge (coarse trigger)": replace(base, method="edge"),
    }


def noise_streams(store_path, template: PulseTemplate, pz, captures, seed=0):
    """(real noise, pure noise) streams from the clean signal of the given captures."""
    rng = np.random.default_rng(seed)
    parts = []
    with CaptureStore(store_path) as s:
        fs = 1.0 / s.dt_s
        pre = Preprocessor(fs, PreprocessConfig(), pz)
        det = Detector(template, replace(DetectConfig(), threshold=3.0))
        for i in captures:
            pre.reset()
            y = pre.process(s.volts(i))[1500:]
            free = pulse_free_mask(y.size, det.detect(y).index, 60, 150)
            w, _ = pulse_free_windows(y, free, 256)
            parts.append(w.ravel())
    real = np.concatenate(parts)
    F = np.fft.rfft(real)
    pure = np.fft.irfft(np.abs(F) * np.exp(1j * rng.uniform(0, 2 * np.pi, F.size)), real.size)
    return real, pure, fs


def inject(noise, template, amps, rng, spacing=400, second=None):
    """Add pulses at random sub-sample times every ``spacing`` samples.
    ``second=(d, a2)`` adds a second pulse d samples after each one."""
    y = noise.copy()
    starts = np.arange(300, noise.size - 400, spacing)
    c = starts + rng.uniform(0, 1, starts.size)
    A = rng.choice(np.asarray(amps), starts.size)
    k = np.arange(-8, 250)
    for ci, a in zip(c, A):
        idx = int(np.floor(ci)) + k
        y[idx] += a * template(idx - ci)
        if second is not None:
            d, a2 = second
            idx2 = int(np.floor(ci + d)) + k
            idx2 = idx2[idx2 < y.size]
            y[idx2] += a2 * template(idx2 - ci - d)
    return y, c, A


def _matched(idx, positions, tol=2):
    """For each true position, is there a detection within tol samples?"""
    j = np.searchsorted(idx, positions)
    lo = np.abs(idx[np.clip(j - 1, 0, max(idx.size - 1, 0))] - positions) <= tol if idx.size else np.zeros(positions.size, bool)
    hi = np.abs(idx[np.clip(j, 0, max(idx.size - 1, 0))] - positions) <= tol if idx.size else np.zeros(positions.size, bool)
    return lo | hi


def _both_found(idx, p1, p2, tol=2):
    """Two *distinct* detections, one within tol of p1 and another within tol of p2."""
    out = np.zeros(p1.size, bool)
    for k, (a, b) in enumerate(zip(p1, p2)):
        near = idx[(idx >= a - tol) & (idx <= b + tol)]
        n1 = near[np.abs(near - a) <= tol]
        n2 = near[np.abs(near - b) <= tol]
        out[k] = n1.size > 0 and n2.size > 0 and np.unique(np.r_[n1, n2]).size >= 2
    return out


def evaluate(store_path, template_path, out_dir, captures=range(0, 100, 5), seed=0, verbose=True) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    T = PulseTemplate.load(template_path)
    pz = design_from_template(T)
    Tp = filter_template(T, pz.b, pz.a)
    real, pure, fs = noise_streams(store_path, Tp, pz, list(captures), seed)
    if verbose:
        print(f"noise streams: {real.size / fs:.2f} s, rms {real.std() * 1e3:.2f} mV")
    res = {}
    for name, cfg in default_methods().items():
        det = Detector(Tp, cfg)
        r = {}
        d = det.detect(pure)
        r["sigma_mv"] = d.sigma * 1e3
        r["threshold_mv"] = d.threshold_v * 1e3
        r["false_per_s"] = d.index.size / (pure.size / fs)
        y, c, A = inject(real, Tp, AMPLITUDES_V, rng)
        idx = det.detect(y).index
        pos = np.round(c).astype(np.int64)
        found = _matched(idx, pos)
        # extra detections: in (pos+3 .. pos+250) or (pos-60 .. pos-3), not matching any injection
        extra = np.array([np.sum(((idx > p + 2) & (idx < p + 250)) | ((idx < p - 2) & (idx > p - 60))) for p in pos])
        r["efficiency"] = {f"{a * 1e3:g}": float(found[A == a].mean()) for a in AMPLITUDES_V}
        r["extra_per_pulse"] = {f"{a * 1e3:g}": float(extra[A == a].mean()) for a in AMPLITUDES_V}
        r["pairs"] = {}
        for pname, (a1, a2) in PAIRS_V.items():
            row = {}
            for dsep in SEPARATIONS:
                y, c, _ = inject(real, Tp, [a1], rng, spacing=500, second=(dsep, a2))
                idx = det.detect(y).index
                p1 = np.round(c).astype(np.int64)
                p2 = np.round(c + dsep).astype(np.int64)
                row[str(dsep)] = float(np.mean(_both_found(idx, p1, p2)))
            r["pairs"][pname] = row
        res[name] = r
        if verbose:
            print(f"  {name:34s} sigma {r['sigma_mv']:.3f} mV, false {r['false_per_s']:.1f}/s, "
                  f"eff@2mV {r['efficiency']['2']:.2f}, @4mV {r['efficiency']['4']:.2f}, "
                  f"extra@1V {r['extra_per_pulse']['1000']:.3f}")
    (out / "detection_summary.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    _plot(res, out)
    return res


def _plot(res, out):
    names = list(res)
    fig, ax = P.figure(1, 2, height_per_row=3.6)
    amps = [float(a) for a in res[names[0]]["efficiency"]]
    for k, n in enumerate(names):
        ax[0, 0].plot(amps, list(res[n]["efficiency"].values()), color=P.SERIES[k], marker="o",
                      markersize=4, linewidth=1.3, label=n)
        ax[0, 1].plot(amps, list(res[n]["extra_per_pulse"].values()), color=P.SERIES[k], marker="o",
                      markersize=4, linewidth=1.3, label=n)
    for a in ax[0]:
        a.set_xscale("log")
        a.set_xlabel("injected pulse amplitude (mV)")
        a.legend()
    ax[0, 0].set_title("Detection efficiency (real noise)")
    ax[0, 1].set_title("Extra detections per injected pulse (artifacts)")
    ax[0, 1].set_yscale("symlog", linthresh=0.01)
    P.save(fig, out / "detection_efficiency.png")

    pairs = list(res[names[0]]["pairs"])
    fig, ax = P.figure(1, len(pairs), height_per_row=3.4)
    for j, pn in enumerate(pairs):
        for k, n in enumerate(names):
            row = res[n]["pairs"][pn]
            ax[0, j].plot([int(d) for d in row], list(row.values()), color=P.SERIES[k], marker="o",
                          markersize=4, linewidth=1.3, label=n)
        ax[0, j].set_title(f"Both pulses found: {pn}")
        ax[0, j].set_xlabel("separation (samples)")
        ax[0, j].set_ylim(-0.05, 1.05)
    ax[0, 0].legend()
    P.save(fig, out / "detection_pairs.png")
