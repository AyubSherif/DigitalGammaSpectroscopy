"""Evaluate waveform pre-processing (line removal, pole-zero, baseline tracker) on real data.

Every configuration is applied to the same captures, and the *same pulses* (selected
on the raw signal) are measured in each, so differences come from the processing only.

Metrics per configuration:
* noise sigma of the processed signal;
* pre-pulse baseline: level just before each pulse (should be 0 and narrow);
* tail: average normalised waveform after isolated pulses (should return to 0);
* spectrum: template-fit amplitudes, Gaussian fits of the main peaks (FWHM %).

Run:  python -m gamma_spectrometry evaluate-restoration STORE TEMPLATE OUT_DIR
"""

from __future__ import annotations

import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks as _peaks

from . import plotting as P
from .baseline import TrackerConfig, design_from_template, filter_template
from .preprocess import PreprocessConfig, Preprocessor
from .signal import find_pulses, isolated, noise_sigma, to_positive
from .spectrum import fit_peaks
from .store import CaptureStore
from .template import PulseTemplate, extract_windows, fit_pulses


def default_configs() -> dict[str, PreprocessConfig]:
    off = PreprocessConfig(line_removal=False, pole_zero=False, baseline_tracker=False)
    return {
        "raw": off,
        "lines": replace(off, line_removal=True),
        "tracker only": replace(off, baseline_tracker=True, tracker=TrackerConfig(window=151)),
        "pole-zero": replace(off, pole_zero=True),
        "lines + pole-zero": replace(off, line_removal=True, pole_zero=True),
        "lines + pole-zero + tracker": replace(off, line_removal=True, pole_zero=True, baseline_tracker=True),
    }


def _auto_peaks(a: np.ndarray, n: int, lo: float = 0.1, bin_width: float = 0.002):
    """Starting guesses for the n most prominent spectral peaks above ``lo`` volts."""
    edges = np.arange(lo, np.quantile(a, 0.999), bin_width)
    h, _ = np.histogram(a, edges)
    hs = gaussian_filter1d(h.astype(float), 2)
    idx, prop = _peaks(hs, prominence=0.05 * hs.max())
    idx = idx[np.argsort(prop["prominences"])[::-1][:n]]
    return np.sort(0.5 * (edges[idx] + edges[idx + 1]))


def evaluate(store_path, template_path, out_dir, configs: dict | None = None, n_peaks: int = 2,
             limit: int | None = None, verbose: bool = True) -> dict:
    configs = configs or default_configs()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tmpl = PulseTemplate.load(template_path)
    pz = design_from_template(tmpl)
    tmpl_pz = filter_template(tmpl, pz.b, pz.a)
    t0 = time.perf_counter()

    pre, post = 12, 50                       # spectrum windows
    tail_pre, tail_post = 20, 600            # tail windows (isolated pulses)
    acc = {name: dict(win=[], base=[], sigma=[], tail=[]) for name in configs}
    with CaptureStore(store_path) as store:
        fs = 1.0 / store.dt_s
        procs = {name: Preprocessor(fs, cfg, pz if cfg.pole_zero else None) for name, cfg in configs.items()}
        n_cap = len(store) if limit is None else min(limit, len(store))
        for i in range(n_cap):
            volts = store.volts(i)
            x0 = to_positive(volts)
            x0 = x0 - x0.mean()
            sig0 = noise_sigma(x0)
            peaks = find_pulses(x0, sig0, 5.0)
            n = x0.size
            keep = (peaks > 1500) & isolated(peaks, 25, post, n) & isolated(peaks, pre, post, n)
            sel = peaks[keep]
            amp0 = x0[sel]
            prev_gap = np.diff(peaks, prepend=-10**9)[keep]
            tsel = peaks[isolated(peaks, 300, tail_post, n) & (peaks > 1500) & (x0[peaks] > 30 * sig0)]
            for name, proc in procs.items():
                proc.reset()
                y = proc.process(volts)
                acc[name]["win"].append(extract_windows(y, sel, pre, post).astype(np.float32))
                if sel.size:  # mean of samples 25..8 before the peak, and the raw peak height
                    lvl = extract_windows(y, sel - 8, 17, 0).mean(axis=1)
                    acc[name]["base"].append(np.stack([lvl, amp0, prev_gap], 1))
                acc[name]["sigma"].append(proc.info.sigma_out_v)
                if tsel.size:
                    w = extract_windows(y, tsel, tail_pre, tail_post)
                    lvl = w[:, :tail_pre - 5].mean(axis=1, keepdims=True)
                    acc[name]["tail"].append(((w - lvl) / (y[tsel] - lvl[:, 0])[:, None]).astype(np.float32))
            if verbose and (i + 1) % max(1, n_cap // 10) == 0:
                print(f"  processed {i + 1}/{n_cap} captures", flush=True)

    results = {}
    spectra = {}
    for name, cfg in configs.items():
        a = acc[name]
        W = np.concatenate(a["win"]).astype(np.float64)
        base = np.concatenate(a["base"])
        T = tmpl_pz if cfg.pole_zero else tmpl
        restored = cfg.pole_zero or cfg.baseline_tracker
        fr = (-10, 40) if cfg.pole_zero else (-10, post)
        fit_free = fit_pulses(W, T, pre, fit_range=fr, fit_baseline=True)
        fit_zero = fit_pulses(W, T, pre, fit_range=fr, fit_baseline=False)
        tail = np.concatenate(a["tail"]).mean(axis=0) if a["tail"] else np.zeros(tail_pre + tail_post + 1)
        clear = base[:, 2] >= 60          # previous pulse's main lobe is over
        bl = base[clear, 0] * 1e3
        r = dict(
            sigma_mv=float(np.median(a["sigma"]) * 1e3),
            n_pulses=int(W.shape[0]),
            pre_pulse_baseline_mv=dict(median=float(np.median(bl)),
                                       robust_sigma=float(1.4826 * np.median(np.abs(bl - np.median(bl)))),
                                       p01=float(np.percentile(bl, 1)), p05=float(np.percentile(bl, 5)),
                                       p95=float(np.percentile(bl, 95)), p99=float(np.percentile(bl, 99)),
                                       frac_beyond_5mv=float(np.mean(np.abs(bl) > 5))),
            tail_at={str(k): float(tail[tail_pre + k]) for k in (40, 60, 100, 200, 400)},
            tail_rms_40_600=float(np.sqrt(np.mean(tail[tail_pre + 40:] ** 2))),
            restored=restored,
        )
        for label, f in (("baseline_free", fit_free), ("baseline_fixed", fit_zero)):
            amp = f.amplitude
            try:
                guess = _auto_peaks(amp, n_peaks)
                pf = fit_peaks(amp, guess, guess.min() * 0.85, guess.max() * 1.12, 0.002)
                r[label] = dict(centers_v=pf.centers.tolist(), fwhm_pct=(100 * pf.resolution).tolist(),
                                fwhm_pct_err=(100 * pf.resolution_err).tolist(), areas=pf.areas.tolist(),
                                chi2_ndf=pf.chi2_ndf)
            except Exception as exc:  # noqa: BLE001 - report, don't abort the whole evaluation
                r[label] = dict(error=str(exc))
        cats = dict(close=(base[:, 2] >= 25) & (base[:, 2] < 60),
                    on_tail=(base[:, 2] >= 60) & (base[:, 2] < 300),
                    clean=base[:, 2] >= 300)
        r["by_previous_gap"] = {}
        for cname, m in cats.items():
            row = dict(n=int(m.sum()))
            for label, f in (("baseline_free", fit_free), ("baseline_fixed", fit_zero)):
                try:
                    g = r[label]["centers_v"]
                    pf = fit_peaks(f.amplitude[m], g, min(g) * 0.85, max(g) * 1.12, 0.002)
                    row[label] = dict(centers_v=pf.centers.tolist(), fwhm_pct=(100 * pf.resolution).tolist())
                except Exception as exc:  # noqa: BLE001
                    row[label] = dict(error=str(exc))
            r["by_previous_gap"][cname] = row
        results[name] = r
        spectra[name] = (fit_free.amplitude, fit_zero.amplitude, base, tail)

    summary = dict(store=str(store_path), template=str(template_path), pole_zero=pz.describe(),
                   pz_b=pz.b.tolist(), pz_a=pz.a.tolist(), captures=n_cap, configs=results,
                   elapsed_s=time.perf_counter() - t0)
    (out / "restoration_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _plots(out, configs, spectra, tmpl, tmpl_pz, store_path, pz, tail_pre)
    if verbose:
        print(f"done in {summary['elapsed_s']:.1f} s -> {out}")
    return summary


# ----------------------------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------------------------

def _plots(out, configs, spectra, tmpl, tmpl_pz, store_path, pz, tail_pre):
    names = list(configs)
    # 1. template before/after P/Z
    fig, ax = P.figure(1, 2, height_per_row=3.4)
    t = tmpl.t
    for k, (T, lab) in enumerate([(tmpl, "measured pulse"), (tmpl_pz, "after pole-zero correction")]):
        ax[0, 0].plot(t[(t > -5) & (t < 40)], T(t[(t > -5) & (t < 40)]), color=P.SERIES[k], label=lab)
        ax[0, 1].plot(t, T(t), color=P.SERIES[k], label=lab)
    ax[0, 0].set_title("Pulse shape: main lobe")
    ax[0, 1].set_title("Pulse shape: tail (zoomed)")
    ax[0, 1].set_ylim(-0.05, 0.03)
    for a in ax[0]:
        a.set_xlabel("time from peak (samples)")
        a.legend()
    P.save(fig, out / "pole_zero_template.png")

    # 2. example waveform
    with CaptureStore(store_path) as store:
        volts = store.volts(min(1, len(store) - 1))
        fs = 1.0 / store.dt_s
    x = to_positive(volts)
    x = x - x.mean()
    best = names[-1]
    y = Preprocessor(fs, configs[best], pz).process(volts)
    sig = noise_sigma(x)
    pk = find_pulses(x, sig, 5.0)
    # a busy 3000-sample stretch
    starts = np.arange(5000, x.size - 3000, 500)
    counts = [np.sum((pk >= s) & (pk < s + 3000)) for s in starts]
    s0 = int(starts[int(np.argmax(counts))])
    seg = slice(s0, s0 + 3000)
    fig, ax = P.figure(2, 1, height_per_row=3.0)
    tt = np.arange(3000)
    ax[0, 0].plot(tt, x[seg] * 1e3, color=P.SERIES[0], linewidth=0.8, label="raw (polarity flipped, offset removed)")
    ax[0, 0].plot(tt, y[seg] * 1e3, color=P.SERIES[1], linewidth=0.8, label=best)
    ax[0, 0].set_ylabel("mV")
    ax[0, 0].set_title("Example waveform (busiest 3000 samples of capture 1)")
    ax[0, 0].legend()
    ax[1, 0].plot(tt, x[seg] * 1e3, color=P.SERIES[0], linewidth=0.8, label="raw")
    ax[1, 0].plot(tt, y[seg] * 1e3, color=P.SERIES[1], linewidth=0.8, label=best)
    ax[1, 0].set_ylim(-25, 25)
    ax[1, 0].set_ylabel("mV")
    ax[1, 0].set_title("Same, zoomed on the baseline")
    ax[1, 0].set_xlabel("sample")
    ax[1, 0].legend()
    P.save(fig, out / "waveform_example.png")

    # 2b. noise spectrum before/after line removal (same capture, pulse-free windows)
    from .noise import noise_spectrum, pulse_free_windows
    from .signal import pulse_free_mask
    lines_cfg = replace(configs[best], pole_zero=False, baseline_tracker=False, line_removal=True)
    yl = Preprocessor(fs, lines_cfg).process(volts)
    free = pulse_free_mask(x.size, pk, 20, 1500)
    fig, ax = P.figure(1, 1, height_per_row=3.6)
    for k, (sig_, lab) in enumerate([(x, "raw"), (yl, "after line removal")]):
        w, _ = pulse_free_windows(sig_, free, 512)
        if w.shape[0]:
            ns = noise_spectrum(w, fs)
            ax[0, 0].loglog(ns.freq_hz[1:] / 1e3, ns.psd[1:], color=P.SERIES[k], linewidth=1.1,
                            label=f"{lab}: σ = {ns.sigma_v * 1e3:.2f} mV")
    ax[0, 0].set_title("Noise power spectrum before and after line removal (capture 1)")
    ax[0, 0].set_xlabel("frequency (kHz)")
    ax[0, 0].set_ylabel("PSD (V²/Hz)")
    ax[0, 0].legend()
    P.save(fig, out / "line_removal_psd.png")

    # 3. pre-pulse baseline distributions and tails
    fig, ax = P.figure(1, 2, height_per_row=3.6)
    allb = np.concatenate([spectra[n][2][:, 0] for n in names]) * 1e3
    lo, hi = np.percentile(allb, [0.5, 99.5])
    edges = np.linspace(lo, hi, 200)
    for k, n in enumerate(names):
        ax[0, 0].stairs(np.histogram(spectra[n][2][:, 0] * 1e3, edges)[0], edges, color=P.SERIES[k],
                        linewidth=1.1, label=n)
        tail = spectra[n][3]
        tx = np.arange(tail.size) - tail_pre
        ax[0, 1].plot(tx[tx > 20], tail[tx > 20], color=P.SERIES[k], linewidth=1.1, label=n)
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_title("Baseline just before each pulse")
    ax[0, 0].set_xlabel("mV")
    ax[0, 0].legend()
    ax[0, 1].axhline(0, color=P.AXIS, linewidth=0.8)
    ax[0, 1].set_title("Average tail after isolated pulses")
    ax[0, 1].set_xlabel("samples after peak")
    ax[0, 1].set_ylabel("fraction of peak")
    ax[0, 1].legend()
    P.save(fig, out / "baseline_and_tail.png")

    # 4. spectra raw vs best (both baseline-free fits), zoom on the peaks
    fig, ax = P.figure(1, 1, height_per_row=4.0)
    a_raw = spectra[names[0]][0]
    a_best = spectra[best][1]
    hi = np.quantile(a_raw, 0.995) * 1.3
    edges = np.linspace(0, hi, 700)
    ax[0, 0].stairs(np.histogram(a_raw, edges)[0], edges, color=P.SERIES[0], linewidth=1.0,
                    label="raw, template fit with free baseline")
    ax[0, 0].stairs(np.histogram(a_best, edges)[0], edges, color=P.SERIES[1], linewidth=1.1,
                    label=f"{best}, baseline fixed at 0")
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_xlabel("pulse amplitude (V)")
    ax[0, 0].set_ylabel("counts per bin")
    ax[0, 0].set_title("Amplitude spectrum, same pulses (no pile-up correction yet)")
    ax[0, 0].legend()
    P.save(fig, out / "spectrum_comparison.png")
