"""Detector characterisation: pulse template, analytic pulse model, noise, stability, raw spectrum.

Run from the command line:  python -m gammaspec characterize STORE.h5 OUT_DIR
Outputs (in OUT_DIR): template.h5, noise.h5, summary.json and diagnostic PNGs.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import h5py
import numpy as np

from . import plotting as P
from .noise import NoiseSpectrum, noise_spectrum, pulse_free_windows
from .signal import (DEFAULT_POLARITY, find_pulses, isolated, noise_sigma, pulse_free_mask,
                     to_positive)
from .store import CaptureStore
from .template import (PulseFit, PulseTemplate, build_template, extract_windows, fit_multiexp,
                       fit_pulses, _recentre, _stack)


@dataclass
class CharacterizeConfig:
    polarity: int = DEFAULT_POLARITY
    trigger_threshold: float = 5.0       # coarse trigger, in sigma of the edge signal
    # template pulses
    template_pre: int = 40               # samples kept before the peak
    template_post: int = 1500            # samples kept after the peak (covers the undershoot)
    isolation_pre: int = 1500            # no other pulse this many samples before ...
    isolation_post: int = 1500           # ... or after a template pulse
    min_amplitude_sigma: float = 30.0    # template pulses must exceed this (peak sample, in noise sigma)
    max_amplitude_v: float = 4.0         # ... and stay below this (avoid saturation / non-linearity)
    resolution: int = 16                 # template grid points per sample
    # spectrum pulses (short isolation: only the main lobe must be clean)
    spectrum_isolation_pre: int = 50
    spectrum_isolation_post: int = 60
    # noise
    noise_window: int = 512
    noise_exclude_pre: int = 20
    noise_exclude_post: int = 1500
    max_noise_windows_per_capture: int = 100
    amplitude_bin_edges: tuple = (0.03, 0.1, 0.2, 0.3, 0.45, 0.7, 1.2, 4.0)  # volts, shape-vs-amplitude check
    limit: int | None = None             # only the first N captures


# ----------------------------------------------------------------------------------------------
# 1. Scan the store
# ----------------------------------------------------------------------------------------------

def scan(store: CaptureStore, cfg: CharacterizeConfig, verbose=True) -> dict:
    n_cap = len(store) if cfg.limit is None else min(cfg.limit, len(store))
    n = store.n_samples
    tw, sw, sw_naive, nw, nw_cap = [], [], [], [], []
    per = {k: [] for k in ("time_s", "baseline_v", "sigma_v", "n_pulses", "min_code", "n_at_min_code")}
    rng = np.random.default_rng(0)
    s_pre, s_post = 12, cfg.spectrum_isolation_post
    for i in range(n_cap):
        codes = store.codes(i)
        x = to_positive(codes * store.lsb_volts, cfg.polarity)
        base = float(np.median(x))
        x = x - base
        sig = noise_sigma(x)
        peaks = find_pulses(x, sig, cfg.trigger_threshold)

        per["time_s"].append(store.time_s[i])
        per["baseline_v"].append(base)
        per["sigma_v"].append(sig)
        per["n_pulses"].append(peaks.size)
        per["min_code"].append(int(codes.min()))
        per["n_at_min_code"].append(int(np.sum(codes == codes.min())))

        # template candidates: long isolation, sensible amplitude
        ok = isolated(peaks, cfg.isolation_pre, cfg.isolation_post, n)
        ok &= isolated(peaks, cfg.template_pre, cfg.template_post, n)
        amp = x[peaks]
        ok &= (amp > cfg.min_amplitude_sigma * sig) & (amp < cfg.max_amplitude_v)
        if ok.any():
            tw.append(extract_windows(x, peaks[ok], cfg.template_pre, cfg.template_post).astype(np.float32))

        # spectrum candidates: main lobe clean
        ok = isolated(peaks, cfg.spectrum_isolation_pre, cfg.spectrum_isolation_post, n)
        ok &= isolated(peaks, s_pre, s_post, n)
        if ok.any():
            sw.append(extract_windows(x, peaks[ok], s_pre, s_post).astype(np.float32))
            sw_naive.append(amp[ok].astype(np.float32))

        # noise windows
        free = pulse_free_mask(n, peaks, cfg.noise_exclude_pre, cfg.noise_exclude_post)
        win, _ = pulse_free_windows(x, free, cfg.noise_window)
        if win.shape[0] > cfg.max_noise_windows_per_capture:
            win = win[np.sort(rng.choice(win.shape[0], cfg.max_noise_windows_per_capture, replace=False))]
        nw.append(win.astype(np.float32))
        nw_cap.append(np.full(win.shape[0], i))
        if verbose and (i + 1) % max(1, n_cap // 10) == 0:
            print(f"  scanned {i + 1}/{n_cap} captures", flush=True)

    cat = lambda L, w: np.concatenate(L) if L else np.empty((0, w), np.float32)  # noqa: E731
    return dict(
        n_captures=n_cap,
        template_windows=cat(tw, cfg.template_pre + cfg.template_post + 1),
        spectrum_windows=cat(sw, s_pre + s_post + 1),
        spectrum_pre=s_pre,
        naive_amplitude=np.concatenate(sw_naive) if sw_naive else np.empty(0),
        noise_windows=cat(nw, cfg.noise_window),
        noise_capture=np.concatenate(nw_cap) if nw_cap else np.empty(0, int),
        per_capture={k: np.asarray(v) for k, v in per.items()},
    )


# ----------------------------------------------------------------------------------------------
# 2. Analyses
# ----------------------------------------------------------------------------------------------

def shape_by_amplitude(windows, pre, fit: PulseFit, resolution, edges, min_pulses=30):
    """Templates for amplitude bins, plus a reference stacked from all pulses with the
    same per-pulse fits, so differences reflect amplitude dependence only.

    The bin templates are NOT recentred: they share the reference's time axis, so a
    shape that peaks later for large pulses shows up as a difference.
    """
    def stack(m):
        sub = fit.subset(m)
        t, s = _stack(windows[m], pre, sub, resolution, sub.amplitude ** 2, 40.0)
        return PulseTemplate(resolution, float(t[0]), s)

    ref = stack(np.ones(fit.amplitude.size, dtype=bool))
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (fit.amplitude >= lo) & (fit.amplitude < hi)
        if m.sum() >= min_pulses:
            out.append((float(lo), float(hi), int(m.sum()), stack(m)))
    return out, ref


def characterize(store_path, out_dir, cfg: CharacterizeConfig | None = None,
                 compare: str | None = None, verbose=True) -> dict:
    cfg = cfg or CharacterizeConfig()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    with CaptureStore(store_path) as store:
        fs = 1.0 / store.dt_s
        capture_s = store.n_samples * store.dt_s
        if verbose:
            print(f"scanning {store_path}")
        sc = scan(store, cfg, verbose)
    pc = sc["per_capture"]
    live = sc["n_captures"] * capture_s

    # --- template -------------------------------------------------------------------------
    tw = sc["template_windows"].astype(np.float64)
    if verbose:
        print(f"building template from {tw.shape[0]} isolated pulses")
    if tw.shape[0] < 50:
        raise RuntimeError(f"only {tw.shape[0]} isolated pulses: relax the isolation settings")
    tmpl, fit, used = build_template(tw, cfg.template_pre, cfg.resolution, verbose=verbose)
    model3 = fit_multiexp(tmpl)
    model4 = fit_multiexp(tmpl, taus0=tuple(model3.taus) + (5 * model3.taus[-1],),
                          amps0=tuple(model3.amps[1:]) + (0.005,))
    model = model4 if model4.rms < 0.7 * model3.rms else model3
    tmpl.meta.update(source_store=str(store_path), polarity=cfg.polarity,
                     model_t0=model.t0, model_sigma=model.sigma,
                     model_taus=model.taus, model_amps=model.amps)
    tmpl.save(out / "template.h5")
    bins, bins_ref = shape_by_amplitude(tw[used], cfg.template_pre, fit.subset(used), cfg.resolution,
                                        np.asarray(cfg.amplitude_bin_edges))

    # --- noise ----------------------------------------------------------------------------
    ns = noise_spectrum(sc["noise_windows"], fs)
    sig = pc["sigma_v"]
    p10, p90 = np.percentile(sig, [10, 90])
    split = p90 / p10 > 1.5
    noise_thr = float(np.sqrt(p10 * p90)) if split else float("inf")
    quiet_cap = sig < noise_thr
    ns_quiet = noise_spectrum(sc["noise_windows"][quiet_cap[sc["noise_capture"]]], fs) if split else None
    ns_noisy = noise_spectrum(sc["noise_windows"][~quiet_cap[sc["noise_capture"]]], fs) if split else None
    rate = pc["n_pulses"] / capture_s
    with h5py.File(out / "noise.h5", "w") as f:
        for k in ("freq_hz", "psd", "autocorr"):
            f.create_dataset(k, data=getattr(ns, k))
        f.attrs.update(sigma_v=ns.sigma_v, n_windows=ns.n_windows, window_length=ns.window_length,
                       white_level_v2_per_hz=ns.white_level)

    # --- preliminary spectrum ---------------------------------------------------------------
    sw = sc["spectrum_windows"].astype(np.float64)
    sfit = fit_pulses(sw, tmpl, sc["spectrum_pre"], fit_range=(-10, 60), steps_per_sample=cfg.resolution)
    good = sfit.rms < 4 * np.median(sfit.rms) + 3 * ns.sigma_v

    # --- summary --------------------------------------------------------------------------
    core = (tmpl.t >= -5) & (tmpl.t <= 60)
    bin_dev = [float(np.max(np.abs(b[3](tmpl.t[core]) - bins_ref(tmpl.t[core])))) for b in bins]
    under_min = float(tmpl.shape.min())
    summary = dict(
        store=str(store_path),
        captures=sc["n_captures"],
        live_time_s=live,
        count_rate_cps=float(np.sum(pc["n_pulses"]) / live),
        baseline_v=dict(mean=float(pc["baseline_v"].mean()), std=float(pc["baseline_v"].std())),
        noise=dict(sigma_v_diff=float(np.median(pc["sigma_v"])), sigma_v_windows=ns.sigma_v,
                   white_psd_v2_per_hz=ns.white_level, n_windows=ns.n_windows,
                   autocorr_lag1=float(ns.autocorr[1]), autocorr_lag2=float(ns.autocorr[2]),
                   intermittent=bool(split), sigma_threshold_v=noise_thr,
                   fraction_noisy_captures=float(np.mean(~quiet_cap)),
                   sigma_v_quiet=float(np.median(sig[quiet_cap])) if quiet_cap.any() else None,
                   sigma_v_noisy=float(np.median(sig[~quiet_cap])) if (~quiet_cap).any() else None,
                   rate_cps_quiet=float(np.mean(rate[quiet_cap])) if quiet_cap.any() else None,
                   rate_cps_noisy=float(np.mean(rate[~quiet_cap])) if (~quiet_cap).any() else None,
                   lines_khz=strongest_lines(ns_noisy or ns, 8)),
        clipping=dict(min_code=int(pc["min_code"].min()),
                      max_samples_at_min_code=int(pc["n_at_min_code"].max())),
        template=dict(n_pulses_isolated=int(tw.shape[0]), n_pulses_used=int(used.sum()),
                      resolution=cfg.resolution, window=[cfg.template_pre, cfg.template_post],
                      undershoot_min=under_min,
                      undershoot_min_t=float(tmpl.t[np.argmin(tmpl.shape)]),
                      area_samples=float(np.sum(tmpl.shape) / cfg.resolution),
                      area_positive_lobe=float(np.sum(np.clip(tmpl.shape, 0, None)) / cfg.resolution),
                      iterations=int(tmpl.meta["iterations"]),
                      shift_std=float(np.std(fit.shift[used]))),
        model=dict(components=len(model.taus), t0=model.t0, sigma=model.sigma,
                   taus_samples=model.taus.tolist(), taus_us=(model.taus * store_dt_us(fs)).tolist(),
                   amps=model.amps.tolist(), rms=model.rms, max_abs=model.max_abs,
                   rms_3comp=model3.rms, rms_4comp=model4.rms),
        shape_vs_amplitude=[dict(amp_lo=b[0], amp_hi=b[1], n=b[2], max_dev_core=d)
                            for b, d in zip(bins, bin_dev)],
        spectrum=dict(n_isolated=int(sw.shape[0]), n_good_fit=int(good.sum())),
        elapsed_s=time.perf_counter() - t0,
        config=asdict(cfg),
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    # --- plots ----------------------------------------------------------------------------
    ref = PulseTemplate.load(compare) if compare else None
    plot_template(tmpl, model, out / "template.png", ref)
    plot_shape_vs_amplitude(bins_ref, bins, out / "shape_vs_amplitude.png")
    plot_noise(ns, out / "noise.png", ns_quiet, ns_noisy)
    plot_stability(pc, out / "stability.png")
    plot_spectrum(sc["naive_amplitude"], sfit.amplitude[good], out / "spectrum_raw.png")
    if verbose:
        print(f"done in {summary['elapsed_s']:.1f} s -> {out}")
    return summary


def store_dt_us(fs):
    return 1e6 / fs


# ----------------------------------------------------------------------------------------------
# 3. Plots
# ----------------------------------------------------------------------------------------------

def plot_template(tmpl: PulseTemplate, model, path, ref: PulseTemplate | None = None):
    fig, ax = P.figure(3, 1, height_per_row=3.0)
    t = tmpl.t
    core = (t >= -6) & (t <= 40)
    a = ax[0, 0]
    a.plot(t[core], tmpl.shape[core], label="template (data)", color=P.SERIES[0])
    a.plot(t[core], model(t[core]), label=f"{len(model.taus)}-exponential model", color=P.SERIES[1],
           linestyle="--", linewidth=1.2)
    if ref is not None:
        a.plot(t[core], ref(t[core]), label="reference template", color=P.SERIES[2], linewidth=1.2)
    a.set_title("Pulse template: main lobe (peak normalised to 1)")
    a.set_xlabel("time from peak (samples, 1 sample = 0.556 µs)")
    a.legend()
    b = ax[1, 0]
    b.plot(t, tmpl.shape, color=P.SERIES[0], label="template (data)")
    b.plot(t, model(t), color=P.SERIES[1], linestyle="--", linewidth=1.2, label="model")
    if ref is not None:
        b.plot(t, ref(t), color=P.SERIES[2], linewidth=1.2, label="reference template")
    b.set_ylim(1.6 * min(tmpl.shape.min(), -0.01), 0.04)
    b.set_title("Undershoot (zoomed)")
    b.set_xlabel("time from peak (samples)")
    b.legend()
    c = ax[2, 0]
    c.plot(t, tmpl.shape - model(t), color=P.SERIES[0], linewidth=1.0, label="template − model")
    if ref is not None:
        c.plot(t, tmpl.shape - ref(t), color=P.SERIES[2], linewidth=1.0, label="template − reference")
    c.set_xscale("symlog", linthresh=20)
    c.set_title("Residuals (fraction of peak)")
    c.set_xlabel("time from peak (samples, symmetric-log scale)")
    c.legend()
    P.save(fig, path)


def plot_shape_vs_amplitude(tmpl, bins, path):
    fig, ax = P.figure(1, 2, height_per_row=3.6)
    t = tmpl.t
    core = (t >= -4) & (t <= 40)
    for k, (lo, hi, n, b) in enumerate(bins):
        lab = f"{lo * 1e3:.0f}–{hi * 1e3:.0f} mV (n={n})"
        ax[0, 0].plot(t[core], b(t[core]) - tmpl.shape[core], color=P.SERIES[k], linewidth=1.0, label=lab)
        tail = t > 20
        ax[0, 1].plot(t[tail], b(t[tail]) - tmpl.shape[tail], color=P.SERIES[k], linewidth=1.0, label=lab)
    ax[0, 0].set_title("Shape difference vs amplitude: main lobe")
    ax[0, 0].set_xlabel("time from peak (samples)")
    ax[0, 0].set_ylabel("bin template − overall (fraction of peak)")
    ax[0, 0].legend()
    ax[0, 1].set_title("Shape difference vs amplitude: tail")
    ax[0, 1].set_xlabel("time from peak (samples)")
    ax[0, 1].legend()
    P.save(fig, path)


def strongest_lines(ns: NoiseSpectrum, k: int) -> list:
    """Frequencies (kHz) of the k strongest narrow lines above the white floor."""
    ratio = ns.psd / ns.white_level
    loc = np.flatnonzero((ratio[1:-1] > ratio[:-2]) & (ratio[1:-1] >= ratio[2:]) & (ratio[1:-1] > 5)) + 1
    loc = loc[np.argsort(ratio[loc])[::-1][:k]]
    return [dict(khz=round(float(ns.freq_hz[i]) / 1e3, 1), above_floor=round(float(ratio[i]), 1))
            for i in sorted(loc)]


def plot_noise(ns: NoiseSpectrum, path, quiet=None, noisy=None):
    fig, ax = P.figure(1, 2, height_per_row=3.4)
    a = ax[0, 0]
    b = ax[0, 1]
    lags = np.arange(min(60, ns.autocorr.size))
    sets = [(ns, "all captures")] if quiet is None else [(quiet, "quiet captures"), (noisy, "noisy captures")]
    for k, (s, lab) in enumerate(sets):
        lab = f"{lab}: σ = {s.sigma_v * 1e3:.2f} mV ({s.n_windows} windows)"
        a.loglog(s.freq_hz[1:] / 1e3, s.psd[1:], color=P.SERIES[k], linewidth=1.1, label=lab)
        b.plot(lags, s.autocorr[lags], color=P.SERIES[k], marker="o", markersize=3, linewidth=1.0, label=lab)
    a.set_title("Noise power spectrum")
    a.set_xlabel("frequency (kHz)")
    a.set_ylabel("PSD (V²/Hz)")
    a.legend()
    b.legend()
    b.axhline(0, color=P.AXIS, linewidth=0.8)
    b.set_title("Noise autocorrelation")
    b.set_xlabel("lag (samples)")
    P.save(fig, path)


def plot_stability(pc, path):
    fig, ax = P.figure(3, 1, height_per_row=2.4)
    t = (pc["time_s"] - pc["time_s"][0]) / 60
    ax[0, 0].plot(t, pc["baseline_v"] * 1e3, color=P.SERIES[0], marker="o", markersize=2.5, linewidth=0.8)
    ax[0, 0].set_title("Baseline (capture median, positive-pulse polarity)")
    ax[0, 0].set_ylabel("mV")
    ax[1, 0].plot(t, pc["sigma_v"] * 1e3, color=P.SERIES[0], marker="o", markersize=2.5, linewidth=0.8)
    ax[1, 0].set_title("Noise σ per capture")
    ax[1, 0].set_ylabel("mV")
    dur = 524288 / 1.8e6
    ax[2, 0].plot(t, pc["n_pulses"] / dur, color=P.SERIES[0], marker="o", markersize=2.5, linewidth=0.8)
    ax[2, 0].set_title("Trigger rate per capture")
    ax[2, 0].set_ylabel("counts/s")
    ax[2, 0].set_xlabel("time since first capture (min)")
    P.save(fig, path)


def plot_spectrum(naive, fitted, path):
    fig, ax = P.figure(1, 1, height_per_row=4.0)
    a = ax[0, 0]
    hi = np.quantile(fitted, 0.995) * 1.3 if fitted.size else 1.0
    edges = np.linspace(0, hi, 600)
    a.stairs(np.histogram(naive, edges)[0], edges, color=P.SERIES[1], linewidth=1.0,
             label="peak sample − capture median (naive)")
    a.stairs(np.histogram(fitted, edges)[0], edges, color=P.SERIES[0], linewidth=1.2,
             label="template fit (isolated pulses only)")
    a.set_yscale("log")
    a.set_title("Preliminary amplitude spectrum (no pile-up correction, no baseline restoration)")
    a.set_xlabel("pulse amplitude (V)")
    a.set_ylabel("counts per bin")
    a.legend()
    P.save(fig, path)
