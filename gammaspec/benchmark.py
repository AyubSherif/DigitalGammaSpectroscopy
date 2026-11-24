"""Benchmark: the gammaspec pipeline vs a naive peak search vs an (emulated) analog chain.

All three chains see the **same recorded pulses** and are calibrated by the same
automatic procedure (``calibration.calibrate``), so differences come from the processing.

* **naive**: ``scipy.signal.find_peaks`` on the raw waveform (polarity flipped, capture
  median subtracted, nothing else): height > ``k`` x noise sigma, peaks at least
  ``distance`` samples apart. The pulse height is the peak sample.
* **analog** (emulated): the textbook NaI(Tl) chain applied digitally to the waveform:
  shaping amplifier = pole-zero-cancelled CR differentiator + 4 RC integrators
  (CR-RC^4, shaping time tau), baseline restorer, peak-sensing ADC (peak of the shaped
  pulse, parabolic interpolation like a peak-hold), lower-level discriminator, pile-up
  rejector (optional) and a non-extending ADC dead time after each conversion, with an
  MCA-style live-time clock.
* **gammaspec**: the event file of the full pipeline (good events).

Metrics: FWHM of the 1173 / 1332 keV peaks (% and keV), valley-to-peak ratio between them,
full-energy-peak rates (counts per live second), pile-up continuum 1450-2450 keV, total
rate, low-energy threshold, dead-time / rejected fraction.

Run:  gammaspec benchmark STORE EVENTS OUT_DIR [--first N0 --n N]
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, lfilter

from . import plotting as P
from .analysis import NetSpectrum, energy_spectrum, fit_line_groups
from .baseline import TrackerConfig, track_baseline
from .calibration import calibrate
from .detect import negative_side_sigma
from .events import EventFile
from .signal import noise_sigma, pulse_free_mask, to_positive
from .spectrum import fit_peaks
from .store import CaptureStore


@dataclass
class ChainResult:
    name: str
    amplitudes: np.ndarray
    live_time_s: float
    real_time_s: float
    n_triggers: int            # pulses seen by the chain (before rejection / dead time)
    n_recorded: int            # pulses that made it into the spectrum
    note: str = ""


# ----------------------------------------------------------------------------------------------
# naive
# ----------------------------------------------------------------------------------------------

def naive_chain(store_path, captures, k: float = 5.0, distance: int = 5, polarity: int = -1) -> ChainResult:
    amps, live, n = [], 0.0, 0
    with CaptureStore(store_path) as s:
        for i in captures:
            x = to_positive(s.volts(i), polarity)
            x = x - np.median(x)
            p, prop = find_peaks(x, height=k * noise_sigma(x), distance=distance)
            amps.append(prop["peak_heights"])
            live += x.size * s.dt_s
            n += p.size
    a = np.concatenate(amps)
    return ChainResult("naive find_peaks", a, live, live, n, n)


# ----------------------------------------------------------------------------------------------
# emulated analog chain
# ----------------------------------------------------------------------------------------------

def shaping_amplifier(x: np.ndarray, dt: float, tau_shaping: float = 1e-6, tau_pz: float = 4.1e-6,
                      n_integrators: int = 4) -> np.ndarray:
    """CR-RC^n shaping with pole-zero cancellation of the input decay tau_pz.

    PZ-CR: (1 - a_pz z^-1) / (1 - a_s z^-1) turns an exponential with time constant tau_pz
    into one with tau_shaping; each RC integrator is (1 - a_s) / (1 - a_s z^-1). The result
    is a semi-Gaussian pulse peaking at about n * tau_shaping.
    """
    a_s = np.exp(-dt / tau_shaping)
    a_pz = np.exp(-dt / tau_pz)
    y = lfilter([1.0, -a_pz], [1.0, -a_s], x)
    for _ in range(n_integrators):
        y = lfilter([1.0 - a_s], [1.0, -a_s], y)
    return y


def analog_chain(store_path, captures, tau_shaping: float = 1e-6, tau_pz: float = 4.1e-6,
                 lld_sigma: float = 5.0, dead_time: float = 5e-6, pur: bool = True,
                 polarity: int = -1, warmup: int = 1500) -> ChainResult:
    amps, live, real, n_trig, n_rec = [], 0.0, 0.0, 0, 0
    with CaptureStore(store_path) as s:
        dt = s.dt_s
        for i in captures:
            x = to_positive(s.volts(i), polarity)
            x = x - x.mean()
            y = shaping_amplifier(x, dt, tau_shaping, tau_pz)[warmup:]
            # the shaped signal is smooth, so sample-to-sample differences underestimate its
            # noise: measure it on the output itself (negative half, robust to pulses)
            sig = negative_side_sigma(y - np.median(y))
            pk, _ = find_peaks(y, height=lld_sigma * max(sig, 1e-9), prominence=lld_sigma * max(sig, 1e-9))
            # baseline restorer: slow average of the samples away from pulses
            peaking = int(round(4 * tau_shaping / dt))
            width = int(round(10 * tau_shaping / dt))         # semi-Gaussian back to ~baseline
            free = pulse_free_mask(y.size, pk, peaking + 2, width)
            y = y - track_baseline(y, free, TrackerConfig(window=401, min_free=20))
            pk = pk[(pk > 1) & (pk < y.size - 2)]
            # peak-hold: parabolic interpolation of the maximum
            ym, y0, yp = y[pk - 1], y[pk], y[pk + 1]
            den = ym - 2 * y0 + yp
            off = np.where(den < 0, 0.5 * (ym - yp) / np.where(den < 0, den, -1), 0.0)
            amp = y0 - 0.25 * (ym - yp) * off
            n_trig += pk.size
            # pile-up rejection: reject pulses with another trigger within the pulse width
            keep = np.ones(pk.size, bool)
            if pur and pk.size > 1:
                g = np.diff(pk)
                close = g < width
                keep[1:] &= ~close
                keep[:-1] &= ~close
            # ADC dead time (non-extending), counted from each accepted conversion
            nd = int(round(dead_time / dt))
            last = -10 ** 9
            dead_samples = 0
            accepted = []
            for j in np.flatnonzero(keep):
                if pk[j] - last < nd:
                    continue
                accepted.append(j)
                last = pk[j]
                dead_samples += nd
            accepted = np.array(accepted, dtype=int)
            amps.append(amp[accepted])
            n_rec += accepted.size
            # MCA live-time clock: stops during conversions; rejected pile-up time is also
            # lost time (the rejector's live-time correction)
            rejected = int(np.sum(~keep))
            lost = dead_samples + rejected * width / 2
            live += max(y.size - lost, 0) * dt
            real += y.size * dt
    name = f"analog CR-RC4 {tau_shaping * 1e6:g} us, " + ("with" if pur else "no") + " pile-up rejection"
    return ChainResult(name, np.concatenate(amps), live, real, n_trig, n_rec)


# ----------------------------------------------------------------------------------------------
# gammaspec
# ----------------------------------------------------------------------------------------------

def pipeline_chain(event_path, captures=None) -> ChainResult:
    ef = EventFile(event_path)
    ev = ef.good()
    caps = ef.captures
    if captures is not None:
        sel = np.isin(ev["capture"], list(captures))
        ev = ev[sel]
        caps = caps[np.isin(caps["capture"], list(captures))]
    valid = ef.events if captures is None else ef.events[np.isin(ef.events["capture"], list(captures))]
    live = float(caps["live_time_s"].sum())
    return ChainResult("gammaspec pipeline", ev["amplitude"], live, live, int(valid.size), int(ev.size))


# ----------------------------------------------------------------------------------------------
# metrics
# ----------------------------------------------------------------------------------------------

def metrics(chain: ChainResult, bin_kev: float = 2.0) -> dict:
    a = chain.amplitudes
    cal = calibrate(a, "source")
    e = cal.energy(a)
    out = dict(name=chain.name, live_time_s=chain.live_time_s, real_time_s=chain.real_time_s,
               triggers=chain.n_triggers, recorded=chain.n_recorded,
               recorded_fraction=chain.n_recorded / max(chain.n_triggers, 1),
               live_fraction=chain.live_time_s / max(chain.real_time_s, 1e-12),
               calibration=list(map(float, cal.coeffs)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pf = fit_peaks(e, [1173.2, 1332.5], 1040, 1440, bin_kev, sigma0=[32.4, 33.6], center_tol=0.03,
                       sigma_range=(0.4, 2.5))
    out["fwhm_pct"] = (100 * pf.resolution).tolist()
    out["fwhm_kev"] = pf.fwhm.tolist()
    out["fwhm_pct_err"] = (100 * pf.resolution_err).tolist()
    out["centres_kev"] = pf.centers.tolist()
    # valley-to-peak: histogram minimum between the peaks / mean of the two peak maxima (5 keV smoothing)
    edges = np.arange(1100, 1400.1, 5.0)
    h, _ = np.histogram(e, edges)
    c = 0.5 * (edges[1:] + edges[:-1])
    p1 = h[(c > 1150) & (c < 1200)].max()
    p2 = h[(c > 1310) & (c < 1360)].max()
    v = h[(c > 1230) & (c < 1280)].min()
    out["valley_to_peak"] = float(v / (0.5 * (p1 + p2)))
    S = energy_spectrum(e, chain.live_time_s, np.arange(0, 3000.1, 5.0))
    lines = fit_line_groups(NetSpectrum(S.edges_kev, S.rate, S.rate_err))
    out["peak_rates_cps"] = {f"{r['label']} {r['energy_kev']:.0f}": [r["rate_cps"], r["rate_err_cps"]] for r in lines}
    pu = (e > 1450) & (e < 2450)
    out["pileup_continuum_cps"] = float(pu.sum() / chain.live_time_s)
    out["total_rate_cps_above_120kev"] = float(np.sum(e > 120) / chain.live_time_s)
    # threshold: 1st percentile of the energies of the recorded pulses below 300 keV
    low = e[(e > 0) & (e < 300)]
    out["low_energy_threshold_kev"] = float(np.percentile(low, 1)) if low.size else float("nan")
    out["_energies"] = e
    return out


# Typical values for a 3"x3" NaI(Tl) with a conventional analog chain (for context; another
# detector): FWHM about 7 % at 662 keV (6.5-8 %), about 4.8-5.5 % at 1332 keV; shaping 0.5-1 us.
LITERATURE = dict(fwhm_662_pct="6.5-8 (typically ~7)", fwhm_1332_pct="~4.8-5.5", shaping_us="0.5-1")


def run_benchmark(store_path, event_path, out_dir, captures=None, verbose=True) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with CaptureStore(store_path) as s:
        caps = list(range(len(s))) if captures is None else list(captures)
    chains = [
        naive_chain(store_path, caps),
        analog_chain(store_path, caps, tau_shaping=1e-6, pur=False),
        analog_chain(store_path, caps, tau_shaping=1e-6, pur=True),
        analog_chain(store_path, caps, tau_shaping=0.5e-6, pur=True),
        pipeline_chain(event_path, caps),
    ]
    res = []
    for ch in chains:
        if verbose:
            print(f"  {ch.name}: {ch.n_recorded} recorded of {ch.n_triggers} triggers", flush=True)
        res.append(metrics(ch))
    spectra = {r["name"]: r.pop("_energies") for r in res}
    summary = dict(store=str(store_path), events=str(event_path), captures=len(caps), chains=res,
                   literature_3x3_nai_analog=LITERATURE)
    (out / "benchmark_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _plot(out, res, spectra)
    return summary


def _plot(out, res, spectra):
    fig, ax = P.figure(2, 1, height_per_row=3.6)
    edges = np.arange(0, 2800.1, 5.0)
    for k, r in enumerate(res):
        e = spectra[r["name"]]
        h, _ = np.histogram(e, edges)
        ax[0, 0].stairs(h / r["live_time_s"] / 5.0, edges, color=P.SERIES[k], linewidth=1.0, label=r["name"])
        z = (edges[:-1] > 950) & (edges[:-1] < 1550)
        ez = edges[np.r_[z, False] | np.r_[False, z]]
        ax[1, 0].stairs((h / r["live_time_s"] / 5.0)[z], ez, color=P.SERIES[k], linewidth=1.2, label=r["name"])
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_title("Same pulses, three processing chains (each calibrated on its own ⁶⁰Co lines)")
    ax[0, 0].set_xlabel("energy (keV)")
    ax[0, 0].set_ylabel("counts / live s / keV")
    ax[0, 0].legend()
    ax[1, 0].set_title("⁶⁰Co region")
    ax[1, 0].set_xlabel("energy (keV)")
    ax[1, 0].set_ylabel("counts / live s / keV")
    ax[1, 0].legend()
    P.save(fig, out / "benchmark_spectra.png")
