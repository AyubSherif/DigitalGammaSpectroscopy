"""Signal simulator: a synthetic NaI(Tl) PMT waveform, continuous and chunk by chunk (step 10).

The simulator is a test source that is not limited to the recorded data. Every pulse it
makes is known (the *truth*), so the pipeline can be scored against it
(``compare_to_truth``).

Model, per decay or background event:

1. **Emission.** A nuclide emits its gamma lines: ⁶⁰Co emits 1173 and 1332 keV *in
   coincidence* (same decay), ¹³⁷Cs emits 662 keV. The room background is ⁴⁰K and ²⁰⁸Tl lines
   plus a falling continuum.
2. **Interaction.** Each gamma deposits energy with probability ``efficiency``. Of the
   interacting gammas, a fraction P/T(E) (peak-to-total, falling with energy) deposit
   everything (photopeak), ``backscatter_fraction`` deposit the backscatter energy
   E / (1 + 2E/mc²), and the rest Compton-scatter once: the electron energy is drawn from the
   Klein–Nishina distribution up to the Compton edge.
3. **Summing.** Coincident deposits of one decay add up to one pulse (true-coincidence
   summing, e.g. the 2506 keV ⁶⁰Co sum peak). Chance coincidences (pile-up) happen by
   themselves, from the Poisson arrival times.
4. **Detector resolution.** The summed energy is smeared by a Gaussian with
   FWHM(E)/E = fwhm_pct_662 · √(662/E).
5. **Waveform.** Each pulse is amplitude · pulse_shape(t − t_arrival), with amplitude = E / gain
   and arrival times at sub-sample precision. The shape is the measured pulse (default: its
   analytic model, or a template file), including the AC-coupling undershoot. Then come
   Gaussian noise (white, or shaped like a measured noise spectrum), interference
   sinusoids, a DC offset, the polarity of the PMT signal (negative) and ADC quantisation.

``Simulator.read(n)`` returns the next n samples. The waveform is continuous across calls
whatever their length: pulses straddling a chunk boundary are completed in the next chunk,
sinusoids keep their phase, and shaped noise carries its filter tail.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import fftconvolve

from .template import PulseTemplate, multiexp_model

ME_KEV = 510.999

# The analytic pulse model fitted to the 2021 template (characterize; ADR-0004):
# s(t) = sum a_k exp(-(t-t0)/tau_k) * Gaussian(sigma), t in samples (1.8 MS/s).
MODEL_2021 = dict(taus=[0.44158403, 7.51267126, 75.93783709, 317.66654735],
                  amps=[-1.27583098, 1.34773809, -0.05439575, -0.01751136],
                  sigma=0.0643781094904544, t0=-1.1869147489802732)

NUCLIDES = {
    # name: list of gamma cascades; each cascade is [(energy keV, emission probability), ...]
    # whose gammas come from the same decay (coincident)
    "co60": [[(1173.228, 0.9985), (1332.492, 0.999826)]],
    "cs137": [[(661.657, 0.851)]],
    "k40": [[(1460.82, 1.0)]],
    "none": [],
}


@dataclass
class SimConfig:
    nuclide: str = "co60"                 # source: "co60" | "cs137" | "k40" | "none" (background only)
    rate_cps: float = 3900.0              # source pulses per second (any energy deposit; 2021 run: ~3900)
    background_cps: float = 600.0         # room background pulses per second (2021 background run: ~600)
    background_k40_fraction: float = 0.07  # of the background: 40K 1461 keV gammas
    background_tl208_fraction: float = 0.02  # of the background: 208Tl 2615 keV gammas
    background_continuum_kev: float = 250.0  # the rest: deposits falling as exp(-E / this)
    efficiency: float = 0.1               # probability that an emitted gamma deposits energy (sets summing)
    peak_to_total_662: float = 0.45       # photopeak fraction of interactions at 662 keV ...
    peak_to_total_slope: float = 0.6      # ... scaling as (E / 662)^-slope
    backscatter_fraction: float = 0.04    # interactions depositing the backscatter energy
    fwhm_pct_662: float = 8.0             # intrinsic resolution at 662 keV (scales as 1/sqrt(E))
    gain_kev_per_v: float = 3760.0        # pulse amplitude (V) = energy / gain
    fs: float = 1.8e6                     # sample rate (Hz)
    polarity: int = -1                    # -1: negative-going PMT pulses (as recorded)
    offset_v: float = 0.004               # DC offset of the signal
    noise_v: float = 0.0006               # rms of the Gaussian noise
    noise_psd: str = ""                   # noise.h5 from characterize: shape the noise like it ("" = white)
    lines_hz: tuple = (9000.0, 18000.0, 27000.0, 257000.0)   # interference sinusoids: frequencies ...
    lines_v: tuple = (0.0004, 0.0003, 0.0002, 0.0005)        # ... and amplitudes (V)
    lsb_v: float = 8.4771697e-5           # ADC step (0 = no quantisation)
    pulse_template: str = ""              # pulse shape: template.h5 ("" = the 2021 analytic model)
    seed: int = 1                         # random seed: the same seed gives the same signal


def analytic_template(model: dict = MODEL_2021, resolution: int = 16, t_min: float = -40.0,
                      t_max: float = 1500.0) -> PulseTemplate:
    """A PulseTemplate from the analytic model, peak = 1 at t = 0, with the model stored as in
    ``characterize`` output (so the pipeline can design its pole-zero filter from it)."""
    taus = np.asarray(model["taus"], float)
    amps = np.asarray(model["amps"], float)
    t = t_min + np.arange(int((t_max - t_min) * resolution) + 1) / resolution
    y = multiexp_model(t, model["t0"], model["sigma"], taus, amps[1:])
    # find the peak on a fine grid and put it at t = 0 with value 1
    tf = np.linspace(t[np.argmax(y)] - 1, t[np.argmax(y)] + 1, 2001)
    yf = multiexp_model(tf, model["t0"], model["sigma"], taus, amps[1:])
    tp, peak = tf[np.argmax(yf)], yf.max()
    y = multiexp_model(t + tp, model["t0"], model["sigma"], taus, amps[1:]) / peak
    meta = dict(model_taus=taus, model_amps=amps / peak, model_sigma=float(model["sigma"]),
                model_t0=float(model["t0"] - tp), polarity=-1, source_store="analytic model (simulate.MODEL_2021)")
    return PulseTemplate(resolution, t_min, y, meta)


def load_pulse(cfg: SimConfig) -> PulseTemplate:
    return PulseTemplate.load(cfg.pulse_template) if cfg.pulse_template else analytic_template()


# ----------------------------------------------------------------------------------------------
# Energy deposits
# ----------------------------------------------------------------------------------------------

def klein_nishina_electron_cdf(e_kev: float, n: int = 512) -> tuple[np.ndarray, np.ndarray]:
    """Grid of Compton electron energies T in [0, T_max] and the CDF of the Klein-Nishina
    distribution dσ/dT ∝ 2 + s²/(ε²(1-s)²) + s/(1-s)·(s - 2/ε), s = T/E, ε = E/mc²."""
    eps = e_kev / ME_KEV
    t_max = e_kev * 2 * eps / (1 + 2 * eps)
    T = np.linspace(0, t_max, n)
    s = np.minimum(T / e_kev, 1 - 1e-12)
    pdf = 2 + s ** 2 / (eps ** 2 * (1 - s) ** 2) + s / (1 - s) * (s - 2 / eps)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(T))])
    return T, cdf / cdf[-1]


class DepositModel:
    """Energy deposited by one gamma (0 if it does not interact)."""

    def __init__(self, cfg: SimConfig):
        self.cfg = cfg
        self._kn: dict[float, tuple] = {}

    def deposit(self, e_kev: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        c = self.cfg
        e = np.asarray(e_kev, float)
        out = np.zeros(e.size)
        hit = rng.random(e.size) < c.efficiency
        u = rng.random(e.size)
        pt = np.clip(c.peak_to_total_662 * (e / 661.657) ** -c.peak_to_total_slope, 0, 1)
        peak = hit & (u < pt)
        back = hit & (u >= pt) & (u < pt + c.backscatter_fraction)
        compton = hit & ~peak & ~back
        out[peak] = e[peak]
        out[back] = e[back] / (1 + 2 * e[back] / ME_KEV)
        for en in np.unique(e[compton]):
            sel = compton & (e == en)
            T, cdf = self._kn.setdefault(float(en), klein_nishina_electron_cdf(float(en)))
            out[sel] = np.interp(rng.random(int(sel.sum())), cdf, T)
        return out


# ----------------------------------------------------------------------------------------------
# The simulator
# ----------------------------------------------------------------------------------------------

class Simulator:
    """Continuous synthetic waveform. ``read(n)`` -> the next n samples (volts, as recorded)."""

    def __init__(self, cfg: SimConfig | None = None, t0_s: float = 0.0):
        self.cfg = c = cfg or SimConfig()
        if c.nuclide not in NUCLIDES:
            raise ValueError(f"unknown nuclide {c.nuclide!r}; choose from {sorted(NUCLIDES)}")
        self.fs, self.dt, self.t0_s = c.fs, 1.0 / c.fs, t0_s
        self.pulse = load_pulse(c)
        self._lo = int(np.floor(self.pulse.t_start))            # pulse support, samples rel. to arrival
        self._hi = int(np.ceil(self.pulse.t_end))
        # the pulse at integer sample offsets lo..hi, for _PHASES sub-sample arrival phases:
        # rendering is then a table lookup (timing quantised to 1/_PHASES sample)
        j = np.arange(self._lo, self._hi + 1, dtype=np.float64)
        self._table = np.stack([self.pulse(j - ph / _PHASES) for ph in range(_PHASES)])
        ss = np.random.SeedSequence(c.seed)
        self._rng_ev, self._rng_noise = (np.random.default_rng(s) for s in ss.spawn(2))
        self._deposits = DepositModel(c)
        # source decays per second giving the requested pulse rate
        self._cascades = NUCLIDES[c.nuclide]
        per_decay = self._pulses_per_decay() if self._cascades else 0.0
        self.decay_rate = c.rate_cps / per_decay if per_decay > 0 else 0.0
        self.pos = 0                                            # next sample to produce
        self._pending = np.zeros((0, 3))                       # (time [samples], amplitude, energy)
        self._t_gen = 0.0                                       # events generated up to this time
        self.truth: list[np.ndarray] = []                       # rows: time_samples, amplitude_v, energy_kev, kind
        self._noise_fir = self._noise_filter()
        self._noise_tail = np.zeros(0)

    # -- events --------------------------------------------------------------------------------
    def _pulses_per_decay(self, n: int = 20000) -> float:
        rng = np.random.default_rng(12345)
        e = self._decay_energies(n, rng)
        return float(np.mean(e > 0))

    def _decay_energies(self, n: int, rng) -> np.ndarray:
        """Deposited energy (keV, before resolution) of n source decays (cascade summing)."""
        total = np.zeros(n)
        k = rng.integers(len(self._cascades), size=n)
        for ci, cascade in enumerate(self._cascades):
            sel = k == ci
            for e_kev, prob in cascade:
                emitted = sel & (rng.random(n) < prob)
                total[emitted] += self._deposits.deposit(np.full(int(emitted.sum()), e_kev), rng)
        return total

    def _background_energies(self, n: int, rng) -> np.ndarray:
        c = self.cfg
        u = rng.random(n)
        out = rng.exponential(c.background_continuum_kev, n)
        for frac_lo, frac_hi, e_kev in ((0.0, c.background_k40_fraction, 1460.82),
                                        (c.background_k40_fraction,
                                         c.background_k40_fraction + c.background_tl208_fraction, 2614.511)):
            sel = (u >= frac_lo) & (u < frac_hi)
            d = np.zeros(int(sel.sum()))
            while True:                                         # these are pulses: redraw non-interacting
                miss = d <= 0
                if not miss.any():
                    break
                d[miss] = self._deposits.deposit(np.full(int(miss.sum()), e_kev), rng)
            out[sel] = d
        return out

    def _smear(self, e: np.ndarray, rng) -> np.ndarray:
        c = self.cfg
        fwhm = c.fwhm_pct_662 / 100 * np.sqrt(661.657 * np.maximum(e, 1e-9))
        return np.maximum(e + rng.standard_normal(e.size) * fwhm / 2.3548, 0.0)

    def _generate_until(self, t_end: float) -> None:
        """Create the events with arrival time < t_end (samples)."""
        c, rng = self.cfg, self._rng_ev
        if t_end <= self._t_gen:
            return
        # in fixed blocks, so the random sequence (and the waveform) does not depend on how
        # the caller splits the reads
        t_end = np.ceil(t_end / _GEN_BLOCK) * _GEN_BLOCK
        rows = []
        while self._t_gen < t_end:
            rows += self._generate_block(self._t_gen, self._t_gen + _GEN_BLOCK)
            self._t_gen += _GEN_BLOCK
        if rows:
            new = np.concatenate(rows)
            new = new[np.argsort(new[:, 0])]
            self.truth.append(new)
            self._pending = np.concatenate([self._pending, new[:, :3]])

    def _generate_block(self, t_a: float, t_b: float) -> list:
        c, rng = self.cfg, self._rng_ev
        span = t_b - t_a
        rows = []
        for rate, kind in ((self.decay_rate, 0), (c.background_cps, 1)):
            lam = rate * span * self.dt
            if lam <= 0:
                continue
            n = rng.poisson(lam)
            t = t_a + rng.random(n) * span
            e = self._decay_energies(n, rng) if kind == 0 else self._background_energies(n, rng)
            keep = e > 0
            t, e = t[keep], self._smear(e[keep], rng)
            rows.append(np.column_stack([t, e / c.gain_kev_per_v, e, np.full(t.size, kind)]))
        return rows

    # -- waveform ------------------------------------------------------------------------------
    def _noise_filter(self) -> np.ndarray | None:
        if not self.cfg.noise_psd:
            return None
        with h5py.File(self.cfg.noise_psd, "r") as f:
            fr, psd = f["freq_hz"][:], f["psd"][:]
        nfft = 2 * (fr.size - 1)
        amp = np.sqrt(np.interp(np.fft.rfftfreq(nfft, self.dt), fr, psd))
        h = np.fft.fftshift(np.fft.irfft(amp, nfft)) * np.hanning(nfft)
        return h / np.sqrt(np.sum(h ** 2))                     # unit gain for white noise

    def _noise(self, n: int) -> np.ndarray:
        c = self.cfg
        if self._noise_fir is None:
            return c.noise_v * self._rng_noise.standard_normal(n)
        y = fftconvolve(self._rng_noise.standard_normal(n), self._noise_fir)   # length n + m - 1
        y[: self._noise_tail.size] += self._noise_tail
        self._noise_tail = y[n:].copy()
        return c.noise_v * y[:n]

    def read(self, n: int) -> np.ndarray:
        """The next n samples, in volts as the ADC records them."""
        c = self.cfg
        s0, s1 = self.pos, self.pos + n
        self._generate_until(s1 - self._lo + 1)                # pulses starting before the chunk end
        x = np.zeros(n)
        P = self._pending
        if P.size:
            act = (P[:, 0] + self._hi >= s0) & (P[:, 0] + self._lo < s1)
            A = P[act]
            if A.size:
                base = np.floor(A[:, 0])
                ph = np.rint((A[:, 0] - base) * _PHASES).astype(np.int64)
                base = base.astype(np.int64) + (ph == _PHASES)
                ph[ph == _PHASES] = 0
                k = (base + self._lo)[:, None] + np.arange(self._table.shape[1])[None, :]
                val = self._table[ph] * A[:, 1:2]
                ok = (k >= s0) & (k < s1)
                x += np.bincount(k[ok] - s0, weights=val[ok], minlength=n)
            self._pending = P[P[:, 0] + self._hi >= s1]        # still contributing after this chunk
        t = (s0 + np.arange(n)) * self.dt
        for f_hz, a_v in zip(c.lines_hz, c.lines_v):
            x += a_v * np.sin(2 * np.pi * f_hz * t + f_hz * 1e-3)
        x += self._noise(n)
        v = c.polarity * x + c.offset_v
        if c.lsb_v > 0:
            v = np.round(v / c.lsb_v) * c.lsb_v
        self.pos = s1
        return v

    def truth_table(self) -> np.ndarray:
        """Events produced so far: time_s, sample (float), amplitude_v, energy_kev, kind (0 source, 1 bkg)."""
        T = np.concatenate(self.truth) if self.truth else np.zeros((0, 4))
        T = T[T[:, 0] < self.pos]
        out = np.zeros(T.shape[0], TRUTH_DTYPE)
        out["sample"], out["amplitude_v"], out["energy_kev"], out["kind"] = T[:, 0], T[:, 1], T[:, 2], T[:, 3]
        out["time_s"] = self.t0_s + T[:, 0] * self.dt
        return out


_PHASES = 128                                              # sub-sample phases of the pulse table
_GEN_BLOCK = 65536.0                                       # samples per event-generation block

TRUTH_DTYPE = np.dtype([("time_s", "f8"), ("sample", "f8"), ("amplitude_v", "f4"), ("energy_kev", "f4"),
                        ("kind", "u1")])


def save_truth(path, truth: np.ndarray, cfg: SimConfig, extra: dict | None = None) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for name in TRUTH_DTYPE.names:
            f.create_dataset(name, data=truth[name], compression="gzip", compression_opts=4, shuffle=True)
        f.attrs["sim_config"] = json.dumps({k: (list(v) if isinstance(v, tuple) else v)
                                            for k, v in cfg.__dict__.items()})
        for k, v in (extra or {}).items():
            f.attrs[k] = v
    return path


def load_truth(path) -> np.ndarray:
    with h5py.File(path, "r") as f:
        out = np.zeros(f["time_s"].shape[0], TRUTH_DTYPE)
        for name in TRUTH_DTYPE.names:
            out[name] = f[name][:]
    return out


# ----------------------------------------------------------------------------------------------
# Simulated capture store (like the recorded data: 0.29 s captures)
# ----------------------------------------------------------------------------------------------

def simulate_store(out_path, cfg: SimConfig, n_captures: int, n_samples: int = 524288,
                   gap_s: float = 0.0, overwrite: bool = False) -> tuple[Path, Path]:
    """Write simulated captures to an HDF5 capture store (+ truth file next to it).

    With gap_s = 0 the captures are consecutive pieces of one continuous waveform; with a gap,
    the simulator runs on through the gap (nothing is recorded there), like the digitiser.
    """
    from .store import write_store

    if cfg.lsb_v <= 0:
        raise ValueError("a capture store holds integer ADC codes: set lsb_v > 0")
    sim = Simulator(cfg)
    gap = int(round(gap_s * cfg.fs))
    times, truths = [], []

    def captures():
        for i in range(n_captures):
            times.append(sim.pos * sim.dt)
            v = sim.read(n_samples)
            yield np.rint(v / cfg.lsb_v).astype(np.int32)
            if gap:
                sim._generate_until(sim.pos + gap)              # events during the gap are not recorded
                sim._pending = sim._pending[sim._pending[:, 0] + sim._hi >= sim.pos + gap]
                sim.pos += gap

    out_path = Path(out_path)
    write_store(out_path, captures(), n_captures, n_samples, 1.0 / cfg.fs, cfg.lsb_v,
                times_s=times, source=f"simulation ({cfg.nuclide}, seed {cfg.seed})", overwrite=overwrite)
    # truth inside the recorded windows, with the capture and the sample within it
    T = sim.truth_table()
    cap_start = np.array(times) / sim.dt
    ci = np.searchsorted(cap_start, T["sample"], side="right") - 1
    inside = (ci >= 0) & (T["sample"] - cap_start[np.maximum(ci, 0)] < n_samples)
    T = T[inside]
    truth_path = save_truth(out_path.with_name(out_path.stem + ".truth.h5"), T, cfg,
                            dict(capture_start_samples=cap_start, n_samples=n_samples))
    return out_path, truth_path


# ----------------------------------------------------------------------------------------------
# Scoring the pipeline against the truth
# ----------------------------------------------------------------------------------------------

def compare_to_truth(event_path, truth_path, tolerance_samples: float = 2.0,
                     energy_bins=(30, 100, 300, 600, 1000, 1500, 3000)) -> dict:
    """Match detected events to true pulses by time (±tolerance) and score the pipeline.

    Uses event ``time_s`` and truth ``time_s`` on the same clock (both from the simulator's
    t0). Returns detection efficiency by true energy, the false-detection rate, the amplitude
    error of matched good events, and the fraction of close true pairs that are resolved.
    """
    from .events import EventFile

    ef = EventFile(event_path)
    dt = float(ef.attrs["dt_s"])
    ev = ef.events
    T = load_truth(truth_path)
    good_mask = _good_mask(ev)
    tt = T["time_s"]
    order = np.argsort(tt)
    T, tt = T[order], tt[order]
    # only truth inside processed (live) regions: event file captures cover them
    caps = ef.captures
    starts = np.where(np.isfinite(caps["live_start_s"]), caps["live_start_s"], caps["time_s"])
    covered = np.zeros(tt.size, bool)
    for t0, live in zip(starts, caps["live_time_s"]):
        covered |= (tt >= t0) & (tt < t0 + live)
    te = ev["time_s"]
    j = np.searchsorted(tt, te)
    jl, jr = np.clip(j - 1, 0, max(tt.size - 1, 0)), np.clip(j, 0, max(tt.size - 1, 0))
    if tt.size:
        dl, dr = np.abs(te - tt[jl]), np.abs(te - tt[jr])
        nearest = np.where(dl <= dr, jl, jr)
        dmin = np.minimum(dl, dr) / dt
    else:
        nearest, dmin = np.zeros(te.size, int), np.full(te.size, np.inf)
    matched_ev = dmin <= tolerance_samples
    found = np.zeros(tt.size, bool)
    found[nearest[matched_ev]] = True
    ev_in = np.zeros(te.size, bool)
    for t0, live in zip(starts, caps["live_time_s"]):
        ev_in |= (te >= t0) & (te < t0 + live)
    live = float(caps["live_time_s"].sum())
    E = T["energy_kev"]
    eff = []
    for lo, hi in zip(energy_bins[:-1], energy_bins[1:]):
        sel = covered & (E >= lo) & (E < hi)
        eff.append(dict(energy_kev=[lo, hi], true=int(sel.sum()),
                        efficiency=float(found[sel].mean()) if sel.any() else None))
    false = ev_in & ~matched_ev
    m = matched_ev & good_mask & ev_in
    rel = (ev["amplitude"][m] - T["amplitude_v"][nearest[m]]) / T["amplitude_v"][nearest[m]]
    big = T["amplitude_v"][nearest[m]] > 0.05
    # close true pairs (both inside, separation 2-60 samples): both found as distinct events?
    gaps = np.diff(T["sample"]) if T.size > 1 else np.zeros(0)
    pair = covered[:-1] & covered[1:] & (gaps >= 2) & (gaps < 60) & (E[:-1] > 100) & (E[1:] > 100)
    both = found[:-1] & found[1:]
    return dict(
        live_time_s=live, true_pulses=int(covered.sum()), true_rate_cps=float(covered.sum() / live) if live else 0.0,
        detected=int(ev_in.sum()), detected_rate_cps=float(ev_in.sum() / live) if live else 0.0,
        efficiency_by_energy=eff,
        false_rate_cps=float(false.sum() / live) if live else 0.0,
        amplitude_bias_pct=float(100 * np.median(rel[big])) if big.any() else None,
        amplitude_spread_pct=float(100 * 1.4826 * np.median(np.abs(rel[big] - np.median(rel[big])))) if big.any() else None,
        close_pairs=int(pair.sum()), close_pairs_resolved=float(both[pair].mean()) if pair.any() else None,
        good_fraction_of_matched=float((matched_ev & good_mask & ev_in).sum() / max((matched_ev & ev_in).sum(), 1)),
    )


def _good_mask(ev: np.ndarray) -> np.ndarray:
    from .events import INFO_FLAGS

    return ((ev["flags"] & np.uint16(~np.uint16(INFO_FLAGS))) == 0) & np.isfinite(ev["amplitude"])
