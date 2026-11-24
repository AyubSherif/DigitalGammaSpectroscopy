"""Baseline restoration: pole-zero (P/Z) correction of the AC-coupling undershoot, and a
slow baseline tracker.

Pole-zero correction
--------------------
The sampled pulse is (to a very good approximation) a sum of exponentials, so its
z-transform is

    S(z) = N(z) / prod_i (1 - a_i z^-1),    a_i = exp(-1/tau_i)

The AC coupling appears in S twice: as slow *poles* a_i (the undershoot exponentials)
and as *zeros* of N(z) close to z = 1 (they force the zero net area). The exact
correction swaps them:

    G(z) = prod_slow (1 - a_i z^-1) / prod_near1 (1 - r_k z^-1)

Zeros of G cancel the undershoot poles; poles of G sit on N's near-1 zeros. G·S is
then the pulse without coupling: same main lobe, no undershoot. A near-1 root at or
above 1 (net area very slightly negative) would be an unstable pole, so it is clamped
to 1 - 1/leak_samples, which keeps the filter stable and any DC drift bounded.

The filter is a short IIR (2 zeros, 2 poles here), runs in O(1) per sample and keeps
its state between chunks, so it works on streams as well as files.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.signal import lfilter, lfilter_zi

from .template import PulseTemplate


# ----------------------------------------------------------------------------------------------
# Pole-zero correction
# ----------------------------------------------------------------------------------------------

@dataclass
class PoleZeroDesign:
    b: np.ndarray                 # numerator coefficients (in z^-1)
    a: np.ndarray                 # denominator coefficients (in z^-1)
    zero_taus: np.ndarray         # samples: undershoot time constants cancelled
    pole_taus: np.ndarray         # samples: replacement poles (after clamping)

    def describe(self) -> str:
        z = ", ".join(f"{t:.1f}" for t in self.zero_taus)
        p = ", ".join(f"{t:.1f}" for t in self.pole_taus)
        return f"P/Z: cancels tau = [{z}] samples, replacement poles tau = [{p}] samples"


def design_pole_zero(taus, amps, slow_threshold: float = 30.0, leak_samples: float = 20_000.0,
                     near_one: float = 0.05) -> PoleZeroDesign:
    """Design the P/Z filter from a multi-exponential pulse model (taus in samples).

    ``slow_threshold``: exponentials slower than this many samples are treated as
    coupling (undershoot) and cancelled. ``near_one``: N(z) roots within this distance
    of 1 are the coupling zeros that become the filter's poles.
    """
    taus = np.asarray(taus, dtype=float)
    amps = np.asarray(amps, dtype=float)
    a = np.exp(-1.0 / taus)
    N = np.zeros(len(a))
    for i in range(len(a)):
        p = np.array([1.0])
        for j in range(len(a)):
            if j != i:
                p = np.convolve(p, [1.0, -a[j]])
        N = N + amps[i] * p
    roots = np.roots(N)
    near = roots[np.abs(roots - 1.0) < near_one]
    near = np.real(near[np.abs(np.imag(near)) < 1e-9])
    slow = a[taus > slow_threshold]
    if near.size != slow.size:
        raise ValueError(f"found {slow.size} slow exponentials but {near.size} near-DC zeros; "
                         "cannot pair them for pole-zero correction")
    lam = 1.0 - 1.0 / leak_samples
    poles = np.minimum(near, lam)
    b = np.poly(slow)
    den = np.poly(poles)
    return PoleZeroDesign(b, den, taus[taus > slow_threshold], -1.0 / np.log(poles))


def design_from_template(template: PulseTemplate, **kw) -> PoleZeroDesign:
    """Design from the analytic model stored with a template (``characterize`` output)."""
    m = template.meta
    if "model_taus" not in m:
        raise ValueError("template has no fitted model (model_taus/model_amps)")
    return design_pole_zero(m["model_taus"], m["model_amps"], **kw)


class PoleZero:
    """Optional pipeline stage: pole-zero correction. Keeps filter state between calls."""

    name = "pole_zero"

    def __init__(self, design: PoleZeroDesign, initial_level: float = 0.0):
        self.design = design
        self.initial_level = initial_level
        self._zi = None

    def process(self, x: np.ndarray) -> np.ndarray:
        # Start at rest on the input's mean level (0 after offset removal). Starting on the
        # first sample instead would multiply that single noisy value by the filter's large
        # DC gain and leave a slowly decaying offset.
        if self._zi is None:
            self._zi = lfilter_zi(self.design.b, self.design.a) * self.initial_level
        y, self._zi = lfilter(self.design.b, self.design.a, x, zi=self._zi)
        return y

    def reset(self) -> None:
        self._zi = None


def filter_template(template: PulseTemplate, b, a, pre_samples: int = 200) -> PulseTemplate:
    """Apply a sample-rate IIR filter to a fine-grid template.

    The template is evaluated on each of the ``resolution`` sub-sample phases, each
    phase is filtered at the sample rate, and the results are interleaved again.
    The result is renormalised to peak 1 (the peak position stays at t = 0).
    """
    R = template.resolution
    k0 = int(np.floor(template.t_start)) - pre_samples
    k1 = int(np.ceil(template.t_end))
    k = np.arange(k0, k1 + 1, dtype=float)
    out_t, out_s = [], []
    for j in range(R):
        ph = j / R
        y = lfilter(b, a, template(k + ph))
        out_t.append(k + ph)
        out_s.append(y)
    t = np.stack(out_t, 1).ravel()
    s = np.stack(out_s, 1).ravel()
    keep = (t >= template.t_start) & (t <= template.t_end)
    t, s = t[keep], s[keep]
    i = int(np.argmax(s))
    meta = dict(template.meta)
    meta["filtered"] = True
    return PulseTemplate(R, float(t[0]), s / s[i], meta)


# ----------------------------------------------------------------------------------------------
# Slow baseline tracker
# ----------------------------------------------------------------------------------------------

@dataclass
class TrackerConfig:
    window: int = 1001            # samples in the moving average of pulse-free samples
    min_free: int = 50            # need at least this many free samples in a window
    mask_pre: int = 10            # samples masked before each pulse peak
    mask_post: int = 50           # ... and after it (main lobe after P/Z)


def track_baseline(x: np.ndarray, free: np.ndarray, cfg: TrackerConfig) -> np.ndarray:
    """Centred moving average of the pulse-free samples; gaps are bridged by interpolation."""
    w = cfg.window
    xs = np.concatenate([[0.0], np.cumsum(np.where(free, x, 0.0))])
    ns = np.concatenate([[0], np.cumsum(free.astype(np.int64))])
    i = np.arange(x.size)
    lo = np.clip(i - w // 2, 0, x.size)
    hi = np.clip(i + w // 2 + 1, 0, x.size)
    cnt = ns[hi] - ns[lo]
    est = np.where(cnt >= cfg.min_free, (xs[hi] - xs[lo]) / np.maximum(cnt, 1), np.nan)
    good = np.isfinite(est)
    if not good.any():
        return np.zeros_like(x)
    return np.interp(i, i[good], est[good])


class BaselineTracker:
    """Optional pipeline stage: subtract a slowly varying baseline estimated from pulse-free samples.

    Offline (whole capture) it is centred. A streaming version will use the same
    estimate with a fixed look-ahead delay.
    """

    name = "baseline_tracker"

    def __init__(self, config: TrackerConfig | None = None):
        self.cfg = config or TrackerConfig()
        self.last_baseline: np.ndarray | None = None

    def process(self, x: np.ndarray, free: np.ndarray) -> np.ndarray:
        self.last_baseline = track_baseline(x, free, self.cfg)
        return x - self.last_baseline

    def reset(self) -> None:
        self.last_baseline = None
