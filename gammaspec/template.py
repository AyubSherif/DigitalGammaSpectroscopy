"""Pulse template: construction from data, per-pulse amplitude/timing fit, analytic model.

Conventions
-----------
* Time is measured in **samples**. A template is a function s(t) on a fine grid with
  ``resolution`` points per sample, normalised so that its continuous maximum is 1 at t = 0.
* A pulse window is cut around a coarse peak index p: window column k is sample p - pre + k,
  i.e. offset j = k - pre relative to p.
* The pulse model in a window is  w(j) = A * s(j - delta) + b,  where A is the true
  (phase-independent) peak amplitude, delta the sub-sample position of the continuous peak
  relative to sample p, and b a constant baseline.

Because pulses arrive at random sampling phases, stacking many pulses aligned by their
fitted delta samples the shape on a grid much finer than the ADC's ("equivalent-time
sampling").
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import least_squares
from scipy.special import erfc, erfcx


# ----------------------------------------------------------------------------------------------
# Template container
# ----------------------------------------------------------------------------------------------

@dataclass
class PulseTemplate:
    resolution: int          # grid points per sample
    t_start: float           # time (samples) of the first grid point
    shape: np.ndarray        # normalised shape, max 1 at t = 0
    meta: dict = field(default_factory=dict)

    @property
    def t(self) -> np.ndarray:
        return self.t_start + np.arange(self.shape.size) / self.resolution

    @property
    def t_end(self) -> float:
        return self.t_start + (self.shape.size - 1) / self.resolution

    def __call__(self, t: np.ndarray) -> np.ndarray:
        """Evaluate at arbitrary times (samples). Zero before the start and after the end."""
        return np.interp(t, self.t, self.shape, left=0.0, right=0.0)

    def save(self, path: str | Path) -> None:
        with h5py.File(path, "w") as f:
            f.create_dataset("shape", data=self.shape)
            f.attrs["resolution"] = self.resolution
            f.attrs["t_start"] = self.t_start
            for k, v in self.meta.items():
                f.attrs[k] = v

    @classmethod
    def load(cls, path: str | Path) -> "PulseTemplate":
        with h5py.File(path, "r") as f:
            meta = {k: v for k, v in f.attrs.items() if k not in ("resolution", "t_start")}
            return cls(int(f.attrs["resolution"]), float(f.attrs["t_start"]), f["shape"][:], meta)


# ----------------------------------------------------------------------------------------------
# Windows and per-pulse fit
# ----------------------------------------------------------------------------------------------

def extract_windows(x: np.ndarray, peaks: np.ndarray, pre: int, post: int) -> np.ndarray:
    """Rows x[p-pre : p+post+1] for each peak (peaks must leave room at both ends)."""
    idx = np.asarray(peaks)[:, None] + np.arange(-pre, post + 1)[None, :]
    return x[idx]


@dataclass
class PulseFit:
    amplitude: np.ndarray   # A, volts
    baseline: np.ndarray    # b, volts
    shift: np.ndarray       # delta, samples (continuous peak position relative to the window's peak sample)
    rms: np.ndarray         # rms residual of the fit, volts

    def subset(self, mask: np.ndarray) -> "PulseFit":
        return PulseFit(self.amplitude[mask], self.baseline[mask], self.shift[mask], self.rms[mask])


def _linear_fit(Y: np.ndarray, S: np.ndarray):
    """Row-wise least squares of Y ≈ A*S + b (S broadcastable to Y). Returns A, b, sse."""
    n = Y.shape[-1]
    ss = (S * S).sum(-1)
    s1 = S.sum(-1)
    sy = (S * Y).sum(-1)
    y1 = Y.sum(-1)
    yy = (Y * Y).sum(-1)
    det = ss * n - s1 * s1
    A = (n * sy - s1 * y1) / det
    b = (ss * y1 - s1 * sy) / det
    sse = yy - A * sy - b * y1
    return A, b, np.maximum(sse, 0.0)


def fit_pulses(windows: np.ndarray, template: PulseTemplate, pre: int,
               fit_range: tuple[int, int] = (-10, 60), max_shift: float = 1.0,
               steps_per_sample: int = 16, fit_baseline: bool = True) -> PulseFit:
    """Fit A, b and the sub-sample shift of every window against the template.

    Grid search over shifts in [-max_shift, max_shift], parabolic refinement of the
    best shift, then an exact linear solve at the refined shift. With
    ``fit_baseline=False`` the baseline is fixed at 0 (use after baseline restoration).
    """
    lo, hi = fit_range
    Y = windows[:, pre + lo: pre + hi + 1].astype(np.float64)
    j = np.arange(lo, hi + 1, dtype=np.float64)
    deltas = np.arange(-max_shift, max_shift + 1e-9, 1.0 / steps_per_sample)
    S = template(j[None, :] - deltas[:, None])                     # (G, J)
    n = j.size
    ss, s1 = (S * S).sum(1), S.sum(1)
    SY = Y @ S.T                                                   # (N, G)
    y1 = Y.sum(1)[:, None]
    yy = (Y * Y).sum(1)[:, None]
    if fit_baseline:
        det = ss * n - s1 * s1
        A = (n * SY - s1 * y1) / det
        b = (ss * y1 - s1 * SY) / det
        sse = yy - A * SY - b * y1                                 # (N, G)
    else:
        A = SY / ss
        b = np.zeros_like(A)
        sse = yy - A * SY

    g = np.argmin(sse, axis=1)
    gi = np.clip(g, 1, deltas.size - 2)
    rows = np.arange(Y.shape[0])
    c0, cm, cp = sse[rows, gi], sse[rows, gi - 1], sse[rows, gi + 1]
    denom = cm - 2 * c0 + cp
    off = np.where(denom > 0, 0.5 * (cm - cp) / np.where(denom > 0, denom, 1), 0.0)
    delta = deltas[gi] + np.clip(off, -1, 1) / steps_per_sample

    Si = template(j[None, :] - delta[:, None])                     # (N, J)
    if fit_baseline:
        A, b, sse = _linear_fit(Y, Si)
    else:
        sy, ss_i = (Si * Y).sum(-1), (Si * Si).sum(-1)
        A, b = sy / ss_i, np.zeros(Y.shape[0])
        sse = np.maximum((Y * Y).sum(-1) - A * sy, 0.0)
    return PulseFit(A, b, delta, np.sqrt(sse / n))


# ----------------------------------------------------------------------------------------------
# Template construction
# ----------------------------------------------------------------------------------------------

def _stack(windows, pre, fit: PulseFit, resolution, weights, smooth_from):
    """Average (w - b)/A on a fine grid of t = j - delta. Returns (t_start, shape)."""
    N, L = windows.shape
    j = np.arange(L) - pre
    tau = j[None, :] - fit.shift[:, None]
    t_start = -pre - 1.0
    idx = np.rint((tau - t_start) * resolution).astype(np.int64)
    nbins = int(np.rint((L + 1) * resolution)) + 1
    norm = (windows - fit.baseline[:, None]) / fit.amplitude[:, None]
    w = np.broadcast_to(weights[:, None], norm.shape)
    num = np.bincount(idx.ravel(), (norm * w).ravel(), nbins)
    den = np.bincount(idx.ravel(), w.ravel(), nbins)
    t = t_start + np.arange(nbins) / resolution
    good = den > 0
    shape = np.interp(t, t[good], num[good] / den[good])
    # Smoothing that never blurs the rising edge: almost none within +-3 samples of the
    # peak, ~1/8 sample elsewhere in the main lobe, ~1 sample in the slow tail.
    edge = np.abs(t) <= 3.0
    tail = t > smooth_from
    smoothed = _gauss_smooth(shape, resolution / 8.0)
    smoothed = np.where(tail, _gauss_smooth(shape, float(resolution)), smoothed)
    shape = np.where(edge, _gauss_smooth(shape, 0.5), smoothed)
    # The pre-pulse level is zero by definition. Without this the template and the
    # per-pulse baselines could trade a constant offset between iterations.
    pre_region = (t >= -pre + 1) & (t <= -4)
    if pre_region.any():
        shape = shape - shape[pre_region].mean()
    # drop the edges, where only some sampling phases contribute
    keep = (t >= -pre + 1) & (t <= L - pre - 2)
    return t[keep], shape[keep]


def _gauss_smooth(y: np.ndarray, sigma_bins: float) -> np.ndarray:
    if sigma_bins <= 0:
        return y
    h = int(np.ceil(4 * sigma_bins))
    k = np.exp(-0.5 * (np.arange(-h, h + 1) / sigma_bins) ** 2)
    k /= k.sum()
    return np.convolve(np.pad(y, h, mode="edge"), k, mode="valid")


def _recentre(t, shape, resolution):
    """Shift so the continuous maximum sits at t=0 and has value 1; resample onto a grid with a point at 0."""
    i = int(np.argmax(shape))
    if 0 < i < shape.size - 1:
        ym, y0, yp = shape[i - 1], shape[i], shape[i + 1]
        den = ym - 2 * y0 + yp
        off = 0.5 * (ym - yp) / den if den < 0 else 0.0
        peak_t = t[i] + off / resolution
        peak_v = y0 - 0.25 * (ym - yp) * off
    else:
        peak_t, peak_v = t[i], shape[i]
    k0 = np.ceil((t[0] - peak_t) * resolution)
    k1 = np.floor((t[-1] - peak_t) * resolution)
    new_t = np.arange(k0, k1 + 1) / resolution
    new_shape = np.interp(new_t + peak_t, t, shape) / peak_v
    return float(new_t[0]), new_shape, float(peak_t)


def equalize_phase(shift: np.ndarray) -> np.ndarray:
    """Make the fractional sampling phases of the shifts uniformly distributed.

    Pulses arrive at random times, so their true sampling phases are uniform on [0, 1).
    Stacking pulses by fitted shifts has an ambiguity: a template warped within each
    sample period, together with correspondingly warped shifts, fits the data almost as
    well. Mapping the fitted phases through their empirical CDF onto a uniform
    distribution removes the ambiguity ("code density" / histogram equalisation).
    """
    phi = np.mod(shift, 1.0)
    order = np.argsort(phi, kind="stable")
    target = np.empty_like(phi)
    target[order] = (np.arange(phi.size) + 0.5) / phi.size
    return shift + (target - phi)


def build_template(windows: np.ndarray, pre: int, resolution: int = 16, n_iter: int = 8,
                   fit_range: tuple[int, int] = (-10, 60), smooth_from: float = 40.0,
                   outlier_rms: float = 3.0, tol: float = 5e-4, equalize: bool = True,
                   verbose: bool = False):
    """Iteratively build a fine-grid template from isolated pulse windows.

    1. start from peak-sample alignment and peak-sample amplitudes;
    2. stack normalised pulses on the fine grid, recentre, renormalise;
    3. refit every pulse (A, b, delta) against the new template; drop poor fits;
       equalise the sampling phases (see ``equalize_phase``);
    4. repeat until the template stops changing (rms change < tol, fraction of peak).

    Pulses are weighted by A² (inverse variance of the normalised samples).
    Returns (template, fit, used_mask).
    """
    windows = np.asarray(windows, dtype=np.float64)
    N = windows.shape[0]
    b0 = windows[:, : max(1, pre - 5)].mean(axis=1)
    fit = PulseFit(windows[:, pre] - b0, b0, np.zeros(N), np.zeros(N))
    used = fit.amplitude > 0
    prev = None
    history = []
    for it in range(n_iter):
        sub = fit.subset(used)
        t, shape = _stack(windows[used], pre, sub, resolution, sub.amplitude ** 2, smooth_from)
        t_start, shape, peak_t = _recentre(t, shape, resolution)
        tmpl = PulseTemplate(resolution, t_start, shape)
        fit = fit_pulses(windows, tmpl, pre, fit_range, steps_per_sample=resolution)
        med_rms = np.median(fit.rms)
        used = (fit.amplitude > 0) & (fit.rms < outlier_rms * med_rms)
        if equalize and used.sum() >= 100:
            fit.shift[used] = equalize_phase(fit.shift[used])
        change = np.inf if prev is None else float(np.sqrt(np.mean((tmpl(prev.t) - prev.shape) ** 2)))
        history.append(change)
        if verbose:
            print(f"  iter {it}: used {used.sum()}/{N}, median rms {med_rms:.3g} V, max change {change:.2e}")
        if change < tol:
            break
        prev = tmpl
    tmpl.meta.update(n_pulses=int(used.sum()), resolution=resolution, iterations=it + 1,
                     final_change=history[-1])
    return tmpl, fit, used


# ----------------------------------------------------------------------------------------------
# Analytic model: sum of exponentials convolved with a Gaussian (finite front-end bandwidth)
# ----------------------------------------------------------------------------------------------

def _exp_gauss(x: np.ndarray, tau: float, sigma: float) -> np.ndarray:
    """exp(-x/tau)·H(x) convolved with a unit-area Gaussian of width sigma (stable form)."""
    z = (sigma / tau - x / sigma) / np.sqrt(2.0)
    out = np.empty_like(x, dtype=np.float64)
    pos = z >= 0
    out[pos] = 0.5 * np.exp(-0.5 * (x[pos] / sigma) ** 2) * erfcx(z[pos])
    neg = ~pos
    out[neg] = 0.5 * np.exp(0.5 * (sigma / tau) ** 2 - x[neg] / tau) * erfc(z[neg])
    return out


def multiexp_model(t: np.ndarray, t0: float, sigma: float, taus, amps_rest) -> np.ndarray:
    """s(t) = sum_k a_k exp(-(t-t0)/tau_k) for t > t0, smoothed by a Gaussian of width sigma.

    The first amplitude is fixed by sum(a_k) = 0, so the pulse starts from zero (the PMT
    current is integrated, so the voltage starts continuously from the baseline).
    """
    amps = np.concatenate([[-np.sum(amps_rest)], amps_rest])
    x = np.asarray(t, dtype=np.float64) - t0
    return sum(a * _exp_gauss(x, tau, sigma) for a, tau in zip(amps, taus))


@dataclass
class ModelFit:
    t0: float
    sigma: float
    taus: np.ndarray        # samples
    amps: np.ndarray        # all amplitudes (first one implied)
    rms: float              # rms residual vs the template (fraction of peak)
    max_abs: float          # max |residual|

    def __call__(self, t):
        return multiexp_model(t, self.t0, self.sigma, self.taus, self.amps[1:])


def fit_multiexp(template: PulseTemplate, taus0=(0.4, 6.5, 460.0), amps0=(1.2, -0.04),
                 t_min: float | None = None) -> ModelFit:
    """Fit the analytic model to a template. Taus and amplitudes are given in samples
    and peak units; len(amps0) = len(taus0) - 1."""
    t = template.t
    y = template.shape
    if t_min is not None:
        sel = t >= t_min
        t, y = t[sel], y[sel]
    nt = len(taus0)
    p0 = np.concatenate([[-1.5, 0.3], np.log(taus0), amps0])

    def resid(p):
        return multiexp_model(t, p[0], p[1], np.exp(p[2:2 + nt]), p[2 + nt:]) - y

    lo = np.concatenate([[-10, 0.01], np.full(nt, np.log(0.05)), np.full(nt - 1, -50)])
    hi = np.concatenate([[2, 5], np.full(nt, np.log(1e5)), np.full(nt - 1, 50)])
    r = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac", max_nfev=5000)
    p = r.x
    taus = np.exp(p[2:2 + nt])
    rest = p[2 + nt:]
    amps = np.concatenate([[-np.sum(rest)], rest])
    res = resid(p)
    order = np.argsort(taus)
    return ModelFit(float(p[0]), float(p[1]), taus[order], amps[order],
                    float(np.sqrt(np.mean(res ** 2))), float(np.max(np.abs(res))))
