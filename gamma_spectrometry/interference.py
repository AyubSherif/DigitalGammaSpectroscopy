"""Interference removal: detect narrow spectral lines and subtract them as slowly varying sinusoids.

The 2021 data carries narrow-band pickup: a family of lines at multiples of 9 kHz
(= fs/200, from the acquisition electronics) in both runs, plus an intermittent,
slightly frequency-modulated line near 257 kHz in the background run.

Method (per capture or chunk):

1. **Detect** lines in the full-length power spectrum of the pulse-free samples
   (pulses masked to zero): local maxima standing ``ratio_threshold`` times above the
   local median. Refine each frequency by parabolic interpolation, then merge lines
   closer than one block's frequency resolution (the block-wise fit cannot separate
   them, and block-wise amplitude/phase absorbs small frequency wobble).
2. **Subtract**: in overlapping blocks (sin² windows, 50 % overlap) fit cos/sin
   amplitudes of every line by least squares *using only pulse-free samples*, and
   subtract the overlap-added model from every sample. Pulses therefore never shape the
   fit and are not distorted; inside a pulse the interference is interpolated from its
   surroundings.

The masked least squares is solved by filling masked samples with the current model
and refitting with one fixed design matrix for all blocks (an EM iteration). This is
exact at convergence and needs only matrix products.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from functools import lru_cache

import numpy as np
import scipy.fft as sfft
from scipy.ndimage import maximum_filter, median_filter


@lru_cache(maxsize=4)
def _hann(n: int) -> np.ndarray:
    return np.hanning(n)


@dataclass
class LineRemovalConfig:
    block: int = 1024                # samples per fit block (shorter tracks drifting lines better)
    ratio_threshold: float = 100.0   # line must exceed the local spectral median by this factor
    max_lines: int = 40
    merge_resolutions: float = 2.0   # lines closer than this many block resolutions (fs/block) merge
    min_freq_hz: float = 2_000.0     # ignore lines below this (pulse tails live there)
    median_width: int = 1001         # FFT bins for the local-median spectral floor
    iterations: int = 3              # fill-and-refit passes for masked samples
    ridge: float = 1e-6              # relative regularisation of the normal equations


@dataclass
class LineRemovalInfo:
    freqs_hz: np.ndarray = field(default_factory=lambda: np.empty(0))
    removed_rms_v: float = 0.0       # rms of the subtracted interference


def detect_lines(x: np.ndarray, free: np.ndarray, fs: float, cfg: LineRemovalConfig) -> np.ndarray:
    """Frequencies (Hz) of narrow lines in the pulse-free part of ``x``, merged to block resolution."""
    n = x.size
    xm = np.where(free, x - x[free].mean(), 0.0)
    P = np.abs(sfft.rfft(xm * _hann(n), workers=-1)) ** 2
    f = np.fft.rfftfreq(n, 1.0 / fs)
    floor = median_filter(P, size=cfg.median_width, mode="nearest")
    ratio = P / np.maximum(floor, np.finfo(float).tiny)
    # A line must be the strongest peak within +-merge distance. This rejects the
    # sidebands that masking out pulses can create around strong lines; such
    # near-duplicate frequencies would leave the fit unconstrained inside pulse gaps.
    merge = cfg.merge_resolutions * fs / cfg.block
    half = max(4, int(np.ceil(merge / (fs / n))))
    ok = (ratio > cfg.ratio_threshold) & (P == maximum_filter(P, size=2 * half + 1)) \
        & (f > cfg.min_freq_hz) & (f < fs / 2 - cfg.min_freq_hz)
    c = np.flatnonzero(ok)
    if c.size == 0:
        return np.empty(0)
    c = c[np.argsort(P[c])[::-1]]                     # strongest first
    lp = np.log(P)
    den = lp[c - 1] - 2 * lp[c] + lp[c + 1]
    d = np.where(den < 0, 0.5 * (lp[c - 1] - lp[c + 1]) / np.where(den < 0, den, -1), 0.0)
    freqs = (c + np.clip(d, -0.5, 0.5)) * fs / n
    # merge lines the block fit cannot resolve; keep the strongest of each group
    kept: list[float] = []
    for fr in freqs:
        if all(abs(fr - k) > merge for k in kept):
            kept.append(fr)
        if len(kept) == cfg.max_lines:
            break
    return np.sort(np.array(kept))


def subtract_lines(x: np.ndarray, free: np.ndarray, freqs: np.ndarray, fs: float,
                   cfg: LineRemovalConfig) -> np.ndarray:
    """The interference model (same length as x) for the given line frequencies."""
    n, L = x.size, cfg.block
    if freqs.size == 0:
        return np.zeros(n)
    hop = L // 2
    starts = np.arange(-hop, n, hop)
    nb = starts.size
    # one design matrix in block-local time; each block has its own coefficients
    tau = np.arange(L)
    ph = 2 * np.pi * np.outer(tau, freqs) / fs
    X = np.hstack([np.cos(ph), np.sin(ph)])           # (L, 2K)
    G = X.T @ X
    G += cfg.ridge * np.trace(G) / G.shape[0] * np.eye(G.shape[0])
    Ginv_Xt = np.linalg.solve(G, X.T)                 # (2K, L)
    win = np.sin(np.pi * (tau + 0.5) / L) ** 2        # sums to 1 at 50 % overlap

    # padded copies so every block is full length
    pad_lo, pad_hi = hop, (starts[-1] + L) - n
    xp = np.pad(x, (pad_lo, max(pad_hi, 0)))
    fp = np.pad(free, (pad_lo, max(pad_hi, 0)))
    idx = (starts + pad_lo)[:, None] + tau[None, :]   # (nb, L)
    Y, M = xp[idx], fp[idx]
    model_blocks = np.zeros_like(Y)
    for _ in range(max(1, cfg.iterations)):
        Yf = np.where(M, Y, model_blocks)              # fill masked samples with the model
        coef = Yf @ Ginv_Xt.T                          # (nb, 2K)
        model_blocks = coef @ X.T                      # (nb, L)
    # overlap-add with the sin² window
    out = np.zeros(xp.size)
    np.add.at(out, idx.ravel(), (model_blocks * win[None, :]).ravel())
    return out[pad_lo: pad_lo + n]


class LineRemoval:
    """Optional pipeline stage: remove narrow-band interference.

    ``process(x, free)`` takes a positive-pulse signal and a mask of pulse-free samples,
    and returns the cleaned signal. Lines are re-detected for every chunk (they switch
    on and off and drift between captures). Details of the last call are in ``info``.
    """

    name = "line_removal"

    def __init__(self, fs: float, config: LineRemovalConfig | None = None):
        self.fs = fs
        self.cfg = config or LineRemovalConfig()
        self.info = LineRemovalInfo()

    def process(self, x: np.ndarray, free: np.ndarray) -> np.ndarray:
        freqs = detect_lines(x, free, self.fs, self.cfg)
        model = subtract_lines(x, free, freqs, self.fs, self.cfg)
        self.info = LineRemovalInfo(freqs, float(np.sqrt(np.mean(model ** 2))))
        return x - model

    def reset(self) -> None:
        self.info = LineRemovalInfo()
