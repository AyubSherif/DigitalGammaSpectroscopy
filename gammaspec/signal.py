"""Basic signal utilities: polarity, robust noise estimates and a coarse leading-edge trigger.

The coarse trigger is for characterisation (template building, finding pulse-free
noise windows). Production event detection will use a matched filter on the
baseline-restored signal (roadmap step 4).
"""

from __future__ import annotations

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import find_peaks

#: Pulses in the 2021 data are negative-going. Internally pulses are made positive.
DEFAULT_POLARITY = -1

_MAD_TO_SIGMA = 1.4826


def to_positive(volts: np.ndarray, polarity: int = DEFAULT_POLARITY) -> np.ndarray:
    """Return the signal with pulses pointing up."""
    if polarity not in (1, -1):
        raise ValueError("polarity must be +1 or -1")
    return volts if polarity == 1 else -volts


def robust_sigma(x: np.ndarray) -> float:
    """Gaussian-equivalent standard deviation from the median absolute deviation."""
    x = np.asarray(x)
    return float(_MAD_TO_SIGMA * np.median(np.abs(x - np.median(x))))


def noise_sigma(x: np.ndarray) -> float:
    """White-noise sigma from first differences. Insensitive to pulses and slow baseline changes."""
    return robust_sigma(np.diff(np.asarray(x, dtype=np.float64))) / np.sqrt(2.0)


def edge_signal(x: np.ndarray, gap: int = 2, avg: int = 3) -> np.ndarray:
    """d[n] = x[n] - mean(x[n-gap-avg+1 : n-gap+1]): the rise relative to a short pre-edge average.

    With a 1-2 sample rise, ``gap=2`` keeps a half-risen sample out of the reference.
    The first gap+avg-1 samples are set to 0.
    """
    x = np.asarray(x, dtype=np.float64)
    lag = gap + avg - 1
    d = np.zeros_like(x)
    ref = sliding_window_view(x, avg).mean(axis=1)  # ref[i] = mean(x[i:i+avg])
    d[lag:] = x[lag:] - ref[: x.size - lag]
    return d


def find_pulses(x: np.ndarray, sigma: float, threshold: float = 8.0,
                gap: int = 2, avg: int = 3, peak_search: int = 4,
                min_separation: int = 3) -> np.ndarray:
    """Coarse leading-edge trigger on a positive-pulse signal.

    Triggers where the edge signal exceeds ``threshold`` times its noise. Returns the
    index of the highest sample within ``peak_search`` samples after each trigger.
    """
    d = edge_signal(x, gap, avg)
    sigma_d = sigma * np.sqrt(1.0 + 1.0 / avg)
    trig, _ = find_peaks(d, height=threshold * sigma_d, distance=min_separation)
    if trig.size == 0:
        return trig
    trig = trig[trig + peak_search < x.size]
    win = sliding_window_view(x, peak_search + 1)[trig - 1 if peak_search else trig]
    peaks = trig - 1 + np.argmax(win, axis=1)
    return np.unique(peaks)


def isolated(peaks: np.ndarray, pre: int, post: int, n: int | None = None) -> np.ndarray:
    """Boolean mask: no other pulse within ``pre`` samples before or ``post`` after,
    and (if ``n`` is given) the window [p-pre, p+post] fits inside the capture."""
    peaks = np.asarray(peaks)
    ok = np.ones(peaks.size, dtype=bool)
    if peaks.size > 1:
        gaps = np.diff(peaks)
        ok[1:] &= gaps > pre
        ok[:-1] &= gaps > post
    if n is not None:
        ok &= (peaks - pre >= 0) & (peaks + post < n)
    return ok


def pulse_free_mask(n: int, peaks: np.ndarray, pre: int, post: int) -> np.ndarray:
    """Boolean mask of samples at least ``pre`` before and ``post`` after every pulse."""
    delta = np.zeros(n + 1, dtype=np.int32)
    lo = np.clip(np.asarray(peaks) - pre, 0, n)
    hi = np.clip(np.asarray(peaks) + post + 1, 0, n)
    np.add.at(delta, lo, 1)
    np.add.at(delta, hi, -1)
    return np.cumsum(delta[:-1]) == 0
