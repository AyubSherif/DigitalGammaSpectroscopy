"""Noise characterisation from pulse-free stretches of the waveform."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def pulse_free_windows(x: np.ndarray, free: np.ndarray, length: int) -> np.ndarray:
    """Cut non-overlapping windows of ``length`` samples lying entirely in the ``free`` mask.

    Returns (n_windows, length) and the start index of each window.
    """
    free = np.asarray(free, dtype=bool)
    padded = np.concatenate([[False], free, [False]])
    edges = np.flatnonzero(np.diff(padded.astype(np.int8)))
    starts, stops = edges[::2], edges[1::2]
    win_starts = [np.arange(s, e - length + 1, length) for s, e in zip(starts, stops) if e - s >= length]
    if not win_starts:
        return np.empty((0, length)), np.empty(0, dtype=np.int64)
    ws = np.concatenate(win_starts)
    return x[ws[:, None] + np.arange(length)[None, :]], ws


@dataclass
class NoiseSpectrum:
    freq_hz: np.ndarray       # one-sided frequencies
    psd: np.ndarray           # one-sided power spectral density, V²/Hz
    sigma_v: float            # rms noise within a window (window mean removed)
    autocorr: np.ndarray      # normalised autocorrelation, lags 0..length-1 (samples)
    n_windows: int
    window_length: int

    @property
    def white_level(self) -> float:
        """Median PSD over the upper half of the band: the white-noise floor, V²/Hz."""
        return float(np.median(self.psd[self.psd.size // 2:]))


def noise_spectrum(windows: np.ndarray, fs: float) -> NoiseSpectrum:
    """Averaged periodogram (Hann window, mean removed) of pulse-free windows."""
    w = np.asarray(windows, dtype=np.float64)
    w = w - w.mean(axis=1, keepdims=True)
    n = w.shape[1]
    win = np.hanning(n)
    spec = np.fft.rfft(w * win, axis=1)
    psd = (np.abs(spec) ** 2).mean(axis=0) / (fs * np.sum(win ** 2))
    psd[1:-1] *= 2.0  # one-sided
    freq = np.fft.rfftfreq(n, 1.0 / fs)
    # autocorrelation from the unwindowed power spectrum (zero-padded: linear, not circular)
    full = np.fft.rfft(w, 2 * n, axis=1)
    ac = np.fft.irfft((np.abs(full) ** 2).mean(axis=0))[:n] / (n - np.arange(n))  # unbiased
    ac /= ac[0]
    return NoiseSpectrum(freq, psd, float(np.sqrt(np.mean(w ** 2))), ac, w.shape[0], n)
