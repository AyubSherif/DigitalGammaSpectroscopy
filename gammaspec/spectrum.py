"""Spectrum utilities: histogramming and Gaussian peak fits (resolution measurement)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit

FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


@dataclass
class PeakFit:
    centers: np.ndarray
    sigmas: np.ndarray
    areas: np.ndarray
    center_errs: np.ndarray
    sigma_errs: np.ndarray
    chi2_ndf: float

    @property
    def fwhm(self) -> np.ndarray:
        return FWHM_PER_SIGMA * self.sigmas

    @property
    def resolution(self) -> np.ndarray:
        """FWHM / centre (fraction)."""
        return self.fwhm / self.centers

    @property
    def resolution_err(self) -> np.ndarray:
        return self.resolution * np.hypot(self.sigma_errs / self.sigmas, self.center_errs / self.centers)


def _model(x, *p):
    n = (len(p) - 2) // 3
    y = p[-2] + p[-1] * x
    for k in range(n):
        area, mu, sig = p[3 * k: 3 * k + 3]
        y = y + area / (np.sqrt(2 * np.pi) * sig) * np.exp(-0.5 * ((x - mu) / sig) ** 2)
    return y


def fit_peaks(values: np.ndarray, centers, lo: float, hi: float, bin_width: float,
              sigma0: float | None = None) -> PeakFit:
    """Fit N Gaussians on a linear background to the histogram of ``values`` in [lo, hi).

    ``centers`` are starting guesses. Poisson weights; the background may be negative
    in slope. Returns centres, sigmas, areas (counts) and their errors.
    """
    edges = np.arange(lo, hi + bin_width / 2, bin_width)
    h, _ = np.histogram(values, edges)
    x = 0.5 * (edges[1:] + edges[:-1])
    centers = np.asarray(centers, dtype=float)
    sigma0 = sigma0 or 0.03 * float(np.mean(centers))
    p0 = []
    for c in centers:
        near = np.abs(x - c) < sigma0
        peak = h[near].max() if near.any() else h.max()
        p0 += [max(peak, 1) * np.sqrt(2 * np.pi) * sigma0, c, sigma0]
    p0 += [float(np.median(h)), 0.0]
    err = np.sqrt(np.maximum(h, 1))
    popt, pcov = curve_fit(_model, x, h, p0=p0, sigma=err, absolute_sigma=True, maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    n = len(centers)
    chi2 = float(np.sum(((h - _model(x, *popt)) / err) ** 2) / max(1, x.size - len(popt)))
    return PeakFit(popt[1:3 * n:3], np.abs(popt[2:3 * n:3]), popt[0:3 * n:3] / bin_width,
                   perr[1:3 * n:3], perr[2:3 * n:3], chi2)
