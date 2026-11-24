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
    area_errs: np.ndarray = None

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
              sigma0=None, center_tol: float | None = None, sigma_range=None) -> PeakFit:
    """Fit N Gaussians on a linear background to the histogram of ``values`` in [lo, hi).

    ``centers`` are starting guesses. Poisson weights; the background may be negative
    in slope. Returns centres, sigmas, areas (counts) and their errors.

    Optional bounds (used for calibration lines): ``center_tol`` keeps each centre within
    that fraction of its guess; ``sigma_range=(lo, hi)`` keeps each sigma within that
    factor range of its ``sigma0``; with any bound set, areas are kept >= 0. ``sigma0``
    may be one value or one per peak.
    """
    edges = np.arange(lo, hi + bin_width / 2, bin_width)
    h, _ = np.histogram(values, edges)
    x = 0.5 * (edges[1:] + edges[:-1])
    centers = np.asarray(centers, dtype=float)
    s0 = np.broadcast_to(np.asarray(sigma0 if sigma0 is not None else 0.03 * float(np.mean(centers)),
                                    dtype=float), centers.shape)
    p0, lb, ub = [], [], []
    for c, s in zip(centers, s0):
        near = np.abs(x - c) < s
        peak = h[near].max() if near.any() else h.max()
        p0 += [max(peak, 1) * np.sqrt(2 * np.pi) * s, c, s]
        ct = center_tol if center_tol is not None else np.inf
        sr = sigma_range if sigma_range is not None else (0.0, np.inf)
        lb += [0.0 if (center_tol is not None or sigma_range is not None) else -np.inf, c * (1 - ct), s * sr[0]]
        ub += [np.inf, c * (1 + ct), s * sr[1]]
    p0 += [float(np.median(h)), 0.0]
    lb += [-np.inf, -np.inf]
    ub += [np.inf, np.inf]
    err = np.sqrt(np.maximum(h, 1))
    bounded = center_tol is not None or sigma_range is not None
    if bounded:
        p0 = list(np.clip(p0, np.nextafter(lb, np.inf), np.nextafter(ub, -np.inf)))
        popt, pcov = curve_fit(_model, x, h, p0=p0, sigma=err, absolute_sigma=True, maxfev=20000,
                               bounds=(lb, ub))
    else:
        popt, pcov = curve_fit(_model, x, h, p0=p0, sigma=err, absolute_sigma=True, maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    n = len(centers)
    chi2 = float(np.sum(((h - _model(x, *popt)) / err) ** 2) / max(1, x.size - len(popt)))
    return PeakFit(popt[1:3 * n:3], np.abs(popt[2:3 * n:3]), popt[0:3 * n:3] / bin_width,
                   perr[1:3 * n:3], perr[2:3 * n:3], chi2, perr[0:3 * n:3] / bin_width)
