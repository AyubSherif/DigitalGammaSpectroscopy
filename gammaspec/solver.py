"""Pile-up solver (step 6): weighted least squares with sub-sample timing and a residual search.

For one cluster (clean waveform y, detected positions):

1. **Fit** y(t) = sum_k A_k s(t - t_k) by least squares weighted with
   C = C_noise + (eps * A_max)^2 I (noise plus amplitude-proportional template mismatch,
   ADR-0008). The arrival times t_k are refined jointly by Gauss-Newton
   (A s(t - t_k - d) ~ A s(t - t_k) - (A d) s'(t - t_k)), within +-max_shift samples.
2. **Residual search** for pulses the detector merged (closer than 3 samples) or missed:
   for every candidate time tau, the whitened residual r_w and whitened template s_w(tau)
   give the amplitude and significance of one added pulse,
       A* = s_w.r_w / s_w.s_w,   z = s_w.r_w / sqrt(s_w.s_w)   (noise sigma = 1 after whitening),
   z being the improvement in chi^2 expressed in sigma. The best candidate is added if
   z > ``search_threshold``, A* > 0 and it is at least ``min_separation`` samples from the
   existing pulses; then everything is refitted, and the new pulse is kept only if it is
   still significant. Repeated up to ``max_added`` times.
3. **Quality**: rms of the (unweighted) residual and the whitened chi^2 per sample.

The covariance is factorised once per cluster (Cholesky), so each iteration costs only
triangular solves.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky, solve_triangular, toeplitz

from .template import PulseTemplate


@dataclass
class SolverConfig:
    enabled: bool = True               # stage: resolve pile-up clusters (needs pole_zero + baseline_tracker)
    model_error: float = 0.0035        # template mismatch, fraction of amplitude (measured: 0.35 %)
    iterations: int = 6                # max Gauss-Newton iterations (stops earlier on convergence)
    tolerance: float = 0.005           # stop when every time moved less than this (samples)
    max_shift: float = 1.25            # each time stays within +-max_shift of its start value
    residual_search: bool = True       # search the fit residual for hidden pulses
    search_threshold: float = 8.0      # z (sigma) needed to add a pulse (tuned: 0.3 % false additions)
    min_separation: float = 1.0        # added pulses: at least this far (samples) from existing ones
    max_added: int = 2                 # at most this many pulses added per cluster
    max_length: int = 1024             # longer clusters are fitted without whitening
    acov_lags: int = 512               # noise autocovariance lags used for whitening


@dataclass
class SolverResult:
    amplitude: np.ndarray              # V, one per pulse (original pulses first, in input order)
    time: np.ndarray                   # samples, relative to the start of y
    added: np.ndarray                  # bool: found by the residual search
    rms: float                         # rms of y - model (V)
    chi2: float                        # whitened chi^2 per sample (~1 for a good fit)
    whitened: bool


def taper_autocovariance(acov: np.ndarray, lags: int) -> np.ndarray:
    """First ``lags`` values of an autocovariance, smoothly tapered to zero (Tukey half-window).

    Truncating an estimated autocovariance abruptly can make the Toeplitz matrix indefinite;
    the taper keeps it positive (semi)definite in practice.
    """
    r = np.zeros(lags)
    m = min(lags, acov.size)
    r[:m] = acov[:m]
    k = np.arange(lags)
    edge = int(0.75 * lags)
    w = np.where(k < edge, 1.0, 0.5 * (1 + np.cos(np.pi * (k - edge) / max(lags - edge, 1))))
    return r * w


class PileupSolver:
    """Resolve clusters. Construct per capture (the noise autocovariance can change)."""

    name = "pileup_solver"

    def __init__(self, template: PulseTemplate, noise_acov: np.ndarray | None,
                 config: SolverConfig | None = None):
        self.template = template
        self.cfg = config or SolverConfig()
        self.acov = None if noise_acov is None else taper_autocovariance(noise_acov, self.cfg.acov_lags)
        self.noise_var = None if noise_acov is None else float(noise_acov[0])

    # -- helpers -----------------------------------------------------------------------------
    def _s(self, x):
        return self.template(x)

    def _ds(self, x):
        return self.template.derivative(x)

    def _whitener(self, n: int, extra_var: float):
        """Whitening matrix W = L^-1 (C = L L^T, n x n), or None if whitening is not possible.

        W is formed once per cluster, so every later whitening is one matrix product."""
        if self.acov is None or n > self.cfg.max_length:
            return None
        r = np.zeros(n)
        m = min(n, self.acov.size)
        r[:m] = self.acov[:m]
        C = toeplitz(r)
        C[np.diag_indices_from(C)] += extra_var + 1e-6 * r[0]
        try:
            L = cholesky(C, lower=True, check_finite=False)
        except np.linalg.LinAlgError:
            return None
        return solve_triangular(L, np.eye(n), lower=True, check_finite=False)

    # -- core ----------------------------------------------------------------------------------
    def _fit(self, y, t0, L, yw):
        """Weighted LS with Gauss-Newton timing; returns (A, t). L is the whitening matrix (or None)."""
        c = self.cfg
        n = y.size
        tt = np.arange(n, dtype=float)[:, None]
        t = np.asarray(t0, dtype=float).copy()
        k = t.size

        def wsolve(X):
            # normal equations of the (whitened) problem: small (k x k), fast
            Xw = X if L is None else L @ X
            G = Xw.T @ Xw
            try:
                return np.linalg.solve(G + 1e-12 * np.trace(G) * np.eye(G.shape[0]), Xw.T @ yw)
            except np.linalg.LinAlgError:
                return np.linalg.lstsq(Xw, yw, rcond=None)[0]

        for _ in range(c.iterations):
            S = self._s(tt - t)
            D = self._ds(tt - t)
            coef = wsolve(np.hstack([S, -D]))
            A, Ad = coef[:k], coef[k:]
            delta = np.where(np.abs(A) > 1e-12, Ad / np.where(np.abs(A) > 1e-12, A, 1.0), 0.0)
            delta = np.clip(delta, -0.5, 0.5)
            t_new = np.clip(t + delta, t0 - c.max_shift, t0 + c.max_shift)
            moved = np.max(np.abs(t_new - t)) if k else 0.0
            t = t_new
            if moved < c.tolerance:
                break
        A = wsolve(self._s(tt - t))
        return A, t

    def solve(self, y: np.ndarray, positions) -> SolverResult:
        c = self.cfg
        y = np.asarray(y, dtype=float)
        n = y.size
        tt = np.arange(n, dtype=float)
        t0 = np.asarray(positions, dtype=float)
        # first pass without weights gives A_max for the mismatch term
        S0 = self._s(tt[:, None] - t0[None, :])
        A0, *_ = np.linalg.lstsq(S0, y, rcond=None)
        extra = (c.model_error * float(np.max(np.abs(A0)) if A0.size else 0.0)) ** 2
        L = self._whitener(n, extra)
        yw = y if L is None else L @ y
        A, t = self._fit(y, t0, L, yw)
        added = np.zeros(t.size, bool)

        if c.residual_search and L is not None:
            cand = np.arange(max(0.0, t.min() - 5), min(n - 2.0, t.max() + 15) + 1e-9, 0.5)
            for _ in range(c.max_added):
                model = self._s(tt[:, None] - t[None, :]) @ A
                rw = L @ (y - model)
                Sc = L @ self._s(tt[:, None] - cand[None, :])
                num = Sc.T @ rw
                den = np.einsum("ij,ij->j", Sc, Sc)
                z = num / np.sqrt(np.maximum(den, 1e-30))
                far = np.min(np.abs(cand[:, None] - t[None, :]), axis=1) >= c.min_separation
                z = np.where(far & (num > 0), z, -np.inf)
                j = int(np.argmax(z))
                if not np.isfinite(z[j]) or z[j] < c.search_threshold:
                    break
                t_try0 = np.r_[t, cand[j]]
                A_try, t_try = self._fit(y, t_try0, L, yw)
                # keep the new pulse only if it is still significant after the joint refit
                Sw_new = L @ self._s(tt - t_try[-1])
                sigma_new = 1.0 / np.sqrt(np.dot(Sw_new, Sw_new))      # amplitude sigma (whitened units)
                if A_try[-1] <= 0 or A_try[-1] / sigma_new < c.search_threshold:
                    break
                A, t = A_try, t_try
                added = np.r_[added, True]

        model = self._s(tt[:, None] - t[None, :]) @ A
        res = y - model
        rw = res if L is None else L @ res
        chi2 = float(np.mean(rw ** 2)) if L is not None else float(np.mean(res ** 2) / (self.noise_var or 1.0))
        return SolverResult(A, t, added, float(np.sqrt(np.mean(res ** 2))), chi2, L is not None)
