"""Pile-up resolution: estimate the amplitudes of overlapping pulses in one cluster.

All methods take the clean waveform of a cluster ``y`` (baseline 0, positive pulses), the
detected peak positions ``pos`` (samples, relative to the start of ``y``, sorted) and a
context with the pulse template. They return one amplitude per position.

* ``tm``    - template matching as in Mohammadian-Behbahani & Saramad (2020), after
              Vencelj et al.: measured amplitudes b_p = y[pos_p]; b = M a with
              M_pq = s(pos_p - pos_q) for q <= p (lower triangular, unit diagonal).
* ``hyper`` - HYPER (Wong & Li): assumes mono-exponential pulses e^(-lambda t).
              D_p = y(t) + lambda * integral_{pos_p}^{t} y is constant between pulses p and
              p+1; A_p = D_p - D_{p-1} e^(-lambda (pos_p - pos_{p-1})).
* ``ls``    - least squares (the paper's MLE for white noise): y = sum_k A_k s(t - t_k)
              over the whole cluster, linear in A_k. Positions as detected (integer).
* ``ls_refine``  - LS with joint sub-sample timing: Gauss-Newton on t_k using the template
              derivative (A s(t - t - d) ~ A s(t - t) - A d s'(t - t)).
* ``gls_refine`` - as ``ls_refine`` but weighted by the inverse noise covariance
              (prewhitening with a Cholesky factor): the maximum-likelihood estimate in
              coloured Gaussian noise *if the template were exact*.
* ``wls_refine`` - as ``gls_refine``, but the covariance also contains the template
              mismatch, which grows with the amplitude: C = C_noise + (eps * A_max)^2 I.
              Small pulses (noise-dominated) are weighted like GLS, large pulses
              (mismatch-dominated) like ordinary LS.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky, solve_triangular, toeplitz

from .template import PulseTemplate

METHODS = ("tm", "hyper", "ls", "ls_refine", "gls_refine", "wls_refine")


@dataclass
class PileupContext:
    template: PulseTemplate
    tau: float                       # main decay constant of the template (samples), for HYPER
    noise_acov: np.ndarray | None = None   # noise autocovariance r[0..], for GLS
    model_error: float = 0.0035      # template mismatch, fraction of amplitude (rms per sample), for WLS
    refine_iterations: int = 4
    max_shift: float = 0.75          # timing refinement limited to +-max_shift samples per pulse

    def derivative(self, t: np.ndarray, h: float = 1.0 / 16) -> np.ndarray:
        return (self.template(t + h) - self.template(t - h)) / (2 * h)


# ----------------------------------------------------------------------------------------------
# Template matching (paper)
# ----------------------------------------------------------------------------------------------

def tm(y: np.ndarray, pos: np.ndarray, ctx: PileupContext) -> np.ndarray:
    pos = np.asarray(pos, dtype=float)
    b = y[np.round(pos).astype(int)]
    d = pos[:, None] - pos[None, :]                       # pos_p - pos_q
    M = np.where(d >= 0, ctx.template(d), 0.0)            # only earlier pulses affect pulse p
    M[np.diag_indices_from(M)] = ctx.template(np.zeros(1))[0]
    return solve_triangular(M, b, lower=True)


# ----------------------------------------------------------------------------------------------
# HYPER
# ----------------------------------------------------------------------------------------------

def hyper(y: np.ndarray, pos: np.ndarray, ctx: PileupContext, max_len: int | None = None) -> np.ndarray:
    lam = 1.0 / ctx.tau
    p = np.round(np.asarray(pos)).astype(int)
    max_len = max_len or int(3 * ctx.tau)
    D = np.zeros(p.size)
    A = np.zeros(p.size)
    for k in range(p.size):
        a = p[k]
        b = min(p[k + 1] if k + 1 < p.size else y.size, a + max_len, y.size)
        seg = y[a:b]
        if seg.size < 2:
            D[k] = seg[0] if seg.size else 0.0
        else:
            # trapezoidal running integral from the peak; D(t) is ideally constant on the segment
            integ = np.concatenate([[0.0], np.cumsum(0.5 * (seg[1:] + seg[:-1]))])
            D[k] = float(np.mean(seg + lam * integ))
        A[k] = D[k] - (D[k - 1] * np.exp(-lam * (p[k] - p[k - 1])) if k else 0.0)
    return A


# ----------------------------------------------------------------------------------------------
# Least squares family
# ----------------------------------------------------------------------------------------------

def _whitener(ctx: PileupContext, n: int, extra_var: float = 0.0):
    """Cholesky factor of the noise covariance (Toeplitz, from ctx.noise_acov) plus
    ``extra_var`` on the diagonal (white model error)."""
    if ctx.noise_acov is None:
        raise ValueError("gls/wls need ctx.noise_acov")
    r = np.zeros(n)
    m = min(n, ctx.noise_acov.size)
    r[:m] = ctx.noise_acov[:m]
    C = toeplitz(r)
    C[np.diag_indices_from(C)] = C[np.diag_indices_from(C)] * (1.0 + 1e-6) + extra_var
    return cholesky(C, lower=True)


def _ls(y, pos, ctx, refine: bool, whiten: bool, extra_var: float = 0.0):
    t = np.arange(y.size, dtype=float)
    tk = np.asarray(pos, dtype=float).copy()
    L = _whitener(ctx, y.size, extra_var) if whiten else None

    def solve(X, target):
        if L is not None:
            X = solve_triangular(L, X, lower=True)
            target = solve_triangular(L, target, lower=True)
        coef, *_ = np.linalg.lstsq(X, target, rcond=None)
        return coef

    S = ctx.template(t[:, None] - tk[None, :])
    A = solve(S, y)
    if refine:
        t0 = tk.copy()
        for _ in range(ctx.refine_iterations):
            S = ctx.template(t[:, None] - tk[None, :])
            Dv = ctx.derivative(t[:, None] - tk[None, :])
            coef = solve(np.hstack([S, -Dv]), y)          # y = S A + (-s') (A delta)
            A, Ad = coef[: tk.size], coef[tk.size:]
            delta = np.where(np.abs(A) > 1e-12, Ad / np.where(np.abs(A) > 1e-12, A, 1.0), 0.0)
            tk = np.clip(tk + np.clip(delta, -0.5, 0.5), t0 - ctx.max_shift, t0 + ctx.max_shift)
        A = solve(ctx.template(t[:, None] - tk[None, :]), y)
    return A, tk


def ls(y, pos, ctx):
    return _ls(y, pos, ctx, refine=False, whiten=False)[0]


def ls_refine(y, pos, ctx):
    return _ls(y, pos, ctx, refine=True, whiten=False)[0]


def gls_refine(y, pos, ctx):
    return _ls(y, pos, ctx, refine=True, whiten=True)[0]


def wls_refine(y, pos, ctx):
    A0, _ = _ls(y, pos, ctx, refine=False, whiten=False)
    extra = (ctx.model_error * float(np.max(np.abs(A0)))) ** 2
    return _ls(y, pos, ctx, refine=True, whiten=True, extra_var=extra)[0]


FUNCS = {"tm": tm, "hyper": hyper, "ls": ls, "ls_refine": ls_refine, "gls_refine": gls_refine,
         "wls_refine": wls_refine}


def resolve(method: str, y: np.ndarray, pos: np.ndarray, ctx: PileupContext) -> np.ndarray:
    if method not in FUNCS:
        raise ValueError(f"unknown pile-up method {method!r}; choose from {METHODS}")
    return FUNCS[method](np.asarray(y, dtype=float), np.asarray(pos), ctx)
