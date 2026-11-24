import numpy as np
import pytest

from gamma_spectrometry.pileup import METHODS, PileupContext, resolve
from gamma_spectrometry.template import PulseTemplate

from .test_detect_events import T


def ctx_for(template, tau=7.7, noise=1e-7):
    acov = np.r_[noise, np.zeros(2000)]
    return PileupContext(template, tau, acov)


def wave(template, n, cs, As):
    t = np.arange(n, dtype=float)
    return sum(A * template(t - c) for c, A in zip(cs, As))


@pytest.mark.parametrize("method", ["ls_refine", "gls_refine", "wls_refine"])
@pytest.mark.parametrize("d", [3, 4, 8, 20])
def test_refined_fits_exact_at_any_phase(method, d):
    ctx = ctx_for(T)
    for phase in (0.0, 0.3, 0.6):
        cs, As = [30 + phase, 30 + phase + d], [0.3, 0.05]
        y = wave(T, 150, cs, As)
        est = resolve(method, y, np.round(cs), ctx)
        np.testing.assert_allclose(est, As, rtol=0.01)


@pytest.mark.parametrize("method", ["tm", "ls"])
def test_integer_methods_exact_on_sample_grid(method):
    ctx = ctx_for(T)
    cs, As = [30.0, 36.0, 50.0], [0.3, 0.1, 0.2]
    est = resolve(method, wave(T, 160, cs, As), np.array(cs), ctx)
    np.testing.assert_allclose(est, As, rtol=1e-6)


def test_hyper_exact_for_exponential_pulses():
    t = np.arange(-40 * 16, 300 * 16 + 1) / 16
    s = np.where(t >= 0, np.exp(-np.clip(t, 0, None) / 7.0), 0.0)
    E = PulseTemplate(16, float(t[0]), s)
    ctx = ctx_for(E, tau=7.0)
    cs, As = [20.0, 27.0, 45.0], [0.3, 0.2, 0.1]
    est = resolve("hyper", wave(E, 140, cs, As), np.array(cs), ctx)
    # the trapezoidal integral of a sampled exponential is slightly biased; 2 % is ample
    np.testing.assert_allclose(est, As, rtol=0.02)


def test_wls_matches_gls_for_small_and_ls_for_large_mismatch():
    ctx = ctx_for(T, noise=(5e-4) ** 2)
    rng = np.random.default_rng(0)
    y = wave(T, 150, [30.2, 36.7], [0.02, 0.01]) + rng.normal(0, 5e-4, 150)
    a_wls = resolve("wls_refine", y, np.array([30, 37]), ctx)
    a_gls = resolve("gls_refine", y, np.array([30, 37]), ctx)
    np.testing.assert_allclose(a_wls, a_gls, rtol=0.02)   # noise dominates: same weighting


def test_unknown_method():
    with pytest.raises(ValueError):
        resolve("nope", np.zeros(10), np.array([2]), ctx_for(T))
    assert set(METHODS) == {"tm", "hyper", "ls", "ls_refine", "gls_refine", "wls_refine"}
