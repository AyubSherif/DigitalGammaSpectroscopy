import numpy as np
import pytest

from gammaspec.noise import noise_spectrum, pulse_free_windows
from gammaspec.signal import find_pulses, isolated, noise_sigma, pulse_free_mask, robust_sigma
from gammaspec.template import (PulseTemplate, build_template, extract_windows, fit_multiexp,
                                fit_pulses, multiexp_model)

# A pulse shape like the real detector's (samples): fast rise, 7.4-sample decay, small undershoot
TRUE = dict(t0=-1.1, sigma=0.29, taus=(0.4, 7.4, 115.0), amps_rest=(1.34, -0.06))


def true_shape(t):
    s = multiexp_model(t, TRUE["t0"], TRUE["sigma"], TRUE["taus"], np.array(TRUE["amps_rest"]))
    tt = np.linspace(-3, 5, 8001)
    peak = multiexp_model(tt, TRUE["t0"], TRUE["sigma"], TRUE["taus"], np.array(TRUE["amps_rest"]))
    return s / peak.max(), tt[np.argmax(peak)]


def synthetic_trace(n=400_000, n_pulses=300, noise=0.002, seed=0):
    """Isolated pulses at random sub-sample times, amplitudes 0.05-0.5 V, white noise."""
    rng = np.random.default_rng(seed)
    centres = np.sort(rng.choice(np.arange(2000, n - 2000, 1200), n_pulses, replace=False)
                      + rng.uniform(0, 1, n_pulses))
    amps = rng.uniform(0.05, 0.5, n_pulses)
    x = rng.normal(0, noise, n)
    _, tpk = true_shape(np.zeros(1))
    for c, a in zip(centres, amps):
        idx = np.arange(int(c) - 50, int(c) + 900)
        s, _ = true_shape(idx - c + tpk)  # continuous peak exactly at c
        x[idx] += a * s
    return x, centres, amps


def test_noise_sigma_and_trigger():
    x, centres, amps = synthetic_trace()
    sig = noise_sigma(x)
    assert sig == pytest.approx(0.002, rel=0.05)
    assert robust_sigma(np.random.default_rng(1).normal(0, 3, 100_000)) == pytest.approx(3, rel=0.02)
    peaks = find_pulses(x, sig, threshold=8)
    assert peaks.size == centres.size
    assert np.all(np.abs(peaks - centres) < 1.5)


def test_isolation_and_free_mask():
    peaks = np.array([100, 150, 1000])
    np.testing.assert_array_equal(isolated(peaks, pre=40, post=40), [True, True, True])
    np.testing.assert_array_equal(isolated(peaks, pre=60, post=60), [False, False, True])
    np.testing.assert_array_equal(isolated(peaks, pre=60, post=60, n=1030), [False, False, False])
    m = pulse_free_mask(1200, peaks, pre=5, post=10)
    assert m[94] and not m[95:111].any() and m[111:145].all() and not m[145:161].any() and m[161]
    assert not m[995:1011].any() and m[1011]


def test_build_template_and_fit_recover_truth():
    x, centres, amps = synthetic_trace()
    peaks = find_pulses(x, noise_sigma(x), threshold=8)
    pre, post = 40, 800
    w = extract_windows(x, peaks, pre, post)
    tmpl, fit, used = build_template(w, pre, resolution=16)
    t = np.linspace(-1.5, 300, 5000)
    truth, _ = true_shape(t + true_shape(np.zeros(1))[1])
    # 300 pulses: the residual error is statistical and shrinks with more pulses
    assert np.max(np.abs(tmpl(t) - truth)) < 0.035         # shape within 3.5 % of peak, worst point
    assert np.sqrt(np.mean((tmpl(t) - truth) ** 2)) < 0.003
    assert used.sum() >= 0.95 * peaks.size
    np.testing.assert_allclose(fit.amplitude, amps, rtol=0.02, atol=0.002)
    true_delta = centres - peaks
    assert np.std(fit.shift - true_delta) < 0.035            # sub-sample timing to ~1/30 sample


def test_phase_equalisation_removes_warp():
    from gammaspec.template import equalize_phase
    rng = np.random.default_rng(0)
    true = rng.uniform(-0.5, 0.5, 5000)
    warped = true + 0.04 * np.sin(2 * np.pi * true)          # phase-dependent bias
    fixed = equalize_phase(warped)
    # equalisation is defined up to a constant offset (fixed later by recentring)
    resid = fixed - true
    assert np.std(resid - resid.mean()) < 0.003 < np.std(warped - true)


def test_fit_pulses_with_exact_template():
    t = np.arange(-40, 801) / 1.0
    tt = np.arange(-40 * 16, 800 * 16 + 1) / 16
    s, tpk = true_shape(tt + true_shape(np.zeros(1))[1])
    tmpl = PulseTemplate(16, float(tt[0]), s)
    rng = np.random.default_rng(3)
    deltas = rng.uniform(-0.5, 0.5, 200)
    A = rng.uniform(0.1, 1, 200)
    w = A[:, None] * tmpl(t[None, :] - deltas[:, None]) + 0.01 + rng.normal(0, 1e-4, (200, t.size))
    f = fit_pulses(w, tmpl, pre=40)
    np.testing.assert_allclose(f.amplitude, A, rtol=1e-3)
    np.testing.assert_allclose(f.shift, deltas, atol=5e-3)
    np.testing.assert_allclose(f.baseline, 0.01, atol=1e-4)


def test_multiexp_fit_recovers_time_constants():
    tt = np.arange(-30 * 16, 1000 * 16 + 1) / 16
    s, tpk = true_shape(tt + true_shape(np.zeros(1))[1])
    m = fit_multiexp(PulseTemplate(16, float(tt[0]), s))
    np.testing.assert_allclose(m.taus, TRUE["taus"], rtol=0.03)
    assert m.rms < 1e-3


def test_noise_spectrum_white_and_line():
    fs, n, sigma = 1.8e6, 512, 0.003
    rng = np.random.default_rng(0)
    w = rng.normal(0, sigma, (2000, n))
    ns = noise_spectrum(w, fs)
    assert ns.sigma_v == pytest.approx(sigma, rel=0.02)
    assert ns.white_level == pytest.approx(2 * sigma ** 2 / fs, rel=0.1)  # one-sided PSD
    assert abs(ns.autocorr[1]) < 0.02
    # add a 257 kHz line: it must stand out
    tone = 0.01 * np.sin(2 * np.pi * 257e3 * np.arange(n) / fs + rng.uniform(0, 6.3, (2000, 1)))
    ns2 = noise_spectrum(w + tone, fs)
    assert abs(ns2.freq_hz[np.argmax(ns2.psd)] - 257e3) < fs / n


def test_pulse_free_windows():
    x = np.arange(100.0)
    free = np.ones(100, bool)
    free[30:40] = False
    win, starts = pulse_free_windows(x, free, 10)
    assert list(starts) == [0, 10, 20, 40, 50, 60, 70, 80, 90]
    np.testing.assert_array_equal(win[3], np.arange(40, 50))
