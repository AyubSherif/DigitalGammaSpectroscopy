import warnings

import numpy as np
import pytest
from scipy.signal import lfilter, residuez

from gamma_spectrometry.baseline import (PoleZero, TrackerConfig, design_pole_zero, filter_template,
                                track_baseline)
from gamma_spectrometry.interference import LineRemoval
from gamma_spectrometry.preprocess import PreprocessConfig, Preprocessor
from gamma_spectrometry.signal import find_pulses, noise_sigma, pulse_free_mask
from gamma_spectrometry.spectrum import fit_peaks
from gamma_spectrometry.template import PulseTemplate, extract_windows, fit_pulses

FS = 1.8e6


def uncoupled_pulse(n=4000, tau_s=0.4, tau_d=7.4):
    k = np.arange(n)
    return np.exp(-k / tau_d) - np.exp(-k / tau_s)


def coupled(x, tau_c):
    """Digital CR coupling: zero at DC, pole at exp(-1/tau_c)."""
    ac = np.exp(-1 / tau_c)
    return lfilter([1, -1], [1, -ac], x)


def test_pole_zero_design_recovers_uncoupled_pulse():
    x = uncoupled_pulse()
    y = coupled(x, 80.0)
    assert y.min() < -0.02 * y.max()                          # visible undershoot
    # exact discrete model of y as a sum of exponentials:
    # X(z) = 1/(1-a_d z^-1) - 1/(1-a_s z^-1) = (a_d - a_s) z^-1 / ((1-a_d z^-1)(1-a_s z^-1))
    num = np.convolve([1, -1], np.array([0.0, np.exp(-1 / 7.4) - np.exp(-1 / 0.4)]))
    den = np.poly([np.exp(-1 / 0.4), np.exp(-1 / 7.4), np.exp(-1 / 80.0)])
    r, p, _ = residuez(num, den)
    taus = -1 / np.log(np.real(p))
    np.testing.assert_allclose(lfilter(num, den, np.r_[1, np.zeros(3999)]), y, atol=1e-12)
    d = design_pole_zero(taus, np.real(r), leak_samples=1e6)
    z = lfilter(d.b, d.a, y)
    assert np.max(np.abs(z - x)) < 2e-3 * x.max()             # undershoot gone, shape kept
    assert np.max(np.abs(z[100:])) < 1e-3 * x.max()


def test_pole_zero_stage_is_streamable():
    rng = np.random.default_rng(0)
    y = coupled(rng.normal(size=10_000), 80.0)
    d = design_pole_zero([0.4, 7.4, 80.0], [-1.0, 1.05, -0.05])
    whole = PoleZero(d).process(y)
    st = PoleZero(d)
    parts = np.concatenate([st.process(c) for c in np.array_split(y, 7)])
    np.testing.assert_allclose(parts, whole, atol=1e-12)


def test_filter_template_identity():
    t = np.arange(-40 * 16, 300 * 16 + 1) / 16
    s = np.exp(-np.clip(t, 0, None) / 7.4) * (t >= 0)
    T = PulseTemplate(16, float(t[0]), s / s.max())
    U = filter_template(T, [1.0], [1.0])
    np.testing.assert_allclose(U(T.t), T.shape, atol=1e-12)


def _trace_with_lines(seed=0, n=200_000, n_pulses=150, lines=True):
    rng = np.random.default_rng(seed)
    x = rng.normal(0, 0.0005, n)
    # random (Poisson-like) arrival times, as in real data, at least 500 samples apart
    centres = np.sort(rng.choice(np.arange(3000, n - 3000, 500), n_pulses, replace=False)
                      + rng.integers(0, 300, n_pulses))
    amps = rng.uniform(0.05, 0.5, n_pulses)
    p = uncoupled_pulse(200)
    for c, a in zip(centres, amps):
        x[c:c + 200] += a * p / p.max()
    if lines:
        k = np.arange(n)
        x += 0.004 * np.sin(2 * np.pi * 257.4e3 * k / FS + 0.3)
        x += 0.002 * np.sin(2 * np.pi * 378.0e3 * k / FS + 1.1) * (1 + 0.2 * np.sin(2 * np.pi * 30 * k / FS))
    return x, centres, amps


def test_line_removal_reaches_white_noise_and_keeps_pulses():
    clean, centres, amps = _trace_with_lines(lines=False)
    noisy, _, _ = _trace_with_lines(lines=True)
    sig = noise_sigma(noisy)
    free = pulse_free_mask(noisy.size, find_pulses(noisy, sig, 5), 20, 60)
    lr = LineRemoval(FS)
    out = lr.process(noisy, free)
    assert noise_sigma(out) < 1.1 * noise_sigma(clean)
    assert any(abs(f - 257.4e3) < 500 for f in lr.info.freqs_hz)
    assert any(abs(f - 378.0e3) < 500 for f in lr.info.freqs_hz)
    # pulse heights unchanged: compare samples at the pulse maxima
    pk = centres + np.argmax(uncoupled_pulse(200))
    np.testing.assert_allclose(out[pk], clean[pk], atol=0.002)


def test_tracker_removes_slow_baseline():
    rng = np.random.default_rng(1)
    n = 100_000
    k = np.arange(n)
    drift = 0.01 * np.sin(2 * np.pi * k / 40_000)
    x = drift + rng.normal(0, 0.0005, n)
    peaks = np.arange(1000, n - 1000, 997)
    p = uncoupled_pulse(60)
    for c in peaks:
        x[c - 2:c + 58] += 0.2 * p / p.max()
    free = pulse_free_mask(n, peaks, 10, 60)
    base = track_baseline(x, free, TrackerConfig(window=1001))
    assert np.sqrt(np.mean((base - drift)[2000:-2000] ** 2)) < 2e-4


def test_preprocessor_flags():
    rng = np.random.default_rng(2)
    v = -(rng.normal(0, 0.001, 50_000) + 0.01)
    off = PreprocessConfig(line_removal=False, pole_zero=False, baseline_tracker=False)
    y = Preprocessor(FS, off).process(v)
    np.testing.assert_allclose(y, -v - np.mean(-v))
    with pytest.raises(ValueError):
        Preprocessor(FS, PreprocessConfig())                 # pole_zero needs a design
    d = design_pole_zero([0.4, 7.4, 80.0], [-1.0, 1.05, -0.05])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        Preprocessor(FS, PreprocessConfig(baseline_tracker=False), d)
        assert any("baseline_tracker" in str(m.message) for m in w)


def test_fit_pulses_fixed_baseline():
    t = np.arange(-40 * 16, 300 * 16 + 1) / 16
    s = np.where(t >= -1, np.exp(-np.clip(t + 1, 0, None) / 7.4) - np.exp(-np.clip(t + 1, 0, None) / 0.4), 0)
    T = PulseTemplate(16, float(t[0]), s / s.max())
    j = np.arange(-12, 51)
    A = np.array([0.1, 0.3, 0.7])
    w = A[:, None] * T(j[None, :] - 0.2)
    f = fit_pulses(w, T, 12, fit_range=(-10, 40), fit_baseline=False)
    np.testing.assert_allclose(f.amplitude, A, rtol=1e-3)
    np.testing.assert_allclose(f.baseline, 0)


def test_fit_peaks_doublet():
    rng = np.random.default_rng(3)
    v = np.concatenate([rng.normal(0.311, 0.0087, 20_000), rng.normal(0.352, 0.0089, 18_000),
                        rng.uniform(0.25, 0.42, 20_000)])
    pf = fit_peaks(v, [0.31, 0.35], 0.26, 0.40, 0.002)
    np.testing.assert_allclose(pf.centers, [0.311, 0.352], atol=5e-4)
    np.testing.assert_allclose(pf.sigmas, [0.0087, 0.0089], rtol=0.05)
