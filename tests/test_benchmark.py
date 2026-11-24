import numpy as np
import pytest

from gamma_spectrometry.benchmark import shaping_amplifier

DT = 1 / 1.8e6


def exp_pulse(n=400, tau=4.1e-6, start=50):
    k = np.arange(n)
    return np.where(k >= start, np.exp(-(k - start) * DT / tau), 0.0)


@pytest.mark.parametrize("tau_s", [0.5e-6, 1e-6, 2e-6])
def test_shaper_peaks_near_n_tau_and_returns_to_baseline(tau_s):
    y = shaping_amplifier(exp_pulse(), DT, tau_shaping=tau_s, tau_pz=4.1e-6, n_integrators=4)
    # analog CR-RC^n peaks at n*tau; each zero-delay digital stage leads by ~dt/2,
    # so the sampled peak sits ~(n+1)/2 samples earlier
    expected = 4 * tau_s / DT - (4 + 1) / 2
    assert abs((np.argmax(y) - 50) - expected) <= 1.0
    assert abs(y[-1]) < 1e-3 * y.max()                         # pole-zero: no long tail / undershoot
    assert y.min() > -1e-3 * y.max()


def test_shaper_is_linear():
    a = shaping_amplifier(exp_pulse(), DT)
    b = shaping_amplifier(3.0 * exp_pulse(), DT)
    np.testing.assert_allclose(b, 3 * a, rtol=1e-12, atol=1e-15)


def test_wrong_pole_zero_leaves_a_tail():
    y_ok = shaping_amplifier(exp_pulse(), DT, tau_pz=4.1e-6)
    y_bad = shaping_amplifier(exp_pulse(), DT, tau_pz=2.0e-6)  # under-compensated
    assert abs(y_bad[300]) > 10 * abs(y_ok[300])
