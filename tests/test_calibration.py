import numpy as np
import pytest

from gamma_spectrometry.analysis import EnergySpectrum, NetSpectrum, energy_spectrum, fit_line_groups, subtract
from gamma_spectrometry.calibration import (Calibration, LinePoint, anchor_gain, calibrate, fit_calibration,
                                   gain_match)


def synthetic_run(gain_kev_per_v, offset_kev, lines_kev_counts, n_continuum, seed=0):
    """Amplitudes (V) of a run with Gaussian lines (NaI-like widths) on a falling continuum."""
    rng = np.random.default_rng(seed)
    e = [rng.exponential(700.0, n_continuum)]
    for ek, n in lines_kev_counts:
        sig = ek * 0.065 / 2.3548 * np.sqrt(1173 / ek)
        e.append(rng.normal(ek, sig, n))
    e = np.concatenate(e)
    e = e[(e > 100) & (e < 3000)]
    return (e - offset_kev) / gain_kev_per_v


CO = [(1173.228, 15000), (1332.492, 13500), (2505.692, 600)]


def test_anchor_finds_co60_pair():
    a = synthetic_run(3760.0, -12.0, CO, 200_000)
    assert anchor_gain(a, "source") == pytest.approx(3760 * (1 - 12 / 1250), rel=0.02)


def test_calibrate_source_recovers_linear_calibration():
    a = synthetic_run(3760.0, -12.0, CO, 200_000)
    cal = calibrate(a, "source", lines=[("Co-60", 1173.228), ("Co-60", 1332.492), ("Co-60 sum", 2505.692)])
    assert sum(p.accepted for p in cal.points) == 3
    for _, _, r in cal.residuals_kev():
        assert abs(r) < 3
    assert cal.energy(0.4) == pytest.approx(-12 + 3760 * 0.4, abs=4)


def test_gain_match_uses_reference_shape():
    ref = Calibration([-12.0, 3760.0])
    # a run whose amplitudes are 7 % smaller: E = ref(k A) with k = 1/0.93
    k = 1 / 0.93
    a_k40 = (1460.82 + 12.0) / 3760.0 / k
    pts = [LinePoint("K-40", 1460.82, a_k40, 1e-4, 6.0, 900, True)]
    cal = gain_match(pts, ref)
    assert cal.energy(a_k40) == pytest.approx(1460.82, abs=0.01)
    assert cal.coeffs[1] == pytest.approx(3760.0 * k, rel=1e-6)


def test_fit_calibration_needs_two_lines():
    with pytest.raises(ValueError):
        fit_calibration([LinePoint("K-40", 1460.8, 0.39, 1e-4, 6, 100, True)])


def test_subtraction_and_net_lines():
    edges = np.arange(0, 3000.1, 5.0)
    rng = np.random.default_rng(1)
    bkg_e = np.concatenate([rng.exponential(500, 40_000), rng.normal(1460.82, 36, 3000)])
    src_e = np.concatenate([rng.exponential(500, 40_000), rng.normal(1460.82, 36, 3000),
                            rng.normal(1173.2, 32, 20_000), rng.normal(1332.5, 34, 18_000)])
    S = energy_spectrum(src_e, 100.0, edges)
    B = energy_spectrum(bkg_e, 100.0, edges)
    N = subtract(S, B)
    assert np.allclose(N.rate, S.rate - B.rate)
    assert np.allclose(N.rate_err ** 2, np.maximum(S.counts, 1) / 1e4 + np.maximum(B.counts, 1) / 1e4)
    lines = {(r["label"], round(r["energy_kev"])): r for r in fit_line_groups(N)}
    assert lines[("Co-60", 1173)]["rate_cps"] == pytest.approx(200, rel=0.05)
    assert lines[("Co-60", 1332)]["rate_cps"] == pytest.approx(180, rel=0.05)
    k40 = lines[("K-40", 1461)]
    assert abs(k40["rate_cps"]) < 3 * k40["rate_err_cps"]          # background line vanishes
    with pytest.raises(ValueError):
        subtract(S, energy_spectrum(bkg_e, 100.0, np.arange(0, 3000.1, 10.0)))
