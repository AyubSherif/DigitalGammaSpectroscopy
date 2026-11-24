import numpy as np
import pytest

from gamma_spectrometry.baseline import design_from_template
from gamma_spectrometry.simulate import (SimConfig, Simulator, analytic_template, compare_to_truth,
                                klein_nishina_electron_cdf, load_truth, simulate_store)
from gamma_spectrometry.store import CaptureStore


def test_analytic_template_peaks_at_zero_and_has_a_model():
    T = analytic_template()
    assert T(np.array([0.0]))[0] == pytest.approx(1.0, abs=1e-3)
    assert T.shape.max() == pytest.approx(1.0, abs=1e-3)
    assert T.shape.min() < -0.02                                # the AC-coupling undershoot
    design_from_template(T)                                     # the pipeline can design its P/Z filter


@pytest.mark.parametrize("noise", [0.0, 0.0006])
def test_waveform_does_not_depend_on_read_sizes(noise):
    a, b = Simulator(SimConfig(noise_v=noise)), Simulator(SimConfig(noise_v=noise))
    xa = np.concatenate([a.read(n) for n in [1000, 77, 50000, 3, 98920]])
    xb = np.concatenate([b.read(75000), b.read(75000)])
    np.testing.assert_array_equal(xa, xb)
    np.testing.assert_array_equal(a.truth_table(), b.truth_table())


def test_rates_and_cobalt_lines():
    cfg = SimConfig(noise_v=0.0, lines_v=(0, 0, 0, 0))
    s = Simulator(cfg)
    n = int(1.8e6 * 3)
    s.read(n)
    T = s.truth_table()
    live = n / cfg.fs
    expected = cfg.rate_cps + cfg.background_cps
    assert abs(T.size - expected * live) < 5 * np.sqrt(expected * live)
    src = T[T["kind"] == 0]["energy_kev"]
    sigma = lambda e: 0.08 / 2.3548 * np.sqrt(661.657 * e)
    for e, side_dir in ((1173.2, -1), (1332.5, 1), (2505.7, 1)):    # compare with the continuum beside it,
        peak = np.sum(np.abs(src - e) < 2 * sigma(e))                  # away from the other line
        side = np.sum(np.abs(src - e - side_dir * 6 * sigma(e)) < 2 * sigma(e))
        assert peak > 2 * side + 10, e
    np.testing.assert_allclose(T["amplitude_v"], T["energy_kev"] / cfg.gain_kev_per_v, rtol=1e-5)


def test_signal_is_negative_quantised_with_offset():
    cfg = SimConfig()
    x = Simulator(cfg).read(200_000)
    assert x.min() < -0.1                                       # negative pulses
    np.testing.assert_allclose(x / cfg.lsb_v, np.rint(x / cfg.lsb_v), atol=1e-6)
    # without pulses the baseline is the offset (with pulses it is lifted by their undershoots)
    quiet = Simulator(SimConfig(rate_cps=0, background_cps=0)).read(200_000)
    assert np.median(quiet) == pytest.approx(cfg.offset_v, abs=3e-4)


def test_klein_nishina_ends_at_the_compton_edge():
    T, cdf = klein_nishina_electron_cdf(661.657)
    assert T[-1] == pytest.approx(477.3, abs=0.5)               # Cs-137 Compton edge
    assert np.all(np.diff(cdf) >= 0) and cdf[0] == 0 and cdf[-1] == pytest.approx(1.0)


def test_store_and_truth_and_scoring(tmp_path):
    from gamma_spectrometry.pipeline import PipelineConfig, run

    cfg = SimConfig(rate_cps=1500, background_cps=200)
    store, truth = simulate_store(tmp_path / "s.h5", cfg, n_captures=3, n_samples=120_000, gap_s=0.01)
    with CaptureStore(store) as s:
        assert len(s) == 3 and s.n_samples == 120_000
        assert np.diff(s.time_s) == pytest.approx(120_000 / 1.8e6 + 0.01)
    T = load_truth(truth)
    assert T.size > 0
    analytic_template().save(tmp_path / "t.h5")
    run(store, tmp_path / "t.h5", tmp_path / "ev.h5", PipelineConfig(), workers=1, verbose=False)
    c = compare_to_truth(tmp_path / "ev.h5", truth)
    for e in c["efficiency_by_energy"][1:]:
        if e["true"] > 20:
            assert e["efficiency"] > 0.97, e
    assert abs(c["amplitude_bias_pct"]) < 1.5
    assert c["false_rate_cps"] < 0.03 * c["true_rate_cps"]
