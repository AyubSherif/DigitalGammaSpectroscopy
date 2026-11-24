import numpy as np
import pytest

from gammaspec.detect import DetectConfig
from gammaspec.events import FLAG_FIT_ADDED, FLAG_PILEUP
from gammaspec.pipeline import Pipeline, PipelineConfig
from gammaspec.preprocess import PreprocessConfig
from gammaspec.solver import PileupSolver, SolverConfig, taper_autocovariance

from .test_detect_events import T, place

NOISE = 5e-4


def white_acov(sigma=NOISE, n=512):
    return np.r_[sigma ** 2, np.zeros(n - 1)]


def wave(n, cs, As, seed=0):
    t = np.arange(n, dtype=float)
    rng = np.random.default_rng(seed)
    return sum(A * T(t - c) for c, A in zip(cs, As)) + rng.normal(0, NOISE, n)


def test_taper_keeps_head_and_goes_to_zero():
    r = np.ones(100)
    w = taper_autocovariance(r, 64)
    assert w.size == 64 and w[0] == 1 and np.all(w[:48] == 1) and w[-1] < 0.01


def test_solver_fits_pair_with_timing():
    sv = PileupSolver(T, white_acov(), SolverConfig(residual_search=False))
    res = sv.solve(wave(80, [10.4, 17.7], [0.3, 0.05]), [10, 18])
    np.testing.assert_allclose(res.amplitude, [0.3, 0.05], rtol=0.03, atol=0.002)
    np.testing.assert_allclose(res.time, [10.4, 17.7], atol=0.1)
    assert res.whitened and res.chi2 < 1.5 and not res.added.any()


def test_residual_search_finds_merged_pulse():
    sv = PileupSolver(T, white_acov(), SolverConfig())
    res = sv.solve(wave(80, [10.2, 12.4], [0.3, 0.1]), [10])          # detector merged them
    assert res.added.sum() == 1
    np.testing.assert_allclose(np.sort(res.time), [10.2, 12.4], atol=0.25)
    np.testing.assert_allclose(res.amplitude[np.argsort(res.time)], [0.3, 0.1], rtol=0.08)


def test_residual_search_adds_nothing_to_a_single_pulse():
    sv = PileupSolver(T, white_acov(), SolverConfig())
    added = [sv.solve(wave(60, [10.0 + 0.1 * k], [0.2], seed=k), [round(10.0 + 0.1 * k)]).added.any() for k in range(30)]
    assert not any(added)


def test_solver_without_noise_model_falls_back_to_ls():
    sv = PileupSolver(T, None, SolverConfig())
    res = sv.solve(wave(80, [10.0, 20.0], [0.2, 0.1]), [10, 20])
    assert not res.whitened
    np.testing.assert_allclose(res.amplitude, [0.2, 0.1], rtol=0.03)


def test_pipeline_capture_with_solver():
    """Pile-up pulses get solver amplitudes; a merged pair is split by the residual search."""
    rng = np.random.default_rng(3)
    n = 60_000
    pos = list(np.arange(3000, n - 500, 700) + rng.uniform(0, 1))
    truth = {round(p): 0.2 for p in pos}
    extra = [pos[3] + 20.3, pos[6] + 2.4]                 # a normal pile-up and a 2-sample one
    amps = [0.2] * len(pos) + [0.08, 0.1]
    x = place(n, pos + extra, amps, noise=NOISE)
    cfg = PipelineConfig(preprocess=PreprocessConfig(line_removal=False, pole_zero=False,
                                                     baseline_tracker=True, warmup_samples=0),
                         detect=DetectConfig())
    with pytest.warns(UserWarning):
        pipe = Pipeline(T, 1.8e6, cfg)
    pipe.solver_on = True                                # synthetic pulses have no undershoot
    ev, row, sn = pipe.process_capture(-x, 0, 0.0)
    good = np.isfinite(ev["amplitude"])
    assert good.all()
    pu = ev[(ev["flags"] & FLAG_PILEUP) != 0]
    assert pu.size == 4
    np.testing.assert_allclose(np.sort(pu["amplitude"]), [0.08, 0.1, 0.2, 0.2], rtol=0.05)
    assert np.sum((ev["flags"] & FLAG_FIT_ADDED) != 0) == 1
    assert row[0]["n_events"] == ev.size == len(pos) + 2
