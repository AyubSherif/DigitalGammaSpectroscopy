import datetime as dt

import numpy as np
import pytest

from gamma_spectrometry.detect import DetectConfig
from gamma_spectrometry.events import FLAG_PILEUP, FLAG_WARMUP, EventFile
from gamma_spectrometry.pipeline import EventConfig, PipelineConfig, run, summarize
from gamma_spectrometry.preprocess import PreprocessConfig
from gamma_spectrometry.rawio import DEFAULT_LSB_VOLTS, build_index
from gamma_spectrometry.store import convert

from .conftest import write_capture
from .test_detect_events import T, place


def test_end_to_end_synthetic(tmp_path):
    """ASCII captures -> store -> pipeline -> events: counts, amplitudes, pile-up flags."""
    rng = np.random.default_rng(7)
    n = 60_000
    truth = []
    t0 = dt.datetime(2021, 8, 26, 13, 0, 0)
    for cap in range(3):
        pos = np.arange(2000, n - 500, 900).astype(float) + rng.uniform(0, 1)
        pos = np.r_[pos, pos[::5] + 25.0]            # every 5th pulse gets a neighbour 25 samples later
        pos.sort()
        amps = rng.uniform(0.05, 0.4, pos.size)
        x = place(n, pos, amps, noise=0.0005, seed=cap)
        codes = np.rint(-x / DEFAULT_LSB_VOLTS).astype(np.int64)   # negative pulses, as recorded
        write_capture(tmp_path / f"raw_data_{cap}.txt", codes, t0 + dt.timedelta(seconds=2.3 * cap))
        truth.append((pos, amps))
    store = convert(build_index(tmp_path), tmp_path / "s.h5", workers=1)
    tpath = tmp_path / "t.h5"
    T.save(tpath)
    cfg = PipelineConfig(preprocess=PreprocessConfig(line_removal=False, pole_zero=False),
                         detect=DetectConfig(), events=EventConfig())
    r = run(store, tpath, tmp_path / "ev.h5", cfg, workers=1, verbose=False)
    ef = EventFile(tmp_path / "ev.h5")
    assert r["events"] == sum(p.size for p, _ in truth)
    ev = ef.events
    for cap, (pos, amps) in enumerate(truth):
        e = ev[ev["capture"] == cap]
        # peak sample: either neighbour of the true (sub-sample) peak time is correct
        assert np.all(np.abs(e["index"] - pos) <= 1.0)
        iso = (e["flags"] & (FLAG_PILEUP | FLAG_WARMUP)) == 0
        np.testing.assert_allclose(e["amplitude"][iso], amps[iso], rtol=0.03, atol=0.003)
        paired = np.r_[np.diff(pos) < 60, False] | np.r_[False, np.diff(pos) < 60]
        np.testing.assert_array_equal((e["flags"] & FLAG_PILEUP) != 0, paired)
    s = summarize(tmp_path / "ev.h5")
    assert s["snippets"] > 0 and s["events"] == ev.size
    assert ef.captures.size == 3 and np.all(ef.captures["n_events"] > 0)
