import numpy as np
import pytest

from gamma_spectrometry.detect import (DetectConfig, Detector, apply_filter, local_sigma, negative_side_sigma,
                              noise_autocovariance, optimal_kernel, sharpen, successive_cancellation)
from gamma_spectrometry.events import (FLAG_PILEUP, FLAG_WARMUP, EventFile, EventWriter, EVENT_DTYPE,
                              CAPTURE_DTYPE)
from gamma_spectrometry.template import PulseTemplate


def pz_like_template():
    """A pulse like the pole-zero-corrected one: ~1 sample rise, 7.7-sample decay."""
    t = np.arange(-40 * 16, 400 * 16 + 1) / 16
    s = np.where(t >= -1.2, np.exp(-np.clip(t + 1.2, 0, None) / 7.7) - np.exp(-np.clip(t + 1.2, 0, None) / 0.4), 0)
    return PulseTemplate(16, float(t[0]), s / s.max())


T = pz_like_template()


def place(n, positions, amps, noise=0.0003, seed=0, colored=False):
    rng = np.random.default_rng(seed)
    x = rng.normal(0, noise, n)
    if colored:  # add low-frequency noise, as in the real clean signal
        x += np.convolve(rng.normal(0, noise, n), np.ones(40) / 6, mode="same")
    k = np.arange(-10, 300)
    for c, a in zip(positions, amps):
        idx = int(np.floor(c)) + k
        idx = idx[(idx >= 0) & (idx < n)]
        x[idx] += a * T(idx - c)
    return x


def test_apply_filter_gives_amplitude():
    s = T(np.arange(-6, 61.0))
    x = place(2000, [1000.0], [0.25], noise=0)
    y = apply_filter(x, s, s, 6)
    assert np.argmax(y) == 1000 and y[1000] == pytest.approx(0.25, rel=1e-9)
    h = s * np.linspace(1, 0.5, s.size)          # any kernel is normalised by h.s
    assert apply_filter(x, h, s, 6)[1000] == pytest.approx(0.25, rel=1e-9)


def test_noise_autocovariance_white_and_gapped():
    rng = np.random.default_rng(1)
    x = rng.normal(0, 2.0, 200_000)
    free = rng.random(x.size) > 0.3
    r = noise_autocovariance(x, free, 5)
    assert r[0] == pytest.approx(4.0, rel=0.03)
    assert np.all(np.abs(r[1:]) < 0.06)


def test_optimal_kernel_beats_matched_in_colored_noise():
    s = T(np.arange(-6, 61.0))
    noise = place(300_000, [], [], colored=True)
    r = noise_autocovariance(noise, np.ones(noise.size, bool), s.size)
    h = optimal_kernel(s, r)
    sig_opt = negative_side_sigma(apply_filter(noise, h, s, 6))
    sig_mf = negative_side_sigma(apply_filter(noise, s, s, 6))
    assert sig_opt < 0.8 * sig_mf


def test_sharpen_removes_exponential():
    k = np.arange(200)
    x = np.exp(-k / 7.7)
    assert np.max(np.abs(sharpen(x, 7.7)[1:])) < 1e-12


def test_local_sigma_tracks_bursts():
    rng = np.random.default_rng(2)
    y = rng.normal(0, 1, 40_000)
    y[20_000:] *= 4
    ls = local_sigma(y, 4096)
    assert ls[5000] == pytest.approx(1, rel=0.1) and ls[35_000] == pytest.approx(4, rel=0.1)


def test_detector_single_pulses_no_extras():
    rng = np.random.default_rng(3)
    pos = np.arange(2000, 98_000, 700) + rng.uniform(0, 1, 138)
    amps = rng.choice([0.01, 0.05, 0.3, 1.0], pos.size)
    x = place(100_000, pos, amps)
    d = Detector(T).detect(x)
    assert d.index.size == pos.size                       # no artifacts, nothing missed
    assert np.all(np.abs(d.index - np.round(pos)) <= 1)
    np.testing.assert_allclose(d.filter_amplitude, amps, rtol=0.05)


def test_detector_resolves_close_equal_pair():
    pos = []
    for base in range(2000, 90_000, 1000):
        pos += [base + 0.3, base + 4.3]
    x = place(100_000, pos, [0.3] * len(pos))
    d = Detector(T).detect(x)
    assert d.index.size == len(pos)


def test_successive_cancellation_removes_ripple():
    # response with negative dips and a positive bump before the peak
    span = 20
    tau = np.arange(-span, span + 1)
    resp = np.exp(-0.5 * (tau / 1.2) ** 2) - 0.2 * np.exp(-0.5 * ((tau + 6) / 2) ** 2) \
        + 0.05 * np.exp(-0.5 * ((tau + 12) / 1) ** 2)
    R = np.tile(resp, (16, 1))
    y = np.zeros(400)
    y[200 - span: 200 + span + 1] += 1.0 * resp
    idx, amp, _ = successive_cancellation(y, R, 0.01, 0.0)
    assert list(idx) == [200] and amp[0] == pytest.approx(1.0)


def test_event_file_roundtrip(tmp_path):
    ev = np.zeros(3, EVENT_DTYPE)
    ev["index"] = [10, 20, 30]
    ev["flags"] = [0, FLAG_PILEUP, FLAG_WARMUP]
    ev["amplitude"] = [0.1, np.nan, np.nan]
    ev["cluster"] = [0, 1, 1]
    cap = np.zeros(1, CAPTURE_DTYPE)
    cap["live_time_s"] = 0.25
    path = tmp_path / "e.h5"
    with EventWriter(path, {"config": {"a": 1}}) as w:
        w.append(ev, cap, [(1, 0, 15, np.arange(5.0))])
        w.append(ev, cap, None)
    ef = EventFile(path)
    assert ef.events.size == 6 and ef.live_time_s == pytest.approx(0.5)
    assert ef.clean().size == 2 and ef.attrs["config"] == {"a": 1}
    cap_i, start, w = ef.snippet(1)
    assert start == 15 and list(w) == [0, 1, 2, 3, 4]
    with pytest.raises(FileExistsError):
        EventWriter(path, {})


def test_detect_config_validation():
    with pytest.raises(ValueError):
        Detector(T, DetectConfig(method="nope"))
