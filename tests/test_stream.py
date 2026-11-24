import numpy as np
import pytest

from gammaspec import config as C
from gammaspec.events import EventFile
from gammaspec.pipeline import PipelineConfig, run
from gammaspec.simulate import SimConfig, Simulator, analytic_template, simulate_store
from gammaspec.sources import PipeSource, SimulationSource, StoreSource
from gammaspec.stream import StreamConfig, StreamProcessor, check_stream_config, process_stream

SMALL = dict(frame=32768, pre=2048, post=1536, update_s=0)


@pytest.fixture(scope="module")
def template():
    return analytic_template()


def _events(sp):
    sp.close()
    return np.concatenate([r[0] for r in sp.results]) if sp.results else np.zeros(0)


def test_frame_geometry_is_checked():
    cfg = PipelineConfig()
    assert check_stream_config(StreamConfig(), cfg) == []
    errs = check_stream_config(StreamConfig(pre=100, post=10, frame=100), cfg)
    assert len(errs) == 3


def test_stream_with_capture_sized_frames_equals_batch(tmp_path):
    store, _ = simulate_store(tmp_path / "s.h5", SimConfig(rate_cps=2000), n_captures=3, n_samples=60_000)
    analytic_template().save(tmp_path / "t.h5")
    run(store, tmp_path / "t.h5", tmp_path / "b.h5", PipelineConfig(), workers=1, verbose=False)
    process_stream(StoreSource(store, chunk=7000), tmp_path / "t.h5", tmp_path / "s_ev.h5", PipelineConfig(),
                   StreamConfig(frame=65536, update_s=0), workers=1, verbose=False)
    B, S = EventFile(tmp_path / "b.h5"), EventFile(tmp_path / "s_ev.h5")
    for name in ("capture", "index", "flags", "cluster_size"):
        np.testing.assert_array_equal(B.events[name], S.events[name])
    np.testing.assert_allclose(B.events["amplitude"], S.events["amplitude"], equal_nan=True)
    assert B.live_time_s == pytest.approx(S.live_time_s)
    np.testing.assert_allclose(B.captures["live_start_s"], S.captures["live_start_s"])


def test_results_do_not_depend_on_chunk_sizes(template):
    x = Simulator(SimConfig()).read(150_000)
    out = []
    for chunk in (1000, 150_000):
        sp = StreamProcessor(template, 1.8e6, PipelineConfig(), StreamConfig(**SMALL), workers=1)
        for a in range(0, x.size, chunk):
            sp.feed(x[a:a + chunk])
        out.append(_events(sp))
    np.testing.assert_array_equal(out[0]["time_s"], out[1]["time_s"])
    np.testing.assert_array_equal(out[0]["amplitude"], out[1]["amplitude"])


def test_overlapping_frames_match_one_long_frame(template):
    x = Simulator(SimConfig()).read(200_000)
    runs = {}
    for frame in (32768, 200_000):
        sp = StreamProcessor(template, 1.8e6, PipelineConfig(), StreamConfig(**{**SMALL, "frame": frame}), workers=1)
        sp.feed(x)
        e = _events(sp)
        runs[frame] = e[np.isfinite(e["amplitude"]) & (e["flags"] & 0b111110 == 0)]   # good events
        assert np.all(np.diff(e["time_s"]) >= 0)                # written in time order
        assert sp.status.live_time_s == pytest.approx((200_000 - 1500 - 41) / 1.8e6)
    long, short = runs[200_000], runs[32768]
    j = np.clip(np.searchsorted(short["time_s"], long["time_s"]), 1, short.size - 1)
    d = np.minimum(np.abs(short["time_s"][j] - long["time_s"]), np.abs(short["time_s"][j - 1] - long["time_s"]))
    assert np.mean(d * 1.8e6 < 0.5) > 0.98
    assert abs(short.size - long.size) < 0.01 * long.size
    # no event counted twice at a frame boundary
    assert np.min(np.diff(np.sort(short["time_s"]))) * 1.8e6 > 0.5


def test_live_source_drops_frames_and_counts_them_as_dead_time(template):
    x = Simulator(SimConfig()).read(32768 * 12 + 4000)
    sp = StreamProcessor(template, 1.8e6, PipelineConfig(), StreamConfig(**{**SMALL, "max_backlog": 1}),
                         workers=2, live=True)
    for a in range(0, x.size, 32768):
        sp.feed(x[a:a + 32768])
    st = sp.close()
    assert st.frames + st.dropped_frames == sp._frame_id
    # each processed frame contributes its core (the first frame minus the warm-up, the last minus the end margin)
    assert st.live_time_s <= st.frames * 32768 / 1.8e6 + 1e-9


def test_pipe_source_reads_binary_samples(tmp_path):
    v = np.linspace(-0.5, 0.5, 10_001)
    (tmp_path / "f4.bin").write_bytes(v.astype("<f4").tobytes())
    got = np.concatenate([c.data for c in PipeSource(str(tmp_path / "f4.bin"), chunk=999).chunks()])
    np.testing.assert_allclose(got, v, atol=1e-7)
    codes = np.arange(-300, 300, dtype="<i2")
    (tmp_path / "i2.bin").write_bytes(codes.tobytes())
    got = np.concatenate([c.data for c in PipeSource(str(tmp_path / "i2.bin"), dtype="i2", scale=1e-3, chunk=64).chunks()])
    np.testing.assert_allclose(got, codes * 1e-3)


def test_simulation_source_is_one_continuous_segment():
    chunks = list(SimulationSource(SimConfig(), duration_s=0.05, chunk=20_000).chunks())
    assert chunks[0].new_segment and not any(c.new_segment for c in chunks[1:])
    assert sum(c.data.size for c in chunks) == 90_000


def test_config_for_streamed_kinds(tmp_path):
    cfg = C.load(None, ["input.kind='simulation'"])
    assert C.check(cfg) == []                                   # needs no paths; simulates its own background
    C.check(C.load(None, ["input.kind='pipe'", "input.source='-'", f"template.path='{(tmp_path / 't.h5').as_posix()}'"]),
            paths=False)
    with pytest.raises(C.ConfigError, match="stream.pre"):
        C.check(C.load(None, ["input.kind='simulation'", "stream.pre=10"]))
    with pytest.raises(C.ConfigError, match="lines_hz has 1 values"):
        C.check(C.load(None, ["input.kind='simulation'", "simulation.signal.lines_hz=[9000.0]"]))
    with pytest.raises(C.ConfigError, match="not one of"):
        C.load(None, ["simulation.signal.nuclide='u238'"])


def test_runner_simulation_with_truth_and_reuse(tmp_path, capsys):
    from gammaspec import runner

    sets = [f"output='{(tmp_path / 'out').as_posix()}'", "input.kind='simulation'", "simulation.duration_s=0.4",
            "simulation.background_duration_s=0", "spectrum.enabled=false", "workers=1", "stream.update_s=0",
            "stream.frame=65536"]
    s = runner.run(C.load(None, sets), verbose=True)
    check = (tmp_path / "out" / "simulation_check.json")
    assert check.exists() and (tmp_path / "out" / "truth" / "source.h5").exists()
    assert EventFile(s["events"]["source"]).events.size > 1000
    capsys.readouterr()
    runner.run(C.load(None, sets), verbose=True)
    assert "[process] source: reusing" in capsys.readouterr().out
