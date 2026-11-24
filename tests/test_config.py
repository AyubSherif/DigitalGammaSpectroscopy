import datetime as dt
import json
import tomllib

import numpy as np
import pytest

from gamma_spectrometry import config as C
from gamma_spectrometry.pipeline import PipelineConfig
from gamma_spectrometry.rawio import DEFAULT_LSB_VOLTS

from .conftest import write_capture
from .test_detect_events import T, place


def test_defaults_round_trip_through_toml_and_json(tmp_path):
    cfg = C.RunConfig()
    assert C.from_dict(tomllib.loads(C.dumps(cfg))) == cfg
    assert C.load(C.save(cfg, tmp_path / "r.json")) == cfg
    assert C.load(C.save(cfg, tmp_path / "r.toml")) == cfg
    p = PipelineConfig()
    assert C.from_dict(tomllib.loads(C.dumps(p))["pipeline"], PipelineConfig) == p


def test_default_file_documents_the_settings():
    text = C.dumps(C.RunConfig())
    assert "[pipeline.preprocess]" in text and "[pipeline.solver]" in text
    assert "pole_zero = true  # stage:" in text
    assert "# limit = (not set)" in text                   # optional values are commented out


def test_partial_file_keeps_defaults(tmp_path):
    f = tmp_path / "r.toml"
    f.write_text('output = "o"\n[pipeline.detect]\nthreshold = 6\n', encoding="utf-8")
    cfg = C.load(f)
    assert cfg.output == "o" and cfg.pipeline.detect.threshold == 6.0
    assert isinstance(cfg.pipeline.detect.threshold, float)     # int accepted for a float
    assert cfg.pipeline.detect.method == "optimal" and cfg.pipeline.solver.enabled


def test_all_errors_are_reported_together():
    with pytest.raises(C.ConfigError) as e:
        C.from_dict({"input": {"kind": "raws", "limt": 3},
                     "pipeline": {"detect": {"threshold": "high"}, "solver": {"search_threshold": -1},
                                  "preprocess": {"pole_zero": 1}}})
    msg = str(e.value)
    for part in ("input.limt: unknown setting (did you mean 'limit'?)", "input.kind: 'raws' is not one of",
                 "pipeline.detect.threshold: expected a number", "search_threshold: must be > 0",
                 "pipeline.preprocess.pole_zero: expected true or false"):
        assert part in msg


def test_unsupported_schema_version():
    with pytest.raises(C.ConfigError, match="schema_version"):
        C.from_dict({"schema_version": 99})


@pytest.mark.parametrize("text,value", [("6", 6), ("1e-3", 1e-3), ("true", True), ('"edge"', "edge"),
                                        ("edge", "edge"), ("[1, 2]", [1, 2])])
def test_set_key_parses_toml_values(text, value):
    d = {}
    C.set_key(d, f"pipeline.detect.x={text}")
    assert d == {"pipeline": {"detect": {"x": value}}}


def test_overrides_are_validated(tmp_path):
    cfg = C.load(None, ["pipeline.detect.method=matched", "input.limit=5"])
    assert cfg.pipeline.detect.method == "matched" and cfg.input.limit == 5
    with pytest.raises(C.ConfigError, match="did you mean 'threshold'"):
        C.load(None, ["pipeline.detect.treshold=5"])
    with pytest.raises(C.ConfigError, match="KEY=VALUE"):
        C.load(None, ["pipeline.detect.threshold"])


def test_check_combinations(tmp_path):
    store = tmp_path / "s.h5"
    store.write_bytes(b"")
    cfg = C.load(None, [f"input.source='{store.as_posix()}'", "template.build_from='source'",
                        "pipeline.preprocess.baseline_tracker=false"])
    w = C.check(cfg)
    assert any("pole_zero without baseline_tracker" in x for x in w)
    assert any("solver needs" in x for x in w)
    assert any("not background-subtracted" in x for x in w)
    with pytest.raises(C.ConfigError, match="input.source: not set"):
        C.check(C.RunConfig())
    with pytest.raises(C.ConfigError, match="build_from = 'background' but there is no"):
        C.check(C.load(None, [f"input.source='{store.as_posix()}'"]))
    with pytest.raises(C.ConfigError, match="is a folder"):
        C.check(C.load(None, [f"input.source='{tmp_path.as_posix()}'", "template.build_from='source'"]))
    with pytest.raises(C.ConfigError, match="not a folder"):
        C.check(C.load(None, [f"input.source='{store.as_posix()}'", "input.kind='raw'",
                              "template.build_from='source'"]))


def _synthetic_raw(folder, n_cap=2, n=40_000):
    folder.mkdir()
    t0 = dt.datetime(2021, 8, 26, 13, 0, 0)
    for cap in range(n_cap):
        pos = np.arange(2000, n - 500, 700).astype(float) + 0.3
        x = place(n, pos, np.full(pos.size, 0.2), noise=0.0005, seed=cap)
        write_capture(folder / f"raw_data_{cap}.txt", np.rint(-x / DEFAULT_LSB_VOLTS).astype(np.int64),
                      t0 + dt.timedelta(seconds=2 * cap))


def test_runner_raw_to_events_with_reuse(tmp_path, capsys):
    from gamma_spectrometry import runner
    from gamma_spectrometry.events import EventFile

    _synthetic_raw(tmp_path / "raw")
    T.save(tmp_path / "t.h5")
    sets = [f"output='{(tmp_path / 'out').as_posix()}'", "input.kind='raw'",
            f"input.source='{(tmp_path / 'raw').as_posix()}'", f"template.path='{(tmp_path / 't.h5').as_posix()}'",
            "pipeline.preprocess.line_removal=false", "pipeline.preprocess.pole_zero=false",
            "spectrum.enabled=false", "workers=1"]
    s = runner.run(C.load(None, sets), verbose=True)
    ev = EventFile(s["events"]["source"])
    assert ev.captures.size == 2 and ev.events.size > 100
    assert C.load(tmp_path / "out" / "config.toml") == C.load(None, sets)        # effective config written
    assert json.loads((tmp_path / "out" / "run_summary.json").read_text())["steps"]["process_source"]
    assert C.pipeline_from_event_file(s["events"]["source"]) == C.load(None, sets).pipeline
    capsys.readouterr()

    runner.run(C.load(None, sets), verbose=True)                                    # same: all reused
    out = capsys.readouterr().out
    assert "[input] source: reusing" in out and "[process] source: reusing" in out

    runner.run(C.load(None, sets + ["pipeline.detect.threshold=7"]), verbose=True)  # only processing reruns
    out = capsys.readouterr().out
    assert "[input] source: reusing" in out and "settings changed, recomputing" in out
    assert C.pipeline_from_event_file(s["events"]["source"]).detect.threshold == 7.0
