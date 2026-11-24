import time
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402

APP = str(Path(__file__).resolve().parents[1] / "gamma_spectrometry" / "app.py")


def _app():
    at = AppTest.from_file(APP, default_timeout=60)
    at.run()
    assert not at.exception, at.exception
    return at


def test_page_renders_with_defaults():
    at = _app()
    assert at.title[0].value == "Gamma Spectrometry App"
    # default input kind "store" with no source: the run is invalid, Run is disabled
    assert any("input.source: not set" in e.value for e in at.error)
    run = [b for b in at.button if b.label == "Run"][0]
    assert run.disabled


def test_widgets_edit_the_run_file():
    at = _app()
    at.radio(key="w.input.kind").set_value("simulation").run()
    assert not at.exception
    cfg = at.session_state["cfg"]
    assert cfg.input.kind == "simulation"
    at.number_input(key="w.sim.dur").set_value(0.5).run()
    at.checkbox(key="w.pp.line").set_value(False).run()
    cfg = at.session_state["cfg"]
    assert cfg.simulation.duration_s == 0.5 and cfg.pipeline.preprocess.line_removal is False
    assert not at.error                                  # a simulation needs no paths: valid
    assert not [b for b in at.button if b.label == "Run"][0].disabled


def test_toml_editor_applies_and_reports_errors():
    at = _app()
    text = at.text_area[0].value.replace('kind = "store"', 'kind = "simulation"')
    at.text_area[0].set_value(text)
    [b for b in at.button if b.label == "Apply"][0].click().run()
    assert at.session_state["cfg"].input.kind == "simulation"
    assert at.radio(key="w.input.kind").value == "simulation"   # widgets follow the text
    bad = at.text_area[0].value.replace('kind = "simulation"', 'kind = "nope"')
    at.text_area[0].set_value(bad)
    [b for b in at.button if b.label == "Apply"][0].click().run()
    assert any("'nope' is not one of" in e.value for e in at.error)


def test_run_button_starts_a_run_and_shows_results(tmp_path):
    at = _app()
    at.radio(key="w.input.kind").set_value("simulation").run()
    at.text_input(key="w.output").set_value(str(tmp_path / "out")).run()
    at.number_input(key="w.sim.dur").set_value(0.4).run()
    at.number_input(key="w.sim.bdur").set_value(0.0).run()
    at.checkbox(key="w.spec.on").set_value(False).run()
    at.number_input(key="w.workers").set_value(1).run()
    [b for b in at.button if b.label == "Run"][0].click().run()
    proc = at.session_state["proc"]
    assert proc.wait(timeout=300) == 0, (tmp_path / "out" / "ui_run.log").read_text()
    assert (tmp_path / "out" / "ui_run.toml").exists()
    assert (tmp_path / "out" / "simulation_check.json").exists()
    time.sleep(0.1)
    at.run()
    assert not at.exception
    assert any("against the truth" in m.value for m in at.markdown)
