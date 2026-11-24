"""gammaspec app: a browser UI for the whole chain (step 11).

    gammaspec app                 # opens http://localhost:8501
    gammaspec app --port 8600

The UI edits a run file (``docs/configuration.md``):
* common choices as widgets: input, stages, simulation, output;
* every other setting in the "All settings" TOML editor.

**Run** saves the file to ``<output>/ui_run.toml`` and starts ``gammaspec run`` on it as a
separate process, so the page stays responsive and **Stop** can end it. While the run is
going, the page shows its log and the live spectrum (``live/<run>/live_spectrum.json``,
written by the stream processor). Afterwards it shows the results.

Nothing here processes data itself: the UI is a front end for ``config`` and ``runner``.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
import tomllib
from pathlib import Path

import numpy as np
import streamlit as st

from gammaspec import __version__
from gammaspec import config as C

KINDS = {
    "raw": "Raw ADC files (convert, then process)",
    "store": "Prepared HDF5 store (.h5)",
    "replay": "Recorded data as a stream (replay)",
    "simulation": "Simulated signal",
    "pipe": "Live stream from a file / FIFO",
}
STEPS_HELP = {
    "line_removal": "Removes narrow interference lines (sinusoids from power supplies etc.).",
    "pole_zero": "Removes the undershoot after each pulse, so a pulse on another's tail is measured correctly. Needs the baseline tracker.",
    "baseline_tracker": "Follows the slow baseline and subtracts it.",
    "solver": "Separates overlapping pulses (pile-up) and measures each one.",
    "residual_search": "Looks for pulses hidden inside pile-up clusters.",
}


# ----------------------------------------------------------------------------------------------
# State
# ----------------------------------------------------------------------------------------------

def _cfg() -> C.RunConfig:
    if "cfg" not in st.session_state:
        st.session_state.cfg = C.load(None, ["output='results/ui_run'"])
    return st.session_state.cfg


def _replace_cfg(cfg: C.RunConfig) -> None:
    """New configuration (loaded or edited as text): widgets must show its values."""
    st.session_state.cfg = cfg
    for k in [k for k in st.session_state if str(k).startswith("w.")]:
        del st.session_state[k]


def _w(widget, label: str, obj, attr: str, key: str, **kw):
    """A widget bound to a config field (kept in the dataclass on every rerun)."""
    k = "w." + key
    if k not in st.session_state:
        st.session_state[k] = getattr(obj, attr)
    setattr(obj, attr, widget(label, key=k, **kw))
    return getattr(obj, attr)


def _path_note(path: str, want: str) -> None:
    if not path:
        return
    q = Path(path)
    if not q.exists():
        st.caption(f":red[not found: {q.resolve()}]")
    elif want == "dir" and not q.is_dir():
        st.caption(":red[expected a folder]")
    elif want == "file" and q.is_dir():
        st.caption(":red[expected a file]")
    else:
        st.caption(f":green[found] {'folder' if q.is_dir() else f'{q.stat().st_size / 1e6:.1f} MB'}")


# ----------------------------------------------------------------------------------------------
# Settings panels
# ----------------------------------------------------------------------------------------------

def input_panel(cfg: C.RunConfig) -> None:
    i = cfg.input
    kind = _w(st.radio, "Where does the signal come from?", i, "kind", "input.kind",
              options=list(KINDS), format_func=KINDS.get, horizontal=False)
    if kind in ("raw", "store", "replay"):
        want = {"raw": "dir", "store": "file", "replay": "any"}[kind]
        hint = {"raw": "folder of ASCII capture files", "store": ".h5 file from 'gammaspec convert'",
                "replay": "raw folder or .h5 store"}[kind]
        c1, c2 = st.columns(2)
        with c1:
            _w(st.text_input, f"Source run ({hint})", i, "source", "input.source")
            _path_note(i.source, want)
        with c2:
            _w(st.text_input, "Background run (optional, same kind)", i, "background", "input.background")
            _path_note(i.background, want)
        _w(st.number_input, "Only the first N captures of each run (0 = all)", i, "limit", "input.limit",
           min_value=0, step=10)
    elif kind == "pipe":
        _w(st.text_input, "File or FIFO with binary samples (a digitiser writes into it)", i, "source", "input.source")
        c1, c2, c3 = st.columns(3)
        with c1:
            _w(st.number_input, "Sample rate (Hz)", i, "pipe_fs", "input.pipe_fs", min_value=1.0, format="%.0f")
        with c2:
            _w(st.selectbox, "Sample format", i, "pipe_dtype", "input.pipe_dtype", options=["f4", "f8", "i2", "i4"])
        with c3:
            _w(st.number_input, "Volts per unit (ADC step for integer codes)", i, "pipe_scale", "input.pipe_scale",
               min_value=0.0, format="%.6g")
        _w(st.text_input, "Background run (raw folder or .h5, optional)", i, "background", "input.background")
    else:
        sim, sig = cfg.simulation, cfg.simulation.signal
        c1, c2, c3 = st.columns(3)
        with c1:
            _w(st.selectbox, "Source nuclide", sig, "nuclide", "sim.nuclide", options=["co60", "cs137", "k40", "none"])
            _w(st.number_input, "Source pulses per second", sig, "rate_cps", "sim.rate", min_value=0.0, step=100.0)
            _w(st.number_input, "Background pulses per second", sig, "background_cps", "sim.bkg", min_value=0.0, step=50.0)
        with c2:
            _w(st.number_input, "Source run (s of signal)", sim, "duration_s", "sim.dur", min_value=0.1, step=10.0)
            _w(st.number_input, "Background run (s, 0 = none)", sim, "background_duration_s", "sim.bdur", min_value=0.0, step=10.0)
            _w(st.number_input, "Random seed", sig, "seed", "sim.seed", min_value=0, step=1)
        with c3:
            _w(st.number_input, "Resolution: FWHM at 662 keV (%)", sig, "fwhm_pct_662", "sim.fwhm", min_value=0.1, step=0.5)
            _w(st.number_input, "Noise (V rms)", sig, "noise_v", "sim.noise", min_value=0.0, format="%.5f", step=0.0001)
            _w(st.number_input, "Gain (keV per V)", sig, "gain_kev_per_v", "sim.gain", min_value=1.0, step=10.0)
    if kind in ("replay", "simulation"):
        c1, c2 = st.columns(2)
        with c1:
            _w(st.checkbox, "Real time (a live stream: frames are dropped if the PC falls behind)", cfg.stream,
               "realtime", "stream.realtime")
        with c2:
            _w(st.number_input, "Speed (× real time)", cfg.stream, "speed", "stream.speed", min_value=0.1, step=0.5)


def processing_panel(cfg: C.RunConfig) -> None:
    p = cfg.pipeline
    st.markdown("**Template (the expected pulse shape)**")
    c1, c2 = st.columns(2)
    with c1:
        _w(st.text_input, "Template file (.h5; empty = build one from the data)", cfg.template, "path", "template.path")
        _path_note(cfg.template.path, "file")
    with c2:
        if cfg.input.kind != "simulation":
            _w(st.selectbox, "Build it from", cfg.template, "build_from", "template.build_from",
               options=["background", "source"], disabled=bool(cfg.template.path))
        else:
            st.caption("Simulation: the simulator's own pulse shape is used.")
    st.markdown("**Stages** (switch off to see what each one contributes)")
    c1, c2, c3 = st.columns(3)
    with c1:
        _w(st.checkbox, "Interference removal", p.preprocess, "line_removal", "pp.line", help=STEPS_HELP["line_removal"])
        _w(st.checkbox, "Pole-zero correction", p.preprocess, "pole_zero", "pp.pz", help=STEPS_HELP["pole_zero"])
        _w(st.checkbox, "Baseline tracker", p.preprocess, "baseline_tracker", "pp.bt", help=STEPS_HELP["baseline_tracker"])
    with c2:
        _w(st.selectbox, "Detection filter", p.detect, "method", "det.method",
           options=["optimal", "matched", "sharpened", "edge"])
        _w(st.number_input, "Detection threshold (σ)", p.detect, "threshold", "det.thr", min_value=1.0, step=0.5)
    with c3:
        _w(st.checkbox, "Pile-up solver", p.solver, "enabled", "sol.on", help=STEPS_HELP["solver"])
        _w(st.checkbox, "Search for hidden pulses", p.solver, "residual_search", "sol.rs", help=STEPS_HELP["residual_search"],
           disabled=not p.solver.enabled)
    st.markdown("**After processing**")
    c1, c2, c3 = st.columns(3)
    with c1:
        _w(st.checkbox, "Calibrate and subtract background", cfg.spectrum, "enabled", "spec.on")
    with c2:
        _w(st.checkbox, "Benchmark (naive count, analog chain)", cfg.benchmark, "enabled", "bench.on",
           disabled=cfg.input.kind in ("simulation", "pipe"))
    with c3:
        _w(st.number_input, "Spectrum bin width (keV)", cfg.spectrum, "bin_kev", "spec.bin", min_value=0.5, step=1.0)


def output_panel(cfg: C.RunConfig) -> None:
    c1, c2, c3 = st.columns([3, 1, 1])
    with c1:
        _w(st.text_input, "Output folder", cfg, "output", "output")
    with c2:
        _w(st.number_input, "Worker processes (0 = auto)", cfg, "workers", "workers", min_value=0, max_value=64)
    with c3:
        _w(st.checkbox, "Recompute everything", cfg, "overwrite", "overwrite")


def advanced_panel(cfg: C.RunConfig) -> None:
    st.caption("The complete run file. Edit any setting and press Apply; the widgets above follow. "
               "Meanings: docs/configuration.md (and the comment on each line).")
    current = C.dumps(cfg)
    text = st.text_area("Run file (TOML)", current, height=420, key=f"toml.{hash(current)}")
    c1, c2, c3 = st.columns([1, 1, 3])
    with c1:
        if st.button("Apply"):
            try:
                _replace_cfg(C.from_dict(tomllib.loads(text)))
                st.rerun()
            except (C.ConfigError, tomllib.TOMLDecodeError) as exc:
                st.error(str(exc))
    with c2:
        st.download_button("Download", C.dumps(cfg), file_name="run.toml", mime="text/plain")
    with c3:
        up = st.file_uploader("Load a run file", type=["toml", "json"], label_visibility="collapsed")
        if up is not None and st.session_state.get("loaded") != up.file_id:
            try:
                d = json.loads(up.getvalue()) if up.name.endswith(".json") else tomllib.loads(up.getvalue().decode())
                st.session_state.loaded = up.file_id
                _replace_cfg(C.from_dict(d))
                st.rerun()
            except (C.ConfigError, ValueError) as exc:
                st.error(str(exc))


# ----------------------------------------------------------------------------------------------
# Running
# ----------------------------------------------------------------------------------------------

def _proc():
    return st.session_state.get("proc")


def _running() -> bool:
    p = _proc()
    return p is not None and p.poll() is None


def start_run(cfg: C.RunConfig) -> None:
    out = Path(cfg.output)
    out.mkdir(parents=True, exist_ok=True)
    path = C.save(cfg, out / "ui_run.toml")
    log = open(out / "ui_run.log", "w", encoding="utf-8")
    env = dict(os.environ, PYTHONUTF8="1", PYTHONUNBUFFERED="1", MPLBACKEND="Agg")
    flags = subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0
    st.session_state.proc = subprocess.Popen([sys.executable, "-m", "gammaspec", "run", str(path)],
                                             stdout=log, stderr=subprocess.STDOUT, env=env, creationflags=flags)
    st.session_state.run_out = str(out)
    st.session_state.run_started = time.time()


def stop_run() -> None:
    p = _proc()
    if p is None or p.poll() is not None:
        return
    if os.name == "nt":                          # also end the worker processes
        subprocess.run(["taskkill", "/T", "/F", "/PID", str(p.pid)], capture_output=True)
    else:
        p.terminate()


def _plot_live(snap: dict, title: str):
    import matplotlib.pyplot as plt

    from gammaspec import plotting as P

    counts = np.asarray(snap["counts"], float)
    if "edges_kev" in snap:
        x, unit = np.asarray(snap["edges_kev"]), "keV"
    else:
        x, unit = np.asarray(snap["edges_v"]) * 1e3, "mV"
    y = np.maximum(counts / max(snap["live_time_s"], 1e-9) / np.diff(x), 1e-6)
    fig, ax = plt.subplots(figsize=(9, 3.4))
    ax.stairs(y, x, color=P.SERIES[0], linewidth=1.0)
    ax.set_yscale("log")
    ax.set_xlabel(f"energy ({unit})" if unit == "keV" else "pulse amplitude (mV)")
    ax.set_ylabel(f"counts / live s / {unit}")
    ax.set_title(f"{title}: {snap['good']:,} good events in {snap['live_time_s']:.1f} s live"
                 + ("" if unit == "keV" else " (not calibrated yet)"), fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return fig


@st.fragment(run_every=2.0)
def progress_view() -> None:
    out = Path(st.session_state.get("run_out", ""))
    p = _proc()
    if p is None:
        return
    rc = p.poll()
    elapsed = time.time() - st.session_state.get("run_started", time.time())
    if rc is None:
        st.info(f"Running… {elapsed:.0f} s")
    elif rc == 0:
        st.success(f"Finished in {elapsed:.0f} s. Results below.")
    else:
        st.error(f"The run stopped (exit code {rc}). See the log.")
    for run in ("source", "background"):
        f = out / "live" / run / "live_spectrum.json"
        if f.exists():
            try:
                snap = json.loads(f.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                continue                          # being rewritten: next refresh
            st.pyplot(_plot_live(snap, f"live spectrum, {run} run"), clear_figure=True)
    log = out / "ui_run.log"
    if log.exists():
        lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
        st.code("\n".join(lines[-25:]) or "(starting…)", language=None)
    if rc is not None and not st.session_state.get("shown_done"):
        st.session_state.shown_done = True
        st.rerun()                                 # redraw the whole page with the results


def results_view(out: Path) -> None:
    summ = out / "run_summary.json"
    if not summ.exists():
        st.caption("No results in this output folder yet.")
        return
    s = json.loads(summ.read_text(encoding="utf-8"))
    for w in s.get("warnings", []):
        st.warning(w)
    st.markdown("**Steps**")
    st.dataframe([{"step": k, "seconds": v["seconds"], "output": v["output"]} for k, v in s["steps"].items()],
                 hide_index=True, use_container_width=True)
    spec = out / "spectrum" / "spectrum_summary.json"
    if spec.exists():
        sp = json.loads(spec.read_text(encoding="utf-8"))
        st.markdown("**Spectrum**")
        cal = sp["calibration_source"]["coeffs"]
        st.write(f"Calibration (source run): E = {cal[0]:.1f} keV + {cal[1]:.1f} keV/V × A")
        rows = []
        for kind, lines in sp["lines"].items():
            for r in lines:
                rows.append({"run": kind, "line": f"{r['label']} {r['energy_kev']:.0f} keV",
                             "counts/s": round(r["rate_cps"], 2), "±": round(r["rate_err_cps"], 2),
                             "FWHM %": round(r.get("fwhm_pct") or float("nan"), 2)})
        st.dataframe(rows, hide_index=True, use_container_width=True)
        for img in ("spectrum.png", "spectrum_zoom.png", "calibration.png"):
            if (out / "spectrum" / img).exists():
                st.image(str(out / "spectrum" / img))
    chk = out / "simulation_check.json"
    if chk.exists():
        st.markdown("**Simulation: the result against the truth**")
        for run, c in json.loads(chk.read_text(encoding="utf-8")).items():
            st.write(f"{run}: {c['true_pulses']:,} true pulses, {c['detected']:,} detected; false "
                     f"{c['false_rate_cps']:.1f} /s; amplitude bias {c['amplitude_bias_pct']:+.2f} %; "
                     f"close pairs resolved {100 * (c['close_pairs_resolved'] or 0):.1f} %")
            st.dataframe([{"energy (keV)": f"{e['energy_kev'][0]}–{e['energy_kev'][1]}", "true pulses": e["true"],
                           "found %": None if e["efficiency"] is None else round(100 * e["efficiency"], 1)}
                          for e in c["efficiency_by_energy"]], hide_index=True)
    bench = out / "benchmark" / "benchmark_summary.json"
    if bench.exists():
        b = json.loads(bench.read_text(encoding="utf-8"))
        st.markdown("**Benchmark: the same pulses through other chains**")
        st.dataframe([{"chain": r["name"], "FWHM 1173 %": round(r["fwhm_pct"][0], 2), "FWHM 1332 %": round(r["fwhm_pct"][1], 2),
                       "valley/peak": round(r["valley_to_peak"], 3), "pile-up /s": round(r["pileup_continuum_cps"], 1)}
                      for r in b["chains"]], hide_index=True, use_container_width=True)
        if (out / "benchmark" / "benchmark_spectra.png").exists():
            st.image(str(out / "benchmark" / "benchmark_spectra.png"))


# ----------------------------------------------------------------------------------------------
# Page
# ----------------------------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(page_title="gammaspec", layout="wide")
    st.title("gammaspec")
    st.caption(f"Digital gamma-ray spectroscopy, v{__version__}: from ADC samples to a background-subtracted "
               "spectrum, in software.")
    cfg = _cfg()
    t_in, t_proc, t_all = st.tabs(["1 · Input", "2 · Processing", "All settings"])
    with t_in:
        input_panel(cfg)
    with t_proc:
        processing_panel(cfg)
    with t_all:
        advanced_panel(cfg)
    st.divider()
    st.subheader("3 · Run")
    output_panel(cfg)
    try:
        cfg = C.from_dict(C.to_dict(cfg))        # full validation of the edited values
        warns = C.check(cfg)
        valid = True
    except C.ConfigError as exc:
        st.error(str(exc))
        warns, valid = [], False
    for w in warns:
        st.warning(w)
    c1, c2, _ = st.columns([1, 1, 4])
    running = _running()
    with c1:
        if st.button("Run", type="primary", disabled=running or not valid):
            st.session_state.shown_done = False
            start_run(cfg)
            st.rerun()
    with c2:
        if st.button("Stop", disabled=not running):
            stop_run()
            st.rerun()
    if _proc() is not None:
        progress_view()
    st.divider()
    st.subheader("4 · Results")
    results_view(Path(st.session_state.get("run_out", cfg.output)))


def launch(port: int | None = None, open_browser: bool = True, address: str = "localhost") -> int:
    """Start the Streamlit server on this file (``gammaspec app``).

    It listens on this computer only (``localhost``) by default: the app reads files and starts
    processes, so it should not be reachable from the network unless that is intended.
    """
    cmd = [sys.executable, "-m", "streamlit", "run", str(Path(__file__).resolve()),
           "--server.headless", "false" if open_browser else "true", "--browser.gatherUsageStats", "false",
           "--server.address", address]
    if port:
        cmd += ["--server.port", str(port)]
    return subprocess.call(cmd)


if __name__ == "__main__" or st.runtime.exists():
    main()
