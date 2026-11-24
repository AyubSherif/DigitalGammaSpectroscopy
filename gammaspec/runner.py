"""Run the whole chain from a run configuration (``gammaspec run run.toml``).

Steps, each selected by the configuration:

    1. input      raw folder -> HDF5 store (input.kind = "raw"), or an existing store; streamed
                  kinds ("replay", "simulation", "pipe") need no store
    2. template   template.path, or built by the characterisation step (simulation: the
                  simulator's own pulse shape)
    3. process    -> event file, for the source and (if given) the background run: batch
                  (pipeline.run) for "raw" / "store", the stream processor otherwise
    4. spectrum   calibration, spectra, background subtraction   (spectrum.enabled)
    5. benchmark  naive / analog / gammaspec comparison            (benchmark.enabled)
    +  simulation the events are scored against the simulator's truth (simulation_check.json)

Outputs go to ``output``:

    config.toml        the effective configuration (what was actually run)
    run_summary.json   steps, timings, outputs
    stores/  template/  events/  spectrum/  benchmark/

Each output records what it was made from (``<output>.stamp.json``: inputs and settings).
On a rerun it is reused when that is unchanged and recomputed when it changed, so after
editing one setting only the affected steps run again. ``overwrite = true`` recomputes
everything.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from . import config as C
from .store import CaptureStore


class StepError(RuntimeError):
    """A step could not run with these inputs (message says which and why)."""


def _say(verbose, *a):
    if verbose:
        print(*a, flush=True)


def _stamp_path(path: Path) -> Path:
    return path.with_name(path.name + ".stamp.json")


def _fresh(path: Path, stamp: dict, cfg) -> bool:
    """True if ``path`` exists and was made from exactly ``stamp`` (and overwrite is off)."""
    sp = _stamp_path(path)
    if cfg.overwrite or not path.exists() or not sp.exists():
        return False
    try:
        return json.loads(sp.read_text(encoding="utf-8")) == json.loads(json.dumps(stamp))
    except json.JSONDecodeError:
        return False


def _write_stamp(path: Path, stamp: dict) -> None:
    _stamp_path(path).write_text(json.dumps(stamp, indent=2), encoding="utf-8")


def _path_kind(path: str) -> str:
    return "raw" if Path(path).is_dir() else "store"


def _store_for(run: str, src: str, cfg: C.RunConfig, out: Path, verbose: bool, kind: str | None = None) -> Path:
    """Step 1: the capture store of one run (converted from raw if needed)."""
    i = cfg.input
    kind = kind or i.kind
    if kind == "store":
        return Path(src)
    from .rawio import build_index
    from .store import convert

    headers = build_index(src, i.pattern)
    if i.limit:
        headers = headers[: i.limit]
    dst = out / "stores" / f"{run}.h5"
    stamp = dict(source=str(Path(src).resolve()), pattern=i.pattern, files=[h.path.name for h in headers])
    if i.reuse_converted and _fresh(dst, stamp, cfg):
        _say(verbose, f"[input] {run}: reusing {dst} ({len(headers)} captures)")
        return dst
    _say(verbose, f"[input] {run}: converting {len(headers)} captures from {src}")
    convert(headers, dst, workers=cfg.workers or None, overwrite=True)
    _write_stamp(dst, stamp)
    return dst


def _template(cfg: C.RunConfig, stores: dict, out: Path, verbose: bool) -> Path:
    """Step 2: the pulse template."""
    if cfg.template.path:
        return Path(cfg.template.path)
    from .characterize import characterize

    dst = out / "template" / "template.h5"
    src = stores[cfg.template.build_from]
    stamp = dict(store=_identity(src), characterize=C.to_dict(cfg.template.characterize))
    if _fresh(dst, stamp, cfg):
        _say(verbose, f"[template] reusing {dst}")
        return dst
    _say(verbose, f"[template] building from the {cfg.template.build_from} run")
    try:
        characterize(src, dst.parent, cfg.template.characterize, verbose=verbose)
    except ValueError as exc:
        raise StepError(f"template: could not be built from the {cfg.template.build_from} run ({exc}); "
                        "use more captures or give template.path") from exc
    _write_stamp(dst, stamp)
    return dst


def _identity(path: Path) -> dict:
    """What identifies a store's content: path, size, modification time."""
    st = Path(path).stat()
    return dict(path=str(Path(path).resolve()), bytes=st.st_size, mtime=st.st_mtime)


def _events(run: str, store: Path, template: Path, cfg: C.RunConfig, out: Path, verbose: bool) -> Path:
    """Step 3: store -> event file."""
    from .pipeline import run as run_pipeline

    dst = out / "events" / f"{run}.h5"
    dst.parent.mkdir(parents=True, exist_ok=True)
    limit = cfg.input.limit or None
    stamp = dict(store=_identity(store), template=_identity(template), limit=limit,
                 pipeline=C.to_dict(cfg.pipeline))
    if _fresh(dst, stamp, cfg):
        _say(verbose, f"[process] {run}: reusing {dst} (same inputs and settings)")
        return dst
    if dst.exists():
        _say(verbose, f"[process] {run}: inputs or settings changed, recomputing {dst}")
    _say(verbose, f"[process] {run}: {store} -> {dst}")
    r = run_pipeline(store, template, dst, cfg.pipeline, limit=limit, overwrite=True, verbose=verbose,
                     workers=cfg.workers or None)
    _write_stamp(dst, stamp)
    _say(verbose, f"[process] {run}: {r['events']} events from {r['captures']} captures in {r['elapsed_s']:.0f} s")
    return dst


def _stream_events(run: str, make_source, stamp_input, template: Path, cfg: C.RunConfig, out: Path,
                   verbose: bool, on_update=None, after=None) -> Path:
    """Step 3 for streamed input: source -> stream processor -> event file."""
    from .stream import process_stream

    dst = out / "events" / f"{run}.h5"
    dst.parent.mkdir(parents=True, exist_ok=True)
    s = cfg.stream
    reusable = stamp_input is not None and not s.realtime
    stamp = dict(input=stamp_input, template=_identity(template), pipeline=C.to_dict(cfg.pipeline),
                 frames=dict(frame=s.frame, pre=s.pre, post=s.post))
    if reusable and _fresh(dst, stamp, cfg):
        _say(verbose, f"[process] {run}: reusing {dst} (same inputs and settings)")
        return dst
    source = make_source()
    _say(verbose, f"[process] {run}: streaming {source.description} -> {dst}")
    r = process_stream(source, template, dst, cfg.pipeline, s, workers=cfg.workers or None, verbose=verbose,
                       live_dir=out / "live" / run, on_update=on_update, kind=run, overwrite=True)
    if after is not None:
        after(source)
    if reusable:
        _write_stamp(dst, stamp)
    _say(verbose, f"[process] {run}: {r['events']} events in {r['stream_s']:.1f} s of signal "
                  f"({r['speed']:.2f}x real time, {r['dropped_frames']} frames dropped)")
    return dst


def _sim_template(cfg: C.RunConfig, out: Path, verbose: bool) -> Path:
    """Simulation without template.path: the simulator's own pulse shape is the template."""
    from .simulate import load_pulse

    if cfg.template.path:
        return Path(cfg.template.path)
    dst = out / "template" / "template.h5"
    dst.parent.mkdir(parents=True, exist_ok=True)
    stamp = dict(pulse=cfg.simulation.signal.pulse_template or "analytic")
    if not _fresh(dst, stamp, cfg):
        load_pulse(cfg.simulation.signal).save(dst)
        _write_stamp(dst, stamp)
    _say(verbose, f"[template] the simulator's pulse shape: {dst}")
    return dst


def run(cfg: C.RunConfig, verbose: bool = True, on_update=None) -> dict:
    """Run every selected step. Returns (and writes) the run summary.

    ``on_update(status, live_spectrum)`` is called periodically while a stream is processed
    (e.g. a live plot or a UI)."""
    warnings_ = C.check(cfg)
    for w in warnings_:
        _say(verbose, f"warning: {w}")
    out = Path(cfg.output)
    out.mkdir(parents=True, exist_ok=True)
    C.save(cfg, out / "config.toml")
    summary: dict = dict(config=str(out / "config.toml"), warnings=warnings_, steps={})
    t0 = time.perf_counter()

    def step(name, fn):
        t = time.perf_counter()
        try:
            r = fn()
        except StepError:
            raise
        except ValueError as exc:
            hint = (" With few captures (input.limit) a run may have too few counts to calibrate; "
                    "use more captures." if name in ("spectrum", "benchmark") else "")
            raise StepError(f"{name}: {exc}.{hint}") from exc
        summary["steps"][name] = dict(seconds=round(time.perf_counter() - t, 1), output=str(r) if r else None)
        return r

    i, kind = cfg.input, cfg.input.kind
    runs = {"source": i.source}
    if i.background and kind != "simulation":
        runs["background"] = i.background
    stores: dict = {}
    events: dict = {}
    if kind in ("raw", "store"):
        stores = {k: step(f"input_{k}", lambda k=k, v=v: _store_for(k, v, cfg, out, verbose)) for k, v in runs.items()}
        template = step("template", lambda: _template(cfg, stores, out, verbose))
        events = {k: step(f"process_{k}", lambda k=k: _events(k, stores[k], template, cfg, out, verbose))
                  for k in stores}
    elif kind == "simulation":
        events = _run_simulation(cfg, out, verbose, step, summary, on_update)
        template = Path(summary["template"])
    else:
        from .sources import open_source

        sopt = dict(limit=i.limit or None, pattern=i.pattern, chunk=cfg.stream.chunk,
                    realtime=cfg.stream.realtime, speed=cfg.stream.speed)
        # stores only where a step needs one: building the template, the benchmark, a batch background
        need = set()
        if not cfg.template.path:
            need.add(cfg.template.build_from)
        if cfg.benchmark.enabled and kind == "replay":
            need.add("source")
        if kind == "pipe" and "background" in runs:
            need.add("background")
        for k in sorted(need & set(runs)):
            stores[k] = step(f"input_{k}", lambda k=k: _store_for(k, runs[k], cfg, out, verbose, _path_kind(runs[k])))
        template = step("template", lambda: _template(cfg, stores, out, verbose))
        for k, v in runs.items():
            if kind == "pipe" and k == "source":
                events[k] = step(f"process_{k}", lambda v=v: _stream_events(
                    "source", lambda: open_source("pipe", v or "-", fs=i.pipe_fs, pipe_dtype=i.pipe_dtype,
                                                  pipe_scale=i.pipe_scale, chunk=cfg.stream.chunk),
                    None, template, cfg, out, verbose, on_update))
            elif kind == "pipe":                              # the background run is a recording: batch
                events[k] = step(f"process_{k}", lambda k=k: _events(k, stores[k], template, cfg, out, verbose))
            else:                                             # replay
                ident = (_identity(v) if _path_kind(v) == "store" else
                         dict(folder=str(Path(v).resolve()), pattern=i.pattern,
                              files=sorted(q.name for q in Path(v).glob(i.pattern))))
                events[k] = step(f"process_{k}", lambda k=k, v=v, ident=ident: _stream_events(
                    k, lambda: open_source("replay", v, **sopt), dict(path=ident, limit=i.limit),
                    template, cfg, out, verbose, on_update))

    if cfg.spectrum.enabled:
        from .analysis import analyze

        s = cfg.spectrum
        _say(verbose, "[spectrum] calibration, spectra" + (", background subtraction" if "background" in events else ""))
        step("spectrum", lambda: (analyze(events["source"], out / "spectrum", events.get("background"),
                                          bin_kev=s.bin_kev, e_max_kev=s.e_max_kev,
                                          reliable_above_kev=s.reliable_above_kev, verbose=verbose),
                                  out / "spectrum")[1])
    if cfg.benchmark.enabled and "source" in stores:
        from .benchmark import run_benchmark

        with CaptureStore(stores["source"]) as st:
            n_done = len(st) if not cfg.input.limit else min(cfg.input.limit, len(st))
        n = n_done if not cfg.benchmark.captures else min(cfg.benchmark.captures, n_done)
        _say(verbose, f"[benchmark] {n} captures")
        step("benchmark", lambda: (run_benchmark(stores["source"], events["source"], out / "benchmark",
                                                 captures=range(n), verbose=verbose), out / "benchmark")[1])

    summary.update(stores={k: str(v) for k, v in stores.items()}, template=str(template),
                   events={k: str(v) for k, v in events.items()}, elapsed_s=round(time.perf_counter() - t0, 1))
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _say(verbose, f"done in {summary['elapsed_s']:.0f} s; outputs in {out}")
    return summary


def _run_simulation(cfg: C.RunConfig, out: Path, verbose: bool, step, summary: dict, on_update) -> dict:
    """Simulated source (and background) runs, streamed; truth saved and compared."""
    from dataclasses import replace

    from .simulate import compare_to_truth, save_truth
    from .sources import Paced, SimulationSource

    sc = cfg.simulation
    template = step("template", lambda: _sim_template(cfg, out, verbose))
    summary["template"] = str(template)
    runs = {"source": (sc.signal, sc.duration_s)}
    if sc.background_duration_s > 0:
        # same detector and room, no source; its own random sequence
        runs["background"] = (replace(sc.signal, nuclide="none", seed=sc.signal.seed + 1000), sc.background_duration_s)
    events, checks = {}, {}
    for k, (sig, dur) in runs.items():
        truth_path = out / "truth" / f"{k}.h5"

        def make(sig=sig, dur=dur):
            src = SimulationSource(sig, dur, cfg.stream.chunk)
            return Paced(src, cfg.stream.speed) if cfg.stream.realtime else src

        def keep_truth(source, sig=sig, truth_path=truth_path):
            sim = (source.source if isinstance(source, Paced) else source).sim
            save_truth(truth_path, sim.truth_table(), sig)

        stamp = dict(simulation=C.to_dict(sig), duration_s=dur)
        if not truth_path.exists():                        # events without their truth: redo them
            _stamp_path(out / "events" / f"{k}.h5").unlink(missing_ok=True)
        events[k] = step(f"process_{k}", lambda k=k, make=make, stamp=stamp, keep_truth=keep_truth: _stream_events(
            k, make, stamp, template, cfg, out, verbose, on_update, after=keep_truth))
        checks[k] = compare_to_truth(events[k], truth_path)
    (out / "simulation_check.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
    summary["simulation_check"] = str(out / "simulation_check.json")
    if verbose:
        for k, c in checks.items():
            eff = ", ".join(f"{e['energy_kev'][0]}-{e['energy_kev'][1]} keV {100 * e['efficiency']:.1f} %"
                            for e in c["efficiency_by_energy"] if e["efficiency"] is not None)
            print(f"[truth] {k}: {c['true_pulses']} true pulses ({c['true_rate_cps']:.0f} /s), {c['detected']} detected; "
                  f"found: {eff}; false {c['false_rate_cps']:.1f} /s; amplitude bias {c['amplitude_bias_pct']:+.2f} % "
                  f"(spread {c['amplitude_spread_pct']:.2f} %); close pairs resolved "
                  f"{100 * (c['close_pairs_resolved'] or 0):.1f} %", flush=True)
    return events
