"""Run the whole chain from a run configuration (``gammaspec run run.toml``).

Steps, each selected by the configuration:

    1. input      raw folder -> HDF5 store (input.kind = "raw"), or an existing store
    2. template   template.path, or built by the characterisation step
    3. process    store -> event file, for the source and (if given) the background run
    4. spectrum   calibration, spectra, background subtraction   (spectrum.enabled)
    5. benchmark  naive / analog / gammaspec comparison            (benchmark.enabled)

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


def _store_for(run: str, src: str, cfg: C.RunConfig, out: Path, verbose: bool) -> Path:
    """Step 1: the capture store of one run (converted from raw if needed)."""
    i = cfg.input
    if i.kind == "store":
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


def run(cfg: C.RunConfig, verbose: bool = True) -> dict:
    """Run every selected step. Returns (and writes) the run summary."""
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

    runs = {"source": cfg.input.source}
    if cfg.input.background:
        runs["background"] = cfg.input.background
    stores = {k: step(f"input_{k}", lambda k=k, v=v: _store_for(k, v, cfg, out, verbose)) for k, v in runs.items()}
    template = step("template", lambda: _template(cfg, stores, out, verbose))
    events = {k: step(f"process_{k}", lambda k=k: _events(k, stores[k], template, cfg, out, verbose))
              for k in stores}

    if cfg.spectrum.enabled:
        from .analysis import analyze

        s = cfg.spectrum
        _say(verbose, "[spectrum] calibration, spectra" + (", background subtraction" if "background" in events else ""))
        step("spectrum", lambda: (analyze(events["source"], out / "spectrum", events.get("background"),
                                          bin_kev=s.bin_kev, e_max_kev=s.e_max_kev,
                                          reliable_above_kev=s.reliable_above_kev, verbose=verbose),
                                  out / "spectrum")[1])
    if cfg.benchmark.enabled:
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
