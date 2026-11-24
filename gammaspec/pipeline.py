"""The processing pipeline: capture store -> pre-processing -> detection -> event list.

    for each capture:
        clean = Preprocessor(...).process(volts)          # [lines] [pole-zero] [tracker]
        det   = Detector(...).detect(clean)               # optimal filter by default
        events: clusters (pile-up), flags, template fit of isolated pulses, times
        snippets of pile-up clusters (for the pile-up solver)
    -> EventWriter (HDF5)

Every stage is configurable or optional (``PipelineConfig``). The configuration is
stored in the event file, so each spectrum can be traced back to its processing.
"""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .baseline import design_from_template, filter_template
from .detect import DetectConfig, Detector
from .events import (CAPTURE_DTYPE, EVENT_DTYPE, FLAG_AFTER_LARGE, FLAG_BAD_FIT, FLAG_EDGE,
                     FLAG_LARGE, FLAG_PILEUP, FLAG_WARMUP, EventWriter)
from .preprocess import PreprocessConfig, Preprocessor
from .store import CaptureStore
from .template import PulseTemplate, extract_windows, fit_pulses


@dataclass
class EventConfig:
    cluster_gap: int = 60             # pulses closer than this share a fit window -> pile-up cluster
    fit_pre: int = 10                 # template fit window: samples before the peak ...
    fit_post: int = 40                # ... and after it
    snippets: str = "pileup"          # store clean waveforms of: "none" | "pileup" | "all" clusters
    snippet_pre: int = 20
    snippet_post: int = 80
    large_v: float = 1.5              # above this: FLAG_LARGE (cosmic muons; shape not verified)
    large_shadow: int = 200           # events this many samples after a large pulse: FLAG_AFTER_LARGE
    bad_fit_sigma: float = 4.0        # bad fit if rms > bad_fit_sigma * noise + bad_fit_frac * amplitude
    bad_fit_frac: float = 0.02


@dataclass
class PipelineConfig:
    preprocess: PreprocessConfig = field(default_factory=PreprocessConfig)
    detect: DetectConfig = field(default_factory=DetectConfig)
    events: EventConfig = field(default_factory=EventConfig)

    def to_dict(self) -> dict:
        return asdict(self)


class Pipeline:
    def __init__(self, template: PulseTemplate, fs: float, config: PipelineConfig | None = None):
        self.cfg = config or PipelineConfig()
        p = self.cfg.preprocess
        self.fs = fs
        self.pz = design_from_template(template) if p.pole_zero else None
        # the template as it looks after pre-processing (pole-zero changes the shape)
        self.template = filter_template(template, self.pz.b, self.pz.a) if self.pz else template
        self.fit_baseline = not (p.pole_zero and p.baseline_tracker)
        self.pre = Preprocessor(fs, p, self.pz)
        self.detector = Detector(self.template, self.cfg.detect)
        self._cluster0 = 0

    def process_capture(self, volts: np.ndarray, capture: int, t0_s: float):
        """One capture -> (events, capture_row, snippets)."""
        c, ev = self.cfg, self.cfg.events
        dt = 1.0 / self.fs
        self.pre.reset()
        y = self.pre.process(volts)
        det = self.detector.detect(y)
        idx = det.index
        n, ne = y.size, idx.size

        e = np.zeros(ne, EVENT_DTYPE)
        e["capture"] = capture
        e["index"] = idx
        e["filter_amplitude"] = det.filter_amplitude
        e["amplitude"] = np.nan
        e["fit_rms"] = np.nan
        gaps = np.diff(idx)
        e["gap_prev"] = np.r_[-1, gaps] if ne else []
        e["gap_next"] = np.r_[gaps, -1] if ne else []

        # clusters: a new cluster starts when the gap to the previous pulse is large enough
        new = np.r_[True, gaps >= ev.cluster_gap] if ne else np.zeros(0, bool)
        cl = np.cumsum(new) - 1
        sizes = np.bincount(cl) if ne else np.zeros(0, int)
        e["cluster"] = self._cluster0 + cl
        e["cluster_size"] = sizes[cl] if ne else []
        self._cluster0 += int(cl.max() + 1) if ne else 0

        flags = np.zeros(ne, np.uint16)
        flags[e["cluster_size"] > 1] |= FLAG_PILEUP
        flags[idx < c.preprocess.warmup_samples] |= FLAG_WARMUP
        fit_lo = ev.fit_pre + 2
        flags[(idx - fit_lo < 0) | (idx + ev.fit_post >= n)] |= FLAG_EDGE

        # template fit of isolated pulses in the valid region
        iso = (flags & (FLAG_PILEUP | FLAG_EDGE)) == 0
        if iso.any():
            W = extract_windows(y, idx[iso], fit_lo, ev.fit_post)
            f = fit_pulses(W, self.template, fit_lo, fit_range=(-ev.fit_pre, ev.fit_post),
                           fit_baseline=self.fit_baseline)
            e["amplitude"][iso] = f.amplitude
            e["shift"][iso] = f.shift
            e["fit_rms"][iso] = f.rms
            noise = self.pre.info.sigma_out_v
            bad = f.rms > ev.bad_fit_sigma * noise + ev.bad_fit_frac * np.abs(f.amplitude)
            flags[np.flatnonzero(iso)[bad]] |= FLAG_BAD_FIT

        # very large pulses and their shadow
        big = np.maximum(np.nan_to_num(e["amplitude"]), e["filter_amplitude"]) > ev.large_v
        flags[big] |= FLAG_LARGE
        for i in np.flatnonzero(big):
            shadow = (idx > idx[i]) & (idx <= idx[i] + ev.large_shadow)
            flags[shadow] |= FLAG_AFTER_LARGE
        e["flags"] = flags
        e["time_s"] = t0_s + (idx + e["shift"]) * dt

        # snippets for the pile-up solver
        snippets = []
        if ev.snippets != "none" and ne:
            want = sizes > 1 if ev.snippets == "pileup" else np.ones(sizes.size, bool)
            first = np.r_[0, np.cumsum(sizes)[:-1]]
            for k in np.flatnonzero(want):
                a = int(idx[first[k]]) - ev.snippet_pre
                b = int(idx[first[k] + sizes[k] - 1]) + ev.snippet_post
                a, b = max(a, 0), min(b, n)
                snippets.append((int(e["cluster"][first[k]]), capture, a, y[a:b]))

        row = np.zeros(1, CAPTURE_DTYPE)
        info = self.pre.info
        row[0] = (capture, t0_s, (n - c.preprocess.warmup_samples - ev.fit_post - 1) * dt, info.offset_v,
                  info.sigma_in_v, info.sigma_out_v, det.sigma, det.threshold_v, ne, info.line_freqs_hz.size)
        return e, row, snippets


_WORKER: dict = {}


def _worker_init(store_path, template_path, cfg):
    _WORKER["store"] = CaptureStore(store_path)
    _WORKER["pipe"] = Pipeline(PulseTemplate.load(template_path), 1.0 / _WORKER["store"].dt_s, cfg)


def _worker_run(i):
    s, p = _WORKER["store"], _WORKER["pipe"]
    p._cluster0 = 0                       # cluster ids are made global by the writer
    return p.process_capture(s.volts(i), i, float(s.time_s[i]))


def run(store_path, template_path, out_path, config: PipelineConfig | None = None,
        limit: int | None = None, overwrite: bool = False, verbose: bool = True,
        workers: int | None = None) -> dict:
    """Process a capture store into an event file. Returns a summary dict.

    Captures are independent, so they are processed in parallel worker processes
    (``workers``; default: CPU cores - 1, at most 8; 1 = in this process). Results are
    written in capture order, and cluster ids are made unique across the file.
    """
    import os
    from concurrent.futures import ProcessPoolExecutor

    cfg = config or PipelineConfig()
    template = PulseTemplate.load(template_path)
    t_start = time.perf_counter()
    with CaptureStore(store_path) as store:
        fs = 1.0 / store.dt_s
        pipe = Pipeline(template, fs, cfg)
        n_cap = len(store) if limit is None else min(limit, len(store))
        n_samples = store.n_samples
        attrs = dict(config=cfg.to_dict(), template=str(template_path), source_store=str(store_path),
                     lsb_volts=store.lsb_volts, dt_s=store.dt_s, fit_baseline=pipe.fit_baseline,
                     pole_zero=None if pipe.pz is None else dict(b=pipe.pz.b.tolist(), a=pipe.pz.a.tolist()))
        workers = workers or max(1, min(8, (os.cpu_count() or 2) - 1))
        if workers == 1:
            _WORKER.update(store=store, pipe=pipe)
            results = map(_worker_run, range(n_cap))
            pool = None
        else:
            pool = ProcessPoolExecutor(workers, initializer=_worker_init,
                                       initargs=(str(store_path), str(template_path), cfg))
            results = pool.map(_worker_run, range(n_cap), chunksize=2)
        n_events, cluster0 = 0, 0
        try:
            with EventWriter(out_path, attrs, overwrite=overwrite) as w:
                for i, (ev, row, sn) in enumerate(results):
                    ev["cluster"] += cluster0
                    sn = [(cl + cluster0, cap, a, wv) for cl, cap, a, wv in sn]
                    cluster0 = int(ev["cluster"].max()) + 1 if ev.size else cluster0
                    w.append(ev, row, sn)
                    n_events += ev.size
                    if verbose and (i + 1) % max(1, n_cap // 10) == 0:
                        el = time.perf_counter() - t_start
                        print(f"  {i + 1}/{n_cap} captures, {n_events} events, "
                              f"{(i + 1) * n_samples / el / 1e6:.1f} MS/s", flush=True)
        finally:
            if pool:
                pool.shutdown(cancel_futures=True)
            _WORKER.clear()
        raw_samples = n_cap * n_samples
    elapsed = time.perf_counter() - t_start
    size = Path(out_path).stat().st_size
    return dict(captures=n_cap, events=n_events, elapsed_s=elapsed,
                samples_per_s=raw_samples / elapsed, file_bytes=size,
                text_bytes=raw_samples * 14)  # 14 bytes per sample line in the ASCII files


def summarize(event_path) -> dict:
    """Counts, rates and flag statistics of an event file."""
    from .events import FLAG_NAMES, EventFile
    ef = EventFile(event_path)
    ev, live = ef.events, ef.live_time_s
    valid = (ev["flags"] & (FLAG_WARMUP | FLAG_EDGE)) == 0
    out = dict(events=int(ev.size), live_time_s=live, rate_cps=float(valid.sum() / live),
               clean=int(np.sum(ev["flags"] == 0)), snippets=ef.n_snippets,
               flag_fractions_of_valid={name: float(np.mean((ev["flags"][valid] & bit) != 0))
                                        for bit, name in FLAG_NAMES.items()},
               median_filter_sigma_mv=float(np.median(ef.captures["filter_sigma_v"]) * 1e3),
               median_threshold_mv=float(np.median(ef.captures["threshold_v"]) * 1e3))
    return out


def plot_event_spectrum(event_path, png_path) -> None:
    from . import plotting as P
    from .events import FLAG_PILEUP, EventFile
    ef = EventFile(event_path)
    ev = ef.events
    clean = ev[ev["flags"] == 0]["amplitude"]
    pile = ev[(ev["flags"] & FLAG_PILEUP) != 0]["filter_amplitude"]
    hi = float(np.quantile(clean, 0.995) * 1.3) if clean.size else 1.0
    edges = np.linspace(0, hi, 700)
    fig, ax = P.figure(1, 1, height_per_row=4.0)
    a = ax[0, 0]
    a.stairs(np.histogram(clean, edges)[0], edges, color=P.SERIES[0], linewidth=1.1,
             label=f"isolated pulses, template fit ({clean.size})")
    a.stairs(np.histogram(pile, edges)[0], edges, color=P.SERIES[1], linewidth=1.0,
             label=f"pile-up pulses, filter estimate only ({pile.size})")
    a.set_yscale("log")
    a.set_xlabel("pulse amplitude (V)")
    a.set_ylabel("counts per bin")
    a.set_title(f"Event spectrum ({ef.live_time_s:.1f} s live time; pile-up not yet resolved)")
    a.legend()
    P.save(fig, png_path)
