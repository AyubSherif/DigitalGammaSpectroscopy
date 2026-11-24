"""The processing pipeline: capture store -> pre-processing -> detection -> event list.

    for each capture:
        clean = Preprocessor(...).process(volts)          # [lines] [pole-zero] [tracker]
        det   = Detector(...).detect(clean)               # optimal filter by default
        events: clusters (pile-up), flags, template fit of isolated pulses
        snippets of pile-up clusters
        [pile-up solver]: amplitudes and times of pulses in clusters, hidden pulses
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
import warnings

from .detect import DetectConfig, Detector, noise_autocovariance
from .events import (CAPTURE_DTYPE, EVENT_DTYPE, FLAG_AFTER_LARGE, FLAG_BAD_FIT, FLAG_EDGE,
                     FLAG_FIT_ADDED, FLAG_LARGE, FLAG_PILEUP, FLAG_WARMUP, EventWriter)
from .signal import pulse_free_mask
from .solver import PileupSolver, SolverConfig
from .preprocess import PreprocessConfig, Preprocessor
from .store import CaptureStore
from .template import PulseTemplate, extract_windows, fit_pulses


@dataclass
class EventConfig:
    cluster_gap: int = 60             # pulses closer than this share a fit window -> pile-up cluster
    fit_pre: int = 10                 # template fit window: samples before the peak ...
    fit_post: int = 40                # ... and after it
    snippets: str = "pileup"          # store clean waveforms of: "none" | "pileup" | "all" clusters
    snippet_pre: int = 20             # snippet: samples kept before the cluster ...
    snippet_post: int = 80            # ... and after it
    large_v: float = 1.5              # above this: FLAG_LARGE (cosmic muons; shape not verified)
    large_shadow: int = 200           # events this many samples after a large pulse: FLAG_AFTER_LARGE
    bad_fit_sigma: float = 4.0        # bad fit if rms > bad_fit_sigma * noise + bad_fit_frac * amplitude
    bad_fit_frac: float = 0.02


@dataclass
class PipelineConfig:
    preprocess: PreprocessConfig = field(default_factory=PreprocessConfig)
    detect: DetectConfig = field(default_factory=DetectConfig)
    events: EventConfig = field(default_factory=EventConfig)
    solver: SolverConfig = field(default_factory=SolverConfig)

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
        self.solver_on = self.cfg.solver.enabled and not self.fit_baseline
        if self.cfg.solver.enabled and self.fit_baseline:
            warnings.warn("pile-up solver needs a restored baseline (pole_zero + baseline_tracker); "
                          "it is disabled for this configuration", stacklevel=2)
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

        e["flags"] = flags

        # snippets of the detected clusters (raw material for re-analysis)
        snippets = []
        if ev.snippets != "none" and ne:
            want = sizes > 1 if ev.snippets == "pileup" else np.ones(sizes.size, bool)
            first = np.r_[0, np.cumsum(sizes)[:-1]]
            for k in np.flatnonzero(want):
                a = int(idx[first[k]]) - ev.snippet_pre
                b = int(idx[first[k] + sizes[k] - 1]) + ev.snippet_post
                a, b = max(a, 0), min(b, n)
                snippets.append((int(e["cluster"][first[k]]), capture, a, y[a:b]))

        # pile-up solver: clusters and poorly fitted isolated pulses
        if self.solver_on and ne:
            e = self._solve_clusters(e, y, det.index)

        # very large pulses and their shadow
        idx = e["index"]
        flags = e["flags"]
        big = np.maximum(np.nan_to_num(e["amplitude"]), np.nan_to_num(e["filter_amplitude"])) > ev.large_v
        flags[big] |= FLAG_LARGE
        for i in np.flatnonzero(big):
            shadow = (idx > idx[i]) & (idx <= idx[i] + ev.large_shadow)
            flags[shadow] |= FLAG_AFTER_LARGE
        e["flags"] = flags
        e["time_s"] = t0_s + (idx + e["shift"]) * dt
        ne = e.size

        row = np.zeros(1, CAPTURE_DTYPE)
        info = self.pre.info
        row[0] = (capture, t0_s, (n - c.preprocess.warmup_samples - ev.fit_post - 1) * dt, info.offset_v,
                  info.sigma_in_v, info.sigma_out_v, det.sigma, det.threshold_v, ne, info.line_freqs_hz.size,
                  t0_s + c.preprocess.warmup_samples * dt)
        return e, row, snippets

    def _solve_clusters(self, e: np.ndarray, y: np.ndarray, detected: np.ndarray) -> np.ndarray:
        """Run the pile-up solver on every cluster of more than one pulse and on every poorly
        fitted isolated pulse; insert pulses found by its residual search; recompute gaps
        and cluster sizes."""
        ev, sc = self.cfg.events, self.cfg.solver
        n = y.size
        free = pulse_free_mask(n, detected, 20, 120)
        solver = PileupSolver(self.template, noise_autocovariance(y, free, sc.acov_lags), sc)
        noise = self.pre.info.sigma_out_v
        flags = e["flags"]
        starts = np.flatnonzero(np.r_[True, np.diff(e["cluster"]) != 0])
        ends = np.r_[starts[1:], e.size]
        new_rows = []
        for a0, b0 in zip(starts, ends):
            rows = np.arange(a0, b0)
            f = flags[rows]
            if np.any(f & FLAG_EDGE):
                continue
            if rows.size == 1 and not (f[0] & FLAG_BAD_FIT):
                continue
            idx = e["index"][rows]
            a = max(int(idx[0]) - ev.fit_pre, 0)
            b = min(int(idx[-1]) + ev.fit_post + 1, n)
            res = solver.solve(y[a:b], idx - a)
            k = rows.size
            e["amplitude"][rows] = res.amplitude[:k]
            e["shift"][rows] = res.time[:k] + a - idx
            e["fit_rms"][rows] = res.rms
            amax = float(np.max(np.abs(res.amplitude)))
            bad = res.rms > ev.bad_fit_sigma * noise + ev.bad_fit_frac * amax
            not_bad = np.uint16(~np.uint16(FLAG_BAD_FIT))
            flags[rows] = (flags[rows] | FLAG_BAD_FIT) if bad else (flags[rows] & not_bad)
            if res.added.any():
                flags[rows] |= FLAG_PILEUP
                for j in np.flatnonzero(res.added):
                    r = e[rows[0]].copy()
                    t_abs = res.time[j] + a
                    r["index"] = int(np.round(t_abs))
                    r["shift"] = t_abs - r["index"]
                    r["amplitude"] = res.amplitude[j]
                    r["filter_amplitude"] = np.nan
                    r["flags"] = (flags[rows[0]] & (FLAG_WARMUP | FLAG_BAD_FIT)) | FLAG_PILEUP | FLAG_FIT_ADDED
                    new_rows.append(r)
        e["flags"] = flags
        if new_rows:
            e = np.concatenate([e, np.array(new_rows, dtype=EVENT_DTYPE)])
            e = e[np.lexsort((e["shift"], e["index"]))]
        # gaps and cluster sizes after insertion
        gaps = np.diff(e["index"])
        e["gap_prev"] = np.r_[-1, gaps]
        e["gap_next"] = np.r_[gaps, -1]
        _, inv, counts = np.unique(e["cluster"], return_inverse=True, return_counts=True)
        e["cluster_size"] = counts[inv]
        return e


_WORKER: dict = {}
_THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


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
            # One BLAS/OpenMP thread per worker: the workers already use every core, and
            # multi-threaded linear algebra inside each of them would oversubscribe the CPU.
            # Worker processes read these variables when they start (spawn).
            saved = {k: os.environ.get(k) for k in _THREAD_VARS}
            os.environ.update({k: "1" for k in _THREAD_VARS})
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
                for k, v in saved.items():          # restore the caller's thread settings
                    if v is None:
                        os.environ.pop(k, None)
                    else:
                        os.environ[k] = v
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


def spectrum_report(event_path, n_peaks: int = 2, guess=None, min_counts: int = 1000) -> dict:
    """Peak fits of the event spectrum by category: isolated pulses (template fit), pulses
    resolved by the pile-up solver (by distance to their neighbour), and all good events.
    Also the throughput: the fraction of valid events that have a usable amplitude."""
    from .events import FLAG_EDGE, FLAG_FIT_ADDED, FLAG_PILEUP, FLAG_WARMUP, EventFile
    from .restoration import _auto_peaks
    from .spectrum import fit_peaks

    ef = EventFile(event_path)
    ev = ef.events
    valid = (ev["flags"] & (FLAG_WARMUP | FLAG_EDGE)) == 0
    good = ef.good(include_pileup=True)
    iso = ef.clean()
    pile = good[(good["flags"] & FLAG_PILEUP) != 0]
    gap = np.minimum(np.where(pile["gap_prev"] > 0, pile["gap_prev"], 10**6),
                     np.where(pile["gap_next"] > 0, pile["gap_next"], 10**6))
    cats = {
        "isolated": iso["amplitude"],
        "pile-up (all)": pile["amplitude"],
        "pile-up, gap < 3 (incl. found by the solver)": pile["amplitude"][gap < 3],
        "pile-up, gap 3-10": pile["amplitude"][(gap >= 3) & (gap < 10)],
        "pile-up, gap 10-30": pile["amplitude"][(gap >= 10) & (gap < 30)],
        "pile-up, gap 30-60": pile["amplitude"][(gap >= 30) & (gap < 60)],
        "all good events": good["amplitude"],
    }
    g = np.asarray(guess) if guess is not None else _auto_peaks(iso["amplitude"], n_peaks)
    out = dict(valid_events=int(valid.sum()), good_events=int(good.size),
               throughput=float(good.size / max(valid.sum(), 1)),
               isolated_fraction=float(iso.size / max(valid.sum(), 1)),
               fit_added=int(np.sum((ev["flags"] & FLAG_FIT_ADDED) != 0)), categories={})
    ref = None
    lo, hi = g.min() * 0.85, g.max() * 1.12
    for name, a in cats.items():
        in_range = int(np.sum((a >= lo) & (a < hi)))
        if in_range < min_counts:
            out["categories"][name] = dict(n=int(a.size), error=f"too few events in the peak region ({in_range})")
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pf = fit_peaks(a, g, lo, hi, 0.002)
            row = dict(n=int(a.size), centers=pf.centers.tolist(), fwhm_pct=(100 * pf.resolution).tolist(),
                       fwhm_err_pct=(100 * pf.resolution_err).tolist())
            if ref is None:
                ref = pf.centers
            row["shift_pct"] = (100 * (pf.centers / ref - 1)).tolist()
        except Exception as exc:  # noqa: BLE001
            row = dict(n=int(a.size), error=str(exc))
        out["categories"][name] = row
    return out


def plot_event_spectrum(event_path, png_path) -> None:
    from . import plotting as P
    from .events import FLAG_PILEUP, EventFile
    ef = EventFile(event_path)
    ev = ef.events
    clean = ef.clean()["amplitude"]
    good = ef.good(include_pileup=True)
    resolved = good[(good["flags"] & FLAG_PILEUP) != 0]["amplitude"]
    hi = float(np.quantile(clean, 0.995) * 1.3) if clean.size else 1.0
    edges = np.linspace(0, hi, 700)
    fig, ax = P.figure(1, 1, height_per_row=4.0)
    a = ax[0, 0]
    if resolved.size:
        a.stairs(np.histogram(good["amplitude"], edges)[0], edges, color=P.SERIES[2], linewidth=1.1,
                 label=f"all good events ({good.size})")
    a.stairs(np.histogram(clean, edges)[0], edges, color=P.SERIES[0], linewidth=1.1,
             label=f"isolated pulses, template fit ({clean.size})")
    if resolved.size:
        a.stairs(np.histogram(resolved, edges)[0], edges, color=P.SERIES[1], linewidth=1.0,
                 label=f"pile-up pulses, resolved by the solver ({resolved.size})")
        title = "pile-up resolved"
    else:
        pile = ev[(ev["flags"] & FLAG_PILEUP) != 0]["filter_amplitude"]
        a.stairs(np.histogram(pile, edges)[0], edges, color=P.SERIES[1], linewidth=1.0,
                 label=f"pile-up pulses, detection estimate only ({pile.size})")
        title = "pile-up not resolved"
    a.set_yscale("log")
    a.set_xlabel("pulse amplitude (V)")
    a.set_ylabel("counts per bin")
    a.set_title(f"Event spectrum ({ef.live_time_s:.1f} s live time; {title})")
    a.legend()
    P.save(fig, png_path)
