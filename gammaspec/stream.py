"""Streaming: continuous input processed in overlapping frames, with a live spectrum (step 10).

A stream is cut into **frames**. Each frame is a *core* of ``frame`` samples plus ``pre``
samples of context before it and ``post`` samples after it::

    segment:  ...|<- pre ->|<-------- core k -------->|<- post ->|...
                            |<- pre ->|<-------- core k+1 -------->|<- post ->|

The unchanged per-capture pipeline (``Pipeline.process_capture``) processes the whole frame.
Only events whose peak falls in the core are kept. Why this works:
* The pre-processing restarts at the frame start, and its warm-up (``warmup_samples``) is
  inside the pre-context.
* The fit windows, the baseline tracker's centred average and pile-up clusters reach into
  the post-context.
* So every kept event is processed exactly as one in the middle of a long recording.

Frames are therefore independent and run in parallel worker processes. The stream state is
just the buffered context. Results are identical whatever chunk sizes the source delivers.
A chunk marked ``new_segment`` (a new capture) closes the previous segment: its last frame
has no post-context, as at the end of a capture. With ``frame`` ≥ the capture length, each
capture is one frame, and the result equals batch processing (``pipeline.run``).

``ADR-0012`` records the choice.
"""

from __future__ import annotations

import os
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from .events import EventWriter
from .pipeline import _THREAD_VARS, Pipeline, PipelineConfig
from .sources import Source
from .template import PulseTemplate


@dataclass
class StreamConfig:
    frame: int = 262144                   # samples per frame core (0.15 s at 1.8 MS/s)
    pre: int = 4096                       # context before each core (> warmup_samples)
    post: int = 2048                      # context after each core (fit windows, tracker, clusters)
    chunk: int = 65536                    # samples per read from the source
    max_backlog: int = 16                 # frames waiting for a worker; a live source drops frames beyond this
    realtime: bool = False                # replay / simulation: release samples at real time (a live stream)
    speed: float = 1.0                    # ... at this multiple of real time
    update_s: float = 2.0                 # status line and live spectrum (JSON) every this many seconds (0 = never)
    png_s: float = 30.0                   # live spectrum picture every this many seconds (0 = only at the end)
    channels: int = 4096                  # live spectrum: channels ...
    v_max: float = 2.0                    # ... spanning 0 to this amplitude (V)


def check_stream_config(scfg: StreamConfig, cfg: PipelineConfig) -> list[str]:
    """Problems with the frame geometry (errors)."""
    p, ev = cfg.preprocess, cfg.events
    need_pre = p.warmup_samples + ev.fit_pre + 2
    need_post = max(ev.fit_post + 1, p.tracker.window // 2 + 1, ev.snippet_post) + ev.cluster_gap
    err = []
    if scfg.pre < need_pre:
        err.append(f"stream.pre = {scfg.pre}: must be >= warmup_samples + fit_pre + 2 = {need_pre}")
    if scfg.post < need_post:
        err.append(f"stream.post = {scfg.post}: must be >= {need_post} (fit window, baseline tracker, cluster gap)")
    if scfg.frame < scfg.pre + scfg.post:
        err.append(f"stream.frame = {scfg.frame}: must be >= pre + post = {scfg.pre + scfg.post}")
    return err


# ----------------------------------------------------------------------------------------------
# One frame
# ----------------------------------------------------------------------------------------------

def process_frame(pipe: Pipeline, frame: np.ndarray, frame_id: int, t_core: float, pre_used: int,
                  core_len: int, post_used: int):
    """Process a frame and keep its core: (events, capture_row, snippets), indices relative
    to the core start."""
    dt = 1.0 / pipe.fs
    c = pipe.cfg
    pipe._cluster0 = 0
    e, row, sn = pipe.process_capture(frame, frame_id, t_core - pre_used * dt)
    keep = (e["index"] >= pre_used) & (e["index"] < pre_used + core_len)
    e = e[keep]
    e["index"] -= pre_used
    kept = set(e["cluster"].tolist())
    sn = [(cl, cap, a - pre_used, w) for cl, cap, a, w in sn if cl in kept]
    warm = max(0, c.preprocess.warmup_samples - pre_used)        # only in a segment's first frame
    tail = max(0, c.events.fit_post + 1 - post_used)             # only in a segment's last frame
    row["time_s"] = t_core
    row["live_start_s"] = t_core + warm * dt
    row["live_time_s"] = max(core_len - warm - tail, 0) * dt
    row["n_events"] = e.size
    return e, row, sn


_W: dict = {}


def _worker_init(template: PulseTemplate, fs: float, cfg: PipelineConfig):
    import warnings

    warnings.simplefilter("ignore")                              # config warnings were shown already
    _W["pipe"] = Pipeline(template, fs, cfg)


def _worker_frame(args):
    return process_frame(_W["pipe"], _decode(*args[0]), *args[1:])


def _encode(frame: np.ndarray, lsb: float | None):
    """Compact form for sending a frame to a worker: exact int32 ADC codes when the samples
    lie on the ADC grid (stores, raw files, quantised simulation), else float32."""
    if lsb:
        codes = np.rint(frame / lsb)
        if np.max(np.abs(codes * lsb - frame), initial=0.0) <= 0.02 * lsb and np.max(np.abs(codes), initial=0) < 2**31:
            return codes.astype(np.int32), lsb
    return frame.astype(np.float32), None


def _decode(data: np.ndarray, lsb: float | None) -> np.ndarray:
    return data.astype(np.float64) * lsb if lsb else data.astype(np.float64)


# ----------------------------------------------------------------------------------------------
# Live spectrum
# ----------------------------------------------------------------------------------------------

class LiveSpectrum:
    """Histogram of good events, updated as frames complete; calibrated once enough counts exist."""

    def __init__(self, channels: int = 4096, v_max: float = 2.0, kind: str = "source", calibration=None,
                 min_counts: int = 20000, keep_max: int = 20_000_000):
        self.edges_v = np.linspace(0.0, v_max, channels + 1)
        self.counts = np.zeros(channels, np.int64)
        self.kind, self.calibration, self.fixed = kind, calibration, calibration is not None
        self.min_counts, self.keep_max = min_counts, keep_max
        self._amps: list[np.ndarray] = []
        self._n_kept = 0
        self.live_time_s = 0.0
        self.n_events = 0
        self._last_cal_n = 0

    def add(self, events: np.ndarray, live_time_s: float) -> None:
        from .simulate import _good_mask

        a = events["amplitude"][_good_mask(events)].astype(np.float32)
        self.counts += np.histogram(a, self.edges_v)[0]
        self.live_time_s += live_time_s
        self.n_events += events.size
        if self._n_kept < self.keep_max:
            self._amps.append(a)
            self._n_kept += a.size

    def try_calibrate(self) -> bool:
        """(Re)calibrate from the amplitudes so far; every doubling of the counts."""
        if self.fixed or self._n_kept < self.min_counts or self._n_kept < 2 * self._last_cal_n:
            return self.calibration is not None
        from .calibration import calibrate

        try:
            self.calibration = calibrate(np.concatenate(self._amps), self.kind)
            self._last_cal_n = self._n_kept
        except (ValueError, RuntimeError):
            pass
        return self.calibration is not None

    def snapshot(self) -> dict:
        d = dict(live_time_s=self.live_time_s, events=self.n_events, counts=self.counts.tolist(),
                 edges_v=self.edges_v.tolist(), good=int(self.counts.sum()))
        if self.calibration is not None:
            d["calibration"] = list(map(float, self.calibration.coeffs))
            d["edges_kev"] = self.calibration.energy(self.edges_v).tolist()
        return d

    def save(self, out_dir, png: bool = True) -> None:
        import json

        from . import plotting as P

        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        tmp = out / "live_spectrum.json.tmp"                      # replace atomically (a UI may be reading)
        tmp.write_text(json.dumps(self.snapshot()), encoding="utf-8")
        os.replace(tmp, out / "live_spectrum.json")
        if not png:
            return
        fig, ax = P.figure(1, 1, height_per_row=3.2)
        x = self.calibration.energy(self.edges_v) if self.calibration is not None else self.edges_v * 1e3
        rate = self.counts / max(self.live_time_s, 1e-9) / np.diff(x)
        ax[0, 0].stairs(np.maximum(rate, 1e-6), x, color=P.SERIES[0], linewidth=1.0)
        ax[0, 0].set_yscale("log")
        ax[0, 0].set_xlabel("energy (keV)" if self.calibration is not None else "amplitude (mV)")
        ax[0, 0].set_ylabel("counts / live s / " + ("keV" if self.calibration is not None else "mV"))
        ax[0, 0].set_title(f"live spectrum: {self.counts.sum()} good events, {self.live_time_s:.1f} s live")
        P.save(fig, out / "live_spectrum.png")


# ----------------------------------------------------------------------------------------------
# The stream processor
# ----------------------------------------------------------------------------------------------

@dataclass
class StreamStatus:
    elapsed_s: float = 0.0          # wall clock since start
    samples: int = 0                # samples read from the source
    stream_s: float = 0.0           # signal time read (samples / fs)
    frames: int = 0                 # frames processed
    dropped_frames: int = 0         # frames skipped because processing fell behind (live sources)
    backlog: int = 0                # frames waiting for a worker
    events: int = 0
    live_time_s: float = 0.0

    @property
    def speed(self) -> float:       # processed signal time per wall-clock second
        return self.stream_s / self.elapsed_s if self.elapsed_s > 0 else 0.0


class StreamProcessor:
    """Feed chunks; frames are cut, processed (in parallel) and written in order."""

    def __init__(self, template: PulseTemplate, fs: float, cfg: PipelineConfig, scfg: StreamConfig,
                 writer: EventWriter | None = None, workers: int | None = None, live: bool = False,
                 spectrum: LiveSpectrum | None = None, lsb_volts: float | None = None,
                 flush_events: int = 200_000):
        errs = check_stream_config(scfg, cfg)
        if errs:
            raise ValueError("; ".join(errs))
        self.fs, self.cfg, self.scfg, self.writer, self.live = fs, cfg, scfg, writer, live
        self.spectrum, self.lsb = spectrum, lsb_volts
        self.status = StreamStatus()
        self._out: list = []                 # results waiting to be written (batched: HDF5 appends are slow)
        self._out_events = 0
        self._flush_events = flush_events
        self._frame_id = 0
        self._cluster0 = 0
        self._pending: deque = deque()
        self._buf = np.zeros(0)
        self._buf_start = 0          # segment sample index of _buf[0]
        self._core = 0               # segment sample index of the next core
        self._seg_t0 = 0.0
        self._in_segment = False
        self.results: list | None = [] if writer is None else None     # kept in memory without a writer
        self.workers = workers or max(1, min(8, (os.cpu_count() or 2) - 1))
        self._pool = None
        self._saved_env: dict = {}
        if self.workers > 1:
            self._saved_env = {k: os.environ.get(k) for k in _THREAD_VARS}
            os.environ.update({k: "1" for k in _THREAD_VARS})
            self._pool = ProcessPoolExecutor(self.workers, initializer=_worker_init, initargs=(template, fs, cfg))
        else:
            self._pipe = Pipeline(template, fs, cfg)
        self._t_start = time.perf_counter()

    # -- input ---------------------------------------------------------------------------------
    def feed(self, data: np.ndarray, new_segment: bool = False, t0_s: float | None = None) -> None:
        if new_segment or not self._in_segment:
            self.end_segment()
            self._in_segment = True
            self._buf, self._buf_start, self._core = np.zeros(0), 0, 0
            self._seg_t0 = 0.0 if t0_s is None else float(t0_s)
        self._buf = np.concatenate([self._buf, np.asarray(data, dtype=np.float64)])
        self.status.samples += data.size
        self.status.stream_s = self.status.samples / self.fs
        s = self.scfg
        while self._buf_start + self._buf.size >= self._core + s.frame + s.post:
            self._emit(s.frame, s.post)
        self._collect(block=False)

    def end_segment(self) -> None:
        """Process what remains of the current segment (its last frame has no post-context)."""
        if not self._in_segment:
            return
        end = self._buf_start + self._buf.size
        if end > self._core:
            self._emit(end - self._core, 0)
        self._in_segment = False

    def _emit(self, core_len: int, post_used: int) -> None:
        s = self.scfg
        a = max(self._core - s.pre, 0)
        b = self._core + core_len + post_used
        frame = self._buf[a - self._buf_start: b - self._buf_start]
        args = (_encode(frame, self.lsb), self._frame_id, self._seg_t0 + self._core / self.fs, self._core - a,
                core_len, post_used)
        self._frame_id += 1
        self._core += core_len
        keep_from = max(self._core - s.pre, 0)                    # context for the next frame
        self._buf = self._buf[keep_from - self._buf_start:]
        self._buf_start = keep_from
        if self._pool is None:
            self._pending.append(process_frame(self._pipe, _decode(*args[0]), *args[1:]))
        else:
            if len(self._pending) >= s.max_backlog:
                if self.live:                                     # cannot wait: drop (dead time)
                    self.status.dropped_frames += 1
                    return
                self._collect(block=True, until=s.max_backlog - 1)
            self._pending.append(self._pool.submit(_worker_frame, args))
        self.status.backlog = len(self._pending)

    # -- output --------------------------------------------------------------------------------
    def _collect(self, block: bool, until: int = 0) -> None:
        while self._pending:
            head = self._pending[0]
            if hasattr(head, "result"):
                if not block and not head.done():
                    break
                if block and len(self._pending) <= until:
                    break
                res = head.result()
            else:
                res = head
            self._pending.popleft()
            self._write(*res)
        self.status.backlog = len(self._pending)

    def _write(self, e: np.ndarray, row: np.ndarray, sn: list) -> None:
        e["cluster"] += self._cluster0
        sn = [(cl + self._cluster0, cap, a, w) for cl, cap, a, w in sn]
        if e.size:
            self._cluster0 = int(e["cluster"].max()) + 1
        if self.writer is not None:
            self._out.append((e, row, sn))
            self._out_events += e.size
            if self._out_events >= self._flush_events or len(self._out) >= 64:
                self.flush()
        else:
            self.results.append((e, row, sn))
        st = self.status
        st.frames += 1
        st.events += e.size
        live = float(row["live_time_s"][0])
        st.live_time_s += live
        if self.spectrum is not None:
            self.spectrum.add(e, live)

    def flush(self) -> None:
        """Write the buffered results to the event file."""
        if self.writer is None or not self._out:
            return
        e = np.concatenate([r[0] for r in self._out])
        rows = np.concatenate([r[1] for r in self._out])
        sn = [s for r in self._out for s in r[2]]
        self.writer.append(e, rows, sn)
        self._out, self._out_events = [], 0

    def close(self) -> StreamStatus:
        self.end_segment()
        try:
            self._collect(block=True)
            self.flush()
        finally:
            if self._pool is not None:
                self._pool.shutdown(cancel_futures=True)
                for k, v in self._saved_env.items():
                    if v is None:
                        os.environ.pop(k, None)
                    else:
                        os.environ[k] = v
        self.status.elapsed_s = time.perf_counter() - self._t_start
        return self.status


def process_stream(source: Source, template_path, out_path, cfg: PipelineConfig | None = None,
                   scfg: StreamConfig | None = None, workers: int | None = None, verbose: bool = True,
                   live_dir=None, on_update: Callable[[StreamStatus, LiveSpectrum], None] | None = None,
                   kind: str = "source", calibration=None, overwrite: bool = False) -> dict:
    """Process a whole source into an event file, with periodic status and live spectrum."""
    cfg = cfg or PipelineConfig()
    scfg = scfg or StreamConfig()
    template = PulseTemplate.load(template_path)
    probe = Pipeline(template, source.fs, cfg)                    # shows config warnings once, here
    attrs = dict(config=cfg.to_dict(), template=str(template_path), source_store=source.description,
                 lsb_volts=source.lsb_volts, dt_s=1.0 / source.fs, fit_baseline=probe.fit_baseline,
                 pole_zero=None if probe.pz is None else dict(b=probe.pz.b.tolist(), a=probe.pz.a.tolist()),
                 stream=scfg.__dict__)
    spec = LiveSpectrum(scfg.channels, scfg.v_max, kind, calibration)
    t_last = time.perf_counter()
    with EventWriter(out_path, attrs, overwrite=overwrite) as w:
        sp = StreamProcessor(template, source.fs, cfg, scfg, w, workers, source.live, spec, source.lsb_volts)
        t_png = time.perf_counter()
        try:
            for ch in source.chunks():
                sp.feed(ch.data, ch.new_segment, ch.t0_s)
                if scfg.update_s and time.perf_counter() - t_last >= scfg.update_s:
                    t_last = time.perf_counter()
                    st = sp.status
                    st.elapsed_s = t_last - sp._t_start
                    spec.try_calibrate()
                    if verbose:
                        print(_status_line(st, spec), flush=True)
                    if live_dir is not None:
                        png = bool(scfg.png_s) and t_last - t_png >= scfg.png_s
                        spec.save(live_dir, png=png)
                        t_png = t_last if png else t_png
                    if on_update is not None:
                        on_update(st, spec)
        finally:
            st = sp.close()
    spec.try_calibrate()
    if live_dir is not None:
        spec.save(live_dir)
    if on_update is not None:
        on_update(st, spec)
    if verbose:
        print(_status_line(st, spec) + "  (done)", flush=True)
    return dict(frames=st.frames, dropped_frames=st.dropped_frames, samples=st.samples, stream_s=st.stream_s,
                events=st.events, live_time_s=st.live_time_s, elapsed_s=st.elapsed_s, speed=st.speed,
                calibration=None if spec.calibration is None else list(map(float, spec.calibration.coeffs)))


def _status_line(st: StreamStatus, spec: LiveSpectrum) -> str:
    cal = ("E = %.1f + %.1f A keV" % tuple(spec.calibration.coeffs[:2])) if spec.calibration is not None else "uncalibrated"
    rate = st.events / st.live_time_s if st.live_time_s else 0.0
    return (f"[stream] {st.stream_s:7.1f} s signal in {st.elapsed_s:6.1f} s ({st.speed:4.2f}x real time) | "
            f"{st.frames} frames, {st.dropped_frames} dropped, backlog {st.backlog} | {st.events} events, "
            f"{rate:.0f} /s | {cal}")
