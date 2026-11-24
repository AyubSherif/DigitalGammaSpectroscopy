"""Signal sources for the stream processor: where the samples come from (step 10).

Every source yields ``Chunk``s of samples in volts. A chunk with ``new_segment = True``
starts a new, separate recording (a capture). Samples within a segment are continuous.

| Source | Segments | Use |
|---|---|---|
| ``StoreSource`` | one per capture of an HDF5 store | replay converted data as a stream |
| ``RawSource`` | one per ASCII capture file (no conversion) | stream straight from the digitiser's files |
| ``SimulationSource`` | one continuous segment | test input, not limited to the recorded data |
| ``PipeSource`` | one continuous segment | binary samples from stdin, a file or a FIFO (a live digitiser, or ``gamma-spectrometry simulate --stdout``) |

``Paced`` releases a source's chunks no faster than real time (×``speed``), which turns a
replay or a simulation into a live stream. A live source cannot wait for the processor:
if processing falls behind, frames are dropped and counted as dead time.
"""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np


@dataclass
class Chunk:
    data: np.ndarray            # samples, volts (as recorded: polarity and offset untouched)
    new_segment: bool = False   # first chunk of a new recording (discontinuity before it)
    t0_s: float | None = None   # time of the first sample of the segment (with new_segment)


class Source:
    fs: float = 1.8e6
    live: bool = False          # True: cannot be paused (drop frames when processing falls behind)
    description: str = ""
    lsb_volts: float | None = None

    def chunks(self) -> Iterator[Chunk]:                    # pragma: no cover - interface
        raise NotImplementedError


def _split(x: np.ndarray, chunk: int, t0: float) -> Iterator[Chunk]:
    for k, a in enumerate(range(0, x.size, chunk)):
        yield Chunk(x[a:a + chunk], new_segment=k == 0, t0_s=t0 if k == 0 else None)


class StoreSource(Source):
    def __init__(self, path, limit: int | None = None, chunk: int = 65536, first: int = 0):
        from .store import CaptureStore

        self.path, self.limit, self.chunk, self.first = Path(path), limit, chunk, first
        with CaptureStore(self.path) as s:
            self.fs, self.lsb_volts, self.n = 1.0 / s.dt_s, s.lsb_volts, len(s)
        self.description = f"store {self.path}"

    def chunks(self) -> Iterator[Chunk]:
        from .store import CaptureStore

        with CaptureStore(self.path) as s:
            last = len(s) if self.limit is None else min(len(s), self.first + self.limit)
            for i in range(self.first, last):
                yield from _split(s.volts(i), self.chunk, float(s.time_s[i]))


class RawSource(Source):
    def __init__(self, folder, pattern: str = "raw_data_*", limit: int | None = None, chunk: int = 65536):
        from .rawio import build_index

        self.headers = build_index(folder, pattern)
        if limit:
            self.headers = self.headers[:limit]
        self.fs = 1.0 / self.headers[0].dt_s
        self.chunk = chunk
        self.description = f"raw files {folder} ({len(self.headers)} captures)"

    def chunks(self) -> Iterator[Chunk]:
        from .rawio import read_raw_file

        for h in self.headers:
            _, v = read_raw_file(h.path, lsb=None)              # volts, no integer conversion
            yield from _split(v, self.chunk, h.time_s)


class SimulationSource(Source):
    def __init__(self, cfg=None, duration_s: float = 10.0, chunk: int = 65536, t0_s: float = 0.0):
        from .simulate import Simulator

        self.sim = Simulator(cfg, t0_s=t0_s)
        self.fs, self.duration_s, self.chunk = self.sim.fs, duration_s, chunk
        self.lsb_volts = self.sim.cfg.lsb_v or None
        self.description = (f"simulation ({self.sim.cfg.nuclide}, {self.sim.cfg.rate_cps:g} + "
                            f"{self.sim.cfg.background_cps:g} /s, seed {self.sim.cfg.seed})")

    def chunks(self) -> Iterator[Chunk]:
        n_total = int(round(self.duration_s * self.fs)) if self.duration_s > 0 else None
        first = True
        while n_total is None or self.sim.pos < n_total:
            n = self.chunk if n_total is None else min(self.chunk, n_total - self.sim.pos)
            yield Chunk(self.sim.read(n), new_segment=first, t0_s=self.sim.t0_s if first else None)
            first = False


PIPE_DTYPES = {"f4": "<f4", "f8": "<f8", "i2": "<i2", "i4": "<i4"}


class PipeSource(Source):
    """Little-endian binary samples (``dtype`` f4/f8/i2/i4) from a file, a FIFO or stdin ("-").

    Integer samples are multiplied by ``scale`` (the ADC step, volts per code) to get volts.
    """

    live = True

    def __init__(self, path: str = "-", fs: float = 1.8e6, dtype: str = "f4", scale: float = 1.0,
                 chunk: int = 65536, t0_s: float | None = None):
        if dtype not in PIPE_DTYPES:
            raise ValueError(f"pipe dtype {dtype!r}: choose from {', '.join(PIPE_DTYPES)}")
        self.path, self.fs, self.dtype, self.scale, self.chunk = path, fs, np.dtype(PIPE_DTYPES[dtype]), scale, chunk
        self.t0_s = time.time() if t0_s is None else t0_s
        self.lsb_volts = scale if dtype in ("i2", "i4") else None
        self.description = f"pipe {'stdin' if path == '-' else path} ({dtype}, {fs:g} Hz)"

    def chunks(self) -> Iterator[Chunk]:
        f = sys.stdin.buffer if self.path == "-" else open(self.path, "rb")
        try:
            first, rest = True, b""
            nbytes = self.chunk * self.dtype.itemsize
            while True:
                b = f.read(nbytes - len(rest))
                if not b:
                    break
                b = rest + b
                usable = len(b) - len(b) % self.dtype.itemsize
                rest = b[usable:]
                if usable == 0:
                    continue
                x = np.frombuffer(b[:usable], self.dtype).astype(np.float64) * self.scale
                yield Chunk(x, new_segment=first, t0_s=self.t0_s if first else None)
                first = False
        finally:
            if f is not sys.stdin.buffer:
                f.close()


class Paced(Source):
    """Release the chunks of ``source`` at ``speed`` × real time (a live stream)."""

    live = True

    def __init__(self, source: Source, speed: float = 1.0):
        self.source, self.speed = source, speed
        self.fs, self.lsb_volts = source.fs, source.lsb_volts
        self.description = f"{source.description}, paced at {speed:g}x real time"

    def chunks(self) -> Iterator[Chunk]:
        t_start, n = time.perf_counter(), 0
        for c in self.source.chunks():
            n += c.data.size
            wait = n / self.fs / self.speed - (time.perf_counter() - t_start)
            if wait > 0:
                time.sleep(wait)
            yield c


def open_source(kind: str, path: str = "", *, limit: int | None = None, pattern: str = "raw_data_*",
                chunk: int = 65536, sim_cfg=None, duration_s: float = 10.0, fs: float = 1.8e6,
                pipe_dtype: str = "f4", pipe_scale: float = 1.0, realtime: bool = False,
                speed: float = 1.0) -> Source:
    """A source by name: "store" | "raw" | "replay" (store or raw folder, by path) | "simulation" | "pipe"."""
    if kind == "replay":
        kind = "raw" if Path(path).is_dir() else "store"
    if kind == "store":
        s: Source = StoreSource(path, limit, chunk)
    elif kind == "raw":
        s = RawSource(path, pattern, limit, chunk)
    elif kind == "simulation":
        s = SimulationSource(sim_cfg, duration_s, chunk)
    elif kind == "pipe":
        s = PipeSource(path or "-", fs, pipe_dtype, pipe_scale, chunk)
    else:
        raise ValueError(f"unknown source kind {kind!r}")
    return Paced(s, speed) if realtime and not s.live else s
