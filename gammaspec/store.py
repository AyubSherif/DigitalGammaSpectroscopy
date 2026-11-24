"""Compressed binary store of ADC captures (HDF5).

Layout of a store file::

    /codes            int32 (n_captures, n_samples)   ADC codes, one row per capture, chronological
    /index/time_s     float64  wall-clock creation time (naive local time encoded as if UTC)
    /index/created    str      ISO creation timestamp from the file header
    /index/file_number int32   number in the original file name
    /index/source     str      original file path
    attrs: lsb_volts, dt_s, n_samples, gammaspec_version, ...

Rows are written in chronological order.  Codes are lossless; volts = codes * lsb_volts.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Callable, Iterable, Iterator

import h5py
import numpy as np

from . import __version__
from .rawio import DEFAULT_LSB_VOLTS, RawHeader, estimate_lsb, read_raw_file

_STR = h5py.string_dtype("utf-8")


def _read_codes(args: tuple[str, float]) -> np.ndarray:
    path, lsb = args
    return read_raw_file(path, lsb=lsb)[1]


def convert(
    headers: list[RawHeader],
    out_path: str | Path,
    lsb: float | None = DEFAULT_LSB_VOLTS,
    workers: int | None = None,
    overwrite: bool = False,
    progress: Callable[[int, int], None] | None = None,
) -> Path:
    """Convert capture files (already in the desired order) to one compressed HDF5 store.

    lsb=None estimates the ADC step from the first capture.
    """
    if not headers:
        raise ValueError("no captures to convert")
    out_path = Path(out_path)
    if out_path.exists() and not overwrite:
        raise FileExistsError(f"{out_path} exists (pass overwrite=True to replace it)")
    n_samples = headers[0].n_samples
    dt_s = headers[0].dt_s
    for h in headers:
        if h.n_samples != n_samples or not np.isclose(h.dt_s, dt_s, rtol=1e-12):
            raise ValueError(f"{h.path}: capture length/sampling differs from {headers[0].path}")
    if lsb is None:
        lsb = estimate_lsb(read_raw_file(headers[0].path, lsb=None)[1])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_name(out_path.name + ".partial")
    n = len(headers)
    workers = workers or max(1, min(8, (os.cpu_count() or 2) - 1))
    with h5py.File(tmp_path, "w") as f:
        codes = f.create_dataset(
            "codes", shape=(n, n_samples), dtype=np.int32, chunks=(1, n_samples),
            compression="gzip", compression_opts=1, shuffle=True,
        )
        idx = f.create_group("index")
        idx.create_dataset("time_s", data=np.array([h.time_s for h in headers]))
        idx.create_dataset("created", data=[h.created.isoformat() for h in headers], dtype=_STR)
        idx.create_dataset("file_number", data=np.array([h.file_number for h in headers], dtype=np.int32))
        idx.create_dataset("source", data=[str(h.path) for h in headers], dtype=_STR)
        f.attrs.update(
            lsb_volts=lsb, dt_s=dt_s, n_samples=n_samples, sample_rate_hz=1.0 / dt_s,
            gammaspec_version=__version__, source_dir=str(headers[0].path.parent),
        )
        jobs = [(str(h.path), lsb) for h in headers]
        if workers == 1:
            results: Iterable[np.ndarray] = map(_read_codes, jobs)
            pool = None
        else:
            pool = ProcessPoolExecutor(max_workers=workers)
            results = pool.map(_read_codes, jobs, chunksize=2)
        try:
            for i, row in enumerate(results):
                codes[i] = row
                if progress:
                    progress(i + 1, n)
        finally:
            if pool:
                pool.shutdown(cancel_futures=True)
    os.replace(tmp_path, out_path)
    return out_path


class CaptureStore:
    """Read access to a converted store.  Use as a context manager or call close()."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self._f = h5py.File(self.path, "r")
        self._codes = self._f["codes"]
        a = self._f.attrs
        self.lsb_volts = float(a["lsb_volts"])
        self.dt_s = float(a["dt_s"])
        self.n_samples = int(a["n_samples"])
        idx = self._f["index"]
        self.time_s = idx["time_s"][:]
        self.file_number = idx["file_number"][:]
        self.created = idx["created"].asstr()[:]
        self.source = idx["source"].asstr()[:]

    def __len__(self) -> int:
        return self._codes.shape[0]

    @property
    def live_time_s(self) -> float:
        """Total recorded time (sum of capture durations), excluding the gaps between captures."""
        return len(self) * self.n_samples * self.dt_s

    def codes(self, i: int) -> np.ndarray:
        return self._codes[i]

    def volts(self, i: int, dtype=np.float64) -> np.ndarray:
        return self._codes[i].astype(dtype) * dtype(self.lsb_volts)

    def iter_volts(self, dtype=np.float64) -> Iterator[tuple[int, np.ndarray]]:
        for i in range(len(self)):
            yield i, self.volts(i, dtype)

    def close(self) -> None:
        self._f.close()

    def __enter__(self) -> "CaptureStore":
        return self

    def __exit__(self, *exc) -> None:
        self.close()
