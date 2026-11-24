"""Event list: the reduced form of the data (one row per detected pulse) in HDF5.

Layout::

    /events/<field>      one 1-D dataset per field (see EVENT_DTYPE), appended per capture
    /captures/<field>    one row per processed capture (see CAPTURE_DTYPE)
    /snippets/...        optional waveform snippets (clean signal) of pile-up clusters:
                         cluster, capture, start, offset, length (index) and data (flat float32)
    attrs                config (JSON), template, source store, lsb, dt, version

Columns are compressed (gzip + shuffle) and resizable, so the writer can append capture
by capture (or chunk by chunk from a stream).
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np

from . import __version__

# flag bits
FLAG_PILEUP = 1          # pulse shares its fit window with another pulse (needs pile-up resolution)
FLAG_WARMUP = 2          # in the capture's warm-up region (baseline not fully restored)
FLAG_EDGE = 4            # too close to the capture end for a complete fit window
FLAG_BAD_FIT = 8         # isolated, but the template fit is poor (likely unresolved pile-up)
FLAG_LARGE = 16          # amplitude above the large-pulse limit (e.g. cosmic muons)
FLAG_AFTER_LARGE = 32    # shortly after a large pulse: may be an artifact of it
FLAG_FIT_ADDED = 64      # found by the pile-up solver's residual search (not by the detector)

FLAG_NAMES = {FLAG_PILEUP: "pileup", FLAG_WARMUP: "warmup", FLAG_EDGE: "edge",
              FLAG_BAD_FIT: "bad_fit", FLAG_LARGE: "large", FLAG_AFTER_LARGE: "after_large",
              FLAG_FIT_ADDED: "fit_added"}

#: flags that only describe how an event was measured; its amplitude is still valid
INFO_FLAGS = FLAG_PILEUP | FLAG_FIT_ADDED

EVENT_DTYPE = np.dtype([
    ("capture", "u4"),          # capture number (row in the store, time order)
    ("index", "i8"),            # peak sample within the capture
    ("time_s", "f8"),           # capture time + (index + shift) * dt  (wall-clock seconds)
    ("amplitude", "f4"),        # template-fit amplitude, V (NaN unless an isolated pulse)
    ("filter_amplitude", "f4"),  # detection filter output at the peak, V
    ("shift", "f4"),            # sub-sample peak position relative to index (fit)
    ("fit_rms", "f4"),          # rms residual of the fit, V
    ("flags", "u2"),
    ("cluster", "i8"),          # pulses closer than cluster_gap share a cluster id
    ("cluster_size", "u2"),
    ("gap_prev", "i4"),         # samples to the previous pulse in the capture (-1: none)
    ("gap_next", "i4"),         # samples to the next pulse (-1: none)
])

CAPTURE_DTYPE = np.dtype([
    ("capture", "u4"),
    ("time_s", "f8"),
    ("live_time_s", "f8"),      # duration of the valid region (after warm-up, before the end margin)
    ("offset_v", "f4"),
    ("sigma_in_v", "f4"),       # noise before / after pre-processing
    ("sigma_out_v", "f4"),
    ("filter_sigma_v", "f4"),   # detection filter output noise
    ("threshold_v", "f4"),      # detection threshold
    ("n_events", "u4"),
    ("n_lines", "u2"),          # interference lines removed
])

SNIPPET_DTYPE = np.dtype([("cluster", "i8"), ("capture", "u4"), ("start", "i8"),
                          ("offset", "i8"), ("length", "i4")])


def _create_columns(group: h5py.Group, dtype: np.dtype, chunk: int):
    for name in dtype.names:
        group.create_dataset(name, shape=(0,), maxshape=(None,), dtype=dtype[name],
                             chunks=(chunk,), compression="gzip", compression_opts=4, shuffle=True)


def _append(group: h5py.Group, rows: np.ndarray):
    if rows.size == 0:
        return
    n0 = group[rows.dtype.names[0]].shape[0]
    for name in rows.dtype.names:
        ds = group[name]
        ds.resize((n0 + rows.size,))
        ds[n0:] = rows[name]


class EventWriter:
    def __init__(self, path: str | Path, attrs: dict, overwrite: bool = False):
        self.path = Path(path)
        if self.path.exists() and not overwrite:
            raise FileExistsError(f"{self.path} exists (use overwrite)")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._f = h5py.File(self.path, "w")
        _create_columns(self._f.create_group("events"), EVENT_DTYPE, 65536)
        _create_columns(self._f.create_group("captures"), CAPTURE_DTYPE, 1024)
        sn = self._f.create_group("snippets")
        _create_columns(sn, SNIPPET_DTYPE, 4096)
        sn.create_dataset("data", shape=(0,), maxshape=(None,), dtype="f4", chunks=(65536,),
                          compression="gzip", compression_opts=4, shuffle=True)
        for k, v in attrs.items():
            self._f.attrs[k] = json.dumps(v) if isinstance(v, (dict, list)) or v is None else v
        self._f.attrs["gammaspec_version"] = __version__

    def append(self, events: np.ndarray, capture_row: np.ndarray,
               snippets: list[tuple[int, int, int, np.ndarray]] | None = None) -> None:
        _append(self._f["events"], events)
        _append(self._f["captures"], capture_row)
        if snippets:
            data = self._f["snippets/data"]
            off = data.shape[0]
            rows = np.zeros(len(snippets), SNIPPET_DTYPE)
            chunks = []
            for i, (cl, cap, start, w) in enumerate(snippets):
                rows[i] = (cl, cap, start, off, w.size)
                off += w.size
                chunks.append(w.astype(np.float32))
            flat = np.concatenate(chunks)
            data.resize((data.shape[0] + flat.size,))
            data[-flat.size:] = flat
            _append(self._f["snippets"], rows)

    def close(self) -> None:
        self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class EventFile:
    """Read access to an event file. ``events`` and ``captures`` are structured arrays."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        with h5py.File(self.path, "r") as f:
            self.events = _read_table(f["events"], EVENT_DTYPE)
            self.captures = _read_table(f["captures"], CAPTURE_DTYPE)
            self._snip = _read_table(f["snippets"], SNIPPET_DTYPE)
            self.attrs = {k: _maybe_json(v) for k, v in f.attrs.items()}
        self._snip_by_cluster = {int(c): i for i, c in enumerate(self._snip["cluster"])}

    @property
    def live_time_s(self) -> float:
        return float(self.captures["live_time_s"].sum())

    def clean(self) -> np.ndarray:
        """Events with no flags (isolated, fitted, inside the valid region)."""
        return self.events[self.events["flags"] == 0]

    def good(self, include_pileup: bool = True) -> np.ndarray:
        """Events with a valid amplitude and no problem flags. With ``include_pileup``, pulses
        resolved by the pile-up solver are included (their PILEUP / FIT_ADDED flags only
        describe how they were measured)."""
        allowed = INFO_FLAGS if include_pileup else 0
        ev = self.events
        ok = ((ev["flags"] & np.uint16(~np.uint16(allowed))) == 0) & np.isfinite(ev["amplitude"])
        return ev[ok]

    def snippet(self, cluster: int) -> tuple[int, int, np.ndarray]:
        """(capture, start sample, waveform) of a stored cluster."""
        i = self._snip_by_cluster[int(cluster)]
        r = self._snip[i]
        with h5py.File(self.path, "r") as f:
            w = f["snippets/data"][r["offset"]: r["offset"] + r["length"]]
        return int(r["capture"]), int(r["start"]), w

    @property
    def n_snippets(self) -> int:
        return int(self._snip.size)


def _read_table(group: h5py.Group, dtype: np.dtype) -> np.ndarray:
    n = group[dtype.names[0]].shape[0]
    out = np.zeros(n, dtype)
    for name in dtype.names:
        out[name] = group[name][:]
    return out


def _maybe_json(v):
    if isinstance(v, (str, bytes)):
        s = v.decode() if isinstance(v, bytes) else v
        if s == "null":
            return None
        if s[:1] in "{[":
            try:
                return json.loads(s)
            except ValueError:
                return s
    return v
