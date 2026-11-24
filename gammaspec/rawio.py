"""Reading the ASCII ADC capture files.

Each capture file looks like::

    Creation: 2021-08-26 13:34:48.617564
    Signals:  524288
    Dt:       5.555555555555556e-04 ms
    +3.390868e-04
    -8.477170e-05
    ...

The sample values are ADC codes multiplied by a fixed step (LSB) and printed
with 7 significant digits.  Converting them back to integer codes is lossless:
every code maps to exactly one printed string.
"""

from __future__ import annotations

import calendar
import csv
import datetime as _dt
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

#: ADC step of the 2021 acquisition system, fitted from large codes (volts).
DEFAULT_LSB_VOLTS = 8.4771697e-05

_TIME_UNITS = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "µs": 1e-6, "ns": 1e-9}
_FILE_NUMBER_RE = re.compile(r"(\d+)(?:\.\w+)?$")
_HEADER_LINES = 3
# Fixed-width sample line: "+d.dddddde-dd\n"
_LINE_WIDTH = 14


class RawFormatError(ValueError):
    """The file does not look like an ADC capture file."""


@dataclass(frozen=True)
class RawHeader:
    path: Path
    file_number: int
    created: _dt.datetime  # naive local wall-clock time of the acquisition PC
    n_samples: int
    dt_s: float
    header_bytes: int

    @property
    def duration_s(self) -> float:
        return self.n_samples * self.dt_s

    @property
    def time_s(self) -> float:
        """Wall-clock creation time as seconds since 1970 (naive local time encoded as if UTC)."""
        return calendar.timegm(self.created.timetuple()) + self.created.microsecond * 1e-6


def _file_number(path: Path) -> int:
    m = _FILE_NUMBER_RE.search(path.name)
    return int(m.group(1)) if m else -1


def read_header(path: str | Path) -> RawHeader:
    path = Path(path)
    with open(path, "rb") as f:
        lines = [f.readline() for _ in range(_HEADER_LINES)]
    header_bytes = sum(len(line) for line in lines)
    try:
        text = [line.decode("ascii").strip() for line in lines]
        key_vals = dict(line.split(":", 1) for line in text)
        created = _dt.datetime.fromisoformat(key_vals["Creation"].strip())
        n_samples = int(key_vals["Signals"].strip())
        dt_value, dt_unit = key_vals["Dt"].split()
        dt_s = float(dt_value) * _TIME_UNITS[dt_unit]
    except (UnicodeDecodeError, ValueError, KeyError) as exc:
        raise RawFormatError(f"{path}: unrecognised header {lines!r}") from exc
    return RawHeader(path, _file_number(path), created, n_samples, dt_s, header_bytes)


def _parse_fixed_width(body: bytes) -> np.ndarray | None:
    """Vectorised parser for the fixed-width "+d.dddddde-dd" format; None if the body doesn't match."""
    if body and not body.endswith(b"\n"):
        body += b"\n"
    if len(body) % _LINE_WIDTH:
        return None
    a = np.frombuffer(body, dtype=np.uint8).reshape(-1, _LINE_WIDTH)
    sign, exp_sign = a[:, 0], a[:, 10]
    digits = a[:, [1, 3, 4, 5, 6, 7, 8, 11, 12]]
    if not (
        np.all(a[:, 13] == ord("\n"))
        and np.all(a[:, 2] == ord("."))
        and np.all(a[:, 9] == ord("e"))
        and np.all((sign == ord("+")) | (sign == ord("-")))
        and np.all((exp_sign == ord("+")) | (exp_sign == ord("-")))
        and np.all((digits >= ord("0")) & (digits <= ord("9")))
    ):
        return None
    d = digits.astype(np.int64) - ord("0")
    mantissa = d[:, :7] @ (10 ** np.arange(6, -1, -1, dtype=np.int64))
    exponent = np.where(exp_sign == ord("-"), -1, 1) * (d[:, 7] * 10 + d[:, 8]) - 6
    values = mantissa * np.power(10.0, exponent)
    return np.where(sign == ord("-"), -values, values)


def parse_samples(body: bytes) -> np.ndarray:
    """Parse the sample section of a capture file into volts (float64)."""
    values = _parse_fixed_width(body)
    if values is None:  # generic (slower) fallback for any whitespace-separated floats
        try:
            values = np.array(body.split(), dtype=np.float64)
        except ValueError as exc:
            raise RawFormatError("sample section contains non-numeric values") from exc
    return values


def volts_to_codes(volts: np.ndarray, lsb: float, tol: float = 0.02) -> np.ndarray:
    """Convert volts to integer ADC codes, checking that every value sits on the ADC grid."""
    x = np.asarray(volts, dtype=np.float64) / lsb
    codes = np.rint(x)
    residual = float(np.max(np.abs(x - codes), initial=0.0))
    if residual > tol:
        raise RawFormatError(
            f"values are not multiples of LSB={lsb:.10e} V (max residual {residual:.3f} LSB)"
        )
    if codes.size and (codes.min() < np.iinfo(np.int32).min or codes.max() > np.iinfo(np.int32).max):
        raise RawFormatError("ADC codes overflow int32")
    return codes.astype(np.int32)


def estimate_lsb(volts: np.ndarray) -> float:
    """Estimate the ADC step from the data itself (for captures from unknown hardware)."""
    u = np.unique(np.asarray(volts, dtype=np.float64))
    u = u[u != 0]
    if u.size < 2:
        raise ValueError("need at least two distinct non-zero values to estimate the ADC step")
    # Near zero the printed values carry full precision, so neighbouring steps give a good first guess.
    near = np.sort(np.abs(u))[:64]
    steps = np.diff(near)
    steps = steps[steps > 0]
    lsb = float(np.median(steps[steps < 1.5 * steps.min()]))
    # Refine by least squares over progressively larger values; each pass keeps rounding unambiguous.
    bound = 50 * lsb
    max_abs = float(np.abs(u).max())
    while True:
        sel = u[np.abs(u) <= bound]
        codes = np.rint(sel / lsb)
        lsb = float(np.dot(sel, codes) / np.dot(codes, codes))
        if bound >= max_abs:
            return lsb
        bound *= 8


def read_raw_file(path: str | Path, lsb: float | None = DEFAULT_LSB_VOLTS) -> tuple[RawHeader, np.ndarray]:
    """Read one capture file.  Returns (header, int32 codes); if lsb is None, returns volts instead."""
    header = read_header(path)
    with open(header.path, "rb") as f:
        f.seek(header.header_bytes)
        body = f.read()
    volts = parse_samples(body)
    if volts.size != header.n_samples:
        raise RawFormatError(f"{header.path}: header says {header.n_samples} samples, found {volts.size}")
    if lsb is None:
        return header, volts
    return header, volts_to_codes(volts, lsb)


def build_index(directory: str | Path, pattern: str = "raw_data_*") -> list[RawHeader]:
    """Headers of every capture file in a directory, sorted chronologically."""
    paths = [p for p in Path(directory).glob(pattern) if p.is_file()]
    if not paths:
        raise FileNotFoundError(f"no files matching {pattern!r} in {directory}")
    headers = [read_header(p) for p in paths]
    headers.sort(key=lambda h: (h.created, h.file_number))
    return headers


def write_index_csv(headers: list[RawHeader], out_path: str | Path) -> None:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    t0 = headers[0].time_s if headers else 0.0
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["order", "file_number", "created", "t_rel_s", "n_samples", "dt_s", "path"])
        for i, h in enumerate(headers):
            w.writerow([i, h.file_number, h.created.isoformat(), f"{h.time_s - t0:.6f}",
                        h.n_samples, repr(h.dt_s), str(h.path)])
