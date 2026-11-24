"""Command line interface:  python -m gammaspec <command> ...

    index    DIR [--csv OUT]                     list captures in time order, report gaps and live time
    convert  DIR OUT.h5 [--limit N] [--workers]  convert the first N captures (chronological) to HDF5
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from .rawio import build_index, write_index_csv
from .store import CaptureStore, convert


def _summarise(headers) -> None:
    t = np.array([h.time_s for h in headers])
    gaps = np.diff(t)
    nums = np.array([h.file_number for h in headers])
    live = sum(h.duration_s for h in headers)
    span = t[-1] - t[0] if t.size > 1 else 0.0
    print(f"captures         : {len(headers)}")
    print(f"first / last     : {headers[0].created}  ->  {headers[-1].created}")
    print(f"wall-clock span  : {span / 3600:.3f} h")
    print(f"live time        : {live:.1f} s  ({100 * live / span:.1f}% duty)" if span else f"live time: {live:.1f} s")
    if gaps.size:
        print(f"gap between caps : median {np.median(gaps):.3f} s, min {gaps.min():.3f} s, max {gaps.max():.3f} s")
        print(f"file numbers out of time order: {int(np.sum(np.diff(nums) < 0))} places")


def _cmd_index(args) -> int:
    headers = build_index(args.directory, args.pattern)
    _summarise(headers)
    if args.csv:
        write_index_csv(headers, args.csv)
        print(f"index written to {args.csv}")
    return 0


def _cmd_convert(args) -> int:
    headers = build_index(args.directory, args.pattern)
    if args.limit:
        headers = headers[: args.limit]
    _summarise(headers)
    t0 = time.perf_counter()

    def progress(i, n):
        if i == n or i % max(1, n // 20) == 0:
            el = time.perf_counter() - t0
            print(f"  {i}/{n}  {el:.1f} s  ({i * headers[0].n_samples / el / 1e6:.1f} MS/s)", flush=True)

    out = convert(headers, args.out, lsb=None if args.estimate_lsb else args.lsb,
                  workers=args.workers, overwrite=args.overwrite, progress=progress)
    with CaptureStore(out) as s:
        size = out.stat().st_size
        print(f"wrote {out}: {len(s)} captures, {size / 1e6:.1f} MB "
              f"({size / len(s) / 1e6:.2f} MB/capture), lsb={s.lsb_volts:.10e} V, live time {s.live_time_s:.1f} s")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="gammaspec", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    pi = sub.add_parser("index", help="list captures in chronological order")
    pi.add_argument("directory")
    pi.add_argument("--pattern", default="raw_data_*")
    pi.add_argument("--csv", help="write the index to this CSV file")
    pi.set_defaults(func=_cmd_index)

    pc = sub.add_parser("convert", help="convert ASCII captures to a compressed HDF5 store")
    pc.add_argument("directory")
    pc.add_argument("out")
    pc.add_argument("--pattern", default="raw_data_*")
    pc.add_argument("--limit", type=int, help="only the first N captures in time order")
    pc.add_argument("--workers", type=int, help="parallel parser processes (default: cores-1, max 8)")
    pc.add_argument("--lsb", type=float, default=None, help="ADC step in volts (default: built-in value)")
    pc.add_argument("--estimate-lsb", action="store_true", help="estimate the ADC step from the data")
    pc.add_argument("--overwrite", action="store_true")
    pc.set_defaults(func=_cmd_convert)

    args = p.parse_args(argv)
    if getattr(args, "lsb", None) is None and hasattr(args, "lsb"):
        from .rawio import DEFAULT_LSB_VOLTS
        args.lsb = DEFAULT_LSB_VOLTS
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
