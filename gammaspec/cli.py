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


def _cmd_characterize(args) -> int:
    from .characterize import CharacterizeConfig, characterize

    cfg = CharacterizeConfig(limit=args.limit)
    for name in ("template_post", "isolation_pre", "isolation_post", "trigger_threshold",
                 "min_amplitude_sigma", "max_amplitude_v", "resolution"):
        v = getattr(args, name)
        if v is not None:
            setattr(cfg, name, v)
    s = characterize(args.store, args.out, cfg, compare=args.compare)
    m, n, tp = s["model"], s["noise"], s["template"]
    print(f"count rate {s['count_rate_cps']:.0f} cps | noise sigma {n['sigma_v_diff'] * 1e3:.2f} mV | "
          f"template from {tp['n_pulses_used']} pulses | undershoot {tp['undershoot_min']:.4f} "
          f"at {tp['undershoot_min_t']:.0f} samples")
    print("model taus (samples): " + ", ".join(f"{t:.3g}" for t in m["taus_samples"])
          + f" | rms {m['rms']:.2e}")
    return 0


def _cmd_evaluate_restoration(args) -> int:
    from .restoration import evaluate

    s = evaluate(args.store, args.template, args.out, n_peaks=args.n_peaks, limit=args.limit)
    print(s["pole_zero"])
    def fw(d):
        return " / ".join(f"{v:.2f}" for v in d.get("fwhm_pct", [])) or "n/a"

    def ce(d):
        return " / ".join(f"{v:.4f}" for v in d.get("centers_v", [])) or "n/a"

    print(f"{'configuration':30s} {'sigma mV':>9s} {'base sigma mV':>14s} {'|base|>5mV':>11s} {'tail@60':>9s} "
          f"{'FWHM % (b free)':>18s} {'FWHM % (b=0)':>16s}")
    for name, r in s["configs"].items():
        bl = r["pre_pulse_baseline_mv"]
        print(f"{name:30s} {r['sigma_mv']:9.2f} {bl['robust_sigma']:14.2f} {bl['frac_beyond_5mv']:11.3f} "
              f"{r['tail_at']['60']:9.4f} {fw(r['baseline_free']):>18s} {fw(r['baseline_fixed']):>16s}")
    print()
    print("pulses ON THE TAIL of an earlier pulse (60-300 samples after it) vs CLEAN (>300):")
    for name, r in s["configs"].items():
        g = r["by_previous_gap"]
        for lab in ("baseline_free", "baseline_fixed"):
            print(f"{name:30s} {lab:15s} on-tail n={g['on_tail']['n']:6d} centres {ce(g['on_tail'][lab]):>16s} "
                  f"FWHM {fw(g['on_tail'][lab]):>12s} | clean centres {ce(g['clean'][lab]):>16s} "
                  f"FWHM {fw(g['clean'][lab]):>12s}")
    return 0


def _cmd_process(args) -> int:
    from pathlib import Path

    from .detect import DetectConfig
    from .pipeline import EventConfig, PipelineConfig, plot_event_spectrum, run, summarize
    from .preprocess import PreprocessConfig

    cfg = PipelineConfig(
        preprocess=PreprocessConfig(line_removal=not args.no_line_removal, pole_zero=not args.no_pole_zero,
                                    baseline_tracker=not args.no_tracker),
        detect=DetectConfig(method=args.method, threshold=args.threshold),
        events=EventConfig(snippets=args.snippets),
    )
    r = run(args.store, args.template, args.out, cfg, limit=args.limit, overwrite=args.overwrite)
    s = summarize(args.out)
    png = Path(args.out).with_suffix(".png")
    plot_event_spectrum(args.out, png)
    print(f"{r['events']} events from {r['captures']} captures in {r['elapsed_s']:.1f} s "
          f"({r['samples_per_s'] / 1e6:.1f} MS/s, {r['samples_per_s'] / 1.8e6:.1f}x real time)")
    print(f"rate {s['rate_cps']:.0f} /s | clean (isolated, fitted) {s['clean']} | snippets {s['snippets']} | "
          f"filter noise {s['median_filter_sigma_mv']:.2f} mV, threshold {s['median_threshold_mv']:.2f} mV")
    print("flags (fraction of valid events): " + ", ".join(
        f"{k} {v:.3f}" for k, v in s["flag_fractions_of_valid"].items()))
    print(f"event file {r['file_bytes'] / 1e6:.2f} MB vs {r['text_bytes'] / 1e6:.0f} MB of ASCII input "
          f"({r['text_bytes'] / r['file_bytes']:.0f}x smaller); spectrum plot {png}")
    return 0


def _cmd_evaluate_detection(args) -> int:
    from .detection_eval import evaluate

    evaluate(args.store, args.template, args.out)
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

    ph = sub.add_parser("characterize", help="build the pulse template and characterise noise")
    ph.add_argument("store", help="HDF5 capture store")
    ph.add_argument("out", help="output directory")
    ph.add_argument("--limit", type=int, help="only the first N captures")
    ph.add_argument("--template-post", type=int, help="samples kept after the peak (default 1500)")
    ph.add_argument("--isolation-pre", type=int, help="pulse-free samples required before (default 1500)")
    ph.add_argument("--isolation-post", type=int, help="pulse-free samples required after (default 1500)")
    ph.add_argument("--trigger-threshold", type=float, help="coarse trigger threshold in sigma (default 5)")
    ph.add_argument("--min-amplitude-sigma", type=float, help="min template pulse height in sigma (default 30)")
    ph.add_argument("--max-amplitude-v", type=float, help="max template pulse height in volts (default 4)")
    ph.add_argument("--resolution", type=int, help="template points per sample (default 16)")
    ph.add_argument("--compare", help="reference template.h5 to overlay")
    ph.set_defaults(func=_cmd_characterize)

    pe = sub.add_parser("evaluate-restoration",
                        help="compare line removal / pole-zero / baseline tracker on real data")
    pe.add_argument("store", help="HDF5 capture store")
    pe.add_argument("template", help="template.h5 with a fitted model (from characterize)")
    pe.add_argument("out", help="output directory")
    pe.add_argument("--n-peaks", type=int, default=2, help="spectral peaks to fit (2 for 60Co, 1 for 40K)")
    pe.add_argument("--limit", type=int, help="only the first N captures")
    pe.set_defaults(func=_cmd_evaluate_restoration)

    pp = sub.add_parser("process", help="run the pipeline: capture store -> event list")
    pp.add_argument("store", help="HDF5 capture store")
    pp.add_argument("out", help="output event file (.h5); a spectrum plot is written next to it")
    pp.add_argument("--template", required=True, help="template.h5 with a fitted model (from characterize)")
    pp.add_argument("--no-line-removal", action="store_true", help="skip interference line removal")
    pp.add_argument("--no-pole-zero", action="store_true", help="skip pole-zero correction")
    pp.add_argument("--no-tracker", action="store_true", help="skip the baseline tracker")
    pp.add_argument("--method", default="optimal", choices=["optimal", "matched", "sharpened", "edge"],
                    help="detection filter (default optimal)")
    pp.add_argument("--threshold", type=float, default=5.5, help="detection threshold in filter sigma")
    pp.add_argument("--snippets", default="pileup", choices=["none", "pileup", "all"],
                    help="which clusters get waveform snippets (default pileup)")
    pp.add_argument("--limit", type=int, help="only the first N captures")
    pp.add_argument("--overwrite", action="store_true")
    pp.set_defaults(func=_cmd_process)

    pd = sub.add_parser("evaluate-detection", help="compare detection methods with injected pulses")
    pd.add_argument("store", help="HDF5 capture store (its clean signal provides the noise)")
    pd.add_argument("template", help="template.h5 with a fitted model")
    pd.add_argument("out", help="output directory")
    pd.set_defaults(func=_cmd_evaluate_detection)

    args = p.parse_args(argv)
    if getattr(args, "lsb", None) is None and hasattr(args, "lsb"):
        from .rawio import DEFAULT_LSB_VOLTS
        args.lsb = DEFAULT_LSB_VOLTS
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
