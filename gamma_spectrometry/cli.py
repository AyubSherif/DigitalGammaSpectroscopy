"""Command line interface:  python -m gamma_spectrometry <command> ...

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

    from . import config as C
    from .pipeline import PipelineConfig, plot_event_spectrum, run, summarize

    d = C.read_dict(args.config) if args.config else {}
    run_keys = {f.name for f in C.dataclasses.fields(C.RunConfig)} | {"schema_version"}
    if "pipeline" in d:                                 # a run file: use its [pipeline] tables
        d = d["pipeline"]
    elif set(d) & (run_keys - {"pipeline"}):            # a run file without pipeline settings
        d = {}
    for s in args.set or []:
        C.set_key(d, s[len("pipeline."):] if s.startswith("pipeline.") else s)
    # the flags below are shortcuts: they override the file when given
    flags = {"preprocess.line_removal": False if args.no_line_removal else None,
             "preprocess.pole_zero": False if args.no_pole_zero else None,
             "preprocess.baseline_tracker": False if args.no_tracker else None,
             "detect.method": args.method, "detect.threshold": args.threshold, "events.snippets": args.snippets,
             "solver.enabled": False if args.no_pileup_solver else None,
             "solver.residual_search": False if args.no_residual_search else None}
    for k, v in flags.items():
        if v is not None:
            a, b = k.split(".")
            d.setdefault(a, {})[b] = v
    try:
        cfg = C.from_dict(d, PipelineConfig)
    except C.ConfigError as exc:
        print(exc)
        return 2
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


def _cmd_evaluate_pileup(args) -> int:
    from .pileup import METHODS
    from .pileup_eval import LABELS, SEPARATIONS, evaluate

    s = evaluate(args.store, args.template, args.events, args.out)
    print("\nA (template pairs, detected) - second pulse: bias / spread in %, by separation")
    for pname, rows in s["template_pairs"].items():
        print(f"  {pname}  (detected as a pair: " + " ".join(
            f"{d}:{rows[str(d)]['detected_as_pair']:.2f}" for d in SEPARATIONS) + ")")
        for m in METHODS:
            print(f"    {LABELS[m]:30s}" + " ".join(
                f"{d}:{100*rows[str(d)]['methods'][m]['second']['bias']:+.1f}/{100*rows[str(d)]['methods'][m]['second']['spread']:.1f}"
                for d in SEPARATIONS))
    print("\nB (real pulse pairs) - first | second pulse: bias / spread in %")
    for m in METHODS:
        print(f"  {LABELS[m]:30s}" + " ".join(
            f"{d}:{100*s['real_pairs'][str(d)][m]['first']['bias']:+.1f}/{100*s['real_pairs'][str(d)][m]['first']['spread']:.1f}"
            f"|{100*s['real_pairs'][str(d)][m]['second']['bias']:+.1f}/{100*s['real_pairs'][str(d)][m]['second']['spread']:.1f}"
            for d in SEPARATIONS))
    print("\nspeed (microseconds per 2-pulse cluster): "
          + ", ".join(f"{m} {v:.0f}" for m, v in s["speed_us_per_cluster"].items()))
    print(f"\nC (real 60Co pile-up events; {s['n_real_isolated']} isolated references): "
          "peak shift vs isolated (%), FWHM (%)")
    for m in METHODS:
        r = s["real_events"][m]
        iso = r.get("isolated", {})
        line = f"  {LABELS[m]:30s} isolated FWHM {iso.get('fwhm_pct', ['?', '?'])} |"
        for lab in ("all pile-up", "gap < 10", "gap 10-30", "gap 30-60"):
            x = r.get(lab, {})
            if "shift_pct" in x:
                line += f" {lab} (n={x['n']}): shift {x['shift_pct'][0]:+.2f}/{x['shift_pct'][1]:+.2f}, FWHM {x['fwhm_pct'][0]:.2f}/{x['fwhm_pct'][1]:.2f} |"
            else:
                line += f" {lab}: {x.get('error', 'n/a')[:40]} |"
        print(line)
    return 0


def _cmd_report(args) -> int:
    import json
    from pathlib import Path

    from .pipeline import plot_event_spectrum, spectrum_report

    guess = [float(v) for v in args.peaks.split(",")] if args.peaks else None
    r = spectrum_report(args.events, n_peaks=args.n_peaks, guess=guess)
    print(f"valid events {r['valid_events']}, with a usable amplitude {r['good_events']} "
          f"(throughput {100 * r['throughput']:.1f} %; isolated only: {100 * r['isolated_fraction']:.1f} %); "
          f"pulses found by the solver: {r['fit_added']}")
    print(f"{'category':46s} {'n':>7s}  {'centres (V)':>18s}  {'shift vs isolated (%)':>22s}  {'FWHM (%)':>16s}")
    for name, c in r["categories"].items():
        if "error" in c:
            print(f"{name:46s} {c['n']:7d}  not fitted: {c['error'][:50]}")
            continue
        print(f"{name:46s} {c['n']:7d}  {' / '.join(f'{x:.4f}' for x in c['centers']):>18s}  "
              f"{' / '.join(f'{x:+.2f}' for x in c['shift_pct']):>22s}  "
              f"{' / '.join(f'{x:.2f}' for x in c['fwhm_pct']):>16s}")
    out = Path(args.events).with_suffix(".report.json")
    out.write_text(json.dumps(r, indent=2), encoding="utf-8")
    plot_event_spectrum(args.events, Path(args.events).with_suffix(".png"))
    print(f"written {out}")
    return 0


def _cmd_benchmark(args) -> int:
    from .benchmark import run_benchmark

    caps = range(args.first, args.first + args.n) if args.n else None
    try:
        s = run_benchmark(args.store, args.events, args.out, captures=caps)
    except ValueError as exc:   # too few counts to calibrate a chain
        print(f"benchmark failed: {exc}. Use more captures (--n; 100 recommended).")
        return 1
    print(f"\n{'chain':46s} {'FWHM 1173/1332 %':>18s} {'keV':>13s} {'valley/peak':>12s} "
          f"{'peaks cps':>13s} {'pile-up cps':>12s} {'live %':>7s}")
    for r in s["chains"]:
        pr = r["peak_rates_cps"]
        print(f"{r['name']:46s} {r['fwhm_pct'][0]:8.2f}/{r['fwhm_pct'][1]:<8.2f} "
              f"{r['fwhm_kev'][0]:6.1f}/{r['fwhm_kev'][1]:<6.1f} {r['valley_to_peak']:12.3f} "
              f"{pr['Co-60 1173'][0]:6.0f}/{pr['Co-60 1332'][0]:<6.0f} {r['pileup_continuum_cps']:12.1f} "
              f"{100 * r['live_fraction']:7.1f}")
    print(f"typical 3x3 NaI(Tl) with an analog chain (another detector): {s['literature_3x3_nai_analog']}")
    print(f"outputs in {args.out}")
    return 0


def _sim_config(args):
    """SimConfig from --config ([simulation.signal]) and --set KEY=VALUE (keys of [simulation.signal])."""
    from . import config as C

    d = C.read_dict(args.config) if args.config else {}
    sig = d.get("simulation", {}).get("signal", d.get("signal", {}))
    full = {"simulation": {"signal": dict(sig)}}
    for s in args.set or []:
        key = s.split("=", 1)[0]
        prefix = "" if key.startswith("simulation.") else ("simulation." if key.startswith("signal.") else "simulation.signal.")
        C.set_key(full, prefix + s)
    return C.from_dict(full).simulation.signal


def _cmd_simulate(args) -> int:
    from . import config as C

    try:
        sig = _sim_config(args)
    except C.ConfigError as exc:
        print(exc)
        return 2
    if args.out == "-":
        import numpy as np

        from .sources import Paced, SimulationSource

        src = SimulationSource(sig, args.duration)
        it = (Paced(src, args.speed) if args.realtime else src).chunks()
        out = sys.stdout.buffer
        dtype = {"f4": "<f4", "i2": "<i2", "i4": "<i4"}[args.dtype]
        try:
            for ch in it:
                x = ch.data if args.dtype == "f4" else np.rint(ch.data / sig.lsb_v)
                out.write(x.astype(dtype).tobytes())
            out.flush()
        except (BrokenPipeError, KeyboardInterrupt):
            pass
        return 0
    from .simulate import simulate_store

    store, truth = simulate_store(args.out, sig, args.captures, args.samples, args.gap, overwrite=args.overwrite)
    from .simulate import load_truth

    T = load_truth(truth)
    live = args.captures * args.samples / sig.fs
    print(f"wrote {store}: {args.captures} captures of {args.samples} samples ({live:.2f} s), "
          f"{T.size} pulses ({T.size / live:.0f} /s); truth: {truth}")
    return 0


def _cmd_stream(args) -> int:
    from pathlib import Path

    from . import config as C
    from .runner import StepError, run

    d = C.read_dict(args.config) if args.config else {}
    src = args.source
    if src == "sim":
        kind = "simulation"
    elif src == "-" or not (Path(src).is_dir() or Path(src).suffix.lower() in (".h5", ".hdf5")):
        kind = "pipe"
    else:
        kind = "replay"
    sets = [f"output={C._toml_value(args.out)}", f"input.kind={C._toml_value(kind)}"]
    if kind != "simulation":
        sets.append(f"input.source={C._toml_value(src)}")
    if args.template:
        sets.append(f"template.path={C._toml_value(args.template)}")
    if args.background:
        sets.append(f"input.background={C._toml_value(args.background)}")
    if args.duration is not None:
        sets.append(f"simulation.duration_s={args.duration}")
    if kind == "simulation" and not args.background_duration:
        sets.append("simulation.background_duration_s=0")
    elif args.background_duration:
        sets.append(f"simulation.background_duration_s={args.background_duration}")
    if args.realtime:
        sets += ["stream.realtime=true", f"stream.speed={args.speed}"]
    try:
        for s in sets + (args.set or []):
            C.set_key(d, s)
        cfg = C.from_dict(d)
        plot = None
        if args.show:
            from .liveplot import LivePlot

            plot = LivePlot()
        run(cfg, on_update=plot.update if plot else None)
        if plot:
            plot.keep_open()
    except (C.ConfigError, StepError) as exc:
        print(exc)
        return 2
    return 0


def _cmd_app(args) -> int:
    try:
        import streamlit  # noqa: F401
    except ImportError:
        print("the app needs Streamlit: pip install \"gamma-spectrometry[app]\" (or: pip install streamlit)")
        return 1
    from .app import launch

    return launch(args.port, open_browser=not args.no_browser, address=args.address)


def _cmd_run(args) -> int:
    from . import config as C
    from .runner import StepError, run

    try:
        cfg = C.load(args.config, args.set)
        run(cfg)
    except (C.ConfigError, StepError) as exc:
        print(exc)
        return 2
    return 0


def _cmd_config(args) -> int:
    from pathlib import Path

    from . import config as C

    try:
        if args.action == "init":
            if Path(args.path).exists() and not args.force:
                print(f"{args.path} exists (use --force to replace it)")
                return 1
            cfg = C.load(None, args.set)
            print(f"wrote {C.save(cfg, args.path)}")
        elif args.action == "check":
            cfg = C.load(args.path, args.set)
            warns = C.check(cfg)
            if not args.quiet:
                print(C.dumps(cfg, comments=False))
            for w in warns:
                print(f"warning: {w}")
            print(f"{args.path}: OK" + (f" ({len(warns)} warning(s))" if warns else ""))
        else:
            print(C.dumps(C.pipeline_from_event_file(args.events)))
    except C.ConfigError as exc:
        print(exc)
        return 2
    return 0


def _cmd_spectrum(args) -> int:
    from .analysis import analyze

    s = analyze(args.source, args.out, background_events=args.background, bin_kev=args.bin_kev,
                e_max_kev=args.emax)
    for lab, e, r in s["calibration_source"]["residuals_kev"]:
        print(f"  source line {lab:10s} {e:8.1f} keV: residual {r:+.2f} keV")
    if "lines" in s:
        print(f"background calibration: {s['calibration_background']['kind']}")
        print(f"{'line':8s} {'keV':>8s} {'source cps':>14s} {'background cps':>16s} {'net cps':>16s}")
        for k in range(len(s["lines"]["net"])):
            ls, lb, ln = (s["lines"][x][k] for x in ("source", "background", "net"))
            print(f"{ln['label']:8s} {ln['energy_kev']:8.1f} {ls['rate_cps']:8.2f}+-{ls['rate_err_cps']:.2f} "
                  f"{lb['rate_cps']:9.2f}+-{lb['rate_err_cps']:.2f} {ln['rate_cps']:9.2f}+-{ln['rate_err_cps']:.2f}")
        t = s["totals_cps"]
        print(f"total rate above {s['reliable_above_kev']:.0f} keV: source {t['source']:.0f} /s, "
              f"background {t['background']:.0f} /s, net {t['net']:.0f} /s")
    print("gain drift (source run, % per chunk): "
          + ", ".join(f"{g['relative_pct']:+.2f}" for g in s["gain_vs_time_source"]))
    print(f"outputs in {args.out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="gamma-spectrometry", description=__doc__,
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
    pp.add_argument("--method", default=None, choices=["optimal", "matched", "sharpened", "edge"],
                    help="detection filter (default: from the config, optimal)")
    pp.add_argument("--threshold", type=float, default=None, help="detection threshold in filter sigma (default 5.5)")
    pp.add_argument("--no-pileup-solver", action="store_true", help="skip the pile-up solver")
    pp.add_argument("--no-residual-search", action="store_true",
                    help="pile-up solver without the search for hidden pulses")
    pp.add_argument("--snippets", default=None, choices=["none", "pileup", "all"],
                    help="which clusters get waveform snippets (default pileup)")
    pp.add_argument("--limit", type=int, help="only the first N captures")
    pp.add_argument("--overwrite", action="store_true")
    pp.add_argument("--config", help="run file (.toml/.json); its [pipeline] tables are used")
    pp.add_argument("--set", action="append", metavar="KEY=VALUE",
                    help="override one setting, e.g. --set detect.threshold=6 (repeatable)")
    pp.set_defaults(func=_cmd_process)

    pm = sub.add_parser("simulate", help="simulated detector signal: a capture store (+ truth), or samples on stdout")
    pm.add_argument("out", help="output store (.h5), or - for binary samples on stdout (to pipe into 'stream -')")
    pm.add_argument("--captures", type=int, default=20, help="store: number of captures (default 20)")
    pm.add_argument("--samples", type=int, default=524288, help="store: samples per capture (default 524288)")
    pm.add_argument("--gap", type=float, default=0.0, help="store: seconds between captures (default 0: continuous)")
    pm.add_argument("--duration", type=float, default=0.0, help="stdout: seconds of signal (default 0: endless)")
    pm.add_argument("--realtime", action="store_true", help="stdout: at real time (a live source)")
    pm.add_argument("--speed", type=float, default=1.0, help="stdout with --realtime: multiple of real time")
    pm.add_argument("--dtype", default="f4", choices=["f4", "i2", "i4"], help="stdout: sample format (default f4 volts)")
    pm.add_argument("--config", help="run file; its [simulation.signal] table is used")
    pm.add_argument("--set", action="append", metavar="KEY=VALUE",
                    help="simulator setting, e.g. --set nuclide=cs137 --set rate_cps=2000 (repeatable)")
    pm.add_argument("--overwrite", action="store_true")
    pm.set_defaults(func=_cmd_simulate)

    pt = sub.add_parser("stream", help="process a stream: stdin (-), a replayed recording, or a simulation (sim)")
    pt.add_argument("source", help="- (binary samples on stdin), a raw folder or .h5 store (replay), a file/FIFO "
                                   "of binary samples, or sim (the simulator)")
    pt.add_argument("out", help="output folder (events, live spectrum, spectrum)")
    pt.add_argument("--template", help="pulse template (.h5); default: built from --background (or the simulator's pulse)")
    pt.add_argument("--background", help="background run (raw folder or .h5 store) for background subtraction")
    pt.add_argument("--duration", type=float, help="sim: seconds of signal (default 30)")
    pt.add_argument("--background-duration", type=float, default=0.0, help="sim: also simulate a background run this long")
    pt.add_argument("--realtime", action="store_true", help="replay / sim at real time (frames are dropped if processing falls behind)")
    pt.add_argument("--speed", type=float, default=1.0, help="with --realtime: multiple of real time")
    pt.add_argument("--show", action="store_true", help="live spectrum window")
    pt.add_argument("--config", help="run file to start from")
    pt.add_argument("--set", action="append", metavar="KEY=VALUE", help="override a run-file setting (repeatable)")
    pt.set_defaults(func=_cmd_stream)

    pa = sub.add_parser("app", help="open the browser UI (choose input, stages and settings; run; live spectrum)")
    pa.add_argument("--port", type=int, help="port (default 8501)")
    pa.add_argument("--no-browser", action="store_true", help="do not open a browser window")
    pa.add_argument("--address", default="localhost",
                    help="network address to listen on (default localhost: this computer only)")
    pa.set_defaults(func=_cmd_app)

    pr = sub.add_parser("run", help="run the configured chain: input -> template -> events -> spectrum -> benchmark")
    pr.add_argument("config", help="run file (.toml or .json), e.g. from 'gamma-spectrometry config init'")
    pr.add_argument("--set", action="append", metavar="KEY=VALUE",
                    help="override one setting, e.g. --set pipeline.detect.threshold=6 (repeatable)")
    pr.set_defaults(func=_cmd_run)

    pg = sub.add_parser("config", help="write, check or show configurations")
    gs = pg.add_subparsers(dest="action", required=True)
    gi = gs.add_parser("init", help="write a run file with every setting at its default, commented")
    gi.add_argument("path", help="output .toml (or .json)")
    gi.add_argument("--set", action="append", metavar="KEY=VALUE", help="start from these values")
    gi.add_argument("--force", action="store_true", help="replace an existing file")
    gc = gs.add_parser("check", help="validate a run file and print the effective configuration")
    gc.add_argument("path")
    gc.add_argument("--set", action="append", metavar="KEY=VALUE")
    gc.add_argument("--quiet", action="store_true", help="only errors and warnings")
    gw = gs.add_parser("show", help="print the pipeline configuration stored in an event file")
    gw.add_argument("events", help="event file (.h5)")
    pg.set_defaults(func=_cmd_config)

    pd = sub.add_parser("evaluate-detection", help="compare detection methods with injected pulses")
    pd.add_argument("store", help="HDF5 capture store (its clean signal provides the noise)")
    pd.add_argument("template", help="template.h5 with a fitted model")
    pd.add_argument("out", help="output directory")
    pd.set_defaults(func=_cmd_evaluate_detection)

    pu = sub.add_parser("evaluate-pileup", help="benchmark the pile-up resolution methods")
    pu.add_argument("store", help="HDF5 capture store (noise and real pulses)")
    pu.add_argument("template", help="template.h5 with a fitted model")
    pu.add_argument("events", help="event file of the same store (pile-up snippets)")
    pu.add_argument("out", help="output directory")
    pu.set_defaults(func=_cmd_evaluate_pileup)

    pr = sub.add_parser("report", help="spectrum peaks by event category, throughput")
    pr.add_argument("events", help="event file")
    pr.add_argument("--n-peaks", type=int, default=2, help="peaks to fit (found automatically)")
    pr.add_argument("--peaks", help="comma-separated starting guesses for the peak centres (V)")
    pr.set_defaults(func=_cmd_report)

    ps = sub.add_parser("spectrum", help="calibrate, build energy spectra, subtract the background")
    ps.add_argument("source", help="event file of the source run")
    ps.add_argument("out", help="output directory")
    ps.add_argument("--background", help="event file of the background run")
    ps.add_argument("--bin-kev", type=float, default=5.0, help="bin width in keV (default 5)")
    ps.add_argument("--emax", type=float, default=3000.0, help="upper end of the spectrum in keV")
    ps.set_defaults(func=_cmd_spectrum)

    pb = sub.add_parser("benchmark", help="Gamma Spectrometry App vs naive find_peaks vs an emulated analog chain")
    pb.add_argument("store", help="HDF5 capture store")
    pb.add_argument("events", help="event file of the same store (the Gamma Spectrometry App chain)")
    pb.add_argument("out", help="output directory")
    pb.add_argument("--first", type=int, default=0, help="first capture (default 0)")
    pb.add_argument("--n", type=int, default=100, help="number of captures (default 100; 0 = all)")
    pb.set_defaults(func=_cmd_benchmark)

    args = p.parse_args(argv)
    if getattr(args, "lsb", None) is None and hasattr(args, "lsb"):
        from .rawio import DEFAULT_LSB_VOLTS
        args.lsb = DEFAULT_LSB_VOLTS
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
