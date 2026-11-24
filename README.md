# gammaspec

Digital pulse processing for gamma-ray spectroscopy: turns raw, high-rate ADC
waveforms from a scintillation detector into a clean, pile-up–corrected energy
spectrum, entirely in software.

Developed against 3 hours of data from a NaI(Tl) detector read out directly from
the photomultiplier (PMT) and digitised at 1.8 MS/s, with a ⁶⁰Co source and a
separate background run.

> **Status: v0.1.** The full offline chain works: ingestion, characterisation,
> pre-processing, detection, pile-up resolution, calibration and background
> subtraction, run on the full 3-hour dataset (6.4 M events in 19 min, 75× smaller
> than the raw text), and benchmarked against a naive peak count and an emulated analog
> chain. Every stage is selectable from a run file (`gammaspec run`). Input can be raw files,
> stores, a live stream (stdin / FIFO, or a replay at real time) or a built-in signal
> simulator with known truth. Next: a UI.
> See [Roadmap](#roadmap).

## Why

Raw ADC captures are large (about 7 MB of text per 0.29 s capture) and mostly
noise. Simple thresholding loses real pulses, because the baseline is rarely at
zero: each pulse is followed by a slow undershoot, and at a few thousand counts
per second almost every pulse starts on the tail of an earlier one. Pulses that
overlap distort the measured energies. `gammaspec` addresses each of these with
a modular processing chain:

```text
raw ASCII captures
  → lossless binary store        (done)
  → interference removal         (done)
  → baseline restoration         (done)
  → event detection              (done)  ← the data reduction step
  → pile-up resolution & energy  (done)
  → histogram & calibration      (done)
  → background subtraction       (done)
```

Each processing stage is optional and configurable from a run file, so users can switch
individual corrections on or off without code (see [docs/configuration.md](docs/configuration.md)
and [ADR-0003](docs/adr/0003-modular-optional-stages.md)).

## Installation

Requires Python ≥ 3.11.

```bash
git clone https://github.com/AyubSherif/GammaSpectroscopySoftwareAddOns.git
cd GammaSpectroscopySoftwareAddOns
python -m venv .venv
# Windows: .venv\Scripts\activate    Linux/macOS: source .venv/bin/activate
pip install -e ".[dev]"
```

## Quick start

The whole chain, configured by a run file:

```bash
gammaspec config init my_run.toml        # every setting, at its default, commented
# edit: [input] kind = "raw", source = "ADC_raw_output", background = "ADC_raw_output_background"
gammaspec config check my_run.toml       # validate
gammaspec run my_run.toml                # convert -> template -> events -> spectrum (-> benchmark)
gammaspec run examples/run_2021_raw.toml # or a ready-made example
gammaspec run examples/run_simulation.toml   # no data needed: simulated source + background, scored against the truth
gammaspec stream sim out/sim --duration 20 --show   # live spectrum of a simulated stream
```

Or step by step:

```bash
# List the captures in a directory in time order; report gaps and live time
python -m gammaspec index ADC_raw_output --csv data/index_source.csv

# Convert the first 100 captures (in time order) to a compressed HDF5 store
python -m gammaspec convert ADC_raw_output data/source_first100.h5 --limit 100
```

```python
from gammaspec.store import CaptureStore

with CaptureStore("data/source_first100.h5") as store:
    print(len(store), "captures,", store.live_time_s, "s live time")
    volts = store.volts(0)          # float64 waveform of the first capture
    codes = store.codes(0)          # the exact int32 ADC codes
```

## Documentation

In reading order: each document with the code it describes. Steps 2–12 follow the
data through the pipeline, so they can be reviewed one by one, side by side with
their modules.

| # | Document | Contents | Code |
|---|---|---|---|
| 1 | [docs/architecture.md](docs/architecture.md) | Overview: pipeline design, module map, how the stages fit together | all (module map) |
| 2 | [docs/data-formats.md](docs/data-formats.md) | Step 1: raw capture format, HDF5 store schema, index CSV, event-file format | [rawio.py](gammaspec/rawio.py), [store.py](gammaspec/store.py), [events.py](gammaspec/events.py) |
| 3 | [docs/characterization.md](docs/characterization.md) | Step 2: measured pulse shape, pulse template, noise and preliminary spectra | [signal.py](gammaspec/signal.py), [template.py](gammaspec/template.py), [noise.py](gammaspec/noise.py), [characterize.py](gammaspec/characterize.py), [spectrum.py](gammaspec/spectrum.py), [plotting.py](gammaspec/plotting.py) |
| 4 | [docs/baseline-restoration.md](docs/baseline-restoration.md) | Step 3: interference removal, pole-zero correction, baseline tracker: method and results | [interference.py](gammaspec/interference.py), [baseline.py](gammaspec/baseline.py), [preprocess.py](gammaspec/preprocess.py), [restoration.py](gammaspec/restoration.py) |
| 5 | [docs/detection.md](docs/detection.md) | Step 4: event detection and the event list: method, evaluation, results | [detect.py](gammaspec/detect.py), [events.py](gammaspec/events.py), [pipeline.py](gammaspec/pipeline.py), [detection_eval.py](gammaspec/detection_eval.py) |
| 6 | [docs/pileup-benchmark.md](docs/pileup-benchmark.md) | Step 5: pile-up methods (TM, HYPER, LS/GLS/WLS) compared on synthetic and real pile-ups | [pileup.py](gammaspec/pileup.py), [pileup_eval.py](gammaspec/pileup_eval.py) |
| 7 | [docs/pileup-solver.md](docs/pileup-solver.md) | Step 6: the pile-up solver in the pipeline: method, tuning, results | [solver.py](gammaspec/solver.py), [pipeline.py](gammaspec/pipeline.py) (`_solve_clusters`, `spectrum_report`) |
| 8 | [docs/spectrum.md](docs/spectrum.md) | Step 7: energy calibration, gain matching, background subtraction; first-100 and full-run results | [calibration.py](gammaspec/calibration.py), [analysis.py](gammaspec/analysis.py), [spectrum.py](gammaspec/spectrum.py) |
| 9 | [docs/benchmark.md](docs/benchmark.md) | Step 8: gammaspec vs naive `find_peaks` vs an emulated analog chain vs typical analog NaI(Tl) systems | [benchmark.py](gammaspec/benchmark.py) |
| 10 | [docs/configuration.md](docs/configuration.md) | Step 9: the run file: input, stage switches, every parameter, validation, reruns | [config.py](gammaspec/config.py), [runner.py](gammaspec/runner.py), [examples/](examples/) |
| 11 | [docs/simulation.md](docs/simulation.md) | Step 10: the signal simulator, its physics, and the pipeline scored against the truth | [simulate.py](gammaspec/simulate.py) |
| 12 | [docs/streaming.md](docs/streaming.md) | Step 10: sources, overlapping frames, live input, live spectrum; validation and speed | [sources.py](gammaspec/sources.py), [stream.py](gammaspec/stream.py), [liveplot.py](gammaspec/liveplot.py) |
| 13 | [docs/cli.md](docs/cli.md) | Command-line reference (every command, all steps) | [cli.py](gammaspec/cli.py), [__main__.py](gammaspec/__main__.py) |
| 14 | [docs/walkthrough.ipynb](docs/walkthrough.ipynb) | Every function in call order, run on the real data, with inputs, outputs and the maths | all; built by [docs/build_walkthrough.py](docs/build_walkthrough.py) |
| 15 | [docs/development.md](docs/development.md) | Environment setup, testing, conventions | [tests/](tests/), [pyproject.toml](pyproject.toml) |
| 16 | [docs/adr/](docs/adr/) | Architecture Decision Records: why things are the way they are (0001–0012, in step order) | — |
| 17 | [CHANGELOG.md](CHANGELOG.md) | Release history | — |

## Results at a glance

Same 100 captures of the ⁶⁰Co run through three chains, each calibrated the same way
([details](docs/benchmark.md)):

| | naive `find_peaks` | analog CR-RC⁴ 0.5 µs + PUR (emulated) | **gammaspec** |
|---|---|---|---|
| FWHM at 1173 / 1332 keV | 8.27 / 7.44 % | 6.95 / 6.21 % | **6.54 / 6.05 %** |
| Valley-to-peak | 0.39 | 0.24 | **0.215** |
| Low-energy threshold | none (noise) | ≈ 44 keV | **≈ 6 keV** |
| Pile-up | merged / split | rejected | **resolved** |

Typical 3″×3″ NaI(Tl) with analog electronics: ≈ 7 % at 662 keV and ≈ 5 % at 1332 keV.
Those are other detectors, and resolution at 1.3 MeV is set mostly by the crystal and the
PMT.

## Roadmap

1. ✅ Ingestion: fast reader, lossless binary store, chronological index
2. ✅ Pulse template and noise characterisation ([results](docs/characterization.md))
3. ✅ Interference removal and baseline restoration: pole-zero correction plus a slow baseline tracker ([results](docs/baseline-restoration.md))
4. ✅ Event detection (optimal filter + successive cancellation) and the event list ([results](docs/detection.md))
5. ✅ Pile-up benchmark: template matching vs HYPER vs least-squares fits ([results](docs/pileup-benchmark.md))
6. ✅ Pile-up solver: weighted least squares with sub-sample timing and a search for hidden pulses ([results](docs/pileup-solver.md))
7. ✅ Energy calibration, gain matching and background subtraction ([results](docs/spectrum.md)); full 3-hour run processed
8. ✅ Benchmark: naive `find_peaks` count and an emulated preamp + shaping amp + MCA chain on the same data, plus typical published analog NaI(Tl) values ([results](docs/benchmark.md))
9. ✅ User-selectable stages: one run file for input, stages and parameters; `gammaspec run` ([docs](docs/configuration.md))
10. ✅ Streaming (overlapping frames, live sources, live spectrum) and a signal simulator with truth ([streaming](docs/streaming.md), [simulation](docs/simulation.md))
11. ⬜ UI (edits the run file; input: raw / store / stream / simulation) and packaging
