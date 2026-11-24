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
> chain. Next: user-selectable stages, streaming.
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

Each processing stage will be optional and configurable, so users can switch
individual corrections on or off (see
[ADR-0003](docs/adr/0003-modular-optional-stages.md)).

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

| Document | Contents |
|---|---|
| [docs/architecture.md](docs/architecture.md) | Pipeline design, module map, how the stages fit together |
| [docs/characterization.md](docs/characterization.md) | Measured pulse shape, noise and preliminary spectra of the 2021 data |
| [docs/baseline-restoration.md](docs/baseline-restoration.md) | Interference removal, pole-zero correction, baseline tracker: method and results |
| [docs/detection.md](docs/detection.md) | Event detection and the event list: method, evaluation, results |
| [docs/pileup-benchmark.md](docs/pileup-benchmark.md) | Pile-up methods compared on synthetic and real pile-ups |
| [docs/pileup-solver.md](docs/pileup-solver.md) | The pile-up solver in the pipeline: method, tuning, results |
| [docs/spectrum.md](docs/spectrum.md) | Energy calibration, gain matching, background subtraction; first-100 and full-run results |
| [docs/benchmark.md](docs/benchmark.md) | gammaspec vs naive `find_peaks` vs an emulated analog chain vs typical analog NaI(Tl) systems |
| [docs/walkthrough.ipynb](docs/walkthrough.ipynb) | Every function in call order, run on the real data, with inputs, outputs and the maths |
| [docs/data-formats.md](docs/data-formats.md) | Raw capture format, HDF5 store schema, index CSV |
| [docs/cli.md](docs/cli.md) | Command-line reference |
| [docs/development.md](docs/development.md) | Environment setup, testing, conventions |
| [docs/adr/](docs/adr/) | Architecture Decision Records: why things are the way they are |
| [CHANGELOG.md](CHANGELOG.md) | Release history |

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
9. ⬜ User-selectable stages (pipeline configuration)
10. ⬜ Streaming input and live spectrum
