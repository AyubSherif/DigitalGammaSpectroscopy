# gammaspec

Digital pulse processing for gamma-ray spectroscopy: turns raw, high-rate ADC
waveforms from a scintillation detector into a clean, pile-up–corrected energy
spectrum, entirely in software.

Developed against 3 hours of data from a NaI(Tl) detector read out directly from
the photomultiplier (PMT) and digitised at 1.8 MS/s, with a ⁶⁰Co source and a
separate background run.

> **Status: early development (v0.1).** Data ingestion is complete. Detection,
> baseline restoration, pile-up correction and spectrum building are in progress.
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
  → baseline restoration         (planned)
  → event detection              (planned)  ← the data reduction step
  → pile-up resolution & energy  (planned)
  → histogram & calibration      (planned)
  → background subtraction       (planned)
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
| [docs/data-formats.md](docs/data-formats.md) | Raw capture format, HDF5 store schema, index CSV |
| [docs/cli.md](docs/cli.md) | Command-line reference |
| [docs/development.md](docs/development.md) | Environment setup, testing, conventions |
| [docs/adr/](docs/adr/) | Architecture Decision Records: why things are the way they are |
| [CHANGELOG.md](CHANGELOG.md) | Release history |

## Roadmap

1. ✅ Ingestion: fast reader, lossless binary store, chronological index
2. ⬜ Pulse template and noise characterisation
3. ⬜ Baseline restoration (pole-zero correction for the undershoot, plus a slow baseline tracker)
4. ⬜ Event detection (matched filter) and the event list
5. ⬜ Pile-up benchmark: template matching vs least-squares fit vs HYPER
6. ⬜ Pile-up solver and energy estimation
7. ⬜ Spectrum, energy calibration (⁶⁰Co, ⁴⁰K, ²⁰⁸Tl) and background subtraction
8. ⬜ User-selectable stages (pipeline configuration)
9. ⬜ Streaming input and live spectrum
