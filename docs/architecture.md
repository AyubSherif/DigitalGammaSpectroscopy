# Architecture

This document explains how `gammaspec` is organised and why. For file formats
see [data-formats.md](data-formats.md); for the reasoning behind individual
choices see the [ADRs](adr/).

## Signal chain and measured properties

```text
NaI(Tl) crystal → PMT → ADC (1.8 MS/s) → ASCII capture files → gammaspec
```

Properties measured from the data, which drive the design:

| Property | Value | Consequence |
|---|---|---|
| Polarity | negative-going pulses | the pipeline flips the sign internally |
| Rise time | ~1 sample (0.56 µs) | the peak is defined by 1–2 samples, so the peak-sample height depends on sampling phase; energies need a fit with sub-sample timing |
| Decay | τ ≈ 6.5 samples (3.6 µs), set by the PMT output circuit | the pulse shape is stable, which suits template methods |
| Undershoot | ≈ 3.4 % of the amplitude, recovers with τ ≈ 460 samples (≈ 260 µs) | the baseline shifts after each pulse; this needs baseline restoration |
| Source count rate | ≈ 4,100 pulses/s (⁶⁰Co run) | average gap between pulses ≈ 440 samples, so almost every pulse sits on an earlier pulse's undershoot; about 10 % overlap another pulse directly |
| Background count rate | ≈ 600 pulses/s | |
| Noise | σ ≈ 2.9 mV (source run), 3.9 mV (background run) | |
| ADC step (LSB) | 84.771697 µV | values are stored as int32 codes, losslessly |
| Acquisition | 0.29 s captures roughly every 2.3 s | live time is about 13 % of the wall-clock time; rates must be normalised by live time |

## Pipeline

```text
                ┌────────────┐   ┌──────────┐   ┌───────────┐   ┌──────────┐   ┌──────────┐   ┌────────────┐
 raw captures ─►│ ingestion  │──►│ baseline │──►│ detection │──►│ pile-up  │──►│ spectrum │──►│ background │
                │  (store)   │   │ restore  │   │ (events)  │   │ & energy │   │ & calib. │   │ subtraction│
                └────────────┘   └──────────┘   └───────────┘   └──────────┘   └──────────┘   └────────────┘
                   done            planned         planned         planned        planned        planned
```

The order differs from the "clean first" approach on purpose:

* **Baseline restoration comes before detection.** If you threshold the raw
  signal, you lose small pulses that start on an undershoot.
* **Detection is where the data shrinks.** Only an *event list* is kept (time,
  energy, fit quality, flags, optionally a short waveform snippet). This reduces
  the data from about 7 MB per capture to tens of kB.
* **Background subtraction is done on the spectrum, not on the waveform.**
  Background events happen at random times, so they can't be subtracted sample
  by sample. The background run goes through the identical pipeline, and its
  spectrum is scaled by live time and subtracted.

### Modularity

Every stage after ingestion is designed as an independent, optional step. A
pipeline configuration chooses which stages run and with what parameters (for
example, baseline restoration on or off). See
[ADR-0003](adr/0003-modular-optional-stages.md). The intended interface is:

```python
class Stage(Protocol):
    name: str
    def process(self, block: Block, state: StageState) -> Block: ...
```

Stages are stateful and work on chunks of data, so the same code can later
process a live stream, not just stored files.

## Module map

| Module | Responsibility |
|---|---|
| `gammaspec/rawio.py` | Parse ASCII capture files: header, fast fixed-width sample parser, volts → ADC code conversion with grid validation, ADC step estimation, chronological index |
| `gammaspec/store.py` | Convert captures to a compressed HDF5 store (parallel parsing); `CaptureStore` gives read access |
| `gammaspec/cli.py`, `__main__.py` | Command-line entry point: `python -m gammaspec …` |
| `tests/` | Unit tests (synthetic captures) and real-data checks that are skipped when the data is absent |

## Performance

On the development machine, ingestion runs at about 9.5 M samples/s with
parallel parsing, about 5× faster than the ADC records. The full 3-hour run
(4,786 captures) converts in about 5 minutes to about 3.4 GB. Later stages are
budgeted to stay faster than real time, so streaming will be possible.
