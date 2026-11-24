# Architecture

This document explains how `gammaspec` is organised and why. For file formats
see [data-formats.md](data-formats.md); for the reasoning behind individual
choices see the [ADRs](adr/).

## Signal chain and measured properties

```text
NaI(Tl) crystal → PMT → ADC (1.8 MS/s) → ASCII capture files → gammaspec
```

Properties measured from the data, which drive the design (details and figures in
[characterization.md](characterization.md)):

| Property | Value | Consequence |
|---|---|---|
| Polarity | negative-going pulses | the pipeline flips the sign internally |
| Rise | τ = 0.245 µs (the NaI(Tl) light decay); about 1 sample | the peak is defined by 1–2 samples, so the peak-sample height depends on sampling phase; energies need a fit with sub-sample timing |
| Decay | τ = 4.1 µs (7.4 samples), set by the PMT output circuit | the pulse shape is stable and independent of amplitude (< 1 %), which suits template methods |
| Undershoot | −4.3 % of the peak at +26 µs; τ = 42 µs and 176 µs; net pulse area ≈ 0 | it comes from AC coupling, so pole-zero correction can remove it; this is the basis of baseline restoration |
| Source count rate | ≈ 4,300 pulses/s (⁶⁰Co run) | average gap between pulses ≈ 420 samples, so almost every pulse sits on an earlier pulse's undershoot; ≈ 21 % lack a clean main lobe |
| Background count rate | ≈ 740 pulses/s (quiet captures) | |
| Noise | white floor ≈ 0.4–0.5 mV; total 1.4–6.4 mV | dominated by intermittent narrow-band interference (257, 380, 847 kHz …), which causes false triggers |
| Gain | ≈ 0.265 mV/keV (source run), ≈ 6 % lower in the background run | gain matching needed before background subtraction |
| ADC step (LSB) | 84.771697 µV | values are stored as int32 codes, losslessly |
| Acquisition | 0.29 s captures roughly every 2.3 s | live time is about 13 % of the wall-clock time; rates must be normalised by live time |

## Pipeline

```text
                ┌────────────┐   ┌──────────┐   ┌───────────┐   ┌──────────┐   ┌──────────┐   ┌────────────┐
 raw captures ─►│ ingestion  │──►│ baseline │──►│ detection │──►│ pile-up  │──►│ spectrum │──►│ background │
                │  (store)   │   │ restore  │   │ (events)  │   │ & energy │   │ & calib. │   │ subtraction│
                └────────────┘   └──────────┘   └───────────┘   └──────────┘   └──────────┘   └────────────┘
                   done             done            done            done           done           done

 baseline restore = [line removal] → [pole-zero] → [baseline tracker], see baseline-restoration.md
```

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
| `gammaspec/signal.py` | Polarity, robust noise estimates, coarse leading-edge trigger, isolation and pulse-free masks |
| `gammaspec/template.py` | `PulseTemplate` (fine-grid, sub-sample evaluation, HDF5 I/O); per-pulse fit of amplitude, baseline and sub-sample time; iterative template construction with phase equalisation; analytic multi-exponential model |
| `gammaspec/noise.py` | Pulse-free noise windows, power spectral density, autocorrelation |
| `gammaspec/characterize.py` | Detector characterisation: template, model, shape vs amplitude, noise (quiet/noisy split, interference lines), stability, preliminary spectrum, plots |
| `gammaspec/interference.py` | `LineRemoval` stage: line detection in the masked spectrum, block-wise sinusoid subtraction on pulse-free samples |
| `gammaspec/baseline.py` | Pole-zero filter design from the pulse model, `PoleZero` stage (streamable IIR), `filter_template`, `BaselineTracker` stage |
| `gammaspec/preprocess.py` | `PreprocessConfig` / `Preprocessor`: the optional chain polarity+offset → lines → pole-zero → tracker |
| `gammaspec/detect.py` | `Detector`: per-capture optimal filter (noise-whitened matched filter), successive cancellation, local thresholds; matched / sharpened / edge alternatives |
| `gammaspec/events.py` | Event-file format: `EventWriter` (append per capture), `EventFile` (read), flag bits |
| `gammaspec/pipeline.py` | `Pipeline` / `run`: store → pre-processing → detection → clusters, flags, template fits, snippets → event file (parallel over captures) |
| `gammaspec/detection_eval.py` | Injection tests of the detection methods in real and phase-scrambled noise |
| `gammaspec/pileup.py` | Pile-up resolution methods for one cluster: TM (paper), HYPER, LS, LS / GLS / WLS with sub-sample timing |
| `gammaspec/solver.py` | `PileupSolver`: WLS + timing fit of a cluster, residual search for hidden pulses (step 6) |
| `gammaspec/pileup_eval.py` | Pile-up benchmark: template pairs in real noise, real-pulse pairs, real pile-up events |
| `gammaspec/calibration.py` | Energy calibration: anchor, bounded line fits, calibration curve, gain matching, gain vs time |
| `gammaspec/analysis.py` | Energy spectra, live-time normalisation, background subtraction, net line rates, outputs (step 7) |
| `gammaspec/spectrum.py` | Gaussian peak fits on a linear background (resolution, centroids) |
| `gammaspec/benchmark.py` | Benchmark chains on the same captures: naive `find_peaks`, emulated analog CR-RC⁴ chain (BLR, peak-hold, PUR, dead time, live-time clock), gammaspec; common metrics |
| `gammaspec/restoration.py` | Evaluation of pre-processing configurations on real data (same pulses for all), plots |
| `gammaspec/plotting.py` | Shared figure style (validated colour palette, recessive axes) |
| `gammaspec/config.py` | Run file: `RunConfig` (input, template, pipeline, spectrum, benchmark), TOML/JSON load and save, strict validation, stage-combination checks |
| `gammaspec/runner.py` | `run(RunConfig)`: input → template → events → spectrum → benchmark, reusing outputs whose inputs and settings are unchanged |
| `gammaspec/simulate.py` | Signal simulator: pulse model, ⁶⁰Co/¹³⁷Cs cascades with summing, Klein–Nishina Compton, room background, noise, interference, ADC; continuous `read(n)`; simulated stores; truth and `compare_to_truth` |
| `gammaspec/sources.py` | Signal sources yielding chunks: `StoreSource`, `RawSource`, `SimulationSource`, `PipeSource` (stdin / FIFO), `Paced` (real time) |
| `gammaspec/stream.py` | `StreamProcessor`: overlapping frames through the capture pipeline (parallel, live frame dropping), `LiveSpectrum`, `process_stream` |
| `gammaspec/liveplot.py` | Live spectrum window for `gammaspec stream --show` |
| `gammaspec/cli.py`, `__main__.py` | Command-line entry point: `python -m gammaspec …` |
| `tests/` | Unit tests (synthetic captures) and real-data checks that are skipped when the data is absent |

## Performance

On the development machine, ingestion runs at about 9.5 M samples/s with
parallel parsing, about 5× faster than the ADC records. The full 3-hour run
(4,786 captures) converts in about 5 minutes to about 3.4 GB. Later stages are
budgeted to stay faster than real time, so streaming will be possible.
