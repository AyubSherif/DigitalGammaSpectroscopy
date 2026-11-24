# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- **Renamed** from "gammaspec" (an existing, unrelated app) to the generic **Gamma
  Spectrometry App**:

  | | Old | New |
  |---|---|---|
  | Python package | `gammaspec` | `gamma_spectrometry` |
  | Command | `gammaspec` | `gamma-spectrometry` |
  | pip distribution | `gammaspec` | `gamma-spectrometry` (`pip install -e ".[app]"`) |
  | Version attribute in HDF5 files | `gammaspec_version` | `app_version` |
  | Jupyter kernel | `gammaspec` | `gamma-spectrometry` |

  Run files and data formats are otherwise unchanged. Entries below keep the old name, as
  it was at the time.

## [0.3.0] - 2026-11-24

### Added
- Step 11, the browser app `gamma_spectrometry.app` (Streamlit) and CLI command `app`:
  * input choice (raw files, store, replay, simulation, live file / FIFO) and stage switches
    as widgets, with path checks and live validation;
  * a TOML editor for every setting, and run-file load / download;
  * Run / Stop of `gamma-spectrometry run` as a separate process, live log and live spectrum;
  * results: steps, calibration, line rates, plots, truth check, benchmark.

  Listens on localhost by default. Optional dependency `gamma-spectrometry[app]`. `docs/app.md`,
  ADR-0013, `tests/test_app.py`.
- Walkthrough Part 0, "The app in pictures": one capture through every step, with 13
  before/after figures and plain-language explanations; Part K (the app).

### Changed
- README rewritten: what the software is, which analog instruments it replaces (shaping
  amplifier, pile-up rejector, MCA) and why (cost, portability, high-rate performance, data
  kept), results at a glance, background.
- The signal chain is described correctly as PMT → preamplifier → ADC → PC (not a raw PMT
  signal). The benchmark's analog emulation therefore uses the real preamplifier output.

## [0.2.0] - 2026-11-24

### Added
- Step 10, streaming: `gamma_spectrometry.sources` (store, raw files, simulation, binary pipe /
  stdin, real-time pacing) and `gamma_spectrometry.stream` (`StreamProcessor`: overlapping frames
  through the unchanged capture pipeline, parallel, exact int32 transfer, batched writes,
  live frame dropping as dead time; `LiveSpectrum` with automatic calibration;
  `process_stream`), `gamma_spectrometry.liveplot`. Run-file kinds `replay`, `simulation`, `pipe`,
  sections `[stream]` and `[simulation]`; CLI `stream`. ADR-0012, `docs/streaming.md`.
- Signal simulator `gamma_spectrometry.simulate`: measured pulse model, ⁶⁰Co cascade summing,
  Klein–Nishina Compton, backscatter, room background, resolution, noise (white or measured
  spectrum), interference, offset, polarity, ADC; continuous and independent of read sizes;
  simulated stores (`store.write_store`); truth files and `compare_to_truth`. CLI
  `simulate` (store, or samples on stdout). `docs/simulation.md`.
- Event files: capture column `live_start_s` (older files read as NaN); attribute `stream`.
- `examples/run_simulation.toml`, `examples/run_replay_realtime.toml`; walkthrough Part J.
- Step 9, run configuration: `gamma_spectrometry.config` (one TOML/JSON run file mirroring the config
  dataclasses; strict validation with suggestions; stage-combination checks; commented
  defaults from the code) and `gamma_spectrometry.runner` (input → template → events → spectrum →
  benchmark, reusing unchanged outputs via stamps). CLI `run`, `config init|check|show`, and
  `process --config/--set` (the old flags override the file). `examples/*.toml`,
  `docs/configuration.md`, ADR-0011; ADR-0003 implemented.
- `gamma_spectrometry.benchmark` and CLI command `benchmark`: the pipeline vs a naive
  `scipy.signal.find_peaks` count vs an emulated analog chain (CR-RC⁴ shaping amplifier with
  pole-zero, baseline restorer, peak-hold ADC, pile-up rejector, dead time, live-time clock)
  on the same captures, with common metrics (FWHM, valley-to-peak, peak rates, pile-up
  continuum, threshold, live fraction); `docs/benchmark.md`, ADR-0010, walkthrough Part H.
- Full 3-hour run results in `docs/spectrum.md` (4,786 captures, 6.38 M events, gain stable
  to 0.11 % rms, net ⁴⁰K consistent with zero).

### Changed
- `calibration.anchor_gain` searches between the median and the 99.5th percentile of the
  amplitudes, so the ⁶⁰Co pair is found even when many small (noise) triggers are present.
- `calibration.fit_lines` drops lines whose centre or width is at a bound before the least
  significant ones, so a real line is not dropped in favour of a degenerate neighbour.

## [0.1.0] - 2026-11-24

The offline pipeline, steps 1–7.

### Fixed
- `PoleZero` started from the first sample's value, which the filter's DC gain turned
  into a slowly decaying offset of up to ~0.2 V (hidden by the baseline tracker). It now
  starts at the input's mean level (0 after offset removal).

### Changed
- `PulseTemplate` evaluation by direct indexing on its uniform grid (identical results,
  much faster); pipeline workers limited to one BLAS thread each (they oversubscribed the
  CPU).
- Faster processing: noise covariance and line detection with fewer and faster FFTs,
  cumulative sums in `edge_signal`, `noise_sigma` on the central 131k samples.

### Added
- `gamma_spectrometry.calibration`: automatic energy calibration (anchor, bounded line fits with
  iterative removal of rejected lines, linear/quadratic curve), gain matching to a
  reference run, gain drift vs time (ADR-0009).
- `gamma_spectrometry.analysis` and CLI command `spectrum`: calibrated spectra, live-time
  normalisation, background subtraction with errors, net line rates in local windows,
  CSV/JSON/plots; `docs/spectrum.md`.
- `spectrum.fit_peaks`: optional bounds (centre, width, non-negative area) and area errors.
- `gamma_spectrometry.solver`: `PileupSolver` (WLS + sub-sample timing with a residual search for
  hidden pulses), run by the pipeline on every pile-up cluster and on badly fitted
  isolated pulses; new flag `FIT_ADDED`; `EventFile.good()`; CLI options
  `--no-pileup-solver` / `--no-residual-search`; CLI command `report` (peaks by category,
  throughput); `docs/pileup-solver.md`.
- `PulseTemplate.derivative()`.
- `gamma_spectrometry.pileup`: pile-up resolution methods (TM as in the reference paper, HYPER,
  least squares with integer times, and LS / GLS / WLS with sub-sample timing).
- `gamma_spectrometry.pileup_eval` and CLI command `evaluate-pileup`: benchmark on template pairs in
  real noise, real-pulse pairs and the real pile-up events (ADR-0008).
- `docs/pileup-benchmark.md`; walkthrough Part E (pile-up) and Part F (evaluation tools).
- `gamma_spectrometry.detect`: `Detector` with a per-capture optimal (noise-whitened) filter,
  successive cancellation with residual re-search, local thresholds estimated from the
  residual; matched, sharpened and edge alternatives (ADR-0007).
- `gamma_spectrometry.events`: event-file format (`EventWriter`, `EventFile`, flag bits).
- `gamma_spectrometry.pipeline` and CLI command `process`: store → pre-processing → detection →
  event file, in parallel over captures, with clusters, flags, template fits and pile-up
  snippets.
- `gamma_spectrometry.detection_eval` and CLI command `evaluate-detection`.
- `docs/detection.md`; `docs/walkthrough.ipynb` (generated by `docs/build_walkthrough.py`).
- `gamma_spectrometry.interference`: `LineRemoval` stage that detects narrow interference lines and
  subtracts them as block-wise sinusoids fitted to pulse-free samples (ADR-0006).
- `gamma_spectrometry.baseline`: pole-zero correction designed from the pulse model
  (`design_pole_zero`, `PoleZero`, `filter_template`) and a `BaselineTracker` stage (ADR-0005).
- `gamma_spectrometry.preprocess`: `Preprocessor` chaining the optional stages from a
  `PreprocessConfig`; warns when pole-zero is used without the tracker.
- `gamma_spectrometry.spectrum`: Gaussian peak fits for resolution measurements.
- `gamma_spectrometry.restoration` and CLI command `evaluate-restoration`: compares stage
  combinations on the same pulses (noise, pre-pulse baseline, tail, FWHM, pulses on tails).
- `fit_pulses(..., fit_baseline=False)` for restored signals.
- `docs/baseline-restoration.md`: method and results.
- `gamma_spectrometry.signal`: polarity handling, robust noise estimates, coarse leading-edge
  trigger, isolation and pulse-free masks.
- `gamma_spectrometry.template`: fine-grid `PulseTemplate`; per-pulse fit of amplitude,
  baseline and sub-sample time; iterative template construction with sampling-phase
  equalisation (ADR-0004); analytic multi-exponential pulse model.
- `gamma_spectrometry.noise`: noise power spectrum and autocorrelation from pulse-free windows.
- `gamma_spectrometry.characterize` and the CLI command `characterize`: template, model,
  shape-vs-amplitude check, noise with quiet/noisy split and interference lines,
  stability over time, preliminary spectrum, diagnostic plots.
- `docs/characterization.md`: results for the 2021 NaI(Tl) data.
- `gamma_spectrometry.rawio`: reader for ASCII ADC capture files. It has a vectorised
  fixed-width parser (about 7 M samples/s on one core) and a generic fallback,
  converts volts to ADC codes losslessly with grid validation, estimates the
  ADC step, and builds a chronological index.
- `gamma_spectrometry.store`: parallel conversion of captures to a compressed HDF5 store
  (int32 codes, gzip+shuffle, about 10× smaller than text), and a
  `CaptureStore` reader.
- CLI commands `index` and `convert` (`python -m gamma_spectrometry`).
- Test suite with synthetic captures and real-data checks.
- Documentation: architecture, data formats, CLI reference, development guide,
  ADRs 0001–0009.
