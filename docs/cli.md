# Command-line reference

```bash
python -m gammaspec <command> [options]
# or, after `pip install -e .`:
gammaspec <command> [options]
```

## `index`

List the capture files in a directory in time order, and report the wall-clock
span, live time, gaps between captures, and how often the file numbers are out
of time order.

```bash
gammaspec index DIRECTORY [--pattern GLOB] [--csv OUT.csv]
```

| Option | Default | Description |
|---|---|---|
| `DIRECTORY` | (required) | folder containing the capture files |
| `--pattern` | `raw_data_*` | glob pattern for capture files |
| `--csv` | none | also write the index to this CSV file ([format](data-formats.md#index-csv)) |

Only the headers are read, so this is fast (about 2 s for 4,786 files).

## `convert`

Convert capture files to a compressed HDF5 store ([format](data-formats.md#hdf5-capture-store-h5)).

```bash
gammaspec convert DIRECTORY OUT.h5 [--limit N] [--workers W] [--lsb VOLTS | --estimate-lsb] [--overwrite]
```

| Option | Default | Description |
|---|---|---|
| `DIRECTORY` | (required) | folder containing the capture files |
| `OUT.h5` | (required) | output file; parent folders are created |
| `--pattern` | `raw_data_*` | glob pattern for capture files |
| `--limit N` | all | convert only the first *N* captures **in time order** |
| `--workers W` | CPU cores − 1, at most 8 | parallel parser processes (`1` = no parallelism) |
| `--lsb VOLTS` | `8.4771697e-05` | ADC step used to convert volts to integer codes |
| `--estimate-lsb` | off | estimate the ADC step from the first capture (use for unknown hardware) |
| `--overwrite` | off | replace an existing output file |

The conversion stops with an error if any value is not on the ADC grid (more
than 0.02 LSB off), if the header's sample count is wrong, or if captures have
different lengths or sampling intervals.

Example output (convert):

```text
captures         : 100
first / last     : 2021-08-26 13:34:48.617564  ->  2021-08-26 13:38:32.695494
wall-clock span  : 0.062 h
live time        : 29.1 s  (13.0% duty)
...
wrote data\source_first100.h5: 100 captures, 70.7 MB (0.71 MB/capture), lsb=8.4771697000e-05 V, live time 29.1 s
```

## `characterize`

Build the pulse template and characterise the noise, stability and raw spectrum of a
store. Results are described in [characterization.md](characterization.md).

```bash
gammaspec characterize STORE.h5 OUT_DIR [--limit N] [--template-post S] [--isolation-pre S]
                       [--isolation-post S] [--trigger-threshold K] [--min-amplitude-sigma K]
                       [--max-amplitude-v V] [--resolution R] [--compare TEMPLATE.h5]
```

| Option | Default | Description |
|---|---|---|
| `STORE.h5` | (required) | capture store from `convert` |
| `OUT_DIR` | (required) | output folder (created) |
| `--limit N` | all | only the first *N* captures |
| `--template-post S` | 1500 | samples of template kept after the peak |
| `--isolation-pre S`, `--isolation-post S` | 1500 | template pulses need no other pulse this close before / after. **Use about 300 for high-rate data** |
| `--trigger-threshold K` | 5 | coarse trigger threshold, in noise σ of the edge signal |
| `--min-amplitude-sigma K` | 30 | template pulses must exceed K × noise σ |
| `--max-amplitude-v V` | 4 | template pulses must stay below this (avoids saturation) |
| `--resolution R` | 16 | template grid points per sample |
| `--compare TEMPLATE.h5` | none | overlay another template (e.g. from another run) |

Outputs in `OUT_DIR`:

| File | Contents |
|---|---|
| `template.h5` | the template (`shape`, `resolution`, `t_start`) with model parameters and provenance as attributes; load with `PulseTemplate.load` |
| `noise.h5` | `freq_hz`, `psd` (V²/Hz, one-sided), `autocorr`; σ and white level as attributes |
| `summary.json` | every number in the report, plus the configuration used |
| `template.png`, `shape_vs_amplitude.png`, `noise.png`, `stability.png`, `spectrum_raw.png` | diagnostic plots |

## `evaluate-restoration`

Compare combinations of the pre-processing stages (line removal, pole-zero correction,
baseline tracker) on real data. Every configuration measures the same pulses. Results
are described in [baseline-restoration.md](baseline-restoration.md).

```bash
gammaspec evaluate-restoration STORE.h5 TEMPLATE.h5 OUT_DIR [--n-peaks N] [--limit N]
```

| Option | Default | Description |
|---|---|---|
| `STORE.h5` | (required) | capture store |
| `TEMPLATE.h5` | (required) | template with a fitted model (from `characterize`, ideally a low-rate run) |
| `OUT_DIR` | (required) | output folder |
| `--n-peaks N` | 2 | spectral peaks to fit (2 for ⁶⁰Co, 1 for ⁴⁰K) |
| `--limit N` | all | only the first *N* captures |

Prints a table per configuration (noise σ, pre-pulse baseline spread, tail, FWHM with
free and fixed baseline) and a comparison of pulses on an earlier pulse's tail against
clean pulses. Writes `restoration_summary.json` and the plots `waveform_example.png`,
`line_removal_psd.png`, `pole_zero_template.png`, `baseline_and_tail.png` and
`spectrum_comparison.png`.

## `run`

Run the chain described by a run file ([configuration.md](configuration.md)): input
(convert raw captures, or use stores), template, events for the source and background runs,
spectrum and, if enabled, the benchmark.

```bash
gammaspec run RUN.toml [--set KEY=VALUE ...]
```

| Option | Description |
|---|---|
| `--set KEY=VALUE` | override one setting, e.g. `--set pipeline.detect.threshold=6`, `--set input.limit=20` (repeatable; value in TOML syntax) |

Writes `config.toml` (effective configuration), `run_summary.json` and the step outputs to
the run's `output` folder. Outputs whose inputs and settings did not change are reused.
Exit code 2 on an invalid configuration or a step that cannot run.

## `config`

```bash
gammaspec config init RUN.toml [--set KEY=VALUE ...] [--force]   # write all settings (defaults), commented
gammaspec config check RUN.toml [--set KEY=VALUE ...] [--quiet]  # validate; print warnings + effective config
gammaspec config show EVENTS.h5                                  # pipeline settings stored in an event file
```

## `process`

Run the pipeline: capture store → pre-processing → detection → event file
([format](data-formats.md#event-file-h5)). A spectrum plot of the clean events is written
next to the output (`OUT.png`).

```bash
gammaspec process STORE.h5 OUT.h5 --template TEMPLATE.h5 [--no-line-removal] [--no-pole-zero]
                  [--no-tracker] [--method optimal|matched|sharpened|edge] [--threshold K]
                  [--snippets pileup|all|none] [--limit N] [--overwrite]
                  [--config RUN.toml] [--set KEY=VALUE ...]
```

| Option | Default | Description |
|---|---|---|
| `--template` | (required) | template with a fitted model (from `characterize`) |
| `--no-line-removal`, `--no-pole-zero`, `--no-tracker` | all on | switch pre-processing stages off |
| `--method` | `optimal` | detection filter ([detection.md](detection.md)) |
| `--threshold K` | 5.5 | detection threshold, in noise sigma of the filter output |
| `--no-pileup-solver` | off | do not resolve pile-up clusters ([pileup-solver.md](pileup-solver.md)) |
| `--no-residual-search` | off | resolve clusters but do not search for hidden pulses |
| `--snippets` | `pileup` | store clean waveforms of pile-up clusters, of all clusters, or none |
| `--limit N` | all | only the first *N* captures |
| `--config FILE` | none | run file (or pipeline-only file); its `[pipeline]` tables are used ([configuration.md](configuration.md)) |
| `--set KEY=VALUE` | none | override one pipeline setting, e.g. `--set detect.threshold=6` (repeatable) |

The flags above are shortcuts that override `--config` / `--set` when given; defaults come from
the configuration.

Prints the event count, processing speed (× real time), rate, flag fractions and the size
reduction.

## `evaluate-detection`

Compare detection methods by injecting template pulses into real noise and into
phase-scrambled noise (same spectrum, no pulses).

```bash
gammaspec evaluate-detection STORE.h5 TEMPLATE.h5 OUT_DIR
```

Writes `detection_summary.json`, `detection_efficiency.png` and `detection_pairs.png`.

## `evaluate-pileup`

Benchmark the pile-up methods ([pileup-benchmark.md](pileup-benchmark.md)).

```bash
gammaspec evaluate-pileup STORE.h5 TEMPLATE.h5 EVENTS.h5 OUT_DIR
```

`EVENTS.h5` is the event file of the same store (its pile-up snippets are test C). Prints the
three test tables and the speed; writes `pileup_summary.json`, `pileup_template_pairs.png`,
`pileup_real_pairs.png` and `pileup_real_events.png`.

## `report`

Fit the main spectral peaks of an event file by category (isolated, resolved pile-up by
gap, all good events) and report the throughput.

```bash
gammaspec report EVENTS.h5 [--n-peaks N] [--peaks C1,C2]
```

| Option | Default | Description |
|---|---|---|
| `--n-peaks N` | 2 | number of peaks, found automatically in the isolated-pulse spectrum |
| `--peaks` | auto | starting guesses for the peak centres (V) |

Categories with fewer than 1000 events in the peak region are not fitted. Writes
`EVENTS.report.json` and the spectrum plot `EVENTS.png`.

## `spectrum`

Calibrate the runs, build energy spectra and subtract the background ([spectrum.md](spectrum.md)).

```bash
gammaspec spectrum SOURCE_EVENTS.h5 OUT_DIR [--background BACKGROUND_EVENTS.h5] [--bin-kev 5] [--emax 3000]
```

| Option | Default | Description |
|---|---|---|
| `--background` | none | event file of the background run (calibrated, or gain-matched to the source run) |
| `--bin-kev` | 5 | bin width (keV) |
| `--emax` | 3000 | upper end of the spectrum (keV) |

Prints the calibrations and their residuals, the ⁶⁰Co and ⁴⁰K rates in the source,
background and net spectra, total rates, and the gain drift. Writes
`calibration_source.json`, `calibration_background.json`, `spectrum.csv`,
`spectrum_summary.json`, `calibration.png`, `spectrum.png` and `spectrum_zoom.png`.

## `benchmark`

Compare gammaspec with a naive `scipy.signal.find_peaks` count and an emulated analog
chain (CR-RC⁴ shaping amplifier, baseline restorer, peak-hold ADC, pile-up rejector, dead
time) on the same captures ([benchmark.md](benchmark.md)).

```bash
gammaspec benchmark STORE.h5 EVENTS.h5 OUT_DIR [--first 0] [--n 100]
```

| Option | Default | Description |
|---|---|---|
| `--first` | 0 | first capture |
| `--n` | 100 | number of captures (0 = all; the analog emulation runs at about real time) |

`EVENTS.h5` must be the event file of `STORE.h5`. Prints FWHM (% and keV), valley-to-peak,
peak rates, pile-up continuum and live fraction for each chain. Writes
`benchmark_summary.json` and `benchmark_spectra.png`.
