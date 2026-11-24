# Configuration: choosing the input, the stages and their parameters (step 9)

A single **run file** (TOML) describes a whole run:
* where the data comes from;
* which processing stages are applied;
* every parameter;
* which analysis steps follow.

`gammaspec run` executes it. No code is needed to switch a correction off or change a
threshold, and the same file reproduces the run later.

Code: `gammaspec/config.py` (schema, loading, validation), `gammaspec/runner.py` (runs the
steps). Decision records: [ADR-0003](adr/0003-modular-optional-stages.md),
[ADR-0011](adr/0011-run-configuration.md).

```bash
gammaspec config init my_run.toml          # every setting at its default, with a comment on each
gammaspec config check my_run.toml         # validate; show warnings and the effective configuration
gammaspec run my_run.toml                  # run it
gammaspec run my_run.toml --set pipeline.detect.threshold=6 --set input.limit=20   # one-off changes
gammaspec config show results/events/source.h5   # the settings an event file was made with
```

Examples: [examples/run_2021_raw.toml](../examples/run_2021_raw.toml) (the whole chain from
the raw files) and [examples/run_from_stores.toml](../examples/run_from_stores.toml)
(existing stores and template, with some stages switched off).

## A minimal run file

Keys that are left out keep their defaults:

```toml
output = "data/runs/2021"

[input]
kind = "raw"                              # start from the ASCII captures
source = "ADC_raw_output"
background = "ADC_raw_output_background"
```

This converts both runs, builds the pulse template from the background run, processes
both runs with every stage on, then calibrates and subtracts the background.

## What the run does

| Step | Runs when | Output (in `output/`) |
|---|---|---|
| 1. **Input** | `input.kind = "raw"`: convert the ASCII captures to an HDF5 store. `"store"`: use the given `.h5` files as they are. | `stores/source.h5`, `stores/background.h5` |
| 2. **Template** | `template.path` empty: build one with the characterisation step from `template.build_from` (`[template.characterize]` settings). Otherwise use the given file. | `template/` |
| 3. **Process** | Always, for the source run and (if `input.background` is set) the background run: pre-processing → detection → pile-up solver (`[pipeline.*]`). | `events/source.h5`, `events/background.h5` |
| 4. **Spectrum** | `spectrum.enabled` (default on): calibration, spectra, background subtraction if there is a background run. | `spectrum/` |
| 5. **Benchmark** | `benchmark.enabled` (default off): naive / analog / gammaspec comparison. | `benchmark/` |

Every run also writes `config.toml`, the **effective** configuration (the file plus the
`--set` changes plus the defaults), and `run_summary.json` (time and output of each step).
The event files store their pipeline settings too.

**Reruns only redo what changed.** Each output has a `.stamp.json` next to it that records
its inputs and settings. When they are unchanged, the output is reused. For example, after
changing `pipeline.detect.threshold`, conversion and template are reused and only processing,
spectrum and benchmark run again. `overwrite = true` recomputes everything.

## The stage switches

| Setting | Default | What it switches |
|---|---|---|
| `input.background` | `""` | Background run; empty = no background subtraction |
| `pipeline.preprocess.line_removal` | `true` | Interference line removal ([baseline-restoration.md](baseline-restoration.md)) |
| `pipeline.preprocess.pole_zero` | `true` | Pole-zero correction (needs `baseline_tracker`) |
| `pipeline.preprocess.baseline_tracker` | `true` | Slow baseline subtraction |
| `pipeline.detect.method` | `"optimal"` | Detection filter: `optimal`, `matched`, `sharpened`, `edge` ([detection.md](detection.md)) |
| `pipeline.detect.cancel` | `true` | Successive cancellation (not used by `edge`) |
| `pipeline.solver.enabled` | `true` | Pile-up solver (needs `pole_zero` + `baseline_tracker`) ([pileup-solver.md](pileup-solver.md)) |
| `pipeline.solver.residual_search` | `true` | Search for hidden pulses in pile-up |
| `pipeline.events.snippets` | `"pileup"` | Which clusters keep waveform snippets: `none`, `pileup`, `all` |
| `spectrum.enabled` | `true` | Calibration, spectra, background subtraction ([spectrum.md](spectrum.md)) |
| `benchmark.enabled` | `false` | Benchmark ([benchmark.md](benchmark.md)) |

The order of the stages is fixed, because each needs the previous one's output: for example,
detection needs a restored baseline to set its threshold. Only their use is optional.

## All sections

| Section | Code (dataclass) | Contents |
|---|---|---|
| top level | `config.RunConfig` | `output`, `workers`, `overwrite` |
| `[input]` | `config.InputConfig` | `kind`, `source`, `background`, `limit`, `pattern`, `reuse_converted` |
| `[template]`, `[template.characterize]` | `config.TemplateConfig`, `characterize.CharacterizeConfig` | template file, or how to build it |
| `[pipeline.preprocess]` (`.lines`, `.tracker`) | `preprocess.PreprocessConfig`, `interference.LineRemovalConfig`, `baseline.TrackerConfig` | pre-processing |
| `[pipeline.detect]` | `detect.DetectConfig` | detection |
| `[pipeline.events]` | `pipeline.EventConfig` | clusters, fit windows, flags, snippets |
| `[pipeline.solver]` | `solver.SolverConfig` | pile-up solver |
| `[spectrum]` | `config.SpectrumStepConfig` | bins, range |
| `[benchmark]` | `config.BenchmarkStepConfig` | on/off, number of captures |

`gammaspec config init` writes all of them, each line with its meaning. The tables mirror
the dataclasses in the code, so a new parameter in the code appears in the file
automatically.

## Validation

Loading is strict, and **all** problems are reported at once:

```text
invalid configuration:
  input.limt: unknown setting (did you mean 'limit'?)
  pipeline.detect.threshold: expected a number, got 'high'
  input.kind: 'raws' is not one of 'raw', 'store'
  pipeline.solver.search_threshold: must be > 0, got -1.0
```

* **Types:** true/false, integers, numbers (an integer is accepted for a number), strings,
  lists.
* **Choices:** for example `input.kind`, `detect.method`, `events.snippets`, `polarity`
  (−1 or 1).
* **Ranges:** numbers must be ≥ 0, and thresholds, window lengths and bin widths must be > 0.
* **Combinations** (`config check`, and before every run):

  | Case | Result |
  |---|---|
  | No `input.source` | error |
  | Paths that do not exist, or are a folder vs a file of the wrong kind | error |
  | A template to be built from a background run that is not given | error |
  | `pole_zero` without `baseline_tracker` | warning |
  | The solver without a restored baseline | warning (the solver is skipped) |
  | `cancel` with `edge` | warning |
  | A spectrum without a background run | warning |

A `--set KEY=VALUE` value is read as TOML (`6`, `1e-3`, `true`, `"edge"`, `[1, 2]`). A bare
word is taken as a string. The `process` command takes `--config FILE` and `--set` as well.
It uses the `[pipeline]` tables of a run file, and its older flags (`--no-pole-zero`, …)
still work as shortcuts that override the file.

## Limits

* With few captures (`input.limit`), a run can have too few counts to calibrate. The
  background run needs ~40 captures (12 s) for ⁴⁰K. The run then stops with a message saying
  so; the outputs of earlier steps are kept and reused next time.
* Choosing a pile-up method other than the WLS solver (TM, HYPER, …) is possible in the
  benchmark only (`gammaspec evaluate-pileup`), not in the pipeline. These methods
  were measured to be worse ([pileup-benchmark.md](pileup-benchmark.md)).
* Live streams and simulated input come in step 10. The UI (step 11) will edit the same
  run file.
