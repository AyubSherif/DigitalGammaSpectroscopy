# Data formats

## Raw capture files (input)

One file per ADC capture, named `raw_data_<N>.txt`. ASCII, `\n` line endings.

```text
Creation: 2021-08-26 13:34:48.617564      ← acquisition PC wall-clock time (local, no timezone)
Signals:  524288                          ← number of samples
Dt:       5.555555555555556e-04 ms        ← sample interval (units: s, ms, us, ns)
+3.390868e-04                             ← one sample per line, volts
-8.477170e-05
...
```

* Samples are printed in the fixed-width format `%+.6e` (13 characters plus a
  newline). `rawio` parses this format with a fast vectorised parser, and falls
  back to a generic float parser for anything else.
* Every value is an integer ADC code multiplied by the ADC step
  **LSB = 8.4771697 × 10⁻⁵ V**, printed to 7 significant digits. Each code maps
  to exactly one printed string, so converting values to integer codes loses
  nothing. Rebuilding volts as `code × LSB` reproduces the text to within 5 × 10⁻⁷
  relative. Observed codes run from about −65,200 to +3,800 (about −5.5 V to
  +0.32 V).
* **The file number does not give the time order.** In the 2021 data the
  numbers are scrambled relative to the `Creation` timestamps. Always order
  captures by timestamp (`build_index` does this).
* There are gaps between captures (about 2.3 s apart for 0.29 s of data). Live
  time is the sum of capture durations, not the wall-clock span.

### 2021 datasets

| Directory | Captures | Wall-clock span | Live time | Notes |
|---|---|---|---|---|
| `ADC_raw_output/` | 4,786 | 3.000 h (2021-08-26 13:34–16:34) | 1,394 s | ⁶⁰Co source |
| `ADC_raw_output_background/` | 100 | 5 min (2021-08-31) | 29.1 s | room background |

## HDF5 capture store (`*.h5`)

Written by `gammaspec convert` / `store.convert`; read with `store.CaptureStore`.

| Path | Type | Shape | Description |
|---|---|---|---|
| `/codes` | int32 | (n_captures, n_samples) | ADC codes, one row per capture, **in time order**. Chunked per row, gzip level 1 + shuffle (≈ 0.7 MB per capture). |
| `/index/time_s` | float64 | (n_captures,) | creation time as seconds since 1970, with the naive local time treated as if it were UTC |
| `/index/created` | UTF-8 string | (n_captures,) | creation timestamp from the header, ISO 8601 |
| `/index/file_number` | int32 | (n_captures,) | the number in the original file name |
| `/index/source` | UTF-8 string | (n_captures,) | path of the original file |

Root attributes:

| Attribute | Description |
|---|---|
| `lsb_volts` | ADC step; `volts = codes * lsb_volts` |
| `dt_s`, `sample_rate_hz` | sample interval and rate |
| `n_samples` | samples per capture |
| `source_dir` | directory the captures came from |
| `gammaspec_version` | version of the writer |

A conversion is written first to `<name>.partial` and renamed when complete, so
an interrupted run never leaves a store that looks valid.

## Index CSV

Written by `gammaspec index --csv`. One row per capture, in time order.

| Column | Description |
|---|---|
| `order` | position in time order (0-based) |
| `file_number` | number in the original file name |
| `created` | ISO timestamp |
| `t_rel_s` | seconds since the first capture |
| `n_samples`, `dt_s` | from the header |
| `path` | file path |

## Event file (`*.h5`)

Written by `gammaspec process` and by the stream processor (`events.EventWriter`); read with
`events.EventFile`. In a streamed file, a "capture" is a frame (its core): `capture` is
the frame number and `index` is the sample within the core.
Every column is a separate compressed, resizable dataset.

`/events/<column>`, one row per detected pulse:

| Column | Type | Description |
|---|---|---|
| `capture` | u4 | capture number (row in the store) |
| `index` | i8 | peak sample within the capture |
| `time_s` | f8 | capture time + (index + shift) × Δt |
| `amplitude` | f4 | amplitude in volts: template fit (isolated) or pile-up solver (clusters); NaN if not measured |
| `filter_amplitude` | f4 | detection-filter amplitude estimate (NaN for pulses added by the solver) |
| `shift` | f4 | sub-sample peak position relative to `index` |
| `fit_rms` | f4 | rms residual of the template fit |
| `flags` | u2 | bit mask, see below; 0 = clean |
| `cluster`, `cluster_size` | i8, u2 | pulses closer than 60 samples share a cluster |
| `gap_prev`, `gap_next` | i4 | samples to the neighbouring pulses (−1: none) |

Flags: 1 `PILEUP` (cluster size > 1), 2 `WARMUP` (first 1500 samples of a capture),
4 `EDGE` (no room for the fit window), 8 `BAD_FIT` (isolated but poorly fitted),
16 `LARGE` (> 1.5 V), 32 `AFTER_LARGE` (within 200 samples after a large pulse),
64 `FIT_ADDED` (found by the pile-up solver's residual search). `PILEUP` and `FIT_ADDED`
only describe how a pulse was measured: with the solver on, those pulses have valid
amplitudes. `EventFile.good()` returns every event with a usable amplitude;
`EventFile.clean()` only the unflagged (isolated) ones.

`/captures/<column>`, one row per capture (or stream frame): `capture`, `time_s`,
`live_time_s` (duration of the valid region), `live_start_s` (its start: after the warm-up;
added in v0.2, NaN when reading older files), `offset_v`, `sigma_in_v`, `sigma_out_v` (noise
before / after pre-processing), `filter_sigma_v`, `threshold_v`, `n_events`, `n_lines`.

`/snippets/`: `cluster`, `capture`, `start`, `offset`, `length` and `data` (flat float32):
the clean waveform of each stored cluster, from 20 samples before its first pulse to 80
after its last.

Attributes: `config` (the full `PipelineConfig` as JSON), `template`, `source_store` (for a
stream: a description of the source), `lsb_volts`, `dt_s`, `fit_baseline`, `pole_zero`
(filter coefficients), `stream` (the `StreamConfig`, streamed files only), `gammaspec_version`.

## Simulation truth (`*.truth.h5`, `truth/<run>.h5`)

Written by `gammaspec simulate` (next to the store) and by simulated runs; read with
`simulate.load_truth`. One row per true pulse (one per decay or background event that
deposited energy):

| Column | Type | Description |
|---|---|---|
| `time_s` | f8 | arrival time, on the same clock as the event file's `time_s` |
| `sample` | f8 | arrival, in samples since the start of the simulation |
| `amplitude_v` | f4 | true pulse amplitude (= energy / gain) |
| `energy_kev` | f4 | deposited energy after detector resolution |
| `kind` | u1 | 0 source, 1 room background |

Attributes: `sim_config` (the `SimConfig` as JSON); for a store also `capture_start_samples`
and `n_samples`.

## Live spectrum (`live/<run>/live_spectrum.json`)

Rewritten (atomically) every `stream.update_s` seconds while a stream is processed:
`live_time_s`, `events`, `good`, `counts` (per channel), `edges_v` and, once calibrated,
`calibration` (coefficients) and `edges_kev`. `live_spectrum.png` is the same as a picture.
