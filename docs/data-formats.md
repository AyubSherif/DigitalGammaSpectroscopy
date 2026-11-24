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
