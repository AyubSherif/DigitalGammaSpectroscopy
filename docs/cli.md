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
