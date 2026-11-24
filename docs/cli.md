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

Example output:

```text
captures         : 100
first / last     : 2021-08-26 13:34:48.617564  ->  2021-08-26 13:38:32.695494
wall-clock span  : 0.062 h
live time        : 29.1 s  (13.0% duty)
...
wrote data\source_first100.h5: 100 captures, 70.7 MB (0.71 MB/capture), lsb=8.4771697000e-05 V, live time 29.1 s
```
