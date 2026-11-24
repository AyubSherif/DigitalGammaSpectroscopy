# ADR-0001: Store captures as lossless integer ADC codes in HDF5

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

The raw data is ASCII: about 7.3 MB and 524,288 samples per capture, and about
33 GB for the 2021 run. Parsing text is slow, and every later processing step
would pay for it again. The values are quantised: each one is an integer ADC
code times a fixed step (84.771697 µV), printed to 7 significant digits.

## Decision

* Convert each capture once to **int32 ADC codes** and store them in a single
  **HDF5** file per dataset. Each row is one capture, compressed with
  gzip level 1 + shuffle.
* Check every value against the ADC grid while converting. Abort if any value
  is more than 0.02 LSB off the grid.
* Keep the ADC step, sample interval and per-capture metadata (timestamp, file
  number, source path) in the file.

## Alternatives considered

* **float32 volts:** not exactly lossless, and the same size or larger.
* **One `.npy` file per capture:** simple, but no compression and thousands of
  files.
* **Parquet/Arrow:** well suited to tables (it may be used for the event list),
  but awkward for fixed-length waveform rows.
* **Higher gzip levels / LZF:** level 4 saves only about 7 % more for about
  45 % more write time; LZF is 35 % larger. Level 1 is the balance.

## Consequences

* About 10× smaller than text (0.7 MB per capture); reading a capture takes
  about 10 ms instead of about 75 ms of parsing.
* `volts = codes × lsb_volts` reproduces the original text to within 5 × 10⁻⁷
  relative, and no information is lost.
* Data from hardware with a different ADC step needs `--lsb` or `--estimate-lsb`.
* HDF5 is not written concurrently: parsing runs in parallel worker processes,
  and a single process writes.
