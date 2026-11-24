# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- `gammaspec.rawio`: reader for ASCII ADC capture files. It has a vectorised
  fixed-width parser (about 7 M samples/s on one core) and a generic fallback,
  converts volts to ADC codes losslessly with grid validation, estimates the
  ADC step, and builds a chronological index.
- `gammaspec.store`: parallel conversion of captures to a compressed HDF5 store
  (int32 codes, gzip+shuffle, about 10× smaller than text), and a
  `CaptureStore` reader.
- CLI commands `index` and `convert` (`python -m gammaspec`).
- Test suite with synthetic captures and real-data checks.
- Documentation: architecture, data formats, CLI reference, development guide,
  ADRs 0001–0003.
