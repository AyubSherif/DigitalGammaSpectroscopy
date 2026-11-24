# ADR-0009: Calibrate each run on its own lines; gain-match short runs to a reference

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

The source run (26 Aug) and the background run (31 Aug) were recorded five days apart. The
⁴⁰K line shows that the gain differs by about 7 %. Background subtraction must compare the
two runs on the same energy axis. The background run is only 29 s long, and only ⁴⁰K is
significant in it.

Options: (a) subtract in volts after scaling one run by a single factor; (b) calibrate each
run independently to keV; (c) calibrate the source run and transfer its calibration to the
background run with one gain factor from ⁴⁰K.

## Decision

* Calibrate every run against its own lines (automatic anchor, bounded line fits with
  iterative removal of rejected lines, linear or quadratic curve). Spectra are compared in
  keV.
* When a run has only one usable line, **gain-match** it: use the shape (offset, curvature)
  of a reference run's calibration and fit one gain factor, E(A) = E_ref(k·A). This is (c)
  as a fallback within (b), and it is equivalent to (a) plus a correct treatment of the
  offset.
* Normalise each spectrum by its own live time; propagate Poisson errors of both runs.

## Consequences

* The background run's ⁴⁰K lands at 1460.8 keV by construction; the source run's
  independent calibration predicts its ⁴⁰K within 0.6 mV (≈ 2 keV).
* Gain matching assumes the offset and non-linearity did not change between runs. With a
  longer background run (two or more lines) this assumption is no longer needed.
* Gain drift *within* a run is not corrected yet. It is measured (`gain_vs_time`) and was
  ≤ 0.12 % over the first 3.7 min; the full 3-hour run will show whether a time-dependent
  calibration is needed.
