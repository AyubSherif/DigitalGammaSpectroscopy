# ADR-0003: Processing stages are modular and individually optional

* **Status:** Accepted. Stage selection, parameters, checks and recorded configuration are implemented
  by [ADR-0011](0011-run-configuration.md) (step 9); streams by [ADR-0012](0012-streaming-overlapping-frames.md) (step 10).
* **Date:** 2026-11-24

## Context

The pipeline applies several corrections: baseline restoration, pile-up
resolution, energy calibration, background subtraction. Users need to compare
results with and without a given correction, apply only some of them, and swap
algorithms (for example template matching vs a least-squares fit for pile-up).
The software is also meant to become an app with a user interface, and later to
process live streams.

## Decision

* Each processing step is a self-contained **stage** with a common interface:
  it takes a block of data plus its own state, and returns a block.
* A **pipeline configuration** (a file or object; a UI later) lists which stages
  run, in what order, and with what parameters. Any stage that isn't required
  can be disabled, e.g. `baseline_restoration: {enabled: false}`.
* Algorithm choices inside a stage are configuration too, e.g.
  `pileup: {method: least_squares | template_matching | hyper | reject}`.
* Stages keep their state between chunks, so the same pipeline runs on stored
  files and on live streams.
* Outputs record the configuration used (stored with the results), so every
  spectrum can be traced to its processing.

## Consequences

* Each stage is tested on its own, and stages can be benchmarked against each
  other on the same input.
* Stages must declare what they need from earlier stages (for example, pile-up
  resolution needs detected events). The pipeline checks this at configuration
  time and fails with a clear message if a requirement is missing.
* The stages are developed first as plain functions. The configuration layer is
  added once the stages are settled (roadmap step 8), so the interface isn't
  designed before its users exist.
