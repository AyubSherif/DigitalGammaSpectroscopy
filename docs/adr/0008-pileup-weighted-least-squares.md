# ADR-0008: Resolve pile-up with weighted least squares and sub-sample timing

* **Status:** Accepted; implemented in `gamma_spectrometry/solver.py` (step 6, [pileup-solver.md](../pileup-solver.md))
* **Date:** 2026-11-24

## Context

26 % of ⁶⁰Co events are in pile-up clusters. The reference paper (Mohammadian-Behbahani &
Saramad 2020) compares SER, template matching (TM), MLE and HYPER, and finds TM the most
robust at low sampling rates. TM had also worked well in earlier, similar work. Our pulses
rise within about one sample at 1.8 MS/s, so the sampling phase strongly affects any method
that reads individual samples.

The benchmark ([pileup-benchmark.md](../pileup-benchmark.md)) compared TM, HYPER, least
squares with integer times, and least squares with sub-sample timing in three weightings
(ordinary, noise-whitened, noise + template-mismatch), on template pairs in real noise
(found by the real detector), real-pulse pairs, and the real pile-up events.

## Decision

Resolve pile-up clusters by **weighted least squares with sub-sample timing** (`wls_refine`):

* model: y(t) = Σ A_k s(t − t_k) over the cluster window, with the pole-zero-corrected
  data template s;
* t_k refined jointly by Gauss–Newton (template derivative), limited to a window around the detected positions (±0.75 samples in the benchmark; ±1.25 in the solver, so a detection one sample off can still reach the true peak instead of being "explained" by a spurious extra pulse);
* weights: the inverse of C = C_noise + (ε A_max)² I, with C_noise the measured noise
  covariance and ε = 0.35 % the measured template mismatch.

The other methods stay available (`pileup.resolve(method, ...)`, ADR-0003).

## Consequences

* Small pulses next to large ones are measured to ~1 % (TM: 20–76 % spread). Pile-up
  ⁶⁰Co peaks have the same width as isolated ones and shift by only +0.15 / +0.22 %.
* Cost ~2.4 ms per cluster (one core), so ~1.3 cores at the ⁶⁰Co rate. Step 6 will
  optimise it.
* Needs the noise covariance (already measured per capture by the detector) and ε (from
  the isolated-pulse fit residuals, re-measured for new hardware).
* Pulses closer than 3 samples are a single detection; step 6 will add a residual search
  to find such hidden pulses.
