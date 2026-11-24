# ADR-0006: Remove interference by subtracting sinusoids fitted to pulse-free samples

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

The noise is dominated by narrow lines: a family at multiples of 9 kHz (= fs/200) in
both runs, and an intermittent, slightly frequency-modulated line near 257 kHz in the
background run. The white-noise floor is only about 0.45 mV, while the lines bring the
total to 1.4–6.6 mV.

Options considered:

* **Notch filters:** simple and streaming-friendly, but every notch rings after every
  pulse, which distorts the pulse shape in a way that depends on amplitude. A dozen
  notches add up.
* **Folding at period 200 (synchronous averaging):** the 9 kHz family is almost locked
  to the ADC clock, but not exactly; the pattern drifts within a capture (half-to-half
  difference larger than the pattern). It helped only from 6.6 to 5.2 mV.
* **Sinusoid subtraction on pulse-free samples:** detect lines in the capture's own
  spectrum, fit cos/sin amplitudes block by block using only samples away from pulses,
  subtract everywhere.

## Decision

* Use sinusoid subtraction (`interference.LineRemoval`): lines detected per capture
  (they switch on and off), more than 100× above the local spectral median, each the
  strongest within ±3.5 kHz, at most 40. Fitted in 1024-sample blocks with 50 % overlap
  and sin² cross-fade.
* Solve the masked least squares by filling masked samples with the current model and
  refitting with one fixed design matrix for all blocks (three passes). This is exact at
  convergence and fast (≈ 0.09 s per capture).
* Line removal is optional and runs before pole-zero correction.

## Consequences

* Quiet and source captures reach the white floor (≈ 0.45–0.5 mV); pulses are not
  distorted (synthetic test: peak samples unchanged within noise).
* Broadband noise (noisy background captures, σ ≈ 4 mV after removal) cannot be
  removed this way.
* Strictly periodic pulse masks create sidebands around strong lines. In real data,
  arrivals are random, and the ±3.5 kHz uniqueness rule rejects such sidebands anyway.
