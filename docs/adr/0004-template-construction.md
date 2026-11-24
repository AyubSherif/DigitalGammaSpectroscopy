# ADR-0004: Pulse template built from data on a sub-sample grid with phase equalisation

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

The pulse rises in about one ADC sample (0.56 µs), so the shape of the rising edge
can't be seen from any single pulse. Energy and pile-up fits need a template that
can be evaluated at any sub-sample time. The options were an analytic model (sum
of exponentials) or a template measured from the data.

## Decision

* **Measured template (primary).** Isolated pulses are aligned by their fitted
  sub-sample peak time and averaged on a grid of 1/16 sample. Because pulses arrive
  at random sampling phases, together they sample the shape much more finely than
  any single pulse (equivalent-time sampling).
* **Iterative fit.** Each pulse is fitted as `A·s(t − δ) + b` (a grid search over δ,
  parabolic refinement, then an exact linear solve for A and b). The template is
  rebuilt from the fits until it changes by less than 5 × 10⁻⁴ (rms, fraction of peak).
* **Phase equalisation.** On each iteration the fitted sampling phases δ mod 1 are
  mapped onto a uniform distribution. Without this, a template warped within each
  sample period and correspondingly warped δ values fit the data almost equally
  well. In tests this left a ±0.04-sample timing bias and a 5 % rising-edge
  error that did *not* shrink with more pulses. With equalisation, the error is
  statistical and shrinks as more pulses are used.
* **Constraints:** pre-pulse level fixed at 0; peak fixed at t = 0 with value 1;
  light smoothing (no blurring within ±3 samples of the peak, about 1 sample in the
  slow tail).
* **Analytic model (secondary).** A sum of exponentials, smoothed by a Gaussian,
  with the pulse starting from zero. It gives physical time constants (light decay,
  RC, coupling) for stages that need them, such as pole-zero correction.

## Consequences

* On synthetic data, 300 pulses reproduce the shape to within 0.1 % rms (under 3.5 % at
  the worst point on the rising edge), timing to about 1/40 sample, and amplitudes to
  within noise.
* The template depends on the detector, its high voltage and the electronics. It must
  be rebuilt when any of these change, and is stored with its provenance in
  `template.h5`.
* For high-rate data, the isolation window limits how much of the tail can be
  measured. The long undershoot tail is best measured on low-rate data (the background
  run) or modelled.
