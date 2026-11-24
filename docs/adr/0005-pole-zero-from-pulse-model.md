# ADR-0005: Pole-zero correction derived from the pulse model, always followed by a baseline tracker

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

The pulses have an undershoot from AC coupling (net pulse area ≈ 0). The analytic fit
needs two slow exponentials for it (τ = 76 and 318 samples). Three approaches were
tried on the measured template:

1. **Two ideal CR inversions** (zero at each slow pole, integrator at DC): over-restores
   (the tail settles at +8 %). Two ideal CR stages would also force the first moment to
   zero, which the data does not show, so the second slow component is not a second
   ideal CR stage.
2. **Fitting free zeros and poles for a flat tail:** under-determined. The fitted values
   jump with the fit range (22/176 vs 110/370 samples). Leaving the extra pole free
   made the filter sharpen pulses (high-frequency gain ×2, which doubles white noise).
3. **Exact pole-zero swap from the model:** write the sampled model as
   N(z)/Π(1 − aᵢz⁻¹). Put zeros at the slow poles aᵢ and poles at the roots of N(z)
   near z = 1. This is determined by the model alone, keeps unit high-frequency gain,
   and flattens the tail to about 0.1–0.3 % of the peak.

## Decision

* Use approach 3 (`baseline.design_pole_zero`, `design_from_template`). A root at or
  above 1 (net area very slightly negative) is clamped to a 20,000-sample leak so the
  filter stays stable.
* The filter has a large low-frequency gain (≈ ×13 at 90 Hz), so **pole-zero is always
  paired with the baseline tracker** (moving average of pulse-free samples, 1001
  samples). The pipeline warns if pole-zero is enabled without it.
* The filter keeps its state between chunks (streaming-ready). Templates used with
  corrected data are filtered the same way (`filter_template`).

## Consequences

* On the ⁶⁰Co data, pulses on an earlier pulse's tail are measured identically to clean
  pulses (peak position within 0.03 %; FWHM unchanged), and the pre-pulse baseline
  narrows from 2.1 to 0.5 mV robust σ.
* The design depends on the model fitted to a long template (low-rate data). The
  coupling constants should be re-measured whenever the electronics change.
* The first ~1500 samples of each capture are not fully corrected (tails of pulses that
  occurred before the capture started).
