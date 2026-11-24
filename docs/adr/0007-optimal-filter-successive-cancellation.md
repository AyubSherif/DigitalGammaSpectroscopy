# ADR-0007: Detect pulses with a per-capture optimal filter and successive cancellation

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

After pre-processing, the signal has a flat baseline, but its noise is not white: the
residual has excess low-frequency content (≈ 0.7 mV rms, of which the white floor is
≈ 0.45 mV), and in some background captures there are bursts of broadband interference.
Detection must find small pulses (low threshold), separate close pulses (pile-up), and not
invent pulses (false triggers, artifacts of large pulses).

Measured with pulses injected into real noise (`gammaspec evaluate-detection`):

| Filter | Output noise | Smallest pulse found reliably (≥ 96 %) | Two pulses resolved (300/300 mV; 300/30 mV) | Extra detections per 1 V pulse |
|---|---|---|---|---|
| Leading edge (coarse trigger) | 0.49 mV | 6 mV | equal: ≥ 5; unequal: ≥ 5–15 samples | 0.0 (0.2 at 2 V) |
| Matched filter | 0.96 mV | 8 mV | equal: ≥ 6; unequal: ≥ 30 samples | 0.0 |
| Sharpened (decay removed) | 0.55 mV | 6 mV | equal: ≥ 3; unequal: ≥ 4 samples | 0.24 (2.9 at 2 V) |
| Optimal filter, peak search only | 0.25 mV | 3 mV | equal: ≥ 5; unequal: ≥ 15 samples | 3.0 |
| **Optimal filter + successive cancellation** | **0.25 mV** | **2 mV (96 %), 3 mV (100 %)** | **equal: ≥ 4; unequal: ≥ 3 samples** | **0.01** |

The plain matched filter is optimal only for white noise and here collects the
low-frequency excess. The optimal filter whitens the noise, which also makes its response
narrow (good pile-up separation). But its response to a pulse has negative dips and small
bumps up to ~50 samples before the peak, which a simple peak search (even with a
prominence condition) reports as extra pulses.

## Decision

* **Optimal filter per capture:** the noise autocovariance is estimated from that
  capture's pulse-free samples, and the kernel h = C⁻¹s is solved with the Levinson
  algorithm. It adapts when interference changes between captures.
* **Successive cancellation (CLEAN):** take the largest remaining peak, accept it if
  it exceeds the local threshold plus a margin (3 % of any accepted pulse whose response
  overlaps), subtract its predicted response (computed exactly from the template at the
  interpolated sub-sample phase), and repeat.
* **Threshold:** 5.5 × the local noise sigma of the filter output (blocks of 4096
  samples), estimated robustly from the negative half of the *residual after
  cancellation*. Real pulses' negative dips would otherwise inflate it by ~50 %.
* Other filters stay available (`method=`), per ADR-0003.

## Consequences

* ⁶⁰Co run: threshold 1.53 mV (≈ 6 keV), 4,545 events/s. For pulse pairs where both
  exceed 5 mV, the number at each separation ≥ 4 samples matches random (Poisson)
  arrivals: no spurious pile-up.
* Nothing is accepted within ±2 samples of an accepted pulse: on real data that region
  holds template-mismatch residuals (≈ 4–5 % of the pulse), which a residual search
  would otherwise report as ~9 % extra events. Real pulses that close become one event;
  the pile-up solver can recognise them from a poor fit.
* Below ~5 mV (≈ 20 keV), broad low-frequency noise bumps are sometimes detected and
  split into two events; they are flagged as pile-up and stay out of the clean spectrum.
* Heavy-tailed interference bursts in the background run still cause ≈ 300 /s false
  triggers between 10 and 30 mV in noisy captures; above 30 mV quiet and noisy captures
  agree. Later stages should apply an energy threshold or select captures by noise
  (`captures.sigma_in_v`).
* Pulses above ~1.5 V (cosmic muons) are flagged, as is the region after them.
