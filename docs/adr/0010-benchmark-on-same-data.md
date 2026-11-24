# ADR-0010: Benchmark against a naive count and an emulated analog chain on the same data

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

We want to know how much the processing improves the spectrum. We need a lower reference:
a naive peak count with no transformation. We also need a reference representing standard
practice: an NaI(Tl) detector with preamplifier, shaping amplifier and MCA.

There is no analog chain for this detector, and published resolutions come from other
crystals and PMTs. At 1.3 MeV, NaI(Tl) resolution is dominated by the crystal and the PMT,
so comparing our numbers with literature numbers mixes detector and processing effects.

## Decision

* Emulate the analog chain **digitally, on the recorded waveforms**:
  * a pole-zero-cancelled CR-RC⁴ shaper (0.5 and 1 µs) and a baseline restorer;
  * a peak-hold ADC and a lower-level discriminator;
  * an optional pile-up rejector;
  * a 5 µs non-extending dead time with a live-time clock.
* Run all chains (naive, analog variants, gammaspec) on the **same captures**. Calibrate each
  chain with the **same automatic calibration**, and measure it with the **same metrics**
  (FWHM, valley-to-peak, peak rates, pile-up continuum, threshold, recorded / live fraction).
* Quote literature values for typical 3″×3″ NaI(Tl) analog systems **as context only**, and
  state that they come from other detectors.

## Consequences

* Differences between chains come from the processing alone.
* The emulated analog chain is optimistic: it has no electronic noise of its own, an ideal
  pole-zero setting and a perfect ADC. The improvement of gammaspec over it is therefore a
  lower bound.
* The emulation runs at about real-time speed in Python, so the benchmark uses 100 captures
  by default (`--n` selects more).
