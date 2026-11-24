# ADR-0012: Stream by overlapping frames through the unchanged capture pipeline; validate with a simulator

* **Status:** Accepted
* **Date:** 2026-11-24
* **Completes:** [ADR-0003](0003-modular-optional-stages.md) (streams), with [ADR-0011](0011-run-configuration.md)

## Context

The pipeline was built and tuned on separate 0.29 s captures. A stream is continuous and
arrives in chunks of arbitrary size. Some stages are not causal: the baseline tracker is a
centred moving average, line removal fits blocks, and the detector estimates its noise
covariance from the whole capture. Others are long-memory filters: the pole-zero correction
has a 20,000-sample leak, and the undershoot lasts ~1,500 samples. Processing must keep up
with 1.8 MS/s at ~4,500 pulses/s.

The recorded data only covers one detector, one source and one rate. Testing streaming,
pile-up and rate behaviour needs input whose true pulses are known.

Options for streaming:
* (a) Rewrite every stage as a causal, stateful filter.
* (b) Process the stream in **overlapping frames**, each through the existing per-capture
  pipeline, and keep only the events in each frame's core.

## Decision

* (b). A frame is `pre` + `frame` + `post` samples (default 4,096 + 262,144 + 2,048).
  * `pre` > the warm-up (1,500 samples), so every kept event is past the restart transient.
  * `post` covers the fit window, the tracker's half-window and a cluster gap.
  * `check_stream_config` enforces these limits.
* Frames are independent and are processed in parallel workers. The state between chunks is
  the buffered context only, so the results do not depend on how the source splits its
  chunks.
* A recording is a sequence of *segments* (captures). A segment's last frame has no
  post-context, exactly like the end of a capture. With `frame` ≥ the capture length, a
  stream reproduces batch processing exactly.
* Live sources (stdin, a paced replay or simulation) cannot wait. Beyond `max_backlog`
  queued frames, a frame is dropped and its time is dead time (not live time), so rates stay
  correct per live second.
* A **simulator** (`simulate.py`) produces a realistic continuous waveform:
  * the measured pulse model, ⁶⁰Co cascades with summing, Klein–Nishina Compton, background;
  * noise, interference, offset, polarity and ADC rounding.

  It records the truth, and `compare_to_truth` scores any run against it.

## Consequences

* Streaming needed no change to the processing stages. Every improvement of the batch
  pipeline applies to streams.
* Overhead: (pre + post) / frame = 2.3 % extra processing.
* Latency: one frame plus its post-context, ≈ 0.15 s (plus the processing time).
* Validation:
  * **Real data:** streamed frames give the same spectrum as batch processing (FWHM
    6.56 / 6.04 % vs 6.54 / 6.05 %; events within 0.03 %).
  * **Simulated data:** the frames are as accurate as whole captures (better, against the
    truth: 0.84 vs 1.16 mV rms amplitude error).
* The speed is that of batch processing. This PC (6 threads) reaches ~0.9× real time at
  4,500 pulses/s, so a live 1.8 MS/s input loses ~10 % of its frames as dead time. Faster
  machines or a faster solver remove that; correctness does not depend on it.
* Clusters that cross a core boundary are measured in both frames, and each frame keeps its
  own events. The cluster ids then differ on the two sides.
