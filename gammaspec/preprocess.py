"""Waveform pre-processing chain with individually optional stages (ADR-0003).

    raw volts ─► polarity + offset ─► [line removal] ─► [pole-zero] ─► [baseline tracker] ─► clean signal

Each bracketed stage can be switched off in ``PreprocessConfig``. The chain works on
one capture (or chunk) at a time; ``reset()`` clears filter state between independent
captures.

Recommended: all three on, then fit pulses with the baseline fixed at 0. Pole-zero
must be paired with the tracker (its low-frequency gain amplifies slow drift); on the
2021 data pole-zero alone is worse than no processing. The first ``warmup_samples``
of a capture are not fully corrected (tails of pulses before the capture started),
and events there should be flagged by later stages.
"""

from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass, field

import numpy as np

from .baseline import BaselineTracker, PoleZero, PoleZeroDesign, TrackerConfig
from .interference import LineRemoval, LineRemovalConfig
from .signal import DEFAULT_POLARITY, find_pulses, noise_sigma, pulse_free_mask, to_positive


@dataclass
class PreprocessConfig:
    polarity: int = DEFAULT_POLARITY     # -1: pulses are negative-going (flipped to positive); 1: positive
    line_removal: bool = True            # stage: remove narrow interference lines (sinusoids)
    pole_zero: bool = True               # stage: pole-zero correction (removes the undershoot); needs baseline_tracker
    baseline_tracker: bool = True        # stage: subtract the slow baseline (moving average of pulse-free samples)
    trigger_threshold: float = 5.0       # for the pulse masks used by the stages
    line_mask_pre: int = 20              # line fit: samples excluded before each pulse ...
    line_mask_post: int = 60             # ... and after it
    warmup_samples: int = 1500           # start of each capture: earlier pulses' tails are not corrected
    lines: LineRemovalConfig = field(default_factory=LineRemovalConfig)
    tracker: TrackerConfig = field(default_factory=TrackerConfig)

    def label(self) -> str:
        on = [n for n in ("line_removal", "pole_zero", "baseline_tracker") if getattr(self, n)]
        return " + ".join(on) if on else "raw"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class PreprocessInfo:
    offset_v: float = 0.0
    sigma_in_v: float = 0.0
    sigma_out_v: float = 0.0
    line_freqs_hz: np.ndarray = field(default_factory=lambda: np.empty(0))
    line_rms_v: float = 0.0


class Preprocessor:
    def __init__(self, fs: float, config: PreprocessConfig | None = None,
                 pole_zero: PoleZeroDesign | None = None):
        self.fs = fs
        self.cfg = config or PreprocessConfig()
        if self.cfg.pole_zero and pole_zero is None:
            raise ValueError("pole_zero is enabled but no PoleZeroDesign was given")
        if self.cfg.pole_zero and not self.cfg.baseline_tracker:
            warnings.warn("pole_zero without baseline_tracker: the P/Z filter amplifies slow baseline "
                          "drift; enable baseline_tracker unless you are studying the filter itself",
                          stacklevel=2)
        self.lines = LineRemoval(fs, self.cfg.lines) if self.cfg.line_removal else None
        self.pz = PoleZero(pole_zero) if self.cfg.pole_zero else None
        self.tracker = BaselineTracker(self.cfg.tracker) if self.cfg.baseline_tracker else None
        self.info = PreprocessInfo()

    def reset(self) -> None:
        for st in (self.lines, self.pz, self.tracker):
            if st is not None:
                st.reset()

    def process(self, volts: np.ndarray) -> np.ndarray:
        """One capture/chunk of raw volts -> clean positive-pulse signal (volts)."""
        c = self.cfg
        x = to_positive(np.asarray(volts, dtype=np.float64), c.polarity)
        # The coupled signal has zero net pulse area, so its mean is the electronic offset.
        offset = float(x.mean())
        x = x - offset
        info = PreprocessInfo(offset_v=offset, sigma_in_v=noise_sigma(x))
        if self.lines is not None:
            peaks = find_pulses(x, info.sigma_in_v, c.trigger_threshold)
            free = pulse_free_mask(x.size, peaks, c.line_mask_pre, c.line_mask_post)
            x = self.lines.process(x, free)
            info.line_freqs_hz = self.lines.info.freqs_hz
            info.line_rms_v = self.lines.info.removed_rms_v
        if self.pz is not None:
            x = self.pz.process(x)
        if self.tracker is not None:
            sig = noise_sigma(x)
            peaks = find_pulses(x, sig, c.trigger_threshold)
            free = pulse_free_mask(x.size, peaks, c.tracker.mask_pre, c.tracker.mask_post)
            x = self.tracker.process(x, free)
        info.sigma_out_v = noise_sigma(x)
        self.info = info
        return x
