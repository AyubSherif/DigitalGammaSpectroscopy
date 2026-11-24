"""A live spectrum window for ``gamma-spectrometry stream --show`` (matplotlib, interactive).

``LivePlot.update(status, spectrum)`` is the ``on_update`` callback of the stream
processor: it redraws the spectrum (in keV once the live spectrum is calibrated, in mV
before that) and the status line. Closing the window does not stop the stream.
"""

from __future__ import annotations

import numpy as np


class LivePlot:
    def __init__(self, title: str = "Gamma Spectrometry App live spectrum"):
        import matplotlib.pyplot as plt

        from . import plotting as P

        self.plt = plt
        plt.ion()
        self.fig, self.ax = plt.subplots(figsize=(9, 4.5))
        self.fig.canvas.manager.set_window_title(title)
        self.line = None
        self.color = P.SERIES[0]

    def update(self, status, spectrum) -> None:
        if not self.plt.fignum_exists(self.fig.number):
            return                                              # window closed: keep streaming
        cal = spectrum.calibration
        x = cal.energy(spectrum.edges_v) if cal is not None else spectrum.edges_v * 1e3
        width = np.diff(x)
        y = np.maximum(spectrum.counts / max(spectrum.live_time_s, 1e-9) / width, 1e-6)
        centres = 0.5 * (x[1:] + x[:-1])
        ax = self.ax
        ax.clear()
        ax.step(centres, y, where="mid", color=self.color, linewidth=1.0)
        ax.set_yscale("log")
        ax.set_xlabel("energy (keV)" if cal is not None else "amplitude (mV)")
        ax.set_ylabel("counts / live s / " + ("keV" if cal is not None else "mV"))
        ax.set_title(f"{status.stream_s:.1f} s of signal ({status.speed:.2f}x real time), "
                     f"{int(spectrum.counts.sum())} good events, {status.dropped_frames} frames dropped",
                     fontsize=10)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        self.fig.tight_layout()
        self.fig.canvas.draw_idle()
        self.plt.pause(0.01)

    def keep_open(self) -> None:
        if self.plt.fignum_exists(self.fig.number):
            self.plt.ioff()
            self.plt.show()
