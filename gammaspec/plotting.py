"""Shared matplotlib style for diagnostic plots (static PNGs)."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

#: Categorical series colours, used in this fixed order (validated palette, light surface).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"

STYLE = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelsize": 9,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "legend.labelcolor": INK_2,
    "lines.linewidth": 1.5,
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
    "axes.prop_cycle": matplotlib.cycler(color=SERIES),
}


def figure(nrows=1, ncols=1, width=11.0, height_per_row=3.2, **kw):
    plt.rcParams.update(STYLE)
    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height_per_row * nrows),
                             constrained_layout=True, squeeze=False, **kw)
    return fig, axes


def save(fig, path, dpi=110):
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
