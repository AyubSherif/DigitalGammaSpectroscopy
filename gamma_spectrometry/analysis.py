"""Energy spectra and background subtraction (step 7).

    source events ──► calibrate (60Co lines) ──► energies ──► spectrum (counts/s per bin)
    background events ► calibrate / gain-match ► energies ► spectrum (counts/s per bin)
    net = source rate - background rate,  sigma^2 = N_s / t_s^2 + N_b / t_b^2

Both runs are calibrated on their own (gain matching through calibration), histogrammed on
the same energy bins and normalised by their own live time. A background line (40K) must
vanish in the net spectrum; that is the check of the gain matching and the normalisation.

Run:  python -m gamma_spectrometry spectrum SOURCE_EVENTS OUT_DIR --background BACKGROUND_EVENTS
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

from . import plotting as P
from .calibration import Calibration, calibrate, gain_vs_time
from .events import EventFile

#: line groups for the net-rate fits: (lines reported, nuisance lines, energy window keV).
#: Each group gets its own local window with a linear continuum; 40K sits on the
#: high-energy flank of the 1332 keV peak, so that peak is fitted with it as a nuisance.
NET_GROUPS = [
    ([("Co-60", 1173.228), ("Co-60", 1332.492)], [], (1040.0, 1440.0)),
    ([("K-40", 1460.820)], [("Co-60", 1332.492)], (1300.0, 1650.0)),
]


@dataclass
class EnergySpectrum:
    edges_kev: np.ndarray
    counts: np.ndarray
    live_time_s: float
    label: str = ""

    @property
    def centers(self) -> np.ndarray:
        return 0.5 * (self.edges_kev[1:] + self.edges_kev[:-1])

    @property
    def rate(self) -> np.ndarray:
        """Counts per second in each bin."""
        return self.counts / self.live_time_s

    @property
    def rate_err(self) -> np.ndarray:
        return np.sqrt(np.maximum(self.counts, 1.0)) / self.live_time_s


@dataclass
class NetSpectrum:
    edges_kev: np.ndarray
    rate: np.ndarray            # counts/s per bin
    rate_err: np.ndarray

    @property
    def centers(self) -> np.ndarray:
        return 0.5 * (self.edges_kev[1:] + self.edges_kev[:-1])


def energy_spectrum(energies_kev, live_time_s, edges_kev, label="") -> EnergySpectrum:
    counts, _ = np.histogram(energies_kev, edges_kev)
    return EnergySpectrum(np.asarray(edges_kev, float), counts.astype(float), float(live_time_s), label)


def subtract(src: EnergySpectrum, bkg: EnergySpectrum) -> NetSpectrum:
    if not np.allclose(src.edges_kev, bkg.edges_kev):
        raise ValueError("spectra must share their energy bins")
    rate = src.rate - bkg.rate
    err = np.sqrt(np.maximum(src.counts, 1.0) / src.live_time_s ** 2 + np.maximum(bkg.counts, 1.0) / bkg.live_time_s ** 2)
    return NetSpectrum(src.edges_kev, rate, err)


def _gauss_lin(x, *p):
    n = (len(p) - 2) // 3
    y = p[-2] + p[-1] * (x - x.mean())
    for k in range(n):
        area, mu, sig = p[3 * k: 3 * k + 3]
        y = y + area / (np.sqrt(2 * np.pi) * sig) * np.exp(-0.5 * ((x - mu) / sig) ** 2)
    return y


def fit_line_groups(net: NetSpectrum, groups=None):
    """Fit every group of NET_GROUPS in its own window; returns the reported lines (rates)."""
    out = []
    for lines, nuisance, window in (groups or NET_GROUPS):
        res, chi2 = fit_net_peaks(net, list(lines) + list(nuisance), window=window,
                                  nonnegative=tuple(range(len(lines), len(lines) + len(nuisance))))
        for r in res[: len(lines)]:
            r["window_kev"] = list(window)
            r["fit_chi2_ndf"] = chi2
            out.append(r)
    return out


def fit_net_peaks(net: NetSpectrum, lines, fwhm_rel_at_1173: float = 0.065, window_sigmas: float = 4.0,
                  window=None, nonnegative=()):
    """Gaussian fits of lines on a linear continuum in a (net) rate spectrum.
    Returns, per line, the peak rate (counts/s) and its uncertainty. Widths and centres
    are bounded around the expected values. Areas may be negative (a line that should
    vanish in the net spectrum), except for the line indices in ``nonnegative``
    (nuisance peaks, which can be absent but not negative)."""
    e = np.array([x for _, x in lines], float)
    sig = e * fwhm_rel_at_1173 / 2.3548 * np.sqrt(1173.0 / e)
    lo, hi = window if window is not None else (e.min() - window_sigmas * sig[0], e.max() + window_sigmas * sig[-1])
    c = net.centers
    m = (c >= lo) & (c < hi)
    x, y, s = c[m], net.rate[m] / np.diff(net.edges_kev)[m], net.rate_err[m] / np.diff(net.edges_kev)[m]
    p0, lb, ub = [], [], []
    for k, (ek, sk) in enumerate(zip(e, sig)):
        near = np.abs(x - ek) < sk
        p0 += [max(float(y[near].max() if near.any() else 0), 1e-6) * np.sqrt(2 * np.pi) * sk, ek, sk]
        lb += [0.0 if k in nonnegative else -np.inf, ek * 0.98, sk * 0.6]
        ub += [np.inf, ek * 1.02, sk * 1.6]
    p0 += [float(np.median(y)), 0.0]
    lb += [-np.inf, -np.inf]
    ub += [np.inf, np.inf]
    p0 = list(np.clip(p0, np.nextafter(lb, np.inf), np.nextafter(ub, -np.inf)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        popt, pcov = curve_fit(_gauss_lin, x, y, p0=p0, sigma=s, absolute_sigma=True, bounds=(lb, ub), maxfev=20000)
    perr = np.sqrt(np.diag(pcov))
    chi2 = float(np.sum(((y - _gauss_lin(x, *popt)) / s) ** 2) / max(1, x.size - len(popt)))
    out = []
    for k, (lab, ek) in enumerate(lines):
        out.append(dict(label=lab, energy_kev=ek, rate_cps=float(popt[3 * k]), rate_err_cps=float(perr[3 * k]),
                        centre_kev=float(popt[3 * k + 1]), fwhm_pct=float(100 * 2.3548 * popt[3 * k + 2] / popt[3 * k + 1])))
    return out, chi2


def analyze(source_events, out_dir, background_events=None, bin_kev: float = 5.0, e_max_kev: float = 3000.0,
            reliable_above_kev: float = 120.0, verbose: bool = True) -> dict:
    """Calibrate, build spectra, subtract the background, fit net peaks, write outputs."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    src_ef = EventFile(source_events)
    src_ev = src_ef.good()
    cal_s = calibrate(src_ev["amplitude"], "source")
    cal_s.save(out / "calibration_source.json")
    edges = np.arange(0.0, e_max_kev + bin_kev / 2, bin_kev)
    S = energy_spectrum(cal_s.energy(src_ev["amplitude"]), src_ef.live_time_s, edges, "source")
    summary = dict(source=str(source_events), live_time_source_s=S.live_time_s,
                   calibration_source=cal_s.to_dict(), gain_vs_time_source=gain_vs_time(src_ev, "source"))
    B = N = cal_b = None
    if background_events is not None:
        bkg_ef = EventFile(background_events)
        bkg_ev = bkg_ef.good()
        cal_b = calibrate(bkg_ev["amplitude"], "background", reference=cal_s)
        cal_b.save(out / "calibration_background.json")
        B = energy_spectrum(cal_b.energy(bkg_ev["amplitude"]), bkg_ef.live_time_s, edges, "background")
        N = subtract(S, B)
        net_lines = fit_line_groups(N)
        src_lines = fit_line_groups(NetSpectrum(edges, S.rate, S.rate_err))
        bkg_lines = fit_line_groups(NetSpectrum(edges, B.rate, B.rate_err))
        summary.update(background=str(background_events), live_time_background_s=B.live_time_s,
                       calibration_background=cal_b.to_dict(),
                       gain_ratio_background_to_source=float(cal_b.gain_kev_per_v / cal_s.gain_kev_per_v),
                       gain_vs_time_background=gain_vs_time(bkg_ev, "background"),
                       lines=dict(source=src_lines, background=bkg_lines, net=net_lines),
                       totals_cps=dict(source=float(S.rate[edges[:-1] >= reliable_above_kev].sum()),
                                       background=float(B.rate[edges[:-1] >= reliable_above_kev].sum()),
                                       net=float(N.rate[edges[:-1] >= reliable_above_kev].sum())),
                       reliable_above_kev=reliable_above_kev)
    _write_csv(out / "spectrum.csv", S, B, N)
    (out / "spectrum_summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
    _plots(out, cal_s, cal_b, S, B, N, reliable_above_kev)
    if verbose:
        print(f"source calibration: E = {cal_s.coeffs[0]:+.2f} keV + {cal_s.coeffs[1]:.2f} keV/V x A"
              + (f" + {cal_s.coeffs[2]:.1f} keV/V^2 x A^2" if len(cal_s.coeffs) > 2 else ""))
        if cal_b is not None:
            print(f"background: {cal_b.kind}; gain ratio background/source (keV per V) "
                  f"{summary['gain_ratio_background_to_source']:.4f}")
    return summary


def _write_csv(path, S, B, N):
    cols = ["e_lo_kev", "e_hi_kev", "source_counts", "source_cps"]
    data = [S.edges_kev[:-1], S.edges_kev[1:], S.counts, S.rate]
    if B is not None:
        cols += ["background_counts", "background_cps", "net_cps", "net_cps_err"]
        data += [B.counts, B.rate, N.rate, N.rate_err]
    arr = np.column_stack(data)
    np.savetxt(path, arr, delimiter=",", header=",".join(cols), comments="", fmt="%.6g")


def _plots(out, cal_s, cal_b, S, B, N, reliable_above):
    # calibration
    fig, ax = P.figure(1, 2, height_per_row=3.4)
    for k, (cal, lab) in enumerate([(cal_s, "source run"), (cal_b, "background run")]):
        if cal is None:
            continue
        acc = [p for p in cal.points if p.accepted]
        A = np.array([p.amplitude_v for p in acc])
        E = np.array([p.energy_kev for p in acc])
        aa = np.linspace(0, max(A.max() * 1.1, 0.8), 100)
        ax[0, 0].plot(aa, cal.energy(aa), color=P.SERIES[k], linewidth=1.2, label=f"{lab}: {cal.kind}")
        ax[0, 0].plot(A, E, "o", color=P.SERIES[k], markersize=5)
        for p in acc:
            ax[0, 0].annotate(p.label, (p.amplitude_v, p.energy_kev), textcoords="offset points", xytext=(5, -10),
                              fontsize=7, color=P.INK_2)
        ax[0, 1].errorbar(E, cal.energy(A) - E, yerr=[p.amplitude_err_v * cal.gain_kev_per_v for p in acc],
                          fmt="o", color=P.SERIES[k], markersize=5, label=lab)
    ax[0, 0].set_xlabel("pulse amplitude (V)")
    ax[0, 0].set_ylabel("energy (keV)")
    ax[0, 0].set_title("Energy calibration")
    ax[0, 0].legend()
    ax[0, 1].axhline(0, color=P.AXIS, linewidth=0.8)
    ax[0, 1].set_xlabel("line energy (keV)")
    ax[0, 1].set_ylabel("calibration residual (keV)")
    ax[0, 1].set_title("Calibration residuals")
    ax[0, 1].legend()
    P.save(fig, out / "calibration.png")

    # spectra (rate density, log)
    w = np.diff(S.edges_kev)
    fig, ax = P.figure(2 if N is not None else 1, 1, height_per_row=3.6)
    a = ax[0, 0]
    a.stairs(S.rate / w, S.edges_kev, color=P.SERIES[0], linewidth=1.1, label=f"source run ({S.live_time_s:.1f} s live)")
    if B is not None:
        a.stairs(B.rate / w, B.edges_kev, color=P.SERIES[1], linewidth=1.1,
                 label=f"background run ({B.live_time_s:.1f} s live)")
    a.axvspan(0, reliable_above, color=P.GRID, alpha=0.6, linewidth=0)
    a.set_yscale("log")
    a.set_xlabel("energy (keV)")
    a.set_ylabel("counts / s / keV")
    a.set_title("Calibrated spectra (shaded: below the reliable range of the background run)")
    for lab, e in [("Co-60", 1173.2), ("Co-60", 1332.5), ("K-40", 1460.8), ("sum", 2505.7), ("Tl-208", 2614.5)]:
        a.axvline(e, color=P.MUTED, linewidth=0.6, linestyle=":")
    a.legend()
    if N is not None:
        b = ax[1, 0]
        b.axhline(0, color=P.AXIS, linewidth=0.8)
        b.axvspan(0, reliable_above, color=P.GRID, alpha=0.6, linewidth=0)
        b.errorbar(N.centers, N.rate / w, yerr=N.rate_err / w, fmt="none", ecolor=P.SERIES[2], elinewidth=0.6)
        b.stairs(N.rate / w, N.edges_kev, color=P.SERIES[2], linewidth=1.0, label="net = source - background")
        b.set_yscale("symlog", linthresh=0.01)
        b.set_xlabel("energy (keV)")
        b.set_ylabel("counts / s / keV")
        b.set_title("Background-subtracted spectrum")
        for lab, e in [("Co-60", 1173.2), ("Co-60", 1332.5), ("K-40", 1460.8), ("sum", 2505.7)]:
            b.axvline(e, color=P.MUTED, linewidth=0.6, linestyle=":")
        b.legend()
    P.save(fig, out / "spectrum.png")

    if N is not None:
        fig, ax = P.figure(1, 1, height_per_row=3.6)
        a = ax[0, 0]
        m = (S.centers > 900) & (S.centers < 1800)
        a.stairs((S.rate / w)[m[:]][:-1] if False else (S.rate / w)[m], S.edges_kev[np.r_[m, False] | np.r_[False, m]],
                 color=P.SERIES[0], linewidth=1.1, label="source")
        a.stairs((B.rate / w)[m], B.edges_kev[np.r_[m, False] | np.r_[False, m]], color=P.SERIES[1], linewidth=1.1,
                 label="background (live-time scaled)")
        a.errorbar(N.centers[m], (N.rate / w)[m], yerr=(N.rate_err / w)[m], fmt="o", markersize=2.5,
                   color=P.SERIES[2], elinewidth=0.6, label="net")
        a.axhline(0, color=P.AXIS, linewidth=0.8)
        for e in (1173.2, 1332.5, 1460.8):
            a.axvline(e, color=P.MUTED, linewidth=0.6, linestyle=":")
        a.set_xlabel("energy (keV)")
        a.set_ylabel("counts / s / keV")
        a.set_title("⁶⁰Co and ⁴⁰K region: ⁴⁰K must vanish after subtraction")
        a.legend()
        P.save(fig, out / "spectrum_zoom.png")
