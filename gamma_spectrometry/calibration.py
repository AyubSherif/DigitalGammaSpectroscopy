"""Energy calibration and gain matching (step 7).

Each run is calibrated on its own against known gamma lines, so two runs with different
gain end up on the same energy axis ("gain matching" by calibration):

1. **Anchor**: a first gain estimate (keV per volt) from a line that is easy to recognise:
   the ⁶⁰Co pair (two strong peaks with centre ratio 1332.5/1173.2 = 1.1357) for a source
   run, or the strongest peak (⁴⁰K, 1460.8 keV) for a background run.
2. **Line search**: every line of the library is predicted at E / gain; lines closer than
   ``group_rel`` are fitted together (Gaussians on a linear background); a fitted peak is
   accepted if its centre is within ``accept_rel`` of the prediction, its FWHM is plausible
   for NaI(Tl) and its area is significant.
3. **Calibration curve**: E = c0 + c1 A (+ c2 A²) by weighted least squares through the
   accepted (amplitude, energy) points. The quadratic term is used only with >= 4 points
   and if it lowers the reduced chi² clearly.

Line energies (keV) from standard nuclear data tables.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .spectrum import fit_peaks

# (label, energy keV)
SOURCE_LINES = [("Co-60", 1173.228), ("Co-60", 1332.492), ("K-40", 1460.820),
                ("Co-60 sum", 2505.692), ("Tl-208", 2614.511)]
BACKGROUND_LINES = [("Pb-214", 351.932), ("Tl-208", 583.187), ("Bi-214", 609.312), ("Ac-228", 911.204),
                    ("K-40", 1460.820), ("Bi-214", 1764.494), ("Tl-208", 2614.511)]
CO60_RATIO = 1332.492 / 1173.228


@dataclass
class LinePoint:
    label: str
    energy_kev: float
    amplitude_v: float
    amplitude_err_v: float
    fwhm_pct: float
    area: float
    accepted: bool
    note: str = ""


@dataclass
class Calibration:
    coeffs: list                      # E(keV) = sum c_i A^i
    points: list = field(default_factory=list)
    chi2_ndf: float = float("nan")
    kind: str = ""

    def energy(self, amplitude_v) -> np.ndarray:
        return np.polynomial.polynomial.polyval(np.asarray(amplitude_v, dtype=float), self.coeffs)

    @property
    def gain_kev_per_v(self) -> float:
        return float(self.coeffs[1])

    def residuals_kev(self) -> list:
        return [(p.label, p.energy_kev, float(self.energy(p.amplitude_v) - p.energy_kev))
                for p in self.points if p.accepted]

    def to_dict(self) -> dict:
        d = asdict(self)
        d["residuals_kev"] = self.residuals_kev()
        return d

    def save(self, path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path) -> "Calibration":
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls(d["coeffs"], [LinePoint(**p) for p in d["points"]], d["chi2_ndf"], d.get("kind", ""))


# ----------------------------------------------------------------------------------------------

def _peaks_in(a: np.ndarray, lo: float, hi: float, bin_width: float, n: int):
    """The n most prominent peaks of the histogram of a in [lo, hi)."""
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    edges = np.arange(lo, hi, bin_width)
    h, _ = np.histogram(a, edges)
    hs = gaussian_filter1d(h.astype(float), 2)
    idx, prop = find_peaks(hs, prominence=0.05 * hs.max())
    order = np.argsort(prop["prominences"])[::-1][:n]
    return np.sort(0.5 * (edges[idx[order]] + edges[idx[order] + 1]))


def anchor_gain(amplitudes: np.ndarray, kind: str) -> float:
    """First estimate of keV per volt."""
    a = np.asarray(amplitudes)
    a = a[np.isfinite(a) & (a > 0)]
    # search between the median and the 99.5th percentile: the full-energy peaks lie above
    # the (low-energy dominated) median, and a few very large pulses (muons) must not
    # stretch the range so far that the peaks fall outside it or into one bin
    lo, top = np.quantile(a, 0.5), np.quantile(a, 0.995)
    bw = (top - lo) / 400
    if kind == "source":
        cand = _peaks_in(a, lo, top, bw, 6)
        best = None
        for i in range(cand.size):
            for j in range(i + 1, cand.size):
                r = cand[j] / cand[i]
                if abs(r / CO60_RATIO - 1) < 0.02 and (best is None or abs(r / CO60_RATIO - 1) < best[0]):
                    best = (abs(r / CO60_RATIO - 1), cand[i], cand[j])
        if best is None:
            raise ValueError("could not find the 60Co pair (two peaks with ratio 1.136)")
        return float(np.mean([1173.228 / best[1], 1332.492 / best[2]]))
    if kind == "background":
        cand = _peaks_in(a, lo, top, bw, 1)
        if cand.size == 0:
            raise ValueError("no peak found for the 40K anchor")
        return float(1460.820 / cand[0])
    raise ValueError(f"unknown kind {kind!r}")


#: NaI(Tl) resolution used for starting widths: FWHM = 6.5 % at 1173 keV, scaling as 1/sqrt(E)
FWHM_AT_1173 = 0.065


def expected_sigma_rel(energy_kev):
    return FWHM_AT_1173 / 2.3548 * np.sqrt(1173.0 / np.asarray(energy_kev, dtype=float))


def fit_lines(amplitudes: np.ndarray, lines, gain: float, group_rel: float = 0.16,
              accept_rel: float = 0.04, fwhm_range=(2.5, 15.0), min_significance: float = 4.0):
    """Fit every library line near its predicted amplitude; returns a list of LinePoint.

    Lines closer than ``group_rel`` are fitted together. Each Gaussian is bounded: centre
    within ``accept_rel`` of its prediction, width within 0.5-2x the expected NaI(Tl)
    resolution, area >= 0. A line is accepted if its area is ``min_significance`` sigma
    above zero and its FWHM is within ``fwhm_range`` (%)."""
    a = np.asarray(amplitudes)
    lines = sorted(lines, key=lambda x: x[1])
    groups, cur = [], [lines[0]]
    for ln in lines[1:]:
        if ln[1] / cur[-1][1] - 1 < group_rel:
            cur.append(ln)
        else:
            groups.append(cur)
            cur = [ln]
    groups.append(cur)
    points = []
    for g in groups:
        # Fit the group; if a line is rejected, drop the least significant rejected line and
        # refit the others (a weak line next to a stronger one makes the joint fit degenerate).
        active = list(g)
        dropped = {}
        while active:
            res = _fit_group(a, active, gain, accept_rel, fwhm_range, min_significance)
            bad = [p for p in res if not p.accepted]
            if not bad or len(active) == 1:
                break
            # A line whose centre or width ran to its bound was pulled onto a neighbour:
            # drop those first; among the rest, the least significant.
            def badness(p):
                at_bound = p.note.startswith("centre") or p.note.startswith("width pinned")
                sig = p.area / max(p.area_err, 1e-30) if np.isfinite(p.area_err) else -1.0
                return (0 if at_bound else 1, sig)

            worst = min(bad, key=badness)
            dropped[(worst.label, worst.energy_kev)] = worst
            active = [ln for ln in active if ln != (worst.label, worst.energy_kev)]
        final = {(p.label, p.energy_kev): p for p in res}
        for ln in g:
            p = final.get(ln) or dropped[ln]
            points.append(LinePoint(p.label, p.energy_kev, p.amplitude_v, p.amplitude_err_v, p.fwhm_pct,
                                    p.area, p.accepted, p.note))
    return points


@dataclass
class _GroupLine:
    label: str
    energy_kev: float
    amplitude_v: float
    amplitude_err_v: float
    fwhm_pct: float
    area: float
    area_err: float
    accepted: bool
    note: str


def _fit_group(a, group, gain, accept_rel, fwhm_range, min_significance):
    energies = np.array([e for _, e in group])
    pred = energies / gain
    s0 = pred * expected_sigma_rel(energies)
    lo, hi = pred.min() - 4 * s0[0], pred.max() + 4 * s0[-1]
    bw = pred.mean() * 0.004
    n_in = int(np.sum((a >= lo) & (a < hi)))
    if n_in < 50:
        return [_GroupLine(lab, e, np.nan, np.nan, np.nan, 0.0, np.inf, False, f"only {n_in} events in the window")
                for lab, e in group]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pf = fit_peaks(a, pred, lo, hi, bw, sigma0=s0, center_tol=accept_rel, sigma_range=(0.5, 2.0))
    except Exception as exc:  # noqa: BLE001
        return [_GroupLine(lab, e, np.nan, np.nan, np.nan, 0.0, np.inf, False, f"fit failed: {exc}")
                for lab, e in group]
    out = []
    for k, (lab, e) in enumerate(group):
        c, ce, fw, area = pf.centers[k], pf.center_errs[k], 100 * pf.resolution[k], pf.areas[k]
        ae = pf.area_errs[k]
        area_err = ae if np.isfinite(ae) and ae > 0 else np.inf
        note, ok = "", True
        if not np.isfinite(c) or abs(c / pred[k] - 1) > 0.98 * accept_rel:
            ok, note = False, f"centre {c:.4f} V at the edge of its allowed range ({100 * (c / pred[k] - 1):+.1f} %)"
        elif not (fwhm_range[0] <= fw <= fwhm_range[1]):
            ok, note = False, f"FWHM {fw:.1f} % outside {fwhm_range}"
        elif not (0.52 * s0[k] < pf.sigmas[k] < 1.96 * s0[k]):
            ok, note = False, f"width pinned at its bound (FWHM {fw:.1f} %): not a measured peak"
        elif area / area_err < min_significance:
            ok, note = False, f"not significant (area {area:.0f} +- {area_err:.0f})"
        elif not np.isfinite(ce) or ce <= 0:
            ok, note = False, "no centre uncertainty"
        out.append(_GroupLine(lab, e, float(c), float(ce), float(fw), float(area), float(area_err), ok, note))
    return out


def fit_calibration(points, max_degree: int = 2) -> Calibration:
    """Weighted polynomial E(A) through the accepted points (and through E=0 at A=0 as a
    weak constraint when only two points exist)."""
    acc = [p for p in points if p.accepted]
    if len(acc) < 2:
        raise ValueError(f"need at least 2 calibration lines, have {len(acc)}")
    A = np.array([p.amplitude_v for p in acc])
    E = np.array([p.energy_kev for p in acc])
    # energy uncertainty from the centre uncertainty (through the local gain), plus 0.1 % systematic
    gain = np.polyfit(A, E, 1)[0]
    sE = np.hypot(np.array([p.amplitude_err_v for p in acc]) * gain, 1e-3 * E)

    def wfit(deg):
        X = np.vander(A, deg + 1, increasing=True)
        if len(acc) == 2 and deg == 1:            # two points + origin pin: keeps an offset possible
            X = np.vstack([X, [1, 0]])
            y = np.r_[E, 0.0]
            w = np.r_[1 / sE, 1 / (0.05 * E.min())]
        else:
            y, w = E, 1 / sE
        c, *_ = np.linalg.lstsq(X * w[:, None], y * w, rcond=None)
        res = (np.vander(A, deg + 1, increasing=True) @ c - E) / sE
        ndf = max(len(acc) - (deg + 1), 1)
        return c, float(np.sum(res ** 2) / ndf)

    c1, chi1 = wfit(1)
    best_c, best_chi = c1, chi1
    if max_degree >= 2 and len(acc) >= 4:
        c2, chi2 = wfit(2)
        if chi2 < 0.5 * chi1:
            best_c, best_chi = c2, chi2
    return Calibration(list(map(float, best_c)), list(points), best_chi)


def gain_match(points, reference: Calibration) -> Calibration:
    """Calibration with the *shape* of ``reference`` and one free gain factor k:
    E(A) = E_ref(k A). Used when a run has too few lines of its own (same detector and
    electronics, different gain). k is the weighted mean, over the accepted lines, of the
    factor that maps each line exactly."""
    from scipy.optimize import brentq
    acc = [p for p in points if p.accepted]
    if not acc:
        raise ValueError("no accepted line to gain-match on")
    ks, ws = [], []
    for p in acc:
        ks.append(brentq(lambda kk: float(reference.energy(kk * p.amplitude_v)) - p.energy_kev, 0.2, 5.0))
        ws.append((p.amplitude_v / p.amplitude_err_v) ** 2)
    k = float(np.average(ks, weights=ws))
    coeffs = [c * k ** i for i, c in enumerate(reference.coeffs)]     # E_ref(k A) in powers of A
    return Calibration(coeffs, list(points), float("nan"), kind="gain-matched")


def calibrate(amplitudes: np.ndarray, kind: str, lines=None, iterations: int = 2,
              reference: Calibration | None = None) -> Calibration:
    """Anchor -> line search -> calibration curve; iterated once with the improved gain.
    With only one usable line, ``reference`` (another run's calibration) is gain-matched."""
    lines = lines or (SOURCE_LINES if kind == "source" else BACKGROUND_LINES)
    gain = anchor_gain(amplitudes, kind)
    cal = None
    for _ in range(iterations):
        pts = fit_lines(amplitudes, lines, gain)
        n_ok = sum(p.accepted for p in pts)
        if n_ok >= 2:
            cal = fit_calibration(pts)
        elif reference is not None and n_ok == 1:
            cal = gain_match(pts, reference)
        else:
            raise ValueError(f"{kind} run: {n_ok} usable calibration line(s)"
                             + ("" if reference is not None else "; give a reference calibration to gain-match"))
        gain = float(np.median([p.energy_kev / p.amplitude_v for p in pts if p.accepted]))
    cal.kind = kind if not cal.kind else f"{kind}, {cal.kind}"
    return cal


def gain_vs_time(events: np.ndarray, kind: str, chunk_captures: int = 20):
    """Relative gain per chunk of captures, from the anchor line(s) fitted in each chunk."""
    lines = [("Co-60", 1173.228), ("Co-60", 1332.492)] if kind == "source" else [("K-40", 1460.820)]
    all_gain = anchor_gain(events["amplitude"], kind)
    caps = np.unique(events["capture"])
    out = []
    for c0 in range(int(caps.min()), int(caps.max()) + 1, chunk_captures):
        sel = (events["capture"] >= c0) & (events["capture"] < c0 + chunk_captures)
        pts = fit_lines(events["amplitude"][sel], lines, all_gain, min_significance=3.0)
        acc = [p for p in pts if p.accepted]
        if not acc:
            continue
        g = np.mean([p.energy_kev / p.amplitude_v for p in acc])
        ge = g * np.sqrt(np.mean([(p.amplitude_err_v / p.amplitude_v) ** 2 for p in acc]) / len(acc))
        t = float(np.median(events["time_s"][sel]))
        out.append(dict(first_capture=c0, time_s=t, gain_kev_per_v=float(g), gain_err=float(ge)))
    if out:
        ref = np.average([o["gain_kev_per_v"] for o in out], weights=[1 / o["gain_err"] ** 2 for o in out])
        for o in out:
            o["relative_pct"] = 100 * (o["gain_kev_per_v"] / ref - 1)
            o["relative_err_pct"] = 100 * o["gain_err"] / ref
    return out
