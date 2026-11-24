"""Build and execute docs/walkthrough.ipynb: every function of the pipeline, in the order
it is called, run on the real 2021 data with its input, output and the maths behind it.

    python docs/build_walkthrough.py            # needs data/ (stores, template) to exist

The executed notebook (with outputs) is what GitHub shows.
"""

from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
cells = []


def md(s):
    cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))


def code(s):
    cells.append(nbf.v4.new_code_cell(s.strip("\n")))


# =============================================================================================
md(r"""
# gammaspec walkthrough: every function, in the order it is called

This notebook runs the pipeline step by step on the real 2021 NaI(Tl) data and shows, for
each function, **what goes in, what comes out, and how it works** (including the maths).

The order follows the data:

| Part | Step | What happens |
|---|---|---|
| A | 1. Ingestion | ASCII capture files → lossless integer codes → HDF5 store |
| B | 2. Characterisation | pulse template, analytic pulse model, noise |
| C | 3. Pre-processing | interference removal → pole-zero correction → baseline tracker |
| D | 4. Detection and events | optimal filter → successive cancellation → event list |
| E | 5–6. Pile-up resolution | TM, HYPER and least-squares fits on a real cluster; the pile-up solver in the pipeline |
| G | 7. Calibration and background subtraction | keV axis for each run, gain matching, net spectrum |
| H | Benchmark | the same pulses through a naive `find_peaks` count, an emulated analog chain and gammaspec |
| I | 9. Configuration | one run file selects input, stages and parameters; the runner executes it |
| F | Evaluation tools | how the numbers in the docs were measured |

Conventions: time is in **samples** (1 sample = 1/1.8 MHz = 0.556 µs) unless stated;
voltages are in volts; the raw pulses are negative and are flipped to positive early on.

*Rebuild with `python docs/build_walkthrough.py` (needs the local `data/` folder).*
""")

code(r"""
import numpy as np, matplotlib.pyplot as plt, json, tempfile
from pathlib import Path
import gammaspec
from gammaspec import rawio, store, signal, template, noise, interference, baseline, preprocess, detect, events, pipeline, spectrum
from gammaspec import plotting as P
%matplotlib inline
np.set_printoptions(precision=5, suppress=False, linewidth=110, threshold=24, edgeitems=4)
ROOT = Path.cwd().parent if Path.cwd().name == "docs" else Path.cwd()
RAW = ROOT / "ADC_raw_output"
SRC_STORE = ROOT / "data/source_first100.h5"
BKG_STORE = ROOT / "data/background.h5"
TEMPLATE = ROOT / "data/characterization/background/template.h5"
TMP = Path(tempfile.mkdtemp())
FS = 1.8e6
print("gammaspec", gammaspec.__version__)
""")

# =============================================================================================
md(r"""
---
# Part A: Ingestion (step 1)

## A1. `rawio.read_header(path)`

**Called by:** `build_index` (for every file) and `read_raw_file`.
**Input:** a capture file path. **Output:** a `RawHeader`.

It reads only the first three lines (`Creation`, `Signals`, `Dt`), splits each at the first
colon, and converts: the timestamp with `datetime.fromisoformat`, the sample count as an
integer, and `Dt` with its unit (`ms` → × 10⁻³ s). It also records how many bytes the header
occupies, so the samples can be read starting right after it.
""")
code(r"""
path = RAW / "raw_data_0.txt"
print("input: first 3 lines of", path.name)
with open(path, "rb") as f:
    for _ in range(3): print("   ", f.readline())
h = rawio.read_header(path)
print("\noutput:", h)
print("duration_s =", h.duration_s, "  time_s =", h.time_s)
""")

md(r"""
## A2. `rawio.parse_samples(body)` → `_parse_fixed_width(body)`

**Input:** the bytes after the header. **Output:** float64 volts.

Every sample line has exactly 14 bytes: `+d.dddddde-dd\n`. Instead of parsing 524,288
strings one at a time, the bytes are viewed as a **matrix with 14 columns** (one row per
line) and decoded column-wise with array arithmetic:

$$\text{mantissa} = \sum_{i=0}^{6} d_i\,10^{6-i},\qquad
\text{exponent} = \pm(10e_1 + e_2) - 6,\qquad v = \pm\,\text{mantissa}\cdot 10^{\text{exponent}}$$

where $d_i$ are the 7 mantissa digits (column 1, then columns 3–8; column 2 is the `.`).
Digits are bytes minus `ord('0')` = 48. The function first checks every fixed column
(`.` at 2, `e` at 9, newline at 13, signs, digits). If anything does not match, it returns
`None`, and `parse_samples` falls back to a generic float parser.
""")
code(r"""
raw = path.read_bytes()
body = raw[h.header_bytes:]
a = np.frombuffer(body, dtype=np.uint8).reshape(-1, 14)
print("input bytes, first 2 lines:", body[:28])
print("as a 14-column uint8 matrix (row 0):", a[0], "->", bytes(a[0]))
d = a[0, [1, 3, 4, 5, 6, 7, 8]].astype(int) - 48
print("mantissa digits:", d, " exponent sign:", chr(a[0, 10]), " exponent digits:", a[0, 11:13] - 48)
v = rawio.parse_samples(body)
print("\noutput volts:", v, " count:", v.size)
""")

md(r"""
## A3. `rawio.volts_to_codes(volts, lsb)`

**Input:** volts and the ADC step (LSB = 84.771697 µV). **Output:** int32 ADC codes.

$$\text{code} = \operatorname{rint}(v/\text{LSB}),\qquad
\text{residual} = \max|v/\text{LSB} - \text{code}|$$

The printed values are codes × LSB rounded to 7 significant digits, so the residual is
tiny (≤ 0.006 LSB). If it exceeds 0.02 LSB the file is not on this ADC grid and the
function raises an error. This check proves the conversion loses nothing.
""")
code(r"""
codes = rawio.volts_to_codes(v, rawio.DEFAULT_LSB_VOLTS)
x = v / rawio.DEFAULT_LSB_VOLTS
print("input  volts :", v[:6])
print("v / LSB      :", x[:6])
print("output codes :", codes[:6], codes.dtype)
print("max residual :", np.max(np.abs(x - codes)), "LSB   (code range", codes.min(), "to", codes.max(), ")")
""")

md(r"""
## A4. `rawio.read_raw_file(path, lsb)`
Combines A1–A3: header → seek past it → parse → check the count → codes (or volts if
`lsb=None`). This is what each conversion worker runs.
""")
code(r"""
h2, c2 = rawio.read_raw_file(path)
print("output:", h2.n_samples, "samples; codes", c2[:5], "; identical to A3:", np.array_equal(c2, codes))
""")

md(r"""
## A5. `rawio.estimate_lsb(volts)` (for unknown hardware)

**Input:** volts. **Output:** the ADC step.

1. Near zero the printed values keep full precision, so the spacing between the smallest
   distinct |values| gives a first guess.
2. Refine with least squares over progressively larger values (bound × 8 each pass):
   $$\text{LSB} \leftarrow \frac{\sum_i v_i k_i}{\sum_i k_i^2},\qquad k_i = \operatorname{rint}(v_i/\text{LSB})$$
   Each pass only includes values where rounding to the right code is still unambiguous.
""")
code(r"""
print("estimated LSB:", rawio.estimate_lsb(v), " built-in:", rawio.DEFAULT_LSB_VOLTS)
""")

md(r"""
## A6. `rawio.build_index(directory)` and `write_index_csv`

**Input:** a folder. **Output:** all headers **sorted by timestamp** (the file numbers in
the 2021 data are scrambled relative to time; see ADR-0002).
""")
code(r"""
idx = rawio.build_index(RAW)
print("captures:", len(idx))
for hh in idx[:5]:
    print(f"   file {hh.file_number:5d}  created {hh.created}")
print("file numbers of the first 8 in time order:", [hh.file_number for hh in idx[:8]])
""")

md(r"""
## A7. `store.convert(headers, out_path)` and `store.CaptureStore`

**Input:** headers (time order). **Output:** one HDF5 file: `/codes` (int32, one row per
capture, gzip + shuffle), `/index/*` (time, file number, path) and attributes (LSB, Δt).
Parsing runs in parallel worker processes (`_read_codes`); one process writes.

`CaptureStore(path)` reads it back: `codes(i)`, `volts(i) = codes × LSB`, `live_time_s`.
(Here we open the store built earlier rather than converting again.)
""")
code(r"""
with store.CaptureStore(SRC_STORE) as s:
    print("captures:", len(s), " samples/capture:", s.n_samples, " LSB:", s.lsb_volts, " dt:", s.dt_s)
    print("live time:", s.live_time_s, "s")
    print("capture 0 source:", s.source[0], " file number:", s.file_number[0])
    print("codes(0)[:5] =", s.codes(0)[:5], " identical to the text file:", np.array_equal(s.codes(0), codes))
    print("volts(0)[:5] =", s.volts(0)[:5])
""")

# =============================================================================================
md(r"""
---
# Part B: Characterisation (step 2)

We use capture 1 of the ⁶⁰Co run as the running example.

## B1. `signal.to_positive(volts, polarity=-1)`
Multiplies by −1, so pulses point up. Then the capture median is subtracted (in
`characterize.scan`).
""")
code(r"""
with store.CaptureStore(SRC_STORE) as s:
    volts1 = s.volts(1)
x1 = signal.to_positive(volts1)
x1 = x1 - np.median(x1)
print("input min/max:", volts1.min(), volts1.max(), " -> output min/max:", x1.min(), x1.max())
""")

md(r"""
## B2. `signal.robust_sigma(x)` and `signal.noise_sigma(x)`

* `robust_sigma`: $\hat\sigma = 1.4826\cdot\operatorname{median}|x-\operatorname{median}(x)|$.
  For Gaussian noise the median absolute deviation is 0.6745 σ (1/0.6745 = 1.4826).
  Outliers such as pulses barely move a median.
* `noise_sigma`: robust σ of the **first differences** divided by √2. For white noise
  $\operatorname{Var}(x_n - x_{n-1}) = 2\sigma^2$. Differencing removes slow baseline changes
  and the undershoot, and pulses are rare outliers. It uses the central 131,072 samples.
""")
code(r"""
print("plain std           :", x1.std())
print("robust_sigma(x)     :", signal.robust_sigma(x1))
print("noise_sigma(x)      :", signal.noise_sigma(x1))
""")

md(r"""
## B3. `signal.edge_signal(x, gap=2, avg=3)`

$$d[n] = x[n] - \tfrac{1}{3}\left(x[n-4] + x[n-3] + x[n-2]\right)$$

This is the rise relative to a short average ending 2 samples earlier. The pulse rises in
1–2 samples, so a half-risen sample stays out of the reference. On a slowly falling tail
d is small or negative, so the trigger does not fire on decays. The average is computed
with a cumulative sum (cs[i+3] − cs[i]) / 3.
""")
code(r"""
sig1 = signal.noise_sigma(x1)
pk = signal.find_pulses(x1, sig1, 8)
p0 = int(pk[10])
e1 = signal.edge_signal(x1)
seg = slice(p0 - 6, p0 + 8)
print("x[n] around a pulse     :", np.round(x1[seg], 4))
print("edge_signal d[n] (same) :", np.round(e1[seg], 4))
""")

md(r"""
## B4. `signal.find_pulses(x, sigma, threshold)`

Fires where the edge signal exceeds `threshold` × its noise. Its noise is
$\sigma\sqrt{1 + 1/3}$ (one sample minus a 3-sample average). Each trigger is moved to
the highest sample within the next 4 samples, the peak.
""")
code(r"""
print("input: sigma =", sig1, " threshold = 8")
print("output: ", pk.size, "pulse peaks; first:", pk[:8])
print("rate:", pk.size / (x1.size / FS), "per second")
""")

md(r"""
## B5. `signal.isolated(peaks, pre, post, n)` and `signal.pulse_free_mask(n, peaks, pre, post)`

* `isolated`: True where no other pulse lies within `pre` samples before or `post` after
  (and the window fits in the capture). Uses the gaps `diff(peaks)`.
* `pulse_free_mask`: True for samples at least `pre` before and `post` after every pulse.
  It marks +1 at each pulse's start of exclusion and −1 at its end, takes a cumulative sum,
  and keeps samples where the sum is 0. That is linear time, with no loops.
""")
code(r"""
iso = signal.isolated(pk, 300, 300, x1.size)
free = signal.pulse_free_mask(x1.size, pk, 20, 1500)
print("isolated (±300):", iso.sum(), "of", pk.size, "   pulse-free samples (20 before, 1500 after):", free.mean().round(3))
""")

md(r"""
## B6. `template.extract_windows(x, peaks, pre, post)`
One row per pulse: `x[p-pre : p+post+1]`, built with a single fancy-index (`peaks[:,None] + arange`).
""")
code(r"""
W = template.extract_windows(x1, pk[iso], 40, 300)
print("output shape:", W.shape, " (pulses × samples); row 0 around the peak:", np.round(W[0, 38:44], 4))
""")

md(r"""
## B7. `template.fit_pulses(windows, template, pre, ...)`

Fits each window as $w(j) = A\,s(j - \delta) + b$ over $j \in [-10, 60]$
(A = amplitude, δ = sub-sample peak position, b = baseline).

For a **fixed δ** the model is linear in A and b. Least squares has a closed form: with
$S_{ss}=\sum s^2,\ S_s=\sum s,\ S_{sy}=\sum s y,\ S_y=\sum y,\ n$ samples,

$$A = \frac{n S_{sy} - S_s S_y}{n S_{ss} - S_s^2},\qquad
b = \frac{S_{ss} S_y - S_s S_{sy}}{n S_{ss} - S_s^2},\qquad
\text{SSE} = \sum y^2 - A S_{sy} - b S_y$$

So: evaluate this for 33 trial δ ∈ [−1, 1] (all pulses at once with one matrix product
$Y S^T$), take the best δ, refine it with a parabola through the three SSE values around
the minimum, then solve exactly at the refined δ. With `fit_baseline=False` (used after
baseline restoration) b = 0 and $A = S_{sy}/S_{ss}$.
""")
code(r"""
T = template.PulseTemplate.load(TEMPLATE)
f = template.fit_pulses(W[:5].astype(float), T, 40)
print("input : 5 windows; peak samples", np.round(W[:5, 40], 4))
print("output: amplitude", np.round(f.amplitude, 4))
print("        shift δ  ", np.round(f.shift, 3), "(samples)")
print("        baseline ", np.round(f.baseline * 1e3, 3), "mV;  rms", np.round(f.rms * 1e3, 3), "mV")
""")

md(r"""
## B8. `template.build_template(windows, pre)`, with `_stack`, `_recentre`, `equalize_phase`

Pulses arrive at random sampling phases. Stacking many pulses at their *true* sub-sample
times samples the shape on a grid much finer than the ADC's (equivalent-time sampling).
Each iteration:

1. **`_stack`**: put every normalised sample $(w-b)/A$ at time $t = j - \delta$ into bins of
   1/16 sample; weighted average (weight $A^2$, the inverse variance of a normalised
   sample); light smoothing (none within ±3 samples of the peak); pre-pulse level set to 0.
2. **`_recentre`**: a parabola through the top three points finds the continuous maximum;
   shift it to t = 0 and scale it to 1.
3. **`fit_pulses`** with the new template; drop fits with rms > 3 × median (hidden pile-up).
4. **`equalize_phase`**: the true phases $\varphi = \delta \bmod 1$ are uniform, so map the
   fitted ones through their empirical CDF, $\varphi_i \to (\text{rank}_i + 0.5)/N$.
   Without this, a template warped inside each sample period plus correspondingly warped δ
   fit almost equally well (ADR-0004).

It stops when the template changes by less than 5 × 10⁻⁴ rms. Here we run it on the
isolated pulses of the background run (long isolation shows the whole undershoot).
""")
code(r"""
with store.CaptureStore(BKG_STORE) as s:
    wins = []
    for i in range(20, 50):
        xb = signal.to_positive(s.volts(i)); xb = xb - np.median(xb)
        sb = signal.noise_sigma(xb); pb = signal.find_pulses(xb, sb, 5)
        ok = signal.isolated(pb, 1500, 1500, xb.size) & (xb[pb] > 30 * sb) & (xb[pb] < 4)
        wins.append(template.extract_windows(xb, pb[ok], 40, 1500))
Wb = np.concatenate(wins)
print("input:", Wb.shape[0], "isolated pulses, windows of", Wb.shape[1], "samples")
Tb, fb, used = template.build_template(Wb, 40, resolution=16, verbose=True)
print("output: template with", Tb.shape.size, "points at 1/16 sample, t from", Tb.t_start, "to", round(Tb.t_end, 2))
print("template at t = -2, -1, 0, 1, 5, 46 samples:", np.round(Tb(np.array([-2., -1, 0, 1, 5, 46])), 4))
""")
code(r"""
phi_before = np.mod(np.array([0.10, 0.15, 0.18, 0.5, 0.52, 0.9]), 1)
print("equalize_phase example: fitted phases", phi_before, "->", np.round(template.equalize_phase(phi_before), 3))
fig, ax = P.figure(1, 2, height_per_row=3.0)
ax[0,0].plot(Tb.t[(Tb.t>-4)&(Tb.t<30)], Tb.shape[(Tb.t>-4)&(Tb.t<30)]); ax[0,0].set_title("Template: main lobe"); ax[0,0].set_xlabel("samples")
ax[0,1].plot(Tb.t, Tb.shape); ax[0,1].set_ylim(-0.06, 0.03); ax[0,1].set_title("Template: undershoot"); ax[0,1].set_xlabel("samples")
plt.show()
""")

md(r"""
## B9. `template.fit_multiexp(template)` with `multiexp_model` and `_exp_gauss`

The analytic model is a sum of exponentials starting at $t_0$, smoothed by a Gaussian of
width σ (the electronics bandwidth):

$$s(t) = \sum_k a_k\, \big[e^{-(t-t_0)/\tau_k} H(t-t_0)\big] * \mathcal N(0,\sigma),\qquad \sum_k a_k = 0$$

(the sum of amplitudes is 0, so the pulse starts from zero). Each smoothed exponential is
closed form (`_exp_gauss`). With $x = t - t_0$ and $z = (\sigma/\tau - x/\sigma)/\sqrt2$:

$$e^{-x/\tau} * \mathcal N = \tfrac12 e^{\sigma^2/2\tau^2 - x/\tau}\operatorname{erfc}(z)
= \tfrac12 e^{-x^2/2\sigma^2}\operatorname{erfcx}(z)$$

The second form is used for $z \ge 0$ to avoid overflow. `fit_multiexp` fits
$t_0, \sigma, \log\tau_k, a_k$ by nonlinear least squares (`scipy.optimize.least_squares`).
""")
code(r"""
m3 = template.fit_multiexp(Tb)
m4 = template.fit_multiexp(Tb, taus0=tuple(m3.taus) + (5 * m3.taus[-1],), amps0=tuple(m3.amps[1:]) + (0.005,))
for name, m in [("3 exponentials", m3), ("4 exponentials", m4)]:
    print(f"{name}: taus = {np.round(m.taus, 2)} samples = {np.round(m.taus / 1.8, 3)} µs; amps = {np.round(m.amps, 4)}; sigma = {m.sigma:.3f}; rms = {m.rms:.2e}")
""")

md(r"""
## B10. `noise.pulse_free_windows(x, free, 512)` and `noise.noise_spectrum(windows, fs)`

* `pulse_free_windows` finds runs of the free mask (edges of `diff(mask)`) and cuts
  non-overlapping 512-sample windows from runs that are long enough.
* `noise_spectrum`: subtract each window's mean, apply a Hann window $w$, FFT, average:
  $$\text{PSD}(f) = \frac{2\,\overline{|X(f)|^2}}{f_s\sum w^2}\quad[\mathrm{V^2/Hz}]$$
  (×2 for a one-sided spectrum). The autocorrelation is the inverse FFT of the zero-padded
  power spectrum, divided by (n − lag) (unbiased) and normalised to 1 at lag 0.
""")
code(r"""
wnd, starts = noise.pulse_free_windows(x1, free, 512)
ns = noise.noise_spectrum(wnd, FS)
print("input:", wnd.shape, "windows;  output: sigma", round(ns.sigma_v * 1e3, 3), "mV, white floor", f"{ns.white_level:.2e}", "V²/Hz")
top = np.argsort(ns.psd)[::-1][:5]
print("strongest PSD bins (kHz):", np.round(np.sort(ns.freq_hz[top]) / 1e3, 1))
print("autocorrelation lags 0-4:", np.round(ns.autocorr[:5], 3))
""")

md(r"""
## B11. `characterize.characterize(store, out_dir)` (the step-2 driver)

Calls `scan` (B1–B6 for every capture), `build_template` (B8), `fit_multiexp` (B9),
`shape_by_amplitude` (templates per amplitude bin with the same alignment, compared to a
reference stacked from all pulses), `noise_spectrum` (B10, split into quiet and noisy
captures), a template fit of all pulses with a clean main lobe (preliminary spectrum), and
writes `template.h5`, `noise.h5`, `summary.json` and the plots. The results are in
`docs/characterization.md`.
""")

# =============================================================================================
md(r"""
---
# Part C: Pre-processing (step 3)

The chain is `preprocess.Preprocessor.process(volts)`:
polarity + offset → **line removal** (C1–C3) → **pole-zero** (C4–C6) → **baseline tracker** (C7).

## C0. Polarity and offset
`x = to_positive(volts)`, then subtract the **mean**. The coupled pulses have zero net
area, so the mean of the signal equals the electronic offset (a median would be biased by
the undershoots).
""")
code(r"""
xs = signal.to_positive(volts1); off = xs.mean(); xs = xs - off
print("offset (mean):", off, " vs median:", np.median(signal.to_positive(volts1)))
""")

md(r"""
## C1. `interference.detect_lines(x, free, fs, cfg)`

1. Zero out the pulses (mask from a coarse trigger, 20 before / 60 after each), apply a Hann
   window, take the FFT of the **whole capture** (resolution 3.4 Hz), power P(f).
2. Local floor: running median over 1001 bins. A line is a bin with P > 100 × floor that is
   the largest within ±3.5 kHz (2 × fs/block).
3. Refine each frequency with a parabola through log P at the three bins around the peak.
4. Merge lines closer than 3.5 kHz (keep the strongest). At most 40.
""")
code(r"""
sg = signal.noise_sigma(xs)
fr_mask = signal.pulse_free_mask(xs.size, signal.find_pulses(xs, sg, 5), 20, 60)
cfgL = interference.LineRemovalConfig()
lines = interference.detect_lines(xs, fr_mask, FS, cfgL)
print("output:", lines.size, "lines (kHz):", np.round(lines / 1e3, 2))
print("multiples of 9 kHz (= fs/200)?", np.round(lines / 9e3, 2))
""")

md(r"""
## C2. `interference.subtract_lines(x, free, freqs, fs, cfg)`

In blocks of L = 1024 samples with 50 % overlap, model the interference as
$\sum_k \alpha_k\cos(2\pi f_k t) + \beta_k\sin(2\pi f_k t)$ (t = time within the block).
This is linear: $y \approx X c$ with a fixed design matrix X (L × 2K). The fit must use only
pulse-free samples. Instead of solving a different masked problem per block, fill the masked
samples with the current model and refit with the same X (an EM iteration):

$$c \leftarrow (X^TX)^{-1}X^T \tilde y,\qquad \tilde y = m\,y + (1-m)\,Xc$$

All blocks are solved at once ($\tilde Y\,[(X^TX)^{-1}X^T]^T$). Three passes converge.
The blocks are overlap-added with $\sin^2$ windows (they sum to 1), and the model is
subtracted from **every** sample. Inside a pulse the interference is interpolated, so the
pulse itself is untouched.
""")
code(r"""
model = interference.subtract_lines(xs, fr_mask, lines, FS, cfgL)
xl = xs - model
print("noise sigma before:", round(signal.noise_sigma(xs) * 1e3, 3), "mV  after:", round(signal.noise_sigma(xl) * 1e3, 3), "mV")
seg = slice(200000, 200300)
fig, ax = P.figure(1, 1, height_per_row=2.8)
ax[0,0].plot(xs[seg] * 1e3, label="before"); ax[0,0].plot(xl[seg] * 1e3, label="after line removal"); ax[0,0].plot(model[seg] * 1e3, label="subtracted model", linewidth=1)
ax[0,0].set_ylim(-8, 8); ax[0,0].set_ylabel("mV"); ax[0,0].legend(); plt.show()
""")

md(r"""
## C3. `interference.LineRemoval.process(x, free)`
C1 then C2 for one capture; `info` holds the lines and the rms removed.
""")

md(r"""
## C4. `baseline.design_pole_zero(taus, amps)` (from `design_from_template`)

The sampled pulse model is $s[n] = \sum_i a_i\, p_i^{\,n}$ with $p_i = e^{-1/\tau_i}$. Its
z-transform is

$$S(z) = \sum_i \frac{a_i}{1-p_i z^{-1}} = \frac{N(z)}{\prod_i (1-p_i z^{-1})},\qquad
N(z) = \sum_i a_i \prod_{j\ne i}(1-p_j z^{-1})$$

The AC coupling appears twice: as slow **poles** $p_i$ (τ > 30 samples: the undershoot)
and as **zeros of N near z = 1** (they make the net area zero). The correction filter swaps
them:

$$G(z) = \frac{\prod_{\text{slow}}(1-p_i z^{-1})}{\prod_{\text{near 1}}(1-r_k z^{-1})}$$

so $G\cdot S$ has no undershoot. A root at or above 1 would be unstable and is clamped to
$1 - 1/20000$ (a slow leak).
""")
code(r"""
pz = baseline.design_from_template(T)
taus, amps = np.asarray(T.meta["model_taus"]), np.asarray(T.meta["model_amps"])
p = np.exp(-1 / taus)
N = np.zeros(len(p))
for i in range(len(p)):
    q = np.array([1.0])
    for j in range(len(p)):
        if j != i: q = np.convolve(q, [1, -p[j]])
    N = N + amps[i] * q
print("model taus (samples):", np.round(taus, 2), " poles p_i:", np.round(p, 5))
print("roots of N(z):", np.round(np.roots(N), 6))
print("output:", pz.describe())
print("   b =", pz.b, "\n   a =", pz.a)
""")

md(r"""
## C5. `baseline.PoleZero.process(x)`
Runs the IIR filter $y = G\,x$ (`scipy.signal.lfilter`) and keeps its internal state `zi`
between calls, so a capture processed in chunks gives exactly the same result as in one
piece (streaming-ready). It starts "at rest" on the first sample's level.
""")
code(r"""
pzs = baseline.PoleZero(pz)
xp = pzs.process(xl)
pz2 = baseline.PoleZero(pz)
chunks = np.concatenate([pz2.process(c) for c in np.array_split(xl, 7)])
print("chunked == whole:", np.allclose(chunks, xp))
p1 = int(pk[30]); seg = slice(p1 - 20, p1 + 400)
fig, ax = P.figure(1, 1, height_per_row=2.8)
ax[0,0].plot(xl[seg] * 1e3, label="after line removal"); ax[0,0].plot(xp[seg] * 1e3, label="after pole-zero")
ax[0,0].set_ylim(-30, 30); ax[0,0].set_ylabel("mV"); ax[0,0].set_title("Undershoot removed (zoomed on the baseline)"); ax[0,0].legend(); plt.show()
""")

md(r"""
## C6. `baseline.filter_template(template, b, a)`
The detector and the fit must use the pulse shape **as it looks after** pole-zero. The
fine-grid template is evaluated on each of the 16 sub-sample phases, each phase sequence
is filtered at the sample rate, and the results are interleaved back and renormalised.
""")
code(r"""
Tp = baseline.filter_template(T, pz.b, pz.a)
tt = np.array([0., 10, 30, 46, 100, 300])
print("t (samples)        :", tt)
print("template           :", np.round(T(tt), 4))
print("pole-zero template :", np.round(Tp(tt), 4))
""")

md(r"""
## C7. `baseline.track_baseline(x, free, cfg)` (inside `BaselineTracker.process`)
A centred moving average of the pulse-free samples (1001 samples; pulses masked 10 before
and 50 after). With cumulative sums of $x\cdot m$ and of $m$:

$$\hat b[i] = \frac{\sum_{j=i-500}^{i+500} x_j m_j}{\sum_{j=i-500}^{i+500} m_j}$$

computed in O(1) per sample. Windows with fewer than 50 free samples are bridged by
interpolation. The estimate is subtracted. Pole-zero boosts slow drift (×13 at 90 Hz), and
this removes it, so the two always go together (ADR-0005).
""")
code(r"""
sp = signal.noise_sigma(xp)
fm = signal.pulse_free_mask(xp.size, signal.find_pulses(xp, sp, 5), 10, 50)
b_est = baseline.track_baseline(xp, fm, baseline.TrackerConfig())
xc = xp - b_est
w = 1500  # warm-up: pole-zero starts without knowing the tails of pulses before the capture
print(f"baseline estimate in the warm-up region (first {w} samples): {b_est[:w].min()*1e3:.1f} to {b_est[:w].max()*1e3:.1f} mV")
print("baseline estimate at samples 1500, 10k, 30k, 100k, 300k (mV):", np.round(b_est[[1500, 10000, 30000, 100000, 300000]] * 1e3, 2))
print("-> pole-zero restores a slowly wandering DC level of a few mV (its low-frequency gain);")
print("   the tracker (1001-sample window) follows it, so the corrected signal stays at zero")
print("median of pulse-free samples before:", round(np.median(xp[w:][fm[w:]]) * 1e3, 3), "mV  after:", round(np.median(xc[w:][fm[w:]]) * 1e3, 3), "mV")
""")

md(r"""
## C8. `preprocess.Preprocessor.process(volts)` (the whole chain)
C0 → C1–C3 → C5 → C7 in one call, with every stage switchable in `PreprocessConfig`.
""")
code(r"""
pre = preprocess.Preprocessor(FS, preprocess.PreprocessConfig(), pz)
y1 = pre.process(volts1)
print("info:", {k: (np.round(v, 5) if np.ndim(v) == 0 else f"{len(v)} lines") for k, v in vars(pre.info).items()})
print("identical to the step-by-step result:", np.allclose(y1, xc))
""")

# =============================================================================================
md(r"""
---
# Part D: Detection and events (step 4)

## D1. `detect.kernel_from_template(T, pre=6, post=60)`
The pole-zero template sampled at integer offsets −6…60: the pulse shape **s** used by all filters.
""")
code(r"""
det = detect.Detector(Tp)
print("s =", np.round(det.s, 3), " length", det.s.size)
""")

md(r"""
## D2. `detect.noise_autocovariance(x, free, n_lags)`
The noise covariance from pulse-free samples only (gapped estimator), with $m$ the mask:

$$r[k] = \frac{\sum_i x_i m_i\, x_{i+k} m_{i+k}}{\sum_i m_i m_{i+k}}$$

Both sums are autocorrelations, computed with FFTs (zero-padded to avoid wrap-around) on
the central 131,072 samples.
""")
code(r"""
coarse = signal.find_pulses(y1, signal.noise_sigma(y1), 5)
freeD = signal.pulse_free_mask(y1.size, coarse, 20, 120)
r = detect.noise_autocovariance(y1, freeD, det.s.size)
print("r[0..5] (V²):", r[:6], "   -> correlation", np.round(r[:6] / r[0], 3))
""")

md(r"""
## D3. `detect.optimal_kernel(s, r)`

The best linear amplitude estimator in coloured noise (generalised least squares) is
$\hat A = h^T x / h^T s$ with

$$h = C^{-1}s,\qquad C_{ij} = r[|i-j|]$$

C is Toeplitz, so `scipy.linalg.solve_toeplitz` (Levinson recursion, O(L²)) solves it
directly. Its noise is $\sigma_{\hat A} = (s^T C^{-1} s)^{-1/2}$, the smallest possible.
""")
code(r"""
from scipy.linalg import toeplitz
h = detect.optimal_kernel(det.s, r)
C = toeplitz(r)
s_ = det.s
print("h (normalised) =", np.round(h / np.abs(h).max(), 2))
# noise of A = h.x / h.s is sqrt(h' C h) / (h.s): optimal -> (s' C^-1 s)^-1/2; matched (h = s) -> sqrt(s' C s) / (s.s)
print("theoretical output noise: optimal", round(1e3 / np.sqrt(s_ @ h), 3), "mV;  plain matched",
      round(1e3 * np.sqrt(s_ @ C @ s_) / (s_ @ s_), 3), "mV")
""")

md(r"""
## D4. `detect.apply_filter(x, h, s, pre)`

$$y[n] = \frac{\sum_k h[k]\, x[n+k-\text{pre}]}{\sum_k h[k]\, s[k]}$$

If x contains A·s with its peak at n₀, then y[n₀] = A exactly, for **any** kernel h. That
is why the matched filter (h = s) and the optimal filter (h = C⁻¹s) are on the same scale.
Implemented as an FFT correlation.
""")
code(r"""
yo = detect.apply_filter(y1, h, det.s, 6)
ym = detect.apply_filter(y1, det.s, det.s, 6)
q = int(pk[30])
print(f"pulse at {q}: peak sample {y1[q]*1e3:.2f} mV; optimal filter {yo[q]*1e3:.2f} mV; matched filter {ym[q]*1e3:.2f} mV")
""")

md(r"""
## D5. `detect.negative_side_sigma(y)` and `detect.local_sigma(y, 4096)`

Pulses make y positive, so the noise is measured from its negative half: for zero-mean
Gaussian noise $\operatorname{median}(|y|\,:\,y<0) = 0.6745\,\sigma$. `local_sigma` does
this in blocks of 4096 samples and interpolates between block centres, so the threshold
follows interference bursts.
""")
code(r"""
print("negative_side_sigma: optimal", detect.negative_side_sigma(yo) * 1e3, "mV;  matched", detect.negative_side_sigma(ym) * 1e3, "mV")
ls = detect.local_sigma(yo, 4096)
print("local sigma p10/p50/p90 (mV):", np.round(np.percentile(ls, [10, 50, 90]) * 1e3, 3))
""")

md(r"""
## D6. `Detector.response_table(h)`
The filter output produced by a unit pulse whose peak lies p/16 of a sample after a sample
(p = 0…15), over ±66 samples: 16 rows × 133 columns. It is computed directly as one matrix
product of sliding windows with h. Cancellation (D7) subtracts these exact responses.
""")
code(r"""
R = det.response_table(h)
print("shape:", R.shape, "  peak of phase 0:", R[0].max().round(4))
span = R.shape[1] // 2
print("response before the peak (phase 0, offsets -30…-2):", np.round(R[0, span-30:span-1], 3))
""")

md(r"""
## D7. `detect.successive_cancellation(y, R, thr, margin)`

The optimal filter's response has negative dips and small bumps (above), which a plain
peak search would report as extra pulses. Successive cancellation (CLEAN):

1. Take candidates (local maxima above the threshold) in **descending** order (`_cancel_pass`).
2. For each: re-check it in the **residual** output (it may have shrunk when a neighbour
   was subtracted). Accept only if residual ≥ threshold + 3 % × (largest accepted pulse
   whose response overlaps here), it is a local maximum, and it is **not within ±2 samples
   of an accepted pulse** (there, only template mismatch of the real pulse remains).
3. Parabola through the three samples → sub-sample position → phase row p of R.
   $A = y_r[j] / R_p[j-\text{base}]$; subtract $A\,R_p$ from the residual.
4. **Search the residual again** (up to 4 passes): a small pulse hidden under a larger
   pulse's negative dip is not a local maximum of the original output, but it appears once
   the larger pulse has been subtracted.
5. Return the accepted positions, amplitudes and the residual.

Below: a real close pile-up from this capture.
""")
code(r"""
thr, sigma_y = det._threshold(yo)
idxc, ampc, resid = detect.successive_cancellation(yo, R, thr, det.cfg.margin)
gaps = np.diff(idxc); k = int(np.flatnonzero((gaps >= 4) & (gaps <= 12))[0])
a0 = idxc[k] - 25
seg = slice(a0, a0 + 110)
print("close pair at", idxc[k], idxc[k + 1], "(gap", gaps[k], "samples), amplitudes", np.round(ampc[k:k+2] * 1e3, 1), "mV")
fig, ax = P.figure(1, 1, height_per_row=3.0)
t_ = np.arange(110)
ax[0,0].plot(t_, y1[seg] * 1e3, label="clean signal"); ax[0,0].plot(t_, yo[seg] * 1e3, label="optimal filter output")
ax[0,0].plot(t_, resid[seg] * 1e3, label="residual after cancellation")
for j in (k, k + 1): ax[0,0].axvline(idxc[j] - a0, color=P.MUTED, linewidth=0.8, linestyle=":")
ax[0,0].set_ylabel("mV"); ax[0,0].legend(); plt.show()
""")

md(r"""
## D8. `Detector.detect(x)` (two passes)
filter (D2–D4) → threshold from the output (D5) → cancellation (D7) → **re-estimate the
noise from the residual** (real pulses' dips had inflated it) → cancellation again with the
lower threshold.
""")
code(r"""
d1 = det.detect(y1)
print(f"output: {d1.index.size} pulses; filter noise {d1.sigma*1e3:.3f} mV (first pass {sigma_y*1e3:.3f}); median threshold {d1.threshold_v*1e3:.2f} mV")
print("first indices:", d1.index[:8], " amplitudes (mV):", np.round(d1.filter_amplitude[:8] * 1e3, 1))
""")

md(r"""
## D9. Alternatives: `detect.sharpen(x, tau)` and `decay_constant(T)`
`decay_constant` fits log s(t) over 4…25 samples (τ ≈ 7.7). `sharpen` computes
$x[n] - a\,x[n-1]$ with $a=e^{-1/\tau}$: an exponential decay becomes 0, so each pulse
collapses onto its rising edge (method `"sharpened"`).
""")
code(r"""
print("decay constant:", detect.decay_constant(Tp))
print("sharpened template around the peak:", np.round(detect.sharpen(Tp(np.arange(-3, 8.)), detect.decay_constant(Tp)), 3))
""")

md(r"""
## D10. `pipeline.Pipeline.process_capture(volts, capture, t0)`

For one capture: pre-process (C8) → detect (D8) → then

* **gaps** to the neighbours; **clusters**: a new cluster starts when the gap ≥ 60 samples
  (a template fit window is −10…+40, so closer pulses share one);
* **flags**: PILEUP (cluster > 1), WARMUP (first 1500 samples), EDGE, then the template fit
  (B7, baseline fixed at 0) of isolated pulses → amplitude, shift, rms; BAD_FIT if
  rms > 4σ + 2 % × A; LARGE (> 1.5 V) and AFTER_LARGE (200 samples after one);
* **snippets**: clean waveform of each pile-up cluster (−20…+80 samples), kept for re-analysis;
* **pile-up solver** (`_solve_clusters`, explained in E9–E10): amplitudes and sub-sample
  times for every pulse in a cluster (and for badly fitted isolated pulses), plus pulses
  found by its residual search (flag FIT_ADDED); gaps and cluster sizes recomputed;
* **time** = capture time + (index + shift) × Δt.
""")
code(r"""
pipe = pipeline.Pipeline(T, FS)
with store.CaptureStore(SRC_STORE) as s:
    ev, row, sn = pipe.process_capture(s.volts(1), 1, float(s.time_s[1]))
names = ["index", "amplitude", "filter_amplitude", "shift", "fit_rms", "flags", "cluster", "cluster_size", "gap_prev"]
print(" ".join(f"{n:>16s}" for n in names))
for e in ev[5:15]:
    print(" ".join(f"{e[n]:16.5g}" for n in names))
print("\ncapture row:", row[0])
print("snippets:", len(sn), "; first:", sn[0][0], "capture", sn[0][1], "start", sn[0][2], "length", sn[0][3].size)
""")

md(r"""
## D11. `pipeline.run(...)`, `events.EventWriter`, `events.EventFile`, `summarize`

`run` processes all captures (in parallel worker processes, `_worker_run`), makes cluster
ids unique, and appends each capture's rows with `EventWriter.append` to resizable,
compressed column datasets. `EventFile` reads everything back (`events`, `captures`,
`snippet(cluster)`, `clean()`, `live_time_s`). `summarize` reports rates and flags.
""")
code(r"""
r = pipeline.run(SRC_STORE, TEMPLATE, TMP / "ev.h5", limit=5, workers=1, verbose=False)
print("run:", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})
ef = events.EventFile(TMP / "ev.h5")
print("events:", ef.events.size, " live time:", round(ef.live_time_s, 3), "s  clean:", ef.clean().size, " snippets:", ef.n_snippets)
cap, start, w = ef.snippet(int(ef.events[(ef.events["flags"] & events.FLAG_PILEUP) != 0]["cluster"][0]))
print("a snippet:", cap, start, w.shape)
print(json.dumps(pipeline.summarize(TMP / "ev.h5"), indent=1)[:600])
""")

md(r"""
## D12. `spectrum.fit_peaks(values, centers, lo, hi, bin_width)`
Histogram, then fit N Gaussians on a linear background by weighted least squares
(`curve_fit`, Poisson weights $1/\sqrt{\max(n,1)}$). Resolution = FWHM / centre with
FWHM = 2√(2 ln 2) σ = 2.355 σ. Here: the ⁶⁰Co lines from the full 100-capture event file.
""")
code(r"""
efull = events.EventFile(ROOT / "data/events/source_first100.h5")
amp = efull.clean()["amplitude"]
pf = spectrum.fit_peaks(amp, [0.316, 0.357], 0.27, 0.40, 0.002)
print("centres (V):", np.round(pf.centers, 4), " FWHM:", np.round(100 * pf.resolution, 2), "%  ±", np.round(100 * pf.resolution_err, 2))
print("ratio of centres:", round(pf.centers[1] / pf.centers[0], 4), "(1332/1173 = 1.1355)")
edges = np.linspace(0, 0.8, 500)
fig, ax = P.figure(1, 1, height_per_row=3.2)
ax[0,0].stairs(np.histogram(amp, edges)[0], edges); ax[0,0].set_yscale("log"); ax[0,0].set_xlabel("amplitude (V)")
ax[0,0].set_title("Clean events, ⁶⁰Co run (29 s live): template-fit amplitudes"); plt.show()
""")

# =============================================================================================
md(r"""
---
# Part E: Pile-up resolution (step 5)

A pile-up cluster is a stretch of clean signal with two or more pulses closer than 60
samples. Every method in `gammaspec.pileup` takes the cluster waveform **y** (baseline 0),
the detected peak positions **pos** (relative to the start of y) and a `PileupContext`, and
returns one amplitude per pulse. We use a real ⁶⁰Co pile-up cluster from the event file.

## E0. A real cluster: `EventFile.snippet(cluster)`
""")
code(r"""
ev_all = efull.events
ok_flags = (ev_all["flags"] & (events.FLAG_WARMUP | events.FLAG_EDGE)) == 0
cand = ev_all[(ev_all["cluster_size"] == 2) & (ev_all["gap_next"] >= 6) & (ev_all["gap_next"] <= 10)
              & (ev_all["filter_amplitude"] > 0.25) & ok_flags]
cl = int(cand["cluster"][0])
cap, start, yc = efull.snippet(cl)
rows = ev_all[ev_all["cluster"] == cl]
pos = (rows["index"] - start).astype(float)
print(f"cluster {cl}: capture {cap}, snippet starts at sample {start}, length {yc.size}")
print("positions in the snippet:", pos, " detector amplitudes (mV):", np.round(rows["filter_amplitude"] * 1e3, 1))
fig, ax = P.figure(1, 1, height_per_row=2.8)
ax[0,0].plot(yc * 1e3, marker="o", markersize=2.5, linewidth=0.9)
for p_ in pos: ax[0,0].axvline(p_, color=P.MUTED, linestyle=":", linewidth=0.8)
ax[0,0].set_ylabel("mV"); ax[0,0].set_title("Clean waveform of the cluster (dotted: detected peaks)"); plt.show()
""")

md(r"""
## E1. `pileup.PileupContext`
Holds what the methods need: the (pole-zero-corrected) template **s**, its decay constant
τ (for HYPER), the noise autocovariance (for GLS/WLS), the template mismatch ε (for WLS),
and the timing-refinement settings. `derivative(t)` is the central difference
$s'(t) \approx [s(t+h) - s(t-h)]/2h$ with h = 1/16 sample. The noise autocovariance here is
measured from capture 1's pulse-free samples (as in D2), now with 4096 lags.
""")
code(r"""
from gammaspec import pileup
r_long = detect.noise_autocovariance(y1, freeD, 4096)
ctx = pileup.PileupContext(Tp, detect.decay_constant(Tp), r_long)
print("tau =", round(ctx.tau, 3), " model_error eps =", ctx.model_error, " noise variance r[0] =", r_long[0])
print("s(t) and s'(t) at t = -1, 0, 1, 5:", np.round(Tp(np.array([-1., 0, 1, 5])), 3), np.round(ctx.derivative(np.array([-1., 0, 1, 5])), 3))
""")

md(r"""
## E2. `pileup.tm(y, pos, ctx)`: template matching (the paper)

Measured amplitude of pulse p = the sample at its peak, $b_p = y[t_p]$. Each $b_p$ is the
pulse's own amplitude plus the tails of the **earlier** pulses at that point:

$$b_p = \sum_{q \le p} M_{pq}\, a_q,\qquad M_{pq} = s(t_p - t_q),\quad M_{pp} = s(0) = 1$$

M is lower triangular with a unit diagonal, so a is found by forward substitution
(`solve_triangular`). It uses **one sample per pulse**, and assumes the peak is exactly on
a sample.
""")
code(r"""
b_meas = yc[np.round(pos).astype(int)]
dd = pos[:, None] - pos[None, :]
M = np.where(dd >= 0, Tp(dd), 0.0)
print("b =", np.round(b_meas * 1e3, 2), "mV\nM =\n", np.round(M, 4))
print("a = M^-1 b =", np.round(pileup.tm(yc, pos, ctx) * 1e3, 2), "mV")
""")

md(r"""
## E3. `pileup.hyper(y, pos, ctx)`: HYPER

For mono-exponential pulses $A_q e^{-\lambda(t-t_q)}$, between pulse p and p+1:

$$D(t) = y(t) + \lambda\int_{t_p}^{t} y\,dt' = \sum_{q\le p} A_q e^{-\lambda(t_p - t_q)} = \text{const} = D_p$$

(differentiate: $\dot y = -\lambda y$ between pulses, so $\dot D = 0$). Hence
$D_p = A_p + D_{p-1}e^{-\lambda(t_p - t_{p-1})}$ and

$$A_p = D_p - D_{p-1}\,e^{-\lambda (t_p - t_{p-1})}$$

D_p is averaged over the segment (trapezoidal integral, from the peak, at most 3τ long).
Our pulses are not purely exponential (finite rise), which biases HYPER.
""")
code(r"""
lam = 1 / ctx.tau
p_ = np.round(pos).astype(int)
for k in range(p_.size):
    a_, b_ = p_[k], min(p_[k + 1] if k + 1 < p_.size else yc.size, p_[k] + int(3 * ctx.tau))
    seg_ = yc[a_:b_]; integ = np.concatenate([[0.0], np.cumsum(0.5 * (seg_[1:] + seg_[:-1]))])
    Dt = seg_ + lam * integ
    print(f"segment {k}: D(t) first values {np.round(Dt[:6] * 1e3, 1)} mV -> mean D_{k} = {Dt.mean() * 1e3:.2f} mV")
print("HYPER amplitudes:", np.round(pileup.hyper(yc, pos, ctx) * 1e3, 2), "mV")
""")

md(r"""
## E4. `pileup.ls(y, pos, ctx)`: least squares with integer times

Model the **whole** cluster: $y(t) = \sum_k A_k\, s(t - t_k) + \text{noise}$. With the
design matrix $S_{tk} = s(t - t_k)$ (one column per pulse) this is linear:

$$\hat A = \arg\min \lVert y - S A\rVert^2 = (S^TS)^{-1}S^Ty$$

(`np.linalg.lstsq`). It uses all samples, so it averages noise, but it still assumes the
peaks are on samples.
""")
code(r"""
t_ = np.arange(yc.size, dtype=float)
S = Tp(t_[:, None] - pos[None, :])
print("design matrix S:", S.shape, "  S^T S =\n", np.round(S.T @ S, 3))
print("LS amplitudes:", np.round(pileup.ls(yc, pos, ctx) * 1e3, 2), "mV")
""")

md(r"""
## E5. `pileup.ls_refine`: least squares with sub-sample timing (`_ls(refine=True)`)

The arrival times are refined jointly with the amplitudes (Gauss–Newton). For a small
time shift δ:

$$A\,s(t - t_k - \delta) \approx A\,s(t - t_k) - (A\delta)\, s'(t - t_k)$$

which is still linear, now in $A_k$ and $(A\delta)_k$: solve
$y \approx [\,S \;\; -S'\,]\,[A;\; A\delta]$, take $\delta_k = (A\delta)_k / A_k$
(at most ±0.5 per step, ±0.75 in total), update $t_k$, repeat 4 times, then solve for A at
the final times. Below, one step by hand, then the function.
""")
code(r"""
Dmat = ctx.derivative(t_[:, None] - pos[None, :])
coef, *_ = np.linalg.lstsq(np.hstack([S, -Dmat]), yc, rcond=None)
A1, Ad1 = coef[:2], coef[2:]
print("one Gauss-Newton step: A =", np.round(A1 * 1e3, 2), "mV, delta =", np.round(Ad1 / A1, 3), "samples")
A_ref, t_ref = pileup._ls(yc, pos, ctx, refine=True, whiten=False)
print("ls_refine: A =", np.round(A_ref * 1e3, 2), "mV at times", np.round(t_ref, 3), "(detected:", pos, ")")
""")

md(r"""
## E6. `pileup.gls_refine`: whitened least squares (`_whitener`)

With noise covariance C (Toeplitz from the autocovariance), the maximum-likelihood estimate
minimises $(y - SA)^T C^{-1} (y - SA)$. With the Cholesky factor $C = LL^T$, multiply
everything by $L^{-1}$ (`solve_triangular`): the noise becomes white and ordinary least
squares applies to $L^{-1}y$ and $L^{-1}S$ ("prewhitening"). The timing refinement is the
same as E5.
""")
code(r"""
from scipy.linalg import solve_triangular
Lc = pileup._whitener(ctx, yc.size)
print("Cholesky factor:", Lc.shape, " L[0:3,0:3] =\n", np.round(Lc[:3, :3] * 1e3, 4), "(mV)")
# whitening real noise: take pulse-free stretches of capture 1 as long as the cluster
wins_, _ = noise.pulse_free_windows(y1, freeD, yc.size)
z = np.array([solve_triangular(Lc, w_, lower=True) for w_ in wins_])
print(f"{len(wins_)} real noise windows: before whitening std {wins_.std()*1e3:.3f} mV, lag-1 correlation {np.corrcoef(wins_[:, :-1].ravel(), wins_[:, 1:].ravel())[0,1]:.3f}")
print(f"after L^-1: std {z.std():.3f} (should be 1), lag-1 correlation {np.corrcoef(z[:, :-1].ravel(), z[:, 1:].ravel())[0,1]:.3f} (should be 0)")
print("gls_refine: A =", np.round(pileup.gls_refine(yc, pos, ctx) * 1e3, 2), "mV")
""")

md(r"""
## E7. `pileup.wls_refine`: noise + template-mismatch weighting (the chosen method)

GLS assumes the template is exact. Real pulses differ from it by a small, amplitude-
proportional amount. The fit residuals of the isolated pulses show it:

$$\text{rms}^2 = \sigma_0^2 + (\varepsilon A)^2$$

`wls_refine` adds that model error to the covariance:
$C = C_\text{noise} + (\varepsilon\, A_\max)^2 I$, with $A_\max$ from a first unweighted fit.
Small pulses (noise-dominated) are weighted like GLS; large pulses (mismatch-dominated) like
ordinary LS. Below, ε is measured from the clean events of the event file.
""")
code(r"""
cl_ev = efull.clean(); a_ = cl_ev["amplitude"]; r_ = cl_ev["fit_rms"]
m_ = (a_ > 0.005) & (a_ < 0.6)
c0, c1 = np.linalg.lstsq(np.vstack([np.ones(m_.sum()), a_[m_] ** 2]).T, r_[m_] ** 2, rcond=None)[0]
print(f"fit of rms^2 = sigma0^2 + (eps A)^2 over {m_.sum()} clean events: sigma0 = {np.sqrt(c0)*1e3:.3f} mV, eps = {np.sqrt(c1)*100:.3f} %")
A0 = pileup.ls(yc, pos, ctx)
print(f"extra variance for this cluster: (eps * A_max)^2 = ({ctx.model_error} * {A0.max()*1e3:.1f} mV)^2 = {(ctx.model_error*A0.max()*1e3)**2:.3f} mV^2  vs noise r[0] = {r_long[0]*1e6:.3f} mV^2")
print("wls_refine: A =", np.round(pileup.wls_refine(yc, pos, ctx) * 1e3, 2), "mV")
""")

md(r"""
## E8. `pileup.resolve(method, y, pos, ctx)`: all methods on this cluster
`resolve` dispatches by name. Below: every method's amplitudes, and the WLS fit with its
residual.
""")
code(r"""
for m in pileup.METHODS:
    print(f"{m:11s}", np.round(pileup.resolve(m, yc, pos, ctx) * 1e3, 2), "mV")
A_w, t_w = pileup._ls(yc, pos, ctx, refine=True, whiten=True,
                      extra_var=(ctx.model_error * pileup.ls(yc, pos, ctx).max()) ** 2)
model_w = Tp(t_[:, None] - t_w[None, :]) @ A_w
fig, ax = P.figure(1, 1, height_per_row=3.0)
ax[0,0].plot(t_, yc * 1e3, marker="o", markersize=2.5, linewidth=0.8, label="data")
ax[0,0].plot(t_, model_w * 1e3, linewidth=1.2, label="WLS fit (sum of two templates)")
ax[0,0].plot(t_, (yc - model_w) * 1e3, linewidth=1.0, label="residual")
ax[0,0].set_ylabel("mV"); ax[0,0].legend(); plt.show()
""")

md(r"""
## E9. `solver.PileupSolver` (step 6): the chosen method inside the pipeline

`PileupSolver(template, noise_acov, SolverConfig)` is built once per capture.

* `taper_autocovariance(acov, 512)`: keeps 512 lags and tapers the last quarter smoothly
  to 0 (half-cosine). An abruptly truncated autocovariance can make the Toeplitz matrix
  indefinite, and then the Cholesky factorisation fails.
* `_whitener(n, extra_var)`: C = Toeplitz(r) + extra_var·I → Cholesky L → **W = L⁻¹**,
  formed once per cluster, so every later whitening is a single matrix product.
* `_fit(y, t0, W, Wy)`: the WLS + timing fit of E5–E7, solving the small normal equations
  $(X_w^T X_w)\,c = X_w^T y_w$ each iteration and stopping when no time moves by more
  than 0.005 samples (times stay within ±1.25 samples of the detected positions).
  `PulseTemplate.derivative(t)` gives s'(t) from a precomputed central-difference array,
  and `PulseTemplate.__call__` interpolates by direct indexing on its uniform grid.
* `solve(y, positions)`: (1) unweighted LS to get $A_\max$ for the mismatch term;
  (2) whitening matrix; (3) fit; (4) **residual search**: for candidate times τ every half
  sample, with whitened residual $r_w$ and whitened template $s_w(\tau)$,

  $$A^*(\tau) = \frac{s_w\cdot r_w}{s_w\cdot s_w},\qquad z(\tau) = \frac{s_w\cdot r_w}{\sqrt{s_w\cdot s_w}}$$

  ($z^2$ is the drop in $\chi^2$ from adding that pulse). If max z > 8, positive, at least
  one sample from the existing pulses: add it, refit everything, and keep it only if its
  amplitude is still > 8σ. (5) rms and whitened χ² per sample.

Below: a real cluster where the detector saw **one** pulse but the solver found two.
""")
code(r"""
from gammaspec import solver
ev_s = efull.events
added_ev = ev_s[(ev_s["flags"] & events.FLAG_FIT_ADDED) != 0]
pick = added_ev[(added_ev["amplitude"] > 0.03)][0]
cap_s, cl_s = int(pick["capture"]), int(pick["cluster"])
rows_s = ev_s[ev_s["cluster"] == cl_s]
detected_rows = rows_s[(rows_s["flags"] & events.FLAG_FIT_ADDED) == 0]
pipe_s = pipeline.Pipeline(T, FS)
with store.CaptureStore(SRC_STORE) as s:
    pipe_s.pre.reset(); ys = pipe_s.pre.process(s.volts(cap_s))
a_s = int(detected_rows["index"][0]) - 10
b_s = int(detected_rows["index"][-1]) + 41
free_s = signal.pulse_free_mask(ys.size, pipe_s.detector.detect(ys).index, 20, 120)
sv = solver.PileupSolver(Tp, detect.noise_autocovariance(ys, free_s, 512), solver.SolverConfig())
res_s = sv.solve(ys[a_s:b_s], detected_rows["index"] - a_s)
print("detected positions:", detected_rows["index"] - a_s, "-> solver times:", np.round(res_s.time, 2),
      "amplitudes (mV):", np.round(res_s.amplitude * 1e3, 1), "added:", res_s.added)
print(f"rms {res_s.rms*1e3:.3f} mV, whitened chi2 per sample {res_s.chi2:.2f}")
no_search = solver.PileupSolver(Tp, detect.noise_autocovariance(ys, free_s, 512), solver.SolverConfig(residual_search=False))
r1 = no_search.solve(ys[a_s:b_s], detected_rows["index"] - a_s)
tt_s = np.arange(b_s - a_s, dtype=float)
m1 = Tp(tt_s[:, None] - r1.time[None, :]) @ r1.amplitude
m2 = Tp(tt_s[:, None] - res_s.time[None, :]) @ res_s.amplitude
print(f"without the search: one pulse of {r1.amplitude[0]*1e3:.1f} mV, rms {r1.rms*1e3:.3f} mV, chi2 {r1.chi2:.1f}")
fig, ax = P.figure(1, 1, height_per_row=3.0)
ax[0,0].plot(tt_s, ys[a_s:b_s] * 1e3, marker="o", markersize=2.5, linewidth=0.8, label="data")
ax[0,0].plot(tt_s, (ys[a_s:b_s] - m1) * 1e3, linewidth=1.0, label="residual, one pulse")
ax[0,0].plot(tt_s, (ys[a_s:b_s] - m2) * 1e3, linewidth=1.0, label="residual, with the found pulse")
ax[0,0].set_ylabel("mV"); ax[0,0].legend(); plt.show()
""")

md(r"""
## E10. `Pipeline._solve_clusters`, `EventFile.good()` and `pipeline.spectrum_report`

* `_solve_clusters(e, y, detected)`: noise autocovariance of the capture (pulse-free mask
  around the detections); for each cluster of > 1 pulse, and each isolated pulse flagged
  BAD_FIT, `solve` the window [first − 10, last + 40]; write amplitudes, shifts
  (= solver time − index) and rms; set or clear BAD_FIT (rms > 4σ + 2 % A_max); append
  found pulses as new rows (PILEUP | FIT_ADDED); sort; recompute gaps and cluster sizes.
* `EventFile.good(include_pileup=True)`: every event with a finite amplitude whose flags are
  only PILEUP / FIT_ADDED (how it was measured), i.e. no warm-up, edge, bad fit or large
  pulse. `clean()` stays "no flags at all" (isolated pulses).
* `spectrum_report(events)`: Gaussian peak fits (`fit_peaks`, D12) of isolated, resolved
  pile-up (by gap) and all good events, with shifts relative to isolated and the
  throughput (good / valid events). Categories with < 1000 events in the peak region are
  not fitted.
""")
code(r"""
rep = pipeline.spectrum_report(ROOT / "data/events/source_first100.h5")
print(f"throughput {100*rep['throughput']:.1f} % (isolated only {100*rep['isolated_fraction']:.1f} %), pulses found by the solver: {rep['fit_added']}")
for name, c in rep["categories"].items():
    if "error" in c:
        print(f"  {name:46s} n={c['n']:6d}  not fitted ({c['error']})")
    else:
        print(f"  {name:46s} n={c['n']:6d}  FWHM {' / '.join(f'{x:.2f}' for x in c['fwhm_pct'])} %   shift {' / '.join(f'{x:+.2f}' for x in c['shift_pct'])} %")
""")

# =============================================================================================
md(r"""
---
# Part G: Energy calibration and background subtraction (step 7)

Input: the good events of the source run and of the background run (event files from
`gammaspec process`). Output: both on a keV axis, live-time-normalised spectra, and the
net spectrum.

## G1. `calibration.anchor_gain(amplitudes, kind)`
A first gain (keV per volt) from an unmistakable line. `_peaks_in` histograms the
amplitudes between their median and 99.5th percentile (400 bins), smooths them, and returns
the most prominent peaks. Starting at the median keeps a crowd of small (noise) triggers
from swamping the search (this matters for the naive chain of Part H).
* **source:** the pair of peaks whose centre ratio is closest to 1332.5 / 1173.2 = 1.1357
  (within 2 %) is ⁶⁰Co: gain = mean(1173.2 / A₁, 1332.5 / A₂);
* **background:** the strongest peak is ⁴⁰K: gain = 1460.8 / A.
""")
code(r"""
from gammaspec import calibration, analysis
src_good = efull.good()
bkg_good = events.EventFile(ROOT / "data/events/background.h5").good()
g_src = calibration.anchor_gain(src_good["amplitude"], "source")
g_bkg = calibration.anchor_gain(bkg_good["amplitude"], "background")
print(f"source anchor gain {g_src:.1f} keV/V  ->  1173.2 keV predicted at {1173.228/g_src:.4f} V")
print(f"background anchor gain {g_bkg:.1f} keV/V ->  1460.8 keV predicted at {1460.82/g_bkg:.4f} V")
""")

md(r"""
## G2. `calibration.fit_lines(amplitudes, lines, gain)` → `_fit_group`

Each library line is predicted at A = E / gain. Lines within 16 % of each other form a
group and are fitted together by `spectrum.fit_peaks` with **bounds**: centre ±4 %, width
0.5–2× the expected NaI resolution
$\sigma(E) = E \cdot \tfrac{0.065}{2.3548}\sqrt{1173/E}$ (`expected_sigma_rel`), area ≥ 0.
A line is accepted if its area is ≥ 4σ, its centre and width are not at a bound, and its
FWHM is 2.5–15 %. If any are rejected, one is dropped and the rest of the group refitted,
until all remaining lines are accepted. Lines whose centre or width is at a bound are
dropped first (they are the degenerate ones), then the least significant.
""")
code(r"""
pts = calibration.fit_lines(src_good["amplitude"], calibration.SOURCE_LINES, g_src)
for p_ in pts:
    print(f"{p_.label:10s} {p_.energy_kev:8.1f} keV  A = {p_.amplitude_v:.4f} V  FWHM {p_.fwhm_pct:5.2f} %  area {p_.area:7.0f}  "
          + ("accepted" if p_.accepted else "rejected: " + p_.note))
""")

md(r"""
## G3. `calibration.fit_calibration(points)`

Weighted least squares of $E = \sum_i c_i A^i$ through the accepted points. The weight of a
point is $1/\sigma_E$ with $\sigma_E^2 = (\sigma_A \cdot \text{gain})^2 + (0.001\,E)^2$ (centre
uncertainty plus a 0.1 % systematic). Linear by default; quadratic only with ≥ 4 points and
if it halves χ²/ndf. With exactly two points the line passes through them, and a weak pin
at (0, 0) keeps the offset reasonable.
""")
code(r"""
cal_src = calibration.fit_calibration(pts)
print("coefficients (keV, keV/V):", np.round(cal_src.coeffs, 3), "  chi2/ndf", round(cal_src.chi2_ndf, 2))
for lab, e_, r_ in cal_src.residuals_kev():
    print(f"  residual {lab:10s} {e_:8.1f} keV: {r_:+.2f} keV")
""")

md(r"""
## G4. `calibration.gain_match(points, reference)`

The background run has only ⁴⁰K. Same detector and electronics, so the calibration
*shape* (offset, curvature) is taken from the source run, with one gain factor k:

$$E_{bkg}(A) = E_{ref}(k A) = \sum_i c_i k^i A^i$$

k is found for each accepted line by solving $E_{ref}(k A_\text{line}) = E_\text{line}$
(Brent's method), and averaged with weights $(A/\sigma_A)^2$.
""")
code(r"""
pts_b = calibration.fit_lines(bkg_good["amplitude"], calibration.BACKGROUND_LINES, g_bkg)
print("accepted background lines:", [(p_.label, p_.energy_kev) for p_ in pts_b if p_.accepted])
cal_bkg = calibration.gain_match(pts_b, cal_src)
k_ = cal_bkg.coeffs[1] / cal_src.coeffs[1]
print(f"k = {k_:.4f}: the background run's gain is {100*(1-1/k_):.1f} % lower;  E_bkg(A) = {cal_bkg.coeffs[0]:.2f} + {cal_bkg.coeffs[1]:.2f} A")
print(f"without gain matching, the background 40K would appear at {cal_src.energy(pts_b[[p_.label for p_ in pts_b].index('K-40')].amplitude_v):.0f} keV")
""")

md(r"""
## G5. `calibration.calibrate(...)` and `calibration.gain_vs_time(events, kind)`
`calibrate` = anchor (G1) → lines (G2) → curve (G3), or gain matching (G4) when only one
line is usable; repeated once with the gain from the accepted lines. `gain_vs_time` refits
the anchor line(s) in chunks of 20 captures; relative gain per chunk (± its uncertainty).
""")
code(r"""
print("source:", np.round(calibration.calibrate(src_good["amplitude"], "source").coeffs, 3))
for g_ in calibration.gain_vs_time(src_good, "source"):
    print(f"  captures {g_['first_capture']:3d}-{g_['first_capture']+19:3d}: gain {g_['relative_pct']:+.2f} % ± {g_['relative_err_pct']:.2f} %")
""")

md(r"""
## G6. `analysis.energy_spectrum(...)` and `analysis.subtract(S, B)`

Each run's events → keV with its own calibration → histogram on common 5 keV bins →
divide by that run's live time. Per bin, with counts N and live times t:

$$r_\text{net} = \frac{N_s}{t_s} - \frac{N_b}{t_b},\qquad
\sigma^2 = \frac{N_s}{t_s^2} + \frac{N_b}{t_b^2}$$

(Poisson; empty bins count as 1 for the error).
""")
code(r"""
edges_ = np.arange(0, 3000.1, 5.0)
S_ = analysis.energy_spectrum(cal_src.energy(src_good["amplitude"]), efull.live_time_s, edges_)
B_ = analysis.energy_spectrum(cal_bkg.energy(bkg_good["amplitude"]), events.EventFile(ROOT / "data/events/background.h5").live_time_s, edges_)
N_ = analysis.subtract(S_, B_)
i0 = np.searchsorted(edges_, 1450)
print("bin (keV)    source counts  bkg counts   source cps  bkg cps   net cps ± err")
for i in range(i0, i0 + 4):
    print(f"{edges_[i]:5.0f}-{edges_[i+1]:5.0f}   {S_.counts[i]:8.0f}     {B_.counts[i]:6.0f}     {S_.rate[i]:7.2f}   {B_.rate[i]:6.2f}   {N_.rate[i]:+6.2f} ± {N_.rate_err[i]:.2f}")
""")

md(r"""
## G7. `analysis.fit_line_groups(spectrum)` → `fit_net_peaks`

Line rates from Gaussian fits on a linear continuum, in **local windows**
(`NET_GROUPS`): the ⁶⁰Co pair over 1040–1440 keV; ⁴⁰K over 1300–1650 keV with the 1332 keV
peak as a **nuisance** (area ≥ 0: it can be absent, not negative), because ⁴⁰K sits on its
flank. Reported lines may be negative (a line that should vanish in the net spectrum).
""")
code(r"""
for name_, sp_ in (("source", analysis.NetSpectrum(edges_, S_.rate, S_.rate_err)),
                   ("background", analysis.NetSpectrum(edges_, B_.rate, B_.rate_err)), ("net", N_)):
    print(name_.ljust(10), "  ".join(f"{r_['label']} {r_['energy_kev']:.0f}: {r_['rate_cps']:7.2f} ± {r_['rate_err_cps']:.2f}" for r_ in analysis.fit_line_groups(sp_)))
""")

md(r"""
## G8. `analysis.analyze(source, out_dir, background)` (the step-7 driver)
Calibrate both runs, build and subtract the spectra, fit the lines, measure the gain drift,
and write `calibration_*.json`, `spectrum.csv`, `spectrum_summary.json` and the plots. CLI:
`gammaspec spectrum SOURCE_EVENTS OUT_DIR --background BACKGROUND_EVENTS`. Results:
`docs/spectrum.md`.
""")
code(r"""
summ = analysis.analyze(ROOT / "data/events/source_first100.h5", TMP / "spec", ROOT / "data/events/background.h5", verbose=False)
print("totals above 120 keV (counts/s):", {k: round(v) for k, v in summ["totals_cps"].items()})
fig, ax = P.figure(1, 1, height_per_row=3.4)
w_ = np.diff(edges_)
ax[0,0].stairs(S_.rate / w_, edges_, label="source"); ax[0,0].stairs(B_.rate / w_, edges_, label="background")
ax[0,0].stairs(np.clip(N_.rate / w_, 1e-4, None), edges_, label="net (clipped at 1e-4)")
ax[0,0].set_yscale("log"); ax[0,0].set_xlabel("energy (keV)"); ax[0,0].set_ylabel("counts / s / keV"); ax[0,0].legend(); plt.show()
""")

# =============================================================================================
md(r"""
---
# Part H: Benchmark (naive count vs analog chain vs gammaspec)

The same captures through three chains; each is calibrated by `calibration.calibrate`
(Part G) and measured with the same metrics. Results on 100 captures: `docs/benchmark.md`.
Here, 5 captures to keep it quick.

## H1. `benchmark.naive_chain(store, captures, k=5, distance=5)`
No transformation: flip polarity, subtract the median, then `scipy.signal.find_peaks` with
height > k·σ (`signal.noise_sigma`) and ≥ `distance` samples between peaks. The pulse
height is the highest sample. Pulses on an undershoot are measured low, noise maxima count
as pulses, and pile-up merges or splits pulses.
""")
code(r"""
from gammaspec import benchmark
caps5 = range(5)
naive = benchmark.naive_chain(SRC_STORE, caps5)
print(f"{naive.n_recorded} 'pulses' in {naive.live_time_s:.2f} s = {naive.n_recorded/naive.live_time_s:.0f} /s;  smallest heights (mV):", np.round(np.sort(naive.amplitudes)[:5]*1e3, 2))
""")

md(r"""
## H2. `benchmark.shaping_amplifier(x, dt, tau_shaping, tau_pz, n_integrators=4)`
The digital version of the analog shaping amplifier, CR-RC⁴ with pole-zero cancellation.
With $a = e^{-\Delta t/\tau}$ for each time constant:

$$H(z) = \underbrace{\frac{1 - a_{pz} z^{-1}}{1 - a_s z^{-1}}}_{\text{PZ-cancelled CR}}\;
\left(\frac{1 - a_s}{1 - a_s z^{-1}}\right)^{4}$$

The CR zero at $a_{pz}$ cancels the pulse's own 4.1 µs decay pole, so the input becomes a
clean $\tau_s$ exponential. Four RC integrators then make it a semi-Gaussian
$\propto t^4 e^{-t/\tau_s}$. In continuous time this peaks at $4\tau_s$; each zero-delay
digital stage runs ~½ sample ahead, so the sampled peak comes ≈ 2.5 samples earlier.
""")
code(r"""
k_ = np.arange(120); xp = np.where(k_ >= 10, np.exp(-(k_ - 10) / FS / 4.1e-6), 0.0)
fig, ax = P.figure(1, 1, height_per_row=2.8)
ax[0,0].plot(k_, xp, color=P.SERIES[7], label="input: 4.1 µs exponential")
for i_, ts in enumerate([0.5e-6, 1e-6, 2e-6]):
    y_ = benchmark.shaping_amplifier(xp, 1 / FS, ts)
    ax[0,0].plot(k_, y_ / y_.max(), color=P.SERIES[i_], label=f"CR-RC⁴ τ = {ts*1e6:g} µs (normalised)")
    print(f"tau {ts*1e6:g} us: peak {np.argmax(y_) - 10} samples after the step, gain {y_.max():.3f}, end value {y_[-1]:.1e}")
ax[0,0].set_xlabel("sample"); ax[0,0].legend(); plt.show()
""")

md(r"""
## H3. `benchmark.analog_chain(store, captures, tau_shaping, tau_pz, lld_sigma=5, dead_time=5e-6, pur=True)`
Per capture:
1. **shape** (H2);
2. **discriminate**: peaks of the shaped signal above 5σ, with σ from the negative side
   (`detect.negative_side_sigma`; differences between samples would underestimate the noise
   of a smooth signal);
3. **baseline restorer**: subtract `track_baseline` of the pulse-free samples;
4. **peak-hold ADC**: parabolic interpolation through the three samples around each maximum,
   $A = y_0 - \tfrac14 (y_- - y_+)\,\delta$ with $\delta = \tfrac{y_- - y_+}{2(y_- - 2y_0 + y_+)}$;
5. **pile-up rejector**: drop both pulses if two triggers are < 10 τ apart;
6. **dead time**: no conversion for 5 µs after each accepted one (non-extending).

The live-time clock stops during conversions and for rejected pile-up.
""")
code(r"""
ana = benchmark.analog_chain(SRC_STORE, caps5, tau_shaping=1e-6, pur=True)
print(f"triggers {ana.n_triggers}, recorded {ana.n_recorded} ({100*ana.n_recorded/ana.n_triggers:.1f} %), live {ana.live_time_s:.3f} of {ana.real_time_s:.3f} s ({100*ana.live_time_s/ana.real_time_s:.1f} %)")
""")

md(r"""
## H4. `benchmark.pipeline_chain(events, captures)` and `benchmark.metrics(chain)`
`pipeline_chain` takes the good events of the event file. `metrics` calibrates a chain and
returns, among others:
* the FWHM of 1173/1332 keV (joint bounded Gaussian fit, 1040–1440 keV, 2 keV bins);
* the valley-to-peak ratio (minimum between the peaks / mean peak maximum);
* the peak rates per live second (`analysis.fit_line_groups`);
* the pile-up continuum (1450–2450 keV);
* the low-energy threshold (1st percentile of recorded energies below 300 keV).

`run_benchmark` does this for all five chains and writes the summary and plot (CLI
`gammaspec benchmark`).
""")
code(r"""
pipe = benchmark.pipeline_chain(ROOT / "data/events/source_first100.h5", caps5)
for ch in (naive, ana, pipe):
    m_ = benchmark.metrics(ch)
    print(f"{m_['name']:42s} FWHM {m_['fwhm_pct'][0]:.2f}/{m_['fwhm_pct'][1]:.2f} %  valley/peak {m_['valley_to_peak']:.3f}  threshold {m_['low_energy_threshold_kev']:.0f} keV")
print("(5 captures only: widths and valleys fluctuate by up to ~1 % here; docs/benchmark.md uses 100)")
""")

# =============================================================================================
md(r"""
---
# Part I: Configuration and the runner (step 9)

Everything above was called from code. A **run file** does the same without code
(`docs/configuration.md`).

## I1. `config.from_dict(d)` / `config.load(path, overrides)`
A nested dict (from TOML or JSON) → the config dataclasses, recursively. Missing keys keep
their defaults. Each value is checked against its field's type annotation; unknown keys,
wrong types, invalid choices (`CHOICES`) and out-of-range numbers (`POSITIVE`, ≥ 0) are
collected and raised together as a `ConfigError`.
""")
code(r"""
from gammaspec import config as C
cfg = C.from_dict({"input": {"kind": "store", "source": str(SRC_STORE)}, "pipeline": {"detect": {"threshold": 6}}})
print(type(cfg.pipeline.detect).__name__, "threshold =", cfg.pipeline.detect.threshold, "| method (default) =", cfg.pipeline.detect.method)
try:
    C.from_dict({"input": {"knd": "raw"}, "pipeline": {"detect": {"method": "fast", "threshold": -2}}})
except C.ConfigError as e:
    print(e)
""")

md(r"""
## I2. `config.set_key(d, "a.b.c=value")` and `config.dumps(cfg)`
`set_key` applies one `--set` override; the value is parsed as TOML (`6` → int, `true` →
bool, `"x"` → str). `dumps` writes TOML, one table per dataclass. The comment of each field
comes from its comment in the source code (`_comments` reads the class source).
""")
code(r"""
d = {}
for s_ in ["pipeline.preprocess.line_removal=false", "pipeline.detect.method=matched", "input.limit=2"]:
    C.set_key(d, s_)
print(d)
print("\n".join(C.dumps(C.from_dict(d)).splitlines()[38:48]))
""")

md(r"""
## I3. `config.check(cfg)`
Checks that need more than one setting: required inputs and existing paths are errors,
combinations of stages that work badly are warnings.
""")
code(r"""
c_ = C.load(None, [f"input.source='{SRC_STORE.as_posix()}'", f"template.path='{TEMPLATE.as_posix()}'",
                   "pipeline.preprocess.baseline_tracker=false"])
for w_ in C.check(c_):
    print("warning:", w_)
""")

md(r"""
## I4. `runner.run(cfg)`
Input → template → events → spectrum → benchmark, each when selected. Every output gets a
stamp (`<file>.stamp.json`: its inputs and settings). A rerun reuses an output whose stamp
is unchanged, so after changing one setting only the affected steps run.
""")
code(r"""
from gammaspec import runner
sets_ = [f"output='{(TMP / 'run').as_posix()}'", f"input.source='{SRC_STORE.as_posix()}'", f"template.path='{TEMPLATE.as_posix()}'",
         "input.limit=2", "spectrum.enabled=false", "workers=1"]
r_ = runner.run(C.load(None, sets_), verbose=False)
print({k: v["seconds"] for k, v in r_["steps"].items()})
r_ = runner.run(C.load(None, sets_), verbose=True)                                   # unchanged: reused
r_ = runner.run(C.load(None, sets_ + ["pipeline.detect.threshold=7"]), verbose=False)  # only processing reruns
print("stored in the event file:", C.pipeline_from_event_file(r_["events"]["source"]).detect.threshold)
""")

# =============================================================================================
md(r"""
---
# Part F: Evaluation tools (not part of the pipeline)

These modules measure how well the pipeline works; they produced the numbers in the docs.

## F1. `restoration.default_configs()` and `restoration.evaluate(...)`
Every combination of pre-processing stages is applied to the **same** pulses (selected on
the raw signal). For each: noise, the baseline just before each pulse, the average tail
after isolated pulses, and Gaussian fits of the ⁶⁰Co peaks with free / fixed baseline, also
split into pulses on an earlier pulse's tail vs clean ones (`docs/baseline-restoration.md`).
""")
code(r"""
from gammaspec import restoration
for name, cfg in restoration.default_configs().items():
    print(f"{name:30s} -> label: {cfg.label()}")
""")

md(r"""
## F2. `detection_eval.noise_streams`, `inject`, `_matched`, `_both_found`

* `noise_streams`: pulse-free stretches of the clean signal joined together (**real
  noise**), and the same with **randomised Fourier phases** (**pure noise**): identical
  power spectrum, but any pulse is smeared into noise, so detections in it are false.
* `inject`: adds template pulses at random sub-sample times (optionally a second one d
  samples later).
* `_matched` / `_both_found`: is a detection within ±2 samples of the truth? For pairs,
  two **distinct** detections are required.
""")
code(r"""
from gammaspec import detection_eval
real_n, pure_n, _ = detection_eval.noise_streams(SRC_STORE, Tp, pz, [0, 50])
ns_r = noise.noise_spectrum(real_n[: real_n.size // 512 * 512].reshape(-1, 512), FS)
ns_p = noise.noise_spectrum(pure_n[: pure_n.size // 512 * 512].reshape(-1, 512), FS)
print(f"real noise {real_n.size} samples, rms {real_n.std()*1e3:.3f} mV; pure noise rms {pure_n.std()*1e3:.3f} mV")
print("PSD ratio pure/real at 50, 200, 500 kHz:", np.round(np.interp([5e4, 2e5, 5e5], ns_r.freq_hz, ns_p.psd / ns_r.psd), 2))
yi, ci, Ai = detection_eval.inject(real_n, Tp, [0.003], np.random.default_rng(0), second=(6, 0.03))
found_ = det.detect(yi).index
print("injected", ci.size, "pairs (3 mV then 30 mV, 6 samples apart): first found", detection_eval._matched(found_, np.round(ci).astype(int)).mean().round(3),
      "; both found as two events", detection_eval._both_found(found_, np.round(ci).astype(int), np.round(ci + 6).astype(int)).mean().round(3))
""")

md(r"""
## F3. `pileup_eval`: the three pile-up tests

* `make_context(template, noise)`: template → pole-zero design → filtered template, decay
  constant, noise autocovariance → `PileupContext`.
* `_robust(errors)`: median (bias) and 1.4826 × MAD (spread) of relative errors.
* `test_template_pairs`: test A. Pairs injected into real noise, **found by the detector**,
  resolved; each estimate divided by the method's estimate of the same noiseless pulse alone.
* `real_isolated_pulses` + `test_real_pairs`: test B. Real isolated pulses added in pairs;
  each estimate divided by the method's estimate of that real pulse alone.
* `test_real_events`: test C. The real pile-up snippets are resolved; the ⁶⁰Co peaks are
  compared with the method's own isolated-pulse peaks.
* `_speed`: time per two-pulse cluster.

A small version of test B:
""")
code(r"""
from gammaspec import pileup_eval
Wr = pileup_eval.real_isolated_pulses(SRC_STORE, Tp, pz, [3, 4])
print("real isolated pulses:", Wr.shape)
Bres, singles = pileup_eval.test_real_pairs(Wr, ctx, np.random.default_rng(0), n_per_d=60)
for d_ in ("4", "10"):
    print(f"separation {d_}: " + ", ".join(f"{m} {100*Bres[d_][m]['second']['bias']:+.1f}/{100*Bres[d_][m]['second']['spread']:.1f}%" for m in pileup.METHODS))
""")

nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata["kernelspec"] = {"name": "gammaspec", "display_name": "gammaspec (.venv)", "language": "python"}
out = HERE / "walkthrough.ipynb"
NotebookClient(nb, timeout=900, kernel_name="gammaspec", resources={"metadata": {"path": str(HERE)}}).execute()
nbf.write(nb, out)
print("wrote", out)
