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

1. Take candidates (local maxima above the threshold) in **descending** order.
2. For each: re-check it in the **residual** output (it may have shrunk when a neighbour
   was subtracted). Require residual ≥ threshold + 3 % × (largest accepted pulse whose
   response overlaps here) and a local maximum.
3. Parabola through the three samples → sub-sample position → phase row p of R.
   $A = y_r[j] / R_p[j-\text{base}]$; subtract $A\,R_p$ from the residual.
4. Return the accepted positions, amplitudes and the residual.

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
* **time** = capture time + (index + shift) × Δt;
* **snippets**: clean waveform of each pile-up cluster (−20…+80 samples) for the solver.
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

md(r"""
---
## Evaluation tools (not part of the pipeline)

* `restoration.evaluate`: runs every combination of the pre-processing stages on the same
  pulses and compares noise, pre-pulse baseline, tails and peak resolution
  (`docs/baseline-restoration.md`).
* `detection_eval.evaluate`: injects template pulses into real noise and into
  phase-scrambled noise (same spectrum, no pulses) to measure efficiency, false triggers,
  artifacts and pair resolution for every detection method (`docs/detection.md`).
""")

nb = nbf.v4.new_notebook()
nb.cells = cells
nb.metadata["kernelspec"] = {"name": "gammaspec", "display_name": "gammaspec (.venv)", "language": "python"}
out = HERE / "walkthrough.ipynb"
NotebookClient(nb, timeout=900, kernel_name="gammaspec", resources={"metadata": {"path": str(HERE)}}).execute()
nbf.write(nb, out)
print("wrote", out)
