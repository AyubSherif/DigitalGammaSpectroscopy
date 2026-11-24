<div align="center">

# ☢️ Gamma Spectrometry App

**Gamma-ray spectroscopy in software: from raw ADC samples to a calibrated, background-subtracted spectrum.**

![Python](https://img.shields.io/badge/python-%E2%89%A5%203.11-3776AB?logo=python&logoColor=white)
![Version](https://img.shields.io/badge/version-v0.3-blue)
![Status](https://img.shields.io/badge/status-experimental-orange)
![Detector](https://img.shields.io/badge/detector-NaI(Tl)-8A2BE2)

[Quick start](#-quick-start) •
[How it works](#%EF%B8%8F-what-it-does-step-by-step) •
[Results](#-results-at-a-glance) •
[Docs](#-documentation) •
[Roadmap](#%EF%B8%8F-roadmap)

</div>

---

> [!IMPORTANT]
> **📌 Origin of this code.** The core of this project (the **pile-up correction**, **baseline
> restoration**, **template matching** and **interference cleaning**) began as software add-ons
> I developed for a contract with **Serva Energy** in 2021. I then asked **Claude** (Anthropic's
> AI assistant) to package that code into the app in this repository.

> [!WARNING]
> **🚧 Not fully tested.** The app is still experimental. The processing chain has been run on
> recorded and simulated data, but the app as a whole has not been fully tested. Check results
> before relying on them.

---

## 💡 Overview

The Gamma Spectrometry App turns the digitised signal of a scintillation detector into a
calibrated, background-subtracted energy spectrum. It does in software, on an ordinary PC, the
work that is classically done by a **shaping amplifier** and a **multichannel analyser (MCA)**:
analog instruments that sit between the detector's preamplifier and the spectrum.

![analog chain vs Gamma Spectrometry App](docs/img/chain_analog_vs_digital.png)

It was developed and validated with **3 hours of 2021 data** from a **NaI(Tl)** detector
(PMT → preamplifier → ADC at 1.8 MS/s → PC) measuring a **⁶⁰Co** source, plus a room-background
run. A built-in signal simulator lets it be tested beyond that data.

### 📍 Status: v0.3

| | |
|---|---|
| ⚙️ **Processing** | The complete chain from raw ADC files to a net spectrum, run on the full dataset (6.4 M pulses; 35 GB of text become a 470 MB event list) |
| 🧪 **Validation** | Benchmarked against a naive peak count and an emulated analog chain, and checked against simulations with known truth |
| 🖥️ **Ways to run it** | A browser app, run files, live streams |

See the [Roadmap](#%EF%B8%8F-roadmap) for what's done and what's next.

---

## 🔄 What it replaces, and why

### 📟 The analog way

In a classic NaI(Tl) spectrometer, the preamplifier's output goes to two more instruments:

1. **Shaping amplifier.** It filters each pulse into a smooth, slower bump, which reduces
   noise. A pole-zero adjustment removes the undershoot after each pulse, a baseline
   restorer keeps the zero level steady, and a pile-up rejector discards pulses that overlap.
2. **Multichannel analyser (MCA).** For each bump, a peak-hold circuit and an ADC measure its
   height, once. The MCA adds one count to the matching channel and keeps a live-time clock
   for the time it was busy.

The output is a histogram. The pulses themselves are gone: to try other settings, you have to
measure again.

### 💻 The digital way (this software)

The preamplifier output is digitised directly by a fast ADC, and the app does the rest:

| Analog function | In this app |
|---|---|
| Shaping filter (noise reduction) | **Optimal filter**, designed for the measured pulse shape and noise ([detection](docs/detection.md)) |
| Pole-zero adjustment | **Pole-zero correction** computed from the measured pulse model ([baseline](docs/baseline-restoration.md)) |
| Baseline restorer | **Baseline tracker**, plus removal of interference lines |
| Pile-up rejector (discards both pulses) | **Pile-up solver**: overlapping pulses are separated and *each one is measured* ([pile-up](docs/pileup-solver.md)) |
| Peak-hold ADC (one sample of the bump) | **Template fit** over the whole pulse, with sub-sample timing |
| MCA histogram + live-time clock | **Event list**: every pulse with its time, amplitude and quality flags, then histograms in any binning; exact live time |
| Manual energy calibration | **Automatic calibration** from known lines; gain matching between runs ([spectrum](docs/spectrum.md)) |

### ✨ Why do it in software

- 💰 **Lower cost.** The shaping amplifier, pile-up rejector, MCA (and often the instrument
  crate that powers them) are separate, specialised instruments. Here they are replaced by one
  ADC board and a PC, and improvements come as software updates, not new hardware.
- 🎒 **More portable.** Detector, preamplifier, ADC board and a laptop are the whole
  spectrometer, which makes field and in-plant measurements practical.
- ⚡ **Better at high count rates.** Overlapping pulses are measured instead of thrown away, and
  the threshold is lower. On the same recorded pulses, the app beats an emulated analog chain
  on resolution, peak separation and threshold, and loses no pulses to pile-up rejection
  (see below).
- 💾 **Nothing is thrown away.** Every pulse is kept with its time stamp. A measurement can be
  re-processed later with other settings or better algorithms, without measuring again.
- 🔍 **Transparent.** Every step is documented and checked against simulations whose true
  content is known.

---

## 📊 Results at a glance

![spectrum before and after](docs/img/tour_spectrum_before_after.png)

*20 recorded captures: a naive count of the raw signal (before) and this app (after), each
calibrated to keV.*

The same 100 recorded captures of the ⁶⁰Co run, processed three ways and calibrated the same
way ([details](docs/benchmark.md)):

| | Naive peak count (`find_peaks`) | Analog chain (emulated: CR-RC⁴ 0.5 µs + pile-up rejection) | 🏆 **Gamma Spectrometry App** |
|---|---|---|---|
| Resolution (FWHM) at 1173 / 1332 keV | 8.27 / 7.44 % | 6.95 / 6.21 % | **6.54 / 6.05 %** |
| Valley between the ⁶⁰Co peaks (lower = better separated) | 0.39 | 0.24 | **0.215** |
| Lowest energy recorded | none (noise is counted) | ≈ 44 keV | **≈ 6 keV** |
| Overlapping pulses | merged or split | rejected | **separated and measured** |

🧪 **On simulated data with known content** ([simulation](docs/simulation.md)):

- ✅ the program finds **≥ 99.3 %** of pulses above 100 keV;
- ✅ it separates **99.6 %** of pulse pairs that are 1–33 µs apart;
- ✅ it adds **no measurable broadening** to the peaks.

> [!NOTE]
> Typical 3″×3″ NaI(Tl) detectors with analog electronics reach ≈ 7 % at 662 keV and ≈ 5 % at
> 1332 keV. At those energies, resolution is set mostly by the crystal and PMT, so those values
> come from other detectors and are context, not a like-for-like comparison.

---

## 📜 Background

I developed the core algorithms (pile-up correction, baseline restoration, template matching
and interference cleaning) in 2021 as software add-ons under a contract with **Serva Energy**,
for its digital data-acquisition system. At the time the system sampled at 1.8 MS/s; it moved
to 3.4 MS/s as that work ended. A year later, the company presented a modular system recording
up to 250 MS/s at 16 bits, with its own real-time analysis, at the 11th ICI conference
(Saskatoon, July 2023). There it reported a slightly better 662 keV resolution and lower dead
time than a commercial digital MCA (Canberra DSA-1000).

🤖 **This repository** takes those add-ons and packages them, with Claude's help, into a
standalone app with a pipeline, a browser interface, a CLI, documentation and validation. It
is not an official Serva Energy product and is not endorsed by Serva Energy.

---

## ⚙️ What it does, step by step

```text
raw ADC capture files
  → 1. lossless binary store        (7 MB of text per 0.29 s capture → 0.7 MB)
  → 2. characterisation             (pulse template, pulse model, noise)
  → 3. interference removal, pole-zero correction, baseline tracking
  → 4. pulse detection              (optimal filter, successive cancellation)
  → 5–6. pile-up resolution         (weighted least squares, sub-sample timing)
  → 7. calibration, spectra, background subtraction
```

🎛️ Every stage can be switched on or off and configured, from the app or a run file.

📥 Input can be raw files, prepared stores, a live stream (a file or FIFO a digitiser writes
to, or a recording replayed at real time), or the built-in simulator.

> [!TIP]
> 👀 **To see each step with before/after pictures, start with
> [Part 0 of the walkthrough](docs/walkthrough.ipynb).**

---

## 📦 Installation

Requires **Python ≥ 3.11**.

```bash
git clone https://github.com/AyubSherif/GammaSpectroscopySoftwareAddOns.git
cd GammaSpectroscopySoftwareAddOns
python -m venv .venv
# Windows: .venv\Scripts\activate    Linux/macOS: source .venv/bin/activate
pip install -e ".[app]"         # the program and the browser app
pip install -e ".[app,dev]"     # plus tests and the walkthrough notebook
```

---

## 🚀 Quick start

### 🌐 The app

Choose the input, switch stages on and off, run, and watch the live spectrum:

```bash
gamma-spectrometry app                   # opens http://localhost:8501
```

### 🧪 Without data

A simulated ⁶⁰Co measurement, processed and checked against the truth:

```bash
gamma-spectrometry run examples/run_simulation.toml
gamma-spectrometry stream sim out/sim --duration 20 --show        # live spectrum window
```

### 📁 With recorded data

From a run file:

```bash
gamma-spectrometry config init my_run.toml      # every setting, at its default, commented
# edit: [input] kind = "raw", source = "ADC_raw_output", background = "ADC_raw_output_background"
gamma-spectrometry run my_run.toml              # convert → template → events → spectrum (→ benchmark)
```

### 📡 As a live stream

```bash
gamma-spectrometry stream data/source_first100.h5 out/replay --template T.h5 --realtime --show   # replay at real time
my_digitiser | gamma-spectrometry stream - out/live --template T.h5                              # samples on stdin
```

### 🧩 Step by step

Individual commands ([CLI reference](docs/cli.md)): `index`, `convert`, `characterize`,
`process`, `spectrum`, `benchmark`, `simulate`, …

Or from Python:

```python
from gamma_spectrometry.store import CaptureStore

with CaptureStore("data/source_first100.h5") as store:
    volts = store.volts(0)          # the first capture, in volts
```

---

## 📚 Documentation

In reading order, each document with the code it describes:

| # | Document | Contents | Code |
|---|---|---|---|
| 1 | 📓 [docs/walkthrough.ipynb](docs/walkthrough.ipynb) | **Part 0:** the app in pictures, before and after every step. **Parts A–K:** every function in call order, on real data, with the maths | all; built by [docs/build_walkthrough.py](docs/build_walkthrough.py) |
| 2 | 🏗️ [docs/architecture.md](docs/architecture.md) | Pipeline design, module map, how the stages fit together | all (module map) |
| 3 | 🗂️ [docs/data-formats.md](docs/data-formats.md) | Step 1: raw capture format, HDF5 store, event file, truth and live-spectrum files | [rawio.py](gamma_spectrometry/rawio.py), [store.py](gamma_spectrometry/store.py), [events.py](gamma_spectrometry/events.py) |
| 4 | 📈 [docs/characterization.md](docs/characterization.md) | Step 2: measured pulse shape, pulse template, noise and preliminary spectra | [signal.py](gamma_spectrometry/signal.py), [template.py](gamma_spectrometry/template.py), [noise.py](gamma_spectrometry/noise.py), [characterize.py](gamma_spectrometry/characterize.py), [spectrum.py](gamma_spectrometry/spectrum.py), [plotting.py](gamma_spectrometry/plotting.py) |
| 5 | 🧹 [docs/baseline-restoration.md](docs/baseline-restoration.md) | Step 3: interference removal, pole-zero correction, baseline tracker | [interference.py](gamma_spectrometry/interference.py), [baseline.py](gamma_spectrometry/baseline.py), [preprocess.py](gamma_spectrometry/preprocess.py), [restoration.py](gamma_spectrometry/restoration.py) |
| 6 | 🎯 [docs/detection.md](docs/detection.md) | Step 4: pulse detection and the event list | [detect.py](gamma_spectrometry/detect.py), [events.py](gamma_spectrometry/events.py), [pipeline.py](gamma_spectrometry/pipeline.py), [detection_eval.py](gamma_spectrometry/detection_eval.py) |
| 7 | ⚖️ [docs/pileup-benchmark.md](docs/pileup-benchmark.md) | Step 5: pile-up methods (TM, HYPER, LS/GLS/WLS) compared | [pileup.py](gamma_spectrometry/pileup.py), [pileup_eval.py](gamma_spectrometry/pileup_eval.py) |
| 8 | 🔀 [docs/pileup-solver.md](docs/pileup-solver.md) | Step 6: the pile-up solver in the pipeline | [solver.py](gamma_spectrometry/solver.py), [pipeline.py](gamma_spectrometry/pipeline.py) |
| 9 | 🌈 [docs/spectrum.md](docs/spectrum.md) | Step 7: calibration, gain matching, background subtraction; full 3-hour results | [calibration.py](gamma_spectrometry/calibration.py), [analysis.py](gamma_spectrometry/analysis.py), [spectrum.py](gamma_spectrometry/spectrum.py) |
| 10 | 🏁 [docs/benchmark.md](docs/benchmark.md) | Step 8: this app vs a naive count vs an emulated analog chain vs typical analog systems | [benchmark.py](gamma_spectrometry/benchmark.py) |
| 11 | 🛠️ [docs/configuration.md](docs/configuration.md) | Step 9: the run file: input, stage switches, every parameter | [config.py](gamma_spectrometry/config.py), [runner.py](gamma_spectrometry/runner.py), [examples/](examples/) |
| 12 | 🎲 [docs/simulation.md](docs/simulation.md) | Step 10: the signal simulator, and the pipeline scored against the truth | [simulate.py](gamma_spectrometry/simulate.py) |
| 13 | 📡 [docs/streaming.md](docs/streaming.md) | Step 10: live input, overlapping frames, live spectrum, speed | [sources.py](gamma_spectrometry/sources.py), [stream.py](gamma_spectrometry/stream.py), [liveplot.py](gamma_spectrometry/liveplot.py) |
| 14 | 🌐 [docs/app.md](docs/app.md) | Step 11: the browser app | [app.py](gamma_spectrometry/app.py) |
| 15 | ⌨️ [docs/cli.md](docs/cli.md) | Command-line reference | [cli.py](gamma_spectrometry/cli.py) |
| 16 | 🧑‍💻 [docs/development.md](docs/development.md) | Environment, tests, conventions | [tests/](tests/), [pyproject.toml](pyproject.toml) |
| 17 | 📝 [docs/adr/](docs/adr/) | Architecture Decision Records: why things are the way they are (0001–0013) | — |
| 18 | 🕒 [CHANGELOG.md](CHANGELOG.md) | Release history | — |

---

## 🗺️ Roadmap

- [x] **1.** Ingestion: fast reader, lossless binary store, chronological index
- [x] **2.** Pulse template and noise characterisation ([results](docs/characterization.md))
- [x] **3.** Interference removal and baseline restoration ([results](docs/baseline-restoration.md))
- [x] **4.** Pulse detection and the event list ([results](docs/detection.md))
- [x] **5.** Pile-up methods compared: template matching vs HYPER vs least squares ([results](docs/pileup-benchmark.md))
- [x] **6.** Pile-up solver in the pipeline ([results](docs/pileup-solver.md))
- [x] **7.** Calibration, gain matching, background subtraction; full 3-hour run ([results](docs/spectrum.md))
- [x] **8.** Benchmark against a naive count and an emulated analog chain ([results](docs/benchmark.md))
- [x] **9.** Every stage selectable from a run file ([docs](docs/configuration.md))
- [x] **10.** Streaming, live spectrum, signal simulator with truth ([streaming](docs/streaming.md), [simulation](docs/simulation.md))
- [x] **11.** Browser app and packaging ([docs](docs/app.md))

### 🔭 Possible next steps

- [ ] 🧪 Full end-to-end testing of the app
- [ ] ⚡ A faster pile-up solver, for real-time 1.8 MS/s on 4–6-core PCs
- [ ] 📉 A Compton-edge continuum model for unbiased FWHM
- [ ] 🔬 Re-characterisation at higher sampling rates
- [ ] 🧬 Isotope identification against nuclide libraries
