# The app (step 11)

`gamma-spectrometry app` opens a page in your browser to run the whole program without the command
line: choose the input, switch stages on and off, run, and watch the live spectrum.

Code: `gamma_spectrometry/app.py`. Decision record: [ADR-0013](adr/0013-browser-app.md).

```bash
pip install -e ".[app]"         # once: adds Streamlit
gamma-spectrometry app                   # opens http://localhost:8501
gamma-spectrometry app --port 8600 --no-browser
```

The page runs on this computer only (`localhost`). `--address 0.0.0.0` makes it reachable
from other computers. Only do that on a trusted network: the page can read files and start
processing runs.

## The page

**1 · Input: where the signal comes from**

| Choice | You give |
|---|---|
| Raw ADC files | the folder of ASCII capture files (they are converted first) |
| Prepared HDF5 store | an `.h5` file from `gamma-spectrometry convert` |
| Recorded data as a stream | a raw folder or `.h5` store, read capture by capture (optionally at real time) |
| Simulated signal | nuclide, rates, durations, resolution, noise, gain, seed |
| Live stream from a file / FIFO | the file a digitiser writes to; sample rate, format, volts per unit |

A background run is optional; with it, the spectrum is background-subtracted. Each path
shows whether it was found.

**2 · Processing**
* **Template:** the pulse template file, or build one from the data.
* **Stages:** interference removal, pole-zero correction, baseline tracker, detection filter
  and threshold, pile-up solver, search for hidden pulses. Each has a one-line
  explanation (the ⓘ icon).
* **After processing:** calibration and background subtraction, and the benchmark.

**All settings.** The complete run file as text, every parameter with its explanation. Edit
it and press **Apply**: the widgets follow. You can also load a run file or download the
current one; the command line runs the same file (`gamma-spectrometry run run.toml`).

**3 · Run**
* **Before running:** the settings are checked as you change them. Errors (a missing file,
  an invalid value) disable **Run**; warnings (stage combinations that work badly) are shown.
* **Run** starts the run in the background. Every 2 s the page shows the latest log lines
  and the live spectrum (in keV once it has calibrated itself).
* **Stop** ends the run, including its worker processes.

**4 · Results**
* the steps and their times;
* the calibration and the line rates (source, background, net);
* the spectrum and calibration plots;
* for a simulation, the comparison with the truth (pulses found by energy, false
  detections, amplitude bias);
* for the benchmark, the comparison table and plot.

## Files

For a run, the app writes `<output>/ui_run.toml` (the run file it started) and
`<output>/ui_run.log` (its output). Everything else is the usual output of `gamma-spectrometry run`
([configuration.md](configuration.md)). To repeat a run from the command line:

```bash
gamma-spectrometry run results/ui_run/ui_run.toml
```

## Packaging

* `pip install -e ".[app]"` installs the program, the `gamma-spectrometry` command and the app.
  `.[app,dev]` adds the tests and the notebook tools.
* A standalone Windows executable (no Python installation) is possible with PyInstaller, but
  not provided yet. NumPy, SciPy, Numba and Streamlit make it a large bundle (~400 MB).

## Tests

`tests/test_app.py` drives the page without a browser (Streamlit's `AppTest`):
* rendering;
* that the widgets edit the run file;
* the TOML editor (apply, and error messages);
* a complete simulated run started with **Run**, with its results shown.
