# Development guide

## Environment

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows (PowerShell: .venv\Scripts\Activate.ps1)
source .venv/bin/activate         # Linux/macOS
pip install -e ".[app,dev]"
```

`pyproject.toml` declares the package and its minimum dependency versions.
`requirements.txt` pins the exact versions used during development, so results
can be reproduced:

```bash
pip install -r requirements.txt && pip install -e . --no-deps
```

## Repository layout

```text
gammaspec/          the package (see docs/architecture.md#module-map)
tests/              pytest suite
docs/               documentation and ADRs
data/               generated stores and indexes (gitignored)
ADC_raw_output*/    raw input data (gitignored, never committed)
```

Raw data, generated stores (`*.h5`, `*.npy`), `.venv/` and the reference PDFs
are excluded by `.gitignore`. Line endings are normalised to LF by `.gitattributes`.

## Tests

```bash
python -m pytest            # whole suite
python -m pytest -k store   # a subset
```

* **Unit tests** build small synthetic capture files in the exact 2021 format
  (`tests/conftest.py::write_capture`), with file numbers deliberately out of
  time order. They need no real data.
* **Real-data tests** are marked `skipif` and run only when
  `ADC_raw_output/` exists next to the repository.
* Every new processing stage must come with (a) a test on synthetic data with a
  known answer and (b) a diagnostic plot or summary for the real data.

## Walkthrough notebook

`docs/walkthrough.ipynb` steps through every function in call order on the real data.
It is generated and executed by `docs/build_walkthrough.py` (needs the local `data/`
folder and the `gammaspec` Jupyter kernel: `python -m ipykernel install --user --name gammaspec`).
Rebuild it after changing a function it shows:

```bash
python docs/build_walkthrough.py
```

## Conventions

* **Units:** volts and seconds in public APIs; ADC codes (int32) for storage.
  State units in docstrings and attribute names (`dt_s`, `lsb_volts`, `time_s`).
* **Time order:** order captures by header timestamp, never by file name.
* **Lossless by default:** a step that discards information (thresholds, event
  extraction) must be an explicit, configurable stage.
* **Fail loudly:** malformed input raises `RawFormatError` / `ValueError` rather
  than being silently fixed.
* **Performance:** stages must keep up with the ADC rate (1.8 MS/s) on one core,
  ideally with a 5× margin. Use NumPy vectorisation first; use Numba for loops
  that can't be vectorised.
* **Docs:** a user-visible change updates `CHANGELOG.md`; a significant design
  choice gets an ADR in `docs/adr/`.

## Commits

Small, focused commits with imperative subject lines ("Add …", "Fix …"), and a
body explaining *why* when it isn't obvious.
