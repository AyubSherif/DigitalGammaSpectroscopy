# ADR-0013: A Streamlit browser app as a front end for the run file, running runs as a separate process

* **Status:** Accepted
* **Date:** 2026-11-24

## Context

Users should be able to configure and run the program without code or the command line.
They should choose the input (raw files, a store, a replay, the simulator, a live stream),
switch stages, and watch a live spectrum. The configuration layer (ADR-0011) already
defines every setting, validates it and runs it; the streaming layer (ADR-0012) already
writes a live spectrum file.

Options:
* (a) A desktop GUI (Qt). Native, but a large dependency, and a lot of code for widgets.
* (b) A web front end with a custom server (FastAPI + JavaScript). Flexible, but two code
  bases.
* (c) **Streamlit**: Python only, runs in the browser, small code.

For running: in the page's own process (threads), or as a separate process.

## Decision

* (c), in one file (`gamma_spectrometry/app.py`), launched by `gamma-spectrometry app`. By default it listens
  on `localhost` only.
* The page edits a `RunConfig`:
  * common settings are widgets, each bound to its dataclass field;
  * every other setting is in a TOML editor showing the same run file.

  Validation is `config.from_dict` + `config.check`, the same as on the command line.
* **Run** writes `<output>/ui_run.toml` and starts `python -m gamma_spectrometry run` on it as a
  **separate process**:
  * the page stays responsive;
  * **Stop** can end the process tree;
  * processing uses its own worker pool, as from the command line;
  * the run is reproducible from the saved file.
* The live view polls `live/<run>/live_spectrum.json` and the log every 2 s (a Streamlit
  fragment). It does not share memory with the run.

## Consequences

* The app adds no processing code. A new setting appears in the TOML editor automatically,
  and it needs a widget only if it should be prominent.
* Streamlit becomes an optional dependency (`pip install "gamma-spectrometry[app]"`). The rest of the
  program does not import it.
* Live updates arrive in steps of 2 s (stream `update_s`), which suits a spectrum display.
* The page is tested without a browser (`streamlit.testing.v1.AppTest`).
* A standalone executable is left for later (PyInstaller bundle, large).
