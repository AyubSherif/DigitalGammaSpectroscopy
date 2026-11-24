# ADR-0011: One TOML run file, mirroring the config dataclasses, run by an incremental runner

* **Status:** Accepted
* **Date:** 2026-11-24
* **Implements:** [ADR-0003](0003-modular-optional-stages.md) (configuration and validation part)

## Context

ADR-0003 asks for user-selectable stages, a configuration that lists them with their
parameters, checks of stage requirements at configuration time, and outputs that record
their configuration. The stages now exist, and their settings already live in nested
dataclasses (`PipelineConfig` → `PreprocessConfig`, `DetectConfig`, `EventConfig`,
`SolverConfig`, …). The configuration will later be edited by a UI (step 11) and drive a
streaming input (step 10).

Options:
* (a) Command-line flags only. They do not scale to ~100 parameters, and a run is not
  reproducible from a file.
* (b) A separate configuration schema (e.g. pydantic or JSON Schema) next to the
  dataclasses. That means two definitions to keep in sync, and a new dependency.
* (c) A file format that **mirrors the dataclasses**: they are the schema.

Format: TOML (readable, comments, read by the standard library since Python 3.11). JSON is
accepted too, which suits a UI.

## Decision

* (c): `gamma_spectrometry/config.py` maps TOML tables to the dataclasses, recursively. A
  `RunConfig` adds the input, template, spectrum and benchmark steps around `PipelineConfig`.
* Strict loading:
  * unknown keys are errors with a "did you mean" suggestion;
  * values are type-checked against the dataclass annotations;
  * choices and ranges are listed in `config.CHOICES` / `config.POSITIVE`;
  * every problem is reported at once.
* `config.check` validates stage combinations before running. Errors are for missing or wrong
  inputs; warnings are for stage combinations that work but are not recommended.
* `config init` writes the defaults. The comment of each field is taken from its code
  comment, so the code stays the single source.
* `gamma_spectrometry/runner.py` runs the steps (input → template → process → spectrum →
  benchmark). It writes the effective configuration and a stamp per output (inputs +
  settings): an output is reused only when its stamp is unchanged.
* The stage order stays fixed. Each stage depends on the output of the one before it, so
  only on/off and the parameters are configurable.

## Consequences

* A new parameter in a dataclass appears in the run file, the validation and (later) the UI
  without further code. A parameter with restricted values needs one line in `CHOICES` or
  `POSITIVE`.
* Every result can be traced: `config.toml` in the output folder, the configuration in
  each event file (`gamma-spectrometry config show`), and the stamps.
* Changing one setting reruns only the steps it affects. Spectrum and benchmark always rerun
  (they are fast compared with processing).
* `schema_version` lets later versions convert old run files.
* Not covered yet: carrying stage state across chunks of a continuous stream (step 10).
  Choosing the pile-up algorithm in the pipeline is not offered either: only the WLS solver
  is wired in, and the alternatives were measured to be worse.
