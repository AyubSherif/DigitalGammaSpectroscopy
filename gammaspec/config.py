"""Run configuration: one TOML file selects the input, the stages and every parameter (ADR-0011).

    gammaspec config init run.toml          # write the defaults, with a comment on each setting
    gammaspec config check run.toml         # validate; print warnings and the effective config
    gammaspec run run.toml --set pipeline.detect.threshold=6

The file mirrors the dataclasses of the code: ``[pipeline.preprocess]`` is
``PreprocessConfig``, ``[pipeline.detect]`` is ``DetectConfig``, and so on. Keys left out
keep their defaults, so a run file can be as short as::

    output = "results/run1"
    [input]
    kind = "raw"
    source = "ADC_raw_output"
    background = "ADC_raw_output_background"

Loading is strict: an unknown key is an error (with a "did you mean"), values are type
checked, choices and ranges are validated, and stage combinations that do not work are
reported (``check``). The same functions serve the command line and, later, a UI.
"""

from __future__ import annotations

import dataclasses
import difflib
import inspect
import json
import re
import tomllib
import types
import typing
from dataclasses import dataclass, field
from pathlib import Path

from .characterize import CharacterizeConfig
from .pipeline import PipelineConfig

SCHEMA_VERSION = 1


# ----------------------------------------------------------------------------------------------
# The run file
# ----------------------------------------------------------------------------------------------

@dataclass
class InputConfig:
    kind: str = "store"                   # "raw": folders of ASCII captures (converted first) | "store": HDF5 capture stores
    source: str = ""                      # source run: raw folder or .h5 store
    background: str = ""                  # background run, same kind ("" = no background subtraction)
    limit: int = 0                        # only the first N captures of each run, in time order (0 = all)
    pattern: str = "raw_data_*"           # raw input: capture file names
    reuse_converted: bool = True          # raw input: reuse a store already converted into the output folder


@dataclass
class TemplateConfig:
    path: str = ""                        # pulse template (.h5); "" = build one with the characterisation step
    build_from: str = "background"        # run to build it from: "background" (fewer pulses, cleaner) | "source"
    characterize: CharacterizeConfig = field(default_factory=CharacterizeConfig)


@dataclass
class SpectrumStepConfig:
    enabled: bool = True                  # calibrate, build spectra, subtract the background
    bin_kev: float = 5.0                  # bin width (keV)
    e_max_kev: float = 3000.0             # upper end of the spectrum (keV)
    reliable_above_kev: float = 120.0     # below this the background run is marked unreliable


@dataclass
class BenchmarkStepConfig:
    enabled: bool = False                 # compare with a naive find_peaks count and an emulated analog chain
    captures: int = 100                   # captures used (0 = all; the analog emulation runs at ~real time)


@dataclass
class RunConfig:
    output: str = "results"               # output folder (stores, template, events, spectrum, benchmark)
    workers: int = 0                      # parallel worker processes (0 = cores - 1, at most 8)
    overwrite: bool = False               # recompute outputs that already exist
    input: InputConfig = field(default_factory=InputConfig)
    template: TemplateConfig = field(default_factory=TemplateConfig)
    pipeline: PipelineConfig = field(default_factory=PipelineConfig)
    spectrum: SpectrumStepConfig = field(default_factory=SpectrumStepConfig)
    benchmark: BenchmarkStepConfig = field(default_factory=BenchmarkStepConfig)


class ConfigError(ValueError):
    """Invalid configuration: every problem found, one per line."""


# ----------------------------------------------------------------------------------------------
# Rules beyond types
# ----------------------------------------------------------------------------------------------

CHOICES: dict[str, tuple] = {
    "input.kind": ("raw", "store"),
    "template.build_from": ("background", "source"),
    "pipeline.preprocess.polarity": (-1, 1),
    "template.characterize.polarity": (-1, 1),
    "pipeline.detect.method": ("optimal", "matched", "sharpened", "edge"),
    "pipeline.events.snippets": ("none", "pileup", "all"),
}

# must be > 0 (all other numbers must be >= 0)
POSITIVE = {
    "pipeline.preprocess.trigger_threshold", "pipeline.detect.threshold", "pipeline.detect.mask_threshold",
    "pipeline.detect.edge_threshold", "pipeline.preprocess.lines.block", "pipeline.preprocess.lines.ratio_threshold",
    "pipeline.preprocess.lines.median_width", "pipeline.preprocess.tracker.window",
    "pipeline.solver.search_threshold", "pipeline.solver.max_shift", "pipeline.solver.acov_lags",
    "pipeline.events.cluster_gap", "spectrum.bin_kev", "spectrum.e_max_kev",
    "template.characterize.trigger_threshold", "template.characterize.resolution",
}


# ----------------------------------------------------------------------------------------------
# dict <-> dataclass
# ----------------------------------------------------------------------------------------------

def _hints(cls) -> dict:
    return typing.get_type_hints(cls)


def _optional(tp):
    """``X | None`` -> X, else None."""
    if typing.get_origin(tp) in (typing.Union, types.UnionType):
        args = [a for a in typing.get_args(tp) if a is not type(None)]
        if len(args) == 1:
            return args[0]
    return None


def _coerce(value, tp, where: str, errors: list[str]):
    inner = _optional(tp)
    if inner is not None:
        if value is None:                                   # JSON null for an optional setting
            return None
        tp = inner
    if dataclasses.is_dataclass(tp):
        if not isinstance(value, dict):
            errors.append(f"{where}: expected a table [{where}], got {value!r}")
            return None
        return _from_dict(tp, value, where, errors)
    if tp is bool:
        if isinstance(value, bool):
            return value
        errors.append(f"{where}: expected true or false, got {value!r}")
    elif tp is int:
        if isinstance(value, int) and not isinstance(value, bool):
            return value
        if isinstance(value, float) and value.is_integer():
            return int(value)
        errors.append(f"{where}: expected an integer, got {value!r}")
    elif tp is float:
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            return float(value)
        errors.append(f"{where}: expected a number, got {value!r}")
    elif tp is str:
        if isinstance(value, str):
            return value
        errors.append(f"{where}: expected a string, got {value!r}")
    elif tp is tuple or typing.get_origin(tp) is tuple:
        if isinstance(value, (list, tuple)) and all(isinstance(v, (int, float)) and not isinstance(v, bool)
                                                   for v in value):
            return tuple(float(v) for v in value)
        errors.append(f"{where}: expected a list of numbers, got {value!r}")
    else:                                                   # pragma: no cover - new field type
        errors.append(f"{where}: unsupported setting type {tp!r}")
    return None


def _from_dict(cls, d: dict, where: str, errors: list[str]):
    hints = _hints(cls)
    names = [f.name for f in dataclasses.fields(cls)]
    kwargs = {}
    for k, v in d.items():
        key = f"{where}.{k}" if where else k
        if k not in names:
            close = difflib.get_close_matches(k, names, n=1)
            errors.append(f"{key}: unknown setting" + (f" (did you mean '{close[0]}'?)" if close else
                                                      f" (valid: {', '.join(names)})"))
            continue
        c = _coerce(v, hints[k], key, errors)
        if c is not None:
            kwargs[k] = c
    try:
        return cls(**kwargs)
    except TypeError as exc:                                # pragma: no cover
        errors.append(f"{where or 'run file'}: {exc}")
        return cls()


def _validate_values(obj, where: str, errors: list[str]) -> None:
    for f in dataclasses.fields(obj):
        key = f"{where}.{f.name}" if where else f.name
        v = getattr(obj, f.name)
        if dataclasses.is_dataclass(v):
            _validate_values(v, key, errors)
            continue
        if key in CHOICES:
            if v not in CHOICES[key]:
                errors.append(f"{key}: {v!r} is not one of {', '.join(map(repr, CHOICES[key]))}")
            continue
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            continue
        if key in POSITIVE and not v > 0:
            errors.append(f"{key}: must be > 0, got {v!r}")
        elif v < 0:
            errors.append(f"{key}: must be >= 0, got {v!r}")


def from_dict(d: dict, cls=RunConfig):
    """Build and validate a config from a (possibly partial) nested dict. Raises ConfigError."""
    d = dict(d)
    version = d.pop("schema_version", SCHEMA_VERSION)
    errors: list[str] = []
    if version != SCHEMA_VERSION:
        errors.append(f"schema_version: {version!r} is not supported (this version reads {SCHEMA_VERSION})")
    cfg = _from_dict(cls, d, "", errors)
    prefix = "pipeline" if cls is PipelineConfig else ""
    _validate_values(cfg, prefix, errors)
    if errors:
        raise ConfigError("invalid configuration:\n  " + "\n  ".join(errors))
    return cfg


def to_dict(cfg) -> dict:
    d = dataclasses.asdict(cfg)
    return {"schema_version": SCHEMA_VERSION, **d} if isinstance(cfg, RunConfig) else d


# ----------------------------------------------------------------------------------------------
# Files
# ----------------------------------------------------------------------------------------------

def _parse_value(text: str):
    """A --set value: TOML syntax (6, 1e-3, true, "x", [1, 2]); anything else is a plain string."""
    try:
        return tomllib.loads(f"v = {text}")["v"]
    except tomllib.TOMLDecodeError:
        return text


def set_key(d: dict, assignment: str) -> None:
    """Apply ``a.b.c=value`` to a nested dict (in place)."""
    if "=" not in assignment:
        raise ConfigError(f"--set {assignment!r}: expected KEY=VALUE, e.g. pipeline.detect.threshold=6")
    key, text = assignment.split("=", 1)
    parts = key.strip().split(".")
    node = d
    for p in parts[:-1]:
        node = node.setdefault(p, {})
        if not isinstance(node, dict):
            raise ConfigError(f"--set {key}: '{p}' is a value, not a table")
    node[parts[-1]] = _parse_value(text.strip())


def read_dict(path) -> dict:
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    try:
        return json.loads(text) if path.suffix.lower() == ".json" else tomllib.loads(text)
    except (tomllib.TOMLDecodeError, json.JSONDecodeError) as exc:
        raise ConfigError(f"{path}: not valid {'JSON' if path.suffix.lower() == '.json' else 'TOML'}: {exc}") from exc


def load(path=None, overrides: list[str] | None = None, cls=RunConfig):
    """Read a TOML (or .json) run file, apply ``KEY=VALUE`` overrides, validate."""
    d = read_dict(path) if path else {}
    for s in overrides or []:
        set_key(d, s)
    return from_dict(d, cls)


_COMMENT = re.compile(r"^\s+(\w+)\s*:[^=#]+=[^#]*#\s*(.+?)\s*$")


def _comments(cls) -> dict[str, str]:
    try:
        src = inspect.getsource(cls)
    except (OSError, TypeError):                            # pragma: no cover
        return {}
    return {m.group(1): m.group(2) for line in src.splitlines() if (m := _COMMENT.match(line))}


def _toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    return json.dumps(str(v), ensure_ascii=False)          # a valid TOML basic string


def dumps(cfg, comments: bool = True) -> str:
    """The config as TOML, one table per dataclass, with each field's code comment."""
    lines: list[str] = []
    if isinstance(cfg, RunConfig):
        lines += ["# gammaspec run configuration (docs/configuration.md)", f"schema_version = {SCHEMA_VERSION}"]

    def table(obj, name: str):
        notes = _comments(type(obj)) if comments else {}
        scalars = [f for f in dataclasses.fields(obj) if not dataclasses.is_dataclass(getattr(obj, f.name))]
        nested = [f for f in dataclasses.fields(obj) if dataclasses.is_dataclass(getattr(obj, f.name))]
        if name and scalars:
            lines.extend(["", f"[{name}]"])
        for f in scalars:
            v = getattr(obj, f.name)
            note = f"  # {notes[f.name]}" if f.name in notes else ""
            lines.append(f"# {f.name} = (not set){note}" if v is None else f"{f.name} = {_toml_value(v)}{note}")
        for f in nested:
            table(getattr(obj, f.name), f"{name}.{f.name}" if name else f.name)

    table(cfg, "pipeline" if isinstance(cfg, PipelineConfig) else "")
    return "\n".join(lines).lstrip("\n") + "\n"


def save(cfg, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".json":
        path.write_text(json.dumps(to_dict(cfg), indent=2), encoding="utf-8")
    else:
        path.write_text(dumps(cfg), encoding="utf-8")
    return path


def pipeline_from_event_file(path) -> PipelineConfig:
    """The pipeline configuration an event file was made with (to reproduce or compare a run)."""
    from .events import EventFile

    d = EventFile(path).attrs.get("config")
    if not isinstance(d, dict):
        raise ConfigError(f"{path}: no configuration stored")
    return from_dict(d, PipelineConfig)


# ----------------------------------------------------------------------------------------------
# Combinations of stages
# ----------------------------------------------------------------------------------------------

def check(cfg: RunConfig, paths: bool = True) -> list[str]:
    """Problems between settings. Errors raise ConfigError; returns the warnings."""
    errors, warn = [], []
    p, i = cfg.pipeline.preprocess, cfg.input
    if not i.source:
        errors.append("input.source: not set (a raw capture folder or an .h5 store)")
    if not cfg.template.path and cfg.template.build_from == "background" and not i.background:
        errors.append("template: no template.path, and template.build_from = 'background' but there is no "
                      "input.background (set one, or build_from = 'source')")
    if paths:
        for key, val in (("input.source", i.source), ("input.background", i.background),
                         ("template.path", cfg.template.path)):
            if not val:
                continue
            q = Path(val)
            want_dir = key != "template.path" and i.kind == "raw"
            if not q.exists():
                errors.append(f"{key}: {val} does not exist")
            elif want_dir and not q.is_dir():
                errors.append(f"{key}: {val} is not a folder (input.kind = 'raw')")
            elif not want_dir and q.is_dir():
                errors.append(f"{key}: {val} is a folder; expected an .h5 file"
                              + (" (raw captures need input.kind = 'raw')" if key != "template.path" else ""))
    if p.pole_zero and not p.baseline_tracker:
        warn.append("pole_zero without baseline_tracker: the pole-zero filter amplifies slow baseline drift "
                    "(on the 2021 data this is worse than no processing)")
    restored = p.pole_zero and p.baseline_tracker
    if cfg.pipeline.solver.enabled and not restored:
        warn.append("pipeline.solver needs pole_zero + baseline_tracker; it will be skipped "
                    "(pile-up clusters keep their detector amplitudes, flagged PILEUP)")
    if cfg.pipeline.detect.method == "edge" and cfg.pipeline.detect.cancel:
        warn.append("detect.cancel has no effect with method = 'edge'")
    if cfg.spectrum.enabled and not i.background:
        warn.append("no input.background: the spectrum is calibrated but not background-subtracted")
    if errors:
        raise ConfigError("invalid configuration:\n  " + "\n  ".join(errors))
    return warn
