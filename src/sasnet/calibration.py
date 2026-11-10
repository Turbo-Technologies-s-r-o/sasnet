"""Calibration packs: optional coefficient sets layered on the first-principles core.

The core never needs a calibration: every element has a default derived in
docs/theory.md. A calibration pack supplies project- or company-specific
coefficients (discharge coefficients, inlet constants, loss multipliers) that
were fitted against test or reference calculations.

A pack is applied when a network is built from a dict/JSON file. It only fills
parameters that the input does NOT set explicitly, so a value written in the
network file always wins. Every value a pack fills is recorded, and the pack's
name and version are printed in the report for traceability.

Two ways to provide a pack:

1. A JSON file with linear rules (open format, see ``JsonCalibration``)::

       {"name": "my-pack", "version": "1.0",
        "rules": {"labyrinth": {"cd": {"c0": 0.80, "per": {"clearance_mm": -0.2},
                                        "min": 0.5, "max": 0.9}}}}

2. A Python package that subclasses ``Calibration`` and registers it under the
   entry-point group ``sasnet.calibrations``; it is then loadable by name.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from importlib import metadata

ENTRY_POINT_GROUP = "sasnet.calibrations"


@dataclass
class Calibration:
    """Base class. Subclasses override ``fill``."""
    name: str = "none"
    version: str = ""
    licensee: str = ""
    applied: list = field(default_factory=list)   # (element, parameter, value)

    def fill(self, kind: str, params: dict) -> dict:
        """Return extra parameters for one element (kind = JSON type)."""
        return {}

    # used by Network.from_dict
    def apply(self, kind: str, params: dict) -> dict:
        extra = {k: v for k, v in self.fill(kind, dict(params)).items() if k not in params}
        for k, v in extra.items():
            self.applied.append((params.get("name", "?"), k, v))
        return {**params, **extra}

    def label(self) -> str:
        s = self.name + (f" {self.version}" if self.version else "")
        return s + (f", licensed to {self.licensee}" if self.licensee else "")


def _derived(kind: str, p: dict) -> dict:
    """Convenience variables available to JSON rules."""
    v = dict(p)
    if "clearance" in p:
        v["clearance_mm"] = p["clearance"] * 1e3
    if kind == "labyrinth" and "clearance" in p and "pitch" in p and p["pitch"]:
        v["c_over_s"] = p["clearance"] / p["pitch"]
    for key in ("d", "h", "length"):
        if key in p:
            v[key + "_mm"] = p[key] * 1e3
    return v


@dataclass
class JsonCalibration(Calibration):
    """Linear rules: value = c0 + sum(coef * variable), clipped to [min, max]."""
    rules: dict = field(default_factory=dict)

    def fill(self, kind: str, params: dict) -> dict:
        out = {}
        v = _derived(kind, params)
        for param, rule in self.rules.get(kind, {}).items():
            if isinstance(rule, (int, float)):
                out[param] = float(rule)
                continue
            val = rule.get("c0", 0.0)
            ok = True
            for var, coef in rule.get("per", {}).items():
                if var not in v:
                    ok = False
                    break
                val += coef * v[var]
            if not ok:
                continue
            if "min" in rule:
                val = max(rule["min"], val)
            if "max" in rule:
                val = min(rule["max"], val)
            out[param] = val
        return out

    @classmethod
    def from_file(cls, path: str) -> "JsonCalibration":
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        return cls(name=d.get("name", path), version=str(d.get("version", "")),
                   licensee=d.get("licensee", ""), rules=d.get("rules", {}))


def available() -> list[str]:
    """Names of installed calibration packs (entry points)."""
    return sorted(ep.name for ep in metadata.entry_points(group=ENTRY_POINT_GROUP))


def load(spec) -> Calibration | None:
    """None/"" -> no calibration; Calibration -> as is; *.json -> JsonCalibration;
    otherwise the name of an installed pack."""
    if spec is None or spec == "" or spec == "none":
        return None
    if isinstance(spec, Calibration):
        return spec
    if isinstance(spec, dict):
        return JsonCalibration(name=spec.get("name", "inline"), version=str(spec.get("version", "")),
                               licensee=spec.get("licensee", ""), rules=spec.get("rules", {}))
    if str(spec).lower().endswith(".json"):
        return JsonCalibration.from_file(spec)
    for ep in metadata.entry_points(group=ENTRY_POINT_GROUP):
        if ep.name == spec:
            obj = ep.load()
            return obj() if isinstance(obj, type) else obj
    raise LookupError(f"calibration pack '{spec}' not installed (available: {available() or 'none'})")
