# Calibration packs

The sasnet core is complete without any calibration: every coefficient has a
default derived in [theory.md](theory.md). Anyone can use their own coefficients for free (template:
[../examples/my_coefficients.json](../examples/my_coefficients.json)). A **calibration pack** adds
coefficients fitted to test data or reference calculations — for a company,
an engine family or a project.

Rules:

- A pack only fills parameters the network file leaves open. A value written
  in the network file always wins.
- Every filled value, and the pack name and version, are listed in the report
  ("Coefficients" line and "Calibrated parameters" table).
- Without a pack the report says "first-principles defaults".

## Using a pack

```bash
python -m sasnet --list-calibrations
python -m sasnet network.json report.md --calibration <pack-name | rules.json>
```

```python
from sasnet import Network
net = Network.from_json("network.json", calibration="rules.json")
```

A network file may also name its pack: `"calibration": "<pack-name>"`.

## Open rules format (JSON)

```json
{"name": "my-pack", "version": "1.0",
 "rules": {
   "labyrinth": {"cd": {"c0": 0.80, "per": {"clearance_mm": -0.2}, "min": 0.5, "max": 0.9}},
   "orifice":   {"cd": 0.62}
 }}
```

`value = c0 + Σ coef · variable`, clipped to `[min, max]`. Keys are element
types (`orifice`, `labyrinth`, `channel`, `duct`, `fixed_zeta`) and their
parameters. Variables are the element's own parameters plus `clearance_mm`,
`c_over_s`, `d_mm`, `h_mm`, `length_mm`. The numbers above are placeholders,
not recommendations.

## Python packs

Subclass `sasnet.Calibration` (or `JsonCalibration`), override `fill(kind, params)`,
and register the class under the entry-point group `sasnet.calibrations`:

```toml
[project.entry-points."sasnet.calibrations"]
my-pack = "my_package:MyCalibration"
```

## Turbo Technologies coefficients

Turbo Technologies s.r.o. offers its calibrated coefficient pack (`turbotech`)
with its validation record:

- **Trial** — paid, up to 2 calculations within 3 days;
- **Project licence** — one project, its full duration;
- **Annual licence** — all projects of the licensee, with updates.

The pack needs a signed licence file (`$SASNET_TT_LICENCE` or
`~/.sasnet/turbotech.lic`); the licensee, plan and validity are printed in
every report. Contact: turbotechnologies@email.cz
