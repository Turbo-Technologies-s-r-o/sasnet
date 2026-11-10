# sasnet

1D secondary air system (SAS) network solver for gas turbines and turbomachinery:
cooling flows, seal leakage, cavity pressures and rotor axial thrust.

By [Turbo Technologies s.r.o.](https://turbo-technologies-s-r-o.github.io), Prague.

## Why

Every loss model in sasnet is derived from first principles or classical
public-domain laws — Bernoulli, Euler, Kirchhoff, Borda–Carnot, Hagen–Poiseuille,
Prandtl, Colebrook, Saint-Venant–Wantzel. The derivations are published in
[docs/theory.md](docs/theory.md), so every number can be traced and argued with.
No handbook tables or proprietary data are in the code.

sasnet has been validated on gas turbine secondary air networks of 200+ branches
and is consistent with calculation tools used in EASA Part 21 and CS-E work.

## What it does

- Network of cavities (nodes) connected by flow elements; boundary nodes with
  fixed total pressure and temperature
- Cavity pressures from mass balance (damped Newton, element-wise Jacobian),
  temperatures from adiabatic mixing or imposed (`fix_T`)
- Elements:
  - `Duct` — round holes, slots, annular gaps: inlet contraction, friction by
    shape, exit (Borda–Carnot with velocity-profile coefficients), short channels,
    bends, heating, rotation (Euler pumping, swirl entry), compressibility
  - `Labyrinth` — Martin's formula with jet carry-over and isentropic last tooth
    (choking)
  - `Orifice` — compressible isentropic flow with discharge coefficient
  - `Channel` — friction plus local losses ΣK
  - `FixedZeta` — known loss coefficient (from CFD or test)
- `axial` — rotor axial thrust from the solved cavity pressures, optional
  core-swirl pressure profile
- Automatic checks: convergence, choking, Mach number
- JSON input, markdown report, command line, Claude skill (`skill/SKILL.md`)

## Quick start

```bash
pip install -e .
python -m sasnet examples/demo_network.json report.md
python -m sasnet examples/demo_network.json report.md --calibration examples/my_coefficients.json
pytest
```

The demo network uses invented numbers and does not describe a real engine.

## Your coefficients or ours

**Your own coefficients — free.** The core runs on first-principles defaults.
Anyone can add their own discharge and loss coefficients as a JSON rules file
(template: [examples/my_coefficients.json](examples/my_coefficients.json)):

```bash
python -m sasnet network.json report.md --calibration my_coefficients.json
```

A calibration only fills values the network file leaves open, and every filled
value is listed in the report. Details: [docs/calibration.md](docs/calibration.md).

**Turbo Technologies coefficients — paid.** Our calibrated coefficient pack,
fitted on gas turbine networks and delivered with its validation record:

- **Trial** — paid, up to 2 calculations within 3 days;
- **Project licence** — one project, its full duration;
- **Annual licence** — all projects of the licensee, with updates.

It installs next to the core and is used the same way
(`--calibration turbotech`). Contact: turbotechnologies@email.cz

## Limitations

No heat transfer along channels yet; rotating cavities only through swirl and
rotation terms; inlet losses at moderate Reynolds numbers and very short slots
are estimates. Loss coefficients for flows along disks and blades are inputs
(`FixedZeta`). See `docs/theory.md`, section on validation.

## Roadmap

- Heat transfer along channels
- Rotating cavities (free/forced vortex, pumping, pre-swirl)
- Rotating-orifice discharge coefficients, windage
- Real-gas properties, coupling to thermal models

Further tools for gas turbines and turbomachinery are in development.

## License

Apache License 2.0 — see [LICENSE](LICENSE).

Contact: turbotechnologies@email.cz
