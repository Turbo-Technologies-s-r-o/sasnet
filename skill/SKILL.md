---
name: gte-secondary-air-network
description: "Hydraulic calculation of gas turbine secondary air system (SAS) networks (cavities, holes, labyrinth seals, ducts, rotor thrust) with sasnet. Build the network, solve, check, report."
---

# Gas turbine secondary air networks (sasnet)

The user is an experienced engine engineer: keep it engineering-level, no textbook tone.
The result is a network input file, the solution, a report and an honest list of limitations.
Judge results by trends relative to a known reference point, not by absolute numbers alone.
A human signs off the calculation.

## 1. Collect the inputs (ask for everything missing in one list)

- System scheme: cavities (nodes) and the connections between them (elements); where air is bled, where it is dumped
- Boundary nodes: total pressure and total temperature (compressor bleeds, main gas path, oil sumps, ambient)
- For every element — type and geometry:
  - orifice: area, Cd (and where it comes from)
  - labyrinth: clearance, radius, number of teeth, pitch, straight or stepped, Cd
  - duct (holes, slots, annular gaps): shape, size, length, count, area ratios at inlet and outlet, rotation (ω, radii), bends, heating
  - channel: area, hydraulic diameter, length, roughness, sum of local losses ΣK with the source of each value
- Operating point: speed, throttle point
- Reference: a test or a previous calculation for this or a similar system

Anything missing is marked in the input file as an ASSUMPTION (field `source`) — never invent values silently.
Handbook coefficients only from a cited source (diagram or page number). Never from memory.

## 2. Build the input file

JSON format (example: `examples/demo_network.json`):

```json
{"nodes": [{"name": "hp_bleed", "boundary": true, "p": 1500000, "T": 720},
           {"name": "front_cavity"}],
 "elements": [{"type": "orifice", "name": "disk_holes", "a": "front_cavity", "b": "rim_cavity",
               "area": 1.6e-4, "cd": 0.62, "source": "Cd from rig test, project X"}]}
```

SI units: Pa, K, m, m². Use meaningful node and element names, as on the user's scheme.
Show the user the network (table of nodes and links, or a diagram) before solving.

## 3. Solve

```bash
python -m sasnet network.json report.md
```

or from Python: `Network.from_json(...).solve()`.

Coefficients: without a calibration the core uses first-principles defaults. The user's own
coefficients go in a JSON rules file (`examples/my_coefficients.json`):
`--calibration my_coefficients.json`. An installed licensed pack is used by name
(`--calibration <pack>`; list with `python -m sasnet --list-calibrations`). A calibration only
fills values the network file leaves open. Always show the user the report's "Coefficients" line
and "Calibrated parameters" table, so it is clear which coefficients the result is based on.
See docs/calibration.md.

## 4. Check before showing any numbers

- Convergence: `converged = yes`, relative mass imbalance ≤ 1e-6 in every cavity
- Every warning in the report's Checks section — explain each one to the user:
  - choked orifice — physics, not an error; check whether it is expected
  - labyrinth last tooth choked or close to it — the pure Martin formula would overpredict the flow
  - Mach > 0.3 in a channel or duct — outside the low-Mach model; use `compressible` or split the branch
  - rotating duct without swirl-entry loss (`c_rot = 0`) — set it from data
- Flow directions: unexpected reverse flow is the first thing to show
- Cavity mixing temperatures look reasonable
- Compare with the reference in trends: more or less flow, by how much, and why

## 5. Deliver

1. The network input file (JSON) with `source` fields for every coefficient
2. The report (markdown or PDF): nodes, flows, checks, comparison with the reference
3. The list of limitations below and which of them affect this calculation

## Rotor axial thrust

`sasnet.axial.Rotor`: cavity faces (r, R, node, side ±) and blade force; pressures come from the
network solution. By default cavity pressure is uniform over radius, as in hand calculations
(`pi_value=3.14` reproduces older tables exactly); the swirl option K integrates the pressure rise
towards the rim. Deliver a table by face and totals per rotor for every operating point.

## Limitations (always state them when they affect the result)

- Sharp-edged inlet constant K_in = 0.5 (Weisbach) is the main uncertainty for holes; first
  principles tend to give somewhat higher hole flows than calibrated references
- Short channel to thin orifice transition over the reattachment length 2·d_h is an order-of-magnitude estimate
- Swirl-entry loss into rotating holes: `c_rot` defaults to 0 (off), set it from data
- No heat transfer along channels (only an imposed dT); cavity temperatures can be imposed (`fix_T`)
- Rotating cavities (free/forced vortex, disk pumping, windage) are not modelled yet
- Labyrinth: Martin + jet carry-over + isentropic last tooth (choking included)
- Cavity pressures are total pressures (velocity in the cavity neglected)
- Cd and ΣK are inputs or come from a calibration; the library does not fit them
