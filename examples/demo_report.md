# examples/demo_network.json

sasnet 0.3.0. Gas: R = 287.05 J/(kg K), gamma = 1.4. Node pressures are total pressures in plenums; elements adiabatic.

Converged: **yes**, max relative mass imbalance 1.9e-09, outer iterations 8.

## Nodes

| Node | Type | p, kPa | T, K |
|---|---|---|---|
| hp_bleed | boundary | 1500.000 | 720.00 |
| lp_bleed | boundary | 650.000 | 460.00 |
| sump_vent | boundary | 160.000 | 400.00 |
| turbine_rim | boundary | 560.000 | 900.00 |
| front_cavity | cavity | 1237.541 | 720.00 |
| rim_cavity | cavity | 615.181 | 669.71 |

## Elements

| Element | From → To | m, g/s | Description |
|---|---|---|---|
| feed_pipe | hp_bleed → front_cavity | 345.346 | Channel A=0.0003 m2, Dh=19.5 mm, L=450 mm, sumK=2.3 |
| bearing_seal | front_cavity → sump_vent | 160.442 | Labyrinth straight N=5, c=0.25 mm, r=85 mm, Cd=0.7, k=1.633 |
| disk_holes | front_cavity → rim_cavity | 184.904 | Orifice A=0.00016 m2, Cd=0.62 |
| lp_inlet | lp_bleed → rim_cavity | 44.339 | Orifice A=0.00012 m2, Cd=0.65 |
| rim_seal | rim_cavity → turbine_rim | 229.243 | Labyrinth straight N=3, c=0.5 mm, r=210 mm, Cd=0.7, k=1.500 |

## Checks

- feed_pipe: Mach 0.36 > 0.3, low-Mach channel model outside validity
- bearing_seal: last tooth choked (p ratio ~0.280 < 0.528)
- disk_holes: choked
