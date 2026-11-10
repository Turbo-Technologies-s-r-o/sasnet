"""Plain-text and markdown reports of a solved network."""
from __future__ import annotations

from .network import Network, Result
from . import __version__


def _calibration_line(net) -> str:
    cal = getattr(net, "calibration", None)
    if cal is None:
        return "Coefficients: first-principles defaults (no calibration pack).\n"
    return f"Coefficients: calibration pack **{cal.label()}**, {len(cal.applied)} values filled.\n"


def markdown(net: Network, r: Result, title: str = "Secondary air system network") -> str:
    lines = [f"# {title}", "",
             f"sasnet {__version__}. Gas: R = {net.gas.R} J/(kg K), gamma = {net.gas.gamma}. "
             "Node pressures are total pressures in plenums; elements adiabatic.", "",
             _calibration_line(net),
             f"Converged: **{'yes' if r.converged else 'NO'}**, max relative mass imbalance "
             f"{r.residual_max:.1e}, outer iterations {r.iterations}.", "",
             "## Nodes", "", "| Node | Type | p, kPa | T, K |", "|---|---|---|---|"]
    for n in net.nodes.values():
        kind = "boundary" if n.boundary else "cavity"
        lines.append(f"| {n.name} | {kind} | {r.p[n.name]/1e3:.3f} | {r.T[n.name]:.2f} |")
    lines += ["", "## Elements", "", "| Element | From → To | m, g/s | Description |",
              "|---|---|---|---|"]
    for e in net.elements:
        m = r.m[e.name]
        a, b = (e.a, e.b) if m >= 0 else (e.b, e.a)
        lines.append(f"| {e.name} | {a} → {b} | {abs(m)*1e3:.3f} | {e.describe()} |")
    cal = getattr(net, "calibration", None)
    if cal is not None and cal.applied:
        lines += ["", "## Calibrated parameters", "", "| Element | Parameter | Value |", "|---|---|---|"]
        lines += [f"| {e} | {k} | {v:.4g} |" for e, k, v in cal.applied]
    lines += ["", "## Checks", ""]
    lines += [f"- {w}" for w in r.warnings] or ["- none"]
    return "\n".join(lines) + "\n"
