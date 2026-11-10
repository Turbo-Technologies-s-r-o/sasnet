"""Rotor axial (thrust) force from cavity pressures.

Each cavity face on a rotor is an annulus r..R loaded by the pressure of a
network node. ``side`` = +1 if that pressure pushes the rotor in the
positive axial direction, -1 otherwise. The blade (aerodynamic) force is an
input with its own sign.

Default: pressure uniform over the annulus (as in hand calculations).
Optional swirl: in a rotating cavity with core swirl ratio K the pressure
varies with radius as dp/dr = rho (K omega)^2 r, i.e.

    p(r) = p_ref + 0.5 rho (K omega)^2 (r^2 - r_ref^2)

where p_ref is the node pressure, taken at radius r_ref. The force is then
integrated exactly over the annulus.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

KGF = 9.80665  # N per kgf


@dataclass
class Face:
    name: str
    r: float            # inner radius, m
    R: float            # outer radius, m
    node: str | None = None
    side: int = +1
    p: float | None = None      # fixed pressure, Pa (instead of a node)
    swirl: float = 0.0          # core swirl ratio K (0 = uniform pressure)
    r_ref: float | None = None  # radius where the node pressure applies (default: R)
    pi_value: float = math.pi   # 3.14 reproduces older hand calculations exactly

    @property
    def area(self) -> float:
        return self.pi_value * (self.R ** 2 - self.r ** 2)

    def force(self, p_node: float, rho: float = 0.0, omega: float = 0.0) -> float:
        p_ref = self.p if self.p is not None else p_node
        F = p_ref * self.area
        if self.swirl and omega and rho:
            c = 0.5 * rho * (self.swirl * omega) ** 2
            r_ref = self.R if self.r_ref is None else self.r_ref
            # integral of 2 pi r c (r^2 - r_ref^2) dr from r to R
            F += self.pi_value * c * (0.5 * (self.R ** 4 - self.r ** 4) - r_ref ** 2 * (self.R ** 2 - self.r ** 2))
        return self.side * F


@dataclass
class Rotor:
    name: str
    faces: list = field(default_factory=list)
    blade_force: float = 0.0    # N, signed
    omega: float = 0.0          # rad/s, for the swirl option

    def face(self, *args, **kw) -> "Rotor":
        self.faces.append(Face(*args, **kw))
        return self

    def breakdown(self, p: dict | None = None, T: dict | None = None, R_gas: float = 287.05) -> list[tuple]:
        rows = [("blade force", None, None, None, self.blade_force)]
        for f in self.faces:
            pn = f.p if f.p is not None else p[f.node]
            rho = pn / (R_gas * T[f.node]) if (f.swirl and T and f.node) else 0.0
            rows.append((f.name, f.r, f.R, pn, f.force(pn, rho, self.omega)))
        return rows

    def total(self, p: dict | None = None, T: dict | None = None) -> float:
        return sum(r[-1] for r in self.breakdown(p, T))


def thrust_table(rotors: list[Rotor], p: dict | None = None, T: dict | None = None, in_kgf: bool = True) -> str:
    """Markdown table of face forces and rotor totals."""
    unit, k = ("kgf", 1 / KGF) if in_kgf else ("N", 1.0)
    out = []
    grand = 0.0
    for rot in rotors:
        out += [f"### {rot.name}", "", f"| Face | r, mm | R, mm | p, kPa | F, {unit} |", "|---|---|---|---|---|"]
        for name, r, R, pn, F in rot.breakdown(p, T):
            rr = "" if r is None else f"{r*1e3:.2f}"
            RR = "" if R is None else f"{R*1e3:.2f}"
            pp = "" if pn is None else f"{pn/1e3:.2f}"
            out.append(f"| {name} | {rr} | {RR} | {pp} | {F*k:.1f} |")
        tot = rot.total(p, T)
        grand += tot
        out += [f"| **total {rot.name}** | | | | **{tot*k:.1f}** |", ""]
    out.append(f"**Sum of all rotors: {grand*k:.1f} {unit}**")
    return "\n".join(out) + "\n"
