"""Network assembly and solution.

Boundary nodes have fixed total pressure and total temperature.
Internal nodes (cavities) have unknown pressure; their total temperature
follows from adiabatic mixing of the inflows.

Solution: Newton-type solve of the mass balance at internal nodes for the
pressures (in log form, so pressures stay positive), wrapped in a
fixed-point loop for the mixing temperatures.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import json
import math

import numpy as np

from .gas import AIR, PerfectGas
from . import elements as el


@dataclass
class Node:
    name: str
    p: float | None = None   # Pa, total; fixed for boundary nodes
    T: float | None = None   # K, total; fixed for boundary nodes
    boundary: bool = False
    fixed_T: bool = False    # cavity with imposed total temperature


@dataclass
class Result:
    converged: bool
    message: str
    p: dict
    T: dict
    m: dict
    residual_max: float
    iterations: int
    warnings: list = field(default_factory=list)


class Network:
    def __init__(self, gas: PerfectGas = AIR):
        self.gas = gas
        self.nodes: dict[str, Node] = {}
        self.elements: list[el.Element] = []

    # --- building -------------------------------------------------------
    def boundary(self, name: str, p: float, T: float) -> "Network":
        self.nodes[name] = Node(name, p, T, boundary=True)
        return self

    def cavity(self, name: str, p_guess: float | None = None, T_guess: float | None = None,
               fix_T: bool = False) -> "Network":
        if fix_T and T_guess is None:
            raise ValueError(f"cavity {name}: fix_T needs a temperature")
        self.nodes[name] = Node(name, p_guess, T_guess, boundary=False, fixed_T=fix_T)
        return self

    def add(self, element: el.Element) -> "Network":
        for n in (element.a, element.b):
            if n not in self.nodes:
                raise KeyError(f"element {element.name}: unknown node '{n}'")
        element.gas = self.gas
        self.elements.append(element)
        return self

    # --- solving --------------------------------------------------------
    def _internal(self):
        return [n for n in self.nodes.values() if not n.boundary]

    def _initial_guess(self):
        bnd = [n for n in self.nodes.values() if n.boundary]
        if not bnd:
            raise ValueError("network needs at least one boundary node")
        p_mean = math.exp(sum(math.log(n.p) for n in bnd) / len(bnd))
        T_mean = sum(n.T for n in bnd) / len(bnd)
        for n in self._internal():
            if n.p is None:
                n.p = p_mean
            if n.T is None:
                n.T = T_mean

    def _flows(self, p, T):
        return {e.name: e.mass_flow(p[e.a], T[e.a], p[e.b], T[e.b]) for e in self.elements}

    def _balance(self, m):
        net = {n: 0.0 for n in self.nodes}
        for e in self.elements:
            net[e.a] -= m[e.name]
            net[e.b] += m[e.name]
        return net

    def _mix_temperatures(self, m, T):
        Tnew = dict(T)
        for n in self._internal():
            if n.fixed_T:
                continue
            num = den = 0.0
            for e in self.elements:
                f = m[e.name]
                if e.b == n.name and f > 0:
                    num += f * T[e.a]; den += f
                elif e.a == n.name and f < 0:
                    num += -f * T[e.b]; den += -f
            if den > 0:
                Tnew[n.name] = num / den
        return Tnew

    def _newton(self, p, T, names, tol, max_iter=200):
        """Damped Newton on internal-node pressures with an element-wise Jacobian."""
        idx = {k: i for i, k in enumerate(names)}
        n = len(names)
        pb = [nd.p for nd in self.nodes.values() if nd.boundary]
        p_lo, p_hi = 0.2 * min(pb), 5.0 * max(pb)

        def residual(pp):
            m = self._flows(pp, T)
            net = self._balance(m)
            scale = max(1e-12, max(abs(v) for v in m.values()))
            return np.array([net[k] for k in names]), m, scale

        if n == 0:
            return p, 0.0, 0
        r, m, scale = residual(p)
        it = 0
        for it in range(1, max_iter + 1):
            if np.max(np.abs(r)) / scale < tol:
                break
            J = np.zeros((n, n))
            for e in self.elements:
                pa, pb_ = p[e.a], p[e.b]
                f0 = m[e.name]
                ha, hb = 1e-7 * pa, 1e-7 * pb_
                da = (e.mass_flow(pa + ha, T[e.a], pb_, T[e.b]) - f0) / ha
                db = (e.mass_flow(pa, T[e.a], pb_ + hb, T[e.b]) - f0) / hb
                ia, ib = idx.get(e.a), idx.get(e.b)
                if ia is not None:
                    J[ia, ia] -= da
                    if ib is not None:
                        J[ia, ib] -= db
                if ib is not None:
                    J[ib, ib] += db
                    if ia is not None:
                        J[ib, ia] += da
            try:
                dx = np.linalg.solve(J, -r)
            except np.linalg.LinAlgError:
                dx = np.linalg.lstsq(J, -r, rcond=None)[0]
            # limit the relative pressure step
            x = np.array([p[k] for k in names])
            lim = 0.3 * x
            fac = min(1.0, float(np.min(lim / np.maximum(np.abs(dx), 1e-30))))
            r0 = np.linalg.norm(r)
            for _ in range(20):
                xn = np.clip(x + fac * dx, p_lo, p_hi)
                pp = dict(p)
                for k, v in zip(names, xn):
                    pp[k] = float(v)
                rn, mn, sn = residual(pp)
                if np.linalg.norm(rn) < r0 or fac < 1e-4:
                    break
                fac *= 0.5
            p, r, m, scale = pp, rn, mn, sn
        return p, float(np.max(np.abs(r)) / scale) if n else 0.0, it

    def solve(self, tol: float = 1e-8, max_outer: int = 50) -> Result:
        self._initial_guess()
        internal = self._internal()
        names = [n.name for n in internal]
        p = {n.name: n.p for n in self.nodes.values()}
        T = {n.name: n.T for n in self.nodes.values()}
        res_max, it_total = 0.0, 0
        for it in range(1, max_outer + 1):
            p, res_max, k = self._newton(p, T, names, tol)
            it_total += k
            m = self._flows(p, T)
            Tnew = self._mix_temperatures(m, T)
            dT = max([abs(Tnew[k] - T[k]) for k in names] or [0.0])
            T = Tnew
            if dT < 1e-6:
                break
        m = self._flows(p, T)
        net = self._balance(m)
        scale = max(1e-12, max(abs(v) for v in m.values()))
        res_max = max([abs(net[k]) / scale for k in names] or [0.0])
        converged = res_max < 1e-6
        for n in internal:
            n.p, n.T = p[n.name], T[n.name]
        res = Result(converged, f"Newton, {it_total} iterations", p, T, m, res_max, it)
        res.warnings = self._checks(res)
        return res

    # --- engineering checks --------------------------------------------
    def _checks(self, r: Result) -> list[str]:
        w = []
        if not r.converged:
            w.append(f"NOT CONVERGED: max relative mass imbalance {r.residual_max:.2e}")
        for e in self.elements:
            pa, pb = r.p[e.a], r.p[e.b]
            if isinstance(e, el.Orifice) and e.is_choked(pa, pb):
                w.append(f"{e.name}: choked")
            if isinstance(e, el.Labyrinth):
                ratio = e.last_tooth_ratio(pa, pb)
                crit = self.gas.critical_pressure_ratio()
                if ratio < crit:
                    w.append(f"{e.name}: last tooth choked (p ratio ~{ratio:.3f} < {crit:.3f})"
                             + ("" if e.choke_model else "; Martin formula over-predicts, set choke_model=True"))
                elif ratio < 1.15 * crit and not e.choke_model:
                    w.append(f"{e.name}: last tooth near choke (p ratio ~{ratio:.3f}); "
                             "Martin formula may over-predict")
            if isinstance(e, el.Duct):
                up = e.a if r.m[e.name] >= 0 else e.b
                st = e.state(r.m[e.name], 0.5 * (pa + pb), r.T[up])
                if st["mach"] > 0.3:
                    w.append(f"{e.name}: Mach {st['mach']:.2f} > 0.3, consider compressible=True")
                if e.omega and not e.c_rot:
                    w.append(f"{e.name}: rotating, swirl-entry loss off (c_rot = 0); set it if the "
                             "flow enters from a non-rotating chamber")
            if isinstance(e, el.Channel):
                up = e.a if r.m[e.name] >= 0 else e.b
                M = e.mach(r.m[e.name], min(pa, pb), r.T[up])
                if M > 0.3:
                    w.append(f"{e.name}: Mach {M:.2f} > 0.3, low-Mach channel model outside validity")
        return w

    # --- input file -----------------------------------------------------
    @classmethod
    def from_dict(cls, d: dict) -> "Network":
        net = cls()
        for n in d["nodes"]:
            if n.get("boundary"):
                net.boundary(n["name"], n["p"], n["T"])
            else:
                net.cavity(n["name"], n.get("p"), n.get("T"), fix_T=n.get("fix_T", False))
        kinds = {"orifice": el.Orifice, "labyrinth": el.Labyrinth, "channel": el.Channel,
                 "duct": el.Duct, "fixed_zeta": el.FixedZeta}
        for e in d["elements"]:
            e = dict(e)
            kind = e.pop("type")
            e.pop("source", None)  # free-text reference (drawing, report, handbook diagram)
            net.add(kinds[kind](**e))
        return net

    @classmethod
    def from_json(cls, path: str) -> "Network":
        with open(path, encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))
