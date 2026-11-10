"""Flow elements of a secondary air system network.

Every element connects node ``a`` to node ``b`` and returns the mass flow
from a to b (negative if the flow runs b -> a) for the given node states.
Elements are adiabatic in v0.1: total temperature is carried through.

Node pressures are treated as total pressures in plenums (cavity velocity
neglected). This is the usual SAS assumption; it is stated in every report.

References for the correlations are given per element. Nothing here is
calibrated: discharge coefficients are inputs. Calibrated values belong in
private project data, not in this code.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

from .gas import PerfectGas, AIR

# Below this pressure difference fraction the flow is linearised around zero
# to keep Newton's method well behaved (sqrt(dp) has an infinite slope at 0).
_DP_LIN = 1e-4


def _signed(flow_mag, pa: float, pb: float) -> float:
    """Apply direction and small-dp linearisation.

    ``flow_mag(p_up, p_dn)`` returns the non-negative flow magnitude for
    upstream pressure p_up > downstream p_dn.
    """
    pmax = max(pa, pb)
    dp = pa - pb
    dlin = _DP_LIN * pmax
    if abs(dp) >= dlin:
        if dp > 0:
            return flow_mag(pa, pb, upstream="a")
        return -flow_mag(pb, pa, upstream="b")
    # linear blend through zero
    m_lin = flow_mag(pmax, pmax - dlin, upstream="a" if pa >= pb else "b")
    return m_lin * dp / dlin


@dataclass
class Element:
    name: str
    a: str
    b: str
    gas: PerfectGas = field(default=AIR, repr=False)

    def mass_flow(self, pa: float, Ta: float, pb: float, Tb: float) -> float:
        raise NotImplementedError

    def describe(self) -> str:
        return self.__class__.__name__


@dataclass
class Orifice(Element):
    """Compressible flow through an orifice or short hole.

    m = Cd * A * p0 / sqrt(R T0) * F(p/p0), F = isentropic flow function,
    choked at the critical pressure ratio.

    Parameters
    ----------
    area : geometric flow area, m^2
    cd : discharge coefficient (input; for rotating holes apply e.g. the
         McGreehan & Schotsch correction before passing it in - planned v0.2)
    """

    area: float = 0.0
    cd: float = 0.6

    def mass_flow(self, pa, Ta, pb, Tb):
        g = self.gas

        def mag(p_up, p_dn, upstream):
            T0 = Ta if upstream == "a" else Tb
            return self.cd * self.area * p_up / math.sqrt(g.R * T0) * g.flow_function(p_dn / p_up)

        return _signed(mag, pa, pb)

    def is_choked(self, pa, pb) -> bool:
        p_up, p_dn = max(pa, pb), min(pa, pb)
        return p_dn / p_up <= self.gas.critical_pressure_ratio()

    def describe(self):
        return f"Orifice A={self.area:.4g} m2, Cd={self.cd:.3g}"


@dataclass
class Labyrinth(Element):
    """Straight-through or stepped labyrinth seal.

    Martin's formula with a kinetic-energy carry-over factor:

        m = cd * k * A * p0 / sqrt(R T0) * sqrt((1 - (pn/p0)^2) / (N - ln(pn/p0)))

    carry-over (Hodkinson, Egli-type) for straight-through seals:

        k = 1 / sqrt(1 - ((N - 1)/N) * (c/s) / ((c/s) + 0.02))

    k = 1 for stepped seals (no carry-over).

    Limitation v0.1: Martin's formula does not model choking of the last
    tooth; near choke it over-predicts. Reports flag the last-tooth pressure
    ratio when it approaches critical.

    Parameters
    ----------
    clearance : radial clearance c, m
    radius : seal radius at the tooth tip, m (A = 2 pi r c)
    teeth : number of teeth N
    pitch : tooth pitch s, m (straight-through only)
    cd : flow coefficient per throttle (input)
    stepped : True for a stepped seal (no carry-over)
    """

    clearance: float = 0.0
    radius: float = 0.0
    teeth: int = 1
    pitch: float = 0.0
    cd: float = 0.7
    stepped: bool = False
    choke_model: bool = True

    @property
    def area(self) -> float:
        return 2.0 * math.pi * self.radius * self.clearance

    @property
    def carry_over(self) -> float:
        N = self.teeth
        if self.stepped or N <= 1 or self.pitch <= 0:
            return 1.0
        cs = self.clearance / self.pitch
        return 1.0 / math.sqrt(1.0 - ((N - 1.0) / N) * cs / (cs + 0.02))   # JET_SPREAD, theory sec. 9

    def _martin(self, p_up, p_dn, T0, n):
        g = self.gas
        r = p_dn / p_up
        val = (1.0 - r * r) / (n - math.log(r))
        return (self.cd * self.carry_over * self.area * p_up / math.sqrt(g.R * T0)
                * math.sqrt(max(val, 0.0)))

    def _series(self, p_up, p_dn, T0):
        """N-1 throttles by Martin + last throttle as an isentropic nozzle
        (captures last-tooth choking; theory.md sec. 9)."""
        g = self.gas
        N = self.teeth
        if not self.choke_model or N < 2:
            return self._martin(p_up, p_dn, T0, N)
        k_last = self.cd * self.carry_over * self.area / math.sqrt(g.R * T0)

        def last(p):
            return k_last * p * g.flow_function(p_dn / p)

        lo, hi = p_dn, p_up
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if self._martin(p_up, mid, T0, N - 1) > last(mid):
                lo = mid   # upstream part passes more: intermediate pressure must rise
            else:
                hi = mid
            if hi - lo < 1e-10 * p_up:
                break
        pm = 0.5 * (lo + hi)
        return last(pm)

    def mass_flow(self, pa, Ta, pb, Tb):
        def mag(p_up, p_dn, upstream):
            return self._series(p_up, p_dn, Ta if upstream == "a" else Tb)

        return _signed(mag, pa, pb)

    def last_tooth_ratio(self, pa, pb) -> float:
        """Approximate pressure ratio across the last throttle.

        Assumes equal mass flow through N identical throttles in series with
        an isothermal ideal-gas density change (Martin's assumption):
        p_{N-1}^2 - p_N^2 = (p0^2 - pN^2) / N (log term neglected).
        """
        p_up, p_dn = max(pa, pb), min(pa, pb)
        N = self.teeth
        p_prev_sq = p_dn ** 2 + (p_up ** 2 - p_dn ** 2) / N
        return p_dn / math.sqrt(p_prev_sq)

    def describe(self):
        kind = "stepped" if self.stepped else "straight"
        return (f"Labyrinth {kind} N={self.teeth}, c={self.clearance*1e3:.3g} mm, "
                f"r={self.radius*1e3:.4g} mm, Cd={self.cd:.3g}, k={self.carry_over:.3f}")


@dataclass
class Channel(Element):
    """Pipe or channel with wall friction and local losses.

    Low-Mach form with density at mean pressure:

        p_a - p_b = (f L/D_h + sum K) * m^2 / (2 rho_m A^2)

    f from physics.friction_factor (Poiseuille / Prandtl / Colebrook,
    docs/theory.md sec. 6). Solved iteratively for m.

    Valid for Mach < ~0.3 in the channel; reports flag higher Mach numbers.
    Local loss coefficients K come from the user; record their origin in
    the project input. Prefer Duct, which derives inlet/exit losses.

    Parameters
    ----------
    area : flow area, m^2
    d_h : hydraulic diameter, m
    length : length, m
    roughness : absolute roughness, m
    k_local : sum of local loss coefficients referred to channel velocity
    """

    area: float = 0.0
    d_h: float = 0.0
    length: float = 0.0
    roughness: float = 0.0
    k_local: float = 0.0

    def friction_factor(self, Re: float) -> float:
        from . import physics as ph
        return ph.friction_factor(Re, "round", 0.0, self.roughness / self.d_h if self.d_h > 0 else 0.0)

    def mass_flow(self, pa, Ta, pb, Tb):
        g = self.gas

        def mag(p_up, p_dn, upstream):
            T = Ta if upstream == "a" else Tb
            rho_m = g.rho(0.5 * (p_up + p_dn), T)
            mu = g.mu(T)
            dp = p_up - p_dn
            f = 0.02
            m = 0.0
            for _ in range(50):
                zeta = f * self.length / self.d_h + self.k_local
                m_new = self.area * math.sqrt(2.0 * rho_m * dp / max(zeta, 1e-12))
                Re = m_new * self.d_h / (self.area * mu)
                f_new = self.friction_factor(Re)
                if abs(m_new - m) <= 1e-10 * max(m_new, 1e-12) and abs(f_new - f) < 1e-10:
                    m = m_new
                    break
                m, f = m_new, f_new
            return m

        return _signed(mag, pa, pb)

    def mach(self, m: float, p: float, T: float) -> float:
        g = self.gas
        v = abs(m) / (g.rho(p, T) * self.area)
        return v / math.sqrt(g.gamma * g.R * T)

    def describe(self):
        return (f"Channel A={self.area:.4g} m2, Dh={self.d_h*1e3:.4g} mm, "
                f"L={self.length*1e3:.4g} mm, sumK={self.k_local:.3g}")


# ---------------------------------------------------------------------------
# v0.2 elements
# ---------------------------------------------------------------------------
from . import physics as ph  # noqa: E402


def _drive_flow(pa, pb, drive_fn, mag_fn, dlin_frac=_DP_LIN):
    """Direction and small-drive linearisation for elements with body forces.

    drive_fn(direction) -> available pressure difference (Pa) for flow
    a->b (direction=+1) or b->a (-1), including centrifugal pumping.
    mag_fn(drive, direction) -> flow magnitude for that drive.
    """
    d_ab = drive_fn(+1)
    dlin = dlin_frac * max(pa, pb)
    if d_ab >= dlin:
        return mag_fn(d_ab, +1)
    d_ba = drive_fn(-1)
    if d_ba >= dlin:
        return -mag_fn(d_ba, -1)
    m_lin = mag_fn(dlin, +1)
    return m_lin * d_ab / dlin


@dataclass
class Duct(Element):
    """Channel model from first principles (docs/theory.md).

    dp_drive = (p_up - p_dn) + 0.5 rho_m omega^2 (r_dn^2 - r_up^2)        (sec. 8)
    dp_drive = zeta * m^2 / (2 rho_m A^2)                                  (sec. 2)
    zeta = [inlet + lambda L/d_h + expansion]  (blended to the thin-orifice
           limit below the reattachment length)                           (sec. 4-6)
           + bend + k_local + swirl entry + heating                       (sec. 7-8)

    Geometry:
      shape 'round'   : d = diameter, count holes           A = n pi d^2/4, d_h = d
      shape 'slot'    : d = width b, h = height             A = n b h,      d_h = 2bh/(b+h)
      shape 'annulus' : d = mean diameter, h = gap          A = pi d h,     d_h = 2h

    f_in  : F1 = A_channel / A_previous (0 = from a plenum, 1 = no inlet loss)
    f_out : F2 = A_channel / A_next     (0 = into a plenum, 1 = no exit loss)
    k_inlet: sharp-inlet constant (Weisbach 0.5); r_edge rounds the inlet
    omega, r_a, r_b: rotation and radii of the ends (centrifugal term)
    c_rot : swirl-entry loss efficiency (theory sec. 8; empirical, default 0)
    dT    : total temperature rise along the channel, K (heating term and
            outlet temperature; adiabatic if 0)
    compressible: (1 + M^2/4) factor on the dynamic head (sec. 2)
    """

    shape: str = "round"
    d: float = 0.0
    h: float = 0.0
    length: float = 0.0
    count: int = 1
    f_in: float = 0.0
    f_out: float = 0.0
    k_inlet: float = ph.K_INLET_SHARP
    r_edge: float = 0.0
    roughness: float = 0.0
    omega: float = 0.0
    r_a: float = 0.0
    r_b: float = 0.0
    c_rot: float = 0.0
    bend_deg: float = 0.0
    bend_r: float = 0.0
    k_local: float = 0.0
    dT: float = 0.0
    compressible: bool = False

    @property
    def area(self) -> float:
        n = max(int(self.count), 1)
        if self.shape == "round":
            return n * math.pi * self.d ** 2 / 4.0
        if self.shape == "slot":
            return n * self.d * self.h
        if self.shape == "annulus":
            return math.pi * self.d * self.h
        raise ValueError(f"unknown shape {self.shape!r}")

    @property
    def d_h(self) -> float:
        if self.shape == "round":
            return self.d
        if self.shape == "slot":
            return 2.0 * self.d * self.h / (self.d + self.h)
        return 2.0 * self.h

    @property
    def aspect(self) -> float:
        if self.shape != "slot" or self.d <= 0 or self.h <= 0:
            return 0.0
        return min(self.d, self.h) / max(self.d, self.h)

    def zeta(self, Re: float, u: float, direction: int = +1, T_in: float = 300.0) -> float:
        dh = self.d_h
        ld = self.length / dh if dh > 0 else 0.0
        F1, F2 = (self.f_in, self.f_out) if direction > 0 else (self.f_out, self.f_in)
        lam = ph.friction_factor(Re, self.shape, self.aspect, self.roughness / dh if dh > 0 else 0.0)
        if ld >= 10.0:
            alpha, beta = ph.profile_coefficients(self.shape, Re)
        else:
            alpha, beta = 1.0, 1.0   # profile not developed (sec. 4)
        z_long = (ph.inlet_loss(F1, ld, self.k_inlet, self.r_edge / dh if dh > 0 else 0.0)
                  + lam * ld + ph.expansion_loss(F2, alpha, beta))
        z = ph.short_channel_loss(F1, F2, ld, z_long)
        if self.bend_deg:
            z += ph.bend_loss(self.bend_deg, self.bend_r / dh if dh > 0 else 0.0)
        if self.omega and self.c_rot and u > 0:
            rc = 0.5 * (self.r_a + self.r_b)
            z += ph.swirl_entry_loss(abs(self.omega) * rc / u, self.c_rot)
        if self.dT:
            z += ph.heating_loss(T_in, T_in + self.dT)
        return z + self.k_local

    def mass_flow(self, pa, Ta, pb, Tb):
        g = self.gas
        A = self.area

        def rho_T(direction):
            T = Ta if direction > 0 else Tb
            return g.rho(0.5 * (pa + pb), T + 0.5 * self.dT), T

        def drive(direction):
            rho, _ = rho_T(direction)
            if direction > 0:
                return pa - pb + ph.centrifugal_rise(rho, self.omega, self.r_a, self.r_b)
            return pb - pa + ph.centrifugal_rise(rho, self.omega, self.r_b, self.r_a)

        def mag(D, direction):
            rho, T = rho_T(direction)
            mu = g.mu(T + 0.5 * self.dT)
            z, m = 2.0, 0.0
            for _ in range(60):
                m_new = A * math.sqrt(2.0 * rho * D / max(z, 1e-9))
                u = m_new / (rho * A)
                Re = m_new * self.d_h / (A * mu)
                z_new = self.zeta(Re, u, direction, T)
                if self.compressible:
                    M = u / math.sqrt(g.gamma * g.R * T)
                    z_new *= 1.0 + 0.25 * M * M
                if abs(m_new - m) <= 1e-10 * max(m_new, 1e-15) and abs(z_new - z) <= 1e-10 * z_new:
                    m = m_new
                    break
                m, z = m_new, z_new
            return m

        return _drive_flow(pa, pb, drive, mag)

    def state(self, m: float, p_mean: float, T: float) -> dict:
        g = self.gas
        rho = g.rho(p_mean, T)
        u = abs(m) / (rho * self.area)
        Re = abs(m) * self.d_h / (self.area * g.mu(T))
        return dict(u=u, Re=Re, zeta=self.zeta(Re, u, 1 if m >= 0 else -1, T),
                    mach=u / math.sqrt(g.gamma * g.R * T))

    def describe(self):
        rot = f", omega={self.omega:.4g} rad/s" if self.omega else ""
        return (f"Duct {self.shape} n={self.count} d={self.d*1e3:.4g} mm"
                + (f" h={self.h*1e3:.4g} mm" if self.shape != "round" else "")
                + f" L={self.length*1e3:.4g} mm F1={self.f_in:.3g} F2={self.f_out:.3g}{rot}")


@dataclass
class FixedZeta(Element):
    """Element with a given resistance coefficient, e.g. from CFD.

    dp = zeta * m^2 / (2 rho_m A^2) (+ centrifugal term if omega given).
    Use for disk-cavity flows, blades and anything where 1D correlations are
    not informative; record the CFD case behind zeta in the input 'source'.
    """

    area: float = 0.0
    zeta: float = 1.0
    omega: float = 0.0
    r_a: float = 0.0
    r_b: float = 0.0

    def mass_flow(self, pa, Ta, pb, Tb):
        g = self.gas

        def rho(direction):
            return g.rho(0.5 * (pa + pb), Ta if direction > 0 else Tb)

        def drive(direction):
            cent = 0.5 * rho(direction) * self.omega ** 2 * (self.r_b ** 2 - self.r_a ** 2)
            return (pa - pb + cent) if direction > 0 else (pb - pa - cent)

        def mag(D, direction):
            return self.area * math.sqrt(2.0 * rho(direction) * D / max(self.zeta, 1e-12))

        return _drive_flow(pa, pb, drive, mag)

    def describe(self):
        return f"Fixed zeta={self.zeta:.4g}, A={self.area:.4g} m2"
