"""Loss models derived from first principles and classical public-domain laws.

Every function refers to a section of docs/theory.md where it is derived.
All loss coefficients are referred to the mean velocity in the channel:
    dp0 = zeta * rho * u^2 / 2.

Constants that are not derivable are named, documented with their classical
source, and can be overridden (calibration belongs to project data, not here).
"""
from __future__ import annotations

import math

lg = math.log10

# --- classical constants (theory.md, section "Constants") -----------------
CC_KIRCHHOFF_2D = math.pi / (math.pi + 2.0)   # free-streamline jet from a sharp 2D slot (Kirchhoff 1869)
K_INLET_SHARP = 0.5                            # sharp-edged flush inlet, Weisbach (1855)
RE_LAMINAR = 2300.0                            # end of laminar flow (Reynolds 1883, Schiller 1922)
RE_TURBULENT = 4000.0                          # fully turbulent from here (blend between)
L_REATTACH = 2.0                               # jet reattachment length in diameters (sec. 5)
JET_SPREAD = 0.02                              # labyrinth jet spreading parameter (sec. 9, Hodkinson 1939)


# --- 3. Contraction coefficient -------------------------------------------
def contraction_coefficient(F: float) -> float:
    """Vena-contracta coefficient for a sharp contraction with area ratio
    F = A_small / A_large (0 = from a plenum). Weisbach's empirical law
    Cc = 0.63 + 0.37 F^3 (theory.md sec. 3)."""
    F = min(max(F, 0.0), 1.0)
    return 0.63 + 0.37 * F ** 3


# --- 4. Borda-Carnot and velocity-profile coefficients --------------------
def profile_coefficients(shape: str, Re: float) -> tuple[float, float]:
    """Kinetic-energy (alpha) and momentum (beta) coefficients of a developed
    velocity profile (theory.md sec. 4).

    Laminar: exact Poiseuille profiles (pipe: 2, 4/3; slot/annulus: 54/35, 6/5).
    Turbulent: power law u ~ y^(1/n) with n = 1/sqrt(lambda) (Prandtl's rule).
    """
    flat = shape != "round"
    if Re < RE_LAMINAR:
        return (54.0 / 35.0, 6.0 / 5.0) if flat else (2.0, 4.0 / 3.0)
    n = 1.0 / math.sqrt(friction_factor(Re, shape))
    if flat:
        a = (n + 1) ** 3 / (n ** 2 * (n + 3))
        b = (n + 1) ** 2 / (n * (n + 2))
    else:
        a = (n + 1) ** 3 * (2 * n + 1) ** 3 / (4 * n ** 4 * (n + 3) * (2 * n + 3))
        b = (n + 1) * (2 * n + 1) ** 2 / (4 * n ** 2 * (n + 2))
    return a, b


def expansion_loss(F2: float, alpha: float = 1.0, beta: float = 1.0) -> float:
    """Sudden expansion from the channel into a section larger by 1/F2
    (F2 = A_channel / A_next, 0 = into a plenum). Momentum balance with a
    non-uniform inflow profile (theory.md sec. 4):
        zeta = alpha - 2 beta F2 + F2^2   (alpha = beta = 1: Borda-Carnot (1-F2)^2)."""
    return alpha - 2.0 * beta * F2 + F2 ** 2


# --- 5. Inlet ---------------------------------------------------------------
def inlet_loss(F1: float, l_d: float, k_inlet: float = K_INLET_SHARP, r_d: float = 0.0) -> float:
    """Inlet into a channel from a section larger by 1/F1 (F1 = A_channel /
    A_previous, 0 = from a plenum), long enough for the jet to reattach.

    zeta = k_inlet * (1 - F1): contraction to the vena contracta and
    Borda-Carnot re-expansion (theory.md sec. 5); k_inlet = 0.5 for a sharp
    edge (Weisbach). Rounded edge: linear reduction to 0.04 at r/d >= 0.15.
    Channels shorter than the reattachment length are handled in
    short_channel_loss().
    """
    if F1 >= 1.0:
        return 0.0
    k = k_inlet
    if r_d > 0:
        k = max(0.04, k_inlet * (1.0 - r_d / 0.15)) if r_d < 0.15 else 0.04
    return k * (1.0 - F1)


def orifice_loss(F1: float, F2: float) -> float:
    """Thin sharp orifice (l/d -> 0): the jet leaves with velocity u/Cc and its
    kinetic energy above the downstream mean is lost (theory.md sec. 5):
        zeta = (1/Cc - F2)^2,  Cc from the upstream area ratio F1."""
    cc = contraction_coefficient(F1)
    return (1.0 / cc - F2) ** 2


def short_channel_loss(F1: float, F2: float, l_d: float, long_loss: float) -> float:
    """Blend between the thin-orifice limit and the reattached long-channel
    loss over the reattachment length (theory.md sec. 5):
        zeta = long + (orifice - long) * (1 - l/(L_r d))^2  for l/d < L_r."""
    if l_d >= L_REATTACH:
        return long_loss
    phi = (1.0 - l_d / L_REATTACH) ** 2
    return long_loss + (orifice_loss(F1, F2) - long_loss) * phi


# --- 6. Friction ------------------------------------------------------------
def laminar_constant(shape: str, aspect: float = 0.0) -> float:
    """Poiseuille number lambda*Re for developed laminar flow (theory.md sec. 6):
    pipe 64 (exact), parallel plates / narrow annulus 96 (exact), rectangle of
    aspect ratio a = h/b <= 1 by the series solution (Shah & London fit)."""
    if shape == "round":
        return 64.0
    if shape == "slot" and 0.0 < aspect <= 1.0:
        a = aspect
        return 96.0 * (1 - 1.3553 * a + 1.9467 * a ** 2 - 1.7012 * a ** 3 + 0.9564 * a ** 4 - 0.2537 * a ** 5)
    return 96.0


def turbulent_smooth(Re: float) -> float:
    """Prandtl's universal law for smooth pipes, 1/sqrt(l) = 2.0 lg(Re sqrt(l)) - 0.8
    (from the logarithmic wall law with Nikuradse's constants), solved by
    fixed-point iteration (theory.md sec. 6)."""
    lam = 0.02
    for _ in range(50):
        new = (2.0 * lg(Re * math.sqrt(lam)) - 0.8) ** -2
        if abs(new - lam) < 1e-12:
            break
        lam = new
    return lam


def turbulent_rough(Re: float, rel_rough: float) -> float:
    """Colebrook's combination of the smooth and fully rough log laws:
    1/sqrt(l) = -2 lg(e/(3.7 d) + 2.51/(Re sqrt(l)))."""
    if rel_rough <= 0:
        return turbulent_smooth(Re)
    lam = 0.02
    for _ in range(50):
        new = (-2.0 * lg(rel_rough / 3.7 + 2.51 / (Re * math.sqrt(lam)))) ** -2
        if abs(new - lam) < 1e-12:
            break
        lam = new
    return lam


def friction_factor(Re: float, shape: str = "round", aspect: float = 0.0, rel_rough: float = 0.0) -> float:
    """Darcy friction factor: laminar (Poiseuille), turbulent (Prandtl/Colebrook),
    linear blend in Re between RE_LAMINAR and RE_TURBULENT. Turbulent flow in
    non-circular ducts uses the hydraulic diameter (theory.md sec. 6)."""
    Re = max(Re, 1.0)
    lam_lam = laminar_constant(shape, aspect) / Re
    if Re <= RE_LAMINAR:
        return lam_lam
    lam_t = turbulent_rough(max(Re, RE_TURBULENT), rel_rough)
    if Re >= RE_TURBULENT:
        return turbulent_rough(Re, rel_rough)
    w = (Re - RE_LAMINAR) / (RE_TURBULENT - RE_LAMINAR)
    return (1 - w) * laminar_constant(shape, aspect) / RE_LAMINAR + w * lam_t


# --- 7. Bend ----------------------------------------------------------------
def bend_loss(theta_deg: float, r_d: float) -> float:
    """Weisbach's bend formula (1855): zeta = [0.131 + 0.163 (d/R)^3.5] * theta/90,
    R = centreline bend radius, r_d = R/d (theory.md sec. 7). For R/d < 0.5
    the value at R/d = 0.5 is used."""
    rd = max(r_d, 0.5)
    return (0.131 + 0.163 * (1.0 / rd) ** 3.5) * theta_deg / 90.0


# --- 8. Heating and rotation -----------------------------------------------
def heating_loss(T_in: float, T_out: float) -> float:
    """Momentum change due to heat addition at constant area (Rayleigh line,
    low Mach): dp = rho u^2 (T_out/T_in - 1)  ->  zeta_T = 4 (T_out - T_in)/(T_out + T_in)
    referred to the mean head (theory.md sec. 8)."""
    return 4.0 * (T_out - T_in) / (T_out + T_in)


def centrifugal_rise(rho: float, omega: float, r_from: float, r_to: float) -> float:
    """Pressure rise of fluid co-rotating with a channel moving from radius
    r_from to r_to: Euler's equation in the rotating frame (theory.md sec. 8),
    dp = 0.5 rho omega^2 (r_to^2 - r_from^2)."""
    return 0.5 * rho * omega ** 2 * (r_to ** 2 - r_from ** 2)


def swirl_entry_loss(k: float, c_rot: float) -> float:
    """Entry of non-swirling fluid into a rotating channel (or the reverse):
    the relative tangential velocity omega*r is lost with efficiency c_rot,
    zeta_rot = c_rot * k^2, k = omega r / u (theory.md sec. 8). c_rot = 1 is
    the upper bound (all relative tangential kinetic energy lost); the actual
    value is empirical and must be set by the user."""
    return c_rot * k * k
