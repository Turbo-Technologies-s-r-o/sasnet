"""Checks of the first-principles models against their exact limits (docs/theory.md)."""
import math

import pytest

from sasnet import AIR, Duct, FixedZeta, Labyrinth, Network
from sasnet import physics as ph
from sasnet.axial import Face, Rotor


def test_borda_carnot_and_profile_limits():
    assert ph.expansion_loss(0.5) == pytest.approx(0.25)            # (1-F)^2
    assert ph.expansion_loss(0.0, 2.0, 4 / 3) == pytest.approx(2.0)  # laminar pipe into plenum: alpha
    a, b = ph.profile_coefficients("round", 1000)
    assert (a, b) == (2.0, pytest.approx(4 / 3))
    a, b = ph.profile_coefficients("slot", 1000)
    assert a == pytest.approx(54 / 35) and b == pytest.approx(1.2)
    a, b = ph.profile_coefficients("round", 1e5)                     # turbulent: close to uniform
    assert 1.0 < b < a < 1.1


def test_power_law_coefficients_match_numerical_integration():
    n = 7.0
    N = 20000
    rs = [(i + 0.5) / N for i in range(N)]
    u = [(1 - r) ** (1 / n) for r in rs]
    U = sum(ui * 2 * r for ui, r in zip(u, rs)) / N
    alpha = sum((ui / U) ** 3 * 2 * r for ui, r in zip(u, rs)) / N
    beta = sum((ui / U) ** 2 * 2 * r for ui, r in zip(u, rs)) / N
    a = (n + 1) ** 3 * (2 * n + 1) ** 3 / (4 * n ** 4 * (n + 3) * (2 * n + 3))
    b = (n + 1) * (2 * n + 1) ** 2 / (4 * n ** 2 * (n + 2))
    assert a == pytest.approx(alpha, rel=1e-4) and b == pytest.approx(beta, rel=1e-4)


def test_poiseuille_constants():
    assert ph.laminar_constant("round") == 64.0
    assert ph.laminar_constant("annulus") == 96.0
    assert ph.laminar_constant("slot", 1.0) == pytest.approx(56.9, rel=0.01)   # square duct


def test_prandtl_law_values():
    # Prandtl universal law (Nikuradse): lambda(1e5) = 0.0180, lambda(1e6) = 0.0116
    assert ph.turbulent_smooth(1e5) == pytest.approx(0.0180, rel=0.01)
    assert ph.turbulent_smooth(1e6) == pytest.approx(0.0116, rel=0.01)


def test_orifice_limit_and_inlet():
    cc = ph.contraction_coefficient(0.0)
    assert cc == pytest.approx(0.63)
    assert ph.orifice_loss(0.0, 0.0) == pytest.approx(1 / cc ** 2)
    assert ph.inlet_loss(0.0, 10.0) == pytest.approx(0.5)
    assert ph.inlet_loss(1.0, 10.0) == 0.0
    # blend: orifice at l/d = 0, long-channel value from reattachment on
    assert ph.short_channel_loss(0, 0, 0.0, 1.5) == pytest.approx(ph.orifice_loss(0, 0))
    assert ph.short_channel_loss(0, 0, 2.0, 1.5) == pytest.approx(1.5)


def test_kirchhoff_constant():
    assert ph.CC_KIRCHHOFF_2D == pytest.approx(0.6110, abs=1e-4)


def test_heating_and_centrifugal():
    assert ph.heating_loss(300, 300) == 0.0
    assert ph.heating_loss(300, 330) == pytest.approx(4 * 30 / 630)
    assert ph.centrifugal_rise(1.2, 1000, 0.1, 0.2) == pytest.approx(0.5 * 1.2 * 1e6 * 0.03)


def test_duct_dp_equals_zeta_head():
    ch = Duct("c", "a", "b", shape="round", d=0.01, length=0.05, count=4)
    p1, p2, T = 300e3, 290e3, 500.0
    m = ch.mass_flow(p1, T, p2, T)
    st = ch.state(m, 0.5 * (p1 + p2), T)
    rho = AIR.rho(0.5 * (p1 + p2), T)
    assert st["zeta"] * rho * st["u"] ** 2 / 2 == pytest.approx(p1 - p2, rel=1e-6)


def test_duct_shapes():
    s = Duct("s", "a", "b", shape="slot", d=0.02, h=0.005, count=3)
    assert s.area == pytest.approx(3 * 0.02 * 0.005) and s.d_h == pytest.approx(0.008)
    a = Duct("n", "a", "b", shape="annulus", d=0.3, h=0.001)
    assert a.area == pytest.approx(math.pi * 0.3 * 0.001) and a.d_h == pytest.approx(0.002)


def test_centrifugal_pumping_and_reverse():
    ch = Duct("r", "a", "b", d=0.005, length=0.05, count=10, omega=3000.0, r_a=0.05, r_b=0.15)
    assert ch.mass_flow(200e3, 600, 200e3, 600) > 0       # pumped outward at equal pressure
    inward = Duct("r2", "a", "b", d=0.005, length=0.05, count=10, omega=3000.0, r_a=0.15, r_b=0.05)
    assert inward.mass_flow(200e3, 600, 200e3, 600) < 0


def test_labyrinth_choke_model_limits_flow():
    kw = dict(clearance=0.3e-3, radius=0.1, teeth=5, pitch=4e-3, cd=0.7)
    choked = Labyrinth("l", "a", "b", **kw)
    martin = Labyrinth("l", "a", "b", choke_model=False, **kw)
    p0, T0 = 500e3, 500.0
    # small pressure ratio: both close
    assert choked.mass_flow(p0, T0, 0.95 * p0, T0) == pytest.approx(martin.mass_flow(p0, T0, 0.95 * p0, T0), rel=0.03)
    # deep expansion: the last tooth chokes and the flow saturates (monotone),
    # while Martin's formula passes a spurious maximum and falls again
    m1, m2 = choked.mass_flow(p0, T0, 0.15 * p0, T0), choked.mass_flow(p0, T0, 0.05 * p0, T0)
    assert m2 >= m1 * (1 - 1e-9) and m2 == pytest.approx(m1, rel=0.01)
    assert martin.mass_flow(p0, T0, 0.05 * p0, T0) < martin.mass_flow(p0, T0, 0.15 * p0, T0)


def test_fixed_zeta_and_fixed_T():
    e = FixedZeta("f", "a", "b", area=1e-3, zeta=2.0)
    m = e.mass_flow(200e3, 400, 199e3, 400)
    rho = AIR.rho(199.5e3, 400)
    assert 2.0 * m ** 2 / (2 * rho * 1e-6) == pytest.approx(1e3, rel=1e-9)
    net = Network().boundary("in", 300e3, 800).boundary("out", 200e3, 300).cavity("c", T_guess=500, fix_T=True)
    net.add(Duct("c1", "in", "c", d=0.01, length=0.02, count=2))
    net.add(Duct("c2", "c", "out", d=0.01, length=0.02, count=2))
    r = net.solve()
    assert r.converged and r.T["c"] == 500


def test_axial_uniform_and_swirl():
    rot = Rotor("R", blade_force=100.0)
    rot.face("front", 0.05, 0.10, p=2e5, side=+1)
    rot.face("back", 0.05, 0.10, p=1e5, side=-1)
    A = math.pi * (0.10 ** 2 - 0.05 ** 2)
    assert rot.total() == pytest.approx(100.0 + 1e5 * A)
    f = Face("s", 0.05, 0.10, p=1e5, swirl=0.5, r_ref=0.10)
    rho, om = 2.0, 1000.0
    num = sum(2 * math.pi * r * (1e5 + 0.5 * rho * (0.5 * om) ** 2 * (r * r - 0.01)) * 1e-6
              for r in [0.05 + (i + 0.5) * 1e-6 for i in range(50000)])
    assert f.force(1e5, rho, om) == pytest.approx(num, rel=1e-6)
