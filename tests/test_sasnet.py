import math

import pytest

from sasnet import AIR, Channel, Labyrinth, Network, Orifice


def orifice_hand(cd, A, p0, T0, p):
    g, R = AIR.gamma, AIR.R
    pr = max(p / p0, (2 / (g + 1)) ** (g / (g - 1)))
    F = math.sqrt(2 * g / (g - 1) * (pr ** (2 / g) - pr ** ((g + 1) / g)))
    return cd * A * p0 / math.sqrt(R * T0) * F


def test_choked_flow_function_value():
    # textbook choked mass flow parameter for gamma = 1.4
    assert AIR.flow_function(0.1) == pytest.approx(0.6847, abs=2e-4)


def test_single_orifice_matches_hand_calc():
    net = Network().boundary("in", 200e3, 300.0).boundary("out", 150e3, 300.0)
    net.add(Orifice("o1", "in", "out", area=1e-4, cd=0.6))
    r = net.solve()
    assert r.m["o1"] == pytest.approx(orifice_hand(0.6, 1e-4, 200e3, 300.0, 150e3), rel=1e-9)


def test_orifice_choked_flag_and_value():
    net = Network().boundary("in", 300e3, 300.0).boundary("out", 100e3, 300.0)
    o = Orifice("o1", "in", "out", area=1e-4, cd=0.8)
    net.add(o)
    r = net.solve()
    assert any("choked" in w for w in r.warnings)
    expected = 0.8 * 1e-4 * 300e3 / math.sqrt(AIR.R * 300.0) * 0.68473
    assert r.m["o1"] == pytest.approx(expected, rel=1e-3)


def test_reverse_flow_sign():
    o = Orifice("o", "a", "b", area=1e-4, cd=0.6)
    assert o.mass_flow(100e3, 300, 150e3, 300) < 0
    assert o.mass_flow(150e3, 300, 100e3, 300) == pytest.approx(-o.mass_flow(100e3, 300, 150e3, 300))


def test_series_orifices_conserve_mass_and_split_pressure():
    net = Network().boundary("in", 101e3, 300).boundary("out", 100e3, 300).cavity("mid")
    net.add(Orifice("o1", "in", "mid", area=1e-4, cd=0.6))
    net.add(Orifice("o2", "mid", "out", area=1e-4, cd=0.6))
    r = net.solve()
    assert r.converged
    assert r.m["o1"] == pytest.approx(r.m["o2"], rel=1e-8)
    # identical orifices, small dp: intermediate pressure near the midpoint
    assert r.p["mid"] == pytest.approx(100.5e3, rel=2e-4)


def test_labyrinth_martin_hand_calc_and_teeth_trend():
    lab = Labyrinth("l", "a", "b", clearance=0.3e-3, radius=0.1, teeth=4, pitch=4e-3, cd=0.7, choke_model=False)
    p0, pn, T0 = 300e3, 150e3, 500.0
    cs = 0.3e-3 / 4e-3
    k = 1 / math.sqrt(1 - (3 / 4) * cs / (cs + 0.02))
    A = 2 * math.pi * 0.1 * 0.3e-3
    hand = 0.7 * k * A * p0 / math.sqrt(AIR.R * T0) * math.sqrt((1 - (pn / p0) ** 2) / (4 - math.log(pn / p0)))
    assert lab.mass_flow(p0, T0, pn, T0) == pytest.approx(hand, rel=1e-12)
    more = Labyrinth("l8", "a", "b", clearance=0.3e-3, radius=0.1, teeth=8, pitch=4e-3, cd=0.7, choke_model=False)
    assert more.mass_flow(p0, T0, pn, T0) < lab.mass_flow(p0, T0, pn, T0)


def test_stepped_labyrinth_has_no_carry_over():
    lab = Labyrinth("l", "a", "b", clearance=0.3e-3, radius=0.1, teeth=4, pitch=4e-3, stepped=True)
    assert lab.carry_over == 1.0


def test_channel_friction_close_to_colebrook():
    ch = Channel("c", "a", "b", area=1e-4, d_h=0.01, length=0.1, roughness=0.0)
    # smooth pipe, Re = 1e5: Colebrook f = 0.01799
    assert ch.friction_factor(1e5) == pytest.approx(0.01799, rel=0.02)
    assert ch.friction_factor(1000) == pytest.approx(0.064)


def test_channel_pressure_drop_consistent():
    ch = Channel("c", "a", "b", area=2e-4, d_h=0.016, length=0.5, roughness=1e-6, k_local=1.5)
    p1, p2, T = 200e3, 198e3, 400.0
    m = ch.mass_flow(p1, T, p2, T)
    rho = AIR.rho(0.5 * (p1 + p2), T)
    Re = m * ch.d_h / (ch.area * AIR.mu(T))
    zeta = ch.friction_factor(Re) * ch.length / ch.d_h + ch.k_local
    dp = zeta * m ** 2 / (2 * rho * ch.area ** 2)
    assert dp == pytest.approx(p1 - p2, rel=1e-6)


def test_mixing_temperature_in_cavity():
    net = (Network().boundary("hot", 300e3, 800).boundary("cold", 300e3, 400)
           .boundary("sink", 200e3, 300).cavity("c"))
    net.add(Orifice("oh", "hot", "c", area=1e-4, cd=0.6))
    net.add(Orifice("oc", "cold", "c", area=1e-4, cd=0.6))
    net.add(Orifice("os", "c", "sink", area=2e-4, cd=0.6))
    r = net.solve()
    assert r.converged
    mh, mc = r.m["oh"], r.m["oc"]
    assert r.T["c"] == pytest.approx((mh * 800 + mc * 400) / (mh + mc), rel=1e-6)
    assert r.m["os"] == pytest.approx(mh + mc, rel=1e-8)


def test_demo_network_balances():
    import os
    path = os.path.join(os.path.dirname(__file__), "..", "examples", "demo_network.json")
    net = Network.from_json(path)
    r = net.solve()
    assert r.converged
    for node in ("front_cavity", "rim_cavity"):
        inflow = sum(r.m[e.name] for e in net.elements if e.b == node) - sum(
            r.m[e.name] for e in net.elements if e.a == node)
        assert abs(inflow) < 1e-6 * max(abs(v) for v in r.m.values())
