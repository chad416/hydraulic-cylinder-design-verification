"""Tests: each model against an independent closed form or a published worked example, then the design
decisions the report relies on."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from cylver import acceptance, aoa, buckling, cnc_bridge, dynamics, fatigue, seals, stress
from cylver.model import build, load_order

ORD = Path(__file__).resolve().parent.parent / "data" / "orders"


def design(name="HD-160", rev=0):
    n, raw = load_order(ORD / f"{name}.json")[rev]
    return build(raw, n)


# ---------------------------------------------------------------- stress
def test_lame_thin_wall_limit_and_equilibrium():
    p, ri = 10e6, 0.1
    ro = ri * 1.001
    sr, st, sz = stress.lame(p, ri, ro, ri)
    t = ro - ri
    assert st == pytest.approx(p * (ri + ro) / 2 / t, rel=2e-3)
    assert sr == pytest.approx(-p)
    # radial equilibrium d(r sr)/dr = st on a thick tube
    ro = 0.18
    r = np.linspace(ri, ro, 2001)
    srr = np.array([stress.lame(p, ri, ro, x)[0] for x in r])
    stt = np.array([stress.lame(p, ri, ro, x)[1] for x in r])
    assert np.allclose(np.gradient(r * srr, r)[5:-5], stt[5:-5], rtol=1e-4)


def test_junction_edge_moment_and_zero_point():
    p, ri, ro, nu = 10e6, 0.1, 0.105, 0.3
    j0 = stress.junction_bending(p, ri, ro, nu, 0.0)
    rm, t = (ri + ro) / 2, ro - ri
    # thin-shell clamped edge with closed-end growth: sigma_b = 3 (1 - nu/2) p rm / (t sqrt(3 (1 - nu^2)))
    assert j0["bending"] == pytest.approx(3 * (1 - nu / 2) * p * rm / (t * math.sqrt(3 * (1 - nu**2))), rel=1e-9)
    xz = math.pi / (4 * j0["beta"])
    assert stress.junction_bending(p, ri, ro, nu, xz)["bending"] == pytest.approx(0.0, abs=1e-6 * j0["bending"])


def test_thread_stress_area_iso898():
    # ISO 898-1 tabulated stress area for M20 coarse (P 2.5) is 245 mm^2
    assert stress.thread_stress_area(0.020, 0.0025) * 1e6 == pytest.approx(245, abs=1)


# ---------------------------------------------------------------- buckling FE
def _uniform(L, EI, rear, front, n=81):
    x = np.linspace(0, L, n)
    return buckling.Beam(x, np.full(n - 1, EI), np.zeros(n - 1), rear, front)


@pytest.mark.parametrize("rear,front,K", [("pinned", "pinned", 1.0), ("fixed", "free", 2.0), ("fixed", "pinned", 0.6992), ("fixed", "fixed", 0.5)])
def test_fe_matches_euler(rear, front, K):
    L, EI = 3.0, 2e6
    P, _ = buckling.critical_load(_uniform(L, EI, rear, front))
    assert P == pytest.approx(math.pi**2 * EI / (K * L) ** 2, rel=2e-3)


def test_stepped_column_matches_exact_equation():
    EI1, a, EI2, b = 8e6, 1.4, 1.5e6, 1.6
    n = 151  # spacing 20 mm puts a node exactly on the step
    x = np.linspace(0, a + b, n)
    xm = 0.5 * (x[1:] + x[:-1])
    beam = buckling.Beam(x, np.where(xm < a, EI1, EI2), np.zeros(n - 1), "pinned", "pinned")
    P, _ = buckling.critical_load(beam)
    assert P == pytest.approx(buckling.two_segment_pinned(EI1, a, EI2, b), rel=2e-3)


def test_second_order_secant_and_amplification():
    L, EI, e = 3.0, 2e6, 0.002
    beam = _uniform(L, EI, "pinned", "pinned", 121)
    Pcr = math.pi**2 * EI / L**2
    P = 0.6 * Pcr
    w, _ = buckling.second_order(beam, P, np.zeros(121), (P * e, P * e))
    k = math.sqrt(P / EI)
    assert abs(w[60]) == pytest.approx(e * (1 / math.cos(k * L / 2) - 1), rel=1e-2)
    a0 = 0.003
    w0 = a0 * np.sin(np.pi * beam.x / L)
    w, _ = buckling.second_order(beam, P, w0)
    assert abs(w[60] + w0[60]) == pytest.approx(a0 / (1 - P / Pcr), rel=1e-2)


def test_ld80_buckling_iteration():
    a = buckling.analyse(design("LD-80", 0))
    b = buckling.analyse(design("LD-80", 1))
    assert a["sf_yield"] < 1.0  # rev A: rod 50 yields below the relief push
    assert b["sf_yield"] >= 1.5  # rev B: rod 63 + stop tube meets R-07


# ---------------------------------------------------------------- fatigue
def test_rainflow_astm_e1049_example():
    cyc = fatigue.rainflow(np.array([-2, 1, -3, 5, -1, 3, -4, 4, -2], dtype=float))
    tally: dict[float, float] = {}
    for r, _, c in cyc:
        tally[r] = tally.get(r, 0) + c
    assert tally == {3.0: 0.5, 4.0: 1.5, 6.0: 0.5, 8.0: 1.0, 9.0: 0.5}


def test_weld_sn_shape_and_miner():
    sn = fatigue.SN("weld", 71e6, 1.0)
    assert sn.cycles(71e6, 0) == pytest.approx(2e6)
    sd = 71e6 * (2 / 5) ** (1 / 3)
    assert sn.cycles(sd, 0) == pytest.approx(5e6)
    assert sn.cycles(0.5 * sd, 0) == math.inf  # below the cut-off (0.549 x knee)
    y = np.tile([0.0, 100e6], 500)
    D = fatigue.miner(fatigue.rainflow(y), sn)
    assert D == pytest.approx(499.5 / sn.cycles(100e6, 50e6), rel=1e-6)


# ---------------------------------------------------------------- fluid and dynamics
def test_bulk_modulus_limits():
    assert dynamics.k_eff(100e5, 1.6e9, 0.0, 1.0, 1e21) == pytest.approx(1.6e9, rel=1e-9)
    soft = dynamics.k_eff(5e5, 1.6e9, 0.01, 1.0, 1e21)
    stiff = dynamics.k_eff(200e5, 1.6e9, 0.01, 1.0, 1e21)
    assert soft < 0.5e9 < stiff


def test_viscosity_vg46():
    d = design()
    assert d.oil.viscosity(40) / d.oil.rho * 1e6 == pytest.approx(46.0, rel=1e-6)
    assert d.oil.viscosity(0) > 8 * d.oil.viscosity(40)


def test_frequency_formula_matches_eigenvalue():
    d = design()
    pl = dynamics.Plant(d, coulomb_on=False)
    x = 0.5 * d.stroke
    y = pl.hold_state(x)
    wn, _ = dynamics.modal(dynamics.linearise(pl, y))
    assert wn == pytest.approx(dynamics.omega_formula(pl, x, y[2], y[3]), rel=0.02)


def test_loop_gain_tends_to_merritt_as_valve_gets_fast():
    # Merritt's K_v < 2 zeta omega_h assumes an ideal valve and a symmetric actuator. As the valve gets
    # fast the model's limit must converge, and land near Merritt's value (the asymmetric areas of a
    # single-rod cylinder leave a 5-10 % difference, reported in the README).
    d = design("LD-80", 2)
    ratio = []
    for tau in (0.01, 0.001, 0.0001):
        d.g["valve_tau_s"] = tau
        lp = dynamics.loop_gain_limit(dynamics.Plant(d), 0.5 * d.stroke)
        ratio.append(lp["Kv_loop_max"] / lp["merritt_2zw"])
    assert abs(ratio[2] - ratio[1]) < 0.01 * ratio[2]
    assert 0.88 < ratio[2] < 1.12


def test_cushion_trade_and_no_cushion():
    d = design("LD-80", 2)
    small = dynamics.cushion_case(d, "A", needle_scale=0.25)
    big = dynamics.cushion_case(d, "A", needle_scale=4.0)
    none = dynamics.cushion_case(d, "A", cushion_on=False)
    assert small["p_peak_bar"] > big["p_peak_bar"]
    assert small["v_impact"] < big["v_impact"] <= none["v_impact"] + 1e-9


def test_relief_does_not_protect_cushion_chamber():
    d = design("HD-160", 0)
    c = dynamics.cushion_case(d, "A", T_C=0.0)
    assert c["p_peak_bar"] > d.p_relief / 1e5


# ---------------------------------------------------------------- seals
def test_stick_slip_simulation_agrees_with_linear_criterion():
    d = design()
    p = d.o["load_extend_kN"] * 1e3 / d.A_A
    vc = seals.stick_slip_linear(d, p)
    pl = dynamics.Plant(d)
    y = pl.hold_state(0.5 * d.stroke)
    k = d.A_A**2 / ((pl.V0 + d.A_A * 0.5 * d.stroke) / pl.K_ch(y[2]) + pl.V_line / pl.K_ln(y[2]))
    assert seals.stick_slip_sim(d, 0.5 * vc, p, k)["growth"] > 1.0
    assert seals.stick_slip_sim(d, 2.0 * vc, p, k)["growth"] < 1.0


def test_thermal_drift_hand_calc():
    d = design()
    r = seals.drift(d, 1e7, dT_cool=5.0)
    V = d.A_A * d.stroke + d.g["dead_volume_cm3"] * 1e-6
    assert r["dx_thermal_mm"] == pytest.approx(V * 7e-4 * 5 / d.A_A * 1e3, rel=1e-9)


# ---------------------------------------------------------------- design decisions
def test_hd160_rev_b_fixes_fatigue_and_cushion():
    out = []
    for rev in (0, 1):
        d = design("HD-160", rev)
        du = dynamics.run_duty(d)
        fa = fatigue.assess(d, du["pA"], du["pB"])
        out.append(min(f["stress_factor"] for f in fa))
    assert out[0] < 1.0 <= out[1]


def test_acceptance_flags_a_sticky_seal_set():
    d = design("LD-80", 2)
    res = {"p_hold_bar": d.o["load_extend_kN"] * 1e3 / d.A_A / 1e5}
    good = acceptance.evaluate(d, acceptance.bench_data(d, seed=3), res)
    bad = acceptance.evaluate(d, acceptance.bench_data(d, seed=3, scale_rs=3.0), res)
    assert bad["rs_breakaway"] > good["rs_breakaway"]
    assert not bad["F2"]["pass"]


def test_aoa_sensitivity_is_reported():
    t = aoa.coating_trade()
    r = t["result"]
    assert len(r["rank"]) == 4 and abs(sum(r["win_share"]) - 1) < 1e-9


def test_cnc_rows_are_accepted_by_project_07():
    d = design("LD-80", 2)
    rows = cnc_bridge.rows(d)
    if not (cnc_bridge.DEFAULT_07 / "cncgen").exists():
        pytest.skip("project 07 not next to this one")
    res = cnc_bridge.run_07(rows, Path(".pytest_tmp_cnc"))
    assert all(r["result"] == "PASS" for r in res), res
