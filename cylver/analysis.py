"""Run every analysis for one design and collect the numbers the requirements, FMEA and report use."""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm

from . import buckling, dynamics, fatigue, seals, stress
from .model import Design, load_materials


def reliability_index(capacity: float, demand: float, cov_c: float, cov_d: float) -> tuple[float, float]:
    """Stress-strength interference with normal capacity and demand (Wasson 2016, Principle 34.21):
    beta = (muC - muD)/sqrt(sC^2 + sD^2), Pf = Phi(-beta)."""
    b = (capacity - demand) / math.sqrt((cov_c * capacity) ** 2 + (cov_d * demand) ** 2)
    return b, float(norm.cdf(-b))


def run(des: Design, fast: bool = False) -> dict:
    mats = load_materials()
    res: dict = {"design": des}
    res["static"] = stress.static_checks(des)
    res["gap"] = stress.extrusion_gap(des, des.p_relief)
    res["junction"] = stress.junction_bending(1.0, des.ri, des.ro, des.tube.nu, des.mm("base_weld_offset_mm"))
    res["buckling"] = buckling.analyse(des)

    # duty cycle -> pressure history -> fatigue
    duty = dynamics.run_duty(des)
    res["duty"] = duty
    res["fatigue"] = fatigue.assess(des, duty["pA"], duty["pB"], mats)
    gov = min(res["fatigue"], key=lambda r: r["stress_factor"])
    res["fatigue_gov"] = gov
    pch = duty["pA"] if gov["chamber"] == "A" else duty["pB"]
    res["test_pressure_1e6_bar"] = fatigue.equivalent_test_pressure(des, gov, pch, 1e6, mats) / 1e5
    res["test_pressure_2e6_bar"] = fatigue.equivalent_test_pressure(des, gov, pch, 2e6, mats) / 1e5

    # cushions: both ends at 40 C and cold oil; needle sweep at the rod end
    cush = []
    for side in ("B", "A"):
        for T in (40.0, 0.0):
            c = dynamics.cushion_case(des, side, T_C=T)
            cush.append({k: v for k, v in c.items() if k != "r"} | {"r": c["r"] if T == 40.0 else None})
    res["cushion"] = cush
    if not fast:
        sweep = []
        for s in (0.25, 0.5, 1.0, 2.0, 4.0):
            c = dynamics.cushion_case(des, "A", needle_scale=s)
            sweep.append({"needle_scale": s, "needle_mm2": s * des.g["cushion"]["needle_mm2"], "p_peak_bar": c["p_peak_bar"], "v_impact": c["v_impact"]})
        res["cushion_sweep"] = sweep
        c0 = dynamics.cushion_case(des, "A", cushion_on=False)
        res["no_cushion"] = {"p_peak_bar": c0["p_peak_bar"], "v_impact": c0["v_impact"]}

    # friction, stick-slip, drift
    p_hold = des.o["load_extend_kN"] * 1e3 / des.A_A
    res["p_hold_bar"] = p_hold / 1e5
    res["v_stick_linear"] = seals.stick_slip_linear(des, p_hold)
    pl = dynamics.Plant(des)
    x_mid = 0.5 * des.stroke
    y = pl.hold_state(x_mid)
    k_oil = (des.A_A**2) / ((pl.V0 + des.A_A * x_mid) / pl.K_ch(y[2]) + pl.V_line / pl.K_ln(y[2]))
    res["k_oil"] = k_oil
    vmin = des.o["speed_min_smooth_mm_s"] * 1e-3
    res["stick_sim"] = {v: seals.stick_slip_sim(des, v, p_hold, k_oil) for v in sorted({0.5 * vmin, vmin, 2 * vmin, max(3e-3, 2 * res["v_stick_linear"])})}
    rod, pis = seals.seals_of(des)
    res["friction_at_vmax_N"] = seals.total_friction(des, des.o["speed_max_mm_s"] * 1e-3, p_hold, 0.0)
    res["breakaway_N"] = rod.coulomb(0.0) * rod.p["rs"] + pis.coulomb(p_hold) * pis.p["rs"]
    res["gap_allow_mm"] = pis.gap_allow(des.p_relief)
    res["drift"] = seals.drift(des, p_hold)

    # natural frequency
    res["freq"] = dynamics.frequency_vs_stroke(des)
    res["loop"] = dynamics.loop_gain_limit(dynamics.Plant(des), x_mid)
    pl2 = dynamics.Plant(des)
    yh = pl2.hold_state(x_mid)
    res["omega_formula"] = dynamics.omega_formula(pl2, x_mid, yh[2], yh[3])

    # reliability indices for the static-type margins
    b = res["buckling"]
    res["beta_buckling"] = reliability_index(b["P_yield_kN"], b["F_op_kN"], 0.07, 0.10)
    return res


def speed_achieved(res: dict) -> float:
    d = res["duty"]
    ph = dict((n, (a, b)) for n, a, b in d["phases"])
    a, b = ph["extend 2"]
    m = (d["t"] > a) & (d["t"] < b)
    return float(np.max(d["v"][m]))
