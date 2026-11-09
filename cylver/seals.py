"""Seal friction (Stribeck), stick-slip onset, extrusion gap limits and drift.

Friction model per seal (rod seal and piston seal), per unit circumference:
    F(v, p) = pi d [ (f0 + fp p) (1 + (rs - 1) exp(-(v/vs)^2)) sign(v) + fv v ]
i.e. a Coulomb level that grows with the sealed pressure, a static peak rs times the Coulomb level that
decays over the Stribeck velocity vs, and a viscous term (Stribeck form as used by Armstrong-Helouvry
et al. 1994, Automatica 30(7); pressure-proportional Coulomb level as in Mueller & Nau, Fluid Sealing
Technology). The parameter values below are ASSUMED orders of magnitude, not catalogue data. They are
the quantities that the ISO 10100 module F friction test (static friction and friction at test speed)
measures, so the acceptance test identifies them; see acceptance.py.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

# f0 N/mm, fp N/mm per MPa, rs = static/Coulomb, vs m/s, fv N/mm per (m/s); gap_mm = allowable radial
# extrusion gap at 16 / 25 / 40 MPa. ASSUMED values (see module docstring).
SEALS = {
    "PU_ucup": {"f0": 0.40, "fp": 0.10, "rs": 1.8, "vs": 0.012, "fv": 3.0, "gap_mm": (0.30, 0.20, 0.10), "desc": "polyurethane U-cup"},
    "PU_compact": {"f0": 0.40, "fp": 0.09, "rs": 1.7, "vs": 0.012, "fv": 3.0, "gap_mm": (0.30, 0.20, 0.10), "desc": "polyurethane compact piston seal"},
    "PTFE_step+PU": {"f0": 0.25, "fp": 0.05, "rs": 1.25, "vs": 0.006, "fv": 1.5, "gap_mm": (0.40, 0.30, 0.20), "desc": "PTFE step seal + PU U-cup tandem"},
    "PTFE_compact": {"f0": 0.20, "fp": 0.04, "rs": 1.2, "vs": 0.006, "fv": 1.2, "gap_mm": (0.40, 0.30, 0.20), "desc": "PTFE glide ring with O-ring energiser"},
}


@dataclass
class Seal:
    kind: str
    d: float  # sealed diameter (m)

    @property
    def p(self) -> dict:
        return SEALS[self.kind]

    def coulomb(self, p: float) -> float:
        s = self.p
        return math.pi * self.d * 1e3 * (s["f0"] + s["fp"] * max(p, 0.0) / 1e6)

    def force(self, v: float, p: float, v_eps: float = 1e-5) -> float:
        """Smooth friction force opposing v (tanh regularisation of sign(v) below v_eps)."""
        s = self.p
        fc = self.coulomb(p)
        stri = fc * (1 + (s["rs"] - 1) * math.exp(-((v / s["vs"]) ** 2)))
        return stri * math.tanh(v / v_eps) + math.pi * self.d * 1e3 * s["fv"] * v

    def slope(self, v: float, p: float) -> float:
        """dF/dv for v > 0 (sliding branch)."""
        s = self.p
        fc = self.coulomb(p)
        return -fc * (s["rs"] - 1) * 2 * v / s["vs"] ** 2 * math.exp(-((v / s["vs"]) ** 2)) + math.pi * self.d * 1e3 * s["fv"]

    def gap_allow(self, p: float) -> float:
        """Allowable radial extrusion gap (mm) at pressure p, linear between the table points."""
        return float(np.interp(p / 1e6, [16, 25, 40], self.p["gap_mm"]))


def seals_of(des) -> tuple[Seal, Seal]:
    return Seal(des.g["rod_seal"], des.d), Seal(des.g["piston_seal"], des.D)


def total_friction(des, v: float, pA: float, pB: float, v_eps: float = 1e-5) -> float:
    rod, pis = seals_of(des)
    return rod.force(v, pB, v_eps) + pis.force(v, abs(pA - pB), v_eps)


def stick_slip_linear(des, p_hold: float) -> float:
    """Lowest steady sliding speed at which the linearised drive (oil column + mass + friction) is stable:
    m s^2 + c s + k with c = dF/dv. Unstable (stick-slip) while c < 0. Returns the critical speed (m/s)."""
    rod, pis = seals_of(des)
    vs = np.geomspace(1e-5, 0.5, 4000)
    c = np.array([rod.slope(v, 0.0) + pis.slope(v, p_hold) for v in vs])
    neg = np.nonzero(c < 0)[0]
    return 0.0 if len(neg) == 0 else float(vs[neg[-1] + 1])


def stick_slip_sim(des, v0: float, p_hold: float, k_oil: float, t_end: float | None = None) -> dict:
    """Constant inflow A v0 into a compliant oil column (stiffness k_oil) driving the moving mass against a
    constant load p_hold * A_A. Starts in steady sliding with a 5 % speed kick, so the result shows whether
    sliding at v0 is stable (the kick dies out) or turns into stick-slip (it grows into a limit cycle).
    Returns the trace, 'growth' = late / early oscillation amplitude, and 'stick' (speed falls < 0.1 v0)."""
    m = des.mass
    A = des.A_A
    F_load = p_hold * A
    rod, pis = seals_of(des)

    def rhs(t, y):
        x, v, p = y
        F = rod.force(v, 0.0, 2e-5) + pis.force(v, p, 2e-5)
        return [v, (p * A - F_load - F) / m, k_oil / A**2 * (A * v0 - A * v)]

    period = 2 * math.pi * math.sqrt(m / k_oil)
    t_end = t_end or max(3.0, 60 * period, 0.05 / v0)
    p_ss = (F_load + rod.force(v0, 0.0) + pis.force(v0, p_hold)) / A
    sol = solve_ivp(rhs, (0, t_end), [0.0, 1.05 * v0, p_ss], method="LSODA", max_step=min(period / 40, t_end / 2000), rtol=1e-7, atol=[1e-9, 1e-8, 1.0])
    t, v = sol.t, sol.y[1]
    early = (t < 0.2 * t_end)
    late = (t > 0.8 * t_end)
    amp_e = np.ptp(v[early]) / 2
    amp_l = np.ptp(v[late]) / 2
    return {"t": t, "v": v, "growth": float(amp_l / max(amp_e, 1e-12)), "osc": float(2 * amp_l / v0), "stick": bool(np.min(v[late]) < 0.1 * v0)}


def drift(des, p_hold: float, dT_cool: float = 5.0, t: float = 600.0) -> dict:
    """Load holding at full extension with the cap chamber locked:
    - thermal: oil cooling dT contracts the locked volume, dx = V alpha dT / A_A (pressure stays at the load level)
    - the leak flow budget left for the piston seal so that total drift in t stays within the requirement."""
    V = des.A_A * des.stroke + des.g["dead_volume_cm3"] * 1e-6
    dx_T = V * des.oil.alpha_v * dT_cool / des.A_A
    allow = des.o["drift_max_mm_10min"] * 1e-3
    dT_max = allow * des.A_A / (V * des.oil.alpha_v)
    q_budget = max(0.0, allow - dx_T) * des.A_A / t
    return {
        "V_cap_L": V * 1e3,
        "dx_thermal_mm": dx_T * 1e3,
        "dT_cool_K": dT_cool,
        "dT_max_K": dT_max,
        "q_leak_budget_ml_min": q_budget * 6e7,
        "allow_mm": allow * 1e3,
    }
