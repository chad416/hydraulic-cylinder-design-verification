"""Seal friction (Stribeck), stick-slip onset, extrusion gap limits and drift.

Friction model of one seal element, per unit circumference (N/mm), at sealed pressure p (MPa):
    F(v, p) = pi d [ (f0 + fp p + fs(p) exp(-(v/vs)^2)) sign(v) + fv v ]
a Coulomb level that grows with pressure, a static excess fs that dies out over the Stribeck velocity vs,
and a viscous term (Stribeck form as in Armstrong-Helouvry et al. 1994, Automatica 30(7)). A seal system
is a list of elements; the secondary seal of a tandem rod seal sees (almost) no pressure.

CALIBRATION (checked, open access, CC BY):
- F. A. Sarbu, F. Arnaut, A. Deaconescu, T. Deaconescu, "Theoretical and experimental research concerning
  the friction forces developed in hydraulic cylinder coaxial sealing systems made from polymers",
  Polymers 16(1):157, 2024, doi:10.3390/polym16010157, Tables 3 and 5: total friction of two identical
  coaxial piston seals vs pressure (0-160 bar) at two low speeds.
- Same group, "Experimental research on polymer-based coaxial sealing systems of hydraulic cylinders for
  small displacement velocities", Polymers 14(2):290, 2022, doi:10.3390/polym14020290: seal diameter
  60 mm, width 4 mm; stick-slip for Virgin PTFE below 0.2 mm/s (eigen-vibrations to 0.3 mm/s) at pressures
  < 4 MPa; for H-PU 55D stick-slip up to 1.4 mm/s.
Derived per seal (circumference pi x 60 = 188.5 mm):
- Virgin PTFE: 29.8 N at 0 bar, slope 6.43 N/MPa (1.6 mm/s) -> f0 0.158, fp 0.034. Low-speed excess
  (0.125 vs 1.6 mm/s) 3-12 N per seal at 40-160 bar -> fs ~0.04 N/mm with vs 0.2 mm/s.
- H-PU 55D: 83.7 N at 0 bar, slope 8.95 N/MPa (1.75 mm/s) -> f0 0.444, fp 0.0475. Excess (0.7 vs
  1.75 mm/s) 5-27 N per seal -> fs 0.045 N/mm at 0 bar rising to 0.20 N/mm from 2 MPa, vs 1.0 mm/s.
EXTRUSION GAP LIMITS (checked, seal maker's catalogues, Trelleborg Sealing Solutions, Edition January 2011):
- "Hydraulic Seals - Linear, Piston Seals", Turcon Glyd Ring, Table X (installation dimensions): maximum radial
  clearance S max at 10 / 20 / 40 MPa, standard application, by bore DN: 80-132.9 mm 0.50 / 0.30 / 0.20;
  133-329.9 mm 0.60 / 0.35 / 0.25 (40-79.9 mm 0.40 / 0.25 / 0.20).
- same catalogue, Zurcon U-Cup PUA (polyurethane 93 Shore A piston U-cup), Table XLVI: S max at 5 / 10 / 20 /
  30 / 40 MPa = 0.50 / 0.40 / 0.30 / 0.20 / 0.15 mm for DN > 60 mm (0.40 / 0.30 / 0.20 / 0.15 / 0.10 for
  DN < 60 mm), to be reduced by 30 % above 80 C (not reached: oil 0-70 C).
- "Hydraulic Seals - Linear, Rod Seals", Turcon Stepseal 2K, Table XII: S max at 10 / 20 / 40 MPa, standard
  application, rod 38-199.9 mm: 0.70 / 0.40 / 0.25 (19-37.9 mm: 0.50 / 0.30 / 0.20).
The PU rod U-cup uses the PUA piston table as the closest published polyurethane values.
STILL ASSUMED: the viscous coefficients fv (only matter at high speed) and the use of PU coaxial-seal data
for a PU U-cup (closest published data found). test_seal_model_reproduces_measured_onsets checks the
calibrated model against the measured stick-slip onsets.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

PTFE = {"f0": 0.158, "fp": 0.034, "fs0": 0.040, "fs1": 0.040, "vs": 0.2e-3, "fv": 1.2}
PU = {"f0": 0.444, "fp": 0.0475, "fs0": 0.045, "fs1": 0.20, "vs": 1.0e-3, "fv": 3.0}
# elements: (data, share of the sealed pressure it sees); gap = allowable radial extrusion gap S max (mm) from the
# seal maker's installation tables (see the module docstring): pressures (MPa) and, per diameter band (upper
# limit of the band in mm), the S max values at those pressures.
_GAP_PUA = {"p": (5, 10, 20, 30, 40), "bands": ((60.0, (0.40, 0.30, 0.20, 0.15, 0.10)), (1e9, (0.50, 0.40, 0.30, 0.20, 0.15)))}
_GAP_GLYD = {"p": (10, 20, 40), "bands": ((79.9, (0.40, 0.25, 0.20)), (132.9, (0.50, 0.30, 0.20)), (329.9, (0.60, 0.35, 0.25)),
                                         (669.9, (0.60, 0.35, 0.25)))}
_GAP_STEP2K = {"p": (10, 20, 40), "bands": ((37.9, (0.50, 0.30, 0.20)), (199.9, (0.70, 0.40, 0.25)), (255.9, (0.80, 0.60, 0.35)))}
SEALS = {
    "PU_ucup": {"elements": [(PU, 1.0)], "gap": _GAP_PUA, "desc": "polyurethane U-cup (PU coaxial data as proxy)"},
    "PU_compact": {"elements": [(PU, 1.0)], "gap": _GAP_PUA, "desc": "polyurethane compact piston seal"},
    "PTFE_step+PU": {"elements": [(PTFE, 1.0), (PU, 0.0)], "gap": _GAP_STEP2K, "desc": "PTFE step seal + PU U-cup tandem"},
    "PTFE_compact": {"elements": [(PTFE, 1.0)], "gap": _GAP_GLYD, "desc": "PTFE glide ring with O-ring energiser"},
}
FS_SCALE = {"k": 1.0, "v": 1.0}  # multiply every static excess / Stribeck speed (test and what-if use only)


def _fs(e: dict, p_mpa: float) -> float:
    return FS_SCALE["k"] * (e["fs0"] + (e["fs1"] - e["fs0"]) * min(1.0, p_mpa / 2.0))


@dataclass
class Seal:
    kind: str
    d: float  # sealed diameter (m)

    @property
    def p(self) -> dict:
        return SEALS[self.kind]

    def _c(self) -> float:
        return math.pi * self.d * 1e3  # circumference in mm

    def coulomb(self, p: float) -> float:
        return self._c() * sum(e["f0"] + e["fp"] * sh * max(p, 0.0) / 1e6 for e, sh in self.p["elements"])

    def static_excess(self, p: float) -> float:
        return self._c() * sum(_fs(e, sh * max(p, 0.0) / 1e6) for e, sh in self.p["elements"])

    def breakaway(self, p: float) -> float:
        return self.coulomb(p) + self.static_excess(p)

    def force(self, v: float, p: float, v_eps: float = 1e-5) -> float:
        """Smooth friction force opposing v (tanh regularisation of sign(v) below v_eps)."""
        c = self._c()
        f = 0.0
        for e, sh in self.p["elements"]:
            pm = sh * max(p, 0.0) / 1e6
            f += c * (e["f0"] + e["fp"] * pm + _fs(e, pm) * math.exp(-((v / (e["vs"] * FS_SCALE["v"])) ** 2))) * math.tanh(v / v_eps) + c * e["fv"] * v
        return f

    def slope(self, v: float, p: float) -> float:
        """dF/dv for v > 0 (sliding branch)."""
        c = self._c()
        out = 0.0
        for e, sh in self.p["elements"]:
            pm = sh * max(p, 0.0) / 1e6
            vs = e["vs"] * FS_SCALE["v"]
            out += -c * _fs(e, pm) * 2 * v / vs**2 * math.exp(-((v / vs) ** 2)) + c * e["fv"]
        return out

    def gap_allow(self, p: float) -> float:
        """Allowable radial extrusion gap S max (mm) at pressure p (Pa) for this seal's diameter band, linear
        between the catalogue's pressure points (held constant outside them)."""
        g = self.p["gap"]
        d_mm = self.d * 1e3
        vals = next(v for dmax, v in g["bands"] if d_mm <= dmax)
        return float(np.interp(p / 1e6, g["p"], vals))


def seals_of(des) -> tuple[Seal, Seal]:
    return Seal(des.g["rod_seal"], des.d), Seal(des.g["piston_seal"], des.D)


def total_friction(des, v: float, pA: float, pB: float, v_eps: float = 1e-5) -> float:
    rod, pis = seals_of(des)
    return rod.force(v, pB, v_eps) + pis.force(v, abs(pA - pB), v_eps)


def stick_slip_linear(des, p_hold: float) -> float:
    """Lowest steady sliding speed at which the linearised drive (oil column + mass + friction) is stable:
    m s^2 + c s + k with c = dF/dv. Unstable (stick-slip) while c < 0. Returns the critical speed (m/s)."""
    rod, pis = seals_of(des)
    vs = np.geomspace(1e-6, 0.5, 6000)
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
        F = rod.force(v, 0.0, 0.02 * v0) + pis.force(v, p, 0.02 * v0)
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
