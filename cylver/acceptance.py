"""Acceptance test procedure with pass limits derived from the design analysis, plus an evaluator for
bench data and identification of the friction model from it.

Structure: ISO 10100:2020 (table of contents checked: module L basic leakage tests 8.2 low test
pressure, 8.3 proof/external leakage; module P piston seal leakage, optional; module F friction force,
optional, sinusoidal 10.4.1 or constant speed 10.4.2, static friction 10.6.2). The test pressures and
durations of the standard itself were NOT available to this project; where the procedure needs them it
says 'per ISO 10100' instead of inventing values.

Principle used for every limit: a test limit is the value the design analysis relied on (Wasson 2016
ch. 13 and 33: tests verify the product and validate the model). If the bench shows more friction or
leakage than the analysis assumed, the requirements that depend on it are no longer shown to be met.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import curve_fit

from . import seals
from .model import Design
from .seals import FS_SCALE, Seal


def onset_from_fit(Fc: float, dFs: float, vs: float, fv: float) -> float:
    """Stick-slip onset (m/s) of a fitted Stribeck curve: highest speed with dF/dv < 0."""
    v = np.geomspace(1e-6, 0.5, 6000)
    slope = -dFs * 2 * v / vs**2 * np.exp(-((v / vs) ** 2)) + fv
    neg = np.nonzero(slope < 0)[0]
    return 0.0 if len(neg) == 0 else float(v[min(neg[-1] + 1, len(v) - 1)])


def procedure(res: dict) -> list[dict]:
    d = res["design"]
    p_hold = res["p_hold_bar"] * 1e5
    p_proof = d.proof_factor * max(d.p_rated, d.p_relief) / 1e5
    vm_proof = next(r for r in res["static"] if r["item"].startswith("tube bore von Mises, proof"))
    rod, pis = seals.seals_of(d)
    F_s = res["breakaway_N"]
    gov = res["fatigue_gov"]
    return [
        {"step": "I-1", "module": "inspection", "what": "piston guide and gland bearing diametral clearance",
         "limit": f"<= {d.g['piston_clear_mm']:.2f} / {d.g['gland_clear_mm']:.2f} mm", "why": "sets the buckling tilt (R-07) and the extrusion gap (R-09)"},
        {"step": "L-1", "module": "ISO 10100 sec. 8.2", "what": "leakage at low test pressure, sight test", "limit": "no external leakage (test pressure per ISO 10100)", "why": "standard module"},
        {"step": "L-2", "module": "ISO 10100 sec. 8.3", "what": "proof / external leakage test", "limit": f"no external leakage, no permanent set; analysis at {p_proof:.0f} bar gives {vm_proof['stress_MPa']:.0f} MPa < Re {vm_proof['allow_MPa']:.0f} MPa",
         "why": "R-03; the proof pressure value must be taken from ISO 10100 (here 1.5 x max(rated, relief) is assumed)"},
        {"step": "P-1", "module": "ISO 10100 sec. 9", "what": "piston seal leakage at rated pressure, cap side",
         "limit": (f"<= {res['drift']['q_leak_budget_ml_min']:.2f} ml/min" if res["drift"]["q_leak_budget_ml_min"] > 0 else "no budget left: thermal contraction alone uses the R-13 allowance (open with the customer)"),
         "why": "the share of the R-13 drift budget left after thermal contraction"},
        {"step": "F-1", "module": "ISO 10100 sec. 10.6.2", "what": "static (breakaway) friction, unloaded",
         "limit": f"<= {1.3 * (rod.breakaway(0) + pis.breakaway(0)):.0f} N (1.3 x model)", "why": "the friction the speed and smoothness analyses assumed (R-01, R-12)"},
        {"step": "F-2", "module": "ISO 10100 sec. 10.4.2", "what": "friction at constant speed, both directions, 0.1-100 mm/s; fit the Stribeck model",
         "limit": f"stick-slip onset of the fitted curve <= {d.o['speed_min_smooth_mm_s']:g} mm/s (model: {res['v_stick_linear'] * 1e3:.1f} mm/s)",
         "why": "R-12; and it validates the seal model. An unloaded rig sees only the low-pressure static excess: run it with back pressure at the holding pressure if the rig allows"},
        {"step": "C-1", "module": "project test (not in ISO 10100)", "what": "run into both ends at full speed, 0 C and 40 C oil; record the cushion chamber pressure",
         "limit": f"peak <= {p_proof:.0f} bar; impact <= 0.10 m/s", "why": "R-10, R-11 (needs a transducer in the cushion chamber, not at the port)"},
        {"step": "Q-1", "module": "ISO 10771-1 Annex B (type test, once per design)", "what": "cyclic pressure test 0 -> p_test, 1e6 cycles",
         "limit": (f"p_test >= {res['test_pressure_1e6_bar']:.0f} bar on the {gov['chamber']} side, no failure" if res["test_pressure_1e6_bar"] == res["test_pressure_1e6_bar"] else "n/a: the governing location exceeds static strength at the duty peak (redesign first)"),
         "why": f"same Miner damage at {gov['loc']} as the design life (R-08)" + (f"; equal-damage pressure for 2e6 cycles: {res['test_pressure_2e6_bar']:.0f} bar" if res["test_pressure_2e6_bar"] == res["test_pressure_2e6_bar"] else "")},
    ]


def stribeck(v, Fc, dFs, vs, fv):
    return Fc + dFs * np.exp(-((v / vs) ** 2)) + fv * v


def bench_data(des: Design, seed: int = 1, scale_fs: float = 1.0, scale_vs: float = 1.0) -> dict:
    """SYNTHETIC module F data: the calibrated seal model, unloaded, 0.1-100 mm/s, with 1 % + 5 N noise.
    scale_fs / scale_vs multiply the static excess / Stribeck speed so a test can show a 'sticky' seal set."""
    rng = np.random.default_rng(seed)
    v = np.array([0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 5, 10, 30, 100]) * 1e-3
    rod, pis = seals.seals_of(des)
    FS_SCALE["k"], FS_SCALE["v"] = scale_fs, scale_vs
    try:
        F = np.array([rod.force(vi, 0.0, 1e-9) + pis.force(vi, 0.0, 1e-9) for vi in v])
        Fs = rod.breakaway(0.0) + pis.breakaway(0.0)
    finally:
        FS_SCALE["k"], FS_SCALE["v"] = 1.0, 1.0
    F_meas = F * (1 + 0.01 * rng.standard_normal(len(F))) + 5 * rng.standard_normal(len(F))
    return {"v": v, "F": F_meas, "F_true": F, "breakaway": Fs * (1 + 0.01 * rng.standard_normal()), "synthetic": True}


def evaluate(des: Design, data: dict, res: dict) -> dict:
    """Fit the Stribeck model to the bench data and judge F-1 / F-2."""
    p0 = [data["F"][-3], 0.2 * data["F"][-3], 5e-4, 100.0]
    popt, _ = curve_fit(stribeck, data["v"], data["F"], p0=p0, bounds=([0, 0, 1e-5, 0], [1e6, 1e5, 0.05, 1e6]), maxfev=20000)
    Fc, dFs, vs, fv = popt
    rod, pis = seals.seals_of(des)
    F1_lim = 1.3 * (rod.breakaway(0) + pis.breakaway(0))
    onset = onset_from_fit(Fc, dFs, vs, fv) * 1e3
    vmin = des.o["speed_min_smooth_mm_s"]
    return {
        "Fc_N": Fc,
        "dFs_N": dFs,
        "vs_mm_s": vs * 1e3,
        "fv_Ns_per_m": fv,
        "F1": {"value": data["breakaway"], "limit": F1_lim, "pass": data["breakaway"] <= F1_lim},
        "F2": {"value": onset, "limit": vmin, "pass": onset <= vmin},
        "synthetic": data.get("synthetic", False),
    }
