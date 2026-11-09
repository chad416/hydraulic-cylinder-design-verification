"""FMEA of the cylinder (Wasson 2016 ch. 34: RPN = Severity x Occurrence x Detection, eq. 34.41 after
Dhillon). Occurrence is not guessed: it is read from the margin the analysis computed, with the fixed
rule below. Severity and detection are judgement, stated per row.

Occurrence rule (this project's rule, applied to 'margin' = capacity / demand at the design point):
    margin >= 2.0 -> 2,  1.5-2.0 -> 3,  1.2-1.5 -> 5,  1.0-1.2 -> 7,  < 1.0 -> 9
"""

from __future__ import annotations


def occurrence(margin: float) -> int:
    if margin >= 2.0:
        return 2
    if margin >= 1.5:
        return 3
    if margin >= 1.2:
        return 5
    if margin >= 1.0:
        return 7
    return 9


def table(res: dict) -> list[dict]:
    d = res["design"]
    b = res["buckling"]
    fat = {f["loc"]: f for f in res["fatigue"]}
    p_proof = d.proof_factor * max(d.p_rated, d.p_relief) / 1e5
    cush = max(c["p_peak_bar"] for c in res["cushion"])
    vmin = d.o["speed_min_smooth_mm_s"] * 1e-3
    vcrit = max(res["v_stick_linear"], 1e-6)
    dr = res["drift"]
    rows = [
        ("rod / column buckling", "push above usable load at full extension; worn bearings increase the tilt", "bent rod, uncontrolled load drop", 10,
         b["sf_yield"], "usable load / push at relief", 6, "analysis only (R-07); clearance inspection I-1"),
        ("piston seal extrusion", "extrusion gap grows with clearance wear and bore dilation", "internal leak, drift, seal debris", 6,
         res["gap_allow_mm"] / res["gap"]["gap_mm"], "allowable gap / gap at relief", 5, "I-1 clearances; P-1 leakage"),
        ("tube-to-base weld fatigue crack", "pressure cycles x edge bending at the weld root; root defect", "external leak, possible base separation", 10,
         fat["tube-to-base weld"]["stress_factor"], "fatigue stress factor", 7, "Q-1 type test; weld root NDT"),
        ("port cross-bore fatigue crack", "pressure cycles x hole stress concentration", "external leak at the port", 7,
         min(fat["port A cross-bore"]["stress_factor"], fat["port B cross-bore"]["stress_factor"]), "fatigue stress factor", 7, "Q-1 type test"),
        ("rod thread / piston nut fatigue", "pull cycles; preload loss raises the alternating share (phi_n -> 1)", "piston separates from rod, loss of load", 10,
         fat["rod thread at piston nut"]["stress_factor"], "fatigue stress factor (with preload)", 8, "torque check at assembly; Q-1"),
        ("cushion over-pressure", "trapped cushion oil is not protected by the port relief; cold oil", "seal extrusion, fatigue damage, tube yield", 8,
         p_proof / cush, "proof pressure / cushion peak", 4, "C-1 cushion test with chamber transducer"),
        ("stick-slip / jerky motion", "static/Coulomb friction ratio too high; seal swell; low speed", "poor positioning, noise, wear", 5,
         vmin / vcrit, "smooth-speed requirement / stick-slip onset", 4, "F-1/F-2 friction test"),
        ("drift under load", "piston seal leakage; oil cooling in the locked chamber", "load creeps down", 7,
         dr["allow_mm"] / max(dr["dx_thermal_mm"], 1e-9), "drift allowance / thermal drift", 5, "P-1 leakage test (thermal part only by analysis)"),
        ("aerated oil / spongy motion", "air entrainment, poor bleeding", "low natural frequency, overshoot, cavitation erosion", 5,
         2.0, "not quantified as a margin (natural frequency reported as TPM)", 5, "commissioning bleed procedure"),
    ]
    out = []
    for mode, cause, effect, S, margin, mdef, D, ver in rows:
        O = occurrence(margin)
        out.append({"mode": mode, "cause": cause, "effect": effect, "S": S, "margin": margin, "margin_def": mdef, "O": O, "D": D, "RPN": S * O * D, "verification": ver})
    return sorted(out, key=lambda r: -r["RPN"])
