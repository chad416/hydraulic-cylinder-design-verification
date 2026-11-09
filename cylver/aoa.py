"""Analysis of alternatives (Wasson 2016 ch. 32, methodology of Fig. 32.7: criteria and weights, utility
functions, candidate scores, ranking, sensitivity analysis).

Criteria scores are either COMPUTED by this project's models (marked 'model') or JUDGEMENT with a stated
rationale (marked 'judgement'). Judgement scores are the weakest part of any trade study and are kept
visible so a reviewer can change them. Sensitivity: each weight moved by +-10 % (Wasson's definition,
p. 682) and 2000 random weight sets (Dirichlet around the base weights).
"""

from __future__ import annotations

import copy
import math

import numpy as np

from . import seals
from .model import build


def _util(x, lo, hi, higher_better=True):
    u = (x - lo) / (hi - lo)
    u = min(1.0, max(0.0, u))
    return 10 * (u if higher_better else 1 - u)


def score(trade: dict) -> dict:
    alts = trade["alternatives"]
    w = np.array([c["weight"] for c in trade["criteria"]], dtype=float)
    w = w / w.sum()
    S = np.array([[a["scores"][c["name"]] for c in trade["criteria"]] for a in alts])
    total = S @ w
    # Wasson (2016, p. 683): only VIABLE alternatives (those meeting the specification) are ranked
    viable = np.array([a.get("viable", True) for a in alts])
    Sv = np.where(viable[:, None], S, -1e9)
    order = [i for i in np.argsort(-(S @ w)) if viable[i]]
    base_best = int(order[0])
    flips = []
    for i, c in enumerate(trade["criteria"]):
        for f in (0.9, 1.1):
            w2 = w.copy()
            w2[i] *= f
            w2 /= w2.sum()
            if int(np.argmax(Sv @ w2)) != base_best:
                flips.append(f"{c['name']} x{f:g}")
    rng = np.random.default_rng(7)
    W = rng.dirichlet(w * 40, size=2000)
    wins = np.bincount(np.argmax(Sv @ W.T, axis=0), minlength=len(alts)) / 2000
    rank = [alts[i]["name"] for i in order] + [alts[i]["name"] + " (not viable: " + alts[i].get("why_not", "") + ")" for i in range(len(alts)) if not viable[i]]
    return {"total": total, "rank": rank, "flips_10pct": flips, "win_share": wins, "S": S, "w": w}


def seal_trade(raw: dict) -> dict:
    """Rod + piston seal system for this order. Computed: stick-slip onset, friction at speed, extrusion
    margin. Judgement: dry-rod tightness, contamination tolerance, cost."""
    alts = []
    for name, rod_k, pis_k in (("PU U-cup + PU compact", "PU_ucup", "PU_compact"), ("PTFE step seal + PU tandem, PTFE piston", "PTFE_step+PU", "PTFE_compact")):
        r = copy.deepcopy(raw)
        r["design"]["rod_seal"], r["design"]["piston_seal"] = rod_k, pis_k
        d = build(r)
        p_hold = d.o["load_extend_kN"] * 1e3 / d.A_A
        vcrit = seals.stick_slip_linear(d, p_hold) * 1e3
        fr = seals.total_friction(d, d.o["speed_max_mm_s"] * 1e-3, p_hold, 0.0) / d.thrust * 100
        from .stress import extrusion_gap

        gap = extrusion_gap(d, d.p_relief)["gap_mm"]
        allow = seals.Seal(pis_k, d.D).gap_allow(d.p_relief)
        vmin = d.o["speed_min_smooth_mm_s"]
        why = []
        if allow / gap < 1.0:
            why.append(f"extrusion margin {allow / gap:.2f} < 1 (R-09)")
        if vcrit > vmin:
            why.append(f"stick-slip onset {vcrit:.1f} > {vmin:g} mm/s (R-12)")
        alts.append(
            {
                "name": name,
                "viable": not why,
                "why_not": "; ".join(why),
                "raw": {"stick-slip onset mm/s": vcrit, "friction % of thrust": fr, "extrusion margin": allow / gap},
                "scores": {
                    "smooth motion (model)": _util(math.log10(max(vcrit, 1e-3)), math.log10(0.1), math.log10(10.0), False),
                    "friction loss (model)": _util(fr, 0.5, 4.0, False),
                    "extrusion margin (model)": _util(allow / gap, 1.0, 3.0),
                    "dry rod / static tightness (judgement)": 8 if rod_k.startswith("PU") else 7,
                    "contamination tolerance (judgement)": 7 if rod_k.startswith("PU") else 6,
                    "cost (judgement)": 8 if rod_k.startswith("PU") else 5,
                },
            }
        )
    crit = [
        {"name": "smooth motion (model)", "weight": 25},
        {"name": "friction loss (model)", "weight": 10},
        {"name": "extrusion margin (model)", "weight": 20},
        {"name": "dry rod / static tightness (judgement)", "weight": 20},
        {"name": "contamination tolerance (judgement)", "weight": 15},
        {"name": "cost (judgement)", "weight": 10},
    ]
    t = {"title": "Seal system", "criteria": crit, "alternatives": alts, "judgement_rationale": {
        "dry rod / static tightness": "PU lips seal well at low pressure; PTFE step seals rely on the secondary PU seal for a dry rod",
        "contamination tolerance": "PU is the more abrasion-tolerant lip material in dirty mobile service",
        "cost": "PTFE systems need more parts (energiser, back-up) and tighter groove tolerances"}}
    t["result"] = score(t)
    return t


def coating_trade() -> dict:
    """Rod surface: all criteria are judgement except the regulatory status (checked source)."""
    crit = [
        {"name": "corrosion resistance", "weight": 25},
        {"name": "wear / impact resistance", "weight": 20},
        {"name": "regulatory risk", "weight": 20},
        {"name": "cost", "weight": 15},
        {"name": "long-rod capability", "weight": 10},
        {"name": "field repairability", "weight": 10},
    ]
    alts = [
        {"name": "hard chrome (Cr VI electroplate)", "scores": {"corrosion resistance": 6, "wear / impact resistance": 6, "regulatory risk": 2, "cost": 9, "long-rod capability": 9, "field repairability": 8}},
        {"name": "HVOF WC-CoCr thermal spray", "scores": {"corrosion resistance": 7, "wear / impact resistance": 9, "regulatory risk": 9, "cost": 4, "long-rod capability": 6, "field repairability": 5}},
        {"name": "laser cladding (Ni-base alloy)", "scores": {"corrosion resistance": 9, "wear / impact resistance": 7, "regulatory risk": 9, "cost": 3, "long-rod capability": 6, "field repairability": 4}},
        {"name": "nitrocarburise + post-oxidise", "scores": {"corrosion resistance": 6, "wear / impact resistance": 5, "regulatory risk": 9, "cost": 8, "long-rod capability": 7, "field repairability": 3}},
    ]
    t = {"title": "Rod surface", "criteria": crit, "alternatives": alts, "judgement_rationale": {
        "regulatory risk": "Cr(VI) for hard chrome needs REACH authorisation today; ECHA submitted a restriction dossier on 11 Apr 2025, consultation ran to 18 Dec 2025, RAC/SEAC opinions pending in May 2026 (Freyr, 25 May 2026). Score 2 = real chance the process must change during the product life.",
        "corrosion / wear / cost / capability / repairability": "engineering judgement from published comparisons, NOT measured here; a neutral salt-spray (ISO 9227) and an abrasion test on coupons would replace these scores"}}
    t["result"] = score(t)
    return t


def tube_trade(raw: dict) -> dict:
    """Tube material. Computed: wall thickness, mass per metre, port fatigue stress factor."""
    from . import dynamics, fatigue

    alts = []
    for mat in ("E355+SR", "25CrMo4+QT"):
        r = copy.deepcopy(raw)
        r["design"]["tube_material"] = mat
        d = build(r)
        mpm = d.tube.rho * np.pi * (d.ro**2 - d.ri**2)
        alts.append({"name": mat, "design": d, "raw": {"wall mm": d.wall * 1e3, "kg/m": mpm}})
    duty = dynamics.run_duty(alts[0]["design"])
    for a in alts:
        fa = fatigue.assess(a["design"], duty["pA"], duty["pB"])
        a["raw"]["port fatigue factor"] = min(x["stress_factor"] for x in fa if "port" in x["loc"])
        a["raw"]["weld fatigue factor"] = next(x["stress_factor"] for x in fa if "weld" in x["loc"])
        wf = a["raw"]["weld fatigue factor"]
        a["scores"] = {
            "mass (model)": _util(a["raw"]["kg/m"], 10, 80, False),
            "weld fatigue (model)": _util(wf, 0.5, 2.0),
            "port fatigue (model)": _util(a["raw"]["port fatigue factor"], 0.5, 3.0),
            "weldability (judgement)": 9 if a["name"].startswith("E355") else 5,
            "availability honed (judgement)": 9 if a["name"].startswith("E355") else 5,
            "cost (judgement)": 8 if a["name"].startswith("E355") else 5,
        }
        del a["design"]
    crit = [
        {"name": "mass (model)", "weight": 10},
        {"name": "weld fatigue (model)", "weight": 25},
        {"name": "port fatigue (model)", "weight": 15},
        {"name": "weldability (judgement)", "weight": 20},
        {"name": "availability honed (judgement)", "weight": 15},
        {"name": "cost (judgement)", "weight": 15},
    ]
    t = {"title": "Tube material", "criteria": crit, "alternatives": alts, "judgement_rationale": {
        "weldability": "CrMo grades need preheat and controlled cooling at the base and head welds",
        "availability / cost": "E355 honed tube is the stock cylinder tube; CrMo tube is a special order"}}
    t["result"] = score(t)
    return t
