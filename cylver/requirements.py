"""Requirements derived from the order and the standards, each with its source, verification method
(Wasson 2016 Principle 13.10: Inspection, Examination, Analysis, Demonstration, Test) and the evidence
the analysis produced. Statement form follows Wasson Principle 22.5 (subject, shall, action, outcome,
bound). Targets marked 'project target' are design choices of this project, not taken from a standard.
"""

from __future__ import annotations

from .analysis import speed_achieved


def derive(res: dict) -> list[dict]:
    d = res["design"]
    o = d.o
    b = res["buckling"]
    st = {r["item"]: r for r in res["static"]}
    gov = res["fatigue_gov"]
    cush = res["cushion"]
    p_proof = d.proof_factor * max(d.p_rated, d.p_relief) / 1e5
    cush_peak = max(c["p_peak_bar"] for c in cush)
    cush_worst = max(cush, key=lambda c: c["p_peak_bar"])
    v_imp = max(c["v_impact"] for c in cush)
    dr = res["drift"]
    vmin = o["speed_min_smooth_mm_s"]
    vcrit = res["v_stick_linear"] * 1e3
    f_min = min(min(v) for k, v in res["freq"].items() if k != "x" and k == o["air_fraction_atm"])
    R = []

    def add(rid, text, source, method, value, limit, ok, note=""):
        R.append({"id": rid, "text": text, "source": source, "method": method, "value": value, "limit": limit, "status": "met" if ok else "NOT MET", "note": note})

    v_ach = speed_achieved(res) * 1e3
    add("R-01", f"The cylinder shall extend at >= {0.95 * o['speed_max_mm_s']:.0f} mm/s against {o['load_extend_kN']:g} kN with full valve command.",
        "order: speed, load", "A (duty simulation); T (bench)", f"{v_ach:.0f} mm/s", f">= {0.95 * o['speed_max_mm_s']:.0f} mm/s", v_ach >= 0.95 * o["speed_max_mm_s"])
    for key, rid, src in (
        ("tube bore von Mises, relief pressure", "R-02", "Lame; sf from order design rules"),
        ("tube bore von Mises, proof pressure", "R-03", "ISO 10100 sec. 8.3 proof test (factor = project assumption)"),
        ("base plate centre (simply supported bound), relief", "R-04", "plate theory, primary bending"),
        ("tube wall at the base junction, membrane + edge bending, 0 -> relief", "R-05", "shell edge solution; shakedown limit"),
        ("rod thread at piston, pull at relief (no preload credit)", "R-06", "ISO 898-1 stress area"),
    ):
        r = st[key]
        add(rid, f"Stress '{key}' shall not exceed {r['allow_MPa']:.0f} MPa ({r['category']}).", src, "A", f"{r['stress_MPa']:.0f} MPa", f"<= {r['allow_MPa']:.0f} MPa", r["util"] <= 1.0)
    add("R-07", "The fully extended cylinder shall reach first yield (with clearance tilt, bows, 1 mm pin eccentricity and self weight) only above 1.5 x the push at relief pressure.",
        "ISO/TS 13725 principles; 1.5 = project target", "A", f"{b['sf_yield']:.2f} x ({b['P_yield_kN']:.0f} kN)", ">= 1.50", b["sf_yield"] >= 1.5,
        f"elastic P_cr {b['P_cr_kN']:.0f} kN; rod-only Euler {b['P_euler_rod_kN']:.0f} kN")
    add("R-08", f"Every fatigue-critical location shall survive {gov['required']:.3g} design duty cycles with the partial factors applied (stress factor >= 1).",
        "order: duty x years; EN 1993-1-9 / FKM / VDI 2230 S-N", "A; T (ISO 10771-1 Annex B type test)", f"{gov['loc']}: {gov['stress_factor']:.2f}", ">= 1.00", gov["stress_factor"] >= 1.0)
    add("R-09", "The radial extrusion gap of the piston seal at relief pressure shall not exceed the seal's allowable gap.",
        "seal maker gap chart (assumed values)", "A; I (clearances)", f"{res['gap']['gap_mm']:.3f} mm", f"<= {res['gap_allow_mm']:.2f} mm", res["gap"]["gap_mm"] <= res["gap_allow_mm"])
    add("R-10", f"The cushion chamber pressure shall not exceed the proof pressure ({p_proof:.0f} bar) in any end-of-stroke case, 0-40 C oil.",
        "derived: cushion chamber is not protected by the port relief", "A; T (cushion test)", f"{cush_peak:.0f} bar ({cush_worst['side']} end, {cush_worst['T_C']:.0f} C)", f"<= {p_proof:.0f} bar", cush_peak <= p_proof)
    add("R-11", "The piston shall reach the end stop at <= 0.10 m/s.", "project target (end-stop impact)", "A; T", f"{v_imp:.3f} m/s", "<= 0.10 m/s", v_imp <= 0.10)
    add("R-12", f"Motion shall be free of stick-slip down to {vmin:g} mm/s.", "order: smooth speed", "A; T (ISO 10100 module F)", f"onset {vcrit:.1f} mm/s", f"<= {vmin:g} mm/s", vcrit <= vmin)
    add("R-13", f"Load-holding drift shall be <= {o['drift_max_mm_10min']:g} mm in 10 min.", "order: drift", "A; T (ISO 10100 module P)",
        f"thermal {dr['dx_thermal_mm']:.1f} mm for {dr['dT_cool_K']:g} K cooling", f"<= {dr['allow_mm']:g} mm", dr["dx_thermal_mm"] <= dr["allow_mm"],
        f"the requirement holds only while the locked oil cools < {dr['dT_max_K']:.1f} K in 10 min; the order must state the thermal condition")
    add("R-14", "The lowest hydraulic natural frequency over the stroke shall be reported to the machine builder (TPM).", "Wasson ch. 21.8 TPM", "A", f"{f_min:.1f} Hz", "report", True)
    return R
