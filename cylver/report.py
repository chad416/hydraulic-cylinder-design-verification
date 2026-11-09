"""Calculation report per order: revision history, then every check for each revision, figures for the
first and last revision, AoA and CNC hand-off for the last revision."""

from __future__ import annotations

import json
import math
from pathlib import Path

from . import acceptance, analysis, aoa, cnc_bridge, fmea, plots, requirements
from .model import build, load_order

ROOT = Path(__file__).resolve().parent.parent


def _t(rows: list[list], head: list[str]) -> str:
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _f(x, n=2):
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and (math.isinf(x) or math.isnan(x))):
        return "inf"
    if abs(x) >= 1e5 or (0 < abs(x) < 1e-3):
        return f"{x:.3g}"
    return f"{x:.{n}f}"


def revision_block(res: dict, img: str | None) -> tuple[str, dict]:
    d = res["design"]
    o = d.o
    b = res["buckling"]
    reqs = requirements.derive(res)
    fm = fmea.table(res)
    atp = acceptance.procedure(res)
    ev = acceptance.evaluate(d, acceptance.bench_data(d), res)
    L = []
    L.append(f"### {d.variant}\n")
    if d.raw.get("variant_why"):
        L.append(f"Why this revision: {d.raw['variant_why']}\n")
    L.append(
        _t(
            [
                ["bore / rod / stroke", f"{d.D * 1e3:g} / {d.d * 1e3:g} / {d.stroke * 1e3:g} mm"],
                ["rated / relief / supply", f"{o['p_rated_bar']} / {o['p_relief_bar']} / {o['p_supply_bar']} bar"],
                ["tube (sized)", f"{d.tube.name}, OD {d.ro * 2e3:.1f} mm, wall {d.wall * 1e3:.2f} mm"],
                ["base (sized)", f"{d.t_base * 1e3:.0f} mm"],
                ["rod", f"{d.rod.name}, Re {d.rod.Re / 1e6:.0f} MPa (ruling section {d.d * 1e3:g} mm)"],
                ["pin-to-pin closed / extended", f"{d.closed_length * 1e3:.0f} / {d.extended_length * 1e3:.0f} mm"],
                ["ports", "in the base and head blocks" if d.g.get("port_location") == "block" else "through the tube wall"],
                ["girth weld from base face", f"{d.g['base_weld_offset_mm']} mm"],
                ["retract command limit", f"{d.g.get('retract_command', 1.0):g}"],
                ["seals rod / piston", f"{d.g['rod_seal']} / {d.g['piston_seal']}"],
            ],
            ["item", "value"],
        )
    )
    L.append("\n**Requirements**\n")
    L.append(_t([[r["id"], r["text"], r["method"], r["value"], r["limit"], r["status"] + (f" ({r['note']})" if r["note"] else "")] for r in reqs], ["ID", "requirement", "method", "result", "limit", "status"]))
    L.append("\n**Static stresses**\n")
    L.append(_t([[r["item"], r["category"], f"{r['stress_MPa']:.0f}", f"{r['allow_MPa']:.0f}", f"{r['util']:.2f}", r["basis"]] for r in res["static"]], ["check", "category", "MPa", "allow MPa", "util", "basis"]))
    L.append("\n**Buckling (fully extended)**\n")
    L.append(
        _t(
            [
                ["pin-to-pin length", f"{b['L_ext_m'] * 1e3:.0f} mm"],
                ["guide length (piston to gland bearing)", f"{b['guide_m'] * 1e3:.0f} mm, clearance tilt {b['kink_mrad']:.2f} mrad"],
                ["push at relief pressure", f"{b['F_op_kN']:.0f} kN"],
                ["Euler, rod only over pin-to-pin length", f"{b['P_euler_rod_kN']:.0f} kN"],
                ["elastic critical load, stepped FE", f"{b['P_cr_kN']:.0f} kN ({b['sf_elastic']:.2f} x push)"],
                ["usable load (first yield with imperfections)", f"{b['P_yield_kN']:.0f} kN ({b['sf_yield']:.2f} x push)"],
                ["max deflection at push", f"{b['w_max_at_op_mm']:.1f} mm"],
                ["reliability index (Re CoV 7 %, load CoV 10 %)", f"beta {res['beta_buckling'][0]:.2f}, Pf {res['beta_buckling'][1]:.1e}"],
            ],
            ["item", "value"],
        )
    )
    L.append("\n**Fatigue (duty-cycle simulation -> rainflow -> Miner)**\n")
    L.append(
        _t(
            [[f["loc"], f["basis"], f"{f['max_stress_MPa']:.0f}", _f(f["damage_per_duty"]), (_f(f["life_duty_cycles"], 0) if f["life_duty_cycles"] >= 1 else "< 1 (static)"), _f(f["required"]), f"{f['stress_factor']:.2f}"] for f in res["fatigue"]],
            ["location", "basis", "max stress MPa", "damage / duty", "life (duty cycles)", "required", "stress factor"],
        )
    )
    tp = res["test_pressure_1e6_bar"]
    if tp == tp:
        L.append(f"\nEqual-damage constant-amplitude test pressure for the governing location ({res['fatigue_gov']['loc']}): {tp:.0f} bar for 1e6 cycles, {res['test_pressure_2e6_bar']:.0f} bar for 2e6 cycles.\n")
    else:
        L.append(f"\nThe governing location ({res['fatigue_gov']['loc']}) exceeds static strength at the duty peak with the elastic Kt (in reality it yields locally). A life estimate and an equal-damage test pressure mean nothing until the design changes.\n")
    L.append("\n**Cushioning**\n")
    L.append(_t([[c["side"] + (" (rod end, extending)" if c["side"] == "B" else " (cap end, retracting)"), f"{c['T_C']:.0f}", f"{c['v_entry']:.3f}", f"{c['energy_kJ']:.2f}", f"{c['p_peak_bar']:.0f}", f"{c['v_impact']:.3f}"] for c in res["cushion"]], ["end", "oil C", "entry m/s", "energy kJ", "peak bar", "impact m/s"]))
    if res.get("no_cushion"):
        L.append(f"\nWithout cushions the cap-end impact speed would be {res['no_cushion']['v_impact']:.3f} m/s.\n")
    dr = res["drift"]
    lp = res["loop"]
    L.append("\n**Friction, smoothness, drift, natural frequency**\n")
    L.append(
        _t(
            [
                ["holding pressure (load / A_A)", f"{res['p_hold_bar']:.0f} bar"],
                ["breakaway friction at holding pressure (model)", f"{res['breakaway_N']:.0f} N"],
                ["stick-slip onset, linear criterion dF/dv < 0", f"{res['v_stick_linear'] * 1e3:.1f} mm/s"],
                ["stick-slip simulation", ", ".join(f"{v * 1e3:.1f} mm/s: {'stick-slip' if s['stick'] else ('growing' if s['growth'] > 1.05 else 'smooth')}" for v, s in sorted(res["stick_sim"].items()))],
                ["thermal drift, locked cap chamber, 5 K cooling", f"{dr['dx_thermal_mm']:.1f} mm (requirement holds only for < {dr['dT_max_K']:.1f} K)"],
                ["piston seal leakage budget left", f"{dr['q_leak_budget_ml_min']:.2f} ml/min"],
                ["hydraulic natural frequency at mid stroke: formula / model eigenvalue", f"{res['omega_formula'] / 2 / math.pi:.2f} / {lp['wn'] / 2 / math.pi:.2f} Hz"],
                ["damping ratio at null (critical-centre valve, no leakage)", f"{lp['zeta']:.4f}"],
                ["position-loop velocity gain limit: model / Merritt 2 zeta omega", f"{lp['Kv_loop_max']:.2f} / {lp['merritt_2zw']:.2f} 1/s"],
            ],
            ["item", "value"],
        )
    )
    L.append("\n**FMEA** (occurrence from the computed margin, rule in fmea.py)\n")
    L.append(_t([[r["mode"], r["cause"], r["effect"], r["S"], f"{r['margin']:.2f} ({r['margin_def']})", r["O"], r["D"], r["RPN"], r["verification"]] for r in fm], ["failure mode", "cause", "effect", "S", "margin", "O", "D", "RPN", "verification"]))
    L.append("\n**Acceptance test procedure (limits derived from this analysis)**\n")
    L.append(_t([[a["step"], a["module"], a["what"], a["limit"], a["why"]] for a in atp], ["step", "module", "test", "pass limit", "why"]))
    L.append(
        f"\nBench evaluation demo on SYNTHETIC data (model + 1 % + 10 N noise; not a real test): fitted Coulomb {ev['Fc_N']:.0f} N, Stribeck speed {ev['vs_mm_s']:.1f} mm/s, "
        f"static/Coulomb {ev['rs_breakaway']:.2f}. F-1 {'PASS' if ev['F1']['pass'] else 'FAIL'} ({ev['F1']['value']:.0f} <= {ev['F1']['limit']:.0f} N); "
        f"F-2 {'PASS' if ev['F2']['pass'] else 'FAIL'} ({ev['F2']['value']:.2f} <= {ev['F2']['limit']:.2f}).\n"
    )
    if not ev["F2"]["pass"]:
        L.append("F-2 fails even on data generated from the nominal (assumed) seal model: with these seals R-12 cannot be met, the same finding as R-12 above. The limit states what the seal set would have to achieve.\n")
    if img:
        for k, cap in (("buckling", "Buckling"), ("junction", "Base junction"), ("duty", "Duty cycle"), ("fatigue", "Fatigue"), ("cushion", "Cushion"), ("friction", "Friction and stick-slip"), ("frequency", "Natural frequency")):
            L.append(f"\n![{cap}](../docs/img/{img}_{k}.png)\n")
    summary = {
        "variant": d.variant,
        "tube OD": f"{d.ro * 2e3:.1f}",
        "rod": f"{d.d * 1e3:g}",
        "buckling": f"{b['sf_yield']:.2f}",
        "fatigue (gov)": f"{res['fatigue_gov']['stress_factor']:.2f} {res['fatigue_gov']['loc']}",
        "cushion peak bar": f"{max(c['p_peak_bar'] for c in res['cushion']):.0f}",
        "not met": ", ".join(r["id"] for r in reqs if r["status"] != "met") or "none",
    }
    return "\n".join(L), {"summary": summary, "reqs": reqs, "fmea": fm}


def order_report(path: Path, out: Path, img_dir: Path, fast: bool = False, figures: bool = True) -> dict:
    variants = load_order(path)
    blocks, sums, last_res = [], [], None
    oid = None
    for i, (name, raw) in enumerate(variants):
        des = build(raw, name)
        oid = des.oid
        res = analysis.run(des, fast=fast)
        tag = f"{des.oid}_{name.replace(' ', '')}"
        want_img = figures and (i == 0 or i == len(variants) - 1)
        if want_img:
            for k, fn in (("buckling", plots.fig_buckling), ("junction", plots.fig_junction), ("duty", plots.fig_duty), ("fatigue", plots.fig_fatigue), ("cushion", plots.fig_cushion), ("friction", plots.fig_friction), ("frequency", plots.fig_frequency)):
                fn(des, res, img_dir / f"{tag}_{k}.png")
        txt, info = revision_block(res, tag if want_img else None)
        blocks.append(txt)
        sums.append(info["summary"])
        last_res, last_info, last_raw, last_des = res, info, raw, des
    # AoA and FMEA figure for the final revision
    trades = [aoa.seal_trade(last_raw), aoa.coating_trade()] + ([] if fast else [aoa.tube_trade(last_raw)])
    if figures:
        plots.fig_aoa(trades, img_dir / f"{oid}_aoa.png")
        plots.fig_fmea(last_info["fmea"], img_dir / f"{oid}_fmea.png")
    # CNC hand-off
    rows = cnc_bridge.rows(last_des)
    issues = cnc_bridge.design_side_checks(rows)
    cnc = cnc_bridge.run_07(rows, out.parent / "cnc" / oid) if not fast else []
    H = [f"# Calculation report {oid}: {variants[0][1]['title']}\n", "Generated by `python -m cylver report`. Every number below is computed from `data/orders/" + path.name + "` and `data/materials.json`. Sources and their check status: docs/DOSSIER.md section 9.\n"]
    H.append("## Revision history\n")
    H.append(_t([list(s.values()) for s in sums], list(sums[0].keys())))
    H.append("\n## Analyses per revision\n")
    H += blocks
    H.append(f"\n## Trade studies (final revision {last_des.variant})\n")
    for t in trades:
        r = t["result"]
        H.append(f"\n**{t['title']}** - rank: {' > '.join(r['rank'])}. Weight changes of +-10 % that flip the winner: {', '.join(r['flips_10pct']) or 'none'}. Share of 2000 random weight sets won: " + ", ".join(f"{a['name']} {w * 100:.0f} %" for a, w in zip(t["alternatives"], r["win_share"])) + ".\n")
        H.append(_t([[a["name"]] + [f"{a['scores'][c['name']]:.1f}" for c in t["criteria"]] + [f"{tot:.2f}"] for a, tot in zip(t["alternatives"], r["total"])], ["alternative"] + [f"{c['name']} ({c['weight']})" for c in t["criteria"]] + ["total"]))
        if any("raw" in a for a in t["alternatives"]):
            H.append("\nComputed inputs: " + "; ".join(f"{a['name']}: " + ", ".join(f"{k} {_f(v)}" for k, v in a.get("raw", {}).items()) for a in t["alternatives"]) + "\n")
        H.append("\nJudgement rationale: " + " ".join(f"*{k}*: {v}." for k, v in t["judgement_rationale"].items()) + "\n")
    if figures:
        H.append(f"\n![AoA](../docs/img/{oid}_aoa.png)\n\n![FMEA](../docs/img/{oid}_fmea.png)\n")
    H.append(f"\n## Manufacturing hand-off (final revision) to the CNC part-family generator\n")
    H.append("Rows derived from the verified design (groove sizes illustrative, by size band):\n\n```json\n" + json.dumps(rows, indent=1) + "\n```\n")
    if issues:
        H.append("Design-side checks: " + "; ".join(issues) + "\n")
    if cnc:
        H.append(_t([[c["part"], c["result"], c.get("lines", ""), c.get("cycle_time_s", "")] for c in cnc], ["part", "07 generate + verify", "lines", "cycle time s"]))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(H) + "\n", encoding="utf-8")
    return {"oid": oid, "summary": sums, "final_reqs": last_info["reqs"], "cnc": cnc}
