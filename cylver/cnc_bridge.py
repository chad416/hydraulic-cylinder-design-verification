"""Manufacturing end of the chain: turn a verified design into rows for the CNC part-family generator
(github.com/chad416/cnc-part-family-generator, project 07) and, when that project is available, generate
and verify the lathe programs with it. Project 07 is imported, not copied or changed.

Groove sizes follow 07's illustrative rows by size band (seal makers step their cross-sections by
diameter range); they are NOT seal-catalogue values.
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

from .model import Design

DEFAULT_07 = Path(__file__).resolve().parents[2] / "07_cnc-part-family-generator"

# reference rows of project 07 (data/*.csv), used by size band
PISTON_REF = [  # (max bore, seal_w, seal depth, guide_w, guide depth, guide_inset, length/bore)
    (50, 8.1, 5.0, 4.0, 2.5, 2.5),
    (70, 12.5, 7.5, 4.2, 2.5, 3.0),
    (100, 16.0, 9.0, 5.6, 2.5, 4.0),
]
GLAND_REF = [  # (max rod, oring_w, oring depth, wiper_w, wiper depth, seal_w, seal depth, guide_w, guide depth)
    (30, 3.2, 2.35, 4.0, 3.0, 5.0, 4.0, 4.0, 1.5),
    (50, 4.0, 2.8, 5.3, 4.0, 6.3, 5.0, 5.6, 2.5),
    (80, 4.8, 3.5, 6.3, 5.0, 8.1, 6.0, 5.6, 2.5),
]


def _band(table, size, scale_above=1.25):
    for row in table:
        if size <= row[0]:
            return row[1:], 1.0
    return table[-1][1:], scale_above


def rows(des: Design) -> dict[str, dict]:
    D, d = des.D * 1e3, des.d * 1e3
    th = des.g["rod_thread"]
    dt, p = th["d_mm"], th["pitch_mm"]
    wp = des.g["piston_width_mm"]
    tag = f"{des.oid}-{des.variant.replace(' ', '')}"
    spigot = dt + 3 if dt + 3 < d - 2 else d - 3
    rod_end = {
        "id": f"RE-{tag}",
        "rod_d": d,
        "thread_d": dt,
        "pitch": p,
        "thread_len": round(0.85 * dt),
        "relief_w": math.ceil(2.5 * p),
        "spigot_d": spigot,
        "spigot_len": round(0.76 * wp),
        "chamfer": max(1.5, p),
        "material": "42CrMo4",
    }
    (sw, sdp, gw, gdp, gi), k = _band(PISTON_REF, D, 1.5)
    piston = {
        "id": f"PI-{tag}",
        "bore_d": D,
        "od_clear": 0.3,
        "inner_d": spigot,
        "length": wp,
        "seal_w": round(sw * k, 1),
        "seal_bottom_d": round(D - 2 * sdp * k, 2),
        "guide_w": round(gw * k, 1),
        "guide_bottom_d": round(D - 2 * gdp * k, 2),
        "guide_inset": gi,
        "chamfer": 1,
        "stock_d": D + 5,
        "material": "C45",
    }
    (ow, odp, ww, wdp, rsw, rsdp, gw2, gdp2), k2 = _band(GLAND_REF, d, 1.25)
    L = des.g["gland_length_mm"]
    gland = {
        "id": f"GL-{tag}",
        "rod_d": d,
        "id_clear": 0.2,
        "od": D - 0.05,
        "length": L,
        "oring_w": round(ow * k2, 1),
        "oring_bottom_d": round(D - 0.05 - 2 * odp * k2, 2),
        "oring_z": round(0.68 * L),
        "wiper_w": round(ww * k2, 1),
        "wiper_bottom_d": round(d + 2 * wdp * k2, 2),
        "wiper_z": 2,
        "seal_w": round(rsw * k2, 1),
        "seal_bottom_d": round(d + 2 * rsdp * k2, 2),
        "seal_z": round(0.25 * L),
        "guide_w": round(gw2 * k2, 1),
        "guide_bottom_d": round(d + 2 * gdp2 * k2, 2),
        "guide_z": round(0.5 * L),
        "chamfer": 1,
        "stock_d": D + 5,
        "material": "C45",
    }
    return {"rod_end": rod_end, "piston": piston, "gland": gland}


def design_side_checks(r: dict[str, dict]) -> list[str]:
    """Checks that belong to the design, not to the lathe: wall left in the gland and piston."""
    out = []
    g = r["gland"]
    wall = (g["od"] - g["seal_bottom_d"]) / 2
    if wall < 3.0:
        out.append(f"gland: {wall:.1f} mm wall between the rod-seal groove and the OD (< 3 mm)")
    pi = r["piston"]
    wall_p = (pi["seal_bottom_d"] - pi["inner_d"]) / 2
    if wall_p < 3.0:
        out.append(f"piston: {wall_p:.1f} mm wall between the seal groove and the rod seat (< 3 mm)")
    return out


def run_07(r: dict[str, dict], out_dir: Path, root07: Path | None = None) -> list[dict]:
    """Generate and verify one program per part (FANUC canned-cycle style) with project 07.
    Returns one result row per part; 'available': False when 07 is not found."""
    root = Path(os.environ.get("CNCGEN_ROOT", root07 or DEFAULT_07))
    if not (root / "cncgen").exists():
        return [{"part": v["id"], "available": False, "result": f"project 07 not found at {root}"} for v in r.values()]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from cncgen.interp import run  # type: ignore
    from cncgen.parts import build_part  # type: ignore
    from cncgen.planner import plan  # type: ignore
    from cncgen.post import write_program  # type: ignore
    from cncgen.timing import cycle_time  # type: ignore
    from cncgen.tooling import load_machine, load_tools  # type: ignore
    from cncgen.verify import verify  # type: ignore

    tools = load_tools(root / "data" / "tools.json")
    m = load_machine(root / "data" / "machine.json")
    out_dir.mkdir(parents=True, exist_ok=True)
    res = []
    for fam, row in r.items():
        srow = {k: str(v) for k, v in row.items()}
        try:
            part = build_part(fam, srow, tools, m)
        except ValueError as e:
            res.append({"part": row["id"], "available": True, "result": f"rejected by 07: {e}"})
            continue
        ops = plan(part, tools, m)
        text = write_program(part, ops, "fanuc_cycles", m, prog_no=2000 + len(res))
        (out_dir / f"{row['id']}.nc").write_text(text, encoding="utf-8", newline="\n")
        prog = run(text, "fanuc_cycles", tools, m)
        v = verify(part, prog, tools, m, h=0.1)
        tm = cycle_time(prog, m)
        n_err = len(v.errors)
        res.append(
            {
                "part": row["id"],
                "available": True,
                "result": "PASS" if n_err == 0 else f"{n_err} verifier errors: " + "; ".join(i.msg for i in v.errors[:3]),
                "lines": len(text.splitlines()),
                "cycle_time_s": round(tm.total_s, 1),
            }
        )
        del v
    return res
