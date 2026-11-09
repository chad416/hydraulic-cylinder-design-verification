"""Static stresses: thick-wall tube (Lame), flat base, tube-to-base junction bending, rod, rod thread,
and the seal extrusion gap. Each function returns plain floats so the report can show the numbers.

Formula sources (stated in docs/DOSSIER.md section 9, references):
- Lame thick cylinder: Budynas & Nisbett, Shigley's Mechanical Engineering Design, 'Pressurized cylinders'.
- Clamped / simply supported circular plate: Timoshenko & Woinowsky-Krieger, Theory of Plates and Shells, ch. 3.
- Edge bending of a cylinder joined to a rigid end: same book, ch. 15 (beam on elastic foundation).
- Thread stress area As = pi/4 ((d2+d3)/2)^2: ISO 898-1.
All quoted from memory, not re-checked against the books in this project; each is cross-checked by a test.
"""

from __future__ import annotations

import math

from .model import Design


def lame(p: float, ri: float, ro: float, r: float, closed: bool = True) -> tuple[float, float, float]:
    """(sigma_r, sigma_theta, sigma_z) at radius r for internal pressure p."""
    k = ri**2 / (ro**2 - ri**2)
    sr = p * k * (1 - ro**2 / r**2)
    st = p * k * (1 + ro**2 / r**2)
    sz = p * k if closed else 0.0
    return sr, st, sz


def von_mises(s1: float, s2: float, s3: float) -> float:
    return math.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s3 - s1) ** 2))


def bore_vm(p: float, ri: float, ro: float) -> float:
    return von_mises(*lame(p, ri, ro, ri))


def bore_dilation(p: float, ri: float, ro: float, E: float, nu: float, closed: bool = True) -> float:
    """Radial growth of the bore u(ri) = ri/E (s_t - nu (s_r + s_z))."""
    sr, st, sz = lame(p, ri, ro, ri, closed)
    return ri / E * (st - nu * (sr + sz))


def plate_edge_stress(p: float, a: float, t: float, clamped: bool = True, nu: float = 0.3) -> float:
    """Maximum radial bending stress of a circular plate under uniform pressure:
    clamped edge 3 p a^2/(4 t^2) (at the edge); simply supported 3 (3+nu) p a^2 / (8 t^2) (at the centre)."""
    if clamped:
        return 3 * p * a**2 / (4 * t**2)
    return 3 * (3 + nu) * p * a**2 / (8 * t**2)


def junction_bending(p: float, ri: float, ro: float, nu: float, x: float = 0.0) -> dict:
    """Longitudinal stress in the tube wall at distance x from a rigid base (closed-end cylinder).

    Thin-shell edge solution: the base stops the radial growth delta = p rm^2 (1 - nu/2)/(E t); the wall
    bends with beta^4 = 3 (1 - nu^2)/(rm^2 t^2) and M(x) = M0 e^{-bx} (cos bx - sin bx),
    M0 = p rm t (1 - nu/2) / (2 sqrt(3 (1 - nu^2))). Valid for rm/t well above ~5; heavy tubes are near
    that limit, which the report states.
    """
    t = ro - ri
    rm = (ro + ri) / 2
    beta = (3 * (1 - nu**2)) ** 0.25 / math.sqrt(rm * t)
    m0 = p * rm * t * (1 - nu / 2) / (2 * math.sqrt(3 * (1 - nu**2)))
    bx = beta * x
    m = m0 * math.exp(-bx) * (math.cos(bx) - math.sin(bx))
    s_mem = p * ri**2 / (ro**2 - ri**2)  # thick-wall axial membrane stress
    s_bend = 6 * abs(m) / t**2
    return {
        "beta": beta,
        "rm_over_t": rm / t,
        "M0": m0,
        "membrane": s_mem,
        "bending": s_bend,
        "total": s_mem + s_bend,
        "decay_length": math.pi / beta,
    }


def thread_stress_area(d: float, pitch: float) -> float:
    """ISO 898-1 stress area: d2 = d - 0.649519 P, d3 = d - 1.226869 P."""
    d2 = d - 0.649519 * pitch
    d3 = d - 1.226869 * pitch
    return math.pi / 4 * ((d2 + d3) / 2) ** 2


def static_checks(des: Design) -> list[dict]:
    """The static checks at rated, relief and proof pressure. 'util' = stress / allowable."""
    ri, ro = des.ri, des.ro
    tube, rod = des.tube, des.rod
    p_r, p_rel = des.p_rated, des.p_relief
    p_pr = des.proof_factor * max(p_r, p_rel)
    sf = des.sf_static
    rows = []

    def add(item, cat, stress, allow, basis):
        rows.append({"item": item, "category": cat, "stress_MPa": stress / 1e6, "allow_MPa": allow / 1e6, "util": stress / allow, "basis": basis})

    # Stress categories after the design-by-analysis split of ASME BPVC VIII-2 Part 5 (from memory):
    # primary membrane Pm <= Re/sf; primary bending Pb <= 1.5 Re/sf; primary + secondary (PL+Pb+Q)
    # range <= 2 Re (shakedown). Secondary stresses are self-limiting and are covered by fatigue.
    add("tube bore von Mises, relief pressure", "Pm", bore_vm(p_rel, ri, ro), tube.Re / sf, f"Lame closed end; Re/{sf:g}")
    add("tube bore von Mises, proof pressure", "proof", bore_vm(p_pr, ri, ro), tube.Re, "Lame closed end; no yield at proof")
    add("base plate centre (simply supported bound), relief", "Pb", plate_edge_stress(p_rel, ri, des.t_base, False, tube.nu), 1.5 * tube.Re / sf, "simply supported plate; 1.5 Re/sf")
    add("base plate edge (clamped bound), relief", "Pb", plate_edge_stress(p_rel, ri, des.t_base, True), 1.5 * tube.Re / sf, "clamped plate; the real base lies between the bounds")
    jb = junction_bending(p_rel, ri, ro, tube.nu)
    add("tube wall at the base junction, membrane + edge bending, 0 -> relief", "PL+Pb+Q", jb["total"], 2 * tube.Re, f"shell edge solution, rm/t = {jb['rm_over_t']:.1f}; shakedown 2 Re")
    add("rod axial stress, push at relief", "Pm", p_rel * des.A_A / des.A_r, rod.Re / sf, "F/A (column stability in buckling)")
    th = des.g["rod_thread"]
    a_s = thread_stress_area(th["d_mm"] * 1e-3, th["pitch_mm"] * 1e-3)
    f_pull = p_rel * des.A_B
    add("rod thread at piston, pull at relief (no preload credit)", "Pm", f_pull / a_s, rod.Re / sf, "F/As, ISO 898-1 stress area")
    return rows


def extrusion_gap(des: Design, p: float) -> dict:
    """Radial gap the piston seal must bridge at pressure p: half the guide clearance plus the bore dilation."""
    c = des.mm("piston_clear_mm") / 2
    u = bore_dilation(p, des.ri, des.ro, des.tube.E, des.tube.nu)
    return {"clearance_mm": c * 1e3, "dilation_mm": u * 1e3, "gap_mm": (c + u) * 1e3}
