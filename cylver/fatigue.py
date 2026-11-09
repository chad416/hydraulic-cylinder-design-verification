"""Fatigue under the duty-cycle pressure history: rainflow counting, S-N curves, Miner's rule.

- Rainflow: the three-point method of ASTM E1049-85 (5.4.4), with the residue counted as half cycles.
- Welds: EN 1993-1-9 S-N shape (m=3 to 5e6 cycles, m=5 to 1e8, cut-off below), detail category from
  the order file, partial factor gamma_Mf on the stress range. Mean stress is ignored for welds
  (as-welded residual stress), as that standard does.
- Machined base metal (port cross-bore): Basquin line through fW*Rm at ND with slope k, amplitude
  corrected for mean stress by Goodman, safety factor j on stress.
- Rod thread at the piston nut: endurance sigma_ASV = 0.85 (150/d + 45) MPa at 2e6 cycles, slope 3
  (VDI 2230-1, threads rolled before heat treatment), with the preloaded-joint load factor phi_n.
Sources and their check status are listed in data/materials.json and DOSSIER.md section 13.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .model import Design, load_materials
from .stress import junction_bending, lame, thread_stress_area


def reversals(y: np.ndarray, tol: float = 0.0) -> np.ndarray:
    """Peaks and valleys of a signal (plateaus collapsed, changes below tol ignored)."""
    y = np.asarray(y, dtype=float)
    out = [float(y[0])]
    direction = 0
    for v in y[1:]:
        v = float(v)
        if direction == 0:
            if abs(v - out[0]) > tol:
                direction = 1 if v > out[0] else -1
                out.append(v)
            continue
        if (v - out[-1]) * direction >= 0:
            out[-1] = v  # the current run continues: move its extreme
        elif abs(v - out[-1]) > tol:
            direction = -direction
            out.append(v)
    return np.array(out)


def rainflow(y: np.ndarray, tol: float = 0.0) -> list[tuple[float, float, float]]:
    """ASTM E1049 rainflow. Returns (range, mean, count) with count 1.0 or 0.5."""
    stack: list[float] = []
    cycles = []
    for pt in reversals(y, tol):
        stack.append(pt)
        while len(stack) >= 3:
            x = abs(stack[-1] - stack[-2])
            yv = abs(stack[-2] - stack[-3])
            if x < yv:
                break
            if len(stack) == 3:
                cycles.append((yv, 0.5 * (stack[-2] + stack[-3]), 0.5))
                stack.pop(0)
            else:
                cycles.append((yv, 0.5 * (stack[-2] + stack[-3]), 1.0))
                last = stack.pop()
                stack.pop()
                stack.pop()
                stack.append(last)
    for i in range(len(stack) - 1):
        cycles.append((abs(stack[i + 1] - stack[i]), 0.5 * (stack[i + 1] + stack[i]), 0.5))
    return cycles


@dataclass
class SN:
    kind: str  # "weld", "base", "thread"
    ref: float  # stress (Pa): weld = detail category at 2e6 (range); base/thread = endurance amplitude at ND
    gamma: float
    Rm: float = 0.0
    ND: float = 2e6
    k: float = 3.0

    def cycles(self, rng: float, mean: float) -> float:
        """Cycles to failure for one stress range (Pa) with mean (Pa)."""
        if rng <= 0:
            return math.inf
        if self.kind == "weld":
            s = self.gamma * rng
            sc = self.ref
            sd = sc * (2 / 5) ** (1 / 3)  # range at 5e6
            sl = sd * (5 / 100) ** (1 / 5)  # cut-off at 1e8
            if s >= sd:
                return 2e6 * (sc / s) ** 3
            if s >= sl:
                return 5e6 * (sd / s) ** 5
            return math.inf
        amp = self.gamma * rng / 2
        m = self.gamma * mean
        if self.Rm and amp + abs(m) >= self.Rm:
            return 1.0  # beyond static strength: not a fatigue problem any more
        if self.kind == "base" and m > 0:
            amp = amp / max(1e-6, 1 - m / self.Rm)  # Goodman
        if amp <= self.ref:
            return math.inf
        return self.ND * (self.ref / amp) ** self.k


def miner(cycles: list[tuple[float, float, float]], sn: SN, scale: float = 1.0) -> float:
    return sum(c / sn.cycles(scale * r, scale * m) for r, m, c in cycles)


def locations(des: Design, mats: dict | None = None) -> list[dict]:
    """Fatigue-critical locations: stress per unit pressure and the chamber that loads them."""
    mats = mats or load_materials()
    ff = mats["fatigue_factors"]
    ri, ro = des.ri, des.ro
    _, st, _ = lame(1.0, ri, ro, ri)
    xw = des.mm("base_weld_offset_mm")
    jb = junction_bending(1.0, ri, ro, des.tube.nu, xw)
    from .axisym import weld_stress_per_pa

    k_weld = weld_stress_per_pa(des)[0] if des.g.get("junction_model", "fe") == "fe" else jb["total"]
    th = des.g["rod_thread"]
    d_th = th["d_mm"] * 1e-3
    a_s = thread_stress_area(d_th, th["pitch_mm"] * 1e-3)
    s_asv = 0.85 * (150 / th["d_mm"] + 45) * 1e6
    cat = des.g["weld_category"] * 1e6
    base_ref = ff["fW_sigma"] * des.tube.Rm
    # ports drilled through the tube wall see Kt x the Lame hoop stress at the bore; ports in a solid
    # base / head block see roughly Kt x the pressure (thick-wall limit of the hoop factor, assumed)
    in_tube = des.g.get("port_location", "tube") == "tube"
    port_k = des.g["port_Kt"] * (st if in_tube else 1.0)
    where = "through the tube wall" if in_tube else "in the end block"
    return [
        {
            "loc": "tube-to-base weld",
            "chamber": "A",
            "k_per_Pa": k_weld,
            "sn": SN("weld", cat, ff["gamma_Mf_weld"]),
            "basis": f"structural stress at {xw * 1e3:.0f} mm from the base by axisymmetric FE ({k_weld:.2f} x p; rigid-base shell solution {jb['total']:.2f} x p), EN 1993-1-9 cat {cat / 1e6:.0f}",
        },
        {
            "loc": "port A cross-bore",
            "chamber": "A",
            "k_per_Pa": port_k,
            "sn": SN("base", base_ref, ff["j_base_metal"], des.tube.Rm, ff["base_metal_ND"], ff["base_metal_k"]),
            "basis": f"port {where}, Kt {des.g['port_Kt']:g}",
        },
        {
            "loc": "port B cross-bore",
            "chamber": "B",
            "k_per_Pa": port_k,
            "sn": SN("base", base_ref, ff["j_base_metal"], des.tube.Rm, ff["base_metal_ND"], ff["base_metal_k"]),
            "basis": f"port {where}, Kt {des.g['port_Kt']:g}",
        },
        {
            "loc": "rod thread at piston nut",
            "chamber": "B",
            "k_per_Pa": th["phi_n"] * des.A_B / a_s,
            "sn": SN("thread", s_asv, 1.0, des.rod.Rm, ff["thread_ND"], ff["thread_k"]),
            "basis": f"phi_n {th['phi_n']:g} x pull / As, M{th['d_mm']:g}x{th['pitch_mm']:g}, VDI 2230 endurance {s_asv / 1e6:.0f} MPa",
        },
    ]


def required_cycles(des: Design) -> float:
    d = des.o["duty"]
    return d["cycles_per_hour"] * d["hours_per_year"] * d["years"]


def assess(des: Design, pA: np.ndarray, pB: np.ndarray, mats: dict | None = None) -> list[dict]:
    """Damage per duty cycle at each location, life in cycles, and the stress factor that would use up
    the required life (>1 means margin)."""
    n_req = required_cycles(des)
    out = []
    for L in locations(des, mats):
        p = pA if L["chamber"] == "A" else pB
        sig = L["k_per_Pa"] * np.clip(p, 0, None)
        cyc = rainflow(sig, tol=0.002 * float(np.max(sig)))
        D = miner(cyc, L["sn"])
        life = math.inf if D == 0 else 1 / D
        # stress factor f with n_req * D(f) = 1
        lo, hi = 0.05, 50.0
        if n_req * miner(cyc, L["sn"], hi) < 1:
            f = hi
        else:
            for _ in range(60):
                mid = math.sqrt(lo * hi)
                if n_req * miner(cyc, L["sn"], mid) > 1:
                    hi = mid
                else:
                    lo = mid
            f = lo
        out.append(
            {
                "loc": L["loc"],
                "chamber": L["chamber"],
                "basis": L["basis"],
                "max_stress_MPa": float(np.max(sig)) / 1e6,
                "n_cycles_counted": len(cyc),
                "damage_per_duty": D,
                "life_duty_cycles": life,
                "required": n_req,
                "life_ratio": life / n_req,
                "stress_factor": f,
                "cycles": cyc,
            }
        )
    return out


def equivalent_test_pressure(des: Design, loc: dict, p: np.ndarray, n_test: float, mats: dict | None = None) -> float:
    """Constant-amplitude 0 -> p_test pressure which, applied n_test times (ISO 10771-1 style type test),
    does the same damage at `loc` as the required life of duty cycles."""
    L = next(x for x in locations(des, mats) if x["loc"] == loc["loc"])
    D_life = loc["damage_per_duty"] * loc["required"]
    if D_life == 0:
        return 0.0
    if any(L["sn"].cycles(r, m) <= 1.0 for r, m, _ in loc["cycles"]):
        return float("nan")  # the peak already exceeds static strength: no fatigue test can be equivalent
    lo, hi = 1e5, 2e9
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        s = L["k_per_Pa"] * mid
        d = n_test / L["sn"].cycles(s, s / 2)
        if d > D_life:
            hi = mid
        else:
            lo = mid
    return lo
