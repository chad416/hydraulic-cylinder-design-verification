"""Buckling of the whole cylinder as a stepped beam-column.

Method (follows the published principles of ISO/TS 13725:2021: the cylinder is modelled geometrically,
fully extended, with off-axis loading and its own weight; the equations of the standard itself were not
available to this project and are not reproduced):

1. Euler-Bernoulli beam finite elements along the pin-to-pin axis. Sections: tube only (base to piston),
   tube + rod in the overlap (piston to gland, both bend together between the two bearings), rod only
   (gland to rod-eye pin).
2. Elastic critical load: generalised eigenproblem K phi = P K_G phi (geometric stiffness).
3. Second-order (P-delta) response with imperfections: the tilt of the rod in the tube allowed by the
   piston and gland clearances (a kink at the overlap), rod and tube bow, the load eccentricity at the
   pins, and self weight across the axis. Usable load = the load at which the most stressed fibre reaches
   yield (rod: axial + bending; tube: bending + hoop, von Mises).
4. Reference values: Euler for the rod alone over the pin-to-pin length (a common conservative shortcut),
   Johnson parabola for intermediate slenderness, and the exact two-segment characteristic equation
   (Timoshenko & Gere, Theory of Elastic Stability, stepped bars; from memory) used as a test oracle.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh, solve_banded
from scipy.optimize import brentq

from .model import G, Design

ENDS = {"pinned": (True, False), "fixed": (True, True), "free": (False, False), "guided": (False, True)}


@dataclass
class Beam:
    x: np.ndarray  # node positions (m)
    EI: np.ndarray  # per element
    q: np.ndarray  # transverse line load per element (N/m), + in the deflection direction
    rear: str
    front: str


def _ke(ei: float, le: float) -> np.ndarray:
    l2 = le * le
    return ei / le**3 * np.array([[12, 6 * le, -12, 6 * le], [6 * le, 4 * l2, -6 * le, 2 * l2], [-12, -6 * le, 12, -6 * le], [6 * le, 2 * l2, -6 * le, 4 * l2]])


def _kg(le: float) -> np.ndarray:
    l2 = le * le
    return 1 / (30 * le) * np.array([[36, 3 * le, -36, 3 * le], [3 * le, 4 * l2, -3 * le, -l2], [-36, -3 * le, 36, -3 * le], [3 * le, -l2, -3 * le, 4 * l2]])


def assemble(b: Beam) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = len(b.x)
    K = np.zeros((2 * n, 2 * n))
    KG = np.zeros_like(K)
    f = np.zeros(2 * n)
    for e in range(n - 1):
        le = b.x[e + 1] - b.x[e]
        idx = np.arange(2 * e, 2 * e + 4)
        K[np.ix_(idx, idx)] += _ke(b.EI[e], le)
        KG[np.ix_(idx, idx)] += _kg(le)
        q = b.q[e]
        f[idx] += q * np.array([le / 2, le * le / 12, le / 2, -le * le / 12])
    fixed = []
    w0, t0 = ENDS[b.rear]
    w1, t1 = ENDS[b.front]
    if w0:
        fixed.append(0)
    if t0:
        fixed.append(1)
    if w1:
        fixed.append(2 * n - 2)
    if t1:
        fixed.append(2 * n - 1)
    free = np.setdiff1d(np.arange(2 * n), fixed)
    return K, KG, f, free


def _solve_banded(A: np.ndarray, b: np.ndarray, bw: int) -> np.ndarray:
    """Solve the (banded) beam system without a dense factorisation."""
    n = A.shape[0]
    ab = np.zeros((2 * bw + 1, n))
    for k in range(-bw, bw + 1):
        diag = np.diagonal(A, k)
        if k >= 0:
            ab[bw - k, k:] = diag
        else:
            ab[bw - k, : n + k] = diag
    return solve_banded((bw, bw), ab, b)


def critical_load(b: Beam) -> tuple[float, np.ndarray]:
    K, KG, _, free = assemble(b)
    Kf, KGf = K[np.ix_(free, free)], KG[np.ix_(free, free)]
    # largest mu of KG phi = mu K phi gives P = 1/mu (K is positive definite once supported)
    mu, vec = eigh(KGf, Kf)
    i = int(np.argmax(mu))
    phi = np.zeros(K.shape[0])
    phi[free] = vec[:, i]
    w = phi[0::2]
    return 1.0 / mu[i], w / np.max(np.abs(w))


def second_order(b: Beam, P: float, w0: np.ndarray, M_ends: tuple[float, float] = (0.0, 0.0), asm=None) -> tuple[np.ndarray, np.ndarray]:
    """Elastic deflection w (from the imperfect shape w0) and bending moment at the nodes under axial P.
    (K - P KG) w = P KG w0 + f_q + end moments."""
    K, KG, f, free = asm or assemble(b)
    u0 = np.zeros(K.shape[0])
    u0[0::2] = w0
    # slopes of the initial shape (needed by KG)
    u0[1::2] = np.gradient(w0, b.x)
    rhs = P * KG @ u0 + f
    rhs[1] += M_ends[0]
    rhs[-1] -= M_ends[1]
    A = K - P * KG
    u = np.zeros_like(u0)
    u[free] = _solve_banded(A[np.ix_(free, free)], rhs[free], 3)
    # bending moment from element curvature (Hermite cubic), averaged at shared nodes
    n = len(b.x)
    M = np.zeros(n)
    cnt = np.zeros(n)
    for e in range(n - 1):
        le = b.x[e + 1] - b.x[e]
        ue = u[2 * e : 2 * e + 4]
        for j, xi in ((0, 0.0), (1, 1.0)):
            B = np.array([(-6 + 12 * xi) / le**2, (-4 + 6 * xi) / le, (6 - 12 * xi) / le**2, (-2 + 6 * xi) / le])
            M[e + j] += b.EI[e] * B @ ue
            cnt[e + j] += 1
    return u[0::2], M / cnt


@dataclass
class Column:
    beam: Beam
    region: np.ndarray  # per node: 0 tube only, 1 overlap, 2 rod only
    w0: np.ndarray
    e: float
    EI_t: float
    EI_r: float
    lay: dict
    asm: tuple | None = None


def build_column(des: Design, n: int = 161, s: float | None = None, imperfect: bool = True) -> Column:
    lay = des.layout(s)
    L = lay["L"]
    x = np.linspace(0.0, L, n)
    E_t, E_r = des.tube.E, des.rod.E
    I_t = math.pi / 4 * (des.ro**4 - des.ri**4)
    I_r = math.pi / 64 * des.d**4
    EI_t, EI_r = E_t * I_t, E_r * I_r
    xm = 0.5 * (x[:-1] + x[1:])
    mpl = des.mass_per_length()
    EI = np.empty(n - 1)
    q = np.empty(n - 1)
    reg_e = np.empty(n - 1, dtype=int)
    for i, xc in enumerate(xm):
        if xc < lay["x_piston_c"]:
            EI[i], q[i], reg_e[i] = EI_t, mpl["tube_cap"], 0
        elif xc < lay["x_gland_out"]:
            EI[i], q[i], reg_e[i] = EI_t + EI_r, mpl["tube_ann"] + mpl["rod"], 1
        else:
            EI[i], q[i], reg_e[i] = EI_r, mpl["rod"], 2
    q = q * G * math.cos(des.angle)
    m = des.o["mounting"]
    beam = Beam(x, EI, q if imperfect else 0 * q, m["rear"], m["front"])
    region = np.concatenate([reg_e, reg_e[-1:]])
    # imperfections, all taken in the same direction (worst case)
    w0 = np.zeros(n)
    if imperfect:
        theta = (des.mm("piston_clear_mm") + des.mm("gland_clear_mm")) / (2 * lay["guide"])
        xk = 0.5 * (lay["x_piston_c"] + lay["x_gland_c"])
        if m["front"] == "free":
            w0 += np.where(x > xk, theta * (x - xk), 0.0)
        else:
            delta = theta * xk * (L - xk) / L
            w0 += np.where(x <= xk, delta * x / xk, delta * (L - x) / (L - xk))
        lr = L - lay["x_piston_c"]
        lt = lay["x_gland_out"]
        rb = des.g["rod_bow_mm_per_m"] * 1e-3 * lr
        tb = des.g["tube_bow_mm_per_m"] * 1e-3 * lt
        w0 += np.where(x >= lay["x_piston_c"], rb * np.sin(np.pi * np.clip((x - lay["x_piston_c"]) / lr, 0, 1)), 0.0)
        w0 += np.where(x <= lt, tb * np.sin(np.pi * np.clip(x / lt, 0, 1)), 0.0)
    e = m["pin_eccentricity_mm"] * 1e-3 if imperfect else 0.0
    return Column(beam, region, w0, e, EI_t, EI_r, lay, assemble(beam))


def fibre_stress(des: Design, col: Column, P: float) -> dict:
    """Max utilisation (stress / Re) along the column at axial load P, cap-side pressure P/A_A."""
    w, M = second_order(col.beam, P, col.w0, (P * col.e, P * col.e), col.asm)
    p = P / des.A_A
    W_r = math.pi / 32 * des.d**3
    W_t = math.pi / 32 * (des.ro**4 - des.ri**4) / des.ro
    share_r = col.EI_r / (col.EI_t + col.EI_r)
    reg = col.region
    s_rod = np.where(reg == 2, P / des.A_r + np.abs(M) / W_r, np.where(reg == 1, P / des.A_r + np.abs(M) * share_r / W_r, 0.0))
    # tube: bending + hoop at the outer fibre (cap chamber pressure only acts where reg==0)
    hoop_o = 2 * p * des.ri**2 / (des.ro**2 - des.ri**2)
    sb_t = np.where(reg == 0, np.abs(M) / W_t, np.where(reg == 1, np.abs(M) * (1 - share_r) / W_t, 0.0))
    hoop = np.where(reg == 0, hoop_o, 0.0)
    s_tube = np.sqrt(sb_t**2 + hoop**2 - sb_t * hoop)
    u_rod = s_rod / des.rod.Re
    u_tube = s_tube / des.tube.Re
    i_r, i_t = int(np.argmax(u_rod)), int(np.argmax(u_tube))
    return {
        "w": w,
        "M": M,
        "w_max_mm": float(np.max(np.abs(w + col.w0))) * 1e3,
        "util_rod": float(u_rod[i_r]),
        "x_rod": float(col.beam.x[i_r]),
        "util_tube": float(u_tube[i_t]),
        "x_tube": float(col.beam.x[i_t]),
        "util": float(max(u_rod[i_r], u_tube[i_t])),
    }


def euler_rod_only(des: Design) -> float:
    L = des.extended_length
    k = {("pinned", "pinned"): 1.0, ("fixed", "pinned"): 0.699, ("fixed", "free"): 2.0, ("fixed", "guided"): 1.0, ("fixed", "fixed"): 0.5}
    m = des.o["mounting"]
    K = k.get((m["rear"], m["front"]), 1.0)
    I_r = math.pi / 64 * des.d**4
    return math.pi**2 * des.rod.E * I_r / (K * L) ** 2


def johnson_rod_only(des: Design) -> float:
    """Johnson parabola for the rod alone (same effective length); equals Euler above the transition slenderness."""
    pe = euler_rod_only(des)
    A = des.A_r
    s_e = pe / A
    if s_e >= des.rod.Re / 2:
        return A * des.rod.Re * (1 - des.rod.Re / (4 * s_e))
    return pe


def two_segment_pinned(EI1: float, a: float, EI2: float, b: float) -> float:
    """Exact lowest critical load of a pinned-pinned bar with two segments (a: EI1, b: EI2):
    k1 cos(k1 a) sin(k2 b) + k2 sin(k1 a) cos(k2 b) = 0, k_i^2 = P/EI_i."""

    def f(P):
        k1, k2 = math.sqrt(P / EI1), math.sqrt(P / EI2)
        return k1 * math.cos(k1 * a) * math.sin(k2 * b) + k2 * math.sin(k1 * a) * math.cos(k2 * b)

    lo = 1e-6 * min(EI1, EI2) / (a + b) ** 2
    hi_cap = math.pi**2 * max(EI1, EI2) / (a + b) ** 2 * 1.0001
    Ps = np.linspace(lo, hi_cap, 4000)
    vals = [f(P) for P in Ps]
    for i in range(len(Ps) - 1):
        if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
            return brentq(f, Ps[i], Ps[i + 1])
    raise ValueError("no root")


def analyse(des: Design) -> dict:
    """Buckling verdict for the fully extended cylinder at its push load."""
    col = build_column(des)
    col0 = build_column(des, imperfect=False)
    Pcr, mode = critical_load(col0.beam)
    F_op = des.p_relief * des.A_A  # highest push the circuit allows
    # usable load: first yield with imperfections (bisection on utilisation = 1)
    hi = 0.999 * Pcr
    if fibre_stress(des, col, hi)["util"] < 1.0:
        Py = hi
    else:
        Py = brentq(lambda P: fibre_stress(des, col, P)["util"] - 1.0, 1.0, hi, xtol=10.0)
    at_op = fibre_stress(des, col, F_op)
    lay = col.lay
    theta = (des.mm("piston_clear_mm") + des.mm("gland_clear_mm")) / (2 * lay["guide"])
    return {
        "L_ext_m": lay["L"],
        "guide_m": lay["guide"],
        "kink_mrad": theta * 1e3,
        "F_op_kN": F_op / 1e3,
        "P_cr_kN": Pcr / 1e3,
        "P_euler_rod_kN": euler_rod_only(des) / 1e3,
        "P_johnson_rod_kN": johnson_rod_only(des) / 1e3,
        "P_yield_kN": Py / 1e3,
        "sf_elastic": Pcr / F_op,
        "sf_yield": Py / F_op,
        "util_at_op": at_op["util"],
        "w_max_at_op_mm": at_op["w_max_mm"],
        "mode": mode,
        "x": col.beam.x,
        "col": col,
    }


def load_deflection(des: Design, n: int = 30) -> dict:
    col = build_column(des)
    col0 = build_column(des, imperfect=False)
    Pcr, _ = critical_load(col0.beam)
    Ps = np.linspace(0.02, 0.97, n) * Pcr
    w = [fibre_stress(des, col, P)["w_max_mm"] for P in Ps]
    u = [fibre_stress(des, col, P)["util"] for P in Ps]
    return {"P": Ps, "w_mm": np.array(w), "util": np.array(u), "Pcr": Pcr}
