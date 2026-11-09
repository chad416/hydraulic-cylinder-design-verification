"""Axisymmetric finite-element check of the tube-to-base junction.

The weld finding of this project rests on the thin-shell edge solution with a RIGID base
(stress.junction_bending). This module checks it with a continuum model that has neither assumption:
4-node axisymmetric elements (bilinear, 2x2 Gauss), linear elastic, one material, the base as a real
plate of thickness t_base that can rotate, the tube wall meshed through its thickness.

Loads: internal pressure on the bore and on the base inner face; at the far end of the tube the
closed-end axial stress p ri^2/(ro^2 - ri^2) (the gland's pull). Supports: u = 0 on the axis, one node
fixed axially (the load set is self-equilibrated).

Output along the tube: the through-wall LINEARISED longitudinal stress (membrane + bending), i.e. the
structural stress that the fatigue assessment uses. The sharp re-entrant corner itself is singular in
a linear model and is not used; a real weld has a toe radius and is covered by the weld detail category.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from .model import Design

GP = np.array([-1, 1]) / math.sqrt(3)


def _D(E: float, nu: float) -> np.ndarray:
    c = E / ((1 + nu) * (1 - 2 * nu))
    return c * np.array([[1 - nu, nu, nu, 0], [nu, 1 - nu, nu, 0], [nu, nu, 1 - nu, 0], [0, 0, 0, (1 - 2 * nu) / 2]])


def _B(xy: np.ndarray, xi: float, eta: float) -> tuple[np.ndarray, float, float]:
    """Strain-displacement matrix [er, eth, ez, grz] for a Q4 element at (xi, eta); returns B, detJ, r."""
    dN = 0.25 * np.array([[-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)], [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]])
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta), (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    J = dN @ xy
    detJ = np.linalg.det(J)
    dNx = np.linalg.solve(J, dN)
    r = N @ xy[:, 0]
    B = np.zeros((4, 8))
    B[0, 0::2] = dNx[0]
    B[1, 0::2] = N / max(r, 1e-12)
    B[2, 1::2] = dNx[1]
    B[3, 0::2] = dNx[1]
    B[3, 1::2] = dNx[0]
    return B, detJ, r


def junction_fe(des: Design, p: float = 1e6, nt: int = 8, nb: int = 24, nzb: int = 8, n_decay: float = 6.0, nz: int = 160) -> dict:
    ri, ro, tb = des.ri, des.ro, des.t_base
    E, nu = des.tube.E, des.tube.nu
    t = ro - ri
    rm = 0.5 * (ri + ro)
    beta = (3 * (1 - nu**2)) ** 0.25 / math.sqrt(rm * t)
    L = n_decay * math.pi / beta
    R = np.concatenate([np.linspace(0, ri, nb + 1), np.linspace(ri, ro, nt + 1)[1:]])
    zt = L * (np.linspace(0, 1, nz + 1) ** 1.6)  # graded: fine near the base
    Z = np.concatenate([np.linspace(-tb, 0, nzb + 1), zt[1:]])
    nr_, nz_ = len(R), len(Z)
    node_id = -np.ones((nr_, nz_), dtype=int)
    elems = []
    for j in range(nz_ - 1):
        zc = 0.5 * (Z[j] + Z[j + 1])
        for i in range(nr_ - 1):
            rc = 0.5 * (R[i] + R[i + 1])
            if zc > 0 and rc < ri:
                continue  # oil space
            elems.append((i, j))
    for i, j in elems:
        for a, b in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)):
            if node_id[a, b] < 0:
                node_id[a, b] = 0
    ids = np.argwhere(node_id >= 0)
    for k, (a, b) in enumerate(ids):
        node_id[a, b] = k
    nn = len(ids)
    xy_all = np.array([[R[a], Z[b]] for a, b in ids])
    D = _D(E, nu)
    rows, cols, vals = [], [], []
    for i, j in elems:
        conn = [node_id[i, j], node_id[i + 1, j], node_id[i + 1, j + 1], node_id[i, j + 1]]
        xy = xy_all[conn]
        ke = np.zeros((8, 8))
        for xi in GP:
            for eta in GP:
                B, dJ, r = _B(xy, xi, eta)
                ke += B.T @ D @ B * dJ * 2 * math.pi * r
        dofs = np.array([[2 * n, 2 * n + 1] for n in conn]).ravel()
        rr, cc = np.meshgrid(dofs, dofs, indexing="ij")
        rows.append(rr.ravel())
        cols.append(cc.ravel())
        vals.append(ke.ravel())
    K = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * nn, 2 * nn)).tocsr()
    f = np.zeros(2 * nn)

    def edge_load(n1, n2, comp, trac):
        """Consistent nodal loads of a uniform traction on a straight edge (axisymmetric, 2 pi r ds)."""
        (r1, z1), (r2, z2) = xy_all[n1], xy_all[n2]
        ds = math.hypot(r2 - r1, z2 - z1)
        # integral of N_a * 2 pi r ds for linear N: (2 r1 + r2)/6 and (r1 + 2 r2)/6
        f[2 * n1 + comp] += trac * 2 * math.pi * ds * (2 * r1 + r2) / 6
        f[2 * n2 + comp] += trac * 2 * math.pi * ds * (r1 + 2 * r2) / 6

    ib = int(np.argmin(abs(R - ri)))
    jb0 = nzb  # z = 0
    for j in range(jb0, nz_ - 1):  # bore surface: pressure pushes +r
        edge_load(node_id[ib, j], node_id[ib, j + 1], 0, p)
    for i in range(ib):  # base inner face: pressure pushes -z
        edge_load(node_id[i, jb0], node_id[i + 1, jb0], 1, -p)
    s_end = p * ri**2 / (ro**2 - ri**2)
    for i in range(ib, nr_ - 1):  # far end: closed-end axial stress (+z)
        edge_load(node_id[i, nz_ - 1], node_id[i + 1, nz_ - 1], 1, s_end)
    fixed = [2 * node_id[0, j] for j in range(nzb + 1)]  # axis: u = 0
    fixed.append(2 * node_id[0, 0] + 1)  # one node axially
    free = np.setdiff1d(np.arange(2 * nn), fixed)
    u = np.zeros(2 * nn)
    u[free] = spsolve(K[free][:, free], f[free])
    # longitudinal stress at element centroids through the wall, linearised per z station
    zs, mem, ben = [], [], []
    for j in range(jb0, nz_ - 1):
        s_col, r_col, w_col = [], [], []
        for i in range(ib, nr_ - 1):
            conn = [node_id[i, j], node_id[i + 1, j], node_id[i + 1, j + 1], node_id[i, j + 1]]
            B, _, r = _B(xy_all[conn], 0.0, 0.0)
            ue = np.array([[u[2 * n], u[2 * n + 1]] for n in conn]).ravel()
            s_col.append((D @ B @ ue)[2])
            r_col.append(r)
            w_col.append(R[i + 1] - R[i])
        s_col, r_col, w_col = map(np.array, (s_col, r_col, w_col))
        m = np.sum(s_col * w_col) / t
        b = 6 / t**2 * np.sum(s_col * (r_col - rm) * w_col)
        zs.append(0.5 * (Z[j] + Z[j + 1]))
        mem.append(m)
        ben.append(b)
    zs, mem, ben = map(np.array, (zs, mem, ben))
    sgn = np.sign(ben[0]) if ben[0] != 0 else 1.0
    cross = np.nonzero(np.diff(np.sign(ben)) != 0)[0]
    z0 = float(np.interp(0, [ben[cross[0]] * sgn, ben[cross[0] + 1] * sgn][::-1], [zs[cross[0]], zs[cross[0] + 1]][::-1])) if len(cross) else float("nan")
    return {"z": zs, "membrane": mem / p, "bending": ben / p, "total": (mem + np.abs(ben)) / p, "z_zero": z0, "beta": beta, "n_elem": len(elems), "n_dof": 2 * nn}


def compare(des: Design) -> dict:
    """FE vs rigid-base shell solution at the weld position and at the zero-moment point."""
    from .stress import junction_bending

    fe = junction_fe(des)
    xw = des.mm("base_weld_offset_mm")
    sh_w = junction_bending(1.0, des.ri, des.ro, des.tube.nu, xw)
    fe_w = float(np.interp(xw, fe["z"], fe["total"]))
    sh0 = math.pi / (4 * sh_w["beta"])
    x_peak_fe = float(fe["z"][2])
    return {
        "weld_mm": xw * 1e3,
        "shell_total_per_p": sh_w["total"],
        "fe_total_per_p": fe_w,
        "ratio_fe_shell": fe_w / sh_w["total"],
        "zero_shell_mm": sh0 * 1e3,
        "zero_fe_mm": fe["z_zero"] * 1e3,
        "fe": fe,
        "near_base_fe_per_p": float(np.interp(x_peak_fe, fe["z"], fe["total"])),
    }


_CACHE: dict = {}


def weld_stress_per_pa(des: Design) -> tuple[float, float]:
    """(FE structural stress per Pa at the weld position, shell value) - cached per geometry."""
    from .stress import junction_bending

    key = (round(des.ri, 6), round(des.ro, 6), round(des.t_base, 6), des.tube.E, des.tube.nu, round(des.mm("base_weld_offset_mm"), 6))
    if key not in _CACHE:
        fe = junction_fe(des)
        xw = des.mm("base_weld_offset_mm")
        _CACHE[key] = (float(np.interp(xw, fe["z"], fe["total"])), junction_bending(1.0, des.ri, des.ro, des.tube.nu, xw)["total"])
    return _CACHE[key]
