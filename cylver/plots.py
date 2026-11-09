"""PNG figures for the calculation report. Palette: the validated default categorical order
(blue, orange, aqua, yellow, magenta); limits drawn in the 'critical' status red with a text label."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import buckling, dynamics, seals, stress  # noqa: E402

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
LIMIT = "#d03b3b"
INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update(
    {
        "figure.facecolor": SURF,
        "axes.facecolor": SURF,
        "axes.edgecolor": "#bdbcb7",
        "axes.labelcolor": INK2,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "text.color": INK,
        "axes.grid": True,
        "grid.color": "#e6e5e1",
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.6,
        "font.size": 9,
        "axes.titlesize": 10,
        "legend.frameon": False,
        "legend.fontsize": 8,
    }
)


def _limit(ax, y, text, horizontal=True):
    if horizontal:
        ax.axhline(y, color=LIMIT, lw=1.0, ls="--")
        ax.annotate(text, (1.0, y), xycoords=("axes fraction", "data"), xytext=(-4, 3), textcoords="offset points", ha="right", va="bottom", color=LIMIT, fontsize=8)
    else:
        ax.axvline(y, color=LIMIT, lw=1.0, ls="--")
        ax.annotate(text, (y, 1.0), xycoords=("data", "axes fraction"), xytext=(3, -12), textcoords="offset points", ha="left", color=LIMIT, fontsize=8)


def _save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def fig_buckling(des, res, path: Path):
    b = res["buckling"]
    col = b["col"]
    ld = buckling.load_deflection(des, 24)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    F = b["F_op_kN"] * 1e3
    fs = buckling.fibre_stress(des, col, F)
    x = col.beam.x * 1e3
    ax[0].plot(x, col.w0 * 1e3, color=C[1], ls=":", label="initial imperfection (clearance tilt + bows)")
    ax[0].plot(x, (fs["w"] + col.w0) * 1e3, color=C[0], label=f"deflection at max push {b['F_op_kN']:.0f} kN")
    for k, (a, c, yy, ha, dx) in {"piston": (b["col"].lay["x_piston_c"], C[2], 0.97, "right", -3), "gland bearing": (b["col"].lay["x_gland_c"], C[3], 0.90, "left", 3)}.items():
        ax[0].axvline(a * 1e3, color=c, lw=0.8)
        ax[0].annotate(k, (a * 1e3, yy), xycoords=("data", "axes fraction"), xytext=(dx, 0), textcoords="offset points", ha=ha, color=INK2, fontsize=8)
    ax[0].set_xlabel("position from rear pin (mm)")
    ax[0].set_ylabel("lateral deflection (mm)")
    ax[0].set_title("Fully extended column, second-order deflection")
    ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    ax[1].plot(ld["P"] / 1e3, ld["util"], color=C[0], label="max fibre stress / Re")
    ax[1].set_xlabel("axial push (kN)")
    ax[1].set_ylabel("utilisation")
    _limit(ax[1], 1.0, "first yield")
    ax[1].axvline(b["F_op_kN"], color=INK2, lw=1, ls="-.")
    ax[1].annotate("push at relief", (b["F_op_kN"], 0.05), xytext=(3, 0), textcoords="offset points", color=INK2, fontsize=8)
    for val, name, c in ((b["P_cr_kN"], "elastic P_cr (FE)", C[1]), (b["P_euler_rod_kN"], "Euler, rod only", C[4]), (b["P_yield_kN"], "usable load", C[2])):
        ax[1].axvline(val, color=c, lw=1)
        ax[1].annotate(f"{name} {val:.0f}", (val, 0.6 if c != C[2] else 0.3), xytext=(3, 0), textcoords="offset points", color=INK2, fontsize=8, rotation=90)
    ax[1].set_ylim(0, 1.6)
    ax[1].set_title(f"Usable load = {b['sf_yield']:.2f} x push at relief")
    _save(fig, path)


def fig_junction(des, res, path: Path):
    p = des.p_relief
    xs = np.linspace(0, 0.12, 200)
    tot = [stress.junction_bending(p, des.ri, des.ro, des.tube.nu, x)["total"] / 1e6 for x in xs]
    jb = res["junction"]
    xw = des.mm("base_weld_offset_mm")
    x0 = math.pi / (4 * jb["beta"])
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    ax.plot(xs * 1e3, tot, color=C[0], label="shell theory, rigid base (membrane + edge bending)")
    from .axisym import junction_fe

    fe = junction_fe(des)
    m = fe["z"] <= 0.12
    ax.plot(fe["z"][m] * 1e3, fe["total"][m] * p / 1e6, color=C[1], ls="--", label="axisymmetric FE, flexible base (linearised)")
    ax.axhline(jb["membrane"] * p / 1e6, color=C[2], ls=":", label="membrane only")
    ax.axvline(xw * 1e3, color=C[1], lw=1.2)
    ax.annotate(f"weld at {xw * 1e3:.0f} mm", (xw * 1e3, 0.92), xycoords=("data", "axes fraction"), xytext=(4, 0), textcoords="offset points", color=INK2)
    ax.axvline(x0 * 1e3, color=INK2, lw=0.8, ls="-.")
    ax.annotate(f"zero moment, beta x = pi/4: {x0 * 1e3:.0f} mm", (x0 * 1e3, 0.80), xycoords=("data", "axes fraction"), xytext=(4, 0), textcoords="offset points", color=INK2)
    ax.set_xlabel("distance from the base inner face (mm)")
    ax.set_ylabel("stress at relief pressure (MPa)")
    ax.set_title(f"Tube wall next to the base, rm/t = {jb['rm_over_t']:.1f}")
    ax.legend(loc="upper right", bbox_to_anchor=(1, 0.7))
    _save(fig, path)


def fig_duty(des, res, path: Path):
    d = res["duty"]
    fig, ax = plt.subplots(2, 1, figsize=(10, 5.2), sharex=True)
    ax[0].plot(d["t"], d["x"] * 1e3, color=C[0], label="piston position")
    ax[0].set_ylabel("position (mm)")
    ax0b = ax[0]
    for i, (n, a, b) in enumerate(d["phases"]):
        ax0b.axvspan(a, b, color="#efeee9" if i % 2 else SURF, zorder=0)
        ax0b.annotate(n, ((a + b) / 2, 1.0), xycoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=7, color=INK2)
    ax[1].plot(d["t"], d["pA"] / 1e5, color=C[0], label="cap chamber pA")
    ax[1].plot(d["t"], d["pB"] / 1e5, color=C[1], ls="--", label="rod chamber pB")
    _limit(ax[1], des.p_relief / 1e5, f"relief {des.p_relief / 1e5:.0f} bar")
    ax[1].set_ylabel("pressure (bar)")
    ax[1].set_xlabel("time (s)")
    ax[1].legend(loc="upper left")
    ax[0].set_title("Design duty cycle (simulated): the pressure history that feeds the fatigue check", pad=14)
    _save(fig, path)


def fig_fatigue(des, res, path: Path):
    fa = res["fatigue"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    names = [f["loc"] for f in fa]
    vals = [min(f["stress_factor"], 5) for f in fa]
    ax[0].set_axisbelow(True)
    ax[0].barh(names, vals, color=C[0], height=0.5)
    for i, f in enumerate(fa):
        ax[0].annotate(f"{f['stress_factor']:.2f}", (vals[i], i), xytext=(4, 0), textcoords="offset points", va="center", color=INK)
    _limit(ax[0], 1.0, "required life", horizontal=False)
    ax[0].set_xlabel("stress factor to use up the required life (>1 = margin, capped at 5)")
    ax[0].set_title("Fatigue margin by location")
    ax[0].invert_yaxis()
    g = res["fatigue_gov"]
    cyc = g["cycles"]
    r = np.array([c[0] for c in cyc]) / 1e6
    w = np.array([c[2] for c in cyc])
    ax[1].hist(r, bins=20, weights=w, color=C[1], rwidth=0.9)
    ax[1].set_xlabel(f"stress range at {g['loc']} (MPa)")
    ax[1].set_ylabel("cycles per duty cycle")
    ax[1].set_title("Rainflow count, governing location")
    _save(fig, path)


def fig_cushion(des, res, path: Path):
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    c40 = next(c for c in res["cushion"] if c["side"] == "A" and c["T_C"] == 40.0)
    r = c40["r"]
    lc = des.g["cushion"]["length_mm"]
    m = r["x"] * 1e3 < lc * 2.5
    ax[0].plot(r["x"][m] * 1e3, r["pA"][m] / 1e5, color=C[0], label="cushion chamber pressure (40 C)")
    ax2 = ax[0].twinx()
    ax2.plot(r["x"][m] * 1e3, np.abs(r["v"][m]), color=C[1], ls="--", label="piston speed")
    ax2.set_ylabel("speed (m/s)")
    ax2.grid(False)
    ax[0].axvline(lc, color=INK2, lw=0.8)
    ax[0].annotate("cushion engages", (lc, 0.5), xycoords=("data", "axes fraction"), xytext=(3, 0), textcoords="offset points", color=INK2, fontsize=8)
    _limit(ax[0], des.p_relief / 1e5, "port relief setting")
    ax[0].invert_xaxis()
    ax[0].set_xlabel("distance to the end stop (mm)")
    ax[0].set_ylabel("pressure (bar)")
    ax[0].set_title("Retracting into the cap-end cushion")
    h1, l1 = ax[0].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax[0].legend(h1 + h2, l1 + l2, loc="lower left")
    sw = res.get("cushion_sweep")
    if sw:
        a = [s["needle_mm2"] for s in sw]
        ax[1].plot(a, [s["p_peak_bar"] for s in sw], marker="o", ms=5, color=C[0], label="peak cushion pressure")
        ax3 = ax[1].twinx()
        ax3.plot(a, [s["v_impact"] for s in sw], marker="s", ms=5, ls="--", color=C[1], label="impact speed at the stop")
        ax3.set_ylabel("impact speed (m/s)")
        ax3.grid(False)
        ax[1].set_xscale("log")
        ax[1].set_xticks(a)
        ax[1].set_xticklabels([f"{v:g}" for v in a])
        ax[1].minorticks_off()
        ax[1].set_xlabel("needle valve area (mm^2)")
        ax[1].set_ylabel("peak pressure (bar)")
        _limit(ax[1], des.proof_factor * max(des.p_rated, des.p_relief) / 1e5, "proof pressure")
        ax[1].set_title("Needle setting: pressure vs impact trade")
        h1, l1 = ax[1].get_legend_handles_labels()
        h2, l2 = ax3.get_legend_handles_labels()
        ax[1].legend(h1 + h2, l1 + l2, loc="center right")
    _save(fig, path)


def fig_friction(des, res, path: Path):
    p = res["p_hold_bar"] * 1e5
    v = np.geomspace(1e-5, 0.5, 400)
    F = [seals.total_friction(des, vi, p, 0.0, 1e-6) for vi in v]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].semilogx(v * 1e3, np.array(F) / 1e3, color=C[0], label=f"rod + piston seal friction at {res['p_hold_bar']:.0f} bar")
    vc = res["v_stick_linear"] * 1e3
    ax[0].axvspan(1e-2, vc, color="#f6dede", zorder=0)
    ax[0].annotate(f"falling branch: stick-slip below {vc:.1f} mm/s", (0.012, 0.92), xycoords=("data", "axes fraction"), color=INK2, fontsize=8)
    ax[0].axvline(des.o["speed_min_smooth_mm_s"], color=LIMIT, ls="--", lw=1)
    ax[0].annotate("required smooth speed", (des.o["speed_min_smooth_mm_s"], 0.45), xycoords=("data", "axes fraction"), xytext=(3, 0), textcoords="offset points", color=LIMIT, fontsize=8)
    ax[0].set_xlabel("speed (mm/s)")
    ax[0].set_ylabel("friction (kN)")
    ax[0].set_title("Seal friction, calibrated to measured data (Sarbu et al. 2024)")
    ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.18))
    for i, (vv, s) in enumerate(sorted(res["stick_sim"].items())):
        ax[1].plot(s["t"], s["v"] / vv, color=C[i % len(C)], lw=1.0, label=f"{vv * 1e3:.1f} mm/s: {'stick-slip' if s['stick'] else ('growing' if s['growth'] > 1.05 else 'smooth')}")
    ax[1].set_xlim(0, min(s["t"][-1] for s in res["stick_sim"].values()))
    ax[1].set_xlabel("time (s)")
    ax[1].set_ylabel("speed / commanded speed")
    ax[1].set_title("Constant-flow drive, kicked from steady sliding")
    ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
    _save(fig, path)


def fig_frequency(des, res, path: Path):
    f = res["freq"]
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    for i, k in enumerate([k for k in f if k != "x"]):
        ax.plot(f["x"] * 1e3, f[k], color=C[i], ls=["-", "--", ":"][i % 3], label=f"free air {k * 100:.1f} % at 1 bar")
    ax.set_xlabel("piston position (mm)")
    ax.set_ylabel("hydraulic natural frequency (Hz)")
    lp = res["loop"]
    ax.set_title(f"Oil-column frequency, valve closed; damping ratio {lp['zeta']:.4f}")
    ax.legend(loc="upper center")
    _save(fig, path)


def fig_aoa(trades: list[dict], path: Path):
    fig, axes = plt.subplots(1, len(trades), figsize=(4.2 * len(trades), 3.4))
    for ax, t in zip(np.atleast_1d(axes), trades):
        r = t["result"]
        names = [a["name"] + ("" if a.get("viable", True) else "\n(not viable)") for a in t["alternatives"]]
        ax.barh(names, r["total"], color=C[0], height=0.5)
        for i, (tot, ws) in enumerate(zip(r["total"], r["win_share"])):
            ax.annotate(f"{tot:.1f} ({ws * 100:.0f} %)", (tot, i), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8)
        ax.set_xlim(0, 12)
        ax.set_axisbelow(True)
        ax.invert_yaxis()
        ax.set_title(t["title"])
        ax.set_xlabel("score 0-10 (weight sets won)")
    _save(fig, path)


def fig_fmea(rows: list[dict], path: Path):
    fig, ax = plt.subplots(figsize=(8, 3.8))
    names = [r["mode"] for r in rows]
    ax.barh(names, [r["RPN"] for r in rows], color=C[1], height=0.55)
    for i, r in enumerate(rows):
        ax.annotate(f"{r['RPN']}  (S{r['S']} O{r['O']} D{r['D']}, margin {r['margin']:.2f})", (r["RPN"], i), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, max(r["RPN"] for r in rows) * 1.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("RPN = S x O x D (O from the computed margin)")
    ax.set_title("FMEA ranking")
    _save(fig, path)
