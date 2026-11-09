"""Valve-controlled cylinder: one nonlinear model used for cushioning, pressure peaks, the duty-cycle
pressure history (input to fatigue), and the hydraulic natural frequency.

States: x (m, 0 = fully retracted), v, pA, pB (cylinder chambers), pAl, pBl (hose/line volumes between
valve and ports), u (valve spool, -1..1, first-order lag). Pressures are gauge (Pa).

Physics and sources (see docs/DOSSIER.md section 9):
- Chamber continuity dp/dt = (K_e / V)(sum Q - dV/dt): Merritt, Hydraulic Control Systems (1967) ch. 2-5;
  Jelali & Kroll, Hydraulic Servo-systems (2003) ch. 3.
- Effective bulk modulus of oil + entrained air + container in series (Merritt's form, discussed by
  Gholizadeh, Burton & Schoenau 2011, IJFP 12(3)); the free air volume follows the gas law (polytropic
  exponent kappa), which ignores air dissolving at higher pressure (Ruan & Burton 2006 per the same
  survey), so the model is conservative (too soft) at high pressure.
- Valve edges: turbulent orifice q = K u sqrt(dp), smoothed below 0.1 bar (laminar transition).
- Cushion: outflow from the closing chamber through the tapered spear gap (laminar annulus,
  dp = 12 mu Q integral(dz / (pi d h^3))) in series with the gap entrance (orifice), in parallel with the
  needle valve; free inflow through the check valve on reversal.
- Port relief valves on both lines (linear above the setting). Pump: ideal pressure source at p_supply
  (flow not limited) - stated limit.
- Friction: seals.py Stribeck model.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .model import G, Design
from .seals import SEALS, total_friction
from .stress import bore_dilation

P_ATM = 1.013e5
P_TR = 1e4  # laminar transition of the orifice law (Pa)


def k_eff(p: float, K_l: float, eps0: float, kappa: float, K_c: float) -> float:
    """Effective bulk modulus: liquid, free air (fraction eps0 at atmospheric pressure) and container."""
    pa = max(p + P_ATM, 0.02 * P_ATM)
    vg = eps0 * (P_ATM / pa) ** (1 / kappa)
    x = vg / (1 - eps0 + vg)
    inv = (1 - x) / K_l + x / (kappa * pa) + 1 / K_c
    return 1 / inv


def orifice(K: float, dp: float) -> float:
    return K * dp / math.sqrt(math.sqrt(dp * dp + P_TR * P_TR))


@dataclass
class Plant:
    des: Design
    T_C: float = 40.0
    eps0: float | None = None
    K_hose: float | None = None
    Kv: float = 0.0  # valve edge coefficient (m^3/s per sqrt(Pa))
    needle_scale: float = 1.0
    cushion_on: bool = True
    coulomb_on: bool = True
    log: dict = field(default_factory=dict)

    def __post_init__(self):
        d = self.des
        g = d.g
        self.mu = d.oil.viscosity(self.T_C)
        self.eps0 = d.o["air_fraction_atm"] if self.eps0 is None else self.eps0
        self.K_hose = g["hose"]["K_hose_GPa"] * 1e9 if self.K_hose is None else self.K_hose
        u = bore_dilation(1e6, d.ri, d.ro, d.tube.E, d.tube.nu)
        self.K_tube = 1e6 * d.ri / (2 * u)
        h = g["hose"]
        self.V_line = math.pi / 4 * (h["id_mm"] * 1e-3) ** 2 * h["length_m"]
        self.V0 = g["dead_volume_cm3"] * 1e-6
        c = g["cushion"]
        self.lc = c["length_mm"] * 1e-3
        self.h_tip = c["gap_tip_mm"] * 1e-3 / 2  # diametral -> radial
        self.h_root = c["gap_root_mm"] * 1e-3 / 2
        self.A_n = {"A": c["needle_mm2"] * 1e-6, "B": c.get("needle_rod_mm2", c["needle_mm2"]) * 1e-6}
        self.ds = {"A": c["spear_d_cap_mm"] * 1e-3, "B": c["spear_d_rod_mm"] * 1e-3}
        self.K_port = 0.7 * math.pi / 4 * (g["port_d_mm"] * 1e-3) ** 2 * math.sqrt(2 / d.oil.rho)
        self.k_stop = 2e9
        self.c_stop = 2e6
        self.F_load = d.o["load_extend_kN"] * 1e3
        if not self.Kv:
            self.Kv = size_valve(self)

    # --- pieces
    def K_ch(self, p: float) -> float:
        return k_eff(p, self.des.oil.K, self.eps0, self.des.oil.kappa, self.K_tube)

    def K_ln(self, p: float) -> float:
        return k_eff(p, self.des.oil.K, self.eps0, self.des.oil.kappa, self.K_hose)

    def engaged(self, side: str, x: float) -> float:
        s = self.des.stroke
        return max(0.0, self.lc - x) if side == "A" else max(0.0, x - (s - self.lc))

    def cushion_out(self, side: str, xi: float, dp: float) -> float:
        """Outflow (m^3/s) from the cushion chamber for pressure drop dp > 0 at engagement xi."""
        if dp <= 0:
            return 0.0
        a = (self.h_tip - self.h_root) / self.lc
        hx = self.h_tip - a * xi
        if a > 1e-12:
            integ = (1 / hx**2 - 1 / self.h_tip**2) / (2 * a)
        else:
            integ = xi / self.h_tip**3
        d = self.ds[side]
        R = 12 * self.mu * integ / (math.pi * d)
        A_ent = math.pi * d * hx
        k2 = self.des.oil.rho / (2 * (0.7 * A_ent) ** 2)
        q_gap = (-R + math.sqrt(R * R + 4 * k2 * dp)) / (2 * k2) if R > 0 else math.sqrt(dp / k2)
        q_nd = 0.7 * self.A_n[side] * self.needle_scale * math.sqrt(2 * dp / self.des.oil.rho)
        return q_gap + q_nd

    def port_flow(self, side: str, x: float, p_line: float, p_ch: float) -> float:
        """Flow line -> chamber (m^3/s)."""
        xi = self.engaged(side, x) if self.cushion_on else 0.0
        if xi > 0 and p_ch > p_line:
            return -self.cushion_out(side, xi, p_ch - p_line)
        return orifice(self.K_port, p_line - p_ch)

    # --- the model
    def rhs(self, t: float, y: np.ndarray, u_cmd) -> list[float]:
        d = self.des
        x, v, pA, pB, pAl, pBl, u = y
        ps, pt = d.p_supply, 0.0
        uc = u_cmd(t) if callable(u_cmd) else u_cmd
        Kv = self.Kv
        up, un = max(u, 0.0), max(-u, 0.0)
        q_PA = Kv * up * orifice(1.0, ps - pAl) - Kv * un * orifice(1.0, pAl - pt)
        q_PB = Kv * un * orifice(1.0, ps - pBl) - Kv * up * orifice(1.0, pBl - pt)
        k_rel = 1e-9  # m^3/s per Pa above the setting
        q_relA = k_rel * max(0.0, pAl - d.p_relief)
        q_relB = k_rel * max(0.0, pBl - d.p_relief)
        qA = self.port_flow("A", x, pAl, pA)
        qB = self.port_flow("B", x, pBl, pB)
        VA = self.V0 + d.A_A * x
        VB = self.V0 + d.A_B * (d.stroke - x)
        dpA = self.K_ch(pA) / VA * (qA - d.A_A * v)
        dpB = self.K_ch(pB) / VB * (qB + d.A_B * v)
        dpAl = self.K_ln(pAl) / (self.V_line + 1e-6) * (q_PA - qA - q_relA)
        dpBl = self.K_ln(pBl) / (self.V_line + 1e-6) * (q_PB - qB - q_relB)
        F_stop = 0.0
        if x < 0:
            F_stop = -self.k_stop * x - self.c_stop * min(v, 0.0)
        elif x > d.stroke:
            F_stop = -self.k_stop * (x - d.stroke) - self.c_stop * max(v, 0.0)
        if self.coulomb_on:
            F_f = total_friction(d, v, pA, pB, 1e-4)
        else:  # viscous part only (used for the small-signal model around standstill)
            fv = lambda k: sum(e["fv"] for e, _ in SEALS[k]["elements"])  # noqa: E731
            F_f = math.pi * 1e3 * (d.d * fv(d.g["rod_seal"]) + d.D * fv(d.g["piston_seal"])) * v
        F_g = self.F_load  # constant load opposing extension (gravity-type, includes the payload)
        dv = (pA * d.A_A - pB * d.A_B - F_f - F_g + F_stop) / d.mass
        du = (uc - u) / d.g["valve_tau_s"]
        return [v, dv, dpA, dpB, dpAl, dpBl, du]

    def hold_state(self, x: float) -> np.ndarray:
        """Static equilibrium with the valve closed: rod side at the pressure its trapped oil had when the
        load stopped (taken as 5 bar), cap side carrying the load."""
        d = self.des
        pB = 5e5
        pA = (self.F_load + pB * d.A_B) / d.A_A
        return np.array([x, 0.0, pA, pB, pA, pB, 0.0])


def size_valve(pl: Plant) -> float:
    """Valve edge coefficient so that full opening extends at speed_max against the load (steady state)."""
    d = pl.des
    v = d.o["speed_max_mm_s"] * 1e-3

    def f(K):
        pA = d.p_supply - (d.A_A * v / K) ** 2
        pB = (d.A_B * v / K) ** 2
        return pA * d.A_A - pB * d.A_B - total_friction(d, v, pA, pB) - pl.F_load

    lo = d.A_A * v / math.sqrt(d.p_supply) * 1.0001
    return brentq(f, lo, lo * 1e3)


def duty_command(des: Design) -> tuple[callable, float, list[tuple[str, float, float]]]:
    """The design duty cycle as a valve command u(t): extend with one mid-stroke stop, hold, retract, hold.
    Returns (u(t), duration, phase list)."""
    v = des.o["speed_max_mm_s"] * 1e-3
    s = des.stroke
    t_ext = s / v
    ramp, stop_ramp = 0.15, 0.05
    ph = []
    t = 0.0
    ph.append(("extend 1", t, t + 0.5 * t_ext))
    t += 0.5 * t_ext
    ph.append(("mid-stroke stop", t, t + 0.6))
    t += 0.6
    ph.append(("extend 2", t, t + 0.5 * t_ext + 1.2))
    t += 0.5 * t_ext + 1.2
    ph.append(("hold extended", t, t + 1.0))
    t += 1.0
    t_ret = 1.3 * t_ext / max(des.g.get("retract_command", 1.0), 0.3) + 1.0
    ph.append(("retract", t, t + t_ret))
    t += t_ret
    ph.append(("hold retracted", t, t + 1.0))
    t += 1.0
    ur = des.g.get("retract_command", 1.0)
    level = {"extend 1": 1.0, "mid-stroke stop": 0.0, "extend 2": 1.0, "hold extended": 0.0, "retract": -ur, "hold retracted": 0.0}

    def u(tt: float) -> float:
        prev = 0.0
        for name, a, b in ph:
            if a <= tt < b:
                r = stop_ramp if level[name] == 0.0 else ramp
                frac = min(1.0, (tt - a) / r)
                return prev + (level[name] - prev) * frac
            prev = level[name]
        return 0.0

    return u, t, ph


def simulate(pl: Plant, u_cmd, t_end: float, y0: np.ndarray | None = None, dt_out: float = 1e-3) -> dict:
    y0 = pl.hold_state(0.0) if y0 is None else y0
    t_eval = np.arange(0.0, t_end, dt_out)
    sol = solve_ivp(
        lambda t, y: pl.rhs(t, y, u_cmd),
        (0.0, t_end),
        y0,
        method="LSODA",
        t_eval=t_eval,
        max_step=2e-3,
        rtol=1e-6,
        atol=[1e-7, 1e-6, 50.0, 50.0, 50.0, 50.0, 1e-6],
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    x, v, pA, pB, pAl, pBl, u = sol.y
    return {"t": sol.t, "x": x, "v": v, "pA": pA, "pB": pB, "pAl": pAl, "pBl": pBl, "u": u}


def run_duty(des: Design, **kw) -> dict:
    pl = Plant(des, **kw)
    u, T, ph = duty_command(des)
    r = simulate(pl, u, T)
    r["phases"] = ph
    r["Kv"] = pl.Kv
    return r


# ---------------------------------------------------------------- natural frequency
def omega_formula(pl: Plant, x: float, pA: float, pB: float) -> float:
    """omega_h^2 = (A_A^2 / C_A + A_B^2 / C_B) / m with C = V_ch/K_ch + V_line/K_line (valve closed)."""
    d = pl.des
    VA = pl.V0 + d.A_A * x
    VB = pl.V0 + d.A_B * (d.stroke - x)
    CA = VA / pl.K_ch(pA) + pl.V_line / pl.K_ln(pA)
    CB = VB / pl.K_ch(pB) + pl.V_line / pl.K_ln(pB)
    return math.sqrt((d.A_A**2 / CA + d.A_B**2 / CB) / d.mass)


def linearise(pl: Plant, y: np.ndarray, u_cmd: float = 0.0, h: np.ndarray | None = None) -> np.ndarray:
    """Numerical Jacobian of the state equations at y (central differences)."""
    n = len(y)
    h = h if h is not None else np.array([1e-6, 1e-6, 1e3, 1e3, 1e3, 1e3, 1e-6])
    J = np.zeros((n, n))
    for i in range(n):
        yp, ym = y.copy(), y.copy()
        yp[i] += h[i]
        ym[i] -= h[i]
        J[:, i] = (np.array(pl.rhs(0, yp, u_cmd)) - np.array(pl.rhs(0, ym, u_cmd))) / (2 * h[i])
    return J


def modal(J: np.ndarray) -> tuple[float, float]:
    """Lowest oscillatory pole pair: (natural frequency rad/s, damping ratio)."""
    ev = np.linalg.eigvals(J)
    osc = [e for e in ev if e.imag > 1e-3]
    e = min(osc, key=lambda z: abs(z))
    wn = abs(e)
    return wn, -e.real / wn


def frequency_vs_stroke(des: Design, eps_list=(0.0, 0.005, 0.02), n: int = 41) -> dict:
    out = {"x": np.linspace(0.02, 0.98, n) * des.stroke}
    for eps in eps_list:
        pl = Plant(des, eps0=eps)
        ws = []
        for x in out["x"]:
            y = pl.hold_state(x)
            ws.append(omega_formula(pl, x, y[2], y[3]) / (2 * math.pi))
        out[eps] = np.array(ws)
    return out


def loop_gain_limit(pl: Plant, x: float, Kp_max: float = 200.0) -> dict:
    """Position loop u = Kp (x_ref - x) around the model at x, linearised. Returns the largest stable Kp,
    the corresponding velocity gain Kv_loop = Kp q_gain / A_A (1/s), and Merritt's estimate 2 zeta omega."""
    d = pl.des
    pl.coulomb_on = False
    y = pl.hold_state(x)
    J0 = linearise(pl, y, 0.0)
    # steady velocity per unit spool opening of the linear model (x held, all other states settled)
    z = np.linalg.solve(J0[1:6, 1:6], -J0[1:6, 6])
    v_per_u = abs(z[0])
    wn, zeta = modal(J0)

    def stable(Kp):
        J = J0.copy()
        # du/dt = (Kp (xr - x) - u)/tau -> d(du)/dx = -Kp/tau
        J[6, 0] += -Kp / d.g["valve_tau_s"]
        ev = np.linalg.eigvals(J)
        # the oscillatory (hydraulic) pair decides; the real modes near zero are numerical integrators
        return np.max(ev[np.abs(ev.imag) > 1.0].real) < 0

    lo, hi = 1e-4, Kp_max
    if stable(hi):
        kp = hi
    else:
        for _ in range(50):
            mid = math.sqrt(lo * hi)
            lo, hi = (mid, hi) if stable(mid) else (lo, mid)
        kp = lo
    return {"wn": wn, "zeta": zeta, "Kp_max": kp, "Kv_loop_max": kp * v_per_u, "merritt_2zw": 2 * zeta * wn}


# ---------------------------------------------------------------- cushion study
def cushion_case(des: Design, side: str, T_C: float = 40.0, needle_scale: float = 1.0, cushion_on: bool = True) -> dict:
    """Run into one end at full speed and record the cushion chamber pressure and speed."""
    pl = Plant(des, T_C=T_C, needle_scale=needle_scale, cushion_on=cushion_on)
    v = des.o["speed_max_mm_s"] * 1e-3
    s = des.stroke
    if side == "B":  # extending into the rod-end cushion
        x0 = s - pl.lc - 0.25 * s
        u = 1.0
    else:
        x0 = pl.lc + 0.25 * s
        u = -des.g.get("retract_command", 1.0)
    y0 = pl.hold_state(x0)
    t_end = 0.25 * s / v + 1.5
    r = simulate(pl, lambda t: u * min(1.0, t / 0.15), t_end, y0, 5e-4)
    pc = r["pB"] if side == "B" else r["pA"]
    eng = (r["x"] > s - pl.lc) if side == "B" else (r["x"] < pl.lc)
    i_in = int(np.argmax(eng)) if eng.any() else len(r["t"]) - 1
    v_in = abs(r["v"][i_in])
    # impact speed: speed when the end stop is first touched
    hit = (r["x"] >= s) if side == "B" else (r["x"] <= 0)
    # impact speed = speed at the last sample before the stop is touched (the stop is stiff, so the
    # first sample after contact is already decelerating)
    ih = int(np.argmax(hit))
    v_hit = float(abs(r["v"][max(ih - 1, 0)])) if hit.any() else 0.0
    return {
        "side": side,
        "T_C": T_C,
        "needle_scale": needle_scale,
        "p_peak_bar": float(np.max(pc)) / 1e5,
        "v_entry": v_in,
        "v_impact": v_hit,
        "energy_kJ": 0.5 * des.mass * v_in**2 / 1e3,
        "r": r,
    }
