"""Order and design data: load an order file, size the tube and base, and derive the geometry.

Internal units are SI (m, Pa, N, kg, s). Order files use mm, bar, kN.
"""

from __future__ import annotations

import copy
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
G = 9.81
BAR = 1e5


def load_materials(path: Path | None = None) -> dict:
    return json.loads(Path(path or DATA / "materials.json").read_text(encoding="utf-8"))


def _merge(a: dict, b: dict) -> dict:
    out = copy.deepcopy(a)
    for k, v in b.items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out


def load_order(path: Path) -> list[tuple[str, dict]]:
    """Return [(variant name, merged order file)]: the base design first, then each variant."""
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    out = [("rev A", raw)]
    for v in raw.get("variants", []):
        merged = _merge(raw, {"order": v.get("order", {}), "design": v.get("design", {})})
        merged["variant_why"] = v.get("why", "")
        out.append((v["name"], merged))
    return out


@dataclass
class Steel:
    name: str
    E: float
    nu: float
    rho: float
    Re: float
    Rm: float
    checked: bool
    source: str


def steel(mats: dict, name: str, d: float | None = None) -> Steel:
    s = mats["steels"][name]
    if "Re_by_d" in s:
        dmm = (d or 0.0) * 1e3
        row = next((r for r in s["Re_by_d"] if dmm <= r[0]), s["Re_by_d"][-1])
        re, rm = row[1], row[2]
    else:
        re, rm = s["Re_MPa"], s["Rm_MPa"]
    return Steel(name, s["E_GPa"] * 1e9, s["nu"], s["rho"], re * 1e6, rm * 1e6, s["checked"], s["source"])


@dataclass
class Oil:
    rho: float
    nu40: float
    nu100: float
    K: float
    alpha_v: float
    kappa: float

    def viscosity(self, T_C: float) -> float:
        """Dynamic viscosity (Pa s) from the Walther form log log(nu+0.7) = A - B log T (ASTM D341)."""
        t1, t2 = 313.15, 373.15
        y1, y2 = math.log10(math.log10(self.nu40 + 0.7)), math.log10(math.log10(self.nu100 + 0.7))
        b = (y1 - y2) / (math.log10(t2) - math.log10(t1))
        a = y1 + b * math.log10(t1)
        nu = 10 ** (10 ** (a - b * math.log10(T_C + 273.15))) - 0.7
        return nu * 1e-6 * self.rho


def oil(mats: dict) -> Oil:
    o = mats["oil"]
    return Oil(o["rho"], o["nu40_cSt"], o["nu100_cSt"], o["K_liquid_GPa"] * 1e9, o["alpha_v_per_K"], o["kappa_air"])


# ISO 3320 preferred bores and rod diameters (mm), quoted from memory, not re-checked.
ISO3320_BORES = [25, 32, 40, 50, 63, 80, 100, 125, 140, 160, 180, 200, 220, 250]
ISO3320_RODS = [12, 14, 16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 63, 70, 80, 90, 100, 110, 125, 140, 160]


@dataclass
class Design:
    """A sized cylinder. Lengths in m, pressures in Pa."""

    oid: str
    variant: str
    title: str
    D: float
    d: float
    stroke: float
    p_rated: float
    p_relief: float
    p_supply: float
    proof_factor: float
    sf_static: float
    tube: Steel
    rod: Steel
    oil: Oil
    ro: float = 0.0
    t_base: float = 0.0
    raw: dict = field(default_factory=dict)

    # ---- order quantities
    @property
    def o(self) -> dict:
        return self.raw["order"]

    @property
    def g(self) -> dict:
        return self.raw["design"]

    @property
    def ri(self) -> float:
        return self.D / 2

    @property
    def A_A(self) -> float:
        return math.pi / 4 * self.D**2

    @property
    def A_r(self) -> float:
        return math.pi / 4 * self.d**2

    @property
    def A_B(self) -> float:
        return self.A_A - self.A_r

    @property
    def wall(self) -> float:
        return self.ro - self.ri

    @property
    def thrust(self) -> float:
        """Push force at rated pressure, rod side vented."""
        return self.p_rated * self.A_A

    @property
    def angle(self) -> float:
        return math.radians(self.o["mounting"]["angle_deg"])

    @property
    def mass(self) -> float:
        return self.o["moving_mass_kg"]

    def mm(self, key: str) -> float:
        return self.g[key] * 1e-3

    # ---- axial layout, measured from the rear pin, cylinder fully extended
    def layout(self, s: float | None = None) -> dict:
        s = self.stroke if s is None else s
        x_base = self.mm("rear_pin_to_base_mm")
        x_pist = x_base + s
        w_p = self.mm("piston_width_mm")
        x_gl_in = x_base + self.stroke + w_p + self.mm("stop_tube_mm")
        x_gl_out = x_gl_in + self.mm("gland_length_mm")
        # the tube is fixed; the rod end (and its pin) moves with the piston
        # rod_end_to_pin_mm = gland outer face to front pin with the cylinder closed
        x_end = x_gl_out + self.mm("rod_end_to_pin_mm") + s
        x_gl_c = x_gl_in + self.mm("gland_bearing_centre_mm")
        guide = x_gl_c - (x_pist + w_p / 2)  # piston bearing centre to gland bearing centre
        return {
            "x_base": x_base,
            "x_piston": x_pist,
            "x_piston_c": x_pist + w_p / 2,
            "x_gland_in": x_gl_in,
            "x_gland_c": x_gl_in + self.mm("gland_bearing_centre_mm"),
            "x_gland_out": x_gl_out,
            "x_end": x_end,
            "guide": guide,
            "L": x_end,
        }

    @property
    def closed_length(self) -> float:
        return self.layout(0.0)["L"]

    @property
    def extended_length(self) -> float:
        return self.layout()["L"]

    # ---- masses per unit length (kg/m) for self weight
    def mass_per_length(self) -> dict:
        steel_t = self.tube.rho * math.pi * (self.ro**2 - self.ri**2)
        steel_r = self.rod.rho * self.A_r
        oil_cap = self.oil.rho * self.A_A
        oil_ann = self.oil.rho * self.A_B
        return {"tube_cap": steel_t + oil_cap, "tube_ann": steel_t + oil_ann, "rod": steel_r}


def size_tube(p_rated: float, ri: float, Re: float, sf: float, proof: float, step: float = 2.5e-3) -> float:
    """Smallest outer radius (on a 2.5 mm OD grid) that keeps von Mises at the bore <= Re/sf at rated
    pressure and <= Re at proof pressure. Closed-end Lame: sigma_vM(bore) = sqrt(3) p ro^2/(ro^2-ri^2)."""
    od = 2 * ri + 2 * step
    while True:
        ro = od / 2
        k = math.sqrt(3) * ro**2 / (ro**2 - ri**2)
        if k * p_rated <= Re / sf and k * p_rated * proof <= Re:
            return ro
        od += step


def size_base(p: float, ri: float, Re: float, sf: float, proof: float, nu: float = 0.3) -> float:
    """Flat base sized on the conservative bound (simply supported plate, centre stress
    3 (3+nu) p a^2 / (8 t^2)) as primary bending (allowable 1.5 Re/sf, rectangular-section shape factor),
    and with no yield at proof pressure; rounded up to 5 mm."""
    c = 3 * (3 + nu) / 8
    t = max(math.sqrt(c * p * ri**2 / (1.5 * Re / sf)), math.sqrt(c * p * proof * ri**2 / (1.5 * Re)))
    return math.ceil(t * 1e3 / 5.0) * 5e-3


def build(raw: dict, variant: str = "rev A", mats: dict | None = None) -> Design:
    mats = mats or load_materials()
    o, g = raw["order"], raw["design"]
    D, d = o["bore_mm"] * 1e-3, o["rod_mm"] * 1e-3
    tube = steel(mats, g["tube_material"])
    rod = steel(mats, g["rod_material"], d)
    des = Design(
        raw["id"],
        variant,
        raw["title"],
        D,
        d,
        o["stroke_mm"] * 1e-3,
        o["p_rated_bar"] * BAR,
        o["p_relief_bar"] * BAR,
        o["p_supply_bar"] * BAR,
        g["proof_factor"],
        g["sf_static"],
        tube,
        rod,
        oil(mats),
        raw=raw,
    )
    # the tube and base are sized at the highest static pressure either chamber can hold: the relief setting
    p_design = max(des.p_rated, des.p_relief)
    des.ro = g.get("tube_od_mm", 0) * 5e-4 or size_tube(p_design, des.ri, tube.Re, des.sf_static, des.proof_factor)
    des.t_base = size_base(p_design, des.ri, tube.Re, des.sf_static, des.proof_factor)
    return des
