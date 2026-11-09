# Hydraulic cylinder: from order specification to verified design to acceptance test

**Systems engineering dossier.** Structure after C. S. Wasson, *System Engineering Analysis, Design, and Development*, 2nd ed. (Wiley, 2016); page numbers are the book's printed pages. All engineering numbers come from `python -m cylver report` (reports in `reports/`, figures in `docs/img/`). Section 12 lists what is assumed. Section 13 lists every source and whether it was checked.

Author: Chander Vinayak Sharma, 9 Nov 2026. Status: analysis and simulation only. No hardware was built or tested.

---

## 0. The question and the answer

**Question.** A customer orders a custom double-acting cylinder: bore, rod, stroke, pressure, mounting, duty. Which design is verified to be safe and to last, and where are the real failure margins?

**Answer in one paragraph.** For both example orders, a design sized by the usual static rules (Lamé wall, Euler rod) passes every static check and still fails in service. The failure margins are not where the static rules look:

1. **Fatigue at the tube-to-base weld.** The flat base stops the tube from expanding, which bends the wall next to it. At the original weld position, 10 mm from the base face, this raises the longitudinal stress to 2.2 × the membrane value (HD-160).
    - Rev A fatigue stress factor: **0.38** for the required 4 × 10^5 duty cycles.
    - Moving the girth weld to the point where the edge moment changes sign (βx = π/4, 19 mm from the base face) gives **1.09** with the same tube.
    - Thickening the tube instead works poorly. In a design-exploration run (not in the report), a 25 mm wall, 2.2 × thicker, still gave only 0.81 with the weld left at 10 mm.
2. **The cushion chamber is not protected by the port relief valve.** With an overrunning load and cold oil, the cap cushion reached **446 bar** against a 280 bar relief and a 420 bar proof pressure (HD-160 rev A).
    - That peak is also what drives the weld and port fatigue damage.
    - It is fixed at system level: an interface requirement limits the retract command to 0.55 because the load overruns on retraction, and the needle opens to 8 mm².
3. **Buckling of the long loader cylinder (LD-80).** The ordered 50 mm rod reaches first yield at **0.96 ×** the push at relief, once the rod tilt allowed by the bearing clearances, bow, pin eccentricity and self weight are included.
    - The rod-only Euler shortcut says 0.82×; the elastic stepped-column load says 1.28×.
    - A 56 mm rod cannot reach the 1.5× target at any stop-tube length: it peaks at 1.33 near 100 mm and falls again for longer tubes.
    - A 63 mm rod with a 50 mm stop tube reaches **1.65×**.
4. **One order requirement cannot be met as written.** Drift ≤ 5 mm in 10 min (HD-160): oil cooling of only 4.4 K in the locked cap chamber already uses the whole allowance. It goes back to the customer as a specification issue (Wasson Principle 19.4, p. 400: a specification must bound a *feasible* solution space).
5. **Seal friction is calibrated to published measurements**, not assumed (section 7). A first version of this dossier used assumed seal data and concluded that smooth motion at 10 mm/s was infeasible (stick-slip onset 12 mm/s). The measured data put the onset at 1.8 mm/s (PTFE system) and 2.3 mm/s (PU), so R-12 is met. That reversal is the reason to calibrate before concluding.

The verified designs meet every other requirement: HD-160 rev B misses only the drift specification issue, and LD-80 rev C meets all fourteen. They also show **one design method for both orders**:

- move the base weld to the zero-moment point;
- put the ports in the end blocks, not through the tube wall;
- limit the retract command, because the load overruns;
- set each cushion needle for its own end.

LD-80 rev C also takes the PTFE seal system that its own trade study picked in 100 % of weight sets.

The parts are handed to the CNC part-family generator (project 07), which generates and verifies the lathe programs. For LD-80 rev C, 07 accepts and verifies all three parts. For HD-160 rev B, it rejects the gland: its internal grooving tool reaches 6 mm and the wiper groove is 6.15 mm deep.

---

## 1. Need, problem space and concept of operations (Wasson ch. 4-6)

**Problem statement** (Principle 4.15, p. 91). A make-to-order cylinder maker turns each order line into a design, a set of drawings, NC programs and a test record, often within days. Whether the design is right is usually settled by experience, by catalogue tables and by a final pressure test. The pressure test is ISO 10100 style: leakage, proof, optionally friction (ISO 10100:2020 contents, section 13). That test cannot reveal fatigue, buckling margins or cushion peaks at the customer's duty. Those failures appear in the field, months later.

**Stakeholders and their needs** (Principle 4.12, p. 89, operational needs):

| Stakeholder | Need |
|---|---|
| Machine builder (customer) | force and speed at the duty; no field failures; known natural frequency for its control |
| Operator | smooth, predictable motion; load stays where it was put |
| Maintenance | known wear parts, reseal interval, no surprises |
| Cylinder maker: design | a repeatable route from order line to verified design |
| Cylinder maker: production | part-family dimensions that the lathes and their programs can make |
| Cylinder maker: test bay | pass limits that mean something, traceable to the design |

**Use cases** (Principles 5.12-5.23, pp. 120-122). Each use case has a trigger, a main success scenario and the scenarios that go wrong; the last column is the analysis that covers them.

| UC | Title (outcome) | Trigger | Main scenario | When it goes wrong | Covered by |
|---|---|---|---|---|---|
| UC-1 | Lift the load | operator commands extend | extend at v_max against the load, stop at mid stroke, continue | stall, pressure peaks at the stop | duty simulation |
| UC-2 | Hold the load | valve centred | load stays within the drift limit | seal leakage, oil cooling | drift analysis |
| UC-3 | Lower the load | operator commands retract | controlled retraction | **load overruns**, high speed, cushion peak | duty + cushion simulation |
| UC-4 | Reach the end of stroke | piston enters cushion | decelerate to a gentle stop | over-pressure in the cushion chamber, hard impact | cushion study |
| UC-5 | Position slowly | fine command | smooth slow motion | stick-slip | friction model |
| UC-6 | Accept the cylinder | finished unit on the test bench | leakage, proof, friction tests pass | limits not traceable to the design | ATP (section 9) |
| UC-7 | Manufacture the parts | released design | NC programs generated and verified | part outside tool limits | CNC hand-off (section 10) |

**Concept of operations** (Principle 6.6, p. 139). The order file (`data/orders/*.json`) is the single source of the order (Principle 27.4, p. 579, single information source). `cylver` sizes the tube and base, runs every analysis, derives the requirements status, the FMEA, the acceptance limits and the part-family rows, and writes a calculation report. A design change is a new revision in the same order file, with its reason. The report shows the revision history. Each revision is a decision with its evidence.

---

## 2. Phases, modes and states (Wasson ch. 7)

Wasson separates *phases* of the life cycle (Section 7.3, p. 150), *modes* the user selects (Principle 7.4, p. 157) and *states* (operational, environmental and dynamic; Principles 7.3, 7.8 and 7.9, pp. 156-166). Each mode and state that loads the cylinder becomes a **load case**: modes and states are where the load cases come from.

| Kind | Name | What it does to the cylinder | Load case in `cylver` |
|---|---|---|---|
| phase | manufacture | weld, machine, assemble | weld detail category, CNC hand-off |
| phase | acceptance test | proof pressure, leakage, friction | proof stress check, ATP |
| phase | operation | the duty cycle | duty simulation |
| phase | maintenance | reseal; clearances grow with wear | clearance tilt (buckling), extrusion gap |
| mode | extend | cap side pressurised, push | buckling at full extension; duty |
| mode | hold | both ports blocked | drift; oil cooling |
| mode | retract | rod side pressurised; load may overrun | cushion at the cap end; duty |
| mode | float / off | no pressure | (not a load case) |
| dynamic state | cushioning | trapped oil decelerates the mass | cushion study, 40 °C and 0 °C oil |
| dynamic state | stop at mid stroke | valve closes in 50 ms | duty pressure peaks |
| dynamic state | stall at relief | pressure at the relief setting | static checks at relief pressure |
| environmental state | cold oil (0 °C) | viscosity about 10 × higher | cushion peak, friction |
| environmental state | aerated oil | effective bulk modulus drops | natural frequency |

**Allowable and prohibited actions** (Principle 7.10, p. 169). The analysis produced one prohibited action that only the machine can enforce: *full-speed retraction under an overrunning load.* It becomes the interface requirement IR-1 in section 4.

---

## 3. Requirements derived from the order and the standards (Wasson ch. 19-22)

**How they are derived** (Section 21.4, p. 436). Every requirement comes from one of three places:

- an order clause (force, speed, pressure, duty, drift, smooth speed);
- a standard's method (ISO 10100 acceptance modules, the ISO/TS 13725 approach to buckling, EN 1993-1-9 weld fatigue);
- a derived need found by analysis, for example "the cushion chamber is not protected by the relief, so its peak needs its own requirement".

Each requirement follows Principle 22.5 (p. 449): subject, *shall*, action, outcome, bound. Each names its verification method (Principle 22.4, p. 448). The methods are Wasson's Inspection, Examination, Analysis, Demonstration and Test (Principle 13.10, p. 280). Where the order is silent on conditions (oil temperature, the thermal state during holding), Principle 19.20 (p. 407) on environmental conditions applies and the gap is reported.

| ID | Requirement (short) | Source | Method | HD-160 rev B | LD-80 rev C |
|---|---|---|---|---|---|
| R-01 | extend at ≥ 0.95 v_max against the load | order | A; T | met | met |
| R-02 | tube von Mises ≤ Re/2 at relief (Pm) | Lamé; design rule | A | met | met |
| R-03 | no yield at proof pressure | ISO 10100 §8.3 (factor assumed 1.5) | A; T | met | met |
| R-04 | base plate ≤ 1.5 Re/2 (Pb, conservative bound) | plate theory | A | met | met |
| R-05 | base junction ≤ 2 Re (PL+Pb+Q, shakedown) | ASME VIII-2 Part 5 categories | A | met | met |
| R-06 | rod thread ≤ Re/2 at relief pull | ISO 898-1 | A | met | met |
| R-07 | first yield at full extension ≥ 1.5 × push at relief | ISO/TS 13725 approach; 1.5 = project target | A | met (1.65) | met (1.65) |
| R-08 | every fatigue location survives the required duty cycles | order duty; S-N sources | A; T (ISO 10771-1) | met (1.09, weld) | met (1.15, weld) |
| R-09 | piston-seal extrusion gap ≤ allowable | seal makers' gap charts (values assumed, tables not retrieved) | A; I | met | met |
| R-10 | cushion peak ≤ proof pressure, 0-40 °C | derived (relief does not protect the cushion) | A; T | met (324 / 420 bar) | met (273 / 330 bar) |
| R-11 | end-stop impact ≤ 0.10 m/s | project target | A; T | met (0.080) | met (0.090) |
| R-12 | no stick-slip down to v_min | order; seal data from Sârbu et al. 2022/2024 | A; T (ISO 10100 module F) | met (onset 1.8 mm/s) | met (onset 1.8 mm/s) |
| R-13 | drift ≤ 5 mm / 10 min | order | A; T (ISO 10100 module P) | **not met** (4.4 K cooling uses it up) | met only below 5.9 K cooling |
| R-14 | lowest natural frequency reported | TPM (Section 21.8, p. 442) | A | 11.0 Hz | 11.1 Hz |

The full statements, results and limits for every revision are in `reports/HD-160.md` and `reports/LD-80.md`. Traceability runs both ways (Principle 19.13, p. 402; Section 21.7, p. 439):

- each requirement names its source;
- the report gives the analysis that verified it;
- the acceptance procedure (section 9) names the requirement each test limit protects.

**Specification issues found** (Principle 19.4, p. 400, feasibility):

- **R-12 (withdrawn as an issue).** With assumed seal data the model put the stick-slip onset at 12 mm/s, above the 10 mm/s requirement. Calibrated to measured friction (section 7), the onset is 1.8 mm/s, so the requirement is met. It stays verified by test (F-2), because seal friction varies between batches and with wear.
- **R-13.** Drift is defined without a thermal condition. A locked 32 L cap chamber that cools by 4.4 K in 10 minutes uses the whole 5 mm allowance (oil expansion coefficient assumed 7 × 10^-4 /K). The order must say "at constant oil temperature" or give the cooling rate.

---

## 4. Architecture and interfaces (Wasson ch. 26-27)

**Physical architecture** (Principle 26.5, p. 551):

```
rear pin ─ base (welded) ─ tube ─ head/gland (bearing, rod seal, wiper) ─ rod ─ rod eye ─ front pin
                │                         │
            cap cushion             rod cushion spear
            port A (tube or block)  port B (tube or block)
          piston (guide rings, piston seal) on rod thread + nut
```

**Interfaces** (Principles 27.1, p. 575; 27.6, p. 583). Each interface has an owner (Principle 27.3, p. 577) and becomes either a requirement or a load case.

| ID | Interface | Across it | Owner | Becomes |
|---|---|---|---|---|
| IF-1 | rear and front pins | axial force; eccentricity ≤ 1 mm; no moment | machine builder | buckling eccentricity, pin sizes |
| IF-2 | ports A/B ↔ hoses ↔ valve | flow, pressure; hose compliance | machine builder | dynamic model; port location |
| IF-3 | port relief valves | relief setting | machine builder | static checks at relief; *not* the cushion |
| IF-4 | valve command | command profile | machine builder (control) | **IR-1: the retract command shall be limited so the overrunning load does not exceed v_max (0.55 for HD-160, 0.6 for LD-80)** |
| IF-5 | oil | temperature 0-70 °C, air content | machine builder | viscosity, bulk modulus |
| IF-6 | design → production | part-family rows | cylinder maker | CNC hand-off (section 10) |
| IF-7 | design → test bay | pass limits | cylinder maker | ATP (section 9) |

IF-3 and IF-4 are where the HD-160 rev A cushion failure came from: the relief valve on the machine side does not protect a volume inside the cylinder. FMEA is applied to interfaces as well as parts (Principle 27.9, p. 596).

**Software architecture of the verification tool** (`cylver/`). One order file goes in. Each module owns one physics question, and all of them share the `Design` object. The dynamic model is the one source of pressure histories, so cushion peaks, fatigue spectra, natural frequency and stick-slip all come from the same equations.

| Module | Question | Feeds |
|---|---|---|
| `model.py` | geometry, sizing of tube and base, materials | all |
| `stress.py` | Lamé, plate, edge bending, thread | requirements, fatigue |
| `buckling.py` | stepped beam-column FE, P_cr, first yield with imperfections | R-07, FMEA |
| `dynamics.py` | valve + hoses + chambers + cushions + friction + stops | duty history, cushion, frequency, loop gain |
| `fatigue.py` | rainflow, S-N, Miner, equal-damage test pressure | R-08, ATP Q-1 |
| `seals.py` | Stribeck friction, stick-slip, drift | R-12, R-13, ATP F/P |
| `requirements.py`, `fmea.py`, `acceptance.py`, `aoa.py` | systems-engineering layers | report |
| `cnc_bridge.py` | part-family rows → project 07 | manufacturing |

---

## 5. Budgets and margins (Wasson ch. 31)

A design safety margin is "an additive or multiplicative factor that is applied to a Measure of Performance… to accommodate variations and uncertainty" (Section 31.2, p. 667). Every capability needs a budget and a margin (Principle 31.1, p. 667). The margins used here, and where each came from:

| Margin | Value | Basis |
|---|---|---|
| static primary membrane | Re / 2 | project design rule (common practice); ASME VIII-2 categories for the split |
| static primary bending | 1.5 Re / 2 | rectangular-section shape factor 1.5 |
| primary + secondary (junction) | 2 Re | shakedown (ASME VIII-2 Part 5), from memory |
| buckling | first yield ≥ 1.5 × push at **relief** | project target, applied to a *usable* load, not to Euler |
| fatigue, welds | γ_Mf = 1.35 on stress range | EN 1993-1-9 safe-life, high consequence |
| fatigue, machined | j = 1.5 on stress | FKM order of magnitude |
| cushion | peak ≤ proof pressure | derived: the proof test is the highest pressure the unit has seen without damage |

**Pressure budget.** Supply ≤ relief < proof, and the cushion peak must stay ≤ proof. HD-160: 250 / 280 / 420 bar, with the cushion peak 446 bar in rev A and 324 bar in rev B.

**Drift budget**, apportioned (Principle 31.4, p. 668): 5 mm = thermal contraction + piston-seal leakage over 10 min. The thermal share is computed and the remainder is the leak budget, which becomes the ISO 10100 module P pass limit. For HD-160 the remainder is zero: that is how R-13 was found.

**Reliability index for buckling** (stress-strength interference, Principle 34.21, p. 752). Re CoV 7 % and load CoV 10 % (both assumed) give β = 4.25 for HD-160 rev B, about 1 × 10^-5 per extreme push.

---

## 6. Analysis of alternatives (Wasson ch. 32)

Method after Wasson Fig. 32.7 (p. 690):

1. Criteria with weights.
2. Utility functions (linear, 0-10).
3. Two to four viable candidates.
4. Scoring.
5. Ranking.
6. **Sensitivity**: each weight ±10 % (Wasson's definition, p. 682) and 2000 random weight sets.

Each score is marked **model** (computed by `cylver`) or **judgement** (rationale written next to it). An undocumented trade study "is nothing more than personal opinion" (Principle 32.5, p. 697), so the judgement scores stay visible and editable in `cylver/aoa.py`.

| Trade | Result (HD-160 rev B) | Robust? | What decides it |
|---|---|---|---|
| Seal system | PTFE step seal + PU tandem, PTFE piston seal | HD-160: the PU system is **not viable** (extrusion margin 0.90 < 1, R-09), so it is screened out before scoring (Wasson's viable-alternative rule, p. 683). LD-80: PTFE wins 81 % of weight sets; no ±10 % weight change flips it | **model (measured-data seals):** stick-slip onset 1.8 vs 2.3 mm/s, so both are smooth at 10 mm/s; lower friction and a larger extrusion margin for PTFE. These outweigh the judgement criteria (dry rod, contamination, cost) that favour PU. |
| Rod surface | HVOF WC-CoCr (65 % of weight sets) ahead of laser cladding (24 %); hard chrome last (6 %) | no ±10 % flip, but **all scores are judgement**, so the ranking is only as good as the scores | the regulatory criterion: Cr(VI) restriction dossier submitted 11 Apr 2025, opinions pending (news source, section 13). Coupon tests (salt spray ISO 9227, abrasion) would replace the judgement scores. |
| Tube material | E355 +SR | wins 100 % of weight sets | **model:** the stronger 25CrMo4 tube is sized thinner (7.5 mm instead of 11.25 mm), so its *weld* fatigue factor drops from 1.09 to 0.53. Higher strength does not help a weld. |

The tube trade is the useful lesson. A static-strength argument picks the stronger steel. The fatigue model picks the stock tube, because the governing location is a weld and weld S-N curves do not improve with parent-metal strength (EN 1993-1-9 detail categories are independent of grade).

---

## 7. Modelling and simulation (Wasson ch. 33)

"All models are wrong but some are useful" (Principle 33.5, p. 718). The fidelity of each model is stated (Principle 33.1, p. 706), and each was **verified** against an independent result. None is **validated** against hardware yet. The acceptance tests in section 9 are designed to do that.

| Model | Fidelity | Verified against (test in `tests/test_cylver.py`) | Not modelled |
|---|---|---|---|
| Lamé, plates | closed form, linear elastic | thin-wall limit; radial equilibrium | plasticity, autofrettage |
| base-junction edge bending | thin-shell beam-on-elastic-foundation | closed-form edge moment; zero at βx = π/4 | base rotation (real base is not rigid), weld geometry; rm/t 4-9 is near the theory's limit |
| buckling | Euler-Bernoulli FE, 160 elements, geometric stiffness, P-δ with imperfections | Euler for 4 end conditions (0.2 %); exact two-segment equation (0.2 %); secant formula and 1/(1-P/Pcr) amplification (1 %) | the ISO/TS 13725 equations themselves; shear; spherical-bearing friction; rod-in-tube contact beyond the two bearings |
| rainflow / Miner | ASTM E1049 three-point | the standard's worked example, exact counts | sequence effects; mean stress for welds (deliberately ignored, as in EN 1993-1-9) |
| dynamic model | lumped: 2 chambers + 2 lines + valve lag, orifice edges, cushions, Stribeck friction, stiff stops | natural frequency = oil-column formula (2 %); loop-gain limit converges to Merritt's K_v < 2ζω_h as the valve gets fast (to 0.90-0.95 ×, the remainder from the asymmetric areas); cushion monotonic in needle size | line dynamics (waves), pump flow limit, valve spool forces, temperature rise, cavitation details, dissolved-air release |
| bulk modulus | Merritt series form, free air by gas law | liquid limit; air softens at low pressure | air dissolving at high pressure (Gholizadeh et al. 2011): conservative |
| seal friction | Stribeck per seal element (Coulomb level growing with pressure, static excess, Stribeck speed, viscous); tandem rod seal = two elements | **calibrated** to Sârbu et al. (Polymers 2024, Tables 3 and 5; 60 mm PTFE and PU coaxial seals, 0-160 bar). Predicted stick-slip onset 0.52 mm/s (PTFE) and 2.3 mm/s (PU) against the measured 0.2-0.3 and 1.4 mm/s (Polymers 2022): conservative by about 2 ×, consistent with the damping the model leaves out | the viscous coefficient (assumed; matters only at speed); PU coaxial data used for a PU U-cup; batch scatter |
| stick-slip | 1-DOF, constant-flow drive, Stribeck | linear criterion dF/dv < 0 agrees with the simulation (stick-slip at 0.5 × onset, smooth at 2 × onset) | frictional memory and lag |

**Two model-based findings that only the simulation could give:**

1. **Overrunning retraction.** With a critical-centre valve, the overrunning load sets the retraction speed:
    - HD-160 retracted at 0.41 m/s against a 0.20 m/s order speed.
    - LD-80 retracted at 0.52 m/s against 0.30 m/s.

    That speed is what overloads the cap cushion: 446 and 463 bar with 0 °C oil. Limiting the command (IR-1) brings the entry speed to 0.21 and 0.28 m/s and the peaks to 324 and 273 bar.
2. **Natural frequency and damping.** The oil-column frequency at mid stroke is 11.7 Hz (HD-160) and 12.5 Hz (LD-80 rev C). The formula and the model's eigenvalue agree to 0.1 Hz. The lowest values over the stroke are 11.0 and 11.1 Hz.
    - With a critical-centre valve and tight seals, the damping ratio at null is only 0.0014 (HD-160) and 0.0028 (LD-80). A position loop could then use a velocity gain of only 0.2-0.4 1/s (Merritt's 2ζω_h, which the model reproduces for a fast valve).
    - With the 30 ms proportional valve assumed here, the model allows 1.2-2.5 1/s, because the slow valve attenuates the loop at the resonance. That is a fragile kind of stability, not a design margin.
    - Anyone closing a loop around this cylinder needs added damping: a bleed orifice, pressure feedback, or a valve with underlap. This goes to the machine builder as a TPM (R-14).

---

## 8. Reliability and FMEA (Wasson ch. 34)

FMEA with RPN = S × O × D (Wasson eq. 34.41, p. 765, after Dhillon). Occurrence is **read from the computed margin** with a fixed rule (`cylver/fmea.py`), so the ranking moves when the design changes. Severity and detection are judgement, written per row. Ranking for HD-160 rev A, the "as ordered" design:

| Failure mode | Margin | O | RPN |
|---|---|---|---|
| tube-to-base weld fatigue | 0.38 | 9 | 630 |
| port cross-bore fatigue (ports in the tube wall) | 0.33 | 9 | 441 |
| drift under load | 0.89 | 9 | 315 |
| cushion over-pressure | 0.94 | 9 | 288 |
| column buckling | 1.65 | 3 | 180 |
| piston-nut thread fatigue | 6.5 | 2 | 160 |
| piston-seal extrusion | 1.41 | 5 | 150 |

For rev B, the occurrence drops for three rows. The weld goes from O = 9 to 7 (margin 1.09), the port from 9 to 2 (margin 3.25) and the cushion from 9 to 5 (margin 1.30). The top of the rev B ranking is:

| Failure mode | Margin | O | RPN |
|---|---|---|---|
| base weld fatigue | 1.09 | 7 | 490 |
| drift | 0.89 | 9 | 315 |
| buckling | 1.65 | 3 | 180 |
| cushion over-pressure | 1.30 | 5 | 160 |

The weld stays at the top because its margin is thin, and that margin depends on an inspected root (FAT 71). With an uninspected root (FAT 36, IIW) the same weld would fail. **Weld-root inspection is therefore a derived manufacturing requirement**, not a quality nicety.

Other reliability points:

- The strength-stress view (Principle 34.21, p. 752) is applied to buckling (β above).
- Wear-out (Principle 34.20, p. 751) enters through the clearances: the buckling tilt and the extrusion gap both grow with wear. The inspection step I-1 records them so that a reseal can restore the analysed state.

---

## 9. Verification and validation (Wasson ch. 13 and 28)

Verification asks "did we build it right?" and validation asks "did we build the right thing?" (Principles 13.3 and 13.4, p. 275). The verification methods are listed per requirement (section 3). The acceptance test procedure below applies three principles:

- one test verifies as many requirements as practical (Principle 28.3, p. 608);
- the test equipment is calibrated (Principle 28.9, p. 618);
- **every pass limit is a number the design analysis relied on**, so a test that passes also validates the model.

| Step | Module | Test | Pass limit (HD-160 rev B; values in the report) | Protects |
|---|---|---|---|---|
| I-1 | inspection | piston and gland bearing clearances | ≤ design values | R-07 (tilt), R-09 (gap) |
| L-1 | ISO 10100 §8.2 | leakage at low test pressure | no external leakage | standard module |
| L-2 | ISO 10100 §8.3 | proof / external leakage | no leakage, no permanent set; analysis: 314 MPa < 450 MPa at 420 bar | R-03 |
| P-1 | ISO 10100 §9 | piston-seal leakage | leak budget from R-13 (zero for HD-160: open item) | R-13 |
| F-1 | ISO 10100 §10.6.2 | static friction | ≤ 1.3 × model breakaway | R-01, R-12 |
| F-2 | ISO 10100 §10.4.2 | friction at constant speed 0.1-100 mm/s, fit the Stribeck model | stick-slip onset of the fitted curve ≤ 10 mm/s | R-12; validates the seal model |
| C-1 | project test | cushion at 0 °C and 40 °C, transducer **in the cushion chamber** | ≤ proof; impact ≤ 0.10 m/s | R-10, R-11 |
| Q-1 | ISO 10771-1 Annex B (type test) | cyclic pressure 0 → p_test | p_test = equal-damage pressure for 10^6 cycles (report) | R-08 |

**What is and is not checked in the standards.** ISO 10100's contents (modules, clause numbers) were checked. Its test pressures and durations were not available, so the procedure says "per ISO 10100" instead of inventing them.

The ISO 10771-1 scope says the user sets the test pressure, so this project computes it: the constant-amplitude pressure that does the same Miner damage at the governing location as the design life.

- For HD-160 rev B this is **183 bar for 10^6 cycles**. For LD-80 rev C it is **150 bar**.
- Both are *below* the rated pressure, because most of the duty damage comes from a few peaks above the S-N knee and the rest falls under the cut-off.
- So a type test at rated pressure for 10^6 cycles is more severe than the design life at this location. Passing it shows the life with margin; failing it would not by itself prove the design wrong.
- For the rev A designs the equal-damage pressure is undefined. The port peak exceeds static strength with the elastic Kt, and the report says "n/a, redesign first" instead of printing a number.

**The evaluator** (`acceptance.evaluate`) fits the Stribeck model to module F data and judges F-1 and F-2. It is demonstrated on **synthetic** data (model + noise), clearly labelled; it shows the code works, not that a cylinder passes. A test with a "sticky" seal set (static/Coulomb × 3) is caught by F-2.

Integration is incremental (Principle 28.1, p. 605): inspection first, then pressure tests, then function tests, then the type test. A failure is isolated by the step it appears in (Principle 28.2, p. 606).

---

## 10. Manufacturing hand-off (the CNC part-family generator)

The verified design produces the rows that project 07 (`cnc-part-family-generator`) consumes: rod end, piston and gland. `cylver/cnc_bridge.py` imports 07 unchanged, generates the FANUC canned-cycle program and runs 07's material-removal verifier on it. The part-family table is derived from the design, so a design revision regenerates the NC programs too.

| Part | HD-160 rev B | LD-80 rev C |
|---|---|---|
| rod end | PASS | PASS |
| piston | PASS | PASS |
| gland | rejected: wiper groove 6.15 mm > 6 mm internal-grooving tool limit | PASS |

The rejection is useful feedback. Either a deeper internal grooving bar goes into the tool library, or the large-rod gland gets a different seal cross-section.

---

## 11. Results summary

| | HD-160 rev A (as ordered) | HD-160 rev B | LD-80 rev A (as ordered) | LD-80 rev B | LD-80 rev C |
|---|---|---|---|---|---|
| tube OD / rod (mm) | 182.5 / 100 | 182.5 / 100 | 90 / 50 | 90 / 63 | 90 / 63 |
| buckling (first yield / push at relief) | 1.65 | 1.65 | **0.96** | 1.65 | 1.65 |
| governing fatigue stress factor | **0.33** (port in tube) | 1.09 (weld) | **0.35** (port in tube) | **0.39** (port in tube) | 1.15 (weld) |
| cushion peak (bar), 0-40 °C | **446** | 324 | **463** | **411** | 273 |
| requirements not met | R-08, R-10, R-13 | R-13 | R-07, R-08, R-10 | R-08, R-10, R-11 | none |
| CNC hand-off (project 07) | - | gland rejected (tool depth) | - | - | 3 / 3 PASS |

---

## 12. Limits and assumptions (what is assumed vs sourced)

- **Not validated on hardware.** Every result is analysis and simulation. The acceptance procedure is the plan to validate the models.
- **Seal friction is calibrated to one published test series** (two 60 mm coaxial seals, PTFE and PU, 0-160 bar, Sârbu et al. 2022/2024). The viscous coefficients are still assumed, and the PU U-cup uses PU coaxial-seal data. **The extrusion-gap limits are still assumed:** the seal makers' tables could not be retrieved. Module F on the real cylinder would replace all of these.
- **ISO/TS 13725 equations were not available.** The buckling model follows the standard's published principles (geometric model, off-axis load, own weight, fully extended) with a classical beam-column FE. The first-yield criterion and the 1.5 target are this project's choices.
- **Material values are from memory or secondary summaries** (marked `checked: false` in `data/materials.json`). E355 +SR uses 450 / 580 MPa from a mill catalogue summary; a recalled 520 / 640 was not confirmed and is not used.
- **The edge-bending solution assumes a rigid base** and thin-shell theory. At rm/t = 4-9 that is at the edge of its validity. A real base rotates a little, which lowers the bending. The zero-moment weld position should be confirmed by an axisymmetric FE model or strain gauges on the type-test unit.
- **The port Kt of 3** is a typical cross-bore value (assumed). Ports in end blocks are modelled as Kt × pressure. Peaks above Rm with the elastic Kt are reported as static failures, not given a fatigue life.
- **The duty cycle is a design duty cycle** (extend with a mid-stroke stop, hold, retract, hold) invented for the example. The real spectrum must come from the customer or from measurement.
- **The pump is an ideal pressure source.** Line wave dynamics are not modelled, and air is treated as free (conservative at high pressure).
- **Proof pressure factor 1.5** is assumed until checked in ISO 10100 (a vendor guide, S25, gives the same typical value; the ISO text is paid and was not available).
- **The groove sizes in the CNC rows are illustrative**, not seal-catalogue values.

---

## 13. References and how far each was checked

"Checked" means the cited text was read in this project (page or section stated). "Not re-checked" means the formula or value is quoted from memory or from a search summary. Every not-re-checked item is covered by a test against an independent closed form, or is listed in section 12 as an assumption.

| # | Source | Used for | Status |
|---|---|---|---|
| W | C. S. Wasson, *System Engineering Analysis, Design, and Development*, 2nd ed., Wiley 2016 | dossier structure; the principles cited per section (book page numbers) | checked: chapters 4-7, 13, 19-22, 26-28, 31-34 extracted and the cited principles read |
| S1 | ISO 10100:2020 *Hydraulic fluid power - Cylinders - Acceptance tests* | ATP structure: modules L (8.2, 8.3), P (9), F (10.4.1, 10.4.2, 10.6.2) | table of contents checked (DIN preview); test pressures and durations NOT available |
| S2 | ISO 10771-1:2015 *Fatigue pressure testing of metal pressure-containing envelopes* | type test Q-1 | scope checked via store pages: test pressure set by the user, ISO/TR 10771-2 rating method, Annex B = cylinders; full text not available |
| S3 | ISO/TS 13725:2021 *Method for evaluating the buckling load of a hydraulic cylinder* | buckling principles | scope checked via store pages: geometric model, off-axis load, own weight, fully extended; based on Hoblit and NFPA T3.6.37; not for internal rod buckling; spherical-bearing friction not included. The equations of the standard were NOT available and are not reproduced |
| S4 | H. Gholizadeh, R. Burton, G. Schoenau, "Fluid bulk modulus: a literature survey", Int. J. Fluid Power 12(3), 2011, pp. 5-15 | effective bulk modulus, entrained vs dissolved air, critical pressure | checked (paper read) |
| S5 | H. E. Merritt, *Hydraulic Control Systems*, Wiley 1967 | chamber continuity, hydraulic natural frequency, K_v < 2 zeta omega_h | not re-checked; natural frequency and gain rule reproduced by tests against the model's own eigenvalues |
| S6 | M. Jelali, A. Kroll, *Hydraulic Servo-systems*, Springer 2003 | valve-cylinder modelling | not re-checked |
| S7 | R. G. Budynas, J. K. Nisbett, *Shigley's Mechanical Engineering Design*, "pressurized cylinders" | Lame | not re-checked; test: thin-wall limit and radial equilibrium |
| S8 | S. Timoshenko, S. Woinowsky-Krieger, *Theory of Plates and Shells*, 2nd ed. | circular plate; cylinder edge bending (beam on elastic foundation) | not re-checked; test: edge moment closed form and zero-moment point |
| S9 | S. Timoshenko, J. Gere, *Theory of Elastic Stability* | stepped bars, secant formula, amplification 1/(1-P/Pcr) | not re-checked; tests: FE vs Euler (4 end conditions), vs the exact two-segment equation, vs secant formula and amplification |
| S10 | ASTM E1049-85 *Standard Practices for Cycle Counting in Fatigue Analysis* | rainflow | not re-checked; test: the standard's worked example (counts 3:0.5, 4:1.5, 6:0.5, 8:1, 9:0.5) |
| S11 | EN 1993-1-9 *Fatigue* | weld S-N shape, gamma_Mf | not re-checked |
| S12 | IIW recommendations (Hobbacher, 2nd ed. 2016) | FAT 71 / 36 for one-sided butt welds with / without root inspection | not re-checked |
| S13 | FKM guideline *Analytical strength assessment* | machined steel endurance order of magnitude | not re-checked |
| S14 | VDI 2230 Part 1 (2015) | thread endurance sigma_ASV, load factor phi_n | not re-checked |
| S15 | ISO 898-1 | thread stress area | test against the tabulated M20 value (245 mm^2) |
| S16 | B. Armstrong-Helouvry, P. Dupont, C. Canudas de Wit, "A survey of models, analysis tools and compensation methods for the control of machines with friction", Automatica 30(7), 1994 | Stribeck model, stick-slip | not re-checked; test: simulation agrees with the linear criterion |
| S17 | H. K. Mueller, B. S. Nau, *Fluid Sealing Technology* | pressure-dependent seal friction, extrusion | not re-checked |
| S23 | F. A. Sârbu, F. Arnăuţ, A. Deaconescu, T. Deaconescu, "Theoretical and experimental research concerning the friction forces developed in hydraulic cylinder coaxial sealing systems made from polymers", *Polymers* 16(1):157, 2024, doi:10.3390/polym16010157 (CC BY) | seal friction vs pressure and speed (Tables 3, 5) | checked: full text read, tables used |
| S24 | same group, "Experimental research on polymer-based coaxial sealing systems of hydraulic cylinders for small displacement velocities", *Polymers* 14(2):290, 2022, doi:10.3390/polym14020290 (CC BY) | seal size (60 × 4 mm); measured stick-slip onsets (PTFE < 0.2 mm/s, PU < 1.4 mm/s) | checked: full text read |
| S25 | hydraulic-cylinders.net, "Hydraulic cylinder pressure testing: ISO standards and quality control procedures" (vendor guide) | proof pressure "typically 1.5 × the maximum working pressure" | secondary, vendor source; ISO 10100 §8.3 itself not available |
| S18 | EN 10305-1 (E355 +SR tube), EN 10083-2/-3 (C45, 42CrMo4, 25CrMo4) | material values | not re-checked; E355 +SR 450/580 MPa from a mill catalogue seen as a search summary (my recalled 520/640 was NOT confirmed and is not used) |
| S19 | ECHA Cr(VI) restriction process: Annex XV dossier 11 Apr 2025; consultation 18 Jun - 18 Dec 2025; RAC/SEAC opinions pending (Freyr regulatory news, 25 May 2026) | coating trade, regulatory criterion | checked via search summary of the news item; final adoption date unknown |
| S20 | ASME BPVC Section VIII Division 2, Part 5 (design by analysis) | stress categories Pm / Pb / PL+Pb+Q and the 2 Re shakedown limit | not re-checked |
| S21 | ASTM D341 | viscosity-temperature (Walther) | not re-checked; test: VG 46 = 46 cSt at 40 C by construction |
| S22 | ISO 3320 | preferred bores and rods | not re-checked (values only used as a list) |
