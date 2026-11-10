# Open items

Everything that software or simulation could close has been closed (see the decision log in [DOSSIER.md](DOSSIER.md), section 14). The items below need something this project cannot provide: a paid document, a supplier, hardware, or a customer's answer.

| # | Item | Why it is open | What closes it | Affects |
|---|---|---|---|---|
| O-1 | ISO 10100:2020 test pressures and durations (modules L, P, F) | paid standard; only the table of contents was available | buy ISO 10100:2020 and replace the "per ISO 10100" entries and the assumed proof factor 1.5 | R-03, ATP L-1, L-2 |
| O-2 | ISO/TS 13725:2021 equations | paid standard; the buckling model follows its published principles with a classical FE instead | buy ISO/TS 13725 and compare its buckling load with `buckling.py` for both orders | R-07 |
| O-3 | ISO 10771-1:2015 Annex B (cylinder fatigue test details) | paid standard; only the scope was available | buy ISO 10771-1 and fix the type-test set-up around the computed equal-damage pressure | R-08, ATP Q-1 |
| O-4 | ISO 6022 / ISO 6020-2 mounting dimension tables | paid standards | check the order's pin dimensions against the tables (R-15 is "open" until then) | R-15 |
| ~~O-5~~ | ~~Seal makers' allowable extrusion-gap tables~~ | **closed 10 Nov 2026**: Trelleborg installation tables (Glyd Ring, U-Cup PUA, Stepseal 2K) are in `cylver/seals.py` (DOSSIER S26). They made PU viable and moved the seal trade; see D-15 and O-12 | - | R-09 |
| O-6 | Seal friction for the actual seal profiles and sizes | calibrated to one published test series (60 mm coaxial seals); PU U-cup uses PU coaxial data | ISO 10100 module F on the first unit; the evaluator in `acceptance.py` fits and judges it | R-01, R-12 |
| O-7 | Hardware validation of every model | no cylinder built | run the acceptance procedure on the first unit, with a transducer in the cushion chamber (C-1) and strain gauges near the base weld | all |
| O-8 | The real duty spectrum | the duty cycle is a design assumption | customer data or a measured pressure history on the machine | R-08 |
| O-9 | Drift requirement wording (HD-160) | the order gives no thermal condition; 4.4 K of oil cooling uses the whole 5 mm | customer states "at constant oil temperature" or a cooling rate | R-13 |
| O-10 | Weld procedure and root inspection | the fatigue result assumes a full-penetration weld with an inspected root (detail category 71) | welding engineer confirms the joint design, WPS and NDT of the root | R-08 |
| O-11 | Material certificates | material values are marked `checked: false` in `data/materials.json` | mill certificates for the tube and rod heats | R-02 to R-08 |
| O-12 | Seal system choice (PTFE kept, PU now ranked first on judgement scores) | the trade is decided by judgement criteria (dry rod, contamination, cost), not by the models | supplier quote, the customer's contamination class and leak-tightness requirement | D-15, R-01, R-12 |
