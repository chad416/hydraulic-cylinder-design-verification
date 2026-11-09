# Progress notes (for resuming after an interruption)

Task (9 Nov 2026, from the fair-follow-up session): flagship cylinder project. Spine = Wasson 2nd ed.
Wasson extracts: scratchpad wasson/ch*.txt (book page = PDF page - 30).

Checked sources so far:
- ISO 10100:2020 table of contents (DIN preview): modules L (8.2 low-pressure leakage, 8.3 proof/external), P (9, piston seal leakage, optional), F (10, friction, optional; 10.4.1 sinusoidal, 10.4.2 constant speed; 10.6.2 static friction).
- ISO 10771-1:2015 scope: test pressure set by user, ISO/TR 10771-2 rating method; Annex B = cylinders (search summary of store pages).
- ISO/TS 13725:2021 scope: geometric model, off-axis load, own weight, fully extended, Hoblit + NFPA T3.6.37; not for internal rod buckling; spherical-bearing friction not included.
- Gholizadeh, Burton, Schoenau 2011, IJFP 12(3) 5-15 (read PDF).
- REACH Cr(VI) restriction: dossier 11 Apr 2025, consultation to 18 Dec 2025, RAC/SEAC opinions pending (Freyr 25 May 2026).
- E355+SR: Interpipe catalogue 450/580 MPa (search summary). My memory of 520/640 NOT confirmed.

Steps: [ ] package  [ ] tests  [ ] reports+PNGs  [ ] DOSSIER + PDF  [ ] git commit  [ ] PLAN.md update

## 9 Nov evening status
- Package modules written: model, stress, buckling, fatigue, seals, dynamics, analysis, requirements, acceptance, aoa, fmea, cnc_bridge.
- Design iterations recorded in data/orders/*.json (HD-160 rev A/B; LD-80 rev A/B/C).
- Key findings so far: weld at zero-moment point (beta x = pi/4) lifts HD weld fatigue 0.50 -> 1.09; cushion chamber not protected by relief (445 bar at 0 C, HD rev A); overrunning retract; LD-80 rev A buckling 0.96; stop tube optimum ~100 mm; Merritt Kv<2 zeta w matches model within 1 % for a fast valve; drift spec needs thermal condition; stick-slip onset 12 mm/s with assumed PTFE data.
- Next: plots.py, report.py, cli.py, tests, run, DOSSIER.md + PDF, git init + commit, PLAN.md.

## 9 Nov night: done
- [x] package [x] 24 tests [x] reports + PNGs [x] DOSSIER.md (+ PDF via tools/export_pdf.py) [x] PLAN.md, memory, shared note
- Final: HD-160 rev B (not met: R-12, R-13 = spec issues), LD-80 rev C (not met: R-12). CNC via 07: LD-80 3/3 PASS, HD-160 gland rejected (tool depth).
- Open: seal data assumed; ISO/TS 13725 equations, ISO 10100 values, ISO 10771 text not available; no hardware validation; CI never run (not pushed).
