"""Command line.

    python -m cylver report            # all orders in data/orders -> reports/*.md, docs/img/*.png, reports/cnc/
    python -m cylver report --fast     # skip the slow sweeps and the CNC run (CI smoke run)
    python -m cylver check data/orders/HD-160.json   # print the requirement status of every revision
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .report import ROOT, order_report


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="cylver")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("report")
    r.add_argument("--orders", default=str(ROOT / "data" / "orders"))
    r.add_argument("--out", default=str(ROOT / "reports"))
    r.add_argument("--img", default=str(ROOT / "docs" / "img"))
    r.add_argument("--fast", action="store_true")
    r.add_argument("--no-figures", action="store_true")
    c = sub.add_parser("check")
    c.add_argument("order")
    a = ap.parse_args(argv)
    if a.cmd == "report":
        for p in sorted(Path(a.orders).glob("*.json")):
            info = order_report(p, Path(a.out) / f"{p.stem}.md", Path(a.img), fast=a.fast, figures=not a.no_figures)
            print(f"{info['oid']}:")
            for s in info["summary"]:
                print("   " + ", ".join(f"{k} {v}" for k, v in s.items()))
            for cc in info["cnc"]:
                print(f"   CNC {cc['part']}: {cc['result']}")
        return 0
    from . import analysis, requirements
    from .model import build, load_order

    for name, raw in load_order(Path(a.order)):
        res = analysis.run(build(raw, name), fast=True)
        print(name)
        for q in requirements.derive(res):
            print(f"  {q['id']} {q['status']:8} {q['value']}  (limit {q['limit']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
