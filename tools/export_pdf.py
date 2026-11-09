"""Export DOSSIER.md to DOSSIER.pdf (Markdown -> HTML -> headless Edge/Chrome print).

    pip install markdown
    python tools/export_pdf.py [--browser PATH]
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent.parent
CSS = """
body { font-family: 'Segoe UI', Arial, sans-serif; font-size: 10pt; line-height: 1.4; color: #0b0b0b; max-width: 180mm; margin: auto; }
h1 { font-size: 18pt; } h2 { font-size: 13.5pt; border-bottom: 1px solid #ccc; margin-top: 1.4em; } h3 { font-size: 11pt; }
table { border-collapse: collapse; width: 100%; font-size: 8.3pt; margin: 0.6em 0; page-break-inside: auto; }
th, td { border: 1px solid #c8c8c8; padding: 2px 4px; vertical-align: top; } th { background: #f0efec; }
code, pre { font-size: 8.5pt; background: #f6f6f4; } pre { padding: 6px; white-space: pre-wrap; }
img { max-width: 100%; } blockquote { color: #52514e; border-left: 3px solid #ccc; margin-left: 0; padding-left: 8px; }
@page { size: A4; margin: 14mm 12mm; }
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(ROOT / "DOSSIER.md"))
    ap.add_argument("--browser", default=None)
    a = ap.parse_args()
    src = Path(a.src)
    html = markdown.markdown(src.read_text(encoding="utf-8"), extensions=["tables", "fenced_code", "toc"])
    page = f"<!doctype html><html><head><meta charset='utf-8'><title>{src.stem}</title><style>{CSS}</style></head><body>{html}</body></html>"
    out_html = src.with_suffix(".html")
    out_html.write_text(page, encoding="utf-8")
    browser = a.browser or shutil.which("msedge") or r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"
    out_pdf = src.with_suffix(".pdf")
    if out_pdf.exists():
        out_pdf.unlink()  # so an old file cannot pass for a new one
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as prof:  # own profile: a running browser cannot swallow the job
        subprocess.run([browser, "--headless=new", "--disable-gpu", f"--user-data-dir={prof}", "--no-pdf-header-footer", f"--print-to-pdf={out_pdf}", out_html.resolve().as_uri()], check=True, timeout=180)
        for _ in range(30):
            if out_pdf.exists():
                break
            time.sleep(1)
    out_html.unlink()
    if not out_pdf.exists():
        raise SystemExit("no PDF written")
    print(out_pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
