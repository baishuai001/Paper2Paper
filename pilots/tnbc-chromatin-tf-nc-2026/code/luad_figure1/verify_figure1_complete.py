#!/usr/bin/env python3
"""Independent verification for the corrected and cross-platform Figure 1 package."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import fitz


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    processed = root / "data" / "processed"
    tables = root / "results" / "tables"
    figures = root / "results" / "figures"
    audit = root / "audit" / "figure1_complete"
    checks: list[dict[str, object]] = []

    def require(name: str, condition: bool) -> None:
        checks.append({"check": name, "passed": bool(condition)})

    for filename, expected in (("tcga_manifest.tsv", {"LUAD": 516, "LUSC": 501}),
                               ("gse81089_manifest.tsv", {"LUAD": 108, "LUSC": 67}),
                               ("gse41271_manifest.tsv", {"LUAD": 183, "LUSC": 80})):
        manifest = read_tsv(processed / filename)
        observed = {group: sum(row["group"] == group for row in manifest) for group in expected}
        require(f"{filename} frozen group counts", observed == expected)

    discovery = read_tsv(tables / "tcga_specific_TFs.tsv")
    require("158 LUAD TFs remain the Figure 2 input",
            sum(row["discovery_category"] == "LUAD" for row in discovery) == 158)
    cross = read_tsv(tables / "figure1_cross_platform_replication.tsv")
    require("cross-platform table covers every discovered TF once",
            len(cross) == len(discovery) and len({row["TF"] for row in cross}) == len(cross))

    expected_pdfs = {"Figure1_complete.pdf": 1, "SupplementaryFigure1_complete.pdf": 1,
                     "SupplementaryFigure2_complete.pdf": 2}
    for name, pages in expected_pdfs.items():
        path = figures / name
        require(f"{name} exists", path.exists() and path.stat().st_size > 1000)
        if path.exists():
            pdf = fitz.open(path)
            require(f"{name} page count", len(pdf) == pages)
            require(f"{name} pages have nonzero geometry",
                    all(page.rect.width > 0 and page.rect.height > 0 for page in pdf))
            pdf.close()

    receipt = {"verification_status": "passed" if all(x["passed"] for x in checks) else "failed",
               "checks": checks}
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "figure1_complete_independent_verification.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt))
    return 0 if receipt["verification_status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())

