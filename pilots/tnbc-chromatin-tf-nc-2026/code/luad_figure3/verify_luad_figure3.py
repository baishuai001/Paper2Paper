#!/usr/bin/env python3
"""Independent atomic checks for the LUAD Figure 3 transfer."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import fitz


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    tables = root / "results" / "tables"
    figures = root / "results" / "figures"
    audit = root / "audit"
    checks: list[dict[str, object]] = []

    def require(name: str, value: bool) -> None:
        checks.append({"check": name, "passed": bool(value)})

    composition = rows(tables / "Figure3A_regulon_composition.tsv")
    hcs = {row["TF"] for row in composition}
    require("31 HC-TFs in regulon composition", len(hcs) == 31)
    require("four signed/shared components per HC-TF", len(composition) == 31 * 4)
    require("nonnegative target counts", all(int(row["targets"]) >= 0 for row in composition))

    skew = rows(tables / "Figure3F_TCGA_LUAD_HC_TF_skewness.tsv")
    require("31 HC-TF skewness rows", len(skew) == 31)
    require("516 LUAD observations per skew test", all(int(row["n"]) == 516 for row in skew))
    require("finite skewness statistics", all(math.isfinite(float(row["skewness"])) for row in skew))

    motif = {row["TF"] for row in skew if row["motif_three_system"].lower() == "true"}
    require("three-system motif set remains FOXA3/NFATC4/XBP1",
            motif == {"FOXA3", "NFATC4", "XBP1"})

    modules = rows(tables / "Figure3B_HC_TF_modules.tsv")
    require("all HC-TFs assigned once to a module",
            len(modules) == 31 and len({row["TF"] for row in modules}) == 31)

    required = ("Figure3_complete.pdf", "Figure3_complete.png",
                "SupplementaryFigure6_complete.pdf", "SupplementaryFigure6_complete.png")
    for name in required:
        path = figures / name
        require(f"{name} exists and is nonempty", path.exists() and path.stat().st_size > 1000)
        if path.suffix == ".pdf" and path.exists():
            pdf = fitz.open(path)
            require(f"{name} has one page", len(pdf) == 1)
            pdf.close()

    receipt = {
        "verification_status": "passed" if all(x["passed"] for x in checks) else "failed",
        "checks": checks,
    }
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "figure3_independent_verification.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt))
    return 0 if receipt["verification_status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())

