#!/usr/bin/env python3
"""Independent checks for Figure 4 Cox, permutation and figure outputs."""

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

    def require(name: str, condition: bool) -> None:
        checks.append({"check": name, "passed": bool(condition)})

    models = rows(tables / "Figure4_all_Cox_models.tsv")
    require("31 TF x 2 cohorts x 2 endpoints x 2 model types", len(models) == 31 * 2 * 2 * 2)
    require("all Cox estimates finite", all(math.isfinite(float(row[key])) for row in models
                                             for key in ("HR", "lower95", "upper95", "p_value")))
    require("all Cox analyses have at least 100 complete cases",
            all(int(row["n"]) >= 100 for row in models))
    require("all Cox analyses have at least 20 ten-year events",
            all(int(row["events"]) >= 20 for row in models))

    permutation = rows(tables / "SupplementaryFigure7_permutation_summary.tsv")
    require("eight independent permutation analyses", len(permutation) == 8)
    require("all use exactly 5,000 permutations", all(int(row["permutations"]) == 5000 for row in permutation))
    require(
        "all permutation null distributions are non-degenerate",
        all(float(row["null_sd"]) > 0 for row in permutation),
    )
    require("finite empirical p and fold enrichment", all(
        math.isfinite(float(row["empirical_p"])) and math.isfinite(float(row["fold_enrichment"]))
        for row in permutation
    ))

    km = rows(tables / "Figure4FH_KM_receipts.tsv")
    require("two representatives x two endpoints x two cohorts", len(km) == 8)
    require("KM groups contain at least 10 patients", all(
        int(row["low_n"]) >= 10 and int(row["high_n"]) >= 10 for row in km
    ))

    for name in ("Figure4_complete.pdf", "SupplementaryFigure7_complete.pdf"):
        path = figures / name
        require(f"{name} exists", path.exists() and path.stat().st_size > 1000)
        if path.exists():
            pdf = fitz.open(path)
            require(f"{name} is one page", len(pdf) == 1)
            if name == "Figure4_complete.pdf" and len(pdf) == 1:
                composite_text = pdf[0].get_text()
                require(
                    "Figure4 composite contains all eight KM panels",
                    composite_text.count("Number at risk") >= 8
                    and all(token in composite_text for token in ("ETV1", "ZNF444", "CREBRF", "PHC2")),
                )
            pdf.close()

    receipt = {"verification_status": "passed" if all(x["passed"] for x in checks) else "failed",
               "checks": checks}
    (audit / "figure4_independent_verification.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt))
    return 0 if receipt["verification_status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
