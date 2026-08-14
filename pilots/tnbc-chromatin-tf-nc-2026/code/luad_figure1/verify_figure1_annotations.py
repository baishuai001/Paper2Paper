#!/usr/bin/env python3
"""Independent checks for the Figure 1 annotation/supplement-only run."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cloud_run_root", type=Path)
    args = parser.parse_args()
    root = args.cloud_run_root.resolve()
    annotations = root / "data" / "processed" / "annotations"
    figures = root / "results" / "figures"
    tables = root / "results" / "tables"
    audit = root / "audit" / "figure1_annotations"
    audit.mkdir(parents=True, exist_ok=True)

    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    expected_rows = {
        "tcga_figure1_annotations.tsv": 1017,
        "gse81089_figure1_annotations.tsv": 173,
        "pdmr_figure1_annotations.tsv": 56,
        "depmap_figure1_annotations.tsv": 103,
    }
    loaded: dict[str, list[dict[str, str]]] = {}
    for filename, expected in expected_rows.items():
        path = annotations / filename
        rows = read_tsv(path) if path.exists() else []
        loaded[filename] = rows
        check(f"row_count:{filename}", len(rows) == expected, {"observed": len(rows), "expected": expected})
        sample_ids = [row.get("sample_id", "") for row in rows]
        check(f"unique_sample_id:{filename}", len(set(sample_ids)) == len(sample_ids) and "" not in sample_ids, len(set(sample_ids)))

    tcga = loaded["tcga_figure1_annotations.tsv"]
    subtype_n = sum(row.get("published_expression_subtype", "") not in {"", "NA"} for row in tcga)
    check("published_TCGA_expression_subtype_coverage", subtype_n >= 380, subtype_n)

    pdmr = loaded["pdmr_figure1_annotations.tsv"]
    pdmr_ok = sum(row.get("parse_status") == "ok" for row in pdmr)
    check("PDMR_public_page_parse_coverage", pdmr_ok >= 45, {"parsed_ok": pdmr_ok, "total": len(pdmr)})

    sha_diff = audit / "frozen_numeric_sha256.diff"
    diff_text = sha_diff.read_text(encoding="utf-8") if sha_diff.exists() else "missing diff file"
    check("frozen_Figure1_numeric_outputs_unchanged", diff_text.strip() == "", diff_text[:2000])

    required_figures = [
        "Figure1C_TCGA_TF_activity_annotated.pdf",
        "Figure1D_PDMR_PDX_TF_activity_annotated.pdf",
        "Figure1E_DepMap_22Q2_TF_activity_annotated.pdf",
        "SupplementaryFigure1_inclusion_flow.pdf",
        "SupplementaryFigure2A_GSE81089_logFC_msviper.pdf",
        "SupplementaryFigure2B_TF_overlap.pdf",
        "SupplementaryFigure2C_TCGA_sample_correlations.pdf",
        "SupplementaryFigure2D_PDMR_PDX_sample_correlations.pdf",
        "SupplementaryFigure2E_DepMap_22Q2_sample_correlations.pdf",
        "SupplementaryFigure2F_replicated_LUAD_TF_effects.pdf",
    ]
    for filename in required_figures:
        path = figures / filename
        check(f"figure_exists:{filename}", path.exists() and path.stat().st_size > 1000, path.stat().st_size if path.exists() else 0)

    required_tables = [
        "SupplementaryFigure1_inclusion_counts.tsv",
        "SupplementaryFigure2B_TF_overlap_counts.tsv",
        "SupplementaryFigure2F_replicated_LUAD_TF_effect_zscores.tsv",
        "SupplementaryFigures1_2_core_index.tsv",
    ]
    for filename in required_tables:
        path = tables / filename
        check(f"table_exists:{filename}", path.exists() and path.stat().st_size > 20, path.stat().st_size if path.exists() else 0)

    passed = all(item["passed"] for item in checks)
    receipt = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if passed else "FAIL",
        "scope": "Figure 1 annotation-only replot plus Supplementary Figures 1-2 core",
        "checks": checks,
    }
    output = audit / "figure1_annotations_independent_verification.json"
    output.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": receipt["status"], "checks": len(checks), "failed": [x["check"] for x in checks if not x["passed"]]}, ensure_ascii=False))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
