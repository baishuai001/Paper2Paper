#!/usr/bin/env python3
"""Independent checks for the GSE41271 cross-platform Figure 1 layer."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(checks: list[dict[str, object]], condition: bool, name: str) -> None:
    checks.append({"check": name, "passed": bool(condition)})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    processed = root / "data" / "processed"
    results = root / "results"
    tables = results / "tables"
    figures = results / "figures"
    audit = root / "audit" / "gse41271"
    checks: list[dict[str, object]] = []

    manifest = read_tsv(processed / "gse41271_manifest.tsv")
    counts = {
        group: sum(row["group"] == group for row in manifest)
        for group in ("LUAD", "LUSC")
    }
    require(checks, counts == {"LUAD": 183, "LUSC": 80}, "official 183/80 histology counts")
    require(checks, len({row["geo_accession"] for row in manifest}) == 263,
            "263 unique GEO accessions")

    with (processed / "gse41271_aracne_expression.tsv").open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        first = handle.readline().rstrip("\n").split("\t")
    require(checks, len(header) == 264 and len(first) == 264,
            "gene-by-263-sample expression matrix is rectangular")
    require(checks, header[1:] == [row["geo_accession"] for row in manifest],
            "expression columns equal manifest order")

    viper_receipt = json.loads(
        (results / "gse41271" / "GSE41271_viper_receipt.json").read_text(encoding="utf-8")
    )
    require(checks, viper_receipt["status"] == "passed", "VIPER receipt passed")
    require(checks, viper_receipt["samples"] == 263, "VIPER used 263 samples")
    require(checks, viper_receipt["group_counts"] == {"LUAD": 183, "LUSC": 80},
            "VIPER group counts match manifest")

    table = read_tsv(tables / "gse41271_limma_msviper.tsv")
    require(checks, bool(table), "GSE41271 limma/msVIPER table is nonempty")
    require(
        checks,
        all(math.isfinite(float(row[column])) for row in table
            for column in ("NES", "p_value", "FDR", "logFC", "adj.P.Val")),
        "all core statistics are finite",
    )

    discovery = read_tsv(tables / "figure1b_tcga_limma_msviper.tsv")
    replication = read_tsv(tables / "figure1_cross_platform_replication.tsv")
    gse_by_tf = {row["TF"]: row for row in table}
    expected = {
        (row["TF"], row["category"])
        for row in discovery
        if row["category"] in {"LUAD", "LUSC"}
        and row["TF"] in gse_by_tf
        and gse_by_tf[row["TF"]]["category"] == row["category"]
    }
    observed = {
        (row["TF"], row["discovery_category"])
        for row in replication
        if row["replicated_GSE41271"].lower() == "true"
    }
    require(checks, expected == observed, "cross-platform TF replication recomputed independently")

    required_figures = (
        "SupplementaryFigure2G_GSE41271_logFC_msviper.pdf",
        "SupplementaryFigure2G_GSE41271_logFC_msviper.png",
        "SupplementaryFigure2H_patient_cohort_TF_overlap.pdf",
        "SupplementaryFigure2H_patient_cohort_TF_overlap.png",
    )
    for name in required_figures:
        path = figures / name
        require(checks, path.exists() and path.stat().st_size > 0, f"{name} exists")

    receipt = {
        "verification_status": "passed" if all(item["passed"] for item in checks) else "failed",
        "manifest_counts": counts,
        "checks": checks,
    }
    audit.mkdir(parents=True, exist_ok=True)
    (audit / "gse41271_independent_verification.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt))
    return 0 if receipt["verification_status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())

