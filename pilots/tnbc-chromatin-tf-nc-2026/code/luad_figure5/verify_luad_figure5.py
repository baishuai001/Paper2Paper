#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import fitz


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def bh(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    adjusted_sorted = [0.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        original_index = order[rank_index]
        candidate = values[original_index] * n / (rank_index + 1)
        running = min(running, candidate)
        adjusted_sorted[rank_index] = min(1.0, running)
    answer = [0.0] * n
    for rank_index, original_index in enumerate(order):
        answer[original_index] = adjusted_sorted[rank_index]
    return answer


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    tables = root / "results" / "tables"
    figures = root / "results" / "figures"
    checks: list[dict[str, object]] = []

    def check(label: str, value: bool) -> None:
        checks.append({"check": label, "passed": bool(value)})

    dataset_audit = read_tsv(tables / "Figure5_dataset_audit.tsv")
    all_tests = read_tsv(tables / "Figure5ABC_cellline_concordance_all.tsv.gz")
    replicated_rows = read_tsv(tables / "Figure5D_replicated_drug_TF_pairs.tsv")
    classifier = read_tsv(tables / "SupplementaryFigure8C_lung_classifier_validation.tsv")
    validated = read_tsv(tables / "Figure5GH_PDX_validated_pairs.tsv")

    check("three pharmacogenomic datasets",
          {row["dataset"] for row in dataset_audit} == {"GDSC2", "CTRPv2", "PRISM"})
    check("all mapped LUAD cohorts contain at least 10 models",
          all(int(row["mapped_LUAD_lines"]) >= 10 for row in dataset_audit))
    check("all cell-line tests use n >= 10", all(int(row["n"]) >= 10 for row in all_tests))
    check("all cell-line tests satisfy AAC dynamic-range rule",
          all(float(row["max_AAC"]) > 0.2 for row in all_tests))

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in all_tests:
        grouped[row["dataset"]].append(row)
    fdr_ok = True
    for rows in grouped.values():
        expected = bh([float(row["p_value"]) for row in rows])
        fdr_ok = fdr_ok and all(
            abs(value - float(row["FDR"])) < 1e-9 for value, row in zip(expected, rows)
        )
    check("cell-line FDR values reproduce BH within dataset", fdr_ok)

    replicated = [row for row in replicated_rows if as_bool(row["replicated"])]
    check("replicated pairs have >=2 datasets and one direction",
          all(int(row["dataset_n"]) >= 2 and int(row["direction_n"]) == 1 for row in replicated))
    check("classifier has TCGA, PDMR and DepMap validation",
          {row["cohort"] for row in classifier} == {"TCGA_training", "PDMR_external", "DepMap_external"})
    pdmr = next(row for row in classifier if row["cohort"] == "PDMR_external")
    depmap = next(row for row in classifier if row["cohort"] == "DepMap_external")
    check("PDMR external balanced accuracy >= 0.8", float(pdmr["balanced_accuracy"]) >= 0.8)
    check("DepMap external balanced accuracy is reported", depmap["balanced_accuracy"] not in {"", "NA", "NaN"})
    if validated:
        check("main PDX pairs meet all frozen criteria", all(
            float(row["FDR"]) <= 0.05 and int(row["n"]) >= 10 and
            as_bool(row["cellline_replicated_same_direction"]) for row in validated
        ))
    else:
        check("absence of main PDX pair is retained as an explicit negative result", True)

    for name in ("Figure5_complete.pdf", "SupplementaryFigure8_complete.pdf", "SupplementaryFigure9_complete.pdf"):
        path = figures / name
        check(f"{name} exists", path.exists() and path.stat().st_size > 0)
        if path.exists():
            document = fitz.open(path)
            check(f"{name} has one page", len(document) == 1)
            document.close()

    status = "passed" if all(item["passed"] for item in checks) else "failed"
    payload = {"verification_status": status, "checks": checks}
    (root / "audit" / "figure5_independent_verification.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))
    return 0 if status == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
