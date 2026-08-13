#!/usr/bin/env python3
"""Independently verify LUAD Figure 1 manifests, TF sets, and hard verdict."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__, reverse=True)
    adjusted = [1.0] * len(values)
    running = 1.0
    total = len(values)
    for index in order:
        rank = sum(values[other] <= values[index] for other in range(total))
        running = min(running, values[index] * total / rank)
        adjusted[index] = min(running, 1.0)
    return adjusted


def count_groups(manifest: list[dict[str, str]]) -> dict[str, int]:
    return {group: sum(row["group"] == group for row in manifest) for group in ("LUAD", "LUSC")}


def require(condition: bool, message: str, checks: list[dict[str, object]]) -> None:
    checks.append({"check": message, "passed": bool(condition)})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    processed = root / "data" / "processed"
    if not processed.is_dir() and (root / "manifests").is_dir():
        processed = root / "manifests"
    tables = root / "results" / "tables"
    figures = root / "results" / "figures"
    audit = root / "audit"
    checks: list[dict[str, object]] = []

    manifests = {
        "TCGA": rows(processed / "tcga_manifest.tsv"),
        "GSE81089": rows(processed / "gse81089_manifest.tsv"),
        "PDMR_PDX": rows(processed / "pdmr_manifest.tsv"),
        "DepMap_22Q2": rows(processed / "depmap_manifest.tsv"),
    }
    expected_counts = {
        "TCGA": {"LUAD": 516, "LUSC": 501},
        "GSE81089": {"LUAD": 106, "LUSC": 67},
        "PDMR_PDX": {"LUAD": 22, "LUSC": 34},
        "DepMap_22Q2": {"LUAD": 76, "LUSC": 27},
    }
    for cohort, manifest in manifests.items():
        require(count_groups(manifest) == expected_counts[cohort], f"{cohort} frozen group counts", checks)
        sample_ids = [row["sample_id"] for row in manifest]
        require(len(sample_ids) == len(set(sample_ids)), f"{cohort} unique sample IDs", checks)
    require(
        len({row["patient_id"] for row in manifests["TCGA"]}) == len(manifests["TCGA"]),
        "TCGA one primary tumor per patient",
        checks,
    )
    require(
        len({row["patient_id"] for row in manifests["PDMR_PDX"]}) == len(manifests["PDMR_PDX"]),
        "PDMR one PDX per donor",
        checks,
    )

    selected = rows(tables / "tcga_specific_TFs.tsv")
    replicated_file = rows(tables / "externally_replicated_TFs.tsv")
    replicated_from_selected = {row["TF"] for row in selected if row["externally_replicated"].lower() == "true"}
    require(
        replicated_from_selected == {row["TF"] for row in replicated_file},
        "replicated TF file equals flagged discovery TFs",
        checks,
    )
    tcga = rows(tables / "figure1b_tcga_limma_msviper.tsv")
    gse = rows(tables / "gse81089_limma_msviper.tsv")
    expected_replicated = {
        row["TF"]
        for row in tcga
        if row["category"] in {"LUAD", "LUSC"}
        for validation in gse
        if validation["TF"] == row["TF"] and validation["category"] == row["category"]
    }
    require(replicated_from_selected == expected_replicated, "independent-patient TF intersection", checks)
    for cohort in ("TCGA", "GSE81089"):
        ms_rows = rows(root / "results" / cohort.lower() / f"{cohort}_msviper.tsv")
        require(
            bool(ms_rows) and all(
                math.isfinite(float(row[column]))
                for row in ms_rows
                for column in ("NES", "p_value", "FDR", "regulon_size")
            ),
            f"{cohort} msVIPER values finite",
            checks,
        )

    program = rows(tables / "LUAD_program_cross_system_statistics.tsv")
    p_values = [float(row["p_value"]) for row in program]
    independently_adjusted = bh(p_values)
    require(
        all(abs(float(row["FDR"]) - value) < 1e-8 for row, value in zip(program, independently_adjusted)),
        "program BH adjustment independently reproduced",
        checks,
    )
    by_cohort = {row["cohort"]: row for row in program}
    replicated_luad = sum(row["discovery_category"] == "LUAD" for row in replicated_file)
    model_pass = all(
        float(by_cohort[cohort]["mean_difference"]) > 0
        and float(by_cohort[cohort]["FDR"]) <= 0.05
        and float(by_cohort[cohort]["hedges_g"]) >= 0.5
        for cohort in ("PDMR_PDX", "DepMap_22Q2")
    )
    biological_gate = replicated_luad >= 3 and model_pass
    result_receipt = json.loads((audit / "figure1_result_receipt.json").read_text(encoding="utf-8"))
    execution = result_receipt["method_execution"]
    require(
        execution == {
            "ttestNull": "viper::ttestNull",
            "permutations": 1000,
            "replacement": True,
            "seed": 1,
            "cores": 32,
            "rng_kind": "L'Ecuyer-CMRG",
            "adaptation_timing": "frozen after serial runtime measurement at 2 percent and before any msVIPER/TF result",
            "scientific_definition_changed": False,
            "viper_1_38_aREA_single_target_dimension_patch": True,
            "compatibility_patch_effect": "shape restoration only; minsize=1 and numeric definition unchanged",
        },
        "parallel null adaptation fully disclosed",
        checks,
    )
    for cohort in ("TCGA", "GSE81089"):
        viper_receipt = json.loads(
            (root / "results" / cohort.lower() / f"{cohort}_viper_receipt.json").read_text(encoding="utf-8")
        )
        require(
            viper_receipt["area_single_target_dimension_patch"] is True
            and viper_receipt["area_patch_effect"]
            == "shape restoration only; minsize and numeric definition unchanged",
            f"{cohort} viper 1.38.0 compatibility patch disclosed",
            checks,
        )
    require(result_receipt["status"] == ("passed" if biological_gate else "failed"), "hard verdict independently reproduced", checks)

    required_panels = [f"Figure1{panel}_" for panel in "ABCDE"]
    figure_names = [path.name for path in figures.iterdir() if path.is_file() and path.stat().st_size > 0]
    for prefix in required_panels:
        require(any(name.startswith(prefix) and name.endswith(".pdf") for name in figure_names), f"{prefix} PDF exists", checks)
        require(any(name.startswith(prefix) and name.endswith(".png") for name in figure_names), f"{prefix} PNG exists", checks)

    receipt = {
        "verification_status": "passed" if all(check["passed"] for check in checks) else "failed",
        "biological_gate": "PASS" if biological_gate else "FAIL",
        "checks": checks,
        "manifest_counts": {cohort: count_groups(manifest) for cohort, manifest in manifests.items()},
        "replicated_LUAD_TFs": replicated_luad,
    }
    output = audit / "figure1_independent_verification.json"
    output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt))
    return 0 if receipt["verification_status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
