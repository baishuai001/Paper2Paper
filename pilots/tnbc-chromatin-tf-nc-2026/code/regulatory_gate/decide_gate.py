#!/usr/bin/env python3
"""Combine immutable receipts into the frozen PASS/FAIL/INDETERMINATE verdict."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from common import write_json


def load(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-audit", required=True, type=Path)
    parser.add_argument("--pango", required=True, type=Path)
    parser.add_argument("--tcga", required=True, type=Path)
    parser.add_argument("--aracne-input", required=True, type=Path)
    parser.add_argument("--aracne-run", required=True, type=Path)
    parser.add_argument("--pseudobulk", required=True, type=Path)
    parser.add_argument("--viper", required=True, type=Path)
    parser.add_argument("--statistics", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    audit = load(args.input_audit)
    pango = load(args.pango)
    tcga = load(args.tcga)
    aracne_input = load(args.aracne_input)
    aracne_run = load(args.aracne_run)
    pseudobulk = load(args.pseudobulk)
    viper = load(args.viper)
    statistics = load(args.statistics)
    primary_counts = pseudobulk["group_counts_by_cell_threshold"]["50"]

    data_conditions = {
        "immutable_inputs_passed": audit["status"] == "passed",
        "pango_discrepancy_audited_and_2059_unique_symbols": (
            pango["methods_claimed_symbols"] == 2139
            and pango["annotation_records"] == 2138
            and pango["unique_nonempty_symbols"] == 2059
        ),
        "tcga_participants_ge_550": tcga["participants"] >= 550,
        "tcga_genes_ge_10000": tcga["retained_genes"] >= 10_000,
        "tcga_duplicate_symbols_use_anchor_mean": (
            tcga["processing"]["duplicate_symbol_aggregation"]
            == "mean (anchor-code behavior)"
        ),
        "tcga_has_no_added_low_expression_filter": (
            tcga["processing"]["low_expression_filter"]
            == "none (anchor-code behavior)"
        ),
        "aracne_regulators_ge_1500": aracne_input["aracne_regulators"] >= 1_500,
        "anchor_code_style_aracne_single_subnetwork_passed": (
            aracne_run["status"] == "passed"
            and aracne_run["subnetwork_files"] == 1
            and aracne_run["regulon_source"] == "subnets/subnet1_crc.tsv"
        ),
        "author_subnetwork_regulators_ge_500": (
            viper["network"]["author_subnetwork_regulators"] >= 500
        ),
        "measured_regulons_ge1_target_ge_400": (
            viper["network"]["regulons_ge1_measured_target"] >= 400
        ),
        "primary_case_ge_30": primary_counts["case"] >= 30,
        "primary_control_ge_80": primary_counts["control"] >= 80,
        **statistics["cohort_conditions"],
    }
    scientific_conditions = dict(statistics["scientific_conditions"])

    if not all(data_conditions.values()):
        verdict = "INDETERMINATE"
        failed = [name for name, passed in data_conditions.items() if not passed]
        interpretation = "Data, network or execution minimums were not all met; this is not a biological negative result."
    elif all(scientific_conditions.values()):
        verdict = "PASS"
        failed = []
        interpretation = "The predefined M-vs-rest contrast has a reproducible CRC Cancer-cell TF program and passes LODO validation."
    else:
        verdict = "FAIL"
        failed = [name for name, passed in scientific_conditions.items() if not passed]
        interpretation = "The data were adequate, but M-vs-rest did not satisfy every frozen TF-program and LODO requirement."

    receipt = {
        "verdict": verdict,
        "stop_now": True,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "route": "crc_atlas_M_vs_rest__TCGA_CRC_ARACNe3__VIPER_msVIPER",
        "data_conditions": data_conditions,
        "scientific_conditions": scientific_conditions,
        "failed_conditions": failed,
        "interpretation": interpretation,
        "next_action": (
            "Freeze a different CRC subtype manifest and reuse the CRC ARACNe3 regulon plus patient VIPER matrix."
            if verdict == "FAIL"
            else "Stop this gate and review the report before any downstream chromatin or therapy work."
        ),
        "claim_ceiling": "Cross-dataset association; no causal, treatment-selection or clinical-utility claim.",
    }
    write_json(args.output, receipt)
    print(json.dumps(receipt, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
