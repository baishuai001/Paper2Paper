#!/usr/bin/env python3
"""Compile TNBC replay, old LUAD, and harmonized LUAD motif diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from parse_harmonized_motifs import (
    LINEAGE_CONTROLS,
    build_token_map,
    parse_known_results,
    read_tsv,
    write_tsv,
)


SYSTEMS = ("patient", "PDX", "cell_line")


def load_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def bool_value(value: object) -> bool:
    return str(value).upper() == "TRUE"


def summarize_old_controls(root: Path) -> list[dict[str, object]]:
    inventory = read_tsv(root / "results/motif_gate/hc_tf_motif_inventory.tsv")
    hc_tfs = {row["TF"] for row in inventory if row["motif_testable"] == "TRUE"}
    token_map, _ = build_token_map(root, hc_tfs)
    sample_sets: dict[str, set[str]] = defaultdict(set)
    positive: dict[tuple[str, str], set[str]] = defaultdict(set)
    pattern_root = root / "results/motif_gate/homer_per_sample"
    for known in pattern_root.glob("*/*/*/knownResults.txt"):
        system = known.parts[-3]
        sample = known.parts[-2]
        sample_sets[system].add(sample)
        for row in parse_known_results(known, token_map):
            tf = str(row["TF"])
            if tf in LINEAGE_CONTROLS and row["motif_enriched_q1e-5"] == "TRUE":
                positive[(system, tf)].add(sample)
    rows: list[dict[str, object]] = []
    for system in SYSTEMS:
        total = len(sample_sets[system])
        cutoff = math.ceil(total / 2)
        for tf in LINEAGE_CONTROLS:
            count = len(positive[(system, tf)])
            rows.append(
                {
                    "analysis": "unharmonized_original",
                    "system": system,
                    "TF": tf,
                    "n_positive": count,
                    "denominator": total,
                    "half_sample_cutoff": cutoff,
                    "supported": str(count >= cutoff).upper(),
                    "denominator_type": "completed_samples_in_old_run",
                }
            )
    return rows


def harmonized_control_rows(path: Path, analysis: str, denominator_mode: str) -> list[dict[str, object]]:
    source = read_tsv(path)
    rows: list[dict[str, object]] = []
    for row in source:
        if denominator_mode == "fixed":
            denominator = int(row["n_all_samples"])
            cutoff = int(row["half_threshold_all_samples"])
            supported = row["support_fixed_original_denominator"]
        else:
            denominator = int(row["n_peak_count_matchable_samples"])
            cutoff = int(row["half_threshold_matchable_samples"])
            supported = row["support_matchable_denominator"]
        rows.append(
            {
                "analysis": analysis,
                "system": row["system"],
                "TF": row["TF"],
                "n_positive": int(row["n_motif_enriched"]),
                "denominator": denominator,
                "half_sample_cutoff": cutoff,
                "supported": supported,
                "denominator_type": denominator_mode,
            }
        )
    return rows


def support_count(rows: list[dict[str, object]], analysis: str, system: str, denominator_type: str) -> int:
    return sum(
        row["analysis"] == analysis
        and row["system"] == system
        and row["denominator_type"] == denominator_type
        and bool_value(row["supported"])
        for row in rows
    )


def load_hc_support(path: Path) -> dict[str, dict[str, str]]:
    return {row["TF"]: row for row in read_tsv(path)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("luad_root", type=Path)
    parser.add_argument("tnbc_root", type=Path)
    args = parser.parse_args()
    luad = args.luad_root.resolve()
    tnbc = args.tnbc_root.resolve()
    out = luad / "results/motif_equivalence/compiled_audit"
    out.mkdir(parents=True, exist_ok=True)

    tnbc_receipt = load_json(tnbc / "audit/tnbc_author_motif_replay_receipt.json")
    anchor_candidate_receipt = load_json(luad / "audit/motif_equivalence/candidate_anchor_equivalent_motif_receipt.json")
    combined_candidate_receipt = load_json(luad / "audit/motif_equivalence/candidate_database_motif_receipt.json")
    full_receipt = load_json(luad / "audit/motif_equivalence/harmonized_motif_receipt.json")
    rows = summarize_old_controls(luad)
    anchor_candidate_controls = luad / "results/motif_equivalence/candidate_anchor_equivalent/luad_lineage_positive_control_summary.tsv"
    combined_candidate_controls = luad / "results/motif_equivalence/candidate_database/luad_lineage_positive_control_summary.tsv"
    full_controls = luad / "results/motif_equivalence/luad_lineage_positive_control_summary.tsv"
    for mode in ("fixed", "matchable"):
        rows.extend(harmonized_control_rows(anchor_candidate_controls, "anchor_equivalent_candidate", mode))
        rows.extend(harmonized_control_rows(combined_candidate_controls, "combined_candidate_database", mode))
        rows.extend(harmonized_control_rows(full_controls, "harmonized_full_database", mode))
    write_tsv(
        out / "lineage_positive_control_comparison.tsv",
        rows,
        ["analysis", "system", "TF", "n_positive", "denominator", "half_sample_cutoff", "supported", "denominator_type"],
    )

    support_sources = {
        "anchor": load_hc_support(luad / "results/motif_equivalence/candidate_anchor_equivalent/harmonized_hc_tf_triple_system_summary.tsv"),
        "combined47": load_hc_support(luad / "results/motif_equivalence/candidate_database/harmonized_hc_tf_triple_system_summary.tsv"),
        "full1944": load_hc_support(luad / "results/motif_equivalence/harmonized_hc_tf_triple_system_summary.tsv"),
    }
    support_fields = (
        "patient_fixed_support", "PDX_fixed_support", "cell_line_fixed_support",
        "triple_system_fixed_original_denominator",
    )
    support_comparison: list[dict[str, object]] = []
    for tf in sorted(set().union(*(set(rows) for rows in support_sources.values()))):
        row: dict[str, object] = {"TF": tf}
        for analysis, source_rows in support_sources.items():
            for field in support_fields:
                row[f"{analysis}_{field}"] = source_rows[tf][field]
        support_comparison.append(row)
    write_tsv(out / "hc_tf_support_comparison.tsv", support_comparison)
    combined_full_identical = all(
        all(support_sources["combined47"][tf][field] == support_sources["full1944"][tf][field] for field in support_fields)
        for tf in support_sources["combined47"]
    )
    anchor_combined_changed_tfs = sorted(
        tf
        for tf in support_sources["anchor"]
        if any(support_sources["anchor"][tf][field] != support_sources["combined47"][tf][field] for field in support_fields)
    )

    anchor_candidate_fixed_pass = all(
        support_count(rows, "anchor_equivalent_candidate", system, "fixed") > 0
        for system in ("PDX", "cell_line")
    )
    anchor_candidate_matchable_pass = all(
        support_count(rows, "anchor_equivalent_candidate", system, "matchable") > 0
        for system in ("PDX", "cell_line")
    )
    combined_candidate_fixed_pass = all(
        support_count(rows, "combined_candidate_database", system, "fixed") > 0
        for system in ("PDX", "cell_line")
    )
    combined_candidate_matchable_pass = all(
        support_count(rows, "combined_candidate_database", system, "matchable") > 0
        for system in ("PDX", "cell_line")
    )
    full_fixed_pass = all(
        support_count(rows, "harmonized_full_database", system, "fixed") > 0
        for system in ("PDX", "cell_line")
    )
    full_matchable_pass = all(
        support_count(rows, "harmonized_full_database", system, "matchable") > 0
        for system in ("PDX", "cell_line")
    )
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tnbc_author_result_recovered": bool(tnbc_receipt["full_published_three_system_result_recovered"]),
        "tnbc_replayed_supported_any_system": int(tnbc_receipt["replayed_supported_in_at_least_one_system_count"]),
        "tnbc_replayed_triple_system": int(tnbc_receipt["replayed_triple_system_count"]),
        "old_unharmonized_positive_controls_supported": {
            system: support_count(rows, "unharmonized_original", system, "completed_samples_in_old_run")
            for system in SYSTEMS
        },
        "anchor_equivalent_candidate_positive_controls_supported_fixed_denominator": {
            system: support_count(rows, "anchor_equivalent_candidate", system, "fixed") for system in SYSTEMS
        },
        "anchor_equivalent_candidate_positive_controls_supported_matchable_denominator": {
            system: support_count(rows, "anchor_equivalent_candidate", system, "matchable") for system in SYSTEMS
        },
        "combined_candidate_positive_controls_supported_fixed_denominator": {
            system: support_count(rows, "combined_candidate_database", system, "fixed") for system in SYSTEMS
        },
        "combined_candidate_positive_controls_supported_matchable_denominator": {
            system: support_count(rows, "combined_candidate_database", system, "matchable") for system in SYSTEMS
        },
        "full_database_positive_controls_supported_fixed_denominator": {
            system: support_count(rows, "harmonized_full_database", system, "fixed") for system in SYSTEMS
        },
        "full_database_positive_controls_supported_matchable_denominator": {
            system: support_count(rows, "harmonized_full_database", system, "matchable") for system in SYSTEMS
        },
        "anchor_equivalent_candidate_calibration_pass_fixed_denominator": anchor_candidate_fixed_pass,
        "anchor_equivalent_candidate_calibration_pass_matchable_denominator": anchor_candidate_matchable_pass,
        "combined_candidate_calibration_pass_fixed_denominator": combined_candidate_fixed_pass,
        "combined_candidate_calibration_pass_matchable_denominator": combined_candidate_matchable_pass,
        "full_database_calibration_pass_fixed_denominator": full_fixed_pass,
        "full_database_calibration_pass_matchable_denominator": full_matchable_pass,
        "anchor_equivalent_candidate_hc_tf_triple_fixed_denominator": int(anchor_candidate_receipt["harmonized_triple_hc_tf_count_fixed_original_denominator"]),
        "anchor_equivalent_candidate_hc_tf_triple_matchable_denominator": int(anchor_candidate_receipt["harmonized_triple_hc_tf_count_matchable_denominator"]),
        "combined_candidate_hc_tf_triple_fixed_denominator": int(combined_candidate_receipt["harmonized_triple_hc_tf_count_fixed_original_denominator"]),
        "combined_candidate_hc_tf_triple_matchable_denominator": int(combined_candidate_receipt["harmonized_triple_hc_tf_count_matchable_denominator"]),
        "full_database_hc_tf_triple_fixed_denominator": int(full_receipt["harmonized_triple_hc_tf_count_fixed_original_denominator"]),
        "full_database_hc_tf_triple_matchable_denominator": int(full_receipt["harmonized_triple_hc_tf_count_matchable_denominator"]),
        "combined47_and_full1944_hc_support_matrices_identical": combined_full_identical,
        "anchor_vs_combined47_changed_hc_tfs": anchor_combined_changed_tfs,
        "interpretation": (
            "anchor-equivalent source-split/JASPAR-priority positive controls calibrate both model arms; HC-TF motif results are technically interpretable under the exact author selection rule"
            if anchor_candidate_fixed_pass
            else (
                "the peak-to-HOMER chain detects lineage controls in both model arms when all catalogued candidate motifs are allowed, but the exact author JASPAR-priority selection does not; any zero under the exact rule is motif-model-selection-sensitive rather than a clean biological absence"
                if combined_candidate_fixed_pass
                else "at least one model arm fails both exact-author and all-candidate positive-control calibration; HC-TF motif absence remains technically unresolved"
            )
        ),
        "combined_candidate_role": "more stringent sensitivity analysis because JASPAR and CIS-BP q values are corrected jointly rather than independently",
        "full_database_role": "additional stringent 1,944-model multiple-testing sensitivity analysis; it does not replace the source-split anchor-equivalent candidate diagnostic",
    }
    receipt_path = luad / "audit/motif_equivalence/final_motif_equivalence_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
