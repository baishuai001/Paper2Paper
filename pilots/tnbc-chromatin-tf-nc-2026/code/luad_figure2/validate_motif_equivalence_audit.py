#!/usr/bin/env python3
"""Independently validate the completed TNBC/LUAD motif-equivalence audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path


SYSTEMS = ("patient", "PDX", "cell_line")
EXPECTED_ALL = {"patient": 22, "PDX": 13, "cell_line": 19}
EXPECTED_MATCHABLE = {"patient": 22, "PDX": 10, "cell_line": 19}
EXPECTED_TARGET_PEAKS = 3584
EXPECTED_BACKGROUND_PEAKS = 28672
EXPECTED_WIDTH = 200
EXPECTED_CANDIDATE_MOTIFS = 47
EXPECTED_FULL_MOTIFS = 1944


def read_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_motifs(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(line.startswith(">") for line in handle)


def bed_profile(path: Path) -> tuple[int, set[int]]:
    count = 0
    widths: set[int] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise AssertionError(f"Malformed BED row in {path}: {line.rstrip()}")
            widths.add(int(fields[2]) - int(fields[1]))
            count += 1
    return count, widths


def truth(value: object) -> bool:
    return str(value).upper() == "TRUE"


def add_check(checks: list[dict[str, object]], name: str, observed: object, expected: object) -> None:
    checks.append(
        {
            "check": name,
            "observed": json.dumps(observed, sort_keys=True) if isinstance(observed, (dict, list)) else observed,
            "expected": json.dumps(expected, sort_keys=True) if isinstance(expected, (dict, list)) else expected,
            "pass": str(observed == expected).upper(),
        }
    )


def validate_result_set(
    root: Path,
    label: str,
    receipt_path: Path,
    result_dir: Path,
    expected_motifs: int,
    checks: list[dict[str, object]],
) -> None:
    receipt = read_json(receipt_path)
    add_check(checks, f"{label}.motif_count", int(receipt["unified_motif_database_motif_count"]), expected_motifs)
    add_check(
        checks,
        f"{label}.completed_samples",
        receipt["completed_or_explicitly_not_testable_samples_by_system"],
        EXPECTED_ALL,
    )
    add_check(
        checks,
        f"{label}.matchable_samples",
        receipt["peak_count_matchable_samples_by_system"],
        EXPECTED_MATCHABLE,
    )
    best = read_tsv(result_dir / "harmonized_sample_tf_best_motif.tsv")
    summary = read_tsv(result_dir / "harmonized_system_tf_motif_summary.tsv")
    triple = read_tsv(result_dir / "harmonized_hc_tf_triple_system_summary.tsv")
    controls = read_tsv(result_dir / "luad_lineage_positive_control_summary.tsv")
    hc_count = int(receipt["motif_testable_luad_hc_tf_count"])
    all_tf_count = len({row["TF"] for row in summary})
    add_check(checks, f"{label}.best_table_rows", len(best), sum(EXPECTED_ALL.values()) * all_tf_count)
    add_check(checks, f"{label}.system_summary_rows", len(summary), len(SYSTEMS) * all_tf_count)
    add_check(checks, f"{label}.hc_triple_rows", len(triple), hc_count)
    add_check(checks, f"{label}.lineage_control_rows", len(controls), len(SYSTEMS) * 8)
    add_check(
        checks,
        f"{label}.triple_fixed_receipt_matches_table",
        int(receipt["harmonized_triple_hc_tf_count_fixed_original_denominator"]),
        sum(truth(row["triple_system_fixed_original_denominator"]) for row in triple),
    )
    add_check(
        checks,
        f"{label}.triple_matchable_receipt_matches_table",
        int(receipt["harmonized_triple_hc_tf_count_matchable_denominator"]),
        sum(truth(row["triple_system_matchable_denominator"]) for row in triple),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("luad_root", type=Path)
    parser.add_argument("tnbc_root", type=Path)
    args = parser.parse_args()
    luad = args.luad_root.resolve()
    tnbc = args.tnbc_root.resolve()
    checks: list[dict[str, object]] = []

    tnbc_receipt = read_json(tnbc / "audit/tnbc_author_motif_replay_receipt.json")
    add_check(checks, "tnbc.full_three_system_result_recovered", tnbc_receipt["full_published_three_system_result_recovered"], True)
    add_check(checks, "tnbc.triple_tf_count", int(tnbc_receipt["replayed_triple_system_count"]), 31)
    add_check(checks, "tnbc.intersection_mismatch_count", int(tnbc_receipt["intersection_mismatch_count"]), 0)
    add_check(checks, "tnbc.lor_mismatch_count", int(tnbc_receipt["lor_summary_mismatch_count_at_tolerance_1e-12"]), 0)

    harmonized = luad / "data/processed/motif_equivalence/harmonized_200bp"
    manifest = read_tsv(harmonized / "harmonized_motif_manifest.tsv")
    all_counts = {system: sum(row["system"] == system for row in manifest) for system in SYSTEMS}
    matchable_counts = {
        system: sum(row["system"] == system and truth(row["motif_matchable"]) for row in manifest)
        for system in SYSTEMS
    }
    add_check(checks, "inputs.all_samples", all_counts, EXPECTED_ALL)
    add_check(checks, "inputs.matchable_samples", matchable_counts, EXPECTED_MATCHABLE)
    for row in manifest:
        path = Path(row["matched_target_path"])
        count, widths = bed_profile(path)
        expected_count = EXPECTED_TARGET_PEAKS if truth(row["motif_matchable"]) else 0
        add_check(checks, f"target.{row['system']}.{row['sample_slug']}.count", count, expected_count)
        add_check(
            checks,
            f"target.{row['system']}.{row['sample_slug']}.width",
            sorted(widths),
            [EXPECTED_WIDTH] if expected_count else [],
        )
    background_count, background_widths = bed_profile(harmonized / "LUSC_accessible_common_GC_matched_background.bed")
    add_check(checks, "background.count", background_count, EXPECTED_BACKGROUND_PEAKS)
    add_check(checks, "background.width", sorted(background_widths), [EXPECTED_WIDTH])

    candidate_db = luad / "reference/motifs/candidate_luad_hc_and_lineage_controls.homer"
    candidate_jaspar = luad / "reference/motifs/candidate_luad_hc_and_lineage_controls.JASPAR2024.homer"
    candidate_cisbp = luad / "reference/motifs/candidate_luad_hc_and_lineage_controls.CIS-BP2.00.homer"
    full_db = luad / "reference/motifs/full_reference_both_databases.homer"
    add_check(checks, "candidate_database.motif_count", count_motifs(candidate_db), EXPECTED_CANDIDATE_MOTIFS)
    add_check(checks, "candidate_database.jaspar_motif_count", count_motifs(candidate_jaspar), 21)
    add_check(checks, "candidate_database.cisbp_motif_count", count_motifs(candidate_cisbp), 26)
    add_check(checks, "full_database.motif_count", count_motifs(full_db), EXPECTED_FULL_MOTIFS)

    validate_result_set(
        luad,
        "anchor_candidate",
        luad / "audit/motif_equivalence/candidate_anchor_equivalent_motif_receipt.json",
        luad / "results/motif_equivalence/candidate_anchor_equivalent",
        EXPECTED_CANDIDATE_MOTIFS,
        checks,
    )
    anchor_receipt = read_json(luad / "audit/motif_equivalence/candidate_anchor_equivalent_motif_receipt.json")
    add_check(
        checks,
        "anchor_candidate.database_selection_rule",
        anchor_receipt["tf_database_selection_rule"],
        "JASPAR2024_if_available_else_CIS-BP2.00",
    )
    add_check(checks, "anchor_candidate.homer_result_root_count", len(anchor_receipt["homer_result_roots"]), 2)
    token_rows = read_tsv(luad / "audit/motif_equivalence/candidate_motif_token_map.tsv")
    sources_by_tf: dict[str, set[str]] = {}
    for row in token_rows:
        sources_by_tf.setdefault(row["TF"], set()).add(row["motif_token"].split("|", 1)[0])
    expected_source = {
        tf: "JASPAR2024" if "JASPAR2024" in sources else "CIS-BP2.00"
        for tf, sources in sources_by_tf.items()
    }
    anchor_best = read_tsv(luad / "results/motif_equivalence/candidate_anchor_equivalent/harmonized_sample_tf_best_motif.tsv")
    source_violations = [
        row
        for row in anchor_best
        if row["motif_token"] != "NA"
        and not row["motif_token"].startswith(f"{expected_source[row['TF']]}|")
    ]
    add_check(checks, "anchor_candidate.preferred_source_token_violations", len(source_violations), 0)
    validate_result_set(
        luad,
        "combined_candidate",
        luad / "audit/motif_equivalence/candidate_database_motif_receipt.json",
        luad / "results/motif_equivalence/candidate_database",
        EXPECTED_CANDIDATE_MOTIFS,
        checks,
    )
    validate_result_set(
        luad,
        "full",
        luad / "audit/motif_equivalence/harmonized_motif_receipt.json",
        luad / "results/motif_equivalence",
        EXPECTED_FULL_MOTIFS,
        checks,
    )

    final_receipt = read_json(luad / "audit/motif_equivalence/final_motif_equivalence_receipt.json")
    add_check(checks, "compiled.tnbc_recovered", final_receipt["tnbc_author_result_recovered"], True)
    failed = [row for row in checks if row["pass"] != "TRUE"]
    out_dir = luad / "audit/motif_equivalence"
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "motif_equivalence_independent_checks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "observed", "expected", "pass"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "check_count": len(checks),
        "passed_check_count": len(checks) - len(failed),
        "failed_check_count": len(failed),
        "all_checks_pass": not failed,
        "failed_checks": [row["check"] for row in failed],
        "scope": "independent structural and arithmetic validation; no biological pass/fail threshold added",
    }
    (out_dir / "motif_equivalence_independent_validation_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
