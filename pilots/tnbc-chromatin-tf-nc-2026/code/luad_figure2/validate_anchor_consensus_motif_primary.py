#!/usr/bin/env python3
"""Independent integrity checks for the primary anchor-consensus motif run."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_SYSTEMS = {"patient": 22, "PDX": 13, "cell_line": 19}
EXPECTED_CONTROLS = {"NKX2-1", "FOXA1", "FOXA2", "CEBPA", "CEBPB", "GATA6", "HNF4A", "ELF3"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "pass", "detail"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bed_stats(path: Path) -> tuple[int, bool, bool, bool]:
    count = 0
    widths_ok = True
    canonical_ok = True
    unique_ok = True
    seen: set[tuple[str, int, int]] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            count += 1
            widths_ok &= end - start == 200
            canonical_ok &= chrom == "chrX" or (
                chrom.startswith("chr") and chrom[3:].isdigit() and 1 <= int(chrom[3:]) <= 22
            )
            key = (chrom, start, end)
            if key in seen:
                unique_ok = False
            seen.add(key)
    return count, widths_ok, canonical_ok, unique_ok


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    input_root = root / "data/processed/motif_primary/anchor_consensus_200bp"
    manifest_path = input_root / "anchor_consensus_motif_manifest.tsv"
    background_path = input_root / "lung_accessible_consensus_background.bed"
    input_receipt_path = root / "audit/motif_primary/anchor_consensus_input_receipt.json"
    analysis_receipt_path = root / "audit/motif_primary/anchor_consensus_author_receipt.json"
    result_root = root / "results/motif_primary/anchor_consensus_author"
    triple_path = result_root / "anchor_consensus_hc_tf_triple_system_summary.tsv"
    system_path = result_root / "anchor_consensus_system_tf_motif_summary.tsv"
    best_path = result_root / "anchor_consensus_sample_tf_best_motif.tsv"
    control_path = result_root / "luad_lineage_positive_control_summary.tsv"
    required = [
        manifest_path, background_path, input_receipt_path, analysis_receipt_path,
        triple_path, system_path, best_path, control_path,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing primary motif artifacts: " + "; ".join(missing))

    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"check": name, "pass": str(bool(condition)).upper(), "detail": detail})

    manifest = read_tsv(manifest_path)
    system_counts = Counter(row["system"] for row in manifest)
    check("manifest_total_54", len(manifest) == 54, f"observed={len(manifest)}")
    for system, expected in EXPECTED_SYSTEMS.items():
        check(f"manifest_{system}_{expected}", system_counts[system] == expected, f"observed={system_counts[system]}")
    check(
        "all_samples_keep_all_usable_peaks",
        all(int(row["matched_peak_count"]) == int(row["harmonized_200bp_unique_peaks"]) for row in manifest),
        "no peak-count downsampling allowed in primary analysis",
    )
    check(
        "testability_tracks_nonempty_canonical_peaks",
        all((row["motif_matchable"] == "TRUE") == (int(row["matched_peak_count"]) > 0) for row in manifest),
        "zero-peak samples remain explicit in the fixed denominator",
    )

    for row in manifest:
        target = Path(row["matched_target_path"])
        count, widths_ok, canonical_ok, unique_ok = bed_stats(target)
        expected = int(row["harmonized_200bp_unique_peaks"])
        label = f"{row['system']}:{row['sample_id']}"
        check(f"target_count:{label}", count == expected, f"observed={count}; expected={expected}")
        check(f"target_width_200:{label}", widths_ok, str(target))
        check(f"target_canonical:{label}", canonical_ok, str(target))
        check(f"target_unique:{label}", unique_ok, str(target))

    background_count, background_widths, background_canonical, background_unique = bed_stats(background_path)
    check("background_nonempty", background_count > 0, f"observed={background_count}")
    check("background_width_200", background_widths, str(background_path))
    check("background_canonical", background_canonical, str(background_path))
    check("background_unique", background_unique, str(background_path))

    input_receipt = json.loads(input_receipt_path.read_text(encoding="utf-8"))
    check("input_role_primary", input_receipt.get("analysis_role") == "primary_TNBC_anchor_analog_all_peaks_shared_consensus_background", str(input_receipt.get("analysis_role")))
    check("input_no_manual_peak_matching", input_receipt.get("manual_target_peak_count_matching") is False, str(input_receipt.get("manual_target_peak_count_matching")))
    check("input_background_count_matches", input_receipt.get("shared_background_peak_count") == background_count, f"receipt={input_receipt.get('shared_background_peak_count')}; bed={background_count}")
    check("input_same_background_all_systems", input_receipt.get("same_background_file_reused_for_all_three_systems") is True, str(input_receipt.get("same_background_file_reused_for_all_three_systems")))

    testable_samples = sum(row["motif_matchable"] == "TRUE" for row in manifest)
    not_testable_samples = len(manifest) - testable_samples
    for source in ("JASPAR2024", "CIS-BP2.00"):
        source_root = root / "results/motif_primary/anchor_consensus_homer" / source
        known = list(source_root.glob("*/*/knownResults.txt"))
        not_testable = list(source_root.glob("*/*/.not_testable"))
        check(f"{source}_known_results_testable_samples", len(known) == testable_samples, f"observed={len(known)}; expected={testable_samples}")
        check(f"{source}_explicit_not_testable_samples", len(not_testable) == not_testable_samples, f"observed={len(not_testable)}; expected={not_testable_samples}")

    analysis_receipt = json.loads(analysis_receipt_path.read_text(encoding="utf-8"))
    check("author_tf_source_priority", analysis_receipt.get("tf_database_selection_rule") == "JASPAR2024_if_available_else_CIS-BP2.00", str(analysis_receipt.get("tf_database_selection_rule")))
    check("full_separate_result_roots", len(analysis_receipt.get("homer_result_roots", [])) == 2, str(analysis_receipt.get("homer_result_roots")))
    check("motif_testable_hc_tf_20", analysis_receipt.get("motif_testable_luad_hc_tf_count") == 20, str(analysis_receipt.get("motif_testable_luad_hc_tf_count")))

    triple = read_tsv(triple_path)
    check("triple_table_20_hc_tfs", len(triple) == 20 and len({row["TF"] for row in triple}) == 20, f"rows={len(triple)}")
    system_rows = read_tsv(system_path)
    check("system_table_three_systems", {row["system"] for row in system_rows} == set(EXPECTED_SYSTEMS), str(sorted({row["system"] for row in system_rows})))
    controls = read_tsv(control_path)
    check("eight_controls_three_systems", len(controls) == 24 and {row["TF"] for row in controls} == EXPECTED_CONTROLS, f"rows={len(controls)}; controls={sorted({row['TF'] for row in controls})}")

    best = read_tsv(best_path)
    sources_by_tf: dict[str, set[str]] = defaultdict(set)
    for row in best:
        token = row["motif_token"]
        if token != "NA":
            sources_by_tf[row["TF"]].add(token.split("|", 1)[0])
    mixed = {tf: sorted(sources) for tf, sources in sources_by_tf.items() if len(sources) > 1}
    check("one_preferred_database_source_per_tf", not mixed, json.dumps(mixed, sort_keys=True))

    output_checks = root / "audit/motif_primary/anchor_consensus_independent_checks.tsv"
    write_tsv(output_checks, checks)
    failures = [row for row in checks if row["pass"] != "TRUE"]
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if not failures else "FAIL",
        "checks_total": len(checks),
        "checks_passed": len(checks) - len(failures),
        "checks_failed": len(failures),
        "failed_checks": [row["check"] for row in failures],
        "triple_hc_tf_count": sum(row["triple_system_fixed_original_denominator"] == "TRUE" for row in triple),
        "triple_hc_tfs": [row["TF"] for row in triple if row["triple_system_fixed_original_denominator"] == "TRUE"],
        "positive_control_support_by_system": {
            system: sum(row["system"] == system and row["support_fixed_original_denominator"] == "TRUE" for row in controls)
            for system in EXPECTED_SYSTEMS
        },
    }
    receipt_path = root / "audit/motif_primary/anchor_consensus_independent_validation_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
