#!/usr/bin/env python3
"""Verify the atomic Figure 2 + Supplementary Figure 3 + 4A-C delivery."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


TF_SHA = "6f449276a97b8c66cbf93bd2a4ca6d51095fa426015dc9046a7c084431d5037b"
TF_COUNT = 158
PRIMARY = "PRIMARY_anchor_promoter_tcga_cpm1_both"


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path):
    digest = hashlib.sha256(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    project = args.project_root.resolve()
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, observed: object, expected: object):
        checks.append({"check": name, "status": "PASS" if condition else "FAIL", "observed": observed, "expected": expected})

    tf_path = project / "pilots/tnbc-chromatin-tf-nc-2026/execution/luad-figure1/results/tables/tcga_specific_TFs.tsv"
    check("frozen_tf_parent_sha", sha256(tf_path) == TF_SHA, sha256(tf_path), TF_SHA)
    frozen_luad_path = root / "audit/promoter_gate/frozen_figure2_luad_tf_input.tsv"
    frozen_luad_rows = read_tsv(frozen_luad_path)
    promoter_receipt = json.loads((root / "audit/promoter_gate/promoter_gate_receipt.json").read_text(encoding="utf-8"))
    check("frozen_luad_tf_sha", sha256(frozen_luad_path) == promoter_receipt["frozen_tf_sha256"],
          sha256(frozen_luad_path), promoter_receipt["frozen_tf_sha256"])
    check("frozen_luad_tf_count", len(frozen_luad_rows) == TF_COUNT, len(frozen_luad_rows), TF_COUNT)
    check("frozen_input_is_all_luad_discovery",
          all(row["discovery_category"] == "LUAD" for row in frozen_luad_rows),
          sum(row["discovery_category"] == "LUAD" for row in frozen_luad_rows), TF_COUNT)
    manifest = read_tsv(root / "audit/manifests/figure2_atomic/figure2_atomic_sample_manifest.tsv")
    counts = {system: sum(row["system"] == system for row in manifest) for system in ("patient", "PDX", "cell_line")}
    raw_receipt = json.loads((root / "audit/raw_atac_qc/figure2_raw_atac_qc_receipt.json").read_text(encoding="utf-8"))
    expected_counts = {
        "patient": 22,
        "PDX": int(raw_receipt["systems"]["PDX"]["analysis_included"]),
        "cell_line": int(raw_receipt["systems"]["cell_line"]["analysis_included"]),
    }
    check("raw_data_completion", raw_receipt["raw_data_completion"] == "COMPLETE",
          raw_receipt["raw_data_completion"], "COMPLETE")
    check("atomic_manifest_counts", counts == expected_counts, counts, expected_counts)
    total_n = sum(counts.values())
    check("unique_biological_samples", len({(row["system"], row["sample_id"]) for row in manifest}) == total_n,
          len({(row["system"], row["sample_id"]) for row in manifest}), total_n)

    promoter_sample = read_tsv(root / "results/promoter_gate/sample_tf_promoter_accessibility.tsv")
    primary_rows = [row for row in promoter_sample if row["analysis_definition"] == PRIMARY]
    check("primary_promoter_matrix_rows", len(primary_rows) == total_n * TF_COUNT, len(primary_rows), total_n * TF_COUNT)
    check("primary_promoter_unique_cells", len({(row["system"], row["sample_id"], row["TF"]) for row in primary_rows}) == total_n * TF_COUNT,
          len({(row["system"], row["sample_id"], row["TF"]) for row in primary_rows}), total_n * TF_COUNT)
    promoter_tf = [
        row for row in read_tsv(root / "results/promoter_gate/tf_promoter_gate_summary.tsv")
        if row["analysis_definition"] == PRIMARY
    ]
    hc_rows = [row for row in promoter_tf if row["HC_TF_promoter_activity_definition"] == "TRUE"]
    check("hc_rule_is_triple_open_not_all_three_negative",
          all(row["triple_system_promoter_accessible"] == "TRUE" and
              row["all_three_system_mean_NES_negative"] == "FALSE" for row in hc_rows),
          len(hc_rows), "every HC-TF satisfies the anchor rule")

    inventory = read_tsv(root / "results/motif_gate/hc_tf_motif_inventory.tsv")
    testable = sum(row["motif_testable"] == "TRUE" for row in inventory)
    motif_best = read_tsv(root / "results/motif_gate/sample_tf_best_motif_results.tsv")
    check("motif_best_matrix_rows", len(motif_best) == total_n * testable, len(motif_best), total_n * testable)
    check("motif_best_unique_cells", len({(row["system"], row["sample_id"], row["TF"]) for row in motif_best}) == total_n * testable,
          len({(row["system"], row["sample_id"], row["TF"]) for row in motif_best}), total_n * testable)

    sat = read_tsv(root / "results/supplementary_figure3/saturation_1000_permutations.tsv")
    expected_sat = 1000 * total_n
    check("saturation_rows", len(sat) == expected_sat, len(sat), expected_sat)
    annotations = read_tsv(root / "results/supplementary_figure3/peak_genomic_annotation.tsv")
    check("genomic_annotation_rows", len(annotations) == total_n * 4, len(annotations), total_n * 4)

    ataqv_rows = read_tsv(root / "audit/ataqv/figure2_ataqv_metrics.tsv")
    expected_raw_n = expected_counts["PDX"] + expected_counts["cell_line"]
    check("ataqv_metrics_rows", len(ataqv_rows) == expected_raw_n, len(ataqv_rows), expected_raw_n)
    check("ataqv_unique_raw_samples",
          len({(row["system"], row["sample_id"]) for row in ataqv_rows}) == len(ataqv_rows),
          len({(row["system"], row["sample_id"]) for row in ataqv_rows}), len(ataqv_rows))
    check("ataqv_all_json_parsed", all(row["parse_status"] == "PASS" for row in ataqv_rows),
          sum(row["parse_status"] == "PASS" for row in ataqv_rows), len(ataqv_rows))
    tn5_rows = read_tsv(root / "audit/tn5_tss/figure2_tn5_tss_enrichment.tsv")
    check("tn5_tss_metrics_rows", len(tn5_rows) == expected_raw_n, len(tn5_rows), expected_raw_n)
    check("tn5_tss_all_profiles_complete",
          all(row["status"] == "PASS" and int(row["profile_points"]) == 2001 for row in tn5_rows),
          sum(row["status"] == "PASS" and int(row["profile_points"]) == 2001 for row in tn5_rows), expected_raw_n)

    required_stems = (
        "Figure2_complete",
        "SupplementaryFigure3_complete",
        "SupplementaryFigure4A-C_complete",
    )
    for stem in required_stems:
        for extension in ("pdf", "png"):
            path = root / f"results/figures/{stem}.{extension}"
            check(f"artifact_{stem}_{extension}", path.is_file() and path.stat().st_size > 10_000,
                  path.stat().st_size if path.is_file() else 0, ">10000 bytes")

    final = json.loads((root / "audit/final_gate/figure2_final_verdict.json").read_text(encoding="utf-8"))
    check("final_status", final["final_verdict"] == "ORIGINAL_STYLE_ANALYSIS_COMPLETE",
          final["final_verdict"], "ORIGINAL_STYLE_ANALYSIS_COMPLETE")
    check("no_invented_biological_threshold", final["biological_pass_fail_threshold_applied"] is False,
          final["biological_pass_fail_threshold_applied"], False)
    check("no_automatic_Figure3_decision", final.get("automatic_Figure3_decision") is None,
          final.get("automatic_Figure3_decision"), None)
    overall = "PASS" if all(row["status"] == "PASS" for row in checks) else "FAIL"
    out_dir = root / "audit/final_gate"
    with (out_dir / "figure2_atomic_verification.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verification": overall,
        "checks_passed": sum(row["status"] == "PASS" for row in checks),
        "checks_total": len(checks),
        "final_gate_verdict": final["final_verdict"],
    }
    (out_dir / "figure2_atomic_verification_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if overall != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
