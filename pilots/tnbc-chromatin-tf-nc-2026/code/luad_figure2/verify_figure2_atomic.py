#!/usr/bin/env python3
"""Verify the atomic Figure 2 + Supplementary Figure 3 + 4A-C delivery."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


TF_SHA = "eace4fdf6509d01e8a9648b34b4dc3fdfb5955627c2c7a526c6ce1495a073844"
LUAD_TF_SHA = "64b4aa4e226c6e9142ecd560ec51ec1f49e0c9823d07e332757deeea3aa3c69b"
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

    tf_path = project / "pilots/tnbc-chromatin-tf-nc-2026/execution/luad-figure1/results/tables/externally_replicated_TFs.tsv"
    check("frozen_tf_parent_sha", sha256(tf_path) == TF_SHA, sha256(tf_path), TF_SHA)
    frozen_luad_path = root / "audit/promoter_gate/frozen_figure2_luad_tf_input.tsv"
    frozen_luad_rows = read_tsv(frozen_luad_path)
    check("frozen_luad_tf_sha", sha256(frozen_luad_path) == LUAD_TF_SHA, sha256(frozen_luad_path), LUAD_TF_SHA)
    check("frozen_luad_tf_count", len(frozen_luad_rows) == 97, len(frozen_luad_rows), 97)
    manifest = read_tsv(root / "audit/manifests/figure2_atomic/figure2_atomic_sample_manifest.tsv")
    counts = {system: sum(row["system"] == system for row in manifest) for system in ("patient", "PDX", "cell_line")}
    raw_receipt = json.loads((root / "audit/raw_atac_qc/figure2_raw_atac_qc_receipt.json").read_text(encoding="utf-8"))
    expected_counts = {
        "patient": 22,
        "PDX": int(raw_receipt["systems"]["PDX"]["qualified"]),
        "cell_line": int(raw_receipt["systems"]["cell_line"]["qualified"]),
    }
    check("raw_data_gate", raw_receipt["raw_data_gate"] == "PASS", raw_receipt["raw_data_gate"], "PASS")
    check("atomic_manifest_counts", counts == expected_counts, counts, expected_counts)
    total_n = sum(counts.values())
    check("unique_biological_samples", len({(row["system"], row["sample_id"]) for row in manifest}) == total_n,
          len({(row["system"], row["sample_id"]) for row in manifest}), total_n)

    promoter_sample = read_tsv(root / "results/promoter_gate/sample_tf_promoter_accessibility.tsv")
    primary_rows = [row for row in promoter_sample if row["analysis_definition"] == PRIMARY]
    check("primary_promoter_matrix_rows", len(primary_rows) == total_n * 97, len(primary_rows), total_n * 97)
    check("primary_promoter_unique_cells", len({(row["system"], row["sample_id"], row["TF"]) for row in primary_rows}) == total_n * 97,
          len({(row["system"], row["sample_id"], row["TF"]) for row in primary_rows}), total_n * 97)

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
    check("final_verdict_allowed", final["final_verdict"] in {"PASS", "FAIL_DATA", "FAIL_SIGNAL", "ANCHOR_PASS_ROBUSTNESS_FAIL"},
          final["final_verdict"], "allowed frozen status")
    check("explicit_stop", final["stop_after_this_gate"] is True, final["stop_after_this_gate"], True)
    check("no_automatic_Figure3_execution", final.get("continue_to_Figure3") is False,
          final.get("continue_to_Figure3"), False)
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
