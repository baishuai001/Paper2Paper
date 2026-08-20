#!/usr/bin/env python3
"""Validate and summarize B0-B3 shared-background motif sensitivity."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


SYSTEMS = ("patient", "PDX", "cell_line")
VARIANTS = ("B0", "B1", "B2", "B3")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    input_receipt_path = root / "audit/motif_background_sensitivity/background_input_receipt.json"
    input_receipt = json.loads(input_receipt_path.read_text(encoding="utf-8"))

    result_roots = {
        "B0": root / "results/motif_primary/anchor_consensus_author",
        **{
            variant: root / f"results/motif_background_sensitivity/{variant}/parsed"
            for variant in ("B1", "B2", "B3")
        },
    }
    prefix = {"B0": "anchor_consensus", "B1": "B1", "B2": "B2", "B3": "B3"}
    background_info = {
        "B0": input_receipt["immutable_primary_background"],
        **input_receipt["variants"],
    }
    variant_rows: list[dict[str, object]] = []
    triple_sets: dict[str, set[str]] = {}
    all_hc_tfs: set[str] = set()
    validation_checks: list[dict[str, object]] = []

    for variant in VARIANTS:
        result_root = result_roots[variant]
        triple_path = result_root / f"{prefix[variant]}_hc_tf_triple_system_summary.tsv"
        system_path = result_root / f"{prefix[variant]}_system_tf_motif_summary.tsv"
        control_path = result_root / "luad_lineage_positive_control_summary.tsv"
        for path in (triple_path, system_path, control_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        triple = read_tsv(triple_path)
        system = read_tsv(system_path)
        controls = read_tsv(control_path)
        all_hc_tfs.update(row["TF"] for row in triple)
        supported = {
            row["TF"] for row in triple
            if row["triple_system_fixed_original_denominator"] == "TRUE"
        }
        triple_sets[variant] = supported
        all_systems_ok = {row["system"] for row in system} == set(SYSTEMS)
        control_systems_ok = {row["system"] for row in controls} == set(SYSTEMS)
        validation_checks.extend([
            {"check": f"{variant}_20_motif_testable_hc_tfs", "pass": len(triple) == 20, "detail": len(triple)},
            {"check": f"{variant}_three_systems", "pass": all_systems_ok, "detail": sorted({row['system'] for row in system})},
            {"check": f"{variant}_controls_three_systems", "pass": control_systems_ok, "detail": sorted({row['system'] for row in controls})},
        ])
        info = background_info[variant]
        row: dict[str, object] = {
            "variant": variant,
            "background_definition": info["definition"],
            "background_peak_count": info["background_peak_count"],
            "background_sha256": info["background_sha256"],
            "triple_system_hc_tf_count": len(supported),
            "triple_system_hc_tfs": ";".join(sorted(supported)),
        }
        for arm in SYSTEMS:
            row[f"{arm}_hc_tf_support_count"] = sum(
                item["system"] == arm
                and item["TF_role"] in {"LUAD_HC_TF", "LUAD_HC_TF_AND_LINEAGE_POSITIVE_CONTROL"}
                and item["support_fixed_original_denominator"] == "TRUE"
                for item in system
            )
            row[f"{arm}_lineage_control_support_count"] = sum(
                item["system"] == arm and item["support_fixed_original_denominator"] == "TRUE"
                for item in controls
            )
        variant_rows.append(row)

    stability_rows = []
    for tf in sorted(all_hc_tfs):
        flags = {variant: tf in triple_sets[variant] for variant in VARIANTS}
        stability_rows.append({
            "TF": tf,
            **{f"{variant}_triple_support": str(flags[variant]).upper() for variant in VARIANTS},
            "backgrounds_with_triple_support": sum(flags.values()),
            "stable_all_four_backgrounds": str(all(flags.values())).upper(),
            "primary_B0_supported": str(flags["B0"]).upper(),
        })

    out_dir = root / "results/motif_background_sensitivity"
    write_tsv(
        out_dir / "background_variant_summary.tsv", variant_rows,
        [
            "variant", "background_definition", "background_peak_count", "background_sha256",
            "triple_system_hc_tf_count", "triple_system_hc_tfs",
            "patient_hc_tf_support_count", "patient_lineage_control_support_count",
            "PDX_hc_tf_support_count", "PDX_lineage_control_support_count",
            "cell_line_hc_tf_support_count", "cell_line_lineage_control_support_count",
        ],
    )
    write_tsv(
        out_dir / "hc_tf_background_stability.tsv", stability_rows,
        [
            "TF", "B0_triple_support", "B1_triple_support", "B2_triple_support",
            "B3_triple_support", "backgrounds_with_triple_support",
            "stable_all_four_backgrounds", "primary_B0_supported",
        ],
    )
    failures = [row for row in validation_checks if not row["pass"]]
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS" if not failures else "FAIL",
        "analysis_role": "sensitivity_only_primary_B0_unchanged",
        "variants": {
            row["variant"]: {
                "background_peak_count": row["background_peak_count"],
                "triple_system_hc_tf_count": row["triple_system_hc_tf_count"],
                "triple_system_hc_tfs": row["triple_system_hc_tfs"].split(";") if row["triple_system_hc_tfs"] else [],
            }
            for row in variant_rows
        },
        "hc_tfs_stable_in_all_four_backgrounds": [
            row["TF"] for row in stability_rows if row["stable_all_four_backgrounds"] == "TRUE"
        ],
        "checks": validation_checks,
        "failed_checks": failures,
        "interpretation_rule": "B1-B3 quantify background dependence and cannot replace the immutable B0 primary result.",
    }
    receipt_path = root / "audit/motif_background_sensitivity/final_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

