#!/usr/bin/env python3
"""Combine all predeclared, completely processed samples into the Figure 2 manifest."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    rows: list[dict[str, str]] = []

    patient_rows = read_tsv(root / "audit/tcga_luad_accessibility/tcga_luad_peak_file_manifest.tsv")
    for row in patient_rows:
        if row["definition"] != "cpm1_both_reps":
            continue
        rows.append(
            {
                "system": "patient",
                "sample_id": row["sample_id"],
                "sample_slug": row["sample_id"],
                "peak_path": row["peak_path"],
                "bam_path": "NA",
                "peak_definition": "GDC_open_counts_CPM1_both_technical_replicates",
                "idr_status": "not_rerun_open_fixed_peak_matrix",
            }
        )

    raw_rows = read_tsv(root / "audit/raw_atac_qc/figure2_analysis_raw_atac_samples.tsv")
    for row in raw_rows:
        rows.append(
            {
                "system": row["system"],
                "sample_id": row["sample_id"],
                "sample_slug": row["sample_slug"],
                "peak_path": row["narrowpeak"],
                "bam_path": row["filtered_bam"],
                "peak_definition": "MACS2_q0.01_anchor_aligned",
                "idr_status": "not_applicable_single_public_library",
            }
        )
    rows.sort(key=lambda row: ({"patient": 0, "PDX": 1, "cell_line": 2}[row["system"]], row["sample_id"].lower()))
    counts = {system: sum(row["system"] == system for row in rows) for system in ("patient", "PDX", "cell_line")}
    expected = {"patient": 22, "PDX": 13, "cell_line": 19}
    if counts != expected:
        raise RuntimeError(f"Figure 2 processing is incomplete: observed {counts}, expected {expected}")
    for row in rows:
        path = Path(row["peak_path"])
        if not path.is_file():
            raise RuntimeError(f"Missing peak file: {path}")

    out_dir = root / "audit/manifests/figure2_atomic"
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / "figure2_atomic_sample_manifest.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "counts": counts,
        "predeclared_counts": expected,
        "post_sequencing_read_or_peak_threshold_used_for_exclusion": False,
        "zero_peak_samples_retained": [
            row["sample_id"] for row in rows if Path(row["peak_path"]).stat().st_size == 0
        ],
        "one_row_per_independent_biological_sample": True,
        "manifest": str(output),
        "shared_by": ["Figure2", "SupplementaryFigure3", "SupplementaryFigure4A-C"],
    }
    (out_dir / "figure2_atomic_sample_manifest_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
