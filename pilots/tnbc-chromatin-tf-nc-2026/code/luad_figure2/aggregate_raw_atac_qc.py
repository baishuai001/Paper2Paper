#!/usr/bin/env python3
"""Freeze raw-ATAC processing completeness and report sample QC diagnostics.

The TNBC anchor did not discard sequenced libraries with post hoc minimum-read,
minimum-peak, or cohort-failure-fraction rules.  Accordingly, the LUAD primary
analysis includes every predeclared biological model whose expected artifacts
are complete and manifest-consistent.  Read/peak/FRiP/TSS metrics remain in the
audit tables, but do not act as biological stopping rules.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED = {"PDX": 13, "cell_line": 19}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    manifest = root / "audit/manifests/raw_atac/figure2_raw_atac_runs.tsv"
    rows = read_tsv(manifest)
    if len(rows) != sum(EXPECTED.values()):
        raise RuntimeError(f"Expected 32 raw runs, observed {len(rows)}")

    output_rows: list[dict[str, object]] = []
    for row in rows:
        system = row["system"]
        slug = row["sample_slug"]
        if system not in EXPECTED:
            raise RuntimeError(f"Unexpected system in manifest: {system}")
        qc_path = root / "audit/raw_atac_qc" / system / f"{slug}.tsv"
        out_dir = root / "data/processed/raw_atac" / system / slug
        bam_path = out_dir / f"{slug}.filtered.bam"
        peak_path = out_dir / "peaks" / f"{slug}_peaks.narrowPeak"
        complete_path = out_dir / ".complete"
        qc_record: dict[str, str] = {}
        parse_error = ""
        if qc_path.is_file():
            parsed = read_tsv(qc_path)
            if len(parsed) == 1:
                qc_record = parsed[0]
            else:
                parse_error = f"qc_rows={len(parsed)}"
        else:
            parse_error = "missing_qc"

        artifacts_ok = all(path.is_file() and path.stat().st_size > 0 for path in (bam_path, peak_path, complete_path))
        manifest_match = bool(qc_record) and all(
            qc_record.get(key, "") == row[manifest_key]
            for key, manifest_key in (
                ("system", "system"),
                ("sample_id", "sample_id"),
                ("sample_slug", "sample_slug"),
                ("run", "run"),
                ("library_layout", "layout"),
            )
        )
        attempted = bool(qc_record)
        analysis_included = attempted and not parse_error and artifacts_ok and manifest_match
        exclusion_reason = ""
        if not analysis_included:
            if parse_error:
                exclusion_reason = parse_error
            elif not artifacts_ok:
                exclusion_reason = "missing_or_empty_final_artifact"
            elif not manifest_match:
                exclusion_reason = "qc_manifest_mismatch"

        diagnostic_qc = qc_record.get("diagnostic_qc_status", qc_record.get("hard_qc_status", "MISSING"))

        output_rows.append(
            {
                **row,
                "filtered_human_reads": qc_record.get("filtered_human_reads", ""),
                "filtered_human_units": qc_record.get("filtered_human_units", ""),
                "peak_count": qc_record.get("peak_count", ""),
                "frip": qc_record.get("frip", ""),
                "worker_diagnostic_qc_status": diagnostic_qc,
                "processing_status": "COMPLETED" if analysis_included else "INCOMPLETE",
                "final_artifacts_ok": str(artifacts_ok).upper(),
                "qc_manifest_match": str(manifest_match).upper(),
                "analysis_included": str(analysis_included).upper(),
                "exclusion_reason": exclusion_reason,
                "filtered_bam": str(bam_path),
                "narrowpeak": str(peak_path),
                "qc_record": str(qc_path),
            }
        )

    output_rows.sort(key=lambda item: (str(item["system"]), str(item["sample_id"]).lower()))
    columns = list(output_rows[0])
    audit_dir = root / "audit/raw_atac_qc"
    all_path = audit_dir / "figure2_raw_atac_qc.tsv"
    analysis_path = audit_dir / "figure2_analysis_raw_atac_samples.tsv"
    write_tsv(all_path, output_rows, columns)
    analysis_rows = [row for row in output_rows if row["analysis_included"] == "TRUE"]
    write_tsv(analysis_path, analysis_rows, columns)

    counts: dict[str, dict[str, object]] = {}
    completion_errors: list[str] = []
    for system in EXPECTED:
        system_rows = [row for row in output_rows if row["system"] == system]
        included_rows = [row for row in system_rows if row["analysis_included"] == "TRUE"]
        diagnostic_failures = [
            row for row in included_rows if row["worker_diagnostic_qc_status"] == "BELOW_DIAGNOSTIC_REFERENCE"
            or row["worker_diagnostic_qc_status"] == "FAIL"
        ]
        complete = len(included_rows) == EXPECTED[system]
        evaluation_status = "COMPLETE" if complete else "INCOMPLETE"
        counts[system] = {
            "expected": EXPECTED[system],
            "manifest_count": len(system_rows),
            "processed_with_complete_artifacts": len(included_rows),
            "analysis_included": len(included_rows),
            "diagnostic_reference_failures": len(diagnostic_failures),
            "evaluation_status": evaluation_status,
        }
        if len(system_rows) != EXPECTED[system]:
            completion_errors.append(f"{system}:manifest_count")
        if not complete:
            completion_errors.append(f"{system}:processing_incomplete")

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest),
        "manifest_sha256": sha256(manifest),
        "sample_qc_table": str(all_path),
        "analysis_sample_table": str(analysis_path),
        "systems": counts,
        "raw_data_completion": "COMPLETE" if not completion_errors else "INCOMPLETE",
        "completion_errors": completion_errors,
        "diagnostic_reference_only": {
            "filtered_human_units": 1_000_000,
            "MACS2_peak_count": 10_000,
            "used_to_exclude_samples": False,
        },
        "cohort_failure_fraction_rule_removed": True,
        "patient_open_matrix_gate_is_evaluated_separately": True,
    }
    receipt_path = audit_dir / "figure2_raw_atac_qc_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if completion_errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
