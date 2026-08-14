#!/usr/bin/env python3
"""Freeze sample-level PDX/cell-line ATAC QC and evaluate the raw-data gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


EXPECTED = {"PDX": 13, "cell_line": 19}
MIN_QUALIFIED = {"PDX": 10, "cell_line": 10}
MAX_FAILURE_FRACTION = 0.20


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
    parser.add_argument(
        "--early-stop-locked-system",
        choices=tuple(EXPECTED),
        default=None,
        help=(
            "A fully evaluated system that already makes FAIL_DATA irreversible. "
            "Unattempted samples in the other system are then recorded as not run, "
            "rather than misclassified as QC failures."
        ),
    )
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
        early_stop_path = root / "audit/raw_atac_qc" / system / f"{slug}.early_stop.json"
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
        worker_pass = qc_record.get("hard_qc_status") == "PASS"
        attempted = bool(qc_record)
        aborted_after_lock = False
        if early_stop_path.is_file():
            early_stop_record = json.loads(early_stop_path.read_text(encoding="utf-8"))
            aborted_after_lock = (
                early_stop_record.get("status") == "ABORTED_AFTER_LOCKED_DATA_GATE"
                and early_stop_record.get("counted_as_qc_failure") is False
            )
        qualified = attempted and not parse_error and artifacts_ok and manifest_match and worker_pass
        failure_reason = ""
        if not qualified:
            if aborted_after_lock:
                failure_reason = "aborted_after_locked_data_gate"
            elif not attempted and args.early_stop_locked_system is not None and system != args.early_stop_locked_system:
                failure_reason = "not_run_after_locked_data_gate"
            elif parse_error:
                failure_reason = parse_error
            elif not artifacts_ok:
                failure_reason = "missing_or_empty_final_artifact"
            elif not manifest_match:
                failure_reason = "qc_manifest_mismatch"
            else:
                failure_reason = "worker_hard_qc_fail"

        output_rows.append(
            {
                **row,
                "filtered_human_reads": qc_record.get("filtered_human_reads", ""),
                "filtered_human_units": qc_record.get("filtered_human_units", ""),
                "peak_count": qc_record.get("peak_count", ""),
                "frip": qc_record.get("frip", ""),
                "worker_hard_qc_status": qc_record.get("hard_qc_status", "MISSING"),
                "processing_status": (
                    "COMPLETED"
                    if attempted
                    else "ABORTED_AFTER_LOCKED_DATA_GATE"
                    if aborted_after_lock
                    else "NOT_RUN_AFTER_LOCKED_DATA_GATE"
                    if args.early_stop_locked_system is not None and system != args.early_stop_locked_system
                    else "MISSING_UNEXPECTEDLY"
                ),
                "final_artifacts_ok": str(artifacts_ok).upper(),
                "qc_manifest_match": str(manifest_match).upper(),
                "figure2_qualified": str(qualified).upper(),
                "failure_reason": failure_reason,
                "filtered_bam": str(bam_path),
                "narrowpeak": str(peak_path),
                "qc_record": str(qc_path),
                "early_stop_receipt": str(early_stop_path) if early_stop_path.is_file() else "",
            }
        )

    output_rows.sort(key=lambda item: (str(item["system"]), str(item["sample_id"]).lower()))
    columns = list(output_rows[0])
    audit_dir = root / "audit/raw_atac_qc"
    all_path = audit_dir / "figure2_raw_atac_qc.tsv"
    qualified_path = audit_dir / "figure2_qualified_raw_atac_samples.tsv"
    write_tsv(all_path, output_rows, columns)
    qualified_rows = [row for row in output_rows if row["figure2_qualified"] == "TRUE"]
    write_tsv(qualified_path, qualified_rows, columns)

    counts: dict[str, dict[str, object]] = {}
    reasons: list[str] = []
    for system in EXPECTED:
        system_rows = [row for row in output_rows if row["system"] == system]
        attempted_rows = [row for row in system_rows if row["processing_status"] == "COMPLETED"]
        aborted_rows = [row for row in system_rows if row["processing_status"] == "ABORTED_AFTER_LOCKED_DATA_GATE"]
        not_run_rows = [row for row in system_rows if row["processing_status"] == "NOT_RUN_AFTER_LOCKED_DATA_GATE"]
        n_qualified = sum(row["figure2_qualified"] == "TRUE" for row in system_rows)
        n_failed = len(attempted_rows) - n_qualified
        failure_fraction = n_failed / len(attempted_rows) if attempted_rows else None
        complete = len(attempted_rows) == EXPECTED[system]
        if complete:
            evaluation_status = "COMPLETE"
        elif args.early_stop_locked_system is not None and system != args.early_stop_locked_system:
            evaluation_status = "NOT_FULLY_EVALUATED_AFTER_EARLY_STOP"
        else:
            evaluation_status = "INCOMPLETE_UNEXPECTEDLY"
        counts[system] = {
            "expected": EXPECTED[system],
            "manifest_count": len(system_rows),
            "observed": len(attempted_rows),
            "attempted": len(attempted_rows),
            "aborted_after_locked_gate": len(aborted_rows),
            "not_run": len(not_run_rows),
            "qualified": n_qualified,
            "failed": n_failed,
            "failure_fraction": failure_fraction,
            "evaluation_status": evaluation_status,
            "minimum_qualified": MIN_QUALIFIED[system],
            "maximum_failure_fraction": MAX_FAILURE_FRACTION,
        }
        if len(system_rows) != EXPECTED[system]:
            reasons.append(f"{system}:manifest_count")
        if complete:
            if n_qualified < MIN_QUALIFIED[system]:
                reasons.append(f"{system}:qualified_below_minimum")
            if failure_fraction is not None and failure_fraction > MAX_FAILURE_FRACTION:
                reasons.append(f"{system}:failure_fraction_above_0.20")
        elif evaluation_status == "INCOMPLETE_UNEXPECTEDLY":
            reasons.append(f"{system}:processing_incomplete")

    if args.early_stop_locked_system is not None:
        locked = counts[args.early_stop_locked_system]
        locked_failure = (
            locked["evaluation_status"] == "COMPLETE"
            and (
                int(locked["qualified"]) < MIN_QUALIFIED[args.early_stop_locked_system]
                or (
                    locked["failure_fraction"] is not None
                    and float(locked["failure_fraction"]) > MAX_FAILURE_FRACTION
                )
            )
        )
        if not locked_failure:
            raise RuntimeError(
                "--early-stop-locked-system was supplied, but that system does not "
                "have a complete, irreversible data-gate failure"
            )

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest),
        "manifest_sha256": sha256(manifest),
        "sample_qc_table": str(all_path),
        "qualified_sample_table": str(qualified_path),
        "systems": counts,
        "early_stop_locked_system": args.early_stop_locked_system,
        "unattempted_samples_are_not_counted_as_qc_failures": True,
        "raw_data_gate": "PASS" if not reasons else "FAIL_DATA",
        "failure_reasons": reasons,
        "patient_open_matrix_gate_is_evaluated_separately": True,
    }
    receipt_path = audit_dir / "figure2_raw_atac_qc_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
