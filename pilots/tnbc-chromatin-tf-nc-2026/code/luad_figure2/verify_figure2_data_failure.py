#!/usr/bin/env python3
"""Independently verify the terminal FAIL_DATA path and atomic non-generation."""

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
    raw = json.loads(
        (root / "audit/raw_atac_qc/figure2_raw_atac_qc_receipt.json").read_text(encoding="utf-8")
    )
    final = json.loads(
        (root / "audit/final_gate/figure2_final_verdict.json").read_text(encoding="utf-8")
    )
    bundle = json.loads(
        (root / "audit/final_gate/figure2_atomic_bundle_status.json").read_text(encoding="utf-8")
    )
    manifest = read_tsv(root / "audit/manifests/raw_atac/figure2_raw_atac_runs.tsv")
    qc = read_tsv(root / "audit/raw_atac_qc/figure2_raw_atac_qc.tsv")
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, observed: object, expected: object) -> None:
        checks.append(
            {
                "check": name,
                "status": "PASS" if condition else "FAIL",
                "observed": observed,
                "expected": expected,
            }
        )

    manifest_counts = {
        system: sum(row["system"] == system for row in manifest)
        for system in ("PDX", "cell_line")
    }
    check("raw_manifest_32", len(manifest) == 32, len(manifest), 32)
    check("raw_manifest_system_counts", manifest_counts == {"PDX": 13, "cell_line": 19}, manifest_counts, {"PDX": 13, "cell_line": 19})
    check("raw_gate_fail_data", raw["raw_data_gate"] == "FAIL_DATA", raw["raw_data_gate"], "FAIL_DATA")
    check("locked_system", raw.get("early_stop_locked_system") == "cell_line", raw.get("early_stop_locked_system"), "cell_line")

    cell = raw["systems"]["cell_line"]
    check("cell_line_evaluation_complete", cell["evaluation_status"] == "COMPLETE", cell["evaluation_status"], "COMPLETE")
    check("cell_line_all_19_attempted", int(cell["attempted"]) == 19, cell["attempted"], 19)
    check("cell_line_actual_failure_count", int(cell["failed"]) >= 4, cell["failed"], ">=4")
    check("cell_line_failure_fraction", float(cell["failure_fraction"]) > 0.20, cell["failure_fraction"], ">0.20")
    check(
        "cell_line_failure_reason",
        "cell_line:failure_fraction_above_0.20" in raw["failure_reasons"],
        raw["failure_reasons"],
        "contains cell_line:failure_fraction_above_0.20",
    )

    pdx = raw["systems"]["PDX"]
    check(
        "PDX_explicit_early_stop_status",
        pdx["evaluation_status"] == "NOT_FULLY_EVALUATED_AFTER_EARLY_STOP",
        pdx["evaluation_status"],
        "NOT_FULLY_EVALUATED_AFTER_EARLY_STOP",
    )
    check(
        "PDX_accounting",
        int(pdx["attempted"]) + int(pdx["aborted_after_locked_gate"]) + int(pdx["not_run"]) == 13,
        {
            "completed": pdx["attempted"],
            "aborted_after_locked_gate": pdx["aborted_after_locked_gate"],
            "not_run": pdx["not_run"],
        },
        "sum=13",
    )
    not_run_rows = [row for row in qc if row["processing_status"] == "NOT_RUN_AFTER_LOCKED_DATA_GATE"]
    check("not_run_row_count", len(not_run_rows) == int(pdx["not_run"]), len(not_run_rows), pdx["not_run"])
    check(
        "not_run_not_mislabelled_as_observed_qc_failure",
        all(row["failure_reason"] == "not_run_after_locked_data_gate" and row["worker_hard_qc_status"] == "MISSING" for row in not_run_rows),
        sum(row["failure_reason"] == "not_run_after_locked_data_gate" for row in not_run_rows),
        len(not_run_rows),
    )
    aborted_rows = [row for row in qc if row["processing_status"] == "ABORTED_AFTER_LOCKED_DATA_GATE"]
    check(
        "aborted_after_lock_row_count",
        len(aborted_rows) == int(pdx["aborted_after_locked_gate"]) == 1,
        len(aborted_rows),
        1,
    )
    check(
        "aborted_after_lock_not_mislabelled_as_qc_failure",
        all(
            row["failure_reason"] == "aborted_after_locked_data_gate"
            and row["worker_hard_qc_status"] == "MISSING"
            and bool(row["early_stop_receipt"])
            for row in aborted_rows
        ),
        sum(row["failure_reason"] == "aborted_after_locked_data_gate" for row in aborted_rows),
        len(aborted_rows),
    )

    check("final_verdict_fail_data", final["final_verdict"] == "FAIL_DATA", final["final_verdict"], "FAIL_DATA")
    check("explicit_stop", final["stop_after_this_gate"] is True, final["stop_after_this_gate"], True)
    check("no_figure3", final["continue_to_Figure3"] is False, final["continue_to_Figure3"], False)
    check("not_eligible_for_figure3", final["eligible_to_plan_Figure3"] is False, final["eligible_to_plan_Figure3"], False)
    check("atomic_bundle_status", bundle["status"] == "NOT_GENERATED_DUE_FAIL_DATA", bundle["status"], "NOT_GENERATED_DUE_FAIL_DATA")

    forbidden_stems = (
        "Figure2_complete",
        "SupplementaryFigure3_complete",
        "SupplementaryFigure4A-C_complete",
    )
    forbidden_artifacts = [
        root / f"results/figures/{stem}.{extension}"
        for stem in forbidden_stems
        for extension in ("pdf", "png")
    ]
    present_forbidden = [str(path) for path in forbidden_artifacts if path.exists()]
    check("atomic_signal_artifacts_absent", not present_forbidden, present_forbidden, [])
    forbidden_tables = (
        root / "results/promoter_gate/tf_promoter_gate_summary.tsv",
        root / "results/motif_gate/tf_motif_gate_summary.tsv",
    )
    present_tables = [str(path) for path in forbidden_tables if path.exists()]
    check("downstream_signal_tables_absent", not present_tables, present_tables, [])
    report_path = root / "reports/luad-figure2-gate-report.md"
    report_text = report_path.read_text(encoding="utf-8") if report_path.is_file() else ""
    check("terminal_report_exists", bool(report_text), len(report_text), ">0 characters")
    check("terminal_report_marks_atomic_non_generation", "NOT_GENERATED_DUE_FAIL_DATA" in report_text, "NOT_GENERATED_DUE_FAIL_DATA" in report_text, True)

    overall = "PASS" if all(row["status"] == "PASS" for row in checks) else "FAIL"
    out = root / "audit/final_gate"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "figure2_data_failure_verification.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(checks)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "verification": overall,
        "checks_passed": sum(row["status"] == "PASS" for row in checks),
        "checks_total": len(checks),
        "final_gate_verdict": final["final_verdict"],
        "atomic_bundle_status": bundle["status"],
    }
    (out / "figure2_data_failure_verification_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if overall != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
