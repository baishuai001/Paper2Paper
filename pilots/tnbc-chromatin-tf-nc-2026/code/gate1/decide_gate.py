#!/usr/bin/env python3
"""Fail-closed decision engine for the frozen first-level gate."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from gate1_common import sha256_file, write_json


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def run(audit_path: Path, phenotype_path: Path, output_dir: Path, git_commit: str) -> dict:
    audit = read(audit_path)
    phenotype = read(phenotype_path)
    if audit.get("status") != "passed":
        verdict = "UNINTERPRETABLE"
        reason = "The frozen minimum data/identity conditions were not met."
    elif phenotype.get("status") != "passed":
        verdict = "FAIL"
        reason = "Data were sufficient, but at least one pre-specified phenotype-reconstruction criterion failed."
    else:
        verdict = "PENDING_DOWNSTREAM"
        reason = "The phenotype gate passed; TF and validation stages are required before a final PASS/FAIL decision."
    stop = verdict in {"UNINTERPRETABLE", "FAIL"}
    payload = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit,
        "verdict": verdict,
        "stop_now": stop,
        "reason": reason,
        "data_audit_status": audit.get("status"),
        "phenotype_reconstruction_status": phenotype.get("status"),
        "failed_phenotype_conditions": [key for key, value in phenotype.get("conditions", {}).items() if not value],
        "downstream_status": "NOT_RUN_DUE_TO_PRE_SPECIFIED_STOP" if stop else "REQUIRED",
        "prohibited_after_stop": ["Cancer-cell pseudobulk", "TF activity", "TCGA context network", "drug", "ATAC", "spatial", "wet-lab inference"],
        "inputs": {
            audit_path.name: {"path": str(audit_path.resolve()), "sha256": sha256_file(audit_path)},
            phenotype_path.name: {"path": str(phenotype_path.resolve()), "sha256": sha256_file(phenotype_path)},
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "G0_gate1_decision.json", payload)
    print(payload)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--phenotype", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--git-commit", required=True)
    args = parser.parse_args()
    result = run(args.audit, args.phenotype, args.output_dir, args.git_commit)
    return 2 if result["stop_now"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
