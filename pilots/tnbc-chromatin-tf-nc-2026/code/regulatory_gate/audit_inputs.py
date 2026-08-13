#!/usr/bin/env python3
"""Fail-closed audit of immutable inputs and primary software identities."""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py

from common import output_manifest, sha256_file, write_json


REQUIRED_PIPELINE_FILES = [
    "audit_inputs.py",
    "audit_aracne_run.py",
    "build_pseudobulk.py",
    "common.py",
    "decide_gate.py",
    "extract_pango_regulators.py",
    "phenotype.py",
    "prepare_aracne_inputs.py",
    "prepare_tcga_crc.R",
    "run_aracne3.sh",
    "run_cloud.sh",
    "statistics.py",
    "viper_msviper.R",
]


def command_output(command: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(command, cwd=cwd, text=True, stderr=subprocess.STDOUT).strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--expected-h5ad-bytes", required=True, type=int)
    parser.add_argument("--expected-h5ad-sha256", required=True)
    parser.add_argument("--aracne-repo", required=True, type=Path)
    parser.add_argument("--expected-aracne-commit", required=True)
    parser.add_argument("--anchor-code-repo", required=True, type=Path)
    parser.add_argument("--expected-anchor-code-commit", required=True)
    parser.add_argument("--supplement", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--pipeline-code-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--skip-h5ad-hash", action="store_true")
    args = parser.parse_args()

    if not args.h5ad.is_file() or not args.supplement.is_file() or not args.manifest.is_file():
        raise FileNotFoundError("Required H5AD, supplement or phenotype manifest is absent")
    pipeline_paths = [args.pipeline_code_dir / name for name in REQUIRED_PIPELINE_FILES]
    missing_pipeline_files = [path.name for path in pipeline_paths if not path.is_file()]
    if missing_pipeline_files:
        raise FileNotFoundError(f"Missing pipeline files: {missing_pipeline_files}")
    actual_bytes = args.h5ad.stat().st_size
    actual_hash = None if args.skip_h5ad_hash else sha256_file(args.h5ad)
    commit = command_output(["git", "rev-parse", "HEAD"], cwd=args.aracne_repo)
    anchor_code_commit = command_output(
        ["git", "rev-parse", "HEAD"], cwd=args.anchor_code_repo
    )
    with h5py.File(args.h5ad, "r") as handle:
        shape = [int(handle["X"].attrs.get("shape", [len(handle["obs/_index"]), len(handle["var/_index"])])[0]),
                 int(handle["X"].attrs.get("shape", [len(handle["obs/_index"]), len(handle["var/_index"])])[1])]
        obs_fields = sorted(handle["obs"].keys())

    conditions = {
        "h5ad_size_matches": actual_bytes == args.expected_h5ad_bytes,
        "h5ad_sha256_matches": args.skip_h5ad_hash or actual_hash == args.expected_h5ad_sha256.upper(),
        "aracne_commit_matches": commit == args.expected_aracne_commit,
        "anchor_code_commit_matches": (
            anchor_code_commit == args.expected_anchor_code_commit
        ),
        "supplement_nonempty": args.supplement.stat().st_size > 0,
        "phenotype_manifest_nonempty": args.manifest.stat().st_size > 0,
        "pipeline_code_complete": not missing_pipeline_files,
    }
    receipt = {
        "status": "passed" if all(conditions.values()) else "failed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "conditions": conditions,
        "h5ad": {"path": str(args.h5ad), "bytes": actual_bytes, "sha256": actual_hash, "shape": shape},
        "h5ad_obs_field_count": len(obs_fields),
        "aracne3": {"path": str(args.aracne_repo), "commit": commit},
        "anchor_code": {
            "path": str(args.anchor_code_repo),
            "commit": anchor_code_commit,
            "tag": "v1.0",
            "source": "https://git.codeocean.com/capsule-7227095.git",
        },
        "supplement": {
            "path": str(args.supplement),
            "bytes": args.supplement.stat().st_size,
            "sha256": sha256_file(args.supplement),
        },
        "phenotype_manifest": {
            "path": str(args.manifest),
            "bytes": args.manifest.stat().st_size,
            "sha256": sha256_file(args.manifest),
        },
        "pipeline_code": {
            "path": str(args.pipeline_code_dir),
            "files": output_manifest(pipeline_paths),
        },
        "environment": {
            "python_executable": sys.executable,
            "python": platform.python_version(),
            "platform": platform.platform(),
            "git": command_output(["git", "--version"]),
        },
    }
    write_json(args.output, receipt)
    print(json.dumps(receipt, ensure_ascii=False))
    return 0 if receipt["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
