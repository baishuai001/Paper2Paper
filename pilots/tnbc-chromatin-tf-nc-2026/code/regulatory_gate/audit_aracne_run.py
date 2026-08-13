#!/usr/bin/env python3
"""Create a compact, fail-closed receipt for a completed ARACNe3 run."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from common import output_manifest, sha256_file, write_json


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--regulators", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--expected-subnetworks", type=int, default=100)
    parser.add_argument("--threads", type=int, default=24)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    log = args.run_dir / "log_crc.txt"
    network = args.run_dir / "consolidated-net_crc.tsv"
    author_network = args.run_dir / "subnets" / "subnet1_crc.tsv"
    binary = args.repo / "build-system" / "ARACNe3_app_release"
    subnetworks = sorted((args.run_dir / "subnets").glob("subnet*_crc.tsv"))
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=args.repo, text=True
    ).strip()
    log_text = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
    match = re.search(r"Total subnetworks generated:\s*(\d+)", log_text)
    reported = int(match.group(1)) if match else None
    with network.open(encoding="utf-8", errors="strict") as stream:
        edge_count = sum(1 for _ in stream) - 1 if network.is_file() else -1
    with author_network.open(encoding="utf-8", errors="strict") as stream:
        author_edge_count = (
            sum(1 for _ in stream) - 1 if author_network.is_file() else -1
        )
    conditions = {
        "official_commit": commit == "3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14",
        "log_success": "SUCCESS!" in log_text,
        "subnetwork_files": len(subnetworks) == args.expected_subnetworks,
        "reported_subnetworks": reported == args.expected_subnetworks,
        "network_nonempty": edge_count > 0,
        "author_subnetwork_nonempty": author_edge_count > 0,
    }
    receipt = {
        "status": "passed" if all(conditions.values()) else "failed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "conditions": conditions,
        "commit": commit,
        "parameters": {
            "subnetworks": args.expected_subnetworks,
            "subsample": 0.63212,
            "alpha": 0.05,
            "multiple_testing": "FDR",
            "maximum_entropy_pruning": True,
            "seed": 1729,
            "threads": args.threads,
        },
        "subnetwork_files": len(subnetworks),
        "consolidated_edges": edge_count,
        "author_subnetwork_edges": author_edge_count,
        "regulon_source": "subnets/subnet1_crc.tsv",
        "input": {
            "expression_bytes": args.expression.stat().st_size,
            "expression_sha256": sha256_file(args.expression),
            "regulator_count": len(args.regulators.read_text(encoding="utf-8").splitlines()),
            "regulators_sha256": sha256_file(args.regulators),
        },
        "outputs": output_manifest([binary, log, author_network, network]),
    }
    write_json(args.output, receipt)
    print(json.dumps(receipt, ensure_ascii=False))
    return 0 if receipt["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
