#!/usr/bin/env python3
"""Download and Salmon-quantify the frozen DRA001846 baseline LUAD panel.

The driver is resumable and keeps the SRA archive while deleting only the
explicit per-run temporary FASTQ directory after a successful quantification.
It is intended to run on the cloud server under nohup.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_command(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] {' '.join(command)}\n")
        log.flush()
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    if completed.returncode:
        raise RuntimeError(f"Command failed ({completed.returncode}); see {log_path}")


def locate_sra(root: Path, run: str) -> Path:
    candidates = [root / run / f"{run}.sra", root / run / run, root / f"{run}.sra"]
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size:
            return candidate
    found = list(root.rglob(f"{run}*"))
    found = [path for path in found if path.is_file() and path.stat().st_size]
    if len(found) == 1:
        return found[0]
    raise FileNotFoundError(f"Cannot resolve SRA file for {run}; candidates={found}")


def quantify(row: dict[str, str], root: Path, index: Path, threads: int, manifest_hash: str) -> dict:
    run = row["run"]
    cell_line = row["cell_line"]
    slug = cell_line.replace("/", "_")
    sample_root = root / "data/processed/dra001846_quant" / slug
    quant_dir = sample_root / "salmon"
    quant_file = quant_dir / "quant.sf"
    receipt_path = root / "audit/dra001846_quant" / f"{slug}.json"
    log_path = root / "logs/dra001846_quant" / f"{slug}.log"
    if quant_file.is_file() and quant_file.stat().st_size and receipt_path.is_file():
        return {"cell_line": cell_line, "run": run, "status": "RESUMED_COMPLETE"}

    sra_root = root / "data/raw/dra001846_sra"
    sra_root.mkdir(parents=True, exist_ok=True)
    run_command(["prefetch", "--max-size", "u", "-O", str(sra_root), run], log_path)
    sra_path = locate_sra(sra_root, run)

    fastq_dir = sample_root / "fastq_tmp"
    temp_dir = sample_root / "fasterq_tmp"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            "fasterq-dump", "--split-files", "--threads", str(threads),
            "--temp", str(temp_dir), "--outdir", str(fastq_dir), str(sra_path),
        ],
        log_path,
    )
    read1 = fastq_dir / f"{run}_1.fastq"
    read2 = fastq_dir / f"{run}_2.fastq"
    if not read1.is_file() or not read2.is_file():
        raise FileNotFoundError(f"Expected paired FASTQ files missing for {run}")
    run_command(
        [
            "salmon", "quant", "-i", str(index), "-l", "A", "-1", str(read1),
            "-2", str(read2), "--validateMappings", "--seqBias", "--gcBias",
            "-p", str(threads), "-o", str(quant_dir),
        ],
        log_path,
    )
    if not quant_file.is_file() or not quant_file.stat().st_size:
        raise RuntimeError(f"Salmon did not produce {quant_file}")

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cell_line": cell_line,
        "run": run,
        "submission": row["submission"],
        "experiment": row["experiment"],
        "biosample": row["biosample"],
        "input_manifest_sha256": manifest_hash,
        "sra_path": str(sra_path),
        "sra_bytes": sra_path.stat().st_size,
        "salmon_index": str(index),
        "salmon_quant": str(quant_file),
        "salmon_quant_sha256": sha256(quant_file),
        "temporary_fastq_deleted_after_success": True,
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # Only the explicit temporary directories under this sample are removed.
    shutil.rmtree(fastq_dir)
    shutil.rmtree(temp_dir)
    return {"cell_line": cell_line, "run": run, "status": "COMPLETE"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--salmon-index", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--threads-per-sample", type=int, default=6)
    args = parser.parse_args()
    root = args.run_root.resolve()
    manifest = args.manifest.resolve()
    index = args.salmon_index.resolve()
    if args.workers < 1 or args.threads_per_sample < 1:
        raise ValueError("workers and threads-per-sample must be positive")
    if not manifest.is_file():
        raise FileNotFoundError(manifest)
    # A detached job can be launched before the reference job finishes.
    deadline = time.time() + 6 * 60 * 60
    while not (index / "versionInfo.json").is_file():
        if time.time() > deadline:
            raise TimeoutError(f"Salmon index was not ready after 6 hours: {index}")
        time.sleep(30)

    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 19 or len({row["cell_line"] for row in rows}) != 19:
        raise RuntimeError("Frozen DRA001846 manifest must contain 19 unique cell lines")
    manifest_hash = sha256(manifest)
    statuses: list[dict] = []
    errors: list[dict] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                quantify, row, root, index, args.threads_per_sample, manifest_hash
            ): row
            for row in rows
        }
        for future in as_completed(futures):
            row = futures[future]
            try:
                status = future.result()
                statuses.append(status)
                print(json.dumps(status, sort_keys=True), flush=True)
            except Exception as exc:  # each model is recorded; other models continue
                error = {"cell_line": row["cell_line"], "run": row["run"], "error": repr(exc)}
                errors.append(error)
                print(json.dumps(error, sort_keys=True), flush=True)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest),
        "manifest_sha256": manifest_hash,
        "salmon_index": str(index),
        "workers": args.workers,
        "threads_per_sample": args.threads_per_sample,
        "completed": sorted(statuses, key=lambda row: row["cell_line"]),
        "errors": sorted(errors, key=lambda row: row["cell_line"]),
    }
    summary_path = root / "audit/dra001846_quant/dra001846_quantification_receipt.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
