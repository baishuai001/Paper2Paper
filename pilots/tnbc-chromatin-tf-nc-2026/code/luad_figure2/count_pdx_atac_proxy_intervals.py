#!/usr/bin/env python3
"""Count PDX human ATAC alignments in the frozen TCGA-LUAD intervals.

The script is resumable, verifies that the TCGA proxy validation passed, and
never substitutes sample-specific peaks for the frozen 139,135-coordinate
universe.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import io
import json
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_frozen_bed(zip_path: Path, output: Path) -> int:
    if output.is_file() and output.stat().st_size:
        with output.open(encoding="utf-8") as handle:
            rows = sum(1 for _ in handle)
        if rows != 139_135:
            raise RuntimeError(f"Existing frozen BED has {rows} rows, expected 139135")
        return rows
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".partial")
    with zipfile.ZipFile(zip_path) as archive:
        member = "LUAD_raw_counts.txt"
        if member not in archive.namelist():
            raise RuntimeError(f"{member} is absent from {zip_path}")
        with archive.open(member) as raw, io.TextIOWrapper(raw, encoding="utf-8") as text_in:
            reader = csv.DictReader(text_in, delimiter="\t")
            required = {"seqnames", "start", "end", "name"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise RuntimeError("Unexpected TCGA-LUAD raw-count matrix schema")
            rows = 0
            with temporary.open("w", encoding="utf-8", newline="") as text_out:
                writer = csv.writer(text_out, delimiter="\t", lineterminator="\n")
                for row in reader:
                    writer.writerow([row["seqnames"], row["start"], row["end"], row["name"]])
                    rows += 1
    if rows != 139_135:
        raise RuntimeError(f"Extracted {rows} intervals, expected 139135")
    temporary.replace(output)
    return rows


def run_text(command: list[str]) -> str:
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}"
        )
    return result.stdout.strip()


def count_model(model: str, bam: Path, bed: Path, root: Path, bed_hash: str) -> dict:
    output_dir = root / "data/processed/pdx_atac_proxy_counts"
    audit_dir = root / "audit/pdx_atac_proxy_counts"
    log_dir = root / "logs/pdx_atac_proxy_counts"
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"{model}.tcga_luad_intervals.counts.tsv.gz"
    receipt_path = audit_dir / f"{model}.json"
    log_path = log_dir / f"{model}.log"
    if output.is_file() and output.stat().st_size and receipt_path.is_file():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("frozen_interval_bed_sha256") == bed_hash:
            return {"model": model, "status": "RESUMED_COMPLETE", **receipt}

    temporary = output.with_suffix(output.suffix + ".partial")
    command = ["bedtools", "coverage", "-counts", "-a", str(bed), "-b", str(bam)]
    row_count = 0
    count_sum = 0
    with log_path.open("a", encoding="utf-8") as log:
        log.write(f"[{datetime.now(timezone.utc).isoformat()}] {' '.join(command)}\n")
        log.flush()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
        assert process.stdout is not None
        with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
            handle.write("chromosome\tstart\tend\tpeak_id\tcount\n")
            for raw_line in io.TextIOWrapper(process.stdout, encoding="utf-8"):
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) != 5:
                    process.kill()
                    raise RuntimeError(f"Unexpected bedtools output for {model}: {raw_line[:200]}")
                value = int(fields[4])
                handle.write(raw_line)
                row_count += 1
                count_sum += value
        return_code = process.wait()
    if return_code:
        raise RuntimeError(f"bedtools coverage failed for {model}; see {log_path}")
    if row_count != 139_135 or count_sum <= 0:
        raise RuntimeError(
            f"Invalid interval counts for {model}: rows={row_count}, sum={count_sum}"
        )
    temporary.replace(output)
    bam_alignments = int(run_text(["samtools", "view", "-c", str(bam)]))
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "model": model,
        "status": "COMPLETE",
        "bam": str(bam),
        "bam_bytes": bam.stat().st_size,
        "human_filtered_bam_alignments": bam_alignments,
        "frozen_interval_bed": str(bed),
        "frozen_interval_bed_sha256": bed_hash,
        "frozen_interval_count": row_count,
        "count_sum_over_frozen_intervals": count_sum,
        "cpm_library_size_definition": "count_sum_over_frozen_intervals",
        "output": str(output),
        "output_sha256": sha256(output),
        "counting_command": command,
    }
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state_audit_root", type=Path)
    parser.add_argument("figure2_root", type=Path)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    state_root = args.state_audit_root.resolve()
    figure2_root = args.figure2_root.resolve()
    if args.workers < 1:
        raise ValueError("--workers must be positive")

    validation_path = (
        state_root
        / "audit/atac_state_proxy/tcga_atac_state_proxy_validation_receipt.json"
    )
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if not validation["tru_vs_rest"].get("pdx_transfer_authorized", False):
        raise RuntimeError("Frozen TCGA ATAC proxy validation did not authorize PDX transfer")

    raw_zip = (
        figure2_root
        / "data/raw/tcga_atac_gdc/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip"
    )
    bed = state_root / "reference/atac_proxy/tcga_luad_139135_consensus_intervals.bed"
    interval_count = write_frozen_bed(raw_zip, bed)
    bed_hash = sha256(bed)

    pdx_root = figure2_root / "data/processed/raw_atac/PDX"
    models = sorted(path.name for path in pdx_root.iterdir() if path.is_dir())
    if len(models) != 13:
        raise RuntimeError(f"Frozen PDX set must contain 13 model directories, found {len(models)}")
    bams: dict[str, Path] = {}
    for model in models:
        bam = pdx_root / model / f"{model}.filtered.bam"
        if not bam.is_file() or not bam.stat().st_size:
            raise FileNotFoundError(bam)
        bams[model] = bam

    completed: list[dict] = []
    errors: list[dict] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(count_model, model, bams[model], bed, state_root, bed_hash): model
            for model in models
        }
        for future in concurrent.futures.as_completed(futures):
            model = futures[future]
            try:
                result = future.result()
                completed.append(result)
                print(json.dumps({"model": model, "status": result["status"]}), flush=True)
            except Exception as error:  # preserve all successful models and record failures
                errors.append({"model": model, "error": repr(error)})
                print(json.dumps(errors[-1], sort_keys=True), flush=True)

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tcga_proxy_validation_sha256": sha256(validation_path),
        "pdx_transfer_authorized": True,
        "frozen_interval_bed": str(bed),
        "frozen_interval_bed_sha256": bed_hash,
        "frozen_interval_count": interval_count,
        "workers": args.workers,
        "completed_models": sorted(completed, key=lambda row: row["model"]),
        "errors": sorted(errors, key=lambda row: row["model"]),
    }
    aggregate_path = (
        state_root / "audit/pdx_atac_proxy_counts/pdx_atac_proxy_counting_receipt.json"
    )
    aggregate_path.parent.mkdir(parents=True, exist_ok=True)
    aggregate_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

