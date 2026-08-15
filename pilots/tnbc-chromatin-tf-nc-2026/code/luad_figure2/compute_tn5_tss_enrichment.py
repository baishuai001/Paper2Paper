#!/usr/bin/env python3
"""Compute a layout-independent Tn5 insertion-site TSS enrichment QC profile.

ataqv defines HQAA using properly paired reads, so its TSS score is undefined
for the public single-end DRA cell-line libraries.  This report-only metric
uses the same filtered BAMs and the same Tn5 cut-site convention for both
single- and paired-end samples.  It does not exclude an otherwise completely
processed sample from the anchor-style analysis.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path


EXTENSION = 1000
FLANK_WIDTH = 100
CENTRAL_HALF_WIDTH = 50


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_tss(tss_path: Path, blacklist_path: Path) -> tuple[dict[str, list[int]], dict[str, list[str]], int]:
    command = ["bedtools", "intersect", "-v", "-a", str(tss_path), "-b", str(blacklist_path)]
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    records: dict[str, list[tuple[int, str]]] = {}
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        chrom, position, strand = fields[0], int(fields[1]), fields[5]
        if chrom == "chrX" or chrom.startswith("chr") and chrom[3:].isdigit() and 1 <= int(chrom[3:]) <= 22:
            records.setdefault(chrom, []).append((position, strand))
    positions: dict[str, list[int]] = {}
    strands: dict[str, list[str]] = {}
    for chrom, values in records.items():
        ordered = sorted(values)
        positions[chrom] = [item[0] for item in ordered]
        strands[chrom] = [item[1] for item in ordered]
    return positions, strands, sum(len(values) for values in positions.values())


def process_sample(
    root: Path,
    sample: dict[str, str],
    tss_positions: dict[str, list[int]],
    tss_strands: dict[str, list[str]],
    tss_count: int,
) -> dict[str, object]:
    system = sample["system"]
    sample_id = sample["sample_id"]
    slug = sample["sample_slug"]
    bam = Path(sample["filtered_bam"])
    if not bam.is_file() or bam.stat().st_size == 0:
        raise FileNotFoundError(f"Missing filtered BAM for {system}/{sample_id}: {bam}")

    profile = [0] * (2 * EXTENSION + 1)
    insertion_count = 0
    overlap_events = 0
    process = subprocess.Popen(
        ["bedtools", "bamtobed", "-i", str(bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        chrom = fields[0]
        if chrom not in tss_positions:
            continue
        start, end, strand = int(fields[1]), int(fields[2]), fields[5]
        cut = start + 4 if strand == "+" else max(0, end - 5)
        insertion_count += 1
        coordinates = tss_positions[chrom]
        orientations = tss_strands[chrom]
        left = bisect.bisect_left(coordinates, cut - EXTENSION)
        right = bisect.bisect_right(coordinates, cut + EXTENSION)
        for index in range(left, right):
            tss = coordinates[index]
            distance = cut - tss if orientations[index] == "+" else tss - cut
            if -EXTENSION <= distance <= EXTENSION:
                profile[distance + EXTENSION] += 1
                overlap_events += 1
    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"bedtools bamtobed failed for {system}/{sample_id}: {stderr.strip()}")
    if insertion_count == 0:
        raise RuntimeError(f"No Tn5 insertion sites were found for {system}/{sample_id}")

    flank_counts = profile[:FLANK_WIDTH] + profile[-FLANK_WIDTH:]
    flank_mean = sum(flank_counts) / len(flank_counts)
    if flank_mean <= 0:
        raise RuntimeError(f"Zero TSS flank coverage for {system}/{sample_id}")
    normalized = [count / flank_mean for count in profile]
    center_index = EXTENSION
    central = normalized[center_index - CENTRAL_HALF_WIDTH : center_index + CENTRAL_HALF_WIDTH + 1]

    profile_path = root / "audit/tn5_tss" / system / slug / "tn5_tss_profile.tsv"
    profile_rows = [
        {
            "distance_from_tss": distance,
            "tn5_insertion_count": profile[distance + EXTENSION],
            "flank_normalized_enrichment": normalized[distance + EXTENSION],
        }
        for distance in range(-EXTENSION, EXTENSION + 1)
    ]
    write_tsv(profile_path, profile_rows, list(profile_rows[0]))
    return {
        "system": system,
        "sample_id": sample_id,
        "sample_slug": slug,
        "run": sample["run"],
        "library_layout": sample["layout"],
        "filtered_bam": str(bam),
        "filtered_bam_sha256": sha256(bam),
        "protein_coding_tss_after_blacklist": tss_count,
        "tn5_insertions": insertion_count,
        "tn5_tss_overlap_events": overlap_events,
        "flank_mean_insertions_per_position": flank_mean,
        "tn5_tss_enrichment_exact_center": normalized[center_index],
        "tn5_tss_enrichment_max_central_101bp": max(central),
        "tn5_tss_enrichment_mean_central_101bp": sum(central) / len(central),
        "profile_points": len(profile_rows),
        "profile_tsv": str(profile_path),
        "profile_sha256": sha256(profile_path),
        "status": "PASS",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    manifest_path = root / "audit/raw_atac_qc/figure2_analysis_raw_atac_samples.tsv"
    tss_path = root / "reference/figure2_features/gencode.v47.protein_coding.gene_tss.bed"
    blacklist_path = root / "reference/hg38-blacklist.v2.bed"
    manifest = read_tsv(manifest_path)
    if not manifest:
        raise RuntimeError("No analysis samples available for Tn5 TSS enrichment")
    tss_positions, tss_strands, tss_count = load_tss(tss_path, blacklist_path)
    if tss_count == 0:
        raise RuntimeError("No protein-coding TSS remained after blacklist exclusion")

    workers = max(1, int(os.environ.get("TSS_JOBS", "4")))
    rows: list[dict[str, object]] = []
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(process_sample, root, sample, tss_positions, tss_strands, tss_count): sample
            for sample in manifest
        }
        for future in as_completed(futures):
            sample = futures[future]
            try:
                rows.append(future.result())
                print(f"TN5_TSS_COMPLETE {sample['system']} {sample['sample_id']}", flush=True)
            except Exception as exc:
                errors.append(f"{sample['system']}/{sample['sample_id']}:{type(exc).__name__}:{exc}")
                print(f"TN5_TSS_FAIL {errors[-1]}", flush=True)

    rows.sort(key=lambda item: (str(item["system"]), str(item["sample_id"]).lower()))
    out_dir = root / "audit/tn5_tss"
    summary_path = out_dir / "figure2_tn5_tss_enrichment.tsv"
    if rows:
        write_tsv(summary_path, rows, list(rows[0]))
    bedtools_version = subprocess.run(["bedtools", "--version"], check=True, text=True, capture_output=True).stdout.strip()
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Tn5 insertion sites (+4/-5 bp), strand-oriented +/-1000 bp around GENCODE v47 protein-coding TSS; normalized by mean of outermost 100 bp on each side",
        "score_fields": {
            "exact_center": "normalized insertion count at distance 0",
            "max_central_101bp": "maximum normalized insertion count from -50 to +50 bp",
            "mean_central_101bp": "mean normalized insertion count from -50 to +50 bp",
        },
        "ataqv_hqaa_requires_proper_pairs": True,
        "this_layout_independent_metric_is_report_only_not_a_gate": True,
        "ataqv_source": "https://github.com/ParkerLab/ataqv/blob/master/src/cpp/Metrics.cpp",
        "bedtools_version": bedtools_version,
        "analysis_manifest": str(manifest_path),
        "analysis_manifest_sha256": sha256(manifest_path),
        "tss_bed": str(tss_path),
        "tss_bed_sha256": sha256(tss_path),
        "blacklist_bed": str(blacklist_path),
        "blacklist_bed_sha256": sha256(blacklist_path),
        "protein_coding_tss_after_blacklist": tss_count,
        "expected_samples": len(manifest),
        "completed_samples": len(rows),
        "profile_points_total": sum(int(row["profile_points"]) for row in rows),
        "status": "PASS" if not errors and len(rows) == len(manifest) else "FAIL",
        "errors": errors,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = out_dir / "figure2_tn5_tss_enrichment_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
