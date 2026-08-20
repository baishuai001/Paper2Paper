#!/usr/bin/env python3
"""Build the primary TNBC-analog LUAD motif inputs.

Unlike the peak-count-matched diagnostic analysis, this primary input keeps
every usable peak from every pre-defined sample.  Query peaks and the shared
lung accessible consensus background are all represented as canonical,
blacklist-free 200-bp intervals.  The same background is reused for patients,
PDX, and cell lines; HOMER performs its own oligonucleotide normalization.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pysam

from harmonize_motif_inputs import (
    WIDTH,
    Blacklist,
    canonical_chromosome,
    gc_bin_for_interval,
    normalize_peak_lines,
    read_peak_file,
    read_tsv,
    sha256,
    write_bed,
    write_tsv,
)


SYSTEM_ORDER = ("patient", "PDX", "cell_line")
TCGA_TYPES = ("LUAD", "LUSC")


def summarize_gc(peaks: list[tuple[str, int, int, float, int]]) -> tuple[float, float]:
    values = [peak[3] for peak in peaks]
    if not values:
        return math.nan, math.nan
    return statistics.mean(values), statistics.pstdev(values)


def read_tcga_type_peaks(
    archive: Path,
    cancer_type: str,
    chrom_sizes: dict[str, int],
    blacklist: Blacklist,
    fasta: pysam.FastaFile,
    gc_cache: dict[tuple[str, int, int], tuple[float, int] | None],
) -> tuple[list[tuple[str, int, int, float, int]], dict[str, int]]:
    with zipfile.ZipFile(archive) as zipped:
        member = next(name for name in zipped.namelist() if name.endswith(f"{cancer_type}_raw_counts.txt"))
        with zipped.open(member) as raw:
            lines = (line.decode("utf-8", errors="replace") for line in raw)
            next(lines, None)
            return normalize_peak_lines(lines, chrom_sizes, blacklist, fasta, gc_cache)


def merge_intervals_once(
    coordinates: set[tuple[str, int, int]], chrom_sizes: dict[str, int]
) -> set[tuple[str, int, int]]:
    """Match bedtools merge semantics, then recenter each consensus CRE to 200 bp."""
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in coordinates:
        by_chrom[chrom].append((start, end))
    merged_fixed: set[tuple[str, int, int]] = set()
    for chrom, values in by_chrom.items():
        values.sort()
        current_start, current_end = values[0]
        merged: list[tuple[int, int]] = []
        for start, end in values[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append((current_start, current_end))
                current_start, current_end = start, end
        merged.append((current_start, current_end))
        chrom_size = chrom_sizes[chrom]
        for start, end in merged:
            center = (start + end) // 2
            fixed_start = max(0, min(center - WIDTH // 2, chrom_size - WIDTH))
            merged_fixed.add((chrom, fixed_start, fixed_start + WIDTH))
    return merged_fixed


def build_consensus(
    peaks: list[tuple[str, int, int, float, int]],
    chrom_sizes: dict[str, int],
    blacklist: Blacklist,
    fasta: pysam.FastaFile,
    gc_cache: dict[tuple[str, int, int], tuple[float, int] | None],
) -> tuple[list[tuple[str, int, int, float, int]], int]:
    coordinates = {(chrom, start, end) for chrom, start, end, _, _ in peaks}
    iterations = 0
    while True:
        iterations += 1
        merged = merge_intervals_once(coordinates, chrom_sizes)
        if merged == coordinates:
            break
        coordinates = merged
        if iterations > 20:
            raise RuntimeError("Consensus interval merging did not converge")
    consensus: list[tuple[str, int, int, float, int]] = []
    for interval in sorted(coordinates, key=lambda value: (value[0], value[1], value[2])):
        chrom, start, end = interval
        if blacklist.overlaps(chrom, start, end):
            continue
        gc_result = gc_bin_for_interval(fasta, gc_cache, interval)
        if gc_result is None:
            continue
        gc, gc_bin = gc_result
        consensus.append((chrom, start, end, gc, gc_bin))
    return consensus, iterations


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    manifest_path = root / "audit/manifests/figure2_atomic/figure2_atomic_sample_manifest.tsv"
    blacklist_path = root / "reference/hg38-blacklist.v2.bed"
    fasta_path = root / "reference/hg38.fa"
    fai_path = root / "reference/hg38.fa.fai"
    tcga_archive = root / "data/raw/tcga_atac_gdc/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip"
    for required in (manifest_path, blacklist_path, fasta_path, fai_path, tcga_archive):
        if not required.is_file():
            raise FileNotFoundError(required)

    chrom_sizes: dict[str, int] = {}
    with fai_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields and canonical_chromosome(fields[0]):
                chrom_sizes[fields[0]] = int(fields[1])
    blacklist = Blacklist(blacklist_path)
    fasta = pysam.FastaFile(str(fasta_path))
    gc_cache: dict[tuple[str, int, int], tuple[float, int] | None] = {}
    manifest = read_tsv(manifest_path)

    sample_peaks: dict[tuple[str, str], list[tuple[str, int, int, float, int]]] = {}
    sample_stats: dict[tuple[str, str], dict[str, int]] = {}
    for row in manifest:
        key = (row["system"], row["sample_id"])
        peaks, stats = read_peak_file(Path(row["peak_path"]), chrom_sizes, blacklist, fasta, gc_cache)
        sample_peaks[key] = peaks
        sample_stats[key] = stats

    tcga_peaks: dict[str, list[tuple[str, int, int, float, int]]] = {}
    tcga_stats: dict[str, dict[str, int]] = {}
    for cancer_type in TCGA_TYPES:
        peaks, stats = read_tcga_type_peaks(
            tcga_archive, cancer_type, chrom_sizes, blacklist, fasta, gc_cache
        )
        tcga_peaks[cancer_type] = peaks
        tcga_stats[cancer_type] = stats

    pooled = [peak for values in tcga_peaks.values() for peak in values]
    pooled.extend(peak for values in sample_peaks.values() for peak in values)
    consensus, merge_iterations = build_consensus(
        pooled, chrom_sizes, blacklist, fasta, gc_cache
    )
    if not consensus:
        raise RuntimeError("Shared lung consensus background is empty")

    out_root = root / "data/processed/motif_primary/anchor_consensus_200bp"
    target_dir = out_root / "targets"
    background_path = out_root / "lung_accessible_consensus_background.bed"
    write_bed(background_path, consensus)

    output_rows: list[dict[str, object]] = []
    for row in manifest:
        key = (row["system"], row["sample_id"])
        peaks = sample_peaks[key]
        stats = sample_stats[key]
        target_path = target_dir / f"{row['system']}.{row['sample_slug']}.bed"
        write_bed(target_path, peaks)
        gc_mean, gc_sd = summarize_gc(peaks)
        testable = bool(peaks)
        output_rows.append(
            {
                "system": row["system"],
                "sample_id": row["sample_id"],
                "sample_slug": row["sample_slug"],
                "source_peak_path": row["peak_path"],
                "source_input_rows": stats.get("input_rows", 0),
                "noncanonical_or_invalid_removed": stats.get("noncanonical_or_invalid", 0),
                "blacklist_removed": stats.get("blacklist_removed", 0),
                "non_acgt_removed": stats.get("non_acgt_removed", 0),
                "harmonized_200bp_unique_peaks": len(peaks),
                "cell_line_arm_matching_floor": "NOT_APPLICABLE_PRIMARY_ALL_PEAKS",
                "motif_matchable": str(testable).upper(),
                "nonmatchable_reason": (
                    "ALL_USABLE_200BP_PEAKS_INCLUDED"
                    if testable
                    else "NO_CANONICAL_200BP_PEAKS_AT_MACS2_Q0.01"
                ),
                "matched_peak_count": len(peaks),
                "original_gc_mean": gc_mean,
                "original_gc_sd": gc_sd,
                "matched_gc_mean": gc_mean,
                "matched_gc_sd": gc_sd,
                "matched_target_path": str(target_path),
                "matched_target_sha256": sha256(target_path),
            }
        )

    manifest_out = out_root / "anchor_consensus_motif_manifest.tsv"
    columns = [
        "system", "sample_id", "sample_slug", "source_peak_path", "source_input_rows",
        "noncanonical_or_invalid_removed", "blacklist_removed", "non_acgt_removed",
        "harmonized_200bp_unique_peaks", "cell_line_arm_matching_floor", "motif_matchable",
        "nonmatchable_reason", "matched_peak_count", "original_gc_mean", "original_gc_sd",
        "matched_gc_mean", "matched_gc_sd", "matched_target_path", "matched_target_sha256",
    ]
    write_tsv(manifest_out, output_rows, columns)
    background_gc_mean, background_gc_sd = summarize_gc(consensus)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "primary_TNBC_anchor_analog_all_peaks_shared_consensus_background",
        "author_code_basis": {
            "script": "Figure2/05B-HOMER.sh",
            "query_rule": "one sample peak BED; HOMER -size 200",
            "background_rule": "one shared consensus accessible-peak BED reused across samples",
            "manual_peak_count_matching_in_author_script": False,
        },
        "fixed_peak_width_bp": WIDTH,
        "manual_target_peak_count_matching": False,
        "manual_target_gc_matching": False,
        "gc_control": "same 200-bp background for all samples plus HOMER oligonucleotide normalization",
        "blacklist": str(blacklist_path),
        "blacklist_sha256": sha256(blacklist_path),
        "fasta_fai_sha256": sha256(fai_path),
        "input_manifest_sha256": sha256(manifest_path),
        "all_samples_by_system": {
            system: sum(row["system"] == system for row in output_rows) for system in SYSTEM_ORDER
        },
        "motif_testable_samples_by_system": {
            system: sum(row["system"] == system and row["motif_matchable"] == "TRUE" for row in output_rows)
            for system in SYSTEM_ORDER
        },
        "explicitly_not_testable_samples": [
            {"system": row["system"], "sample_id": row["sample_id"], "reason": row["nonmatchable_reason"]}
            for row in output_rows
            if row["motif_matchable"] != "TRUE"
        ],
        "target_peak_count_by_system": {
            system: {
                "min": min(int(row["matched_peak_count"]) for row in output_rows if row["system"] == system),
                "max": max(int(row["matched_peak_count"]) for row in output_rows if row["system"] == system),
                "median": statistics.median(int(row["matched_peak_count"]) for row in output_rows if row["system"] == system),
            }
            for system in SYSTEM_ORDER
        },
        "background_source": "union of TCGA-LUAD and TCGA-LUSC cancer-type accessible peaks plus all predefined LUAD patient, PDX, and cell-line target peaks; overlapping intervals merged to consensus CREs and recentered to 200 bp",
        "tcga_source_stats": tcga_stats,
        "pooled_background_input_peak_records": len(pooled),
        "consensus_merge_iterations": merge_iterations,
        "shared_background_peak_count": len(consensus),
        "background_gc_mean": background_gc_mean,
        "background_gc_sd": background_gc_sd,
        "background_sha256": sha256(background_path),
        "target_and_background_lengths_are_exactly_200bp": True,
        "same_background_file_reused_for_all_three_systems": True,
        "anchor_analogy_limitation": "the TNBC capsule names but does not provide construction code for its consensus background; this reproducible LUAD consensus uses the same shared-accessible-universe role rather than claiming byte-level reconstruction",
    }
    receipt_path = root / "audit/motif_primary/anchor_consensus_input_receipt.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
