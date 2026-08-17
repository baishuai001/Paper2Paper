#!/usr/bin/env python3
"""Build fixed-width, peak-count- and GC-matched inputs for the LUAD motif audit.

This is a diagnostic analysis.  It does not replace the anchor-style Figure 2
result.  Every target is resized to 200 bp, filtered against one hg38 blacklist,
and sampled to the largest common GC-bin quota supported by every technically
matchable sample.  A single LUSC accessible-region background is built with the
same 200-bp width and GC-bin distribution and is reused for every system.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
import zipfile
from bisect import bisect_left
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pysam


WIDTH = 200
GC_BIN_WIDTH = 0.05
MAX_BACKGROUND_MULTIPLIER = 10
SYSTEM_ORDER = ("patient", "PDX", "cell_line")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_key(seed: int, namespace: str, peak: tuple[str, int, int, float, int]) -> str:
    chrom, start, end, _, _ = peak
    return hashlib.sha256(f"{seed}|{namespace}|{chrom}:{start}-{end}".encode()).hexdigest()


def canonical_chromosome(chrom: str) -> bool:
    return chrom == "chrX" or (
        chrom.startswith("chr") and chrom[3:].isdigit() and 1 <= int(chrom[3:]) <= 22
    )


class Blacklist:
    def __init__(self, path: Path) -> None:
        intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) >= 3 and canonical_chromosome(fields[0]):
                    intervals[fields[0]].append((int(fields[1]), int(fields[2])))
        self.intervals = {chrom: sorted(values) for chrom, values in intervals.items()}
        self.starts = {chrom: [start for start, _ in values] for chrom, values in self.intervals.items()}

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        values = self.intervals.get(chrom, [])
        starts = self.starts.get(chrom, [])
        if not values:
            return False
        index = bisect_left(starts, end)
        if index > 0 and values[index - 1][1] > start:
            return True
        return index < len(values) and values[index][0] < end and values[index][1] > start


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", errors="replace"
    )


def fixed_interval(fields: list[str], chrom_sizes: dict[str, int]) -> tuple[str, int, int] | None:
    if len(fields) < 3 or not canonical_chromosome(fields[0]):
        return None
    chrom = fields[0]
    if chrom not in chrom_sizes:
        return None
    try:
        original_start, original_end = int(fields[1]), int(fields[2])
    except ValueError:
        return None
    if original_end <= original_start:
        return None
    summit = None
    if len(fields) >= 10:
        try:
            offset = int(float(fields[9]))
            if offset >= 0:
                summit = original_start + offset
        except ValueError:
            summit = None
    center = summit if summit is not None else (original_start + original_end) // 2
    chrom_size = chrom_sizes[chrom]
    if chrom_size < WIDTH:
        return None
    start = max(0, min(center - WIDTH // 2, chrom_size - WIDTH))
    return chrom, start, start + WIDTH


def gc_bin_for_interval(
    fasta: pysam.FastaFile,
    cache: dict[tuple[str, int, int], tuple[float, int] | None],
    interval: tuple[str, int, int],
) -> tuple[float, int] | None:
    if interval in cache:
        return cache[interval]
    sequence = fasta.fetch(*interval).upper()
    if len(sequence) != WIDTH or any(base not in "ACGT" for base in sequence):
        cache[interval] = None
        return None
    gc = (sequence.count("G") + sequence.count("C")) / WIDTH
    gc_bin = min(int(gc / GC_BIN_WIDTH), int(1 / GC_BIN_WIDTH) - 1)
    cache[interval] = (gc, gc_bin)
    return cache[interval]


def normalize_peak_lines(
    lines,
    chrom_sizes: dict[str, int],
    blacklist: Blacklist,
    fasta: pysam.FastaFile,
    gc_cache: dict[tuple[str, int, int], tuple[float, int] | None],
) -> tuple[list[tuple[str, int, int, float, int]], dict[str, int]]:
    peaks: dict[tuple[str, int, int], tuple[str, int, int, float, int]] = {}
    stats = defaultdict(int)
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        stats["input_rows"] += 1
        fields = line.rstrip().split("\t")
        interval = fixed_interval(fields, chrom_sizes)
        if interval is None:
            stats["noncanonical_or_invalid"] += 1
            continue
        chrom, start, end = interval
        if blacklist.overlaps(chrom, start, end):
            stats["blacklist_removed"] += 1
            continue
        gc_result = gc_bin_for_interval(fasta, gc_cache, interval)
        if gc_result is None:
            stats["non_acgt_removed"] += 1
            continue
        gc, gc_bin = gc_result
        peaks.setdefault(interval, (chrom, start, end, gc, gc_bin))
    stats["deduplicated_usable"] = len(peaks)
    return list(peaks.values()), dict(stats)


def read_peak_file(
    path: Path,
    chrom_sizes: dict[str, int],
    blacklist: Blacklist,
    fasta: pysam.FastaFile,
    gc_cache: dict[tuple[str, int, int], tuple[float, int] | None],
) -> tuple[list[tuple[str, int, int, float, int]], dict[str, int]]:
    with open_text(path) as handle:
        return normalize_peak_lines(handle, chrom_sizes, blacklist, fasta, gc_cache)


def read_lusc_background(
    archive: Path,
    chrom_sizes: dict[str, int],
    blacklist: Blacklist,
    fasta: pysam.FastaFile,
    gc_cache: dict[tuple[str, int, int], tuple[float, int] | None],
) -> tuple[list[tuple[str, int, int, float, int]], dict[str, int]]:
    with zipfile.ZipFile(archive) as zipped:
        member = next(name for name in zipped.namelist() if name.endswith("LUSC_raw_counts.txt"))
        with zipped.open(member) as raw:
            lines = (line.decode("utf-8", errors="replace") for line in raw)
            next(lines, None)
            return normalize_peak_lines(lines, chrom_sizes, blacklist, fasta, gc_cache)


def write_bed(path: Path, peaks: list[tuple[str, int, int, float, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for chrom, start, end, _, _ in sorted(peaks, key=lambda value: (value[0], value[1], value[2])):
            handle.write(f"{chrom}\t{start}\t{end}\n")


def summarize_gc(peaks: list[tuple[str, int, int, float, int]]) -> tuple[float, float]:
    values = [peak[3] for peak in peaks]
    if not values:
        return math.nan, math.nan
    return statistics.mean(values), statistics.pstdev(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--seed", type=int, default=20260817)
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

    cell_counts = [len(sample_peaks[(row["system"], row["sample_id"])]) for row in manifest if row["system"] == "cell_line"]
    if not cell_counts or min(cell_counts) <= 0:
        raise RuntimeError("Cannot define matching floor from the complete cell-line arm")
    match_floor = min(cell_counts)
    eligible_rows = [
        row
        for row in manifest
        if len(sample_peaks[(row["system"], row["sample_id"])]) >= match_floor
    ]
    if not eligible_rows:
        raise RuntimeError("No motif-matchable samples")

    n_bins = int(1 / GC_BIN_WIDTH)
    per_sample_bins: dict[tuple[str, str], dict[int, list[tuple[str, int, int, float, int]]]] = {}
    for row in eligible_rows:
        key = (row["system"], row["sample_id"])
        bins: dict[int, list[tuple[str, int, int, float, int]]] = defaultdict(list)
        for peak in sample_peaks[key]:
            bins[peak[4]].append(peak)
        per_sample_bins[key] = bins
    quotas = {
        gc_bin: min(len(per_sample_bins[(row["system"], row["sample_id"])].get(gc_bin, [])) for row in eligible_rows)
        for gc_bin in range(n_bins)
    }
    matched_peak_count = sum(quotas.values())
    if matched_peak_count <= 0:
        raise RuntimeError("The common GC-stratified peak quota is empty")

    out_root = root / "data/processed/motif_equivalence/harmonized_200bp"
    target_dir = out_root / "targets"
    selected_by_sample: dict[tuple[str, str], list[tuple[str, int, int, float, int]]] = {}
    for row in eligible_rows:
        key = (row["system"], row["sample_id"])
        selected: list[tuple[str, int, int, float, int]] = []
        for gc_bin, quota in quotas.items():
            candidates = sorted(
                per_sample_bins[key].get(gc_bin, []),
                key=lambda peak: stable_key(args.seed, f"{key[0]}|{key[1]}|bin{gc_bin}", peak),
            )
            selected.extend(candidates[:quota])
        if len(selected) != matched_peak_count:
            raise RuntimeError(f"GC quota mismatch for {key}: {len(selected)} != {matched_peak_count}")
        selected_by_sample[key] = selected
        write_bed(target_dir / f"{row['system']}.{row['sample_slug']}.bed", selected)

    lusc_peaks, lusc_stats = read_lusc_background(tcga_archive, chrom_sizes, blacklist, fasta, gc_cache)
    selected_coordinates = {
        (peak[0], peak[1], peak[2]) for peaks in selected_by_sample.values() for peak in peaks
    }
    background_bins: dict[int, list[tuple[str, int, int, float, int]]] = defaultdict(list)
    for peak in lusc_peaks:
        if (peak[0], peak[1], peak[2]) not in selected_coordinates:
            background_bins[peak[4]].append(peak)
    possible_multipliers = [
        len(background_bins[gc_bin]) // quota for gc_bin, quota in quotas.items() if quota > 0
    ]
    background_multiplier = min(MAX_BACKGROUND_MULTIPLIER, min(possible_multipliers))
    if background_multiplier <= 0:
        raise RuntimeError("LUSC accessible background cannot satisfy the shared GC-bin quota")
    background: list[tuple[str, int, int, float, int]] = []
    for gc_bin, quota in quotas.items():
        candidates = sorted(
            background_bins[gc_bin],
            key=lambda peak: stable_key(args.seed, f"LUSC_background|bin{gc_bin}", peak),
        )
        background.extend(candidates[: quota * background_multiplier])
    background_path = out_root / "LUSC_accessible_common_GC_matched_background.bed"
    write_bed(background_path, background)

    output_rows: list[dict[str, object]] = []
    for row in manifest:
        key = (row["system"], row["sample_id"])
        eligible = key in selected_by_sample
        target_path = target_dir / f"{row['system']}.{row['sample_slug']}.bed"
        if not eligible:
            write_bed(target_path, [])
        selected = selected_by_sample.get(key, [])
        original_gc_mean, original_gc_sd = summarize_gc(sample_peaks[key])
        selected_gc_mean, selected_gc_sd = summarize_gc(selected)
        stats = sample_stats[key]
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
                "harmonized_200bp_unique_peaks": len(sample_peaks[key]),
                "cell_line_arm_matching_floor": match_floor,
                "motif_matchable": str(eligible).upper(),
                # Keep this field non-empty because Bash treats tab as IFS
                # whitespace and otherwise collapses the empty column while
                # reading the manifest.
                "nonmatchable_reason": "MATCHABLE" if eligible else f"FEWER_THAN_COMPLETE_CELL_LINE_ARM_MINIMUM_{match_floor}_PEAKS",
                "matched_peak_count": len(selected),
                "original_gc_mean": original_gc_mean,
                "original_gc_sd": original_gc_sd,
                "matched_gc_mean": selected_gc_mean,
                "matched_gc_sd": selected_gc_sd,
                "matched_target_path": str(target_path),
                "matched_target_sha256": sha256(target_path),
            }
        )
    manifest_out = out_root / "harmonized_motif_manifest.tsv"
    write_tsv(
        manifest_out,
        output_rows,
        [
            "system", "sample_id", "sample_slug", "source_peak_path", "source_input_rows",
            "noncanonical_or_invalid_removed", "blacklist_removed", "non_acgt_removed",
            "harmonized_200bp_unique_peaks", "cell_line_arm_matching_floor", "motif_matchable",
            "nonmatchable_reason", "matched_peak_count", "original_gc_mean", "original_gc_sd",
            "matched_gc_mean", "matched_gc_sd", "matched_target_path", "matched_target_sha256",
        ],
    )
    quota_rows = [
        {
            "gc_bin": gc_bin,
            "gc_lower_inclusive": gc_bin * GC_BIN_WIDTH,
            "gc_upper_exclusive": (gc_bin + 1) * GC_BIN_WIDTH,
            "target_peaks_per_sample": quota,
            "background_peaks": quota * background_multiplier,
        }
        for gc_bin, quota in quotas.items()
    ]
    quota_path = out_root / "shared_gc_bin_quotas.tsv"
    write_tsv(
        quota_path,
        quota_rows,
        ["gc_bin", "gc_lower_inclusive", "gc_upper_exclusive", "target_peaks_per_sample", "background_peaks"],
    )
    background_gc_mean, background_gc_sd = summarize_gc(background)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "diagnostic_method_equivalence_audit_not_primary_biological_gate",
        "seed": args.seed,
        "fixed_peak_width_bp": WIDTH,
        "gc_bin_width": GC_BIN_WIDTH,
        "blacklist": str(blacklist_path),
        "blacklist_sha256": sha256(blacklist_path),
        "fasta_fai_sha256": sha256(fai_path),
        "input_manifest_sha256": sha256(manifest_path),
        "tcga_lusc_background_source_sha256": sha256(tcga_archive),
        "motif_matchability_rule": "at least the minimum harmonized peak count observed across the complete 19-cell-line arm",
        "cell_line_arm_matching_floor": match_floor,
        "eligible_samples_by_system": {
            system: sum(row["system"] == system and row["motif_matchable"] == "TRUE" for row in output_rows)
            for system in SYSTEM_ORDER
        },
        "all_samples_by_system": {
            system: sum(row["system"] == system for row in output_rows) for system in SYSTEM_ORDER
        },
        "shared_matched_target_peak_count": matched_peak_count,
        "background_source": "TCGA-LUSC accessible regions from GDC cancer-type raw-count matrix",
        "background_multiplier": background_multiplier,
        "shared_background_peak_count": len(background),
        "background_gc_mean": background_gc_mean,
        "background_gc_sd": background_gc_sd,
        "background_sha256": sha256(background_path),
        "lusc_background_preselection_stats": lusc_stats,
        "target_and_background_lengths_are_exactly_200bp": True,
        "same_background_file_reused_for_all_three_systems": True,
        "exact_target_gc_bin_quotas_reused_for_every_matchable_sample": True,
    }
    receipt_path = root / "audit/motif_equivalence/harmonized_input_receipt.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
