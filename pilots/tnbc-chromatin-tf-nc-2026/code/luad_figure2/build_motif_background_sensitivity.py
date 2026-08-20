#!/usr/bin/env python3
"""Build the three prespecified shared-background sensitivity variants.

B0 remains the immutable completed primary analysis. B1-B3 reuse byte-exact
B0 target BEDs and vary only the shared 200-bp accessible background.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import pysam

from build_anchor_consensus_motif_inputs import (
    build_consensus,
    read_tcga_type_peaks,
    summarize_gc,
)
from harmonize_motif_inputs import (
    Blacklist,
    canonical_chromosome,
    read_peak_file,
    read_tsv,
    sha256,
    write_bed,
)


VARIANTS = {
    "B1": "TCGA_LUAD_plus_all_three_LUAD_target_systems",
    "B2": "TCGA_LUAD_only",
    "B3": "all_three_LUAD_target_systems_only",
}


def bed_count_and_width(path: Path) -> tuple[int, bool]:
    count = 0
    width_ok = True
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            count += 1
            width_ok &= int(fields[2]) - int(fields[1]) == 200
    return count, width_ok


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    b0_root = root / "data/processed/motif_primary/anchor_consensus_200bp"
    manifest_path = b0_root / "anchor_consensus_motif_manifest.tsv"
    b0_background = b0_root / "lung_accessible_consensus_background.bed"
    blacklist_path = root / "reference/hg38-blacklist.v2.bed"
    fasta_path = root / "reference/hg38.fa"
    fai_path = root / "reference/hg38.fa.fai"
    tcga_archive = root / "data/raw/tcga_atac_gdc/TCGA-ATAC_Cancer_Type-specific_Count_Matrices_raw_counts.zip"
    for path in (manifest_path, b0_background, blacklist_path, fasta_path, fai_path, tcga_archive):
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest = read_tsv(manifest_path)
    if len(manifest) != 54:
        raise RuntimeError(f"Immutable B0 target manifest must have 54 samples, found {len(manifest)}")
    chrom_sizes: dict[str, int] = {}
    with fai_path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if fields and canonical_chromosome(fields[0]):
                chrom_sizes[fields[0]] = int(fields[1])
    blacklist = Blacklist(blacklist_path)
    fasta = pysam.FastaFile(str(fasta_path))
    gc_cache: dict[tuple[str, int, int], tuple[float, int] | None] = {}

    target_peaks = []
    target_hashes: dict[str, str] = {}
    for row in manifest:
        path = Path(row["matched_target_path"])
        if not path.is_file():
            raise FileNotFoundError(path)
        observed_hash = sha256(path)
        if observed_hash != row["matched_target_sha256"]:
            raise RuntimeError(f"Immutable B0 target hash changed: {row['system']} {row['sample_id']}")
        peaks, _ = read_peak_file(path, chrom_sizes, blacklist, fasta, gc_cache)
        if len(peaks) != int(row["matched_peak_count"]):
            raise RuntimeError(f"B0 target count changed: {row['system']} {row['sample_id']}")
        target_peaks.extend(peaks)
        target_hashes[f"{row['system']}:{row['sample_id']}"] = observed_hash

    tcga_luad, tcga_stats = read_tcga_type_peaks(
        tcga_archive, "LUAD", chrom_sizes, blacklist, fasta, gc_cache
    )
    pools = {
        "B1": tcga_luad + target_peaks,
        "B2": tcga_luad,
        "B3": target_peaks,
    }
    out_root = root / "data/processed/motif_background_sensitivity"
    out_root.mkdir(parents=True, exist_ok=True)
    variant_receipts: dict[str, object] = {}
    for variant, pool in pools.items():
        consensus, merge_iterations = build_consensus(
            pool, chrom_sizes, blacklist, fasta, gc_cache
        )
        if not consensus:
            raise RuntimeError(f"{variant} background is empty")
        path = out_root / variant / f"{variant}_shared_accessible_background_200bp.bed"
        write_bed(path, consensus)
        count, width_ok = bed_count_and_width(path)
        if count != len(consensus) or not width_ok:
            raise RuntimeError(f"{variant} BED validation failed")
        gc_mean, gc_sd = summarize_gc(consensus)
        variant_receipts[variant] = {
            "definition": VARIANTS[variant],
            "pooled_input_peak_records": len(pool),
            "consensus_merge_iterations": merge_iterations,
            "background_peak_count": count,
            "background_gc_mean": gc_mean,
            "background_gc_sd": gc_sd,
            "background": str(path),
            "background_sha256": sha256(path),
            "all_intervals_exactly_200bp": width_ok,
        }

    b0_count, b0_width_ok = bed_count_and_width(b0_background)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "prespecified_shared_background_sensitivity_only",
        "immutable_primary_background": {
            "variant": "B0",
            "definition": "TCGA_LUAD_plus_TCGA_LUSC_plus_all_three_LUAD_target_systems",
            "background": str(b0_background),
            "background_sha256": sha256(b0_background),
            "background_peak_count": b0_count,
            "all_intervals_exactly_200bp": b0_width_ok,
        },
        "target_manifest": str(manifest_path),
        "target_manifest_sha256": sha256(manifest_path),
        "target_files_reused_without_rewriting": True,
        "target_count": len(manifest),
        "target_hashes": target_hashes,
        "tcga_luad_source_stats": tcga_stats,
        "blacklist_sha256": sha256(blacklist_path),
        "fasta_fai_sha256": sha256(fai_path),
        "variants": variant_receipts,
        "interpretation": "B1-B3 cannot replace B0; they quantify dependence of motif support on the shared accessible-universe construction.",
    }
    receipt_path = root / "audit/motif_background_sensitivity/background_input_receipt.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

