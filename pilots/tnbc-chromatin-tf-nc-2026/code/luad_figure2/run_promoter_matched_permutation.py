#!/usr/bin/env python3
"""Matched all-protein-coding-gene permutation for the three-system promoter overlap."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import random
import re
import statistics
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pysam

from compute_promoter_gate_and_annotations import (
    IntervalIndex,
    merge_interval_map,
    parse_attributes,
    read_peaks,
    read_tsv,
    write_tsv,
)


PERMUTATIONS = 10_000
SEED = 20260814
MATCH_POOL = 100
PRIMARY_DEFINITION = "PRIMARY_anchor_promoter_tcga_cpm1_both"


def parse_all_promoters(gtf_path: Path):
    values: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(lambda: defaultdict(list))
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            chrom, start_text, end_text, strand = fields[0], fields[3], fields[4], fields[6]
            if chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX"}:
                continue
            attrs = parse_attributes(fields[8])
            if (attrs.get("gene_type") or attrs.get("gene_biotype")) != "protein_coding":
                continue
            symbol = attrs.get("gene_name", "")
            if not symbol:
                continue
            start = int(start_text) - 1
            end = int(end_text)
            tss = start if strand == "+" else end - 1
            promoter = (max(0, tss - 2500), tss + 1000) if strand == "+" else (max(0, tss - 1000), tss + 2500)
            values[symbol][chrom].append(promoter)
    return {gene: merge_interval_map(chrom_map) for gene, chrom_map in values.items()}


def load_samples(root: Path):
    patient_rows = read_tsv(root / "audit/tcga_luad_accessibility/tcga_luad_peak_file_manifest.tsv")
    samples = [
        ("patient", row["sample_id"], Path(row["peak_path"]))
        for row in patient_rows if row["definition"] == "cpm1_both_reps"
    ]
    raw_rows = read_tsv(root / "audit/raw_atac_qc/figure2_qualified_raw_atac_samples.tsv")
    samples.extend((row["system"], row["sample_id"], Path(row["narrowpeak"])) for row in raw_rows)
    observed = {system: sum(item[0] == system for item in samples) for system in ("patient", "PDX", "cell_line")}
    if observed["patient"] != 22 or observed["PDX"] < 10 or observed["cell_line"] < 10:
        raise RuntimeError(f"Matched permutation sample count below frozen minimum: {observed}")
    return samples


def all_gene_accessibility(promoters, samples):
    genes = sorted(promoters)
    counts = {system: {gene: 0 for gene in genes} for system in ("patient", "PDX", "cell_line")}
    sample_counts = {system: 0 for system in counts}
    for system, sample_id, peak_path in samples:
        peaks_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for chrom, start, end in read_peaks(peak_path):
            peaks_by_chrom[chrom].append((start, end))
        index = IntervalIndex(merge_interval_map(peaks_by_chrom))
        sample_counts[system] += 1
        for gene in genes:
            if any(index.overlaps(chrom, start, end) for chrom, values in promoters[gene].items() for start, end in values):
                counts[system][gene] += 1
    thresholds = {system: math.ceil(value / 2) for system, value in sample_counts.items()}
    return counts, sample_counts, thresholds


def promoter_metrics(promoters, fasta_path: Path):
    fasta = pysam.FastaFile(str(fasta_path))
    metrics: dict[str, dict[str, float]] = {}
    for gene, chrom_map in promoters.items():
        length = 0
        gc = 0
        valid = 0
        for chrom, intervals in chrom_map.items():
            for start, end in intervals:
                sequence = fasta.fetch(chrom, start, end).upper()
                length += len(sequence)
                gc += sequence.count("G") + sequence.count("C")
                valid += sum(sequence.count(base) for base in "ACGT")
        metrics[gene] = {
            "promoter_length": float(length),
            "promoter_gc": gc / valid if valid else math.nan,
        }
    fasta.close()
    return metrics


def tcga_mean_expression(path: Path, wanted: set[str]):
    library_sizes: np.ndarray | None = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        n_samples = len(header) - 1
        library_sizes = np.zeros(n_samples, dtype=np.float64)
        for line in handle:
            _gene, values = line.rstrip("\n").split("\t", 1)
            library_sizes += np.fromstring(values, dtype=np.float64, sep="\t")
    if np.any(library_sizes <= 0):
        raise RuntimeError("TCGA LUAD expression matrix has invalid library sizes")
    output: dict[str, float] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            gene, values = line.rstrip("\n").split("\t", 1)
            if gene not in wanted:
                continue
            counts = np.fromstring(values, dtype=np.float64, sep="\t")
            output[gene] = float(np.mean(np.log2(counts / library_sizes * 1e6 + 1.0)))
    return output


def zscore(values: dict[str, float]):
    numbers = list(values.values())
    mean = statistics.mean(numbers)
    sd = statistics.stdev(numbers)
    return {key: (value - mean) / sd for key, value in values.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    project_root = args.project_root.resolve()
    promoters = parse_all_promoters(root / "reference/downloads/gencode.v47.basic.annotation.gtf.gz")
    samples = load_samples(root)
    counts, sample_counts, thresholds = all_gene_accessibility(promoters, samples)
    metrics = promoter_metrics(promoters, root / "reference/hg38.fa")
    expression_path = project_root / "tmp/tnbc-chromatin-tf-nc-2026/luad-figure1/data/processed/tcga_counts_symbols.tsv.gz"
    expression = tcga_mean_expression(expression_path, set(promoters))
    eligible = sorted(
        gene for gene in promoters
        if gene in expression
        and metrics[gene]["promoter_length"] > 0
        and math.isfinite(metrics[gene]["promoter_gc"])
    )
    if len(eligible) < 10_000:
        raise RuntimeError(f"Unexpectedly small matched-gene universe: {len(eligible)}")

    length_z = zscore({gene: math.log2(metrics[gene]["promoter_length"]) for gene in eligible})
    gc_z = zscore({gene: metrics[gene]["promoter_gc"] for gene in eligible})
    expression_z = zscore({gene: expression[gene] for gene in eligible})
    triple_open = {
        gene: all(counts[system][gene] >= thresholds[system] for system in thresholds)
        for gene in eligible
    }
    summary_rows = [
        {
            "gene": gene,
            "promoter_length": metrics[gene]["promoter_length"],
            "promoter_gc": metrics[gene]["promoter_gc"],
            "TCGA_LUAD_mean_log2_CPM_plus1": expression[gene],
            "patient_accessible_n": counts["patient"][gene],
            "PDX_accessible_n": counts["PDX"][gene],
            "cell_line_accessible_n": counts["cell_line"][gene],
            "triple_system_promoter_accessible": str(triple_open[gene]).upper(),
        }
        for gene in eligible
    ]
    write_tsv(root / "audit/promoter_gate/all_gene_matched_universe.tsv", summary_rows)

    promoter_gate = read_tsv(root / "results/promoter_gate/tf_promoter_gate_summary.tsv")
    frozen = sorted(row["TF"] for row in promoter_gate if row["analysis_definition"] == PRIMARY_DEFINITION)
    missing = [gene for gene in frozen if gene not in eligible]
    if missing:
        raise RuntimeError(f"Frozen TFs missing from matched universe: {missing}")
    frozen_set = set(frozen)
    candidates = [gene for gene in eligible if gene not in frozen_set]
    coordinates = np.array([[length_z[gene], gc_z[gene], expression_z[gene]] for gene in candidates])
    match_pools: dict[str, list[tuple[str, float]]] = {}
    match_rows: list[dict[str, object]] = []
    for tf in frozen:
        target = np.array([length_z[tf], gc_z[tf], expression_z[tf]])
        distances = np.sqrt(((coordinates - target[None, :]) ** 2).sum(axis=1))
        indexes = np.argpartition(distances, MATCH_POOL)[:MATCH_POOL]
        ordered = sorted(((candidates[index], float(distances[index])) for index in indexes), key=lambda item: item[1])
        match_pools[tf] = ordered
        for rank, (gene, distance) in enumerate(ordered, start=1):
            match_rows.append({"frozen_TF": tf, "rank": rank, "matched_gene": gene, "standardized_euclidean_distance": distance})
    write_tsv(root / "audit/promoter_gate/frozen_tf_matched_gene_pools.tsv", match_rows)

    observed = sum(triple_open[tf] for tf in frozen)
    rng = random.Random(SEED)
    null_counts: list[int] = []
    match_distance_means: list[float] = []
    for _ in range(PERMUTATIONS):
        selected: set[str] = set()
        distances_selected: list[float] = []
        count = 0
        order = frozen.copy()
        rng.shuffle(order)
        for tf in order:
            available = [(gene, distance) for gene, distance in match_pools[tf] if gene not in selected]
            if not available:
                raise RuntimeError(f"Matched pool exhausted for {tf}")
            top = available[: min(25, len(available))]
            gene, distance = rng.choice(top)
            selected.add(gene)
            distances_selected.append(distance)
            count += int(triple_open[gene])
        null_counts.append(count)
        match_distance_means.append(statistics.mean(distances_selected))
    empirical_p = (1 + sum(value >= observed for value in null_counts)) / (PERMUTATIONS + 1)
    write_tsv(
        root / "audit/promoter_gate/promoter_matched_permutation_null.tsv",
        [
            {"permutation": index + 1, "triple_system_promoter_genes": value, "mean_match_distance": match_distance_means[index]}
            for index, value in enumerate(null_counts)
        ],
        ["permutation", "triple_system_promoter_genes", "mean_match_distance"],
    )
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "observed_frozen_TFs": len(frozen),
        "observed_triple_system_promoter_TFs": observed,
        "eligible_protein_coding_gene_universe": len(eligible),
        "permutations": PERMUTATIONS,
        "seed": SEED,
        "matching_variables": ["promoter GC", "log2 promoter union length", "TCGA-LUAD mean log2(CPM+1)"],
        "match_pool_per_TF": MATCH_POOL,
        "draw_rule": "without replacement; random among closest 25 still available",
        "null_mean": statistics.mean(null_counts),
        "null_sd": statistics.stdev(null_counts),
        "null_max": max(null_counts),
        "mean_standardized_match_distance": statistics.mean(match_distance_means),
        "empirical_p_ge_observed": empirical_p,
        "robustness_pass": empirical_p < 0.05,
    }
    (root / "audit/promoter_gate/promoter_matched_permutation_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
