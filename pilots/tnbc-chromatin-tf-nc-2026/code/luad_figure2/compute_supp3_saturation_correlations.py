#!/usr/bin/env python3
"""Compute Figure S3 peak saturation and within-system accessibility correlations."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_presence(path: Path) -> tuple[list[str], np.ndarray]:
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        sample_ids = header[3:]
        values = np.loadtxt(handle, delimiter="\t", usecols=range(3, len(header)), dtype=np.uint8)
    if values.ndim == 1:
        values = values[:, None]
    return sample_ids, values.astype(bool, copy=False)


def saturation(system: str, sample_ids: list[str], presence: np.ndarray, permutations: int, seed: int):
    rng = np.random.default_rng(seed)
    raw_rows: list[dict[str, object]] = []
    sums = np.zeros(len(sample_ids), dtype=np.float64)
    sums_squared = np.zeros(len(sample_ids), dtype=np.float64)
    for permutation_index in range(1, permutations + 1):
        order = rng.permutation(len(sample_ids))
        covered = np.zeros(presence.shape[0], dtype=bool)
        for k, column in enumerate(order, start=1):
            np.logical_or(covered, presence[:, column], out=covered)
            count = int(covered.sum())
            sums[k - 1] += count
            sums_squared[k - 1] += count * count
            raw_rows.append(
                {
                    "system": system,
                    "permutation": permutation_index,
                    "sample_number": k,
                    "cumulative_consensus_peaks": count,
                }
            )
    means = sums / permutations
    variances = np.maximum(0, (sums_squared - permutations * means * means) / max(1, permutations - 1))
    summary_rows = [
        {
            "system": system,
            "sample_number": index + 1,
            "mean_cumulative_consensus_peaks": means[index],
            "sem_cumulative_consensus_peaks": math.sqrt(variances[index] / permutations),
            "observed_full_consensus_peaks": presence.shape[0],
            "permutations": permutations,
            "seed": seed,
        }
        for index in range(len(sample_ids))
    ]
    return raw_rows, summary_rows


def raw_count_correlations(system: str, count_path: Path, sample_ids: list[str]):
    values = np.loadtxt(count_path, delimiter="\t", usecols=range(3, 3 + len(sample_ids)), dtype=np.float64)
    if values.ndim == 1:
        values = values[:, None]
    library_sizes = values.sum(axis=0)
    if np.any(library_sizes <= 0):
        raise RuntimeError(f"{system} consensus count matrix contains zero library")
    log_cpm = np.log2(values / library_sizes[None, :] * 1e6 + 1.0)
    correlation = np.corrcoef(log_cpm, rowvar=False)
    return [
        {"system": system, "sample_1": sample_ids[i], "sample_2": sample_ids[j], "pearson_r": correlation[i, j]}
        for i in range(len(sample_ids))
        for j in range(len(sample_ids))
    ], library_sizes


def patient_correlations(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_ids = header[1:]
        matrix = np.array([[float(value) for value in row[1:]] for row in reader], dtype=np.float64)
    if matrix.shape != (len(sample_ids), len(sample_ids)):
        raise RuntimeError("Unexpected TCGA correlation matrix dimensions")
    return [
        {"system": "patient", "sample_1": sample_ids[i], "sample_2": sample_ids[j], "pearson_r": matrix[i, j]}
        for i in range(len(sample_ids))
        for j in range(len(sample_ids))
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260814)
    args = parser.parse_args()
    root = args.run_root.resolve()
    processed = root / "data/processed/supplementary_figure3"
    out = root / "results/supplementary_figure3"
    out.mkdir(parents=True, exist_ok=True)

    all_raw: list[dict[str, object]] = []
    all_summary: list[dict[str, object]] = []
    all_correlations = patient_correlations(
        root / "data/processed/tcga_luad/tcga_luad_sample_pearson_correlation.tsv.gz"
    )
    system_receipts: dict[str, dict[str, object]] = {}
    for system in ("patient", "PDX", "cell_line"):
        sample_ids, presence = load_presence(processed / f"{system}_consensus_peak_presence.tsv")
        raw, summary = saturation(system, sample_ids, presence, args.permutations, args.seed)
        all_raw.extend(raw)
        all_summary.extend(summary)
        receipt = {"samples": len(sample_ids), "consensus_peaks": int(presence.shape[0])}
        if system != "patient":
            correlations, library_sizes = raw_count_correlations(
                system,
                processed / f"{system}_consensus_peak_read_counts.tsv",
                sample_ids,
            )
            all_correlations.extend(correlations)
            receipt["consensus_read_count_min"] = float(library_sizes.min())
            receipt["consensus_read_count_median"] = float(np.median(library_sizes))
            receipt["consensus_read_count_max"] = float(library_sizes.max())
        system_receipts[system] = receipt

    write_tsv(
        out / "saturation_1000_permutations.tsv",
        all_raw,
        ["system", "permutation", "sample_number", "cumulative_consensus_peaks"],
    )
    write_tsv(
        out / "saturation_summary.tsv",
        all_summary,
        ["system", "sample_number", "mean_cumulative_consensus_peaks", "sem_cumulative_consensus_peaks", "observed_full_consensus_peaks", "permutations", "seed"],
    )
    write_tsv(
        out / "sample_accessibility_pearson_correlations.tsv",
        all_correlations,
        ["system", "sample_1", "sample_2", "pearson_r"],
    )
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "permutations": args.permutations,
        "seed": args.seed,
        "saturation_definition": "cumulative full-cohort consensus peaks overlapped by at least one added sample",
        "raw_system_correlation_definition": "Pearson correlation of log2(consensus-peak CPM + 1)",
        "patient_correlation_definition": "Pearson correlation of official GDC log2(merged technical-replicate CPM + 1)",
        "systems": system_receipts,
    }
    (root / "audit/supplementary_figure3/saturation_correlation_receipt.json").parent.mkdir(parents=True, exist_ok=True)
    (root / "audit/supplementary_figure3/saturation_correlation_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
