#!/usr/bin/env python3
"""Audit batch mixing and biological retention in a published CRC embedding.

This is deliberately an audit of the released integrated object.  Without a
pre-integration embedding generated from the same cells it cannot establish
that integration improved the data, so the receipt remains diagnostic-only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from join_crc_h5ad_table_s1 import clean_value, read_obs_column


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def read_dense_obsm(handle: h5py.File, key: str, rows: np.ndarray) -> np.ndarray:
    path = f"obsm/{key}"
    if path not in handle:
        raise ValueError(f"H5AD missing embedding: {path}")
    item = handle[path]
    if not isinstance(item, h5py.Dataset) or item.ndim != 2:
        raise ValueError(f"embedding must be a dense two-dimensional dataset: {path}")
    values = np.asarray(item[rows], dtype=np.float32)
    if not np.isfinite(values).all():
        raise ValueError(f"embedding contains non-finite values: {path}")
    return values


def stratified_indices(strata: pd.DataFrame, max_cells: int, seed: int) -> np.ndarray:
    if max_cells < 100:
        raise ValueError("max_cells must be at least 100")
    rng = np.random.default_rng(seed)
    groups = strata.groupby(list(strata.columns), sort=True, observed=True, dropna=False).indices
    if not groups:
        raise ValueError("no non-empty strata")
    cap = max(1, max_cells // len(groups))
    selected: list[np.ndarray] = []
    for indices in groups.values():
        values = np.asarray(indices, dtype=np.int64)
        if len(values) > cap:
            values = rng.choice(values, cap, replace=False)
        selected.append(values)
    merged = np.concatenate(selected)
    if len(merged) > max_cells:
        merged = rng.choice(merged, max_cells, replace=False)
    return np.sort(np.unique(merged))


def normalized_neighbor_entropy(values: np.ndarray, neighbor_indices: np.ndarray) -> np.ndarray:
    categories, encoded = np.unique(values, return_inverse=True)
    denominator = np.log(min(neighbor_indices.shape[1], len(categories)))
    if denominator <= 0:
        return np.zeros(len(values), dtype=float)
    out = np.zeros(len(values), dtype=float)
    for row, neighbors in enumerate(neighbor_indices):
        counts = np.bincount(encoded[neighbors], minlength=len(categories))
        probabilities = counts[counts > 0] / len(neighbors)
        out[row] = float(-(probabilities * np.log(probabilities)).sum() / denominator)
    return out


def safe_silhouette(embedding: np.ndarray, labels: np.ndarray, seed: int) -> float:
    _, counts = np.unique(labels, return_counts=True)
    if len(counts) < 2 or counts.min() < 2:
        return float("nan")
    sample_size = min(20000, len(labels))
    return float(silhouette_score(embedding, labels, sample_size=sample_size, random_state=seed))


def audit_integration(
    h5ad: Path,
    output_dir: Path,
    embedding_key: str = "X_scANVI",
    batch_key: str = "dataset",
    label_key: str = "cell_type_coarse",
    patient_key: str = "patient_id",
    max_cells: int = 50000,
    neighbors: int = 30,
    seed: int = 20260810,
    expected_obs: int | None = None,
) -> dict[str, object]:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD is missing or empty: {h5ad}")
    with h5py.File(h5ad, "r") as handle:
        if "obs" not in handle:
            raise ValueError("H5AD missing obs")
        obs = handle["obs"]
        raw = {
            batch_key: read_obs_column(obs, batch_key),
            label_key: read_obs_column(obs, label_key),
            patient_key: read_obs_column(obs, patient_key),
        }
        lengths = {len(values) for values in raw.values()}
        if len(lengths) != 1:
            raise ValueError(f"obs columns have inconsistent lengths: {lengths}")
        n_obs = lengths.pop()
        if expected_obs is not None and n_obs != expected_obs:
            raise ValueError(f"H5AD cell count mismatch: expected {expected_obs}, got {n_obs}")
        metadata = pd.DataFrame({key: pd.Series(values).map(clean_value) for key, values in raw.items()})
        if (metadata == "").any().any():
            missing = {column: int((metadata[column] == "").sum()) for column in metadata if (metadata[column] == "").any()}
            raise ValueError(f"integration audit metadata contain missing values: {missing}")
        selected = stratified_indices(metadata[[batch_key, label_key]], max_cells=max_cells, seed=seed)
        embedding = read_dense_obsm(handle, embedding_key, selected)

    sampled = metadata.iloc[selected].reset_index(drop=True)
    if sampled[batch_key].nunique() < 2 or sampled[label_key].nunique() < 2:
        raise ValueError("integration audit requires at least two batches and two biological labels")
    k = min(neighbors + 1, len(sampled))
    if k < 2:
        raise ValueError("too few sampled cells for neighbor audit")
    model = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    neighbor_indices = model.fit(embedding).kneighbors(embedding, return_distance=False)[:, 1:]
    batch = sampled[batch_key].to_numpy(dtype=str)
    label = sampled[label_key].to_numpy(dtype=str)
    patient = sampled[patient_key].to_numpy(dtype=str)
    neighbor_batch = batch[neighbor_indices]
    neighbor_label = label[neighbor_indices]
    neighbor_patient = patient[neighbor_indices]
    per_cell = pd.DataFrame(
        {
            "batch_different_neighbor_fraction": (neighbor_batch != batch[:, None]).mean(axis=1),
            "batch_neighbor_entropy": normalized_neighbor_entropy(batch, neighbor_indices),
            "label_same_neighbor_fraction": (neighbor_label == label[:, None]).mean(axis=1),
            "patient_different_neighbor_fraction": (neighbor_patient != patient[:, None]).mean(axis=1),
        }
    )
    metrics = [
        {"scope": "global", "metric": column, "value": float(per_cell[column].mean()), "n": len(per_cell)}
        for column in per_cell.columns
    ]
    metrics.extend(
        [
            {"scope": "global", "metric": "batch_silhouette", "value": safe_silhouette(embedding, batch, seed), "n": len(sampled)},
            {"scope": "global", "metric": "label_silhouette", "value": safe_silhouette(embedding, label, seed), "n": len(sampled)},
        ]
    )
    coverage = (
        sampled.groupby([batch_key, label_key], observed=True, sort=True)
        .agg(sampled_cells=(patient_key, "size"), patients=(patient_key, "nunique"))
        .reset_index()
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    metric_path = output_dir / "P0_integration_metrics.tsv"
    coverage_path = output_dir / "P0_integration_sampled_coverage.tsv"
    pd.DataFrame(metrics).to_csv(metric_path, sep="\t", index=False, lineterminator="\n")
    coverage.to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "diagnostic_only_no_preintegration_comparator",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "h5ad_sha256": sha256_file(h5ad),
        "n_obs": n_obs,
        "sampled_cells": len(sampled),
        "embedding_key": embedding_key,
        "batch_key": batch_key,
        "label_key": label_key,
        "patient_key": patient_key,
        "neighbors": neighbor_indices.shape[1],
        "seed": seed,
        "claim_boundary": "Metrics describe the released integrated embedding only; they do not prove improvement over unintegrated counts or successful target-project retraining.",
        "outputs": {},
    }
    for path in (metric_path, coverage_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (output_dir / "P0_integration_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--embedding-key", default="X_scANVI")
    parser.add_argument("--batch-key", default="dataset")
    parser.add_argument("--label-key", default="cell_type_coarse_crc_atlas")
    parser.add_argument("--patient-key", default="donor_id")
    parser.add_argument("--max-cells", type=int, default=50000)
    parser.add_argument("--neighbors", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--expected-obs", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = audit_integration(
        args.h5ad,
        args.output_dir,
        args.embedding_key,
        args.batch_key,
        args.label_key,
        args.patient_key,
        args.max_cells,
        args.neighbors,
        args.seed,
        args.expected_obs,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
