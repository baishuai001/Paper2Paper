#!/usr/bin/env python3
"""Audit author labels against author-provided marker genes across donors/datasets."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def marker_group_for_label(label: str, marker_groups: set[str]) -> str:
    if label in marker_groups:
        return label
    if label.startswith("Tumor ") or label.startswith("Cancer "):
        return "Cancer cell" if "Cancer cell" in marker_groups else ""
    if label.startswith("Endothelial ") and "Pan-Endothelial" in marker_groups:
        return "Pan-Endothelial"
    if label.startswith("Fibroblast ") and "Pan-Fibroblast" in marker_groups:
        return "Pan-Fibroblast"
    if "cycling" in label.lower() and "Dividing" in marker_groups:
        return "Dividing"
    return ""


def sample_balanced_indices(obs: pd.DataFrame, patient_key: str, label_key: str, cap: int, seed: int) -> np.ndarray:
    if cap < 1:
        raise ValueError("per_patient_label_cap must be positive")
    rng = np.random.default_rng(seed)
    groups = obs.groupby([patient_key, label_key], observed=True, sort=True, dropna=False).indices
    selected: list[np.ndarray] = []
    for values in groups.values():
        indices = np.asarray(values, dtype=np.int64)
        if len(indices) > cap:
            indices = rng.choice(indices, cap, replace=False)
        selected.append(indices)
    if not selected:
        raise ValueError("no patient-label groups available for marker audit")
    return np.sort(np.concatenate(selected))


def gene_symbols(adata: ad.AnnData, symbol_key: str) -> pd.Series:
    if symbol_key in adata.var.columns:
        symbols = adata.var[symbol_key].astype(str)
    else:
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var_names)
    symbols = symbols.str.strip()
    return symbols


def resolve_marker_features(
    adata: ad.AnnData,
    symbols: pd.Series,
    marker_genes: list[str],
) -> tuple[dict[str, int], pd.DataFrame]:
    if "n_cells" in adata.var.columns:
        feature_score = pd.to_numeric(adata.var["n_cells"], errors="coerce").fillna(-1).to_numpy()
        criterion = "maximum_var_n_cells_then_lowest_feature_index"
    else:
        feature_score = np.zeros(adata.n_vars, dtype=float)
        criterion = "lowest_feature_index_no_var_n_cells_available"
    symbol_values = symbols.to_numpy(dtype=str)
    resolved: dict[str, int] = {}
    rows: list[dict[str, object]] = []
    for gene in sorted(set(marker_genes)):
        candidates = np.flatnonzero(symbol_values == gene)
        if len(candidates) == 0:
            continue
        order = sorted(candidates.tolist(), key=lambda index: (-feature_score[index], index))
        selected = int(order[0])
        resolved[gene] = selected
        rows.append(
            {
                "gene": gene,
                "candidate_features": len(candidates),
                "selected_feature_index": selected,
                "selected_var_name": str(adata.var_names[selected]),
                "selected_n_cells": float(feature_score[selected]),
                "selection_rule": criterion,
            }
        )
    return resolved, pd.DataFrame(rows)


def extract_expression(
    adata: ad.AnnData,
    row_indices: np.ndarray,
    gene_indices: np.ndarray,
    row_chunk_size: int = 256,
) -> np.ndarray:
    """Read selected genes without issuing two HDF5 fancy indices at once."""
    if row_chunk_size < 1:
        raise ValueError("row_chunk_size must be positive")
    matrix = adata.X
    output = np.empty((len(row_indices), len(gene_indices)), dtype=np.float32)
    for start in range(0, len(row_indices), row_chunk_size):
        stop = min(start + row_chunk_size, len(row_indices))
        rows = row_indices[start:stop]
        block = matrix[rows, :]
        if sparse.issparse(block):
            selected = block[:, gene_indices].toarray()
        else:
            selected = np.asarray(block)[:, gene_indices]
        output[start:stop] = np.asarray(selected, dtype=np.float32)
    return output


def audit_markers(
    h5ad: Path,
    marker_table: Path,
    output_dir: Path,
    label_key: str = "cell_type_fine",
    patient_key: str = "patient_id",
    dataset_key: str = "dataset",
    symbol_key: str = "symbol",
    per_patient_label_cap: int = 30,
    seed: int = 20260810,
    expected_obs: int | None = None,
) -> dict[str, object]:
    for path in (h5ad, marker_table):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"input missing or empty: {path}")
    markers = pd.read_csv(marker_table, dtype=str, keep_default_na=False)
    missing_marker_columns = sorted({"cell_type", "symbol"} - set(markers.columns))
    if missing_marker_columns:
        raise ValueError(f"marker table missing columns: {missing_marker_columns}")
    markers["cell_type"] = markers["cell_type"].str.strip()
    markers["symbol"] = markers["symbol"].str.strip()
    markers = markers[(markers["cell_type"] != "") & (markers["symbol"] != "")].drop_duplicates(
        ["cell_type", "symbol"]
    )
    if markers.empty:
        raise ValueError("marker table contains no usable marker rows")

    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        if expected_obs is not None and adata.n_obs != expected_obs:
            raise ValueError(f"H5AD cell count mismatch: expected {expected_obs}, got {adata.n_obs}")
        missing_obs = sorted({label_key, patient_key, dataset_key} - set(adata.obs.columns))
        if missing_obs:
            raise ValueError(f"H5AD obs missing marker-audit fields: {missing_obs}")
        obs = adata.obs[[label_key, patient_key, dataset_key]].astype(str).apply(lambda column: column.str.strip())
        if (obs == "").any().any():
            counts = {column: int((obs[column] == "").sum()) for column in obs if (obs[column] == "").any()}
            raise ValueError(f"marker audit metadata contain missing values: {counts}")
        selected = sample_balanced_indices(obs, patient_key, label_key, per_patient_label_cap, seed)
        sampled_obs = obs.iloc[selected].reset_index(drop=True)
        marker_groups = set(markers["cell_type"])
        sampled_obs["marker_group"] = sampled_obs[label_key].map(lambda value: marker_group_for_label(value, marker_groups))
        symbols = gene_symbols(adata, symbol_key)
        symbol_to_index, gene_map = resolve_marker_features(adata, symbols, markers["symbol"].tolist())
        available = markers[markers["symbol"].isin(symbol_to_index)].copy()
        missing_genes = sorted(set(markers["symbol"]) - set(available["symbol"]))
        gene_order = sorted(available["symbol"].unique())
        if not gene_order:
            raise ValueError("none of the marker genes are present in H5AD var")
        gene_indices = np.array([symbol_to_index[gene] for gene in gene_order], dtype=np.int64)
        expression = extract_expression(adata, selected, gene_indices)
        if expression.shape != (len(selected), len(gene_order)):
            raise ValueError(f"unexpected marker expression shape: {expression.shape}")
        if not np.isfinite(expression).all() or (expression < 0).any():
            raise ValueError("marker audit expression must be finite and nonnegative")
    finally:
        adata.file.close()

    gene_to_col = {gene: column for column, gene in enumerate(gene_order)}
    rows: list[dict[str, object]] = []
    for marker_group, group_markers in available.groupby("cell_type", sort=True, observed=True):
        target = sampled_obs["marker_group"].to_numpy() == marker_group
        if not target.any():
            continue
        background = ~target
        target_patients = sampled_obs.loc[target, patient_key].nunique()
        target_datasets = sampled_obs.loc[target, dataset_key].nunique()
        for gene in sorted(group_markers["symbol"].unique()):
            values = expression[:, gene_to_col[gene]]
            target_mean = float(values[target].mean())
            background_mean = float(values[background].mean()) if background.any() else float("nan")
            target_detection = float((values[target] > 0).mean())
            background_detection = float((values[background] > 0).mean()) if background.any() else float("nan")
            rows.append(
                {
                    "marker_group": marker_group,
                    "gene": gene,
                    "target_cells": int(target.sum()),
                    "target_patients": int(target_patients),
                    "target_datasets": int(target_datasets),
                    "target_mean_expression": target_mean,
                    "background_mean_expression": background_mean,
                    "target_detection_fraction": target_detection,
                    "background_detection_fraction": background_detection,
                    "log2_mean_ratio": float(np.log2((target_mean + 1e-3) / (background_mean + 1e-3))) if background.any() else float("nan"),
                    "multi_patient": bool(target_patients >= 2),
                    "multi_dataset": bool(target_datasets >= 2),
                }
            )
    evidence = pd.DataFrame(rows)
    if evidence.empty:
        raise ValueError("no author labels could be mapped to available marker groups")
    label_coverage = (
        sampled_obs.groupby([label_key, "marker_group"], observed=True, sort=True)
        .agg(sampled_cells=(patient_key, "size"), patients=(patient_key, "nunique"), datasets=(dataset_key, "nunique"))
        .reset_index()
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "P0_annotation_evidence.tsv"
    coverage_path = output_dir / "P0_annotation_sampled_coverage.tsv"
    missing_path = output_dir / "P0_annotation_missing_markers.tsv"
    gene_map_path = output_dir / "P0_annotation_gene_map.tsv"
    evidence.to_csv(evidence_path, sep="\t", index=False, lineterminator="\n")
    label_coverage.to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"missing_marker_gene": missing_genes}).to_csv(missing_path, sep="\t", index=False, lineterminator="\n")
    gene_map.to_csv(gene_map_path, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "diagnostic_only_author_labels_not_overwritten",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "h5ad_sha256": sha256_file(h5ad),
        "marker_table": str(marker_table.resolve()),
        "marker_table_sha256": sha256_file(marker_table),
        "n_obs": int(expected_obs if expected_obs is not None else len(obs)),
        "sampled_cells": len(sampled_obs),
        "sampled_patients": int(sampled_obs[patient_key].nunique()),
        "sampled_datasets": int(sampled_obs[dataset_key].nunique()),
        "available_marker_genes": len(gene_order),
        "missing_marker_genes": len(missing_genes),
        "duplicate_marker_symbols_resolved": int((gene_map["candidate_features"] > 1).sum()),
        "seed": seed,
        "claim_boundary": "This cross-patient marker audit supports or challenges author labels; it does not independently establish malignancy and does not replace author_label with reviewed_label.",
        "outputs": {},
    }
    for path in (evidence_path, coverage_path, missing_path, gene_map_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (output_dir / "P0_annotation_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--marker-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--label-key", default="cell_type_fine_crc_atlas")
    parser.add_argument("--patient-key", default="donor_id")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--symbol-key", default="GeneSymbol")
    parser.add_argument("--per-patient-label-cap", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--expected-obs", type=int)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = audit_markers(
        args.h5ad,
        args.marker_table,
        args.output_dir,
        args.label_key,
        args.patient_key,
        args.dataset_key,
        args.symbol_key,
        args.per_patient_label_cap,
        args.seed,
        args.expected_obs,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
