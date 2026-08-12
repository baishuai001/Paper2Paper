#!/usr/bin/env python3
"""Prepare a donor/dataset-balanced real CRC cNMF diagnostic input."""

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

from join_crc_h5ad_table_s1 import canonical_patient
from prepare_crc_cna_inputs import unique_gene_features


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def prepare_cnmf_input(
    h5ad: Path,
    output_dir: Path,
    patient_key: str = "donor_id",
    dataset_key: str = "dataset",
    sample_key: str = "sample_id",
    state_key: str = "sample_type",
    tissue_key: str = "tissue",
    coarse_label_key: str = "cell_type_coarse_crc_atlas",
    cancer_label: str = "Cancer cell",
    symbol_key: str = "GeneSymbol",
    max_datasets: int = 5,
    max_patients_per_dataset: int = 4,
    max_cells_per_patient: int = 100,
    min_cells_per_patient: int = 50,
    seed: int = 20260810,
) -> dict[str, object]:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD missing or empty: {h5ad}")
    for value, label in (
        (max_datasets, "max_datasets"),
        (max_patients_per_dataset, "max_patients_per_dataset"),
        (max_cells_per_patient, "max_cells_per_patient"),
        (min_cells_per_patient, "min_cells_per_patient"),
    ):
        if value < 1:
            raise ValueError(f"{label} must be positive")
    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        obs_fields = [patient_key, dataset_key, sample_key, state_key, tissue_key, coarse_label_key]
        missing = sorted(set(obs_fields) - set(adata.obs.columns))
        if missing:
            raise ValueError(f"H5AD obs missing cNMF fields: {missing}")
        if adata.raw is None:
            raise ValueError("H5AD raw is missing; cNMF input must not use normalized X")
        obs = adata.obs[obs_fields].astype(str).copy()
        obs["analysis_patient_id"] = obs[patient_key].map(canonical_patient)
        cancer = obs[coarse_label_key].eq(cancer_label)
        candidate_obs = obs.loc[cancer].copy()
        candidate_obs["row_index"] = np.flatnonzero(cancer.to_numpy())
        units = (
            candidate_obs.groupby([dataset_key, "analysis_patient_id"], observed=True, sort=True)
            .size()
            .rename("cancer_cells")
            .reset_index()
        )
        eligible = units[units["cancer_cells"] >= min_cells_per_patient].copy()
        if eligible.empty:
            raise ValueError("no dataset-patient unit meets the minimum cancer-cell count")
        dataset_rank = (
            eligible.groupby(dataset_key, observed=True)
            .agg(eligible_patients=("analysis_patient_id", "nunique"), eligible_cells=("cancer_cells", "sum"))
            .reset_index()
            .sort_values(["eligible_patients", "eligible_cells", dataset_key], ascending=[False, False, True])
        )
        selected_datasets = dataset_rank.head(max_datasets)[dataset_key].tolist()
        rng = np.random.default_rng(seed)
        selected_units: list[tuple[str, str]] = []
        for dataset in selected_datasets:
            subset = eligible[eligible[dataset_key] == dataset].sort_values("analysis_patient_id")
            patients = subset["analysis_patient_id"].to_numpy(dtype=str)
            if len(patients) > max_patients_per_dataset:
                patients = np.sort(rng.choice(patients, max_patients_per_dataset, replace=False))
            selected_units.extend((dataset, patient) for patient in patients)
        selected_rows: list[np.ndarray] = []
        for dataset, patient in selected_units:
            indices = candidate_obs.loc[
                (candidate_obs[dataset_key] == dataset)
                & (candidate_obs["analysis_patient_id"] == patient),
                "row_index",
            ].to_numpy(dtype=np.int64)
            if len(indices) > max_cells_per_patient:
                indices = rng.choice(indices, max_cells_per_patient, replace=False)
            selected_rows.append(np.sort(indices))
        rows = np.sort(np.concatenate(selected_rows))
        selected_obs = obs.iloc[rows].copy()
        selected_obs.insert(0, "source_cell_id", adata.obs_names[rows].astype(str))
        selected_obs.index = selected_obs["source_cell_id"].astype(str)
        raw_matrix = adata.raw.X[rows, :]
        if not sparse.issparse(raw_matrix):
            raw_matrix = sparse.csr_matrix(np.asarray(raw_matrix))
        raw_matrix = raw_matrix.tocsr()
        feature_map = unique_gene_features(adata.raw.var.copy(), symbol_key)
        feature_indices = feature_map["feature_index"].to_numpy(dtype=np.int64)
        counts = raw_matrix[:, feature_indices].tocsr()
        if counts.nnz and (counts.data < 0).any():
            raise ValueError("raw/X contains negative values")
        if counts.nnz and not np.allclose(counts.data, np.rint(counts.data), atol=1e-6):
            raise ValueError("raw/X contains non-integer values; refusing cNMF input")
        counts.data = np.rint(counts.data).astype(np.float32)
        counts.eliminate_zeros()
        var = pd.DataFrame(index=feature_map["gene"].astype(str).to_numpy())
        var.index.name = "gene"
        cnmf_input = ad.AnnData(X=counts, obs=selected_obs, var=var)
    finally:
        adata.file.close()

    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = output_dir / "crc_author_cancer_balanced_counts.h5ad"
    unit_path = output_dir / "P0_cnmf_selected_units.tsv"
    feature_path = output_dir / "P0_cnmf_gene_feature_map.tsv"
    cnmf_input.write_h5ad(input_path, compression="gzip")
    selected_unit_table = (
        selected_obs.groupby([dataset_key, "analysis_patient_id"], observed=True, sort=True)
        .agg(
            selected_cells=("source_cell_id", "size"),
            samples=(sample_key, "nunique"),
            states=(state_key, lambda values: "|".join(sorted(set(values)))),
            tissues=(tissue_key, lambda values: "|".join(sorted(set(values)))),
        )
        .reset_index()
    )
    selected_unit_table.to_csv(unit_path, sep="\t", index=False, lineterminator="\n")
    feature_map.to_csv(feature_path, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "diagnostic_author_cancer_input_pending_cna",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "raw_matrix": "raw/X",
        "cells": int(cnmf_input.n_obs),
        "genes": int(cnmf_input.n_vars),
        "datasets": int(selected_obs[dataset_key].nunique()),
        "patients": int(selected_obs["analysis_patient_id"].nunique()),
        "selected_units": len(selected_unit_table),
        "seed": seed,
        "selection_rule": "datasets ranked by eligible-patient coverage; patients sampled outcome-blind; equal cell cap per dataset-patient unit",
        "claim_boundary": "This diagnostic input uses the atlas author Cancer cell label. It must be replaced by P0.5 high-confidence malignant calls before scientific P0.6 acceptance.",
        "outputs": {},
    }
    for path in (input_path, unit_path, feature_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (output_dir / "P0_cnmf_input_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--patient-key", default="donor_id")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--state-key", default="sample_type")
    parser.add_argument("--tissue-key", default="tissue")
    parser.add_argument("--coarse-label-key", default="cell_type_coarse_crc_atlas")
    parser.add_argument("--cancer-label", default="Cancer cell")
    parser.add_argument("--symbol-key", default="GeneSymbol")
    parser.add_argument("--max-datasets", type=int, default=5)
    parser.add_argument("--max-patients-per-dataset", type=int, default=4)
    parser.add_argument("--max-cells-per-patient", type=int, default=100)
    parser.add_argument("--min-cells-per-patient", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260810)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = prepare_cnmf_input(
        args.h5ad,
        args.output_dir,
        args.patient_key,
        args.dataset_key,
        args.sample_key,
        args.state_key,
        args.tissue_key,
        args.coarse_label_key,
        args.cancer_label,
        args.symbol_key,
        args.max_datasets,
        args.max_patients_per_dataset,
        args.max_cells_per_patient,
        args.min_cells_per_patient,
        args.seed,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
