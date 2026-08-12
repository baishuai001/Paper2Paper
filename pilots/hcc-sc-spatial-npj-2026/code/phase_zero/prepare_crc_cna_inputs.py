#!/usr/bin/env python3
"""Prepare one patient/sample raw-count input for CopyKAT and SCEVAN."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from join_crc_h5ad_table_s1 import canonical_patient


DEFAULT_REFERENCE_LABELS = [
    "T cell",
    "B cell",
    "Plasma cell",
    "Myeloid cell",
    "Neutrophil",
    "NK",
    "Mast cell",
    "ILC",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def choose_indices(indices: np.ndarray, cap: int, rng: np.random.Generator) -> np.ndarray:
    values = np.asarray(indices, dtype=np.int64)
    if len(values) > cap:
        values = rng.choice(values, cap, replace=False)
    return np.sort(values)


def unique_gene_features(var: pd.DataFrame, symbol_key: str) -> pd.DataFrame:
    if symbol_key not in var.columns:
        raise ValueError(f"raw.var missing gene symbol field: {symbol_key}")
    frame = pd.DataFrame(
        {
            "feature_index": np.arange(len(var), dtype=np.int64),
            "var_name": var.index.astype(str),
            "gene": var[symbol_key].astype(str).str.strip().to_numpy(),
        }
    )
    if "n_cells" in var.columns:
        frame["n_cells"] = pd.to_numeric(var["n_cells"], errors="coerce").fillna(-1).to_numpy()
        rule = "maximum_var_n_cells_then_lowest_feature_index"
    else:
        frame["n_cells"] = -1.0
        rule = "lowest_feature_index_no_var_n_cells_available"
    frame = frame[frame["gene"] != ""].copy()
    frame = frame.sort_values(["gene", "n_cells", "feature_index"], ascending=[True, False, True])
    frame["candidate_features"] = frame.groupby("gene", observed=True)["gene"].transform("size")
    selected = frame.drop_duplicates("gene", keep="first").copy()
    selected["selection_rule"] = rule
    return selected.sort_values("feature_index").reset_index(drop=True)


def prepare_cna_input(
    h5ad: Path,
    output_dir: Path,
    sample_id: str,
    patient_id: str,
    sample_key: str = "sample_id",
    patient_key: str = "donor_id",
    coarse_label_key: str = "cell_type_coarse_crc_atlas",
    symbol_key: str = "GeneSymbol",
    reference_labels: list[str] | None = None,
    cancer_label: str = "Cancer cell",
    epithelial_label: str = "Epithelial cell",
    max_cancer: int = 300,
    max_epithelial: int = 100,
    max_reference: int = 300,
    min_cancer: int = 50,
    min_reference: int = 50,
    seed: int = 20260810,
) -> dict[str, object]:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD missing or empty: {h5ad}")
    reference_labels = reference_labels or DEFAULT_REFERENCE_LABELS
    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        required_obs = {sample_key, patient_key, coarse_label_key}
        missing = sorted(required_obs - set(adata.obs.columns))
        if missing:
            raise ValueError(f"H5AD obs missing CNA fields: {missing}")
        if adata.raw is None:
            raise ValueError("H5AD raw is missing; CNA input must not use normalized X")
        obs = adata.obs[[sample_key, patient_key, coarse_label_key]].astype(str)
        canonical_patients = obs[patient_key].map(canonical_patient)
        target = (obs[sample_key] == sample_id) & (canonical_patients == canonical_patient(patient_id))
        if not target.any():
            raise ValueError(f"no cells found for sample={sample_id}, patient={patient_id}")
        labels = obs[coarse_label_key]
        cancer_indices = np.flatnonzero(target.to_numpy() & labels.eq(cancer_label).to_numpy())
        epithelial_indices = np.flatnonzero(target.to_numpy() & labels.eq(epithelial_label).to_numpy())
        reference_indices = np.flatnonzero(target.to_numpy() & labels.isin(reference_labels).to_numpy())
        if len(cancer_indices) < min_cancer:
            raise ValueError(f"too few author cancer cells: {len(cancer_indices)} < {min_cancer}")
        if len(reference_indices) < min_reference:
            raise ValueError(f"too few known-normal reference cells: {len(reference_indices)} < {min_reference}")
        rng = np.random.default_rng(seed)
        selected_by_role = {
            "author_cancer": choose_indices(cancer_indices, max_cancer, rng),
            "other_epithelial": choose_indices(epithelial_indices, max_epithelial, rng),
            "known_normal_reference": choose_indices(reference_indices, max_reference, rng),
        }
        role_by_index = {
            int(index): role for role, indices in selected_by_role.items() for index in indices
        }
        selected_rows = np.array(sorted(role_by_index), dtype=np.int64)
        if len(selected_rows) != len(role_by_index):
            raise ValueError("CNA role selections overlap unexpectedly")
        raw_matrix = adata.raw.X[selected_rows, :]
        if not sparse.issparse(raw_matrix):
            raw_matrix = sparse.csr_matrix(np.asarray(raw_matrix))
        raw_matrix = raw_matrix.tocsr()
        feature_map = unique_gene_features(adata.raw.var.copy(), symbol_key)
        selected_features = feature_map["feature_index"].to_numpy(dtype=np.int64)
        counts = raw_matrix[:, selected_features].tocsr()
        if counts.nnz and (counts.data < 0).any():
            raise ValueError("raw/X contains negative values")
        if counts.nnz and not np.allclose(counts.data, np.rint(counts.data), atol=1e-6):
            raise ValueError("raw/X contains non-integer values; refusing CNA input")
        counts.data = np.rint(counts.data).astype(np.int64)
        counts.eliminate_zeros()
        cell_metadata = obs.iloc[selected_rows].copy()
        cell_metadata.insert(0, "cell_id", adata.obs_names[selected_rows].astype(str))
        cell_metadata["analysis_patient_id"] = canonical_patients.iloc[selected_rows].to_numpy()
        cell_metadata["cna_input_role"] = [role_by_index[int(index)] for index in selected_rows]
        cell_metadata["is_known_normal_reference"] = cell_metadata["cna_input_role"].eq(
            "known_normal_reference"
        )
    finally:
        adata.file.close()

    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = output_dir / "counts_genes_by_cells.mtx.gz"
    genes_path = output_dir / "genes.tsv"
    cells_path = output_dir / "cells.tsv"
    metadata_path = output_dir / "cell_metadata.tsv"
    feature_map_path = output_dir / "gene_feature_map.tsv"
    with gzip.open(matrix_path, "wb") as handle:
        mmwrite(handle, counts.T.tocoo(), field="integer", symmetry="general")
    feature_map[["gene"]].to_csv(genes_path, sep="\t", index=False, header=False, lineterminator="\n")
    cell_metadata[["cell_id"]].to_csv(cells_path, sep="\t", index=False, header=False, lineterminator="\n")
    cell_metadata.to_csv(metadata_path, sep="\t", index=False, lineterminator="\n")
    feature_map.to_csv(feature_map_path, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "prepared_real_sample_input",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "sample_id": sample_id,
        "patient_id": canonical_patient(patient_id),
        "raw_matrix": "raw/X",
        "genes": int(counts.shape[1]),
        "cells": int(counts.shape[0]),
        "nonzero_counts": int(counts.nnz),
        "role_counts": cell_metadata["cna_input_role"].value_counts().sort_index().to_dict(),
        "reference_labels": reference_labels,
        "seed": seed,
        "duplicate_gene_symbols_resolved": int((feature_map["candidate_features"] > 1).sum()),
        "outputs": {},
    }
    for path in (matrix_path, genes_path, cells_path, metadata_path, feature_map_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (output_dir / "P0_cna_input_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--patient-id", required=True)
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--patient-key", default="donor_id")
    parser.add_argument("--coarse-label-key", default="cell_type_coarse_crc_atlas")
    parser.add_argument("--symbol-key", default="GeneSymbol")
    parser.add_argument("--reference-label", action="append", dest="reference_labels")
    parser.add_argument("--max-cancer", type=int, default=300)
    parser.add_argument("--max-epithelial", type=int, default=100)
    parser.add_argument("--max-reference", type=int, default=300)
    parser.add_argument("--min-cancer", type=int, default=50)
    parser.add_argument("--min-reference", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260810)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = prepare_cna_input(
        args.h5ad,
        args.output_dir,
        args.sample_id,
        args.patient_id,
        args.sample_key,
        args.patient_key,
        args.coarse_label_key,
        args.symbol_key,
        args.reference_labels,
        max_cancer=args.max_cancer,
        max_epithelial=args.max_epithelial,
        max_reference=args.max_reference,
        min_cancer=args.min_cancer,
        min_reference=args.min_reference,
        seed=args.seed,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
