#!/usr/bin/env python3
"""Build a balanced cNMF input from two-method CNA-supported CRC cells."""

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


REQUIRED_EVIDENCE_COLUMNS = {
    "cell_id",
    "panel_order",
    "panel_dataset",
    "panel_sample_id",
    "panel_patient_id",
    "panel_sample_type",
    "panel_tissue",
    "cna_input_role",
    "copykat.pred",
    "scevan_class",
    "cna_selection_class",
}
SELECTED_CLASS = "author_cancer_two_method_malignancy_support"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def load_supported_cells(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"CNA cell-evidence table missing or empty: {path}")
    table = pd.read_csv(path, sep="\t", dtype=str)
    missing = REQUIRED_EVIDENCE_COLUMNS - set(table.columns)
    if missing:
        raise ValueError(f"CNA cell-evidence fields missing: {sorted(missing)}")
    if table["cell_id"].isna().any() or table["cell_id"].duplicated().any():
        raise ValueError("CNA cell-evidence IDs are missing or duplicated")
    table["panel_order"] = pd.to_numeric(table["panel_order"], errors="raise").astype(int)
    references = table[table["cna_input_role"].eq("known_normal_reference")].copy()
    if references.empty:
        raise ValueError("CNA panel evidence has no known-normal reference cells")
    references["unexpected_malignancy_support"] = (
        references["copykat.pred"].eq("aneuploid") | references["scevan_class"].eq("tumor")
    )
    reference_qc = (
        references.groupby("panel_order", observed=True, sort=True)
        .agg(
            unit_normal_reference_cells=("cell_id", "size"),
            unit_normal_reference_unexpected_support=("unexpected_malignancy_support", "sum"),
        )
        .reset_index()
    )
    reference_qc["unit_normal_reference_unexpected_support_fraction"] = (
        reference_qc["unit_normal_reference_unexpected_support"]
        / reference_qc["unit_normal_reference_cells"]
    )
    selected = table[table["cna_selection_class"].eq(SELECTED_CLASS)].copy()
    if selected.empty:
        raise ValueError("no two-method CNA-supported author Cancer cells are available")
    invalid = selected[
        ~selected["cna_input_role"].eq("author_cancer")
        | ~selected["copykat.pred"].eq("aneuploid")
        | ~selected["scevan_class"].eq("tumor")
    ]
    if not invalid.empty:
        raise ValueError("CNA selection class disagrees with role or method calls")
    selected = selected.merge(reference_qc, on="panel_order", how="left", validate="many_to_one")
    if selected["unit_normal_reference_cells"].isna().any():
        raise ValueError("a CNA-supported unit has no matched normal-reference QC evidence")
    return selected


def balance_supported_cells(
    selected: pd.DataFrame,
    min_cells_per_unit: int,
    max_cells_per_unit: int,
    seed: int,
    required_states: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if min_cells_per_unit < 1 or max_cells_per_unit < min_cells_per_unit:
        raise ValueError("cNMF per-unit limits are invalid")
    unit_columns = [
        "panel_order",
        "panel_dataset",
        "panel_sample_id",
        "panel_patient_id",
        "panel_sample_type",
        "panel_tissue",
    ]
    counts = selected.groupby(unit_columns, observed=True, sort=True).size().rename("supported_cells").reset_index()
    eligible = counts[counts["supported_cells"] >= min_cells_per_unit].copy()
    if eligible.empty:
        raise ValueError("no CNA panel unit meets the minimum two-method-supported cell count")
    if required_states:
        missing_states = sorted(set(required_states) - set(eligible["panel_sample_type"].astype(str)))
        if missing_states:
            raise ValueError(f"required CRC states lack eligible CNA-supported units: {missing_states}")
    eligible_orders = set(eligible["panel_order"].astype(int))
    rng = np.random.default_rng(seed)
    chosen: list[pd.DataFrame] = []
    for order, group in selected[selected["panel_order"].isin(eligible_orders)].groupby(
        "panel_order", observed=True, sort=True
    ):
        group = group.sort_values("cell_id")
        if len(group) > max_cells_per_unit:
            positions = np.sort(rng.choice(len(group), max_cells_per_unit, replace=False))
            group = group.iloc[positions]
        chosen.append(group)
    balanced = pd.concat(chosen, ignore_index=True).sort_values(["panel_order", "cell_id"])
    selection = (
        balanced.groupby(unit_columns, observed=True, sort=True)
        .size()
        .rename("selected_cells")
        .reset_index()
        .merge(counts, on=unit_columns, how="left", validate="one_to_one")
    )
    qc_columns = [
        "unit_normal_reference_cells",
        "unit_normal_reference_unexpected_support",
        "unit_normal_reference_unexpected_support_fraction",
    ]
    if set(qc_columns).issubset(balanced.columns):
        unit_qc = balanced.groupby("panel_order", observed=True, sort=True)[qc_columns].first().reset_index()
        selection = selection.merge(unit_qc, on="panel_order", how="left", validate="one_to_one")
    return balanced, selection


def prepare_cnmf_from_cna(
    h5ad: Path,
    evidence_path: Path,
    output_dir: Path,
    min_cells_per_unit: int = 30,
    max_cells_per_unit: int = 150,
    required_states: list[str] | None = None,
    sample_key: str = "sample_id",
    patient_key: str = "donor_id",
    dataset_key: str = "dataset",
    state_key: str = "sample_type",
    tissue_key: str = "tissue",
    symbol_key: str = "GeneSymbol",
    seed: int = 20260810,
) -> dict[str, object]:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD missing or empty: {h5ad}")
    selected = load_supported_cells(evidence_path)
    balanced, unit_table = balance_supported_cells(
        selected,
        min_cells_per_unit,
        max_cells_per_unit,
        seed,
        required_states,
    )

    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        required_obs = {sample_key, patient_key, dataset_key, state_key, tissue_key}
        missing = sorted(required_obs - set(adata.obs.columns))
        if missing:
            raise ValueError(f"H5AD obs missing CNA-derived cNMF fields: {missing}")
        if adata.raw is None:
            raise ValueError("H5AD raw is missing; cNMF input must not use normalized X")
        if not adata.obs_names.is_unique:
            raise ValueError("H5AD cell IDs are not unique")
        rows = adata.obs_names.get_indexer(balanced["cell_id"].astype(str))
        if (rows < 0).any():
            missing_ids = balanced.loc[rows < 0, "cell_id"].head(3).tolist()
            raise ValueError(f"CNA-supported cells are absent from H5AD: {missing_ids}")
        # Backed sparse matrices require monotonic row access on some anndata/
        # h5py combinations.  Reorder metadata by source row before slicing;
        # this changes no selection and makes the extraction deterministic.
        source_order = np.argsort(rows, kind="stable")
        rows = rows[source_order]
        balanced = balanced.iloc[source_order].reset_index(drop=True)
        source_obs = adata.obs.iloc[rows][[sample_key, patient_key, dataset_key, state_key, tissue_key]].astype(str)
        checks = {
            sample_key: balanced["panel_sample_id"].astype(str).to_numpy(),
            dataset_key: balanced["panel_dataset"].astype(str).to_numpy(),
            state_key: balanced["panel_sample_type"].astype(str).to_numpy(),
            tissue_key: balanced["panel_tissue"].astype(str).to_numpy(),
        }
        for column, expected in checks.items():
            observed = source_obs[column].astype(str).to_numpy()
            if not np.array_equal(observed, expected):
                raise ValueError(f"CNA panel identity and H5AD {column} disagree")
        observed_patients = source_obs[patient_key].map(canonical_patient).to_numpy()
        expected_patients = balanced["panel_patient_id"].map(canonical_patient).to_numpy()
        if not np.array_equal(observed_patients, expected_patients):
            raise ValueError("CNA panel identity and H5AD patient IDs disagree")

        raw_matrix = adata.raw.X[rows, :]
        if not sparse.issparse(raw_matrix):
            raw_matrix = sparse.csr_matrix(np.asarray(raw_matrix))
        raw_matrix = raw_matrix.tocsr()
        feature_map = unique_gene_features(adata.raw.var.copy(), symbol_key)
        feature_indices = feature_map["feature_index"].to_numpy(dtype=np.int64)
        counts = raw_matrix[:, feature_indices].tocsr()
        if counts.nnz and (not np.isfinite(counts.data).all() or (counts.data < 0).any()):
            raise ValueError("raw/X contains non-finite or negative values")
        if counts.nnz and not np.allclose(counts.data, np.rint(counts.data), atol=1e-6):
            raise ValueError("raw/X contains non-integer values; refusing cNMF input")
        counts.data = np.rint(counts.data).astype(np.float32)
        counts.eliminate_zeros()
        totals = np.asarray(counts.sum(axis=1)).ravel()
        if (totals <= 0).any():
            raise ValueError("CNA-supported cNMF input contains zero-count cells")

        obs = balanced.copy()
        obs.index = obs["cell_id"].astype(str)
        obs.index.name = None
        obs["analysis_patient_id"] = obs["panel_patient_id"].map(canonical_patient)
        obs["dataset"] = obs["panel_dataset"].astype(str)
        obs["sample_id"] = obs["panel_sample_id"].astype(str)
        obs["sample_type"] = obs["panel_sample_type"].astype(str)
        obs["tissue"] = obs["panel_tissue"].astype(str)
        var = pd.DataFrame(index=feature_map["gene"].astype(str).to_numpy())
        var.index.name = "gene"
        cnmf_input = ad.AnnData(X=counts, obs=obs, var=var)
    finally:
        adata.file.close()

    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix with existing CNA-derived cNMF input: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = output_dir / "crc_two_method_cna_supported_counts.h5ad"
    unit_path = output_dir / "P0_cnmf_cna_selected_units.tsv"
    feature_path = output_dir / "P0_cnmf_cna_gene_feature_map.tsv"
    cnmf_input.write_h5ad(input_path, compression="gzip")
    unit_table.to_csv(unit_path, sep="\t", index=False, lineterminator="\n")
    feature_map.to_csv(feature_path, sep="\t", index=False, lineterminator="\n")

    units_with_reference_warning = int(
        unit_table["unit_normal_reference_unexpected_support"].gt(0).sum()
    )
    receipt: dict[str, object] = {
        "status": "diagnostic_two_method_cna_supported_cnmf_input",
        "scientific_gate": (
            "not_passed_unexpected_normal_reference_support"
            if units_with_reference_warning
            else "not_passed_inferred_CNA_without_DNA_truth"
        ),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "h5ad": str(h5ad.resolve()),
        "cna_evidence_sha256": sha256_file(evidence_path),
        "raw_matrix": "raw/X",
        "selection_class": SELECTED_CLASS,
        "cells": int(cnmf_input.n_obs),
        "genes": int(cnmf_input.n_vars),
        "units": int(len(unit_table)),
        "datasets": int(cnmf_input.obs["dataset"].nunique()),
        "patients": int(cnmf_input.obs["analysis_patient_id"].nunique()),
        "states": sorted(cnmf_input.obs["sample_type"].astype(str).unique().tolist()),
        "units_with_unexpected_normal_reference_support": units_with_reference_warning,
        "maximum_unit_normal_reference_unexpected_support_fraction": float(
            unit_table["unit_normal_reference_unexpected_support_fraction"].max()
        ),
        "min_cells_per_unit": min_cells_per_unit,
        "max_cells_per_unit": max_cells_per_unit,
        "seed": seed,
        "selection_rule": (
            "retain author Cancer candidates supported as aneuploid by CopyKAT and tumor by SCEVAN; "
            "exclude units below the frozen minimum; sample cells within unit without reading genes or outcomes"
        ),
        "claim_boundary": (
            "The input is more conservative than the author Cancer label but remains an inferred, sampled "
            "CNA-supported candidate set rather than DNA-validated malignancy truth. It represents selected "
            "panel units only; sample-level unexpected malignant calls among predefined normal references "
            "are retained in obs rather than hidden. It is suitable for diagnostic cNMF sensitivity analysis, "
            "not yet atlas-wide program freezing."
        ),
        "outputs": {},
    }
    for path in (input_path, unit_path, feature_path):
        receipt["outputs"][path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    receipt_path = output_dir / "P0_cnmf_cna_input_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--cna-evidence", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--min-cells-per-unit", type=int, default=30)
    parser.add_argument("--max-cells-per-unit", type=int, default=150)
    parser.add_argument("--required-state", action="append", dest="required_states")
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--patient-key", default="donor_id")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--state-key", default="sample_type")
    parser.add_argument("--tissue-key", default="tissue")
    parser.add_argument("--symbol-key", default="GeneSymbol")
    parser.add_argument("--seed", type=int, default=20260810)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = prepare_cnmf_from_cna(
        args.h5ad,
        args.cna_evidence,
        args.output_dir,
        args.min_cells_per_unit,
        args.max_cells_per_unit,
        args.required_states,
        args.sample_key,
        args.patient_key,
        args.dataset_key,
        args.state_key,
        args.tissue_key,
        args.symbol_key,
        args.seed,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
