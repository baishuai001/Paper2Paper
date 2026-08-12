#!/usr/bin/env python3
"""Create patient- and dataset-disjoint, state-stratified cNMF holdout folds."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


REQUIRED_OBS = {
    "panel_order",
    "dataset",
    "analysis_patient_id",
    "sample_id",
    "sample_type",
    "cna_selection_class",
}
EXPECTED_SELECTION_CLASS = "author_cancer_two_method_malignancy_support"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def assign_folds(obs: pd.DataFrame, seed: int) -> pd.DataFrame:
    missing = REQUIRED_OBS - set(obs.columns)
    if missing:
        raise ValueError(f"cNMF holdout metadata fields missing: {sorted(missing)}")
    if set(obs["cna_selection_class"].astype(str)) != {EXPECTED_SELECTION_CLASS}:
        raise ValueError("holdout input is not restricted to the two-method CNA-supported class")
    unit_columns = ["panel_order", "dataset", "analysis_patient_id", "sample_id", "sample_type"]
    units = (
        obs.groupby(unit_columns, observed=True, sort=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    if units["panel_order"].duplicated().any():
        raise ValueError("one panel_order maps to multiple holdout identities")
    if units["dataset"].duplicated().any() or units["analysis_patient_id"].duplicated().any():
        raise ValueError("holdout units are not dataset- and patient-disjoint before folding")
    state_counts = units["sample_type"].value_counts()
    insufficient = state_counts[state_counts < 2]
    if not insufficient.empty:
        raise ValueError(f"states need at least two independent units for holdout: {insufficient.to_dict()}")

    rng = np.random.default_rng(seed)
    assignments: list[dict[str, object]] = []
    for state, group in units.groupby("sample_type", observed=True, sort=True):
        group = group.sort_values("panel_order").reset_index(drop=True)
        positions = rng.permutation(len(group))
        randomized = group.iloc[positions].reset_index(drop=True)
        randomized["holdout_fold"] = ["A" if index % 2 == 0 else "B" for index in range(len(randomized))]
        assignments.extend(randomized.to_dict(orient="records"))
    manifest = pd.DataFrame(assignments).sort_values(["sample_type", "holdout_fold", "panel_order"])
    for state, group in manifest.groupby("sample_type", observed=True):
        if set(group["holdout_fold"]) != {"A", "B"}:
            raise ValueError(f"state {state!r} is absent from one holdout fold")
    return manifest


def split_holdouts(input_h5ad: Path, output_dir: Path, seed: int = 20260810) -> dict[str, object]:
    if not input_h5ad.is_file() or input_h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"cNMF input missing or empty: {input_h5ad}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix with existing holdout output: {output_dir}")
    data = ad.read_h5ad(input_h5ad)
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("cNMF holdout input cell or gene identifiers are not unique")
    manifest = assign_folds(data.obs, seed)
    fold_by_order = manifest.set_index("panel_order")["holdout_fold"].to_dict()
    cell_orders = pd.to_numeric(data.obs["panel_order"], errors="raise").astype(int)
    cell_folds = cell_orders.map(fold_by_order)
    if cell_folds.isna().any():
        raise ValueError("some cNMF cells were not assigned to a holdout fold")

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "P0_cnmf_holdout_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False, lineterminator="\n")
    fold_paths: dict[str, Path] = {}
    for fold in ("A", "B"):
        subset = data[cell_folds.eq(fold).to_numpy()].copy()
        if subset.n_obs == 0:
            raise ValueError(f"cNMF holdout fold {fold} is empty")
        fold_paths[fold] = output_dir / f"crc_two_method_cna_supported_holdout_{fold}.h5ad"
        subset.write_h5ad(fold_paths[fold], compression="gzip")

    a_patients = set(data.obs.loc[cell_folds.eq("A"), "analysis_patient_id"].astype(str))
    b_patients = set(data.obs.loc[cell_folds.eq("B"), "analysis_patient_id"].astype(str))
    a_datasets = set(data.obs.loc[cell_folds.eq("A"), "dataset"].astype(str))
    b_datasets = set(data.obs.loc[cell_folds.eq("B"), "dataset"].astype(str))
    if a_patients & b_patients or a_datasets & b_datasets:
        raise ValueError("holdout folds leak patients or datasets")

    receipt: dict[str, object] = {
        "status": "diagnostic_patient_dataset_disjoint_holdouts_prepared",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "input_sha256": sha256_file(input_h5ad),
        "seed": seed,
        "assignment_rule": (
            "within each state, shuffle independent units using the frozen seed and alternate them "
            "between folds; do not read genes, program usage or outcomes"
        ),
        "folds": {},
        "claim_boundary": (
            "The folds test replication across disjoint selected patients and datasets. They are not "
            "randomized disease cohorts, and state remains partly confounded with source datasets."
        ),
        "outputs": {},
    }
    for fold, path in fold_paths.items():
        subset_obs = data.obs.loc[cell_folds.eq(fold)]
        receipt["folds"][fold] = {
            "cells": int(len(subset_obs)),
            "patients": int(subset_obs["analysis_patient_id"].astype(str).nunique()),
            "datasets": int(subset_obs["dataset"].astype(str).nunique()),
            "states": sorted(subset_obs["sample_type"].astype(str).unique().tolist()),
        }
        receipt["outputs"][path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    receipt["outputs"][manifest_path.name] = {
        "bytes": manifest_path.stat().st_size,
        "sha256": sha256_file(manifest_path),
    }
    receipt_path = output_dir / "P0_cnmf_holdout_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260810)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = split_holdouts(args.input_h5ad, args.output_dir, args.seed)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
