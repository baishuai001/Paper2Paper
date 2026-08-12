#!/usr/bin/env python3
"""Create paired-patient-disjoint Liu cNMF validation folds."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad

from formal_cnmf_design import assign_patient_folds


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def split_patient_holdouts(input_h5ad: Path, output_dir: Path, seed: int) -> dict[str, object]:
    if not input_h5ad.is_file() or input_h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"formal cNMF input missing or empty: {input_h5ad}")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix with existing holdout output: {output_dir}")
    data = ad.read_h5ad(input_h5ad)
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("formal cNMF holdout input IDs are not unique")
    if "formal_analysis_role" not in data.obs.columns or set(
        data.obs["formal_analysis_role"].astype(str)
    ) != {"liu_primary_discovery"}:
        raise ValueError("patient holdouts require the frozen Liu primary-discovery input")
    manifest = assign_patient_folds(data.obs, seed)
    fold_by_patient = manifest.set_index("analysis_patient_id")["holdout_fold"].to_dict()
    cell_folds = data.obs["analysis_patient_id"].astype(str).map(fold_by_patient)
    if cell_folds.isna().any():
        raise ValueError("some cNMF cells were not assigned to a patient holdout")

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "P0_formal_cnmf_patient_holdout_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False, lineterminator="\n")
    fold_paths: dict[str, Path] = {}
    fold_details: dict[str, object] = {}
    patient_sets: dict[str, set[str]] = {}
    for fold in ("A", "B"):
        mask = cell_folds.eq(fold).to_numpy()
        subset = data[mask].copy()
        if subset.n_obs == 0:
            raise ValueError(f"formal cNMF patient holdout {fold} is empty")
        patients = set(subset.obs["analysis_patient_id"].astype(str))
        patient_sets[fold] = patients
        if set(subset.obs["sample_type"].astype(str)) != {"tumor", "metastasis"}:
            raise ValueError(f"holdout {fold} lacks a required disease state")
        fold_paths[fold] = output_dir / f"crc_liu_patient_holdout_{fold}_counts.h5ad"
        subset.write_h5ad(fold_paths[fold], compression="gzip")
        fold_details[fold] = {
            "cells": int(subset.n_obs),
            "patients": int(len(patients)),
            "states": sorted(subset.obs["sample_type"].astype(str).unique().tolist()),
        }
    if patient_sets["A"] & patient_sets["B"]:
        raise ValueError("patient IDs leak across formal cNMF holdouts")

    receipt: dict[str, object] = {
        "status": "formal_liu_patient_holdouts_frozen_pending_cnmf",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "design_helper_sha256": sha256_file(Path(__file__).with_name("formal_cnmf_design.py")),
        "input": str(input_h5ad.resolve()),
        "input_sha256": sha256_file(input_h5ad),
        "seed": seed,
        "assignment_rule": (
            "shuffle paired Liu patients with the frozen seed and alternate whole patients between "
            "folds; primary and metastasis cells from one patient never cross folds"
        ),
        "folds": fold_details,
        "claim_boundary": (
            "These folds are patient-disjoint but share the Liu dataset. They test internal patient "
            "replication only; Che is reserved for external cohort replication after K is frozen."
        ),
        "outputs": {},
    }
    for path in (*fold_paths.values(), manifest_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    receipt_path = output_dir / "P0_formal_cnmf_patient_holdout_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = split_patient_holdouts(args.input_h5ad, args.output_dir, args.seed)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
