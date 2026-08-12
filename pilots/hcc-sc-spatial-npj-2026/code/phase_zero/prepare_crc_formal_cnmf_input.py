#!/usr/bin/env python3
"""Build a frozen patient-state-balanced Liu discovery input for formal cNMF."""

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

from formal_cnmf_design import (
    ONE_METHOD_CLASS,
    SELECTED_CLASS,
    balance_patient_state_cells,
    load_evidence_tables,
    select_cohort_cells,
)
from join_crc_h5ad_table_s1 import canonical_patient
from prepare_crc_cna_inputs import unique_gene_features


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def prepare_formal_input(
    h5ad: Path,
    evidence_paths: list[Path],
    output_dir: Path,
    h5ad_sha256: str,
    cohort_dataset: str,
    analysis_role: str,
    expected_patients: int,
    expected_samples: int,
    expected_dual_cells: int,
    max_cells_per_patient_state: int,
    seed: int,
    sample_key: str = "sample_id",
    patient_key: str = "donor_id",
    dataset_key: str = "dataset",
    state_key: str = "sample_type",
    tissue_key: str = "tissue",
    symbol_key: str = "GeneSymbol",
    selection_mode: str = "dual_only",
    excluded_sample_ids: tuple[str, ...] = (),
    require_paired_states: bool = True,
    expected_patient_state_units: int | None = None,
    input_filename: str = "crc_liu_dual_cna_patient_state_balanced_counts.h5ad",
) -> dict[str, object]:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD missing or empty: {h5ad}")
    normalized_h5ad_sha256 = h5ad_sha256.strip().upper()
    if len(normalized_h5ad_sha256) != 64 or any(
        character not in "0123456789ABCDEF" for character in normalized_h5ad_sha256
    ):
        raise ValueError("h5ad_sha256 must be one hexadecimal SHA256")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix with existing formal cNMF input: {output_dir}")

    if selection_mode == "dual_only":
        selected_classes = (SELECTED_CLASS,)
    elif selection_mode == "dual_plus_one_method":
        selected_classes = (SELECTED_CLASS, ONE_METHOD_CLASS)
    else:
        raise ValueError(f"unsupported cNMF sensitivity selection mode: {selection_mode}")
    if Path(input_filename).name != input_filename or not input_filename.endswith(".h5ad"):
        raise ValueError("input_filename must be one local .h5ad filename")

    evidence = load_evidence_tables(evidence_paths)
    selected, sample_qc = select_cohort_cells(
        evidence,
        cohort_dataset,
        expected_patients=expected_patients,
        expected_samples=expected_samples,
        expected_dual_cells=expected_dual_cells,
        selected_classes=selected_classes,
        excluded_sample_ids=excluded_sample_ids,
        require_paired_states=require_paired_states,
    )
    balanced, patient_state, sample_composition = balance_patient_state_cells(
        selected, max_cells_per_patient_state, seed
    )

    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        required_obs = {sample_key, patient_key, dataset_key, state_key, tissue_key}
        missing = sorted(required_obs - set(adata.obs.columns))
        if missing:
            raise ValueError(f"H5AD obs missing formal cNMF fields: {missing}")
        if adata.raw is None:
            raise ValueError("H5AD raw is missing; formal cNMF must use raw integer counts")
        if not adata.obs_names.is_unique:
            raise ValueError("H5AD cell IDs are not unique")
        rows = adata.obs_names.get_indexer(balanced["cell_id"].astype(str))
        if (rows < 0).any():
            missing_ids = balanced.loc[rows < 0, "cell_id"].head(3).tolist()
            raise ValueError(f"selected cells are absent from H5AD: {missing_ids}")
        source_order = np.argsort(rows, kind="stable")
        rows = rows[source_order]
        balanced = balanced.iloc[source_order].reset_index(drop=True)
        source_obs = adata.obs.iloc[rows][
            [sample_key, patient_key, dataset_key, state_key, tissue_key]
        ].astype(str)
        checks = {
            sample_key: balanced["panel_sample_id"].astype(str).to_numpy(),
            dataset_key: balanced["panel_dataset"].astype(str).to_numpy(),
            state_key: balanced["panel_sample_type"].astype(str).to_numpy(),
            tissue_key: balanced["panel_tissue"].astype(str).to_numpy(),
        }
        for column, expected in checks.items():
            if not np.array_equal(source_obs[column].astype(str).to_numpy(), expected):
                raise ValueError(f"CNA evidence and H5AD {column} disagree")
        observed_patients = source_obs[patient_key].map(canonical_patient).to_numpy()
        expected_patient_ids = balanced["panel_patient_id"].map(canonical_patient).to_numpy()
        if not np.array_equal(observed_patients, expected_patient_ids):
            raise ValueError("CNA evidence and H5AD patient IDs disagree")

        raw_matrix = adata.raw.X[rows, :]
        if not sparse.issparse(raw_matrix):
            raw_matrix = sparse.csr_matrix(np.asarray(raw_matrix))
        feature_map = unique_gene_features(adata.raw.var.copy(), symbol_key)
        counts = raw_matrix.tocsr()[:, feature_map["feature_index"].to_numpy(dtype=np.int64)].tocsr()
        if counts.nnz and (not np.isfinite(counts.data).all() or (counts.data < 0).any()):
            raise ValueError("raw/X contains non-finite or negative values")
        if counts.nnz and not np.allclose(counts.data, np.rint(counts.data), atol=1e-6):
            raise ValueError("raw/X contains non-integer values")
        counts.data = np.rint(counts.data).astype(np.float32)
        counts.eliminate_zeros()
        if (np.asarray(counts.sum(axis=1)).ravel() <= 0).any():
            raise ValueError("formal cNMF input contains zero-count cells")

        obs = balanced.copy()
        obs.index = obs["cell_id"].astype(str)
        obs.index.name = None
        obs["analysis_patient_id"] = obs["panel_patient_id"].map(canonical_patient)
        obs["dataset"] = obs["panel_dataset"].astype(str)
        obs["sample_id"] = obs["panel_sample_id"].astype(str)
        obs["sample_type"] = obs["panel_sample_type"].astype(str)
        obs["tissue"] = obs["panel_tissue"].astype(str)
        obs["formal_analysis_role"] = analysis_role
        obs["balance_unit"] = "patient_state"
        obs["cell_sampling_seed"] = seed
        var = pd.DataFrame(index=feature_map["gene"].astype(str).to_numpy())
        var.index.name = "gene"
        cnmf_input = ad.AnnData(X=counts, obs=obs, var=var)
    finally:
        adata.file.close()

    observed_patient_states = cnmf_input.obs[
        ["analysis_patient_id", "sample_type"]
    ].drop_duplicates()
    if len(cnmf_input.obs["analysis_patient_id"].unique()) != expected_patients:
        raise ValueError("patient loss occurred during H5AD extraction")
    expected_units = (
        expected_patients * 2 if expected_patient_state_units is None else expected_patient_state_units
    )
    if len(observed_patient_states) != expected_units:
        raise ValueError(
            f"expected {expected_units} patient-state units after H5AD extraction, "
            f"observed {len(observed_patient_states)}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    input_path = output_dir / input_filename
    patient_state_path = output_dir / "P0_formal_cnmf_patient_state_units.tsv"
    sample_path = output_dir / "P0_formal_cnmf_sample_composition.tsv"
    selected_path = output_dir / "P0_formal_cnmf_selected_cells.tsv"
    reference_path = output_dir / "P0_formal_cnmf_reference_qc.tsv"
    feature_path = output_dir / "P0_formal_cnmf_gene_feature_map.tsv"
    cnmf_input.write_h5ad(input_path, compression="gzip")
    patient_state.to_csv(patient_state_path, sep="\t", index=False, lineterminator="\n")
    sample_composition.to_csv(sample_path, sep="\t", index=False, lineterminator="\n")
    balanced.to_csv(selected_path, sep="\t", index=False, lineterminator="\n")
    sample_qc.to_csv(reference_path, sep="\t", index=False, lineterminator="\n")
    feature_map.to_csv(feature_path, sep="\t", index=False, lineterminator="\n")

    high_reference_unresolved = sample_qc.loc[
        sample_qc["normal_reference_method_unresolved_fraction"].gt(0.5),
        "panel_sample_id",
    ].astype(str).tolist()
    retained_sample_ids = set(selected["panel_sample_id"].astype(str))
    retained_high_reference_unresolved = sorted(
        set(high_reference_unresolved) & retained_sample_ids
    )
    sensitivity = analysis_role != "liu_primary_discovery"
    if sensitivity:
        status = "formal_liu_sensitivity_input_frozen_pending_cnmf"
        claim_boundary = (
            "This is one frozen Liu sensitivity input at candidate K=9-11 only. It may test "
            "whether the primary program solution changes, but cannot freeze K by itself, name "
            "programs, support a primary-to-metastasis contrast when patient states are incomplete, "
            "or replace Che external replication."
        )
    else:
        status = "formal_liu_discovery_input_frozen_pending_cnmf"
        claim_boundary = (
            "This is a patient-state-balanced expression-count input from two-method CNA-supported "
            "author Cancer candidates, not DNA-validated malignancy. It is frozen for Liu program "
            "discovery; Che remains external replication and no K or program name is yet accepted."
        )
    receipt: dict[str, object] = {
        "status": status,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "design_helper_sha256": sha256_file(Path(__file__).with_name("formal_cnmf_design.py")),
        "h5ad": str(h5ad.resolve()),
        "h5ad_sha256": normalized_h5ad_sha256,
        "h5ad_sha256_verification": (
            "verified once by the execution supervisor immediately before both deterministic inputs"
        ),
        "evidence": [
            {"path": str(path.resolve()), "sha256": sha256_file(path)} for path in evidence_paths
        ],
        "cohort_dataset": cohort_dataset,
        "analysis_role": analysis_role,
        "selection_classes": list(selected_classes),
        "selection_mode": selection_mode,
        "available_selected_cells": int(len(selected)),
        "available_dual_supported_cells": int(
            selected["cna_selection_class"].eq(SELECTED_CLASS).sum()
        ),
        "available_one_method_supported_cells": int(
            selected["cna_selection_class"].eq(ONE_METHOD_CLASS).sum()
        ),
        "selected_cells": int(cnmf_input.n_obs),
        "genes": int(cnmf_input.n_vars),
        "patients": int(cnmf_input.obs["analysis_patient_id"].nunique()),
        "samples": int(cnmf_input.obs["sample_id"].nunique()),
        "states": sorted(cnmf_input.obs["sample_type"].astype(str).unique().tolist()),
        "patient_state_units": int(len(observed_patient_states)),
        "excluded_sample_ids": list(excluded_sample_ids),
        "paired_states_required": require_paired_states,
        "balance_unit": "patient_state",
        "max_cells_per_patient_state": max_cells_per_patient_state,
        "cell_sampling_seed": seed,
        "high_reference_unresolved_samples_gt_0_5": high_reference_unresolved,
        "retained_high_reference_unresolved_samples_gt_0_5": (
            retained_high_reference_unresolved
        ),
        "selection_rule": (
            f"within the frozen Liu cohort, retain {selection_mode} author Cancer candidates; "
            f"exclude only the predeclared samples {list(excluded_sample_ids)}; pool replicate "
            "samples inside each patient-state and cap that patient-state without reading expression "
            "values or downstream outcomes"
        ),
        "claim_boundary": claim_boundary,
        "outputs": {},
    }
    for path in (
        input_path,
        patient_state_path,
        sample_path,
        selected_path,
        reference_path,
        feature_path,
    ):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    receipt_path = output_dir / "P0_formal_cnmf_input_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--cna-evidence", required=True, action="append", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--h5ad-sha256", required=True)
    parser.add_argument("--cohort-dataset", required=True)
    parser.add_argument("--analysis-role", default="liu_primary_discovery")
    parser.add_argument("--expected-patients", required=True, type=int)
    parser.add_argument("--expected-samples", required=True, type=int)
    parser.add_argument("--expected-dual-cells", required=True, type=int)
    parser.add_argument("--max-cells-per-patient-state", type=int, default=500)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--patient-key", default="donor_id")
    parser.add_argument("--dataset-key", default="dataset")
    parser.add_argument("--state-key", default="sample_type")
    parser.add_argument("--tissue-key", default="tissue")
    parser.add_argument("--symbol-key", default="GeneSymbol")
    parser.add_argument(
        "--selection-mode",
        choices=["dual_only", "dual_plus_one_method"],
        default="dual_only",
    )
    parser.add_argument("--exclude-sample-id", action="append", default=[])
    parser.add_argument("--allow-incomplete-patient-states", action="store_true")
    parser.add_argument("--expected-patient-state-units", type=int)
    parser.add_argument(
        "--input-filename",
        default="crc_liu_dual_cna_patient_state_balanced_counts.h5ad",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = prepare_formal_input(
        args.h5ad,
        args.cna_evidence,
        args.output_dir,
        args.h5ad_sha256,
        args.cohort_dataset,
        args.analysis_role,
        args.expected_patients,
        args.expected_samples,
        args.expected_dual_cells,
        args.max_cells_per_patient_state,
        args.seed,
        args.sample_key,
        args.patient_key,
        args.dataset_key,
        args.state_key,
        args.tissue_key,
        args.symbol_key,
        args.selection_mode,
        tuple(args.exclude_sample_id),
        not args.allow_incomplete_patient_states,
        args.expected_patient_state_units,
        args.input_filename,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
