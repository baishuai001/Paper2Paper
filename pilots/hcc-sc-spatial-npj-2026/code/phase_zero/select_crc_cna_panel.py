#!/usr/bin/env python3
"""Select an outcome-blind, dataset-balanced CRC panel for CNA qualification."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd

from join_crc_h5ad_table_s1 import canonical_patient
from prepare_crc_cna_inputs import DEFAULT_REFERENCE_LABELS


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def verify_upstream_h5ad_receipt(h5ad: Path, receipt_path: Path) -> dict[str, object]:
    if not receipt_path.is_file() or receipt_path.stat().st_size == 0:
        raise FileNotFoundError(f"upstream H5AD receipt missing or empty: {receipt_path}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("status") != "passed":
        raise ValueError("upstream H5AD receipt is not passed")
    try:
        recorded_path = Path(str(receipt["h5ad"])).resolve()
        recorded_bytes = int(receipt["bytes"])
        recorded_sha256 = str(receipt["sha256"]).upper()
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("upstream H5AD receipt lacks path, bytes, or SHA256") from exc
    if recorded_path != h5ad.resolve():
        raise ValueError("upstream H5AD receipt points to a different file")
    if h5ad.stat().st_size != recorded_bytes:
        raise ValueError("H5AD byte size changed since the upstream receipt")
    if not re.fullmatch(r"[0-9A-F]{64}", recorded_sha256):
        raise ValueError("upstream H5AD receipt contains an invalid SHA256")
    return {
        "receipt": str(receipt_path.resolve()),
        "bytes": recorded_bytes,
        "sha256": recorded_sha256,
        "verification": (
            "content SHA256 inherited from passed P0.0 receipt; downstream step rechecked "
            "resolved path and byte size without rehashing the 31 GB object"
        ),
    }


def build_eligible_table(
    h5ad: Path,
    min_cancer: int,
    min_reference: int,
    sample_key: str = "sample_id",
    patient_key: str = "donor_id",
    dataset_key: str = "dataset",
    state_key: str = "sample_type",
    tissue_key: str = "tissue",
    label_key: str = "cell_type_coarse_crc_atlas",
) -> pd.DataFrame:
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD missing or empty: {h5ad}")
    data = ad.read_h5ad(h5ad, backed="r")
    try:
        fields = [sample_key, patient_key, dataset_key, state_key, tissue_key, label_key]
        missing = sorted(set(fields) - set(data.obs.columns))
        if missing:
            raise ValueError(f"H5AD obs missing CNA-panel fields: {missing}")
        obs = data.obs[fields].astype(str).copy()
    finally:
        data.file.close()
    obs["analysis_patient_id"] = obs[patient_key].map(canonical_patient)
    obs["is_cancer"] = obs[label_key].eq("Cancer cell")
    obs["is_reference"] = obs[label_key].isin(DEFAULT_REFERENCE_LABELS)
    obs["is_other_epithelial"] = obs[label_key].eq("Epithelial cell")
    grouped = (
        obs.groupby(
            [dataset_key, sample_key, "analysis_patient_id", state_key, tissue_key],
            observed=True,
            sort=True,
        )
        .agg(
            cells=(label_key, "size"),
            cancer_cells=("is_cancer", "sum"),
            reference_cells=("is_reference", "sum"),
            other_epithelial_cells=("is_other_epithelial", "sum"),
        )
        .reset_index()
        .rename(
            columns={
                dataset_key: "dataset",
                sample_key: "sample_id",
                state_key: "sample_type",
                tissue_key: "tissue",
            }
        )
    )
    grouped["eligible"] = (
        (grouped["cancer_cells"] >= min_cancer)
        & (grouped["reference_cells"] >= min_reference)
    )
    grouped["limiting_role_cells"] = grouped[["cancer_cells", "reference_cells"]].min(axis=1)
    grouped["analysis_unit_id"] = grouped["sample_id"] + "::" + grouped["analysis_patient_id"]
    if grouped["analysis_unit_id"].duplicated().any():
        raise ValueError("sample-patient CNA analysis units are not unique")
    return grouped


def select_panel(
    eligible_table: pd.DataFrame,
    states: list[str],
    max_datasets_per_state: int,
) -> pd.DataFrame:
    if max_datasets_per_state < 1:
        raise ValueError("max_datasets_per_state must be positive")
    required = {
        "dataset",
        "sample_id",
        "analysis_patient_id",
        "sample_type",
        "tissue",
        "cancer_cells",
        "reference_cells",
        "limiting_role_cells",
        "eligible",
    }
    missing = required - set(eligible_table.columns)
    if missing:
        raise ValueError(f"eligible table fields missing: {sorted(missing)}")
    selected_rows: list[pd.Series] = []
    for state in states:
        candidates = eligible_table[
            eligible_table["eligible"].astype(bool) & eligible_table["sample_type"].eq(state)
        ].copy()
        if candidates.empty:
            raise ValueError(f"no CNA-eligible samples for requested state: {state}")
        dataset_rank = (
            candidates.groupby("dataset", observed=True)
            .agg(
                eligible_samples=("sample_id", "nunique"),
                eligible_patients=("analysis_patient_id", "nunique"),
                limiting_role_cells=("limiting_role_cells", "sum"),
            )
            .reset_index()
            .sort_values(
                ["eligible_patients", "eligible_samples", "limiting_role_cells", "dataset"],
                ascending=[False, False, False, True],
            )
        )
        chosen_datasets = dataset_rank.head(max_datasets_per_state)["dataset"].tolist()
        for dataset_rank_index, dataset in enumerate(chosen_datasets, start=1):
            subset = candidates[candidates["dataset"].eq(dataset)].sort_values(
                ["limiting_role_cells", "cancer_cells", "reference_cells", "sample_id", "analysis_patient_id"],
                ascending=[False, False, False, True, True],
            )
            row = subset.iloc[0].copy()
            row["state_dataset_rank"] = dataset_rank_index
            row["selection_rule"] = (
                "rank datasets by eligible-patient/sample coverage and limiting-role cells; "
                "within dataset select the sample maximizing min(cancer, reference), outcome-blind"
            )
            selected_rows.append(row)
    selected = pd.DataFrame(selected_rows).reset_index(drop=True)
    selected.insert(0, "panel_order", range(1, len(selected) + 1))
    selected["patient_selected_units"] = selected.groupby("analysis_patient_id")[
        "analysis_patient_id"
    ].transform("size")
    selected["paired_patient_in_panel"] = selected["patient_selected_units"] > 1
    if selected.duplicated(["sample_type", "dataset"]).any():
        raise ValueError("selected panel contains duplicate dataset within a state")
    return selected


def run_selection(
    h5ad: Path,
    h5ad_receipt: Path,
    output_dir: Path,
    states: list[str],
    max_datasets_per_state: int,
    min_cancer: int,
    min_reference: int,
) -> dict[str, object]:
    resource = verify_upstream_h5ad_receipt(h5ad, h5ad_receipt)
    eligible = build_eligible_table(h5ad, min_cancer, min_reference)
    selected = select_panel(eligible, states, max_datasets_per_state)
    output_dir.mkdir(parents=True, exist_ok=True)
    eligible_path = output_dir / "P0_cna_eligible_samples.tsv"
    selected_path = output_dir / "P0_cna_selected_panel.tsv"
    eligible.to_csv(eligible_path, sep="\t", index=False, lineterminator="\n")
    selected.to_csv(selected_path, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "real_data_cna_panel_selected_not_yet_executed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "h5ad_resource": resource,
        "min_cancer": min_cancer,
        "min_reference": min_reference,
        "states": states,
        "max_datasets_per_state": max_datasets_per_state,
        "eligible_samples": int(eligible["eligible"].sum()),
        "eligible_datasets": int(eligible.loc[eligible["eligible"], "dataset"].nunique()),
        "eligible_patients": int(eligible.loc[eligible["eligible"], "analysis_patient_id"].nunique()),
        "selected_samples": int(len(selected)),
        "selected_datasets": int(selected["dataset"].nunique()),
        "selected_patients": int(selected["analysis_patient_id"].nunique()),
        "selection_is_outcome_blind": True,
        "claim_boundary": (
            "This table selects a diagnostic CNA panel from author labels and cell-count coverage. "
            "It does not establish malignancy or make state comparisons exchangeable across datasets."
        ),
        "outputs": {},
    }
    for path in (eligible_path, selected_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (output_dir / "P0_cna_panel_selection_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--h5ad-receipt", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--state", action="append", dest="states", default=None)
    parser.add_argument("--max-datasets-per-state", type=int, default=4)
    parser.add_argument("--min-cancer", type=int, default=50)
    parser.add_argument("--min-reference", type=int, default=50)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = run_selection(
        args.h5ad,
        args.h5ad_receipt,
        args.output_dir,
        args.states or ["polyp", "tumor", "metastasis"],
        args.max_datasets_per_state,
        args.min_cancer,
        args.min_reference,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
