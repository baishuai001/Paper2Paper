#!/usr/bin/env python3
"""Audit CRC-atlas metadata needed by the frozen TF gate without reading expression."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from crc_cohort import eligible_patient_table
from gate1_common import clean_text, read_obs_column, sha256_file, write_json


FIELDS = [
    "donor_id",
    "sample_id",
    "dataset",
    "study_id",
    "medical_condition",
    "sample_type",
    "tumor_source",
    "anatomic_region",
    "anatomic_location",
    "treatment_status_before_resection",
    "treatment_drug",
    "enrichment_cell_types",
    "tissue_cell_state",
    "cell_type_coarse_crc_atlas",
    "cell_type_middle_crc_atlas",
    "cell_type_fine_crc_atlas",
    "immune_infiltration_type",
    "CMS_type",
    "microsatellite_status",
    "tumor_stage",
    "age",
    "sex",
]


def normalized(series: pd.Series) -> pd.Series:
    return series.map(clean_text).str.lower()


def run(h5ad: Path, output_dir: Path, expected_bytes: int) -> dict[str, object]:
    if not h5ad.is_file():
        raise FileNotFoundError(h5ad)
    if h5ad.stat().st_size != expected_bytes:
        raise ValueError(f"H5AD byte mismatch: expected {expected_bytes}, got {h5ad.stat().st_size}")
    with h5py.File(h5ad, "r") as handle:
        if "obs" not in handle:
            raise ValueError("H5AD has no obs")
        obs_group = handle["obs"]
        missing = sorted(set(FIELDS) - set(obs_group))
        if missing:
            raise ValueError(f"H5AD missing gate fields: {missing}")
        data = {field: read_obs_column(obs_group, field) for field in FIELDS}
    lengths = {len(values) for values in data.values()}
    if len(lengths) != 1:
        raise ValueError(f"inconsistent obs lengths: {lengths}")
    cells = pd.DataFrame(data)
    for field in FIELDS:
        cells[field] = cells[field].map(clean_text)

    sample_type = normalized(cells["sample_type"])
    tumor_source = normalized(cells["tumor_source"])
    medical = normalized(cells["medical_condition"])
    cells["candidate_primary_tumor"] = (
        sample_type.str.contains("tumor", regex=False)
        & ~sample_type.str.contains("metasta", regex=False)
        & ~tumor_source.str.contains("metasta", regex=False)
        & ~medical.str.contains("normal", regex=False)
        & ~medical.str.contains("polyp", regex=False)
    )
    cells["is_author_cancer"] = cells["cell_type_coarse_crc_atlas"].eq("Cancer cell")
    cells["has_immune_type"] = cells["immune_infiltration_type"].ne("")

    output_dir.mkdir(parents=True, exist_ok=True)
    level_frames: list[pd.DataFrame] = []
    for field in FIELDS:
        cell_n = cells[field].value_counts(dropna=False, sort=False).rename("cells")
        patient_n = cells[[field, "donor_id"]].drop_duplicates().groupby(field, dropna=False).size().rename("patients")
        sample_n = cells[[field, "sample_id"]].drop_duplicates().groupby(field, dropna=False).size().rename("samples")
        summary = pd.concat([cell_n, patient_n, sample_n], axis=1).fillna(0).astype(int).rename_axis("value").reset_index()
        summary.insert(0, "field", field)
        level_frames.append(summary)
    levels = pd.concat(level_frames, ignore_index=True)
    levels_path = output_dir / "G1_metadata_levels.tsv"
    levels.to_csv(levels_path, sep="\t", index=False, lineterminator="\n")

    patient = cells.groupby("donor_id", observed=True, sort=True).size().rename("all_cells").to_frame()
    for field in [
        "immune_infiltration_type",
        "dataset",
        "study_id",
        "sample_id",
        "sample_type",
        "tumor_source",
        "treatment_status_before_resection",
        "treatment_drug",
        "enrichment_cell_types",
        "CMS_type",
        "microsatellite_status",
        "tumor_stage",
        "anatomic_location",
        "age",
        "sex",
    ]:
        pairs = cells[["donor_id", field]].drop_duplicates()
        pairs = pairs[pairs[field].ne("")]
        patient[field] = pairs.groupby("donor_id", observed=True)[field].agg(lambda x: "|".join(sorted(x))).reindex(patient.index).fillna("")
        patient[f"{field}_n"] = pairs.groupby("donor_id", observed=True).size().reindex(patient.index).fillna(0).astype(int)
    grouped = cells.groupby("donor_id", observed=True, sort=True)
    patient["primary_cells"] = grouped["candidate_primary_tumor"].sum().astype(int)
    cells["candidate_primary_cancer"] = cells["candidate_primary_tumor"] & cells["is_author_cancer"]
    patient["primary_cancer_cells"] = grouped["candidate_primary_cancer"].sum().astype(int)
    patient["immune_type_cells"] = grouped["has_immune_type"].sum().astype(int)
    patient = patient.reset_index()
    patient["immune_label_unique"] = patient["immune_infiltration_type_n"].eq(1)
    patient_path = output_dir / "G1_patient_metadata_audit.tsv"
    patient.to_csv(patient_path, sep="\t", index=False, lineterminator="\n")

    scoped_patient, eligible_mask = eligible_patient_table(cells)
    typed = scoped_patient[scoped_patient["immune_label_n"].eq(1)].copy()
    coverage = (
        typed.groupby(["immune_label", "dataset"], observed=True, sort=True)
        .agg(
            patients=("donor_id", "nunique"),
            patients_ge20_cancer=("cancer_cells", lambda x: int((x >= 20).sum())),
            patients_ge50_cancer=("cancer_cells", lambda x: int((x >= 50).sum())),
            patients_ge100_cancer=("cancer_cells", lambda x: int((x >= 100).sum())),
            primary_cancer_cells=("cancer_cells", "sum"),
        )
        .reset_index()
    )
    coverage_path = output_dir / "G1_label_dataset_coverage.tsv"
    coverage.to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")

    contradictions = scoped_patient[scoped_patient["immune_label_n"].ne(1)]
    main = typed[typed["cancer_cells"].ge(50)].copy()
    m_counts = main.loc[main["immune_label"].eq("M"), "dataset"].value_counts()
    cross = pd.crosstab(main["dataset"], main["immune_label"])
    cross["nonM"] = cross.drop(columns=["M"], errors="ignore").sum(axis=1)
    informative = cross[(cross.get("M", 0) >= 3) & (cross["nonM"] >= 3)]
    n_m = int(main["immune_label"].eq("M").sum())
    n_nonm = int(main["immune_label"].ne("M").sum())
    max_m_share = float(m_counts.max() / n_m) if n_m else 1.0
    minimums = {
        "unique_labels": bool(contradictions.empty),
        "M_ge_15": n_m >= 15,
        "nonM_ge_30": n_nonm >= 30,
        "informative_datasets_ge_3": int(len(informative)) >= 3,
        "max_M_dataset_share_le_0_60": max_m_share <= 0.60,
    }
    receipt = {
        "status": "passed" if all(minimums.values()) else "failed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "h5ad_bytes": h5ad.stat().st_size,
        "n_obs": int(len(cells)),
        "patients": int(cells["donor_id"].nunique()),
        "samples": int(cells["sample_id"].nunique()),
        "datasets": int(cells["dataset"].nunique()),
        "typed_patients": int(typed["donor_id"].nunique()),
        "eligible_cells": int(eligible_mask.sum()),
        "typed_patients_ge50_primary_cancer": int(len(main)),
        "M_patients_ge50": n_m,
        "nonM_patients_ge50": n_nonm,
        "informative_datasets_ge3_each": int(len(informative)),
        "max_M_dataset_share": max_m_share,
        "minimum_conditions": minimums,
        "contradictory_patient_labels": int(len(contradictions)),
        "outputs": {},
    }
    for path in (levels_path, patient_path, coverage_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    receipt_path = output_dir / "G1_input_audit_receipt.json"
    write_json(receipt_path, receipt)
    print(receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-bytes", type=int, default=30_875_155_333)
    args = parser.parse_args()
    receipt = run(args.h5ad, args.output_dir, args.expected_bytes)
    return 0 if receipt["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
