#!/usr/bin/env python3
"""Frozen CRC-atlas cohort definitions used by every gate-1 module."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from gate1_common import clean_text, read_obs_column


OBS_FIELDS = [
    "donor_id",
    "sample_id",
    "dataset",
    "study_id",
    "sample_type",
    "tumor_source",
    "medical_condition",
    "treatment_status_before_resection",
    "enrichment_cell_types",
    "cell_type_coarse_crc_atlas",
    "cell_type_middle_crc_atlas",
    "immune_infiltration_type",
    "CMS_type",
    "microsatellite_status",
    "anatomic_location",
    "tumor_stage",
    "age",
    "sex",
]


def load_obs(h5ad: Path, fields: list[str] | None = None) -> pd.DataFrame:
    selected = fields or OBS_FIELDS
    with h5py.File(h5ad, "r") as handle:
        missing = sorted(set(selected) - set(handle["obs"]))
        if missing:
            raise ValueError(f"H5AD missing required obs fields: {missing}")
        data = {field: read_obs_column(handle["obs"], field) for field in selected}
    frame = pd.DataFrame(data)
    for field in selected:
        frame[field] = frame[field].map(clean_text)
    return frame


def eligible_primary_mask(obs: pd.DataFrame) -> np.ndarray:
    """Reproduce the CRC paper's primary, untreated, non-enriched cohort."""
    sample = obs["sample_type"].str.lower()
    source = obs["tumor_source"].str.lower()
    medical = obs["medical_condition"].str.lower()
    treatment = obs["treatment_status_before_resection"].str.lower()
    enrichment = obs["enrichment_cell_types"].str.lower()
    return np.asarray(
        sample.eq("tumor")
        & ~source.str.contains("normal", regex=False)
        & ~source.str.contains("metasta", regex=False)
        & ~medical.str.contains("normal", regex=False)
        & ~medical.str.contains("polyp", regex=False)
        & treatment.eq("naive")
        & enrichment.eq("naive"),
        dtype=bool,
    )


def eligible_patient_table(obs: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    eligible = eligible_primary_mask(obs)
    scoped = obs.loc[eligible].copy()
    typed = scoped.loc[scoped["immune_infiltration_type"].ne("")].copy()
    label_n = typed.groupby("donor_id", observed=True)["immune_infiltration_type"].nunique()
    contradictory = set(label_n[label_n.ne(1)].index)

    rows: list[dict[str, object]] = []
    for donor, group in scoped.groupby("donor_id", observed=True, sort=True):
        labels = sorted(set(group["immune_infiltration_type"]) - {""})
        datasets = sorted(set(group["dataset"]) - {""})
        studies = sorted(set(group["study_id"]) - {""})
        if not labels:
            continue
        dataset_counts = group.loc[group["dataset"].ne(""), "dataset"].value_counts()
        study_counts = group.loc[group["study_id"].ne(""), "study_id"].value_counts()
        modal_dataset = sorted(dataset_counts[dataset_counts.eq(dataset_counts.max())].index)[0]
        modal_study = sorted(study_counts[study_counts.eq(study_counts.max())].index)[0]
        rows.append(
            {
                "donor_id": donor,
                "immune_label": "|".join(labels),
                "immune_label_n": len(labels),
                "dataset": modal_dataset,
                "dataset_all": "|".join(datasets),
                "dataset_n": len(datasets),
                "study_id": modal_study,
                "study_id_all": "|".join(studies),
                "study_id_n": len(studies),
                "eligible_cells": int(len(group)),
                "cancer_cells": int(group["cell_type_coarse_crc_atlas"].eq("Cancer cell").sum()),
                "contradictory": donor in contradictory,
            }
        )
    patient = pd.DataFrame(rows)
    if patient.empty:
        raise ValueError("No author-labelled patients in the eligible cohort")
    patient["M"] = patient["immune_label"].eq("M").astype(int)
    return patient, eligible
