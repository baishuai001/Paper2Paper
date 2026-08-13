#!/usr/bin/env python3
"""Manifest-driven CRC-atlas cohort construction."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from common import clean_text, read_h5_column


BASE_FIELDS = [
    "donor_id",
    "sample_id",
    "dataset",
    "study_id",
    "sample_type",
    "tumor_source",
    "medical_condition",
    "treatment_status_before_resection",
    "enrichment_cell_types",
]


@dataclass(frozen=True)
class PhenotypeManifest:
    raw: dict[str, object]

    @classmethod
    def load(cls, path: Path) -> "PhenotypeManifest":
        with path.open(encoding="utf-8") as stream:
            raw = json.load(stream)
        required = {
            "phenotype_id",
            "label_column",
            "case_values",
            "control_values",
            "analysis_cell_column",
            "analysis_cell_values",
            "patient_column",
            "dataset_column",
            "scope",
            "primary_min_cells",
        }
        missing = sorted(required - set(raw))
        if missing:
            raise ValueError(f"Phenotype manifest missing fields: {missing}")
        overlap = set(raw["case_values"]) & set(raw["control_values"])
        if overlap:
            raise ValueError(f"Case/control values overlap: {sorted(overlap)}")
        return cls(raw)

    def __getitem__(self, key: str):
        return self.raw[key]

    @property
    def required_obs_fields(self) -> list[str]:
        return sorted(
            set(BASE_FIELDS)
            | {
                str(self["label_column"]),
                str(self["analysis_cell_column"]),
                str(self["patient_column"]),
                str(self["dataset_column"]),
            }
        )


def load_obs(h5ad: Path, fields: list[str]) -> pd.DataFrame:
    with h5py.File(h5ad, "r") as handle:
        missing = sorted(set(fields) - set(handle["obs"]))
        if missing:
            raise ValueError(f"H5AD missing required obs fields: {missing}")
        data = {field: read_h5_column(handle["obs"], field) for field in fields}
    frame = pd.DataFrame(data)
    for field in fields:
        frame[field] = frame[field].map(clean_text)
    return frame


def scope_mask(obs: pd.DataFrame, manifest: PhenotypeManifest) -> np.ndarray:
    scope = manifest["scope"]
    mask = np.ones(len(obs), dtype=bool)
    equal_rules = {
        "sample_type_equals": "sample_type",
        "treatment_status_before_resection_equals": "treatment_status_before_resection",
        "enrichment_cell_types_equals": "enrichment_cell_types",
    }
    for rule, column in equal_rules.items():
        if rule in scope:
            mask &= obs[column].str.casefold().eq(str(scope[rule]).casefold()).to_numpy()
    contains_rules = {
        "exclude_tumor_source_contains": "tumor_source",
        "exclude_medical_condition_contains": "medical_condition",
    }
    for rule, column in contains_rules.items():
        for token in scope.get(rule, []):
            mask &= ~obs[column].str.casefold().str.contains(str(token).casefold(), regex=False).to_numpy()
    return mask


def patient_table(obs: pd.DataFrame, manifest: PhenotypeManifest) -> tuple[pd.DataFrame, np.ndarray]:
    eligible = scope_mask(obs, manifest)
    label_column = str(manifest["label_column"])
    patient_column = str(manifest["patient_column"])
    dataset_column = str(manifest["dataset_column"])
    allowed = set(manifest["case_values"]) | set(manifest["control_values"])
    scoped = obs.loc[eligible & obs[label_column].isin(allowed)].copy()
    if scoped.empty:
        raise ValueError("No cells in the manifest-defined labelled cohort")

    rows: list[dict[str, object]] = []
    for patient, group in scoped.groupby(patient_column, observed=True, sort=True):
        labels = sorted(set(group[label_column]) & allowed)
        dataset_counts = group.loc[group[dataset_column].ne(""), dataset_column].value_counts()
        if dataset_counts.empty:
            dataset = ""
        else:
            dataset = sorted(dataset_counts[dataset_counts.eq(dataset_counts.max())].index)[0]
        cell_mask = group[str(manifest["analysis_cell_column"])].isin(manifest["analysis_cell_values"])
        rows.append(
            {
                "donor_id": patient,
                "phenotype_id": manifest["phenotype_id"],
                "label": "|".join(labels),
                "label_n": len(labels),
                "group": (
                    "case"
                    if len(labels) == 1 and labels[0] in manifest["case_values"]
                    else "control" if len(labels) == 1 and labels[0] in manifest["control_values"] else "conflict"
                ),
                "dataset": dataset,
                "dataset_n": int(group[dataset_column].replace("", np.nan).nunique()),
                "eligible_cells": int(len(group)),
                "analysis_cells": int(cell_mask.sum()),
            }
        )
    patients = pd.DataFrame(rows).sort_values("donor_id").reset_index(drop=True)
    return patients, eligible
