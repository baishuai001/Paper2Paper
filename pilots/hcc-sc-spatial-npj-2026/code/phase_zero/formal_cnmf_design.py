#!/usr/bin/env python3
"""Pure table helpers for the frozen CRC formal cNMF design."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


SELECTED_CLASS = "author_cancer_two_method_malignancy_support"
ONE_METHOD_CLASS = "author_cancer_one_method_malignancy_support"
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
    "method_unresolved",
    "cna_selection_class",
}


def parse_bool(series: pd.Series, field: str) -> pd.Series:
    normalized = series.astype(str).str.strip().str.lower()
    mapping = {"true": True, "false": False, "1": True, "0": False}
    invalid = sorted(set(normalized) - set(mapping))
    if invalid:
        raise ValueError(f"{field} contains invalid booleans: {invalid[:3]}")
    return normalized.map(mapping).astype(bool)


def load_evidence_tables(paths: list[Path]) -> pd.DataFrame:
    if not paths:
        raise ValueError("at least one CNA cell-evidence table is required")
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"CNA cell-evidence table missing or empty: {path}")
        frame = pd.read_csv(path, sep="\t", dtype=str)
        missing = sorted(REQUIRED_EVIDENCE_COLUMNS - set(frame.columns))
        if missing:
            raise ValueError(f"CNA cell-evidence fields missing in {path}: {missing}")
        frame["evidence_source_file"] = path.name
        frames.append(frame)
    table = pd.concat(frames, ignore_index=True)
    if table["cell_id"].isna().any() or table["cell_id"].duplicated().any():
        duplicates = table.loc[table["cell_id"].duplicated(False), "cell_id"].head(3).tolist()
        raise ValueError(f"CNA cell-evidence IDs are missing or duplicated: {duplicates}")
    table["panel_order"] = pd.to_numeric(table["panel_order"], errors="raise").astype(int)
    table["method_unresolved_bool"] = parse_bool(table["method_unresolved"], "method_unresolved")
    return table


def select_cohort_cells(
    table: pd.DataFrame,
    cohort_dataset: str,
    required_states: tuple[str, ...] = ("tumor", "metastasis"),
    expected_patients: int | None = None,
    expected_samples: int | None = None,
    expected_dual_cells: int | None = None,
    selected_classes: tuple[str, ...] = (SELECTED_CLASS,),
    excluded_sample_ids: tuple[str, ...] = (),
    require_paired_states: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    allowed_classes = {SELECTED_CLASS, ONE_METHOD_CLASS}
    if not selected_classes or not set(selected_classes).issubset(allowed_classes):
        raise ValueError(f"unsupported cNMF selection classes: {selected_classes}")
    if SELECTED_CLASS not in selected_classes:
        raise ValueError("cNMF selection must retain the two-method-supported core")
    cohort = table.loc[table["panel_dataset"].eq(cohort_dataset)].copy()
    if cohort.empty:
        raise ValueError(f"cohort dataset not found in CNA evidence: {cohort_dataset}")
    cohort = cohort.loc[cohort["panel_sample_type"].isin(required_states)].copy()
    observed_states = set(cohort["panel_sample_type"].astype(str))
    missing_states = sorted(set(required_states) - observed_states)
    if missing_states:
        raise ValueError(f"cohort lacks required states: {missing_states}")

    references = cohort.loc[cohort["cna_input_role"].eq("known_normal_reference")].copy()
    if references.empty:
        raise ValueError("cohort evidence has no known-normal reference cells")
    references["unexpected_malignancy_support"] = (
        references["copykat.pred"].eq("aneuploid")
        | references["scevan_class"].eq("tumor")
    )
    sample_qc = (
        references.groupby(
            ["panel_order", "panel_sample_id", "panel_patient_id", "panel_sample_type"],
            observed=True,
            sort=True,
        )
        .agg(
            normal_reference_cells=("cell_id", "size"),
            normal_reference_unexpected_support=("unexpected_malignancy_support", "sum"),
            normal_reference_method_unresolved=("method_unresolved_bool", "sum"),
        )
        .reset_index()
    )
    sample_qc["normal_reference_unexpected_support_fraction"] = (
        sample_qc["normal_reference_unexpected_support"] / sample_qc["normal_reference_cells"]
    )
    sample_qc["normal_reference_method_unresolved_fraction"] = (
        sample_qc["normal_reference_method_unresolved"] / sample_qc["normal_reference_cells"]
    )

    selected = cohort.loc[cohort["cna_selection_class"].isin(selected_classes)].copy()
    if excluded_sample_ids:
        missing_exclusions = sorted(set(excluded_sample_ids) - set(cohort["panel_sample_id"]))
        if missing_exclusions:
            raise ValueError(f"excluded sensitivity samples not found: {missing_exclusions}")
        selected = selected.loc[~selected["panel_sample_id"].isin(excluded_sample_ids)].copy()
    copykat_support = selected["copykat.pred"].eq("aneuploid")
    scevan_support = selected["scevan_class"].eq("tumor")
    dual_rows = selected["cna_selection_class"].eq(SELECTED_CLASS)
    one_method_rows = selected["cna_selection_class"].eq(ONE_METHOD_CLASS)
    invalid = selected.loc[
        ~selected["cna_input_role"].eq("author_cancer")
        | (dual_rows & ~(copykat_support & scevan_support))
        | (one_method_rows & ~(copykat_support ^ scevan_support))
    ]
    if selected.empty or not invalid.empty:
        raise ValueError("selected sensitivity cells disagree with their recorded CNA support class")
    selected = selected.merge(
        sample_qc,
        on=["panel_order", "panel_sample_id", "panel_patient_id", "panel_sample_type"],
        how="left",
        validate="many_to_one",
    )
    if selected["normal_reference_cells"].isna().any():
        raise ValueError("a selected Cancer candidate lacks same-sample reference QC")

    observed_patients = selected["panel_patient_id"].nunique()
    observed_samples = selected["panel_sample_id"].nunique()
    if expected_patients is not None and observed_patients != expected_patients:
        raise ValueError(f"expected {expected_patients} patients, observed {observed_patients}")
    if expected_samples is not None and observed_samples != expected_samples:
        raise ValueError(f"expected {expected_samples} samples, observed {observed_samples}")
    if expected_dual_cells is not None and len(selected) != expected_dual_cells:
        raise ValueError(f"expected {expected_dual_cells} selected cells, observed {len(selected)}")
    patient_states = selected[["panel_patient_id", "panel_sample_type"]].drop_duplicates()
    state_counts = patient_states.groupby("panel_patient_id")["panel_sample_type"].nunique()
    incomplete = state_counts.loc[state_counts.ne(len(required_states))]
    if require_paired_states and not incomplete.empty:
        raise ValueError(f"patients lack paired states: {incomplete.index.tolist()}")
    return selected, sample_qc


def _group_rng(seed: int, patient: str, state: str) -> np.random.Generator:
    token = f"{seed}|{patient}|{state}".encode("utf-8")
    derived = int.from_bytes(hashlib.sha256(token).digest()[:8], "big", signed=False)
    return np.random.default_rng(derived)


def balance_patient_state_cells(
    selected: pd.DataFrame,
    max_cells_per_patient_state: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if max_cells_per_patient_state < 30:
        raise ValueError("patient-state cap must be at least 30")
    chosen: list[pd.DataFrame] = []
    for (patient, state), group in selected.groupby(
        ["panel_patient_id", "panel_sample_type"], observed=True, sort=True
    ):
        group = group.sort_values("cell_id", kind="mergesort")
        if len(group) > max_cells_per_patient_state:
            rng = _group_rng(seed, str(patient), str(state))
            positions = np.sort(rng.choice(len(group), max_cells_per_patient_state, replace=False))
            group = group.iloc[positions]
        chosen.append(group)
    balanced = pd.concat(chosen, ignore_index=True).sort_values(
        ["panel_patient_id", "panel_sample_type", "panel_sample_id", "cell_id"],
        kind="mergesort",
    )
    if balanced["cell_id"].duplicated().any():
        raise ValueError("balanced formal cNMF selection contains duplicate cells")

    all_counts = (
        selected.groupby(["panel_patient_id", "panel_sample_type"], observed=True, sort=True)
        .size()
        .rename("available_dual_supported_cells")
        .reset_index()
    )
    chosen_counts = (
        balanced.groupby(["panel_patient_id", "panel_sample_type"], observed=True, sort=True)
        .size()
        .rename("selected_cells")
        .reset_index()
    )
    patient_state = all_counts.merge(
        chosen_counts,
        on=["panel_patient_id", "panel_sample_type"],
        how="left",
        validate="one_to_one",
    )
    patient_state["cap"] = max_cells_per_patient_state
    patient_state["cell_sampling_seed"] = seed

    available_by_sample = (
        selected.groupby(
            ["panel_patient_id", "panel_sample_type", "panel_sample_id"],
            observed=True,
            sort=True,
        )
        .size()
        .rename("available_dual_supported_cells")
        .reset_index()
    )
    selected_by_sample = (
        balanced.groupby(
            ["panel_patient_id", "panel_sample_type", "panel_sample_id"],
            observed=True,
            sort=True,
        )
        .size()
        .rename("selected_cells")
        .reset_index()
    )
    sample_composition = available_by_sample.merge(
        selected_by_sample,
        on=["panel_patient_id", "panel_sample_type", "panel_sample_id"],
        how="left",
        validate="one_to_one",
    ).fillna({"selected_cells": 0})
    sample_composition["selected_cells"] = sample_composition["selected_cells"].astype(int)
    return balanced.reset_index(drop=True), patient_state, sample_composition


def assign_patient_folds(
    obs: pd.DataFrame,
    seed: int,
    required_states: tuple[str, ...] = ("tumor", "metastasis"),
) -> pd.DataFrame:
    required = {"analysis_patient_id", "sample_type"}
    missing = sorted(required - set(obs.columns))
    if missing:
        raise ValueError(f"formal cNMF obs lacks fold fields: {missing}")
    units = (
        obs.groupby(["analysis_patient_id", "sample_type"], observed=True, sort=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    state_counts = units.groupby("analysis_patient_id")["sample_type"].nunique()
    incomplete = state_counts.loc[state_counts.ne(len(required_states))]
    if not incomplete.empty:
        raise ValueError(f"patient holdout requires paired states: {incomplete.index.tolist()}")
    patients = sorted(units["analysis_patient_id"].astype(str).unique())
    if len(patients) < 4:
        raise ValueError("patient holdout requires at least four paired patients")
    rng = np.random.default_rng(seed)
    shuffled = np.asarray(patients, dtype=object)
    rng.shuffle(shuffled)
    assignments = {str(patient): ("A" if i % 2 == 0 else "B") for i, patient in enumerate(shuffled)}
    manifest = (
        units.pivot(index="analysis_patient_id", columns="sample_type", values="cells")
        .fillna(0)
        .reset_index()
    )
    manifest.columns.name = None
    manifest["holdout_fold"] = manifest["analysis_patient_id"].astype(str).map(assignments)
    manifest["total_cells"] = manifest[list(required_states)].sum(axis=1).astype(int)
    manifest["assignment_seed"] = seed
    manifest = manifest.sort_values(["holdout_fold", "analysis_patient_id"], kind="mergesort")
    if set(manifest["holdout_fold"]) != {"A", "B"}:
        raise ValueError("patient holdout did not create both folds")
    return manifest.reset_index(drop=True)
