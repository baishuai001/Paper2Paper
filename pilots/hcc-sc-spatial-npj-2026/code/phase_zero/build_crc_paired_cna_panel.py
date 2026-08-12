#!/usr/bin/env python3
"""Build the prespecified paired CRC CNA analysis panel.

This does not redesign the CRC atlas.  It selects only within-patient,
within-dataset contrasts that can support stage comparisons, while preserving
all eligible samples for a paired patient when technical replicates exist.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {
    "dataset",
    "sample_id",
    "analysis_patient_id",
    "sample_type",
    "tissue",
    "cells",
    "cancer_cells",
    "reference_cells",
    "other_epithelial_cells",
    "eligible",
    "analysis_unit_id",
}

COHORTS = {
    "Liu_2024_mixCD45PosCD45Neg": {
        "analysis_role": "primary_discovery",
        "contrast": "primary_vs_liver_metastasis",
        "expected_pairs": 11,
    },
    "Che_2021": {
        "analysis_role": "independent_replication",
        "contrast": "primary_vs_liver_metastasis",
        "expected_pairs": 5,
    },
    "Ji_2024_scopeV2": {
        "analysis_role": "supporting_case",
        "contrast": "primary_vs_liver_metastasis",
        "expected_pairs": 1,
    },
    "Zheng_2022": {
        "analysis_role": "early_stage_exploration",
        "contrast": "adenoma_vs_carcinoma",
        "expected_pairs": 3,
    },
}


def _truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def build_panel(eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = REQUIRED_COLUMNS - set(eligible.columns)
    if missing:
        raise ValueError(f"eligible-sample table is missing fields: {sorted(missing)}")
    eligible = eligible[_truthy(eligible["eligible"])].copy()
    eligible = eligible[eligible["dataset"].isin(COHORTS)].copy()
    if eligible.empty:
        raise ValueError("no eligible rows from the prespecified paired cohorts")

    selected_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []
    for dataset, contract in COHORTS.items():
        cohort = eligible[eligible["dataset"].eq(dataset)].copy()
        if dataset == "Zheng_2022":
            cohort = cohort[
                cohort["sample_type"].eq("polyp")
                | (
                    cohort["sample_type"].eq("tumor")
                    & cohort["sample_id"].str.contains("carcinoma", case=False, na=False)
                )
            ].copy()
            required_states = {"polyp", "tumor"}
            member_map = {"polyp": "adenoma", "tumor": "carcinoma"}
        else:
            cohort = cohort[cohort["sample_type"].isin(["tumor", "metastasis"])].copy()
            required_states = {"tumor", "metastasis"}
            member_map = {"tumor": "primary", "metastasis": "liver_metastasis"}

        states = cohort.groupby("analysis_patient_id", observed=True)["sample_type"].agg(set)
        paired_patients = states[states.map(lambda value: required_states.issubset(value))].index
        paired = cohort[cohort["analysis_patient_id"].isin(paired_patients)].copy()
        observed_pairs = int(paired["analysis_patient_id"].nunique())
        if observed_pairs != contract["expected_pairs"]:
            raise ValueError(
                f"{dataset}: expected {contract['expected_pairs']} paired patients, "
                f"found {observed_pairs}"
            )
        if paired.empty:
            raise ValueError(f"{dataset}: paired selection is empty")

        paired["paired_group"] = dataset + "::" + paired["analysis_patient_id"].astype(str)
        paired["contrast"] = contract["contrast"]
        paired["analysis_role"] = contract["analysis_role"]
        paired["pair_member"] = paired["sample_type"].map(member_map)
        paired["selection_rule"] = (
            "eligible same-dataset patient with both prespecified biological states; "
            "retain all eligible state samples for that patient"
        )
        selected_parts.append(paired)
        summary_rows.append(
            {
                "dataset": dataset,
                "analysis_role": contract["analysis_role"],
                "contrast": contract["contrast"],
                "paired_patients": observed_pairs,
                "selected_samples": int(len(paired)),
                "primary_or_carcinoma_samples": int(
                    paired["pair_member"].isin(["primary", "carcinoma"]).sum()
                ),
                "comparison_samples": int(
                    paired["pair_member"].isin(["liver_metastasis", "adenoma"]).sum()
                ),
                "author_cancer_cells": int(paired["cancer_cells"].sum()),
                "known_normal_reference_cells": int(paired["reference_cells"].sum()),
            }
        )

    panel = pd.concat(selected_parts, ignore_index=True)
    role_order = {
        "primary_discovery": 1,
        "independent_replication": 2,
        "supporting_case": 3,
        "early_stage_exploration": 4,
    }
    member_order = {
        "primary": 1,
        "liver_metastasis": 2,
        "adenoma": 1,
        "carcinoma": 2,
    }
    panel["_role_order"] = panel["analysis_role"].map(role_order)
    panel["_member_order"] = panel["pair_member"].map(member_order)
    panel = panel.sort_values(
        ["_role_order", "analysis_patient_id", "_member_order", "sample_id"]
    ).drop(columns=["_role_order", "_member_order"])
    panel.insert(0, "panel_order", range(1, len(panel) + 1))
    if panel.duplicated(["sample_id", "analysis_patient_id"]).any():
        raise ValueError("paired panel contains duplicate sample-patient analysis units")
    if panel["pair_member"].isna().any():
        raise ValueError("paired panel contains an unmapped biological state")

    summary = pd.DataFrame(summary_rows)
    return panel.reset_index(drop=True), summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eligible", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.eligible, sep="\t")
    panel, summary = build_panel(table)
    for path in (args.panel, args.summary, args.receipt):
        path.parent.mkdir(parents=True, exist_ok=True)
    panel.to_csv(args.panel, sep="\t", index=False, lineterminator="\n")
    summary.to_csv(args.summary, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "paired_crc_cna_panel_built",
        "selected_samples": int(len(panel)),
        "paired_patients": int(panel["paired_group"].nunique()),
        "cohorts": int(panel["dataset"].nunique()),
        "contrasts": sorted(panel["contrast"].unique()),
        "panel": str(args.panel.resolve()),
        "summary": str(args.summary.resolve()),
        "claim_boundary": (
            "This manifest defines within-dataset paired comparisons. It does not redesign the "
            "full CRC atlas, and unpaired cohorts remain available for later program recurrence "
            "or external validation rather than primary stage-effect estimation."
        ),
    }
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
