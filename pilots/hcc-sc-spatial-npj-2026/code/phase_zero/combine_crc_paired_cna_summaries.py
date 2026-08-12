#!/usr/bin/env python3
"""Combine Liu and non-Liu CNA shards into paired-sample review tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


PAIR_COLUMNS = [
    "panel_order",
    "paired_group",
    "contrast",
    "analysis_role",
    "pair_member",
]


def _read(path: Path, required: set[str]) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"required table is missing or empty: {path}")
    table = pd.read_csv(path, sep="\t")
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{path.name} missing fields: {sorted(missing)}")
    return table


def _load_shards(
    summary_dirs: list[Path],
    filename: str,
    required: set[str],
    *,
    unique_panel_order_per_row: bool,
) -> pd.DataFrame:
    """Load disjoint panel shards without mistaking cell rows for shard overlap.

    Unit summaries must contain one row per ``panel_order``. Cell-evidence tables
    legitimately contain many rows per sample, so for those tables we compare the
    *set* of panel orders between shards while allowing repeats within one shard.
    """

    tables: list[pd.DataFrame] = []
    seen_orders: set[object] = set()
    for directory in summary_dirs:
        table = _read(directory / filename, required)
        if "panel_order" in table:
            if unique_panel_order_per_row and table["panel_order"].duplicated().any():
                raise ValueError(f"{filename}: panel order repeats within {directory}")
            shard_orders = set(table["panel_order"])
            overlap = seen_orders & shard_orders
            if overlap:
                raise ValueError(
                    f"{filename}: panel orders overlap across shards: {sorted(overlap)}"
                )
            seen_orders.update(shard_orders)
        tables.append(table)
    return pd.concat(tables, ignore_index=True)


def classify_samples(units: pd.DataFrame) -> pd.DataFrame:
    units = units.copy()
    units["normal_reference_unexpected_support_fraction"] = (
        units["normal_reference_unexpected_malignancy_support"]
        / units["normal_reference_cells"]
    )
    units["author_cancer_method_unresolved_fraction"] = (
        units["author_cancer_method_unresolved"] / units["author_cancer_cells"]
    )
    units["reference_review"] = "pass_le_5pct"
    units.loc[
        units["normal_reference_unexpected_support_fraction"].gt(0.05)
        & units["normal_reference_unexpected_support_fraction"].le(0.10),
        "reference_review",
    ] = "warning_gt_5pct_le_10pct"
    units.loc[
        units["normal_reference_unexpected_support_fraction"].gt(0.10),
        "reference_review",
    ] = "sensitivity_required_gt_10pct"
    units["two_method_cell_review"] = "adequate_ge_30"
    units.loc[
        units["author_cancer_two_method_support"].lt(30), "two_method_cell_review"
    ] = "low_lt_30"
    units["method_stability_review"] = "pass_le_50pct_unresolved"
    units.loc[
        units["author_cancer_method_unresolved_fraction"].gt(0.50),
        "method_stability_review",
    ] = "warning_gt_50pct_unresolved"

    units["next_action"] = "complete_full_cancer_identification"
    units.loc[
        units["reference_review"].eq("warning_gt_5pct_le_10pct"), "next_action"
    ] = "complete_full_identification_with_reference_warning"
    units.loc[
        units["reference_review"].eq("sensitivity_required_gt_10pct"), "next_action"
    ] = "run_one_lymphoid_reference_sensitivity_before_full_identification"
    return units


def summarize_pairs(units: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group_id, group in units.groupby("paired_group", observed=True, sort=True):
        members = set(group["pair_member"])
        contrast = str(group["contrast"].iloc[0])
        required = (
            {"primary", "liver_metastasis"}
            if contrast == "primary_vs_liver_metastasis"
            else {"adenoma", "carcinoma"}
        )
        if not required.issubset(members):
            raise ValueError(f"paired group lost a required member: {group_id}, {members}")
        acceptable = group[~group["reference_review"].eq("sensitivity_required_gt_10pct")]
        acceptable_members = set(acceptable["pair_member"])
        rows.append(
            {
                "paired_group": group_id,
                "dataset": str(group["dataset"].iloc[0]),
                "analysis_patient_id": str(group["analysis_patient_id"].iloc[0]),
                "analysis_role": str(group["analysis_role"].iloc[0]),
                "contrast": contrast,
                "samples": int(len(group)),
                "members": ";".join(sorted(members)),
                "samples_requiring_reference_sensitivity": int(
                    group["reference_review"].eq("sensitivity_required_gt_10pct").sum()
                ),
                "samples_with_low_two_method_cells": int(
                    group["two_method_cell_review"].eq("low_lt_30").sum()
                ),
                "two_method_supported_cells": int(
                    group["author_cancer_two_method_support"].sum()
                ),
                "pair_screen_status": (
                    "ready_for_full_identification"
                    if required.issubset(acceptable_members)
                    else "pending_reference_sensitivity"
                ),
            }
        )
    return pd.DataFrame(rows)


def combine(
    panel_path: Path,
    summary_dirs: list[Path],
    output_dir: Path,
    stage: str,
) -> dict[str, object]:
    panel = _read(panel_path, set(PAIR_COLUMNS) | {"dataset", "sample_id", "analysis_patient_id"})
    units = _load_shards(
        summary_dirs,
        "P0_cna_panel_unit_summary.tsv",
        {
            "panel_order",
            "dataset",
            "sample_id",
            "analysis_patient_id",
            "author_cancer_cells",
            "author_cancer_two_method_support",
            "author_cancer_method_unresolved",
            "normal_reference_cells",
            "normal_reference_unexpected_malignancy_support",
        },
        unique_panel_order_per_row=True,
    )
    if set(panel["panel_order"]) != set(units["panel_order"]):
        raise ValueError("paired manifest and shard summaries have different panel orders")
    units = units.merge(panel[PAIR_COLUMNS], on="panel_order", validate="one_to_one")
    units = classify_samples(units).sort_values("panel_order")
    pairs = summarize_pairs(units)

    cells = _load_shards(
        summary_dirs,
        "P0_cna_panel_cell_evidence.tsv",
        {
            "panel_order",
            "cna_input_role",
            "either_malignancy_support",
            "cell_type_coarse_crc_atlas",
        },
        unique_panel_order_per_row=False,
    )
    references = cells[cells["cna_input_role"].eq("known_normal_reference")].copy()
    references["unexpected"] = references["either_malignancy_support"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    reference_labels = (
        references.groupby(["panel_order", "cell_type_coarse_crc_atlas"], observed=True)
        .agg(reference_cells=("unexpected", "size"), unexpected_support=("unexpected", "sum"))
        .reset_index()
    )
    reference_labels["unexpected_support_fraction"] = (
        reference_labels["unexpected_support"] / reference_labels["reference_cells"]
    )
    reference_labels = reference_labels.merge(
        panel[["panel_order", "dataset", "sample_id", "analysis_patient_id"]],
        on="panel_order",
        validate="many_to_one",
    )

    cohort = (
        pairs.groupby(["dataset", "analysis_role", "contrast"], observed=True)
        .agg(
            paired_patients=("paired_group", "size"),
            pairs_ready_for_full_identification=(
                "pair_screen_status",
                lambda value: int((value == "ready_for_full_identification").sum()),
            ),
            pairs_pending_reference_sensitivity=(
                "pair_screen_status",
                lambda value: int((value == "pending_reference_sensitivity").sum()),
            ),
            two_method_supported_cells=("two_method_supported_cells", "sum"),
        )
        .reset_index()
    )

    output_dir.mkdir(parents=True, exist_ok=False)
    outputs = {
        "P0_paired_cna_sample_review.tsv": units,
        "P0_paired_cna_pair_review.tsv": pairs,
        "P0_paired_cna_cohort_review.tsv": cohort,
        "P0_paired_cna_reference_label_review.tsv": reference_labels,
    }
    for name, table in outputs.items():
        table.to_csv(output_dir / name, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": f"paired_cna_{stage}_review_tables_completed",
        "stage": stage,
        "samples": int(len(units)),
        "paired_patients": int(len(pairs)),
        "samples_requiring_reference_sensitivity": int(
            units["reference_review"].eq("sensitivity_required_gt_10pct").sum()
        ),
        "pairs_pending_reference_sensitivity": int(
            pairs["pair_screen_status"].eq("pending_reference_sensitivity").sum()
        ),
        "samples_with_at_least_30_two_method_cells": int(
            units["author_cancer_two_method_support"].ge(30).sum()
        ),
        "formal_cnmf_started": False,
        "claim_boundary": (
            "These are expression-derived CNA support and reference-integrity review tables, not "
            "DNA truth. The script does not execute cNMF."
        ),
    }
    (output_dir / "P0_paired_cna_review_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--summary-dir", required=True, type=Path, action="append")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--stage", choices=["screen", "full"], required=True)
    args = parser.parse_args()
    combine(args.panel, args.summary_dir, args.output_dir, args.stage)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
