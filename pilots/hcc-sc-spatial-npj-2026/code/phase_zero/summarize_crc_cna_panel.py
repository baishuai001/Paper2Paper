#!/usr/bin/env python3
"""Aggregate a frozen CRC CNA panel without treating inferred calls as DNA truth."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from run_crc_cna_panel import load_panel
from summarize_crc_cna_concordance import build_cell_evidence


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def read_completed_state(run_root: Path) -> dict[str, object]:
    state_path = run_root / "P0_cna_panel_run_receipt.json"
    if not state_path.is_file() or state_path.stat().st_size == 0:
        raise FileNotFoundError(f"completed CNA panel receipt missing or empty: {state_path}")
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if state.get("status") != "completed" or state.get("units_failed") != 0:
        raise ValueError(
            f"CNA panel is not a complete successful run: status={state.get('status')!r}, "
            f"units_failed={state.get('units_failed')!r}"
        )
    units = state.get("units")
    if not isinstance(units, list) or not units or any(unit.get("status") != "completed" for unit in units):
        raise ValueError("CNA panel receipt does not contain only completed analysis units")
    if state.get("units_total") != len(units) or state.get("units_completed") != len(units):
        raise ValueError("CNA panel receipt unit totals disagree")
    return state


def classify_cell_evidence(table: pd.DataFrame) -> pd.Series:
    role = table["cna_input_role"]
    both = table["copykat_malignancy_support"] & table["scevan_malignancy_support"]
    either = table["either_malignancy_support"]
    unresolved = table["method_unresolved"]
    classes = pd.Series("unclassified", index=table.index, dtype="object")
    author = role.eq("author_cancer")
    reference = role.eq("known_normal_reference")
    other = ~(author | reference)
    classes.loc[author & both] = "author_cancer_two_method_malignancy_support"
    classes.loc[author & either & ~both] = "author_cancer_one_method_malignancy_support"
    classes.loc[author & ~either & unresolved] = "author_cancer_method_unresolved"
    classes.loc[author & ~either & ~unresolved] = "author_cancer_without_cna_support_not_nonmalignancy_proof"
    classes.loc[reference & either] = "normal_reference_with_unexpected_malignancy_support"
    classes.loc[reference & ~either & unresolved] = "normal_reference_method_unresolved"
    classes.loc[reference & ~either & ~unresolved] = "normal_reference_without_malignancy_support_not_truth"
    classes.loc[other] = "non_author_cancer_non_reference_context_only"
    if classes.eq("unclassified").any():
        raise ValueError("some CNA cells could not be classified")
    return classes


def _unit_summary(joined: pd.DataFrame, row: object) -> dict[str, object]:
    author = joined[joined["cna_input_role"].eq("author_cancer")]
    reference = joined[joined["cna_input_role"].eq("known_normal_reference")]
    both = author["copykat_malignancy_support"] & author["scevan_malignancy_support"]
    one_method = author["either_malignancy_support"] & ~both
    unresolved_without_support = ~author["either_malignancy_support"] & author["method_unresolved"]
    defined_without_support = ~author["either_malignancy_support"] & ~author["method_unresolved"]
    defined = author[author["both_defined"]]
    return {
        "panel_order": int(row.panel_order),
        "dataset": str(row.dataset),
        "sample_id": str(row.sample_id),
        "analysis_patient_id": str(row.analysis_patient_id),
        "sample_type": str(row.sample_type),
        "tissue": str(row.tissue),
        "author_cancer_cells": int(len(author)),
        "author_cancer_two_method_support": int(both.sum()),
        "author_cancer_one_method_support": int(one_method.sum()),
        "author_cancer_unresolved_without_support": int(unresolved_without_support.sum()),
        "author_cancer_defined_without_support_not_nonmalignancy_proof": int(
            defined_without_support.sum()
        ),
        "author_cancer_one_or_two_method_support": int(author["either_malignancy_support"].sum()),
        "author_cancer_method_unresolved": int(author["method_unresolved"].sum()),
        "author_cancer_two_method_support_fraction": float(both.mean()) if len(author) else None,
        "author_cancer_defined_binary_agreement_fraction": (
            float(defined["defined_binary_agreement"].mean()) if len(defined) else None
        ),
        "normal_reference_cells": int(len(reference)),
        "normal_reference_copykat_aneuploid": int(
            reference["copykat_malignancy_support"].sum()
        ),
        "normal_reference_scevan_tumor": int(
            reference["scevan_malignancy_support"].sum()
        ),
        "normal_reference_unexpected_malignancy_support": int(
            reference["either_malignancy_support"].sum()
        ),
        "normal_reference_method_unresolved": int(reference["method_unresolved"].sum()),
    }


def summarize_panel(panel_path: Path, run_root: Path, output_dir: Path) -> dict[str, object]:
    panel = load_panel(panel_path)
    state = read_completed_state(run_root)
    state_units = {int(unit["panel_order"]): unit for unit in state["units"]}
    if set(state_units) != set(panel["panel_order"].astype(int)):
        raise ValueError("frozen panel orders and completed run units disagree")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix with existing CNA panel summary: {output_dir}")

    cell_tables: list[pd.DataFrame] = []
    unit_rows: list[dict[str, object]] = []
    for row in panel.itertuples(index=False):
        order = int(row.panel_order)
        state_unit = state_units[order]
        expected_identity = {
            "dataset": str(row.dataset),
            "sample_id": str(row.sample_id),
            "patient_id": str(row.analysis_patient_id),
            "sample_type": str(row.sample_type),
            "tissue": str(row.tissue),
        }
        observed_identity = {key: str(state_unit.get(key)) for key in expected_identity}
        if observed_identity != expected_identity:
            raise ValueError(
                f"CNA unit {order} identity differs from frozen panel: "
                f"expected={expected_identity}, observed={observed_identity}"
            )
        unit_dir = run_root / f"unit_{order:02d}"
        joined = build_cell_evidence(
            unit_dir / "input" / "cell_metadata.tsv",
            unit_dir / "copykat" / "copykat_prediction.tsv",
            unit_dir / "scevan" / "scevan_prediction.tsv",
        )
        joined.insert(0, "panel_order", order)
        joined.insert(1, "panel_dataset", str(row.dataset))
        joined.insert(2, "panel_sample_id", str(row.sample_id))
        joined.insert(3, "panel_patient_id", str(row.analysis_patient_id))
        joined.insert(4, "panel_sample_type", str(row.sample_type))
        joined.insert(5, "panel_tissue", str(row.tissue))
        joined["cna_selection_class"] = classify_cell_evidence(joined)
        cell_tables.append(joined)
        unit_rows.append(_unit_summary(joined, row))

    cells = pd.concat(cell_tables, ignore_index=True)
    if cells["cell_id"].duplicated().any():
        duplicates = cells.loc[cells["cell_id"].duplicated(False), "cell_id"].head(3).tolist()
        raise ValueError(f"cell IDs recur across CNA panel units: {duplicates}")
    units = pd.DataFrame(unit_rows).sort_values("panel_order")
    state_summary = (
        units.groupby("sample_type", observed=True, sort=True)
        .agg(
            samples=("sample_id", "size"),
            datasets=("dataset", "nunique"),
            patients=("analysis_patient_id", "nunique"),
            author_cancer_cells=("author_cancer_cells", "sum"),
            author_cancer_two_method_support=("author_cancer_two_method_support", "sum"),
            author_cancer_one_method_support=("author_cancer_one_method_support", "sum"),
            author_cancer_unresolved_without_support=(
                "author_cancer_unresolved_without_support",
                "sum",
            ),
            author_cancer_defined_without_support_not_nonmalignancy_proof=(
                "author_cancer_defined_without_support_not_nonmalignancy_proof",
                "sum",
            ),
            author_cancer_one_or_two_method_support=("author_cancer_one_or_two_method_support", "sum"),
            author_cancer_method_unresolved=("author_cancer_method_unresolved", "sum"),
            normal_reference_cells=("normal_reference_cells", "sum"),
            normal_reference_copykat_aneuploid=("normal_reference_copykat_aneuploid", "sum"),
            normal_reference_scevan_tumor=("normal_reference_scevan_tumor", "sum"),
            normal_reference_unexpected_malignancy_support=(
                "normal_reference_unexpected_malignancy_support",
                "sum",
            ),
        )
        .reset_index()
    )
    state_summary["author_cancer_two_method_support_fraction"] = (
        state_summary["author_cancer_two_method_support"] / state_summary["author_cancer_cells"]
    )
    state_summary["normal_reference_unexpected_support_fraction"] = (
        state_summary["normal_reference_unexpected_malignancy_support"]
        / state_summary["normal_reference_cells"]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    cell_path = output_dir / "P0_cna_panel_cell_evidence.tsv"
    unit_path = output_dir / "P0_cna_panel_unit_summary.tsv"
    state_path = output_dir / "P0_cna_panel_state_summary.tsv"
    cells.to_csv(cell_path, sep="\t", index=False, lineterminator="\n")
    units.to_csv(unit_path, sep="\t", index=False, lineterminator="\n")
    state_summary.to_csv(state_path, sep="\t", index=False, lineterminator="\n")

    receipt: dict[str, object] = {
        "status": "diagnostic_panel_completed_cna_supported_cell_set_available",
        "scientific_gate": "not_passed_diagnostic_panel_without_dna_truth",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "panel_sha256": sha256_file(panel_path),
        "run_receipt_sha256": sha256_file(run_root / "P0_cna_panel_run_receipt.json"),
        "samples": int(len(units)),
        "datasets": int(units["dataset"].nunique()),
        "patients": int(units["analysis_patient_id"].nunique()),
        "states": int(units["sample_type"].nunique()),
        "cells_evaluated": int(len(cells)),
        "author_cancer_two_method_support": int(units["author_cancer_two_method_support"].sum()),
        "units_with_at_least_30_two_method_supported_cells": int(
            units["author_cancer_two_method_support"].ge(30).sum()
        ),
        "normal_reference_unexpected_malignancy_support": int(
            units["normal_reference_unexpected_malignancy_support"].sum()
        ),
        "normal_reference_copykat_aneuploid": int(
            units["normal_reference_copykat_aneuploid"].sum()
        ),
        "normal_reference_scevan_tumor": int(
            units["normal_reference_scevan_tumor"].sum()
        ),
        "units_with_unexpected_normal_reference_support": int(
            units["normal_reference_unexpected_malignancy_support"].gt(0).sum()
        ),
        "maximum_unit_normal_reference_unexpected_support_fraction": float(
            (
                units["normal_reference_unexpected_malignancy_support"]
                / units["normal_reference_cells"]
            ).max()
        ),
        "selection_for_cnmf": "author Cancer candidates supported as aneuploid by CopyKAT and tumor by SCEVAN",
        "claim_boundary": (
            "The two-method intersection is a conservative CNA-supported candidate set, not DNA-validated "
            "malignancy truth. The panel covers selected samples across states and datasets, not every atlas "
            "sample; absence of a CNA call is not evidence of nonmalignancy. State comparisons remain "
            "observational and dataset-confounded."
        ),
        "outputs": {},
    }
    for path in (cell_path, unit_path, state_path):
        receipt["outputs"][path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    receipt_path = output_dir / "P0_cna_panel_summary_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = summarize_panel(args.panel, args.run_root, args.output_dir)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
