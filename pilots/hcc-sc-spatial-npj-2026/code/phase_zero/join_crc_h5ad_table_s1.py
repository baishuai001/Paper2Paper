#!/usr/bin/env python3
"""Build a sample-level, bidirectional H5AD↔Table S1 identity receipt."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


REQUIRED_LOGICAL = ["sample_id", "patient_id", "dataset", "study_id", "sample_type", "sample_tissue"]
COMPARE_COLUMNS = ["patient_id", "dataset", "study_id", "sample_type", "sample_tissue"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in "SUO":
        return np.array(
            [value.decode("utf-8", "replace") if isinstance(value, bytes) else str(value) for value in values],
            dtype=object,
        )
    return values.astype(object)


def read_obs_column(obs: h5py.Group, column: str) -> np.ndarray:
    if column not in obs:
        raise ValueError(f"H5AD obs missing required column: {column}")
    item = obs[column]
    if isinstance(item, h5py.Dataset):
        return decode(item[:])
    if not isinstance(item, h5py.Group) or "codes" not in item or "categories" not in item:
        raise ValueError(f"unsupported H5AD obs encoding for {column}: {type(item).__name__}")
    codes = np.asarray(item["codes"][:], dtype=np.int64)
    categories = decode(item["categories"][:])
    values = np.empty(codes.shape[0], dtype=object)
    values[:] = None
    valid = codes >= 0
    if valid.any():
        if int(codes[valid].max()) >= len(categories):
            raise ValueError(f"categorical code exceeds categories for {column}")
        values[valid] = categories[codes[valid]]
    return values


def clean_value(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def canonical_patient(value: object) -> str:
    text = clean_value(value)
    for suffix in ("-smartseq-smartseq",):
        if text.endswith(suffix):
            return text[: -len(suffix)]
    return text


def canonical_tissue(value: object) -> str:
    text = clean_value(value).lower()
    if not text:
        return ""
    if text == "blood":
        return "blood"
    if text == "liver":
        return "liver"
    if "lymph node" in text:
        return "mesenteric lymph nodes"
    colorectal_terms = (
        "colon",
        "rectum",
        "rectosigmoid",
        "colorect",
        "caecum",
        "cecum",
        "mucosa",
    )
    if any(term in text for term in colorectal_terms):
        return "colon"
    return text


def one_value(series: pd.Series) -> str:
    values = sorted({clean_value(value) for value in series if clean_value(value)})
    return values[0] if len(values) == 1 else "|".join(values)


def value_count(series: pd.Series) -> int:
    return len({clean_value(value) for value in series if clean_value(value)})


def read_sample_metadata(
    h5ad: Path,
    field_map: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    field_map = field_map or {column: column for column in REQUIRED_LOGICAL}
    missing_map = sorted(set(REQUIRED_LOGICAL) - set(field_map))
    if missing_map:
        raise ValueError(f"field map missing logical fields: {missing_map}")
    with h5py.File(h5ad, "r") as handle:
        if "obs" not in handle or not isinstance(handle["obs"], h5py.Group):
            raise ValueError("H5AD missing obs dataframe")
        obs = handle["obs"]
        data = {logical: read_obs_column(obs, field_map[logical]) for logical in REQUIRED_LOGICAL}
        lengths = {len(values) for values in data.values()}
        if len(lengths) != 1:
            raise ValueError(f"obs columns have inconsistent lengths: {lengths}")
        n_obs = lengths.pop()

    cells = pd.DataFrame(data)
    for column in REQUIRED_LOGICAL:
        cells[column] = cells[column].map(clean_value)
    cells["patient_id_raw"] = cells["patient_id"]
    cells["sample_tissue_raw"] = cells["sample_tissue"]
    cells["patient_id"] = cells["patient_id"].map(canonical_patient)
    cells["sample_tissue"] = cells["sample_tissue"].map(canonical_tissue)
    required_missing = {
        column: int((cells[column] == "").sum()) for column in REQUIRED_LOGICAL if (cells[column] == "").any()
    }
    if required_missing:
        raise ValueError(f"H5AD contains cells with missing identity fields: {required_missing}")
    grouped = cells.groupby(["sample_id", "patient_id"], sort=True, observed=True)
    summary = grouped.size().rename("cell_count").to_frame()
    check_columns = ["dataset", "study_id", "sample_type", "sample_tissue"]
    for column in check_columns:
        summary[column] = grouped[column].agg(one_value)
        summary[f"{column}_value_count"] = grouped[column].agg(value_count)
    summary["patient_id_raw"] = grouped["patient_id_raw"].agg(one_value)
    summary["sample_tissue_raw"] = grouped["sample_tissue_raw"].agg(one_value)
    analysis_units = summary.reset_index()
    sample_grouped = analysis_units.groupby("sample_id", sort=True, observed=True)
    sample_summary = sample_grouped["cell_count"].sum().to_frame()
    sample_summary["analysis_units"] = sample_grouped.size()
    sample_summary["patient_count"] = sample_grouped["patient_id"].nunique()
    for column in ["dataset", "study_id", "sample_type", "sample_tissue"]:
        sample_summary[column] = sample_grouped[column].agg(one_value)
        sample_summary[f"{column}_value_count"] = sample_grouped[column].agg(value_count)
    sample_summary["pooled_sample"] = sample_summary["patient_count"] > 1
    return analysis_units, sample_summary.reset_index(), n_obs


def build_join(
    h5ad: Path,
    table_samples: Path,
    output_dir: Path,
    expected_obs: int | None = None,
    field_map: dict[str, str] | None = None,
) -> dict[str, object]:
    for path in (h5ad, table_samples):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"input missing or empty: {path}")
    field_map = field_map or {column: column for column in REQUIRED_LOGICAL}
    h5_units, h5_samples, n_obs = read_sample_metadata(h5ad, field_map)
    if expected_obs is not None and n_obs != expected_obs:
        raise ValueError(f"H5AD cell count mismatch: expected {expected_obs}, got {n_obs}")

    conflict_columns = [f"{column}_value_count" for column in ["dataset", "study_id", "sample_type", "sample_tissue"]]
    internal_conflicts = h5_units[(h5_units[conflict_columns] > 1).any(axis=1)].copy()

    table = pd.read_csv(table_samples, sep="\t", dtype=str, keep_default_na=False)
    missing = sorted(set(["sample_id", *COMPARE_COLUMNS]) - set(table.columns))
    if missing:
        raise ValueError(f"Table S1 sample table missing columns: {missing}")
    if table["sample_id"].duplicated().any():
        raise ValueError("Table S1 sample table contains duplicate sample_id")
    table = table[["sample_id", *COMPARE_COLUMNS]].copy()
    table["patient_id"] = table["patient_id"].map(canonical_patient)
    table["sample_tissue"] = table["sample_tissue"].map(canonical_tissue)
    table.columns = ["sample_id", *[f"table_{column}" for column in COMPARE_COLUMNS]]
    h5_keep = h5_units[["sample_id", "cell_count", "patient_id_raw", "sample_tissue_raw", *COMPARE_COLUMNS]].copy()
    h5_keep = h5_keep.rename(columns={column: f"h5ad_{column}" for column in COMPARE_COLUMNS})
    joined = table.merge(h5_keep, on="sample_id", how="outer", indicator=True, validate="one_to_many")
    joined["join_status"] = joined["_merge"].map(
        {"both": "in_both", "left_only": "table_only", "right_only": "h5ad_only"}
    ).astype(str)
    joined = joined.drop(columns="_merge")
    pooled_sample_ids = set(h5_samples.loc[h5_samples["pooled_sample"], "sample_id"].astype(str))
    for column in COMPARE_COLUMNS:
        left = joined[f"table_{column}"].fillna("").astype(str).str.strip()
        right = joined[f"h5ad_{column}"].fillna("").astype(str).str.strip()
        comparison = np.where(
            joined["join_status"] != "in_both",
            "not_comparable",
            np.where(left == right, "true", "false"),
        )
        if column == "patient_id":
            comparison = np.where(
                joined["sample_id"].isin(pooled_sample_ids) & (joined["join_status"] == "in_both"),
                "not_comparable_pooled_sample",
                comparison,
            )
        joined[f"{column}_match"] = comparison

    both = joined["join_status"] == "in_both"
    mismatch_counts = {
        column: int((both & (joined[f"{column}_match"] == "false")).sum()) for column in COMPARE_COLUMNS
    }
    h5ad_only = sorted(set(joined.loc[joined["join_status"] == "h5ad_only", "sample_id"].astype(str)))
    table_only = sorted(set(joined.loc[joined["join_status"] == "table_only", "sample_id"].astype(str)))

    output_dir.mkdir(parents=True, exist_ok=True)
    h5_path = output_dir / "P0_h5ad_sample_summary.tsv"
    unit_path = output_dir / "P0_h5ad_analysis_units.tsv"
    join_path = output_dir / "P0_h5ad_table_s1_join.tsv"
    conflict_path = output_dir / "P0_h5ad_internal_conflicts.tsv"
    summary_path = output_dir / "P0_h5ad_table_s1_join_summary.tsv"
    receipt_path = output_dir / "P0_h5ad_table_s1_join_receipt.json"
    h5_samples.to_csv(h5_path, sep="\t", index=False, lineterminator="\n")
    h5_units.to_csv(unit_path, sep="\t", index=False, lineterminator="\n")
    joined.sort_values(["join_status", "sample_id"]).to_csv(join_path, sep="\t", index=False, na_rep="", lineterminator="\n")
    internal_conflicts.to_csv(conflict_path, sep="\t", index=False, lineterminator="\n")
    summary_rows = [
        {"metric": "h5ad_cells", "value": n_obs},
        {"metric": "h5ad_samples", "value": len(h5_samples)},
        {"metric": "h5ad_analysis_units", "value": len(h5_units)},
        {"metric": "pooled_sample_ids", "value": len(pooled_sample_ids)},
        {"metric": "table_samples", "value": len(table)},
        {"metric": "samples_in_both", "value": int(joined.loc[joined["join_status"] == "in_both", "sample_id"].nunique())},
        {"metric": "table_only_samples", "value": len(table_only)},
        {"metric": "h5ad_only_samples", "value": len(h5ad_only)},
        {"metric": "h5ad_internal_conflict_samples", "value": len(internal_conflicts)},
        *[
            {"metric": f"{column}_mismatch_samples", "value": count}
            for column, count in mismatch_counts.items()
        ],
    ]
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False, lineterminator="\n")

    failures: list[str] = []
    limitations: list[str] = []
    if len(internal_conflicts):
        failures.append(f"{len(internal_conflicts)} H5AD samples contain internally conflicting metadata")
    if h5ad_only:
        failures.append(f"{len(h5ad_only)} H5AD samples are absent from Table S1")
    for column in ["patient_id", "dataset", "study_id", "sample_type"]:
        count = mismatch_counts[column]
        if count:
            failures.append(f"{count} joined analysis units disagree on {column}")
    if mismatch_counts["sample_tissue"]:
        limitations.append(
            f"{mismatch_counts['sample_tissue']} joined analysis units disagree on coarse tissue after explicit ontology mapping"
        )
    if pooled_sample_ids:
        limitations.append(
            f"{len(pooled_sample_ids)} sample IDs pool multiple donor IDs; patient identity is retained from H5AD donor_id and Table S1 patient comparison is not applicable"
        )
    receipt: dict[str, object] = {
        "status": "failed" if failures else ("passed_with_limitations" if limitations else "passed"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "h5ad_sha256": sha256_file(h5ad),
        "table_samples": str(table_samples.resolve()),
        "table_samples_sha256": sha256_file(table_samples),
        "n_obs": n_obs,
        "h5ad_samples": len(h5_samples),
        "h5ad_analysis_units": len(h5_units),
        "table_sample_rows": len(table),
        "samples_in_both": int(joined.loc[joined["join_status"] == "in_both", "sample_id"].nunique()),
        "pooled_sample_ids": sorted(pooled_sample_ids),
        "table_only_sample_ids": table_only,
        "h5ad_only_sample_ids": h5ad_only,
        "metadata_mismatch_counts": mismatch_counts,
        "h5ad_field_map": field_map,
        "failures": failures,
        "limitations": limitations,
        "outputs": {},
    }
    for path in (h5_path, unit_path, join_path, conflict_path, summary_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if failures:
        raise ValueError("; ".join(failures))
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--table-samples", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-obs", type=int)
    parser.add_argument("--patient-field", default="donor_id")
    parser.add_argument("--tissue-field", default="tissue")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    field_map = {
        "sample_id": "sample_id",
        "patient_id": args.patient_field,
        "dataset": "dataset",
        "study_id": "study_id",
        "sample_type": "sample_type",
        "sample_tissue": args.tissue_field,
    }
    receipt = build_join(args.h5ad, args.table_samples, args.output_dir, args.expected_obs, field_map)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
