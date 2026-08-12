#!/usr/bin/env python3
"""Audit post-QC metrics and doublet evidence in the published CRC atlas object."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from join_crc_h5ad_table_s1 import canonical_patient, canonical_tissue, clean_value, read_obs_column


IDENTITY_COLUMNS = ["sample_id", "patient_id", "dataset", "sample_type", "sample_tissue"]
COUNT_CANDIDATES = ["n_counts", "total_counts"]
GENE_CANDIDATES = ["n_genes_by_counts", "n_genes"]
MITO_CANDIDATES = ["pct_counts_mito"]
DOUBLET_COLUMNS = ["SOLO_doublet_status", "SOLO_doublet_prob", "SOLO_singlet_prob"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def first_available(obs: h5py.Group, candidates: list[str], label: str) -> str:
    for candidate in candidates:
        if candidate in obs:
            return candidate
    raise ValueError(f"H5AD obs has no supported {label} field; tried {candidates}")


def numeric(values: np.ndarray, label: str) -> np.ndarray:
    converted = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    if np.isfinite(converted).sum() == 0:
        raise ValueError(f"{label} has no numeric values")
    return converted


def safe_fraction(mask: pd.Series | np.ndarray) -> float:
    array = np.asarray(mask, dtype=bool)
    return float(array.mean()) if array.size else float("nan")


def summarize_group(group: pd.DataFrame, count_field: str, gene_field: str, mito_field: str) -> pd.Series:
    status = group["SOLO_doublet_status_clean"]
    doublet_prob = group["SOLO_doublet_prob_numeric"]
    singlet_prob = group["SOLO_singlet_prob_numeric"]
    probability_available = np.isfinite(doublet_prob) | np.isfinite(singlet_prob)
    status_available = status != ""
    if probability_available.any():
        evidence = "probability_available"
    elif status_available.any():
        evidence = "status_without_probability"
    else:
        evidence = "not_available"
    return pd.Series(
        {
            "cell_count": len(group),
            "dataset": "|".join(sorted(set(group["dataset_clean"]) - {""})),
            "sample_type": "|".join(sorted(set(group["sample_type_clean"]) - {""})),
            "sample_tissue": "|".join(sorted(set(group["sample_tissue_clean"]) - {""})),
            "count_field": count_field,
            "count_median": float(np.nanmedian(group[count_field])),
            "count_min": float(np.nanmin(group[count_field])),
            "count_below_or_equal_400": int(np.sum(group[count_field] <= 400)),
            "gene_field": gene_field,
            "genes_median": float(np.nanmedian(group[gene_field])),
            "genes_min": float(np.nanmin(group[gene_field])),
            "genes_below_or_equal_100": int(np.sum(group[gene_field] <= 100)),
            "mito_field": mito_field,
            "mito_median": float(np.nanmedian(group[mito_field])),
            "mito_max": float(np.nanmax(group[mito_field])),
            "mito_greater_or_equal_50": int(np.sum(group[mito_field] >= 50)),
            "doublet_evidence": evidence,
            "doublet_status_missing_fraction": safe_fraction(~status_available),
            "doublet_probability_missing_fraction": safe_fraction(~probability_available),
            "doublet_label_count": int(status.str.lower().str.contains("doublet", na=False).sum()),
            "singlet_label_count": int(status.str.lower().str.contains("singlet", na=False).sum()),
        }
    )


def audit_qc(
    h5ad: Path,
    output_dir: Path,
    expected_obs: int | None = None,
    identity_field_map: dict[str, str] | None = None,
) -> dict[str, object]:
    identity_field_map = identity_field_map or {column: column for column in IDENTITY_COLUMNS}
    missing_map = sorted(set(IDENTITY_COLUMNS) - set(identity_field_map))
    if missing_map:
        raise ValueError(f"identity field map missing logical fields: {missing_map}")
    if not h5ad.is_file() or h5ad.stat().st_size == 0:
        raise FileNotFoundError(f"H5AD is missing or empty: {h5ad}")
    with h5py.File(h5ad, "r") as handle:
        if "obs" not in handle or not isinstance(handle["obs"], h5py.Group):
            raise ValueError("H5AD missing obs dataframe")
        obs = handle["obs"]
        for logical, actual in identity_field_map.items():
            if actual not in obs:
                raise ValueError(f"H5AD obs missing identity column for {logical}: {actual}")
        count_field = first_available(obs, COUNT_CANDIDATES, "count")
        gene_field = first_available(obs, GENE_CANDIDATES, "detected-gene")
        mito_field = first_available(obs, MITO_CANDIDATES, "mitochondrial percentage")
        data: dict[str, np.ndarray] = {
            logical: read_obs_column(obs, actual) for logical, actual in identity_field_map.items()
        }
        data[count_field] = numeric(read_obs_column(obs, count_field), count_field)
        data[gene_field] = numeric(read_obs_column(obs, gene_field), gene_field)
        data[mito_field] = numeric(read_obs_column(obs, mito_field), mito_field)
        for column in DOUBLET_COLUMNS:
            data[column] = read_obs_column(obs, column) if column in obs else np.full(len(data["sample_id"]), None)

    lengths = {len(values) for values in data.values()}
    if len(lengths) != 1:
        raise ValueError(f"obs columns have inconsistent lengths: {lengths}")
    n_obs = lengths.pop()
    if expected_obs is not None and n_obs != expected_obs:
        raise ValueError(f"H5AD cell count mismatch: expected {expected_obs}, got {n_obs}")
    frame = pd.DataFrame(data)
    for column in IDENTITY_COLUMNS:
        frame[f"{column}_clean"] = frame[column].map(clean_value)
    frame["patient_id_clean"] = frame["patient_id_clean"].map(canonical_patient)
    frame["sample_tissue_clean"] = frame["sample_tissue_clean"].map(canonical_tissue)
    if (frame["sample_id_clean"] == "").any():
        raise ValueError("H5AD contains cells without sample_id")
    frame["SOLO_doublet_status_clean"] = frame["SOLO_doublet_status"].map(clean_value)
    frame["SOLO_doublet_prob_numeric"] = pd.to_numeric(frame["SOLO_doublet_prob"], errors="coerce")
    frame["SOLO_singlet_prob_numeric"] = pd.to_numeric(frame["SOLO_singlet_prob"], errors="coerce")

    sample_summary = (
        frame.groupby(["sample_id_clean", "patient_id_clean"], sort=True, observed=True)
        .apply(summarize_group, count_field=count_field, gene_field=gene_field, mito_field=mito_field, include_groups=False)
        .reset_index()
        .rename(columns={"sample_id_clean": "sample_id", "patient_id_clean": "analysis_patient_id"})
    )
    sample_summary["patient_id"] = sample_summary["analysis_patient_id"]
    identity_conflicts = sample_summary[
        sample_summary[["dataset", "sample_type", "sample_tissue"]].astype(str).apply(
            lambda column: column.str.contains("\\|", regex=True)
        ).any(axis=1)
    ]
    if len(identity_conflicts):
        raise ValueError(f"{len(identity_conflicts)} samples contain internally conflicting identity metadata")

    dataset_summary = (
        sample_summary.groupby("dataset", sort=True, observed=True)
        .agg(
            samples=("sample_id", "size"),
            cells=("cell_count", "sum"),
            samples_probability_available=("doublet_evidence", lambda values: int((values == "probability_available").sum())),
            samples_status_without_probability=("doublet_evidence", lambda values: int((values == "status_without_probability").sum())),
            samples_doublet_not_available=("doublet_evidence", lambda values: int((values == "not_available").sum())),
            cells_below_count_threshold=("count_below_or_equal_400", "sum"),
            cells_below_gene_threshold=("genes_below_or_equal_100", "sum"),
            cells_above_mito_threshold=("mito_greater_or_equal_50", "sum"),
        )
        .reset_index()
    )
    doublet_receipt = sample_summary[
        [
            "sample_id",
            "analysis_patient_id",
            "patient_id",
            "dataset",
            "cell_count",
            "doublet_evidence",
            "doublet_status_missing_fraction",
            "doublet_probability_missing_fraction",
            "doublet_label_count",
            "singlet_label_count",
        ]
    ].copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    sample_path = output_dir / "P0_qc_by_sample.tsv"
    dataset_path = output_dir / "P0_qc_by_dataset.tsv"
    doublet_path = output_dir / "P0_doublet_receipt.tsv"
    receipt_path = output_dir / "P0_qc_receipt.json"
    sample_summary.to_csv(sample_path, sep="\t", index=False, lineterminator="\n")
    dataset_summary.to_csv(dataset_path, sep="\t", index=False, lineterminator="\n")
    doublet_receipt.to_csv(doublet_path, sep="\t", index=False, lineterminator="\n")

    unsupported_samples = int((sample_summary["doublet_evidence"] != "probability_available").sum())
    receipt: dict[str, object] = {
        "status": "passed_with_limitations" if unsupported_samples else "passed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "h5ad": str(h5ad.resolve()),
        "h5ad_sha256": sha256_file(h5ad),
        "n_obs": n_obs,
        "samples": int(sample_summary["sample_id"].nunique()),
        "sample_patient_analysis_units": len(sample_summary),
        "datasets": int(sample_summary["dataset"].nunique()),
        "selected_fields": {"counts": count_field, "genes": gene_field, "mitochondrial": mito_field},
        "identity_field_map": identity_field_map,
        "published_threshold_violations": {
            "count_below_or_equal_400": int(sample_summary["count_below_or_equal_400"].sum()),
            "genes_below_or_equal_100": int(sample_summary["genes_below_or_equal_100"].sum()),
            "mito_greater_or_equal_50": int(sample_summary["mito_greater_or_equal_50"].sum()),
        },
        "samples_without_probability_evidence": unsupported_samples,
        "claim_boundary": "Post-QC audit only; it cannot reconstruct removed cells or prove ambient-RNA/doublet processing for inputs that were not retained.",
        "outputs": {},
    }
    for path in (sample_path, dataset_path, doublet_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-obs", type=int)
    parser.add_argument("--patient-field", default="donor_id")
    parser.add_argument("--tissue-field", default="tissue")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    identity_field_map = {
        "sample_id": "sample_id",
        "patient_id": args.patient_field,
        "dataset": "dataset",
        "sample_type": "sample_type",
        "sample_tissue": args.tissue_field,
    }
    receipt = audit_qc(args.h5ad, args.output_dir, args.expected_obs, identity_field_map)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
