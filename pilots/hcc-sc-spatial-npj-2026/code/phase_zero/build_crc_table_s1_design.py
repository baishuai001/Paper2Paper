#!/usr/bin/env python3
"""Build the CRC Table S1 patient-sample design receipt for P0.0/P0.1."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


PATIENT_SHEET = "C Patient metadata"
SAMPLE_SHEET = "D Sample metadata"
PATIENT_REQUIRED = {
    "patient_id",
    "dataset",
    "study_id",
    "treatment_status_before_resection",
    "treatment_response",
    "RECIST",
}
SAMPLE_REQUIRED = {
    "sample_id",
    "patient_id",
    "dataset",
    "study_id",
    "sample_type",
    "sample_tissue",
}
PATIENT_AUGMENT = [
    "patient_id",
    "platform",
    "tissue_processing_lab",
    "hospital_location",
    "country",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def load_sheet(workbook: Path, sheet: str) -> pd.DataFrame:
    frame = pd.read_excel(
        workbook,
        sheet_name=sheet,
        header=1,
        dtype=object,
        engine="openpyxl",
    )
    frame = frame.dropna(how="all").copy()
    frame.columns = [str(value).strip() for value in frame.columns]
    return frame


def require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}")


def clean_identifier(series: pd.Series, label: str) -> pd.Series:
    if series.isna().any():
        raise ValueError(f"{label} contains missing identifiers")
    cleaned = series.astype(str).str.strip()
    if (cleaned == "").any():
        raise ValueError(f"{label} contains blank identifiers")
    duplicates = cleaned[cleaned.duplicated()].unique().tolist()
    if duplicates:
        raise ValueError(f"{label} contains duplicate identifiers: {duplicates[:10]}")
    return cleaned


def append_counts(rows: list[dict[str, object]], frame: pd.DataFrame, column: str, prefix: str) -> None:
    values = frame[column].fillna("<missing>").astype(str)
    for group, count in values.value_counts(dropna=False).items():
        rows.append({"metric": prefix, "group": group, "value": int(count)})


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, na_rep="", lineterminator="\n")


def build_design(workbook: Path, output_dir: Path, expected_sha256: str = "") -> dict[str, object]:
    if not workbook.is_file() or workbook.stat().st_size == 0:
        raise FileNotFoundError(f"workbook is missing or empty: {workbook}")
    actual_sha256 = sha256_file(workbook)
    if expected_sha256 and actual_sha256 != expected_sha256.upper():
        raise ValueError(
            f"workbook SHA256 mismatch: expected {expected_sha256.upper()}, got {actual_sha256}"
        )

    patients = load_sheet(workbook, PATIENT_SHEET)
    samples = load_sheet(workbook, SAMPLE_SHEET)
    require_columns(patients, PATIENT_REQUIRED, PATIENT_SHEET)
    require_columns(samples, SAMPLE_REQUIRED, SAMPLE_SHEET)
    patients["patient_id"] = clean_identifier(patients["patient_id"], "patient_id")
    samples["sample_id"] = clean_identifier(samples["sample_id"], "sample_id")
    if samples["patient_id"].isna().any() or (samples["patient_id"].astype(str).str.strip() == "").any():
        raise ValueError("sample metadata contains missing patient_id")
    samples["patient_id"] = samples["patient_id"].astype(str).str.strip()

    patient_ids = set(patients["patient_id"])
    sample_patient_ids = set(samples["patient_id"])
    orphan_sample_patients = sorted(sample_patient_ids - patient_ids)
    if orphan_sample_patients:
        raise ValueError(f"sample patient IDs absent from patient table: {orphan_sample_patients[:10]}")
    patients_without_samples = sorted(patient_ids - sample_patient_ids)

    augment = [column for column in PATIENT_AUGMENT if column in patients.columns]
    patient_context = patients[augment].rename(
        columns={column: f"patient_{column}" for column in augment if column != "patient_id"}
    )
    design = samples.merge(patient_context, on="patient_id", how="left", validate="many_to_one")

    output_dir.mkdir(parents=True, exist_ok=True)
    patient_path = output_dir / "P0_table_s1_patients.tsv"
    sample_path = output_dir / "P0_table_s1_samples.tsv"
    design_path = output_dir / "P0_patient_sample_design_table_s1.tsv"
    summary_path = output_dir / "P0_table_s1_summary.tsv"
    receipt_path = output_dir / "P0_table_s1_receipt.json"
    write_tsv(patients, patient_path)
    write_tsv(samples, sample_path)
    write_tsv(design, design_path)

    summary_rows: list[dict[str, object]] = [
        {"metric": "patient_rows", "group": "all", "value": len(patients)},
        {"metric": "sample_rows", "group": "all", "value": len(samples)},
        {"metric": "unique_sample_patient_ids", "group": "all", "value": len(sample_patient_ids)},
        {"metric": "patients_without_samples", "group": "all", "value": len(patients_without_samples)},
        {"metric": "orphan_sample_patient_ids", "group": "all", "value": 0},
        {
            "metric": "response_nonmissing_patients",
            "group": "all",
            "value": int(patients["treatment_response"].notna().sum()),
        },
        {
            "metric": "RECIST_nonmissing_patients",
            "group": "all",
            "value": int(patients["RECIST"].notna().sum()),
        },
    ]
    append_counts(summary_rows, samples, "sample_type", "sample_type_samples")
    append_counts(summary_rows, samples, "sample_tissue", "sample_tissue_samples")
    append_counts(
        summary_rows,
        patients,
        "treatment_status_before_resection",
        "treatment_status_patients",
    )
    for column, metric in (("sample_type", "sample_type_patients"), ("sample_tissue", "sample_tissue_patients")):
        grouped = samples.assign(_group=samples[column].fillna("<missing>").astype(str)).groupby("_group")[
            "patient_id"
        ].nunique()
        for group, count in grouped.sort_values(ascending=False).items():
            summary_rows.append({"metric": metric, "group": group, "value": int(count)})
    summary = pd.DataFrame(summary_rows, columns=["metric", "group", "value"])
    write_tsv(summary, summary_path)

    receipt: dict[str, object] = {
        "status": "passed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "workbook": str(workbook.resolve()),
        "workbook_bytes": workbook.stat().st_size,
        "workbook_sha256": actual_sha256,
        "patient_rows": len(patients),
        "patient_ids_unique": True,
        "sample_rows": len(samples),
        "sample_ids_unique": True,
        "unique_sample_patient_ids": len(sample_patient_ids),
        "orphan_sample_patient_ids": orphan_sample_patients,
        "patients_without_samples_count": len(patients_without_samples),
        "patients_without_samples": patients_without_samples,
        "outputs": {},
    }
    for path in (patient_path, sample_path, design_path, summary_path):
        receipt["outputs"][path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-sha256", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = build_design(args.workbook, args.output_dir, args.expected_sha256)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
