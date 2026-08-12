#!/usr/bin/env python3
"""Join CopyKAT calls to frozen input roles and write aggregate evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


VALID_CALLS = {"aneuploid", "diploid", "not.defined"}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def summarize(metadata_path: Path, prediction_path: Path, output_dir: Path) -> dict[str, object]:
    for path in (metadata_path, prediction_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"CopyKAT evidence missing or empty: {path}")
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    predictions = pd.read_csv(prediction_path, sep="\t", dtype=str)
    required_metadata = {"cell_id", "cna_input_role", "is_known_normal_reference"}
    required_predictions = {"cell.names", "copykat.pred"}
    if not required_metadata.issubset(metadata.columns):
        raise ValueError(f"metadata fields missing: {sorted(required_metadata - set(metadata.columns))}")
    if not required_predictions.issubset(predictions.columns):
        raise ValueError(f"prediction fields missing: {sorted(required_predictions - set(predictions.columns))}")
    if metadata["cell_id"].duplicated().any() or predictions["cell.names"].duplicated().any():
        raise ValueError("CopyKAT metadata or predictions contain duplicate cell IDs")
    unexpected = sorted(set(predictions["copykat.pred"].dropna()) - VALID_CALLS)
    if unexpected:
        raise ValueError(f"unexpected CopyKAT calls: {unexpected}")
    joined = metadata.merge(
        predictions,
        left_on="cell_id",
        right_on="cell.names",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    if not joined["_merge"].eq("both").all():
        counts = joined["_merge"].value_counts().to_dict()
        raise ValueError(f"CopyKAT cell identity join is incomplete: {counts}")
    joined["cna_evidence_class"] = joined["copykat.pred"].map(
        {
            "aneuploid": "aneuploid_supports_malignancy",
            "diploid": "diploid_does_not_exclude_malignancy",
            "not.defined": "unresolved",
        }
    )
    summary = (
        joined.groupby(["cna_input_role", "copykat.pred", "cna_evidence_class"], observed=True)
        .size()
        .rename("cells")
        .reset_index()
    )
    role_totals = summary.groupby("cna_input_role")["cells"].transform("sum")
    summary["fraction_within_input_role"] = summary["cells"] / role_totals
    output_dir.mkdir(parents=True, exist_ok=True)
    source_path = output_dir / "P0_copykat_role_prediction.tsv"
    summary.to_csv(source_path, sep="\t", index=False, lineterminator="\n")
    author = joined[joined["cna_input_role"] == "author_cancer"]
    normals = joined[joined["cna_input_role"] == "known_normal_reference"]
    receipt = {
        "status": "diagnostic_one_sample_copykat_completed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "metadata_sha256": sha256_file(metadata_path),
        "prediction_sha256": sha256_file(prediction_path),
        "cells": int(len(joined)),
        "author_cancer_cells": int(len(author)),
        "author_cancer_aneuploid": int(author["copykat.pred"].eq("aneuploid").sum()),
        "author_cancer_diploid": int(author["copykat.pred"].eq("diploid").sum()),
        "author_cancer_unresolved": int(author["copykat.pred"].eq("not.defined").sum()),
        "known_normal_cells": int(len(normals)),
        "known_normal_aneuploid": int(normals["copykat.pred"].eq("aneuploid").sum()),
        "claim_boundary": "Aneuploid calls support malignancy, but diploid calls do not prove nonmalignancy. This is one-sample CopyKAT evidence without SCEVAN sensitivity or DNA truth and cannot qualify the full atlas malignant set.",
        "outputs": {
            source_path.name: {"bytes": source_path.stat().st_size, "sha256": sha256_file(source_path)}
        },
    }
    (output_dir / "P0_copykat_summary_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--prediction", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(json.dumps(summarize(args.metadata, args.prediction, args.output_dir), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
