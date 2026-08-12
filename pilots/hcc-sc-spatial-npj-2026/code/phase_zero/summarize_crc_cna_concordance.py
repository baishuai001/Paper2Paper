#!/usr/bin/env python3
"""Compare CopyKAT and SCEVAN calls without treating either as DNA truth."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


VALID_COPYKAT = {"aneuploid", "diploid", "not.defined"}
VALID_SCEVAN = {"tumor", "normal", "filtered"}
TRUE_STRINGS = {"true", "1", "yes"}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def _read_required(path: Path, columns: set[str], label: str) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} evidence missing or empty: {path}")
    table = pd.read_csv(path, sep="\t", dtype=str)
    missing = columns - set(table.columns)
    if missing:
        raise ValueError(f"{label} fields missing: {sorted(missing)}")
    return table


def _evidence_class(copykat_call: str, scevan_call: str) -> str:
    lookup = {
        ("aneuploid", "tumor"): "concordant_malignancy_support",
        ("aneuploid", "normal"): "discordant_aneuploid_vs_scevan_normal",
        ("aneuploid", "filtered"): "copykat_only_support_scevan_filtered",
        ("diploid", "tumor"): "scevan_only_support_copykat_diploid_nonexclusive",
        ("diploid", "normal"): "concordant_nonaneuploid_classification_not_nonmalignancy_proof",
        ("diploid", "filtered"): "unresolved_scevan_filtered_copykat_diploid_nonexclusive",
        ("not.defined", "tumor"): "scevan_only_support_copykat_unresolved",
        ("not.defined", "normal"): "unresolved_copykat_not_defined",
        ("not.defined", "filtered"): "unresolved_both_methods",
    }
    try:
        return lookup[(copykat_call, scevan_call)]
    except KeyError as exc:
        raise ValueError(f"unsupported CNA call pair: {copykat_call!r}, {scevan_call!r}") from exc


def build_cell_evidence(
    metadata_path: Path,
    copykat_path: Path,
    scevan_path: Path,
) -> pd.DataFrame:
    """Identity-join two CNA methods and retain every cell-level evidence state."""
    metadata = _read_required(
        metadata_path,
        {"cell_id", "cna_input_role", "is_known_normal_reference"},
        "CNA metadata",
    )
    copykat = _read_required(copykat_path, {"cell.names", "copykat.pred"}, "CopyKAT")
    scevan = _read_required(scevan_path, {"cell.names", "class"}, "SCEVAN")
    for table, column, label in (
        (metadata, "cell_id", "metadata"),
        (copykat, "cell.names", "CopyKAT"),
        (scevan, "cell.names", "SCEVAN"),
    ):
        if table[column].isna().any() or table[column].duplicated().any():
            raise ValueError(f"{label} cell IDs are missing or duplicated")
    unexpected_copykat = sorted(set(copykat["copykat.pred"].dropna()) - VALID_COPYKAT)
    unexpected_scevan = sorted(set(scevan["class"].dropna()) - VALID_SCEVAN)
    if unexpected_copykat:
        raise ValueError(f"unexpected CopyKAT calls: {unexpected_copykat}")
    if unexpected_scevan:
        raise ValueError(f"unexpected SCEVAN calls: {unexpected_scevan}")

    reference_flag = metadata["is_known_normal_reference"].str.lower().isin(TRUE_STRINGS)
    role_is_reference = metadata["cna_input_role"].eq("known_normal_reference")
    if not reference_flag.equals(role_is_reference):
        raise ValueError("known-normal reference flag and CNA input role disagree")

    joined = metadata.merge(
        copykat[["cell.names", "copykat.pred"]],
        left_on="cell_id",
        right_on="cell.names",
        how="outer",
        validate="one_to_one",
        indicator="copykat_join",
    )
    if not joined["copykat_join"].eq("both").all():
        raise ValueError(f"CopyKAT identity join is incomplete: {joined['copykat_join'].value_counts().to_dict()}")
    joined = joined.drop(columns=["cell.names", "copykat_join"]).merge(
        scevan[["cell.names", "class"]].rename(columns={"class": "scevan_class"}),
        left_on="cell_id",
        right_on="cell.names",
        how="outer",
        validate="one_to_one",
        indicator="scevan_join",
    )
    if not joined["scevan_join"].eq("both").all():
        raise ValueError(f"SCEVAN identity join is incomplete: {joined['scevan_join'].value_counts().to_dict()}")
    joined = joined.drop(columns=["cell.names", "scevan_join"])
    joined["evidence_class"] = [
        _evidence_class(copykat_call, scevan_call)
        for copykat_call, scevan_call in zip(joined["copykat.pred"], joined["scevan_class"])
    ]
    joined["copykat_malignancy_support"] = joined["copykat.pred"].eq("aneuploid")
    joined["scevan_malignancy_support"] = joined["scevan_class"].eq("tumor")
    joined["either_malignancy_support"] = (
        joined["copykat_malignancy_support"] | joined["scevan_malignancy_support"]
    )
    joined["method_unresolved"] = (
        joined["copykat.pred"].eq("not.defined") | joined["scevan_class"].eq("filtered")
    )
    joined["both_defined"] = (
        ~joined["copykat.pred"].eq("not.defined") & ~joined["scevan_class"].eq("filtered")
    )
    joined["defined_binary_agreement"] = joined["both_defined"] & (
        joined["copykat_malignancy_support"] == joined["scevan_malignancy_support"]
    )
    return joined


def summarize(
    metadata_path: Path,
    copykat_path: Path,
    scevan_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    joined = build_cell_evidence(metadata_path, copykat_path, scevan_path)

    detail = (
        joined.groupby(
            ["cna_input_role", "copykat.pred", "scevan_class", "evidence_class"],
            observed=True,
        )
        .size()
        .rename("cells")
        .reset_index()
    )
    detail["fraction_within_input_role"] = detail["cells"] / detail.groupby("cna_input_role")[
        "cells"
    ].transform("sum")

    role_rows: list[dict[str, object]] = []
    for role, group in joined.groupby("cna_input_role", observed=True):
        defined = group[group["both_defined"]]
        role_rows.append(
            {
                "cna_input_role": role,
                "cells": len(group),
                "copykat_aneuploid": int(group["copykat_malignancy_support"].sum()),
                "scevan_tumor": int(group["scevan_malignancy_support"].sum()),
                "both_malignancy_support": int(
                    (group["copykat_malignancy_support"] & group["scevan_malignancy_support"]).sum()
                ),
                "either_malignancy_support": int(group["either_malignancy_support"].sum()),
                "method_unresolved": int(group["method_unresolved"].sum()),
                "both_defined": int(group["both_defined"].sum()),
                "defined_binary_agreement_fraction": (
                    float(defined["defined_binary_agreement"].mean()) if len(defined) else None
                ),
                "aneuploid_vs_scevan_normal": int(
                    (
                        group["copykat.pred"].eq("aneuploid")
                        & group["scevan_class"].eq("normal")
                    ).sum()
                ),
                "diploid_vs_scevan_tumor": int(
                    (
                        group["copykat.pred"].eq("diploid")
                        & group["scevan_class"].eq("tumor")
                    ).sum()
                ),
            }
        )
    role_summary = pd.DataFrame(role_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    detail_path = output_dir / "P0_cna_method_pair_by_role.tsv"
    role_path = output_dir / "P0_cna_role_summary.tsv"
    detail.to_csv(detail_path, sep="\t", index=False, lineterminator="\n")
    role_summary.to_csv(role_path, sep="\t", index=False, lineterminator="\n")

    author = joined[joined["cna_input_role"].eq("author_cancer")]
    references = joined[joined["cna_input_role"].eq("known_normal_reference")]
    receipt: dict[str, object] = {
        "status": "diagnostic_one_sample_two_method_cna_completed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "metadata_sha256": sha256_file(metadata_path),
        "copykat_prediction_sha256": sha256_file(copykat_path),
        "scevan_prediction_sha256": sha256_file(scevan_path),
        "cells": int(len(joined)),
        "author_cancer_cells": int(len(author)),
        "author_cancer_both_malignancy_support": int(
            (author["copykat_malignancy_support"] & author["scevan_malignancy_support"]).sum()
        ),
        "author_cancer_either_malignancy_support": int(author["either_malignancy_support"].sum()),
        "author_cancer_method_unresolved": int(author["method_unresolved"].sum()),
        "known_normal_reference_either_malignancy_support": int(
            references["either_malignancy_support"].sum()
        ),
        "claim_boundary": (
            "This is a one-sample diagnostic comparison using author labels and predefined normal "
            "references, not DNA truth. Aneuploid/tumor calls support malignancy; diploid/normal "
            "calls do not prove nonmalignancy; filtered/not.defined remain unresolved. The result "
            "cannot qualify malignant cells across the CRC atlas or freeze the cNMF discovery set."
        ),
        "outputs": {
            detail_path.name: {"bytes": detail_path.stat().st_size, "sha256": sha256_file(detail_path)},
            role_path.name: {"bytes": role_path.stat().st_size, "sha256": sha256_file(role_path)},
        },
    }
    (output_dir / "P0_cna_concordance_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--copykat", required=True, type=Path)
    parser.add_argument("--scevan", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(
        json.dumps(
            summarize(args.metadata, args.copykat, args.scevan, args.output_dir),
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
