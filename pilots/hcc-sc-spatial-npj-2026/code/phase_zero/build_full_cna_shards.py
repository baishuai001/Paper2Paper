#!/usr/bin/env python3
"""Freeze allowed full-CNA samples into deterministic workload-balanced shards."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ALLOWED_REFERENCE_REVIEWS = {"pass_le_5pct", "warning_gt_5pct_le_10pct"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    return fields, rows


def freeze_shards(
    panel_path: Path,
    review_path: Path,
    output_dir: Path,
    shard_count: int = 4,
) -> dict[str, object]:
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    fields, panel = read_tsv(panel_path)
    _, review = read_tsv(review_path)
    required_panel = {"panel_order", "sample_id", "analysis_patient_id", "cancer_cells"}
    required_review = {"panel_order", "sample_id", "reference_review"}
    if required_panel - set(fields):
        raise ValueError(f"panel missing fields: {sorted(required_panel - set(fields))}")
    if not review or required_review - set(review[0]):
        raise ValueError("sample review is empty or missing required fields")

    review_by_order = {int(row["panel_order"]): row for row in review}
    if len(review_by_order) != len(review):
        raise ValueError("sample review contains duplicate panel_order values")
    allowed: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    for row in panel:
        order = int(row["panel_order"])
        reviewed = review_by_order.get(order)
        if reviewed is None or reviewed["sample_id"] != row["sample_id"]:
            raise ValueError(f"panel/review identity mismatch at panel_order {order}")
        if reviewed["reference_review"] in ALLOWED_REFERENCE_REVIEWS:
            allowed.append(row)
        else:
            excluded.append(row)
    if set(review_by_order) != {int(row["panel_order"]) for row in panel}:
        raise ValueError("panel and review have different panel_order values")

    bins: list[list[dict[str, str]]] = [[] for _ in range(shard_count)]
    loads = [0] * shard_count
    for row in sorted(
        allowed,
        key=lambda value: (-int(value["cancer_cells"]), int(value["panel_order"])),
    ):
        target = min(range(shard_count), key=lambda index: (loads[index], index))
        bins[target].append(row)
        loads[target] += int(row["cancer_cells"])

    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite frozen full-CNA manifest: {output_dir}")
    output_dir.mkdir(parents=True)
    allowed_path = output_dir / "P0_paired_cna_full_allowed.tsv"
    with allowed_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(allowed, key=lambda row: int(row["panel_order"])))

    shard_rows: list[dict[str, object]] = []
    observed_orders: list[int] = []
    for index, rows in enumerate(bins, start=1):
        path = output_dir / f"P0_paired_cna_full_shard_{index:02d}.tsv"
        rows = sorted(rows, key=lambda row: int(row["panel_order"]))
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        observed_orders.extend(int(row["panel_order"]) for row in rows)
        shard_rows.append(
            {
                "shard": index,
                "samples": len(rows),
                "author_cancer_cells": sum(int(row["cancer_cells"]) for row in rows),
                "panel_orders": [int(row["panel_order"]) for row in rows],
                "path": str(path.resolve()),
                "sha256": sha256_file(path),
            }
        )
    if sorted(observed_orders) != sorted(int(row["panel_order"]) for row in allowed):
        raise RuntimeError("full-CNA shards do not partition allowed samples exactly once")

    receipt: dict[str, object] = {
        "status": "full_cna_allowed_shards_frozen",
        "panel": str(panel_path.resolve()),
        "panel_sha256": sha256_file(panel_path),
        "screen_review": str(review_path.resolve()),
        "screen_review_sha256": sha256_file(review_path),
        "allowed_reference_reviews": sorted(ALLOWED_REFERENCE_REVIEWS),
        "samples_total": len(panel),
        "samples_allowed": len(allowed),
        "samples_excluded": len(excluded),
        "author_cancer_cells_allowed": sum(int(row["cancer_cells"]) for row in allowed),
        "shards": shard_rows,
        "formal_cnmf_started": False,
    }
    (output_dir / "P0_paired_cna_full_shards_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--sample-review", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--shards", type=int, default=4)
    args = parser.parse_args()
    print(
        json.dumps(
            freeze_shards(args.panel, args.sample_review, args.output_dir, args.shards),
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
