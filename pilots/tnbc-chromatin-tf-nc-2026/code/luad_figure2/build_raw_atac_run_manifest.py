#!/usr/bin/env python3
"""Combine the frozen LUAD PDX and cell-line baseline run manifests."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    if not slug:
        raise RuntimeError(f"Cannot make filesystem-safe slug from {value!r}")
    return slug


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdx_tsv", type=Path)
    parser.add_argument("cell_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    rows: list[dict] = []
    for item in read_tsv(args.pdx_tsv):
        if item.get("figure2_primary_include") != "TRUE":
            continue
        rows.append(
            {
                "system": "PDX",
                "sample_id": item["model_id"],
                "sample_slug": slugify(item["model_id"]),
                "run": item["run"],
                "layout": item["library_layout"],
                "host_depletion": "mouse_then_human",
                "source_accession": "GSE269746",
                "independent_biological_unit": item["model_id"],
                "figure2_primary_include": "TRUE",
            }
        )
    for item in read_tsv(args.cell_tsv):
        rows.append(
            {
                "system": "cell_line",
                "sample_id": item["cell_line"],
                "sample_slug": slugify(item["cell_line"]),
                "run": item["run"],
                "layout": item["library_layout"],
                "host_depletion": "none",
                "source_accession": item["submission"],
                "independent_biological_unit": item["cell_line"],
                "figure2_primary_include": "TRUE",
            }
        )
    rows.sort(key=lambda row: (row["system"], row["sample_id"].lower()))
    if sum(row["system"] == "PDX" for row in rows) != 13:
        raise RuntimeError("Expected 13 LUAD PDX raw runs")
    if sum(row["system"] == "cell_line" for row in rows) != 19:
        raise RuntimeError("Expected 19 strict LUAD cell-line raw runs")
    if len({row["run"] for row in rows}) != 32:
        raise RuntimeError("Run accessions must be unique")
    if len({(row["system"], row["independent_biological_unit"]) for row in rows}) != 32:
        raise RuntimeError("Independent biological units must be unique within system")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "figure2_raw_atac_runs.tsv"
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    receipt = {
        "total_raw_runs": len(rows),
        "pdx_luad_runs": 13,
        "strict_luad_cell_line_runs": 19,
        "pdx_layouts": sorted({row["layout"] for row in rows if row["system"] == "PDX"}),
        "cell_line_layouts": sorted({row["layout"] for row in rows if row["system"] == "cell_line"}),
        "one_run_per_independent_model": True,
        "output": str(output),
    }
    (args.output_dir / "figure2_raw_atac_runs_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
