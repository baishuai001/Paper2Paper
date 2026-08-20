#!/usr/bin/env python3
"""Fetch TCGA-LUAD smoking status from the official GDC cases API."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import requests


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_root", type=Path)
    args = parser.parse_args()
    output = args.output_root.resolve() / "data"
    output.mkdir(parents=True, exist_ok=True)
    raw_path = output / "gdc_tcga_luad_smoking_cases.json"
    table_path = output / "tcga_luad_smoking.tsv"

    filters = {"op": "in", "content": {"field": "project.project_id", "value": ["TCGA-LUAD"]}}
    parameters = {
        "filters": json.dumps(filters, separators=(",", ":")),
        "fields": "submitter_id,exposures.tobacco_smoking_status,exposures.pack_years_smoked",
        "format": "JSON",
        "size": "600",
    }
    response = requests.get("https://api.gdc.cancer.gov/cases", params=parameters, timeout=120)
    response.raise_for_status()
    raw_path.write_bytes(response.content)
    payload = response.json()
    hits = payload["data"]["hits"]
    if len(hits) < 500:
        raise RuntimeError(f"Expected at least 500 TCGA-LUAD cases, received {len(hits)}")

    rows: list[dict[str, object]] = []
    for hit in hits:
        exposures = hit.get("exposures") or []
        statuses = [str(item.get("tobacco_smoking_status", "")).strip() for item in exposures]
        statuses = [value for value in statuses if value and value.lower() not in {"not reported", "unknown"}]
        status = statuses[0] if statuses else ""
        if status:
            smoking_model = "Never" if "non-smoker" in status.lower() else "Ever"
        else:
            smoking_model = ""
        packs = [item.get("pack_years_smoked") for item in exposures if item.get("pack_years_smoked") is not None]
        rows.append({
            "patient_id": hit["submitter_id"],
            "tobacco_smoking_status": status,
            "smoking_model": smoking_model,
            "gdc_pack_years_smoked": packs[0] if packs else "",
        })
    if len({row["patient_id"] for row in rows}) != len(rows):
        raise RuntimeError("GDC returned duplicate TCGA-LUAD patient IDs")
    rows.sort(key=lambda row: str(row["patient_id"]))
    with table_path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    receipt = {
        "status": "passed",
        "endpoint": "https://api.gdc.cancer.gov/cases",
        "project": "TCGA-LUAD",
        "cases": len(rows),
        "smoking_available": sum(bool(row["smoking_model"]) for row in rows),
        "never": sum(row["smoking_model"] == "Never" for row in rows),
        "ever": sum(row["smoking_model"] == "Ever" for row in rows),
        "raw_sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
    }
    (output / "tcga_luad_smoking_receipt.json").write_text(
        json.dumps(receipt, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

