#!/usr/bin/env python3
"""Freeze one ancestral-named NCI PDMR PDX RNA-seq sample per patient donor."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re


STRICT_LABELS = {"LUAD": "LUAD", "LUSC": "LUSC"}


def lineage_depth(sample: str) -> tuple[int, int, int, str]:
    """Rank PDX names deterministically toward an ancestral unqualified sample.

    NCI PDMR names append two or three alphanumeric characters at each passage.
    The public RNA-seq manifest does not expose a numeric passage column, so
    this name-based rule prevents repeat mice/descendants from being counted as
    independent donors; it is not presented as a measured passage number.
    """

    qualifiers = len(re.findall(r"(?:_RG-|_AL-|POOL|ORIGINATOR)", sample, flags=re.I))
    punctuation = len(re.findall(r"[^A-Za-z0-9]", sample))
    return qualifiers, punctuation, len(sample), sample


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("receipt_json", type=Path)
    args = parser.parse_args()

    rows = json.loads(args.input_json.read_text(encoding="utf-8"))
    eligible = [
        row
        for row in rows
        if row.get("pdm_type") == "PDX"
        and row.get("oncotree") in STRICT_LABELS
        and row.get("rsem_gene_href")
        and row.get("sample", "").upper() != "ORIGINATOR"
    ]
    by_donor: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in eligible:
        key = (str(row["patient"]), str(row["oncotree"]))
        by_donor.setdefault(key, []).append(row)

    selected = []
    for (patient, group), candidates in sorted(by_donor.items()):
        winner = min(candidates, key=lambda row: lineage_depth(str(row["sample"])))
        selected.append(
            {
                "sample_id": f"PDMR_{patient}",
                "patient_id": patient,
                "specimen_id": winner["specimen"],
                "pdmr_sample_id": winner["sample"],
                "group": group,
                "diagnosis": winner["diagnosis"],
                "pdm_type": winner["pdm_type"],
                "rsem_version": winner["rsem_version"],
                "rsem_gene_url": winner["rsem_gene_href"],
                "candidate_samples_for_donor": len(candidates),
                "selection_rule": "deterministic_shortest_unqualified_pdx_name",
            }
        )

    counts = {group: sum(row["group"] == group for row in selected) for group in STRICT_LABELS}
    if min(counts.values()) < 10:
        raise SystemExit(f"PDX donor group below frozen n=10 minimum: {counts}")
    if len({row["patient_id"] for row in selected}) != len(selected):
        raise SystemExit("A PDMR donor was selected more than once")

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=selected[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    receipt = {
        "input_rows": len(rows),
        "eligible_pdx_rsem_rows": len(eligible),
        "unique_selected_donors": len(selected),
        "group_counts": counts,
        "selection_unit": "unique NCI PDMR patient donor",
        "selection_rule": "one deterministic shortest unqualified PDX name per donor; no originator, organoid, culture, or repeat mouse; numeric passage unavailable",
    }
    args.receipt_json.parent.mkdir(parents=True, exist_ok=True)
    args.receipt_json.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
