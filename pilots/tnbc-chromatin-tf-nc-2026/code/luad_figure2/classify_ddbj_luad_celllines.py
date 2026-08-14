#!/usr/bin/env python3
"""Freeze the histology-qualified DDBJ cell-line set for LUAD Figure 2.

The DDBJ perturbation paper calls the 23-model panel lung cancer cells.  That
does not make every model a lung adenocarcinoma.  This script joins the frozen
DepMap 22Q2 annotation used in Figure 1 and records explicit Cellosaurus
exceptions for models absent from that matrix.  Cellosaurus disease takes
precedence over a broader/legacy DepMap subtype when they disagree.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


SPECIAL = {
    "H1299": {
        "cellosaurus_id": "CVCL_0060",
        "cellosaurus_disease": "Lung large cell carcinoma",
        "cellosaurus_issue": "",
        "evidence_url": "https://www.cellosaurus.org/CVCL_0060",
        "group": "NSCLC_large_cell",
    },
    "II18": {
        "cellosaurus_id": "CVCL_6659",
        "cellosaurus_disease": "Lung adenocarcinoma",
        "cellosaurus_issue": "",
        "evidence_url": "https://www.cellosaurus.org/CVCL_6659",
        "group": "LUAD",
    },
    "LC2AD": {
        "cellosaurus_id": "CVCL_1373",
        "cellosaurus_disease": "Lung adenocarcinoma",
        "cellosaurus_issue": "",
        "evidence_url": "https://www.cellosaurus.org/CVCL_1373",
        "group": "LUAD",
    },
    "RERFLCOK": {
        "cellosaurus_id": "CVCL_3154",
        "cellosaurus_disease": "Astrocytoma",
        "cellosaurus_issue": "Contaminated; Marcus derivative; ICLAC-00344",
        "evidence_url": "https://www.cellosaurus.org/CVCL_3154",
        "group": "MISIDENTIFIED",
    },
    "VMRCLCD": {
        "cellosaurus_id": "CVCL_1787",
        "cellosaurus_disease": "Lung adenocarcinoma",
        "cellosaurus_issue": "",
        "evidence_url": "https://www.cellosaurus.org/CVCL_1787",
        "group": "LUAD",
    },
}


def normalize(name: str) -> str:
    value = re.sub(r"[^A-Za-z0-9]", "", name).upper()
    if value.startswith("NCIH"):
        value = value[3:]
    return value


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("identity_tsv", type=Path)
    parser.add_argument("depmap_manifest_tsv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    identity = read_tsv(args.identity_tsv)
    canonical = [row for row in identity if row.get("canonical_submission") == "TRUE"]
    if len(canonical) != 23 or len({row["cell_line"] for row in canonical}) != 23:
        raise RuntimeError("Expected exactly 23 unique canonical Dataset 2 cell lines")

    depmap = read_tsv(args.depmap_manifest_tsv)
    depmap_by_name: dict[str, dict] = {}
    for row in depmap:
        for field in ("cell_line_name", "stripped_cell_line_name", "cell_line"):
            if row.get(field):
                depmap_by_name.setdefault(normalize(row[field]), row)

    rows: list[dict] = []
    unresolved: list[str] = []
    for item in sorted(canonical, key=lambda row: row["cell_line"].lower()):
        cell_line = item["cell_line"]
        key = normalize(cell_line)
        dep = depmap_by_name.get(key)
        source = "DepMap 22Q2 model annotation + embedded Cellosaurus fields"
        evidence_url = "https://depmap.org/portal/"
        if dep is not None:
            disease = dep.get("Cellosaurus_NCIt_disease", "")
            issue = dep.get("Cellosaurus_issues", "")
            cellosaurus_id = dep.get("RRID", "")
            group = dep.get("group", "")
            depmap_id = dep.get("depmap_id", "")
            depmap_subtype = dep.get("subtype_disease", "")
            if cellosaurus_id:
                evidence_url = f"https://www.cellosaurus.org/{cellosaurus_id}"
        elif key in SPECIAL:
            special = SPECIAL[key]
            disease = special["cellosaurus_disease"]
            issue = special["cellosaurus_issue"]
            cellosaurus_id = special["cellosaurus_id"]
            group = special["group"]
            depmap_id = ""
            depmap_subtype = ""
            source = "Cellosaurus direct record"
            evidence_url = special["evidence_url"]
        else:
            unresolved.append(cell_line)
            continue

        disease_lower = disease.lower()
        is_luad = "lung adenocarcinoma" in disease_lower
        is_large_cell = "lung large cell carcinoma" in disease_lower
        is_lusc = "lung squamous cell carcinoma" in disease_lower or group == "LUSC"
        is_misidentified = group == "MISIDENTIFIED" or bool(issue and "contamin" in issue.lower())
        primary = is_luad and not is_misidentified
        broad_nsclc = (is_luad or is_large_cell) and not is_lusc and not is_misidentified

        if primary:
            decision = "include_strict_LUAD"
        elif is_misidentified:
            decision = "exclude_misidentified"
        elif is_lusc:
            decision = "exclude_LUSC"
        elif is_large_cell:
            decision = "exclude_primary_include_broad_NSCLC_sensitivity"
        else:
            decision = "exclude_non_LUAD"

        rows.append(
            {
                "cell_line": cell_line,
                "canonical_atac_submission": item["submission"],
                "depmap_id": depmap_id,
                "depmap_subtype": depmap_subtype,
                "cellosaurus_id": cellosaurus_id,
                "cellosaurus_disease": disease,
                "cellosaurus_issue": issue,
                "classification_source": source,
                "evidence_url": evidence_url,
                "figure2_primary_include": str(primary).upper(),
                "broad_nsclc_sensitivity_include": str(broad_nsclc).upper(),
                "author_panel_reconstruction_include": "TRUE",
                "decision": decision,
            }
        )

    if unresolved:
        raise RuntimeError(f"Unresolved cell-line identities: {unresolved}")

    columns = [
        "cell_line",
        "canonical_atac_submission",
        "depmap_id",
        "depmap_subtype",
        "cellosaurus_id",
        "cellosaurus_disease",
        "cellosaurus_issue",
        "classification_source",
        "evidence_url",
        "figure2_primary_include",
        "broad_nsclc_sensitivity_include",
        "author_panel_reconstruction_include",
        "decision",
    ]
    output = args.output_dir / "ddbj_cell_line_histology_audit.tsv"
    write_tsv(output, rows, columns)

    primary_n = sum(row["figure2_primary_include"] == "TRUE" for row in rows)
    broad_n = sum(row["broad_nsclc_sensitivity_include"] == "TRUE" for row in rows)
    decisions: dict[str, int] = {}
    for row in rows:
        decisions[row["decision"]] = decisions.get(row["decision"], 0) + 1
    if primary_n != 19:
        raise RuntimeError(f"Frozen strict LUAD count changed: expected 19, observed {primary_n}")
    receipt = {
        "author_panel_unique_models": len(rows),
        "strict_luad_primary_models": primary_n,
        "broad_nsclc_sensitivity_models": broad_n,
        "decision_counts": decisions,
        "precedence_rule": "Cellosaurus disease/identity issue overrides broader legacy subtype labels",
        "strict_rule": "Cellosaurus disease contains Lung adenocarcinoma and model is not misidentified",
        "output": str(output),
    }
    (args.output_dir / "ddbj_cell_line_histology_audit_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
