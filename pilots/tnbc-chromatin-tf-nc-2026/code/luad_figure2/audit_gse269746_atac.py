#!/usr/bin/env python3
"""Create an independent-model ATAC manifest for GSE269746 from SRA XML."""

from __future__ import annotations

import argparse
import csv
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path


TITLE_RE = re.compile(r"ATAC-(MGH[0-9]+)-(adeno|transformedSCLC|denovoSCLC)", re.I)
STATE = {
    "adeno": "LUAD",
    "transformedsclc": "transformed_SCLC",
    "denovosclc": "de_novo_SCLC",
}


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sra_xml", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    root = ET.parse(args.sra_xml).getroot()
    rows: list[dict] = []
    for package in root.findall(".//EXPERIMENT_PACKAGE"):
        experiment = package.find("EXPERIMENT")
        sample = package.find("SAMPLE")
        runs = package.findall("./RUN_SET/RUN")
        if experiment is None or sample is None or len(runs) != 1:
            raise RuntimeError("Every ATAC experiment must have one SAMPLE and exactly one RUN")
        title = experiment.findtext("TITLE") or ""
        match = TITLE_RE.search(title)
        if not match:
            raise RuntimeError(f"Cannot classify experiment title: {title}")
        model, raw_state = match.groups()
        state = STATE[raw_state.lower()]
        run = runs[0]
        external = {
            item.findtext("DB"): item.findtext("LABEL") or item.findtext("ID")
            for item in package.findall("./EXPERIMENT/IDENTIFIERS/EXTERNAL_ID")
        }
        gsm_match = re.search(r"(GSM[0-9]+)", title)
        row = {
            "model_id": model,
            "disease_state": state,
            "figure2_primary_include": str(state == "LUAD").upper(),
            "independent_biological_model": "TRUE",
            "experiment": experiment.attrib.get("accession", ""),
            "run": run.attrib.get("accession", ""),
            "biosample": next(
                (
                    x.findtext("LABEL") or x.findtext("ID") or ""
                    for x in sample.findall("./IDENTIFIERS/EXTERNAL_ID")
                    if x.findtext("DB") == "BioSample"
                ),
                "",
            ),
            "sra_sample": sample.attrib.get("accession", ""),
            "gsm": gsm_match.group(1) if gsm_match else external.get("GEO", ""),
            "library_strategy": experiment.findtext("./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY") or "",
            "library_layout": (
                list(experiment.find("./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT"))[0].tag
                if experiment.find("./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT") is not None
                else ""
            ),
            "instrument_model": experiment.findtext("./PLATFORM/ILLUMINA/INSTRUMENT_MODEL") or "",
            "total_spots": run.attrib.get("total_spots", ""),
            "total_bases": run.attrib.get("total_bases", ""),
            "run_size_bytes": run.attrib.get("size", ""),
            "title": title,
            "source_series": "GSE269746",
            "source_bioproject": "PRJNA1123604",
        }
        rows.append(row)

    rows.sort(key=lambda row: (row["disease_state"], row["model_id"]))
    model_ids = [row["model_id"] for row in rows]
    if len(rows) != 24 or len(set(model_ids)) != 24:
        raise RuntimeError(f"Expected 24 independent ATAC PDX models; observed {len(rows)}/{len(set(model_ids))}")
    primary = [row for row in rows if row["figure2_primary_include"] == "TRUE"]
    if len(primary) != 13:
        raise RuntimeError(f"Expected 13 LUAD PDX models; observed {len(primary)}")

    columns = list(rows[0])
    write_tsv(args.output_dir / "gse269746_atac_all_models.tsv", rows, columns)
    write_tsv(args.output_dir / "gse269746_atac_luad_models.tsv", primary, columns)
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["disease_state"]] = counts.get(row["disease_state"], 0) + 1
    receipt = {
        "series_reported_total_pdx": 26,
        "atac_independent_models": len(rows),
        "atac_state_counts": counts,
        "figure2_luad_pdx_models": len(primary),
        "run_per_model": 1,
        "primary_rule": "experiment title suffix adeno; one MGH model counted once",
        "excluded_from_luad_primary": "transformed SCLC and de novo SCLC",
        "source_xml": str(args.sra_xml),
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "gse269746_atac_manifest_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
