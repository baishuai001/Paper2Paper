#!/usr/bin/env python3
"""Fetch public PDMR patient/model annotations for the frozen Figure 1 cohort.

All network access is intentionally performed on the cloud run host.  The raw
CSV and HTML pages are retained there as an audit trail; only compact parsed
tables and receipts are intended for version control.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import requests
from bs4 import BeautifulSoup


PDMR_EXPORT = "https://pdmdb.cancer.gov/web/apex/r/dctd_01/pdm/41?request=CSV"
PDMR_LOT = "https://pdmdb.cancer.gov/web/apex/r/dctd_01/pdm/26"
PDMR_PATIENT = "https://pdmdb.cancer.gov/web/apex/r/dctd_01/pdm/3"


def clean(value: object) -> str:
    if value is None:
        return ""
    return " ".join(str(value).replace("\xa0", " ").split())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def element_value(soup: BeautifulSoup, element_id: str) -> str:
    element = soup.find(id=element_id)
    if element is None:
        return ""
    if element.name in {"input", "textarea", "select"}:
        return clean(element.get("value", ""))
    return clean(element.get_text(" ", strip=True))


def selected_radio(soup: BeautifulSoup, prefix: str) -> str:
    for element in soup.find_all("input", id=re.compile(rf"^{re.escape(prefix)}_\d+$")):
        if element.has_attr("checked"):
            return clean(element.get("data-display") or element.get("value"))
    return ""


def smallest_table(soup: BeautifulSoup, required_header: str):
    candidates = []
    for table in soup.find_all("table"):
        headers = [clean(cell.get_text(" ", strip=True)) for cell in table.find_all("th")]
        if required_header in headers:
            candidates.append((len(headers), table))
    if not candidates:
        return None
    return min(candidates, key=lambda item: item[0])[1]


def table_rows(soup: BeautifulSoup, required_header: str) -> list[dict[str, str]]:
    table = smallest_table(soup, required_header)
    if table is None:
        return []
    header_row = None
    for row in table.find_all("tr"):
        values = [clean(cell.get_text(" ", strip=True)) for cell in row.find_all(["th", "td"], recursive=False)]
        if required_header in values:
            header_row = values
            break
    if not header_row:
        return []
    records: list[dict[str, str]] = []
    for row in table.find_all("tr"):
        values = [clean(cell.get_text(" ", strip=True)) for cell in row.find_all("td", recursive=False)]
        if len(values) == len(header_row) and any(values):
            records.append(dict(zip(header_row, values)))
    return records


def choose_export_row(rows: list[dict[str, str]], patient_id: str, specimen_id: str) -> dict[str, str] | None:
    patient_rows = [row for row in rows if clean(row.get("Patient ID")) == patient_id]
    exact = [row for row in patient_rows if clean(row.get("Specimen ID")) == specimen_id]
    return (exact or patient_rows or [None])[0]


def parse_patient_page(soup: BeautifulSoup, specimen_id: str) -> dict[str, str]:
    social_rows = table_rows(soup, "Has Smoked 100 Cigarettes")
    social = social_rows[0] if social_rows else {}
    specimen_rows = table_rows(soup, "Biopsy Site")
    exact_specimen = [row for row in specimen_rows if clean(row.get("Specimen ID")) == specimen_id]
    specimen = (exact_specimen or specimen_rows or [{}])[0]

    gender_code = element_value(soup, "P3_GENDER_HIDDENVALUE")
    gender = {"F": "Female", "M": "Male", "U": "Unknown"}.get(gender_code, gender_code)
    molecular_text = element_value(soup, "P3_KNOWNGENETICMUTATIONS_DISPLAY")
    if not molecular_text:
        molecular_text = element_value(soup, "P3_KNOWNGENETICMUTATIONS")

    return {
        "sex": gender,
        "age_at_diagnosis": element_value(soup, "P3_AGEATDIAGNOSISRANGE"),
        "known_metastatic_disease": selected_radio(soup, "P3_HASMETASTATICDISEASESEQNBR"),
        "known_molecular_ihc_data": molecular_text,
        "has_smoked_100_cigarettes": clean(social.get("Has Smoked 100 Cigarettes")),
        "total_pack_years": clean(social.get("Total Pack Years")),
        "tobacco_use_history": clean(social.get("Tobacco Use History")),
        "biopsy_site": clean(specimen.get("Biopsy Site")),
        "tissue_type": clean(specimen.get("Tissue Type")),
        "age_at_sampling": clean(specimen.get("Age at Sampling")),
        "collection_date": clean(specimen.get("Collection Date")),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cloud_run_root", type=Path)
    parser.add_argument("--delay-seconds", type=float, default=0.15)
    args = parser.parse_args()

    root = args.cloud_run_root.resolve()
    manifest_path = root / "data" / "processed" / "pdmr_manifest.tsv"
    raw_dir = root / "data" / "raw" / "annotations" / "pdmr"
    out_dir = root / "data" / "processed" / "annotations"
    audit_dir = root / "audit" / "figure1_annotations"
    for directory in (raw_dir, out_dir, audit_dir):
        directory.mkdir(parents=True, exist_ok=True)

    with manifest_path.open(encoding="utf-8", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))

    session = requests.Session()
    session.headers.update({"User-Agent": "Paper2Paper/1.0 (public-data reproducibility audit)"})
    response = session.get(PDMR_EXPORT, timeout=90)
    response.raise_for_status()
    export_path = raw_dir / "pdmr_public_models.csv"
    export_path.write_bytes(response.content)
    export_rows = list(csv.DictReader(io.StringIO(response.text)))

    output: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    for index, item in enumerate(manifest, start=1):
        patient_id = clean(item.get("patient_id"))
        specimen_id = clean(item.get("specimen_id"))
        sample_id = clean(item.get("sample_id"))
        export_row = choose_export_row(export_rows, patient_id, specimen_id)
        record = {
            "sample_id": sample_id,
            "patient_id": patient_id,
            "specimen_id": specimen_id,
            "group": clean(item.get("group")),
            "parse_status": "unmapped_export",
            "pdmr_view_id": "",
            "pdmr_patient_sequence": "",
            "disease_body_location": "",
            "oncotree_code": "",
            "diagnosis_subtype": "",
            "max_passage": "",
            "sex": "",
            "age_at_diagnosis": "",
            "known_metastatic_disease": "",
            "known_molecular_ihc_data": "",
            "has_smoked_100_cigarettes": "",
            "total_pack_years": "",
            "tobacco_use_history": "",
            "biopsy_site": "",
            "tissue_type": "",
            "age_at_sampling": "",
            "collection_date": "",
            "source": "NCI PDMR public database",
        }
        if export_row is None:
            errors.append({"sample_id": sample_id, "error": "patient/specimen not found in public CSV"})
            output.append(record)
            continue

        view_id = clean(export_row.get("View"))
        record.update(
            {
                "pdmr_view_id": view_id,
                "disease_body_location": clean(export_row.get("Disease BodyLocation")),
                "oncotree_code": clean(export_row.get("OncoTreeCode")),
                "diagnosis_subtype": clean(export_row.get("DiagnosisSubtype")),
                "max_passage": clean(export_row.get("Max.Passage")),
            }
        )
        try:
            lot_response = session.get(
                PDMR_LOT,
                params={"p26_distributionlotseqnbr": view_id, "clear": "26"},
                timeout=90,
            )
            lot_response.raise_for_status()
            lot_path = raw_dir / f"{sample_id}.distribution_lot.html"
            lot_path.write_bytes(lot_response.content)
            match = re.search(r"p3_patientseqnbr(?:%3D|=)([0-9]+)", lot_response.text, re.I)
            if match is None:
                raise RuntimeError("patient sequence was not present on distribution-lot page")
            patient_sequence = match.group(1)
            record["pdmr_patient_sequence"] = patient_sequence

            patient_response = session.get(
                PDMR_PATIENT,
                params={"p3_patientseqnbr": patient_sequence, "clear": "3"},
                timeout=90,
            )
            patient_response.raise_for_status()
            patient_path = raw_dir / f"{sample_id}.patient.html"
            patient_path.write_bytes(patient_response.content)
            record.update(parse_patient_page(BeautifulSoup(patient_response.text, "lxml"), specimen_id))
            record["parse_status"] = "ok"
        except Exception as exc:  # keep partial public metadata and report every failure
            record["parse_status"] = "page_error"
            errors.append({"sample_id": sample_id, "error": clean(exc)})
        output.append(record)
        if index != len(manifest):
            time.sleep(args.delay_seconds)

    fieldnames = list(output[0].keys()) if output else []
    output_path = out_dir / "pdmr_figure1_annotations.tsv"
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

    receipt = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_urls": [PDMR_EXPORT, PDMR_LOT, PDMR_PATIENT],
        "manifest_samples": len(manifest),
        "export_rows": len(export_rows),
        "parsed_ok": sum(row["parse_status"] == "ok" for row in output),
        "partial_or_failed": sum(row["parse_status"] != "ok" for row in output),
        "raw_export_sha256": sha256(export_path),
        "parsed_table_sha256": sha256(output_path),
        "errors": errors,
    }
    (audit_dir / "pdmr_annotation_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
