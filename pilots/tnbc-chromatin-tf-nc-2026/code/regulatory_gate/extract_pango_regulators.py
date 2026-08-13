#!/usr/bin/env python3
"""Extract the exact PAN-GO regulator symbols from TNBC Supplementary Data 2."""

from __future__ import annotations

import argparse
import json
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET

from common import output_manifest, sha256_file, write_json


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def column_index(cell_ref: str) -> int:
    match = re.match(r"[A-Z]+", cell_ref)
    if not match:
        raise ValueError(f"Invalid XLSX cell reference: {cell_ref}")
    value = 0
    for letter in match.group(0):
        value = value * 26 + ord(letter) - 64
    return value - 1


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    values: list[str] = []
    with archive.open("xl/sharedStrings.xml") as stream:
        for _, element in ET.iterparse(stream, events=("end",)):
            if element.tag == f"{{{MAIN_NS}}}si":
                values.append("".join(node.text or "" for node in element.iter(f"{{{MAIN_NS}}}t")))
                element.clear()
    return values


def target_sheet(archive: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in relationships.findall(f"{{{PKG_REL_NS}}}Relationship")
    }
    for sheet in workbook.findall(f".//{{{MAIN_NS}}}sheet"):
        if sheet.attrib["name"] == sheet_name:
            path = targets[sheet.attrib[f"{{{REL_NS}}}id"]].lstrip("/")
            return path if path.startswith("xl/") else f"xl/{path}"
    raise ValueError(f"Sheet not found: {sheet_name}")


def cell_value(cell: ET.Element, strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.iter(f"{{{MAIN_NS}}}t"))
    value = cell.find(f"{{{MAIN_NS}}}v")
    if value is None or value.text is None:
        return ""
    return strings[int(value.text)] if cell_type == "s" else value.text


def extract(path: Path, sheet_name: str) -> list[str]:
    with zipfile.ZipFile(path) as archive:
        strings = shared_strings(archive)
        sheet_path = target_sheet(archive, sheet_name)
        symbol_column: int | None = None
        symbols: list[str] = []
        with archive.open(sheet_path) as stream:
            for _, row in ET.iterparse(stream, events=("end",)):
                if row.tag != f"{{{MAIN_NS}}}row":
                    continue
                values = {
                    column_index(cell.attrib.get("r", "A1")): cell_value(cell, strings).strip()
                    for cell in row.findall(f"{{{MAIN_NS}}}c")
                }
                if symbol_column is None:
                    matches = [index for index, value in values.items() if value == "gene_symbol"]
                    if len(matches) != 1:
                        raise ValueError("Expected one gene_symbol header")
                    symbol_column = matches[0]
                else:
                    symbol = values.get(symbol_column, "")
                    if symbol:
                        symbols.append(symbol)
                row.clear()
    return symbols


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sheet", default="PANGO_TR_RegulationOfTranscript")
    parser.add_argument("--methods-claimed", type=int, default=2139)
    parser.add_argument("--expected-records", type=int, default=2138)
    parser.add_argument("--expected-unique", type=int, default=2059)
    args = parser.parse_args()
    records = extract(args.xlsx, args.sheet)
    symbols = sorted(set(records))
    if len(records) != args.expected_records:
        raise ValueError(f"PAN-GO record count {len(records)} != frozen expectation {args.expected_records}")
    if len(symbols) != args.expected_unique:
        raise ValueError(f"PAN-GO unique symbol count {len(symbols)} != frozen expectation {args.expected_unique}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "pango_regulators.txt"
    output.write_text("\n".join(symbols) + "\n", encoding="utf-8", newline="\n")
    receipt = {
        "status": "passed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {"name": args.xlsx.name, "bytes": args.xlsx.stat().st_size, "sha256": sha256_file(args.xlsx)},
        "sheet": args.sheet,
        "methods_claimed_symbols": args.methods_claimed,
        "annotation_records": len(records),
        "expected_annotation_records": args.expected_records,
        "unique_nonempty_symbols": len(symbols),
        "expected_unique_symbols": args.expected_unique,
        "duplicate_symbol_records": len(records) - len(symbols),
        "methods_minus_observed_records": args.methods_claimed - len(records),
        "output": output_manifest([output]),
    }
    write_json(args.output_dir / "pango_receipt.json", receipt)
    print(json.dumps(receipt, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
