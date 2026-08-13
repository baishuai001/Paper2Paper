#!/usr/bin/env python3
"""Validate TCGA expression and intersect the frozen PAN-GO regulator list."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

from common import output_manifest, write_json


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--pango", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    genes: list[str] = []
    samples: list[str] = []
    with args.expression.open(encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if len(header) < 551 or header[0] != "gene":
            raise ValueError("TCGA expression header is not gene + >=550 participants")
        samples = header[1:]
        if len(samples) != len(set(samples)):
            raise ValueError("Duplicate TCGA participant columns")
        expected_fields = len(header)
        for line_number, row in enumerate(reader, start=2):
            if len(row) != expected_fields:
                raise ValueError(f"TCGA expression row {line_number} has {len(row)} fields; expected {expected_fields}")
            genes.append(row[0])
    if len(genes) < 10_000 or len(genes) != len(set(genes)) or "" in genes:
        raise ValueError("TCGA expression gene universe fails uniqueness or size checks")

    pango = [line.strip() for line in args.pango.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(pango) != 2_059 or len(pango) != len(set(pango)):
        raise ValueError("PAN-GO list is not the frozen 2,059-symbol unique list")
    gene_set = set(genes)
    retained = sorted(set(pango) & gene_set)
    missing = sorted(set(pango) - gene_set)
    if len(retained) < 1_500:
        raise ValueError(f"Only {len(retained)} PAN-GO regulators are present in TCGA expression")

    regulator_path = args.output_dir / "aracne_regulators.txt"
    missing_path = args.output_dir / "pango_missing_from_tcga.txt"
    regulator_path.write_text("\n".join(retained) + "\n", encoding="utf-8", newline="\n")
    missing_path.write_text("\n".join(missing) + ("\n" if missing else ""), encoding="utf-8", newline="\n")
    receipt = {
        "status": "passed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "tcga_participants": len(samples),
        "tcga_genes": len(genes),
        "pango_symbols": len(pango),
        "aracne_regulators": len(retained),
        "pango_missing_from_tcga": len(missing),
        "outputs": output_manifest([regulator_path, missing_path]),
    }
    write_json(args.output_dir / "aracne_input_receipt.json", receipt)
    print(json.dumps(receipt, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
