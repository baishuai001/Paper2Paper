#!/usr/bin/env python3
"""Create a deduplicated BED6 protein-coding gene TSS set for ataqv."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


def attributes(value: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)"', value))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_gz", type=Path)
    parser.add_argument("output_bed", type=Path)
    args = parser.parse_args()
    tss: set[tuple[str, int, int, str, str]] = set()
    with gzip.open(args.gtf_gz, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chrom, start_text, end_text, strand = fields[0], fields[3], fields[4], fields[6]
            if chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX"}:
                continue
            attrs = attributes(fields[8])
            if (attrs.get("gene_type") or attrs.get("gene_biotype")) != "protein_coding":
                continue
            symbol = attrs.get("gene_name", attrs.get("gene_id", "unknown"))
            start = int(start_text) - 1
            end = int(end_text)
            position = start if strand == "+" else end - 1
            tss.add((chrom, position, position + 1, symbol, strand))
    chrom_order = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23}
    ordered = sorted(tss, key=lambda row: (chrom_order[row[0]], row[1], row[3]))
    args.output_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.output_bed.open("w", encoding="utf-8") as handle:
        for chrom, start, end, symbol, strand in ordered:
            handle.write(f"{chrom}\t{start}\t{end}\t{symbol}\t0\t{strand}\n")
    print(f"GENCODE_GENE_TSS_COMPLETE\t{len(ordered)}\t{args.output_bed}")


if __name__ == "__main__":
    main()
