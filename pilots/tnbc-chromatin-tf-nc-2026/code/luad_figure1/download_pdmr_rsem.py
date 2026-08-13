#!/usr/bin/env python3
"""Download frozen NCI PDMR RSEM gene-result files and build a TPM matrix."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import time
import urllib.error
import urllib.request


def download(url: str, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": "Paper2Paper/1.0"})
    for attempt in range(6):
        try:
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as handle:
                while block := response.read(1024 * 1024):
                    handle.write(block)
            partial.replace(destination)
            return
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            if attempt == 5:
                raise
            time.sleep(3 * (attempt + 1))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tpm(path: Path) -> dict[str, float]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "gene_id" not in reader.fieldnames or "TPM" not in reader.fieldnames:
            raise ValueError(f"Unexpected RSEM columns in {path}: {reader.fieldnames}")
        result = {}
        for row in reader:
            # PDMR's RSEM files use gene symbols/UCSC knownGene labels in the
            # gene_id field, not Ensembl IDs. Preserve the identifier exactly;
            # splitting on a period would corrupt valid symbols.
            gene = row["gene_id"].strip()
            value = float(row["TPM"])
            if gene in result:
                result[gene] = (result[gene] + value) / 2
            else:
                result[gene] = value
        return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest_tsv", type=Path)
    parser.add_argument("raw_dir", type=Path)
    parser.add_argument("matrix_tsv", type=Path)
    parser.add_argument("receipt_json", type=Path)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()

    with args.manifest_tsv.open(encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))

    if not 1 <= args.workers <= 12:
        raise SystemExit("--workers must be between 1 and 12")

    def fetch(row: dict[str, str]) -> tuple[str, dict[str, float], dict[str, object]]:
        destination = args.raw_dir / f"{row['sample_id']}.RSEM.genes.results"
        download(row["rsem_gene_url"], destination)
        return (
            row["sample_id"],
            read_tpm(destination),
            {
                "sample_id": row["sample_id"],
                "url": row["rsem_gene_url"],
                "path": str(destination),
                "size": destination.stat().st_size,
                "sha256": sha256(destination),
            },
        )

    matrices: dict[str, dict[str, float]] = {}
    file_receipts = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for sample_id, matrix, file_receipt in executor.map(fetch, manifest):
            matrices[sample_id] = matrix
            file_receipts.append(file_receipt)

    common_genes = sorted(set.intersection(*(set(matrix) for matrix in matrices.values())))
    if len(common_genes) < 15000:
        raise SystemExit(f"Only {len(common_genes)} genes shared by PDMR RSEM files")
    args.matrix_tsv.parent.mkdir(parents=True, exist_ok=True)
    sample_ids = [row["sample_id"] for row in manifest]
    with args.matrix_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene"] + sample_ids)
        for gene in common_genes:
            writer.writerow([gene] + [matrices[sample][gene] for sample in sample_ids])

    receipt = {
        "samples": len(sample_ids),
        "shared_gene_symbols": len(common_genes),
        "matrix_path": str(args.matrix_tsv),
        "matrix_sha256": sha256(args.matrix_tsv),
        "download_workers": args.workers,
        "files": file_receipts,
    }
    args.receipt_json.parent.mkdir(parents=True, exist_ok=True)
    args.receipt_json.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in receipt.items() if k != "files"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
