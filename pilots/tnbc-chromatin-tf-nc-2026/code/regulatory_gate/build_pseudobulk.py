#!/usr/bin/env python3
"""Build manifest-defined patient Cancer-cell pseudobulks from H5AD raw counts."""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from common import output_manifest, read_h5_column, sha256_file, write_json
from phenotype import PhenotypeManifest, load_obs, patient_table


def raw_gene_symbols(h5ad: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(h5ad, "r") as handle:
        if "raw" not in handle or "X" not in handle["raw"] or "var" not in handle["raw"]:
            raise ValueError("H5AD has no raw/X and raw/var")
        group = handle["raw/var"]
        if "GeneSymbol" not in group:
            raise ValueError("H5AD raw/var has no GeneSymbol")
        symbols = read_h5_column(group, "GeneSymbol").astype(str)
        n_cells = (
            np.asarray(group["n_cells"][:], dtype=float)
            if "n_cells" in group
            else np.zeros(len(symbols), dtype=float)
        )
    frame = pd.DataFrame(
        {"symbol": symbols, "n_cells": n_cells, "column": np.arange(len(symbols), dtype=int)}
    )
    frame = frame[frame["symbol"].ne("")].sort_values(
        ["symbol", "n_cells", "column"], ascending=[True, False, True]
    )
    selected = frame.drop_duplicates("symbol", keep="first").sort_values("column")
    return selected["symbol"].to_numpy(str), selected["column"].to_numpy(int)


def write_expression(path: Path, matrix: np.ndarray, genes: np.ndarray, patients: list[str]) -> None:
    frame = pd.DataFrame(matrix.T, index=genes, columns=patients)
    frame.index.name = "gene"
    with gzip.open(path, "wt", encoding="utf-8", newline="") as stream:
        frame.to_csv(stream, sep="\t", lineterminator="\n", float_format="%.8g")


def run(h5ad: Path, manifest_path: Path, output_dir: Path, block_size: int) -> dict[str, object]:
    manifest = PhenotypeManifest.load(manifest_path)
    obs = load_obs(h5ad, manifest.required_obs_fields)
    patients, eligible = patient_table(obs, manifest)
    if patients["label_n"].ne(1).any() or patients["group"].eq("conflict").any():
        conflicts = patients.loc[patients["label_n"].ne(1), "donor_id"].tolist()
        raise ValueError(f"Contradictory phenotype labels in scoped cohort: {conflicts[:20]}")
    if patients["dataset"].eq("").any() or patients["dataset_n"].ne(1).any():
        conflicts = patients.loc[(patients["dataset"].eq("")) | patients["dataset_n"].ne(1), "donor_id"].tolist()
        raise ValueError(f"Missing or multiple dataset assignments: {conflicts[:20]}")
    patients = patients.sort_values("donor_id").reset_index(drop=True)
    patient_to_row = {patient: index for index, patient in enumerate(patients["donor_id"])}

    patient_column = str(manifest["patient_column"])
    cell_column = str(manifest["analysis_cell_column"])
    selected_cells = np.flatnonzero(
        eligible
        & obs[patient_column].isin(patient_to_row).to_numpy()
        & obs[cell_column].isin(manifest["analysis_cell_values"]).to_numpy()
    )
    if len(selected_cells) == 0:
        raise ValueError("No analysis cells in manifest-defined cohort")

    genes, gene_columns = raw_gene_symbols(h5ad)
    counts = np.zeros((len(patients), len(genes)), dtype=np.float64)
    cell_counts = np.zeros(len(patients), dtype=np.int64)
    backed = ad.read_h5ad(h5ad, backed="r")
    try:
        for start in range(0, len(selected_cells), block_size):
            rows = selected_cells[start : start + block_size]
            matrix = backed.raw.X[rows, :]
            if not sparse.issparse(matrix):
                matrix = sparse.csr_matrix(matrix)
            matrix = matrix[:, gene_columns]
            codes = np.fromiter(
                (patient_to_row[value] for value in obs.iloc[rows][patient_column]),
                dtype=np.int64,
                count=len(rows),
            )
            grouping = sparse.csr_matrix(
                (np.ones(len(rows)), (codes, np.arange(len(rows)))),
                shape=(len(patients), len(rows)),
            )
            counts += (grouping @ matrix).toarray()
            cell_counts += np.bincount(codes, minlength=len(patients))
    finally:
        backed.file.close()

    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("raw/X pseudobulk contains non-finite or negative values")
    if not np.allclose(counts, np.rint(counts), atol=1e-6):
        raise ValueError("raw/X pseudobulk is not integer count data")
    if not np.array_equal(cell_counts, patients["analysis_cells"].to_numpy(int)):
        raise ValueError("Aggregated cell counts disagree with metadata-derived counts")

    library = counts.sum(axis=1)
    if (library <= 0).any():
        raise ValueError("At least one patient has a zero pseudobulk library")
    cpm = counts / library[:, None] * 1e6
    keep = (cpm >= 1.0).mean(axis=0) >= 0.10
    if keep.sum() < 10_000:
        raise ValueError(f"Only {int(keep.sum())} genes survived the frozen pseudobulk filter")

    patients["analysis_cells_aggregated"] = cell_counts
    patients["raw_library_size"] = library
    patients["primary_eligible"] = patients["analysis_cells_aggregated"].ge(
        int(manifest["primary_min_cells"])
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = output_dir / "patient_metadata.tsv"
    counts_path = output_dir / "pseudobulk_counts.tsv.gz"
    log_path = output_dir / "log2cpm.tsv.gz"
    filter_path = output_dir / "gene_filter.tsv.gz"
    patients.to_csv(metadata_path, sep="\t", index=False, lineterminator="\n")
    write_expression(counts_path, counts[:, keep], genes[keep], patients["donor_id"].tolist())
    write_expression(log_path, np.log2(cpm[:, keep] + 1.0), genes[keep], patients["donor_id"].tolist())
    with gzip.open(filter_path, "wt", encoding="utf-8", newline="") as stream:
        pd.DataFrame({"gene": genes, "retained": keep}).to_csv(
            stream, sep="\t", index=False, lineterminator="\n"
        )

    thresholds = [int(manifest["primary_min_cells"])] + [int(x) for x in manifest.raw.get("sensitivity_min_cells", [])]
    group_counts: dict[str, dict[str, int]] = {}
    for threshold in sorted(set(thresholds)):
        subset = patients[patients["analysis_cells_aggregated"].ge(threshold)]
        group_counts[str(threshold)] = {
            "patients": int(len(subset)),
            "case": int(subset["group"].eq("case").sum()),
            "control": int(subset["group"].eq("control").sum()),
            "datasets": int(subset["dataset"].nunique()),
        }

    outputs = [metadata_path, counts_path, log_path, filter_path]
    receipt = {
        "status": "passed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phenotype_id": manifest["phenotype_id"],
        "manifest_sha256": sha256_file(manifest_path),
        "patients_all": int(len(patients)),
        "analysis_cells": int(len(selected_cells)),
        "unique_gene_symbols": int(len(genes)),
        "retained_genes": int(keep.sum()),
        "raw_integer_validation": True,
        "group_counts_by_cell_threshold": group_counts,
        "outputs": output_manifest(outputs),
    }
    write_json(output_dir / "pseudobulk_receipt.json", receipt)
    print(json.dumps(receipt, ensure_ascii=False))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--block-size", type=int, default=2048)
    args = parser.parse_args()
    run(args.h5ad, args.manifest, args.output_dir, args.block_size)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
