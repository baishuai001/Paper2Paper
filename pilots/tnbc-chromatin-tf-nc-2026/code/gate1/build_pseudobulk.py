#!/usr/bin/env python3
"""Build author-Cancer-cell patient pseudobulks from integer raw/X on cloud."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from crc_cohort import eligible_patient_table, load_obs
from gate1_common import read_obs_column, sha256_file, write_json


def raw_gene_symbols(h5ad: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(h5ad, "r") as handle:
        group = handle["raw/var"]
        symbols = read_obs_column(group, "GeneSymbol").astype(str)
        n_cells = np.asarray(group["n_cells"][:], dtype=float)
    frame = pd.DataFrame({"symbol": symbols, "n_cells": n_cells, "column": np.arange(len(symbols))})
    frame = frame[frame["symbol"].ne("")].sort_values(["symbol", "n_cells", "column"], ascending=[True, False, True])
    selected = frame.drop_duplicates("symbol", keep="first").sort_values("column")
    return selected["symbol"].to_numpy(str), selected["column"].to_numpy(int)


def run(h5ad: Path, output_dir: Path, block_size: int = 2048) -> dict[str, object]:
    fields = [
        "donor_id", "sample_id", "dataset", "study_id", "sample_type", "tumor_source", "medical_condition",
        "treatment_status_before_resection", "enrichment_cell_types", "cell_type_coarse_crc_atlas",
        "cell_type_middle_crc_atlas", "immune_infiltration_type", "CMS_type", "microsatellite_status",
        "anatomic_location", "tumor_stage", "age", "sex",
    ]
    obs = load_obs(h5ad, fields)
    patient, eligible = eligible_patient_table(obs)
    patient = patient[patient["immune_label_n"].eq(1)].copy().sort_values("donor_id").reset_index(drop=True)
    donor_to_row = {donor: i for i, donor in enumerate(patient["donor_id"])}
    selected_cells = np.flatnonzero(
        eligible & obs["donor_id"].isin(donor_to_row) & obs["cell_type_coarse_crc_atlas"].eq("Cancer cell")
    )
    if not len(selected_cells):
        raise ValueError("No eligible author-labelled Cancer cells")

    genes, gene_columns = raw_gene_symbols(h5ad)
    totals = np.zeros((len(patient), len(genes)), dtype=np.float64)
    cell_counts = np.zeros(len(patient), dtype=np.int64)
    data = ad.read_h5ad(h5ad, backed="r")
    try:
        for start in range(0, len(selected_cells), block_size):
            rows = selected_cells[start : start + block_size]
            matrix = data.raw.X[rows, :]
            if not sparse.issparse(matrix):
                matrix = sparse.csr_matrix(matrix)
            matrix = matrix[:, gene_columns]
            codes = np.fromiter((donor_to_row[x] for x in obs.iloc[rows]["donor_id"]), dtype=np.int64, count=len(rows))
            grouping = sparse.csr_matrix((np.ones(len(rows)), (codes, np.arange(len(rows)))), shape=(len(patient), len(rows)))
            totals += (grouping @ matrix).toarray()
            cell_counts += np.bincount(codes, minlength=len(patient))
    finally:
        data.file.close()

    if not np.isfinite(totals).all() or (totals < 0).any() or not np.allclose(totals, np.rint(totals), atol=1e-6):
        raise ValueError("raw/X pseudobulk is not finite, nonnegative integer count data")
    patient["cancer_cells_aggregated"] = cell_counts
    if not np.array_equal(patient["cancer_cells"].to_numpy(int), cell_counts):
        raise ValueError("Cancer-cell aggregation count disagrees with metadata audit")

    library = totals.sum(axis=1)
    cpm = totals / np.maximum(library[:, None], 1.0) * 1e6
    keep = (cpm >= 1.0).mean(axis=0) >= 0.10
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_dir / "G3_pseudobulk_counts.npz", counts=totals[:, keep], genes=genes[keep])
    np.savez_compressed(output_dir / "G3_log2cpm.npz", log2cpm=np.log2(cpm[:, keep] + 1.0), genes=genes[keep])
    patient["raw_library_size"] = library
    patient.to_csv(output_dir / "G3_patient_metadata.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"gene": genes, "retained": keep}).to_csv(
        output_dir / "G3_gene_filter.tsv", sep="\t", index=False, lineterminator="\n"
    )
    receipt = {
        "status": "passed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "patients": int(len(patient)),
        "Cancer_cells": int(len(selected_cells)),
        "unique_gene_symbols": int(len(genes)),
        "retained_genes": int(keep.sum()),
        "raw_integer_validation": True,
        "threshold_counts": {str(n): int((cell_counts >= n).sum()) for n in (20, 50, 100, 200)},
        "outputs": {},
    }
    for path in sorted(output_dir.glob("G3_*")):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    write_json(output_dir / "G3_pseudobulk_receipt.json", receipt)
    print(receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--block-size", type=int, default=2048)
    args = parser.parse_]u×kh‘éì¶»§q«^uÍÑ½À%Á¡•¹½ÑåÁ•}Ù…±¥‘…Ñ¥½¸%É¥Ñ¥…°%%¹‘•Á•¹‘•¹Ğ™½ÕÈµ±…ÍÌI$É¥Ñ•É¥„™…¥±•%ÈÉ••¥ÁĞè€À¸ÄÜÌ¡¥•É…É¡¥…°…¹€À¸ÈÌÈµ•‘¥…¸,µµ•…¹Ì%Q9	µÍÑå±”‘½İ¹ÍÑÉ•…´ÑÉ…¹Í™•È¥Ì¹½Ğ…ÕÑ¡½É¥é•%ÁÁ±ä™É½é•¸ÍÑ½ÀÉÕ±”%ÕÑ½µ…Ñ¥ŒÙ•É‘¥Ğ%0°ÍÑ½Á}¹½ÜõÑÉÕ”%MÑ½À‘½İ¹ÍÑÉ•…´É•…°µ‘…Ñ„…¹…±åÍ¥Ì%±½Í•%ÑÉÕ”%	¥¹…Éä4Í¥¹…°É•µ…¥¹ÌÍÑÉ½¹œ‰ÕĞ…¹¹½Ğ½Ù•ÉÉ¥‘”½¹©Õ¹Ñ¥Ù”ÉÕ±”4