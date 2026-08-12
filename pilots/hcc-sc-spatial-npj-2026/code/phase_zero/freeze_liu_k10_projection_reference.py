#!/usr/bin/env python3
"""Freeze a Liu-only K=10 rank-score projection reference before reading Che results."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr


K = 10
PROGRAMS = list(range(1, K + 1))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_signatures(top_gene_path: Path, top_n: int = 50) -> dict[int, list[str]]:
    table = pd.read_csv(top_gene_path, sep="\t")
    table = table.loc[(table["k"].astype(int) == K) & (table["rank"].astype(int) <= top_n)].copy()
    table["program"] = table["program"].astype(int)
    table["rank"] = table["rank"].astype(int)
    result: dict[int, list[str]] = {}
    for program in PROGRAMS:
        genes = table.loc[table["program"] == program].sort_values("rank")["gene"].astype(str).tolist()
        if len(genes) != top_n or len(set(genes)) != top_n:
            raise ValueError(f"Program {program} does not have {top_n} unique genes")
        result[program] = genes
    return result


def rank_scores(
    matrix: object,
    var_names: pd.Index,
    signatures: dict[int, list[str]],
    max_rank: int = 1500,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(np.asarray(matrix))
    matrix = matrix.tocsr()
    gene_index = pd.Index(var_names.astype(str))
    if not gene_index.is_unique:
        raise ValueError("Expression gene identifiers are not unique")
    union_genes = sorted({gene for genes in signatures.values() for gene in genes})
    union_locations = gene_index.get_indexer(union_genes)
    ordered_present_genes = [gene for gene, location in zip(union_genes, union_locations) if location >= 0]
    union_pos = {gene: i for i, gene in enumerate(ordered_present_genes)}
    signature_union_positions: dict[int, np.ndarray] = {}
    coverage_rows: list[dict[str, object]] = []
    for program, genes in signatures.items():
        present = [gene for gene in genes if gene in union_pos]
        missing = [gene for gene in genes if gene not in union_pos]
        signature_union_positions[program] = np.asarray([union_pos[gene] for gene in present], dtype=int)
        coverage_rows.append(
            {
                "program_id": program,
                "signature_genes": len(genes),
                "present_genes": len(present),
                "coverage_fraction": len(present) / len(genes),
                "missing_genes": ";".join(missing),
            }
        )
    coverage = pd.DataFrame(coverage_rows)
    if (coverage["present_genes"] < 45).any():
        bad = coverage.loc[coverage["present_genes"] < 45, "program_id"].tolist()
        raise ValueError(f"Signature coverage below 45/50 for programs {bad}")

    valid_union_locations = np.asarray([gene_index.get_loc(gene) for gene in ordered_present_genes])
    gene_to_union = np.full(matrix.shape[1], -1, dtype=np.int32)
    gene_to_union[valid_union_locations] = np.arange(len(ordered_present_genes), dtype=np.int32)
    score_values = np.zeros((matrix.shape[0], K), dtype=np.float64)
    union_ranks = np.full(len(ordered_present_genes), max_rank + 1.0, dtype=np.float64)
    for row_index in range(matrix.shape[0]):
        start, stop = matrix.indptr[row_index], matrix.indptr[row_index + 1]
        indices = matrix.indices[start:stop]
        values = matrix.data[start:stop]
        if len(values):
            ranks = np.minimum(rankdata(-values, method="average"), max_rank + 1.0)
            positions = gene_to_union[indices]
            keep = positions >= 0
            union_ranks[positions[keep]] = ranks[keep]
        for program in PROGRAMS:
            selected = signature_union_positions[program]
            n_genes = len(selected)
            u_stat = float(union_ranks[selected].sum() - n_genes * (n_genes + 1) / 2.0)
            score_values[row_index, program - 1] = np.clip(
                1.0 - u_stat / (n_genes * max_rank), 0.0, 1.0
            )
        if len(values):
            union_ranks[positions[keep]] = max_rank + 1.0
    return pd.DataFrame(score_values, columns=PROGRAMS), coverage


def paired_differences(
    scores: pd.DataFrame, obs: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    joined = scores.copy()
    joined.index = obs.index.astype(str)
    joined = joined.join(obs[["analysis_patient_id", "sample_type"]].astype(str))
    sample = (
        joined.groupby(["analysis_patient_id", "sample_type"], observed=True)[PROGRAMS]
        .mean()
        .reset_index()
    )
    wide = sample.pivot(index="analysis_patient_id", columns="sample_type", values=PROGRAMS)
    required_states = {"tumor", "metastasis"}
    if set(wide.columns.get_level_values(1)) != required_states:
        raise ValueError("Liu direction reference requires tumor and metastasis for every program")
    difference_rows: list[dict[str, object]] = []
    direction_rows: list[dict[str, object]] = []
    for program in PROGRAMS:
        delta = wide[(program, "metastasis")] - wide[(program, "tumor")]
        if delta.isna().any() or (delta == 0).all():
            raise ValueError(f"Program {program} has incomplete or zero paired differences")
        median_delta = float(np.median(delta))
        direction = "higher_in_liver_metastasis" if median_delta > 0 else "lower_in_liver_metastasis"
        for patient, value in delta.items():
            difference_rows.append(
                {"program_id": program, "analysis_patient_id": patient, "delta_lm_minus_primary": float(value)}
            )
        direction_rows.append(
            {
                "program_id": program,
                "liu_median_delta_lm_minus_primary": median_delta,
                "liu_frozen_direction": direction,
                "pairs_in_frozen_direction": int((delta > 0).sum() if median_delta > 0 else (delta < 0).sum()),
                "pairs_opposite_frozen_direction": int((delta < 0).sum() if median_delta > 0 else (delta > 0).sum()),
                "paired_patients": int(len(delta)),
            }
        )
    return pd.DataFrame(difference_rows), pd.DataFrame(direction_rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--liu-h5ad", type=Path, required=True)
    parser.add_argument("--liu-usages", type=Path, required=True)
    parser.add_argument("--top-genes", type=Path, required=True)
    parser.add_argument("--eligibility-table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    inputs = [args.liu_h5ad, args.liu_usages, args.top_genes, args.eligibility_table]
    for path in inputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    signatures = load_signatures(args.top_genes)
    data = ad.read_h5ad(args.liu_h5ad)
    scores, coverage = rank_scores(data.X, data.var_names, signatures)
    scores.index = data.obs_names.astype(str)
    usage = pd.read_csv(args.liu_usages, sep="\t", index_col=0)
    usage.index = usage.index.astype(str)
    usage.columns = [int(value) for value in usage.columns]
    if set(usage.index) != set(scores.index):
        raise ValueError("Liu usage and rank-score cells do not match")
    usage = usage.loc[scores.index, PROGRAMS].astype(float)
    usage = usage.div(usage.sum(axis=1), axis=0)

    corr_rows: list[dict[str, object]] = []
    calibration_rows: list[dict[str, object]] = []
    for signature_program in PROGRAMS:
        correlations: dict[int, float] = {}
        for usage_program in PROGRAMS:
            rho = float(spearmanr(scores[signature_program], usage[usage_program]).statistic)
            correlations[usage_program] = rho
            corr_rows.append(
                {
                    "signature_program": signature_program,
                    "usage_program": usage_program,
                    "spearman_rho": rho,
                }
            )
        ranked = sorted(correlations, key=lambda p: (-correlations[p], p))
        own_rank = ranked.index(signature_program) + 1
        own_rho = correlations[signature_program]
        calibration_rows.append(
            {
                "program_id": signature_program,
                "own_usage_spearman_rho": own_rho,
                "own_usage_rank_among_10": own_rank,
                "projection_calibration_pass": own_rho >= 0.30 and own_rank == 1,
                "calibration_rule": "own rho >= 0.30 and own usage is the highest-correlated K10 usage",
            }
        )
    calibration = pd.DataFrame(calibration_rows)
    differences, directions = paired_differences(scores, data.obs)
    eligibility = pd.read_csv(args.eligibility_table, sep="\t")[["program_id", "che_eligibility", "eligibility_use"]]
    reference = eligibility.merge(calibration, on="program_id", validate="one_to_one").merge(
        directions, on="program_id", validate="one_to_one"
    )
    eligible = reference["che_eligibility"].eq("eligible_for_confirmatory_che")
    if not reference.loc[eligible, "projection_calibration_pass"].all():
        failed = reference.loc[eligible & ~reference["projection_calibration_pass"], "program_id"].tolist()
        raise RuntimeError(f"Eligible programs fail Liu projection calibration: {failed}")

    cell_path = args.output_dir / "P0_liu_k10_rank_projection_cell_scores.tsv.gz"
    corr_path = args.output_dir / "P0_liu_k10_rank_projection_cross_correlation.tsv"
    pair_path = args.output_dir / "P0_liu_k10_rank_projection_pair_differences.tsv"
    ref_path = args.output_dir / "P0_liu_k10_projection_reference.tsv"
    coverage_path = args.output_dir / "P0_liu_k10_projection_gene_coverage.tsv"
    cell_table = scores.copy()
    cell_table.insert(0, "cell_id", cell_table.index)
    cell_table.to_csv(cell_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(corr_rows).to_csv(corr_path, sep="\t", index=False)
    differences.to_csv(pair_path, sep="\t", index=False)
    reference.to_csv(ref_path, sep="\t", index=False)
    coverage.to_csv(coverage_path, sep="\t", index=False)
    receipt = {
        "status": "LIU_K10_PROJECTION_REFERENCE_FROZEN_BEFORE_CHE",
        "k": K,
        "top_n": 50,
        "rank_score_max_rank": 1500,
        "sample_summary": "mean cell rank score within patient-state",
        "paired_effect": "liver metastasis minus primary",
        "eligible_programs": reference.loc[eligible, "program_id"].astype(int).tolist(),
        "eligible_program_directions": {
            str(int(row.program_id)): row.liu_frozen_direction
            for row in reference.loc[eligible].itertuples(index=False)
        },
        "che_data_read_or_scored": False,
        "inputs_sha256": {str(path): sha256_file(path) for path in inputs},
        "outputs": [str(cell_path), str(corr_path), str(pair_path), str(ref_path), str(coverage_path)],
    }
    receipt_path = args.output_dir / "P0_liu_k10_projection_reference_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, ensure_ascii=False))
    data.file.close() if getattr(data, "isbacked", False) else None
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
