#!/usr/bin/env python3
"""Run the frozen Che external projection without re-learning Liu programs."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from formal_cnmf_design import (
    ONE_METHOD_CLASS,
    SELECTED_CLASS,
    balance_patient_state_cells,
    load_evidence_tables,
    select_cohort_cells,
)
from freeze_liu_k10_projection_reference import PROGRAMS, load_signatures, rank_scores, sha256_file
from join_crc_h5ad_table_s1 import canonical_patient
from prepare_crc_cna_inputs import unique_gene_features


ELIGIBLE = {1, 7, 8}


def finalize_replication_status(
    *,
    primary_status: str,
    eligible_for_confirmatory_che: bool,
    sensitivity_direction_stable: bool,
) -> str:
    """Apply sensitivity only after the primary analysis supports Liu's direction.

    A program that already fails the primary external check remains
    ``not_replicated_or_opposite``.  It must not be relabelled as merely
    sensitivity-unstable because that wording would conceal the primary
    non-replication.
    """
    if (
        eligible_for_confirmatory_che
        and primary_status != "not_replicated_or_opposite"
        and not sensitivity_direction_stable
    ):
        return "sensitivity_unstable"
    return primary_status


def exact_sign_flip_pvalue(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values != 0)]
    if not len(values):
        return 1.0
    observed = abs(values.mean())
    statistics = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(values)):
        statistics.append(abs(np.mean(values * np.asarray(signs))))
    return float(np.mean(np.asarray(statistics) >= observed - 1e-15))


def paired_hodges_lehmann(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    walsh = [(values[i] + values[j]) / 2.0 for i in range(len(values)) for j in range(i, len(values))]
    return float(np.median(walsh))


def holm_adjust(pvalues: pd.Series) -> pd.Series:
    values = pvalues.to_numpy(dtype=float)
    order = np.argsort(values, kind="stable")
    adjusted = np.empty_like(values)
    running = 0.0
    m = len(values)
    for rank, position in enumerate(order):
        running = max(running, (m - rank) * values[position])
        adjusted[position] = min(1.0, running)
    return pd.Series(adjusted, index=pvalues.index)


def extract_counts(
    h5ad_path: Path,
    selected: pd.DataFrame,
    *,
    symbol_key: str,
) -> tuple[ad.AnnData, pd.DataFrame]:
    source = ad.read_h5ad(h5ad_path, backed="r")
    try:
        if source.raw is None:
            raise ValueError("Che projection requires raw integer counts")
        rows = source.obs_names.get_indexer(selected["cell_id"].astype(str))
        if (rows < 0).any():
            raise ValueError(f"Selected Che cells missing from H5AD: {selected.loc[rows < 0, 'cell_id'].head().tolist()}")
        order = np.argsort(rows, kind="stable")
        rows = rows[order]
        selected = selected.iloc[order].reset_index(drop=True)
        checks = {
            "sample_id": selected["panel_sample_id"].astype(str).to_numpy(),
            "dataset": selected["panel_dataset"].astype(str).to_numpy(),
            "sample_type": selected["panel_sample_type"].astype(str).to_numpy(),
            "tissue": selected["panel_tissue"].astype(str).to_numpy(),
        }
        for column, expected in checks.items():
            observed = source.obs.iloc[rows][column].astype(str).to_numpy()
            if not np.array_equal(observed, expected):
                raise ValueError(f"Che CNA evidence and H5AD disagree for {column}")
        observed_patients = source.obs.iloc[rows]["donor_id"].map(canonical_patient).to_numpy()
        expected_patients = selected["panel_patient_id"].map(canonical_patient).to_numpy()
        if not np.array_equal(observed_patients, expected_patients):
            raise ValueError("Che CNA evidence and H5AD patient IDs disagree")
        matrix = source.raw.X[rows, :]
        if not sparse.issparse(matrix):
            matrix = sparse.csr_matrix(np.asarray(matrix))
        feature_map = unique_gene_features(source.raw.var.copy(), symbol_key)
        matrix = matrix.tocsr()[:, feature_map["feature_index"].to_numpy(dtype=np.int64)].tocsr()
        if matrix.nnz and (
            not np.isfinite(matrix.data).all()
            or (matrix.data < 0).any()
            or not np.allclose(matrix.data, np.rint(matrix.data), atol=1e-6)
        ):
            raise ValueError("Che raw expression is not finite nonnegative integer counts")
        matrix.data = np.rint(matrix.data).astype(np.float32)
        matrix.eliminate_zeros()
        obs = selected.copy()
        obs.index = obs["cell_id"].astype(str)
        obs.index.name = None
        obs["analysis_patient_id"] = obs["panel_patient_id"].map(canonical_patient)
        obs["sample_type"] = obs["panel_sample_type"].astype(str)
        obs["sample_id"] = obs["panel_sample_id"].astype(str)
        obs["dataset"] = obs["panel_dataset"].astype(str)
        var = pd.DataFrame(index=feature_map["gene"].astype(str).to_numpy())
        var.index.name = "gene"
        data = ad.AnnData(X=matrix, obs=obs, var=var)
    finally:
        source.file.close()
    return data, feature_map


def summarize_scores(
    scores: pd.DataFrame,
    obs: pd.DataFrame,
    liu_reference: pd.DataFrame,
    analysis_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cell = scores.copy()
    cell.index = obs.index.astype(str)
    cell = cell.join(obs[["analysis_patient_id", "sample_type", "sample_id"]].astype(str))
    patient_state = (
        cell.groupby(["analysis_patient_id", "sample_type"], observed=True)[PROGRAMS]
        .mean()
        .reset_index()
    )
    wide = patient_state.pivot(index="analysis_patient_id", columns="sample_type", values=PROGRAMS)
    if set(wide.columns.get_level_values(1)) != {"tumor", "metastasis"}:
        raise ValueError("Che projection lost a paired state")
    pair_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    direction_map = liu_reference.set_index("program_id")["liu_frozen_direction"].to_dict()
    eligibility_map = liu_reference.set_index("program_id")["che_eligibility"].to_dict()
    use_map = liu_reference.set_index("program_id")["eligibility_use"].to_dict()
    for program in PROGRAMS:
        primary = wide[(program, "tumor")]
        metastasis = wide[(program, "metastasis")]
        delta = metastasis - primary
        frozen_direction = direction_map[int(program)]
        expected_positive = frozen_direction == "higher_in_liver_metastasis"
        same_direction = int((delta > 0).sum() if expected_positive else (delta < 0).sum())
        median_delta = float(np.median(delta))
        direction_consistent = median_delta > 0 if expected_positive else median_delta < 0
        for patient in delta.index:
            pair_rows.append(
                {
                    "analysis": analysis_name,
                    "program_id": program,
                    "analysis_patient_id": patient,
                    "primary_mean_score": float(primary.loc[patient]),
                    "liver_metastasis_mean_score": float(metastasis.loc[patient]),
                    "delta_lm_minus_primary": float(delta.loc[patient]),
                    "consistent_with_liu_direction": bool(
                        delta.loc[patient] > 0 if expected_positive else delta.loc[patient] < 0
                    ),
                }
            )
        if direction_consistent and same_direction >= 4:
            status = "directionally_replicated_underpowered"
        elif direction_consistent and same_direction == 3:
            status = "directionally_consistent_weak"
        else:
            status = "not_replicated_or_opposite"
        summary_rows.append(
            {
                "analysis": analysis_name,
                "program_id": program,
                "che_eligibility": eligibility_map[int(program)],
                "eligibility_use": use_map[int(program)],
                "liu_frozen_direction": frozen_direction,
                "che_median_delta_lm_minus_primary": median_delta,
                "che_hodges_lehmann_delta": paired_hodges_lehmann(delta.to_numpy()),
                "che_exact_sign_flip_p_two_sided": exact_sign_flip_pvalue(delta.to_numpy()),
                "pairs_in_liu_direction": same_direction,
                "pairs_opposite_liu_direction": int(len(delta) - same_direction),
                "paired_patients": int(len(delta)),
                "direction_consistent_with_liu": direction_consistent,
                "replication_status_before_sensitivity": status,
            }
        )
    summary = pd.DataFrame(summary_rows)
    return patient_state, pd.DataFrame(pair_rows), summary


def run_analysis(
    *,
    name: str,
    evidence: pd.DataFrame,
    h5ad_path: Path,
    signatures: dict[int, list[str]],
    liu_reference: pd.DataFrame,
    selection_classes: tuple[str, ...],
    cap: int | None,
    seed: int,
    symbol_key: str,
    output_dir: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    expected_cells = 6115 if selection_classes == (SELECTED_CLASS,) else 6984
    selected, sample_qc = select_cohort_cells(
        evidence,
        "Che_2021",
        expected_patients=5,
        expected_samples=10,
        expected_dual_cells=expected_cells,
        selected_classes=selection_classes,
    )
    if cap is None:
        chosen = selected.sort_values(["panel_patient_id", "panel_sample_type", "cell_id"], kind="mergesort")
        patient_state = (
            selected.groupby(["panel_patient_id", "panel_sample_type"], observed=True)
            .size().rename("available_and_selected_cells").reset_index()
        )
    else:
        chosen, patient_state, _sample_comp = balance_patient_state_cells(selected, cap, seed)
    data, feature_map = extract_counts(h5ad_path, chosen, symbol_key=symbol_key)
    scores, coverage = rank_scores(data.X, data.var_names, signatures, max_rank=1500)
    scores.index = data.obs_names.astype(str)
    patient_scores, pair_scores, summary = summarize_scores(scores, data.obs, liu_reference, name)

    analysis_dir = output_dir / name
    analysis_dir.mkdir(parents=True, exist_ok=True)
    cell_table = scores.copy()
    cell_table.insert(0, "cell_id", cell_table.index)
    cell_table = cell_table.join(
        data.obs[["analysis_patient_id", "sample_type", "sample_id", "cna_selection_class"]]
    )
    cell_table.to_csv(analysis_dir / "P0_che_k10_cell_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    patient_scores.to_csv(analysis_dir / "P0_che_k10_patient_state_scores.tsv", sep="\t", index=False)
    pair_scores.to_csv(analysis_dir / "P0_che_k10_pair_differences.tsv", sep="\t", index=False)
    summary.to_csv(analysis_dir / "P0_che_k10_program_summary.tsv", sep="\t", index=False)
    chosen.to_csv(analysis_dir / "P0_che_k10_selected_cells.tsv.gz", sep="\t", index=False, compression="gzip")
    patient_state.to_csv(analysis_dir / "P0_che_k10_patient_state_cell_counts.tsv", sep="\t", index=False)
    sample_qc.to_csv(analysis_dir / "P0_che_k10_reference_qc.tsv", sep="\t", index=False)
    coverage.to_csv(analysis_dir / "P0_che_k10_signature_coverage.tsv", sep="\t", index=False)
    feature_map.to_csv(analysis_dir / "P0_che_k10_gene_feature_map.tsv", sep="\t", index=False)
    receipt = {
        "analysis": name,
        "selection_classes": list(selection_classes),
        "patient_state_cap": cap,
        "seed": seed,
        "selected_cells": int(data.n_obs),
        "patients": int(data.obs["analysis_patient_id"].nunique()),
        "samples": int(data.obs["sample_id"].nunique()),
        "genes": int(data.n_vars),
        "minimum_signature_coverage": int(coverage["present_genes"].min()),
    }
    return summary, receipt


def plot_confirmatory(pair_table: pd.DataFrame, final_summary: pd.DataFrame, output_path: Path) -> None:
    programs = [1, 7, 8]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8), sharey=False)
    for axis, program in zip(axes, programs):
        subset = pair_table[(pair_table["analysis"] == "primary_dual_balanced") & (pair_table["program_id"] == program)]
        for row in subset.itertuples(index=False):
            axis.plot([0, 1], [row.primary_mean_score, row.liver_metastasis_mean_score], color="#888888", alpha=0.75)
            axis.scatter([0, 1], [row.primary_mean_score, row.liver_metastasis_mean_score], s=28, color=["#4472C4", "#C55A11"])
        status = final_summary.loc[final_summary["program_id"] == program, "final_replication_status"].iloc[0]
        axis.set_title(f"P{program}\n{status}", fontsize=9)
        axis.set_xticks([0, 1], ["Primary", "Liver met."])
        axis.set_ylabel("Mean rank score" if program == 1 else "")
        axis.grid(axis="y", alpha=0.2)
    fig.suptitle("Che external projection: frozen Liu K=10 confirmatory programs")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--h5ad-sha256", required=True)
    parser.add_argument("--cna-evidence", type=Path, action="append", required=True)
    parser.add_argument("--top-genes", type=Path, required=True)
    parser.add_argument("--liu-reference", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--cap", type=int, default=500)
    parser.add_argument("--symbol-key", default="GeneSymbol")
    args = parser.parse_args()
    inputs = [args.h5ad, args.top_genes, args.liu_reference, args.contract, *args.cna_evidence]
    for path in inputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    expected_sha = args.h5ad_sha256.upper()
    if sha256_file(args.h5ad).upper() != expected_sha:
        raise RuntimeError("CRC atlas H5AD SHA256 mismatch")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    evidence = load_evidence_tables(args.cna_evidence)
    signatures = load_signatures(args.top_genes)
    liu_reference = pd.read_csv(args.liu_reference, sep="\t")
    if set(liu_reference.loc[liu_reference["che_eligibility"] == "eligible_for_confirmatory_che", "program_id"]) != ELIGIBLE:
        raise ValueError("Frozen eligible programs are not P1/P7/P8")
    if not liu_reference.loc[liu_reference["program_id"].isin(ELIGIBLE), "projection_calibration_pass"].all():
        raise ValueError("An eligible program lacks Liu projection calibration")

    specifications = [
        ("primary_dual_balanced", (SELECTED_CLASS,), args.cap),
        ("sensitivity_dual_all_cells", (SELECTED_CLASS,), None),
        ("sensitivity_dual_plus_one_method_balanced", (SELECTED_CLASS, ONE_METHOD_CLASS), args.cap),
    ]
    summaries: list[pd.DataFrame] = []
    receipts: list[dict[str, object]] = []
    for name, classes, cap in specifications:
        summary, receipt = run_analysis(
            name=name,
            evidence=evidence,
            h5ad_path=args.h5ad,
            signatures=signatures,
            liu_reference=liu_reference,
            selection_classes=classes,
            cap=cap,
            seed=args.seed,
            symbol_key=args.symbol_key,
            output_dir=args.output_dir,
        )
        summaries.append(summary)
        receipts.append(receipt)
    combined = pd.concat(summaries, ignore_index=True)
    combined.to_csv(args.output_dir / "P0_che_k10_all_analysis_summary.tsv", sep="\t", index=False)
    primary = combined[combined["analysis"] == "primary_dual_balanced"].copy()
    sensitivity = combined[combined["analysis"] != "primary_dual_balanced"]
    final_rows: list[dict[str, object]] = []
    for row in primary.itertuples(index=False):
        program_sensitivity = sensitivity[sensitivity["program_id"] == row.program_id]
        sensitivity_direction_stable = bool(program_sensitivity["direction_consistent_with_liu"].all())
        final_status = finalize_replication_status(
            primary_status=row.replication_status_before_sensitivity,
            eligible_for_confirmatory_che=row.program_id in ELIGIBLE,
            sensitivity_direction_stable=sensitivity_direction_stable,
        )
        final_rows.append(
            {
                **row._asdict(),
                "sensitivity_direction_stable": sensitivity_direction_stable,
                "final_replication_status": final_status,
            }
        )
    final = pd.DataFrame(final_rows)
    confirmatory_mask = final["program_id"].isin(ELIGIBLE)
    final["confirmatory_holm_p"] = np.nan
    final.loc[confirmatory_mask, "confirmatory_holm_p"] = holm_adjust(
        final.loc[confirmatory_mask, "che_exact_sign_flip_p_two_sided"]
    )
    final.to_csv(args.output_dir / "P0_che_k10_confirmatory_and_exploratory_results.tsv", sep="\t", index=False)
    pair_tables = [
        pd.read_csv(args.output_dir / name / "P0_che_k10_pair_differences.tsv", sep="\t")
        for name, _, _ in specifications
    ]
    all_pairs = pd.concat(pair_tables, ignore_index=True)
    all_pairs.to_csv(args.output_dir / "P0_che_k10_all_pair_differences.tsv", sep="\t", index=False)
    plot_confirmatory(all_pairs, final, args.output_dir / "P0_che_k10_confirmatory_paired_plot.png")
    receipt = {
        "status": "READY_FOR_CHE_EXTERNAL_REVIEW",
        "frozen_eligible_programs": sorted(ELIGIBLE),
        "programs_scored_blind": PROGRAMS,
        "liu_directions": {
            str(int(row.program_id)): row.liu_frozen_direction
            for row in liu_reference[liu_reference["program_id"].isin(ELIGIBLE)].itertuples(index=False)
        },
        "analysis_receipts": receipts,
        "h5ad": str(args.h5ad),
        "h5ad_sha256": expected_sha,
        "input_sha256": {str(path): sha256_file(path) for path in inputs},
        "che_results_used_to_change_k_or_eligibility": False,
        "biological_naming_executed": False,
        "final_confirmatory_status": {
            str(int(row.program_id)): row.final_replication_status
            for row in final[final["program_id"].isin(ELIGIBLE)].itertuples(index=False)
        },
    }
    (args.output_dir / "P0_che_k10_external_projection_receipt.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
