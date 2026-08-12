#!/usr/bin/env python3
"""Run one pinned cNMF discovery replicate and write auditable diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


EXPECTED_CNMF_VERSION = "1.7.1"
REQUIRED_OBS = ["analysis_patient_id", "dataset"]
FORMAL_RUN_ROLES = {
    "liu_primary_seed",
    "liu_cell_resample",
    "liu_patient_holdout",
}
SENSITIVITY_RUN_ROLES = {
    "liu_reference_qc_exclusion_sensitivity",
    "liu_one_method_expansion_sensitivity",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def validate_counts_input(path: Path) -> dict[str, object]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"cNMF counts H5AD missing or empty: {path}")
    data = ad.read_h5ad(path, backed="r")
    try:
        missing = sorted(set(REQUIRED_OBS) - set(data.obs.columns))
        if missing:
            raise ValueError(f"cNMF obs fields missing: {missing}")
        if not data.obs_names.is_unique or not data.var_names.is_unique:
            raise ValueError("cNMF cell and gene identifiers must be unique")
        matrix = data.X
        if not sparse.issparse(matrix):
            if hasattr(matrix, "to_memory"):
                matrix = matrix.to_memory()
            else:
                matrix = np.asarray(matrix[:, :])
        values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
        if values.size and (not np.isfinite(values).all() or (values < 0).any()):
            raise ValueError("cNMF input contains non-finite or negative values")
        if values.size and not np.allclose(values, np.rint(values), atol=1e-6):
            raise ValueError("cNMF input must contain integer-valued counts")
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        if (totals <= 0).any():
            raise ValueError("cNMF input contains zero-count cells")
        selection_basis = "author_cancer_label_pending_cna"
        selection_classes: list[str] = []
        if "cna_selection_class" in data.obs.columns:
            selection_classes = sorted(data.obs["cna_selection_class"].astype(str).unique().tolist())
            allowed = {
                "author_cancer_two_method_malignancy_support",
                "author_cancer_one_method_malignancy_support",
            }
            if (
                "author_cancer_two_method_malignancy_support" not in selection_classes
                or not set(selection_classes).issubset(allowed)
            ):
                raise ValueError(
                    "cNMF CNA selection contains unexpected evidence classes: "
                    f"{selection_classes}"
                )
            selection_basis = (
                "copykat_or_scevan_malignancy_support_expanded"
                if "author_cancer_one_method_malignancy_support" in selection_classes
                else "copykat_aneuploid_and_scevan_tumor_intersection"
            )
        formal_analysis_roles = (
            sorted(data.obs["formal_analysis_role"].astype(str).unique().tolist())
            if "formal_analysis_role" in data.obs.columns
            else []
        )
        cell_sampling_seeds = (
            sorted(pd.to_numeric(data.obs["cell_sampling_seed"], errors="raise").astype(int).unique().tolist())
            if "cell_sampling_seed" in data.obs.columns
            else []
        )
        return {
            "cells": int(data.n_obs),
            "genes": int(data.n_vars),
            "patients": int(data.obs["analysis_patient_id"].astype(str).nunique()),
            "datasets": int(data.obs["dataset"].astype(str).nunique()),
            "states": (
                sorted(data.obs["sample_type"].astype(str).unique().tolist())
                if "sample_type" in data.obs.columns
                else []
            ),
            "selection_basis": selection_basis,
            "selection_classes": selection_classes,
            "formal_analysis_roles": formal_analysis_roles,
            "cell_sampling_seeds": cell_sampling_seeds,
            "sha256": sha256_file(path),
        }
    finally:
        data.file.close()


def load_npz_frame(path: Path) -> pd.DataFrame:
    with np.load(path, allow_pickle=True) as archive:
        return pd.DataFrame(data=archive["data"], index=archive["index"], columns=archive["columns"])


def normalized_entropy(weights: pd.Series) -> tuple[float, float, int]:
    values = weights.to_numpy(dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if not len(values):
        return 0.0, 0.0, 0
    probabilities = values / values.sum()
    entropy = float(-(probabilities * np.log(probabilities)).sum())
    normalized = float(entropy / np.log(len(probabilities))) if len(probabilities) > 1 else 0.0
    contributors = int((probabilities >= 0.05).sum())
    return normalized, float(np.exp(entropy)), contributors


def usage_mixing_table(usage: pd.DataFrame, obs: pd.DataFrame, k: int) -> pd.DataFrame:
    missing_cells = sorted(set(usage.index.astype(str)) - set(obs.index.astype(str)))
    if missing_cells:
        raise ValueError(f"cNMF usage has cells missing from input metadata: {missing_cells[:3]}")
    metadata = obs.copy()
    metadata.index = metadata.index.astype(str)
    joined = usage.copy()
    joined.index = joined.index.astype(str)
    joined = joined.join(metadata[REQUIRED_OBS], how="left")
    rows: list[dict[str, object]] = []
    for program in usage.columns:
        for unit in REQUIRED_OBS:
            unit_usage = joined.groupby(unit, observed=True)[program].mean()
            entropy, effective, contributors = normalized_entropy(unit_usage)
            rows.append(
                {
                    "k": k,
                    "program": str(program),
                    "unit": unit,
                    "units_total": int(len(unit_usage)),
                    "normalized_shannon_entropy": entropy,
                    "effective_unit_number": effective,
                    "contributors_ge_5pct": contributors,
                    "calculation": "mean cell usage within unit, then normalize contributions across units",
                }
            )
    return pd.DataFrame(rows)


def maximum_program_redundancy(spectra: pd.DataFrame) -> float:
    values = spectra.to_numpy(dtype=float)
    norms = np.linalg.norm(values, axis=0)
    normalized = values / np.where(norms == 0, 1.0, norms)
    similarities = normalized.T @ normalized
    np.fill_diagonal(similarities, np.nan)
    return float(np.nanmax(similarities)) if similarities.shape[0] > 1 else 0.0


def relabel_k_selection_density_metadata(
    k_stats: pd.DataFrame, final_consensus_density_threshold: float
) -> pd.DataFrame:
    """Make cNMF 1.7.1's no-filter K-statistics label unambiguous."""

    if "local_density_threshold" not in k_stats.columns:
        raise ValueError("cNMF K-selection statistics are missing local_density_threshold")
    observed = sorted(
        pd.to_numeric(k_stats["local_density_threshold"], errors="raise").unique().tolist()
    )
    if observed != [0.5]:
        raise ValueError(f"unexpected cNMF 1.7.1 K-selection threshold metadata: {observed}")
    result = k_stats.rename(
        columns={"local_density_threshold": "cnmf_reported_default_threshold_no_filter"}
    ).copy()
    result.insert(2, "k_selection_filtering", "disabled_by_cnmf_skip_density_stats_path")
    result.insert(3, "final_consensus_density_threshold", final_consensus_density_threshold)
    return result


def factorize_worker(output_dir: str, run_name: str, worker_i: int, total_workers: int) -> int:
    """Run the officially supported cNMF task partition in an isolated process."""

    from cnmf import cNMF

    worker = cNMF(output_dir=output_dir, name=run_name)
    worker.factorize(worker_i=worker_i, total_workers=total_workers)
    return worker_i


def run_cnmf(
    counts_h5ad: Path,
    output_dir: Path,
    run_name: str,
    components: list[int],
    n_iter: int,
    seed: int,
    density_threshold: float,
    num_highvar_genes: int,
    max_nmf_iter: int,
    total_workers: int = 1,
    run_role: str = "diagnostic",
) -> dict[str, object]:
    if sorted(set(components)) != sorted(components) or any(k < 2 for k in components):
        raise ValueError("components must be unique, sorted integers >=2")
    if n_iter < 5:
        raise ValueError("n_iter must be at least 5 so the cNMF consensus neighborhood is nonzero")
    if total_workers < 1 or total_workers > n_iter * len(components):
        raise ValueError("total_workers must be between 1 and the number of NMF tasks")
    if run_role != "diagnostic" and run_role not in FORMAL_RUN_ROLES | SENSITIVITY_RUN_ROLES:
        raise ValueError(f"unsupported cNMF run role: {run_role}")
    installed = importlib.metadata.version("cnmf")
    if installed != EXPECTED_CNMF_VERSION:
        raise RuntimeError(f"cNMF version {installed}; expected {EXPECTED_CNMF_VERSION}")
    input_receipt = validate_counts_input(counts_h5ad)
    formal_criteria = (
        run_role in FORMAL_RUN_ROLES
        and input_receipt["selection_basis"] == "copykat_aneuploid_and_scevan_tumor_intersection"
        and input_receipt["formal_analysis_roles"] == ["liu_primary_discovery"]
        and input_receipt["states"] == ["metastasis", "tumor"]
        and components == list(range(5, 21))
        and n_iter >= 100
    )
    sensitivity_criteria = (
        run_role in SENSITIVITY_RUN_ROLES
        and input_receipt["states"] == ["metastasis", "tumor"]
        and components == [9, 10, 11]
        and n_iter >= 100
        and (
            (
                run_role == "liu_reference_qc_exclusion_sensitivity"
                and input_receipt["selection_basis"]
                == "copykat_aneuploid_and_scevan_tumor_intersection"
                and input_receipt["formal_analysis_roles"]
                == ["liu_reference_qc_exclusion_sensitivity"]
            )
            or (
                run_role == "liu_one_method_expansion_sensitivity"
                and input_receipt["selection_basis"]
                == "copykat_or_scevan_malignancy_support_expanded"
                and input_receipt["formal_analysis_roles"]
                == ["liu_one_method_expansion_sensitivity"]
            )
        )
    )
    if run_role in FORMAL_RUN_ROLES and not formal_criteria:
        raise ValueError(
            "formal cNMF run does not satisfy the frozen Liu input, state, K-range and initialization contract"
        )
    if run_role in SENSITIVITY_RUN_ROLES and not sensitivity_criteria:
        raise ValueError(
            "cNMF sensitivity run does not satisfy the frozen role, K=9-11 and initialization contract"
        )
    run_dir = output_dir / run_name
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"refusing to mix with existing cNMF run: {run_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    from cnmf import cNMF

    cnmf_object = cNMF(output_dir=str(output_dir), name=run_name)
    cnmf_object.prepare(
        counts_fn=str(counts_h5ad),
        components=np.asarray(components, dtype=int),
        n_iter=n_iter,
        seed=seed,
        num_highvar_genes=num_highvar_genes,
        max_NMF_iter=max_nmf_iter,
    )
    if total_workers == 1:
        cnmf_object.factorize(worker_i=0, total_workers=1)
    else:
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=total_workers, mp_context=context) as pool:
            futures = {
                pool.submit(
                    factorize_worker,
                    str(output_dir),
                    run_name,
                    worker_i,
                    total_workers,
                ): worker_i
                for worker_i in range(total_workers)
            }
            completed: list[int] = []
            for future in as_completed(futures):
                completed.append(future.result())
        if sorted(completed) != list(range(total_workers)):
            raise RuntimeError(f"cNMF worker completion mismatch: {sorted(completed)}")
    cnmf_object.combine()

    # cNMF 1.7.1's k_selection_plot() always computes its statistics through
    # consensus(..., skip_density_and_return_after_stats=True).  In that code
    # path no spectra are filtered, but the returned table still records the
    # consensus() default (0.5) as ``local_density_threshold``.  Relabel that
    # value with the actual interpretation instead of allowing reviewers to
    # mistake it for the frozen 0.1 threshold used below for final consensus.
    cnmf_object.k_selection_plot(close_fig=True)

    input_data = ad.read_h5ad(counts_h5ad)
    k_stats = load_npz_frame(Path(cnmf_object.paths["k_selection_stats"])).reset_index(drop=True)
    k_stats = relabel_k_selection_density_metadata(k_stats, density_threshold)
    mixing_frames: list[pd.DataFrame] = []
    top_gene_frames: list[pd.DataFrame] = []
    redundancies: dict[int, float] = {}
    for k in components:
        cnmf_object.consensus(
            k=k,
            density_threshold=density_threshold,
            show_clustering=False,
            close_clustergram_fig=True,
        )
        usage, spectra_scores, _spectra_tpm, top_genes = cnmf_object.load_results(
            K=k,
            density_threshold=density_threshold,
            n_top_genes=100,
        )
        mixing_frames.append(usage_mixing_table(usage, input_data.obs, k))
        long_top = top_genes.reset_index(names="rank").melt(id_vars="rank", var_name="program", value_name="gene")
        long_top.insert(0, "k", k)
        long_top["rank"] = long_top["rank"].astype(int) + 1
        top_gene_frames.append(long_top)
        redundancies[k] = maximum_program_redundancy(spectra_scores)

    k_stats["maximum_within_k_program_cosine"] = k_stats["k"].astype(int).map(redundancies)
    k_stats_path = output_dir / "P0_cnmf_k_selection.tsv"
    mixing_path = output_dir / "P0_cnmf_patient_dataset_mixing.tsv"
    top_path = output_dir / "P0_cnmf_top_genes.tsv"
    k_stats.to_csv(k_stats_path, sep="\t", index=False, lineterminator="\n")
    pd.concat(mixing_frames, ignore_index=True).to_csv(mixing_path, sep="\t", index=False, lineterminator="\n")
    pd.concat(top_gene_frames, ignore_index=True).to_csv(top_path, sep="\t", index=False, lineterminator="\n")

    if formal_criteria:
        status = "formal_liu_candidate_run_completed_pending_k_and_replication_review"
        claim_boundary = (
            "This run satisfies the frozen Liu K=5-20 and >=100-initialization candidate contract. "
            "It does not itself freeze K or name a program. Cross-NMF-seed, cell-resampling and paired-"
            "patient-holdout comparisons must agree, followed by Che external projection."
        )
    elif sensitivity_criteria:
        status = "formal_liu_sensitivity_run_completed_pending_k_review"
        claim_boundary = (
            "This run is one predeclared Liu sensitivity at K=9-11 with 100 initializations. "
            "It can only be compared with the frozen primary solution before K review; it does not "
            "freeze K, name programs, or authorize Che or downstream analyses."
        )
    elif input_receipt["selection_basis"] == "copykat_aneuploid_and_scevan_tumor_intersection":
        status = "diagnostic_real_data_smoke_not_scientific_acceptance"
        claim_boundary = (
            "The input is a sampled two-method CNA-supported candidate set, not DNA-validated "
            "malignancy truth or an atlas-wide census. Patient/dataset holdouts, balanced resampling "
            "and 100-200 initializations are still required before freezing K or malignant meta-programs."
        )
    else:
        status = "diagnostic_real_data_smoke_not_scientific_acceptance"
        claim_boundary = (
            "The input uses author Cancer cell labels, not the final P0.5 CNA-supported set. This run "
            "lacks patient/dataset holdouts and 100-200 initializations, so it tests executability and "
            "produces diagnostics but cannot freeze K or malignant meta-programs."
        )
    receipt = {
        "status": status,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "cnmf_version": installed,
        "counts_h5ad": str(counts_h5ad.resolve()),
        "input": input_receipt,
        "run_name": run_name,
        "run_role": run_role,
        "components": components,
        "n_iter": n_iter,
        "seed": seed,
        "density_threshold": density_threshold,
        "num_highvar_genes": num_highvar_genes,
        "max_nmf_iter": max_nmf_iter,
        "total_workers": total_workers,
        "factorization_execution": (
            "official cNMF task partition across isolated worker processes"
            if total_workers > 1
            else "single official cNMF worker"
        ),
        "input_processing": "official cNMF TPM normalization, HVG selection and nonnegative variance scaling; no cross-dataset batch correction in this diagnostic run",
        "claim_boundary": claim_boundary,
        "outputs": {},
    }
    for path in (k_stats_path, mixing_path, top_path, Path(cnmf_object.paths["k_selection_plot"])):
        receipt["outputs"][str(path.relative_to(output_dir))] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    receipt_path = output_dir / "P0_cnmf_run_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts-h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--components", type=int, nargs="+", default=[5, 6, 7])
    parser.add_argument("--n-iter", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--density-threshold", type=float, default=0.1)
    parser.add_argument("--num-highvar-genes", type=int, default=2000)
    parser.add_argument("--max-nmf-iter", type=int, default=1000)
    parser.add_argument("--total-workers", type=int, default=1)
    parser.add_argument(
        "--run-role",
        choices=["diagnostic", *sorted(FORMAL_RUN_ROLES | SENSITIVITY_RUN_ROLES)],
        default="diagnostic",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = run_cnmf(
        args.counts_h5ad,
        args.output_dir,
        args.run_name,
        args.components,
        args.n_iter,
        args.seed,
        args.density_threshold,
        args.num_highvar_genes,
        args.max_nmf_iter,
        args.total_workers,
        args.run_role,
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
