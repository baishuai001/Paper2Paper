#!/usr/bin/env python3
"""Compare cNMF programs across seeds and adjacent K with an empirical null."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


SPECTRA_PATTERN = re.compile(r"\.gene_spectra_score\.k_(\d+)\.dt_[^.]+\.txt$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes())
    return digest.hexdigest().upper()


def cosine_matrix(left: pd.DataFrame, right: pd.DataFrame) -> np.ndarray:
    genes = left.columns.intersection(right.columns)
    if len(genes) < 2:
        raise ValueError("fewer than two shared genes between cNMF spectra")
    a = left.loc[:, genes].to_numpy(dtype=float)
    b = right.loc[:, genes].to_numpy(dtype=float)
    a /= np.where(np.linalg.norm(a, axis=1, keepdims=True) == 0, 1.0, np.linalg.norm(a, axis=1, keepdims=True))
    b /= np.where(np.linalg.norm(b, axis=1, keepdims=True) == 0, 1.0, np.linalg.norm(b, axis=1, keepdims=True))
    return a @ b.T


def top_gene_sets(spectra: pd.DataFrame, top_n: int) -> list[set[str]]:
    return [set(row.nlargest(min(top_n, len(row))).index.astype(str)) for _, row in spectra.iterrows()]


def match_programs(left: pd.DataFrame, right: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    similarities = cosine_matrix(left, right)
    left_index, right_index = linear_sum_assignment(-similarities)
    left_top = top_gene_sets(left, top_n)
    right_top = top_gene_sets(right, top_n)
    rows = []
    for i, j in zip(left_index, right_index):
        union = left_top[i] | right_top[j]
        rows.append(
            {
                "left_program": str(left.index[i]),
                "right_program": str(right.index[j]),
                "cosine": float(similarities[i, j]),
                "top_gene_jaccard": float(len(left_top[i] & right_top[j]) / len(union)) if union else 0.0,
            }
        )
    return pd.DataFrame(rows)


def empirical_thresholds(
    left: pd.DataFrame,
    right: pd.DataFrame,
    top_n: int,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    genes = left.columns.intersection(right.columns)
    left_values = left.loc[:, genes].to_numpy(dtype=float)
    right_values = right.loc[:, genes].to_numpy(dtype=float)
    left_norms = np.linalg.norm(left_values, axis=1, keepdims=True)
    right_norms = np.linalg.norm(right_values, axis=1, keepdims=True)
    left_normalized = left_values / np.where(left_norms == 0, 1.0, left_norms)
    right_norms = np.where(right_norms == 0, 1.0, right_norms)
    effective_top_n = min(top_n, len(genes))
    left_top = [
        set(np.argpartition(-row, effective_top_n - 1)[:effective_top_n])
        for row in left_values
    ]
    cosine_null: list[float] = []
    jaccard_null: list[float] = []
    for _ in range(permutations):
        permuted = np.vstack([rng.permutation(row) for row in right_values])
        similarities = left_normalized @ (permuted / right_norms).T
        left_index, right_index = linear_sum_assignment(-similarities)
        right_top = [
            set(np.argpartition(-row, effective_top_n - 1)[:effective_top_n])
            for row in permuted
        ]
        for i, j in zip(left_index, right_index):
            intersection = len(left_top[i] & right_top[j])
            union = len(left_top[i] | right_top[j])
            cosine_null.append(float(similarities[i, j]))
            jaccard_null.append(float(intersection / union) if union else 0.0)
    return float(np.quantile(cosine_null, 0.99)), float(np.quantile(jaccard_null, 0.99))


def discover_spectra(run_root: Path, run_name: str) -> dict[int, pd.DataFrame]:
    run_dir = run_root / run_name
    if not run_dir.is_dir():
        raise FileNotFoundError(f"cNMF run directory missing: {run_dir}")
    discovered: dict[int, pd.DataFrame] = {}
    for path in sorted(run_dir.glob(f"{run_name}.gene_spectra_score.k_*.txt")):
        match = SPECTRA_PATTERN.search(path.name)
        if match:
            k = int(match.group(1))
            spectra = pd.read_csv(path, sep="\t", index_col=0)
            # cNMF 1.7.1 writes gene_spectra_score as K program rows by gene
            # columns.  A second transpose turns genes into apparent programs,
            # makes the Hungarian problem enormous, and invalidates matching.
            if spectra.shape[0] != k:
                raise ValueError(
                    f"cNMF spectra orientation/shape mismatch for K={k}: "
                    f"expected {k} program rows, got {spectra.shape}"
                )
            if spectra.shape[1] <= k:
                raise ValueError(f"cNMF spectra for K={k} contain too few gene columns: {spectra.shape}")
            if spectra.index.has_duplicates or spectra.columns.has_duplicates:
                raise ValueError(f"cNMF spectra for K={k} have duplicate program or gene identifiers")
            numeric = spectra.apply(pd.to_numeric, errors="coerce")
            if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
                raise ValueError(f"cNMF spectra for K={k} contain nonnumeric or nonfinite values")
            discovered[k] = numeric
    if not discovered:
        raise FileNotFoundError(f"no consensus spectra found under {run_dir}")
    return discovered


def cross_run_comparison_type(left_label: str, right_label: str) -> str:
    labels = {left_label, right_label}
    if all(label.startswith("seed_") for label in labels):
        return "same_k_cross_seed"
    if labels == {"patient_holdout_A", "patient_holdout_B"}:
        return "same_k_independent_patient_holdout"
    if any(label.startswith("cell_resample_") for label in labels) and any(
        label.startswith("seed_") for label in labels
    ):
        return "same_k_cell_resample"
    if labels == {"holdout_A", "holdout_B"}:
        return "same_k_independent_patient_dataset_holdout"
    return "same_k_cross_input"


def evaluation_status(has_patient_holdout: bool, has_cell_resample: bool = False) -> str:
    if has_patient_holdout and has_cell_resample:
        return "candidate_seed_cell_resample_and_patient_holdout_comparison"
    if has_patient_holdout:
        return "diagnostic_seed_adjacent_k_and_independent_holdout_comparison"
    return "diagnostic_seed_and_adjacent_k_comparison_only"


def summarize_independent_holdout_pairs(
    matches: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep direct, same-K A/B evidence without inventing cross-K families.

    Transitive graph union across adjacent K can merge two distinct programs from
    the same run and K through intermediate nodes.  Until K is frozen, the valid
    unit of holdout replication is therefore the one-to-one matched pair at each
    K, not a graph-derived meta-program family.
    """

    holdout_types = {
        "same_k_independent_patient_holdout",
        "same_k_independent_patient_dataset_holdout",
    }
    pairs = matches.loc[matches["comparison_type"].isin(holdout_types)].copy()
    if pairs.empty:
        return pairs, pd.DataFrame(
            columns=[
                "k",
                "matched_pairs",
                "pairs_exceeding_both_nulls",
                "replication_fraction",
                "median_cosine",
                "median_top_gene_jaccard",
                "cosine_null_q99",
                "jaccard_null_q99",
            ]
        )
    if not (pairs["left_k"].astype(int) == pairs["right_k"].astype(int)).all():
        raise ValueError("independent holdout replication must compare the same K")
    if pairs["comparison_type"].nunique() != 1:
        raise ValueError("cannot merge patient-only and patient-dataset holdout designs")
    pairs["k"] = pairs["left_k"].astype(int)
    pairs = pairs.sort_values(
        ["k", "left_program", "right_program"], kind="mergesort"
    ).reset_index(drop=True)
    within_k_order = pairs.groupby("k", sort=False).cumcount() + 1
    pairs.insert(
        0,
        "holdout_pair_id",
        [f"K{k}-M{order:02d}" for k, order in zip(pairs["k"], within_k_order)],
    )
    by_k = (
        pairs.groupby("k", as_index=False, sort=True)
        .agg(
            matched_pairs=("exceeds_both_nulls", "size"),
            pairs_exceeding_both_nulls=("exceeds_both_nulls", "sum"),
            median_cosine=("cosine", "median"),
            median_top_gene_jaccard=("top_gene_jaccard", "median"),
            cosine_null_q99=("cosine_null_q99", "first"),
            jaccard_null_q99=("jaccard_null_q99", "first"),
        )
    )
    by_k["pairs_exceeding_both_nulls"] = by_k[
        "pairs_exceeding_both_nulls"
    ].astype(int)
    by_k.insert(
        3,
        "replication_fraction",
        by_k["pairs_exceeding_both_nulls"] / by_k["matched_pairs"],
    )
    return pairs, by_k


def evaluate(
    run_specs: list[tuple[str, Path, str]],
    output_dir: Path,
    top_n: int,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    if len(run_specs) < 2:
        raise ValueError("at least two cNMF runs are required")
    runs = {label: discover_spectra(root, name) for label, root, name in run_specs}
    rng = np.random.default_rng(seed)
    comparisons: list[pd.DataFrame] = []
    threshold_rows: list[dict[str, object]] = []
    labels = list(runs)
    comparison_specs: list[tuple[str, int, str, int, str]] = []
    for i, left_label in enumerate(labels):
        for right_label in labels[i + 1 :]:
            for k in sorted(set(runs[left_label]) & set(runs[right_label])):
                comparison_specs.append(
                    (
                        left_label,
                        k,
                        right_label,
                        k,
                        cross_run_comparison_type(left_label, right_label),
                    )
                )
    for label in labels:
        ks = sorted(runs[label])
        for left_k, right_k in zip(ks, ks[1:]):
            if right_k - left_k == 1:
                comparison_specs.append((label, left_k, label, right_k, "adjacent_k_same_seed"))
    for left_label, left_k, right_label, right_k, comparison_type in comparison_specs:
        left = runs[left_label][left_k]
        right = runs[right_label][right_k]
        cosine_threshold, jaccard_threshold = empirical_thresholds(
            left, right, top_n, permutations, rng
        )
        matched = match_programs(left, right, top_n=top_n)
        matched.insert(0, "comparison_type", comparison_type)
        matched.insert(1, "left_run", left_label)
        matched.insert(2, "left_k", left_k)
        matched.insert(3, "right_run", right_label)
        matched.insert(4, "right_k", right_k)
        matched["cosine_null_q99"] = cosine_threshold
        matched["jaccard_null_q99"] = jaccard_threshold
        matched["exceeds_both_nulls"] = (
            (matched["cosine"] > cosine_threshold)
            & (matched["top_gene_jaccard"] > jaccard_threshold)
        )
        comparisons.append(matched)
        threshold_rows.append(
            {
                "comparison_type": comparison_type,
                "left_run": left_label,
                "left_k": left_k,
                "right_run": right_label,
                "right_k": right_k,
                "permutations": permutations,
                "cosine_null_q99": cosine_threshold,
                "jaccard_null_q99": jaccard_threshold,
            }
        )
    matches = pd.concat(comparisons, ignore_index=True)
    thresholds = pd.DataFrame(threshold_rows)

    holdout_comparisons, holdout_by_k = summarize_independent_holdout_pairs(matches)

    output_dir.mkdir(parents=True, exist_ok=True)
    match_path = output_dir / "P0_cnmf_stability.tsv"
    threshold_path = output_dir / "P0_cnmf_null_thresholds.tsv"
    holdout_pairs_path = output_dir / "P0_cnmf_holdout_replication_pairs.tsv"
    holdout_by_k_path = output_dir / "P0_cnmf_holdout_replication_by_k.tsv"
    matches.to_csv(match_path, sep="\t", index=False, lineterminator="\n")
    thresholds.to_csv(threshold_path, sep="\t", index=False, lineterminator="\n")
    holdout_comparisons.to_csv(
        holdout_pairs_path, sep="\t", index=False, lineterminator="\n"
    )
    holdout_by_k.to_csv(
        holdout_by_k_path, sep="\t", index=False, lineterminator="\n"
    )
    has_patient_holdout = bool(
        matches["comparison_type"].eq("same_k_independent_patient_holdout").any()
    )
    has_patient_dataset_holdout = bool(
        matches["comparison_type"].eq("same_k_independent_patient_dataset_holdout").any()
    )
    has_any_holdout = has_patient_holdout or has_patient_dataset_holdout
    has_cell_resample = bool(
        matches["comparison_type"].eq("same_k_cell_resample").any()
    )
    if has_patient_holdout:
        claim_boundary = (
            "The formal Liu comparison includes paired-patient-disjoint A/B folds and a separately "
            "sampled cell input. Both folds share the Liu dataset, so they test patient replication, "
            "not external cohort transport. Evidence is reported as direct one-to-one pairs at each K. "
            "No cross-K family or program name is accepted until K review; Che external projection, "
            "CNA-set sensitivity and technical-program review remain required."
        )
    elif has_patient_dataset_holdout:
        claim_boundary = (
            "The run includes patient- and dataset-disjoint state-stratified A/B folds, but the folds "
            "remain selected observational CRC samples and are not randomized cohorts. Evidence is "
            "reported as direct one-to-one A/B pairs separately at each K. No cross-K program family is "
            "counted until K is frozen; sufficient initializations, balanced resampling, CNA-set "
            "sensitivity and technical-program review are also still required."
        )
    else:
        claim_boundary = (
            "These are same-input seed and adjacent-K diagnostics. Patient holdout, dataset holdout, "
            "balanced resampling and P0.5 malignant-set sensitivity have not been tested; program "
            "families are provisional and K is not frozen."
        )
    receipt = {
        "status": evaluation_status(has_any_holdout, has_cell_resample),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "runs": [{"label": label, "root": str(root.resolve()), "name": name} for label, root, name in run_specs],
        "top_n": top_n,
        "permutations": permutations,
        "seed": seed,
        "comparisons": int(len(thresholds)),
        "matched_program_pairs": int(len(matches)),
        "pairs_exceeding_both_nulls": int(matches["exceeds_both_nulls"].sum()),
        "independent_holdout_matched_pairs": int(len(holdout_comparisons)),
        "independent_holdout_pairs_exceeding_both_nulls": int(
            holdout_comparisons["exceeds_both_nulls"].sum()
        ),
        "independent_holdout_k_values_with_at_least_one_pass": int(
            (holdout_by_k["pairs_exceeding_both_nulls"] > 0).sum()
        ),
        "patient_only_holdout": has_patient_holdout,
        "patient_and_dataset_holdout": has_patient_dataset_holdout,
        "cell_resample_comparison": has_cell_resample,
        "program_family_freeze_status": "not_assessed_until_k_is_frozen",
        "claim_boundary": claim_boundary,
        "outputs": {},
    }
    for path in (match_path, threshold_path, holdout_pairs_path, holdout_by_k_path):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    (output_dir / "P0_cnmf_stability_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return receipt


def parse_run(value: str) -> tuple[str, Path, str]:
    parts = value.split("=", 1)
    if len(parts) != 2 or "/" not in parts[1].replace("\\", "/"):
        raise argparse.ArgumentTypeError("run must be LABEL=OUTPUT_ROOT/RUN_NAME")
    label, full = parts
    path = Path(full)
    return label, path.parent, path.name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="append", required=True, type=parse_run)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260810)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = evaluate(args.run, args.output_dir, args.top_n, args.permutations, args.seed)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
