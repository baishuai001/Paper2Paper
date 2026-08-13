#!/usr/bin/env python3
"""Study-stratified TF meta-analysis and leakage-safe LOSO validation."""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from common import bh_fdr, output_manifest, write_json
from phenotype import PhenotypeManifest


def read_activity(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", compression="infer")
    if table.columns[0] != "TF":
        raise ValueError("VIPER activity first column must be TF")
    if table["TF"].duplicated().any() or table.columns[1:].duplicated().any():
        raise ValueError("Duplicate TF or patient in VIPER activity")
    return table.set_index("TF").T


def hedges_g(case: np.ndarray, control: np.ndarray) -> tuple[float, float]:
    case = np.asarray(case, dtype=float)
    control = np.asarray(control, dtype=float)
    n1, n0 = len(case), len(control)
    pooled_variance = (
        (n1 - 1) * case.var(ddof=1) + (n0 - 1) * control.var(ddof=1)
    ) / (n1 + n0 - 2)
    if not np.isfinite(pooled_variance) or pooled_variance <= 0:
        return np.nan, np.nan
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    effect = correction * (case.mean() - control.mean()) / np.sqrt(pooled_variance)
    variance = (n1 + n0) / (n1 * n0) + effect * effect / (2.0 * (n1 + n0 - 2))
    return float(effect), float(variance)


def reml_meta(effects: np.ndarray, variances: np.ndarray) -> dict[str, float]:
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)
    valid = np.isfinite(effects) & np.isfinite(variances) & (variances > 0)
    effects, variances = effects[valid], variances[valid]
    if len(effects) < 2:
        return {key: np.nan for key in ("effect", "SE", "p", "tau2", "I2", "mKH_scale", "df")}

    fixed_weights = 1.0 / variances
    fixed = np.sum(fixed_weights * effects) / np.sum(fixed_weights)
    q = float(np.sum(fixed_weights * (effects - fixed) ** 2))
    df = len(effects) - 1
    i2 = max(0.0, (q - df) / q) * 100.0 if q > 0 else 0.0

    def negative_reml(tau2: float) -> float:
        total_variance = variances + tau2
        weights = 1.0 / total_variance
        mean = np.sum(weights * effects) / np.sum(weights)
        return 0.5 * (
            np.sum(np.log(total_variance))
            + np.log(np.sum(weights))
            + np.sum(weights * (effects - mean) ** 2)
        )

    upper = max(10.0, float(np.var(effects, ddof=1) * 100.0))
    optimization = minimize_scalar(
        negative_reml, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-10}
    )
    tau2 = max(0.0, float(optimization.x)) if optimization.success else 0.0
    weights = 1.0 / (variances + tau2)
    pooled = float(np.sum(weights * effects) / np.sum(weights))
    raw_hk_scale = float(np.sum(weights * (effects - pooled) ** 2) / (len(effects) - 1))
    modified_hk_scale = max(1.0, raw_hk_scale)
    se = float(np.sqrt(modified_hk_scale / np.sum(weights)))
    meta_df = float(len(effects) - 1)
    p_value = float(2.0 * t.sf(abs(pooled / se), df=meta_df))
    return {
        "effect": pooled,
        "SE": se,
        "p": p_value,
        "tau2": tau2,
        "I2": i2,
        "mKH_scale": modified_hk_scale,
        "df": meta_df,
    }


def informative_datasets(metadata: pd.DataFrame, min_group: int) -> list[str]:
    result = []
    for dataset, group in metadata.groupby("dataset", observed=True, sort=True):
        case = int(group["group"].eq("case").sum())
        control = int(group["group"].eq("control").sum())
        if case >= min_group and control >= min_group:
            result.append(str(dataset))
    return result


def meta_table(
    activity: np.ndarray,
    tfs: np.ndarray,
    metadata: pd.DataFrame,
    min_group: int,
    min_datasets: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    datasets = informative_datasets(metadata, min_group)
    per_dataset_rows: list[dict[str, object]] = []
    meta_rows: list[dict[str, object]] = []
    for tf_index, tf in enumerate(tfs):
        effects: list[float] = []
        variances: list[float] = []
        directions: list[float] = []
        for dataset in datasets:
            index = metadata["dataset"].eq(dataset).to_numpy()
            labels = metadata.loc[index, "group"].eq("case").to_numpy()
            effect, variance = hedges_g(
                activity[index, tf_index][labels], activity[index, tf_index][~labels]
            )
            if not np.isfinite(effect) or not np.isfinite(variance):
                continue
            effects.append(effect)
            variances.append(variance)
            directions.append(float(np.sign(effect)))
            per_dataset_rows.append(
                {
                    "TF": tf,
                    "dataset": dataset,
                    "Hedges_g": effect,
                    "variance": variance,
                    "case_n": int(labels.sum()),
                    "control_n": int((~labels).sum()),
                }
            )
        if len(effects) < 2:
            continue
        meta = reml_meta(np.asarray(effects), np.asarray(variances))
        meta_direction = np.sign(meta["effect"])
        direction_fraction = (
            float(np.mean(np.asarray(directions) == meta_direction)) if meta_direction != 0 else 0.0
        )
        meta_rows.append(
            {
                "TF": tf,
                "meta_effect": meta["effect"],
                "meta_SE": meta["SE"],
                "meta_p": meta["p"],
                "tau2_REML": meta["tau2"],
                "I2": meta["I2"],
                "modified_Knapp_Hartung_scale": meta["mKH_scale"],
                "meta_df": meta["df"],
                "informative_datasets": len(effects),
                "same_direction_fraction": direction_fraction,
            }
        )
    meta_frame = pd.DataFrame(meta_rows)
    if not meta_frame.empty:
        meta_frame["meta_FDR"] = bh_fdr(meta_frame["meta_p"].to_numpy())
        meta_frame["reproducible"] = (
            meta_frame["informative_datasets"].ge(min_datasets)
            & meta_frame["meta_FDR"].le(0.05)
            & meta_frame["meta_effect"].abs().ge(0.50)
            & meta_frame["same_direction_fraction"].ge(0.75)
            & meta_frame["I2"].le(75.0)
        )
        meta_frame = meta_frame.sort_values(
            ["reproducible", "meta_FDR", "meta_effect"], ascending=[False, True, False]
        ).reset_index(drop=True)
    return meta_frame, pd.DataFrame(per_dataset_rows)


def vectorized_hedges(activity: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    case = activity[labels]
    control = activity[~labels]
    n1, n0 = len(case), len(control)
    pooled_variance = (
        (n1 - 1) * case.var(axis=0, ddof=1) + (n0 - 1) * control.var(axis=0, ddof=1)
    ) / (n1 + n0 - 2)
    correction = 1.0 - 3.0 / (4.0 * (n1 + n0) - 9.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        effect = correction * (case.mean(axis=0) - control.mean(axis=0)) / np.sqrt(pooled_variance)
        variance = (n1 + n0) / (n1 * n0) + effect * effect / (2.0 * (n1 + n0 - 2))
    effect[~np.isfinite(effect)] = 0.0
    variance[~np.isfinite(variance) | (variance <= 0)] = np.inf
    return effect, variance


def select_training_features(
    activity: np.ndarray,
    labels: np.ndarray,
    datasets: np.ndarray,
    train: np.ndarray,
    top_n: int,
) -> np.ndarray:
    numerator = np.zeros(activity.shape[1], dtype=float)
    denominator = np.zeros(activity.shape[1], dtype=float)
    for dataset in sorted(set(datasets[train])):
        index = train & (datasets == dataset)
        local = labels[index]
        if local.sum() < 3 or (~local).sum() < 3:
            continue
        effect, variance = vectorized_hedges(activity[index], local)
        weight = np.where(np.isfinite(variance), 1.0 / variance, 0.0)
        numerator += weight * effect
        denominator += weight
    score = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
    order = np.lexsort((np.arange(len(score)), -np.abs(score)))
    return order[: min(top_n, len(order))]


def loso_predictions(
    activity: np.ndarray,
    labels: np.ndarray,
    datasets: np.ndarray,
    held_datasets: list[str],
    top_n: int = 20,
) -> tuple[np.ndarray, pd.DataFrame]:
    predictions = np.full(len(labels), np.nan, dtype=float)
    rows: list[dict[str, object]] = []
    for held in held_datasets:
        test = datasets == held
        train = np.isin(datasets, held_datasets) & ~test
        if labels[train].sum() < 3 or (~labels[train]).sum() < 3:
            continue
        selected = select_training_features(activity, labels, datasets, train, top_n)
        scaler = StandardScaler().fit(activity[train][:, selected])
        model = LogisticRegression(
            class_weight="balanced",
            penalty="l2",
            C=1.0,
            max_iter=5000,
            random_state=1729,
            solver="liblinear",
        )
        model.fit(scaler.transform(activity[train][:, selected]), labels[train])
        predictions[test] = model.predict_proba(scaler.transform(activity[test][:, selected]))[:, 1]
        local_auc = float(roc_auc_score(labels[test], predictions[test]))
        rows.append(
            {
                "held_dataset": held,
                "train_n": int(train.sum()),
                "test_n": int(test.sum()),
                "case_n": int(labels[test].sum()),
                "control_n": int((~labels[test]).sum()),
                "selected_features": int(len(selected)),
                "AUROC": local_auc,
            }
        )
    return predictions, pd.DataFrame(rows)


def stratified_auc(labels: np.ndarray, predictions: np.ndarray, datasets: np.ndarray) -> float:
    numerator = 0.0
    denominator = 0.0
    for dataset in sorted(set(datasets)):
        index = (datasets == dataset) & np.isfinite(predictions)
        local = labels[index]
        if local.sum() == 0 or (~local).sum() == 0:
            continue
        pairs = float(local.sum() * (~local).sum())
        numerator += pairs * roc_auc_score(local, predictions[index])
        denominator += pairs
    return float(numerator / denominator) if denominator > 0 else np.nan


def validate_loso(
    activity: np.ndarray,
    labels: np.ndarray,
    datasets: np.ndarray,
    held_datasets: list[str],
    permutations: int,
    bootstraps: int,
    seed: int = 1729,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    predictions, folds = loso_predictions(activity, labels, datasets, held_datasets)
    observed = stratified_auc(labels, predictions, datasets)
    rng = np.random.default_rng(seed)
    strata = [
        np.flatnonzero((datasets == dataset) & (labels == outcome))
        for dataset in held_datasets
        for outcome in (False, True)
    ]
    bootstrap_values: list[float] = []
    for _ in range(bootstraps):
        sampled = np.concatenate([rng.choice(index, size=len(index), replace=True) for index in strata if len(index)])
        value = stratified_auc(labels[sampled], predictions[sampled], datasets[sampled])
        if np.isfinite(value):
            bootstrap_values.append(value)

    permutation_rows: list[dict[str, float]] = []
    for permutation in range(permutations):
        permuted = labels.copy()
        for dataset in held_datasets:
            index = np.flatnonzero(datasets == dataset)
            permuted[index] = rng.permutation(permuted[index])
        null_predictions, _ = loso_predictions(activity, permuted, datasets, held_datasets)
        value = stratified_auc(permuted, null_predictions, datasets)
        permutation_rows.append({"permutation": permutation + 1, "stratified_AUROC": value})
    permutation_frame = pd.DataFrame(permutation_rows)
    empirical_p = float(
        (1 + permutation_frame["stratified_AUROC"].ge(observed).sum())
        / (1 + len(permutation_frame))
    )
    prediction_frame = pd.DataFrame(
        {"dataset": datasets, "case": labels.astype(int), "prediction": predictions}
    )
    bootstrap_frame = pd.DataFrame({"bootstrap": np.arange(1, len(bootstrap_values) + 1), "stratified_AUROC": bootstrap_values})
    receipt = {
        "stratified_AUROC": observed,
        "bootstrap_95_CI": [
            float(np.quantile(bootstrap_values, 0.025)),
            float(np.quantile(bootstrap_values, 0.975)),
        ],
        "bootstrap_replicates": len(bootstrap_values),
        "permutations": len(permutation_frame),
        "permutation_empirical_p": empirical_p,
    }
    return receipt, prediction_frame, folds, permutation_frame, bootstrap_frame


def run(
    activity_path: Path,
    metadata_path: Path,
    msviper_path: Path,
    manifest_path: Path,
    output_dir: Path,
    permutations: int,
    bootstraps: int,
) -> dict[str, object]:
    manifest = PhenotypeManifest.load(manifest_path)
    activity_frame = read_activity(activity_path)
    metadata = pd.read_csv(metadata_path, sep="\t")
    position = metadata["donor_id"].map({patient: index for index, patient in enumerate(activity_frame.index)})
    if position.isna().any() or position.duplicated().any():
        raise ValueError("Patient metadata does not map one-to-one to VIPER activity")
    activity_frame = activity_frame.loc[metadata["donor_id"]]
    tfs = activity_frame.columns.to_numpy(str)
    all_activity = activity_frame.to_numpy(float)
    if not np.isfinite(all_activity).all():
        raise ValueError("VIPER activity contains non-finite values")

    primary_threshold = int(manifest["primary_min_cells"])
    primary_mask = metadata["analysis_cells_aggregated"].ge(primary_threshold).to_numpy()
    primary_metadata = metadata.loc[primary_mask].reset_index(drop=True)
    primary_activity = all_activity[primary_mask]
    min_group = int(manifest.raw["minimum_group_per_dataset"])
    info = informative_datasets(primary_metadata, min_group)
    labels = primary_metadata["group"].eq("case").to_numpy()
    datasets = primary_metadata["dataset"].astype(str).to_numpy()
    max_case_fraction = float(
        primary_metadata.loc[labels, "dataset"].value_counts(normalize=True).max()
    )
    cohort_conditions = {
        "case_patients_ge_minimum": int(labels.sum()) >= int(manifest.raw["minimum_case_patients"]),
        "control_patients_ge_minimum": int((~labels).sum()) >= int(manifest.raw["minimum_control_patients"]),
        "informative_datasets_ge_minimum": len(info) >= int(manifest.raw["minimum_informative_datasets"]),
        "max_case_dataset_fraction_le_limit": max_case_fraction <= float(manifest.raw["maximum_case_fraction_from_one_dataset"]),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    minimum_informative = int(manifest.raw["minimum_informative_datasets"])
    meta, per_dataset = meta_table(
        primary_activity, tfs, primary_metadata, min_group, minimum_informative
    )
    meta_path = output_dir / "primary_random_effects_meta.tsv.gz"
    effects_path = output_dir / "primary_dataset_effects.tsv.gz"
    meta.to_csv(meta_path, sep="\t", index=False, compression="gzip", lineterminator="\n")
    per_dataset.to_csv(effects_path, sep="\t", index=False, compression="gzip", lineterminator="\n")
    reproducible = meta[meta["reproducible"]].copy()

    ms = pd.read_csv(msviper_path, sep="\t")
    merged = reproducible.merge(ms[["TF", "NES", "FDR"]], on="TF", how="left", suffixes=("", "_msviper"))
    merged["msviper_direction_agrees"] = np.sign(merged["meta_effect"]) == np.sign(merged["NES"])
    merged["msviper_confirmed"] = merged["FDR"].le(0.01) & merged["msviper_direction_agrees"]
    reproducible_path = output_dir / "reproducible_TF_program.tsv"
    merged.to_csv(reproducible_path, sep="\t", index=False, lineterminator="\n")

    loso, predictions, folds, permutation_frame, bootstrap_frame = validate_loso(
        primary_activity,
        labels,
        datasets,
        info,
        permutations=permutations,
        bootstraps=bootstraps,
    )
    predictions.insert(0, "donor_id", primary_metadata["donor_id"].to_numpy())
    predictions.to_csv(output_dir / "loso_predictions.tsv", sep="\t", index=False, lineterminator="\n")
    folds.to_csv(output_dir / "loso_folds.tsv", sep="\t", index=False, lineterminator="\n")
    permutation_frame.to_csv(output_dir / "loso_permutations.tsv.gz", sep="\t", index=False, compression="gzip", lineterminator="\n")
    bootstrap_frame.to_csv(output_dir / "loso_bootstrap.tsv.gz", sep="\t", index=False, compression="gzip", lineterminator="\n")

    sensitivity_rows: list[dict[str, object]] = []
    for threshold in sorted(set(int(x) for x in manifest.raw.get("sensitivity_min_cells", []))):
        mask = metadata["analysis_cells_aggregated"].ge(threshold).to_numpy()
        local_metadata = metadata.loc[mask].reset_index(drop=True)
        local_activity = all_activity[mask]
        local_info = informative_datasets(local_metadata, min_group)
        local_meta, _ = meta_table(
            local_activity, tfs, local_metadata, min_group, minimum_informative
        )
        local_path = output_dir / f"sensitivity_meta_cells_{threshold}.tsv.gz"
        local_meta.to_csv(local_path, sep="\t", index=False, compression="gzip", lineterminator="\n")
        if len(local_info) >= 2:
            local_labels = local_metadata["group"].eq("case").to_numpy()
            local_datasets = local_metadata["dataset"].astype(str).to_numpy()
            local_predictions, _ = loso_predictions(local_activity, local_labels, local_datasets, local_info)
            local_auc = stratified_auc(local_labels, local_predictions, local_datasets)
        else:
            local_auc = np.nan
        local_repro = set(local_meta.loc[local_meta.get("reproducible", False), "TF"]) if not local_meta.empty else set()
        primary_repro = set(reproducible["TF"])
        sensitivity_rows.append(
            {
                "min_cells": threshold,
                "patients": int(mask.sum()),
                "case": int(local_metadata["group"].eq("case").sum()),
                "control": int(local_metadata["group"].eq("control").sum()),
                "informative_datasets": len(local_info),
                "reproducible_TFs": len(local_repro),
                "primary_program_overlap": len(local_repro & primary_repro),
                "observed_stratified_LOSO_AUROC": local_auc,
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(output_dir / "sensitivity_summary.tsv", sep="\t", index=False, lineterminator="\n")

    scientific_conditions = {
        "reproducible_TFs_ge_10": len(reproducible) >= 10,
        "msviper_confirmed_TFs_ge_5": int(merged["msviper_confirmed"].sum()) >= 5,
    }
    validation_conditions = {
        "LOSO_AUROC_ge_0_65": loso["stratified_AUROC"] >= 0.65,
        "LOSO_bootstrap_lower_gt_0_55": loso["bootstrap_95_CI"][0] > 0.55,
        "LOSO_permutation_p_le_0_05": loso["permutation_empirical_p"] <= 0.05,
    }
    output_paths = [path for path in output_dir.iterdir() if path.is_file() and path.name != "statistics_receipt.json"]
    receipt = {
        "status": "completed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "primary_threshold": primary_threshold,
        "patients": int(len(primary_metadata)),
        "case_patients": int(labels.sum()),
        "control_patients": int((~labels).sum()),
        "informative_datasets": info,
        "max_case_fraction_from_one_dataset": max_case_fraction,
        "activity_TFs": int(len(tfs)),
        "meta_TFs_tested": int(len(meta)),
        "reproducible_TFs": int(len(reproducible)),
        "msviper_confirmed_reproducible_TFs": int(merged["msviper_confirmed"].sum()),
        "cohort_conditions": cohort_conditions,
        "scientific_conditions": scientific_conditions,
        "secondary_validation_conditions": validation_conditions,
        "loso": loso,
        "outputs": output_manifest(output_paths),
    }
    write_json(output_dir / "statistics_receipt.json", receipt)
    print(json.dumps(receipt, ensure_ascii=False))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--activity", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--msviper", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--permutations", type=int, default=500)
    parser.add_argument("--bootstraps", type=int, default=1000)
    args = parser.parse_args()
    run(
        args.activity,
        args.metadata,
        args.msviper,
        args.manifest,
        args.output_dir,
        args.permutations,
        args.bootstraps,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
