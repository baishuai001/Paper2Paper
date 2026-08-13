#!/usr/bin/env python3
"""Independent, post-decision arithmetic audit of a completed regulatory gate.

This script does not alter the primary analysis or verdict. It recomputes the
key multiplicity, cohort, reproducibility and LOSO summaries directly from the
saved result tables so that the final report does not rely on one code path.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def bh_fdr(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    p = np.clip(values[finite], 0.0, 1.0)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    result[finite] = restored
    return result


def concordance_auc(case: np.ndarray, control: np.ndarray) -> float:
    comparison = case[:, None] - control[None, :]
    return float((np.sum(comparison > 0) + 0.5 * np.sum(comparison == 0)) / comparison.size)


def stratified_auc(frame: pd.DataFrame) -> float:
    numerator = 0.0
    denominator = 0.0
    for _, group in frame[np.isfinite(frame["prediction"])].groupby("dataset", sort=True):
        case = group.loc[group["case"].eq(1), "prediction"].to_numpy(float)
        control = group.loc[group["case"].eq(0), "prediction"].to_numpy(float)
        if len(case) == 0 or len(control) == 0:
            continue
        pairs = float(len(case) * len(control))
        numerator += pairs * concordance_auc(case, control)
        denominator += pairs
    return float(numerator / denominator)


def close(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return bool(np.isfinite(left) and np.isfinite(right) and abs(left - right) <= tolerance)


def records(frame: pd.DataFrame, columns: list[str], n: int = 10) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in frame.loc[:, columns].head(n).to_dict(orient="records"):
        rows.append(
            {
                key: value.item() if isinstance(value, np.generic) else value
                for key, value in row.items()
            }
        )
    return rows


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="\n", dir=path.parent, delete=False
    ) as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
        stream.write("\n")
        temporary = Path(stream.name)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--statistics-dir", required=True, type=Path)
    parser.add_argument("--msviper", required=True, type=Path)
    parser.add_argument("--viper-receipt", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--decision", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    with args.decision.open(encoding="utf-8") as stream:
        decision = json.load(stream)
    with args.viper_receipt.open(encoding="utf-8") as stream:
        viper_receipt = json.load(stream)
    with (args.statistics_dir / "statistics_receipt.json").open(encoding="utf-8") as stream:
        receipt = json.load(stream)

    meta = pd.read_csv(args.statistics_dir / "primary_random_effects_meta.tsv.gz", sep="\t")
    effects = pd.read_csv(args.statistics_dir / "primary_dataset_effects.tsv.gz", sep="\t")
    program = pd.read_csv(args.statistics_dir / "reproducible_TF_program.tsv", sep="\t")
    sensitivity = pd.read_csv(args.statistics_dir / "sensitivity_summary.tsv", sep="\t")
    predictions = pd.read_csv(args.statistics_dir / "loso_predictions.tsv", sep="\t")
    folds = pd.read_csv(args.statistics_dir / "loso_folds.tsv", sep="\t")
    null = pd.read_csv(args.statistics_dir / "loso_permutations.tsv.gz", sep="\t")
    bootstrap = pd.read_csv(args.statistics_dir / "loso_bootstrap.tsv.gz", sep="\t")
    ms = pd.read_csv(args.msviper, sep="\t")
    metadata = pd.read_csv(args.metadata, sep="\t")

    meta_bh = bh_fdr(meta["meta_p"].to_numpy(float))
    ms_bh = bh_fdr(ms["p.value"].to_numpy(float))
    max_meta_bh_difference = float(np.max(np.abs(meta_bh - meta["meta_FDR"].to_numpy(float))))
    max_ms_bh_difference = float(np.max(np.abs(ms_bh - ms["FDR"].to_numpy(float))))

    recomputed_reproducible = (
        meta["informative_datasets"].ge(int(manifest["minimum_informative_datasets"]))
        & meta["meta_FDR"].le(0.05)
        & meta["meta_effect"].abs().ge(0.50)
        & meta["same_direction_fraction"].ge(0.75)
        & meta["I2"].le(75.0)
    )
    stored_reproducible = meta["reproducible"].astype(str).str.casefold().eq("true")

    threshold = int(manifest["primary_min_cells"])
    primary = metadata[metadata["analysis_cells_aggregated"].ge(threshold)].copy()
    study_counts = (
        primary.groupby(["dataset", "group"], observed=True).size().unstack(fill_value=0)
    )
    for group in ("case", "control"):
        if group not in study_counts:
            study_counts[group] = 0
    minimum_group = int(manifest["minimum_group_per_dataset"])
    informative = sorted(
        study_counts.index[
            study_counts["case"].ge(minimum_group)
            & study_counts["control"].ge(minimum_group)
        ].astype(str)
    )
    case_counts = primary.loc[primary["group"].eq("case"), "dataset"].value_counts()
    max_case_fraction = float(case_counts.max() / case_counts.sum())

    observed_auc = stratified_auc(predictions)
    bootstrap_values = bootstrap["stratified_AUROC"].to_numpy(float)
    bootstrap_ci = np.quantile(bootstrap_values, [0.025, 0.975])
    null_values = null["stratified_AUROC"].to_numpy(float)
    permutation_p = float((1 + np.sum(null_values >= observed_auc)) / (1 + len(null_values)))
    fold_auc = {
        str(study): stratified_auc(predictions[predictions["dataset"].eq(study)])
        for study in informative
    }

    effect_study_counts = effects.groupby("TF")["dataset"].nunique()
    checks = {
        "decision_is_FAIL_and_stop": decision["verdict"] == "FAIL" and decision["stop_now"] is True,
        "all_data_conditions_pass": all(decision["data_conditions"].values()),
        "scientific_conditions_fail": not all(decision["scientific_conditions"].values()),
        "meta_has_2028_unique_TFs": len(meta) == 2028 and meta["TF"].nunique() == 2028,
        "meta_BH_exact": max_meta_bh_difference <= 1e-12,
        "msVIPER_BH_exact": max_ms_bh_difference <= 1e-12,
        "stored_reproducible_flags_exact": bool(np.array_equal(recomputed_reproducible, stored_reproducible)),
        "program_table_count_matches": len(program) == int(recomputed_reproducible.sum()),
        "all_meta_TFs_have_three_studies": bool((effect_study_counts == 3).all()),
        "msVIPER_FDR_0_01_count_matches": int(ms["FDR"].le(0.01).sum()) == int(viper_receipt["msviper"]["TFs_FDR_0_01"]),
        "primary_patient_counts_match": (
            len(primary) == int(receipt["patients"])
            and int(primary["group"].eq("case").sum()) == int(receipt["case_patients"])
            and int(primary["group"].eq("control").sum()) == int(receipt["control_patients"])
        ),
        "informative_studies_match": informative == sorted(receipt["informative_datasets"]),
        "max_case_fraction_matches": close(max_case_fraction, float(receipt["max_case_fraction_from_one_dataset"])),
        "LOSO_AUROC_matches": close(observed_auc, float(receipt["loso"]["stratified_AUROC"])),
        "LOSO_bootstrap_CI_matches": bool(np.allclose(bootstrap_ci, receipt["loso"]["bootstrap_95_CI"], atol=1e-12, rtol=0)),
        "LOSO_permutation_p_matches": close(permutation_p, float(receipt["loso"]["permutation_empirical_p"])),
        "LOSO_fold_AUCs_match": all(
            close(fold_auc[str(row.held_dataset)], float(row.AUROC))
            for row in folds.itertuples(index=False)
        ),
        "all_sensitivities_have_zero_reproducible_TFs": bool(sensitivity["reproducible_TFs"].eq(0).all()),
    }

    criterion_counts = {
        "meta_nominal_p_le_0_05": int(meta["meta_p"].le(0.05).sum()),
        "meta_FDR_le_0_05": int(meta["meta_FDR"].le(0.05).sum()),
        "absolute_meta_effect_ge_0_50": int(meta["meta_effect"].abs().ge(0.50).sum()),
        "three_of_three_same_direction": int(meta["same_direction_fraction"].eq(1.0).sum()),
        "I2_le_75": int(meta["I2"].le(75.0).sum()),
        "all_nonmultiplicity_criteria": int(
            (
                meta["informative_datasets"].ge(3)
                & meta["meta_effect"].abs().ge(0.50)
                & meta["same_direction_fraction"].ge(0.75)
                & meta["I2"].le(75.0)
            ).sum()
        ),
        "all_reproducibility_criteria": int(recomputed_reproducible.sum()),
        "msVIPER_nominal_p_le_0_05": int(ms["p.value"].le(0.05).sum()),
        "msVIPER_FDR_le_0_01": int(ms["FDR"].le(0.01).sum()),
    }

    study_rows = []
    for study in informative:
        row = study_counts.loc[study]
        study_rows.append(
            {
                "study": study,
                "case": int(row["case"]),
                "control": int(row["control"]),
                "LOSO_AUROC": fold_auc[study],
            }
        )

    payload = {
        "status": "passed" if all(checks.values()) else "failed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "role": "post-decision independent arithmetic audit; does not alter the frozen verdict",
        "checks": checks,
        "criterion_counts": criterion_counts,
        "maximum_absolute_BH_difference": {
            "meta": max_meta_bh_difference,
            "msVIPER": max_ms_bh_difference,
        },
        "primary_cohort": {
            "patients": int(len(primary)),
            "case": int(primary["group"].eq("case").sum()),
            "control": int(primary["group"].eq("control").sum()),
            "max_case_fraction_from_one_study": max_case_fraction,
            "informative_studies": study_rows,
        },
        "LOSO": {
            "stratified_AUROC": observed_auc,
            "bootstrap_95_CI": bootstrap_ci.tolist(),
            "permutation_empirical_p": permutation_p,
        },
        "sensitivity": json.loads(sensitivity.to_json(orient="records")),
        "top_meta_by_p": records(
            meta.sort_values("meta_p"),
            ["TF", "meta_effect", "meta_p", "meta_FDR", "I2", "same_direction_fraction"],
        ),
        "top_msVIPER_by_p": records(
            ms.sort_values("p.value"), ["TF", "NES", "size", "p.value", "FDR"]
        ),
    }
    write_json(args.output, payload)
    print(json.dumps(payload, ensure_ascii=False))
    return 0 if payload["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
