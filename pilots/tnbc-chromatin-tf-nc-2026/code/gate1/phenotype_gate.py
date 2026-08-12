#!/usr/bin/env python3
"""Independently reconstruct and validate the four CRC immune phenotypes."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from crc_cohort import eligible_patient_table, load_obs
from gate1_common import sha256_file, write_json


IMMUNE_COARSE = ["B cell", "Mast cell", "Myeloid cell", "NK", "Plasma cell", "T cell"]


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Positive means values in x tend to exceed values in y."""
    return float((np.greater.outer(x, y).sum() - np.less.outer(x, y).sum()) / (len(x) * len(y)))


def residualize_by_dataset(matrix: np.ndarray, dataset: pd.Series) -> np.ndarray:
    design = pd.get_dummies(dataset, drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(dataset)), design.to_numpy()])
    return matrix - X @ np.linalg.lstsq(X, matrix, rcond=None)[0]


def lodo_auc(features: np.ndarray, labels: np.ndarray, datasets: np.ndarray) -> tuple[float, pd.DataFrame]:
    prediction = np.full(len(labels), np.nan)
    rows: list[dict[str, object]] = []
    for held in sorted(set(datasets)):
        test = datasets == held
        train = ~test
        if labels[train].min() == labels[train].max():
            continue
        scale = StandardScaler().fit(features[train])
        model = LogisticRegression(class_weight="balanced", penalty="l2", C=1.0, max_iter=5000, random_state=1729)
        model.fit(scale.transform(features[train]), labels[train])
        prediction[test] = model.predict_proba(scale.transform(features[test]))[:, 1]
        fold_auc = roc_auc_score(labels[test], prediction[test]) if len(np.unique(labels[test])) == 2 else np.nan
        rows.append({"held_out_dataset": held, "patients": int(test.sum()), "M": int(labels[test].sum()), "auc": fold_auc})
    valid = np.isfinite(prediction)
    return float(roc_auc_score(labels[valid], prediction[valid])), pd.DataFrame(rows)


def run(h5ad: Path, output_dir: Path) -> dict[str, object]:
    obs = load_obs(h5ad)
    patient, eligible = eligible_patient_table(obs)
    patient = patient[patient["immune_label_n"].eq(1)].copy().sort_values("donor_id").reset_index(drop=True)
    scoped = obs.loc[eligible & obs["donor_id"].isin(patient["donor_id"])].copy()

    counts = pd.crosstab(scoped["donor_id"], scoped["cell_type_coarse_crc_atlas"]).reindex(
        index=patient["donor_id"], columns=IMMUNE_COARSE, fill_value=0
    )
    # The paper removes the neutrophil fraction as a feature.  Fractions remain
    # fractions of all eligible primary-tumour cells; otherwise the immune-
    # desert class is mathematically erased by closure over immune cells only.
    denominator = scoped.groupby("donor_id", observed=True).size().reindex(counts.index)
    fractions = counts.div(denominator.replace(0, np.nan), axis=0).fillna(0.0)
    patient = patient.set_index("donor_id").join(fractions.add_prefix("fraction__")).reset_index()
    feature_columns = [f"fraction__{x}" for x in IMMUNE_COARSE]
    feature = patient[feature_columns].to_numpy(float)
    corrected = residualize_by_dataset(feature, patient["dataset"])

    distance = pdist(corrected, metric="correlation")
    if not np.isfinite(distance).all():
        raise ValueError("Non-finite correlation distances after batch residualization")
    hierarchical = fcluster(linkage(distance, method="average"), 4, criterion="maxclust")
    truth4 = pd.Categorical(patient["immune_label"], categories=["desert", "B", "T", "M"]).codes
    hierarchical_ari = float(adjusted_rand_score(truth4, hierarchical))

    scaled = StandardScaler().fit_transform(corrected)
    kmeans_ari = np.asarray(
        [adjusted_rand_score(truth4, KMeans(n_clusters=4, n_init=20, random_state=seed).fit_predict(scaled)) for seed in range(50)]
    )

    m_fraction = patient["fraction__Myeloid cell"].to_numpy(float)
    is_m = patient["immune_label"].eq("M").to_numpy()
    mw = mannwhitneyu(m_fraction[is_m], m_fraction[~is_m], alternative="two-sided")
    delta = cliffs_delta(m_fraction[is_m], m_fraction[~is_m])

    pseudocount = 0.5
    clr_counts = counts.to_numpy(float) + pseudocount
    clr = np.log(clr_counts) - np.log(clr_counts).mean(axis=1, keepdims=True)
    auc, folds = lodo_auc(clr, is_m.astype(int), patient["dataset"].to_numpy(str))

    conditions = {
        "M_myeloid_fraction_p_lt_0_01": bool(mw.pvalue < 0.01),
        "M_myeloid_fraction_delta_ge_0_33": bool(delta >= 0.33),
        "hierarchical_ARI_ge_0_50": bool(hierarchical_ari >= 0.50),
        "median_kmeans_ARI_ge_0_40": bool(np.median(kmeans_ari) >= 0.40),
        "LODO_M_AUROC_ge_0_70": bool(auc >= 0.70),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    patient["hierarchical_cluster"] = hierarchical
    patient.to_csv(output_dir / "G2_patient_immune_fractions.tsv", sep="\t", index=False, lineterminator="\n")
    folds.to_csv(output_dir / "G2_phenotype_lodo_folds.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"seed": np.arange(50), "ARI": kmeans_ari}).to_csv(
        output_dir / "G2_kmeans_seed_ARI.tsv", sep="\t", index=False, lineterminator="\n"
    )
    receipt = {
        "status": "passed" if all(conditions.values()) else "failed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "patients": int(len(patient)),
        "M_patients": int(is_m.sum()),
        "nonM_patients": int((~is_m).sum()),
        "datasets": int(patient["dataset"].nunique()),
        "M_myeloid_fraction_median": float(np.median(m_fraction[is_m])),
        "nonM_myeloid_fraction_median": float(np.median(m_fraction[~is_m])),
        "mann_whitney_p_two_sided": float(mw.pvalue),
        "cliffs_delta_M_minus_nonM": delta,
        "hierarchical_ARI": hierarchical_ari,
        "kmeans_ARI_median": float(np.median(kmeans_ari)),
        "kmeans_ARI_range": [float(kmeans_ari.min()), float(kmeans_ari.max())],
        "LODO_M_AUROC": auc,
        "conditions": conditions,
        "implementation_note": "Independent reconstruction from frozen H5AD annotations; no author cluster labels were used during clustering.",
        "outputs": {},
    }
    for path in sorted(output_dir.glob("G2_*.tsv")):
        receipt["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    write_json(output_dir / "G2_phenotype_gate_receipt.json", receipt)
    print(receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    run(args.h5ad, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
