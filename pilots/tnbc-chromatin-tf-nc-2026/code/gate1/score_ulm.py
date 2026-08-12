#!/usr/bin/env python3
"""Score signed TF activities with the univariate-linear-model t statistic."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from gate1_common import sha256_file, write_json


def ulm(expression: np.ndarray, genes: np.ndarray, edges: pd.DataFrame, min_targets: int = 10) -> tuple[np.ndarray, list[str], pd.DataFrame]:
    lookup = {str(g): i for i, g in enumerate(genes)}
    activities: list[np.ndarray] = []
    tf_names: list[str] = []
    rows: list[dict[str, object]] = []
    for tf, group in edges.groupby("source", sort=True):
        group = group[group["target"].isin(lookup)].drop_duplicates("target")
        if len(group) < min_targets:
            continue
        idx = np.asarray([lookup[x] for x in group["target"]], dtype=int)
        weights = group["sign"].to_numpy(float)
        centered_w = weights - weights.mean()
        denom = float(centered_w @ centered_w)
        if denom <= 0:
            continue
        y = expression[:, idx]
        centered_y = y - y.mean(axis=1, keepdims=True)
        slope = centered_y @ centered_w / denom
        residual = centered_y - slope[:, None] * centered_w[None, :]
        dof = len(idx) - 2
        se = np.sqrt(np.maximum((residual * residual).sum(axis=1) / dof / denom, np.finfo(float).tiny))
        activities.append(slope / se)
        tf_names.append(str(tf))
        rows.append({"TF": tf, "targets": int(len(idx)), "positive_targets": int((weights > 0).sum()), "negative_targets": int((weights < 0).sum())})
    return np.column_stack(activities), tf_names, pd.DataFrame(rows)


def run(logcpm: Path, network: Path, output_dir: Path, prefix: str) -> dict[str, object]:
    loaded = np.load(logcpm, allow_pickle=False)
    expression = loaded["log2cpm"].astype(float)
    genes = loaded["genes"].astype(str)
    edges = pd.read_csv(network, sep="\t")
    activity, tfs, coverage = ulm(expression, genes, edges)
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = output_dir / f"{prefix}_activities.npz"
    coverage_path = output_dir / f"{prefix}_coverage.tsv"
    np.savez_compressed(matrix_path, activity=activity, TF=np.asarray(tfs, dtype="U"))
    coverage.to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")
    receipt = {
        "status": "passed", "generated_at": datetime.now(timezone.utc).isoformat(), "patients": int(activity.shape[0]),
        "TFs": int(activity.shape[1]), "min_targets": 10, "score": "ULM slope t statistic with intercept",
        "outputs": {p.name: {"bytes": p.stat().st_size, "sha256": sha256_file(p)} for p in (matrix_path, coverage_path)},
    }
    write_json(output_dir / f"{prefix}_receipt.json", receipt); print(receipt); return receipt


def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument("--logcpm",required=True,type=Path); p.add_argument("--network",required=True,type=Path); p.add_argument("--output-dir",required=True,type=Path); p.add_argument("--prefix",required=True); a=p.parse_args(); run(a.logcpm,a.network,a.output_dir,a.prefix); return 0


if __name__ == "__main__": raise SystemExit(main())
