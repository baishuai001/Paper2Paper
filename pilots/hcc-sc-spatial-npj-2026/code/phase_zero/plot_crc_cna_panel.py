#!/usr/bin/env python3
"""Plot auditable sample-level CNA evidence without implying DNA truth."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REQUIRED_COLUMNS = {
    "panel_order",
    "sample_type",
    "author_cancer_cells",
    "author_cancer_two_method_support",
    "author_cancer_one_method_support",
    "author_cancer_unresolved_without_support",
    "author_cancer_defined_without_support_not_nonmalignancy_proof",
    "normal_reference_cells",
    "normal_reference_copykat_aneuploid",
    "normal_reference_scevan_tumor",
    "normal_reference_unexpected_malignancy_support",
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def build_plot_source(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"CNA unit summary missing or empty: {path}")
    table = pd.read_csv(path, sep="\t")
    missing = REQUIRED_COLUMNS - set(table.columns)
    if missing:
        raise ValueError(f"CNA plot fields missing: {sorted(missing)}")
    if table.empty or table["panel_order"].duplicated().any():
        raise ValueError("CNA plot requires nonempty unique panel orders")
    count_columns = sorted(REQUIRED_COLUMNS - {"panel_order", "sample_type"})
    for column in count_columns:
        table[column] = pd.to_numeric(table[column], errors="raise")
        if (table[column] < 0).any():
            raise ValueError(f"CNA plot count is negative: {column}")
    exclusive = [
        "author_cancer_two_method_support",
        "author_cancer_one_method_support",
        "author_cancer_unresolved_without_support",
        "author_cancer_defined_without_support_not_nonmalignancy_proof",
    ]
    if not table[exclusive].sum(axis=1).eq(table["author_cancer_cells"]).all():
        raise ValueError("exclusive author-Cancer CNA classes do not sum to the input total")
    if (table["normal_reference_cells"] <= 0).any():
        raise ValueError("CNA plot contains a unit without normal references")
    table = table.sort_values("panel_order").reset_index(drop=True)
    for column in exclusive:
        table[column + "_fraction"] = table[column] / table["author_cancer_cells"]
    for column in (
        "normal_reference_copykat_aneuploid",
        "normal_reference_scevan_tumor",
        "normal_reference_unexpected_malignancy_support",
    ):
        table[column + "_fraction"] = table[column] / table["normal_reference_cells"]
    return table


def plot_panel(unit_summary: Path, output_dir: Path) -> dict[str, object]:
    source = build_plot_source(unit_summary)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to mix with existing CNA plot output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    source_path = output_dir / "P0_cna_panel_plot_source.tsv"
    figure_path = output_dir / "P0_cna_panel_diagnostic.png"
    source.to_csv(source_path, sep="\t", index=False, lineterminator="\n")

    x = np.arange(len(source))
    figure, axes = plt.subplots(
        2,
        1,
        figsize=(12, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]},
        constrained_layout=True,
    )
    stack = [
        ("author_cancer_two_method_support_fraction", "Two-method support", "#2468a2"),
        ("author_cancer_one_method_support_fraction", "One-method support", "#f28e2b"),
        ("author_cancer_unresolved_without_support_fraction", "Unresolved, no support", "#bab0ac"),
        (
            "author_cancer_defined_without_support_not_nonmalignancy_proof_fraction",
            "Defined without CNA support (not nonmalignancy proof)",
            "#e6e6e6",
        ),
    ]
    bottom = np.zeros(len(source))
    for column, label, color in stack:
        values = source[column].to_numpy(dtype=float)
        axes[0].bar(x, values, bottom=bottom, label=label, color=color, width=0.82)
        bottom += values
    axes[0].set_ylim(0, 1.02)
    axes[0].set_ylabel("Fraction of author Cancer candidates")
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2, frameon=False)

    unexpected = source["normal_reference_unexpected_malignancy_support_fraction"].to_numpy(dtype=float)
    axes[1].bar(x, unexpected, color="#d62728", alpha=0.45, label="Either method")
    axes[1].scatter(
        x,
        source["normal_reference_copykat_aneuploid_fraction"],
        color="#7f0000",
        marker="o",
        label="CopyKAT aneuploid",
        zorder=3,
    )
    axes[1].scatter(
        x,
        source["normal_reference_scevan_tumor_fraction"],
        color="#9467bd",
        marker="x",
        label="SCEVAN tumor",
        zorder=3,
    )
    axes[1].axhline(0, color="black", linewidth=0.7)
    axes[1].set_ylabel("Unexpected support\nin normal references")
    axes[1].set_xlabel("Frozen CNA panel order")
    axes[1].legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.27))
    axes[1].set_xticks(x, source["panel_order"].astype(str))

    previous_state = None
    start = 0
    for index, state in enumerate(source["sample_type"].astype(str).tolist() + [None]):
        if previous_state is None:
            previous_state = state
        elif state != previous_state:
            midpoint = (start + index - 1) / 2
            axes[0].text(midpoint, 1.035, previous_state, ha="center", va="bottom", fontsize=10)
            if index < len(source):
                for axis in axes:
                    axis.axvline(index - 0.5, color="black", linewidth=0.7, linestyle=":")
            start = index
            previous_state = state
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    receipt: dict[str, object] = {
        "status": "diagnostic_noninteractive_figure_written",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": sha256_file(Path(__file__)),
        "unit_summary_sha256": sha256_file(unit_summary),
        "samples": int(len(source)),
        "scientific_gate": "not_assessed_by_plot",
        "claim_boundary": (
            "The figure visualizes inferred CNA support and calibration warnings. It does not establish "
            "DNA truth, disease evolution, or statistical differences between cross-sectional states."
        ),
        "outputs": {},
    }
    for path in (source_path, figure_path):
        receipt["outputs"][path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    receipt_path = output_dir / "P0_cna_panel_plot_receipt.json"
    receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unit-summary", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    receipt = plot_panel(args.unit_summary, args.output_dir)
    print(json.dumps(receipt, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
