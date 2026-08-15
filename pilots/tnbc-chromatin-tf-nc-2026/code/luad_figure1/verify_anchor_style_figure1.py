#!/usr/bin/env python3

import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


if len(sys.argv) != 2:
    raise SystemExit("Usage: verify_anchor_style_figure1.py CLOUD_RUN_ROOT")

root = Path(sys.argv[1]).resolve()
figures = root / "results" / "figures"
tables = root / "results" / "tables"
processed = root / "data" / "processed"
audit = root / "audit" / "figure1_anchor_style"
audit.mkdir(parents=True, exist_ok=True)

stems = [
    "Figure1A_workflow",
    "Figure1C_TCGA_TF_activity_annotated",
    "SupplementaryFigure1A_TCGA_inclusion",
    "SupplementaryFigure1B_GSE81089_inclusion",
]
checks = []
for stem in stems:
    for suffix in (".pdf", ".png"):
        path = figures / f"{stem}{suffix}"
        require(path.is_file(), f"Missing output: {path}")
        require(path.stat().st_size > 10_000, f"Output is unexpectedly small: {path}")
        if suffix == ".pdf":
            require(path.read_bytes()[:4] == b"%PDF", f"Invalid PDF header: {path}")
        checks.append({"check": f"output:{stem}{suffix}", "status": "passed", "bytes": path.stat().st_size})

with gzip.open(tables / "Figure1C_TCGA_activity_matrix.tsv.gz", "rt") as handle:
    matrix = pd.read_csv(handle, sep="\t", index_col=0)
display = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1, ddof=1), axis=0)
require(display.shape == (351, 1017), f"Unexpected Figure 1C dimensions: {display.shape}")
require(np.isfinite(display.to_numpy()).all(), "Figure 1C display matrix contains non-finite values")
require(float(display.mean(axis=1).abs().max()) < 1e-8, "Figure 1C rows are not centered")
require(float((display.std(axis=1, ddof=1) - 1).abs().max()) < 1e-8, "Figure 1C rows are not unit scaled")
checks.append({"check": "Figure1C:row_scaled_351x1017", "status": "passed"})

selected = pd.read_csv(tables / "tcga_specific_TFs.tsv", sep="\t")
counts = selected.groupby("discovery_category").size().to_dict()
require(counts == {"LUAD": 158, "LUSC": 193}, f"Unexpected TF categories: {counts}")
checks.append({"check": "Figure1C:158_LUAD_193_LUSC_TFs", "status": "passed"})

tcga = pd.read_csv(processed / "tcga_manifest.tsv", sep="\t")
gse = pd.read_csv(processed / "gse81089_manifest.tsv", sep="\t")
require(tcga["group"].value_counts().to_dict() == {"LUAD": 516, "LUSC": 501}, "TCGA flow counts changed")
require(gse["group"].value_counts().to_dict() == {"LUAD": 106, "LUSC": 67}, "GSE81089 flow counts changed")
checks.append({"check": "SupplementaryFigure1:manifest_counts", "status": "passed"})

receipt = pd.read_csv(tables / "Figure1_anchor_style_render_receipt.tsv", sep="\t")
required_flow = {
    "tcga_all": 1141,
    "tcga_primary": 1029,
    "tcga_non_primary": 112,
    "tcga_unique": 1017,
    "tcga_duplicate": 12,
    "gse_all": 218,
    "gse_tumor": 199,
    "gse_normal": 19,
    "gse_target_histology": 175,
    "gse_non_target_histology": 24,
    "gse_target_with_expression": 173,
    "gse_unmatched": 2,
}
observed = dict(zip(receipt["metric"], receipt["value"]))
for key, expected in required_flow.items():
    require(int(round(observed[key])) == expected, f"Flow count mismatch for {key}: {observed[key]}")
checks.append({"check": "SupplementaryFigure1:flow_arithmetic", "status": "passed"})

verification = {
    "status": "passed",
    "outputs": stems,
    "figure1c_shape": list(display.shape),
    "tf_categories": counts,
    "checks": checks,
}
(audit / "verification.json").write_text(json.dumps(verification, indent=2), encoding="utf-8")
print(json.dumps(verification, indent=2))
