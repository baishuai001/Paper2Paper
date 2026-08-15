#!/usr/bin/env python3

import csv
import gzip
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def pdf_nonwhite_fraction(path: Path) -> float | None:
    pdftoppm = shutil.which("pdftoppm")
    if pdftoppm is None:
        return None
    with tempfile.TemporaryDirectory() as directory:
        prefix = Path(directory) / "page"
        subprocess.run(
            [pdftoppm, "-f", "1", "-singlefile", "-r", "36", str(path), str(prefix)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        payload = (prefix.with_suffix(".ppm")).read_bytes()
    match = re.match(br"P6\s+(?:#[^\n]*\s+)*(\d+)\s+(\d+)\s+(\d+)\s", payload)
    require(match is not None, f"Could not parse PDF preview for {path}")
    pixels = payload[match.end():]
    require(int(match.group(3)) == 255, f"Unsupported PDF preview depth for {path}")
    total = len(pixels) // 3
    nonwhite = sum(1 for index in range(0, len(pixels) - 2, 3) if min(pixels[index:index + 3]) < 245)
    return nonwhite / total


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
    "Figure1A_workflow_stacked",
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
            nonwhite_fraction = pdf_nonwhite_fraction(path)
            if nonwhite_fraction is not None:
                require(nonwhite_fraction > 0.005, f"PDF rendered as an empty or nearly empty page: {path}")
        else:
            nonwhite_fraction = None
        checks.append({
            "check": f"output:{stem}{suffix}", "status": "passed",
            "bytes": path.stat().st_size, "pdf_nonwhite_fraction": nonwhite_fraction,
        })

matrix_rows = 0
matrix_columns = None
max_abs_row_mean = 0.0
max_abs_row_sd_minus_1 = 0.0
with gzip.open(tables / "Figure1C_TCGA_activity_matrix.tsv.gz", "rt", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)
    matrix_columns = len(header) - 1
    for row in reader:
        values = [float(value) for value in row[1:]]
        require(len(values) == matrix_columns, f"Ragged Figure 1C matrix row: {row[0]}")
        require(all(math.isfinite(value) for value in values), f"Non-finite Figure 1C row: {row[0]}")
        mean = sum(values) / len(values)
        variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
        sd = math.sqrt(variance)
        max_abs_row_mean = max(max_abs_row_mean, abs(mean))
        max_abs_row_sd_minus_1 = max(max_abs_row_sd_minus_1, abs(sd - 1.0))
        matrix_rows += 1
matrix_shape = (matrix_rows, matrix_columns)
require(matrix_shape == (351, 1017), f"Unexpected Figure 1C dimensions: {matrix_shape}")
require(max_abs_row_mean < 1e-8, "Figure 1C rows are not centered")
require(max_abs_row_sd_minus_1 < 1e-8, "Figure 1C rows are not unit scaled")
checks.append({"check": "Figure1C:row_scaled_351x1017", "status": "passed"})

with (tables / "tcga_specific_TFs.tsv").open(newline="", encoding="utf-8") as handle:
    counts = dict(Counter(row["discovery_category"] for row in csv.DictReader(handle, delimiter="\t")))
require(counts == {"LUAD": 158, "LUSC": 193}, f"Unexpected TF categories: {counts}")
checks.append({"check": "Figure1C:158_LUAD_193_LUSC_TFs", "status": "passed"})

def group_counts(path: Path) -> dict[str, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        return dict(Counter(row["group"] for row in csv.DictReader(handle, delimiter="\t")))


require(group_counts(processed / "tcga_manifest.tsv") == {"LUAD": 516, "LUSC": 501}, "TCGA flow counts changed")
require(group_counts(processed / "gse81089_manifest.tsv") == {"LUAD": 106, "LUSC": 67}, "GSE81089 flow counts changed")
checks.append({"check": "SupplementaryFigure1:manifest_counts", "status": "passed"})

with (tables / "Figure1_anchor_style_render_receipt.tsv").open(newline="", encoding="utf-8") as handle:
    receipt = list(csv.DictReader(handle, delimiter="\t"))
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
observed = {row["metric"]: float(row["value"]) for row in receipt}
for key, expected in required_flow.items():
    require(int(round(observed[key])) == expected, f"Flow count mismatch for {key}: {observed[key]}")
checks.append({"check": "SupplementaryFigure1:flow_arithmetic", "status": "passed"})

verification = {
    "status": "passed",
    "outputs": stems,
    "figure1c_shape": list(matrix_shape),
    "max_abs_row_mean": max_abs_row_mean,
    "max_abs_row_sd_minus_1": max_abs_row_sd_minus_1,
    "tf_categories": counts,
    "checks": checks,
}
(audit / "verification.json").write_text(json.dumps(verification, indent=2), encoding="utf-8")
print(json.dumps(verification, indent=2))
