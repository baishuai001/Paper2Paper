#!/usr/bin/env python3
"""Verify the strict official-code Supplementary Figure 3 cloud run."""

from __future__ import annotations

import csv
import json
import struct
import subprocess
import sys
from pathlib import Path


EXPECTED = {
    "SupplementaryFigure3A": "SupplementaryFigure3A_TCGA_Number_Of_Called_Peaks.pdf",
    "SupplementaryFigure3B": "SupplementaryFigure3B_PDX_Number_Of_Called_Peaks.pdf",
    "SupplementaryFigure3C": "SupplementaryFigure3C_CellLines_Number_Of_Called_Peaks.pdf",
    "SupplementaryFigure3D": "SupplementaryFigure3D_TCGA_peak_saturation.pdf",
    "SupplementaryFigure3E": "SupplementaryFigure3E_PDX_peak_saturation.pdf",
    "SupplementaryFigure3F": "SupplementaryFigure3F_CellLines_peak_saturation.pdf",
    "SupplementaryFigure3G": "SupplementaryFigure3G_TCGA_Pearson_Correlation_Signal_ConsensusPeakSet.pdf",
    "SupplementaryFigure3H": "SupplementaryFigure3H_PDX_Pearson_Correlation_Signal_ConsensusPeakSet.pdf",
    "SupplementaryFigure3I": "SupplementaryFigure3I_CellLines_Pearson_Correlation_Signal_ConsensusPeakSet.pdf",
    "SupplementaryFigure3J": "SupplementaryFigure3J_TCGA_GenomicAnnotation.pdf",
    "SupplementaryFigure3K": "SupplementaryFigure3K_PDX_GenomicAnnotation.pdf",
    "SupplementaryFigure3L": "SupplementaryFigure3L_CellLines_GenomicAnnotation.pdf",
}


def read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def png_dimensions(path: Path):
    with path.open("rb") as fh:
        header = fh.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise RuntimeError(f"invalid PNG: {path}")
    return struct.unpack(">II", header[16:24])


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: verify_official_supp3.py WORK_ROOT PANEL_STATUS_TSV")
    root = Path(sys.argv[1]).resolve()
    status_path = Path(sys.argv[2]).resolve()
    visual_qa_path = status_path.parent / "VISUAL_QA.tsv"
    panels = root / "results" / "panels"
    png_dir = root / "results" / "rendered_png"
    audit = root / "runtime" / "audit"

    statuses = read_tsv(status_path)
    if {r["panel"] for r in statuses} != set(EXPECTED):
        raise RuntimeError("panel status table does not cover A-L exactly")
    visual_qa = read_tsv(visual_qa_path)
    if {r["panel"] for r in visual_qa} != set(EXPECTED):
        raise RuntimeError("visual QA table does not cover A-L exactly")

    materialization = json.loads((audit / "materialization_receipt.json").read_text())
    if materialization["plotting_constructor_rewrites"] != 0:
        raise RuntimeError("plot constructor rewrite detected")
    if materialization["official_sources_verified"] != 5:
        raise RuntimeError("official source lock incomplete")
    substitutions = read_tsv(audit / "runtime_substitution_manifest.tsv")
    if len(substitutions) != 4:
        raise RuntimeError(f"expected four path-only substitutions, found {len(substitutions)}")
    if {r["change_kind"] for r in substitutions} - {"PATH_ONLY", "PATH_ARGUMENT_TYPO_ONLY"}:
        raise RuntimeError("unapproved runtime substitution")

    cohort_receipt = read_tsv(root / "adapter" / "audit" / "cohort_receipt.tsv")
    observed_n = {r["cohort"]: int(r["n_samples"]) for r in cohort_receipt}
    if observed_n != {"TCGA": 22, "PDX": 13, "CellLines": 19}:
        raise RuntimeError(f"cohort counts changed: {observed_n}")
    genomic = read_tsv(root / "adapter" / "audit" / "genomic_annotation_adapter_receipt.tsv")
    nonzero_n = {r["cohort"]: int(r["n_nonzero_samples"]) for r in genomic}
    if nonzero_n != {"TCGA": 22, "PDX": 12, "CellLines": 19}:
        raise RuntimeError(f"nonzero genomic-annotation sample counts changed: {nonzero_n}")
    baseline = read_tsv(root / "adapter" / "audit" / "genome_baseline_receipt.tsv")
    if len(baseline) != 1 or int(baseline[0]["bins"]) != 6062095:
        raise RuntimeError("the real 6,062,095-bin hg38 500-bp baseline was not verified")

    results = []
    for panel, filename in EXPECTED.items():
        pdf = panels / filename
        if not pdf.is_file() or pdf.stat().st_size < 5000:
            raise RuntimeError(f"missing or undersized PDF: {pdf}")
        with pdf.open("rb") as fh:
            if fh.read(4) != b"%PDF":
                raise RuntimeError(f"invalid PDF header: {pdf}")
        info = subprocess.run(["pdfinfo", str(pdf)], check=True, capture_output=True, text=True).stdout
        pages_line = next((line for line in info.splitlines() if line.startswith("Pages:")), "")
        if pages_line.split()[-1:] != ["1"]:
            raise RuntimeError(f"expected one-page PDF: {pdf}")
        png = png_dir / f"{pdf.stem}.png"
        if not png.is_file() or png.stat().st_size < 10000:
            raise RuntimeError(f"missing or undersized render: {png}")
        width, height = png_dimensions(png)
        if min(width, height) < 400:
            raise RuntimeError(f"render resolution too small: {png} {width}x{height}")
        status = next(r["status"] for r in statuses if r["panel"] == panel)
        results.append(
            {
                "panel": panel,
                "status": status,
                "pdf": str(pdf),
                "pdf_bytes": pdf.stat().st_size,
                "render_png": str(png),
                "render_width": width,
                "render_height": height,
            }
        )

    visual_warnings = [r for r in visual_qa if r["qa_status"].startswith("OFFICIAL_")]
    verification = {
        "verdict": "PASS_STRICT_PORT_WITH_DECLARED_INPUT_AND_VISUAL_BOUNDARIES",
        "publication_ready": False,
        "official_plotting_code_used": True,
        "plotting_constructor_rewrites": 0,
        "panels_verified": len(results),
        "partial_panels": ["SupplementaryFigure3A"],
        "idr_not_applicable_panels": ["SupplementaryFigure3B", "SupplementaryFigure3C"],
        "real_hg38_500bp_baseline_used": True,
        "real_hg38_500bp_baseline_bins": 6062095,
        "visual_warnings": visual_warnings,
        "results": results,
    }
    out = root / "results" / "verification_receipt.json"
    out.write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
