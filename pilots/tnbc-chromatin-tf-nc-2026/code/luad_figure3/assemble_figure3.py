#!/usr/bin/env python3
"""Assemble vector Figure 3 and its supplementary heterogeneity panels."""

from __future__ import annotations

import argparse
from pathlib import Path

import fitz


def place(page: fitz.Page, path: Path, box: tuple[float, float, float, float]) -> None:
    source = fitz.open(path)
    page.show_pdf_page(fitz.Rect(*box), source, 0, keep_proportion=True)
    source.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    figures = args.root.resolve() / "results" / "figures"
    width, height = 1728.0, 1152.0
    document = fitz.open()
    page = document.new_page(width=width, height=height)
    placements = (
        ("Figure3A_regulon_composition.pdf", (0, 0, 520, 650)),
        ("Figure3B_module_GO_enrichment.pdf", (520, 0, 1110, 550)),
        ("Figure3C_HC_TF_activity_correlation.pdf", (1110, 0, 1728, 560)),
        ("Figure3D_TF_collaboration_network.pdf", (500, 545, 1085, 1152)),
        ("Figure3E_regulon_size_partners.pdf", (0, 650, 505, 1152)),
        ("Figure3F_HC_TF_activity_skewness.pdf", (1080, 550, 1405, 900)),
        ("Figure3G_representative_distributions.pdf", (1400, 550, 1728, 900)),
    )
    for name, box in placements:
        place(page, figures / name, box)
    output = figures / "Figure3_complete.pdf"
    document.save(output, garbage=4, deflate=True)
    document.close()

    supplement = fitz.open()
    page = supplement.new_page(width=1000, height=720)
    place(page, figures / "SupplementaryFigure6A_MKI67_correlations.pdf", (0, 0, 500, 720))
    place(page, figures / "SupplementaryFigure6B_cross_system_skewness.pdf", (500, 0, 1000, 720))
    supplement.save(figures / "SupplementaryFigure6_complete.pdf", garbage=4, deflate=True)
    supplement.close()

    for path in (figures / "Figure3_complete.pdf", figures / "SupplementaryFigure6_complete.pdf"):
        pdf = fitz.open(path)
        pixmap = pdf[0].get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False)
        pixmap.save(path.with_suffix(".png"))
        pdf.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

