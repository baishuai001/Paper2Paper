#!/usr/bin/env python3
"""Assemble Figure 1 and Supplementary Figures 1-2 without rasterizing PDFs."""

from __future__ import annotations

import argparse
from pathlib import Path

import fitz


def place(page: fitz.Page, figures: Path, name: str, box: tuple[float, float, float, float]) -> None:
    source = fitz.open(figures / name)
    page.show_pdf_page(fitz.Rect(*box), source, 0, keep_proportion=True)
    source.close()


def save_preview(path: Path, page_number: int = 0) -> None:
    pdf = fitz.open(path)
    pix = pdf[page_number].get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False)
    suffix = "" if len(pdf) == 1 else f"_page{page_number + 1}"
    pix.save(path.with_name(path.stem + suffix + ".png"))
    pdf.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    figures = args.root.resolve() / "results" / "figures"

    main = fitz.open()
    page = main.new_page(width=1800, height=1400)
    place(page, figures, "Figure1A_workflow_final.pdf", (0, 0, 710, 730))
    place(page, figures, "Figure1B_TCGA_logFC_msviper.pdf", (710, 0, 1180, 540))
    place(page, figures, "Figure1C_TCGA_TF_activity_annotated.pdf", (1180, 0, 1800, 620))
    place(page, figures, "Figure1D_PDMR_PDX_TF_activity_annotated.pdf", (0, 720, 850, 1400))
    place(page, figures, "Figure1E_DepMap_22Q2_TF_activity_annotated.pdf", (850, 620, 1800, 1400))
    main_path = figures / "Figure1_complete.pdf"
    main.save(main_path, garbage=4, deflate=True)
    main.close()

    supp1 = fitz.open()
    page = supp1.new_page(width=1500, height=780)
    place(page, figures, "SupplementaryFigure1A_TCGA_inclusion.pdf", (0, 0, 500, 780))
    place(page, figures, "SupplementaryFigure1B_GSE81089_inclusion.pdf", (500, 0, 1000, 780))
    place(page, figures, "SupplementaryFigure1C_GSE41271_inclusion.pdf", (1000, 0, 1500, 780))
    supp1_path = figures / "SupplementaryFigure1_complete.pdf"
    supp1.save(supp1_path, garbage=4, deflate=True)
    supp1.close()

    supp2 = fitz.open()
    page1 = supp2.new_page(width=1500, height=1050)
    place(page1, figures, "SupplementaryFigure2A_GSE81089_logFC_msviper.pdf", (0, 0, 500, 525))
    place(page1, figures, "SupplementaryFigure2B_TF_overlap.pdf", (500, 0, 1000, 525))
    place(page1, figures, "SupplementaryFigure2C_TCGA_sample_correlations.pdf", (1000, 0, 1500, 525))
    place(page1, figures, "SupplementaryFigure2D_PDMR_PDX_sample_correlations.pdf", (0, 525, 750, 1050))
    place(page1, figures, "SupplementaryFigure2E_DepMap_22Q2_sample_correlations.pdf", (750, 525, 1500, 1050))
    page2 = supp2.new_page(width=1500, height=720)
    place(page2, figures, "SupplementaryFigure2F_replicated_LUAD_TF_effects.pdf", (0, 0, 500, 720))
    place(page2, figures, "SupplementaryFigure2G_GSE41271_logFC_msviper.pdf", (500, 0, 1000, 720))
    place(page2, figures, "SupplementaryFigure2H_patient_cohort_TF_overlap.pdf", (1000, 0, 1500, 720))
    supp2_path = figures / "SupplementaryFigure2_complete.pdf"
    supp2.save(supp2_path, garbage=4, deflate=True)
    supp2.close()

    save_preview(main_path)
    save_preview(supp1_path)
    save_preview(supp2_path, 0)
    save_preview(supp2_path, 1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

