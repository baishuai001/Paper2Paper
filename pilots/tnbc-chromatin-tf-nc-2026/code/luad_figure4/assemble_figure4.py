#!/usr/bin/env python3
"""Vector assembly for the LUAD Figure 4 survival package."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import fitz


def place(page: fitz.Page, path: Path, box: tuple[float, float, float, float]) -> None:
    source = fitz.open(path)
    # survminer/ggsurvplot writes a nominally blank first page followed by the
    # composed survival plot on page 2.  The blank page can still contain one
    # background drawing, so "first page with drawings" is not sufficient.
    # Select the page with the largest amount of semantic/graphic content.
    def content_score(candidate: fitz.Page) -> tuple[int, int, int]:
        return (
            len(candidate.get_text().strip()),
            len(candidate.get_drawings()),
            len(candidate.get_images(full=True)),
        )

    source_page = max(range(len(source)), key=lambda index: content_score(source[index]))
    page.show_pdf_page(fitz.Rect(*box), source, source_page, keep_proportion=True)
    source.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    figures = root / "results" / "figures"
    tables = root / "results" / "tables"
    with (tables / "Figure4FH_KM_receipts.tsv").open("rt", encoding="utf-8", newline="") as handle:
        km = list(csv.DictReader(handle, delimiter="\t"))

    document = fitz.open()
    page = document.new_page(width=2400, height=1850)
    forests = (
        "Figure4A_GSE41271_LUAD_OS_multivariable_forest.pdf",
        "Figure4B_GSE41271_LUAD_RFS_multivariable_forest.pdf",
        "Figure4C_TCGA_LUAD_OS_multivariable_forest.pdf",
        "Figure4D_TCGA_LUAD_RFS_multivariable_forest.pdf",
    )
    for index, name in enumerate(forests):
        place(page, figures / name, (index * 600, 0, (index + 1) * 600, 760))
    place(page, figures / "Figure4E_GSE41271_endpoint_overlap.pdf", (0, 750, 1200, 1110))
    place(page, figures / "Figure4G_TCGA_endpoint_overlap.pdf", (1200, 750, 2400, 1110))

    for cohort, x0, x1 in (("GSE41271_LUAD", 0, 1200), ("TCGA_LUAD", 1200, 2400)):
        cohort_rows = [row for row in km if row["cohort"] == cohort]
        cohort_rows.sort(key=lambda row: (row["TF"], row["endpoint"]))
        for index, row in enumerate(cohort_rows):
            col = index % 2
            line = index // 2
            width = (x1 - x0) / 2
            box = (x0 + col * width, 1090 + line * 380,
                   x0 + (col + 1) * width, 1090 + (line + 1) * 380)
            place(page, figures / f"{row['filename']}.pdf", box)
    main_path = figures / "Figure4_complete.pdf"
    document.save(main_path, garbage=4, deflate=True)
    document.close()

    supplement = fitz.open()
    page = supplement.new_page(width=1400, height=780)
    place(page, figures / "SupplementaryFigure7_permutation_nulls.pdf", (0, 0, 1400, 780))
    supp_path = figures / "SupplementaryFigure7_complete.pdf"
    supplement.save(supp_path, garbage=4, deflate=True)
    supplement.close()

    for path in (main_path, supp_path):
        pdf = fitz.open(path)
        pix = pdf[0].get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False)
        pix.save(path.with_suffix(".png"))
        pdf.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
