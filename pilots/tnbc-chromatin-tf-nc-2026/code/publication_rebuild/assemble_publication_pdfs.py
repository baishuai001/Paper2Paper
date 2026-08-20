#!/usr/bin/env python3
"""Assemble the LUAD main figures, supplements, and the TNBC comparison atlas.

The script uses PyMuPDF's ``show_pdf_page`` so that vector panels stay vector.
No scientific statistic is recomputed here; every source is a frozen PDF and
is recorded in a machine-readable manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
from typing import Iterable, Sequence

import fitz


PAGE_W = 13.0 * 72
PAGE_H = 10.0 * 72
MARGIN = 18
INK = (0.10, 0.12, 0.14)
MUTED = (0.32, 0.36, 0.39)
LIGHT = (0.93, 0.94, 0.95)
BLUE = (0.16, 0.47, 0.71)
RUST = (0.77, 0.35, 0.22)


def rect(x0: float, y0: float, x1: float, y1: float) -> fitz.Rect:
    return fitz.Rect(x0, y0, x1, y1)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def title(page: fitz.Page, text: str, subtitle: str | None = None) -> None:
    page.insert_text((MARGIN, 19), text, fontsize=11, fontname="helv", color=INK)
    if subtitle:
        page.insert_text((MARGIN, 31), subtitle, fontsize=6.5, fontname="helv", color=MUTED)
    page.draw_line((MARGIN, 35), (PAGE_W - MARGIN, 35), color=(0.78, 0.80, 0.82), width=0.45)


def footer(page: fitz.Page, text: str) -> None:
    page.insert_textbox(
        rect(MARGIN, PAGE_H - 15, PAGE_W - MARGIN, PAGE_H - 3),
        text,
        fontsize=5.7,
        fontname="helv",
        color=MUTED,
        align=fitz.TEXT_ALIGN_RIGHT,
    )


def draw_missing(page: fitz.Page, destination: fitz.Rect, path: Path) -> None:
    page.draw_rect(destination, color=(0.75, 0.30, 0.30), fill=(0.98, 0.94, 0.94), width=0.8)
    page.insert_textbox(
        destination + (6, 6, -6, -6),
        f"MISSING SOURCE\n{path}",
        fontsize=7,
        color=(0.55, 0.10, 0.10),
        align=fitz.TEXT_ALIGN_CENTER,
    )


def normalized_clip(page: fitz.Page, values: Sequence[float]) -> fitz.Rect:
    x0, y0, x1, y1 = values
    r = page.rect
    return fitz.Rect(r.x0 + x0 * r.width, r.y0 + y0 * r.height,
                     r.x0 + x1 * r.width, r.y0 + y1 * r.height)


def place_pdf(
    page: fitz.Page,
    source_path: Path,
    destination: fitz.Rect,
    source_page: int = 0,
    clip_norm: Sequence[float] | None = None,
    border: bool = False,
) -> bool:
    if not source_path.exists():
        draw_missing(page, destination, source_path)
        return False
    if source_path.suffix.lower() in {".png", ".jpg", ".jpeg", ".tif", ".tiff"}:
        try:
            page.insert_image(destination, filename=str(source_path),
                              keep_proportion=True, overlay=True)
        except Exception:
            draw_missing(page, destination, source_path)
            return False
        if border:
            page.draw_rect(destination, color=(0.83, 0.84, 0.85), width=0.45)
        return True
    with fitz.open(source_path) as source:
        if source_page >= source.page_count:
            draw_missing(page, destination, source_path)
            return False
        clip = normalized_clip(source[source_page], clip_norm) if clip_norm else None
        page.show_pdf_page(destination, source, pno=source_page, clip=clip,
                           keep_proportion=True, overlay=True)
    if border:
        page.draw_rect(destination, color=(0.83, 0.84, 0.85), width=0.45)
    return True


def make_page(doc: fitz.Document, heading: str, subtitle: str | None = None) -> fitz.Page:
    page = doc.new_page(width=PAGE_W, height=PAGE_H)
    page.draw_rect(page.rect, fill=(1, 1, 1), color=None)
    title(page, heading, subtitle)
    return page


def save(doc: fitz.Document, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    doc.set_metadata({
        "title": path.stem,
        "author": "Paper2Paper LUAD publication rebuild",
        "subject": "TNBC-anchor-informed LUAD regulatory analysis",
    })
    doc.save(path, garbage=4, deflate=True, clean=True)
    doc.close()


def compose_block(output: Path, heading: str, placements: Sequence[tuple[Path, fitz.Rect]]) -> Path:
    doc = fitz.open()
    page = make_page(doc, heading)
    for src, dst in placements:
        place_pdf(page, src, dst)
    footer(page, "Composite panel; all component paths are recorded in the publication manifest")
    save(doc, output)
    return output


def main_figure(
    output: Path,
    heading: str,
    subtitle: str,
    placements: Sequence[tuple[str, Path, fitz.Rect]],
    manifest_rows: list[dict[str, str]],
) -> None:
    doc = fitz.open()
    page = make_page(doc, heading, subtitle)
    for panel, src, dst in placements:
        ok = place_pdf(page, src, dst)
        manifest_rows.append({
            "artifact": output.name,
            "panel": panel,
            "source": str(src),
            "source_exists": str(ok),
            "destination": ",".join(f"{v:.1f}" for v in dst),
        })
    footer(page, "LUAD transfer; visual hierarchy follows the anchor while colors, geometry, and annotations are independently redrawn")
    save(doc, output)


def append_scaled_pages(doc: fitz.Document, source: Path, heading: str, manifest_rows: list[dict[str, str]]) -> None:
    if not source.exists():
        page = make_page(doc, heading)
        draw_missing(page, rect(24, 48, PAGE_W - 24, PAGE_H - 24), source)
        return
    with fitz.open(source) as src_doc:
        page_count = src_doc.page_count
    for pno in range(page_count):
        page = make_page(doc, heading if page_count == 1 else f"{heading} | page {pno + 1}/{page_count}")
        place_pdf(page, source, rect(22, 44, PAGE_W - 22, PAGE_H - 22), source_page=pno)
        footer(page, "Canonical LUAD supplementary numbering; source panel retained as vector PDF")
        manifest_rows.append({
            "artifact": "LUAD_Supplementary_Figures_1-9.pdf",
            "panel": heading,
            "source": str(source),
            "source_exists": "True",
            "destination": "full-page",
        })


def write_labeled_box(page: fitz.Page, box: fitz.Rect, label: str, text: str, accent: tuple[float, float, float]) -> None:
    page.draw_rect(box, fill=(0.975, 0.978, 0.981), color=(0.84, 0.85, 0.86), width=0.45)
    page.draw_rect(rect(box.x0, box.y0, box.x0 + 4, box.y1), fill=accent, color=None)
    page.insert_text((box.x0 + 10, box.y0 + 13), label, fontsize=7.2, fontname="helv", color=accent)
    page.insert_textbox(rect(box.x0 + 10, box.y0 + 18, box.x1 - 7, box.y1 - 5),
                        text, fontsize=6.8, lineheight=1.22, fontname="helv", color=INK)


def comparison_atlas(
    output: Path,
    anchor_pdf: Path,
    anchor_supp_pdf: Path,
    gap_matrix: Path,
    luad_sources: dict[tuple[str, str], Path],
    anchor_crops: dict[tuple[str, str], tuple[int, Sequence[float]]],
) -> None:
    rows: list[dict[str, str]] = []
    with gap_matrix.open("r", encoding="utf-8-sig", newline="") as handle:
        rows.extend(csv.DictReader(handle, delimiter="\t"))
    doc = fitz.open()
    for index, row in enumerate(rows, 1):
        figure = row["figure"]
        panel = row["panel"]
        key = (figure, panel)
        page = make_page(
            doc,
            f"{figure}{panel} | TNBC anchor versus LUAD transfer",
            f"Panel {index}/{len(rows)} | side-by-side scientific and visual audit",
        )
        left = rect(20, 55, 454, 333)
        right = rect(482, 55, 916, 333)
        page.insert_text((left.x0, 49), "TNBC anchor", fontsize=7.2, color=RUST)
        page.insert_text((right.x0, 49), "LUAD transfer", fontsize=7.2, color=BLUE)
        if key in anchor_crops:
            pno, clip = anchor_crops[key]
            anchor_source = anchor_supp_pdf if figure == "SupplementaryFigure" else anchor_pdf
            place_pdf(page, anchor_source, left, source_page=pno, clip_norm=clip, border=True)
        else:
            draw_missing(page, left, anchor_pdf)
        place_pdf(page, luad_sources.get(key, Path("MISSING")), right, border=True)

        write_labeled_box(page, rect(20, 350, 454, 438), "TNBC biological meaning",
                          row["tnbc_biological_meaning"], RUST)
        write_labeled_box(page, rect(482, 350, 916, 438), "LUAD biological meaning",
                          row["luad_current_meaning"], BLUE)
        write_labeled_box(page, rect(20, 451, 454, 542), "Interpretive difference",
                          row["interpretive_gap"], (0.38, 0.27, 0.55))
        write_labeled_box(page, rect(482, 451, 916, 542), "Visual deficit in the first pass",
                          row["current_visual_deficit"], (0.73, 0.39, 0.10))
        write_labeled_box(page, rect(20, 555, 916, 645), "Required / implemented improvement",
                          row["required_improvement"], (0.10, 0.48, 0.34))
        page.insert_textbox(rect(20, 658, 916, 692),
                            f"Evidence status: {row['evidence_status']}",
                            fontsize=7.2, fontname="helv", color=INK)
        footer(page, "Anchor crop is reproduced only for methodological critique; LUAD panel is independently generated from frozen LUAD results")
    save(doc, output)


def build(args: argparse.Namespace) -> None:
    fig1 = Path(args.fig1)
    fig2 = Path(args.fig2)
    fig3 = Path(args.fig3)
    fig4 = Path(args.fig4)
    fig5 = Path(args.fig5)
    rebuild = Path(args.rebuild)
    atom = rebuild / "results" / "atomic_figures"
    publication = rebuild / "results" / "publication"
    publication.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str]] = []

    f1 = fig1 / "results" / "figures"
    f2 = fig2 / "results" / "figures"
    f3 = fig3 / "results" / "figures"
    f4 = fig4 / "results" / "figures"
    f5 = fig5 / "results" / "figures"

    # Composite KM blocks are useful both in Figure 4 and in the panel atlas.
    km_gse = compose_block(
        atom / "Figure4F_GSE41271_KM_block.pdf",
        "F  GSE41271-LUAD representative Kaplan-Meier analyses",
        [
            (f4 / "Figure4F_GSE41271_LUAD_ETV1_OS_KM.png", rect(18, 42, 460, 365)),
            (f4 / "Figure4F_GSE41271_LUAD_ETV1_RFS_KM.png", rect(476, 42, 918, 365)),
            (f4 / "Figure4F_GSE41271_LUAD_ZNF444_OS_KM.png", rect(18, 380, 460, 700)),
            (f4 / "Figure4F_GSE41271_LUAD_ZNF444_RFS_KM.png", rect(476, 380, 918, 700)),
        ],
    )
    km_tcga = compose_block(
        atom / "Figure4H_TCGA_KM_block.pdf",
        "H  TCGA-LUAD representative Kaplan-Meier analyses",
        [
            (f4 / "Figure4H_TCGA_LUAD_CREBRF_OS_KM.png", rect(18, 42, 460, 365)),
            (f4 / "Figure4H_TCGA_LUAD_CREBRF_RFS_KM.png", rect(476, 42, 918, 365)),
            (f4 / "Figure4H_TCGA_LUAD_PHC2_OS_KM.png", rect(18, 380, 460, 700)),
            (f4 / "Figure4H_TCGA_LUAD_PHC2_RFS_KM.png", rect(476, 380, 918, 700)),
        ],
    )
    supp4_compare = compose_block(
        atom / "SupplementaryFigure4_comparison_block.pdf",
        "Supplementary Figure 4 | HC-TF gate, motif availability, and proliferation",
        [
            (f2 / "SupplementaryFigure4A-C_complete.pdf", rect(18, 42, 575, 700)),
            (atom / "SupplementaryFigure4D_MKI67_four_system.pdf", rect(590, 42, 918, 700)),
        ],
    )
    supp5_compare = compose_block(
        atom / "SupplementaryFigure5_comparison_block.pdf",
        "Supplementary Figure 5 | Heterogeneity and microenvironment confounding",
        [
            (atom / "SupplementaryFigure5A_cross_system_skewness.pdf", rect(18, 42, 460, 700)),
            (atom / "SupplementaryFigure5B_ESTIMATE_correlations.pdf", rect(476, 42, 918, 700)),
        ],
    )
    supp6_compare = compose_block(
        atom / "SupplementaryFigure6_comparison_block.pdf",
        "Supplementary Figure 6 | Complete survival and permutation screens",
        [
            (atom / "SupplementaryFigure6A-H_Cox_volcanoes.pdf", rect(18, 42, 918, 365)),
            (atom / "SupplementaryFigure6I-P_permutation_nulls.pdf", rect(18, 380, 918, 700)),
        ],
    )
    supp7_compare = compose_block(
        atom / "SupplementaryFigure7_comparison_block.pdf",
        "Supplementary Figure 7 | Expanded pharmacogenomic evidence",
        [
            (atom / "SupplementaryFigure7A_drug_overlap.pdf", rect(18, 42, 330, 320)),
            (atom / "SupplementaryFigure7B_replicated_pair_scatters.pdf", rect(342, 42, 918, 470)),
            (atom / "SupplementaryFigure7C_complete_bubble_matrix.pdf", rect(18, 482, 918, 700)),
        ],
    )

    main_figure(
        publication / "Figure1_revised.pdf",
        "Figure 1 | A transferable LUAD TF-activity framework",
        "TCGA discovery, two patient validations, and independent PDX/cell-line projection",
        [
            ("A", f1 / "Figure1A_workflow_final.pdf", rect(20, 45, 290, 315)),
            ("B", f1 / "Figure1B_TCGA_logFC_msviper.pdf", rect(300, 45, 560, 315)),
            ("C", f1 / "Figure1C_TCGA_TF_activity_annotated.pdf", rect(570, 45, 916, 315)),
            ("D", f1 / "Figure1D_PDMR_PDX_TF_activity_annotated.pdf", rect(20, 330, 463, 700)),
            ("E", f1 / "Figure1E_DepMap_22Q2_TF_activity_annotated.pdf", rect(473, 330, 916, 700)),
        ],
        manifest_rows,
    )
    main_figure(
        publication / "Figure2_revised.pdf",
        "Figure 2 | Chromatin prioritization and selective motif conservation",
        "All 158 discovery TFs enter the gate; 31 HC-TFs pass promoter/activity criteria",
        [
            ("A", f2 / "Figure2A_design.pdf", rect(20, 45, 245, 255)),
            ("B", f2 / "Figure2B_promoter_upset.pdf", rect(255, 45, 465, 255)),
            ("D", f2 / "Figure2D_motif_upset.pdf", rect(475, 45, 675, 255)),
            ("E", atom / "Figure2E_cross_system_motif_prevalence.pdf", rect(685, 45, 916, 255)),
            ("C", atom / "Figure2C_compact_HC_TF_evidence.pdf", rect(20, 270, 916, 700)),
        ],
        manifest_rows,
    )
    main_figure(
        publication / "Figure3_revised.pdf",
        "Figure 3 | Network architecture and intertumoral regulatory heterogeneity",
        "All 31 HC-TFs are analyzed; motif support is an annotation rather than a downstream inclusion rule",
        [
            ("A", f3 / "Figure3A_regulon_composition.pdf", rect(20, 45, 285, 355)),
            ("B", atom / "Figure3B_target_gene_community_network.pdf", rect(295, 45, 605, 355)),
            ("C", f3 / "Figure3C_HC_TF_activity_correlation.pdf", rect(615, 45, 916, 355)),
            ("D", f3 / "Figure3D_TF_collaboration_network.pdf", rect(20, 370, 300, 700)),
            ("E", f3 / "Figure3E_regulon_size_partners.pdf", rect(310, 370, 505, 700)),
            ("F", f3 / "Figure3F_HC_TF_activity_skewness.pdf", rect(515, 370, 710, 700)),
            ("G", f3 / "Figure3G_representative_distributions.pdf", rect(720, 370, 916, 700)),
        ],
        manifest_rows,
    )
    main_figure(
        publication / "Figure4_revised.pdf",
        "Figure 4 | Exploratory clinical associations of the LUAD HC-TF network",
        "Nominal and BH-adjusted evidence are visually separated; permutation support is reported in Supplementary Figure 6",
        [
            ("A", atom / "Figure4A_compact_forest.pdf", rect(20, 45, 235, 245)),
            ("B", atom / "Figure4B_compact_forest.pdf", rect(245, 45, 460, 245)),
            ("C", atom / "Figure4C_compact_forest.pdf", rect(470, 45, 685, 245)),
            ("D", atom / "Figure4D_compact_forest.pdf", rect(695, 45, 916, 245)),
            ("E", f4 / "Figure4E_GSE41271_endpoint_overlap.pdf", rect(20, 260, 455, 405)),
            ("G", f4 / "Figure4G_TCGA_endpoint_overlap.pdf", rect(481, 260, 916, 405)),
            ("F", km_gse, rect(20, 418, 455, 700)),
            ("H", km_tcga, rect(481, 418, 916, 700)),
        ],
        manifest_rows,
    )
    main_figure(
        publication / "Figure5_revised.pdf",
        "Figure 5 | Replicated pharmacogenomic associations and the in-vivo validation gap",
        "Five cell-line associations replicate across datasets; no pair satisfies the matched public PDX validation rule",
        [
            ("A", f5 / "Figure5A_GDSC2_volcano.pdf", rect(20, 45, 315, 250)),
            ("B", f5 / "Figure5B_CTRPv2_volcano.pdf", rect(321, 45, 615, 250)),
            ("C", f5 / "Figure5C_PRISM_volcano.pdf", rect(621, 45, 916, 250)),
            ("D", atom / "Figure5D_true_upset.pdf", rect(20, 265, 305, 470)),
            ("E", f5 / "Figure5E_replicated_association_heatmap.pdf", rect(315, 265, 600, 470)),
            ("F", atom / "Figure5F_pathway_TF_bubble.pdf", rect(610, 265, 916, 470)),
            ("G", atom / "Figure5G_PDX_validation_funnel.pdf", rect(20, 485, 460, 700)),
            ("H", atom / "Figure5H_PDX_complete_screen.pdf", rect(475, 485, 916, 700)),
        ],
        manifest_rows,
    )

    # Canonical supplementary figures. Dense figures are allowed to span pages.
    supp = fitz.open()
    append_scaled_pages(supp, f1 / "SupplementaryFigure1_complete.pdf", "Supplementary Figure 1 | Cohort inclusion", manifest_rows)
    append_scaled_pages(supp, f1 / "SupplementaryFigure2_complete.pdf", "Supplementary Figure 2 | Independent TF discovery and transfer", manifest_rows)
    append_scaled_pages(supp, f2 / "SupplementaryFigure3_complete.pdf", "Supplementary Figure 3 | ATAC quality, saturation, and genomic distribution", manifest_rows)

    page = make_page(supp, "Supplementary Figure 4 | HC-TF selection, motif availability, and proliferation")
    place_pdf(page, f2 / "SupplementaryFigure4A-C_complete.pdf", rect(20, 45, 570, 700))
    place_pdf(page, atom / "SupplementaryFigure4D_MKI67_four_system.pdf", rect(580, 45, 916, 700))
    footer(page, "Panels A-C reproduce the frozen gate; panel D tests four anchor-equivalent expression systems")

    page = make_page(supp, "Supplementary Figure 5 | Regulatory heterogeneity and TME confounding")
    place_pdf(page, atom / "SupplementaryFigure5A_cross_system_skewness.pdf", rect(20, 45, 455, 700))
    place_pdf(page, atom / "SupplementaryFigure5B_ESTIMATE_correlations.pdf", rect(465, 45, 916, 700))
    footer(page, "Skewness and ESTIMATE correlations use all 31 HC-TFs with within-analysis BH correction")

    append_scaled_pages(supp, atom / "SupplementaryFigure6A-H_Cox_volcanoes.pdf", "Supplementary Figure 6A-H | Cox screens", manifest_rows)
    append_scaled_pages(supp, atom / "SupplementaryFigure6I-P_permutation_nulls.pdf", "Supplementary Figure 6I-P | Permutation tests", manifest_rows)

    page = make_page(supp, "Supplementary Figure 7A-B | Drug coverage and replicated examples")
    place_pdf(page, atom / "SupplementaryFigure7A_drug_overlap.pdf", rect(20, 45, 350, 330))
    place_pdf(page, atom / "SupplementaryFigure7B_replicated_pair_scatters.pdf", rect(360, 45, 916, 700))
    footer(page, "Every dataset contributing to a replicated pair is shown; no illustrative pair is omitted")
    append_scaled_pages(supp, atom / "SupplementaryFigure7C_complete_bubble_matrix.pdf", "Supplementary Figure 7C | Complete association matrix", manifest_rows)
    append_scaled_pages(supp, f5 / "SupplementaryFigure8_complete.pdf", "Extended Supplementary Figure 8 | Cell-line and classifier audit", manifest_rows)
    append_scaled_pages(supp, f5 / "SupplementaryFigure9_complete.pdf", "Extended Supplementary Figure 9 | PDX coverage and negative screen", manifest_rows)
    save(supp, publication / "LUAD_Supplementary_Figures_1-9.pdf")

    # Panel sources for the 33-page comparison atlas.
    luad_sources = {
        ("Figure1", "A"): f1 / "Figure1A_workflow_final.pdf",
        ("Figure1", "B"): f1 / "Figure1B_TCGA_logFC_msviper.pdf",
        ("Figure1", "C"): f1 / "Figure1C_TCGA_TF_activity_annotated.pdf",
        ("Figure1", "D"): f1 / "Figure1D_PDMR_PDX_TF_activity_annotated.pdf",
        ("Figure1", "E"): f1 / "Figure1E_DepMap_22Q2_TF_activity_annotated.pdf",
        ("Figure2", "A"): f2 / "Figure2A_design.pdf",
        ("Figure2", "B"): f2 / "Figure2B_promoter_upset.pdf",
        ("Figure2", "C"): atom / "Figure2C_compact_HC_TF_evidence.pdf",
        ("Figure2", "D"): f2 / "Figure2D_motif_upset.pdf",
        ("Figure2", "E"): atom / "Figure2E_cross_system_motif_prevalence.pdf",
        ("Figure3", "A"): f3 / "Figure3A_regulon_composition.pdf",
        ("Figure3", "B"): atom / "Figure3B_target_gene_community_network.pdf",
        ("Figure3", "C"): f3 / "Figure3C_HC_TF_activity_correlation.pdf",
        ("Figure3", "D"): f3 / "Figure3D_TF_collaboration_network.pdf",
        ("Figure3", "E"): f3 / "Figure3E_regulon_size_partners.pdf",
        ("Figure3", "F"): f3 / "Figure3F_HC_TF_activity_skewness.pdf",
        ("Figure3", "G"): f3 / "Figure3G_representative_distributions.pdf",
        ("Figure4", "A"): atom / "Figure4A_compact_forest.pdf",
        ("Figure4", "B"): atom / "Figure4B_compact_forest.pdf",
        ("Figure4", "C"): atom / "Figure4C_compact_forest.pdf",
        ("Figure4", "D"): atom / "Figure4D_compact_forest.pdf",
        ("Figure4", "E"): f4 / "Figure4E_GSE41271_endpoint_overlap.pdf",
        ("Figure4", "F"): km_gse,
        ("Figure4", "G"): f4 / "Figure4G_TCGA_endpoint_overlap.pdf",
        ("Figure4", "H"): km_tcga,
        ("Figure5", "A"): f5 / "Figure5A_GDSC2_volcano.pdf",
        ("Figure5", "B"): f5 / "Figure5B_CTRPv2_volcano.pdf",
        ("Figure5", "C"): f5 / "Figure5C_PRISM_volcano.pdf",
        ("Figure5", "D"): atom / "Figure5D_true_upset.pdf",
        ("Figure5", "E"): f5 / "Figure5E_replicated_association_heatmap.pdf",
        ("Figure5", "F"): atom / "Figure5F_pathway_TF_bubble.pdf",
        ("Figure5", "G"): atom / "Figure5G_PDX_validation_funnel.pdf",
        ("Figure5", "H"): atom / "Figure5H_PDX_complete_screen.pdf",
        ("SupplementaryFigure", "1"): f1 / "SupplementaryFigure1_complete.pdf",
        ("SupplementaryFigure", "2"): f1 / "SupplementaryFigure2_complete.pdf",
        ("SupplementaryFigure", "3"): f2 / "SupplementaryFigure3_complete.pdf",
        ("SupplementaryFigure", "4"): supp4_compare,
        ("SupplementaryFigure", "5"): supp5_compare,
        ("SupplementaryFigure", "6"): supp6_compare,
        ("SupplementaryFigure", "7"): supp7_compare,
    }
    anchor_crops = {
        ("Figure1", "A"): (35, (0.02, 0.05, 0.40, 0.44)),
        ("Figure1", "B"): (35, (0.39, 0.05, 0.64, 0.44)),
        ("Figure1", "C"): (35, (0.64, 0.05, 0.98, 0.44)),
        ("Figure1", "D"): (35, (0.02, 0.43, 0.39, 0.95)),
        ("Figure1", "E"): (35, (0.38, 0.43, 0.98, 0.95)),
        ("Figure2", "A"): (36, (0.01, 0.04, 0.32, 0.39)),
        ("Figure2", "B"): (36, (0.01, 0.29, 0.24, 0.56)),
        ("Figure2", "C"): (36, (0.01, 0.48, 0.99, 0.98)),
        ("Figure2", "D"): (36, (0.22, 0.25, 0.46, 0.55)),
        ("Figure2", "E"): (36, (0.44, 0.04, 0.99, 0.51)),
        ("Figure3", "A"): (37, (0.01, 0.04, 0.44, 0.41)),
        ("Figure3", "B"): (37, (0.43, 0.04, 0.99, 0.42)),
        ("Figure3", "C"): (37, (0.01, 0.38, 0.47, 0.73)),
        ("Figure3", "D"): (37, (0.45, 0.38, 0.99, 0.72)),
        ("Figure3", "E"): (37, (0.48, 0.69, 0.99, 0.85)),
        ("Figure3", "F"): (37, (0.01, 0.70, 0.47, 0.97)),
        ("Figure3", "G"): (37, (0.48, 0.82, 0.99, 0.98)),
        ("Figure4", "A"): (38, (0.01, 0.04, 0.25, 0.35)),
        ("Figure4", "B"): (38, (0.25, 0.04, 0.50, 0.35)),
        ("Figure4", "C"): (38, (0.50, 0.04, 0.75, 0.35)),
        ("Figure4", "D"): (38, (0.75, 0.04, 0.99, 0.35)),
        ("Figure4", "E"): (38, (0.01, 0.31, 0.50, 0.50)),
        ("Figure4", "F"): (38, (0.01, 0.45, 0.50, 0.98)),
        ("Figure4", "G"): (38, (0.50, 0.31, 0.99, 0.50)),
        ("Figure4", "H"): (38, (0.50, 0.45, 0.99, 0.98)),
        ("Figure5", "A"): (39, (0.01, 0.04, 0.34, 0.33)),
        ("Figure5", "B"): (39, (0.33, 0.04, 0.67, 0.33)),
        ("Figure5", "C"): (39, (0.66, 0.04, 0.99, 0.33)),
        ("Figure5", "D"): (39, (0.01, 0.29, 0.35, 0.58)),
        ("Figure5", "E"): (39, (0.01, 0.52, 0.35, 0.72)),
        ("Figure5", "F"): (39, (0.33, 0.29, 0.99, 0.72)),
        ("Figure5", "G"): (39, (0.01, 0.68, 0.38, 0.98)),
        ("Figure5", "H"): (39, (0.36, 0.68, 0.99, 0.98)),
        ("SupplementaryFigure", "1"): (0, (0.00, 0.00, 1.00, 1.00)),
        ("SupplementaryFigure", "2"): (2, (0.00, 0.00, 1.00, 1.00)),
        ("SupplementaryFigure", "3"): (4, (0.00, 0.00, 1.00, 1.00)),
        ("SupplementaryFigure", "4"): (6, (0.00, 0.00, 1.00, 1.00)),
        ("SupplementaryFigure", "5"): (8, (0.00, 0.00, 1.00, 1.00)),
        ("SupplementaryFigure", "6"): (10, (0.00, 0.00, 1.00, 1.00)),
        ("SupplementaryFigure", "7"): (12, (0.00, 0.00, 1.00, 1.00)),
    }
    comparison_atlas(
        publication / "TNBC_LUAD_panel_comparison_atlas.pdf",
        Path(args.anchor_pdf),
        Path(args.anchor_supp_pdf),
        Path(args.gap_matrix),
        luad_sources,
        anchor_crops,
    )

    manifest_path = rebuild / "results" / "tables" / "publication_panel_source_manifest.tsv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["artifact", "panel", "source", "source_exists", "destination"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    output_index = rebuild / "results" / "tables" / "publication_pdf_index.tsv"
    with output_index.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["file_name", "pages", "bytes", "md5"])
        for path in sorted(publication.glob("*.pdf")):
            with fitz.open(path) as doc:
                pages = doc.page_count
            writer.writerow([path.name, pages, path.stat().st_size, md5(path)])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("fig1")
    p.add_argument("fig2")
    p.add_argument("fig3")
    p.add_argument("fig4")
    p.add_argument("fig5")
    p.add_argument("rebuild")
    p.add_argument("anchor_pdf")
    p.add_argument("anchor_supp_pdf")
    p.add_argument("gap_matrix")
    return p.parse_args()


if __name__ == "__main__":
    build(parse_args())
