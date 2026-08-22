#!/usr/bin/env python3
"""Assemble the frozen FDR05/HC45 LUAD panels without redrawing them.

All scientific panels are PDFs emitted by the TNBC official-code port (or by the
declared schematic renderer where the capsule did not provide plotting code).
This script only places those vector pages on journal-style sheets, adds panel
letters, and writes a source/hash manifest.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import io
import os
from pathlib import Path
from typing import Iterable, Sequence

from pypdf import PdfReader, PdfWriter, Transformation
from pypdf._page import PageObject
from reportlab.pdfgen import canvas


PAGE_WIDTH = 612.0
PAGE_HEIGHT = 792.0
MARGIN = 14.0
GUTTER = 4.0


def extended(path: Path) -> str:
    """Return a Windows long-path-safe absolute path string."""
    absolute = str(path.resolve())
    if os.name == "nt" and not absolute.startswith("\\\\?\\"):
        return "\\\\?\\" + absolute
    return absolute


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(extended(path), "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reader(path: Path) -> PdfReader:
    return PdfReader(extended(path), strict=False)


def normalized_box(box: Sequence[float]) -> tuple[float, float, float, float]:
    """Convert [x, top, width, height] fractions to PDF points."""
    x, top, width, height = (float(value) for value in box)
    content_width = PAGE_WIDTH - 2 * MARGIN
    content_height = PAGE_HEIGHT - 2 * MARGIN
    return (
        MARGIN + x * content_width,
        MARGIN + (1.0 - top - height) * content_height,
        width * content_width,
        height * content_height,
    )


def grid_boxes(
    box: tuple[float, float, float, float], rows: int, cols: int, count: int
) -> list[tuple[float, float, float, float]]:
    x, y, width, height = box
    cell_width = (width - GUTTER * (cols - 1)) / cols
    cell_height = (height - GUTTER * (rows - 1)) / rows
    result: list[tuple[float, float, float, float]] = []
    for index in range(count):
        row = index // cols
        col = index % cols
        result.append(
            (
                x + col * (cell_width + GUTTER),
                y + (rows - 1 - row) * (cell_height + GUTTER),
                cell_width,
                cell_height,
            )
        )
    return result


def place_pdf(
    target: PageObject,
    source_path: Path,
    destination: tuple[float, float, float, float],
) -> None:
    source = copy.copy(reader(source_path).pages[0])
    left = float(source.mediabox.left)
    bottom = float(source.mediabox.bottom)
    source_width = float(source.mediabox.width)
    source_height = float(source.mediabox.height)
    x, y, width, height = destination
    # Keep a small panel-letter band, while preserving the source aspect ratio.
    inner_x = x + 7.0
    inner_y = y + 1.5
    inner_width = max(1.0, width - 8.5)
    inner_height = max(1.0, height - 5.0)
    scale = min(inner_width / source_width, inner_height / source_height)
    tx = inner_x + (inner_width - source_width * scale) / 2.0
    ty = inner_y + (inner_height - source_height * scale) / 2.0
    transform = (
        Transformation()
        .translate(-left, -bottom)
        .scale(scale, scale)
        .translate(tx, ty)
    )
    target.merge_transformed_page(source, transform, over=True, expand=False)


def add_labels(
    target: PageObject,
    labels: Iterable[tuple[str, tuple[float, float, float, float]]],
) -> None:
    stream = io.BytesIO()
    overlay = canvas.Canvas(stream, pagesize=(PAGE_WIDTH, PAGE_HEIGHT), pageCompression=1)
    overlay.setFillColorRGB(0, 0, 0)
    overlay.setFont("Helvetica-Bold", 11)
    for label, (x, y, _width, height) in labels:
        overlay.drawString(x + 0.8, y + height - 10.5, label)
    overlay.showPage()
    overlay.save()
    stream.seek(0)
    target.merge_page(PdfReader(stream).pages[0], over=True, expand=False)


def assemble(
    output: Path,
    panels: Sequence[dict],
    manifest_rows: list[dict[str, str]],
) -> None:
    page = PageObject.create_blank_page(width=PAGE_WIDTH, height=PAGE_HEIGHT)
    labels: list[tuple[str, tuple[float, float, float, float]]] = []
    for panel in panels:
        outer = normalized_box(panel["box"])
        sources = [Path(item) for item in panel["sources"]]
        rows, cols = panel.get("grid", (1, 1))
        destinations = grid_boxes(outer, rows, cols, len(sources))
        for index, (source_path, destination) in enumerate(zip(sources, destinations), start=1):
            if not os.path.isfile(extended(source_path)):
                raise FileNotFoundError(source_path)
            place_pdf(page, source_path, destination)
            manifest_rows.append(
                {
                    "output": output.name,
                    "panel": panel["label"],
                    "source_index": str(index),
                    "source": str(source_path.resolve()),
                    "source_sha256": sha256(source_path),
                    "embedding": "pypdf_vector_merge",
                }
            )
        if panel.get("show_label", True):
            labels.append((panel["label"], outer))
    add_labels(page, labels)
    writer = PdfWriter()
    writer.add_page(page)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(extended(output), "wb") as handle:
        writer.write(handle)


def assemble_single(
    output: Path,
    source: Path,
    label: str,
    manifest_rows: list[dict[str, str]],
) -> None:
    original = reader(source).pages[0]
    width = float(original.mediabox.width)
    height = float(original.mediabox.height)
    page = PageObject.create_blank_page(width=width, height=height)
    copied = copy.copy(original)
    page.merge_page(copied, over=True, expand=False)
    stream = io.BytesIO()
    overlay = canvas.Canvas(stream, pagesize=(width, height), pageCompression=1)
    overlay.setFont("Helvetica-Bold", 11)
    overlay.drawString(7.0, height - 14.0, label)
    overlay.showPage()
    overlay.save()
    stream.seek(0)
    page.merge_page(PdfReader(stream).pages[0], over=True, expand=False)
    writer = PdfWriter()
    writer.add_page(page)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(extended(output), "wb") as handle:
        writer.write(handle)
    manifest_rows.append(
        {
            "output": output.name,
            "panel": label,
            "source_index": "1",
            "source": str(source.resolve()),
            "source_sha256": sha256(source),
            "embedding": "pypdf_vector_merge_original_size",
        }
    )


def main() -> int:
    project = Path(__file__).resolve().parents[4]
    run = project / "pilots/tnbc-chromatin-tf-nc-2026/execution/luad-fdr05-hc45-primary-v1"
    out = run / "output/pdf"
    fig1 = run / "official-fig12-port-compact/run/adapter_final_corrected/results/visuals/Figure1"
    fig2 = run / "official-fig12-port-compact/run/adapter_final_corrected/results/visuals/Figure2"
    fig3 = run / "results/figure3b-official/results/visuals/Figure3"
    fig4 = run / "results/figure4-official-logrank-v5/results/atomic"
    schematic1 = run / "schematics/figure1-schematics-updated/results/figures"
    schematic = run / "schematics/schematics-updated/anchor_style/atomic"
    manifest_rows: list[dict[str, str]] = []

    assemble(
        out / "Figure1_primary_FDR05.pdf",
        [
            {"label": "A", "sources": [schematic1 / "Figure1A_workflow_final.pdf"], "box": [0.00, 0.00, 0.40, 0.34], "show_label": False},
            {"label": "B", "sources": [fig1 / "Figure1B_log2FC_Viper_TNBC.pdf"], "box": [0.40, 0.00, 0.22, 0.34]},
            {"label": "C", "sources": [fig1 / "Figure1C_TCGA_TRActivity_TNBC150_NonTNBC155_test.pdf"], "box": [0.62, 0.00, 0.38, 0.34]},
            {"label": "D", "sources": [fig1 / "Figure1D_PDX_TRActivity_TNBC150_NonTNBC154.pdf"], "box": [0.00, 0.34, 0.50, 0.66]},
            {"label": "E", "sources": [fig1 / "Figure1E_CellLine_TRActivity_TNBC150_NonTNBC155.pdf"], "box": [0.50, 0.34, 0.50, 0.66]},
        ],
        manifest_rows,
    )

    assemble(
        out / "Figure2_primary_HC45.pdf",
        [
            {"label": "A", "sources": [schematic / "Figure2A_anchor_schematic.pdf"], "box": [0.00, 0.00, 0.22, 0.25], "show_label": False},
            {"label": "B", "sources": [fig2 / "Figure2B_TF_PromoterAccessibility_Intersection_Sample_Groups.pdf"], "box": [0.00, 0.25, 0.22, 0.25]},
            {"label": "D", "sources": [fig2 / "Figure2D_JASPAR_CISBP_TR_MotifEnrichment_Intersection_Sample_Groups.pdf"], "box": [0.22, 0.00, 0.18, 0.50]},
            {"label": "E", "sources": [fig2 / "Figure2E_Complete_HC-TRs.pdf"], "box": [0.40, 0.00, 0.60, 0.50]},
            {"label": "C", "sources": [fig2 / "Figure2C_JASPAR_CISBP_Combined_ChromatinAccessibility_MotifEnrichment_Heatmap.pdf"], "box": [0.00, 0.50, 1.00, 0.50]},
        ],
        manifest_rows,
    )

    assemble_single(
        out / "Figure3B_shared_ge2.pdf",
        fig3 / "Figure3B_Target_Network_Shared_TFs.pdf",
        "B",
        manifest_rows,
    )

    assemble(
        out / "Figure4_primary_unadjusted_logrank.pdf",
        [
            {"label": "A", "sources": [fig4 / "Figure4A.pdf"], "box": [0.00, 0.00, 0.25, 0.36]},
            {"label": "B", "sources": [fig4 / "Figure4B.pdf"], "box": [0.25, 0.00, 0.25, 0.36]},
            {"label": "C", "sources": [fig4 / "Figure4C.pdf"], "box": [0.50, 0.00, 0.25, 0.36]},
            {"label": "D", "sources": [fig4 / "Figure4D.pdf"], "box": [0.75, 0.00, 0.25, 0.36]},
            {"label": "E", "sources": [fig4 / "Figure4E_protective.pdf", fig4 / "Figure4E_risk.pdf"], "box": [0.00, 0.36, 0.50, 0.18], "grid": [1, 2]},
            {"label": "G", "sources": [fig4 / "Figure4G_protective.pdf", fig4 / "Figure4G_risk.pdf"], "box": [0.50, 0.36, 0.50, 0.18], "grid": [1, 2]},
            {"label": "F", "sources": [fig4 / "Figure4F_ZNF540_OS.pdf", fig4 / "Figure4F_ZNF540_RFS.pdf", fig4 / "Figure4F_PHC2_OS.pdf", fig4 / "Figure4F_PHC2_RFS.pdf"], "box": [0.00, 0.54, 0.50, 0.46], "grid": [2, 2]},
            {"label": "H", "sources": [fig4 / "Figure4H_CRY2_OS.pdf", fig4 / "Figure4H_CRY2_DFS.pdf"], "box": [0.50, 0.54, 0.50, 0.46], "grid": [1, 2]},
        ],
        manifest_rows,
    )

    assemble(
        out / "SupplementaryFigure1_primary_FDR05.pdf",
        [
            {"label": "A", "sources": [schematic / "SupplementaryFigure1A_TCGA_inclusion.pdf"], "box": [0.00, 0.00, 1.00, 0.32], "show_label": False},
            {"label": "B", "sources": [schematic / "SupplementaryFigure1B_GSE81089_inclusion.pdf"], "box": [0.00, 0.32, 1.00, 0.34], "show_label": False},
            {"label": "C", "sources": [schematic / "SupplementaryFigure1C_GSE41271_inclusion.pdf"], "box": [0.00, 0.66, 1.00, 0.34], "show_label": False},
        ],
        manifest_rows,
    )

    assemble(
        out / "SupplementaryFigure2_primary_FDR05.pdf",
        [
            {"label": "A", "sources": [fig1 / "SupplementaryFigure2A_METABRIC_log2FC_Viper.pdf"], "box": [0.00, 0.00, 0.34, 0.31]},
            {"label": "B", "sources": [fig1 / "SupplementaryFigure2B_VennDiagram_NonTNBC_TCGA_METABRIC.pdf", fig1 / "SupplementaryFigure2B_VennDiagram_TNBC_TCGA_METABRIC.pdf"], "box": [0.34, 0.00, 0.66, 0.31], "grid": [1, 2]},
            {"label": "C", "sources": [fig1 / "SupplementaryFigure2C_TCGA_TRActivity_TNBC150_NonTNBC155_PearsonCor.pdf"], "box": [0.00, 0.31, 0.34, 0.31]},
            {"label": "D", "sources": [fig1 / "SupplementaryFigure2D_PDX_TRActivity_TNBC150_NonTNBC155_PearsonCor.pdf"], "box": [0.34, 0.31, 0.33, 0.31]},
            {"label": "E", "sources": [fig1 / "SupplementaryFigure2E_Cell_Lines_TRActivity_TNBC150_NonTNBC155_PearsonCor.pdf"], "box": [0.67, 0.31, 0.33, 0.31]},
            {"label": "F", "sources": [fig1 / "SupplementaryFigure2F_TCGA_Basal156TRs.pdf", fig1 / "SupplementaryFigure2F_PDX_Basal156TRs.pdf", fig1 / "SupplementaryFigure2F_CellLines_Basal156TRs.pdf", fig1 / "SupplementaryFigure2F_METABRIC_Basal122TRs.pdf"], "box": [0.00, 0.62, 1.00, 0.38], "grid": [1, 4]},
        ],
        manifest_rows,
    )

    assemble(
        out / "SupplementaryFigure4_HC45.pdf",
        [
            {"label": "A", "sources": [schematic / "SupplementaryFigure4A_anchor_gate_schematic.pdf"], "box": [0.00, 0.00, 1.00, 0.28], "show_label": False},
            {"label": "B", "sources": [fig2 / "SupplementaryFigure4B_NES_Average_Across_CellLines_PDX_TCGA.pdf"], "box": [0.00, 0.28, 0.50, 0.72]},
            {"label": "C", "sources": [fig2 / "JASPAR_vs_CISBP_HC-TRs_VennDiagram.pdf"], "box": [0.50, 0.28, 0.50, 0.72]},
        ],
        manifest_rows,
    )

    assemble(
        out / "SupplementaryFigure6_logrank_sensitivity.pdf",
        [
            {"label": "A", "sources": [fig4 / "SupplementaryFigure6A.pdf"], "box": [0.00, 0.00, 0.25, 0.50]},
            {"label": "B", "sources": [fig4 / "SupplementaryFigure6B.pdf"], "box": [0.25, 0.00, 0.25, 0.50]},
            {"label": "C", "sources": [fig4 / "SupplementaryFigure6C.pdf"], "box": [0.50, 0.00, 0.25, 0.50]},
            {"label": "D", "sources": [fig4 / "SupplementaryFigure6D.pdf"], "box": [0.75, 0.00, 0.25, 0.50]},
            {"label": "E", "sources": [fig4 / "SupplementaryFigure6E.pdf"], "box": [0.00, 0.50, 0.25, 0.50]},
            {"label": "F", "sources": [fig4 / "SupplementaryFigure6F.pdf"], "box": [0.25, 0.50, 0.25, 0.50]},
            {"label": "G", "sources": [fig4 / "SupplementaryFigure6G.pdf"], "box": [0.50, 0.50, 0.25, 0.50]},
            {"label": "H", "sources": [fig4 / "SupplementaryFigure6H.pdf"], "box": [0.75, 0.50, 0.25, 0.50]},
        ],
        manifest_rows,
    )

    manifest = out.parent / "primary_delivery_embedding_manifest.tsv"
    with open(extended(manifest), "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["output", "panel", "source_index", "source", "source_sha256", "embedding"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    checksums = out.parent / "primary_delivery_sha256.tsv"
    with open(extended(checksums), "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["file", "sha256", "pages", "bytes"])
        for pdf in sorted(out.glob("*.pdf")):
            writer.writerow([pdf.name, sha256(pdf), len(reader(pdf).pages), pdf.stat().st_size])

    print(f"Wrote {len(list(out.glob('*.pdf')))} PDFs to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
