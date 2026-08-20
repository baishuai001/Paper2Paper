#!/usr/bin/env python3
"""Fail-closed verification and page rendering for publication rebuild outputs."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import fitz
import numpy as np
from PIL import Image, ImageDraw, ImageFont


EXPECTED = {
    "Figure1_revised.pdf": 1,
    "Figure2_revised.pdf": 1,
    "Figure3_revised.pdf": 1,
    "Figure4_revised.pdf": 1,
    "Figure5_revised.pdf": 1,
    "TNBC_LUAD_panel_comparison_atlas.pdf": 40,
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def render_pdf(path: Path, qa_root: Path, metrics: list[dict[str, object]]) -> list[Path]:
    out_dir = qa_root / "pages" / path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    rendered = []
    with fitz.open(path) as doc:
        for pno, page in enumerate(doc):
            pix = page.get_pixmap(matrix=fitz.Matrix(1.25, 1.25), alpha=False)
            png_path = out_dir / f"page-{pno + 1:03d}.png"
            pix.save(png_path)
            image = Image.open(png_path).convert("RGB")
            arr = np.asarray(image)
            gray = arr.mean(axis=2)
            white_fraction = float((gray > 248).mean())
            ink = gray < 245
            if ink.any():
                ys, xs = np.where(ink)
                bbox_fraction = float(((xs.max() - xs.min() + 1) * (ys.max() - ys.min() + 1)) / ink.size)
            else:
                bbox_fraction = 0.0
            text = page.get_text("text")
            metrics.append({
                "pdf": path.name,
                "page": pno + 1,
                "width_pt": round(page.rect.width, 2),
                "height_pt": round(page.rect.height, 2),
                "white_fraction": round(white_fraction, 5),
                "ink_bbox_fraction": round(bbox_fraction, 5),
                "text_chars": len(text),
                "missing_source_marker": "MISSING SOURCE" in text,
                "status": "PASS" if white_fraction < 0.985 and bbox_fraction > 0.08 and "MISSING SOURCE" not in text else "FAIL",
            })
            rendered.append(png_path)
    return rendered


def contact_sheets(images: list[Path], output_prefix: Path, cols: int = 4) -> list[Path]:
    if not images:
        return []
    thumb_w, thumb_h = 360, 280
    per_sheet = 12
    outputs = []
    font = ImageFont.load_default()
    for block_index in range(math.ceil(len(images) / per_sheet)):
        block = images[block_index * per_sheet:(block_index + 1) * per_sheet]
        rows = math.ceil(len(block) / cols)
        canvas = Image.new("RGB", (cols * thumb_w, rows * (thumb_h + 24)), "white")
        draw = ImageDraw.Draw(canvas)
        for i, image_path in enumerate(block):
            image = Image.open(image_path).convert("RGB")
            # Pillow < 9.1 exposes LANCZOS directly on Image rather than through
            # Image.Resampling.  Keep verification portable across the cloud
            # runtime and the bundled local document runtime.
            resampling = getattr(Image, "Resampling", Image)
            image.thumbnail((thumb_w - 12, thumb_h - 12), resampling.LANCZOS)
            x = (i % cols) * thumb_w + (thumb_w - image.width) // 2
            y = (i // cols) * (thumb_h + 24) + 6
            canvas.paste(image, (x, y))
            draw.text(((i % cols) * thumb_w + 8, (i // cols) * (thumb_h + 24) + thumb_h + 4),
                      image_path.parent.name + "/" + image_path.stem, fill="black", font=font)
        out = output_prefix.with_name(f"{output_prefix.name}-{block_index + 1:02d}.png")
        out.parent.mkdir(parents=True, exist_ok=True)
        canvas.save(out)
        outputs.append(out)
    return outputs


def assertions(root: Path, publication: Path) -> list[dict[str, object]]:
    rows = []

    def check(name: str, observed, expected, ok: bool):
        rows.append({"check": name, "observed": observed, "expected": expected, "status": "PASS" if ok else "FAIL"})

    pdfs = sorted(publication.glob("*.pdf"))
    check("final_pdf_count", len(pdfs), 7, len(pdfs) == 7)
    for name, pages in EXPECTED.items():
        path = publication / name
        observed = 0
        if path.exists():
            with fitz.open(path) as doc:
                observed = doc.page_count
        check(f"pages:{name}", observed, pages, observed == pages)
    supp_path = publication / "LUAD_Supplementary_Figures_1-9.pdf"
    supp_pages = 0
    if supp_path.exists():
        with fitz.open(supp_path) as doc:
            supp_pages = doc.page_count
    check("supplement_page_count", supp_pages, ">=10", supp_pages >= 10)

    receipt = {r["key"]: r["value"] for r in read_tsv(root / "audit" / "publication_rebuild_receipt.tsv")}
    check("HC_TFs", int(receipt["frozen_HC_TFs"]), 31, int(receipt["frozen_HC_TFs"]) == 31)
    check("strict_motif", receipt["strict_three_system_motif_TFs"], "FOXA3;NFATC4;XBP1",
          receipt["strict_three_system_motif_TFs"] == "FOXA3;NFATC4;XBP1")
    check("sensitivity_motif", receipt["sensitivity_q0.05_prev0.5_TFs"],
          "ETV1;FOXA3;NFATC4;XBP1;ZNF75D",
          receipt["sensitivity_q0.05_prev0.5_TFs"] == "ETV1;FOXA3;NFATC4;XBP1;ZNF75D")
    check("replicated_cellline_pairs", int(receipt["replicated_cellline_pairs"]), 5,
          int(receipt["replicated_cellline_pairs"]) == 5)
    check("validated_PDX_pairs", int(receipt["validated_PDX_pairs"]), 0,
          int(receipt["validated_PDX_pairs"]) == 0)
    check("target_network_nodes", int(receipt["target_network_nodes"]), ">0", int(receipt["target_network_nodes"]) > 0)
    check("target_network_edges", int(receipt["target_network_edges"]), ">0", int(receipt["target_network_edges"]) > 0)

    mki = read_tsv(root / "results" / "tables" / "SupplementaryFigure4D_MKI67_correlations.tsv")
    check("MKI67_rows", len(mki), 155, len(mki) == 155)
    est = read_tsv(root / "results" / "tables" / "SupplementaryFigure5B_ESTIMATE_correlations.tsv")
    check("ESTIMATE_rows", len(est), 186, len(est) == 186)

    source_manifest = read_tsv(root / "results" / "tables" / "publication_panel_source_manifest.tsv")
    missing = [r for r in source_manifest if r["source_exists"] != "True"]
    check("main_panel_sources", len(missing), 0, len(missing) == 0)

    # Figure 4F/H were originally exported as blank PDFs by the upstream R
    # graphics device even though the paired high-resolution PNGs were valid.
    # Guard the lower KM region explicitly so a nonblank header cannot mask a
    # missing survival panel in whole-page white-fraction checks.
    figure4_render = root / "qa" / "pages" / "Figure4_revised" / "page-001.png"
    lower_ink = 0.0
    if figure4_render.exists():
        image = np.asarray(Image.open(figure4_render).convert("RGB"))
        gray = image.mean(axis=2)
        y0, y1 = int(gray.shape[0] * 0.58), int(gray.shape[0] * 0.93)
        lower_ink = float((gray[y0:y1, :] < 245).mean())
    check("Figure4_FH_lower_region_ink", round(lower_ink, 5), ">=0.012", lower_ink >= 0.012)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root")
    args = parser.parse_args()
    root = Path(args.root)
    publication = root / "results" / "publication"
    qa_root = root / "qa"
    qa_root.mkdir(parents=True, exist_ok=True)

    page_metrics: list[dict[str, object]] = []
    all_images = []
    for pdf in sorted(publication.glob("*.pdf")):
        all_images.extend(render_pdf(pdf, qa_root, page_metrics))
    contacts = contact_sheets(all_images, qa_root / "contact_sheets" / "all-final-pages")
    checks = assertions(root, publication)

    with (qa_root / "page_render_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(page_metrics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(page_metrics)
    with (qa_root / "verification_checks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(checks[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(checks)

    failed = [r for r in page_metrics if r["status"] != "PASS"] + [r for r in checks if r["status"] != "PASS"]
    print(f"Rendered {len(page_metrics)} pages; contact sheets={len(contacts)}; failed checks={len(failed)}")
    if failed:
        for row in failed:
            print(row)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
