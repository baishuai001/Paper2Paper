#!/usr/bin/env python3
"""Render and sanity-check the frozen primary-delivery PDF set."""

from __future__ import annotations

import csv
import os
from pathlib import Path

import pypdfium2 as pdfium
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageStat
from pypdf import PdfReader


def extended(path: Path) -> str:
    absolute = str(path.resolve())
    if os.name == "nt" and not absolute.startswith("\\\\?\\"):
        return "\\\\?\\" + absolute
    return absolute


def main() -> int:
    project = Path(__file__).resolve().parents[4]
    root = project / "pilots/tnbc-chromatin-tf-nc-2026/execution/luad-fdr05-hc45-primary-v1/output"
    pdf_dir = root / "pdf"
    render_dir = root / "qa_rendered"
    render_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    rendered_images: list[tuple[str, Path]] = []
    failed = False

    for path in sorted(pdf_dir.glob("*.pdf")):
        header_ok = path.read_bytes()[:5] == b"%PDF-"
        parsed_pages = len(PdfReader(extended(path), strict=False).pages)
        document = pdfium.PdfDocument(extended(path))
        rendered_pages = len(document)
        if rendered_pages != 1:
            failed = True
        page = document[0]
        bitmap = page.render(scale=2.0)
        image = bitmap.to_pil().convert("RGB")
        output_png = render_dir / f"{path.stem}.png"
        image.save(extended(output_png), format="PNG", optimize=True)
        rendered_images.append((path.name, output_png))
        gray = image.convert("L")
        difference = ImageChops.difference(gray, gray.point(lambda _value: 255))
        bbox = difference.getbbox()
        mean = ImageStat.Stat(gray).mean[0]
        nonwhite_bbox = "" if bbox is None else ",".join(str(value) for value in bbox)
        status = "PASS" if header_ok and parsed_pages == 1 and rendered_pages == 1 and bbox else "FAIL"
        if status != "PASS":
            failed = True
        rows.append(
            {
                "pdf": path.name,
                "status": status,
                "bytes": str(path.stat().st_size),
                "parsed_pages": str(parsed_pages),
                "rendered_pages": str(rendered_pages),
                "render_width_px": str(image.width),
                "render_height_px": str(image.height),
                "mean_gray": f"{mean:.3f}",
                "nonwhite_bbox_px": nonwhite_bbox,
                "rendered_png": str(output_png.resolve()),
            }
        )
        document.close()

    thumb_width, thumb_height = 560, 740
    label_height = 30
    contact = Image.new("RGB", (thumb_width * 2, (thumb_height + label_height) * 4), "white")
    draw = ImageDraw.Draw(contact)
    font = ImageFont.load_default()
    for index, (name, png) in enumerate(rendered_images):
        row, col = divmod(index, 2)
        with Image.open(extended(png)) as source_image:
            thumbnail = source_image.convert("RGB")
            thumbnail.thumbnail((thumb_width - 12, thumb_height - 12))
            x = col * thumb_width + (thumb_width - thumbnail.width) // 2
            y = row * (thumb_height + label_height) + label_height + (thumb_height - thumbnail.height) // 2
            contact.paste(thumbnail, (x, y))
        draw.text((col * thumb_width + 8, row * (thumb_height + label_height) + 8), name, fill="black", font=font)
    contact_path = root / "primary_delivery_contact_sheet.png"
    contact.save(extended(contact_path), format="PNG", optimize=True)

    receipt = root / "primary_delivery_pdf_qa.tsv"
    with receipt.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Verified {len(rows)} PDFs; failures={sum(row['status'] != 'PASS' for row in rows)}")
    print(receipt)
    print(contact_path)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
