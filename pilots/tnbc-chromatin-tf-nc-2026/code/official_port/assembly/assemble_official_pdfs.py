#!/usr/bin/env python3
"""Assemble strict official-code atomic PDFs without redrawing any panel.

The LUAD plot atoms are imported as PDF content streams with pypdf.  They are
never rasterized, recolored, re-themed, cropped, or otherwise reconstructed.
Only page placement, panel letters, provenance notes, and declared neutral
boundary placeholders are authored here.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import io
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

try:
    from PIL import Image
    from pypdf import PdfReader, PdfWriter, Transformation
    from pypdf._page import PageObject
    from reportlab.lib.colors import HexColor
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfgen import canvas
except ImportError as exc:  # pragma: no cover - exercised on deployment hosts
    raise SystemExit(
        "Missing PDF dependency. Install pypdf, reportlab, and Pillow in the "
        f"active Python environment. Original error: {exc}"
    )


SCRIPT_DIR = Path(__file__).resolve().parent
EXPECTED_FINALS = (
    "Figure1.pdf",
    "Figure2.pdf",
    "Figure3.pdf",
    "Figure4.pdf",
    "Figure5.pdf",
    "Supplementary_Figures.pdf",
    "TNBC_vs_LUAD_comparison_atlas.pdf",
)


@dataclass(frozen=True)
class SourceSlice:
    path: Path
    page_index: int = 0
    crop_top: tuple[float, float, float, float] = (0.0, 0.0, 1.0, 1.0)
    role: str = "LUAD_OFFICIAL_ATOM"
    atom_id: str = ""


@dataclass
class ResolvedAtom:
    atom_id: str
    status: str
    sources: list[SourceSlice] = field(default_factory=list)
    note: str = ""
    boundary: str = ""
    required: bool = False
    root_key: str = ""
    patterns: list[str] = field(default_factory=list)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate_pdf(path: Path) -> tuple[int, str]:
    if not path.is_file():
        raise ValueError(f"PDF does not exist: {path}")
    if path.stat().st_size < 1000:
        raise ValueError(f"PDF is empty or too small: {path}")
    with path.open("rb") as handle:
        if handle.read(5) != b"%PDF-":
            raise ValueError(f"Invalid PDF header: {path}")
    reader = PdfReader(str(path), strict=False)
    if not reader.pages:
        raise ValueError(f"PDF has no pages: {path}")
    return len(reader.pages), sha256_file(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Vector-assemble the seven strict official-code LUAD deliverables."
    )
    parser.add_argument("--fig12-root", required=True, type=Path)
    parser.add_argument("--fig3-root", required=True, type=Path)
    parser.add_argument("--fig45-root", required=True, type=Path)
    parser.add_argument(
        "--supp3-root",
        type=Path,
        help="Root containing official Supplementary Figure 3 A-L atoms; may be supplied later.",
    )
    parser.add_argument(
        "--supp3-verification-receipt",
        type=Path,
        help=(
            "Explicit Supplementary Figure 3 verification_receipt.json. By default, "
            "exactly one receipt is discovered recursively below --supp3-root."
        ),
    )
    parser.add_argument("--tnbc-main-pdf", type=Path)
    parser.add_argument("--tnbc-supplement-pdf", type=Path)
    parser.add_argument(
        "--tnbc-reference-dir",
        type=Path,
        help=(
            "Directory containing s41467-026-76385-8_reference_36.png through _40.png "
            "and 41467_2026_76385_MOESM1_ESM_01.png etc. Used when source PDFs are absent."
        ),
    )
    parser.add_argument("--spec", type=Path, default=SCRIPT_DIR / "assembly_spec.json")
    parser.add_argument("--tnbc-panel-map", type=Path, default=SCRIPT_DIR / "tnbc_panel_map.tsv")
    parser.add_argument("--explanations", type=Path, default=SCRIPT_DIR / "panel_explanations.tsv")
    parser.add_argument(
        "--source-override",
        type=Path,
        help=(
            "Optional TSV with atom_id and path columns. Paths may be absolute or glob patterns; "
            "page_1based and normalized crop columns are optional."
        ),
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--audit-dir",
        type=Path,
        help="Defaults to OUTPUT_DIR/../audit.",
    )
    parser.add_argument(
        "--font-path",
        type=Path,
        help="Optional TTF font for non-ASCII atlas explanation text.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Resolve and validate the frozen input contract, write audit manifests, "
            "and create no final PDFs. Exits non-zero when a required input is absent."
        ),
    )
    parser.add_argument(
        "--require-publication-ready",
        action="store_true",
        help=(
            "Activate the release gate. No final PDFs are generated and the command exits "
            "non-zero unless all parsed upstream release receipts are publication-ready."
        ),
    )
    return parser.parse_args()


class FontBook:
    def __init__(self, font_path: Path | None):
        self.regular = "Helvetica"
        self.bold = "Helvetica-Bold"
        if font_path:
            path = font_path.resolve()
            if not path.is_file():
                raise ValueError(f"Font file not found: {path}")
            pdfmetrics.registerFont(TTFont("AssemblyFont", str(path)))
            font = pdfmetrics.getFont("AssemblyFont")
            coverage_probe = "A1中审阅"
            missing = [
                char
                for char in coverage_probe
                if font.face.charToGlyph.get(ord(char)) is None
            ]
            if missing:
                raise ValueError(
                    "Assembly font lacks required Latin/CJK glyph coverage: "
                    + ",".join(f"U+{ord(char):04X}" for char in missing)
                    + f" ({path})"
                )
            self.regular = "AssemblyFont"
            self.bold = "AssemblyFont"

    def ensure(self, text: str) -> None:
        if self.regular == "Helvetica" and any(ord(ch) > 127 for ch in text):
            raise ValueError(
                "Non-ASCII text requires --font-path pointing to a Unicode TTF font. "
                "The assembly stops instead of emitting missing-glyph squares."
            )


class AtomResolver:
    def __init__(
        self,
        atoms: dict[str, dict[str, Any]],
        roots: dict[str, Path | None],
        overrides: dict[str, list[SourceSlice]],
    ):
        self.atoms = atoms
        self.roots = roots
        self.overrides = overrides
        self.resolved: dict[str, ResolvedAtom] = {}
        self.resolution_rows: list[dict[str, Any]] = []
        self.missing_rows: list[dict[str, Any]] = []
        self.boundary_rows: list[dict[str, Any]] = []
        self._receipt_cache: dict[Path, list[dict[str, str]]] = {}

    @staticmethod
    def fixed_pattern_scope(pattern: str) -> str:
        """Return the non-glob directory prefix used as a provenance scope."""
        normalized = pattern.replace("\\", "/")
        positions = [normalized.find(char) for char in "*?[" if char in normalized]
        cut = min(positions) if positions else len(normalized)
        prefix = normalized[:cut]
        if prefix and not prefix.endswith("/"):
            prefix = prefix.rsplit("/", 1)[0] + "/" if "/" in prefix else ""
        return prefix.lstrip("./")

    def enforce_source_scope(
        self,
        atom_id: str,
        candidates: list[SourceSlice],
        root: Path | None,
        patterns: list[str],
        required: bool,
        root_key: str,
    ) -> list[SourceSlice]:
        """Confine every atom, including overrides, to its declared official root."""
        if not root:
            return candidates
        official_root = root.resolve()
        scopes = [scope for scope in (self.fixed_pattern_scope(p) for p in patterns) if scope]
        accepted: list[SourceSlice] = []
        for source in candidates:
            resolved = source.path.resolve()
            try:
                relative = resolved.relative_to(official_root).as_posix()
            except ValueError:
                relative = ""
            in_declared_root = bool(relative)
            in_pattern_scope = not scopes or any(relative.startswith(scope) for scope in scopes)
            if in_declared_root and in_pattern_scope:
                accepted.append(source)
                continue
            self.missing_rows.append(
                {
                    "atom_id": atom_id,
                    "status": "SOURCE_OUTSIDE_OFFICIAL_SCOPE",
                    "root_key": root_key,
                    "patterns": ";".join(patterns),
                    "detail": str(resolved),
                    "is_required": "YES" if required else "NO",
                }
            )
        return accepted

    def receipts(self, root: Path) -> list[dict[str, str]]:
        root = root.resolve()
        if root in self._receipt_cache:
            return self._receipt_cache[root]
        rows: list[dict[str, str]] = []
        for path in sorted(root.rglob("render_manifest.tsv")):
            try:
                for row in read_tsv(path):
                    row = dict(row)
                    row["_manifest"] = str(path)
                    rows.append(row)
            except Exception:
                continue
        self._receipt_cache[root] = rows
        return rows

    def supported_zero(self, root: Path, panel: str) -> tuple[bool, str]:
        matches = [row for row in self.receipts(root) if row.get("panel") == panel]
        zero = [row for row in matches if row.get("status") == "SUPPORTED_ZERO"]
        if not zero:
            return False, ""
        note = "; ".join(filter(None, (row.get("note", "") for row in zero)))
        return True, note

    def resolve_all(self) -> dict[str, ResolvedAtom]:
        for atom_id in self.atoms:
            self.resolve(atom_id)
        return self.resolved

    def resolve(self, atom_id: str) -> ResolvedAtom:
        if atom_id in self.resolved:
            return self.resolved[atom_id]
        spec = self.atoms[atom_id]
        declared = spec.get("status")
        required = bool(spec.get("required", False))
        boundary = str(spec.get("boundary", ""))
        note = str(spec.get("note", ""))
        root_key = str(spec.get("root", ""))
        patterns = [str(value) for value in spec.get("patterns", [])]

        if declared in {"MANUAL_REQUIRED", "UNSUPPORTED"}:
            result = ResolvedAtom(
                atom_id=atom_id,
                status=declared,
                note=note,
                boundary=boundary,
                required=False,
            )
            self.boundary_rows.append(
                {"atom_id": atom_id, "status": declared, "note": note, "is_missing": "NO"}
            )
            self.resolved[atom_id] = result
            return result

        root = self.roots.get(root_key)
        if atom_id in self.overrides:
            candidates = self.overrides[atom_id]
        else:
            candidates = []
            if root and root.is_dir():
                for pattern in patterns:
                    matches = sorted(path.resolve() for path in root.glob(pattern) if path.is_file())
                    if matches:
                        candidates = [SourceSlice(path=path, atom_id=atom_id) for path in matches]
                        break

        candidates = self.enforce_source_scope(
            atom_id, candidates, root, patterns, required, root_key
        )

        match_mode = str(spec.get("match", "one"))
        if match_mode == "one" and len(candidates) > 1:
            self.missing_rows.append(
                {
                    "atom_id": atom_id,
                    "status": "AMBIGUOUS_INPUT",
                    "root_key": root_key,
                    "patterns": ";".join(patterns),
                    "detail": ";".join(str(item.path) for item in candidates),
                    "is_required": "YES" if required else "NO",
                }
            )
            candidates = []

        valid: list[SourceSlice] = []
        for source in candidates:
            try:
                pages, digest = validate_pdf(source.path)
                if source.page_index < 0 or source.page_index >= pages:
                    raise ValueError(
                        f"Requested page {source.page_index + 1}, but source has {pages} pages"
                    )
                valid.append(source)
                self.resolution_rows.append(
                    {
                        "atom_id": atom_id,
                        "status": "RESOLVED_OFFICIAL_PDF",
                        "root_key": root_key,
                        "path": str(source.path),
                        "page_1based": source.page_index + 1,
                        "page_count": pages,
                        "sha256": digest,
                        "crop_top": ",".join(f"{v:.6g}" for v in source.crop_top),
                        "boundary": boundary,
                    }
                )
            except Exception as exc:
                self.missing_rows.append(
                    {
                        "atom_id": atom_id,
                        "status": "INVALID_INPUT",
                        "root_key": root_key,
                        "patterns": ";".join(patterns),
                        "detail": f"{source.path}: {exc}",
                        "is_required": "YES" if required else "NO",
                    }
                )

        if valid:
            result = ResolvedAtom(
                atom_id=atom_id,
                status="RESOLVED_OFFICIAL_PDF",
                sources=valid,
                note=note,
                boundary=boundary,
                required=required,
                root_key=root_key,
                patterns=patterns,
            )
        else:
            is_zero = False
            zero_note = ""
            if bool(spec.get("zero_allowed", False)) and root and root.is_dir():
                receipt_panel = str(spec.get("receipt_panel", atom_id))
                is_zero, zero_note = self.supported_zero(root, receipt_panel)
            if is_zero:
                result = ResolvedAtom(
                    atom_id=atom_id,
                    status="SUPPORTED_ZERO",
                    note=zero_note or "No eligible rows under the official constructor input rule.",
                    boundary=boundary,
                    required=False,
                    root_key=root_key,
                    patterns=patterns,
                )
                self.boundary_rows.append(
                    {
                        "atom_id": atom_id,
                        "status": "SUPPORTED_ZERO",
                        "note": result.note,
                        "is_missing": "NO",
                    }
                )
            else:
                status = "MISSING_REQUIRED" if required else "MISSING_OPTIONAL"
                result = ResolvedAtom(
                    atom_id=atom_id,
                    status=status,
                    note="No source matched the frozen patterns.",
                    boundary=boundary,
                    required=required,
                    root_key=root_key,
                    patterns=patterns,
                )
                self.missing_rows.append(
                    {
                        "atom_id": atom_id,
                        "status": status,
                        "root_key": root_key,
                        "patterns": ";".join(patterns),
                        "detail": (
                            f"Root unavailable: {root}" if not root or not root.is_dir()
                            else "No source matched the frozen patterns."
                        ),
                        "is_required": "YES" if required else "NO",
                    }
                )
        self.resolved[atom_id] = result
        return result


class PDFAssembler:
    def __init__(
        self,
        spec: dict[str, Any],
        resolver: AtomResolver,
        font_book: FontBook,
        output_dir: Path,
        audit_dir: Path,
        tnbc_map: dict[str, dict[str, str]],
        explanations: dict[str, dict[str, str]],
        tnbc_sources: dict[str, Path | None],
        reference_dir: Path | None,
        release_assessment: dict[str, Any],
        source_identity: dict[str, Any],
    ):
        self.spec = spec
        self.resolver = resolver
        self.fonts = font_book
        self.output_dir = output_dir
        self.audit_dir = audit_dir
        self.tnbc_map = tnbc_map
        self.explanations = explanations
        self.tnbc_sources = tnbc_sources
        self.reference_dir = reference_dir
        self.release_assessment = release_assessment
        self.source_identity = source_identity
        self.release_status = str(release_assessment["release_status"])
        self.publication_ready = bool(release_assessment["publication_ready"])
        self.reader_cache: dict[Path, PdfReader] = {}
        self.embedding_rows: list[dict[str, Any]] = []
        self.atlas_reference_rows: list[dict[str, Any]] = []

    @property
    def defaults(self) -> dict[str, Any]:
        return self.spec["page_defaults"]

    def new_page(self, width: float, height: float) -> PageObject:
        return PageObject.create_blank_page(width=width, height=height)

    def pdf_reader(self, path: Path) -> PdfReader:
        path = path.resolve()
        if path not in self.reader_cache:
            self.reader_cache[path] = PdfReader(str(path), strict=False)
        return self.reader_cache[path]

    def source_size(self, source: SourceSlice) -> tuple[float, float]:
        if source.path.suffix.lower() == ".pdf":
            page = self.pdf_reader(source.path).pages[source.page_index]
            return float(page.mediabox.width), float(page.mediabox.height)
        with Image.open(source.path) as image:
            return float(image.width), float(image.height)

    def merge_pdf_slice(
        self,
        target: PageObject,
        source: SourceSlice,
        dest: tuple[float, float, float, float],
        *,
        output_name: str,
        output_page: int,
        panel: str,
    ) -> None:
        reader = self.pdf_reader(source.path)
        original = reader.pages[source.page_index]
        page = copy.copy(original)
        media_left = float(page.mediabox.left)
        media_bottom = float(page.mediabox.bottom)
        media_width = float(page.mediabox.width)
        media_height = float(page.mediabox.height)
        x0, y0, x1, y1 = source.crop_top
        crop_left = media_left + x0 * media_width
        crop_right = media_left + x1 * media_width
        crop_bottom = media_bottom + (1.0 - y1) * media_height
        crop_top = media_bottom + (1.0 - y0) * media_height
        page.cropbox.lower_left = (crop_left, crop_bottom)
        page.cropbox.upper_right = (crop_right, crop_top)
        source_width = crop_right - crop_left
        source_height = crop_top - crop_bottom
        dx, dy, dw, dh = dest
        scale = min(dw / source_width, dh / source_height)
        tx = dx + (dw - source_width * scale) / 2.0
        ty = dy + (dh - source_height * scale) / 2.0
        transform = (
            Transformation()
            .translate(-crop_left, -crop_bottom)
            .scale(scale, scale)
            .translate(tx, ty)
        )
        target.merge_transformed_page(page, transform, over=True, expand=False)
        self.embedding_rows.append(
            {
                "output": output_name,
                "output_page": output_page,
                "panel": panel,
                "atom_id": source.atom_id,
                "source_role": source.role,
                "source": str(source.path.resolve()),
                "source_page_1based": source.page_index + 1,
                "source_sha256": sha256_file(source.path),
                "source_crop_top": ",".join(f"{value:.6g}" for value in source.crop_top),
                "destination_pt": ",".join(f"{value:.3f}" for value in dest),
                "embedding_method": "pypdf.merge_transformed_page_vector",
            }
        )

    def merge_png_slice(
        self,
        target: PageObject,
        source: SourceSlice,
        dest: tuple[float, float, float, float],
        *,
        output_name: str,
        output_page: int,
        panel: str,
    ) -> None:
        with Image.open(source.path) as original:
            image = original.convert("RGB")
            x0, y0, x1, y1 = source.crop_top
            crop = image.crop(
                (
                    round(x0 * image.width),
                    round(y0 * image.height),
                    round(x1 * image.width),
                    round(y1 * image.height),
                )
            )
            png = io.BytesIO()
            crop.save(png, format="PNG")
            png.seek(0)
            dx, dy, dw, dh = dest
            scale = min(dw / crop.width, dh / crop.height)
            width = crop.width * scale
            height = crop.height * scale
            tx = dx + (dw - width) / 2.0
            ty = dy + (dh - height) / 2.0
            overlay = io.BytesIO()
            cv = canvas.Canvas(overlay, pagesize=(float(target.mediabox.width), float(target.mediabox.height)))
            cv.drawImage(ImageReader(png), tx, ty, width=width, height=height, mask="auto")
            cv.save()
            overlay.seek(0)
            target.merge_page(PdfReader(overlay).pages[0], over=True, expand=False)
        self.embedding_rows.append(
            {
                "output": output_name,
                "output_page": output_page,
                "panel": panel,
                "atom_id": source.atom_id,
                "source_role": source.role,
                "source": str(source.path.resolve()),
                "source_page_1based": source.page_index + 1,
                "source_sha256": sha256_file(source.path),
                "source_crop_top": ",".join(f"{value:.6g}" for value in source.crop_top),
                "destination_pt": ",".join(f"{value:.3f}" for value in dest),
                "embedding_method": "reference_png_crop_no_luad_atom_rasterization",
            }
        )

    def merge_source(
        self,
        target: PageObject,
        source: SourceSlice,
        dest: tuple[float, float, float, float],
        *,
        output_name: str,
        output_page: int,
        panel: str,
    ) -> None:
        if source.path.suffix.lower() == ".pdf":
            self.merge_pdf_slice(
                target, source, dest, output_name=output_name, output_page=output_page, panel=panel
            )
        elif source.role == "TNBC_REFERENCE_CROP" and source.path.suffix.lower() in {".png", ".jpg", ".jpeg"}:
            self.merge_png_slice(
                target, source, dest, output_name=output_name, output_page=output_page, panel=panel
            )
        else:
            raise ValueError(
                f"LUAD atoms must be PDF; only TNBC reference crops may be raster images: {source.path}"
            )

    def overlay_page(self, width: float, height: float, draw_fn: Any) -> PageObject:
        stream = io.BytesIO()
        cv = canvas.Canvas(stream, pagesize=(width, height), pageCompression=1)
        draw_fn(cv)
        # ReportLab may emit a zero-page PDF when an atlas overlay has no
        # labels, placeholders, or boundary text.  Finalise the current page
        # explicitly so pypdf always receives one valid transparent overlay.
        cv.showPage()
        cv.save()
        stream.seek(0)
        return PdfReader(stream).pages[0]

    def page_rect(self, width: float, height: float) -> tuple[float, float, float, float]:
        margin = float(self.defaults["margin_pt"])
        header = float(self.defaults["header_pt"])
        footer = float(self.defaults["footer_pt"])
        return margin, margin + footer, width - 2 * margin, height - 2 * margin - header - footer

    def normalize_box(
        self, box: Sequence[float], content: tuple[float, float, float, float]
    ) -> tuple[float, float, float, float]:
        cx, cy, cw, ch = content
        x, y_top, w, h = (float(value) for value in box)
        return cx + x * cw, cy + (1.0 - y_top - h) * ch, w * cw, h * ch

    def grid_boxes(
        self,
        box: tuple[float, float, float, float],
        count: int,
        requested: Sequence[int] | None,
    ) -> list[tuple[float, float, float, float]]:
        if count <= 0:
            return []
        if requested:
            rows, cols = int(requested[0]), int(requested[1])
            if rows * cols < count:
                rows = math.ceil(count / cols)
        else:
            cols = math.ceil(math.sqrt(count))
            rows = math.ceil(count / cols)
        gutter = float(self.defaults["gutter_pt"])
        x, y, width, height = box
        cell_width = (width - gutter * (cols - 1)) / cols
        cell_height = (height - gutter * (rows - 1)) / rows
        cells: list[tuple[float, float, float, float]] = []
        for index in range(count):
            row = index // cols
            col = index % cols
            cells.append(
                (
                    x + col * (cell_width + gutter),
                    y + (rows - 1 - row) * (cell_height + gutter),
                    cell_width,
                    cell_height,
                )
            )
        return cells

    def units_for_item(self, item: dict[str, Any]) -> list[tuple[str, SourceSlice | None, str, str]]:
        units: list[tuple[str, SourceSlice | None, str, str]] = []
        for atom_id in item["atoms"]:
            result = self.resolver.resolve(atom_id)
            if result.sources:
                units.extend((atom_id, source, result.status, result.note) for source in result.sources)
            else:
                units.append((atom_id, None, result.status, result.note))
        return units

    def placeholder(self, cv: canvas.Canvas, box: tuple[float, float, float, float], status: str, note: str) -> None:
        x, y, width, height = box
        cv.setFillColor(HexColor("#F7F7F7"))
        cv.setStrokeColor(HexColor("#9B9B9B"))
        cv.setLineWidth(0.55)
        cv.rect(x + 1, y + 1, max(0, width - 2), max(0, height - 2), fill=1, stroke=1)
        cv.setFillColor(HexColor("#333333"))
        cv.setFont(self.fonts.bold, min(9.0, max(6.0, width / 35.0)))
        cv.drawCentredString(x + width / 2, y + height / 2 + 4, status)
        if note:
            cv.setFillColor(HexColor("#666666"))
            cv.setFont(self.fonts.regular, min(6.5, max(5.0, width / 55.0)))
            lines = self.wrap_text(note, width - 14, self.fonts.regular, min(6.5, max(5.0, width / 55.0)))
            for index, line in enumerate(lines[:3]):
                cv.drawCentredString(x + width / 2, y + height / 2 - 8 - index * 7, line)

    def wrap_text(self, text: str, width: float, font: str, size: float) -> list[str]:
        self.fonts.ensure(text)
        if not text:
            return []
        tokens = text.split(" ") if " " in text else list(text)
        separator = " " if " " in text else ""
        lines: list[str] = []
        current = ""
        for token in tokens:
            trial = token if not current else current + separator + token
            if pdfmetrics.stringWidth(trial, font, size) <= width or not current:
                current = trial
            else:
                lines.append(current)
                current = token
        if current:
            lines.append(current)
        return lines

    def draw_wrapped(
        self,
        cv: canvas.Canvas,
        text: str,
        x: float,
        top: float,
        width: float,
        font: str,
        size: float,
        leading: float,
        max_lines: int | None = None,
    ) -> float:
        lines = self.wrap_text(text, width, font, size)
        if max_lines is not None:
            lines = lines[:max_lines]
        cv.setFont(font, size)
        for index, line in enumerate(lines):
            cv.drawString(x, top - index * leading, line)
        return top - len(lines) * leading

    def render_item(
        self,
        page: PageObject,
        item: dict[str, Any],
        dest: tuple[float, float, float, float],
        *,
        output_name: str,
        output_page: int,
        draw_label: bool = True,
    ) -> PageObject:
        panel = str(item["panel"])
        label_band = 12.0 if draw_label else 0.0
        x, y, width, height = dest
        atom_box = (x, y, width, max(1.0, height - label_band))
        units = self.units_for_item(item)
        cells = self.grid_boxes(atom_box, len(units), item.get("grid"))
        placeholders: list[tuple[tuple[float, float, float, float], str, str]] = []
        for cell, (atom_id, source, status, note) in zip(cells, units):
            if source is None:
                placeholders.append((cell, status, note))
            else:
                self.merge_source(
                    page,
                    source,
                    cell,
                    output_name=output_name,
                    output_page=output_page,
                    panel=panel,
                )

        boundary = str(item.get("boundary", ""))
        if not boundary:
            boundaries = [self.resolver.resolve(atom_id).boundary for atom_id in item["atoms"]]
            boundary = next((value for value in boundaries if value), "")

        def draw(cv: canvas.Canvas) -> None:
            for cell, status, note in placeholders:
                self.placeholder(cv, cell, status, note)
            if draw_label:
                cv.setFillColor(HexColor("#111111"))
                cv.setFont(self.fonts.bold, float(self.defaults["panel_label_pt"]))
                cv.drawString(x + 1, y + height - 10, str(item.get("label", "")))
            if boundary:
                cv.setFillColor(HexColor("#555555"))
                cv.setFont(self.fonts.regular, 5.8)
                cv.drawRightString(x + width - 1, y + height - 8, boundary)

        page.merge_page(self.overlay_page(float(page.mediabox.width), float(page.mediabox.height), draw), over=True)
        return page

    def page_chrome(
        self,
        page: PageObject,
        title: str,
        page_no: int,
        total: int,
    ) -> None:
        width = float(page.mediabox.width)
        height = float(page.mediabox.height)
        margin = float(self.defaults["margin_pt"])
        footer = (
            f"{self.release_status} | publication_ready="
            f"{str(self.publication_ready).lower()} | official atoms; layout only"
        )

        def draw(cv: canvas.Canvas) -> None:
            self.fonts.ensure(title)
            cv.setFillColor(HexColor("#222222"))
            cv.setFont(self.fonts.regular, 7.2)
            cv.drawString(margin, height - margin + 1, title)
            cv.setStrokeColor(HexColor("#D2D2D2"))
            cv.setLineWidth(0.4)
            cv.line(margin, height - margin - 4, width - margin, height - margin - 4)
            cv.setFillColor(HexColor("#777777"))
            cv.setFont(self.fonts.regular, 5.8)
            cv.drawString(margin, margin - 2, footer)
            cv.drawRightString(width - margin, margin - 2, f"{page_no}/{total}")

        page.merge_page(self.overlay_page(width, height, draw), over=True)

    def atomic_write(self, writer: PdfWriter, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp_dir = path.parent / "tmp" / "pdfs"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        tmp_path = tmp_dir / (path.name + ".partial")
        with tmp_path.open("wb") as handle:
            writer.write(handle)
        validate_pdf(tmp_path)
        os.replace(tmp_path, path)
        try:
            tmp_dir.rmdir()
            tmp_dir.parent.rmdir()
        except OSError:
            pass

    def build_configured_outputs(self) -> list[dict[str, Any]]:
        final_rows: list[dict[str, Any]] = []
        for output in self.spec["outputs"]:
            width = float(output.get("width_pt", self.defaults["width_pt"]))
            height = float(output.get("height_pt", self.defaults["height_pt"]))
            writer = PdfWriter()
            writer.add_metadata(
                {
                    "/Title": str(output["title"]),
                    "/Subject": (
                        "Strict TNBC official-code port to LUAD; layout-only assembly; "
                        f"{self.release_status}; publication_ready="
                        f"{str(self.publication_ready).lower()}"
                    ),
                    "/Creator": "Paper2Paper official-port assembly",
                }
            )
            pages = output["pages"]
            for page_index, page_spec in enumerate(pages, start=1):
                page = self.new_page(width, height)
                content = self.page_rect(width, height)
                for item in page_spec["items"]:
                    dest = self.normalize_box(item["box"], content)
                    self.render_item(
                        page,
                        item,
                        dest,
                        output_name=output["filename"],
                        output_page=page_index,
                    )
                self.page_chrome(
                    page,
                    str(page_spec.get("title", output["title"])),
                    page_index,
                    len(pages),
                )
                writer.add_page(page)
            path = self.output_dir / output["filename"]
            self.atomic_write(writer, path)
            pages_count, digest = validate_pdf(path)
            final_rows.append(
                {
                    "output": str(path.resolve()),
                    "filename": path.name,
                    "pages": pages_count,
                    "bytes": path.stat().st_size,
                    "sha256": digest,
                }
            )
        return final_rows

    def reference_png_path(self, source_key: str, page_1based: int) -> Path | None:
        if not self.reference_dir or not self.reference_dir.is_dir():
            return None
        if source_key == "MAIN_ARTICLE":
            name = f"s41467-026-76385-8_reference_{page_1based:02d}.png"
        elif source_key == "SUPPLEMENT_PDF":
            name = f"41467_2026_76385_MOESM1_ESM_{page_1based:02d}.png"
        else:
            return None
        path = self.reference_dir / name
        return path.resolve() if path.is_file() else None

    def tnbc_slice(self, panel: str) -> SourceSlice | None:
        row = self.tnbc_map.get(panel)
        if not row:
            self.resolver.missing_rows.append(
                {
                    "atom_id": panel,
                    "status": "MISSING_TNBC_CROP_MAP",
                    "root_key": "TNBC_REFERENCE",
                    "patterns": "",
                    "detail": "No audited TNBC crop coordinates.",
                    "is_required": "YES",
                }
            )
            return None
        source_key = row["source_key"]
        page_1based = int(row["page_1based"])
        path = self.tnbc_sources.get(source_key)
        if not path or not path.is_file():
            path = self.reference_png_path(source_key, page_1based)
        if not path or not path.is_file():
            self.resolver.missing_rows.append(
                {
                    "atom_id": panel,
                    "status": "MISSING_TNBC_REFERENCE",
                    "root_key": source_key,
                    "patterns": "",
                    "detail": "Neither source PDF nor reference PNG was available.",
                    "is_required": "YES",
                }
            )
            return None
        crop = tuple(float(row[key]) for key in ("x0_top", "y0_top", "x1_top", "y1_top"))
        source_page = page_1based - 1 if path.suffix.lower() == ".pdf" else 0
        source = SourceSlice(
            path=path.resolve(),
            page_index=source_page,
            crop_top=crop,
            role="TNBC_REFERENCE_CROP",
            atom_id=panel,
        )
        self.atlas_reference_rows.append(
            {
                "panel": panel,
                "source_key": source_key,
                "source": str(path.resolve()),
                "source_page_1based": page_1based,
                "crop_top": ",".join(f"{value:.6g}" for value in crop),
                "source_format": path.suffix.lower().lstrip("."),
            }
        )
        return source

    def atlas_items(self) -> list[dict[str, Any]]:
        ordered: list[dict[str, Any]] = []
        seen: set[str] = set()
        for output in self.spec["outputs"]:
            for page in output["pages"]:
                for item in page["items"]:
                    panel = str(item["panel"])
                    if panel not in seen:
                        ordered.append(item)
                        seen.add(panel)
        return ordered

    def explanation_lines(self, panel: str) -> list[tuple[str, str]]:
        row = self.explanations.get(panel, {})
        return [
            ("TNBC biological role", row.get("biological_role_tnbc", "Explanation not supplied.")),
            ("LUAD biological role", row.get("biological_role_luad", "Explanation not supplied.")),
            ("Interpretation boundary", row.get("interpretation_boundary", "Explanation not supplied.")),
            ("Visual comparison", row.get("visual_comparison", "Explanation not supplied.")),
        ]

    def audit_atlas_inputs(self) -> None:
        """Validate comparison references and prose without producing an atlas."""
        required_explanation_fields = (
            "biological_role_tnbc",
            "biological_role_luad",
            "interpretation_boundary",
            "visual_comparison",
        )
        for item in self.atlas_items():
            panel = str(item["panel"])
            row = self.explanations.get(panel)
            if not row or any(not row.get(field, "").strip() for field in required_explanation_fields):
                self.resolver.missing_rows.append(
                    {
                        "atom_id": panel,
                        "status": "MISSING_PANEL_EXPLANATION",
                        "root_key": "ATLAS_EXPLANATIONS",
                        "patterns": "",
                        "detail": "One or more required explanation fields are absent.",
                        "is_required": "YES",
                    }
                )
            source = self.tnbc_slice(panel)
            if not source:
                continue
            try:
                if source.path.suffix.lower() == ".pdf":
                    pages, _ = validate_pdf(source.path)
                    if source.page_index < 0 or source.page_index >= pages:
                        raise ValueError(
                            f"TNBC source has {pages} pages; requested {source.page_index + 1}"
                        )
                else:
                    with Image.open(source.path) as reference:
                        reference.verify()
            except Exception as exc:
                self.resolver.missing_rows.append(
                    {
                        "atom_id": panel,
                        "status": "INVALID_TNBC_REFERENCE",
                        "root_key": "TNBC_REFERENCE",
                        "patterns": "",
                        "detail": f"{source.path}: {exc}",
                        "is_required": "YES",
                    }
                )

    def build_atlas(self) -> dict[str, Any]:
        width, height = 792.0, 612.0
        writer = PdfWriter()
        writer.add_metadata(
            {
                "/Title": "TNBC vs LUAD panel comparison atlas",
                "/Subject": (
                    "TNBC reference crops versus strict official-code LUAD atomic PDFs; "
                    f"{self.release_status}; publication_ready="
                    f"{str(self.publication_ready).lower()}"
                ),
                "/Creator": "Paper2Paper official-port assembly",
            }
        )
        items = self.atlas_items()
        for page_index, item in enumerate(items, start=1):
            panel = str(item["panel"])
            page = self.new_page(width, height)
            left = (24.0, 188.0, 362.0, 360.0)
            right = (406.0, 188.0, 362.0, 360.0)
            tnbc = self.tnbc_slice(panel)
            tnbc_placeholder = tnbc is None
            if tnbc:
                self.merge_source(
                    page,
                    tnbc,
                    left,
                    output_name="TNBC_vs_LUAD_comparison_atlas.pdf",
                    output_page=page_index,
                    panel=panel,
                )
            atlas_item = dict(item)
            atlas_item["label"] = ""
            self.render_item(
                page,
                atlas_item,
                right,
                output_name="TNBC_vs_LUAD_comparison_atlas.pdf",
                output_page=page_index,
                draw_label=False,
            )
            explanations = self.explanation_lines(panel)

            def draw(cv: canvas.Canvas) -> None:
                cv.setFillColor(HexColor("#222222"))
                cv.setFont(self.fonts.bold, 9.5)
                cv.drawString(24, 585, panel)
                cv.setFont(self.fonts.regular, 7.0)
                cv.setFillColor(HexColor("#555555"))
                cv.drawRightString(768, 585, f"{page_index}/{len(items)}")
                cv.setStrokeColor(HexColor("#D0D0D0"))
                cv.setLineWidth(0.45)
                cv.line(24, 576, 768, 576)
                cv.setFillColor(HexColor("#333333"))
                cv.setFont(self.fonts.bold, 7.0)
                cv.drawString(left[0], 558, "TNBC reference panel")
                cv.drawString(right[0], 558, "LUAD official-code atom")
                if tnbc_placeholder:
                    self.placeholder(cv, left, "MISSING_TNBC_REFERENCE", "See missing_inputs.tsv")
                cv.setStrokeColor(HexColor("#D0D0D0"))
                cv.line(24, 174, 768, 174)
                y = 162.0
                for label, value in explanations:
                    self.fonts.ensure(label + value)
                    cv.setFillColor(HexColor("#333333"))
                    cv.setFont(self.fonts.bold, 6.3)
                    cv.drawString(24, y, label + ":")
                    label_width = pdfmetrics.stringWidth(label + ": ", self.fonts.bold, 6.3)
                    cv.setFont(self.fonts.regular, 6.3)
                    lines = self.wrap_text(value, 744 - label_width, self.fonts.regular, 6.3)
                    if lines:
                        cv.drawString(24 + label_width, y, lines[0])
                        for extra_index, line in enumerate(lines[1:3], start=1):
                            cv.drawString(24, y - extra_index * 8, line)
                        y -= max(1, min(3, len(lines))) * 8 + 2
                    else:
                        y -= 10
                cv.setFillColor(HexColor("#666666"))
                cv.setFont(self.fonts.regular, 5.8)
                cv.drawString(
                    24,
                    10,
                    f"{self.release_status} | publication_ready="
                    f"{str(self.publication_ready).lower()} | comparison for review",
                )

            page.merge_page(self.overlay_page(width, height, draw), over=True)
            writer.add_page(page)

        path = self.output_dir / "TNBC_vs_LUAD_comparison_atlas.pdf"
        self.atomic_write(writer, path)
        pages, digest = validate_pdf(path)
        return {
            "output": str(path.resolve()),
            "filename": path.name,
            "pages": pages,
            "bytes": path.stat().st_size,
            "sha256": digest,
        }

    def write_input_contract(self) -> None:
        usages: dict[str, list[str]] = {atom_id: [] for atom_id in self.spec["atoms"]}
        for output in self.spec["outputs"]:
            for page_index, page in enumerate(output["pages"], start=1):
                for item in page["items"]:
                    for atom_id in item["atoms"]:
                        usages.setdefault(atom_id, []).append(
                            f"{output['filename']}#page={page_index}:{item['panel']}"
                        )
        rows: list[dict[str, Any]] = []
        for atom_id, atom in self.spec["atoms"].items():
            patterns = [str(value) for value in atom.get("patterns", [])]
            scopes = sorted(
                {
                    scope
                    for scope in (
                        self.resolver.fixed_pattern_scope(pattern) for pattern in patterns
                    )
                    if scope
                }
            )
            rows.append(
                {
                    "atom_id": atom_id,
                    "declared_status": atom.get("status", "OFFICIAL_PDF_REQUIRED"),
                    "required": "YES" if atom.get("required", False) else "NO",
                    "root_key": atom.get("root", ""),
                    "path_scope": ";".join(scopes) or atom.get("root", ""),
                    "patterns": ";".join(patterns),
                    "match_mode": atom.get("match", "one"),
                    "publication_boundary": atom.get("boundary", ""),
                    "vector_pdf_only": "YES" if not atom.get("status") else "N/A",
                    "used_by": ";".join(usages.get(atom_id, [])),
                }
            )
        write_tsv(
            self.audit_dir / "input_contract.tsv",
            rows,
            (
                "atom_id", "declared_status", "required", "root_key", "path_scope",
                "patterns", "match_mode", "publication_boundary", "vector_pdf_only",
                "used_by",
            ),
        )

    def write_audit(self, final_rows: list[dict[str, Any]], run_mode: str = "ASSEMBLY") -> None:
        resolution_fields = (
            "atom_id", "status", "root_key", "path", "page_1based", "page_count",
            "sha256", "crop_top", "boundary",
        )
        missing_fields = (
            "atom_id", "status", "root_key", "patterns", "detail", "is_required",
        )
        boundary_fields = ("atom_id", "status", "note", "is_missing")
        embedding_fields = (
            "output", "output_page", "panel", "atom_id", "source_role", "source",
            "source_page_1based", "source_sha256", "source_crop_top", "destination_pt",
            "embedding_method",
        )
        reference_fields = (
            "panel", "source_key", "source", "source_page_1based", "crop_top", "source_format",
        )
        final_fields = ("output", "filename", "pages", "bytes", "sha256")
        write_tsv(self.audit_dir / "resolved_atoms.tsv", self.resolver.resolution_rows, resolution_fields)
        write_tsv(self.audit_dir / "missing_inputs.tsv", self.resolver.missing_rows, missing_fields)
        write_tsv(self.audit_dir / "declared_boundaries.tsv", self.resolver.boundary_rows, boundary_fields)
        write_tsv(self.audit_dir / "embedding_manifest.tsv", self.embedding_rows, embedding_fields)
        write_tsv(self.audit_dir / "atlas_reference_crops.tsv", self.atlas_reference_rows, reference_fields)
        write_tsv(self.audit_dir / "final_output_sha256.tsv", final_rows, final_fields)
        self.write_input_contract()
        (self.audit_dir / "release_status.json").write_text(
            json.dumps(self.release_assessment, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        identity_fields = (
            "root_key", "atom_id", "source", "bytes_observed", "sha256_observed",
            "bytes_expected", "sha256_expected", "identity_status", "commit_evidence",
            "hash_manifest", "receipt_paths", "receipt_sha256s", "detail",
        )
        write_tsv(
            self.audit_dir / "source_identity_manifest.tsv",
            self.source_identity.get("rows", []),
            identity_fields,
        )
        identity_summary = {key: value for key, value in self.source_identity.items() if key != "rows"}
        (self.audit_dir / "source_identity_summary.json").write_text(
            json.dumps(identity_summary, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        summary = {
            "schema_version": 1,
            "run_mode": run_mode,
            "release_status": self.release_status,
            "publication_ready": self.publication_ready,
            "release_blocker_count": len(self.release_assessment.get("release_blockers", [])),
            "source_identity_ready": bool(self.source_identity.get("identity_ready", False)),
            "expected_outputs": list(EXPECTED_FINALS),
            "created_outputs": [row["filename"] for row in final_rows],
            "resolved_atom_pdf_count": len(self.resolver.resolution_rows),
            "declared_boundary_count": len(self.resolver.boundary_rows),
            "missing_required_count": sum(
                1 for row in self.resolver.missing_rows if row.get("is_required") == "YES"
            ),
            "embedding_count": len(self.embedding_rows),
            "luad_embedding_methods": sorted(
                {
                    row["embedding_method"]
                    for row in self.embedding_rows
                    if row["source_role"] == "LUAD_OFFICIAL_ATOM"
                }
            ),
        }
        (self.audit_dir / "assembly_summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )


def load_overrides(path: Path | None) -> dict[str, list[SourceSlice]]:
    if not path:
        return {}
    rows = read_tsv(path)
    overrides: dict[str, list[SourceSlice]] = {}
    for row in rows:
        atom_id = row["atom_id"]
        raw = os.path.expandvars(os.path.expanduser(row["path"]))
        paths: list[Path]
        if any(char in raw for char in "*?["):
            candidate = Path(raw)
            parent = candidate.parent if str(candidate.parent) else Path(".")
            paths = sorted(item.resolve() for item in parent.glob(candidate.name) if item.is_file())
        else:
            paths = [Path(raw).resolve()]
        page = int(row.get("page_1based") or "1") - 1
        crop = tuple(
            float(row.get(key) or default)
            for key, default in zip(
                ("x0_top", "y0_top", "x1_top", "y1_top"),
                (0.0, 0.0, 1.0, 1.0),
            )
        )
        overrides.setdefault(atom_id, []).extend(
            SourceSlice(path=item, page_index=page, crop_top=crop, atom_id=atom_id)
            for item in paths
        )
    return overrides


def keyed_rows(path: Path, key: str) -> dict[str, dict[str, str]]:
    return {row[key]: row for row in read_tsv(path)}


def assess_supp3_release(
    resolver: AtomResolver,
    supp3_root: Path | None,
    explicit_receipt: Path | None,
) -> dict[str, Any]:
    """Bind the Supplementary Figure 3 release verdict to the exact PDF atoms."""
    expected_panels = [f"SupplementaryFigure3{letter}" for letter in "ABCDEFGHIJKL"]
    assessment: dict[str, Any] = {
        "schema_version": 1,
        "release_status": "STRICT_OFFICIAL_PORT_REVIEW_ONLY",
        "publication_ready": False,
        "upstream_component": "Supplementary Figure 3",
        "receipt_path": "",
        "receipt_sha256": "",
        "upstream_verdict": "MISSING_OR_INVALID_RECEIPT",
        "upstream_publication_ready": False,
        "receipt_contract_valid": False,
        "partial_panels": [],
        "visual_warning_panels": [],
        "visual_warnings": [],
        "release_blockers": [],
    }

    candidates: list[Path] = []
    root = supp3_root.resolve() if supp3_root else None
    if explicit_receipt:
        candidate = explicit_receipt.resolve()
        if root:
            try:
                candidate.relative_to(root)
            except ValueError:
                assessment["release_blockers"].append(
                    "Explicit Supplementary Figure 3 receipt is outside --supp3-root."
                )
            else:
                candidates = [candidate]
        else:
            candidates = [candidate]
    elif root and root.is_dir():
        candidates = sorted(path.resolve() for path in root.rglob("verification_receipt.json"))

    if len(candidates) != 1 or not candidates[0].is_file():
        detail = (
            "No verification_receipt.json found below Supplementary Figure 3 root."
            if not candidates
            else "Ambiguous Supplementary Figure 3 receipts: "
            + ";".join(str(path) for path in candidates)
        )
        assessment["release_blockers"].append(detail)
        resolver.missing_rows.append(
            {
                "atom_id": "SupplementaryFigure3_verification_receipt",
                "status": "MISSING_OR_AMBIGUOUS_UPSTREAM_VERIFICATION_RECEIPT",
                "root_key": "supp3",
                "patterns": "**/verification_receipt.json",
                "detail": detail,
                "is_required": "YES",
            }
        )
        return assessment

    receipt_path = candidates[0]
    assessment["receipt_path"] = str(receipt_path)
    assessment["receipt_sha256"] = sha256_file(receipt_path)
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except Exception as exc:
        detail = f"Cannot parse Supplementary Figure 3 receipt: {exc}"
        assessment["release_blockers"].append(detail)
        resolver.missing_rows.append(
            {
                "atom_id": "SupplementaryFigure3_verification_receipt",
                "status": "INVALID_UPSTREAM_VERIFICATION_RECEIPT",
                "root_key": "supp3",
                "patterns": str(receipt_path),
                "detail": detail,
                "is_required": "YES",
            }
        )
        return assessment

    contract_errors: list[str] = []
    required_keys = {
        "verdict",
        "publication_ready",
        "official_plotting_code_used",
        "plotting_constructor_rewrites",
        "panels_verified",
        "partial_panels",
        "visual_warnings",
        "results",
    }
    missing_keys = sorted(required_keys - set(receipt))
    if missing_keys:
        contract_errors.append("Missing receipt keys: " + ",".join(missing_keys))
    if receipt.get("official_plotting_code_used") is not True:
        contract_errors.append("official_plotting_code_used is not true")
    if receipt.get("plotting_constructor_rewrites") != 0:
        contract_errors.append("plotting_constructor_rewrites is not zero")
    if receipt.get("panels_verified") != len(expected_panels):
        contract_errors.append("panels_verified is not 12")

    result_rows = receipt.get("results")
    if not isinstance(result_rows, list):
        result_rows = []
        contract_errors.append("results is not a list")
    result_by_panel = {
        str(row.get("panel")): row for row in result_rows if isinstance(row, dict)
    }
    if set(result_by_panel) != set(expected_panels):
        contract_errors.append("results does not contain exactly SupplementaryFigure3A-L")

    for panel in expected_panels:
        resolved = resolver.resolve(panel)
        row = result_by_panel.get(panel)
        if len(resolved.sources) != 1 or not row:
            continue
        source = resolved.sources[0].path.resolve()
        reported = Path(str(row.get("pdf", "")))
        if reported.name != source.name:
            contract_errors.append(
                f"{panel} receipt PDF basename {reported.name!r} != resolved {source.name!r}"
            )
        try:
            reported_bytes = int(row.get("pdf_bytes"))
        except (TypeError, ValueError):
            contract_errors.append(f"{panel} has invalid pdf_bytes")
        else:
            if reported_bytes != source.stat().st_size:
                contract_errors.append(
                    f"{panel} receipt bytes {reported_bytes} != resolved {source.stat().st_size}"
                )

    partial_panels = receipt.get("partial_panels", [])
    warnings = receipt.get("visual_warnings", [])
    if not isinstance(partial_panels, list):
        partial_panels = []
        contract_errors.append("partial_panels is not a list")
    if not isinstance(warnings, list):
        warnings = []
        contract_errors.append("visual_warnings is not a list")
    warning_panels = [
        str(item.get("panel")) for item in warnings if isinstance(item, dict) and item.get("panel")
    ]

    assessment.update(
        {
            "upstream_verdict": str(receipt.get("verdict", "UNKNOWN")),
            "upstream_publication_ready": receipt.get("publication_ready") is True,
            "receipt_contract_valid": not contract_errors,
            "partial_panels": [str(panel) for panel in partial_panels],
            "visual_warning_panels": warning_panels,
            "visual_warnings": warnings,
        }
    )
    blockers = assessment["release_blockers"]
    blockers.extend(contract_errors)
    if receipt.get("publication_ready") is not True:
        blockers.append("Supplementary Figure 3 upstream receipt sets publication_ready=false.")
    if partial_panels:
        blockers.append("Partial Supplementary Figure 3 panels: " + ",".join(partial_panels))
    if warning_panels:
        blockers.append(
            "Supplementary Figure 3 visual warnings: " + ",".join(warning_panels)
        )

    if contract_errors:
        resolver.missing_rows.append(
            {
                "atom_id": "SupplementaryFigure3_verification_receipt",
                "status": "RECEIPT_DOES_NOT_BIND_TO_RESOLVED_ATOMS",
                "root_key": "supp3",
                "patterns": str(receipt_path),
                "detail": "; ".join(contract_errors),
                "is_required": "YES",
            }
        )

    assessment["publication_ready"] = bool(
        assessment["receipt_contract_valid"]
        and assessment["upstream_publication_ready"]
        and not partial_panels
        and not warning_panels
    )
    assessment["release_status"] = (
        "PUBLICATION_READY"
        if assessment["publication_ready"]
        else "STRICT_OFFICIAL_PORT_REVIEW_ONLY"
    )
    return assessment


def assess_source_identity(
    resolver: AtomResolver,
    roots: dict[str, Path | None],
    supp3_release: dict[str, Any],
) -> dict[str, Any]:
    """Verify official-port receipts and bind resolved PDFs to hash/byte evidence."""
    expected_commit = "edf5314ce5b9ee0e2f88b2310e7c2df5619ad888"
    rows: list[dict[str, Any]] = []
    root_summaries: dict[str, Any] = {}

    def one(root: Path, names: Sequence[str]) -> tuple[Path | None, str]:
        found = sorted({path.resolve() for name in names for path in root.rglob(name)})
        if len(found) == 1:
            return found[0], ""
        return None, f"expected one of {list(names)}, found {len(found)}"

    def hash_evidence(path: Path | None) -> list[dict[str, str]]:
        return read_tsv(path) if path and path.is_file() else []

    def row_value(row: dict[str, str], candidates: Sequence[str]) -> str:
        return next((row.get(key, "") for key in candidates if row.get(key, "")), "")

    def match_hash_row(
        evidence: list[dict[str, str]], source: Path, root: Path
    ) -> dict[str, str] | None:
        relative = source.relative_to(root).as_posix()
        exact: list[dict[str, str]] = []
        basename: list[dict[str, str]] = []
        for item in evidence:
            named = row_value(
                item,
                (
                    "relative_path",
                    "artifact",
                    "path",
                    "output",
                    "file",
                    "pdf",
                    "filename",
                ),
            ).replace("\\", "/")
            if not named:
                continue
            if (
                named == relative
                or named.endswith("/" + relative)
                or relative.endswith("/" + named.lstrip("/"))
            ):
                exact.append(item)
            if Path(named).name == source.name:
                basename.append(item)
        if len(exact) == 1:
            return exact[0]
        if len(basename) == 1:
            return basename[0]
        return None

    for root_key in ("fig12", "fig3", "fig45", "supp3"):
        root = roots.get(root_key)
        root_errors: list[str] = []
        receipt_paths: list[Path] = []
        commit_evidence = ""
        hash_manifest: Path | None = None
        hashes_required = root_key in {"fig12", "fig3", "fig45"}
        expected_by_name: dict[str, dict[str, Any]] = {}
        if not root or not root.is_dir():
            root_errors.append("official-port root is unavailable")
            root_path = root.resolve() if root else None
        else:
            root_path = root.resolve()
            hash_manifest, hash_error = one(
                root_path,
                (
                    "official_output_sha256.tsv",
                    "atomic_output_sha256.tsv",
                    "output_sha256.tsv",
                    "pdf_sha256.tsv",
                    "pdf_sha256_manifest.tsv",
                    "atomic_pdf_sha256.tsv",
                ),
            )
            if hash_error and hashes_required:
                root_errors.append("output hash manifest: " + hash_error)

            if root_key == "fig3":
                receipt, error = one(root_path, ("official_port_receipt.json",))
                source_manifest, source_error = one(root_path, ("official_source_manifest.tsv",))
                if error or not receipt:
                    root_errors.append("official port receipt: " + error)
                else:
                    receipt_paths.append(receipt)
                    data = json.loads(receipt.read_text(encoding="utf-8"))
                    if data.get("status") != "PASS" or data.get("plot_constructors_changed") is not False:
                        root_errors.append("Figure 3 official port receipt is not PASS/unchanged")
                if source_error or not source_manifest:
                    root_errors.append("official source manifest: " + source_error)
                else:
                    receipt_paths.append(source_manifest)
                    source_rows = read_tsv(source_manifest)
                    if len(source_rows) != 8 or any(
                        not re.fullmatch(r"[0-9a-f]{64}", row.get("official_sha256", ""))
                        for row in source_rows
                    ):
                        root_errors.append("Figure 3 source-hash lock is incomplete")
                    else:
                        commit_evidence = "edf5314 source-hash lock (8 official scripts)"

            elif root_key == "fig45":
                material, material_error = one(root_path, ("materialization_receipt.json",))
                verify, verify_error = one(root_path, ("official_port_verification.json",))
                for path, error, label in (
                    (material, material_error, "materialization receipt"),
                    (verify, verify_error, "official port verification"),
                ):
                    if error or not path:
                        root_errors.append(f"{label}: {error}")
                    else:
                        receipt_paths.append(path)
                if material:
                    data = json.loads(material.read_text(encoding="utf-8"))
                    if str(data.get("commit", "")) != "edf5314":
                        root_errors.append("Figure 4/5 materialization commit is not edf5314")
                if verify:
                    data = json.loads(verify.read_text(encoding="utf-8"))
                    if data.get("status") != "PASS" or str(data.get("commit", "")) != "edf5314":
                        root_errors.append("Figure 4/5 verification is not PASS at edf5314")
                    else:
                        commit_evidence = "edf5314"

            elif root_key == "supp3":
                material, material_error = one(root_path, ("materialization_receipt.json",))
                source_manifest, source_error = one(root_path, ("official_source_manifest.tsv",))
                if material_error or not material:
                    root_errors.append("materialization receipt: " + material_error)
                else:
                    receipt_paths.append(material)
                    data = json.loads(material.read_text(encoding="utf-8"))
                    if data.get("plotting_constructor_rewrites") != 0:
                        root_errors.append("Supplementary Figure 3 constructor rewrites are not zero")
                if source_error or not source_manifest:
                    root_errors.append("official source manifest: " + source_error)
                else:
                    receipt_paths.append(source_manifest)
                    source_rows = read_tsv(source_manifest)
                    if len(source_rows) != 5 or any(
                        not re.fullmatch(r"[0-9a-f]{64}", row.get("sha256", ""))
                        for row in source_rows
                    ):
                        root_errors.append("Supplementary Figure 3 source-hash lock is incomplete")
                    else:
                        commit_evidence = "edf5314 source-hash lock (5 official scripts)"
                release_path = Path(str(supp3_release.get("receipt_path", "")))
                if release_path.is_file():
                    receipt_paths.append(release_path.resolve())
                for item in supp3_release.get("visual_warnings", []):
                    if isinstance(item, dict):
                        pass
                try:
                    release_data = json.loads(release_path.read_text(encoding="utf-8"))
                    expected_by_name = {
                        Path(str(item.get("pdf", ""))).name: item
                        for item in release_data.get("results", [])
                        if isinstance(item, dict)
                    }
                except Exception:
                    root_errors.append("Supplementary Figure 3 verification receipt is unreadable")

            elif root_key == "fig12":
                manifests = sorted(root_path.glob("runtime_*corrected/audit/runtime_manifest.json"))
                if len(manifests) != 2:
                    root_errors.append(
                        f"expected final and secondary runtime manifests, found {len(manifests)}"
                    )
                for manifest in manifests:
                    receipt_paths.append(manifest.resolve())
                    data = json.loads(manifest.read_text(encoding="utf-8"))
                    if data.get("official_commit") != expected_commit:
                        root_errors.append(f"{manifest} official_commit mismatch")
                if manifests and not any("runtime_final_corrected" in str(p) for p in manifests):
                    root_errors.append("final corrected runtime manifest is absent")
                if manifests and not any("runtime_secondary_GSE81089_corrected" in str(p) for p in manifests):
                    root_errors.append("secondary corrected runtime manifest is absent")
                commit_evidence = expected_commit if len(manifests) == 2 else ""

        if hash_manifest:
            receipt_paths.append(hash_manifest.resolve())
        receipt_identities = [
            {"path": str(path), "sha256": sha256_file(path)}
            for path in sorted(set(receipt_paths))
            if path.is_file()
        ]
        evidence = hash_evidence(hash_manifest)
        resolved_sources = [
            (atom_id, source.path.resolve())
            for atom_id, resolved in resolver.resolved.items()
            if resolved.root_key == root_key
            for source in resolved.sources
        ]
        for atom_id, source in resolved_sources:
            observed_sha = sha256_file(source)
            observed_bytes = source.stat().st_size
            expected_sha = ""
            expected_bytes = ""
            identity_errors = list(root_errors)
            evidence_row = match_hash_row(evidence, source, root_path) if root_path else None
            if evidence_row:
                expected_sha = row_value(evidence_row, ("sha256", "hash"))
                expected_bytes = row_value(evidence_row, ("bytes", "size"))
                if expected_sha and expected_sha.lower() != observed_sha:
                    identity_errors.append("PDF SHA-256 mismatch")
                if expected_bytes and int(expected_bytes) != observed_bytes:
                    identity_errors.append("PDF byte-count mismatch")
                if evidence_row.get("status") and evidence_row.get("status") != "PASS":
                    identity_errors.append("hash manifest status is not PASS")
            elif hashes_required:
                identity_errors.append("resolved PDF has no unique hash-manifest row")
            elif source.name in expected_by_name:
                expected_bytes = str(expected_by_name[source.name].get("pdf_bytes", ""))
                if not expected_bytes or int(expected_bytes) != observed_bytes:
                    identity_errors.append("verification-receipt byte-count mismatch")
            else:
                identity_errors.append("resolved PDF is absent from upstream verification receipt")

            rows.append(
                {
                    "root_key": root_key,
                    "atom_id": atom_id,
                    "source": str(source),
                    "bytes_observed": observed_bytes,
                    "sha256_observed": observed_sha,
                    "bytes_expected": expected_bytes,
                    "sha256_expected": expected_sha,
                    "identity_status": "VERIFIED" if not identity_errors else "FAILED",
                    "commit_evidence": commit_evidence,
                    "hash_manifest": str(hash_manifest or ""),
                    "receipt_paths": ";".join(str(path) for path in sorted(set(receipt_paths))),
                    "receipt_sha256s": ";".join(
                        f"{item['path']}={item['sha256']}" for item in receipt_identities
                    ),
                    "detail": "; ".join(identity_errors),
                }
            )

        if root_errors or any(row["identity_status"] != "VERIFIED" for row in rows if row["root_key"] == root_key):
            details = list(root_errors)
            details.extend(
                row["detail"] for row in rows
                if row["root_key"] == root_key and row["identity_status"] != "VERIFIED"
            )
            resolver.missing_rows.append(
                {
                    "atom_id": f"{root_key}_official_source_identity",
                    "status": "OFFICIAL_SOURCE_IDENTITY_FAILED",
                    "root_key": root_key,
                    "patterns": "official receipts + output hash manifest",
                    "detail": "; ".join(dict.fromkeys(filter(None, details))),
                    "is_required": "YES",
                }
            )
        root_summaries[root_key] = {
            "root": str(root_path or ""),
            "receipt_paths": [str(path) for path in sorted(set(receipt_paths))],
            "receipts": receipt_identities,
            "commit_evidence": commit_evidence,
            "hash_manifest": str(hash_manifest or ""),
            "resolved_pdf_slices": len(resolved_sources),
            "verified_pdf_slices": sum(
                row["identity_status"] == "VERIFIED" for row in rows if row["root_key"] == root_key
            ),
            "identity_ready": not root_errors and all(
                row["identity_status"] == "VERIFIED" for row in rows if row["root_key"] == root_key
            ) and bool(resolved_sources),
            "errors": list(dict.fromkeys(root_errors)),
        }

    return {
        "schema_version": 1,
        "expected_official_commit": expected_commit,
        "identity_ready": all(item["identity_ready"] for item in root_summaries.values()),
        "roots": root_summaries,
        "rows": rows,
    }


def main() -> int:
    args = parse_args()
    spec = json.loads(args.spec.read_text(encoding="utf-8"))
    if int(spec.get("schema_version", 0)) != 1:
        raise ValueError("Unsupported assembly_spec schema_version")
    output_dir = args.output_dir.resolve()
    audit_dir = (
        args.audit_dir.resolve()
        if args.audit_dir
        else (output_dir.parent / "audit").resolve()
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    roots: dict[str, Path | None] = {
        "fig12": args.fig12_root.resolve(),
        "fig3": args.fig3_root.resolve(),
        "fig45": args.fig45_root.resolve(),
        "supp3": args.supp3_root.resolve() if args.supp3_root else None,
    }
    resolver = AtomResolver(spec["atoms"], roots, load_overrides(args.source_override))
    resolver.resolve_all()
    release_assessment = assess_supp3_release(
        resolver,
        roots["supp3"],
        args.supp3_verification_receipt,
    )
    source_identity = assess_source_identity(resolver, roots, release_assessment)
    tnbc_sources = {
        "MAIN_ARTICLE": args.tnbc_main_pdf.resolve() if args.tnbc_main_pdf else None,
        "SUPPLEMENT_PDF": args.tnbc_supplement_pdf.resolve() if args.tnbc_supplement_pdf else None,
    }
    assembler = PDFAssembler(
        spec=spec,
        resolver=resolver,
        font_book=FontBook(args.font_path),
        output_dir=output_dir,
        audit_dir=audit_dir,
        tnbc_map=keyed_rows(args.tnbc_panel_map, "panel"),
        explanations=keyed_rows(args.explanations, "panel"),
        tnbc_sources=tnbc_sources,
        reference_dir=args.tnbc_reference_dir.resolve() if args.tnbc_reference_dir else None,
        release_assessment=release_assessment,
        source_identity=source_identity,
    )
    if args.dry_run:
        assembler.audit_atlas_inputs()
        assembler.write_audit([], run_mode="DRY_RUN_INPUT_CONTRACT")
        missing_required = sum(
            1 for row in resolver.missing_rows if row.get("is_required") == "YES"
        )
        print(f"Dry-run input contract: {len(spec['atoms'])} atoms checked")
        print(f"Required missing/invalid inputs: {missing_required}")
        print(
            f"Release status: {assembler.release_status}; "
            f"publication_ready={str(assembler.publication_ready).lower()}"
        )
        print(f"Audit manifests: {audit_dir}")
        if missing_required:
            return 4
        if args.require_publication_ready and not assembler.publication_ready:
            print("RELEASE GATE BLOCKED BY UPSTREAM RECEIPT", file=sys.stderr)
            return 5
        return 0

    if args.require_publication_ready:
        assembler.audit_atlas_inputs()
        assembler.write_audit([], run_mode="RELEASE_GATE_PREFLIGHT")
        missing_required = sum(
            1 for row in resolver.missing_rows if row.get("is_required") == "YES"
        )
        if missing_required:
            print(
                f"RELEASE GATE BLOCKED: {missing_required} required inputs are missing/invalid",
                file=sys.stderr,
            )
            return 4
        if not assembler.publication_ready:
            print(
                "RELEASE GATE BLOCKED: STRICT_OFFICIAL_PORT_REVIEW_ONLY; "
                "publication_ready=false",
                file=sys.stderr,
            )
            return 5
        assembler.atlas_reference_rows.clear()

    final_rows = assembler.build_configured_outputs()
    final_rows.append(assembler.build_atlas())
    assembler.write_audit(final_rows, run_mode="ASSEMBLY")
    created = {row["filename"] for row in final_rows}
    if created != set(EXPECTED_FINALS):
        raise RuntimeError(f"Expected seven stable outputs; created {sorted(created)}")
    print(f"Assembled {len(final_rows)} PDFs into {output_dir}")
    print(
        f"Release status: {assembler.release_status}; "
        f"publication_ready={str(assembler.publication_ready).lower()}"
    )
    print(f"Audit manifests: {audit_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ASSEMBLY FAILED: {exc}", file=sys.stderr)
        raise
