#!/usr/bin/env python3
"""Render and verify the seven layout-only official-port PDF deliverables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence

try:
    from PIL import Image, ImageDraw, ImageFont
    from pypdf import PdfReader
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"Missing pypdf/Pillow dependency: {exc}")


# Pillow <9.1 exposes LANCZOS directly on Image; newer releases place it in
# Image.Resampling.  QA must run on the frozen cloud image without changing
# any rendered PDF content.
RESAMPLE_LANCZOS = getattr(getattr(Image, "Resampling", Image), "LANCZOS")


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

EXPECTED_ATOM_COUNT = 75
EXPECTED_PANEL_COUNT = 73
ATLAS_FILENAME = "TNBC_vs_LUAD_comparison_atlas.pdf"
BOUNDARY_TOKENS = (
    "MANUAL_REQUIRED",
    "UNSUPPORTED",
    "SUPPORTED_ZERO",
    "CAPSULE_PRINT_MISMATCH",
    "ANCHOR_BUGFIX",
)
# These two scientific-zero results are frozen by the current Figure 4/5
# official-port receipt.  They are not interchangeable with missing inputs.
FROZEN_SUPPORTED_ZERO = {
    "Figure5F_high": "SUPPORTED_ZERO",
    "SupplementaryFigure7C_high": "SUPPORTED_ZERO",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render and QA strict official-port assembled PDFs.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument("--qa-dir", required=True, type=Path)
    parser.add_argument("--spec", type=Path, default=SCRIPT_DIR / "assembly_spec.json")
    parser.add_argument("--pdftoppm", help="Explicit Poppler pdftoppm executable.")
    parser.add_argument("--dpi", type=int, default=144)
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Fail for missing required atoms/references as well as PDF/render defects. "
            "This is a technical check and does not imply publication readiness."
        ),
    )
    parser.add_argument(
        "--require-publication-ready",
        action="store_true",
        help="Activate the release gate and fail when release_status.json is review-only.",
    )
    return parser.parse_args()


def expected_pages(spec: dict[str, Any]) -> dict[str, int]:
    counts = {output["filename"]: len(output["pages"]) for output in spec["outputs"]}
    seen: set[str] = set()
    atlas_count = 0
    for output in spec["outputs"]:
        for page in output["pages"]:
            for item in page["items"]:
                panel = item["panel"]
                if panel not in seen:
                    seen.add(panel)
                    atlas_count += 1
    counts["TNBC_vs_LUAD_comparison_atlas.pdf"] = atlas_count
    return counts


def ordered_panels(spec: dict[str, Any]) -> list[str]:
    panels: list[str] = []
    seen: set[str] = set()
    for output in spec["outputs"]:
        for page in output["pages"]:
            for item in page["items"]:
                panel = str(item["panel"])
                if panel not in seen:
                    panels.append(panel)
                    seen.add(panel)
    return panels


def expected_boundary_tokens(
    spec: dict[str, Any],
) -> tuple[dict[tuple[str, int], set[str]], dict[str, set[str]]]:
    """Return frozen boundary tokens for configured pages and atlas panels."""
    atom_status = {
        atom_id: str(atom.get("status", ""))
        for atom_id, atom in spec["atoms"].items()
    }
    atom_status.update(FROZEN_SUPPORTED_ZERO)
    configured: dict[tuple[str, int], set[str]] = {}
    atlas: dict[str, set[str]] = {}
    for output in spec["outputs"]:
        for page_index, page in enumerate(output["pages"], start=1):
            page_tokens: set[str] = set()
            for item in page["items"]:
                panel = str(item["panel"])
                panel_tokens: set[str] = set()
                for atom_id in item["atoms"]:
                    status = atom_status.get(str(atom_id), "")
                    if status in BOUNDARY_TOKENS:
                        panel_tokens.add(status)
                    boundary = str(spec["atoms"][str(atom_id)].get("boundary", ""))
                    panel_tokens.update(
                        token for token in BOUNDARY_TOKENS if token in boundary
                    )
                item_boundary = str(item.get("boundary", ""))
                panel_tokens.update(
                    token for token in BOUNDARY_TOKENS if token in item_boundary
                )
                atlas[panel] = panel_tokens
                page_tokens.update(panel_tokens)
            configured[(str(output["filename"]), page_index)] = page_tokens
    return configured, atlas


def parse_crop(value: str) -> tuple[float, ...]:
    if not value:
        return ()
    return tuple(float(item) for item in value.split(","))


def same_path(left: str, right: str) -> bool:
    try:
        return Path(left).resolve() == Path(right).resolve()
    except Exception:
        return left == right


def locate_pdftoppm(explicit: str | None) -> str:
    if explicit:
        path = Path(explicit)
        if path.is_file():
            return str(path.resolve())
        resolved = shutil.which(explicit)
        if resolved:
            return resolved
        raise FileNotFoundError(f"pdftoppm not found: {explicit}")
    resolved = shutil.which("pdftoppm")
    if not resolved:
        raise FileNotFoundError(
            "Poppler pdftoppm is required for final visual QA. Install poppler-utils "
            "or pass --pdftoppm /absolute/path/to/pdftoppm."
        )
    return resolved


def render_pdf(executable: str, pdf: Path, render_dir: Path, dpi: int) -> list[Path]:
    render_dir.mkdir(parents=True, exist_ok=True)
    prefix = render_dir / "page"
    completed = subprocess.run(
        [executable, "-png", "-r", str(dpi), str(pdf), str(prefix)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"pdftoppm failed for {pdf}: {completed.stderr.strip() or completed.stdout.strip()}"
        )
    pages = sorted(render_dir.glob("page-*.png"), key=lambda path: int(path.stem.split("-")[-1]))
    if not pages:
        raise RuntimeError(f"pdftoppm produced no PNG pages for {pdf}")
    return pages


def image_metrics(path: Path) -> dict[str, Any]:
    with Image.open(path) as original:
        image = original.convert("RGB")
        width, height = image.size
        sample = image.copy()
        sample.thumbnail((1000, 1000), RESAMPLE_LANCZOS)
        pixels = list(sample.getdata())
        total = max(1, len(pixels))
        white = sum(1 for red, green, blue in pixels if min(red, green, blue) >= 248)
        dark = sum(1 for red, green, blue in pixels if max(red, green, blue) <= 10)
        ink = sum(1 for red, green, blue in pixels if min(red, green, blue) <= 244)
        return {
            "png": str(path.resolve()),
            "width_px": width,
            "height_px": height,
            "white_fraction": white / total,
            "dark_fraction": dark / total,
            "ink_fraction": ink / total,
        }


def contact_sheet(images: Sequence[Path], destination: Path, title: str) -> None:
    columns = 3
    thumb_width, thumb_height = 230, 300
    caption = 20
    rows = math.ceil(len(images) / columns)
    sheet = Image.new("RGB", (columns * thumb_width, 30 + rows * (thumb_height + caption)), "white")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.load_default()
    draw.text((8, 8), title, fill="black", font=font)
    for index, path in enumerate(images):
        row, col = divmod(index, columns)
        with Image.open(path) as original:
            thumb = original.convert("RGB")
            thumb.thumbnail((thumb_width - 10, thumb_height - 10), RESAMPLE_LANCZOS)
            x = col * thumb_width + (thumb_width - thumb.width) // 2
            y = 30 + row * (thumb_height + caption) + (thumb_height - thumb.height) // 2
            sheet.paste(thumb, (x, y))
        draw.text(
            (col * thumb_width + 5, 30 + row * (thumb_height + caption) + thumb_height),
            f"page {index + 1}",
            fill="black",
            font=font,
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(destination)


def page_content_bytes(page: Any) -> int:
    contents = page.get_contents()
    if contents is None:
        return 0
    try:
        return len(contents.get_data())
    except Exception:
        return 0


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    audit_dir = args.audit_dir.resolve()
    qa_dir = args.qa_dir.resolve()
    qa_dir.mkdir(parents=True, exist_ok=True)
    spec = json.loads(args.spec.read_text(encoding="utf-8"))
    expected = expected_pages(spec)
    atlas_panels = ordered_panels(spec)
    configured_boundaries, atlas_boundaries = expected_boundary_tokens(spec)
    explanation_token_counts: dict[str, dict[str, int]] = {}
    early_explanations_path = args.spec.resolve().parent / "panel_explanations.tsv"
    if early_explanations_path.is_file():
        for row in read_tsv(early_explanations_path):
            panel = row.get("panel", "")
            prose = " ".join(
                row.get(field, "")
                for field in (
                    "biological_role_tnbc",
                    "biological_role_luad",
                    "interpretation_boundary",
                    "visual_comparison",
                )
            )
            explanation_token_counts[panel] = {
                token: prose.count(token) for token in BOUNDARY_TOKENS
            }
    failures: list[str] = []
    pdf_rows: list[dict[str, Any]] = []
    render_rows: list[dict[str, Any]] = []
    page_texts: dict[str, list[str]] = {}
    if len(spec.get("atoms", {})) != EXPECTED_ATOM_COUNT:
        failures.append(
            f"Assembly spec atom count {len(spec.get('atoms', {}))}, "
            f"expected frozen count {EXPECTED_ATOM_COUNT}"
        )
    if len(atlas_panels) != EXPECTED_PANEL_COUNT or len(set(atlas_panels)) != EXPECTED_PANEL_COUNT:
        failures.append(
            f"Assembly spec panel count {len(atlas_panels)} "
            f"({len(set(atlas_panels))} unique), expected {EXPECTED_PANEL_COUNT}"
        )
    if expected.get(ATLAS_FILENAME) != EXPECTED_PANEL_COUNT:
        failures.append(
            f"Atlas page contract {expected.get(ATLAS_FILENAME)}, expected {EXPECTED_PANEL_COUNT}"
        )
    release_status = "UNKNOWN"
    publication_ready = False
    release_path = audit_dir / "release_status.json"
    if not release_path.is_file():
        failures.append("release_status.json was not generated")
    else:
        try:
            release = json.loads(release_path.read_text(encoding="utf-8"))
            release_status = str(release.get("release_status", "UNKNOWN"))
            publication_ready = release.get("publication_ready") is True
            if publication_ready and release_status != "PUBLICATION_READY":
                failures.append(
                    "release_status.json is internally inconsistent: publication_ready=true "
                    f"but release_status={release_status}"
                )
            if not publication_ready and release_status != "STRICT_OFFICIAL_PORT_REVIEW_ONLY":
                failures.append(
                    "release_status.json must label a non-ready assembly as "
                    "STRICT_OFFICIAL_PORT_REVIEW_ONLY"
                )
            if not release.get("receipt_path") or not release.get("receipt_sha256"):
                failures.append("Supplementary Figure 3 upstream receipt provenance is absent")
            if release.get("upstream_component") != "Supplementary Figure 3":
                failures.append("release_status.json does not identify Supplementary Figure 3")
            if release.get("receipt_contract_valid") is not True:
                failures.append("Supplementary Figure 3 receipt contract is not valid")
            if release.get("upstream_publication_ready") is not False:
                failures.append(
                    "Supplementary Figure 3 upstream publication_ready boundary changed"
                )
            expected_partial = {"SupplementaryFigure3A"}
            observed_partial = set(release.get("partial_panels", []))
            if observed_partial != expected_partial:
                failures.append(
                    f"Supplementary Figure 3 partial panels {sorted(observed_partial)}, "
                    f"expected {sorted(expected_partial)}"
                )
            expected_warnings = {
                "SupplementaryFigure3D",
                "SupplementaryFigure3F",
                "SupplementaryFigure3K",
                "SupplementaryFigure3L",
            }
            observed_warnings = set(release.get("visual_warning_panels", []))
            if observed_warnings != expected_warnings:
                failures.append(
                    f"Supplementary Figure 3 visual-warning panels "
                    f"{sorted(observed_warnings)}, expected {sorted(expected_warnings)}"
                )
            receipt_path = Path(str(release.get("receipt_path", "")))
            if receipt_path.is_file():
                if sha256_file(receipt_path) != str(release.get("receipt_sha256", "")):
                    failures.append("Supplementary Figure 3 receipt SHA-256 changed")
            else:
                failures.append(f"Supplementary Figure 3 receipt disappeared: {receipt_path}")
        except Exception as exc:
            failures.append(f"Cannot parse release_status.json: {exc}")

    if args.require_publication_ready and not publication_ready:
        failures.append(
            "Release gate blocked: STRICT_OFFICIAL_PORT_REVIEW_ONLY; publication_ready=false"
        )

    required_suffix = ("execution", "luad-official-code-port", "final")
    if args.strict and tuple(output_dir.parts[-len(required_suffix):]) != required_suffix:
        failures.append(
            "Strict publication output must end in "
            "execution/luad-official-code-port/final; observed " + str(output_dir)
        )

    actual_names = {path.name for path in output_dir.glob("*.pdf")}
    if actual_names != set(EXPECTED_FINALS):
        failures.append(
            f"Final PDF set mismatch: expected {sorted(EXPECTED_FINALS)}, found {sorted(actual_names)}"
        )

    for name in EXPECTED_FINALS:
        path = output_dir / name
        if not path.is_file():
            failures.append(f"Missing final PDF: {path}")
            continue
        try:
            with path.open("rb") as handle:
                header = handle.read(5)
            if header != b"%PDF-":
                failures.append(f"Invalid PDF header: {path}")
            reader = PdfReader(str(path), strict=False)
            pages = len(reader.pages)
            if pages != expected[name]:
                failures.append(f"{name} page count {pages}, expected {expected[name]}")
            content_sizes = [page_content_bytes(page) for page in reader.pages]
            empty_pages = [index + 1 for index, size in enumerate(content_sizes) if size < 50]
            if empty_pages:
                failures.append(f"{name} has empty content streams on pages {empty_pages}")
            texts = [(page.extract_text() or "") for page in reader.pages]
            page_texts[name] = texts
            expected_ready_text = f"publication_ready={str(publication_ready).lower()}"
            for page_index, page_text in enumerate(texts, start=1):
                if name != ATLAS_FILENAME and "official atoms; layout only" not in page_text:
                    failures.append(
                        f"{name} page {page_index} is missing the assembly provenance footer"
                    )
                if release_status != "UNKNOWN" and release_status not in page_text:
                    failures.append(
                        f"{name} page {page_index} is missing visible release status: "
                        f"{release_status}"
                    )
                if expected_ready_text not in page_text:
                    failures.append(
                        f"{name} page {page_index} is missing visible release flag: "
                        f"{expected_ready_text}"
                    )
                if f"{page_index}/{pages}" not in page_text:
                    failures.append(
                        f"{name} page {page_index} is missing the expected page counter "
                        f"{page_index}/{pages}"
                    )
            extracted = "\n".join(texts)
            pdf_rows.append(
                {
                    "filename": name,
                    "path": str(path),
                    "pages": pages,
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "min_content_stream_bytes": min(content_sizes) if content_sizes else 0,
                    "text_characters": len(extracted),
                }
            )
        except Exception as exc:
            failures.append(f"Cannot reopen {name}: {exc}")

    # Boundary checks are page- and panel-specific.  A token appearing once
    # elsewhere in a multi-page PDF cannot satisfy the intended placeholder.
    for (name, page_index), expected_tokens in configured_boundaries.items():
        texts = page_texts.get(name, [])
        if page_index > len(texts):
            continue
        observed_text = texts[page_index - 1]
        observed_tokens = {token for token in BOUNDARY_TOKENS if token in observed_text}
        if observed_tokens != expected_tokens:
            failures.append(
                f"{name} page {page_index} boundary tokens {sorted(observed_tokens)}, "
                f"expected {sorted(expected_tokens)}"
            )

    atlas_texts = page_texts.get(ATLAS_FILENAME, [])
    for page_index, panel in enumerate(atlas_panels, start=1):
        if page_index > len(atlas_texts):
            continue
        observed_text = atlas_texts[page_index - 1]
        if panel not in observed_text:
            failures.append(
                f"{ATLAS_FILENAME} page {page_index} is missing frozen panel name {panel}"
            )
        for required_heading in (
            "TNBC reference panel",
            "LUAD official-code atom",
            "TNBC biological role",
            "LUAD biological role",
            "Interpretation boundary",
            "Visual comparison",
        ):
            if required_heading not in observed_text:
                failures.append(
                    f"{ATLAS_FILENAME} page {page_index} ({panel}) is missing "
                    f"heading {required_heading}"
                )
        expected_tokens = atlas_boundaries.get(panel, set())
        narrative_counts = explanation_token_counts.get(panel, {})
        for token in BOUNDARY_TOKENS:
            observed_count = observed_text.count(token)
            narrative_count = int(narrative_counts.get(token, 0))
            if token in expected_tokens:
                if observed_count < narrative_count + 1:
                    failures.append(
                        f"{ATLAS_FILENAME} page {page_index} ({panel}) lacks the "
                        f"panel-specific {token} boundary in addition to atlas prose"
                    )
            elif observed_count != narrative_count:
                failures.append(
                    f"{ATLAS_FILENAME} page {page_index} ({panel}) has unexpected "
                    f"{token} text outside the frozen atlas explanation"
                )

    final_manifest_path = audit_dir / "final_output_sha256.tsv"
    if not final_manifest_path.is_file():
        failures.append("final_output_sha256.tsv was not generated")
    else:
        final_manifest_rows = read_tsv(final_manifest_path)
        final_manifest_by_name = {
            row.get("filename", ""): row for row in final_manifest_rows
        }
        if (
            len(final_manifest_rows) != len(EXPECTED_FINALS)
            or len(final_manifest_by_name) != len(EXPECTED_FINALS)
            or set(final_manifest_by_name) != set(EXPECTED_FINALS)
        ):
            failures.append(
                "final_output_sha256.tsv must contain exactly one row for each of the seven PDFs"
            )
        for name in EXPECTED_FINALS:
            row = final_manifest_by_name.get(name)
            path = output_dir / name
            if not row or not path.is_file():
                continue
            try:
                if not same_path(row.get("output", ""), str(path)):
                    failures.append(f"Final-output manifest path mismatch for {name}")
                if int(row.get("pages", "-1")) != expected[name]:
                    failures.append(f"Final-output manifest page-count mismatch for {name}")
                if int(row.get("bytes", "-1")) != path.stat().st_size:
                    failures.append(f"Final-output manifest byte-count mismatch for {name}")
                if row.get("sha256", "").lower() != sha256_file(path):
                    failures.append(f"Final-output manifest SHA-256 mismatch for {name}")
            except Exception as exc:
                failures.append(f"Cannot validate final-output manifest row for {name}: {exc}")

    input_contract_path = audit_dir / "input_contract.tsv"
    if not input_contract_path.is_file():
        failures.append("input_contract.tsv was not generated")
    else:
        contract_rows = read_tsv(input_contract_path)
        contract_by_atom = {row.get("atom_id", ""): row for row in contract_rows}
        if (
            len(contract_rows) != EXPECTED_ATOM_COUNT
            or len(contract_by_atom) != EXPECTED_ATOM_COUNT
            or set(contract_by_atom) != set(spec["atoms"])
        ):
            failures.append(
                "input_contract.tsv does not contain the exact frozen 75-atom contract"
            )
        expected_usages: dict[str, list[str]] = {
            atom_id: [] for atom_id in spec["atoms"]
        }
        for output in spec["outputs"]:
            for page_index, page in enumerate(output["pages"], start=1):
                for item in page["items"]:
                    for atom_id in item["atoms"]:
                        expected_usages[str(atom_id)].append(
                            f"{output['filename']}#page={page_index}:{item['panel']}"
                        )
        for atom_id, atom in spec["atoms"].items():
            row = contract_by_atom.get(atom_id)
            if not row:
                continue
            expected_fields = {
                "declared_status": str(atom.get("status", "OFFICIAL_PDF_REQUIRED")),
                "required": "YES" if atom.get("required", False) else "NO",
                "root_key": str(atom.get("root", "")),
                "patterns": ";".join(str(value) for value in atom.get("patterns", [])),
                "match_mode": str(atom.get("match", "one")),
                "publication_boundary": str(atom.get("boundary", "")),
                "vector_pdf_only": "YES" if not atom.get("status") else "N/A",
                "used_by": ";".join(expected_usages[atom_id]),
            }
            for field, expected_value in expected_fields.items():
                if row.get(field, "") != expected_value:
                    failures.append(
                        f"input_contract.tsv mismatch for {atom_id}/{field}: "
                        f"{row.get(field, '')!r} vs {expected_value!r}"
                    )

    declared_path = audit_dir / "declared_boundaries.tsv"
    expected_declared = {
        atom_id: str(atom["status"])
        for atom_id, atom in spec["atoms"].items()
        if atom.get("status") in {"MANUAL_REQUIRED", "UNSUPPORTED"}
    }
    expected_declared.update(FROZEN_SUPPORTED_ZERO)
    if not declared_path.is_file():
        failures.append("declared_boundaries.tsv was not generated")
    else:
        declared_rows = read_tsv(declared_path)
        declared_by_atom = {row.get("atom_id", ""): row for row in declared_rows}
        if (
            len(declared_rows) != len(expected_declared)
            or len(declared_by_atom) != len(expected_declared)
            or set(declared_by_atom) != set(expected_declared)
        ):
            failures.append(
                "declared_boundaries.tsv must contain exactly 5 MANUAL_REQUIRED, "
                "3 UNSUPPORTED, and 2 frozen SUPPORTED_ZERO atoms"
            )
        for atom_id, status in expected_declared.items():
            row = declared_by_atom.get(atom_id)
            if not row:
                continue
            if row.get("status") != status or row.get("is_missing") != "NO":
                failures.append(
                    f"Declared boundary mismatch for {atom_id}: "
                    f"status={row.get('status')}, is_missing={row.get('is_missing')}"
                )
            if not row.get("note", "").strip():
                failures.append(f"Declared boundary note is empty for {atom_id}")

    resolved_path = audit_dir / "resolved_atoms.tsv"
    if not resolved_path.is_file():
        failures.append("resolved_atoms.tsv was not generated")
        resolved_rows: list[dict[str, str]] = []
    else:
        resolved_rows = read_tsv(resolved_path)
        resolved_by_atom: dict[str, list[dict[str, str]]] = {}
        resolution_keys: set[tuple[str, str, str, str]] = set()
        for row in resolved_rows:
            atom_id = row.get("atom_id", "")
            resolved_by_atom.setdefault(atom_id, []).append(row)
            key = (
                atom_id,
                row.get("path", ""),
                row.get("page_1based", ""),
                row.get("crop_top", ""),
            )
            if key in resolution_keys:
                failures.append(f"Duplicate resolved atom/source slice: {key}")
            resolution_keys.add(key)
            if row.get("status") != "RESOLVED_OFFICIAL_PDF":
                failures.append(f"Unexpected resolved-atom status for {atom_id}: {row.get('status')}")
            source = Path(row.get("path", ""))
            if not source.is_file():
                failures.append(f"Resolved atomic source disappeared: {source}")
            else:
                if sha256_file(source) != row.get("sha256"):
                    failures.append(f"Resolved atomic source SHA-256 changed: {source}")
                try:
                    reader = PdfReader(str(source), strict=False)
                    page_1based = int(row.get("page_1based", "0"))
                    if int(row.get("page_count", "-1")) != len(reader.pages):
                        failures.append(f"Resolved atomic page-count changed: {source}")
                    if page_1based < 1 or page_1based > len(reader.pages):
                        failures.append(f"Resolved atomic page is out of range: {source}")
                except Exception as exc:
                    failures.append(f"Cannot reopen resolved atomic source {source}: {exc}")
        for atom_id, atom in spec["atoms"].items():
            if atom_id in expected_declared:
                if atom_id in resolved_by_atom:
                    failures.append(f"Boundary atom unexpectedly resolved to a PDF: {atom_id}")
                continue
            rows = resolved_by_atom.get(atom_id, [])
            if not rows:
                failures.append(f"Official PDF atom has no resolved source slice: {atom_id}")
            elif atom.get("match", "one") == "one" and len(rows) != 1:
                failures.append(
                    f"match:one atom {atom_id} has {len(rows)} resolved source slices"
                )

    static_map_path = args.spec.resolve().parent / "tnbc_panel_map.tsv"
    atlas_crop_path = audit_dir / "atlas_reference_crops.tsv"
    if not static_map_path.is_file():
        failures.append(f"Frozen TNBC panel map is absent: {static_map_path}")
        static_map_rows: list[dict[str, str]] = []
    else:
        static_map_rows = read_tsv(static_map_path)
    static_map_by_panel = {row.get("panel", ""): row for row in static_map_rows}
    if (
        len(static_map_rows) != EXPECTED_PANEL_COUNT
        or len(static_map_by_panel) != EXPECTED_PANEL_COUNT
        or set(static_map_by_panel) != set(atlas_panels)
    ):
        failures.append("tnbc_panel_map.tsv does not cover the exact frozen 73-panel atlas")
    for panel, row in static_map_by_panel.items():
        try:
            crop = tuple(
                float(row[field])
                for field in ("x0_top", "y0_top", "x1_top", "y1_top")
            )
            if not (
                0.0 <= crop[0] < crop[2] <= 1.0
                and 0.0 <= crop[1] < crop[3] <= 1.0
            ):
                failures.append(f"Invalid normalized TNBC crop for {panel}: {crop}")
            if int(row.get("page_1based", "0")) < 1:
                failures.append(f"Invalid TNBC source page for {panel}")
            if row.get("source_key") not in {"MAIN_ARTICLE", "SUPPLEMENT_PDF"}:
                failures.append(f"Invalid TNBC source key for {panel}: {row.get('source_key')}")
        except Exception as exc:
            failures.append(f"Cannot parse frozen TNBC crop for {panel}: {exc}")

    if not atlas_crop_path.is_file():
        failures.append("atlas_reference_crops.tsv was not generated")
        atlas_crop_rows: list[dict[str, str]] = []
    else:
        atlas_crop_rows = read_tsv(atlas_crop_path)
    atlas_crop_by_panel = {row.get("panel", ""): row for row in atlas_crop_rows}
    if (
        len(atlas_crop_rows) != EXPECTED_PANEL_COUNT
        or len(atlas_crop_by_panel) != EXPECTED_PANEL_COUNT
        or set(atlas_crop_by_panel) != set(atlas_panels)
    ):
        failures.append("atlas_reference_crops.tsv does not contain 73 unique frozen panels")
    for panel, expected_row in static_map_by_panel.items():
        observed_row = atlas_crop_by_panel.get(panel)
        if not observed_row:
            continue
        if observed_row.get("source_key") != expected_row.get("source_key"):
            failures.append(f"TNBC atlas source-key mismatch for {panel}")
        if observed_row.get("source_page_1based") != expected_row.get("page_1based"):
            failures.append(f"TNBC atlas source-page mismatch for {panel}")
        try:
            expected_crop = tuple(
                float(expected_row[field])
                for field in ("x0_top", "y0_top", "x1_top", "y1_top")
            )
            observed_crop = parse_crop(observed_row.get("crop_top", ""))
            if len(observed_crop) != 4 or any(
                abs(left - right) > 1e-9
                for left, right in zip(expected_crop, observed_crop)
            ):
                failures.append(f"TNBC atlas crop-coordinate mismatch for {panel}")
        except Exception as exc:
            failures.append(f"Cannot validate TNBC atlas crop for {panel}: {exc}")
        source = Path(observed_row.get("source", ""))
        if not source.is_file():
            failures.append(f"TNBC atlas reference source disappeared for {panel}: {source}")
        expected_format = source.suffix.lower().lstrip(".")
        if observed_row.get("source_format") != expected_format:
            failures.append(f"TNBC atlas source-format mismatch for {panel}")

    explanations_path = args.spec.resolve().parent / "panel_explanations.tsv"
    if not explanations_path.is_file():
        failures.append(f"Panel explanations are absent: {explanations_path}")
    else:
        explanation_rows = read_tsv(explanations_path)
        explanation_by_panel = {row.get("panel", ""): row for row in explanation_rows}
        if (
            len(explanation_rows) != EXPECTED_PANEL_COUNT
            or len(explanation_by_panel) != EXPECTED_PANEL_COUNT
            or set(explanation_by_panel) != set(atlas_panels)
        ):
            failures.append("panel_explanations.tsv does not cover the exact 73-panel atlas")
        for panel, row in explanation_by_panel.items():
            for field in (
                "biological_role_tnbc",
                "biological_role_luad",
                "interpretation_boundary",
                "visual_comparison",
            ):
                if not row.get(field, "").strip():
                    failures.append(f"Panel explanation is empty for {panel}/{field}")

    missing_path = audit_dir / "missing_inputs.tsv"
    if not missing_path.is_file():
        failures.append("missing_inputs.tsv was not generated")
        missing_required: list[dict[str, str]] = []
    else:
        missing_rows = read_tsv(missing_path)
        missing_required = [row for row in missing_rows if row.get("is_required") == "YES"]
        if args.strict and missing_required:
            failures.append(
                f"Strict input audit has {len(missing_required)} required missing/ambiguous entries"
            )

    embedding_path = audit_dir / "embedding_manifest.tsv"
    if not embedding_path.is_file():
        failures.append("embedding_manifest.tsv was not generated")
        embedding_rows: list[dict[str, str]] = []
    else:
        embedding_rows = read_tsv(embedding_path)
        luad_rows = [row for row in embedding_rows if row.get("source_role") == "LUAD_OFFICIAL_ATOM"]
        tnbc_rows = [row for row in embedding_rows if row.get("source_role") == "TNBC_REFERENCE_CROP"]
        unexpected_roles = sorted(
            {
                row.get("source_role", "")
                for row in embedding_rows
                if row.get("source_role") not in {"LUAD_OFFICIAL_ATOM", "TNBC_REFERENCE_CROP"}
            }
        )
        if unexpected_roles:
            failures.append(f"Unexpected embedding source roles: {unexpected_roles}")
        if not luad_rows:
            failures.append("No LUAD official atomic PDF embeddings were recorded")
        for row in luad_rows:
            if row.get("embedding_method") != "pypdf.merge_transformed_page_vector":
                failures.append(
                    f"Non-vector LUAD atom embedding recorded: {row.get('source')} "
                    f"({row.get('embedding_method')})"
                )
            source = Path(row["source"])
            if not source.is_file():
                failures.append(f"Embedded source disappeared: {source}")
            elif sha256_file(source) != row.get("source_sha256"):
                failures.append(f"Embedded source changed after assembly: {source}")
        atom_consumers: dict[str, tuple[str, int, str]] = {}
        for output in spec["outputs"]:
            for page_index, page in enumerate(output["pages"], start=1):
                for item in page["items"]:
                    for atom_id in item["atoms"]:
                        atom_consumers[str(atom_id)] = (
                            str(output["filename"]),
                            page_index,
                            str(item["panel"]),
                        )
        atlas_page_by_panel = {
            panel: page_index for page_index, panel in enumerate(atlas_panels, start=1)
        }

        def matches_resolution(
            embedding: dict[str, str], resolution: dict[str, str]
        ) -> bool:
            if embedding.get("atom_id") != resolution.get("atom_id"):
                return False
            if not same_path(embedding.get("source", ""), resolution.get("path", "")):
                return False
            if embedding.get("source_page_1based") != resolution.get("page_1based"):
                return False
            try:
                left = parse_crop(embedding.get("source_crop_top", ""))
                right = parse_crop(resolution.get("crop_top", ""))
                return len(left) == len(right) and all(
                    abs(a - b) <= 1e-9 for a, b in zip(left, right)
                )
            except Exception:
                return False

        if len(luad_rows) != 2 * len(resolved_rows):
            failures.append(
                f"LUAD embedding count {len(luad_rows)}, expected exactly twice the "
                f"{len(resolved_rows)} resolved source slices"
            )
        for resolution in resolved_rows:
            matching = [row for row in luad_rows if matches_resolution(row, resolution)]
            atom_id = resolution.get("atom_id", "")
            if len(matching) != 2:
                failures.append(
                    f"Resolved source slice {atom_id}/{resolution.get('path')} has "
                    f"{len(matching)} embeddings, expected configured-output plus atlas"
                )
                continue
            configured = [row for row in matching if row.get("output") != ATLAS_FILENAME]
            atlas = [row for row in matching if row.get("output") == ATLAS_FILENAME]
            if len(configured) != 1 or len(atlas) != 1:
                failures.append(
                    f"Resolved source slice {atom_id} is not embedded exactly once in "
                    "a configured output and once in the atlas"
                )
                continue
            expected_consumer = atom_consumers.get(atom_id)
            if expected_consumer:
                expected_output, expected_page, expected_panel = expected_consumer
                row = configured[0]
                if (
                    row.get("output") != expected_output
                    or row.get("output_page") != str(expected_page)
                    or row.get("panel") != expected_panel
                ):
                    failures.append(f"Configured-output embedding consumer mismatch for {atom_id}")
                row = atlas[0]
                expected_atlas_page = atlas_page_by_panel[expected_panel]
                if (
                    row.get("output_page") != str(expected_atlas_page)
                    or row.get("panel") != expected_panel
                ):
                    failures.append(f"Atlas embedding consumer mismatch for {atom_id}")

        if len(tnbc_rows) != EXPECTED_PANEL_COUNT:
            failures.append(
                f"TNBC reference embedding count {len(tnbc_rows)}, expected {EXPECTED_PANEL_COUNT}"
            )
        tnbc_by_panel: dict[str, list[dict[str, str]]] = {}
        for row in tnbc_rows:
            tnbc_by_panel.setdefault(row.get("panel", ""), []).append(row)
            method = row.get("embedding_method")
            if method not in {
                "pypdf.merge_transformed_page_vector",
                "reference_png_crop_no_luad_atom_rasterization",
            }:
                failures.append(
                    f"Invalid TNBC reference embedding method for {row.get('panel')}: {method}"
                )
            source = Path(row.get("source", ""))
            if not source.is_file():
                failures.append(f"TNBC reference source disappeared: {source}")
            elif sha256_file(source) != row.get("source_sha256"):
                failures.append(f"TNBC reference source changed after assembly: {source}")
        if set(tnbc_by_panel) != set(atlas_panels) or any(
            len(rows) != 1 for rows in tnbc_by_panel.values()
        ):
            failures.append("TNBC embeddings must contain exactly one row for every atlas panel")
        for panel, page_index in atlas_page_by_panel.items():
            rows = tnbc_by_panel.get(panel, [])
            if not rows:
                continue
            row = rows[0]
            if row.get("output") != ATLAS_FILENAME or row.get("output_page") != str(page_index):
                failures.append(f"TNBC atlas embedding page mismatch for {panel}")

    identity_path = audit_dir / "source_identity_manifest.tsv"
    identity_summary_path = audit_dir / "source_identity_summary.json"
    if not identity_path.is_file() or not identity_summary_path.is_file():
        failures.append("Official source identity manifests were not generated")
    else:
        identity_rows = read_tsv(identity_path)
        try:
            identity_summary = json.loads(identity_summary_path.read_text(encoding="utf-8"))
        except Exception as exc:
            failures.append(f"Cannot parse source_identity_summary.json: {exc}")
            identity_summary = {}
        if identity_summary.get("identity_ready") is not True:
            failures.append("Official source identity is not ready for all four roots")
        if set(identity_summary.get("roots", {})) != {"fig12", "fig3", "fig45", "supp3"}:
            failures.append("Official source identity summary does not contain the exact four roots")
        if len(identity_rows) != len(resolved_rows):
            failures.append(
                f"Official source identity has {len(identity_rows)} rows for "
                f"{len(resolved_rows)} resolved PDF slices"
            )
        for root_key, root_identity in identity_summary.get("roots", {}).items():
            if root_identity.get("identity_ready") is not True:
                failures.append(f"Official source identity is not ready for root {root_key}")
            if root_identity.get("resolved_pdf_slices") != root_identity.get("verified_pdf_slices"):
                failures.append(f"Official source identity slice counts differ for root {root_key}")
            for receipt in root_identity.get("receipts", []):
                receipt_path = Path(str(receipt.get("path", "")))
                expected_receipt_sha = str(receipt.get("sha256", ""))
                if not receipt_path.is_file():
                    failures.append(f"Identity receipt disappeared for {root_key}: {receipt_path}")
                elif sha256_file(receipt_path) != expected_receipt_sha:
                    failures.append(f"Identity receipt/hash manifest changed for {root_key}: {receipt_path}")
        for row in identity_rows:
            if row.get("identity_status") != "VERIFIED":
                failures.append(
                    f"Official source identity failed: {row.get('root_key')}/"
                    f"{row.get('atom_id')}: {row.get('detail')}"
                )
                continue
            source = Path(row["source"])
            if not source.is_file():
                failures.append(f"Identity-bound source disappeared: {source}")
            elif sha256_file(source) != row.get("sha256_observed"):
                failures.append(f"Identity-bound source hash changed: {source}")

    summary_path = audit_dir / "assembly_summary.json"
    if not summary_path.is_file():
        failures.append("assembly_summary.json was not generated")
    else:
        try:
            assembly_summary = json.loads(summary_path.read_text(encoding="utf-8"))
            if assembly_summary.get("run_mode") != "ASSEMBLY":
                failures.append(
                    f"assembly_summary.json run_mode={assembly_summary.get('run_mode')}, "
                    "expected ASSEMBLY"
                )
            if assembly_summary.get("release_status") != release_status:
                failures.append("assembly_summary.json release status does not match release_status.json")
            if assembly_summary.get("publication_ready") is not publication_ready:
                failures.append("assembly_summary.json publication-ready flag is inconsistent")
            if assembly_summary.get("source_identity_ready") is not True:
                failures.append("assembly_summary.json source identity is not ready")
            if assembly_summary.get("expected_outputs") != list(EXPECTED_FINALS):
                failures.append("assembly_summary.json expected output contract changed")
            if set(assembly_summary.get("created_outputs", [])) != set(EXPECTED_FINALS):
                failures.append("assembly_summary.json does not record exactly seven created PDFs")
            if assembly_summary.get("resolved_atom_pdf_count") != len(resolved_rows):
                failures.append("assembly_summary.json resolved-atom count is inconsistent")
            if assembly_summary.get("declared_boundary_count") != len(expected_declared):
                failures.append("assembly_summary.json declared-boundary count is inconsistent")
            if assembly_summary.get("missing_required_count") != 0:
                failures.append("assembly_summary.json reports required missing inputs")
            if assembly_summary.get("embedding_count") != len(embedding_rows):
                failures.append("assembly_summary.json embedding count is inconsistent")
            if assembly_summary.get("luad_embedding_methods") != [
                "pypdf.merge_transformed_page_vector"
            ]:
                failures.append("assembly_summary.json records a non-vector LUAD embedding method")
        except Exception as exc:
            failures.append(f"Cannot validate assembly_summary.json: {exc}")

    try:
        pdftoppm = locate_pdftoppm(args.pdftoppm)
    except Exception as exc:
        failures.append(str(exc))
        pdftoppm = ""
    if pdftoppm:
        rendered_root = qa_dir / "rendered"
        sheets_root = qa_dir / "contact_sheets"
        for name in EXPECTED_FINALS:
            path = output_dir / name
            if not path.is_file():
                continue
            render_dir = rendered_root / path.stem
            if render_dir.exists():
                for old in render_dir.glob("page-*.png"):
                    old.unlink()
            try:
                pages = render_pdf(pdftoppm, path, render_dir, args.dpi)
                if len(pages) != expected[name]:
                    failures.append(
                        f"Rendered page count mismatch for {name}: {len(pages)} vs {expected[name]}"
                    )
                for page_index, png in enumerate(pages, start=1):
                    metrics = image_metrics(png)
                    metrics.update({"filename": name, "page": page_index})
                    render_rows.append(metrics)
                    if metrics["width_px"] < 300 or metrics["height_px"] < 300:
                        failures.append(f"Rendered page is too small: {name} page {page_index}")
                    if metrics["ink_fraction"] < 0.0002:
                        failures.append(f"Rendered page appears blank: {name} page {page_index}")
                    if metrics["dark_fraction"] > 0.90:
                        failures.append(f"Rendered page appears black: {name} page {page_index}")
                contact_sheet(pages, sheets_root / f"{path.stem}.png", name)
            except Exception as exc:
                failures.append(f"Render failure for {name}: {exc}")

    write_tsv(
        qa_dir / "pdf_integrity.tsv",
        pdf_rows,
        ("filename", "path", "pages", "bytes", "sha256", "min_content_stream_bytes", "text_characters"),
    )
    write_tsv(
        qa_dir / "render_metrics.tsv",
        render_rows,
        (
            "filename", "page", "png", "width_px", "height_px", "white_fraction",
            "dark_fraction", "ink_fraction",
        ),
    )
    summary = {
        "status": "FAIL" if failures else release_status,
        "technical_checks": "FAILED" if failures else "VALIDATED",
        "release_status": release_status,
        "publication_ready": publication_ready,
        "release_gate_required": bool(args.require_publication_ready),
        "strict": bool(args.strict),
        "expected_output_count": len(EXPECTED_FINALS),
        "verified_output_count": len(pdf_rows),
        "rendered_page_count": len(render_rows),
        "required_missing_count": len(missing_required),
        "failure_count": len(failures),
        "failures": failures,
    }
    (qa_dir / "qa_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if failures:
        print("ASSEMBLY QA FAILED", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 2
    print(
        f"ASSEMBLY TECHNICAL QA VALIDATED: {len(pdf_rows)} PDFs, "
        f"{len(render_rows)} rendered pages; release_status={release_status}; "
        f"publication_ready={str(publication_ready).lower()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
