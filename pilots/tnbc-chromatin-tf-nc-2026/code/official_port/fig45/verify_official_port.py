#!/usr/bin/env python3
"""Static/provenance verification for the strict official Figure 4/5 port."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_replacements(spec: str) -> list[tuple[str, str]]:
    if not spec:
        return []
    answer = []
    for item in spec.split(";"):
        old, new = item.split("=>", 1)
        answer.append((old, new))
    return answer


def check(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def pdf_media_box_points(path: Path) -> tuple[float, float, float, float] | None:
    """Read the first explicit PDF MediaBox; Cairo emits it uncompressed."""
    match = re.search(
        rb"/MediaBox\s*\[\s*([-+0-9.]+)\s+([-+0-9.]+)\s+"
        rb"([-+0-9.]+)\s+([-+0-9.]+)\s*\]",
        path.read_bytes(),
    )
    if not match:
        return None
    return tuple(float(value) for value in match.groups())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--code-dir", required=True, type=Path)
    parser.add_argument("--generated-dir", required=True, type=Path)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--rscript", type=Path)
    args = parser.parse_args()

    code_dir = args.code_dir.resolve()
    generated = args.generated_dir.resolve()
    rows = read_tsv(code_dir / "source_manifest.tsv")
    failures: list[str] = []
    upstream_license = code_dir / "UPSTREAM_LICENSE.txt"
    check(upstream_license.is_file(), "missing upstream MIT licence copy", failures)
    if upstream_license.is_file():
        check(
            sha256(upstream_license)
            == "470349c2333c9e595d08f10fab1bfb18a29c8c909a502baac93e12b1fe13488c",
            "upstream MIT licence checksum mismatch",
            failures,
        )
    block_receipt_path = generated / "materialization_receipt.json"
    check(block_receipt_path.is_file(), "missing materialization_receipt.json", failures)
    receipt = json.loads(block_receipt_path.read_text(encoding="utf-8"))
    check(receipt.get("commit") == "edf5314", "unexpected upstream commit", failures)
    check(receipt.get("capsule") == "7227095/v1", "unexpected upstream capsule", failures)
    receipt_by_id = {x["block_id"]: x for x in receipt.get("blocks", [])}
    rows_by_id = {row["block_id"]: row for row in rows}
    check(len(receipt_by_id) == len(rows), "receipt/manifest block count mismatch", failures)
    check(len(rows_by_id) == len(rows), "duplicate block_id in source manifest", failures)

    for row in rows:
        block_id = row["block_id"]
        exact = generated / "exact" / f"{block_id}.R"
        patched = generated / "patched" / f"{block_id}.R"
        diff = generated / "diff" / f"{block_id}.diff"
        check(exact.is_file(), f"{block_id}: missing exact block", failures)
        check(patched.is_file(), f"{block_id}: missing patched block", failures)
        check(diff.is_file(), f"{block_id}: missing unified diff", failures)
        if not exact.is_file() or not patched.is_file():
            continue
        expected_patched = exact.read_text(encoding="utf-8")
        for old, new in parse_replacements(row["allowed_replacements"]):
            check(old in expected_patched, f"{block_id}: declared source literal absent", failures)
            expected_patched = expected_patched.replace(old, new)
        observed_patched = patched.read_text(encoding="utf-8")
        check(
            observed_patched == expected_patched,
            f"{block_id}: patched code has a non-allowlisted edit",
            failures,
        )
        rec = receipt_by_id.get(block_id, {})
        check(
            rec.get("source_sha256_observed", "").lower() == row["sha256"].lower(),
            f"{block_id}: upstream checksum not preserved",
            failures,
        )
        check(
            rec.get("exact_block_sha256") == sha256(exact),
            f"{block_id}: exact block checksum mismatch",
            failures,
        )
        check(
            rec.get("patched_block_sha256") == sha256(patched),
            f"{block_id}: patched block checksum mismatch",
            failures,
        )

    # P0 fidelity locks from the independent audit.
    heat_row = rows_by_id.get("supp5b_heatmap", {})
    heat_exact_path = generated / "exact" / "supp5b_heatmap.R"
    heat_patched_path = generated / "patched" / "supp5b_heatmap.R"
    heat_exact = heat_exact_path.read_text(encoding="utf-8") if heat_exact_path.is_file() else ""
    heat_patched = heat_patched_path.read_text(encoding="utf-8") if heat_patched_path.is_file() else ""
    check(heat_row.get("mode") == "verbatim", "supp5b_heatmap is not manifest-locked as verbatim", failures)
    check(not heat_row.get("allowed_replacements"), "supp5b_heatmap declares a replacement", failures)
    check(heat_exact == heat_patched, "supp5b_heatmap exact/patched blocks differ", failures)
    check("Liberation Sans" not in heat_exact, "supp5b_heatmap contains a non-official font", failures)
    check(heat_exact.count('fontfamily = "Helvetica"') >= 5,
          "supp5b_heatmap Helvetica declarations changed", failures)
    check(
        'grid.text("* FDR<0.05  ** FDR<0.01  *** FDR<0.001"' in heat_exact,
        "supp5b_heatmap is missing the official FDR-star footer",
        failures,
    )
    check(
        "cairo_pdf(out_path, width = pdf_width, height = pdf_height)" in heat_exact,
        "supp5b_heatmap is missing official dynamic device sizing",
        failures,
    )
    font_substitution_rows = [
        row["block_id"] for row in rows
        if "font" in row["mode"].lower() or "Helvetica=>" in row["allowed_replacements"]
    ]
    check(not font_substitution_rows, f"font substitutions remain: {font_substitution_rows}", failures)

    for block_id in ("supp6d_uni_tcga_rfs", "supp6h_multi_tcga_rfs"):
        patched_path = generated / "patched" / f"{block_id}.R"
        text = patched_path.read_text(encoding="utf-8") if patched_path.is_file() else ""
        check('label = "DFS\\nTCGA_LUAD"' in text, f"{block_id}: visible DFS label missing", failures)
        check("TCGA_LUAD_RFS_" in text, f"{block_id}: frozen RFS input path/object was altered", failures)
        check("TCGA_LUAD_DFS_" not in text, f"{block_id}: DFS leaked beyond visible label", failures)

    strict_case_ids = {f"supp7b_case{i}" for i in range(1, 5)}
    check(strict_case_ids.issubset(rows_by_id), "four official Supplementary Figure 7B case blocks are not frozen", failures)
    for block_id in strict_case_ids:
        case_path = generated / "patched" / f"{block_id}.R"
        text = case_path.read_text(encoding="utf-8") if case_path.is_file() else ""
        check(text.count("plot_aac_vs_nes_with_ci(") == 1, f"{block_id}: official plot call missing/duplicated", failures)
        check(text.count("cairo_pdf(") == 1, f"{block_id}: official device call missing/duplicated", failures)
        check("width=" not in text and "width =" not in text and "height=" not in text and "height =" not in text,
              f"{block_id}: official default 7x7 device was overridden", failures)

    # The port must remain isolated from the withdrawn custom renderer.
    forbidden_imports = (
        "render_anchor_" + "figure45_supp",
        "publication_rebuild/" + "tnbc_anchor_theme",
    )
    for path in code_dir.glob("*"):
        if path.is_file() and path.suffix.lower() in {".r", ".py", ".sh"}:
            text = path.read_text(encoding="utf-8")
            check(
                all(marker not in text for marker in forbidden_imports),
                f"{path.name}: imports withdrawn custom renderer",
                failures,
            )

    renderer_text = (code_dir / "render_official_fig45.R").read_text(encoding="utf-8")
    check("Liberation Sans" not in renderer_text, "renderer still overrides the official Helvetica font", failures)
    for title in (
        "Multivariate GSE41271-LUAD OS", "Multivariate GSE41271-LUAD RFS",
        "Multivariate TCGA-LUAD OS", "Multivariate TCGA-LUAD DFS",
    ):
        check(title in renderer_text, f"Figure 4 title lock missing: {title}", failures)
    for snippet in (
        'min_cohorts = 1L', 'col_width_mm = 10', 'row_height_mm = 4',
        'plot_title = "HC-TR activity correlation"',
        'source_block("supp5b_heatmap", heat_env, patched = FALSE)',
    ):
        check(snippet in renderer_text, f"Supplementary Figure 5B wrapper lock missing: {snippet}", failures)
    for block_id in ("fig5g_pdx", "fig5h_waterfall", "supp7d_pdx"):
        call_literal = "source_block(" + chr(34) + block_id + chr(34)
        check(
            call_literal not in renderer_text,
            f"unsupported block is executable in renderer: {block_id}",
            failures,
        )

    r_parse = []
    parse_not_applicable_blocks: set[str] = set()
    if args.rscript:
        rscript = args.rscript.resolve()
        check(rscript.is_file(), f"Rscript not found: {rscript}", failures)
        if rscript.is_file():
            parse_paths = list(sorted(code_dir.glob("*.R")))
            parse_paths += list(sorted((generated / "exact").glob("*.R")))
            parse_paths += list(sorted((generated / "patched").glob("*.R")))
            for path in parse_paths:
                proc = subprocess.run(
                    [str(rscript), "-e", f"parse(file={json.dumps(str(path))})"],
                    text=True,
                    capture_output=True,
                )
                generated_variant = path.parent.name in {"exact", "patched"}
                row = rows_by_id.get(path.stem) if generated_variant else None
                # An unsupported source slice is retained solely for provenance and
                # checksum comparison; it is never sourced by the renderer.  Some
                # official slices (currently Figure 5H) end inside an upstream
                # control-flow scope and therefore are not standalone R programs.
                # We still attempt to parse every slice.  Only a parse failure from
                # a manifest-locked, non-executed unsupported block is N/A; all
                # supported blocks and all complete unsupported slices remain strict.
                if proc.returncode != 0 and row and row.get("mode") == "unsupported":
                    parse_not_applicable_blocks.add(path.stem)
                    r_parse.append({
                        "file": str(path),
                        "returncode": proc.returncode,
                        "status": "R_PARSE_NOT_APPLICABLE_UNSUPPORTED_SOURCE_SLICE",
                        "reason": "provenance-only block is not executed and is not a standalone R program",
                    })
                else:
                    r_parse.append({
                        "file": str(path),
                        "returncode": proc.returncode,
                        "status": "PARSED" if proc.returncode == 0 else "R_PARSE_FAILED",
                    })
                    check(proc.returncode == 0, f"R parse failure: {path.name}: {proc.stderr}", failures)

    result_summary = {}
    if args.results_dir:
        results = args.results_dir.resolve()
        render_manifest_path = results / "render_manifest.tsv"
        unsupported_path = results / "unsupported_panels.tsv"
        executed_path = results / "audit" / "executed_official_blocks.tsv"
        check(render_manifest_path.is_file(), "missing render_manifest.tsv", failures)
        check(unsupported_path.is_file(), "missing unsupported_panels.tsv", failures)
        check(executed_path.is_file(), "missing executed_official_blocks.tsv", failures)
        executed_rows = read_tsv(executed_path) if executed_path.is_file() else []
        executed_ids = [row["block_id"] for row in executed_rows]
        executed_set = set(executed_ids)
        check(
            all(executed_ids.count(block_id) == 1 for block_id in strict_case_ids),
            "each official Supplementary Figure 7B case block must execute exactly once",
            failures,
        )
        check(
            any(row["block_id"] == "supp5b_heatmap" and row["variant"] == "exact" for row in executed_rows),
            "official exact Supplementary Figure 5B heatmap block was not executed",
            failures,
        )
        if render_manifest_path.is_file():
            render_rows = read_tsv(render_manifest_path)
            declared_pdfs: set[Path] = set()
            pdf_hash_rows: list[dict[str, str | int]] = []
            for row in render_rows:
                if row["status"] == "RENDERED_OFFICIAL_BLOCK":
                    output = Path(row["output"]).resolve()
                    declared_pdfs.add(output)
                    check(output.is_file() and output.stat().st_size > 1000, f"missing/empty PDF: {output}", failures)
                    for block_id in row["official_block_ids"].split(";"):
                        check(block_id in executed_set, f"{row['panel']}: claimed block was not executed: {block_id}", failures)
                    if output.is_file():
                        pdf_hash_rows.append({
                            "panel": row["panel"],
                            "filename": output.name,
                            "bytes": output.stat().st_size,
                            "sha256": sha256(output),
                            "official_block_ids": row["official_block_ids"],
                        })
                    if row["panel"].startswith("SupplementaryFigure7B_case") and output.is_file():
                        media_box = pdf_media_box_points(output)
                        check(media_box is not None, f"{row['panel']}: PDF MediaBox missing", failures)
                        if media_box is not None:
                            width_pt = media_box[2] - media_box[0]
                            height_pt = media_box[3] - media_box[1]
                            check(
                                abs(width_pt - 504.0) <= 0.5 and abs(height_pt - 504.0) <= 0.5,
                                f"{row['panel']}: expected official 7x7-inch device, got {width_pt}x{height_pt} pt",
                                failures,
                            )
            observed_pdfs = {
                path.resolve() for path in (results / "atomic").glob("*.pdf")
            }
            check(
                observed_pdfs == declared_pdfs,
                "atomic PDF set differs from render manifest: "
                f"extra={sorted(map(str, observed_pdfs - declared_pdfs))}; "
                f"missing={sorted(map(str, declared_pdfs - observed_pdfs))}",
                failures,
            )
            unsupported = {r["panel"] for r in render_rows if r["status"] == "UNSUPPORTED"}
            check(
                unsupported == {"Figure5G", "Figure5H", "SupplementaryFigure7D"},
                f"unexpected unsupported set: {sorted(unsupported)}",
                failures,
            )
            result_summary["render_rows"] = len(render_rows)
            result_summary["rendered"] = sum(r["status"] == "RENDERED_OFFICIAL_BLOCK" for r in render_rows)
            pdf_hash_path = results / "atomic_pdf_sha256.tsv"
            with pdf_hash_path.open("wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    delimiter="\t",
                    fieldnames=["panel", "filename", "bytes", "sha256", "official_block_ids"],
                )
                writer.writeheader()
                writer.writerows(sorted(pdf_hash_rows, key=lambda row: (str(row["panel"]), str(row["filename"]))))
            result_summary["atomic_pdf_sha256_manifest"] = str(pdf_hash_path)
            result_summary["atomic_pdf_sha256_manifest_sha256"] = sha256(pdf_hash_path)
            result_summary["executed_official_blocks"] = len(executed_rows)
        for forbidden in ("Figure5G.pdf", "Figure5H.pdf", "SupplementaryFigure7D.pdf"):
            check(not (results / "atomic" / forbidden).exists(), f"forbidden substitute exists: {forbidden}", failures)

    report = {
        "status": "PASS" if not failures else "FAIL",
        "capsule": "7227095/v1",
        "commit": "edf5314",
        "manifest_sha256": sha256(code_dir / "source_manifest.tsv"),
        "blocks_checked": len(rows),
        "font_substitution_count": len(font_substitution_rows),
        "supp5b_heatmap_verbatim": heat_exact == heat_patched and not heat_row.get("allowed_replacements"),
        "r_parse_not_applicable_unsupported_blocks": sorted(parse_not_applicable_blocks),
        "generated_blocks_strictly_parseable": len(rows) - len(parse_not_applicable_blocks),
        "r_parse": r_parse,
        "results": result_summary,
        "failures": failures,
    }
    out_path = (
        args.results_dir.resolve() / "audit" / "official_port_verification.json"
        if args.results_dir
        else generated / "official_port_static_verification.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
