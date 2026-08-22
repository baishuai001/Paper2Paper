#!/usr/bin/env python3
"""Verify and compact the completed LUAD Figure 1-2 official-code port."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path
from typing import Iterable


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def truth(value: str) -> bool:
    return value.strip().upper() == "TRUE"


def pdf_pages(path: Path) -> int:
    result = subprocess.run(
        ["pdfinfo", str(path)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=60,
    )
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise RuntimeError(f"pdfinfo did not report a page count for {path}")


def render_first_page(path: Path, scratch: Path, token: str) -> bool:
    prefix = scratch / token
    subprocess.run(
        [
            "pdftoppm",
            "-f",
            "1",
            "-singlefile",
            "-r",
            "18",
            "-png",
            str(path),
            str(prefix),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=120,
    )
    output = prefix.with_suffix(".png")
    return output.exists() and output.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def iter_files(path: Path) -> Iterable[Path]:
    if path.is_file():
        yield path
    elif path.is_dir():
        yield from sorted(p for p in path.rglob("*") if p.is_file())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("official_port_root", type=Path)
    parser.add_argument("--expected-discovery-luad", type=int)
    parser.add_argument("--expected-discovery-lusc", type=int)
    parser.add_argument("--expected-hc-tfs", type=int)
    parser.add_argument("--expected-motif-testable", type=int)
    parser.add_argument("--expected-triple-tfs")
    parser.add_argument("--expected-motif-supported", type=int)
    ns = parser.parse_args()

    run_root = ns.run_root.resolve()
    port_root = ns.official_port_root.resolve()
    audit_root = run_root / "audit"
    audit_root.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []
    failures: list[str] = []

    def check(
        check_id: str,
        observed: object,
        expected: object,
        evidence: Path | str,
        passed: bool | None = None,
    ) -> None:
        is_ok = observed == expected if passed is None else passed
        checks.append(
            {
                "check_id": check_id,
                "status": "PASS" if is_ok else "FAIL",
                "observed": observed,
                "expected": expected,
                "evidence": str(evidence),
            }
        )
        if not is_ok:
            failures.append(check_id)

    full_log = run_root / "full_run.log"
    full_text = full_log.read_text(encoding="utf-8", errors="replace")
    for marker in (
        "Strict replay complete:",
        "Final corrected Figure 1-2 port complete:",
        "Secondary GSE81089 official-block replay complete:",
    ):
        check(f"full_run_marker::{marker}", marker in full_text, True, full_log)

    receipt_specs = (
        (
            "strict_primary_figure1",
            run_root / "runtime_strict_replay/audit/figure1_official_block_receipt.tsv",
            11,
            "GSE41271",
            "CAPSULE_BUG_RETAINED",
        ),
        (
            "final_primary_figure1",
            run_root / "runtime_final_corrected/audit/figure1_official_block_receipt.tsv",
            11,
            "GSE41271",
            "ANCHOR_BUGFIX",
        ),
        (
            "secondary_figure1",
            run_root
            / "runtime_secondary_GSE81089_corrected/audit/figure1_official_block_receipt.tsv",
            11,
            "GSE81089",
            "ANCHOR_BUGFIX",
        ),
    )
    for label, path, expected_rows, dataset, bug_status in receipt_specs:
        rows = read_tsv(path)
        check(f"{label}::row_count", len(rows), expected_rows, path)
        check(
            f"{label}::all_generated",
            all(r["status"] == "GENERATED_OFFICIAL_CONSTRUCTOR" for r in rows),
            True,
            path,
        )
        check(
            f"{label}::validation_dataset",
            sorted({r["validation_dataset"] for r in rows}),
            [dataset],
            path,
        )
        check(
            f"{label}::bug_status",
            sorted({r["anchor_bug_status"] for r in rows}),
            [bug_status],
            path,
        )

    f2_receipt_path = (
        run_root / "runtime_final_corrected/audit/figure2_official_script_receipt.tsv"
    )
    f2_receipt = read_tsv(f2_receipt_path)
    check("final_figure2::row_count", len(f2_receipt), 6, f2_receipt_path)
    check(
        "final_figure2::all_generated",
        all(
            r["status"] == "GENERATED_OFFICIAL_CONSTRUCTOR"
            and r["exit_status"] == "0"
            for r in f2_receipt
        ),
        True,
        f2_receipt_path,
    )
    f2_boundaries = {r["panel"]: r["publication_boundary"] for r in f2_receipt}
    check(
        "figure2c::capsule_print_mismatch",
        f2_boundaries.get("Figure2C_CAPSULE_OUTPUT"),
        "CAPSULE_PRINT_MISMATCH;PAIRED_TCGA_ATAC_RNA_N21_OF_N22",
        f2_receipt_path,
    )
    check(
        "supplementary4c::empty_set_bugfix_boundary",
        f2_boundaries.get("SupplementaryFigure4C"),
        "ANCHOR_EMPTY_SET_BUGFIX",
        f2_receipt_path,
    )

    adapter_audit = run_root / "adapter_final_corrected/audit"
    tf_availability_path = adapter_audit / "figure1_cross_system_heatmap_tf_availability.tsv"
    tf_availability = read_tsv(tf_availability_path)
    displayed = [r for r in tf_availability if truth(r["display_in_cross_system_heatmaps"])]
    missing = [r["TF"] for r in tf_availability if not truth(r["display_in_cross_system_heatmaps"])]
    tcga_stats_path = (
        run_root
        / "adapter_final_corrected/data/Figure1/Tables_ForPlotting/TCGA_TR_activities_NES_viper.tsv"
    )
    tcga_stats = read_tsv(tcga_stats_path)
    discovery = {
        r["tf"].upper() for r in tcga_stats if r["category"] in {"LUAD", "LUSC"}
    }
    discovery_luad = {r["tf"].upper() for r in tcga_stats if r["category"] == "LUAD"}
    discovery_lusc = {r["tf"].upper() for r in tcga_stats if r["category"] == "LUSC"}
    availability_tfs = {r["TF"].upper() for r in tf_availability}
    check("figure1::availability_equals_discovery_set", availability_tfs, discovery,
          tf_availability_path)
    if ns.expected_discovery_luad is not None:
        check("figure1::LUAD_discovery_tf_count", len(discovery_luad),
              ns.expected_discovery_luad, tcga_stats_path)
    if ns.expected_discovery_lusc is not None:
        check("figure1::LUSC_discovery_tf_count", len(discovery_lusc),
              ns.expected_discovery_lusc, tcga_stats_path)
    check(
        "figure1::display_flag_matches_three_system_availability",
        all(
            truth(r["display_in_cross_system_heatmaps"])
            == all(truth(r[k]) for k in ("available_TCGA", "available_PDX", "available_CellLines"))
            for r in tf_availability
        ),
        True,
        tf_availability_path,
    )

    missing_tf_path = adapter_audit / "figure1_model_heatmap_missing_tf_audit.tsv"
    missing_tf = read_tsv(missing_tf_path)
    check("figure1::missing_heatmap_audit_set",
          {r["TF"].upper() for r in missing_tf}, set(missing), missing_tf_path)
    heatmap_matrix_path = (
        run_root
        / "adapter_final_corrected/data/Figure1/Tables_ForPlotting/TCGA_TR_Activity_TNBC150_NonTNBC155.tsv"
    )
    with heatmap_matrix_path.open("r", encoding="utf-8", newline="") as handle:
        matrix_rows = list(csv.reader(handle, delimiter="\t"))[1:]
    matrix_tfs = {r[0].upper() for r in matrix_rows if r}
    check("figure1::heatmap_matrix_equals_displayed_set", matrix_tfs,
          {r["TF"].upper() for r in displayed}, heatmap_matrix_path)

    validation_manifest_path = adapter_audit / "external_validation_slot_manifest.tsv"
    validation_manifest = {r["dataset"]: r for r in read_tsv(validation_manifest_path)}
    check(
        "validation::GSE41271_samples",
        (validation_manifest["GSE41271"]["n_LUAD"], validation_manifest["GSE41271"]["n_LUSC"]),
        ("183", "80"),
        validation_manifest_path,
    )
    check(
        "validation::GSE81089_samples",
        (validation_manifest["GSE81089"]["n_LUAD"], validation_manifest["GSE81089"]["n_LUSC"]),
        ("108", "67"),
        validation_manifest_path,
    )

    def shared_luad_count(path: Path) -> int:
        return len(
            {
                r["TR"].upper()
                for r in read_tsv(path)
                if r["Group"] == "TNBC_METABRIC_TCGA_Shared"
            }
        )

    primary_overlap_path = (
        run_root
        / "adapter_final_corrected/data/Figure1/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv"
    )
    secondary_overlap_path = (
        run_root
        / "adapter_final_corrected/secondary_validation/GSE81089/data/Figure1/Tables_ForPlotting/TNBC_TCGA_METABRIC_comparision.tsv"
    )
    primary_validation_stats_path = (
        run_root
        / "adapter_final_corrected/data/Figure1/Tables_ForPlotting/METABRIC_TR_activities_NES_viper.tsv"
    )
    secondary_validation_stats_path = (
        run_root
        / "adapter_final_corrected/secondary_validation/GSE81089/data/Figure1/Tables_ForPlotting/METABRIC_TR_activities_NES_viper.tsv"
    )
    primary_validation_luad = {
        r["tf"].upper() for r in read_tsv(primary_validation_stats_path)
        if r["category"] == "LUAD"
    }
    secondary_validation_luad = {
        r["tf"].upper() for r in read_tsv(secondary_validation_stats_path)
        if r["category"] == "LUAD"
    }
    check("validation::GSE41271_reproduced_LUAD_TFs",
          shared_luad_count(primary_overlap_path),
          len(discovery_luad & primary_validation_luad), primary_overlap_path)
    check("validation::GSE81089_reproduced_LUAD_TFs",
          shared_luad_count(secondary_overlap_path),
          len(discovery_luad & secondary_validation_luad), secondary_overlap_path)

    pairing_path = adapter_audit / "figure2_tcga_rna_sample_id_mapping_summary.tsv"
    pairing = read_tsv(pairing_path)[0]
    pairing_expected = {
        "n_input_luad_aliquots": "516",
        "n_unique_patient_barcodes": "516",
        "n_duplicate_patient_barcodes": "0",
        "n_promoter_atac_patients": "22",
        "n_promoter_patients_matched_in_activity": "21",
        "n_promoter_patients_missing_rna": "1",
        "missing_rna_patient_ids": "TCGA-44-A47F",
        "paired_21_promoter_activity_exact_set_equal": "TRUE",
        "official_selector_exact_set_equal": "TRUE",
    }
    for field, expected in pairing_expected.items():
        check(f"figure2c_pairing::{field}", pairing[field], expected, pairing_path)

    hc_path = run_root / "adapter_final_corrected/data/Figure2/High_RNA_NES.tsv"
    hc_tfs = [line.strip().upper() for line in hc_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    check("figure2::hc_tf_unique", len(hc_tfs), len(set(hc_tfs)), hc_path)
    if ns.expected_hc_tfs is not None:
        check("figure2::hc_tf_count", len(set(hc_tfs)), ns.expected_hc_tfs, hc_path)

    version_path = adapter_audit / "motif_analysis_version_boundary.tsv"
    version_rows = {r["analysis_slot"]: r for r in read_tsv(version_path)}
    formal = version_rows["FORMAL_PRIMARY_AUTHOR_PRIORITY"]
    candidate = version_rows["CANDIDATE_ANCHOR_VERSION_BOUNDARY"]
    if ns.expected_motif_testable is not None:
        check("motif::formal_testable_hc_tf_count", formal["n_motif_testable_hc_tfs"],
              str(ns.expected_motif_testable), version_path)
    if ns.expected_triple_tfs is not None:
        expected_triple = ";".join(sorted(filter(None, ns.expected_triple_tfs.upper().split(";"))))
        observed_triple = ";".join(sorted(filter(None, formal["triple_system_tfs"].upper().split(";"))))
        check("motif::formal_triple_tfs", observed_triple, expected_triple, version_path)
        check("motif::formal_triple_count", formal["n_triple_system_tfs"],
              str(len(expected_triple.split(";")) if expected_triple else 0), version_path)
    check("motif::formal_database_priority_violations", formal["database_priority_violations"], "0", version_path)
    check("motif::candidate_sensitivity_triple_count", candidate["n_triple_system_tfs"], "0", version_path)
    check("motif::candidate_not_formal", candidate["permitted_use"], "SENSITIVITY_AUDIT_ONLY_NOT_FORMAL_FIGURE", version_path)

    s2f_path = run_root / "runtime_final_corrected/audit/supplementary2f_anchor_bugfix_manifest.tsv"
    s2f = read_tsv(s2f_path)[0]
    check("supplementary2f::bugfix_class", s2f["change_class"], "ANCHOR_BUGFIX", s2f_path)
    check("supplementary2f::changed_line_count", s2f["changed_line_count"], "1", s2f_path)

    s4c_path = run_root / "runtime_final_corrected/audit/supplementary4c_empty_set_bugfix_manifest.tsv"
    s4c = read_tsv(s4c_path)[0]
    adapter_manifest_path = adapter_audit / "adapter_object_manifest.tsv"
    adapter_manifest = read_tsv(adapter_manifest_path)
    inventory_source = Path(next(
        r["source_path"] for r in adapter_manifest
        if r["adapter_object"] == "Figure2_motif_database_map"
    ))
    inventory_rows = [r for r in read_tsv(inventory_source) if truth(r["motif_testable"])]
    inventory_categories = {
        "jaspar_only_count": sum(r["motif_database_category"] == "JASPAR_only" for r in inventory_rows),
        "cisbp_only_count": sum(r["motif_database_category"] == "CIS-BP_only" for r in inventory_rows),
        "both_count": sum(r["motif_database_category"] == "both" for r in inventory_rows),
    }
    for field, expected in {
        "change_class": "ANCHOR_EMPTY_SET_BUGFIX",
        "changed_line_count": "3",
        "jaspar_only_count": str(inventory_categories["jaspar_only_count"]),
        "cisbp_only_count": str(inventory_categories["cisbp_only_count"]),
        "both_count": str(inventory_categories["both_count"]),
        "dummy_tf_added": "False",
        "eulerr_constructor_changed": "False",
        "palette_changed": "False",
    }.items():
        check(f"supplementary4c::{field}", s4c[field], expected, s4c_path)

    figure2de_input_path = adapter_audit / "figure2de_input_receipt.tsv"
    figure2de_input = {
        row["metric"]: row["value"] for row in read_tsv(figure2de_input_path)
    }
    formal_best_table = read_tsv(Path(formal["best_motif_path"]))
    formal_best_rows = sum(
        r.get("motif_test_status", "TESTED") == "TESTED"
        and r.get("TF_role", "LUAD_HC_TF")
        in {"LUAD_HC_TF", "LUAD_HC_TF_AND_LINEAGE_POSITIVE_CONTROL"}
        for r in formal_best_table
    )
    dynamic_figure2de_expected = {
        "motif_testable_hc_tfs": formal["n_motif_testable_hc_tfs"],
        "formal_best_motif_rows": str(formal_best_rows),
        "duplicate_tf_sample_keys": "0",
        "triple_system_supported": formal["n_triple_system_tfs"],
    }
    if ns.expected_motif_supported is not None:
        dynamic_figure2de_expected["supported_in_at_least_one_system"] = str(
            ns.expected_motif_supported
        )
    for metric, expected in dynamic_figure2de_expected.items():
        check(f"figure2de_input::{metric}", figure2de_input[metric], expected,
              figure2de_input_path)

    figure2e_patch_path = (
        run_root
        / "runtime_final_corrected/audit/figure2e_luad_motif_bugfix_manifest.tsv"
    )
    figure2e_patch = read_tsv(figure2e_patch_path)[0]
    for field, expected in {
        "change_class": "LUAD_MOTIF_INPUT_BUGFIX_AND_DECLARED_ZOOM",
        "official_source_unchanged": "True",
        "supported_tf_count": figure2de_input["supported_in_at_least_one_system"],
        "case_tf_count": "6",
        "case_tf_sample_duplicate_keys": "0",
        "etv1_unique_samples": "52",
        "sd_size_direction_fixed": "True",
        "full_observed_range_retained": "True",
        "core_zoom_tfs": "FOXA3;NFATC4;XBP1",
        "core_zoom_reuses_full_range_observations": "True",
        "official_boxpoint_geometry_retained": "True",
        "official_palette_retained": "True",
    }.items():
        check(f"figure2e_patch::{field}", figure2e_patch[field], expected,
              figure2e_patch_path)

    figure2e_runtime_path = (
        run_root
        / "adapter_final_corrected/results/visuals/Figure2/Figure2E_Runtime_Receipt.tsv"
    )
    figure2e_runtime = {
        row["metric"]: row["value"] for row in read_tsv(figure2e_runtime_path)
    }
    for metric, expected in {
        "summary_tf_count": figure2de_input["supported_in_at_least_one_system"],
        "case_tf_count": "6",
        "duplicate_case_tf_sample_keys": "0",
    }.items():
        check(f"figure2e_runtime::{metric}", figure2e_runtime[metric], expected,
              figure2e_runtime_path)

    for required_figure2e_pdf in (
        "Figure2E_Motif_Enrichment_HC-TRs.pdf",
        "Figure2E_Case_Examples_HC-TRs.pdf",
        "Figure2E_Core_Zoom_HC-TRs.pdf",
        "Figure2E_Complete_HC-TRs.pdf",
    ):
        path = figure2e_runtime_path.parent / required_figure2e_pdf
        check(f"figure2e_output::{required_figure2e_pdf}", path.is_file(), True, path)

    substitution_path = run_root / "runtime_final_corrected/audit/runtime_substitution_manifest.tsv"
    substitutions = read_tsv(substitution_path)
    check("runtime::substitution_count", len(substitutions), 68, substitution_path)
    forbidden_classes = sorted(
        {r["class"] for r in substitutions}
        - {"PATH_ONLY", "CANCER_LABEL", "ADAPTER_NAME", "SAMPLE_COUNT"}
    )
    check("runtime::forbidden_substitution_classes", forbidden_classes, [], substitution_path)

    pdf_roots = (
        ("strict_primary", run_root / "adapter_strict_replay/results/visuals"),
        ("final_primary", run_root / "adapter_final_corrected/results/visuals"),
        (
            "secondary_GSE81089",
            run_root
            / "adapter_final_corrected/secondary_validation/GSE81089/results/visuals",
        ),
    )
    pdf_rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="official_fig12_render_") as scratch_name:
        scratch = Path(scratch_name)
        index = 0
        for suite, root in pdf_roots:
            for path in sorted(root.rglob("*.pdf")):
                index += 1
                magic_ok = path.read_bytes()[:5] == b"%PDF-"
                try:
                    pages = pdf_pages(path)
                    render_ok = render_first_page(path, scratch, f"pdf_{index:03d}")
                    render_error = ""
                except Exception as error:  # receipt preserves the exact failure
                    pages = 0
                    render_ok = False
                    render_error = str(error)
                pdf_rows.append(
                    {
                        "suite": suite,
                        "relative_path": path.relative_to(run_root).as_posix(),
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                        "pdf_magic_ok": magic_ok,
                        "pages": pages,
                        "first_page_render_ok": render_ok,
                        "render_error": render_error,
                    }
                )
    pdf_manifest_path = audit_root / "pdf_sha256_manifest.tsv"
    write_tsv(pdf_manifest_path, pdf_rows)
    suite_counts = {
        suite: sum(row["suite"] == suite for row in pdf_rows)
        for suite, _ in pdf_roots
    }
    check("pdf::strict_primary_count", suite_counts["strict_primary"], 14, pdf_manifest_path)
    check("pdf::final_primary_count", suite_counts["final_primary"], 23, pdf_manifest_path)
    check("pdf::secondary_count", suite_counts["secondary_GSE81089"], 14, pdf_manifest_path)
    check("pdf::all_magic_ok", all(bool(r["pdf_magic_ok"]) for r in pdf_rows), True, pdf_manifest_path)
    check("pdf::all_first_pages_render", all(bool(r["first_page_render_ok"]) for r in pdf_rows), True, pdf_manifest_path)

    verification_path = audit_root / "final_verification_receipt.tsv"
    write_tsv(verification_path, checks)
    if failures:
        raise RuntimeError("Verification failed: " + ", ".join(failures))

    package_path = run_root / "official-fig12-port-compact.tar.gz"
    include_roots = [
        full_log,
        run_root / "formal_author_priority_figure2_refresh.log",
        audit_root,
        run_root / "adapter_strict_replay/audit",
        run_root / "adapter_final_corrected/audit",
        run_root / "runtime_strict_replay/audit",
        run_root / "runtime_strict_replay/logs",
        run_root / "runtime_final_corrected/audit",
        run_root / "runtime_final_corrected/logs",
        run_root / "runtime_secondary_GSE81089_corrected/audit",
        run_root / "runtime_secondary_GSE81089_corrected/logs",
        run_root / "adapter_strict_replay/results/visuals",
        run_root / "adapter_final_corrected/results/visuals",
        run_root
        / "adapter_final_corrected/secondary_validation/GSE81089/results/visuals",
        run_root / "sensitivity_candidate_anchor_equivalent",
    ]
    run_files = sorted(
        {
            path.resolve()
            for root in include_roots
            for path in iter_files(root)
            if path.resolve() != package_path.resolve()
            and path.name != "compact_package_sha256.tsv"
        }
    )
    code_files = sorted(
        [p for p in port_root.iterdir() if p.is_file()]
        + list(iter_files(port_root / "upstream_edf5314"))
    )
    package_rows: list[dict[str, object]] = []
    archive_entries: list[tuple[Path, str]] = []
    for path in run_files:
        arcname = "run/" + path.relative_to(run_root).as_posix()
        archive_entries.append((path, arcname))
    for path in code_files:
        if path.suffix == ".pyc" or "__pycache__" in path.parts:
            continue
        arcname = "code/official_port/" + path.relative_to(port_root).as_posix()
        archive_entries.append((path, arcname))
    seen: set[str] = set()
    archive_entries = [
        entry
        for entry in archive_entries
        if not (entry[1] in seen or seen.add(entry[1]))
    ]
    for path, arcname in archive_entries:
        package_rows.append(
            {
                "archive_path": arcname,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    package_manifest_path = audit_root / "compact_package_file_manifest.tsv"
    write_tsv(package_manifest_path, package_rows)
    archive_entries.append(
        (package_manifest_path, "run/audit/compact_package_file_manifest.tsv")
    )

    temporary_package = package_path.with_suffix(".tar.gz.partial")
    with tarfile.open(temporary_package, "w:gz", compresslevel=6) as archive:
        for path, arcname in archive_entries:
            archive.add(path, arcname=arcname, recursive=False)
    temporary_package.replace(package_path)
    package_sha = sha256_file(package_path)
    package_sha_path = audit_root / "compact_package_sha256.tsv"
    write_tsv(
        package_sha_path,
        [
            {
                "package_path": package_path.as_posix(),
                "bytes": package_path.stat().st_size,
                "sha256": package_sha,
                "verification_status": "PASS",
                "pdf_count": len(pdf_rows),
                "file_count": len(archive_entries),
            }
        ],
    )
    print(f"VERIFICATION_PASS checks={len(checks)} pdfs={len(pdf_rows)}")
    print(f"PACKAGE={package_path}")
    print(f"SHA256={package_sha}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
