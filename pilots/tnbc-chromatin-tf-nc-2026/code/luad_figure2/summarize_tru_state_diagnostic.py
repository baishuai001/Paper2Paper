#!/usr/bin/env python3
"""Recount frozen LUAD Figure 2 promoter/motif evidence in predefined states.

This is a diagnostic restriction of the immutable 31 HC-TFs. It does not
perform TRU-specific discovery and it never overwrites the all-LUAD B0 result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


SYSTEMS = ("patient", "PDX", "cell_line")
VARIANTS = ("B0", "B1", "B2", "B3")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_true(value: str | None) -> bool:
    return str(value).strip().upper() == "TRUE"


def support_label(count: int, denominator: int) -> str:
    if denominator == 0:
        return "NOT_EVALUABLE_NO_SELECTED_SAMPLES"
    return "TRUE" if count >= math.ceil(denominator / 2) else "FALSE"


def load_state_membership(state_root: Path) -> tuple[list[dict[str, object]], dict[str, dict[str, set[str]]]]:
    paths = {
        "patient": state_root / "results/state_classifier/tcga_luad_atac_rna_matched_state_labels.tsv",
        "PDX": state_root / "results/state_classifier/gse269746_13_pdx_wilkerson506_atac_proxy_classification.tsv",
        "cell_line": state_root / "results/state_classifier/dra001846_19_cellline_wilkerson506_classification.tsv",
    }
    rows_by_system = {system: read_tsv(path) for system, path in paths.items()}
    expected = {"patient": 22, "PDX": 13, "cell_line": 19}
    for system, rows in rows_by_system.items():
        if len(rows) != expected[system]:
            raise RuntimeError(f"{system} state table has {len(rows)} rows; expected {expected[system]}")

    membership: list[dict[str, object]] = []
    selected: dict[str, dict[str, set[str]]] = {
        "TRU_like": {system: set() for system in SYSTEMS},
        "NKX2_1_high_TRU_concordant": {system: set() for system in SYSTEMS},
    }
    for system, rows in rows_by_system.items():
        for row in rows:
            sample_id = row.get("atac_sample_id", "") if system == "patient" else row.get("sample_id", "")
            subtype = row.get("subtype", "").strip()
            classifiable = bool(subtype)
            tru_like = subtype == "TRU"
            nkx_evaluable = system != "PDX" and "nkx2_1_high_tru_concordant" in row
            nkx_tru = nkx_evaluable and is_true(row.get("nkx2_1_high_tru_concordant"))
            if tru_like:
                selected["TRU_like"][system].add(sample_id)
            if nkx_tru:
                selected["NKX2_1_high_TRU_concordant"][system].add(sample_id)
            membership.append({
                "system": system,
                "sample_id": sample_id,
                "subtype": subtype or "UNCLASSIFIABLE_NO_MATCHED_RNA",
                "classifiable": str(classifiable).upper(),
                "TRU_like_selected": str(tru_like).upper(),
                "NKX2_1_high_TRU_evaluable": str(nkx_evaluable).upper(),
                "NKX2_1_high_TRU_selected": str(nkx_tru).upper() if nkx_evaluable else "NOT_EVALUABLE_NO_MATCHED_RNA",
                "state_source": (
                    "Wilkerson506_RNA" if system != "PDX"
                    else "TCGA_validated_Wilkerson506_ATAC_proxy"
                ),
                "diagnostic_qc_status": row.get("worker_diagnostic_qc_status", "NOT_APPLICABLE"),
            })
    return membership, selected


def summarize_state_counts(membership: list[dict[str, object]]) -> list[dict[str, object]]:
    output = []
    for system in SYSTEMS:
        rows = [row for row in membership if row["system"] == system]
        output.append({
            "system": system,
            "total_samples": len(rows),
            "classifiable_samples": sum(row["classifiable"] == "TRUE" for row in rows),
            "TRU_like_samples": sum(row["TRU_like_selected"] == "TRUE" for row in rows),
            "NKX2_1_high_TRU_evaluable_samples": sum(row["NKX2_1_high_TRU_evaluable"] == "TRUE" for row in rows),
            "NKX2_1_high_TRU_selected_samples": sum(row["NKX2_1_high_TRU_selected"] == "TRUE" for row in rows),
        })
    return output


def promoter_diagnostic(
    figure2_root: Path,
    selected: dict[str, dict[str, set[str]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    gate_path = figure2_root / "results/promoter_gate/tf_promoter_gate_summary.tsv"
    sample_path = figure2_root / "results/promoter_gate/sample_tf_promoter_accessibility.tsv"
    gate_rows = read_tsv(gate_path)
    primary_rows = [
        row for row in gate_rows
        if row["analysis_definition"] == "PRIMARY_anchor_promoter_tcga_cpm1_both"
    ]
    hc_tfs = sorted({row["TF"] for row in primary_rows if is_true(row["HC_TF_promoter_activity_definition"])})
    if len(hc_tfs) != 31:
        raise RuntimeError(f"Immutable primary HC-TF count is {len(hc_tfs)}, expected 31")
    sample_rows = [
        row for row in read_tsv(sample_path)
        if row["analysis_definition"] == "PRIMARY_anchor_promoter_tcga_cpm1_both"
        and row["TF"] in hc_tfs
    ]
    accessible = {
        (row["system"], row["sample_id"], row["TF"]): is_true(row["promoter_accessible"])
        for row in sample_rows
    }
    system_rows: list[dict[str, object]] = []
    triple_rows: list[dict[str, object]] = []
    layer = "TRU_like"
    for tf in hc_tfs:
        flags: dict[str, str] = {}
        for system in SYSTEMS:
            samples = selected[layer][system]
            count = sum(accessible.get((system, sample, tf), False) for sample in samples)
            denominator = len(samples)
            status = support_label(count, denominator)
            flags[system] = status
            system_rows.append({
                "layer": layer,
                "system": system,
                "TF": tf,
                "selected_samples": denominator,
                "half_sample_threshold": math.ceil(denominator / 2) if denominator else "NA",
                "accessible_samples": count,
                "support_at_least_half": status,
            })
        triple_rows.append({
            "layer": layer,
            "TF": tf,
            "patient_support": flags["patient"],
            "PDX_support": flags["PDX"],
            "cell_line_support": flags["cell_line"],
            "triple_system_support": str(all(flags[system] == "TRUE" for system in SYSTEMS)).upper(),
        })
    return system_rows, triple_rows, hc_tfs


def motif_paths(figure2_root: Path, variant: str) -> tuple[Path, Path]:
    if variant == "B0":
        root = figure2_root / "results/motif_primary/anchor_consensus_author"
        return (
            root / "anchor_consensus_sample_tf_best_motif.tsv",
            root / "anchor_consensus_hc_tf_triple_system_summary.tsv",
        )
    root = figure2_root / f"results/motif_background_sensitivity/{variant}/parsed"
    return root / f"{variant}_sample_tf_best_motif.tsv", root / f"{variant}_hc_tf_triple_system_summary.tsv"


def motif_diagnostic(
    figure2_root: Path,
    selected: dict[str, dict[str, set[str]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[str]]:
    b0_sample_path, b0_triple_path = motif_paths(figure2_root, "B0")
    b0_triple = read_tsv(b0_triple_path)
    motif_hc_tfs = sorted(row["TF"] for row in b0_triple)
    if len(motif_hc_tfs) != 20:
        raise RuntimeError(f"Immutable motif-testable HC-TF count is {len(motif_hc_tfs)}, expected 20")

    system_rows: list[dict[str, object]] = []
    triple_rows: list[dict[str, object]] = []
    lineage_rows: list[dict[str, object]] = []
    for variant in VARIANTS:
        sample_path, _ = motif_paths(figure2_root, variant)
        rows = read_tsv(sample_path)
        enriched = {
            (row["system"], row["sample_id"], row["TF"]): is_true(row["motif_enriched_q1e-5"])
            for row in rows
        }
        lineage_tfs = sorted({
            row["TF"] for row in rows if "LINEAGE_POSITIVE_CONTROL" in row.get("TF_role", "")
        })
        for layer in ("TRU_like", "NKX2_1_high_TRU_concordant"):
            for tf in motif_hc_tfs:
                flags: dict[str, str] = {}
                for system in SYSTEMS:
                    samples = selected[layer][system]
                    if layer == "NKX2_1_high_TRU_concordant" and system == "PDX":
                        status = "NOT_EVALUABLE_NO_MATCHED_RNA"
                        count = 0
                        denominator = 0
                    else:
                        count = sum(enriched.get((system, sample, tf), False) for sample in samples)
                        denominator = len(samples)
                        status = support_label(count, denominator)
                    flags[system] = status
                    system_rows.append({
                        "variant": variant,
                        "layer": layer,
                        "system": system,
                        "TF": tf,
                        "selected_samples": denominator,
                        "half_sample_threshold": math.ceil(denominator / 2) if denominator else "NA",
                        "motif_enriched_samples": count,
                        "support_at_least_half": status,
                    })
                triple_evaluable = all(flags[system] in {"TRUE", "FALSE"} for system in SYSTEMS)
                triple_rows.append({
                    "variant": variant,
                    "layer": layer,
                    "TF": tf,
                    "patient_support": flags["patient"],
                    "PDX_support": flags["PDX"],
                    "cell_line_support": flags["cell_line"],
                    "triple_system_evaluable": str(triple_evaluable).upper(),
                    "triple_system_support": (
                        str(all(flags[system] == "TRUE" for system in SYSTEMS)).upper()
                        if triple_evaluable else "NOT_EVALUABLE"
                    ),
                })
            for tf in lineage_tfs:
                for system in SYSTEMS:
                    samples = selected[layer][system]
                    if layer == "NKX2_1_high_TRU_concordant" and system == "PDX":
                        status = "NOT_EVALUABLE_NO_MATCHED_RNA"
                        count = 0
                        denominator = 0
                    else:
                        count = sum(enriched.get((system, sample, tf), False) for sample in samples)
                        denominator = len(samples)
                        status = support_label(count, denominator)
                    lineage_rows.append({
                        "variant": variant,
                        "layer": layer,
                        "system": system,
                        "TF": tf,
                        "selected_samples": denominator,
                        "half_sample_threshold": math.ceil(denominator / 2) if denominator else "NA",
                        "motif_enriched_samples": count,
                        "support_at_least_half": status,
                    })
    return system_rows, triple_rows, lineage_rows, motif_hc_tfs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state_root", type=Path)
    parser.add_argument("figure2_root", type=Path)
    args = parser.parse_args()
    state_root = args.state_root.resolve()
    figure2_root = args.figure2_root.resolve()
    required_receipts = [
        state_root / "audit/state_classifier/state_classifier_receipt.json",
        figure2_root / "audit/motif_background_sensitivity/final_receipt.json",
    ]
    for path in required_receipts:
        if not path.is_file():
            raise FileNotFoundError(path)

    out_dir = state_root / "results/tru_state_diagnostic"
    audit_dir = state_root / "audit/tru_state_diagnostic"
    out_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    membership, selected = load_state_membership(state_root)
    state_counts = summarize_state_counts(membership)
    promoter_system, promoter_triple, hc_tfs = promoter_diagnostic(figure2_root, selected)
    motif_system, motif_triple, lineage_rows, motif_hc_tfs = motif_diagnostic(figure2_root, selected)

    summary_rows: list[dict[str, object]] = []
    for variant in VARIANTS:
        supported = sorted({
            row["TF"] for row in motif_triple
            if row["variant"] == variant and row["layer"] == "TRU_like"
            and row["triple_system_support"] == "TRUE"
        })
        summary_rows.append({
            "variant": variant,
            "layer": "TRU_like",
            "triple_system_evaluable": str(all(selected["TRU_like"][system] for system in SYSTEMS)).upper(),
            "triple_system_supported_hc_tf_count": len(supported),
            "triple_system_supported_hc_tfs": ";".join(supported),
            "interpretation": "diagnostic_recount_of_immutable_31_HC_TFs_not_TRU_specific_discovery",
        })
        patient_cell_supported = sorted({
            row["TF"] for row in motif_triple
            if row["variant"] == variant and row["layer"] == "NKX2_1_high_TRU_concordant"
            and row["patient_support"] == "TRUE" and row["cell_line_support"] == "TRUE"
        })
        summary_rows.append({
            "variant": variant,
            "layer": "NKX2_1_high_TRU_concordant",
            "triple_system_evaluable": "FALSE_PDX_HAS_NO_MATCHED_RNA",
            "triple_system_supported_hc_tf_count": "NA",
            "triple_system_supported_hc_tfs": "NA",
            "patient_cell_supported_hc_tf_count": len(patient_cell_supported),
            "patient_cell_supported_hc_tfs": ";".join(patient_cell_supported),
            "interpretation": "two_arm_diagnostic_only_no_PDX_NKX2_1_expression_label",
        })

    write_tsv(out_dir / "state_membership.tsv", membership, [
        "system", "sample_id", "subtype", "classifiable", "TRU_like_selected",
        "NKX2_1_high_TRU_evaluable", "NKX2_1_high_TRU_selected", "state_source",
        "diagnostic_qc_status",
    ])
    write_tsv(out_dir / "state_counts.tsv", state_counts, [
        "system", "total_samples", "classifiable_samples", "TRU_like_samples",
        "NKX2_1_high_TRU_evaluable_samples", "NKX2_1_high_TRU_selected_samples",
    ])
    write_tsv(out_dir / "tru_promoter_system_tf_support.tsv", promoter_system, [
        "layer", "system", "TF", "selected_samples", "half_sample_threshold",
        "accessible_samples", "support_at_least_half",
    ])
    write_tsv(out_dir / "tru_promoter_triple_system_summary.tsv", promoter_triple, [
        "layer", "TF", "patient_support", "PDX_support", "cell_line_support",
        "triple_system_support",
    ])
    write_tsv(out_dir / "state_restricted_motif_system_tf_support.tsv", motif_system, [
        "variant", "layer", "system", "TF", "selected_samples", "half_sample_threshold",
        "motif_enriched_samples", "support_at_least_half",
    ])
    write_tsv(out_dir / "state_restricted_motif_tf_summary.tsv", motif_triple, [
        "variant", "layer", "TF", "patient_support", "PDX_support", "cell_line_support",
        "triple_system_evaluable", "triple_system_support",
    ])
    write_tsv(out_dir / "state_restricted_lineage_positive_controls.tsv", lineage_rows, [
        "variant", "layer", "system", "TF", "selected_samples", "half_sample_threshold",
        "motif_enriched_samples", "support_at_least_half",
    ])
    write_tsv(out_dir / "state_restricted_diagnostic_summary.tsv", summary_rows, [
        "variant", "layer", "triple_system_evaluable",
        "triple_system_supported_hc_tf_count", "triple_system_supported_hc_tfs",
        "patient_cell_supported_hc_tf_count", "patient_cell_supported_hc_tfs",
        "interpretation",
    ])

    tru_counts = {row["system"]: row["TRU_like_samples"] for row in state_counts}
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETE",
        "analysis_role": "diagnostic_state_restriction_primary_B0_unchanged",
        "immutable_hc_tf_count": len(hc_tfs),
        "immutable_motif_testable_hc_tf_count": len(motif_hc_tfs),
        "tru_like_sample_counts": tru_counts,
        "tru_like_three_systems_evaluable": all(int(tru_counts[system]) > 0 for system in SYSTEMS),
        "nkx2_1_high_tru_three_systems_evaluable": False,
        "nkx2_1_high_tru_limitation": "PDX has no exact matched RNA; ATAC proxy defines TRU-like state but cannot establish NKX2-1-high expression.",
        "variant_summaries": summary_rows,
        "upstream_receipts": {
            str(path): sha256(path) for path in required_receipts
        },
        "interpretation_rule": (
            "This recount can diagnose state mixing but cannot establish a TRU-specific regulatory discovery. "
            "Formal discovery requires a separately frozen Figure 1 TRU-versus-rest rerun."
        ),
    }
    receipt_path = audit_dir / "final_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
