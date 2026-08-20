#!/usr/bin/env python3
"""Parse the harmonized HOMER audit for LUAD HC TFs and lineage controls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


Q_THRESHOLD = 1e-5
SYSTEM_ORDER = ("patient", "PDX", "cell_line")
LINEAGE_CONTROLS = ("NKX2-1", "FOXA1", "FOXA2", "CEBPA", "CEBPB", "GATA6", "HNF4A", "ELF3")


def tf_role(tf: str, hc_tfs: set[str]) -> str:
    if tf in hc_tfs and tf in LINEAGE_CONTROLS:
        return "LUAD_HC_TF_AND_LINEAGE_POSITIVE_CONTROL"
    if tf in hc_tfs:
        return "LUAD_HC_TF"
    return "LUAD_LINEAGE_POSITIVE_CONTROL"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_number(value: str) -> float:
    cleaned = value.strip().replace("%", "").replace(",", "")
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", cleaned)
    return float(match.group(0)) if match else math.nan


def locate(header: list[str], terms: tuple[str, ...], optional: bool = False) -> int | None:
    normalized = [re.sub(r"\s+", " ", value.lower()).strip() for value in header]
    for index, value in enumerate(normalized):
        if all(term in value for term in terms):
            return index
    if optional:
        return None
    raise RuntimeError(f"Missing HOMER column {terms}: {header}")


def parse_meme_gene_map(path: Path, source: str) -> dict[tuple[str, str], set[str]]:
    mapping: dict[tuple[str, str], set[str]] = defaultdict(set)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("MOTIF "):
                continue
            fields = line.strip().split(maxsplit=2)
            if len(fields) < 3:
                continue
            motif_id, gene = fields[1], fields[2]
            mapping[(source, motif_id)].add(gene.upper())
    return mapping


def build_token_map(
    root: Path,
    hc_tfs: set[str],
    motif_database_path: Path | None = None,
) -> tuple[dict[str, set[str]], dict[str, list[str]]]:
    token_map: dict[str, set[str]] = defaultdict(set)
    mapping_rows = read_tsv(root / "results/motif_gate/selected_motif_mapping.tsv")
    for row in mapping_rows:
        if row["TF"] in hc_tfs:
            token_map[row["homer_motif_token"]].add(row["TF"])
    gene_map: dict[tuple[str, str], set[str]] = defaultdict(set)
    for key, values in parse_meme_gene_map(
        root / "reference/motifs/downloads/JASPAR2024_CORE_vertebrates_non-redundant.meme", "JASPAR2024"
    ).items():
        gene_map[key].update(values)
    for key, values in parse_meme_gene_map(
        root / "reference/motifs/CIS-BP_2.00_Homo_sapiens.meme", "CIS-BP2.00"
    ).items():
        gene_map[key].update(values)
    motif_headers: dict[str, list[str]] = defaultdict(list)
    if motif_database_path is None:
        motif_database_path = root / "reference/motifs/full_reference_both_databases.homer"
    with motif_database_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            fields = line[1:].strip().split()
            if len(fields) < 2:
                continue
            token = fields[1]
            match = re.fullmatch(r"(JASPAR2024|CIS-BP2\.00)\|[0-9]{6}\|(.+)", token)
            if not match:
                continue
            source, motif_id = match.groups()
            for gene in gene_map.get((source, motif_id), set()):
                if gene in {control.upper() for control in LINEAGE_CONTROLS}:
                    canonical = next(control for control in LINEAGE_CONTROLS if control.upper() == gene)
                    token_map[token].add(canonical)
                    motif_headers[canonical].append(token)
    return token_map, {tf: sorted(set(tokens)) for tf, tokens in motif_headers.items()}


def parse_known_results(path: Path, token_map: dict[str, set[str]]) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        values = [row for row in reader if row]
    p_index = locate(header, ("p-value",), optional=True)
    if p_index is None:
        p_index = locate(header, ("p value",))
    q_index = locate(header, ("q-value",), optional=True)
    if q_index is None:
        q_index = locate(header, ("benjamini",))
    target_index = locate(header, ("target", "with motif"))
    background_index = locate(header, ("background", "with motif"))

    def denominator(index: int) -> float:
        match = re.search(r"\(\s*of\s+([0-9.eE+-]+)\s*\)", header[index], flags=re.IGNORECASE)
        return float(match.group(1)) if match else math.nan

    target_total = denominator(target_index)
    background_total = denominator(background_index)
    parsed: list[dict[str, object]] = []
    for row in values:
        if max(q_index, target_index, background_index) >= len(row):
            continue
        motif_name = row[0]
        token_match = re.search(r"(JASPAR2024|CIS-BP2\.00)\|[0-9]{6}\|[^\s/()]+", motif_name)
        if not token_match:
            continue
        token = token_match.group(0)
        if token not in token_map:
            continue
        q_value = parse_number(row[q_index])
        p_value = parse_number(row[p_index])
        target_with = parse_number(row[target_index])
        background_with = parse_number(row[background_index])
        target_without = max(0.0, target_total - target_with)
        background_without = max(0.0, background_total - background_with)
        odds = ((target_with + 0.5) * (background_without + 0.5)) / (
            (target_without + 0.5) * (background_with + 0.5)
        )
        for tf in token_map[token]:
            parsed.append(
                {
                    "TF": tf,
                    "motif_token": token,
                    "p_value": p_value,
                    "q_value": q_value,
                    "target_with_motif": target_with,
                    "target_total": target_total,
                    "background_with_motif": background_with,
                    "background_total": background_total,
                    "log2_odds_ratio": math.log2(odds),
                    # The TNBC author code uses an inclusive threshold
                    # (qvalue_Benjamini <= 0.00001).
                    "motif_enriched_q1e-5": str(q_value <= Q_THRESHOLD).upper(),
                }
            )
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument(
        "--homer-results-relative",
        action="append",
        default=None,
        help="May be repeated; q values are retained from each separately run HOMER database before TF-level union.",
    )
    parser.add_argument(
        "--output-relative",
        default="results/motif_equivalence",
    )
    parser.add_argument(
        "--motif-database-relative",
        default="reference/motifs/full_reference_both_databases.homer",
    )
    parser.add_argument(
        "--receipt-relative",
        default="audit/motif_equivalence/harmonized_motif_receipt.json",
    )
    parser.add_argument(
        "--analysis-label",
        default="full JASPAR2024 plus CIS-BP2.00 database sensitivity analysis",
    )
    parser.add_argument(
        "--manifest-relative",
        default="data/processed/motif_equivalence/harmonized_200bp/harmonized_motif_manifest.tsv",
    )
    parser.add_argument(
        "--background-relative",
        default="data/processed/motif_equivalence/harmonized_200bp/LUSC_accessible_common_GC_matched_background.bed",
    )
    parser.add_argument(
        "--output-prefix",
        default="harmonized",
        help="Filename prefix only; defaults preserve the completed diagnostic audit paths.",
    )
    parser.add_argument(
        "--prefer-jaspar-then-cisbp",
        action="store_true",
        help="Match the TNBC author logic: use JASPAR motifs for a TF when available and CIS-BP only otherwise.",
    )
    parser.add_argument(
        "--not-testable-status",
        default="NOT_TESTABLE_PEAK_COUNT_MATCH",
        help="Status recorded for fixed samples that cannot enter HOMER; defaults preserve the sensitivity-analysis vocabulary.",
    )
    args = parser.parse_args()
    root = args.run_root.resolve()
    manifest_path = root / args.manifest_relative
    background_path = root / args.background_relative
    manifest = read_tsv(manifest_path)
    inventory = read_tsv(root / "results/motif_gate/hc_tf_motif_inventory.tsv")
    hc_tfs = {row["TF"] for row in inventory if row["motif_testable"] == "TRUE"}
    motif_database_path = root / args.motif_database_relative
    token_map, control_catalog = build_token_map(root, hc_tfs, motif_database_path)
    preferred_source: dict[str, str] = {}
    if args.prefer_jaspar_then_cisbp:
        sources_by_tf: dict[str, set[str]] = defaultdict(set)
        for token, tfs in token_map.items():
            source = token.split("|", 1)[0]
            for tf in tfs:
                sources_by_tf[tf].add(source)
        preferred_source = {
            tf: "JASPAR2024" if "JASPAR2024" in sources else "CIS-BP2.00"
            for tf, sources in sources_by_tf.items()
        }
    homer_result_relatives = args.homer_results_relative or ["results/motif_equivalence/harmonized_homer"]
    homer_results = [root / relative for relative in homer_result_relatives]
    out = root / args.output_relative
    all_rows: list[dict[str, object]] = []
    best_rows: list[dict[str, object]] = []
    all_tfs = sorted(hc_tfs | set(LINEAGE_CONTROLS))
    for sample in manifest:
        system, sample_id, slug = sample["system"], sample["sample_id"], sample["sample_slug"]
        matchable = sample["motif_matchable"] == "TRUE"
        if matchable:
            parsed = []
            for homer_result_root in homer_results:
                known = homer_result_root / system / slug / "knownResults.txt"
                if not known.is_file():
                    raise FileNotFoundError(known)
                source_rows = parse_known_results(known, token_map)
                if args.prefer_jaspar_then_cisbp:
                    source_rows = [
                        row
                        for row in source_rows
                        if str(row["motif_token"]).startswith(f"{preferred_source[str(row['TF'])]}|")
                    ]
                parsed.extend(source_rows)
            for row in parsed:
                row.update({"system": system, "sample_id": sample_id, "sample_slug": slug})
            all_rows.extend(parsed)
            by_tf: dict[str, list[dict[str, object]]] = defaultdict(list)
            for row in parsed:
                by_tf[str(row["TF"])].append(row)
        else:
            by_tf = defaultdict(list)
        for tf in all_tfs:
            candidates = by_tf.get(tf, [])
            if not matchable:
                status, selected = args.not_testable_status, None
            elif not candidates:
                status, selected = "TESTED_MOTIF_NOT_AVAILABLE_OR_NOT_RETURNED", None
            else:
                significant = [row for row in candidates if row["motif_enriched_q1e-5"] == "TRUE"]
                pool = significant if significant else candidates
                selected = sorted(pool, key=lambda row: (float(row["q_value"]), -float(row["log2_odds_ratio"]), str(row["motif_token"])))[0]
                status = "TESTED"
            best_rows.append(
                {
                    "system": system,
                    "sample_id": sample_id,
                    "sample_slug": slug,
                    "TF": tf,
                    "TF_role": tf_role(tf, hc_tfs),
                    "motif_test_status": status,
                    "motif_token": selected["motif_token"] if selected else "NA",
                    "q_value": selected["q_value"] if selected else 1.0,
                    "log2_odds_ratio": selected["log2_odds_ratio"] if selected else 0.0,
                    "motif_enriched_q1e-5": selected["motif_enriched_q1e-5"] if selected else "FALSE",
                }
            )
    write_tsv(out / f"{args.output_prefix}_all_selected_motif_results.tsv", all_rows)
    write_tsv(out / f"{args.output_prefix}_sample_tf_best_motif.tsv", best_rows)

    system_rows: list[dict[str, object]] = []
    for system in SYSTEM_ORDER:
        total_ids = sorted({row["sample_id"] for row in manifest if row["system"] == system})
        testable_ids = sorted({row["sample_id"] for row in manifest if row["system"] == system and row["motif_matchable"] == "TRUE"})
        for tf in all_tfs:
            rows = [row for row in best_rows if row["system"] == system and row["TF"] == tf]
            enriched = sum(row["motif_enriched_q1e-5"] == "TRUE" for row in rows)
            tested = [row for row in rows if row["motif_test_status"] == "TESTED"]
            lors = [float(row["log2_odds_ratio"]) for row in tested]
            total_threshold = math.ceil(len(total_ids) / 2)
            testable_threshold = math.ceil(len(testable_ids) / 2)
            system_rows.append(
                {
                    "system": system,
                    "TF": tf,
                    "TF_role": tf_role(tf, hc_tfs),
                    "n_all_samples": len(total_ids),
                    "n_peak_count_matchable_samples": len(testable_ids),
                    "n_motif_enriched": enriched,
                    "fraction_enriched_all_samples": enriched / len(total_ids),
                    "fraction_enriched_matchable_samples": enriched / len(testable_ids) if testable_ids else math.nan,
                    "half_threshold_all_samples": total_threshold,
                    "half_threshold_matchable_samples": testable_threshold,
                    "support_fixed_original_denominator": str(enriched >= total_threshold).upper(),
                    "support_matchable_denominator": str(enriched >= testable_threshold).upper(),
                    "mean_log2_odds_ratio_tested": statistics.mean(lors) if lors else math.nan,
                }
            )
    write_tsv(out / f"{args.output_prefix}_system_tf_motif_summary.tsv", system_rows)

    triple_rows: list[dict[str, object]] = []
    for tf in sorted(hc_tfs):
        rows = [row for row in system_rows if row["TF"] == tf]
        triple_rows.append(
            {
                "TF": tf,
                "patient_fixed_support": next(row["support_fixed_original_denominator"] for row in rows if row["system"] == "patient"),
                "PDX_fixed_support": next(row["support_fixed_original_denominator"] for row in rows if row["system"] == "PDX"),
                "cell_line_fixed_support": next(row["support_fixed_original_denominator"] for row in rows if row["system"] == "cell_line"),
                "triple_system_fixed_original_denominator": str(all(row["support_fixed_original_denominator"] == "TRUE" for row in rows)).upper(),
                "patient_matchable_support": next(row["support_matchable_denominator"] for row in rows if row["system"] == "patient"),
                "PDX_matchable_support": next(row["support_matchable_denominator"] for row in rows if row["system"] == "PDX"),
                "cell_line_matchable_support": next(row["support_matchable_denominator"] for row in rows if row["system"] == "cell_line"),
                "triple_system_matchable_denominator": str(all(row["support_matchable_denominator"] == "TRUE" for row in rows)).upper(),
            }
        )
    write_tsv(out / f"{args.output_prefix}_hc_tf_triple_system_summary.tsv", triple_rows)
    # A lineage control remains a positive control even when it is also one of
    # the LUAD HC-TFs (ELF3 in the current frozen candidate set).
    positive_rows = [row for row in system_rows if row["TF"] in LINEAGE_CONTROLS]
    for row in positive_rows:
        row["catalogued_motif_tokens"] = ";".join(control_catalog.get(str(row["TF"]), []))
    write_tsv(out / "luad_lineage_positive_control_summary.tsv", positive_rows)

    previous_path = root / "results/motif_gate/tf_motif_gate_summary.tsv"
    previous = read_tsv(previous_path) if previous_path.is_file() else []
    previous_triple = sum(row.get("triple_system_motif_enriched") == "TRUE" for row in previous)
    motif_count = sum(1 for line in motif_database_path.open(encoding="utf-8") if line.startswith(">"))
    target_denominators = sorted({float(row["target_total"]) for row in all_rows})
    background_denominators = sorted({float(row["background_total"]) for row in all_rows})
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_label": args.analysis_label,
        "tf_database_selection_rule": (
            "JASPAR2024_if_available_else_CIS-BP2.00"
            if args.prefer_jaspar_then_cisbp
            else "all_catalogued_motifs_in_configured_database"
        ),
        "homer_result_roots": [str(path) for path in homer_results],
        "q_threshold": Q_THRESHOLD,
        "unified_motif_database": str(motif_database_path),
        "unified_motif_database_sha256": sha256(motif_database_path),
        "unified_motif_database_motif_count": motif_count,
        "parsed_selected_motif_result_row_count": len(all_rows),
        "homer_target_denominator_values": target_denominators,
        "homer_effective_background_denominator_min": min(background_denominators) if background_denominators else None,
        "homer_effective_background_denominator_max": max(background_denominators) if background_denominators else None,
        "homer_effective_background_note": "the same input BED is reused; HOMER removes target-overlapping background regions and reports the resulting per-sample effective denominator",
        "shared_background_sha256": sha256(background_path),
        "analysis_manifest": str(manifest_path),
        "analysis_manifest_sha256": sha256(manifest_path),
        "analysis_background": str(background_path),
        "harmonized_manifest_sha256": sha256(manifest_path),
        "completed_or_explicitly_not_testable_samples_by_system": {
            system: sum(row["system"] == system for row in manifest) for system in SYSTEM_ORDER
        },
        "peak_count_matchable_samples_by_system": {
            system: sum(row["system"] == system and row["motif_matchable"] == "TRUE" for row in manifest)
            for system in SYSTEM_ORDER
        },
        "motif_testable_luad_hc_tf_count": len(hc_tfs),
        "lineage_positive_controls": list(LINEAGE_CONTROLS),
        "lineage_controls_with_at_least_one_catalogued_motif": sum(bool(control_catalog.get(tf)) for tf in LINEAGE_CONTROLS),
        "previous_variable_width_self_union_background_triple_hc_tf_count": previous_triple,
        "harmonized_triple_hc_tf_count_fixed_original_denominator": sum(row["triple_system_fixed_original_denominator"] == "TRUE" for row in triple_rows),
        "harmonized_triple_hc_tf_count_matchable_denominator": sum(row["triple_system_matchable_denominator"] == "TRUE" for row in triple_rows),
        "positive_control_support_by_system_fixed_denominator": {
            system: sum(row["system"] == system and row["support_fixed_original_denominator"] == "TRUE" for row in positive_rows)
            for system in SYSTEM_ORDER
        },
        "positive_control_support_by_system_matchable_denominator": {
            system: sum(row["system"] == system and row["support_matchable_denominator"] == "TRUE" for row in positive_rows)
            for system in SYSTEM_ORDER
        },
        "interpretation_rule": "failure of known LUAD lineage controls in a model arm indicates an assay, pipeline, or motif-model-selection calibration problem before absence of LUAD HC-TF motifs can be interpreted as clean biology",
    }
    receipt_path = root / args.receipt_relative
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
