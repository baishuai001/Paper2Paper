#!/usr/bin/env python3
"""Parse per-sample HOMER results, summarize motifs, and evaluate the anchor gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


Q_THRESHOLD = 1e-5
PERMUTATIONS = 10_000
SEED = 20260814
PRIMARY_PROMOTER_DEFINITION = "PRIMARY_anchor_promoter_tcga_cpm1_both"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        if not rows:
            raise RuntimeError(f"Cannot infer columns for empty table: {path}")
        columns = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_number(value: str) -> float:
    cleaned = value.strip().replace("%", "").replace(",", "")
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", cleaned)
    return float(match.group(0)) if match else math.nan


def bh_adjust(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=lambda index: values[index])
    output = [1.0] * n
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = n - reverse_rank + 1
        running = min(running, values[index] * n / rank)
        output[index] = min(1.0, running)
    return output


def locate(header: list[str], required_terms: tuple[str, ...], optional: bool = False) -> int | None:
    normalized = [re.sub(r"\s+", " ", item.lower()).strip() for item in header]
    for index, item in enumerate(normalized):
        if all(term in item for term in required_terms):
            return index
    if optional:
        return None
    raise RuntimeError(f"Could not find HOMER column containing {required_terms}: {header}")


def parse_known_results(path: Path, target_total: int, background_total: int) -> list[dict[str, object]]:
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        values = [row for row in reader if row]
    motif_index = 0
    p_index = locate(header, ("p-value",), optional=True)
    if p_index is None:
        p_index = locate(header, ("p value",))
    q_index = locate(header, ("q-value",), optional=True)
    if q_index is None:
        q_index = locate(header, ("benjamini",), optional=True)
    target_index = locate(header, ("target", "with motif"))
    background_index = locate(header, ("background", "with motif"))

    # HOMER can discard a few sequences during genome extraction and reports
    # GC-normalized background motif counts as weighted decimals.  Therefore
    # the input position-file line counts and rounded motif counts are not the
    # correct values for an odds ratio.  Prefer the post-extraction denominator
    # embedded in each knownResults header, retaining the weighted numerator.
    def header_denominator(index: int, fallback: int) -> float:
        match = re.search(r"\(\s*of\s+([0-9.eE+-]+)\s*\)", header[index], flags=re.IGNORECASE)
        if not match:
            return float(fallback)
        value = float(match.group(1))
        return value if math.isfinite(value) and value > 0 else float(fallback)

    homer_target_total = header_denominator(target_index, target_total)
    homer_background_total = header_denominator(background_index, background_total)

    parsed: list[dict[str, object]] = []
    for row in values:
        if max(motif_index, p_index, target_index, background_index) >= len(row):
            continue
        motif_name = row[motif_index]
        match = re.search(r"(JASPAR2024|CIS-BP2\.00)\|([A-Za-z0-9_.-]+)\|([^\s/()]+)", motif_name)
        if not match:
            continue
        source, tf, motif_id = match.groups()
        p_value = parse_number(row[p_index])
        q_value = parse_number(row[q_index]) if q_index is not None and q_index < len(row) else math.nan
        target_with = parse_number(row[target_index])
        background_with = parse_number(row[background_index])
        target_without = max(0.0, homer_target_total - target_with)
        background_without = max(0.0, homer_background_total - background_with)
        odds_ratio = ((target_with + 0.5) * (background_without + 0.5)) / (
            (target_without + 0.5) * (background_with + 0.5)
        )
        parsed.append(
            {
                "database": source,
                "TF": tf,
                "motif_id": motif_id,
                "homer_motif_name": motif_name,
                "p_value": p_value,
                "q_value": q_value,
                "target_with_motif": target_with,
                "target_total": homer_target_total,
                "background_with_motif": background_with,
                "background_total": homer_background_total,
                "odds_ratio": odds_ratio,
                "log2_odds_ratio": math.log2(odds_ratio),
            }
        )
    if parsed and any(not math.isfinite(float(row["q_value"])) for row in parsed):
        adjusted = bh_adjust([float(row["p_value"]) for row in parsed])
        for row, q_value in zip(parsed, adjusted):
            row["q_value"] = q_value
    return parsed


def line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def motif_permutation(system_rows: list[dict[str, object]], inventory: dict[str, dict[str, str]]):
    systems = ("patient", "PDX", "cell_line")
    tfs = sorted(inventory)
    flags = {
        system: {str(row["TF"]): row["system_motif_enriched"] == "TRUE" for row in system_rows if row["system"] == system}
        for system in systems
    }
    observed = sum(all(flags[system].get(tf, False) for system in systems) for tf in tfs)

    def motif_count_bin(value: int) -> str:
        return "1" if value == 1 else ("2-3" if value <= 3 else "4+")

    strata: dict[tuple[str, str], list[str]] = defaultdict(list)
    for tf in tfs:
        record = inventory[tf]
        strata[(record["motif_database_category"], motif_count_bin(int(record["total_tested_motif_count"])))].append(tf)
    rng = random.Random(SEED)
    null_counts: list[int] = []
    for _ in range(PERMUTATIONS):
        permuted: dict[str, dict[str, bool]] = {system: {} for system in systems}
        for system in systems:
            for members in strata.values():
                values = [flags[system].get(tf, False) for tf in members]
                rng.shuffle(values)
                permuted[system].update(dict(zip(members, values)))
        null_counts.append(sum(all(permuted[system][tf] for system in systems) for tf in tfs))
    empirical_p = (1 + sum(value >= observed for value in null_counts)) / (PERMUTATIONS + 1)
    return {
        "observed_triple_system_motif_TFs": observed,
        "permutations": PERMUTATIONS,
        "seed": SEED,
        "matching_strata": "motif database availability category x total tested motif count (1, 2-3, 4+)",
        "null_mean": statistics.mean(null_counts),
        "null_sd": statistics.stdev(null_counts) if len(null_counts) > 1 else 0.0,
        "null_max": max(null_counts),
        "empirical_p_ge_observed": empirical_p,
        "robustness_pass": empirical_p < 0.05,
    }, null_counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    manifest = read_tsv(root / "audit/manifests/figure2_atomic/figure2_atomic_sample_manifest.tsv")
    inventory_rows = read_tsv(root / "results/motif_gate/hc_tf_motif_inventory.tsv")
    inventory = {row["TF"]: row for row in inventory_rows if row["motif_testable"] == "TRUE"}
    background_total = line_count(root / "data/processed/motif_inputs/lung_accessible_CRE_background.homer.pos")
    all_motif_rows: list[dict[str, object]] = []
    best_rows: list[dict[str, object]] = []

    for sample in manifest:
        system = sample["system"]
        slug = sample["sample_slug"]
        target_path = root / f"data/processed/motif_inputs/targets/{system}.{slug}.homer.pos"
        if not target_path.is_file():
            raise RuntimeError(f"Missing HOMER output for {system}/{sample['sample_id']}")
        parsed: list[dict[str, object]] = []
        known_paths: list[Path] = []
        for database in ("JASPAR2024", "CIS-BP2.00"):
            known_path = root / f"results/motif_gate/homer_per_sample/{database}/{system}/{slug}/knownResults.txt"
            known_paths.append(known_path)
            if not known_path.is_file():
                raise RuntimeError(f"Missing independent {database} HOMER output for {system}/{sample['sample_id']}")
            parsed.extend(parse_known_results(known_path, line_count(target_path), background_total))
        if not parsed:
            raise RuntimeError(f"No selected HC-TF motifs parsed from {known_paths}")
        for row in parsed:
            row.update({"system": system, "sample_id": sample["sample_id"], "sample_slug": slug})
            row["motif_enriched_q1e-5"] = str(float(row["q_value"]) < Q_THRESHOLD and float(row["log2_odds_ratio"]) > 0).upper()
        all_motif_rows.extend(parsed)
        by_tf: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in parsed:
            by_tf[str(row["TF"])].append(row)
        for tf in inventory:
            candidates = by_tf.get(tf, [])
            if not candidates:
                best_rows.append(
                    {
                        "system": system,
                        "sample_id": sample["sample_id"],
                        "sample_slug": slug,
                        "TF": tf,
                        "database": inventory[tf]["primary_source"],
                        "motif_id": "NOT_RETURNED",
                        "q_value": 1.0,
                        "log2_odds_ratio": 0.0,
                        "motif_enriched_q1e-5": "FALSE",
                    }
                )
                continue
            enriched_candidates = [
                row for row in candidates
                if float(row["q_value"]) < Q_THRESHOLD and float(row["log2_odds_ratio"]) > 0
            ]
            # A TF is supported when any of its frozen database motifs is
            # significantly enriched.  Prefer the strongest enriched motif;
            # only fall back to the best overall motif when none is enriched.
            selection_pool = enriched_candidates if enriched_candidates else candidates
            selected = sorted(selection_pool, key=lambda row: (float(row["q_value"]), -float(row["log2_odds_ratio"]), str(row["motif_id"])))[0]
            best_rows.append({key: selected[key] for key in (
                "system", "sample_id", "sample_slug", "TF", "database", "motif_id", "q_value", "log2_odds_ratio", "motif_enriched_q1e-5"
            )})

    out = root / "results/motif_gate"
    write_tsv(out / "homer_all_selected_motif_results.tsv", all_motif_rows)
    write_tsv(out / "sample_tf_best_motif_results.tsv", best_rows)
    system_rows: list[dict[str, object]] = []
    for system in ("patient", "PDX", "cell_line"):
        sample_ids = sorted({row["sample_id"] for row in best_rows if row["system"] == system})
        threshold = math.ceil(len(sample_ids) / 2)
        for tf in sorted(inventory):
            rows = [row for row in best_rows if row["system"] == system and row["TF"] == tf]
            enriched = sum(row["motif_enriched_q1e-5"] == "TRUE" for row in rows)
            lors = [float(row["log2_odds_ratio"]) for row in rows]
            system_rows.append(
                {
                    "system": system,
                    "TF": tf,
                    "database": inventory[tf]["primary_source"],
                    "n_samples": len(sample_ids),
                    "half_sample_threshold": threshold,
                    "n_motif_enriched": enriched,
                    "fraction_motif_enriched": enriched / len(sample_ids),
                    "mean_log2_odds_ratio": statistics.mean(lors),
                    "sd_log2_odds_ratio": statistics.stdev(lors) if len(lors) > 1 else 0.0,
                    "system_motif_enriched": str(enriched >= threshold).upper(),
                }
            )
    write_tsv(out / "system_tf_motif_summary.tsv", system_rows)

    triple_rows: list[dict[str, object]] = []
    for tf in sorted(inventory):
        values = [row for row in system_rows if row["TF"] == tf]
        triple = all(row["system_motif_enriched"] == "TRUE" for row in values)
        triple_rows.append(
            {
                "TF": tf,
                "database": inventory[tf]["primary_source"],
                "patient_motif_enriched": next(row["system_motif_enriched"] for row in values if row["system"] == "patient"),
                "PDX_motif_enriched": next(row["system_motif_enriched"] for row in values if row["system"] == "PDX"),
                "cell_line_motif_enriched": next(row["system_motif_enriched"] for row in values if row["system"] == "cell_line"),
                "triple_system_motif_enriched": str(triple).upper(),
            }
        )
    write_tsv(out / "tf_motif_gate_summary.tsv", triple_rows)

    promoter_rows = read_tsv(root / "results/promoter_gate/tf_promoter_gate_summary.tsv")
    primary_hc = [
        row for row in promoter_rows
        if row["analysis_definition"] == PRIMARY_PROMOTER_DEFINITION and row["HC_TF_promoter_activity_definition"] == "TRUE"
    ]
    triple_count = sum(row["triple_system_motif_enriched"] == "TRUE" for row in triple_rows)
    testable_count = len(inventory)
    proportion = triple_count / testable_count if testable_count else 0.0
    anchor_pass = len(primary_hc) >= 10 and triple_count >= 3 and proportion >= 0.10
    motif_robustness, null_counts = motif_permutation(system_rows, inventory)
    write_tsv(
        root / "audit/motif_gate/motif_matched_permutation_null.tsv",
        [{"permutation": index + 1, "triple_system_motif_TFs": value} for index, value in enumerate(null_counts)],
        ["permutation", "triple_system_motif_TFs"],
    )
    (root / "audit/motif_gate/motif_matched_permutation_receipt.json").write_text(
        json.dumps(motif_robustness, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "HOMER_adjusted_q_threshold": Q_THRESHOLD,
        "HC_TFs_after_promoter_activity": len(primary_hc),
        "motif_testable_HC_TFs": testable_count,
        "triple_system_motif_enriched_TFs": triple_count,
        "triple_system_fraction_of_testable": proportion,
        "anchor_requirements": {"minimum_HC_TFs": 10, "minimum_triple_motif_TFs": 3, "minimum_triple_fraction": 0.10},
        "anchor_logic_verdict": "ANCHOR_PASS" if anchor_pass else "FAIL_SIGNAL",
        "motif_matched_permutation": motif_robustness,
        "final_verdict_requires_promoter_matched_permutation": True,
    }
    (root / "audit/motif_gate/motif_gate_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
