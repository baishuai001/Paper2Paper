#!/usr/bin/env python3
"""Summarize fixed-denominator motif-support sensitivity scenarios.

The input must be the harmonized per-sample, per-TF best-motif table created
after applying the frozen JASPAR-first/CIS-BP-fallback rule.  A sample supports
a TF when the selected HOMER enrichment result has q <= the scenario threshold.
The log2 odds ratio is reported but is not used as an extra, non-anchor filter.
Missing or untested sample/TF combinations remain in the frozen system
denominator and therefore count as unsupported.
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


EXPECTED_SYSTEM_SIZES = {"patient": 22, "PDX": 13, "cell_line": 19}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_scenario(value: str) -> tuple[str, float, float]:
    fields = value.split(":")
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("scenario must be label:q_threshold:prevalence")
    label, q_text, prevalence_text = fields
    q_threshold = float(q_text)
    prevalence = float(prevalence_text)
    if not label or not (0 < q_threshold <= 1) or not (0 < prevalence <= 1):
        raise argparse.ArgumentTypeError("invalid scenario label, q threshold, or prevalence")
    return label, q_threshold, prevalence


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--scenario",
        action="append",
        type=parse_scenario,
        default=[],
        help="label:q_threshold:prevalence; may be supplied more than once",
    )
    args = parser.parse_args()

    scenarios = args.scenario or [
        ("primary_anchor", 1e-5, 0.50),
        ("relaxed_q05_half_samples", 0.05, 0.50),
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample_sets: dict[str, set[str]] = defaultdict(set)
    records: list[dict] = []
    motif_testable_hc_tfs: set[str] = set()
    with args.input.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "system",
            "sample_id",
            "TF",
            "TF_role",
            "motif_test_status",
            "q_value",
            "log2_odds_ratio",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"missing input columns: {sorted(missing)}")
        for row in reader:
            system = row["system"]
            sample_sets[system].add(row["sample_id"])
            if "LUAD_HC_TF" not in row["TF_role"]:
                continue
            tested = row["motif_test_status"] == "TESTED"
            if tested:
                motif_testable_hc_tfs.add(row["TF"])
            try:
                q_value = float(row["q_value"]) if tested else math.nan
                log2_odds_ratio = float(row["log2_odds_ratio"]) if tested else math.nan
            except ValueError:
                q_value = math.nan
                log2_odds_ratio = math.nan
            records.append(
                {
                    "system": system,
                    "sample_id": row["sample_id"],
                    "TF": row["TF"],
                    "tested": tested,
                    "q_value": q_value,
                    "log2_odds_ratio": log2_odds_ratio,
                }
            )

    observed_sizes = {system: len(samples) for system, samples in sample_sets.items()}
    if observed_sizes != EXPECTED_SYSTEM_SIZES:
        raise SystemExit(
            f"frozen system denominators changed: observed={observed_sizes}, "
            f"expected={EXPECTED_SYSTEM_SIZES}"
        )
    systems = list(EXPECTED_SYSTEM_SIZES)

    system_rows: list[dict] = []
    triple_rows: list[dict] = []
    scenario_rows: list[dict] = []
    for label, q_threshold, prevalence in scenarios:
        counts: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: {"tested": 0, "enriched": 0})
        for row in records:
            tf = row["TF"]
            if tf not in motif_testable_hc_tfs:
                continue
            key = (row["system"], tf)
            if row["tested"]:
                counts[key]["tested"] += 1
                if row["q_value"] <= q_threshold:
                    counts[key]["enriched"] += 1

        supported_by_system: dict[str, set[str]] = {}
        for system in systems:
            denominator = EXPECTED_SYSTEM_SIZES[system]
            required_count = math.ceil(prevalence * denominator)
            supported_by_system[system] = set()
            for tf in sorted(motif_testable_hc_tfs):
                values = counts[(system, tf)]
                support = values["enriched"] >= required_count
                if support:
                    supported_by_system[system].add(tf)
                system_rows.append(
                    {
                        "scenario": label,
                        "q_threshold": f"{q_threshold:.10g}",
                        "prevalence_threshold": f"{prevalence:.6g}",
                        "system": system,
                        "TF": tf,
                        "n_all_samples": denominator,
                        "n_tested": values["tested"],
                        "n_enriched": values["enriched"],
                        "fraction_enriched_all_samples": f"{values['enriched'] / denominator:.12g}",
                        "support_count_required": required_count,
                        "system_supported": str(support).upper(),
                    }
                )

        triple = set.intersection(*(supported_by_system[system] for system in systems))
        for tf in sorted(triple):
            row = {
                "scenario": label,
                "q_threshold": f"{q_threshold:.10g}",
                "prevalence_threshold": f"{prevalence:.6g}",
                "TF": tf,
            }
            for system in systems:
                values = counts[(system, tf)]
                row[f"{system}_n_enriched"] = values["enriched"]
                row[f"{system}_n_all_samples"] = EXPECTED_SYSTEM_SIZES[system]
            triple_rows.append(row)
        scenario_rows.append(
            {
                "scenario": label,
                "q_threshold": f"{q_threshold:.10g}",
                "prevalence_threshold": f"{prevalence:.6g}",
                "patient_required": math.ceil(prevalence * EXPECTED_SYSTEM_SIZES["patient"]),
                "PDX_required": math.ceil(prevalence * EXPECTED_SYSTEM_SIZES["PDX"]),
                "cell_line_required": math.ceil(prevalence * EXPECTED_SYSTEM_SIZES["cell_line"]),
                "triple_system_hc_tf_count": len(triple),
                "triple_system_hc_tfs": ";".join(sorted(triple)),
            }
        )

    system_path = args.output_dir / "motif_threshold_sensitivity_system_tf.tsv"
    triple_path = args.output_dir / "motif_threshold_sensitivity_triple_system_tfs.tsv"
    scenario_path = args.output_dir / "motif_threshold_sensitivity_scenarios.tsv"
    write_tsv(
        system_path,
        system_rows,
        [
            "scenario",
            "q_threshold",
            "prevalence_threshold",
            "system",
            "TF",
            "n_all_samples",
            "n_tested",
            "n_enriched",
            "fraction_enriched_all_samples",
            "support_count_required",
            "system_supported",
        ],
    )
    triple_fields = ["scenario", "q_threshold", "prevalence_threshold", "TF"]
    for system in systems:
        triple_fields.extend([f"{system}_n_enriched", f"{system}_n_all_samples"])
    write_tsv(triple_path, triple_rows, triple_fields)
    write_tsv(
        scenario_path,
        scenario_rows,
        [
            "scenario",
            "q_threshold",
            "prevalence_threshold",
            "patient_required",
            "PDX_required",
            "cell_line_required",
            "triple_system_hc_tf_count",
            "triple_system_hc_tfs",
        ],
    )

    receipt = {
        "status": "PASS",
        "analysis_role": "supplementary_threshold_sensitivity_primary_unchanged",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(args.input.resolve()),
        "input_sha256": sha256(args.input),
        "selection_rule": "harmonized JASPAR-first/CIS-BP-fallback best motif per sample and TF",
        "enrichment_rule": "HOMER enrichment q_value <= scenario threshold; no extra LOR filter",
        "denominator_rule": "all frozen samples in each system; missing or untested combinations are unsupported",
        "system_sizes": EXPECTED_SYSTEM_SIZES,
        "motif_testable_hc_tf_count": len(motif_testable_hc_tfs),
        "motif_testable_hc_tfs": sorted(motif_testable_hc_tfs),
        "scenarios": scenario_rows,
        "outputs": {
            str(path.name): sha256(path)
            for path in (scenario_path, system_path, triple_path)
        },
    }
    receipt_path = args.output_dir / "motif_threshold_sensitivity_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
