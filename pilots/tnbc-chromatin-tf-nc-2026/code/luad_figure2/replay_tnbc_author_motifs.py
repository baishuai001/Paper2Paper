#!/usr/bin/env python3
"""Replay the TNBC author motif aggregation and compare it with published outputs.

This script starts from the per-sample HOMER ``knownResults.txt`` files after
the author's unmodified 06B R script has combined them.  It mirrors the
thresholding and database-selection logic in Figure2/07A and the log-odds
summary logic in Figure2/08, then compares every TF/system call with the
published Code Ocean reference tables.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


Q_THRESHOLD = 1e-5
SYSTEMS = ("TCGA", "PDX", "CellLines")
EXPECTED_TOTALS = {"TCGA": 7, "PDX": 38, "CellLines": 26}
EXPECTED_CUTOFFS = {system: math.ceil(total / 2) for system, total in EXPECTED_TOTALS.items()}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_float(value: str) -> float:
    value = value.strip().replace("%", "")
    if value in {"", "NA", "NaN", "nan"}:
        return math.nan
    return float(value)


def read_mapping(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle)
        mapping = {row["TR"].strip().upper(): row["Database"].strip() for row in rows}
    if not mapping:
        raise RuntimeError(f"No TF-to-database mappings in {path}")
    return mapping


def read_combined(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) != 12:
            raise RuntimeError(f"Expected 12 HOMER columns in {path}, found {len(header)}")
        for values in reader:
            if not values:
                continue
            if len(values) != 12:
                raise RuntimeError(f"Malformed row with {len(values)} columns in {path}: {values[:3]}")
            rows.append(
                {
                    "TF": values[0].strip().upper(),
                    "Consensus": values[1],
                    "P_value": parse_float(values[2]),
                    "LogPvalue": parse_float(values[3]),
                    "q_value": parse_float(values[4]),
                    "target_with": parse_float(values[5]),
                    "target_pct": parse_float(values[6]) / 100.0,
                    "background_with": parse_float(values[7]),
                    "background_pct": parse_float(values[8]) / 100.0,
                    "Sample": values[9],
                    "Sample_Group": values[10],
                    "Database": values[11],
                }
            )
    return rows


def read_expected_intersection(path: Path) -> dict[str, dict[str, object]]:
    expected: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header[:4] != ["TCGA", "PDX", "CellLines", "hc_group"]:
            raise RuntimeError(f"Unexpected reference intersection header: {header}")
        for values in reader:
            if not values:
                continue
            if len(values) != 5:
                raise RuntimeError(f"Expected row-name plus four fields, got {values}")
            expected[values[0].upper()] = {
                "TCGA": int(values[1]),
                "PDX": int(values[2]),
                "CellLines": int(values[3]),
                "hc_group": values[4],
            }
    return expected


def read_author_summary(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {(row["TR"].upper(), row["Sample_Group"]): row for row in rows}


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sample_sd(values: list[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else math.nan


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("tnbc_root", type=Path)
    args = parser.parse_args()
    root = args.tnbc_root.resolve()
    source = root / "data/codeocean"
    replay_combined = root / "results/author_code_replay/combined"
    out = root / "results/author_code_replay"
    audit = root / "audit"
    out.mkdir(parents=True, exist_ok=True)
    audit.mkdir(parents=True, exist_ok=True)

    mapping_path = source / "JASPAR_CISBP_comparison_long.csv"
    expected_path = source / "JASPAR_TR_MotifEnrichment_Intersection_Sample_Groups.tsv"
    author_summary_path = source / "TR_MotifEnrichment_HC-TRs.tsv"
    mapping = read_mapping(mapping_path)
    expected = read_expected_intersection(expected_path)

    published_paths = {
        database: source / f"{database}_Combined_MotifEnrichment_TNBC.tsv"
        for database in ("JASPAR", "CISBP")
    }
    replay_paths = {
        database: replay_combined / f"{database}_Combined_MotifEnrichment_TNBC.tsv"
        for database in ("JASPAR", "CISBP")
    }
    for path in [*published_paths.values(), *replay_paths.values(), mapping_path, expected_path, author_summary_path]:
        if not path.is_file():
            raise FileNotFoundError(path)

    byte_identity = {
        database: sha256(replay_paths[database]) == sha256(published_paths[database])
        for database in ("JASPAR", "CISBP")
    }
    jaspar = read_combined(replay_paths["JASPAR"])
    cisbp = read_combined(replay_paths["CISBP"])
    cisbp_only = {tf for tf, database in mapping.items() if database == "CISBP_only"}
    selected = [row for row in jaspar if row["TF"] in mapping]
    selected.extend(row for row in cisbp if row["TF"] in cisbp_only)

    observed_samples = {
        database: {
            system: len({str(row["Sample"]) for row in rows if row["Sample_Group"] == system})
            for system in SYSTEMS
        }
        for database, rows in (("JASPAR", jaspar), ("CISBP", cisbp))
    }
    jaspar_sample_totals_match = observed_samples["JASPAR"] == EXPECTED_TOTALS

    significant_samples: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in selected:
        q_value = float(row["q_value"])
        if not math.isnan(q_value) and q_value <= Q_THRESHOLD:
            significant_samples[(str(row["TF"]), str(row["Sample_Group"]))].add(str(row["Sample"]))

    support_rows: list[dict[str, object]] = []
    computed_all: dict[str, dict[str, object]] = {}
    for tf in sorted(mapping):
        record: dict[str, object] = {
            "TF": tf,
            "hc_group": mapping[tf].replace("_", " "),
        }
        any_support = False
        for system in SYSTEMS:
            n_positive = len(significant_samples[(tf, system)])
            supported = n_positive >= EXPECTED_CUTOFFS[system]
            record[f"{system}_n_q_le_1e-5"] = n_positive
            record[f"{system}_cutoff"] = EXPECTED_CUTOFFS[system]
            record[system] = int(supported)
            any_support = any_support or supported
        record["any_system_support"] = str(any_support).upper()
        record["triple_system_support"] = str(all(int(record[system]) == 1 for system in SYSTEMS)).upper()
        support_rows.append(record)
        if any_support:
            computed_all[tf] = record
    write_tsv(out / "tnbc_replayed_all_mapped_tf_support.tsv", support_rows)
    write_tsv(out / "tnbc_replayed_supported_intersection.tsv", list(computed_all.values()))

    mismatches: list[dict[str, object]] = []
    for tf in sorted(set(computed_all) | set(expected)):
        observed = computed_all.get(tf)
        reference = expected.get(tf)
        for field in (*SYSTEMS, "hc_group"):
            observed_value = observed.get(field) if observed else "ABSENT"
            expected_value = reference.get(field) if reference else "ABSENT"
            if observed_value != expected_value:
                mismatches.append(
                    {
                        "TF": tf,
                        "field": field,
                        "replayed": observed_value,
                        "published": expected_value,
                    }
                )
    write_tsv(out / "tnbc_replay_vs_published_intersection_mismatches.tsv", mismatches,
              ["TF", "field", "replayed", "published"])

    supported_tfs = set(computed_all)
    lor_rows: list[dict[str, object]] = []
    total_per_group = {
        system: len({str(row["Sample"]) for row in selected if row["TF"] in supported_tfs and row["Sample_Group"] == system})
        for system in SYSTEMS
    }
    lor_values: dict[tuple[str, str], list[tuple[str, float]]] = defaultdict(list)
    for row in selected:
        tf = str(row["TF"])
        system = str(row["Sample_Group"])
        if tf not in supported_tfs or system not in SYSTEMS:
            continue
        k1 = float(row["target_with"])
        k0 = float(row["background_with"])
        p1 = float(row["target_pct"])
        p0 = float(row["background_pct"])
        n1 = round(k1 / p1) if not math.isnan(p1) and p1 > 0 else k1 + 1
        n0 = round(k0 / p0) if not math.isnan(p0) and p0 > 0 else k0 + 1
        odds = ((k1 + 0.5) / (n1 - k1 + 0.5)) / ((k0 + 0.5) / (n0 - k0 + 0.5))
        lor_values[(tf, system)].append((str(row["Sample"]), math.log2(odds)))
    for tf in sorted(supported_tfs):
        for system in SYSTEMS:
            values = lor_values[(tf, system)]
            samples = {sample for sample, _ in values}
            lors = [value for _, value in values]
            lor_rows.append(
                {
                    "TR": tf,
                    "Sample_Group": system,
                    "Mean_LOR": statistics.mean(lors),
                    "SD_LOR": sample_sd(lors),
                    "n_samples": len(samples),
                    "total": total_per_group[system],
                    "Pct_samples": len(samples) / total_per_group[system] * 100,
                    "sig_cohort": "TRUE" if int(computed_all[tf][system]) == 1 else "FALSE",
                }
            )
    write_tsv(out / "tnbc_replayed_motif_lor_summary.tsv", lor_rows)

    author_summary = read_author_summary(author_summary_path)
    numeric_fields = ("Mean_LOR", "SD_LOR", "Pct_samples")
    integer_fields = ("n_samples", "total")
    max_abs_delta = {field: 0.0 for field in numeric_fields}
    summary_mismatches: list[dict[str, object]] = []
    for row in lor_rows:
        key = (str(row["TR"]), str(row["Sample_Group"]))
        reference = author_summary.get(key)
        if reference is None:
            summary_mismatches.append({"TR": key[0], "Sample_Group": key[1], "field": "row", "replayed": "PRESENT", "published": "ABSENT"})
            continue
        for field in numeric_fields:
            delta = abs(float(row[field]) - float(reference[field]))
            max_abs_delta[field] = max(max_abs_delta[field], delta)
            if delta > 1e-12:
                summary_mismatches.append({"TR": key[0], "Sample_Group": key[1], "field": field, "replayed": row[field], "published": reference[field]})
        for field in integer_fields:
            if int(row[field]) != int(reference[field]):
                summary_mismatches.append({"TR": key[0], "Sample_Group": key[1], "field": field, "replayed": row[field], "published": reference[field]})
        if str(row["sig_cohort"]).upper() != str(reference["sig_cohort"]).upper():
            summary_mismatches.append({"TR": key[0], "Sample_Group": key[1], "field": "sig_cohort", "replayed": row["sig_cohort"], "published": reference["sig_cohort"]})
    for key in sorted(set(author_summary) - {(str(row["TR"]), str(row["Sample_Group"])) for row in lor_rows}):
        summary_mismatches.append({"TR": key[0], "Sample_Group": key[1], "field": "row", "replayed": "ABSENT", "published": "PRESENT"})
    write_tsv(out / "tnbc_replay_vs_published_lor_mismatches.tsv", summary_mismatches,
              ["TR", "Sample_Group", "field", "replayed", "published"])

    triple_tfs = sorted(tf for tf, row in computed_all.items() if all(int(row[system]) == 1 for system in SYSTEMS))
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_scope": "TNBC positive-control replay from per-sample published HOMER knownResults through author 06B plus an exact 07A/08 logic reimplementation",
        "not_claimed": "HOMER was not rerun from ATAC peak BEDs because the Code Ocean capsule exposes per-sample HOMER outputs but not all raw/query peak inputs; PDX raw ATAC is controlled-access",
        "q_threshold_inclusive": Q_THRESHOLD,
        "fixed_system_totals": EXPECTED_TOTALS,
        "half_sample_cutoffs": EXPECTED_CUTOFFS,
        "observed_samples_by_database": observed_samples,
        "jaspar_sample_totals_match_author_constants": jaspar_sample_totals_match,
        "author_06B_combined_tables_byte_identical_to_published": byte_identity,
        "mapped_tf_count": len(mapping),
        "replayed_supported_in_at_least_one_system_count": len(computed_all),
        "published_supported_in_at_least_one_system_count": len(expected),
        "replayed_triple_system_count": len(triple_tfs),
        "published_triple_system_count": sum(all(int(row[system]) == 1 for system in SYSTEMS) for row in expected.values()),
        "replayed_triple_system_tfs": triple_tfs,
        "intersection_mismatch_count": len(mismatches),
        "lor_summary_mismatch_count_at_tolerance_1e-12": len(summary_mismatches),
        "lor_summary_max_absolute_numeric_delta": max_abs_delta,
        "full_published_three_system_result_recovered": bool(
            all(byte_identity.values())
            and jaspar_sample_totals_match
            and len(computed_all) == len(expected) == 43
            and len(triple_tfs) == 31
            and not mismatches
            and not summary_mismatches
        ),
        "input_sha256": {
            "homer_zip": sha256(source / "capsule-7227095-HOMER.zip"),
            "mapping": sha256(mapping_path),
            "published_intersection": sha256(expected_path),
            "published_lor_summary": sha256(author_summary_path),
            **{f"published_{database.lower()}_combined": sha256(path) for database, path in published_paths.items()},
        },
    }
    receipt_path = audit / "tnbc_author_motif_replay_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
