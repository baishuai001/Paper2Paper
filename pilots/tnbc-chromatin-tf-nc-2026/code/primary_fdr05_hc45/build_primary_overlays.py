#!/usr/bin/env python3
"""Build isolated Figure 1/2 inputs for the FDR05 -> HC45 primary branch.

The script never mutates the completed 158/31 run.  It copies only small
tabular/reference inputs and symlinks immutable matrices, manifests, HOMER
results and sensitivity outputs from the source runs.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


PRIMARY = "PRIMARY_anchor_promoter_tcga_cpm1_both"
SCENARIO = "S3_activity_FDR05"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is None:
        if not rows:
            raise RuntimeError(f"Cannot infer an empty TSV schema: {path}")
        columns = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def numeric(row: dict[str, str], name: str) -> float:
    value = row.get(name, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def classify(row: dict[str, str]) -> str:
    fdr = numeric(row, "FDR")
    nes = numeric(row, "NES")
    logfc = numeric(row, "logFC")
    if fdr <= 0.05 and nes > 0 and logfc > 0:
        return "LUAD"
    if fdr <= 0.05 and nes < 0 and logfc < 0:
        return "LUSC"
    return "Non-Significant"


def recategorize(path: Path) -> tuple[list[dict[str, str]], dict[str, set[str]]]:
    rows = read_tsv(path)
    required = {"TF", "NES", "FDR", "logFC", "category"}
    if not rows or not required.issubset(rows[0]):
        raise RuntimeError(f"Unexpected msVIPER table schema: {path}")
    sets = {"LUAD": set(), "LUSC": set()}
    for row in rows:
        row["category"] = classify(row)
        if row["category"] in sets:
            sets[row["category"]].add(row["TF"].upper())
    return rows, sets


def symlink(source: Path, target: Path) -> None:
    if not source.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    os.symlink(source.resolve(), target, target_is_directory=source.is_dir())


def exact_memberships(sets: dict[str, set[str]]) -> list[dict[str, object]]:
    names = list(sets)
    universe = set().union(*sets.values())
    counts: Counter[str] = Counter()
    for tf in universe:
        present = [name for name in names if tf in sets[name]]
        counts["+".join(present)] += 1
    return [{"combination": key, "TFs": counts[key]} for key in sorted(counts)]


def build_figure1(source: Path, output: Path) -> dict[str, object]:
    source_tables = source / "results/tables"
    output_tables = output / "results/tables"
    output_tables.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source_tables, output_tables)
    symlink(source / "data", output / "data")
    for child in (source / "results").iterdir():
        if child.name != "tables":
            symlink(child, output / "results" / child.name)

    table_names = {
        "TCGA": "figure1b_tcga_limma_msviper.tsv",
        "GSE81089": "gse81089_limma_msviper.tsv",
        "GSE41271": "gse41271_limma_msviper.tsv",
    }
    category_sets: dict[str, dict[str, set[str]]] = {}
    for cohort, name in table_names.items():
        rows, sets = recategorize(source_tables / name)
        category_sets[cohort] = sets
        write_tsv(output_tables / name, rows)

    tcga = category_sets["TCGA"]
    gse = category_sets["GSE81089"]
    micro = category_sets["GSE41271"]
    membership_rows: list[dict[str, object]] = []
    replication_rows: list[dict[str, object]] = []
    for direction in ("LUAD", "LUSC"):
        for tf in sorted(tcga[direction]):
            gse_hit = tf in gse[direction]
            micro_hit = tf in micro[direction]
            replication_rows.append(
                {
                    "TF": tf,
                    "discovery_category": direction,
                    "replicated_GSE81089": str(gse_hit).upper(),
                    "replicated_GSE41271": str(micro_hit).upper(),
                    "replicated_both": str(gse_hit and micro_hit).upper(),
                }
            )
            membership_rows.append(
                {
                    "TF": tf,
                    "discovery_category": direction,
                    "externally_replicated": str(gse_hit or micro_hit).upper(),
                }
            )
    write_tsv(output_tables / "tcga_specific_TFs.tsv", membership_rows)
    write_tsv(output_tables / "figure1_cross_platform_replication.tsv", replication_rows)
    write_tsv(
        output_tables / "externally_replicated_TFs.tsv",
        [row for row in replication_rows if row["replicated_GSE81089"] == "TRUE"],
    )
    write_tsv(
        output_tables / "gse41271_replicated_TFs.tsv",
        [row for row in replication_rows if row["replicated_GSE41271"] == "TRUE"],
    )

    summary_rows = []
    three_way_rows = []
    pair_rows = []
    for direction in ("LUAD", "LUSC"):
        t, g, m = tcga[direction], gse[direction], micro[direction]
        summary_rows.append(
            {
                "specificity": direction,
                "TCGA_specific": len(t),
                "GSE81089_specific": len(g),
                "GSE41271_specific": len(m),
                "TCGA_replicated_GSE81089": len(t & g),
                "TCGA_replicated_GSE41271": len(t & m),
                "TCGA_replicated_both": len(t & g & m),
            }
        )
        combos = exact_memberships({"TCGA": t, "GSE81089": g, "GSE41271": m})
        three_way_rows.extend({"specificity": direction, **row} for row in combos)
        pair_rows.extend(
            [
                {"direction": direction, "component": "TCGA only", "n": len(t - g)},
                {"direction": direction, "component": "Overlap", "n": len(t & g)},
                {"direction": direction, "component": "GSE81089 only", "n": len(g - t)},
            ]
        )
    write_tsv(output_tables / "figure1_cross_platform_replication_summary.tsv", summary_rows)
    write_tsv(output_tables / "figure1_patient_cohort_TF_overlap_counts.tsv", three_way_rows)
    write_tsv(output_tables / "SupplementaryFigure2B_TF_overlap_counts.tsv", pair_rows)

    return {
        "rule": "BH FDR<=0.05 and NES/expression-logFC direction concordant",
        "counts": {
            cohort: {direction: len(values) for direction, values in sets.items()}
            for cohort, sets in category_sets.items()
        },
        "TCGA_LUAD_replicated_GSE81089": len(tcga["LUAD"] & gse["LUAD"]),
        "TCGA_LUAD_replicated_GSE41271": len(tcga["LUAD"] & micro["LUAD"]),
        "TCGA_LUAD_replicated_both": len(tcga["LUAD"] & gse["LUAD"] & micro["LUAD"]),
    }


def filter_scenario(path: Path, output: Path) -> tuple[int, set[str]]:
    rows = read_tsv(path)
    selected: list[dict[str, object]] = []
    for row in rows:
        if row.get("analysis_definition") != SCENARIO:
            continue
        row = dict(row)
        row["source_analysis_definition"] = SCENARIO
        row["analysis_definition"] = PRIMARY
        selected.append(row)
    if not selected:
        raise RuntimeError(f"No {SCENARIO} rows found in {path}")
    write_tsv(output, selected)
    hcs = {
        row["TF"]
        for row in selected
        if row.get("HC_TF_promoter_activity_definition", "").upper() == "TRUE"
    }
    return len({row["TF"] for row in selected}), hcs


def build_figure2(source: Path, entry: Path, output: Path) -> dict[str, object]:
    symlink(source / "data", output / "data")
    for name in ("tools", "cache"):
        symlink(source / name, output / name)

    motif_source = source / "reference/motifs"
    motif_output = output / "reference/motifs"
    shutil.copytree(motif_source, motif_output)

    promoter_out = output / "results/promoter_gate"
    promoter_out.mkdir(parents=True, exist_ok=True)
    candidate_n, hcs = filter_scenario(
        entry / "tf_promoter_gate_summary.tsv", promoter_out / "tf_promoter_gate_summary.tsv"
    )
    sample_n, _ = filter_scenario(
        entry / "sample_tf_promoter_accessibility.tsv",
        promoter_out / "sample_tf_promoter_accessibility.tsv",
    )
    system_n, _ = filter_scenario(
        entry / "system_tf_promoter_activity_summary.tsv",
        promoter_out / "system_tf_promoter_activity_summary.tsv",
    )
    if (candidate_n, len(hcs)) != (187, 45):
        raise RuntimeError(f"Expected 187 candidates and 45 HC-TFs; got {candidate_n}/{len(hcs)}")

    symlink(
        source / "results/motif_primary/anchor_consensus_homer",
        output / "results/motif_primary/anchor_consensus_homer",
    )
    symlink(source / "results/motif_equivalence", output / "results/motif_equivalence")
    audit = output / "audit/primary_fdr05_hc45"
    audit.mkdir(parents=True, exist_ok=True)
    write_tsv(audit / "HC45.tsv", [{"TF": tf} for tf in sorted(hcs)])
    write_tsv(
        audit / "entry_counts.tsv",
        [
            {
                "scenario": SCENARIO,
                "candidate_TFs": candidate_n,
                "HC_TFs": len(hcs),
                "sample_table_TFs": sample_n,
                "system_table_TFs": system_n,
            }
        ],
    )
    return {
        "source_scenario": SCENARIO,
        "compatibility_analysis_definition": PRIMARY,
        "candidate_TFs": candidate_n,
        "HC_TFs": len(hcs),
        "HC_TF_names": sorted(hcs),
        "source_full_JASPAR_sha256": sha256(source / "reference/motifs/full_reference_JASPAR2024.homer"),
        "source_full_CISBP_sha256": sha256(source / "reference/motifs/full_reference_CIS-BP2.00.homer"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure1-source", required=True, type=Path)
    parser.add_argument("--figure2-source", required=True, type=Path)
    parser.add_argument("--figure2-entry-audit", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()

    output = args.output_root.resolve()
    if output.exists():
        raise SystemExit(f"Refusing to overwrite audited primary branch: {output}")
    output.mkdir(parents=True)
    f1_out = output / "figure1"
    f2_out = output / "figure2"
    f1_out.mkdir()
    f2_out.mkdir()

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS",
        "analysis_id": "luad-fdr05-hc45-primary-v1",
        "figure1": build_figure1(args.figure1_source.resolve(), f1_out),
        "figure2_entry": build_figure2(
            args.figure2_source.resolve(), args.figure2_entry_audit.resolve(), f2_out
        ),
        "immutability": "source 158/31 runs are symlinked read-only by convention and never overwritten",
    }
    (output / "overlay_build_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
