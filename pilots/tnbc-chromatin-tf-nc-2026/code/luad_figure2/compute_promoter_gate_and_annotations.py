#!/usr/bin/env python3
"""Compute the frozen Figure 2 promoter gate and Supplementary Figure 3 annotations.

The script intentionally stops before motif analysis.  It emits the complete
sample-by-TF accessibility evidence used by Figure 2B/C and Supplementary
Figure 4A/B, plus primary/legend-window and TCGA-threshold sensitivities.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


FROZEN_PARENT_TABLE_SHA256 = "6f449276a97b8c66cbf93bd2a4ca6d51095fa426015dc9046a7c084431d5037b"
EXPECTED_TF_COUNT = 158
SYSTEM_ACTIVITY_COHORT = {
    "patient": "TCGA",
    "PDX": "PDMR_PDX",
    "cell_line": "DepMap_22Q2",
}
SYSTEM_EXPECTED = {"patient": 22, "PDX": 13, "cell_line": 19}
PATIENT_DEFINITIONS = (
    "cpm0.5_both_reps",
    "cpm1_both_reps",
    "cpm2_both_reps",
    "cpm1_merged_reps",
)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_attributes(value: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)"', value))


Interval = tuple[int, int]
IntervalMap = dict[str, list[Interval]]


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def merge_interval_map(values: dict[str, list[Interval]]) -> IntervalMap:
    return {chrom: merge_intervals(intervals) for chrom, intervals in values.items()}


def parse_gencode(
    gtf_path: Path,
    frozen_tfs: set[str],
) -> tuple[
    dict[str, dict[str, IntervalMap]],
    dict[str, IntervalMap],
    dict[str, set[str]],
]:
    promoters: dict[str, dict[str, dict[str, list[Interval]]]] = {
        "anchor_-2500_+1000": defaultdict(lambda: defaultdict(list)),
        "legend_-1000_+100": defaultdict(lambda: defaultdict(list)),
    }
    global_promoter: dict[str, list[Interval]] = defaultdict(list)
    exons: dict[str, list[Interval]] = defaultdict(list)
    genes: dict[str, list[Interval]] = defaultdict(list)
    frozen_gene_types: dict[str, set[str]] = defaultdict(set)

    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, _source, feature, start_text, end_text, _score, strand, _frame, attr_text = fields
            if chrom not in {f"chr{i}" for i in range(1, 23)} | {"chrX"}:
                continue
            attrs = parse_attributes(attr_text)
            gene_type = attrs.get("gene_type") or attrs.get("gene_biotype")
            start = int(start_text) - 1
            end = int(end_text)
            symbol = attrs.get("gene_name", "")

            # The frozen Figure 1 regulator list is PAN-GO-derived and contains
            # three explicitly named pseudogene regulators. The anchor's
            # promoter-accessibility rule is a locus/TSS rule, not a
            # protein-coding filter. Map every exact frozen symbol with a
            # GENCODE transcript, while retaining protein-coding genes as the
            # background feature universe used for peak annotation.
            if feature == "transcript" and symbol in frozen_tfs:
                frozen_gene_types[symbol].add(gene_type or "UNKNOWN")
                tss = start if strand == "+" else end - 1
                if strand == "+":
                    anchor = (max(0, tss - 2500), tss + 1000)
                    legend = (max(0, tss - 1000), tss + 100)
                else:
                    anchor = (max(0, tss - 1000), tss + 2500)
                    legend = (max(0, tss - 100), tss + 1000)
                promoters["anchor_-2500_+1000"][symbol][chrom].append(anchor)
                promoters["legend_-1000_+100"][symbol][chrom].append(legend)

            if gene_type != "protein_coding":
                continue
            if feature == "gene":
                genes[chrom].append((start, end))
            elif feature == "exon":
                exons[chrom].append((start, end))
            elif feature == "transcript":
                tss = start if strand == "+" else end - 1
                if strand == "+":
                    anchor = (max(0, tss - 2500), tss + 1000)
                    legend = (max(0, tss - 1000), tss + 100)
                else:
                    anchor = (max(0, tss - 1000), tss + 2500)
                global_promoter[chrom].append(anchor)

    frozen: dict[str, dict[str, IntervalMap]] = {}
    for definition, by_tf in promoters.items():
        frozen[definition] = {
            tf: merge_interval_map(by_tf.get(tf, {})) for tf in sorted(frozen_tfs)
        }
    features = {
        "promoter": merge_interval_map(global_promoter),
        "exon": merge_interval_map(exons),
        "gene": merge_interval_map(genes),
    }
    return frozen, features, frozen_gene_types


class IntervalIndex:
    def __init__(self, intervals: IntervalMap):
        self.intervals = intervals
        self.starts = {chrom: [item[0] for item in values] for chrom, values in intervals.items()}
        self.prefix_max_ends: dict[str, list[int]] = {}
        for chrom, values in intervals.items():
            maxima: list[int] = []
            current = -1
            for _start, end in values:
                current = max(current, end)
                maxima.append(current)
            self.prefix_max_ends[chrom] = maxima

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        starts = self.starts.get(chrom)
        if not starts:
            return False
        index = bisect.bisect_left(starts, end) - 1
        return index >= 0 and self.prefix_max_ends[chrom][index] > start


def read_peaks(path: Path) -> list[tuple[str, int, int]]:
    peaks: list[tuple[str, int, int]] = []
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            if chrom in {"chrM", "chrY"} or not chrom.startswith("chr"):
                continue
            peaks.append((chrom, int(fields[1]), int(fields[2])))
    return peaks


def promoter_accessibility(peaks: list[tuple[str, int, int]], promoter_maps: dict[str, IntervalMap]) -> dict[str, bool]:
    peak_values: dict[str, list[Interval]] = defaultdict(list)
    for chrom, start, end in peaks:
        peak_values[chrom].append((start, end))
    peak_index = IntervalIndex(merge_interval_map(peak_values))
    output: dict[str, bool] = {}
    for tf, chrom_map in promoter_maps.items():
        output[tf] = any(
            peak_index.overlaps(chrom, start, end)
            for chrom, intervals in chrom_map.items()
            for start, end in intervals
        )
    return output


def annotate_peaks(peaks: list[tuple[str, int, int]], features: dict[str, IntervalMap]) -> dict[str, int]:
    indexes = {name: IntervalIndex(values) for name, values in features.items()}
    counts = {"promoter": 0, "exonic": 0, "intronic": 0, "distal": 0}
    for chrom, start, end in peaks:
        if indexes["promoter"].overlaps(chrom, start, end):
            counts["promoter"] += 1
        elif indexes["exon"].overlaps(chrom, start, end):
            counts["exonic"] += 1
        elif indexes["gene"].overlaps(chrom, start, end):
            counts["intronic"] += 1
        else:
            counts["distal"] += 1
    return counts


def sample_manifests(root: Path) -> tuple[dict[str, dict[str, Path]], list[dict[str, str]]]:
    patient_table = read_tsv(root / "audit/tcga_luad_accessibility/tcga_luad_peak_file_manifest.tsv")
    patients: dict[str, dict[str, Path]] = defaultdict(dict)
    for row in patient_table:
        patients[row["definition"]][row["sample_id"]] = Path(row["peak_path"])
    for definition in PATIENT_DEFINITIONS:
        if len(patients[definition]) != 22:
            raise RuntimeError(f"Patient definition {definition} has {len(patients[definition])}, expected 22")

    raw_rows = read_tsv(root / "audit/raw_atac_qc/figure2_analysis_raw_atac_samples.tsv")
    for system in ("PDX", "cell_line"):
        count = sum(row["system"] == system for row in raw_rows)
        if count != SYSTEM_EXPECTED[system]:
            raise RuntimeError(f"{system} analysis n={count}, expected {SYSTEM_EXPECTED[system]}")
    return patients, raw_rows


def load_activity(project_root: Path, frozen_tfs: list[str]) -> dict[tuple[str, str], float | None]:
    path = project_root / "pilots/tnbc-chromatin-tf-nc-2026/execution/luad-figure1/results/tables/cross_system_TF_effects.tsv.gz"
    rows = read_tsv(path)
    lookup: dict[tuple[str, str], float | None] = {}
    for row in rows:
        key = (row["cohort"], row["TF"])
        value = row.get("mean_LUAD", "")
        try:
            lookup[key] = float(value) if value and value.upper() != "NA" else None
        except ValueError:
            lookup[key] = None
    for cohort in SYSTEM_ACTIVITY_COHORT.values():
        missing = [tf for tf in frozen_tfs if (cohort, tf) not in lookup]
        # The anchor explicitly retained a "Not Present" activity category.
        # Preserve it as missing/failed activity rather than silently imputing.
        for tf in missing:
            lookup[(cohort, tf)] = None
    return lookup


def evaluate_definition(
    definition_name: str,
    promoter_maps: dict[str, IntervalMap],
    patient_paths: dict[str, Path],
    raw_rows: list[dict[str, str]],
    frozen_tfs: list[str],
    activity: dict[tuple[str, str], float | None],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    samples: list[tuple[str, str, Path, bool]] = [
        ("patient", sample_id, path, False) for sample_id, path in sorted(patient_paths.items())
    ]
    samples.extend(
        (
            row["system"],
            row["sample_id"],
            Path(row["narrowpeak"]),
            row.get("peak_call_status") == "NO_SIGNIFICANT_PEAKS_AT_Q0.01",
        )
        for row in raw_rows
    )
    accessibility_rows: list[dict[str, object]] = []
    sample_cache: dict[Path, list[tuple[str, int, int]]] = {}
    for system, sample_id, peak_path, allow_empty in samples:
        if not peak_path.is_file():
            raise RuntimeError(f"Missing peak file for {system}/{sample_id}: {peak_path}")
        if peak_path.stat().st_size == 0 and not allow_empty:
            raise RuntimeError(
                f"Unexpected empty peak file without a q=0.01 zero-peak receipt "
                f"for {system}/{sample_id}: {peak_path}"
            )
        peaks = sample_cache.setdefault(peak_path, read_peaks(peak_path))
        calls = promoter_accessibility(peaks, promoter_maps)
        for tf in frozen_tfs:
            accessibility_rows.append(
                {
                    "analysis_definition": definition_name,
                    "system": system,
                    "sample_id": sample_id,
                    "TF": tf,
                    "promoter_accessible": str(calls[tf]).upper(),
                }
            )

    system_summary: list[dict[str, object]] = []
    for system in ("patient", "PDX", "cell_line"):
        system_samples = sorted({row["sample_id"] for row in accessibility_rows if row["system"] == system})
        threshold = math.ceil(len(system_samples) / 2)
        for tf in frozen_tfs:
            tf_rows = [row for row in accessibility_rows if row["system"] == system and row["TF"] == tf]
            accessible = sum(row["promoter_accessible"] == "TRUE" for row in tf_rows)
            cohort = SYSTEM_ACTIVITY_COHORT[system]
            mean_nes = activity[(cohort, tf)]
            system_summary.append(
                {
                    "analysis_definition": definition_name,
                    "system": system,
                    "activity_cohort": cohort,
                    "TF": tf,
                    "n_samples": len(system_samples),
                    "half_sample_threshold": threshold,
                    "n_promoter_accessible": accessible,
                    "fraction_promoter_accessible": accessible / len(system_samples),
                    "system_promoter_accessible": str(accessible >= threshold).upper(),
                    "mean_LUAD_NES": "" if mean_nes is None else mean_nes,
                    "mean_NES_nonnegative": str(mean_nes is not None and mean_nes >= 0).upper(),
                }
            )

    tf_summary: list[dict[str, object]] = []
    for tf in frozen_tfs:
        values = [row for row in system_summary if row["TF"] == tf]
        triple_open = all(row["system_promoter_accessible"] == "TRUE" for row in values)
        mean_nes_values = [
            None if row["mean_LUAD_NES"] == "" else float(row["mean_LUAD_NES"])
            for row in values
        ]
        all_three_negative = all(value is not None and value < 0 for value in mean_nes_values)
        retained_by_activity = not all_three_negative
        hc_tf = triple_open and retained_by_activity
        tf_summary.append(
            {
                "analysis_definition": definition_name,
                "TF": tf,
                "patient_promoter_accessible": next(row["system_promoter_accessible"] for row in values if row["system"] == "patient"),
                "PDX_promoter_accessible": next(row["system_promoter_accessible"] for row in values if row["system"] == "PDX"),
                "cell_line_promoter_accessible": next(row["system_promoter_accessible"] for row in values if row["system"] == "cell_line"),
                "triple_system_promoter_accessible": str(triple_open).upper(),
                "all_three_system_mean_NES_negative": str(all_three_negative).upper(),
                "retained_by_anchor_activity_rule": str(retained_by_activity).upper(),
                "HC_TF_promoter_activity_definition": str(hc_tf).upper(),
                "exclusion_reason": "" if hc_tf else (
                    "promoter_not_open_all_systems" if not triple_open else "mean_NES_negative_in_patient_PDX_and_cell_line"
                ),
            }
        )
    return accessibility_rows, system_summary, tf_summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("project_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    project_root = args.project_root.resolve()
    out_dir = root / "results/promoter_gate"
    audit_dir = root / "audit/promoter_gate"
    out_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    tf_path = project_root / "pilots/tnbc-chromatin-tf-nc-2026/execution/luad-figure1/results/tables/tcga_specific_TFs.tsv"
    if sha256(tf_path) != FROZEN_PARENT_TABLE_SHA256:
        raise RuntimeError("Frozen Figure 1 discovery-TF parent table SHA256 mismatch")
    tf_rows = read_tsv(tf_path)
    frozen_rows = sorted(
        (
            {
                "TF": row["TF"],
                "discovery_category": row["discovery_category"],
                "externally_replicated": row["externally_replicated"],
            }
            for row in tf_rows
            if row.get("discovery_category") == "LUAD"
        ),
        key=lambda row: row["TF"],
    )
    frozen_input_path = audit_dir / "frozen_figure2_luad_tf_input.tsv"
    write_tsv(
        frozen_input_path,
        frozen_rows,
        ["TF", "discovery_category", "externally_replicated"],
    )
    frozen_tfs = [row["TF"] for row in frozen_rows]
    if len(frozen_tfs) != EXPECTED_TF_COUNT or len(set(frozen_tfs)) != EXPECTED_TF_COUNT:
        raise RuntimeError(f"Expected 158 unique frozen TFs, observed {len(set(frozen_tfs))}")

    gtf_path = root / "reference/downloads/gencode.v47.basic.annotation.gtf.gz"
    promoters, features, frozen_gene_types = parse_gencode(gtf_path, set(frozen_tfs))
    absent = [tf for tf in frozen_tfs if not promoters["anchor_-2500_+1000"].get(tf)]
    mapping_rows = [
        {
            "TF": tf,
            "gencode_v47_transcript_present": str(tf not in absent).upper(),
            "gencode_v47_gene_types": ";".join(sorted(frozen_gene_types.get(tf, set()))),
            "protein_coding_only": str(frozen_gene_types.get(tf, set()) == {"protein_coding"}).upper(),
            "anchor_promoter_interval_count": sum(len(values) for values in promoters["anchor_-2500_+1000"].get(tf, {}).values()),
            "legend_promoter_interval_count": sum(len(values) for values in promoters["legend_-1000_+100"].get(tf, {}).values()),
        }
        for tf in frozen_tfs
    ]
    write_tsv(audit_dir / "frozen_tf_gencode_mapping.tsv", mapping_rows)
    if absent:
        raise RuntimeError(f"Cannot map frozen regulators to any GENCODE v47 transcript: {absent}")
    non_protein_coding = [row for row in mapping_rows if row["protein_coding_only"] != "TRUE"]
    (audit_dir / "frozen_tf_gencode_mapping_receipt.json").write_text(
        json.dumps(
            {
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "frozen_regulators": len(frozen_tfs),
                "mapped_to_gencode_v47_transcript": len(frozen_tfs) - len(absent),
                "non_protein_coding_regulators": [
                    {"TF": row["TF"], "gene_types": row["gencode_v47_gene_types"]}
                    for row in non_protein_coding
                ],
                "rule": (
                    "Use exact-symbol GENCODE transcripts for the anchor promoter locus rule; "
                    "retain protein-coding transcripts as the genome-wide peak-annotation universe"
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    feature_dir = root / "reference/figure2_features"
    feature_dir.mkdir(parents=True, exist_ok=True)
    for definition, by_tf in promoters.items():
        rows: list[dict[str, object]] = []
        for tf, chrom_map in by_tf.items():
            for chrom, intervals in chrom_map.items():
                rows.extend({"chrom": chrom, "start": start, "end": end, "TF": tf} for start, end in intervals)
        write_tsv(feature_dir / f"frozen_tf_promoters_{definition}.tsv", rows, ["chrom", "start", "end", "TF"])

    patients, raw_rows = sample_manifests(root)
    activity = load_activity(project_root, frozen_tfs)
    all_accessibility: list[dict[str, object]] = []
    all_system_summary: list[dict[str, object]] = []
    all_tf_summary: list[dict[str, object]] = []
    configurations = [
        ("PRIMARY_anchor_promoter_tcga_cpm1_both", "anchor_-2500_+1000", "cpm1_both_reps"),
        ("SENS_legend_promoter_tcga_cpm1_both", "legend_-1000_+100", "cpm1_both_reps"),
        ("SENS_anchor_promoter_tcga_cpm0.5_both", "anchor_-2500_+1000", "cpm0.5_both_reps"),
        ("SENS_anchor_promoter_tcga_cpm2_both", "anchor_-2500_+1000", "cpm2_both_reps"),
        ("SENS_anchor_promoter_tcga_cpm1_merged", "anchor_-2500_+1000", "cpm1_merged_reps"),
    ]
    for name, promoter_definition, patient_definition in configurations:
        accessibility, system_summary, tf_summary = evaluate_definition(
            name,
            promoters[promoter_definition],
            patients[patient_definition],
            raw_rows,
            frozen_tfs,
            activity,
        )
        all_accessibility.extend(accessibility)
        all_system_summary.extend(system_summary)
        all_tf_summary.extend(tf_summary)

    write_tsv(out_dir / "sample_tf_promoter_accessibility.tsv", all_accessibility)
    write_tsv(out_dir / "system_tf_promoter_activity_summary.tsv", all_system_summary)
    write_tsv(out_dir / "tf_promoter_gate_summary.tsv", all_tf_summary)

    primary_name = configurations[0][0]
    primary_tf_rows = [row for row in all_tf_summary if row["analysis_definition"] == primary_name]
    sensitivity_rows: list[dict[str, object]] = []
    for name, _promoter, _patient in configurations:
        rows = [row for row in all_tf_summary if row["analysis_definition"] == name]
        sensitivity_rows.append(
            {
                "analysis_definition": name,
                "frozen_TFs": len(rows),
                "triple_system_promoter_accessible_TFs": sum(row["triple_system_promoter_accessible"] == "TRUE" for row in rows),
                "HC_TFs_after_activity_filter": sum(row["HC_TF_promoter_activity_definition"] == "TRUE" for row in rows),
            }
        )
    write_tsv(out_dir / "promoter_gate_sensitivity_summary.tsv", sensitivity_rows)

    annotation_rows: list[dict[str, object]] = []
    primary_samples: list[tuple[str, str, Path]] = [
        ("patient", sample_id, path) for sample_id, path in sorted(patients["cpm1_both_reps"].items())
    ]
    primary_samples.extend(
        (row["system"], row["sample_id"], Path(row["narrowpeak"])) for row in raw_rows
    )
    for system, sample_id, peak_path in primary_samples:
        peaks = read_peaks(peak_path)
        counts = annotate_peaks(peaks, features)
        total = len(peaks)
        for category in ("promoter", "exonic", "intronic", "distal"):
            annotation_rows.append(
                {
                    "system": system,
                    "sample_id": sample_id,
                    "category": category,
                    "peak_count": counts[category],
                    "total_peaks": total,
                    "fraction": counts[category] / total if total else "",
                    "hierarchical_definition": "promoter_then_exonic_then_intronic_then_distal",
                }
            )
    write_tsv(root / "results/supplementary_figure3/peak_genomic_annotation.tsv", annotation_rows)

    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_tf_parent_file": str(tf_path),
        "frozen_tf_parent_sha256": sha256(tf_path),
        "frozen_tf_selection_rule": "discovery_category == LUAD (all Figure 1 discovery TFs)",
        "frozen_tf_file": str(frozen_input_path),
        "frozen_tf_sha256": sha256(frozen_input_path),
        "frozen_tf_count": len(frozen_tfs),
        "gencode": "v47 basic exact-symbol candidate transcripts; protein-coding genome-wide annotation universe",
        "gencode_sha256": sha256(gtf_path),
        "primary_definition": primary_name,
        "patient_n": SYSTEM_EXPECTED["patient"],
        "PDX_n": sum(row["system"] == "PDX" for row in raw_rows),
        "cell_line_n": sum(row["system"] == "cell_line" for row in raw_rows),
        "triple_system_promoter_accessible_TFs": sum(row["triple_system_promoter_accessible"] == "TRUE" for row in primary_tf_rows),
        "HC_TFs_after_promoter_and_activity": sum(row["HC_TF_promoter_activity_definition"] == "TRUE" for row in primary_tf_rows),
        "motif_gate_not_yet_evaluated": True,
        "output_tables": [
            str(out_dir / "sample_tf_promoter_accessibility.tsv"),
            str(out_dir / "system_tf_promoter_activity_summary.tsv"),
            str(out_dir / "tf_promoter_gate_summary.tsv"),
            str(out_dir / "promoter_gate_sensitivity_summary.tsv"),
            str(root / "results/supplementary_figure3/peak_genomic_annotation.tsv"),
        ],
    }
    (audit_dir / "promoter_gate_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
