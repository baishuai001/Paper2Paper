#!/usr/bin/env python3
"""Parse one ataqv record per qualified raw ATAC sample into an audit table."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def load_json_maybe_gzip(path: Path) -> Any:
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw)


def finite_or_na(value: object) -> object:
    if value is None:
        return "NA"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return number if math.isfinite(number) else "NA"


def fraction(numerator: object, denominator: object) -> object:
    try:
        numerator_float = float(numerator)
        denominator_float = float(denominator)
    except (TypeError, ValueError):
        return "NA"
    if denominator_float <= 0:
        return "NA"
    return numerator_float / denominator_float


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    manifest_path = root / "audit/raw_atac_qc/figure2_qualified_raw_atac_samples.tsv"
    manifest = read_tsv(manifest_path)
    if not manifest:
        raise RuntimeError("No qualified raw ATAC samples were available for ataqv parsing")
    tn5_path = root / "audit/tn5_tss/figure2_tn5_tss_enrichment.tsv"
    tn5_rows = read_tsv(tn5_path)
    tn5_by_sample = {(row["system"], row["sample_id"]): row for row in tn5_rows}
    if len(tn5_by_sample) != len(manifest):
        raise RuntimeError(f"Expected {len(manifest)} layout-independent Tn5 TSS records, observed {len(tn5_by_sample)}")

    rows: list[dict[str, object]] = []
    errors: list[str] = []
    for sample in manifest:
        system = sample["system"]
        sample_id = sample["sample_id"]
        slug = sample["sample_slug"]
        json_path = root / "audit/ataqv" / system / slug / "metrics.ataqv.json.gz"
        parse_status = "PASS"
        parse_error = ""
        metrics: dict[str, Any] = {}
        ataqv_version = ""
        if not json_path.is_file() or json_path.stat().st_size == 0:
            parse_status = "FAIL"
            parse_error = "missing_or_empty_json"
        else:
            try:
                payload = load_json_maybe_gzip(json_path)
                if not isinstance(payload, list) or len(payload) != 1 or not isinstance(payload[0], dict):
                    observed = f"{type(payload).__name__}/{len(payload) if isinstance(payload, list) else 'NA'}"
                    raise ValueError(f"expected one top-level record, observed {observed}")
                record = payload[0]
                ataqv_version = str(record.get("ataqv_version", ""))
                candidate = record.get("metrics")
                if not isinstance(candidate, dict):
                    raise ValueError("metrics is not a JSON object")
                metrics = candidate
                if str(metrics.get("name", "")) != sample_id:
                    raise ValueError(f"sample name mismatch: {metrics.get('name', '')!r} != {sample_id!r}")
            except Exception as exc:  # retain a compact failure row for auditability
                parse_status = "FAIL"
                parse_error = f"{type(exc).__name__}:{exc}"
        if parse_status != "PASS":
            errors.append(f"{system}/{sample_id}:{parse_error}")

        total_reads = metrics.get("total_reads")
        duplicate_reads = metrics.get("duplicate_reads")
        mitochondrial_reads = metrics.get("total_mitochondrial_reads")
        hqaa_in_peaks = metrics.get("hqaa_in_peaks")
        hqaa = metrics.get("hqaa")
        ataqv_frip = fraction(hqaa_in_peaks, hqaa)
        worker_frip = finite_or_na(sample.get("frip"))
        frip_difference: object = "NA"
        if isinstance(ataqv_frip, float) and isinstance(worker_frip, float):
            frip_difference = ataqv_frip - worker_frip
        tn5 = tn5_by_sample.get((system, sample_id), {})
        if tn5.get("status") != "PASS":
            parse_status = "FAIL"
            parse_error = (parse_error + ";" if parse_error else "") + "missing_or_failed_tn5_tss_metric"
            errors.append(f"{system}/{sample_id}:missing_or_failed_tn5_tss_metric")

        rows.append(
            {
                "system": system,
                "sample_id": sample_id,
                "sample_slug": slug,
                "run": sample["run"],
                "ataqv_version": ataqv_version,
                "total_reads": total_reads if total_reads is not None else "NA",
                "high_quality_autosomal_alignments": hqaa if hqaa is not None else "NA",
                "high_quality_alignments_in_peaks": hqaa_in_peaks if hqaa_in_peaks is not None else "NA",
                "ataqv_frip": ataqv_frip,
                "ataqv_frip_percent_reported": finite_or_na(metrics.get("hqaa_overlapping_peaks_percent")),
                "worker_frip": worker_frip,
                "ataqv_minus_worker_frip": frip_difference,
                "tss_enrichment": finite_or_na(metrics.get("tss_enrichment")),
                "ataqv_tss_applicable": str(sample["layout"] == "PAIRED").upper(),
                "tn5_tss_enrichment_exact_center": tn5.get("tn5_tss_enrichment_exact_center", "NA"),
                "tn5_tss_enrichment_max_central_101bp": tn5.get("tn5_tss_enrichment_max_central_101bp", "NA"),
                "tn5_tss_enrichment_mean_central_101bp": tn5.get("tn5_tss_enrichment_mean_central_101bp", "NA"),
                "duplicate_reads": duplicate_reads if duplicate_reads is not None else "NA",
                "duplicate_fraction": fraction(duplicate_reads, total_reads),
                "mitochondrial_reads": mitochondrial_reads if mitochondrial_reads is not None else "NA",
                "mitochondrial_fraction": fraction(mitochondrial_reads, total_reads),
                "short_mononucleosomal_ratio": finite_or_na(metrics.get("short_mononucleosomal_ratio")),
                "mean_mapq": finite_or_na(metrics.get("mean_mapq")),
                "median_mapq": finite_or_na(metrics.get("median_mapq")),
                "total_peaks": metrics.get("total_peaks", "NA"),
                "total_peak_territory": metrics.get("total_peak_territory", "NA"),
                "parse_status": parse_status,
                "parse_error": parse_error,
                "metrics_json": str(json_path),
                "metrics_json_sha256": sha256(json_path) if json_path.is_file() else "",
            }
        )

    rows.sort(key=lambda item: (str(item["system"]), str(item["sample_id"]).lower()))
    out_dir = root / "audit/ataqv"
    table_path = out_dir / "figure2_ataqv_metrics.tsv"
    write_tsv(table_path, rows, list(rows[0]))
    tss_reported = sum(row["tss_enrichment"] != "NA" for row in rows)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "qualified_manifest": str(manifest_path),
        "qualified_manifest_sha256": sha256(manifest_path),
        "metrics_table": str(table_path),
        "metrics_table_sha256": sha256(table_path),
        "qualified_samples": len(manifest),
        "parsed_samples": sum(row["parse_status"] == "PASS" for row in rows),
        "tss_enrichment_reported_samples": tss_reported,
        "layout_independent_tn5_tss_reported_samples": sum(
            row["tn5_tss_enrichment_max_central_101bp"] != "NA" for row in rows
        ),
        "tss_enrichment_is_report_only_not_a_gate": True,
        "parse_status": "PASS" if not errors else "FAIL",
        "errors": errors,
    }
    receipt_path = out_dir / "figure2_ataqv_metrics_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
