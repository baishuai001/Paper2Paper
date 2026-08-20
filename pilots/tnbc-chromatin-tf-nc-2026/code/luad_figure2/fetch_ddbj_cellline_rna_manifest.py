#!/usr/bin/env python3
"""Resolve paired baseline RNA-seq libraries for the DDBJ LUAD ATAC panel.

DRA006875--DRA006902 are the RNA-seq submission batches paired by design with
DRA006903--DRA006930.  The five early submissions are a dense perturbation set;
DRA006880--DRA006902 form the canonical 23-cell-line panel.  This script reads
official DDBJ experiment/run XML, selects DMSO baseline libraries, and can
restrict output to the exact Figure 2 cell-line manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


DDBJ_DRA = "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq"
SUBMISSIONS = [f"DRA{i:06d}" for i in range(6880, 6903)]
PRIMARY_BASELINE_SUBMISSION = "DRA001846"
CELL_LINE_RE = re.compile(r"^RNA-seq_(.+?)_([0-9.]+h)_", re.IGNORECASE)
PRIMARY_CELL_LINE_RE = re.compile(r"^(.+?)\s+\[RNAseq\]$", re.IGNORECASE)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, cache: Path, attempts: int = 5) -> bytes:
    if cache.is_file() and cache.stat().st_size:
        return cache.read_bytes()
    cache.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "Paper2Paper-LUAD-state-audit/1.0"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            cache.write_bytes(payload)
            return payload
        except (OSError, urllib.error.URLError) as exc:
            last_error = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


def xml_url(submission: str, suffix: str) -> str:
    return f"{DDBJ_DRA}/{submission[:6]}/{submission}/{submission}.{suffix}.xml"


def read_eligible(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result: set[str] = set()
    for row in rows:
        if row.get("system") and row.get("system") != "cell_line":
            continue
        include = row.get("figure2_primary_include", row.get("primary_include", "TRUE"))
        if include.upper() != "TRUE":
            continue
        name = row.get("cell_line") or row.get("sample_id")
        if name:
            result.add(name)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("--eligible-manifest", type=Path)
    args = parser.parse_args()
    raw = args.run_root / "audit/manifests/dra_rna/xml"
    eligible = read_eligible(args.eligible_manifest)
    rows: list[dict[str, str]] = []

    for submission in SUBMISSIONS:
        experiments = ET.fromstring(
            fetch(xml_url(submission, "experiment"), raw / f"{submission}.experiment.xml")
        )
        runs = ET.fromstring(fetch(xml_url(submission, "run"), raw / f"{submission}.run.xml"))
        run_by_experiment: dict[str, list[str]] = {}
        for run in runs.findall("RUN"):
            ref = run.find("EXPERIMENT_REF")
            if ref is None:
                continue
            experiment = ref.attrib.get("accession", "")
            accession = run.attrib.get("accession", "")
            if experiment and accession:
                run_by_experiment.setdefault(experiment, []).append(accession)

        for experiment in experiments.findall("EXPERIMENT"):
            accession = experiment.attrib.get("accession", "")
            descriptor = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
            if descriptor is None:
                continue
            library_name = descriptor.findtext("LIBRARY_NAME") or ""
            match = CELL_LINE_RE.match(library_name)
            if not match or "dmso" not in library_name.lower():
                continue
            cell_line, duration = match.group(1), match.group(2)
            if eligible is not None and cell_line not in eligible:
                continue
            layout_node = descriptor.find("LIBRARY_LAYOUT")
            layout = list(layout_node)[0].tag if layout_node is not None and list(layout_node) else ""
            sample = experiment.find("./DESIGN/SAMPLE_DESCRIPTOR")
            rows.append(
                {
                    "submission": submission,
                    "experiment": accession,
                    "run": ";".join(sorted(run_by_experiment.get(accession, []))),
                    "biosample": sample.findtext("./IDENTIFIERS/PRIMARY_ID") if sample is not None else "",
                    "sra_sample": sample.attrib.get("accession", "") if sample is not None else "",
                    "cell_line": cell_line,
                    "duration": duration,
                    "library_name": library_name,
                    "library_strategy": descriptor.findtext("LIBRARY_STRATEGY") or "",
                    "library_source": descriptor.findtext("LIBRARY_SOURCE") or "",
                    "library_selection": descriptor.findtext("LIBRARY_SELECTION") or "",
                    "library_layout": layout,
                    "is_dmso_control": "TRUE",
                    "paired_atac_submission": f"DRA{int(submission[3:]) + 28:06d}",
                }
            )

    rows.sort(key=lambda row: row["cell_line"])
    observed = {row["cell_line"] for row in rows}
    dmso_missing = sorted(eligible - observed) if eligible is not None else []
    duplicate = sorted(name for name in observed if sum(row["cell_line"] == name for row in rows) != 1)
    if duplicate:
        raise RuntimeError(f"DMSO RNA manifest has duplicate baselines: {duplicate}")
    columns = [
        "submission", "experiment", "run", "biosample", "sra_sample", "cell_line",
        "duration", "library_name", "library_strategy", "library_source",
        "library_selection", "library_layout", "is_dmso_control", "paired_atac_submission",
    ]
    out = args.run_root / "audit/manifests/dra_rna/ddbj_rna_dmso_control_libraries.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    # The untreated 26-cell-line RNA panel used by the source publication is
    # the complete primary expression source for state classification.
    primary_experiments = ET.fromstring(
        fetch(
            xml_url(PRIMARY_BASELINE_SUBMISSION, "experiment"),
            raw / f"{PRIMARY_BASELINE_SUBMISSION}.experiment.xml",
        )
    )
    primary_runs = ET.fromstring(
        fetch(
            xml_url(PRIMARY_BASELINE_SUBMISSION, "run"),
            raw / f"{PRIMARY_BASELINE_SUBMISSION}.run.xml",
        )
    )
    primary_run_by_experiment: dict[str, list[str]] = {}
    for run in primary_runs.findall("RUN"):
        ref = run.find("EXPERIMENT_REF")
        if ref is None:
            continue
        experiment = ref.attrib.get("accession", "")
        accession = run.attrib.get("accession", "")
        if experiment and accession:
            primary_run_by_experiment.setdefault(experiment, []).append(accession)
    primary_rows: list[dict[str, str]] = []
    for experiment in primary_experiments.findall("EXPERIMENT"):
        accession = experiment.attrib.get("accession", "")
        descriptor = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
        if descriptor is None:
            continue
        library_name = descriptor.findtext("LIBRARY_NAME") or experiment.findtext("TITLE") or ""
        match = PRIMARY_CELL_LINE_RE.match(library_name)
        if not match:
            continue
        cell_line = match.group(1)
        if eligible is not None and cell_line not in eligible:
            continue
        layout_node = descriptor.find("LIBRARY_LAYOUT")
        layout = list(layout_node)[0].tag if layout_node is not None and list(layout_node) else ""
        sample = experiment.find("./DESIGN/SAMPLE_DESCRIPTOR")
        primary_rows.append(
            {
                "submission": PRIMARY_BASELINE_SUBMISSION,
                "experiment": accession,
                "run": ";".join(sorted(primary_run_by_experiment.get(accession, []))),
                "biosample": sample.findtext("./IDENTIFIERS/PRIMARY_ID") if sample is not None else "",
                "sra_sample": sample.attrib.get("accession", "") if sample is not None else "",
                "cell_line": cell_line,
                "duration": "untreated_baseline",
                "library_name": library_name,
                "library_strategy": descriptor.findtext("LIBRARY_STRATEGY") or "",
                "library_source": descriptor.findtext("LIBRARY_SOURCE") or "",
                "library_selection": descriptor.findtext("LIBRARY_SELECTION") or "",
                "library_layout": layout,
                "is_dmso_control": "NOT_APPLICABLE_UNTREATED_BASELINE",
                "paired_atac_submission": "model_identity_match_not_same_experimental_batch",
            }
        )
    primary_rows.sort(key=lambda row: row["cell_line"])
    primary_observed = {row["cell_line"] for row in primary_rows}
    primary_missing = sorted(eligible - primary_observed) if eligible is not None else []
    primary_duplicate = sorted(
        name for name in primary_observed if sum(row["cell_line"] == name for row in primary_rows) != 1
    )
    if primary_missing or primary_duplicate:
        raise RuntimeError(
            f"DRA001846 baseline manifest mismatch: missing={primary_missing}; duplicate={primary_duplicate}"
        )
    primary_out = args.run_root / "audit/manifests/dra_rna/ddbj_dra001846_baseline_libraries.tsv"
    with primary_out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(primary_rows)
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source": "official DDBJ experiment and run XML",
        "submissions": SUBMISSIONS,
        "selection": "canonical Dataset-2 DMSO baseline RNA-seq libraries",
        "eligible_manifest": str(args.eligible_manifest.resolve()) if args.eligible_manifest else None,
        "eligible_count": len(eligible) if eligible is not None else None,
        "selected_library_count": len(rows),
        "selected_cell_lines": [row["cell_line"] for row in rows],
        "dmso_missing_eligible_cell_lines": dmso_missing,
        "manifest_sha256": sha256(out),
        "primary_baseline_submission": PRIMARY_BASELINE_SUBMISSION,
        "primary_baseline_library_count": len(primary_rows),
        "primary_baseline_cell_lines": [row["cell_line"] for row in primary_rows],
        "primary_baseline_manifest_sha256": sha256(primary_out),
    }
    receipt_path = args.run_root / "audit/manifests/dra_rna/ddbj_rna_manifest_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
