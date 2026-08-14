#!/usr/bin/env python3
"""Build an auditable DDBJ ATAC-seq cell-line manifest on the cloud server.

DRA006903--DRA006930 are submission batches, not independent samples.  This
script resolves each batch through the official DDBJ Search API, identifies the
cell line from experiment library names, and optionally expands all experiments
to select untreated/DMSO baseline libraries.  Every API response is cached.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


API = "https://ddbj.nig.ac.jp/search/api/entries"
DDBJ_DRA = "https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq"
SUBMISSIONS = [f"DRA{i:06d}" for i in range(6903, 6931)]
CELL_LINE_RE = re.compile(r"^ATAC-seq_(.+?)_([0-9.]+h)_", re.IGNORECASE)


def fetch_json(url: str, cache_path: Path, attempts: int = 5) -> dict:
    if cache_path.exists() and cache_path.stat().st_size:
        return json.loads(cache_path.read_text(encoding="utf-8"))
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "Paper2Paper-LUAD-Figure2-manifest/1.0"},
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
            parsed = json.loads(payload)
            cache_path.write_bytes(payload)
            return parsed
        except (OSError, urllib.error.URLError, json.JSONDecodeError) as exc:
            last_error = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


def fetch_bytes(url: str, cache_path: Path, attempts: int = 5) -> bytes:
    if cache_path.exists() and cache_path.stat().st_size:
        return cache_path.read_bytes()
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "Paper2Paper-LUAD-Figure2-manifest/1.0"},
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                payload = response.read()
            cache_path.write_bytes(payload)
            return payload
        except (OSError, urllib.error.URLError) as exc:
            last_error = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"Failed to fetch {url}: {last_error}")


def xrefs(record: dict, kind: str) -> list[str]:
    return sorted(
        x["identifier"]
        for x in record.get("dbXrefs", [])
        if x.get("type") == kind and x.get("identifier")
    )


def parse_library_name(name: str | None) -> tuple[str | None, str | None]:
    if not name:
        return None, None
    match = CELL_LINE_RE.search(name)
    if not match:
        return None, None
    return match.group(1), match.group(2)


def is_baseline(name: str | None) -> bool:
    text = (name or "").lower()
    # Dataset 1 names contain "DMSO (Control)" whereas Dataset 2 uses only
    # "DMSO".  Requiring the literal word Control silently drops all 23
    # canonical Dataset 2 baselines.
    return "dmso" in text


def dataset_name(submission: str) -> str:
    return "Dataset1_dense_5_lines" if int(submission[3:]) <= 6907 else "Dataset2_panel_23_lines"


def is_canonical_submission(submission: str) -> bool:
    # The five Dataset 1 lines recur in Dataset 2.  Dataset 2 is therefore the
    # one-row-per-biological-model panel; Dataset 1 is retained as a technical
    # batch sensitivity set and is never counted as five extra cell lines.
    return int(submission[3:]) >= 6908


def dra_xml_url(submission: str, suffix: str) -> str:
    prefix = submission[:6]
    return f"{DDBJ_DRA}/{prefix}/{submission}/{submission}.{suffix}.xml"


def parse_submission_xml(submission: str, raw: Path) -> list[dict]:
    experiment_payload = fetch_bytes(
        dra_xml_url(submission, "experiment"),
        raw / "submission_xml" / f"{submission}.experiment.xml",
    )
    run_payload = fetch_bytes(
        dra_xml_url(submission, "run"),
        raw / "submission_xml" / f"{submission}.run.xml",
    )
    run_by_experiment: dict[str, list[str]] = {}
    for run in ET.fromstring(run_payload).findall("RUN"):
        experiment_ref = run.find("EXPERIMENT_REF")
        if experiment_ref is None:
            continue
        experiment = experiment_ref.attrib.get("accession", "")
        accession = run.attrib.get("accession", "")
        if experiment and accession:
            run_by_experiment.setdefault(experiment, []).append(accession)

    rows: list[dict] = []
    for experiment in ET.fromstring(experiment_payload).findall("EXPERIMENT"):
        accession = experiment.attrib.get("accession", "")
        descriptor = experiment.find("./DESIGN/LIBRARY_DESCRIPTOR")
        sample = experiment.find("./DESIGN/SAMPLE_DESCRIPTOR")
        library_name = descriptor.findtext("LIBRARY_NAME") if descriptor is not None else ""
        cell_line, duration = parse_library_name(library_name)
        layout = ""
        if descriptor is not None:
            layout_node = descriptor.find("LIBRARY_LAYOUT")
            if layout_node is not None and list(layout_node):
                layout = list(layout_node)[0].tag
        platform_node = experiment.find("./PLATFORM")
        platform = ""
        instrument = ""
        if platform_node is not None and list(platform_node):
            platform = list(platform_node)[0].tag
            instrument = list(platform_node)[0].findtext("INSTRUMENT_MODEL") or ""
        biosample = ""
        sra_sample = ""
        if sample is not None:
            sra_sample = sample.attrib.get("accession", "")
            biosample = sample.findtext("./IDENTIFIERS/PRIMARY_ID") or ""
        rows.append(
            {
                "submission": submission,
                "dataset": dataset_name(submission),
                "canonical_submission": str(is_canonical_submission(submission)).upper(),
                "experiment": accession,
                "run": ";".join(sorted(run_by_experiment.get(accession, []))),
                "biosample": biosample,
                "sra_sample": sra_sample,
                "cell_line": cell_line or "",
                "duration": duration or "",
                "library_name": library_name or "",
                "library_strategy": descriptor.findtext("LIBRARY_STRATEGY") if descriptor is not None else "",
                "library_layout": layout,
                "platform": platform,
                "instrument_model": instrument,
                "is_dmso_control": str(is_baseline(library_name)).upper(),
                "status": "public",
                "accessibility": "open",
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def experiment_row(submission: str, accession: str, record: dict) -> dict:
    library_name = record.get("libraryName")
    cell_line, duration = parse_library_name(library_name)
    return {
        "submission": submission,
        "experiment": accession,
        "run": ";".join(xrefs(record, "sra-run")),
        "biosample": ";".join(xrefs(record, "biosample")),
        "sra_sample": ";".join(xrefs(record, "sra-sample")),
        "cell_line": cell_line or "",
        "duration": duration or "",
        "library_name": library_name or "",
        "library_strategy": ";".join(record.get("libraryStrategy", [])),
        "library_layout": record.get("libraryLayout") or "",
        "platform": record.get("platform") or "",
        "instrument_model": ";".join(record.get("instrumentModel", [])),
        "is_dmso_control": str(is_baseline(library_name)).upper(),
        "status": record.get("status") or "",
        "accessibility": record.get("accessibility") or "",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument(
        "--expand-baselines",
        action="store_true",
        help="Fetch all experiments and emit all DMSO/control baseline libraries.",
    )
    parser.add_argument(
        "--include-dataset1-controls",
        action="store_true",
        help="Also emit controls from the five dense Dataset 1 submissions as a sensitivity set.",
    )
    parser.add_argument(
        "--eligible-cell-lines",
        type=Path,
        help="Optional TSV with cell_line and figure2_primary_include columns.",
    )
    args = parser.parse_args()

    audit = args.run_root / "audit" / "manifests" / "dra"
    raw = audit / "api_json"
    identify_rows: list[dict] = []
    baseline_rows: list[dict] = []
    eligible: set[str] | None = None
    if args.eligible_cell_lines:
        with args.eligible_cell_lines.open(encoding="utf-8", newline="") as handle:
            eligible_rows = list(csv.DictReader(handle, delimiter="\t"))
        eligible = {
            row["cell_line"]
            for row in eligible_rows
            if row.get("figure2_primary_include", "").upper() == "TRUE"
        }

    for submission in SUBMISSIONS:
        submission_record = fetch_json(
            f"{API}/sra-submission/{submission}", raw / "submissions" / f"{submission}.json"
        )
        experiments = xrefs(submission_record, "sra-experiment")
        if not experiments:
            raise RuntimeError(f"No experiments linked from {submission}")

        first_accession = experiments[0]
        first_record = fetch_json(
            f"{API}/sra-experiment/{first_accession}",
            raw / "experiments" / f"{first_accession}.json",
        )
        first = experiment_row(submission, first_accession, first_record)
        identify_rows.append(
            {
                "submission": submission,
                "dataset": dataset_name(submission),
                "canonical_submission": str(is_canonical_submission(submission)).upper(),
                "experiment_count": len(experiments),
                "cell_line": first["cell_line"],
                "first_experiment": first_accession,
                "first_run": first["run"],
                "first_biosample": first["biosample"],
                "first_library_name": first["library_name"],
                "first_is_dmso_control": first["is_dmso_control"],
                "status": submission_record.get("status") or "",
                "accessibility": submission_record.get("accessibility") or "",
            }
        )

        if args.expand_baselines:
            if not is_canonical_submission(submission) and not args.include_dataset1_controls:
                continue
            for row in parse_submission_xml(submission, raw):
                if row["is_dmso_control"] != "TRUE":
                    continue
                if eligible is not None and row["cell_line"] not in eligible:
                    continue
                baseline_rows.append(row)

    identify_columns = [
        "submission",
        "dataset",
        "canonical_submission",
        "experiment_count",
        "cell_line",
        "first_experiment",
        "first_run",
        "first_biosample",
        "first_library_name",
        "first_is_dmso_control",
        "status",
        "accessibility",
    ]
    write_tsv(audit / "ddbj_atac_submission_identity.tsv", identify_rows, identify_columns)

    if args.expand_baselines:
        baseline_columns = [
            "submission",
            "dataset",
            "canonical_submission",
            "experiment",
            "run",
            "biosample",
            "sra_sample",
            "cell_line",
            "duration",
            "library_name",
            "library_strategy",
            "library_layout",
            "platform",
            "instrument_model",
            "is_dmso_control",
            "status",
            "accessibility",
        ]
        write_tsv(audit / "ddbj_atac_dmso_control_libraries.tsv", baseline_rows, baseline_columns)

    duplicate_lines = sorted(
        name
        for name in {row["cell_line"] for row in identify_rows}
        if name and sum(row["cell_line"] == name for row in identify_rows) > 1
    )
    receipt = {
        "submission_count": len(identify_rows),
        "unique_cell_line_count": len({r["cell_line"] for r in identify_rows if r["cell_line"]}),
        "unparsed_cell_line_submissions": [r["submission"] for r in identify_rows if not r["cell_line"]],
        "duplicate_cell_lines": duplicate_lines,
        "baseline_library_count": len(baseline_rows) if args.expand_baselines else None,
        "baseline_unique_cell_line_count": (
            len({row["cell_line"] for row in baseline_rows}) if args.expand_baselines else None
        ),
        "dataset1_controls_included": args.include_dataset1_controls,
        "eligible_cell_line_filter": str(args.eligible_cell_lines) if args.eligible_cell_lines else None,
        "source_api": API,
        "source_bulk_xml": DDBJ_DRA,
        "selection_rule": (
            "canonical Dataset 2 submission; library_name contains DMSO (case-insensitive); "
            "one biological model counted once; Dataset 1 duplicates are sensitivity-only"
        ),
    }
    (audit / "ddbj_atac_manifest_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
