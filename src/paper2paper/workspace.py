from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Any

from .schema import (
    MANUSCRIPT_STAGES,
    PILOT_STAGES,
    PROJECT_KEYS,
    PROJECT_OBJECT_KEYS,
    SCHEMA_VERSION,
    STAGES,
    TABLES,
)


TODO = "<!-- TODO -->"
DATA_READY = {"sample_parsed", "downloaded", "checksum_verified"}
CODE_READY = {"smoke_passed", "tested"}
PATH_READY = {"target_environment_passed", "staged_workaround"}
MODEL_READY = {"locked", "verified"}
FIGURE_RANK = {
    "idea": 0,
    "mapped": 1,
    "spike_generated": 2,
    "generated": 3,
    "verified": 4,
    "blocked": -1,
    "dropped": -1,
}
EVIDENCE_STAGE_RANK = {
    "direction_audited": 0,
    "availability_prechecked": 1,
    "minimal_real_run": 2,
    "figure_loop_closed": 3,
}


@dataclass
class ValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    return header, rows


def append_tsv(path: Path, values: list[str]) -> None:
    header, _ = read_tsv(path)
    if len(values) != len(header):
        raise ValueError(
            f"{path.name}: expected {len(header)} values, got {len(values)}"
        )
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(values)


def _write_tsv_rows(path: Path, rows: list[dict[str, str]]) -> None:
    header, _ = read_tsv(path)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_manifest(project_dir: Path) -> dict[str, Any]:
    path = Path(project_dir) / "PROJECT.json"
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("PROJECT.json root must be an object")
    return value


def load_tables(
    project_dir: Path,
) -> tuple[dict[str, list[dict[str, str]]], ValidationReport]:
    project_dir = Path(project_dir)
    report = ValidationReport()
    tables: dict[str, list[dict[str, str]]] = {}

    for name, spec in TABLES.items():
        path = project_dir / spec["path"]
        if not path.exists():
            report.errors.append(f"missing table: {spec['path']}")
            tables[name] = []
            continue
        header, rows = read_tsv(path)
        expected = list(spec["columns"])
        if header != expected:
            report.errors.append(
                f"{spec['path']}: header must exactly match the workspace contract"
            )
        for row_number, row in enumerate(rows, start=2):
            for field_name in spec["nonempty"]:
                if not row.get(field_name, ""):
                    report.errors.append(
                        f"{spec['path']}:{row_number}: empty {field_name}"
                    )
            for field_name, choices in spec["enums"].items():
                value = row.get(field_name, "")
                if value not in choices:
                    report.errors.append(
                        f"{spec['path']}:{row_number}: invalid "
                        f"{field_name}={value!r}"
                    )
        tables[name] = rows
    return tables, report


def _split_ids(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def _unsafe_recorded_path(value: str) -> bool:
    normalized = value.replace("\\", "/")
    posix = PurePosixPath(normalized)
    windows = PureWindowsPath(value)
    return (
        not value
        or posix.is_absolute()
        or windows.is_absolute()
        or bool(windows.drive)
        or normalized.startswith("~")
        or ".." in posix.parts
    )


def _recorded_path_exists(project_dir: Path, value: str) -> bool:
    if _unsafe_recorded_path(value):
        return False
    candidates = [project_dir / value]
    for ancestor in (project_dir, *project_dir.parents):
        if (ancestor / "pyproject.toml").exists():
            candidates.append(ancestor / value)
            break
    for path in candidates:
        try:
            if path.is_file() and path.stat().st_size > 0:
                return True
        except OSError:
            continue
    return False


def _parse_iso_timestamp(value: str) -> datetime | None:
    if not isinstance(value, str):
        return None
    normalized = value.strip()
    if normalized.endswith("Z"):
        normalized = normalized[:-1] + "+00:00"
    # Python 3.10 rejects fractional seconds longer than six digits, while
    # newer runtimes accept and truncate them. Preserve the instant at the
    # precision datetime can represent so recorded runs validate identically.
    match = re.fullmatch(
        r"(?P<prefix>.+[T ]\d{2}:\d{2}:\d{2})\."
        r"(?P<fraction>\d{7,})(?P<offset>[+-]\d{2}:\d{2})?",
        normalized,
    )
    if match:
        normalized = (
            f"{match.group('prefix')}.{match.group('fraction')[:6]}"
            f"{match.group('offset') or ''}"
        )
    try:
        return datetime.fromisoformat(normalized)
    except (TypeError, ValueError):
        return None


def _passed_run_evidence_gaps(
    project_dir: Path, row: dict[str, str]
) -> list[str]:
    """Return evidence defects that prevent a passed run from qualifying."""

    gaps: list[str] = []
    run_kind = row.get("run_kind", "")
    if run_kind not in {
        "real_data", "unit", "synthetic_integration", "environment_smoke"
    }:
        gaps.append("run_kind is missing or invalid")
    for field_name in ("command", "commit", "environment"):
        if not row.get(field_name, "").strip():
            gaps.append(f"{field_name} is missing")
    if not row.get("input_provenance", "").strip():
        gaps.append("input_provenance is missing")
    started_at = _parse_iso_timestamp(row.get("started_at", ""))
    finished_at = _parse_iso_timestamp(row.get("finished_at", ""))
    if started_at is None or finished_at is None:
        gaps.append("valid started_at and finished_at timestamps are required")
    elif started_at.utcoffset() is None or finished_at.utcoffset() is None:
        gaps.append("started_at and finished_at must include timezone offsets")
    elif finished_at < started_at:
        gaps.append("finished_at precedes started_at")
    if row.get("exit_code", "") != "0":
        gaps.append("a passed run must record exit_code=0")
    if run_kind == "real_data" and not _recorded_path_exists(
        project_dir, row.get("data_manifest", "")
    ):
        gaps.append("real-data data_manifest is missing or unsafe")
    if not _recorded_path_exists(project_dir, row.get("log", "")):
        gaps.append("log is missing or unsafe")
    artifacts = _split_ids(row.get("artifacts", ""))
    if not artifacts:
        gaps.append("artifacts are missing")
    else:
        missing = [
            value for value in artifacts
            if not _recorded_path_exists(project_dir, value)
        ]
        if missing:
            gaps.append("declared artifacts are missing or unsafe: " + ", ".join(missing))
    return gaps


def _qualifying_real_data_runs(
    project_dir: Path,
    route_id: str,
    tables: dict[str, list[dict[str, str]]],
) -> list[dict[str, str]]:
    return [
        row
        for row in tables.get("runs", [])
        if row.get("route_id") == route_id
        and row.get("run_kind") == "real_data"
        and row.get("status") == "passed"
        and not _passed_run_evidence_gaps(project_dir, row)
    ]


def _validate_exact_keys(
    value: dict[str, Any], expected: set[str], label: str, report: ValidationReport
) -> None:
    missing = sorted(expected - set(value))
    unexpected = sorted(set(value) - expected)
    if missing:
        report.errors.append(f"{label}: missing keys {', '.join(missing)}")
    if unexpected:
        report.errors.append(
            f"{label}: unexpected keys {', '.join(unexpected)}"
        )


def _validate_manifest(
    manifest: dict[str, Any], report: ValidationReport
) -> None:
    _validate_exact_keys(manifest, PROJECT_KEYS, "PROJECT.json", report)
    if manifest.get("schema_version") != SCHEMA_VERSION:
        report.errors.append(
            f"PROJECT.json: schema_version must be {SCHEMA_VERSION}"
        )
    for key in {"project_id", "title", "stop_rule"}:
        value = manifest.get(key)
        if not isinstance(value, str) or not value.strip():
            report.errors.append(f"PROJECT.json: {key} must be non-empty text")
    workspace_kind = manifest.get("workspace_kind")
    if workspace_kind not in {"pilot", "manuscript_project"}:
        report.errors.append(
            f"PROJECT.json: invalid workspace_kind={workspace_kind!r}"
        )
    allowed_stages = (
        PILOT_STAGES if workspace_kind == "pilot" else MANUSCRIPT_STAGES
    )
    if manifest.get("stage") not in allowed_stages:
        report.errors.append(f"PROJECT.json: invalid stage={manifest.get('stage')!r}")
    selected = manifest.get("selected_route_id")
    if not isinstance(selected, str):
        report.errors.append("PROJECT.json: selected_route_id must be text")

    for object_name, expected_keys in PROJECT_OBJECT_KEYS.items():
        value = manifest.get(object_name)
        if not isinstance(value, dict):
            report.errors.append(f"PROJECT.json: {object_name} must be an object")
            continue
        _validate_exact_keys(
            value, expected_keys, f"PROJECT.json {object_name}", report
        )

    anchor = manifest.get("anchor", {})
    if isinstance(anchor, dict):
        if not isinstance(anchor.get("title"), str) or not anchor.get(
            "title", ""
        ).strip():
            report.errors.append("PROJECT.json anchor: title must be non-empty")
        for key in {"doi", "citation", "source_location"}:
            if not isinstance(anchor.get(key), str):
                report.errors.append(
                    f"PROJECT.json anchor: {key} must be text"
                )

    goal = manifest.get("goal", {})
    if isinstance(goal, dict):
        for key in PROJECT_OBJECT_KEYS["goal"]:
            if not isinstance(goal.get(key), str) or not goal.get(key, "").strip():
                report.errors.append(
                    f"PROJECT.json goal: {key} must be non-empty text"
                )

    constraints = manifest.get("constraints", {})
    if isinstance(constraints, dict):
        if constraints.get("beginner_led") is not True:
            report.errors.append(
                "PROJECT.json constraints: beginner_led must be true"
            )
        for key in {"deadline", "compute", "skills"}:
            if not isinstance(constraints.get(key), str):
                report.errors.append(
                    f"PROJECT.json constraints: {key} must be text"
                )

    provenance = manifest.get("provenance", {})
    if isinstance(provenance, dict):
        for key in PROJECT_OBJECT_KEYS["provenance"]:
            if not isinstance(provenance.get(key), str):
                report.errors.append(
                    f"PROJECT.json provenance: {key} must be text"
                )
        source_values = [
            provenance.get("source_pilot_id", ""),
            provenance.get("source_route_id", ""),
            provenance.get("source_run_ids", ""),
        ]
        if workspace_kind == "pilot" and any(source_values):
            report.errors.append(
                "PROJECT.json pilot provenance must not claim a source Pilot, "
                "route or run"
            )
        if workspace_kind == "manuscript_project" and not all(source_values):
            report.errors.append(
                "PROJECT.json manuscript_project provenance requires source_pilot_id "
                "source_route_id and source_run_ids"
            )


def _index_rows(
    tables: dict[str, list[dict[str, str]]], report: ValidationReport
) -> dict[str, dict[str, dict[str, str]]]:
    indexes: dict[str, dict[str, dict[str, str]]] = {}
    for name, spec in TABLES.items():
        id_field = spec["id_field"]
        index: dict[str, dict[str, str]] = {}
        for row_number, row in enumerate(tables.get(name, []), start=2):
            entity_id = row.get(id_field, "")
            if not entity_id:
                report.errors.append(
                    f"{spec['path']}:{row_number}: empty {id_field}"
                )
                continue
            if entity_id in index:
                report.errors.append(
                    f"{spec['path']}:{row_number}: duplicate {entity_id}"
                )
            index[entity_id] = row
        indexes[name] = index
    return indexes


def _validate_references(
    tables: dict[str, list[dict[str, str]]],
    indexes: dict[str, dict[str, dict[str, str]]],
    report: ValidationReport,
) -> None:
    for name, spec in TABLES.items():
        for row_number, row in enumerate(tables.get(name, []), start=2):
            for reference in spec["references"]:
                field_name = reference["field"]
                values = (
                    _split_ids(row.get(field_name, ""))
                    if reference.get("multi")
                    else [row.get(field_name, "")]
                )
                if not values or not values[0]:
                    if reference.get("optional"):
                        continue
                    report.errors.append(
                        f"{spec['path']}:{row_number}: empty reference "
                        f"{field_name}"
                    )
                    continue
                target = reference["target"]
                for value in values:
                    if value not in indexes.get(target, {}):
                        report.errors.append(
                            f"{spec['path']}:{row_number}: unresolved "
                            f"{field_name}={value}"
                        )


def _validate_route_alignment(
    tables: dict[str, list[dict[str, str]]],
    indexes: dict[str, dict[str, dict[str, str]]],
    report: ValidationReport,
) -> None:
    data_requirements = indexes.get("data_requirements", {})
    data_candidates = indexes.get("data_candidates", {})
    code_requirements = indexes.get("code_requirements", {})
    figures = indexes.get("figures", {})
    runs = indexes.get("runs", {})

    for row in tables.get("data_candidates", []):
        parent = data_requirements.get(row.get("requirement_id", ""), {})
        if parent and parent.get("route_id") != row.get("route_id"):
            report.errors.append(
                f"data candidate {row.get('data_id')} route does not match "
                "its requirement"
            )
    for row in tables.get("cohort_usage", []):
        parent = data_candidates.get(row.get("data_id", ""), {})
        if parent and parent.get("route_id") != row.get("route_id"):
            report.errors.append(
                f"cohort usage {row.get('usage_id')} route does not match "
                "its data candidate"
            )
    for row in tables.get("code_candidates", []):
        parent = code_requirements.get(row.get("module_id", ""), {})
        if parent and parent.get("route_id") != row.get("route_id"):
            report.errors.append(
                f"code candidate {row.get('code_id')} route does not match "
                "its module"
            )
    for row in tables.get("figures", []):
        route_id = row.get("route_id")
        for requirement_id in _split_ids(row.get("data_requirement_ids", "")):
            parent = data_requirements.get(requirement_id, {})
            if parent and parent.get("route_id") != route_id:
                report.errors.append(
                    f"figure {row.get('figure_id')} uses a data requirement "
                    "from another route"
                )
        for module_id in _split_ids(row.get("code_module_ids", "")):
            parent = code_requirements.get(module_id, {})
            if parent and parent.get("route_id") != route_id:
                report.errors.append(
                    f"figure {row.get('figure_id')} uses a code module from "
                    "another route"
                )
    for row in tables.get("results", []):
        route_id = row.get("route_id")
        run = runs.get(row.get("run_id", ""), {})
        figure = figures.get(row.get("figure_id", ""), {})
        if run and run.get("route_id") != route_id:
            report.errors.append(
                f"result {row.get('result_id')} run belongs to another route"
            )
        if figure and figure.get("route_id") != route_id:
            report.errors.append(
                f"result {row.get('result_id')} figure belongs to another route"
            )


def _positive_integer(value: str) -> bool:
    try:
        return int(value) > 0
    except (TypeError, ValueError):
        return False


def _nonnegative_integer(value: str) -> bool:
    try:
        return int(value) >= 0
    except (TypeError, ValueError):
        return False


def _normalized_items(value: str) -> set[str]:
    return {
        item.strip().casefold()
        for item in value.split(";")
        if item.strip()
    }


def _same_contract_value(left: str, right: str) -> bool:
    return " ".join(left.split()).casefold() == " ".join(right.split()).casefold()


def _qualified_license(value: str) -> bool:
    normalized = " ".join(value.split()).casefold()
    if not normalized:
        return False
    disallowed = (
        "not selected", "unknown", "unlicensed", "no license",
        "all rights reserved", "internal project use only",
    )
    return not any(marker in normalized for marker in disallowed)


def _data_candidate_meets_contract(
    requirement: dict[str, str], candidate: dict[str, str]
) -> bool:
    if candidate.get("contract_match") != "true":
        return False
    for field_name in ("disease", "tissue", "modality", "biological_unit"):
        if not _same_contract_value(
            requirement.get(field_name, ""), candidate.get(field_name, "")
        ):
            return False
    minimum = requirement.get("minimum_subjects", "")
    observed = candidate.get("subjects", "")
    if not _positive_integer(minimum) or not _positive_integer(observed):
        return False
    if int(observed) < int(minimum):
        return False
    required_fields = _normalized_items(requirement.get("required_fields", ""))
    checked_fields = _normalized_items(candidate.get("fields_checked", ""))
    if not required_fields.issubset(checked_fields):
        return False
    return True


def _code_candidate_meets_contract(
    module: dict[str, str], candidate: dict[str, str]
) -> bool:
    if candidate.get("contract_match") != "true":
        return False
    if any(
        not candidate.get(field_name)
        for field_name in (
            "version", "environment", "entrypoint",
            "smoke_input", "smoke_output", "tests_passed",
        )
    ) or not _qualified_license(candidate.get("license", "")):
        return False
    required_tests = _normalized_items(module.get("required_tests", ""))
    passed_tests = _normalized_items(candidate.get("tests_passed", ""))
    return required_tests.issubset(passed_tests)


def _data_resources_meet_contract(
    requirement: dict[str, str],
    candidate: dict[str, str],
    tables: dict[str, list[dict[str, str]]],
) -> bool:
    resources = [
        row for row in tables.get("data_resources", [])
        if row.get("data_id") == candidate.get("data_id")
        and row.get("verification") in DATA_READY
    ]
    if not resources:
        return False
    supplied: set[str] = set()
    for resource in resources:
        supplied.update(_normalized_items(resource.get("fields_supplied", "")))
    required = _normalized_items(requirement.get("required_fields", ""))
    return required.issubset(supplied)


def _cohort_usage_gaps(
    route_id: str, tables: dict[str, list[dict[str, str]]]
) -> list[str]:
    gaps: list[str] = []
    uses = [
        row for row in tables.get("cohort_usage", [])
        if row.get("route_id") == route_id
    ]
    if any(row.get("acceptable") == "false" for row in uses):
        gaps.append("cohort usage ledger contains an unacceptable analysis role")
    by_cohort: dict[str, list[dict[str, str]]] = {}
    for row in uses:
        by_cohort.setdefault(row.get("cohort_key", ""), []).append(row)
    for cohort_key, rows in by_cohort.items():
        influenced_development = any(
            row.get(field) == "true"
            for row in rows
            for field in (
                "features_influenced", "parameters_influenced",
                "cutoff_influenced",
            )
        )
        claimed_external = any(
            row.get("role") == "external_validation"
            or row.get("claimed_external_validation") == "true"
            for row in rows
        )
        if influenced_development and claimed_external:
            gaps.append(
                f"cohort {cohort_key} influenced development and cannot be "
                "claimed as external validation"
            )
    return gaps


def _execution_gaps_from_tables(
    project_dir: Path,
    route: dict[str, str],
    tables: dict[str, list[dict[str, str]]],
) -> list[str]:
    route_id = route.get("route_id", "")
    gaps: list[str] = []

    if route.get("decision_status") in {"rejected", "stopped"}:
        return ["route is rejected or stopped"]
    if route.get("scientific_review") != "passed":
        gaps.append("scientific design has not passed review")
    if not _positive_integer(route.get("minimum_main_figures", "")):
        gaps.append("minimum_main_figures is not a positive integer")
    unresolved_burdens = [
        field_name
        for field_name in ("data_burden", "code_burden", "beginner_burden")
        if route.get(field_name) == "unknown"
    ]
    if unresolved_burdens:
        gaps.append(
            "unresolved execution burden: " + ", ".join(unresolved_burdens)
        )

    logs = [
        row for row in tables.get("search_log", [])
        if row.get("route_id") == route_id and row.get("outcome") != "stop"
    ]
    logged_domains = {row.get("domain") for row in logs}
    for domain in ("data", "code", "literature"):
        if domain not in logged_domains:
            gaps.append(f"no recorded {domain} search")

    requirements = [
        row for row in tables.get("data_requirements", [])
        if row.get("route_id") == route_id and row.get("required") == "true"
    ]
    if not requirements:
        gaps.append("no required data contract")
    candidates = tables.get("data_candidates", [])
    for requirement in requirements:
        requirement_id = requirement.get("requirement_id", "")
        usable = [
            row for row in candidates
            if row.get("requirement_id") == requirement_id
            and row.get("route_id") == route_id
            and row.get("decision") == "use"
            and row.get("verification") in DATA_READY
            and row.get("access") not in {"unavailable", "unknown"}
            and _data_candidate_meets_contract(requirement, row)
            and _data_resources_meet_contract(requirement, row, tables)
        ]
        if requirement.get("independence_required") == "true":
            usable = [
                row for row in usable if row.get("independence") == "independent"
            ]
        if not usable:
            gaps.append(
                f"data requirement {requirement_id} lacks a parsed usable "
                "candidate meeting its contract, field-level resource "
                "provenance and minimum subject count"
            )

    used_data_ids = {
        row.get("data_id") for row in candidates
        if row.get("route_id") == route_id and row.get("decision") == "use"
    }
    usage_data_ids = {
        row.get("data_id") for row in tables.get("cohort_usage", [])
        if row.get("route_id") == route_id
    }
    for data_id in sorted(used_data_ids - usage_data_ids):
        gaps.append(f"data candidate {data_id} lacks a cohort usage ledger entry")
    gaps.extend(_cohort_usage_gaps(route_id, tables))

    modules = [
        row for row in tables.get("code_requirements", [])
        if row.get("route_id") == route_id and row.get("required") == "true"
    ]
    if not modules:
        gaps.append("no required code module contract")
    code_candidates = tables.get("code_candidates", [])
    for module in modules:
        module_id = module.get("module_id", "")
        usable = []
        for row in code_candidates:
            if (
                row.get("module_id") != module_id
                or row.get("route_id") != route_id
                or row.get("decision") != "use"
                or row.get("verification") not in CODE_READY
            ):
                continue
            if any(
                not row.get(field)
                for field in ("version", "license", "environment", "entrypoint")
            ):
                continue
            if not _code_candidate_meets_contract(module, row):
                continue
            if row.get("noninteractive") != "true":
                continue
            if row.get("private_inputs") != "false":
                continue
            if row.get("hardcoded_paths") != "false":
                continue
            if row.get("path_portability") not in PATH_READY:
                continue
            usable.append(row)
        if not usable:
            gaps.append(
                f"code module {module_id} lacks a qualified smoke-tested donor"
            )

    if route.get("model_spec_required") == "true":
        model_specs = [
            row for row in tables.get("model_specifications", [])
            if row.get("route_id") == route_id
            and row.get("status") in MODEL_READY
        ]
        if not model_specs:
            gaps.append("route requires a locked computable model specification")

    required_figures = [
        row for row in tables.get("figures", [])
        if row.get("route_id") == route_id and row.get("required") == "true"
    ]
    required_main = [row for row in required_figures if row.get("role") == "main"]
    minimum = int(route.get("minimum_main_figures", "0")) if _positive_integer(
        route.get("minimum_main_figures", "")
    ) else 0
    if len(required_main) < minimum:
        gaps.append(
            f"only {len(required_main)} required main figures; minimum is {minimum}"
        )
    for figure in required_figures:
        if FIGURE_RANK.get(figure.get("status", ""), -1) < 1:
            gaps.append(f"figure {figure.get('figure_id')} is not mapped")
        if FIGURE_RANK.get(figure.get("status", ""), -1) >= 2:
            source_tables = _split_ids(figure.get("source_table", ""))
            missing_sources = [
                value for value in source_tables
                if not _recorded_path_exists(project_dir, value)
            ]
            if not source_tables or missing_sources:
                gaps.append(
                    f"figure {figure.get('figure_id')} lacks existing source-table "
                    "artifacts"
                )
    if required_figures and not any(
        FIGURE_RANK.get(row.get("status", ""), -1) >= 2
        for row in required_figures
    ):
        gaps.append("no representative figure or source table was generated")

    passed_runs = _qualifying_real_data_runs(project_dir, route_id, tables)
    if not passed_runs:
        gaps.append("no passed real-data run with complete evidence is recorded")
    else:
        passed_run_ids = {row.get("run_id", "") for row in passed_runs}
        run_results = [
            row for row in tables.get("results", [])
            if row.get("route_id") == route_id
            and row.get("run_id") in passed_run_ids
            and row.get("status") in {"provisional", "verified"}
            and all(
                _recorded_path_exists(project_dir, value)
                for value in _split_ids(row.get("source_table", ""))
            )
            and bool(_split_ids(row.get("source_table", "")))
        ]
        if not run_results:
            gaps.append(
                "no result with an existing source table is linked to a passed "
                "real-data run"
            )

    literature = [
        row for row in tables.get("literature", [])
        if row.get("route_id") == route_id
    ]
    if not literature:
        gaps.append("no nearest-publication comparison")
    elif any(
        row.get("overlap_level") == "duplicate"
        or row.get("decision") == "stop"
        for row in literature
    ):
        gaps.append("publication comparison contains a blocking duplicate")
    elif not any(
        row.get("decision") in {"continue", "distinguish"}
        for row in literature
    ):
        gaps.append("publication comparison has no proceed decision")

    blocking_issues = [
        row for row in tables.get("issues", [])
        if row.get("route_id") in {"", route_id}
        and row.get("blocking") == "true"
        and row.get("status") in {"open", "accepted_risk"}
    ]
    if blocking_issues:
        gaps.append(
            "open blocking issues: "
            + ", ".join(row.get("issue_id", "") for row in blocking_issues)
        )

    return gaps


def _promotion_evidence_gaps(
    project_dir: Path,
    route: dict[str, str],
    tables: dict[str, list[dict[str, str]]],
) -> list[str]:
    gaps = _execution_gaps_from_tables(project_dir, route, tables)
    if route.get("route_role") != "manuscript_candidate":
        gaps.insert(
            0,
            f"route role {route.get('route_role') or 'missing'} cannot be promoted "
            "to a manuscript project",
        )
    approvals = [
        row for row in tables.get("decisions", [])
        if row.get("route_id") == route.get("route_id")
        and row.get("decision") == "approve_route_evidence"
        and row.get("authority") == "user"
        and row.get("stage") == "pilot_review"
    ]
    if route.get("route_role") == "manuscript_candidate":
        if not approvals:
            gaps.append(
                "route evidence lacks explicit user approval at pilot_review"
            )
        else:
            approval_times = [
                timestamp
                for row in approvals
                if (
                    (timestamp := _parse_iso_timestamp(
                        row.get("decided_at", "")
                    )) is not None
                    and timestamp.utcoffset() is not None
                )
            ]
            if not approval_times:
                gaps.append(
                    "approve_route_evidence must record a full timestamp with a "
                    "timezone offset"
                )
            qualifying_runs = _qualifying_real_data_runs(
                project_dir, route.get("route_id", ""), tables
            )
            finished = [
                timestamp
                for row in qualifying_runs
                if (
                    timestamp := _parse_iso_timestamp(
                        row.get("finished_at", "")
                    )
                ) is not None
            ]
            if finished and approval_times and not any(
                approval >= max(finished) for approval in approval_times
            ):
                gaps.append(
                    "approve_route_evidence must occur on or after the latest "
                    "qualified real-data run"
                )
    return gaps


def _final_promotion_decision_gaps(
    project_dir: Path,
    route: dict[str, str],
    tables: dict[str, list[dict[str, str]]],
) -> list[str]:
    """Require a fresh final user decision after evidence approval and execution."""

    route_id = route.get("route_id", "")
    promote_rows = [
        row
        for row in tables.get("decisions", [])
        if row.get("route_id") == route_id
        and row.get("decision") == "promote"
        and row.get("authority") == "user"
        and row.get("stage") == "pilot_review"
    ]
    if not promote_rows:
        return ["route lacks an explicit user promote decision at pilot_review"]
    promote_times = [
        timestamp
        for row in promote_rows
        if (
            timestamp := _parse_iso_timestamp(row.get("decided_at", ""))
        ) is not None
        and timestamp.utcoffset() is not None
    ]
    if not promote_times:
        return [
            "promote decision must record a full timestamp with a timezone offset"
        ]

    approval_times = [
        timestamp
        for row in tables.get("decisions", [])
        if row.get("route_id") == route_id
        and row.get("decision") == "approve_route_evidence"
        and row.get("authority") == "user"
        and row.get("stage") == "pilot_review"
        and (
            timestamp := _parse_iso_timestamp(row.get("decided_at", ""))
        ) is not None
        and timestamp.utcoffset() is not None
    ]
    qualifying_runs = _qualifying_real_data_runs(project_dir, route_id, tables)
    finished = [
        timestamp
        for row in qualifying_runs
        if (
            timestamp := _parse_iso_timestamp(row.get("finished_at", ""))
        ) is not None
    ]
    if not approval_times or not finished:
        return [
            "promote decision requires a prior timely route-evidence approval "
            "and qualified real-data run"
        ]
    latest_run = max(finished)
    timely_approvals = [value for value in approval_times if value >= latest_run]
    if not timely_approvals or not any(
        promote >= approval and promote >= latest_run
        for promote in promote_times
        for approval in timely_approvals
    ):
        return [
            "promote decision must occur on or after both the latest qualified "
            "real-data run and a timely approve_route_evidence decision"
        ]
    return []


def _figure_loop_gaps(
    project_dir: Path,
    route: dict[str, str],
    tables: dict[str, list[dict[str, str]]],
) -> list[str]:
    gaps = _execution_gaps_from_tables(project_dir, route, tables)
    route_id = route.get("route_id", "")
    required_figures = [
        row for row in tables.get("figures", [])
        if row.get("route_id") == route_id and row.get("required") == "true"
    ]
    passed_runs = {
        row.get("run_id", "")
        for row in _qualifying_real_data_runs(project_dir, route_id, tables)
    }
    verified_figure_ids = {
        row.get("figure_id", "") for row in tables.get("results", [])
        if row.get("route_id") == route_id
        and row.get("run_id") in passed_runs
        and row.get("status") == "verified"
    }
    for figure in required_figures:
        figure_id = figure.get("figure_id", "")
        if FIGURE_RANK.get(figure.get("status", ""), -1) < 3:
            gaps.append(f"figure {figure_id} has not completed the figure loop")
        if figure_id not in verified_figure_ids:
            gaps.append(f"figure {figure_id} lacks a verified result from a passed run")
        source_tables = _split_ids(figure.get("source_table", ""))
        if not source_tables or any(
            not _recorded_path_exists(project_dir, value) for value in source_tables
        ):
            gaps.append(f"figure {figure_id} lacks existing source-table artifacts")
    return list(dict.fromkeys(gaps))


def route_readiness(project_dir: Path) -> list[dict[str, Any]]:
    project_dir = Path(project_dir)
    tables, _ = load_tables(project_dir)
    workspace_errors = validate_workspace(project_dir).errors
    readiness: list[dict[str, Any]] = []
    for route in tables.get("routes", []):
        execution_gaps = _execution_gaps_from_tables(project_dir, route, tables)
        execution_gaps.extend(
            f"workspace validation error: {error}" for error in workspace_errors
        )
        execution_gaps = list(dict.fromkeys(execution_gaps))
        promotion_gaps = _promotion_evidence_gaps(project_dir, route, tables)
        if workspace_errors:
            promotion_gaps.extend(
                f"workspace validation error: {error}" for error in workspace_errors
            )
            promotion_gaps = list(dict.fromkeys(promotion_gaps))
        readiness.append({
            "route_id": route.get("route_id"),
            "title": route.get("title"),
            "decision_status": route.get("decision_status"),
            "evidence_stage": route.get("evidence_stage"),
            "route_role": route.get("route_role"),
            "execution_ready": not execution_gaps,
            "promotion_evidence_complete": not promotion_gaps,
            "execution_gaps": execution_gaps,
            "promotion_gaps": promotion_gaps,
        })
    return readiness


def _document_complete(path: Path, minimum_length: int) -> bool:
    if not path.exists():
        return False
    text = path.read_text(encoding="utf-8")
    return TODO not in text and len(text.strip()) >= minimum_length


def validate_workspace(project_dir: Path) -> ValidationReport:
    project_dir = Path(project_dir)
    report = ValidationReport()
    manifest_path = project_dir / "PROJECT.json"
    if not manifest_path.exists():
        report.errors.append("missing PROJECT.json")
        return report
    try:
        manifest = load_manifest(project_dir)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        report.errors.append(f"invalid PROJECT.json: {exc}")
        return report

    _validate_manifest(manifest, report)
    tables, table_report = load_tables(project_dir)
    report.errors.extend(table_report.errors)
    report.warnings.extend(table_report.warnings)
    indexes = _index_rows(tables, report)
    _validate_references(tables, indexes, report)
    _validate_route_alignment(tables, indexes, report)

    for row in tables.get("search_log", []):
        if not _nonnegative_integer(row.get("result_count", "")):
            report.errors.append(
                f"search {row.get('search_id')} result_count must be a "
                "non-negative integer"
            )
    for row in tables.get("data_requirements", []):
        if not _positive_integer(row.get("minimum_subjects", "")):
            report.errors.append(
                f"data requirement {row.get('requirement_id')} "
                "minimum_subjects must be positive"
            )
    for row in tables.get("data_candidates", []):
        if row.get("decision") != "use":
            continue
        requirement = indexes.get("data_requirements", {}).get(
            row.get("requirement_id", ""), {}
        )
        if (
            row.get("verification") not in DATA_READY
            or row.get("access") in {"unavailable", "unknown"}
            or not requirement
            or not _data_candidate_meets_contract(requirement, row)
            or not _data_resources_meet_contract(requirement, row, tables)
            or (
                requirement.get("independence_required") == "true"
                and row.get("independence") != "independent"
            )
        ):
            report.errors.append(
                f"data candidate {row.get('data_id')} is marked use but does "
                "not satisfy access, parsing, contract, field-level resource "
                "provenance, sample-size and independence requirements"
            )
    for row in tables.get("data_resources", []):
        verification = row.get("verification")
        if verification in {"downloaded", "checksum_verified"} and not _positive_integer(
            row.get("bytes", "")
        ):
            report.errors.append(
                f"data resource {row.get('resource_id')} must record positive bytes"
            )
        if verification == "checksum_verified" and not row.get("checksum"):
            report.errors.append(
                f"data resource {row.get('resource_id')} must record a checksum"
            )
    for row in tables.get("cohort_usage", []):
        if (
            row.get("claimed_external_validation") == "true"
            and row.get("role") != "external_validation"
        ):
            report.errors.append(
                f"cohort usage {row.get('usage_id')} claims external validation "
                "but its role is not external_validation"
            )
    for row in tables.get("code_candidates", []):
        if row.get("decision") != "use":
            continue
        module = indexes.get("code_requirements", {}).get(
            row.get("module_id", ""), {}
        )
        if (
            row.get("verification") not in CODE_READY
            or not module
            or not _code_candidate_meets_contract(module, row)
            or row.get("noninteractive") != "true"
            or row.get("private_inputs") != "false"
            or row.get("hardcoded_paths") != "false"
            or row.get("path_portability") not in PATH_READY
        ):
            report.errors.append(
                f"code candidate {row.get('code_id')} is marked use but does "
                "not satisfy contract, provenance, non-interactive execution "
                "and required-test evidence"
            )

    for row in tables.get("runs", []):
        run_id = row.get("run_id", "")
        status = row.get("status", "")
        started_at = _parse_iso_timestamp(row.get("started_at", ""))
        finished_value = row.get("finished_at", "")
        finished_at = _parse_iso_timestamp(finished_value) if finished_value else None
        if started_at is None:
            report.errors.append(f"run {run_id} started_at must be an ISO timestamp")
        if status in {"passed", "failed", "invalidated"}:
            if finished_at is None:
                report.errors.append(
                    f"terminal run {run_id} requires a valid finished_at timestamp"
                )
            elif started_at is not None and finished_at < started_at:
                report.errors.append(
                    f"run {run_id} finished_at precedes started_at"
                )
        if status != "passed":
            continue
        for gap in _passed_run_evidence_gaps(project_dir, row):
            report.errors.append(f"passed run {run_id}: {gap}")

    for row in tables.get("results", []):
        if row.get("status") == "invalidated":
            continue
        result_id = row.get("result_id", "")
        run = indexes.get("runs", {}).get(row.get("run_id", ""), {})
        if run.get("run_kind") != "real_data":
            report.errors.append(
                f"result {result_id} must be linked to a real_data run"
            )
        source_tables = _split_ids(row.get("source_table", ""))
        if not source_tables or any(
            not _recorded_path_exists(project_dir, value) for value in source_tables
        ):
            report.errors.append(
                f"result {result_id} source_table must reference existing, safe "
                "repository-relative artifacts"
            )

    for row in tables.get("issues", []):
        issue_id = row.get("issue_id", "")
        disposition = row.get("disposition", "")
        candidate_scope = row.get("candidate_scope", "")
        promotion_id = row.get("promotion_id", "")
        affected_capabilities = _split_ids(
            row.get("affected_capability_ids", "")
        )
        if disposition == "promote_to_module":
            if candidate_scope != "module":
                report.errors.append(
                    f"issue {issue_id} promotes to module but candidate_scope "
                    "is not module"
                )
            if not promotion_id or not affected_capabilities:
                report.errors.append(
                    f"issue {issue_id} module promotion requires promotion_id "
                    "and affected_capability_ids"
                )
        elif disposition == "promote_to_core":
            if candidate_scope != "core":
                report.errors.append(
                    f"issue {issue_id} promotes to core but candidate_scope is not core"
                )
            if not promotion_id or not affected_capabilities:
                report.errors.append(
                    f"issue {issue_id} core promotion requires promotion_id and "
                    "affected_capability_ids"
                )
        elif promotion_id:
            report.errors.append(
                f"issue {issue_id} records promotion_id without a promotion disposition"
            )

    for route in tables.get("routes", []):
        if not _positive_integer(route.get("minimum_main_figures", "")):
            report.errors.append(
                f"route {route.get('route_id')} minimum_main_figures must be positive"
            )
        execution_gaps = _execution_gaps_from_tables(project_dir, route, tables)
        evidence_stage = route.get("evidence_stage", "")
        if (
            manifest.get("workspace_kind") == "pilot"
            and evidence_stage == "direction_audited"
            and route.get("decision_status") in {"active", "backup"}
        ):
            report.errors.append(
                f"route {route.get('route_id')} must remain candidate at "
                "direction_audited; active/backup requires the route-specific "
                "availability precheck used for an execution choice"
            )
        if (
            EVIDENCE_STAGE_RANK.get(evidence_stage, -1)
            >= EVIDENCE_STAGE_RANK["minimal_real_run"]
            and execution_gaps
        ):
            report.errors.append(
                f"route {route.get('route_id')} claims {evidence_stage} but has "
                f"execution gaps: {'; '.join(execution_gaps)}"
            )
        if evidence_stage == "figure_loop_closed":
            figure_gaps = _figure_loop_gaps(project_dir, route, tables)
            if figure_gaps:
                report.errors.append(
                    f"route {route.get('route_id')} claims figure_loop_closed but "
                    f"has gaps: {'; '.join(figure_gaps)}"
                )
        if route.get("decision_status") == "promoted":
            promotion_gaps = _promotion_evidence_gaps(project_dir, route, tables)
            promotion_gaps.extend(
                _final_promotion_decision_gaps(project_dir, route, tables)
            )
            if promotion_gaps:
                report.errors.append(
                    f"route {route.get('route_id')} is marked promoted but its "
                    f"promotion evidence is incomplete: {'; '.join(promotion_gaps)}"
                )
        if not execution_gaps and EVIDENCE_STAGE_RANK.get(evidence_stage, -1) < 2:
            report.warnings.append(
                f"route {route.get('route_id')} has no execution gaps and can "
                "record evidence_stage=minimal_real_run"
            )

    selected_route_id = manifest.get("selected_route_id", "")
    workspace_kind = manifest.get("workspace_kind")
    promoted_rows = [
        row for row in tables.get("routes", [])
        if row.get("decision_status") == "promoted"
    ]
    if workspace_kind == "pilot" and any(
        row.get("decision") == "approve_release"
        for row in tables.get("decisions", [])
    ):
        report.errors.append(
            "pilot cannot record approve_release; release approval belongs only "
            "to a completed manuscript project"
        )
    if workspace_kind == "pilot" and len(promoted_rows) > 1:
        report.errors.append("a pilot cannot promote more than one manuscript route")
    if selected_route_id:
        selected = indexes.get("routes", {}).get(selected_route_id)
        if not selected:
            report.errors.append(
                f"PROJECT.json: unresolved selected_route_id={selected_route_id}"
            )
        elif (
            workspace_kind == "pilot"
            and selected.get("decision_status") != "promoted"
        ):
            report.errors.append(
                "PROJECT.json pilot selected_route_id must point to a promoted route"
            )
        elif (
            workspace_kind == "manuscript_project"
            and selected.get("decision_status") not in {"active", "promoted"}
        ):
            report.errors.append(
                "PROJECT.json manuscript selected_route_id must point to an active route"
            )
    if promoted_rows and promoted_rows[0].get("route_id") != selected_route_id:
        report.errors.append(
            "PROJECT.json selected_route_id does not match the promoted route"
        )

    stage = manifest.get("stage")
    if workspace_kind == "pilot" and stage in PILOT_STAGES and stage != "stopped":
        stage_index = PILOT_STAGES.index(stage)
        if stage_index >= PILOT_STAGES.index("route_generation"):
            if not _document_complete(project_dir / "anchor/audit.md", 300):
                report.errors.append(
                    "anchor/audit.md must be completed before route generation"
                )
            if not tables.get("routes"):
                report.errors.append("route generation stage requires routes")
        if stage_index >= PILOT_STAGES.index("pilot_review"):
            if not _document_complete(project_dir / "reports/pilot-outcome.md", 500):
                report.errors.append(
                    "reports/pilot-outcome.md must separate paper-side and "
                    "product-side outcomes before pilot review"
                )
        if stage == "pilot_complete":
            closing_decisions = [
                row for row in tables.get("decisions", [])
                if row.get("decision") in {
                    "promote", "retain_training", "close_pilot"
                }
                and row.get("authority") == "user"
            ]
            if not closing_decisions:
                report.errors.append(
                    "pilot_complete requires an explicit user promote, "
                    "retain_training or close_pilot decision"
                )

    if (
        workspace_kind == "manuscript_project"
        and stage in MANUSCRIPT_STAGES
        and stage != "stopped"
    ):
        if not selected_route_id:
            report.errors.append("manuscript_project requires selected_route_id")
        stage_index = MANUSCRIPT_STAGES.index(stage)
        if stage_index >= MANUSCRIPT_STAGES.index("execution"):
            if not _document_complete(
                project_dir / "analysis/specification.md", 300
            ):
                report.errors.append(
                    "analysis/specification.md must be completed before execution"
                )
        if stage_index >= MANUSCRIPT_STAGES.index("execution"):
            selected_runs = [
                row for row in tables.get("runs", [])
                if row.get("route_id") == selected_route_id
            ]
            if not selected_runs:
                report.errors.append("execution stage requires a recorded run")
        if stage_index >= MANUSCRIPT_STAGES.index("interpretation"):
            passed_runs = {
                row.get("run_id")
                for row in _qualifying_real_data_runs(
                    project_dir, selected_route_id, tables
                )
            }
            verified_results = [
                row for row in tables.get("results", [])
                if row.get("route_id") == selected_route_id
                and row.get("run_id") in passed_runs
                and row.get("status") == "verified"
            ]
            if not verified_results:
                report.errors.append(
                    "interpretation stage requires a verified result from a "
                    "qualified passed real-data run"
                )
        if stage_index >= MANUSCRIPT_STAGES.index("writing"):
            if not _document_complete(project_dir / "manuscript/draft.md", 800):
                report.errors.append(
                    "manuscript/draft.md must contain a substantive draft"
                )
        if stage == "complete":
            selected_route = indexes.get("routes", {}).get(selected_route_id)
            if selected_route is not None:
                completion_gaps = _figure_loop_gaps(
                    project_dir, selected_route, tables
                )
                if completion_gaps:
                    report.errors.append(
                        "complete manuscript project has unresolved execution or "
                        "figure-loop gaps: " + "; ".join(completion_gaps)
                    )
            release = [
                row for row in tables.get("decisions", [])
                if row.get("route_id") == selected_route_id
                and row.get("decision") == "approve_release"
                and row.get("authority") == "user"
                and row.get("stage") == "complete"
            ]
            if not release:
                report.errors.append(
                    "complete project requires a user approve_release decision "
                    "recorded at stage=complete"
                )
            else:
                qualifying_runs = _qualifying_real_data_runs(
                    project_dir, selected_route_id, tables
                )
                finished = [
                    timestamp
                    for row in qualifying_runs
                    if (timestamp := _parse_iso_timestamp(
                        row.get("finished_at", "")
                    )) is not None
                ]
                if finished:
                    latest_run = max(finished)
                    release_times = [
                        decided
                        for row in release
                        if (
                            (decided := _parse_iso_timestamp(
                                row.get("decided_at", "")
                            )) is not None
                            and decided.utcoffset() is not None
                        )
                    ]
                    if not release_times:
                        report.errors.append(
                            "approve_release must record a full timestamp with a "
                            "timezone offset"
                        )
                    elif not any(decided >= latest_run for decided in release_times):
                        report.errors.append(
                            "approve_release must occur on or after the latest "
                            "qualified real-data run"
                        )

    return report


ANCHOR_TEMPLATE = f"""# Anchor audit

{TODO}

## Paper identity

- Full citation:
- DOI:
- Local source-material location:
- Public data and code statements:

## Scientific grammar

- Population and disease:
- Biological unit:
- Exposure or focal object:
- Comparison:
- Primary outcome:
- Central claim:
- Claim ceiling:
- Result that would falsify the central claim:

## Anchor framework assets

Describe the paper-specific question, central relationship, evidence order,
data/cohort roles, analysis-module sequence and Figure narrative that are worth
retaining. Add any framework elements unique to this paper type; do not force
the audit into a fixed substitution checklist.

## Figure-to-evidence map

For every main and supplementary figure record its manuscript role, data,
metadata, method, code, output, statistical unit and unavailable dependencies.

## Module disposition

Classify every module as retain, repair, substitute, extend, drop or blocked.

## Defect-to-repair contracts

For every material flaw record: evidence; affected Figure/claim; unaffected
framework assets; candidate repair; how the target paper would implement it; the minimum
evidence that could verify or falsify the repair; and residual risk/claim
ceiling. A flaw does not by itself justify discarding the anchor framework.

## Reproduction boundary

List private data, wet-lab evidence, unavailable code and undocumented author
choices that prevent exact reproduction.
"""


PILOT_SPECIFICATION_TEMPLATE = f"""# Pilot minimal-run specification

{TODO}

Complete the relevant sections before a minimal real-data run. This document
does not freeze a future manuscript analysis.

## Source framework retained, repaired and transformed

State which anchor-framework assets the route preserves, which defect repairs
it implements, which paper-specific elements it changes, and how those choices
form the target paper. Predefined substitution examples are not a completeness test.

## Cohorts and exclusions

## Statistical and biological units

## Exposure, comparison, outcome and covariates

## Preprocessing and missing data

## Discovery, validation, sensitivity and exploratory analyses

## Models, parameters, software versions and random seeds

## Signature formula, if applicable

## Module input/output contracts and scientific invariants

## Figure and source-table outputs

## Falsifying and claim-limiting results
"""


MANUSCRIPT_SPECIFICATION_TEMPLATE = f"""# Manuscript analysis specification

{TODO}

Freeze this specification before the first full manuscript execution.

## Source pilot, route and inherited evidence

## Cohorts and exclusions

## Statistical and biological units

## Exposure, comparison, outcome and covariates

## Preprocessing and missing data

## Primary, validation, sensitivity and exploratory analyses

## Models, parameters, software versions and random seeds

## Signature formula, if applicable

## Module releases, input/output contracts and scientific invariants

## Figure and source-table outputs

## Falsifying and claim-limiting results
"""


PILOT_OUTCOME_TEMPLATE = f"""# Pilot outcome

{TODO}

## Paper-side outcome

Record what was learned about the anchor, candidate routes, real-data runs,
scientific limitations and the continue/refine/reroute/stop decision.

## Product-side outcome

List reusable findings separately. For each one record whether it remains
pilot-specific, is a module/core promotion candidate, or was rejected.

## Capability and regression contribution

State exactly which capability and test level this pilot exercised. Do not
describe one pilot as validation of unrelated papers or the whole workflow.

## Human decision and next boundary

Record whether a route is promoted to a manuscript project, retained as
training/supporting evidence, or the pilot is closed without promotion.
"""


MANUSCRIPT_TEMPLATE = f"""# Manuscript draft

{TODO}

## Title

## Abstract

## Introduction

## Methods

## Results

## Discussion

## Limitations

## Data and code availability
"""


PILOT_README = """# Paper2Paper pilot workspace

This directory audits one real anchor paper, tests candidate directions and
produces separate paper-side and product-side outcomes. It is not itself a
manuscript project. Large data, credentials and copyrighted PDFs stay outside
Git.
"""


MANUSCRIPT_README = """# Paper2Paper manuscript project

This directory contains one user-approved route promoted from a pilot. Its
goal is a complete figure and manuscript package, not further platform
generalization. The source pilot and route are locked in PROJECT.json.
"""


def init_workspace(
    target: Path,
    project_id: str,
    title: str,
    anchor_title: str,
    doi: str = "",
    workspace_kind: str = "pilot",
    provenance: dict[str, str] | None = None,
) -> Path:
    if workspace_kind not in {"pilot", "manuscript_project"}:
        raise ValueError(f"unsupported workspace_kind: {workspace_kind}")
    target = Path(target)
    if target.exists() and any(target.iterdir()):
        raise ValueError(f"target directory is not empty: {target}")
    folders = ["anchor", "evidence", "analysis", "execution", "reports"]
    if workspace_kind == "manuscript_project":
        folders.append("manuscript")
    for folder in folders:
        (target / folder).mkdir(parents=True, exist_ok=True)

    source = provenance or {
        "source_pilot_id": "",
        "source_route_id": "",
        "source_run_ids": "",
    }
    expected_source_keys = PROJECT_OBJECT_KEYS["provenance"]
    if set(source) != expected_source_keys:
        raise ValueError("provenance must contain the exact workspace contract keys")
    if workspace_kind == "pilot" and any(source.values()):
        raise ValueError("pilot workspaces cannot claim source-pilot provenance")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "workspace_kind": workspace_kind,
        "project_id": project_id,
        "title": title,
        "stage": "anchor_audit" if workspace_kind == "pilot" else "specification",
        "selected_route_id": "",
        "anchor": {
            "title": anchor_title,
            "doi": doi,
            "citation": "",
            "source_location": "outside Git",
        },
        "provenance": source,
        "goal": {
            "paper_type": (
                "anchor-paper pilot"
                if workspace_kind == "pilot"
                else "promoted manuscript project"
            ),
            "target_audience": "beginner-led project with AI assistance",
            "success_definition": (
                "A reviewed pilot outcome with traceable direction, data, code, "
                "minimal real execution and scoped product findings"
                if workspace_kind == "pilot"
                else "A scientifically defensible, non-duplicate manuscript with "
                "traceable data, code, figures, results and limitations"
            ),
        },
        "constraints": {
            "beginner_led": True,
            "deadline": "",
            "compute": "",
            "skills": "",
        },
        "stop_rule": (
            "Stop product work when the pilot can support a route decision and "
            "record its scoped reusable findings"
            if workspace_kind == "pilot"
            else "Freeze non-blocking workflow work when the selected route can be "
            "executed, reviewed and written reliably"
        ),
    }
    (target / "PROJECT.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    readme = PILOT_README if workspace_kind == "pilot" else MANUSCRIPT_README
    (target / "README.md").write_text(readme, encoding="utf-8")
    (target / "anchor/audit.md").write_text(ANCHOR_TEMPLATE, encoding="utf-8")
    specification = (
        PILOT_SPECIFICATION_TEMPLATE
        if workspace_kind == "pilot"
        else MANUSCRIPT_SPECIFICATION_TEMPLATE
    )
    (target / "analysis/specification.md").write_text(
        specification, encoding="utf-8"
    )
    if workspace_kind == "pilot":
        (target / "reports/pilot-outcome.md").write_text(
            PILOT_OUTCOME_TEMPLATE, encoding="utf-8"
        )
    else:
        (target / "manuscript/draft.md").write_text(
            MANUSCRIPT_TEMPLATE, encoding="utf-8"
        )
    for spec in TABLES.values():
        path = target / spec["path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\t".join(spec["columns"]) + "\n", encoding="utf-8")
    write_readiness_report(target)
    return target


def promote_workspace(
    pilot_dir: Path,
    route_id: str,
    target: Path,
    project_id: str,
    title: str,
) -> Path:
    """Create a manuscript workspace from an explicitly user-promoted pilot route."""
    pilot_dir = Path(pilot_dir)
    manifest = load_manifest(pilot_dir)
    if manifest.get("workspace_kind") != "pilot":
        raise ValueError("only a pilot workspace can promote a manuscript route")
    report = validate_workspace(pilot_dir)
    if not report.ok:
        raise ValueError(
            "source pilot must pass structural and contract validation before promotion: "
            + "; ".join(report.errors)
        )
    tables, _ = load_tables(pilot_dir)
    routes = {
        row.get("route_id", ""): row for row in tables.get("routes", [])
    }
    route = routes.get(route_id)
    if route is None:
        raise ValueError(f"unknown pilot route: {route_id}")
    if route.get("decision_status") != "promoted":
        raise ValueError(
            "pilot route must be marked decision_status=promoted before promotion"
        )
    if manifest.get("selected_route_id") != route_id:
        raise ValueError(
            "pilot selected_route_id must identify the promoted route"
        )
    gaps = _promotion_evidence_gaps(pilot_dir, route, tables)
    gaps.extend(_final_promotion_decision_gaps(pilot_dir, route, tables))
    if gaps:
        raise ValueError("route promotion evidence is incomplete: " + "; ".join(gaps))

    source_run_ids = [
        row.get("run_id", "")
        for row in _qualifying_real_data_runs(pilot_dir, route_id, tables)
    ]
    target = init_workspace(
        target=target,
        project_id=project_id,
        title=title,
        anchor_title=manifest.get("anchor", {}).get("title", ""),
        doi=manifest.get("anchor", {}).get("doi", ""),
        workspace_kind="manuscript_project",
        provenance={
            "source_pilot_id": manifest.get("project_id", ""),
            "source_route_id": route_id,
            "source_run_ids": ";".join(source_run_ids),
        },
    )
    target_manifest = load_manifest(target)
    target_manifest["selected_route_id"] = route_id
    target_manifest["anchor"] = dict(manifest.get("anchor", {}))
    (target / "PROJECT.json").write_text(
        json.dumps(target_manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    source_audit = pilot_dir / "anchor/audit.md"
    if source_audit.exists():
        (target / "anchor/audit.md").write_text(
            source_audit.read_text(encoding="utf-8"), encoding="utf-8"
        )

    route_data_ids = {
        row.get("data_id", "") for row in tables.get("data_candidates", [])
        if row.get("route_id") == route_id
    }
    for table_name, spec in TABLES.items():
        if table_name in {"runs", "results"}:
            continue
        source_rows = tables.get(table_name, [])
        if table_name == "routes":
            promoted_route = dict(route)
            promoted_route["decision_status"] = "active"
            promoted_route["evidence_stage"] = "availability_prechecked"
            copied_rows = [promoted_route]
        elif table_name == "data_resources":
            copied_rows = [
                dict(row) for row in source_rows
                if row.get("data_id") in route_data_ids
            ]
        elif table_name == "decisions":
            copied_rows = [
                dict(row) for row in source_rows
                if row.get("route_id") == route_id
                and row.get("decision") != "approve_release"
            ]
        else:
            copied_rows = [
                dict(row) for row in source_rows
                if row.get("route_id") == route_id
            ]

        if table_name == "figures":
            for row in copied_rows:
                if row.get("status") != "dropped":
                    row["status"] = "mapped"
                row["notes"] = (
                    row.get("notes", "")
                    + " Promoted from pilot; manuscript artifacts must be regenerated."
                ).strip()
        elif table_name == "model_specifications":
            for row in copied_rows:
                if row.get("status") == "verified":
                    row["status"] = "locked"
        elif table_name == "code_candidates":
            for row in copied_rows:
                row["decision"] = "candidate"
                row["verification"] = "inspected"
                row["path_portability"] = "not_checked"
                row["path_test"] = "must be rerun in the manuscript workspace"
        _write_tsv_rows(target / spec["path"], copied_rows)

    write_readiness_report(target)
    return target


def status_summary(project_dir: Path) -> dict[str, Any]:
    project_dir = Path(project_dir)
    manifest = load_manifest(project_dir)
    report = validate_workspace(project_dir)
    readiness = route_readiness(project_dir)
    return {
        "project_id": manifest.get("project_id"),
        "title": manifest.get("title"),
        "workspace_kind": manifest.get("workspace_kind"),
        "stage": manifest.get("stage"),
        "selected_route_id": manifest.get("selected_route_id"),
        "routes": readiness,
        "validation_errors": report.errors,
        "validation_warnings": report.warnings,
    }


def next_actions(project_dir: Path) -> list[str]:
    project_dir = Path(project_dir)
    manifest = load_manifest(project_dir)
    tables, _ = load_tables(project_dir)
    stage = manifest.get("stage")
    actions: list[str] = []

    workspace_kind = manifest.get("workspace_kind")
    if workspace_kind == "pilot":
        if not _document_complete(project_dir / "anchor/audit.md", 300):
            actions.append("Complete anchor/audit.md and remove its TODO marker.")
        if not tables.get("routes"):
            actions.append(
                "Generate the bounded route portfolio and record each route's "
                "capabilities and direction-audit evidence."
            )
            return actions

        promotion_ready: list[str] = []
        for route in tables.get("routes", []):
            if route.get("decision_status") in {"rejected", "stopped"}:
                continue
            route_id = route.get("route_id", "")
            execution_gaps = _execution_gaps_from_tables(project_dir, route, tables)
            if execution_gaps:
                for gap in execution_gaps:
                    actions.append(f"{route_id}: {gap}.")
                continue
            promotion_gaps = _promotion_evidence_gaps(project_dir, route, tables)
            if not promotion_gaps:
                promotion_ready.append(route_id)
            elif route.get("route_role") in {"training", "supporting"}:
                actions.append(
                    f"{route_id}: execution is complete; retain scoped "
                    f"{route.get('route_role')} evidence and do not promote it."
                )
        if promotion_ready and not manifest.get("selected_route_id"):
            actions.append(
                "A route has complete promotion evidence. Record an explicit user "
                "promote or backup decision before creating a manuscript project."
            )
        if not _document_complete(project_dir / "reports/pilot-outcome.md", 500):
            actions.append(
                "Complete reports/pilot-outcome.md with separate paper-side and "
                "product-side conclusions."
            )
        if not actions and stage != "pilot_complete":
            actions.append(
                "Ask the user to promote one route, retain the Pilot as training, "
                "or close it; do not start a manuscript draft inside the pilot."
            )
        return actions

    selected = manifest.get("selected_route_id", "")
    if not selected:
        actions.append("Repair manuscript-project provenance and selected_route_id.")
        return actions
    if not _document_complete(project_dir / "analysis/specification.md", 300):
        actions.append("Complete and freeze the manuscript analysis specification.")
    selected_runs = [
        row for row in tables.get("runs", []) if row.get("route_id") == selected
    ]
    if not selected_runs:
        actions.append("Record and execute the first reproducible manuscript run.")
    selected_results = [
        row for row in tables.get("results", [])
        if row.get("route_id") == selected and row.get("status") == "verified"
    ]
    if selected_runs and not selected_results:
        actions.append("Verify run outputs and register results and claim effects.")
    if selected_results and not _document_complete(
        project_dir / "manuscript/draft.md", 800
    ):
        actions.append("Draft the manuscript from verified artifacts.")
    if not actions and stage != "complete":
        actions.append("Advance the manuscript stage after human review.")
    return actions


def write_readiness_report(project_dir: Path) -> Path:
    project_dir = Path(project_dir)
    manifest = load_manifest(project_dir)
    readiness = route_readiness(project_dir)
    tables, _ = load_tables(project_dir)
    routes_by_id = {
        row.get("route_id", ""): row for row in tables.get("routes", [])
    }
    lines = [
        "# Readiness report",
        "",
        f"- Project: `{manifest.get('project_id')}`",
        f"- Workspace kind: `{manifest.get('workspace_kind')}`",
        f"- Stage: `{manifest.get('stage')}`",
        f"- Selected route: `{manifest.get('selected_route_id') or 'none'}`",
        "",
        "## Routes",
        "",
    ]
    if not readiness:
        lines.append("No candidate routes have been recorded.")
    for item in readiness:
        route = routes_by_id.get(item["route_id"], {})
        lines.extend(
            [
                f"### {item['route_id']}: {item['title']}",
                "",
                f"- Execution ready: `{str(item['execution_ready']).lower()}`",
                "- Promotion evidence complete: "
                f"`{str(item['promotion_evidence_complete']).lower()}`",
                f"- Decision status: `{item['decision_status']}`",
                f"- Evidence stage: `{item['evidence_stage']}`",
                f"- Route role: `{item['route_role']}`",
                f"- Data burden: `{route.get('data_burden')}`",
                f"- Code burden: `{route.get('code_burden')}`",
                f"- Beginner burden: `{route.get('beginner_burden')}`",
                f"- Estimated calendar time: {route.get('estimated_calendar_time')}",
            ]
        )
        if item["execution_gaps"]:
            lines.append("- Execution gaps:")
            lines.extend(f"  - {gap}" for gap in item["execution_gaps"])
        if item["promotion_gaps"] and not item["execution_gaps"]:
            lines.append("- Promotion limits:")
            lines.extend(f"  - {gap}" for gap in item["promotion_gaps"])
        lines.append("")
    lines.extend(["## Next actions", ""])
    lines.extend(f"- {action}" for action in next_actions(project_dir))
    path = project_dir / "reports/readiness.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return path
