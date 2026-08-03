from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .schema import (
    PROJECT_KEYS,
    PROJECT_OBJECT_KEYS,
    SCHEMA_VERSION,
    STAGES,
    TABLES,
)


TODO = "<!-- TODO -->"
DATA_READY = {"sample_parsed", "downloaded", "checksum_verified"}
CODE_READY = {"smoke_passed", "tested"}
FIGURE_RANK = {
    "idea": 0,
    "mapped": 1,
    "spike_generated": 2,
    "generated": 3,
    "verified": 4,
    "blocked": -1,
    "dropped": -1,
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
    if manifest.get("stage") not in STAGES:
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
            "version", "license", "environment", "entrypoint",
            "smoke_input", "smoke_output", "tests_passed",
        )
    ):
        return False
    required_tests = _normalized_items(module.get("required_tests", ""))
    passed_tests = _normalized_items(candidate.get("tests_passed", ""))
    return required_tests.issubset(passed_tests)


def _route_gaps_from_tables(
    route: dict[str, str], tables: dict[str, list[dict[str, str]]]
) -> list[str]:
    route_id = route.get("route_id", "")
    gaps: list[str] = []

    if route.get("status") in {"rejected", "stopped"}:
        return ["route is rejected or stopped"]
    if route.get("science_status") != "pass":
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
        ]
        if requirement.get("independence_required") == "true":
            usable = [
                row for row in usable if row.get("independence") == "independent"
            ]
        if not usable:
            gaps.append(
                f"data requirement {requirement_id} lacks a parsed usable "
                "candidate meeting its contract and minimum subject count"
            )

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
            usable.append(row)
        if not usable:
            gaps.append(
                f"code module {module_id} lacks a qualified smoke-tested donor"
            )

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
    if required_figures and not any(
        FIGURE_RANK.get(row.get("status", ""), -1) >= 2
        for row in required_figures
    ):
        gaps.append("no representative figure or source table was generated")

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

    return gaps


def route_readiness(project_dir: Path) -> list[dict[str, Any]]:
    tables, _ = load_tables(Path(project_dir))
    return [
        {
            "route_id": route.get("route_id"),
            "title": route.get("title"),
            "status": route.get("status"),
            "ready": not _route_gaps_from_tables(route, tables),
            "gaps": _route_gaps_from_tables(route, tables),
        }
        for route in tables.get("routes", [])
    ]


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
            or (
                requirement.get("independence_required") == "true"
                and row.get("independence") != "independent"
            )
        ):
            report.errors.append(
                f"data candidate {row.get('data_id')} is marked use but does "
                "not satisfy access, parsing, contract, sample-size and "
                "independence requirements"
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
        ):
            report.errors.append(
                f"code candidate {row.get('code_id')} is marked use but does "
                "not satisfy contract, provenance, non-interactive execution "
                "and required-test evidence"
            )

    for route in tables.get("routes", []):
        if not _positive_integer(route.get("minimum_main_figures", "")):
            report.errors.append(
                f"route {route.get('route_id')} minimum_main_figures must be positive"
            )
        gaps = _route_gaps_from_tables(route, tables)
        if route.get("status") in {"ready", "selected"} and gaps:
            report.errors.append(
                f"route {route.get('route_id')} is marked {route.get('status')} "
                f"but has gaps: {'; '.join(gaps)}"
            )
        if route.get("status") in {"candidate", "verifying"} and not gaps:
            report.warnings.append(
                f"route {route.get('route_id')} has no readiness gaps and can "
                "be marked ready"
            )

    selected_route_id = manifest.get("selected_route_id", "")
    selected_rows = [
        row for row in tables.get("routes", []) if row.get("status") == "selected"
    ]
    if len(selected_rows) > 1:
        report.errors.append("more than one route is marked selected")
    if selected_route_id:
        selected = indexes.get("routes", {}).get(selected_route_id)
        if not selected:
            report.errors.append(
                f"PROJECT.json: unresolved selected_route_id={selected_route_id}"
            )
        elif selected.get("status") != "selected":
            report.errors.append(
                "PROJECT.json: selected_route_id must point to a selected route"
            )
    if selected_rows and selected_rows[0].get("route_id") != selected_route_id:
        report.errors.append(
            "PROJECT.json selected_route_id does not match the selected route"
        )

    stage = manifest.get("stage")
    if stage in STAGES and stage != "stopped":
        stage_index = STAGES.index(stage)
        if stage_index >= STAGES.index("route_generation"):
            if not _document_complete(project_dir / "anchor/audit.md", 300):
                report.errors.append(
                    "anchor/audit.md must be completed before route generation"
                )
            if not tables.get("routes"):
                report.errors.append("route generation stage requires routes")
        if stage_index >= STAGES.index("selection"):
            if not selected_route_id:
                report.errors.append("selection stage requires selected_route_id")
            elif selected_route_id in indexes.get("routes", {}):
                selected = indexes["routes"][selected_route_id]
                gaps = _route_gaps_from_tables(selected, tables)
                if gaps:
                    report.errors.append(
                        "selected route is not execution-ready: " + "; ".join(gaps)
                    )
                decisions = [
                    row for row in tables.get("decisions", [])
                    if row.get("route_id") == selected_route_id
                    and row.get("decision") == "select"
                ]
                if not decisions:
                    report.errors.append(
                        "selected route lacks a recorded user select decision"
                    )
        if stage_index >= STAGES.index("specification"):
            if not _document_complete(
                project_dir / "analysis/specification.md", 300
            ):
                report.errors.append(
                    "analysis/specification.md must be completed before execution"
                )
        if stage_index >= STAGES.index("execution"):
            selected_runs = [
                row for row in tables.get("runs", [])
                if row.get("route_id") == selected_route_id
            ]
            if not selected_runs:
                report.errors.append("execution stage requires a recorded run")
        if stage_index >= STAGES.index("interpretation"):
            passed_runs = {
                row.get("run_id") for row in tables.get("runs", [])
                if row.get("route_id") == selected_route_id
                and row.get("status") == "passed"
            }
            verified_results = [
                row for row in tables.get("results", [])
                if row.get("route_id") == selected_route_id
                and row.get("run_id") in passed_runs
                and row.get("status") == "verified"
            ]
            if not verified_results:
                report.errors.append(
                    "interpretation stage requires a verified result from a passed run"
                )
        if stage_index >= STAGES.index("writing"):
            if not _document_complete(project_dir / "manuscript/draft.md", 800):
                report.errors.append(
                    "manuscript/draft.md must contain a substantive draft"
                )
        if stage == "complete":
            release = [
                row for row in tables.get("decisions", [])
                if row.get("route_id") == selected_route_id
                and row.get("decision") == "approve_release"
            ]
            if not release:
                report.errors.append(
                    "complete project requires an approve_release decision"
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

## Figure-to-evidence map

For every main and supplementary figure record its manuscript role, data,
metadata, method, code, output, statistical unit and unavailable dependencies.

## Module disposition

Classify every module as retain, repair, substitute, extend, drop or blocked.

## Reproduction boundary

List private data, wet-lab evidence, unavailable code and undocumented author
choices that prevent exact reproduction.
"""


SPECIFICATION_TEMPLATE = f"""# Analysis specification

{TODO}

Complete this only after a route is selected.

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


WORKSPACE_README = """# Paper2Paper workspace

This directory is one isolated paper project. Begin with `anchor/audit.md`,
then use `paper2paper next .` to see evidence gaps. Large data, credentials and
copyrighted source PDFs stay outside Git.
"""


def init_workspace(
    target: Path,
    project_id: str,
    title: str,
    anchor_title: str,
    doi: str = "",
) -> Path:
    target = Path(target)
    if target.exists() and any(target.iterdir()):
        raise ValueError(f"target directory is not empty: {target}")
    for folder in ("anchor", "evidence", "analysis", "execution", "manuscript", "reports"):
        (target / folder).mkdir(parents=True, exist_ok=True)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "project_id": project_id,
        "title": title,
        "stage": "anchor_audit",
        "selected_route_id": "",
        "anchor": {
            "title": anchor_title,
            "doi": doi,
            "citation": "",
            "source_location": "outside Git",
        },
        "goal": {
            "paper_type": "anchor-paper adaptation",
            "target_audience": "beginner-led project with AI assistance",
            "success_definition": (
                "A scientifically defensible, non-duplicate manuscript with "
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
            "Stop workflow engineering when one selected route can be executed, "
            "reviewed and written reliably"
        ),
    }
    (target / "PROJECT.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (target / "README.md").write_text(WORKSPACE_README, encoding="utf-8")
    (target / "anchor/audit.md").write_text(ANCHOR_TEMPLATE, encoding="utf-8")
    (target / "analysis/specification.md").write_text(
        SPECIFICATION_TEMPLATE, encoding="utf-8"
    )
    (target / "manuscript/draft.md").write_text(
        MANUSCRIPT_TEMPLATE, encoding="utf-8"
    )
    for spec in TABLES.values():
        path = target / spec["path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\t".join(spec["columns"]) + "\n", encoding="utf-8")
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

    if not _document_complete(project_dir / "anchor/audit.md", 300):
        actions.append("Complete anchor/audit.md and remove its TODO marker.")
    if not tables.get("routes"):
        actions.append(
            "Generate the route portfolio and add explicit rows to evidence/routes.tsv."
        )
        return actions

    ready_routes: list[str] = []
    for route in tables.get("routes", []):
        gaps = _route_gaps_from_tables(route, tables)
        if not gaps:
            ready_routes.append(route.get("route_id", ""))
        elif route.get("status") not in {"rejected", "stopped"}:
            for gap in gaps:
                actions.append(f"{route.get('route_id')}: {gap}.")

    selected = manifest.get("selected_route_id", "")
    if ready_routes and not selected:
        actions.append(
            "Review ready routes, record one select decision, mark that route "
            "selected and set PROJECT.json selected_route_id."
        )
    if selected:
        if not _document_complete(project_dir / "analysis/specification.md", 300):
            actions.append("Complete and freeze analysis/specification.md.")
        selected_runs = [
            row for row in tables.get("runs", []) if row.get("route_id") == selected
        ]
        if not selected_runs:
            actions.append("Record and execute the first reproducible full run.")
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
        actions.append("Advance the project stage after human review of exit conditions.")
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
                f"- Ready for selection: `{str(item['ready']).lower()}`",
                f"- Recorded status: `{item['status']}`",
                f"- Data burden: `{route.get('data_burden')}`",
                f"- Code burden: `{route.get('code_burden')}`",
                f"- Beginner burden: `{route.get('beginner_burden')}`",
                f"- Estimated calendar time: {route.get('estimated_calendar_time')}",
            ]
        )
        if item["gaps"]:
            lines.append("- Gaps:")
            lines.extend(f"  - {gap}" for gap in item["gaps"])
        lines.append("")
    lines.extend(["## Next actions", ""])
    lines.extend(f"- {action}" for action in next_actions(project_dir))
    path = project_dir / "reports/readiness.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return path
