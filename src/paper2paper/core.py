from __future__ import annotations

import csv
import json
from collections import defaultdict, deque
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path
from typing import Any


@dataclass
class ValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def load_contract() -> dict[str, Any]:
    contract_path = files("paper2paper").joinpath(
        "data", "registry-contract.json"
    )
    with contract_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    return header, rows


def _allowed_values(spec: Any, contract: dict[str, Any]) -> set[str]:
    if spec == "@gates":
        return set(contract["gates"])
    return set(spec)


def _split_reference(value: str, multi: bool) -> list[str]:
    if not value:
        return []
    if not multi:
        return [value]
    return [item.strip() for item in value.split(";") if item.strip()]


def _legacy_manifest_paths(value: Any, prefix: str = "") -> list[str]:
    legacy_keys = {
        "quality_axes",
        "treat_repair_as_non_novel_by_default",
        "require_independent_quality_assessments",
    }
    found: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            path = f"{prefix}.{key}" if prefix else key
            if key in legacy_keys:
                found.append(path)
            found.extend(_legacy_manifest_paths(child, path))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(_legacy_manifest_paths(child, f"{prefix}[{index}]"))
    return found


def load_project(project_dir: Path) -> dict[str, Any]:
    manifest_path = project_dir / "PROJECT.json"
    with manifest_path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("PROJECT.json root must be an object")
    return value


def load_registries(
    project_dir: Path, contract: dict[str, Any]
) -> tuple[dict[str, list[dict[str, str]]], ValidationReport]:
    report = ValidationReport()
    registries: dict[str, list[dict[str, str]]] = {}
    registry_dir = project_dir / "registry"

    for name, spec in contract["registries"].items():
        path = registry_dir / spec["file"]
        if not path.exists():
            report.errors.append(f"missing registry: {path}")
            registries[name] = []
            continue

        header, rows = read_tsv(path)
        missing = [
            column
            for column in spec["required_columns"]
            if column not in header
        ]
        if missing:
            report.errors.append(
                f"{spec['file']}: missing columns {', '.join(missing)}"
            )

        for row_number, row in enumerate(rows, start=2):
            for field_name in spec.get("nonempty_fields", []):
                if not row.get(field_name, ""):
                    report.errors.append(
                        f"{spec['file']}:{row_number}: empty required "
                        f"value {field_name}"
                    )
            for field_name, enum_spec in spec.get("enums", {}).items():
                value = row.get(field_name, "")
                allowed = _allowed_values(enum_spec, contract)
                if value not in allowed:
                    report.errors.append(
                        f"{spec['file']}:{row_number}: invalid "
                        f"{field_name}={value!r}; expected one of "
                        f"{sorted(allowed)}"
                    )
        registries[name] = rows

    return registries, report


def expected_execution_priority(assessment: dict[str, str]) -> str:
    """Derive a transparent execution tier without a weighted total score."""
    validity = assessment.get("scientific_validity")
    data = assessment.get("data_readiness")
    code = assessment.get("code_readiness")
    figures = assessment.get("figure_coverage")
    overlap = assessment.get("publication_overlap")
    allowed_overlap = {
        "clear",
        "adjacent",
        "high_overlap_distinguishable",
    }

    if (
        validity in {"fail", "unknown", None}
        or data in {"unverified", "blocked", None}
        or code in {"unusable", "unknown", None}
        or figures in {"not_mapped", "blocked", None}
        or overlap not in allowed_overlap
    ):
        return "P3"
    if (
        validity == "pass"
        and data == "verified"
        and code == "qualified"
        and figures == "complete"
    ):
        return "P0"
    if (
        validity in {"pass", "conditional"}
        and data == "verified"
        and code in {"qualified", "adaptable"}
        and figures in {"complete", "partial"}
    ):
        return "P1"
    if (
        validity in {"pass", "conditional"}
        and data in {"verified", "partial"}
        and code == "rebuild_required"
        and figures in {"complete", "partial"}
    ):
        return "P2"
    return "P3"


def _validate_manifest(
    manifest: dict[str, Any], contract: dict[str, Any], report: ValidationReport
) -> None:
    required_manifest = {
        "schema_version",
        "project_id",
        "title",
        "project_status",
        "current_gate",
        "anchor_paper_id",
        "active_route_id",
        "publication_goal",
        "adaptation_policy",
        "execution_policy",
        "scope_policy",
        "feedback_policy",
    }
    missing_manifest = sorted(required_manifest - set(manifest))
    if missing_manifest:
        report.errors.append(
            "PROJECT.json: missing keys " + ", ".join(missing_manifest)
        )

    if manifest.get("schema_version") != contract["schema_version"]:
        report.errors.append(
            "PROJECT.json: schema_version does not match installed "
            f"contract {contract['schema_version']}"
        )
    if manifest.get("project_status") not in contract["project_statuses"]:
        report.errors.append(
            "PROJECT.json: invalid project_status "
            f"{manifest.get('project_status')!r}"
        )
    if manifest.get("current_gate") not in contract["gates"]:
        report.errors.append(
            f"PROJECT.json: invalid current_gate "
            f"{manifest.get('current_gate')!r}"
        )

    legacy_paths = _legacy_manifest_paths(manifest)
    if legacy_paths:
        report.errors.append(
            "PROJECT.json: legacy PaperRoute policy keys are not allowed: "
            + ", ".join(sorted(legacy_paths))
        )

    publication_goal = manifest.get("publication_goal", {})
    required_goal_keys = {
        "primary_output",
        "audience",
        "approach",
        "success_definition",
    }
    if not isinstance(publication_goal, dict):
        report.errors.append("PROJECT.json publication_goal must be an object")
    else:
        missing = sorted(required_goal_keys - set(publication_goal))
        if missing:
            report.errors.append(
                "PROJECT.json publication_goal: missing keys "
                + ", ".join(missing)
            )
        for key in required_goal_keys:
            value = publication_goal.get(key)
            if not isinstance(value, str) or not value.strip():
                report.errors.append(
                    f"PROJECT.json publication_goal: {key} must be a "
                    "non-empty string"
                )

    adaptation = manifest.get("adaptation_policy", {})
    required_adaptation_keys = {
        "allowed_route_modes",
        "require_explicit_adaptation_map",
        "prioritize_template_preserving_routes",
        "novelty_required",
        "reject_substantive_duplicate",
        "require_human_route_selection",
    }
    if not isinstance(adaptation, dict):
        report.errors.append("PROJECT.json adaptation_policy must be an object")
    else:
        missing = sorted(required_adaptation_keys - set(adaptation))
        if missing:
            report.errors.append(
                "PROJECT.json adaptation_policy: missing keys "
                + ", ".join(missing)
            )
        modes = adaptation.get("allowed_route_modes")
        if not isinstance(modes, list) or not modes:
            report.errors.append(
                "PROJECT.json adaptation_policy: allowed_route_modes must "
                "be a non-empty list"
            )
        else:
            invalid_modes = sorted(set(modes) - set(contract["route_modes"]))
            if invalid_modes:
                report.errors.append(
                    "PROJECT.json adaptation_policy: invalid route modes "
                    + ", ".join(invalid_modes)
                )
        for key in {
            "require_explicit_adaptation_map",
            "prioritize_template_preserving_routes",
            "reject_substantive_duplicate",
            "require_human_route_selection",
        }:
            if adaptation.get(key) is not True:
                report.errors.append(
                    f"PROJECT.json adaptation_policy: {key} must be true"
                )
        if adaptation.get("novelty_required") is not False:
            report.errors.append(
                "PROJECT.json adaptation_policy: novelty_required must be "
                "false; novelty is not a Paper2Paper gate"
            )

    execution = manifest.get("execution_policy", {})
    required_execution_keys = {
        "require_verified_data_for_approval",
        "require_qualified_code_for_approval",
        "require_complete_figure_plan_for_approval",
        "allowed_approval_code_grades",
        "allowed_approval_overlap",
        "minimum_required_figures",
        "priority_order",
    }
    if not isinstance(execution, dict):
        report.errors.append("PROJECT.json execution_policy must be an object")
    else:
        missing = sorted(required_execution_keys - set(execution))
        if missing:
            report.errors.append(
                "PROJECT.json execution_policy: missing keys "
                + ", ".join(missing)
            )
        for key in {
            "require_verified_data_for_approval",
            "require_qualified_code_for_approval",
            "require_complete_figure_plan_for_approval",
        }:
            if execution.get(key) is not True:
                report.errors.append(
                    f"PROJECT.json execution_policy: {key} must be true"
                )
        allowed_grades = execution.get("allowed_approval_code_grades")
        if allowed_grades != ["A", "B"]:
            report.errors.append(
                "PROJECT.json execution_policy: "
                "allowed_approval_code_grades must be ['A', 'B']"
            )
        allowed_overlap = execution.get("allowed_approval_overlap")
        expected_overlap = [
            "clear",
            "adjacent",
            "high_overlap_distinguishable",
        ]
        if allowed_overlap != expected_overlap:
            report.errors.append(
                "PROJECT.json execution_policy: allowed_approval_overlap "
                "must preserve the three non-duplicate states"
            )
        minimum_figures = execution.get("minimum_required_figures")
        if not isinstance(minimum_figures, int) or minimum_figures < 1:
            report.errors.append(
                "PROJECT.json execution_policy: minimum_required_figures "
                "must be a positive integer"
            )
        expected_order = [
            "code_readiness",
            "data_readiness",
            "figure_coverage",
            "scientific_validity",
            "anchor_reuse",
            "beginner_burden",
            "publication_overlap",
        ]
        if execution.get("priority_order") != expected_order:
            report.errors.append(
                "PROJECT.json execution_policy: priority_order must match "
                "the binding execution-first order"
            )

    scope = manifest.get("scope_policy", {})
    required_scope_keys = {
        "allowed_work_reasons",
        "optimization_stop_rule",
        "require_work_item_for_modules",
        "require_manuscript_implication_for_results",
    }
    if not isinstance(scope, dict):
        report.errors.append("PROJECT.json scope_policy must be an object")
    else:
        missing = sorted(required_scope_keys - set(scope))
        if missing:
            report.errors.append(
                "PROJECT.json scope_policy: missing keys " + ", ".join(missing)
            )
        if set(scope.get("allowed_work_reasons", [])) != set(
            contract["work_reasons"]
        ):
            report.errors.append(
                "PROJECT.json scope_policy: allowed_work_reasons must match "
                "the binding project discipline"
            )
        stop_rule = scope.get("optimization_stop_rule")
        if not isinstance(stop_rule, str) or not stop_rule.strip():
            report.errors.append(
                "PROJECT.json scope_policy: optimization_stop_rule must be "
                "a non-empty string"
            )
        for key in {
            "require_work_item_for_modules",
            "require_manuscript_implication_for_results",
        }:
            if scope.get(key) is not True:
                report.errors.append(
                    f"PROJECT.json scope_policy: {key} must be true"
                )

    feedback = manifest.get("feedback_policy", {})
    if not isinstance(feedback, dict):
        report.errors.append("PROJECT.json feedback_policy must be an object")
    else:
        effects = set(feedback.get("effects_requiring_change_request", []))
        if effects != {"refine", "reroute", "stop"}:
            report.errors.append(
                "PROJECT.json feedback_policy: "
                "effects_requiring_change_request must be refine, reroute, stop"
            )
        for key in {
            "invalidate_downstream_on_approved_change",
            "require_review_before_high_impact_change",
        }:
            if feedback.get(key) is not True:
                report.errors.append(
                    f"PROJECT.json feedback_policy: {key} must be true"
                )


def _validate_code_qualification(
    registries: dict[str, list[dict[str, str]]], report: ValidationReport
) -> None:
    for row in registries.get("code_sources", []):
        grade = row.get("qualification_grade")
        name = row.get("code_id", "(unknown)")
        if grade not in {"A", "B"}:
            continue
        required_text = ["version_ref", "license", "environment", "entrypoint"]
        missing = [field for field in required_text if not row.get(field)]
        if missing:
            report.errors.append(
                f"code source {name} grade {grade} lacks: "
                + ", ".join(missing)
            )
        if row.get("noninteractive") != "true":
            report.errors.append(
                f"code source {name} grade {grade} must be noninteractive"
            )
        if row.get("private_dependencies") != "false":
            report.errors.append(
                f"code source {name} grade {grade} cannot depend on private inputs"
            )
        if row.get("install_status") != "passed":
            report.errors.append(
                f"code source {name} grade {grade} requires a passed install"
            )
        if row.get("smoke_test") != "passed":
            report.errors.append(
                f"code source {name} grade {grade} requires a passed smoke test"
            )
        allowed_tests = {"full_pass"} if grade == "A" else {"smoke_pass", "full_pass"}
        if row.get("test_status") not in allowed_tests:
            report.errors.append(
                f"code source {name} grade {grade} has insufficient tests"
            )
        if grade == "A" and row.get("hardcoded_paths") != "false":
            report.errors.append(
                f"code source {name} grade A cannot contain hardcoded paths"
            )


def _validate_route_approval(
    route: dict[str, str],
    assessment: dict[str, str],
    registries: dict[str, list[dict[str, str]]],
    manifest: dict[str, Any],
    report: ValidationReport,
) -> None:
    route_id = route["route_id"]
    failed = []
    if assessment.get("status") != "approved":
        failed.append("assessment_status")
    if assessment.get("scientific_validity") != "pass":
        failed.append("scientific_validity")
    if assessment.get("data_readiness") != "verified":
        failed.append("data_readiness")
    if assessment.get("code_readiness") != "qualified":
        failed.append("code_readiness")
    if assessment.get("figure_coverage") != "complete":
        failed.append("figure_coverage")
    allowed_overlap = set(
        manifest.get("execution_policy", {}).get(
            "allowed_approval_overlap", []
        )
    )
    if assessment.get("publication_overlap") not in allowed_overlap:
        failed.append("publication_overlap")
    if expected_execution_priority(assessment) != "P0":
        failed.append("execution_priority")
    if failed:
        report.errors.append(
            f"approved route {route_id} fails execution gates: "
            + ", ".join(failed)
        )

    datasets = {
        row["dataset_id"]: row
        for row in registries.get("datasets", [])
        if row.get("dataset_id")
    }
    required_dataset_maps = [
        row
        for row in registries.get("dataset_route_map", [])
        if row.get("route_id") == route_id and row.get("required") == "true"
    ]
    if not required_dataset_maps:
        report.errors.append(
            f"approved route {route_id} has no required dataset mapping"
        )
    download_ready = {"sample_verified", "downloaded", "checksum_verified"}
    for mapping in required_dataset_maps:
        dataset = datasets.get(mapping.get("dataset_id", ""), {})
        if mapping.get("availability_verdict") != "verified":
            report.errors.append(
                f"approved route {route_id} has unverified required dataset "
                f"map {mapping.get('dataset_map_id')}"
            )
        if dataset.get("status") != "accepted":
            report.errors.append(
                f"approved route {route_id} uses a required dataset that is "
                f"not accepted: {mapping.get('dataset_id')}"
            )
        if dataset.get("download_status") not in download_ready:
            report.errors.append(
                f"approved route {route_id} lacks a verified sample or "
                f"download for {mapping.get('dataset_id')}"
            )

    code_sources = {
        row["code_id"]: row
        for row in registries.get("code_sources", [])
        if row.get("code_id")
    }
    required_code_maps = [
        row
        for row in registries.get("code_module_map", [])
        if row.get("route_id") == route_id and row.get("required") == "true"
    ]
    if not required_code_maps:
        report.errors.append(
            f"approved route {route_id} has no required code mapping"
        )
    allowed_grades = set(
        manifest.get("execution_policy", {}).get(
            "allowed_approval_code_grades", []
        )
    )
    for mapping in required_code_maps:
        source = code_sources.get(mapping.get("code_id", ""), {})
        if mapping.get("status") != "verified":
            report.errors.append(
                f"approved route {route_id} has unverified code map "
                f"{mapping.get('code_map_id')}"
            )
        if not mapping.get("source_path"):
            report.errors.append(
                f"approved route {route_id} code map "
                f"{mapping.get('code_map_id')} lacks source_path"
            )
        if source.get("status") != "accepted":
            report.errors.append(
                f"approved route {route_id} uses code that is not accepted: "
                f"{mapping.get('code_id')}"
            )
        if source.get("qualification_grade") not in allowed_grades:
            report.errors.append(
                f"approved route {route_id} uses unqualified code: "
                f"{mapping.get('code_id')}"
            )

    required_figures = [
        row
        for row in registries.get("figure_plan", [])
        if row.get("route_id") == route_id and row.get("required") == "true"
    ]
    minimum_figures = manifest.get("execution_policy", {}).get(
        "minimum_required_figures", 1
    )
    if len(required_figures) < minimum_figures:
        report.errors.append(
            f"approved route {route_id} has {len(required_figures)} required "
            f"figures; minimum is {minimum_figures}"
        )
    for figure in required_figures:
        if figure.get("status") not in {"ready", "generated", "verified"}:
            report.errors.append(
                f"approved route {route_id} has unready required figure "
                f"{figure.get('figure_id')}"
            )

    overlap_rows = [
        row
        for row in registries.get("publication_overlap", [])
        if row.get("route_id") == route_id
    ]
    if not overlap_rows:
        report.errors.append(
            f"approved route {route_id} lacks a publication-overlap record"
        )
    if any(
        row.get("overlap_level") == "duplicate"
        or row.get("decision") == "stop"
        for row in overlap_rows
    ):
        report.errors.append(
            f"approved route {route_id} has a blocking publication duplicate"
        )


def validate_project(project_dir: Path) -> ValidationReport:
    project_dir = Path(project_dir)
    report = ValidationReport()
    contract = load_contract()

    manifest_path = project_dir / "PROJECT.json"
    if not manifest_path.exists():
        report.errors.append(f"missing project manifest: {manifest_path}")
        return report

    try:
        manifest = load_project(project_dir)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        report.errors.append(f"invalid PROJECT.json: {exc}")
        return report

    _validate_manifest(manifest, contract, report)
    registries, registry_report = load_registries(project_dir, contract)
    report.errors.extend(registry_report.errors)
    report.warnings.extend(registry_report.warnings)

    ids_by_registry: dict[str, set[str]] = {}
    global_ids: dict[str, str] = {}
    for name, spec in contract["registries"].items():
        id_field = spec.get("id_field", "")
        ids: set[str] = set()
        if id_field:
            for row_number, row in enumerate(
                registries.get(name, []), start=2
            ):
                entity_id = row.get(id_field, "")
                if not entity_id:
                    report.errors.append(
                        f"{spec['file']}:{row_number}: empty {id_field}"
                    )
                    continue
                if entity_id in ids:
                    report.errors.append(
                        f"{spec['file']}:{row_number}: duplicate {entity_id}"
                    )
                if entity_id in global_ids:
                    report.errors.append(
                        f"{spec['file']}:{row_number}: ID {entity_id} "
                        f"already used in {global_ids[entity_id]}"
                    )
                ids.add(entity_id)
                global_ids[entity_id] = name
        ids_by_registry[name] = ids

    for name, spec in contract["registries"].items():
        for row_number, row in enumerate(registries.get(name, []), start=2):
            for reference in spec.get("references", []):
                field_name = reference["field"]
                values = _split_reference(
                    row.get(field_name, ""), reference.get("multi", False)
                )
                if not values and not reference.get("optional", False):
                    report.errors.append(
                        f"{spec['file']}:{row_number}: empty required "
                        f"reference {field_name}"
                    )
                    continue
                for value in values:
                    target = reference["target"]
                    valid = (
                        value in global_ids
                        if target == "*"
                        else value in ids_by_registry.get(target, set())
                    )
                    if not valid:
                        report.errors.append(
                            f"{spec['file']}:{row_number}: unresolved "
                            f"{field_name}={value}"
                        )

    anchor_id = manifest.get("anchor_paper_id", "")
    if anchor_id and anchor_id not in ids_by_registry.get("papers", set()):
        report.errors.append(
            f"PROJECT.json: unresolved anchor_paper_id={anchor_id!r}"
        )
    if manifest.get("project_status") != "intake" and not anchor_id:
        report.errors.append(
            "PROJECT.json: non-intake project requires anchor_paper_id"
        )

    route_rows = {
        row["route_id"]: row
        for row in registries.get("routes", [])
        if row.get("route_id")
    }
    active_route_id = manifest.get("active_route_id", "")
    if active_route_id and active_route_id not in route_rows:
        report.errors.append(
            f"PROJECT.json: unresolved active_route_id={active_route_id!r}"
        )
    if manifest.get("project_status") in {"active", "paused", "completed"}:
        if not active_route_id:
            report.errors.append(
                "PROJECT.json: active, paused, or completed project requires "
                "active_route_id"
            )
        elif route_rows.get(active_route_id, {}).get("status") not in {
            "approved",
            "active",
        }:
            report.errors.append(
                "PROJECT.json: an active, paused, or completed project must "
                "point to an approved or active route"
            )
    if active_route_id in route_rows and route_rows[active_route_id].get(
        "status"
    ) in {"blocked", "superseded", "rejected"}:
        report.errors.append(
            "PROJECT.json: active route is blocked, superseded, or rejected"
        )
    active_route_rows = [
        row
        for row in route_rows.values()
        if row.get("status") == "active"
    ]
    if len(active_route_rows) > 1:
        report.errors.append("routes.tsv: more than one route is active")
    if active_route_rows and active_route_rows[0].get("route_id") != active_route_id:
        report.errors.append(
            "PROJECT.json: active_route_id does not match the route marked active"
        )

    configured_modes = set(
        manifest.get("adaptation_policy", {}).get(
            "allowed_route_modes", []
        )
    )
    adaptations_by_route: dict[str, list[dict[str, str]]] = defaultdict(list)
    for adaptation in registries.get("route_adaptations", []):
        adaptations_by_route[adaptation.get("route_id", "")].append(adaptation)

    assessments_by_route: dict[str, list[dict[str, str]]] = defaultdict(list)
    overlap_by_route: dict[str, list[dict[str, str]]] = defaultdict(list)
    for overlap in registries.get("publication_overlap", []):
        overlap_by_route[overlap.get("route_id", "")].append(overlap)
    overlap_rank = {
        "clear": 0,
        "adjacent": 1,
        "high_overlap_distinguishable": 2,
        "duplicate": 3,
    }
    for assessment in registries.get("route_assessments", []):
        route_id = assessment.get("route_id", "")
        assessments_by_route[route_id].append(assessment)
        expected = expected_execution_priority(assessment)
        if assessment.get("execution_priority") != expected:
            message = (
                f"route assessment {assessment.get('assessment_id')} records "
                f"{assessment.get('execution_priority')} but execution rules "
                f"derive {expected}"
            )
            route = route_rows.get(route_id, {})
            if route.get("status") in {"approved", "active"}:
                report.errors.append(message)
            else:
                report.warnings.append(message)
        route_overlap_rows = overlap_by_route.get(route_id, [])
        if route_overlap_rows:
            derived_overlap = max(
                (
                    row.get("overlap_level", "clear")
                    for row in route_overlap_rows
                ),
                key=lambda value: overlap_rank.get(value, -1),
            )
            if assessment.get("publication_overlap") != derived_overlap:
                message = (
                    f"route assessment {assessment.get('assessment_id')} "
                    f"records publication_overlap="
                    f"{assessment.get('publication_overlap')!r} but overlap "
                    f"records derive {derived_overlap!r}"
                )
                if route_rows.get(route_id, {}).get("status") in {
                    "approved",
                    "active",
                }:
                    report.errors.append(message)
                else:
                    report.warnings.append(message)
        if (
            assessment.get("publication_overlap") == "duplicate"
            and route_rows.get(route_id, {}).get("status")
            not in {"rejected", "blocked", "superseded"}
        ):
            report.errors.append(
                f"route {route_id} is a substantive publication duplicate "
                "and must be rejected or blocked"
            )

    completed_approval_reviews = {
        row["review_id"]: row
        for row in registries.get("reviews", [])
        if row.get("status") == "completed"
        and row.get("outcome") == "approve"
        and row.get("review_id")
    }

    for route_id, route in route_rows.items():
        if route.get("route_mode") not in configured_modes:
            report.errors.append(
                f"route {route_id} uses project-disallowed mode "
                f"{route.get('route_mode')!r}"
            )
        if (
            route.get("route_mode") != "faithful_reproduction"
            and not adaptations_by_route.get(route_id)
        ):
            report.errors.append(
                f"route {route_id} requires at least one explicit adaptation"
            )
        assessments = assessments_by_route.get(route_id, [])
        if route.get("status") in {
            "triaged",
            "spike_passed",
            "approved",
            "active",
        } and len(assessments) != 1:
            report.errors.append(
                f"route {route_id} requires exactly one assessment; "
                f"found {len(assessments)}"
            )
        if route.get("status") in {"approved", "active"}:
            review_id = route.get("approved_review_id", "")
            approval = completed_approval_reviews.get(review_id)
            if not approval or approval.get("target_id") != route_id:
                report.errors.append(
                    f"approved route {route_id} lacks a matching completed "
                    "approval review"
                )
            if len(assessments) == 1:
                _validate_route_approval(
                    route, assessments[0], registries, manifest, report
                )

    _validate_code_qualification(registries, report)

    for overlap in registries.get("publication_overlap", []):
        if (
            overlap.get("overlap_level") == "duplicate"
            and overlap.get("decision") != "stop"
        ):
            report.errors.append(
                f"publication overlap {overlap.get('overlap_id')} is duplicate "
                "and must use decision=stop"
            )
        if (
            overlap.get("overlap_level") == "duplicate"
            and route_rows.get(overlap.get("route_id", ""), {}).get("status")
            not in {"blocked", "rejected", "superseded"}
        ):
            report.errors.append(
                f"route {overlap.get('route_id')} has a substantive duplicate "
                "record and must be blocked or rejected"
            )

    results_needing_change = {
        row["result_id"]
        for row in registries.get("results", [])
        if row.get("route_effect") in {"refine", "reroute", "stop"}
    }
    linked_results = {
        row.get("source_result_id", "")
        for row in registries.get("change_requests", [])
    }
    for result_id in sorted(results_needing_change - linked_results):
        report.errors.append(
            f"result {result_id} has a route effect but no change request"
        )

    for row in registries.get("decisions", []):
        if (
            row.get("status") == "approved"
            and row.get("requires_review") == "true"
            and row.get("approved_review_id")
            not in completed_approval_reviews
        ):
            report.errors.append(
                f"approved decision {row.get('decision_id')} lacks an "
                "approving completed review"
            )
    for row in registries.get("change_requests", []):
        if (
            row.get("status") in {"approved", "implemented", "closed"}
            and row.get("severity") in {"high", "critical"}
            and row.get("approved_review_id")
            not in completed_approval_reviews
        ):
            report.errors.append(
                f"high-impact change {row.get('change_id')} lacks an "
                "approving completed review"
            )

    return report


def init_project(target: Path, project_id: str, title: str) -> Path:
    target = Path(target)
    if target.exists() and any(target.iterdir()):
        raise ValueError(f"target directory is not empty: {target}")

    contract = load_contract()
    target.mkdir(parents=True, exist_ok=True)
    registry_dir = target / "registry"
    registry_dir.mkdir(exist_ok=True)
    for folder in (
        "config",
        "src",
        "tests",
        "runs",
        "outputs",
        "reports",
        "reviews",
    ):
        (target / folder).mkdir(exist_ok=True)

    manifest = {
        "schema_version": contract["schema_version"],
        "project_id": project_id,
        "title": title,
        "project_status": "intake",
        "current_gate": "G0_SCOPE",
        "anchor_paper_id": "",
        "active_route_id": "",
        "publication_goal": {
            "primary_output": "A scientifically defensible manuscript",
            "audience": "Beginner-led project with AI assistance",
            "approach": (
                "Reproducible adaptation of an anchor paper using accessible "
                "data and qualified code"
            ),
            "success_definition": (
                "A non-duplicate manuscript draft whose figures, source tables, "
                "methods, inputs, code, runs, and limitations are traceable"
            ),
        },
        "adaptation_policy": {
            "allowed_route_modes": list(contract["route_modes"]),
            "require_explicit_adaptation_map": True,
            "prioritize_template_preserving_routes": True,
            "novelty_required": False,
            "reject_substantive_duplicate": True,
            "require_human_route_selection": True,
        },
        "execution_policy": {
            "require_verified_data_for_approval": True,
            "require_qualified_code_for_approval": True,
            "require_complete_figure_plan_for_approval": True,
            "allowed_approval_code_grades": ["A", "B"],
            "allowed_approval_overlap": [
                "clear",
                "adjacent",
                "high_overlap_distinguishable",
            ],
            "minimum_required_figures": 1,
            "priority_order": [
                "code_readiness",
                "data_readiness",
                "figure_coverage",
                "scientific_validity",
                "anchor_reuse",
                "beginner_burden",
                "publication_overlap",
            ],
        },
        "scope_policy": {
            "allowed_work_reasons": list(contract["work_reasons"]),
            "optimization_stop_rule": (
                "Stop workflow engineering when the selected manuscript can "
                "be generated, checked, and reviewed reliably"
            ),
            "require_work_item_for_modules": True,
            "require_manuscript_implication_for_results": True,
        },
        "feedback_policy": {
            "effects_requiring_change_request": [
                "refine",
                "reroute",
                "stop",
            ],
            "invalidate_downstream_on_approved_change": True,
            "require_review_before_high_impact_change": True,
        },
    }
    (target / "PROJECT.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    for spec in contract["registries"].values():
        path = registry_dir / spec["file"]
        path.write_text(
            "\t".join(spec["required_columns"]) + "\n",
            encoding="utf-8",
        )

    return target


def downstream_impact(project_dir: Path, changed_id: str) -> list[dict[str, Any]]:
    contract = load_contract()
    registries, report = load_registries(Path(project_dir), contract)
    if report.errors:
        raise ValueError("; ".join(report.errors))

    adjacency: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in registries.get("dependencies", []):
        if row.get("active") == "true":
            adjacency[row["upstream_id"]].append(
                (row["downstream_id"], row["relation"])
            )

    queue: deque[tuple[str, int, str]] = deque()
    for downstream, relation in adjacency.get(changed_id, []):
        queue.append((downstream, 1, relation))

    seen = {changed_id}
    impacted: list[dict[str, Any]] = []
    while queue:
        entity_id, depth, relation = queue.popleft()
        if entity_id in seen:
            continue
        seen.add(entity_id)
        impacted.append(
            {
                "entity_id": entity_id,
                "depth": depth,
                "via_relation": relation,
            }
        )
        for downstream, child_relation in adjacency.get(entity_id, []):
            queue.append((downstream, depth + 1, child_relation))

    return impacted


def route_priorities(project_dir: Path) -> list[dict[str, Any]]:
    contract = load_contract()
    registries, report = load_registries(Path(project_dir), contract)
    if report.errors:
        raise ValueError("; ".join(report.errors))
    return [
        {
            "route_id": row.get("route_id"),
            "recorded": row.get("execution_priority"),
            "expected": expected_execution_priority(row),
            "matches": row.get("execution_priority")
            == expected_execution_priority(row),
        }
        for row in registries.get("route_assessments", [])
    ]


def status_summary(project_dir: Path) -> dict[str, Any]:
    project_dir = Path(project_dir)
    manifest = load_project(project_dir)
    contract = load_contract()
    registries, _ = load_registries(project_dir, contract)
    validation = validate_project(project_dir)

    return {
        "project_id": manifest.get("project_id"),
        "title": manifest.get("title"),
        "project_status": manifest.get("project_status"),
        "current_gate": manifest.get("current_gate"),
        "active_route_id": manifest.get("active_route_id"),
        "primary_output": manifest.get("publication_goal", {}).get(
            "primary_output"
        ),
        "entity_counts": {
            name: len(rows) for name, rows in registries.items()
        },
        "pending_reviews": [
            row.get("review_id")
            for row in registries.get("reviews", [])
            if row.get("status") == "pending"
        ],
        "pending_decisions": [
            row.get("decision_id")
            for row in registries.get("decisions", [])
            if row.get("status") == "proposed"
        ],
        "open_change_requests": [
            row.get("change_id")
            for row in registries.get("change_requests", [])
            if row.get("status") in {"proposed", "approved", "implemented"}
        ],
        "open_work_items": [
            row.get("work_id")
            for row in registries.get("work_items", [])
            if row.get("status")
            in {"proposed", "approved", "in_progress", "blocked"}
        ],
        "route_priorities": route_priorities(project_dir),
        "validation_errors": validation.errors,
        "validation_warnings": validation.warnings,
    }
