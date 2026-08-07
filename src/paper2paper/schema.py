from __future__ import annotations

from typing import Any


SCHEMA_VERSION = "2.0.0"

PILOT_STAGES = (
    "anchor_audit",
    "route_generation",
    "verification",
    "pilot_review",
    "pilot_complete",
    "stopped",
)

MANUSCRIPT_STAGES = (
    "specification",
    "execution",
    "interpretation",
    "writing",
    "complete",
    "stopped",
)

STAGES = tuple(dict.fromkeys((*PILOT_STAGES, *MANUSCRIPT_STAGES)))

ROUTE_MODES = (
    "reproduction",
    "marker",
    "gene_set",
    "cell_type",
    "cancer_type",
    "pan_cancer",
    "signature",
    "combined",
    "custom",
)

ROUTE_ROLES = (
    "manuscript_candidate",
    "training",
    "supporting",
)

PROJECT_KEYS = {
    "schema_version",
    "workspace_kind",
    "project_id",
    "title",
    "stage",
    "selected_route_id",
    "anchor",
    "provenance",
    "goal",
    "constraints",
    "stop_rule",
}

PROJECT_OBJECT_KEYS = {
    "anchor": {"title", "doi", "citation", "source_location"},
    "provenance": {"source_pilot_id", "source_route_id", "source_run_ids"},
    "goal": {"paper_type", "target_audience", "success_definition"},
    "constraints": {"beginner_led", "deadline", "compute", "skills"},
}


def _spec(
    path: str,
    id_field: str,
    columns: tuple[str, ...],
    nonempty: tuple[str, ...],
    enums: dict[str, tuple[str, ...]] | None = None,
    references: tuple[dict[str, Any], ...] = (),
) -> dict[str, Any]:
    return {
        "path": path,
        "id_field": id_field,
        "columns": columns,
        "nonempty": nonempty,
        "enums": enums or {},
        "references": references,
    }


TABLES = {
    "routes": _spec(
        "evidence/routes.tsv",
        "route_id",
        (
            "route_id", "decision_status", "evidence_stage", "route_role",
            "mode", "capability_ids", "title", "question", "target_disease",
            "target_object", "biological_unit", "comparison",
            "primary_outcome", "claim_ceiling", "falsifier",
            "anchor_reuse", "changed_axes", "scientific_review",
            "science_basis", "minimum_main_figures", "data_burden",
            "code_burden", "beginner_burden", "estimated_calendar_time",
            "model_spec_required", "main_risk", "stop_condition",
        ),
        (
            "decision_status", "evidence_stage", "route_role", "mode",
            "capability_ids", "title", "question", "target_disease",
            "target_object", "biological_unit", "comparison",
            "primary_outcome", "claim_ceiling", "falsifier", "anchor_reuse",
            "changed_axes", "scientific_review", "science_basis",
            "minimum_main_figures", "data_burden", "code_burden",
            "beginner_burden", "estimated_calendar_time", "main_risk",
            "model_spec_required", "stop_condition",
        ),
        {
            "decision_status": (
                "active", "backup", "promoted", "rejected", "stopped",
            ),
            "evidence_stage": (
                "direction_audited", "availability_prechecked",
                "minimal_real_run", "figure_loop_closed",
            ),
            "route_role": ROUTE_ROLES,
            "mode": ROUTE_MODES,
            "scientific_review": ("unreviewed", "conditional", "passed", "failed"),
            "data_burden": ("low", "medium", "high", "unknown"),
            "code_burden": ("low", "medium", "high", "unknown"),
            "beginner_burden": ("low", "medium", "high", "unknown"),
            "model_spec_required": ("true", "false"),
        },
    ),
    "search_log": _spec(
        "evidence/search_log.tsv",
        "search_id",
        (
            "search_id", "route_id", "domain", "source", "query",
            "searched_at", "result_count", "outcome", "notes",
        ),
        (
            "route_id", "domain", "source", "query", "searched_at",
            "result_count", "outcome",
        ),
        {
            "domain": ("data", "code", "literature"),
            "outcome": ("continue", "no_result", "stop"),
        },
        ({"field": "route_id", "target": "routes"},),
    ),
    "data_requirements": _spec(
        "evidence/data_requirements.tsv",
        "requirement_id",
        (
            "requirement_id", "route_id", "role", "disease", "tissue",
            "modality", "biological_unit", "comparison_or_outcome",
            "required_fields", "minimum_subjects", "access_limit",
            "independence_required", "required", "notes",
        ),
        (
            "route_id", "role", "disease", "tissue", "modality",
            "biological_unit", "comparison_or_outcome", "required_fields",
            "minimum_subjects", "access_limit", "independence_required",
            "required",
        ),
        {
            "independence_required": ("true", "false"),
            "required": ("true", "false"),
        },
        ({"field": "route_id", "target": "routes"},),
    ),
    "data_candidates": _spec(
        "evidence/data_candidates.tsv",
        "data_id",
        (
            "data_id", "route_id", "requirement_id", "name", "accession",
            "uri", "source", "disease", "tissue", "modality", "subjects",
            "biological_unit", "fields_checked", "contract_match",
            "access", "verification", "independence", "decision",
            "checked_at", "checksum", "notes",
        ),
        (
            "route_id", "requirement_id", "name", "uri", "source",
            "disease", "tissue", "modality", "subjects", "biological_unit",
            "fields_checked", "contract_match", "access",
            "verification", "independence", "decision", "checked_at",
        ),
        {
            "access": (
                "public", "controlled", "request_only", "local",
                "unavailable", "unknown",
            ),
            "contract_match": ("true", "false", "unknown"),
            "verification": (
                "not_checked", "metadata_checked", "sample_parsed",
                "downloaded", "checksum_verified", "blocked",
            ),
            "independence": (
                "independent", "overlap", "unknown", "not_applicable",
            ),
            "decision": ("candidate", "use", "reject", "blocked"),
        },
        (
            {"field": "route_id", "target": "routes"},
            {"field": "requirement_id", "target": "data_requirements"},
        ),
    ),
    "data_resources": _spec(
        "evidence/data_resources.tsv",
        "resource_id",
        (
            "resource_id", "data_id", "name", "role", "uri",
            "source_version", "checked_at", "verification", "local_name",
            "bytes", "checksum", "fields_supplied", "identifier_field",
            "notes",
        ),
        (
            "data_id", "name", "role", "uri", "source_version",
            "checked_at", "verification", "local_name", "fields_supplied",
            "identifier_field",
        ),
        {
            "verification": (
                "not_checked", "metadata_checked", "sample_parsed",
                "downloaded", "checksum_verified", "blocked",
            ),
        },
        ({"field": "data_id", "target": "data_candidates"},),
    ),
    "cohort_usage": _spec(
        "evidence/cohort_usage.tsv",
        "usage_id",
        (
            "usage_id", "route_id", "data_id", "cohort_key",
            "analysis_step", "role", "outcome_used", "features_influenced",
            "parameters_influenced", "cutoff_influenced",
            "claimed_external_validation", "acceptable", "notes",
        ),
        (
            "route_id", "data_id", "cohort_key", "analysis_step", "role",
            "outcome_used", "features_influenced", "parameters_influenced",
            "cutoff_influenced", "claimed_external_validation", "acceptable",
        ),
        {
            "role": (
                "discovery", "feature_screening", "model_fitting",
                "cutoff_selection", "internal_validation",
                "external_validation", "reproduction", "sensitivity",
                "localization", "exploratory", "other",
            ),
            "outcome_used": ("true", "false"),
            "features_influenced": ("true", "false"),
            "parameters_influenced": ("true", "false"),
            "cutoff_influenced": ("true", "false"),
            "claimed_external_validation": ("true", "false"),
            "acceptable": ("true", "false"),
        },
        (
            {"field": "route_id", "target": "routes"},
            {"field": "data_id", "target": "data_candidates"},
        ),
    ),
    "code_requirements": _spec(
        "evidence/code_requirements.tsv",
        "module_id",
        (
            "module_id", "route_id", "capability_id", "name", "purpose",
            "input_contract", "output_contract", "required_tests", "required",
            "notes",
        ),
        (
            "route_id", "capability_id", "name", "purpose", "input_contract",
            "output_contract", "required_tests", "required",
        ),
        {"required": ("true", "false")},
        ({"field": "route_id", "target": "routes"},),
    ),
    "code_candidates": _spec(
        "evidence/code_candidates.tsv",
        "code_id",
        (
            "code_id", "route_id", "module_id", "name", "source_type",
            "uri", "version", "license", "language", "environment",
            "entrypoint", "contract_match", "noninteractive",
            "private_inputs", "hardcoded_paths", "path_portability",
            "path_test", "verification", "decision", "checked_at",
            "reusable_release_id", "smoke_input", "smoke_output",
            "tests_passed", "notes",
        ),
        (
            "route_id", "module_id", "name", "source_type", "uri",
            "language", "noninteractive", "private_inputs", "hardcoded_paths",
            "path_portability", "path_test", "contract_match", "verification",
            "decision", "checked_at",
        ),
        {
            "source_type": (
                "author", "official", "data_paper", "method_paper",
                "related_paper", "reconstructed", "other",
            ),
            "contract_match": ("true", "false", "unknown"),
            "noninteractive": ("true", "false", "unknown"),
            "private_inputs": ("true", "false", "unknown"),
            "hardcoded_paths": ("true", "false", "unknown"),
            "path_portability": (
                "not_checked", "target_environment_passed",
                "staged_workaround", "blocked",
            ),
            "verification": (
                "not_checked", "inspected", "install_passed",
                "smoke_passed", "tested", "blocked",
            ),
            "decision": ("candidate", "use", "reject", "blocked"),
        },
        (
            {"field": "route_id", "target": "routes"},
            {"field": "module_id", "target": "code_requirements"},
        ),
    ),
    "model_specifications": _spec(
        "evidence/model_specifications.tsv",
        "model_id",
        (
            "model_id", "route_id", "name", "status", "feature_order",
            "coefficients", "input_scale", "feature_mapping",
            "duplicate_feature_policy", "missing_feature_policy",
            "normalization_reference", "cutoff_rule", "output_definition",
            "locked_at", "source", "notes",
        ),
        (
            "route_id", "name", "status", "feature_order", "coefficients",
            "input_scale", "feature_mapping", "duplicate_feature_policy",
            "missing_feature_policy", "normalization_reference",
            "cutoff_rule", "output_definition", "locked_at", "source",
        ),
        {"status": ("draft", "locked", "verified", "blocked")},
        ({"field": "route_id", "target": "routes"},),
    ),
    "figures": _spec(
        "evidence/figures.tsv",
        "figure_id",
        (
            "figure_id", "route_id", "role", "title", "question",
            "data_requirement_ids", "code_module_ids", "source_table",
            "acceptance_test", "required", "status", "notes",
        ),
        (
            "route_id", "role", "title", "question",
            "data_requirement_ids", "code_module_ids", "source_table",
            "acceptance_test", "required", "status",
        ),
        {
            "role": ("main", "supplement", "table"),
            "required": ("true", "false"),
            "status": (
                "idea", "mapped", "spike_generated", "generated",
                "verified", "blocked", "dropped",
            ),
        },
        (
            {"field": "route_id", "target": "routes"},
            {
                "field": "data_requirement_ids",
                "target": "data_requirements",
                "multi": True,
            },
            {
                "field": "code_module_ids",
                "target": "code_requirements",
                "multi": True,
            },
        ),
    ),
    "literature": _spec(
        "evidence/literature.tsv",
        "literature_id",
        (
            "literature_id", "route_id", "query", "nearest_paper", "uri",
            "checked_at", "overlap_level", "disease_overlap",
            "object_overlap", "outcome_overlap", "data_overlap",
            "analysis_overlap", "claim_overlap", "decision", "distinction",
            "notes",
        ),
        (
            "route_id", "query", "nearest_paper", "uri", "checked_at",
            "overlap_level", "disease_overlap", "object_overlap",
            "outcome_overlap", "data_overlap", "analysis_overlap",
            "claim_overlap", "decision", "distinction",
        ),
        {
            "overlap_level": ("clear", "adjacent", "high", "duplicate"),
            "disease_overlap": ("none", "partial", "same", "unknown"),
            "object_overlap": ("none", "partial", "same", "unknown"),
            "outcome_overlap": ("none", "partial", "same", "unknown"),
            "data_overlap": ("none", "partial", "same", "unknown"),
            "analysis_overlap": ("none", "partial", "same", "unknown"),
            "claim_overlap": ("none", "partial", "same", "unknown"),
            "decision": ("continue", "distinguish", "stop"),
        },
        ({"field": "route_id", "target": "routes"},),
    ),
    "decisions": _spec(
        "evidence/decisions.tsv",
        "decision_id",
        (
            "decision_id", "route_id", "stage", "decision", "authority",
            "reviewer", "rationale", "decided_at",
        ),
        (
            "route_id", "stage", "decision", "authority", "reviewer",
            "rationale", "decided_at",
        ),
        {
            "stage": STAGES,
            "decision": (
                "approve_route_evidence", "promote", "retain_training",
                "close_pilot", "backup", "reject", "continue", "refine",
                "reroute", "stop", "approve_release",
            ),
            "authority": ("user", "ai", "collaborator"),
        },
        ({"field": "route_id", "target": "routes"},),
    ),
    "runs": _spec(
        "execution/runs.tsv",
        "run_id",
        (
            "run_id", "route_id", "run_kind", "status", "command", "commit",
            "environment", "input_provenance", "data_manifest",
            "module_release_ids", "started_at", "finished_at", "exit_code",
            "log", "artifacts", "notes",
        ),
        (
            "route_id", "run_kind", "status", "command", "commit",
            "environment", "input_provenance", "started_at", "log",
            "artifacts",
        ),
        {
            "run_kind": (
                "real_data", "unit", "synthetic_integration",
                "environment_smoke",
            ),
            "status": ("planned", "running", "passed", "failed", "invalidated"),
        },
        ({"field": "route_id", "target": "routes"},),
    ),
    "results": _spec(
        "execution/results.tsv",
        "result_id",
        (
            "result_id", "route_id", "run_id", "figure_id", "status",
            "summary", "claim_effect", "next_action", "limitation",
            "source_table",
        ),
        (
            "route_id", "run_id", "figure_id", "status", "summary",
            "claim_effect", "next_action", "limitation", "source_table",
        ),
        {
            "status": ("provisional", "verified", "invalidated"),
            "claim_effect": ("supports", "weakens", "contradicts", "inconclusive"),
            "next_action": ("continue", "refine", "reroute", "stop"),
        },
        (
            {"field": "route_id", "target": "routes"},
            {"field": "run_id", "target": "runs"},
            {"field": "figure_id", "target": "figures"},
        ),
    ),
    "issues": _spec(
        "evidence/issues.tsv",
        "issue_id",
        (
            "issue_id", "route_id", "observed_layer", "candidate_scope",
            "issue_type", "stage", "severity", "observation", "evidence",
            "consequence", "proposed_action", "disposition", "status",
            "blocking", "affected_capability_ids", "promotion_id", "notes",
        ),
        (
            "observed_layer", "candidate_scope", "issue_type", "stage",
            "severity", "observation", "evidence", "consequence",
            "proposed_action", "disposition", "status", "blocking",
        ),
        {
            "observed_layer": ("anchor", "pilot_execution", "both"),
            "candidate_scope": ("pilot", "module", "core"),
            "issue_type": (
                "data_access", "data_identity", "code", "contract",
                "statistical", "scientific_result", "license", "usability",
            ),
            "stage": STAGES,
            "severity": ("low", "medium", "high", "critical"),
            "disposition": (
                "instance_only", "promote_to_module", "promote_to_core",
                "already_covered", "defer", "not_applicable",
            ),
            "status": (
                "open", "resolved", "accepted_risk", "superseded", "wont_fix",
            ),
            "blocking": ("true", "false"),
        },
        ({"field": "route_id", "target": "routes", "optional": True},),
    ),
}
