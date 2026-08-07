from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from paper2paper.workspace import (
    init_workspace,
    next_actions,
    route_readiness,
    validate_workspace,
    write_readiness_report,
)


ROOT = Path(__file__).resolve().parents[1]


def add_row(path: Path, values: dict[str, str]) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = csv.DictReader(handle, delimiter="\t").fieldnames or []
    unknown = set(values) - set(header)
    if unknown:
        raise ValueError(f"unknown test fields: {sorted(unknown)}")
    row = {field: "" for field in header}
    row.update(values)
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writerow(row)


def rewrite_rows(path: Path, change) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(row) for row in reader]
    change(rows)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def complete_anchor(project: Path) -> None:
    text = """# Anchor audit

## Scientific grammar

The anchor defines a patient-level comparison, a measurable primary outcome,
and a falsifiable associational claim. The discovery and validation roles are
separated. Cells and repeated samples are not treated as independent patients.

## Figure-to-evidence map

Each main figure is mapped to an input dataset, required metadata, an analysis
module, a source table, and a manuscript claim. Private inputs and wet-lab
evidence are explicitly outside the public reproduction ceiling.

## Module disposition

Preprocessing is retained, the marker is substituted, statistical inference
is repaired to remain patient-level, and decorative analyses may be dropped.
"""
    (project / "anchor/audit.md").write_text(text, encoding="utf-8")


def make_ready_route(project: Path) -> None:
    complete_anchor(project)
    manifest_path = project / "PROJECT.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["stage"] = "route_generation"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    evidence = project / "evidence"
    add_row(
        evidence / "routes.tsv",
        {
            "route_id": "ROUTE-1",
            "status": "ready",
            "route_role": "manuscript_candidate",
            "mode": "marker",
            "title": "Marker substitution",
            "question": "Is marker B associated with outcome Y at patient level?",
            "target_disease": "cancer B",
            "target_object": "marker B",
            "biological_unit": "patient",
            "comparison": "marker-high versus marker-low",
            "primary_outcome": "outcome Y",
            "claim_ceiling": "associational",
            "falsifier": "no reproducible patient-level association",
            "anchor_reuse": "question;figure order;methods",
            "changed_axes": "marker_gene",
            "science_status": "pass",
            "science_basis": "Patient-level design with an explicit comparison.",
            "minimum_main_figures": "1",
            "data_burden": "low",
            "code_burden": "low",
            "beginner_burden": "low",
            "estimated_calendar_time": "two weeks for the minimum route",
            "model_spec_required": "false",
            "main_risk": "marker instability",
            "stop_condition": "required data or code cannot pass the spike",
        },
    )
    for index, domain in enumerate(("data", "code", "literature"), start=1):
        add_row(
            evidence / "search_log.tsv",
            {
                "search_id": f"SEARCH-{index}",
                "route_id": "ROUTE-1",
                "domain": domain,
                "source": "documented source",
                "query": f"route-specific {domain} query",
                "searched_at": "2026-08-03",
                "result_count": "1",
                "outcome": "continue",
                "notes": "fixture",
            },
        )
    add_row(
        evidence / "data_requirements.tsv",
        {
            "requirement_id": "DATAREQ-1",
            "route_id": "ROUTE-1",
            "role": "discovery",
            "disease": "cancer B",
            "tissue": "tumor",
            "modality": "bulk RNA-seq",
            "biological_unit": "patient",
            "comparison_or_outcome": "outcome Y",
            "required_fields": "patient_id;outcome;expression",
            "minimum_subjects": "50",
            "access_limit": "public or locally authorized",
            "independence_required": "false",
            "required": "true",
            "notes": "fixture",
        },
    )
    add_row(
        evidence / "data_candidates.tsv",
        {
            "data_id": "DATA-1",
            "route_id": "ROUTE-1",
            "requirement_id": "DATAREQ-1",
            "name": "Example cohort",
            "accession": "EXAMPLE-1",
            "uri": "https://example.org/data",
            "source": "public repository",
            "disease": "cancer B",
            "tissue": "tumor",
            "modality": "bulk RNA-seq",
            "subjects": "60",
            "biological_unit": "patient",
            "fields_checked": "patient_id;outcome;expression",
            "contract_match": "true",
            "access": "public",
            "verification": "sample_parsed",
            "independence": "not_applicable",
            "decision": "use",
            "checked_at": "2026-08-03",
            "notes": "fixture",
        },
    )
    add_row(
        evidence / "code_requirements.tsv",
        {
            "module_id": "MODULE-1",
            "route_id": "ROUTE-1",
            "name": "patient-level model",
            "purpose": "primary analysis",
            "input_contract": "one row per patient",
            "output_contract": "source table with estimate and uncertainty",
            "required_tests": "unit;smoke;patient-independence",
            "required": "true",
            "notes": "fixture",
        },
    )
    add_row(
        evidence / "data_resources.tsv",
        {
            "resource_id": "RESOURCE-1",
            "data_id": "DATA-1",
            "name": "Example parsed file",
            "role": "expression and outcome",
            "uri": "https://example.org/data.tsv",
            "source_version": "2026-08-03 snapshot",
            "checked_at": "2026-08-03",
            "verification": "sample_parsed",
            "local_name": "data.tsv",
            "fields_supplied": "patient_id;outcome;expression",
            "identifier_field": "patient_id",
            "notes": "fixture",
        },
    )
    add_row(
        evidence / "cohort_usage.tsv",
        {
            "usage_id": "USAGE-1",
            "route_id": "ROUTE-1",
            "data_id": "DATA-1",
            "cohort_key": "EXAMPLE-COHORT",
            "analysis_step": "primary association",
            "role": "discovery",
            "outcome_used": "true",
            "features_influenced": "false",
            "parameters_influenced": "false",
            "cutoff_influenced": "false",
            "claimed_external_validation": "false",
            "acceptable": "true",
            "notes": "fixture",
        },
    )
    add_row(
        evidence / "code_candidates.tsv",
        {
            "code_id": "CODE-1",
            "route_id": "ROUTE-1",
            "module_id": "MODULE-1",
            "name": "Official implementation",
            "source_type": "official",
            "uri": "https://example.org/code",
            "version": "1.0.0",
            "license": "MIT",
            "language": "R",
            "environment": "R 4.4",
            "entrypoint": "Rscript run.R",
            "contract_match": "true",
            "noninteractive": "true",
            "private_inputs": "false",
            "hardcoded_paths": "false",
            "path_portability": "target_environment_passed",
            "path_test": "Windows path with spaces",
            "verification": "smoke_passed",
            "decision": "use",
            "checked_at": "2026-08-03",
            "smoke_input": "fixtures/representative-patients.tsv",
            "smoke_output": "outputs/smoke/model-source-table.tsv",
            "tests_passed": "unit;smoke;patient-independence",
            "notes": "fixture",
        },
    )
    add_row(
        evidence / "figures.tsv",
        {
            "figure_id": "FIG-1",
            "route_id": "ROUTE-1",
            "role": "main",
            "title": "Primary association",
            "question": "Does marker B associate with outcome Y?",
            "data_requirement_ids": "DATAREQ-1",
            "code_module_ids": "MODULE-1",
            "source_table": "outputs/fig1.tsv",
            "acceptance_test": "one row per patient; prespecified direction",
            "required": "true",
            "status": "spike_generated",
            "notes": "fixture",
        },
    )
    add_row(
        evidence / "literature.tsv",
        {
            "literature_id": "LIT-1",
            "route_id": "ROUTE-1",
            "query": "marker B cancer B outcome Y",
            "nearest_paper": "Nearest relevant paper",
            "uri": "https://example.org/paper",
            "checked_at": "2026-08-03",
            "overlap_level": "adjacent",
            "disease_overlap": "same",
            "object_overlap": "partial",
            "outcome_overlap": "same",
            "data_overlap": "none",
            "analysis_overlap": "partial",
            "claim_overlap": "partial",
            "decision": "distinguish",
            "distinction": "Different marker and independent evidence chain.",
            "notes": "fixture",
        },
    )


def select_route(project: Path, record_decision: bool = True) -> None:
    rewrite_rows(
        project / "evidence/routes.tsv",
        lambda rows: rows[0].update({"status": "selected"}),
    )
    manifest_path = project / "PROJECT.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["stage"] = "selection"
    manifest["selected_route_id"] = "ROUTE-1"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if record_decision:
        add_row(
            project / "evidence/decisions.tsv",
            {
                "decision_id": "DECISION-1",
                "route_id": "ROUTE-1",
                "stage": "selection",
                "decision": "select",
                "reviewer": "project owner",
                "rationale": "All execution evidence passed the minimum spike.",
                "decided_at": "2026-08-03",
            },
        )


class Paper2PaperWorkflowTests(unittest.TestCase):
    def test_init_creates_a_valid_intake_workspace(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor", "10.test/x")
            report = validate_workspace(project)
            self.assertTrue(report.ok, report.errors)
            actions = next_actions(project)
            self.assertTrue(any("anchor/audit.md" in item for item in actions))
            self.assertTrue(any("route portfolio" in item for item in actions))

    def test_fully_verified_route_is_ready(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            report = validate_workspace(project)
            self.assertTrue(report.ok, report.errors)
            readiness = route_readiness(project)
            self.assertTrue(
                readiness[0]["execution_ready"],
                readiness[0]["execution_gaps"],
            )
            self.assertTrue(
                readiness[0]["manuscript_eligible"],
                readiness[0]["manuscript_gaps"],
            )

    def test_training_reproduction_cannot_be_selected_as_manuscript_route(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/routes.tsv",
                lambda rows: rows[0].update({"route_role": "training"}),
            )
            readiness = route_readiness(project)
            self.assertTrue(
                readiness[0]["execution_ready"],
                readiness[0]["execution_gaps"],
            )
            self.assertFalse(readiness[0]["manuscript_eligible"])
            report = validate_workspace(project)
            self.assertTrue(report.ok, report.errors)
            select_route(project)
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("not manuscript-eligible" in e for e in report.errors)
            )

    def test_metadata_only_data_cannot_support_a_ready_route(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/data_candidates.tsv",
                lambda rows: rows[0].update({"verification": "metadata_checked"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("parsed usable candidate" in e for e in report.errors))

    def test_install_only_code_cannot_support_a_ready_route(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/code_candidates.tsv",
                lambda rows: rows[0].update({"verification": "install_passed"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("smoke-tested donor" in e for e in report.errors))

    def test_license_placeholder_cannot_qualify_code(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/code_candidates.tsv",
                lambda rows: rows[0].update(
                    {"license": "Repository license not selected; internal project use only"}
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("qualified smoke-tested donor" in e for e in report.errors))

    def test_data_candidate_requires_field_level_resource_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/data_resources.tsv",
                lambda rows: rows[0].update({"fields_supplied": "expression"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("field-level resource" in e for e in report.errors))

    def test_target_environment_path_smoke_is_required(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/code_candidates.tsv",
                lambda rows: rows[0].update(
                    {"path_portability": "not_checked", "path_test": "not run"}
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("smoke-tested donor" in e for e in report.errors))

    def test_development_cohort_cannot_be_claimed_as_external_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            add_row(
                project / "evidence/cohort_usage.tsv",
                {
                    "usage_id": "USAGE-2",
                    "route_id": "ROUTE-1",
                    "data_id": "DATA-1",
                    "cohort_key": "EXAMPLE-COHORT",
                    "analysis_step": "feature screening",
                    "role": "feature_screening",
                    "outcome_used": "true",
                    "features_influenced": "true",
                    "parameters_influenced": "false",
                    "cutoff_influenced": "false",
                    "claimed_external_validation": "false",
                    "acceptable": "true",
                    "notes": "fixture",
                },
            )
            rewrite_rows(
                project / "evidence/cohort_usage.tsv",
                lambda rows: rows[0].update(
                    {
                        "role": "external_validation",
                        "claimed_external_validation": "true",
                    }
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("cannot be claimed" in e for e in report.errors))

    def test_required_signature_needs_a_locked_computable_specification(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/routes.tsv",
                lambda rows: rows[0].update({"model_spec_required": "true"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("computable model specification" in e for e in report.errors))

    def test_open_blocking_issue_prevents_route_readiness(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            add_row(
                project / "evidence/issues.tsv",
                {
                    "issue_id": "ISSUE-1",
                    "route_id": "ROUTE-1",
                    "scope": "product",
                    "stage": "verification",
                    "severity": "critical",
                    "observation": "Required endpoint cannot be parsed.",
                    "evidence": "smoke log",
                    "consequence": "Primary result cannot be generated.",
                    "proposed_action": "Repair parser.",
                    "disposition": "promote_to_core",
                    "status": "open",
                    "blocking": "true",
                    "notes": "fixture",
                },
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("open blocking issues" in e for e in report.errors))
            rewrite_rows(
                project / "evidence/issues.tsv",
                lambda rows: rows[0].update(
                    {"status": "partially_resolved_in_core"}
                ),
            )
            self.assertFalse(validate_workspace(project).ok)
            rewrite_rows(
                project / "evidence/issues.tsv",
                lambda rows: rows[0].update({"status": "verified_in_core"}),
            )
            self.assertTrue(validate_workspace(project).ok)

    def test_too_few_subjects_cannot_support_a_ready_route(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/data_candidates.tsv",
                lambda rows: rows[0].update({"subjects": "10"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("minimum subject count" in error for error in report.errors)
            )

    def test_required_metadata_contract_is_enforced(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/data_candidates.tsv",
                lambda rows: rows[0].update(
                    {"fields_checked": "patient_id;expression"}
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("minimum subject count" in error for error in report.errors)
            )

    def test_code_contract_and_required_tests_are_enforced(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/code_candidates.tsv",
                lambda rows: rows[0].update({"tests_passed": "unit;smoke"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("smoke-tested donor" in error for error in report.errors)
            )

    def test_each_route_requires_recorded_searches(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/search_log.tsv",
                lambda rows: rows.__setitem__(slice(None), [
                    row for row in rows if row["domain"] != "code"
                ]),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("no recorded code search" in e for e in report.errors))

    def test_duplicate_publication_is_a_hard_stop(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            rewrite_rows(
                project / "evidence/literature.tsv",
                lambda rows: rows[0].update(
                    {"overlap_level": "duplicate", "decision": "stop"}
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("blocking duplicate" in e for e in report.errors))

    def test_selection_requires_a_human_decision(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            select_route(project, record_decision=False)
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("select decision" in e for e in report.errors))

    def test_ready_route_can_be_selected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            select_route(project)
            report = validate_workspace(project)
            self.assertTrue(report.ok, report.errors)

    def test_report_is_regenerated_from_current_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "project"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_ready_route(project)
            path = write_readiness_report(project)
            text = path.read_text(encoding="utf-8")
            self.assertIn("ROUTE-1", text)
            self.assertIn("Execution ready: `true`", text)
            self.assertIn("Manuscript eligible: `true`", text)

    def test_pilot_workspaces_validate(self) -> None:
        for project in (
            ROOT / "pilots/spp1-tam-jitc",
            ROOT / "pilots/gastric-nrrs",
        ):
            report = validate_workspace(project)
            self.assertTrue(report.ok, f"{project}: {report.errors}")


if __name__ == "__main__":
    unittest.main()
