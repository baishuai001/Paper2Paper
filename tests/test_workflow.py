from __future__ import annotations

import csv
import contextlib
import io
import json
import tempfile
import unittest
from datetime import timedelta
from pathlib import Path

from paper2paper.cli import _require_manuscript_target, main as cli_main
from paper2paper.workspace import (
    _parse_iso_timestamp,
    init_workspace,
    next_actions,
    promote_workspace,
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


def update_manifest(project: Path, **changes) -> None:
    path = project / "PROJECT.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest.update(changes)
    path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def complete_anchor(project: Path) -> None:
    text = """# Anchor audit

## Scientific grammar

The anchor defines a patient-level comparison, a measurable primary outcome,
and a falsifiable associational claim. Discovery, model development and
validation roles are recorded separately. Cells and repeated samples are not
treated as independent patients.

## Figure-to-evidence map

Every main figure is mapped to input data, required metadata, a code module, a
source table, a statistical unit and a bounded claim. Private inputs and wet-lab
evidence remain outside the public reproduction ceiling.

## Module disposition and reproduction boundary

Preprocessing is retained, the marker is substituted, patient-level inference
is required and decorative analyses may be dropped. Missing author code and
undocumented choices are recorded rather than silently invented.
"""
    (project / "anchor/audit.md").write_text(text, encoding="utf-8")


def complete_pilot_outcome(project: Path) -> None:
    text = """# Pilot outcome

## Paper-side outcome

The minimal real-data run answered the bounded execution question and retained
the route's scientific limitations. It did not establish publication value or
turn the Pilot into a manuscript. The data roles, code contract, result and
claim ceiling can be reviewed independently.

## Product-side outcome

The instance-specific implementation was repaired locally. Any proposal for a
module or core change is recorded separately and needs its own promotion
evidence. A negative scientific result would remain a scientific result rather
than being relabeled as a workflow failure.

## Capability and regression contribution

This Pilot exercises only its declared patient-level capability and one real
dataset. It does not validate unrelated modalities, anchors, diseases or the
whole workflow. A second independent case would be needed for transfer claims.

## Human decision and next boundary

The user may promote a manuscript candidate, retain a training case or close
the Pilot. Formal analysis and writing begin only in a separate manuscript
workspace created after explicit promotion.
"""
    (project / "reports/pilot-outcome.md").write_text(text, encoding="utf-8")


def make_executable_route(project: Path) -> None:
    complete_anchor(project)
    update_manifest(project, stage="verification")
    (project / "config").mkdir(exist_ok=True)
    (project / "outputs").mkdir(exist_ok=True)
    (project / "config/data.tsv").write_text(
        "resource\tversion\nexample\tfixture\n", encoding="utf-8"
    )
    (project / "outputs/run.log").write_text("exit_code=0\n", encoding="utf-8")
    (project / "outputs/fig1.tsv").write_text(
        "estimate\tse\n0.2\t0.1\n", encoding="utf-8"
    )
    (project / "outputs/fig1.png").write_bytes(b"fixture-png")
    evidence = project / "evidence"
    add_row(
        evidence / "routes.tsv",
        {
            "route_id": "ROUTE-1",
            "decision_status": "active",
            "evidence_stage": "minimal_real_run",
            "route_role": "manuscript_candidate",
            "mode": "marker",
            "capability_ids": "CAP-TEST-PATIENT-MARKER",
            "title": "Marker substitution",
            "question": "Is marker B associated with outcome Y at patient level?",
            "target_disease": "cancer B",
            "target_object": "marker B",
            "biological_unit": "patient",
            "comparison": "marker-high versus marker-low",
            "primary_outcome": "outcome Y",
            "claim_ceiling": "retrospective patient-level association",
            "falsifier": "no reproducible patient-level association",
            "anchor_reuse": "question;figure order;methods",
            "changed_axes": "marker_gene",
            "scientific_review": "passed",
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
                "searched_at": "2026-08-07",
                "result_count": "1",
                "outcome": "continue",
            },
        )
    add_row(
        evidence / "data_requirements.tsv",
        {
            "requirement_id": "DATAREQ-1",
            "route_id": "ROUTE-1",
            "role": "minimal-run cohort",
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
            "checked_at": "2026-08-07",
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
            "source_version": "2026-08-07 snapshot",
            "checked_at": "2026-08-07",
            "verification": "sample_parsed",
            "local_name": "data.tsv",
            "fields_supplied": "patient_id;outcome;expression",
            "identifier_field": "patient_id",
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
            "role": "reproduction",
            "outcome_used": "true",
            "features_influenced": "false",
            "parameters_influenced": "false",
            "cutoff_influenced": "false",
            "claimed_external_validation": "false",
            "acceptable": "true",
        },
    )
    add_row(
        evidence / "code_requirements.tsv",
        {
            "module_id": "MODULE-1",
            "route_id": "ROUTE-1",
            "capability_id": "CAP-TEST-PATIENT-MARKER",
            "name": "patient-level model",
            "purpose": "primary analysis",
            "input_contract": "one row per patient",
            "output_contract": "source table with estimate and uncertainty",
            "required_tests": "unit;smoke;patient-independence",
            "required": "true",
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
            "path_test": "path with spaces and Unicode",
            "verification": "smoke_passed",
            "decision": "use",
            "checked_at": "2026-08-07",
            "smoke_input": "fixtures/patients.tsv",
            "smoke_output": "outputs/smoke/source-table.tsv",
            "tests_passed": "unit;smoke;patient-independence",
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
            "acceptance_test": "one row per patient; estimate and uncertainty",
            "required": "true",
            "status": "spike_generated",
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
            "checked_at": "2026-08-07",
            "overlap_level": "adjacent",
            "disease_overlap": "same",
            "object_overlap": "partial",
            "outcome_overlap": "same",
            "data_overlap": "none",
            "analysis_overlap": "partial",
            "claim_overlap": "partial",
            "decision": "distinguish",
            "distinction": "Different marker and independent evidence chain.",
        },
    )
    add_row(
        project / "execution/runs.tsv",
        {
            "run_id": "RUN-1",
            "route_id": "ROUTE-1",
            "run_kind": "real_data",
            "status": "passed",
            "command": "Rscript run.R --manifest config/data.tsv",
            "commit": "test-commit",
            "environment": "R 4.4",
            "input_provenance": "Repository test fixture described by config/data.tsv.",
            "data_manifest": "config/data.tsv",
            "started_at": "2026-08-07T10:00:00+08:00",
            "finished_at": "2026-08-07T10:01:00+08:00",
            "exit_code": "0",
            "log": "outputs/run.log",
            "artifacts": "outputs/fig1.tsv;outputs/fig1.png",
        },
    )
    add_row(
        project / "execution/results.tsv",
        {
            "result_id": "RESULT-1",
            "route_id": "ROUTE-1",
            "run_id": "RUN-1",
            "figure_id": "FIG-1",
            "status": "provisional",
            "summary": "The bounded patient-level run generated the source table.",
            "claim_effect": "supports",
            "next_action": "continue",
            "limitation": "One retrospective example cohort.",
            "source_table": "outputs/fig1.tsv",
        },
    )


def add_user_evidence_approval(project: Path) -> None:
    add_row(
        project / "evidence/decisions.tsv",
        {
            "decision_id": "DEC-APPROVE-EVIDENCE",
            "route_id": "ROUTE-1",
            "stage": "pilot_review",
            "decision": "approve_route_evidence",
            "authority": "user",
            "reviewer": "project owner",
            "rationale": "The recorded evidence is sufficient for a promotion decision.",
            "decided_at": "2026-08-07T10:02:00+08:00",
        },
    )


def mark_pilot_for_promotion(project: Path) -> None:
    complete_pilot_outcome(project)
    add_user_evidence_approval(project)
    add_row(
        project / "evidence/decisions.tsv",
        {
            "decision_id": "DEC-PROMOTE",
            "route_id": "ROUTE-1",
            "stage": "pilot_review",
            "decision": "promote",
            "authority": "user",
            "reviewer": "project owner",
            "rationale": "Create a separate formal manuscript project.",
            "decided_at": "2026-08-07T10:03:00+08:00",
        },
    )
    rewrite_rows(
        project / "evidence/routes.tsv",
        lambda rows: rows[0].update({"decision_status": "promoted"}),
    )
    update_manifest(project, stage="pilot_review", selected_route_id="ROUTE-1")


class Paper2PaperWorkflowTests(unittest.TestCase):
    def test_iso_timestamp_parser_is_stable_on_python_310(self) -> None:
        parsed = _parse_iso_timestamp("2026-08-07T19:43:52.6379127Z")
        self.assertIsNotNone(parsed)
        assert parsed is not None
        self.assertEqual(parsed.microsecond, 637912)
        self.assertEqual(parsed.utcoffset(), timedelta(0))

    def test_init_creates_pilot_not_manuscript(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor", "10.test/x")
            report = validate_workspace(project)
            self.assertTrue(report.ok, report.errors)
            manifest = json.loads((project / "PROJECT.json").read_text())
            self.assertEqual(manifest["workspace_kind"], "pilot")
            self.assertTrue((project / "reports/pilot-outcome.md").exists())
            self.assertFalse((project / "manuscript").exists())

    def test_executable_route_requires_a_recorded_real_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(project / "execution/runs.tsv", lambda rows: rows.clear())
            rewrite_rows(project / "execution/results.tsv", lambda rows: rows.clear())
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("no passed real-data run" in e for e in report.errors))

    def test_passed_run_cannot_point_to_missing_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "execution/runs.tsv",
                lambda rows: rows[0].update(
                    {"log": "outputs/missing.log", "artifacts": "outputs/missing.tsv"}
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("missing or unsafe" in error for error in report.errors)
            )
            readiness = route_readiness(project)[0]
            self.assertFalse(readiness["execution_ready"])
            self.assertTrue(any("complete evidence" in gap for gap in readiness["execution_gaps"]))

    def test_unit_run_cannot_satisfy_real_data_execution(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "execution/runs.tsv",
                lambda rows: rows[0].update(
                    {
                        "run_kind": "unit",
                        "input_provenance": "Embedded synthetic fixtures only.",
                        "data_manifest": "",
                    }
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            readiness = route_readiness(project)[0]
            self.assertFalse(readiness["execution_ready"])
            self.assertTrue(
                any("no passed real-data run" in gap for gap in readiness["execution_gaps"])
            )

    def test_readiness_fails_when_passed_run_metadata_is_invalid(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "execution/runs.tsv",
                lambda rows: rows[0].update({"started_at": "not-a-timestamp"}),
            )
            self.assertFalse(validate_workspace(project).ok)
            readiness = route_readiness(project)[0]
            self.assertFalse(readiness["execution_ready"])
            self.assertTrue(
                any("timestamps" in gap for gap in readiness["execution_gaps"])
            )
            rewrite_rows(
                project / "execution/runs.tsv",
                lambda rows: rows[0].update(
                    {
                        "started_at": "2026-08-07T10:00:00",
                        "finished_at": "2026-08-07T10:01:00",
                    }
                ),
            )
            report = validate_workspace(project)
            self.assertTrue(
                any("timezone offsets" in error for error in report.errors)
            )
            self.assertFalse(route_readiness(project)[0]["execution_ready"])

    def test_directories_or_empty_files_cannot_masquerade_as_evidence(self) -> None:
        for mode in ("directories", "empty_files"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as temp_dir:
                project = Path(temp_dir) / "pilot"
                init_workspace(project, "P2P-TEST", "Test", "Anchor")
                make_executable_route(project)
                if mode == "directories":
                    rewrite_rows(
                        project / "execution/runs.tsv",
                        lambda rows: rows[0].update(
                            {
                                "data_manifest": "config",
                                "log": "outputs",
                                "artifacts": "outputs",
                            }
                        ),
                    )
                    rewrite_rows(
                        project / "execution/results.tsv",
                        lambda rows: rows[0].update({"source_table": "outputs"}),
                    )
                    rewrite_rows(
                        project / "evidence/figures.tsv",
                        lambda rows: rows[0].update({"source_table": "outputs"}),
                    )
                else:
                    for relative in (
                        "config/data.tsv",
                        "outputs/run.log",
                        "outputs/fig1.tsv",
                        "outputs/fig1.png",
                    ):
                        (project / relative).write_bytes(b"")
                report = validate_workspace(project)
                self.assertFalse(report.ok)
                self.assertFalse(route_readiness(project)[0]["execution_ready"])

    def test_execution_and_promotion_evidence_are_separate(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            report = validate_workspace(project)
            self.assertTrue(report.ok, report.errors)
            readiness = route_readiness(project)[0]
            self.assertTrue(readiness["execution_ready"])
            self.assertFalse(readiness["promotion_evidence_complete"])
            self.assertTrue(any("user approval" in x for x in readiness["promotion_gaps"]))
            add_user_evidence_approval(project)
            self.assertTrue(route_readiness(project)[0]["promotion_evidence_complete"])

    def test_training_route_can_run_but_cannot_be_promoted(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "evidence/routes.tsv",
                lambda rows: rows[0].update({"route_role": "training"}),
            )
            readiness = route_readiness(project)[0]
            self.assertTrue(readiness["execution_ready"])
            self.assertFalse(readiness["promotion_evidence_complete"])
            self.assertTrue(any("training" in x for x in readiness["promotion_gaps"]))

    def test_promotion_creates_a_separate_manuscript_workspace(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            pilot = Path(temp_dir) / "pilot"
            manuscript = Path(temp_dir) / "manuscript"
            init_workspace(pilot, "P2P-TEST", "Test", "Anchor")
            make_executable_route(pilot)
            mark_pilot_for_promotion(pilot)
            promote_workspace(
                pilot, "ROUTE-1", manuscript, "P2P-MS-TEST", "Formal project"
            )
            report = validate_workspace(manuscript)
            self.assertTrue(report.ok, report.errors)
            manifest = json.loads((manuscript / "PROJECT.json").read_text())
            self.assertEqual(manifest["workspace_kind"], "manuscript_project")
            self.assertEqual(manifest["provenance"]["source_pilot_id"], "P2P-TEST")
            self.assertEqual(manifest["provenance"]["source_run_ids"], "RUN-1")
            self.assertTrue((manuscript / "manuscript/draft.md").exists())
            self.assertEqual((manuscript / "execution/runs.tsv").read_text().count("\n"), 1)
            self.assertNotIn(
                "approve_release",
                (manuscript / "evidence/decisions.tsv").read_text(encoding="utf-8"),
            )
            self.assertFalse((pilot / "manuscript").exists())

    def test_pilot_cannot_preapprove_manuscript_release(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            pilot = Path(temp_dir) / "pilot"
            init_workspace(pilot, "P2P-TEST", "Test", "Anchor")
            make_executable_route(pilot)
            add_row(
                pilot / "evidence/decisions.tsv",
                {
                    "decision_id": "DEC-EARLY-RELEASE",
                    "route_id": "ROUTE-1",
                    "stage": "pilot_review",
                    "decision": "approve_release",
                    "authority": "user",
                    "reviewer": "project owner",
                    "rationale": "Adversarial early approval that must be rejected.",
                    "decided_at": "2026-08-07T09:00:00+08:00",
                },
            )
            report = validate_workspace(pilot)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("pilot cannot record approve_release" in e for e in report.errors)
            )

    def test_promotion_target_must_be_new_registered_child_with_unique_id(self) -> None:
        allowed = ROOT / "manuscript-projects" / "future-test-project"
        self.assertFalse(allowed.exists())
        _require_manuscript_target(ROOT, allowed, "P2P-MS-FUTURE-UNIQUE")

        with self.assertRaisesRegex(ValueError, "one direct child"):
            _require_manuscript_target(
                ROOT,
                ROOT / "outside-manuscript-projects",
                "P2P-MS-OUTSIDE",
            )
        with self.assertRaisesRegex(ValueError, "one direct child"):
            _require_manuscript_target(
                ROOT,
                ROOT / "manuscript-projects" / "nested" / "project",
                "P2P-MS-NESTED",
            )

        existing_id = json.loads(
            (ROOT / "pilots/gastric-nrrs/PROJECT.json").read_text(encoding="utf-8")
        )["project_id"]
        with self.assertRaisesRegex(ValueError, "project_id already exists"):
            _require_manuscript_target(
                ROOT,
                ROOT / "manuscript-projects" / "duplicate-id-test",
                existing_id,
            )
        for invalid_id in ("", "bad id", "../BAD"):
            with self.subTest(project_id=invalid_id):
                with self.assertRaisesRegex(ValueError, "must start"):
                    _require_manuscript_target(
                        ROOT,
                        ROOT / "manuscript-projects" / "invalid-id-test",
                        invalid_id,
                    )

    def test_promotion_decisions_must_follow_run_and_evidence_review(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            pilot = Path(temp_dir) / "pilot"
            manuscript = Path(temp_dir) / "manuscript"
            init_workspace(pilot, "P2P-TEST", "Test", "Anchor")
            make_executable_route(pilot)
            mark_pilot_for_promotion(pilot)

            rewrite_rows(
                pilot / "evidence/decisions.tsv",
                lambda rows: [
                    row.update({"decided_at": "2026-08-06T09:00:00+08:00"})
                    for row in rows
                    if row.get("decision") in {
                        "approve_route_evidence", "promote"
                    }
                ],
            )
            report = validate_workspace(pilot)
            self.assertFalse(report.ok)
            self.assertFalse(
                route_readiness(pilot)[0]["promotion_evidence_complete"]
            )
            self.assertTrue(
                any("approve_route_evidence must occur" in e for e in report.errors)
            )
            with self.assertRaisesRegex(ValueError, "source pilot must pass"):
                promote_workspace(
                    pilot,
                    "ROUTE-1",
                    manuscript,
                    "P2P-MS-EARLY",
                    "Early manuscript",
                )
            self.assertFalse(manuscript.exists())

            rewrite_rows(
                pilot / "evidence/decisions.tsv",
                lambda rows: [
                    row.update(
                        {
                            "decided_at": (
                                "2026-08-07T10:03:00+08:00"
                                if row.get("decision") == "approve_route_evidence"
                                else "2026-08-07T10:02:00+08:00"
                            )
                        }
                    )
                    for row in rows
                    if row.get("decision") in {
                        "approve_route_evidence", "promote"
                    }
                ],
            )
            report = validate_workspace(pilot)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("promote decision must occur" in e for e in report.errors)
            )

    def test_cli_rejects_unregistered_external_pilot_for_repo_root(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            pilot = Path(temp_dir) / "external-pilot"
            target = Path(temp_dir) / "manuscript"
            init_workspace(pilot, "P2P-EXTERNAL", "External", "Anchor")
            make_executable_route(pilot)
            mark_pilot_for_promotion(pilot)
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                exit_code = cli_main(
                    [
                        "promote",
                        str(pilot),
                        "ROUTE-1",
                        str(target),
                        "--project-id",
                        "P2P-MS-EXTERNAL",
                        "--title",
                        "External manuscript",
                        "--repo-root",
                        str(ROOT),
                    ]
                )
            self.assertEqual(exit_code, 2)
            self.assertIn("must be inside", stderr.getvalue())
            self.assertFalse(target.exists())

    def test_incomplete_manuscript_cannot_validate_as_complete(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            pilot = Path(temp_dir) / "pilot"
            manuscript = Path(temp_dir) / "manuscript"
            init_workspace(pilot, "P2P-TEST", "Test", "Anchor")
            make_executable_route(pilot)
            mark_pilot_for_promotion(pilot)
            promote_workspace(
                pilot, "ROUTE-1", manuscript, "P2P-MS-TEST", "Formal project"
            )
            (manuscript / "config").mkdir(exist_ok=True)
            (manuscript / "outputs").mkdir(exist_ok=True)
            (manuscript / "config/data.tsv").write_text(
                "resource\tversion\nexample\tfixture\n", encoding="utf-8"
            )
            (manuscript / "outputs/run.log").write_text(
                "exit_code=0\n", encoding="utf-8"
            )
            (manuscript / "outputs/fig1.tsv").write_text(
                "estimate\tse\n0.2\t0.1\n", encoding="utf-8"
            )
            (manuscript / "outputs/fig1.png").write_bytes(b"fixture-png")
            (manuscript / "analysis/specification.md").write_text(
                "# Frozen specification\n\n" + "Prespecified patient-level analysis. " * 20,
                encoding="utf-8",
            )
            (manuscript / "manuscript/draft.md").write_text(
                "# Draft\n\n" + "Substantive reviewed manuscript text. " * 30,
                encoding="utf-8",
            )
            add_row(
                manuscript / "execution/runs.tsv",
                {
                    "run_id": "MS-RUN-1",
                    "route_id": "ROUTE-1",
                    "run_kind": "real_data",
                    "status": "passed",
                    "command": "Rscript run.R --manifest config/data.tsv",
                    "commit": "test-commit",
                    "environment": "R 4.4",
                    "input_provenance": "Repository test fixture.",
                    "data_manifest": "config/data.tsv",
                    "started_at": "2026-08-07T11:00:00+08:00",
                    "finished_at": "2026-08-07T11:01:00+08:00",
                    "exit_code": "0",
                    "log": "outputs/run.log",
                    "artifacts": "outputs/fig1.tsv;outputs/fig1.png",
                },
            )
            add_row(
                manuscript / "execution/results.tsv",
                {
                    "result_id": "MS-RESULT-1",
                    "route_id": "ROUTE-1",
                    "run_id": "MS-RUN-1",
                    "figure_id": "FIG-1",
                    "status": "verified",
                    "summary": "A result exists, but the promoted code and figure loop were not reverified.",
                    "claim_effect": "supports",
                    "next_action": "continue",
                    "limitation": "Deliberately incomplete completion fixture.",
                    "source_table": "outputs/fig1.tsv",
                },
            )
            add_row(
                manuscript / "evidence/decisions.tsv",
                {
                    "decision_id": "MS-APPROVE-RELEASE",
                    "route_id": "ROUTE-1",
                    "stage": "complete",
                    "decision": "approve_release",
                    "authority": "user",
                    "reviewer": "project owner",
                    "rationale": "Adversarial fixture that must not override evidence gaps.",
                    "decided_at": "2026-08-07T11:02:00+08:00",
                },
            )
            update_manifest(manuscript, stage="complete")
            report = validate_workspace(manuscript)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("complete manuscript project" in error for error in report.errors)
            )
            rewrite_rows(
                manuscript / "evidence/decisions.tsv",
                lambda rows: next(
                    row for row in rows
                    if row.get("decision") == "approve_release"
                ).update({"stage": "anchor_audit"}),
            )
            report = validate_workspace(manuscript)
            self.assertTrue(
                any("recorded at stage=complete" in error for error in report.errors)
            )
            rewrite_rows(
                manuscript / "evidence/decisions.tsv",
                lambda rows: next(
                    row for row in rows
                    if row.get("decision") == "approve_release"
                ).update(
                    {
                        "stage": "complete",
                        "decided_at": "2026-08-07T10:59:00+08:00",
                    }
                ),
            )
            report = validate_workspace(manuscript)
            self.assertTrue(
                any("occur on or after" in error for error in report.errors)
            )
            rewrite_rows(
                manuscript / "evidence/decisions.tsv",
                lambda rows: next(
                    row for row in rows
                    if row.get("decision") == "approve_release"
                ).update({"decided_at": "2026-08-07"}),
            )
            report = validate_workspace(manuscript)
            self.assertTrue(
                any("full timestamp" in error for error in report.errors)
            )

    def test_pilot_cannot_enter_manuscript_stage(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            update_manifest(project, stage="writing")
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("stage" in e for e in report.errors))

    def test_user_can_complete_a_training_pilot_without_a_manuscript(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            complete_pilot_outcome(project)
            rewrite_rows(
                project / "evidence/routes.tsv",
                lambda rows: rows[0].update({"route_role": "training"}),
            )
            add_row(
                project / "evidence/decisions.tsv",
                {
                    "decision_id": "DEC-RETAIN",
                    "route_id": "ROUTE-1",
                    "stage": "pilot_review",
                    "decision": "retain_training",
                    "authority": "user",
                    "reviewer": "project owner",
                    "rationale": "Keep the bounded run as training evidence only.",
                    "decided_at": "2026-08-07",
                },
            )
            update_manifest(project, stage="pilot_complete")
            report = validate_workspace(project)
            self.assertTrue(report.ok, report.errors)
            self.assertFalse((project / "manuscript").exists())

    def test_metadata_only_data_cannot_support_minimal_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "evidence/data_candidates.tsv",
                lambda rows: rows[0].update({"verification": "metadata_checked"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("parsed usable candidate" in e for e in report.errors))

    def test_install_only_or_unlicensed_code_cannot_support_minimal_run(self) -> None:
        for changes in (
            {"verification": "install_passed"},
            {"license": "Repository license not selected; internal use only"},
            {"path_portability": "not_checked", "path_test": "not run"},
            {"tests_passed": "unit;smoke"},
        ):
            with self.subTest(changes=changes), tempfile.TemporaryDirectory() as temp_dir:
                project = Path(temp_dir) / "pilot"
                init_workspace(project, "P2P-TEST", "Test", "Anchor")
                make_executable_route(project)
                rewrite_rows(
                    project / "evidence/code_candidates.tsv",
                    lambda rows, changes=changes: rows[0].update(changes),
                )
                report = validate_workspace(project)
                self.assertFalse(report.ok)
                self.assertTrue(any("smoke-tested donor" in e for e in report.errors))

    def test_data_contract_and_field_provenance_are_enforced(self) -> None:
        for table, changes in (
            ("data_candidates.tsv", {"subjects": "10"}),
            ("data_candidates.tsv", {"fields_checked": "patient_id;expression"}),
            ("data_resources.tsv", {"fields_supplied": "expression"}),
        ):
            with self.subTest(table=table, changes=changes), tempfile.TemporaryDirectory() as temp_dir:
                project = Path(temp_dir) / "pilot"
                init_workspace(project, "P2P-TEST", "Test", "Anchor")
                make_executable_route(project)
                rewrite_rows(
                    project / "evidence" / table,
                    lambda rows, changes=changes: rows[0].update(changes),
                )
                self.assertFalse(validate_workspace(project).ok)

    def test_development_exposed_cohort_cannot_be_external_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
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
                },
            )
            rewrite_rows(
                project / "evidence/cohort_usage.tsv",
                lambda rows: rows[0].update(
                    {"role": "external_validation", "claimed_external_validation": "true"}
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("cannot be claimed" in e for e in report.errors))

    def test_required_signature_needs_locked_computable_specification(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "evidence/routes.tsv",
                lambda rows: rows[0].update({"model_spec_required": "true"}),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("computable model specification" in e for e in report.errors))

    def test_open_blocking_issue_prevents_readiness(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            add_row(
                project / "evidence/issues.tsv",
                {
                    "issue_id": "ISSUE-1",
                    "route_id": "ROUTE-1",
                    "observed_layer": "pilot_execution",
                    "candidate_scope": "pilot",
                    "issue_type": "data_access",
                    "stage": "verification",
                    "severity": "critical",
                    "observation": "Required endpoint cannot be parsed.",
                    "evidence": "smoke log",
                    "consequence": "Primary result cannot be generated.",
                    "proposed_action": "Repair parser.",
                    "disposition": "instance_only",
                    "status": "open",
                    "blocking": "true",
                },
            )
            self.assertFalse(validate_workspace(project).ok)
            rewrite_rows(
                project / "evidence/issues.tsv",
                lambda rows: rows[0].update({"status": "resolved"}),
            )
            self.assertTrue(validate_workspace(project).ok)

    def test_promotion_disposition_requires_scoped_link(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            add_row(
                project / "evidence/issues.tsv",
                {
                    "issue_id": "ISSUE-1",
                    "route_id": "ROUTE-1",
                    "observed_layer": "pilot_execution",
                    "candidate_scope": "core",
                    "issue_type": "data_identity",
                    "stage": "verification",
                    "severity": "high",
                    "observation": "Identity must fail closed.",
                    "evidence": "test log",
                    "consequence": "Rows may be mismatched.",
                    "proposed_action": "Promote a guarded identity rule.",
                    "disposition": "promote_to_core",
                    "status": "resolved",
                    "blocking": "false",
                },
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("promotion_id" in e for e in report.errors))

    def test_each_route_requires_data_code_and_literature_searches(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "evidence/search_log.tsv",
                lambda rows: rows.__setitem__(
                    slice(None), [row for row in rows if row["domain"] != "code"]
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("no recorded code search" in e for e in report.errors))

    def test_duplicate_publication_blocks_route(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            rewrite_rows(
                project / "evidence/literature.tsv",
                lambda rows: rows[0].update(
                    {"overlap_level": "duplicate", "decision": "stop"}
                ),
            )
            report = validate_workspace(project)
            self.assertFalse(report.ok)
            self.assertTrue(any("blocking duplicate" in e for e in report.errors))

    def test_report_uses_v2_evidence_language(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            path = write_readiness_report(project)
            text = path.read_text(encoding="utf-8")
            self.assertIn("Execution ready: `true`", text)
            self.assertIn("Promotion evidence complete: `false`", text)
            self.assertIn("Evidence stage: `minimal_real_run`", text)
            self.assertNotIn("Manuscript eligible", text)

    def test_next_action_never_starts_manuscript_inside_pilot(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            project = Path(temp_dir) / "pilot"
            init_workspace(project, "P2P-TEST", "Test", "Anchor")
            make_executable_route(project)
            complete_pilot_outcome(project)
            actions = next_actions(project)
            self.assertTrue(any("do not start a manuscript draft" in action for action in actions))

    def test_repository_pilots_validate(self) -> None:
        for project in (
            ROOT / "pilots/spp1-tam-jitc",
            ROOT / "pilots/gastric-nrrs",
        ):
            report = validate_workspace(project)
            self.assertTrue(report.ok, f"{project}: {report.errors}")


if __name__ == "__main__":
    unittest.main()
