from __future__ import annotations

import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from paper2paper.core import (
    downstream_impact,
    expected_execution_priority,
    init_project,
    load_contract,
    validate_project,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SPP1_PROJECT = REPOSITORY_ROOT / "workspaces" / "spp1-tam-jitc"
NRRS_PROJECT = (
    REPOSITORY_ROOT
    / "workspaces"
    / "bmc-cancer-2025-gastric-nerve-model"
)


def append_tsv(path: Path, values: list[str]) -> None:
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(values)


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv_rows(
    path: Path, fieldnames: list[str], rows: list[dict[str, str]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def update_manifest(target: Path, **changes: object) -> None:
    manifest_path = target / "PROJECT.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(changes)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def make_approved_project(target: Path) -> None:
    init_project(target, "P2P-APPROVED", "Approved fixture")
    update_manifest(
        target,
        project_status="active",
        current_gate="G2_SPIKE",
        anchor_paper_id="PAPER-001",
        active_route_id="ROUTE-001",
    )
    append_tsv(
        target / "registry" / "papers.tsv",
        [
            "PAPER-001",
            "Fixture anchor",
            "10.0000/fixture",
            "anchor",
            "accepted",
            "https://example.org/paper",
            "Synthetic test fixture",
        ],
    )
    append_tsv(
        target / "registry" / "routes.tsv",
        [
            "ROUTE-001",
            "1",
            "",
            "manuscript_candidate",
            "faithful_reproduction",
            "approved",
            "Executable route",
            "Is the prespecified patient-level association present?",
            "test cancer",
            "test cell state",
            "patient outcome",
            "design;figure_order;method",
            "patient-level discovery and validation",
            "a null external-validation effect",
            "associational",
            "REVIEW-001",
        ],
    )
    append_tsv(
        target / "registry" / "route_assessments.tsv",
        [
            "ASSESS-001",
            "ROUTE-001",
            "pass",
            "Patient-level design and bounded claim are prespecified.",
            "verified",
            "A representative input was downloaded and parsed.",
            "qualified",
            "Every required module maps to grade A code.",
            "complete",
            "The minimum figure plan is ready.",
            "low",
            "low",
            "high",
            "clear",
            "A dated search found no substantive duplicate.",
            "P0",
            "Small fixture with complete execution path.",
            "none",
            "primary",
            "approved",
        ],
    )
    append_tsv(
        target / "registry" / "datasets.tsv",
        [
            "DATA-001",
            "Fixture dataset",
            "FIXTURE-001",
            "https://example.org/data",
            "bulk RNA-seq",
            "test cancer",
            "100 patients",
            "patient",
            "counts and clinical table",
            "patient_id;outcome;counts",
            "public",
            "sample_verified",
            "abc123",
            "1 MB",
            "FIXTURE",
            "discovery",
            "accepted",
            "2026-08-02",
            "Synthetic test fixture",
        ],
    )
    append_tsv(
        target / "registry" / "dataset_route_map.tsv",
        [
            "DATAMAP-001",
            "ROUTE-001",
            "DATA-001",
            "discovery and validation fixture",
            "true",
            "patient_id;outcome;counts",
            "verified",
            "not_applicable",
            "Synthetic test fixture",
        ],
    )
    append_tsv(
        target / "registry" / "code_sources.tsv",
        [
            "CODE-001",
            "Fixture pipeline",
            "https://example.org/code",
            "reconstructed",
            "abc123",
            "MIT",
            "Python 3.11",
            "requirements.lock",
            "python -m fixture_pipeline",
            "true",
            "full_pass",
            "false",
            "false",
            "passed",
            "passed",
            "A",
            "accepted",
            "2026-08-02",
            "Synthetic test fixture",
        ],
    )
    append_tsv(
        target / "registry" / "code_module_map.tsv",
        [
            "CODEMAP-001",
            "ROUTE-001",
            "primary_analysis",
            "CODE-001",
            "src/fixture_pipeline.py",
            "none",
            "true",
            "unit;smoke;integration;scientific_invariant",
            "verified",
            "Synthetic test fixture",
        ],
    )
    append_tsv(
        target / "registry" / "publication_overlap.tsv",
        [
            "OVERLAP-001",
            "ROUTE-001",
            "",
            "Nearest searched paper",
            "https://example.org/neighbor",
            "test cancer test cell state patient outcome",
            "2026-08-02",
            "clear",
            "none",
            "none",
            "none",
            "none",
            "none",
            "none",
            "none",
            "No central overlap found in fixture search.",
            "proceed",
            "Synthetic test fixture",
        ],
    )
    append_tsv(
        target / "registry" / "figure_plan.tsv",
        [
            "FIGURE-001",
            "ROUTE-001",
            "main",
            "Figure 1",
            "Patient-level primary effect figure",
            "",
            "DATA-001",
            "primary_analysis",
            "CODEMAP-001",
            "outputs/figure_001_source.tsv",
            "One row per patient and prespecified contrast direction.",
            "true",
            "ready",
            "Synthetic test fixture",
        ],
    )
    append_tsv(
        target / "registry" / "reviews.tsv",
        [
            "REVIEW-001",
            "G2_SPIKE",
            "ROUTE-001",
            "completed",
            "approve",
            "fixture-reviewer",
            "All execution gates pass in the synthetic fixture.",
            "2026-08-02",
        ],
    )


class Paper2PaperCoreTests(unittest.TestCase):
    def test_intake_workspaces_validate(self) -> None:
        for project in (SPP1_PROJECT, NRRS_PROJECT):
            report = validate_project(project)
            self.assertTrue(report.ok, report.errors)

    def test_init_creates_valid_intake_without_legacy_policy(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            init_project(target, "P2P-TEST", "Test project")
            report = validate_project(target)
            self.assertTrue(report.ok, report.errors)
            manifest = json.loads(
                (target / "PROJECT.json").read_text(encoding="utf-8")
            )
            self.assertFalse(manifest["adaptation_policy"]["novelty_required"])
            self.assertNotIn("quality_axes", manifest["publication_goal"])
            self.assertTrue(
                (target / "registry" / "datasets.tsv").exists()
            )
            self.assertTrue(
                (target / "registry" / "code_sources.tsv").exists()
            )
            self.assertTrue(
                (target / "registry" / "figure_plan.tsv").exists()
            )

    def test_legacy_policy_keys_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            init_project(target, "P2P-TEST", "Test project")
            manifest_path = target / "PROJECT.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["publication_goal"]["quality_axes"] = ["novelty"]
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            report = validate_project(target)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("legacy PaperRoute policy keys" in e for e in report.errors)
            )

    def test_novelty_cannot_be_reintroduced_as_a_gate(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            init_project(target, "P2P-TEST", "Test project")
            manifest_path = target / "PROJECT.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["adaptation_policy"]["novelty_required"] = True
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            report = validate_project(target)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("novelty_required must be false" in e for e in report.errors)
            )

    def test_marker_substitution_is_an_allowed_route(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            init_project(target, "P2P-TEST", "Test project")
            update_manifest(
                target,
                project_status="route_review",
                current_gate="G0_SCOPE",
                anchor_paper_id="PAPER-001",
            )
            append_tsv(
                target / "registry" / "papers.tsv",
                [
                    "PAPER-001", "Anchor", "", "anchor", "accepted",
                    "https://example.org/anchor", "Fixture",
                ],
            )
            append_tsv(
                target / "registry" / "routes.tsv",
                [
                    "ROUTE-001", "1", "", "manuscript_candidate",
                    "marker_substitution", "proposed", "Replace marker",
                    "Does marker B reproduce the anchor relation?", "cancer",
                    "marker B cells", "patient outcome", "design;figures",
                    "patient-level association", "null validation effect",
                    "associational", "",
                ],
            )
            append_tsv(
                target / "registry" / "route_adaptations.tsv",
                [
                    "ADAPT-001", "ROUTE-001", "marker_gene", "substitute",
                    "marker A", "marker B", "Biological rationale recorded.",
                    "expression and patient outcome", "anchor workflow",
                    "medium", "proposed",
                ],
            )
            report = validate_project(target)
            self.assertTrue(report.ok, report.errors)

    def test_non_reproduction_route_requires_an_adaptation_map(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            init_project(target, "P2P-TEST", "Test project")
            append_tsv(
                target / "registry" / "routes.tsv",
                [
                    "ROUTE-001", "1", "", "manuscript_candidate",
                    "cancer_type_substitution", "proposed", "Replace cancer",
                    "Does the relation hold in cancer B?", "cancer B", "cell",
                    "outcome", "design", "patient evidence", "null result",
                    "associational", "",
                ],
            )
            report = validate_project(target)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("explicit adaptation" in e for e in report.errors)
            )

    def test_priority_rules_are_transparent(self) -> None:
        baseline = {
            "scientific_validity": "pass",
            "data_readiness": "verified",
            "code_readiness": "qualified",
            "figure_coverage": "complete",
            "publication_overlap": "adjacent",
        }
        self.assertEqual(expected_execution_priority(baseline), "P0")
        self.assertEqual(
            expected_execution_priority({**baseline, "code_readiness": "adaptable"}),
            "P1",
        )
        self.assertEqual(
            expected_execution_priority(
                {
                    **baseline,
                    "data_readiness": "partial",
                    "code_readiness": "rebuild_required",
                }
            ),
            "P2",
        )
        self.assertEqual(
            expected_execution_priority(
                {**baseline, "publication_overlap": "duplicate"}
            ),
            "P3",
        )

    def test_fully_qualified_route_can_be_approved(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            make_approved_project(target)
            report = validate_project(target)
            self.assertTrue(report.ok, report.errors)

    def test_approved_route_requires_sample_verified_data(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            make_approved_project(target)
            path = target / "registry" / "datasets.tsv"
            fields, rows = read_tsv_rows(path)
            rows[0]["download_status"] = "metadata_verified"
            write_tsv_rows(path, fields, rows)
            report = validate_project(target)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("verified sample or download" in e for e in report.errors)
            )

    def test_approved_route_requires_grade_a_or_b_code(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            make_approved_project(target)
            path = target / "registry" / "code_sources.tsv"
            fields, rows = read_tsv_rows(path)
            rows[0]["qualification_grade"] = "C"
            write_tsv_rows(path, fields, rows)
            report = validate_project(target)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("uses unqualified code" in e for e in report.errors)
            )

    def test_grade_a_code_has_strict_engineering_requirements(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            make_approved_project(target)
            path = target / "registry" / "code_sources.tsv"
            fields, rows = read_tsv_rows(path)
            rows[0]["hardcoded_paths"] = "true"
            write_tsv_rows(path, fields, rows)
            report = validate_project(target)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("grade A cannot contain hardcoded paths" in e for e in report.errors)
            )

    def test_adjacent_publication_overlap_is_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            make_approved_project(target)
            assessment_path = target / "registry" / "route_assessments.tsv"
            fields, rows = read_tsv_rows(assessment_path)
            rows[0]["publication_overlap"] = "adjacent"
            write_tsv_rows(assessment_path, fields, rows)
            overlap_path = target / "registry" / "publication_overlap.tsv"
            fields, rows = read_tsv_rows(overlap_path)
            rows[0]["overlap_level"] = "adjacent"
            rows[0]["decision"] = "proceed_with_distinction"
            write_tsv_rows(overlap_path, fields, rows)
            report = validate_project(target)
            self.assertTrue(report.ok, report.errors)

    def test_substantive_duplicate_is_a_hard_stop(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            make_approved_project(target)
            assessment_path = target / "registry" / "route_assessments.tsv"
            fields, rows = read_tsv_rows(assessment_path)
            rows[0]["publication_overlap"] = "duplicate"
            rows[0]["execution_priority"] = "P3"
            write_tsv_rows(assessment_path, fields, rows)
            overlap_path = target / "registry" / "publication_overlap.tsv"
            fields, rows = read_tsv_rows(overlap_path)
            rows[0]["overlap_level"] = "duplicate"
            rows[0]["decision"] = "stop"
            write_tsv_rows(overlap_path, fields, rows)
            report = validate_project(target)
            self.assertFalse(report.ok)
            self.assertTrue(
                any("publication duplicate" in e for e in report.errors)
            )

    def test_route_effect_requires_change_request(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            target = Path(temp_dir) / "project"
            make_approved_project(target)
            append_tsv(
                target / "registry" / "claims.tsv",
                [
                    "CLAIM-001", "ROUTE-001", "approved", "Fixture claim",
                    "patient", "patient-level estimate", "associational", "current",
                ],
            )
            append_tsv(
                target / "registry" / "work_items.tsv",
                [
                    "WORK-001", "Test analysis", "completed", "manuscript_claim",
                    "CLAIM-001", "Figure 1 result", "Verified source table",
                    "Stop after the prespecified test", "high", "Fixture",
                ],
            )
            append_tsv(
                target / "registry" / "modules.tsv",
                [
                    "MODULE-001", "WORK-001", "ROUTE-001", "CLAIM-001",
                    "verified", "Primary module", "DATA-001 patient table",
                    "patient", "", "CODEMAP-001", "current",
                ],
            )
            append_tsv(
                target / "registry" / "runs.tsv",
                [
                    "RUN-001", "completed", "abc123", "config", "data",
                    "", "2026-08-02T00:00:00Z", "2026-08-02T00:01:00Z", "Fixture",
                ],
            )
            append_tsv(
                target / "registry" / "results.tsv",
                [
                    "RESULT-001", "RUN-001", "MODULE-001", "CLAIM-001",
                    "verified", "current", "refine", "associational",
                    "Sensitivity result requires a bounded update.",
                    "Revise the limitation and sensitivity figure.",
                    "outputs/result.tsv", "2026-08-02T00:01:00Z",
                ],
            )
            missing = validate_project(target)
            self.assertFalse(missing.ok)
            self.assertTrue(
                any("no change request" in e for e in missing.errors)
            )
            append_tsv(
                target / "registry" / "change_requests.tsv",
                [
                    "CHANGE-001", "WORK-001", "RESULT-001", "claim", "medium",
                    "CLAIM-001", "Add the prespecified sensitivity limitation.",
                    "proposed", "G5_AUDIT", "", "2026-08-02T00:02:00Z",
                ],
            )
            linked = validate_project(target)
            self.assertTrue(linked.ok, linked.errors)

    def test_downstream_impact_is_retained(self) -> None:
        impacted = downstream_impact(SPP1_PROJECT, "PAPER-SPP1-001")
        self.assertEqual([item["entity_id"] for item in impacted], ["WORK-SPP1-001"])

    def test_active_schema_has_no_legacy_direction_registries(self) -> None:
        registries = load_contract()["registries"]
        self.assertNotIn("directions", registries)
        self.assertNotIn("direction_assessments", registries)
        self.assertNotIn("novelty", json.dumps(registries))


if __name__ == "__main__":
    unittest.main()
