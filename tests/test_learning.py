from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from paper2paper.learning import (
    FINDING_COLUMNS,
    learning_summary,
    validate_learning,
    write_learning_report,
)
from paper2paper.workspace import init_workspace


def append_tsv_row(path: Path, values: dict[str, str]) -> None:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = csv.DictReader(handle, delimiter="\t").fieldnames or []
    unknown = set(values) - set(header)
    if unknown:
        raise ValueError(f"unknown fields: {sorted(unknown)}")
    row = {field: "" for field in header}
    row.update(values)
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writerow(row)


def add_workflow_issue(
    pilot: Path,
    issue_id: str,
    *,
    observation: str,
    action: str = "consider",
    source_scope: str = "pilot_execution",
) -> None:
    append_tsv_row(
        pilot / "evidence/issues.tsv",
        {
            "issue_id": issue_id,
            "source_scope": source_scope,
            "issue_type": "usability",
            "stage": "route_generation",
            "severity": "high",
            "observation": observation,
            "evidence": "real Pilot evidence",
            "consequence": "The paper could receive a wrong workflow decision.",
            "proposed_action": "Classify the issue before changing core behavior.",
            "current_resolution": "The current paper is repaired locally first.",
            "workflow_action": action,
            "status": "open" if action != "implemented" else "resolved",
            "blocking": "false",
        },
    )


def prepare_repo(root: Path) -> tuple[Path, Path]:
    (root / "pilots").mkdir(parents=True)
    (root / "manuscript-projects").mkdir()
    (root / "learning").mkdir()
    pilot_one = root / "pilots/one"
    pilot_two = root / "pilots/two"
    init_workspace(pilot_one, "P2P-ONE", "One", "Anchor one")
    init_workspace(pilot_two, "P2P-TWO", "Two", "Anchor two")
    return pilot_one, pilot_two


def write_findings(root: Path, rows: list[dict[str, str]]) -> None:
    path = root / "learning/findings.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(FINDING_COLUMNS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for values in rows:
            row = {field: "" for field in FINDING_COLUMNS}
            row.update(values)
            writer.writerow(row)


def finding(
    finding_id: str,
    source_refs: str,
    *,
    scope: str = "general",
    status: str = "accepted",
    verification_refs: str = "",
) -> dict[str, str]:
    return {
        "finding_id": finding_id,
        "title": "A classified workflow finding",
        "scope": scope,
        "capability_area": "evidence_routing",
        "status": status,
        "source_refs": source_refs,
        "generalized_failure": "The workflow lacks a reusable check for this failure class.",
        "risk_to_paper": "A target paper can receive an incorrect decision.",
        "core_response": "Apply the smallest check and rerun the source Pilot.",
        "verification_refs": verification_refs,
        "updated_at": "2026-08-09",
        "notes": "The source observation remains in its Pilot.",
    }


class Paper2PaperLearningTests(unittest.TestCase):
    def test_aggregation_across_pilots_and_untriaged_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pilot_one, pilot_two = prepare_repo(root)
            add_workflow_issue(
                pilot_one, "ISSUE-ONE", observation="First paper exposed a gap."
            )
            add_workflow_issue(
                pilot_two, "ISSUE-TWO", observation="Second paper confirmed the gap."
            )
            add_workflow_issue(
                pilot_two,
                "ISSUE-LOCAL",
                observation="Second paper also has an unclassified local obstacle.",
            )
            write_findings(
                root,
                [
                    finding(
                        "P2P-L-ONE",
                        "P2P-ONE#ISSUE-ONE;P2P-TWO#ISSUE-TWO",
                        status="verified",
                        verification_refs="tests/test_learning.py::source-pilot-rerun",
                    )
                ],
            )

            summary = learning_summary(root)

            self.assertTrue(summary["ok"], summary["validation_errors"])
            self.assertEqual(summary["pilots_scanned"], 2)
            self.assertEqual(summary["workflow_candidates"], 3)
            self.assertEqual(summary["classified_source_issues"], 2)
            self.assertEqual(
                [row["source_ref"] for row in summary["untriaged_candidates"]],
                ["P2P-TWO#ISSUE-LOCAL"],
            )

    def test_learning_report_separates_raw_candidates_from_product_findings(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pilot_one, pilot_two = prepare_repo(root)
            add_workflow_issue(
                pilot_one, "ISSUE-PRODUCT", observation="Reusable failure observed."
            )
            add_workflow_issue(
                pilot_two, "ISSUE-RAW", observation="Topic-specific obstacle observed."
            )
            write_findings(
                root,
                [finding("P2P-L-PRODUCT", "P2P-ONE#ISSUE-PRODUCT")],
            )

            path = write_learning_report(root)
            text = path.read_text(encoding="utf-8")

            self.assertIn("候选不自动等于 Paper2Paper 产品缺口", text)
            self.assertIn("P2P-L-PRODUCT", text)
            self.assertIn("P2P-TWO#ISSUE-RAW", text)
            self.assertIn("课题特有、论文类型级还是通用", text)

    def test_unknown_source_reference_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pilot_one, _ = prepare_repo(root)
            add_workflow_issue(
                pilot_one, "ISSUE-ONE", observation="A real issue exists."
            )
            write_findings(
                root,
                [finding("P2P-L-BAD", "P2P-NOT-REAL#ISSUE-NOT-REAL")],
            )

            report = validate_learning(root)

            self.assertFalse(report.ok)
            self.assertTrue(
                any("unknown source_ref" in error for error in report.errors)
            )

    def test_anchor_issue_marked_for_workflow_review_is_not_lost(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pilot_one, _ = prepare_repo(root)
            add_workflow_issue(
                pilot_one,
                "ISSUE-ANCHOR",
                observation="An anchor-only limitation may expose a reusable check.",
                source_scope="anchor",
            )
            write_findings(root, [])

            summary = learning_summary(root)

            self.assertTrue(summary["ok"], summary["validation_errors"])
            self.assertEqual(summary["workflow_candidates"], 1)
            self.assertEqual(summary["classified_source_issues"], 0)
            self.assertEqual(
                summary["untriaged_candidates"][0]["source_ref"],
                "P2P-ONE#ISSUE-ANCHOR",
            )

    def test_verified_finding_requires_verification_reference(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pilot_one, _ = prepare_repo(root)
            add_workflow_issue(
                pilot_one, "ISSUE-ONE", observation="A real issue exists."
            )
            write_findings(
                root,
                [
                    finding(
                        "P2P-L-UNVERIFIED",
                        "P2P-ONE#ISSUE-ONE",
                        status="verified",
                    )
                ],
            )

            report = validate_learning(root)

            self.assertFalse(report.ok)
            self.assertTrue(
                any("requires verification_refs" in error for error in report.errors)
            )

    def test_paper_specific_finding_cannot_claim_product_implementation(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            pilot_one, _ = prepare_repo(root)
            add_workflow_issue(
                pilot_one, "ISSUE-ONE", observation="A local-only issue exists."
            )
            write_findings(
                root,
                [
                    finding(
                        "P2P-L-LOCAL",
                        "P2P-ONE#ISSUE-ONE",
                        scope="paper_specific",
                        status="implemented",
                    )
                ],
            )

            report = validate_learning(root)

            self.assertFalse(report.ok)
            self.assertTrue(
                any("paper_specific finding" in error for error in report.errors)
            )


if __name__ == "__main__":
    unittest.main()
