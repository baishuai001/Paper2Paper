from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

from paper2paper.regression import build_regression_plan


ROOT = Path(__file__).resolve().parents[1]
ALL_CASES = {
    "GN-RC-RESOURCE",
    "GN-RC-RESOURCE-ADVERSARIAL",
    "GN-RC-IDENTITY",
    "GN-RC-IDENTITY-ADVERSARIAL",
    "GN-RC-SIGNATURE",
    "GN-RC-SIGNATURE-FREEZE",
    "GN-RC-SIGNATURE-UNIT",
    "GN-RC-SIGNATURE-ADVERSARIAL",
    "GN-RC-SIGNATURE-RESOURCE-ADVERSARIAL",
}
SIGNATURE_CASES = {
    "GN-RC-SIGNATURE",
    "GN-RC-SIGNATURE-FREEZE",
    "GN-RC-SIGNATURE-UNIT",
    "GN-RC-SIGNATURE-ADVERSARIAL",
    "GN-RC-SIGNATURE-RESOURCE-ADVERSARIAL",
}


def clone_registry_world(destination: Path) -> Path:
    root = destination / "repo"
    root.mkdir()
    shutil.copytree(ROOT / "registries", root / "registries")
    shutil.copytree(ROOT / "modules", root / "modules")
    shutil.copytree(
        ROOT / "pilots/gastric-nrrs",
        root / "pilots/gastric-nrrs",
        ignore=shutil.ignore_patterns("raw_data", "__pycache__", "*.pyc"),
    )
    return root


def add_manuscript_workspace(root: Path) -> str:
    workspace = root / "manuscript-projects/signature-manuscript"
    (workspace / "analysis").mkdir(parents=True)
    (workspace / "evidence").mkdir(parents=True)
    (workspace / "PROJECT.json").write_text(
        json.dumps(
            {
                "schema_version": "2.0.0",
                "workspace_kind": "manuscript_project",
                "project_id": "P2P-MS-SIGNATURE",
                "anchor": {"doi": "10.9999/manuscript-anchor"},
            }
        ),
        encoding="utf-8",
    )
    (workspace / "evidence/routes.tsv").write_text(
        "route_id\tdecision_status\tcapability_ids\n"
        "MS-R01\tactive\tCAP-BULK-FIXED-SIGNATURE\n",
        encoding="utf-8",
    )
    (workspace / "analysis/specification.md").write_text(
        "# Formal signature specification\n", encoding="utf-8"
    )
    return "manuscript-projects/signature-manuscript/analysis/specification.md"


class RegressionPlanningTests(unittest.TestCase):
    def test_regular_pilot_semantic_change_selects_workspace_bindings(self) -> None:
        for changed in (
            "pilots/gastric-nrrs/config/gse62254_resources.tsv",
            "pilots/gastric-nrrs/evidence/routes.tsv",
            "pilots/gastric-nrrs/analysis/specification.md",
            "pilots/gastric-nrrs/PROJECT.json",
        ):
            with self.subTest(changed=changed):
                plan = build_regression_plan(ROOT, [changed])
                self.assertEqual(set(plan.case_ids), ALL_CASES)
                self.assertTrue(all(item.reasons for item in plan.selected))

    def test_release_entrypoint_contract_and_environment_select_shared_release_cases(self) -> None:
        for changed in (
            "pilots/gastric-nrrs/code/gse62254_nrrs_spike.py",
            "modules/bulk-fixed-signature/contract.tsv",
            "pilots/gastric-nrrs/code/requirements-gse62254.txt",
            "pilots/gastric-nrrs/code/test_gse62254_nrrs_spike.py",
        ):
            with self.subTest(changed=changed):
                plan = build_regression_plan(ROOT, [changed])
                self.assertEqual(set(plan.case_ids), SIGNATURE_CASES)

    def test_promotion_change_includes_required_and_counterexample_cases(self) -> None:
        plan = build_regression_plan(ROOT, ["registries/promotions.tsv"])
        self.assertEqual(
            set(plan.case_ids),
            ALL_CASES - {"GN-RC-SIGNATURE-FREEZE"},
        )

    def test_core_and_ci_changes_select_all_active_cases(self) -> None:
        for changed in (
            "src/paper2paper/workspace.py",
            ".github/workflows/ci.yml",
            "./.github/workflows/ci.yml",
        ):
            with self.subTest(changed=changed):
                plan = build_regression_plan(ROOT, [changed])
                self.assertEqual(set(plan.case_ids), ALL_CASES)
                self.assertEqual(plan.unmatched_paths, [])
                if changed.startswith("./"):
                    self.assertEqual(plan.changed_paths, [".github/workflows/ci.yml"])

    def test_manuscript_semantic_change_selects_capability_cases(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            changed = add_manuscript_workspace(root)
            plan = build_regression_plan(root, [changed])
            self.assertEqual(set(plan.case_ids), SIGNATURE_CASES)

    def test_documentation_change_does_not_request_real_data(self) -> None:
        for changed in (
            "docs/workflow.md",
            "README.md",
            "pilots/gastric-nrrs/reports/pilot-findings.md",
            "manuscript-projects/README.md",
        ):
            with self.subTest(changed=changed):
                plan = build_regression_plan(ROOT, [changed])
                self.assertEqual(plan.case_ids, [])
                self.assertEqual(plan.unmatched_paths, [])

    def test_unknown_path_is_reported_without_broad_rerun(self) -> None:
        plan = build_regression_plan(ROOT, ["scratch/idea.txt"])
        self.assertEqual(plan.case_ids, [])
        self.assertEqual(plan.unmatched_paths, ["scratch/idea.txt"])

    def test_absolute_or_traversal_path_is_rejected(self) -> None:
        for changed in ("C:\\Users\\example\\file.py", "../outside.py"):
            with self.subTest(changed=changed), self.assertRaises(ValueError):
                build_regression_plan(ROOT, [changed])


if __name__ == "__main__":
    unittest.main()
