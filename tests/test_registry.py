from __future__ import annotations

import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from paper2paper.registry import (
    capability_coverage,
    derive_dataset_key,
    validate_registries,
)


ROOT = Path(__file__).resolve().parents[1]
GASTRIC_PROJECT = "P2P-GASTRIC-NRRS"
RELEASE = "MOD-BULK-FIXED-SIGNATURE@0.1.0"


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], [dict(row) for row in reader]


def write_rows(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def rewrite(path: Path, change) -> None:
    header, rows = read_rows(path)
    change(rows)
    write_rows(path, header, rows)


def clone_registry_world(destination: Path) -> Path:
    root = destination / "repo"
    root.mkdir(parents=True)
    shutil.copytree(ROOT / "registries", root / "registries")
    shutil.copytree(ROOT / "modules", root / "modules")
    shutil.copytree(
        ROOT / "pilots/gastric-nrrs",
        root / "pilots/gastric-nrrs",
        ignore=shutil.ignore_patterns("raw_data", "__pycache__", "*.pyc"),
    )
    return root


def refresh_release_hash(root: Path, field: str, relative: str) -> None:
    payload = (root / relative).read_bytes().replace(b"\r\n", b"\n").replace(
        b"\r", b"\n"
    )
    digest = hashlib.sha256(payload).hexdigest().upper()
    rewrite(
        root / "registries/module_releases.tsv",
        lambda rows: rows[0].update({field: digest}),
    )


def add_transfer_workspace(
    root: Path,
    *,
    container: str = "pilots",
    slug: str = "signature-transfer",
    project_id: str = "P2P-SIGNATURE-TRANSFER",
    doi: str = "10.9999/independent-anchor",
    add_case: bool = True,
) -> str:
    workspace = root / container / slug
    (workspace / "evidence").mkdir(parents=True)
    (workspace / "execution").mkdir(parents=True)
    (workspace / "inputs").mkdir(parents=True)
    (workspace / "tests").mkdir(parents=True)
    (workspace / "analysis").mkdir(parents=True)
    kind = "pilot" if container == "pilots" else "manuscript_project"
    (workspace / "PROJECT.json").write_text(
        json.dumps(
            {
                "schema_version": "2.0.0",
                "workspace_kind": kind,
                "project_id": project_id,
                "anchor": {"doi": doi},
            }
        ),
        encoding="utf-8",
    )
    (workspace / "analysis/specification.md").write_text(
        "# Transfer specification\n", encoding="utf-8"
    )
    (workspace / "evidence/routes.tsv").write_text(
        "route_id\tdecision_status\tcapability_ids\n"
        "TR-R01\tactive\tCAP-BULK-FIXED-SIGNATURE\n",
        encoding="utf-8",
    )
    (workspace / "evidence/code_requirements.tsv").write_text(
        "module_id\troute_id\tcapability_id\n"
        "TR-M01\tTR-R01\tCAP-BULK-FIXED-SIGNATURE\n",
        encoding="utf-8",
    )
    (workspace / "evidence/code_candidates.tsv").write_text(
        "code_id\troute_id\tmodule_id\tdecision\treusable_release_id\n"
        f"TR-C01\tTR-R01\tTR-M01\tuse\t{RELEASE}\n",
        encoding="utf-8",
    )
    (workspace / "evidence/issues.tsv").write_text(
        "issue_id\troute_id\tcandidate_scope\tdisposition\t"
        "affected_capability_ids\tpromotion_id\n",
        encoding="utf-8",
    )
    (workspace / "evidence/data_candidates.tsv").write_text(
        "data_id\troute_id\tdecision\taccession\turi\tsource\tdisease\t"
        "tissue\tmodality\tsubjects\tbiological_unit\tchecksum\n"
        "TR-D01\tTR-R01\tuse\tTR-DATA-001\thttps://example.org/tr-data\t"
        "independent fixture\tindependent disease\ttumor\tbulk expression\t42\t"
        "patient\tSHA256:1111\n",
        encoding="utf-8",
    )
    (workspace / "evidence/data_resources.tsv").write_text(
        "resource_id\tdata_id\turi\tsource_version\tbytes\tchecksum\t"
        "fields_supplied\tidentifier_field\n"
        "TR-RES01\tTR-D01\thttps://example.org/tr-data\tv1\t42\t"
        "SHA256:1111\texpression;outcome\tpatient_id\n",
        encoding="utf-8",
    )
    (workspace / "evidence/cohort_usage.tsv").write_text(
        "usage_id\troute_id\tdata_id\tcohort_key\tanalysis_step\trole\t"
        "outcome_used\tfeatures_influenced\tparameters_influenced\t"
        "cutoff_influenced\n"
        "TR-U01\tTR-R01\tTR-D01\tTR-COHORT\ttransfer check\tvalidation\t"
        "true\tfalse\tfalse\tfalse\n",
        encoding="utf-8",
    )
    manifest_rel = f"{container}/{slug}/inputs/manifest.tsv"
    assertion_rel = f"{container}/{slug}/tests/assertions.py"
    (workspace / "inputs/manifest.tsv").write_text(
        "resource\tsha256\nTR-DATA-001\t1111\n", encoding="utf-8"
    )
    (workspace / "tests/assertions.py").write_text(
        "def test_independent_transfer():\n"
        "    rows = [dict(patient_id='P1', score=0.2), "
        "dict(patient_id='P2', score=-0.1)]\n"
        "    assert len({row['patient_id'] for row in rows}) == len(rows)\n"
        "    assert all(isinstance(row['score'], float) for row in rows)\n",
        encoding="utf-8",
    )
    (workspace / "execution/runs.tsv").write_text(
        "run_id\troute_id\trun_kind\tstatus\tcommand\tdata_manifest\t"
        "module_release_ids\n"
        f"TR-RUN01\tTR-R01\treal_data\tpassed\tpython {assertion_rel}\t"
        f"{manifest_rel}\t{RELEASE}\n",
        encoding="utf-8",
    )
    if not add_case:
        return ""
    dataset_key = derive_dataset_key(root, project_id, "TR-R01", manifest_rel)
    path = root / "registries/regression_cases.tsv"
    header, rows = read_rows(path)
    rows.append(
        {
            "case_id": "TR-RC-SIGNATURE",
            "target_scope": "module",
            "target_id": RELEASE,
            "capability_id": "CAP-BULK-FIXED-SIGNATURE",
            "pilot_id": project_id,
            "route_id": "TR-R01",
            "case_role": "transfer",
            "test_level": "real_data",
            "command": f"python {assertion_rel}",
            "input_manifest": manifest_rel,
            "assertion_path": assertion_rel,
            "covered_test_ids": "",
            "oracle_source": "property",
            "oracle_ref": "independent transfer property",
            "anchor_key": doi,
            "dataset_key": dataset_key,
            "lifecycle": "active",
            "last_outcome": "passed",
            "last_run_id": "TR-RUN01",
            "notes": "independent transfer fixture",
        }
    )
    write_rows(path, header, rows)
    return dataset_key


class RegistryTests(unittest.TestCase):
    def test_repository_registries_are_valid_and_reference_only(self) -> None:
        report = validate_registries(ROOT)
        self.assertTrue(report.ok, report.errors)
        self.assertEqual(
            capability_coverage(ROOT),
            {
                "CAP-RESOURCE-INTEGRITY": "reference_exercised",
                "CAP-PATIENT-ID-LINKAGE": "reference_exercised",
                "CAP-BULK-FIXED-SIGNATURE": "reference_exercised",
            },
        )

    def test_single_pilot_cannot_claim_transfer_or_confirmed_promotion(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "registries/module_releases.tsv",
                lambda rows: rows[0].update({"maturity": "transfer_verified"}),
            )
            rewrite(
                root / "registries/promotions.tsv",
                lambda rows: rows[2].update(
                    {
                        "maturity": "confirmed",
                        "authority": "user",
                        "review_status": "approved",
                    }
                ),
            )
            report = validate_registries(root)
            self.assertFalse(report.ok)
            self.assertTrue(any("transfer_verified requires" in x for x in report.errors))
            self.assertTrue(any("confirmed capability" in x for x in report.errors))

    def test_transfer_independence_is_derived_not_free_text(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            actual_dataset = add_transfer_workspace(root)
            rewrite(
                root / "registries/module_releases.tsv",
                lambda rows: rows[0].update({"maturity": "transfer_verified"}),
            )
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: rows[-1].update(
                    {"anchor_key": "invented-anchor", "dataset_key": "invented-data"}
                ),
            )
            report = validate_registries(root)
            self.assertFalse(report.ok)
            self.assertTrue(any("PROJECT.json DOI" in x for x in report.errors))
            self.assertTrue(any("dataset_key does not match" in x for x in report.errors))

            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: rows[-1].update(
                    {
                        "anchor_key": "10.9999/independent-anchor",
                        "dataset_key": actual_dataset,
                    }
                ),
            )
            rewrite(
                root / "registries/promotions.tsv",
                lambda rows: rows[2].update(
                    {
                        "required_case_ids": rows[2]["required_case_ids"]
                        + ";TR-RC-SIGNATURE",
                        "maturity": "confirmed",
                        "authority": "user",
                        "review_status": "approved",
                    }
                ),
            )
            report = validate_registries(root)
            self.assertTrue(report.ok, report.errors)
            self.assertEqual(
                capability_coverage(root)["CAP-BULK-FIXED-SIGNATURE"],
                "transfer_exercised",
            )

    def test_reviewed_freeze_and_missing_unit_cannot_advance_release(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: [
                    row.update({"oracle_source": "reviewed_freeze"})
                    for row in rows
                    if row["case_id"] == "GN-RC-SIGNATURE"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(
                any("stronger than reviewed_freeze" in x for x in report.errors),
                report.errors,
            )

            root = clone_registry_world(Path(temp) / "second")
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: [
                    row.update({"last_outcome": "not_run", "last_run_id": ""})
                    for row in rows
                    if row["case_id"] in {
                        "GN-RC-SIGNATURE-UNIT",
                        "GN-RC-SIGNATURE-ADVERSARIAL",
                        "GN-RC-SIGNATURE-RESOURCE-ADVERSARIAL",
                    }
                ],
            )
            report = validate_registries(root)
            self.assertTrue(any("lacks a passed qualifying unit case" in x for x in report.errors))

    def test_workspace_capability_release_and_promotion_links_are_bidirectional(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "pilots/gastric-nrrs/evidence/routes.tsv",
                lambda rows: rows[0].update({"capability_ids": "CAP-NOT-REGISTERED"}),
            )
            report = validate_registries(root)
            self.assertTrue(any("unknown capability CAP-NOT-REGISTERED" in x for x in report.errors))
            self.assertTrue(any("release capability is not declared" in x for x in report.errors))

            root = clone_registry_world(Path(temp) / "second")
            rewrite(
                root / "pilots/gastric-nrrs/evidence/issues.tsv",
                lambda rows: [
                    row.update({"promotion_id": ""})
                    for row in rows
                    if row["issue_id"] == "PF-004"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(any("lacks promotion_id" in x for x in report.errors))
            self.assertTrue(any("lacks reciprocal link" in x for x in report.errors))

    def test_promotion_cannot_use_unrelated_capability_or_target_case(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "registries/promotions.tsv",
                lambda rows: rows[2].update(
                    {
                        "required_case_ids": "GN-RC-RESOURCE",
                        "maturity": "confirmed",
                        "authority": "user",
                        "review_status": "approved",
                    }
                ),
            )
            report = validate_registries(root)
            self.assertTrue(any("unrelated capability" in x for x in report.errors))
            self.assertTrue(any("does not target its release" in x for x in report.errors))

    def test_case_outcome_must_bind_a_real_matching_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "pilots/gastric-nrrs/execution/runs.tsv",
                lambda rows: [
                    row.update({"status": "failed"})
                    for row in rows
                    if row["run_id"] == "GN-RUN06"
                ],
            )
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: [
                    row.update({"last_run_id": "GN-RUN06"})
                    for row in rows
                    if row["case_id"] == "GN-RC-IDENTITY-ADVERSARIAL"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(any("points to run status 'failed'" in x for x in report.errors))
            self.assertTrue(any("not_run must not record last_run_id" in x for x in report.errors))

    def test_project_id_is_not_a_directory_slug_and_unknown_workspace_does_not_cascade(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: rows[0].update({"pilot_id": "gastric-nrrs"}),
            )
            report = validate_registries(root)
            self.assertTrue(any("unknown workspace project_id" in x for x in report.errors))
            self.assertFalse(
                any(
                    "GN-RC-RESOURCE" in x and "non-real-data case" in x
                    for x in report.errors
                ),
                report.errors,
            )

    def test_contract_tests_must_exist_and_be_covered_by_passed_cases(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "modules/bulk-fixed-signature/contract.tsv",
                lambda rows: rows[0].update(
                    {"test_id": "test_metadata_path_never_exposes_an_absolute_directory"}
                ),
            )
            refresh_release_hash(
                root,
                "contract_sha256",
                "modules/bulk-fixed-signature/contract.tsv",
            )
            report = validate_registries(root)
            self.assertTrue(any("required contract tests lack passed case coverage" in x for x in report.errors))

            root = clone_registry_world(Path(temp) / "second")
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: [
                    row.update({"covered_test_ids": "test_does_not_exist"})
                    for row in rows
                    if row["case_id"] == "GN-RC-SIGNATURE-UNIT"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(any("unknown test_id 'test_does_not_exist'" in x for x in report.errors))

    def test_passed_unit_command_must_resolve_executable_class_and_test(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: [
                    row.update(
                        {
                            "command": "python -m unittest "
                            "pilots.gastric-nrrs.code.test_gse62254_nrrs_spike."
                            "MissingClass.test_resource_verification_fails_closed_on_change"
                        }
                    )
                    for row in rows
                    if row["case_id"] == "GN-RC-RESOURCE-ADVERSARIAL"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(any("unknown class 'MissingClass'" in x for x in report.errors))

            root = clone_registry_world(Path(temp) / "second")
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: [
                    row.update(
                        {
                            "command": "not-python -m unittest "
                            "pilots/gastric-nrrs/code/test_gse62254_nrrs_spike.py"
                        }
                    )
                    for row in rows
                    if row["case_id"] == "GN-RC-SIGNATURE-UNIT"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(any("unsupported executable" in x for x in report.errors))

    def test_release_hashes_reject_stale_implementation_contract_environment_or_tests(self) -> None:
        for field, relative in (
            ("implementation_sha256", "pilots/gastric-nrrs/code/gse62254_nrrs_spike.py"),
            ("contract_sha256", "modules/bulk-fixed-signature/contract.tsv"),
            ("environment_sha256", "pilots/gastric-nrrs/code/requirements-gse62254.txt"),
            ("unit_test_sha256", "pilots/gastric-nrrs/code/test_gse62254_nrrs_spike.py"),
        ):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temp:
                root = clone_registry_world(Path(temp))
                path = root / relative
                path.write_text(path.read_text(encoding="utf-8") + "\n# drift\n", encoding="utf-8")
                report = validate_registries(root)
                self.assertTrue(any(f"{field} is stale" in x for x in report.errors))

    def test_release_hashes_are_stable_across_text_line_endings(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            _, releases = read_rows(root / "registries/module_releases.tsv")
            release = releases[0]
            for field in (
                "entrypoint", "environment_lock", "contract_path", "unit_test_path"
            ):
                path = root / release[field]
                payload = path.read_bytes().replace(b"\r\n", b"\n").replace(
                    b"\r", b"\n"
                )
                path.write_bytes(payload.replace(b"\n", b"\r\n"))
            report = validate_registries(root)
            self.assertTrue(report.ok, report.errors)

    def test_dataset_keys_are_stable_across_manifest_line_endings(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            manifest = (
                root / "pilots/gastric-nrrs/config/gse62254_resources.tsv"
            )
            before = derive_dataset_key(
                root,
                "P2P-GASTRIC-NRRS",
                "GN-R01",
                "pilots/gastric-nrrs/config/gse62254_resources.tsv",
            )
            payload = manifest.read_bytes().replace(b"\r\n", b"\n").replace(
                b"\r", b"\n"
            )
            manifest.write_bytes(payload.replace(b"\n", b"\r\n"))
            after = derive_dataset_key(
                root,
                "P2P-GASTRIC-NRRS",
                "GN-R01",
                "pilots/gastric-nrrs/config/gse62254_resources.tsv",
            )
            self.assertEqual(before, after)
            report = validate_registries(root)
            self.assertTrue(report.ok, report.errors)

    def test_regression_case_level_must_match_bound_run_kind(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "pilots/gastric-nrrs/execution/runs.tsv",
                lambda rows: [
                    row.update({"run_kind": "unit"})
                    for row in rows
                    if row["run_id"] == "GN-RUN05"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(
                any("real_data case must bind a real_data run" in x for x in report.errors)
            )

            root = clone_registry_world(Path(temp) / "second")
            rewrite(
                root / "pilots/gastric-nrrs/execution/runs.tsv",
                lambda rows: [
                    row.update({"run_kind": "real_data"})
                    for row in rows
                    if row["run_id"] == "GN-RUN06"
                ],
            )
            report = validate_registries(root)
            self.assertTrue(
                any("unit case must bind a unit run" in x for x in report.errors)
            )

    def test_real_data_case_command_must_match_bound_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: rows[0].update(
                    {"command": "python unrelated_assertion.py"}
                ),
            )
            report = validate_registries(root)
            self.assertTrue(
                any("command differs from its bound run" in x for x in report.errors)
            )

    def test_unsafe_paths_and_missing_foreign_keys_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = clone_registry_world(Path(temp))
            rewrite(
                root / "registries/module_releases.tsv",
                lambda rows: rows[0].update({"path": "C:\\Users\\example\\module"}),
            )
            rewrite(
                root / "registries/regression_cases.tsv",
                lambda rows: rows[0].update({"capability_id": "CAP-MISSING"}),
            )
            report = validate_registries(root)
            self.assertTrue(any("unsafe or absolute path" in x for x in report.errors))
            self.assertTrue(any("unknown capability" in x for x in report.errors))


if __name__ == "__main__":
    unittest.main()
