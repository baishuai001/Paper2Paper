from __future__ import annotations

import ast
import csv
import hashlib
import json
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath, PureWindowsPath
from typing import Iterable

from .schema import SCHEMA_VERSION


REGISTRY_HEADERS: dict[str, tuple[str, ...]] = {
    "capabilities": (
        "capability_id",
        "family",
        "name",
        "input_modalities",
        "biological_unit",
        "output_kind",
        "claim_ceiling",
        "out_of_scope",
        "lifecycle",
        "notes",
    ),
    "module_releases": (
        "release_id",
        "module_id",
        "capability_id",
        "name",
        "version",
        "maturity",
        "lifecycle",
        "path",
        "entrypoint",
        "implementation_sha256",
        "language",
        "environment_lock",
        "environment_sha256",
        "license",
        "contract_path",
        "contract_sha256",
        "unit_test_command",
        "unit_test_path",
        "unit_test_sha256",
        "source_provenance",
        "supersedes",
        "notes",
    ),
    "promotions": (
        "promotion_id",
        "source_refs",
        "target_scope",
        "target_id",
        "affected_capability_ids",
        "change_statement",
        "risk_class",
        "evidence_refs",
        "required_case_ids",
        "counterexample_case_ids",
        "maturity",
        "authority",
        "review_status",
        "reviewer",
        "decided_at",
        "implementation_ref",
        "notes",
    ),
    "regression_cases": (
        "case_id",
        "target_scope",
        "target_id",
        "capability_id",
        "pilot_id",
        "route_id",
        "case_role",
        "test_level",
        "command",
        "input_manifest",
        "assertion_path",
        "covered_test_ids",
        "oracle_source",
        "oracle_ref",
        "anchor_key",
        "dataset_key",
        "lifecycle",
        "last_outcome",
        "last_run_id",
        "notes",
    ),
}

REGISTRY_FILES = {
    name: f"registries/{name}.tsv" for name in REGISTRY_HEADERS
}

ENUMS: dict[tuple[str, str], set[str]] = {
    ("capabilities", "lifecycle"): {"active", "deprecated"},
    ("module_releases", "maturity"): {
        "draft",
        "unit_verified",
        "reference_verified",
        "transfer_verified",
    },
    ("module_releases", "lifecycle"): {"active", "blocked", "deprecated"},
    ("promotions", "target_scope"): {"module", "core"},
    ("promotions", "risk_class"): {
        "safety",
        "scientific",
        "engineering",
        "usability",
    },
    ("promotions", "maturity"): {
        "observed",
        "provisional",
        "confirmed",
        "rejected",
        "deprecated",
    },
    ("promotions", "authority"): {"user", "ai", "collaborator"},
    ("promotions", "review_status"): {"pending", "approved", "rejected"},
    ("regression_cases", "target_scope"): {"capability", "module", "core"},
    ("regression_cases", "case_role"): {
        "reference",
        "transfer",
        "adversarial",
        "safety",
    },
    ("regression_cases", "test_level"): {
        "structure",
        "unit",
        "real_data",
        "scientific_review",
    },
    ("regression_cases", "oracle_source"): {
        "published",
        "independent_implementation",
        "reviewed_freeze",
        "property",
        "failure_expected",
    },
    ("regression_cases", "lifecycle"): {"active", "quarantined", "retired"},
    ("regression_cases", "last_outcome"): {
        "not_run",
        "passed",
        "failed",
        "blocked",
    },
}

CONTRACT_HEADER = (
    "item_id",
    "kind",
    "name",
    "specification",
    "required",
    "test_id",
    "failure_action",
    "notes",
)
CONTRACT_KINDS = {"input", "output", "invariant", "assumption", "failure_condition"}
CONTRACT_ACTIONS = {"fail_closed", "limit_claim", "warn"}
BOOLEAN_VALUES = {"true", "false"}
NON_QUALIFYING_ORACLES = {"reviewed_freeze"}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@-]*$")
SEMVER = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+(?:[-+][0-9A-Za-z.-]+)?$")
SHA256 = re.compile(r"^[0-9A-Fa-f]{64}$")
SOURCE_REF = re.compile(
    r"^(?P<project>[A-Za-z0-9][A-Za-z0-9._-]*)#"
    r"(?P<issue>[A-Za-z0-9][A-Za-z0-9._-]*)$"
)


@dataclass
class RegistryValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


@dataclass(frozen=True)
class WorkspaceRecord:
    project_id: str
    workspace_kind: str
    path: Path
    relative_path: str
    manifest: dict[str, object]


def _split_ids(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows: list[dict[str, str]] = []
        for raw in reader:
            row = {
                key: (value or "").strip()
                for key, value in raw.items()
                if key is not None
            }
            extra = raw.get(None)
            if extra:
                row["__extra_columns__"] = "\t".join(extra)
            rows.append(row)
    return header, rows


def _workspace_table(record: WorkspaceRecord, relative: str) -> list[dict[str, str]]:
    path = record.path / relative
    if not path.exists():
        return []
    return _read_tsv(path)[1]


def load_registries(repo_root: Path) -> dict[str, list[dict[str, str]]]:
    """Load the cross-paper registries without claiming that they are valid."""

    root = Path(repo_root)
    tables: dict[str, list[dict[str, str]]] = {}
    for name, relative in REGISTRY_FILES.items():
        path = root / relative
        tables[name] = _read_tsv(path)[1] if path.exists() else []
    return tables


def _unsafe_relative_path(value: str) -> bool:
    if not value:
        return False
    normalized = value.replace("\\", "/")
    posix = PurePosixPath(normalized)
    windows = PureWindowsPath(value)
    return (
        posix.is_absolute()
        or windows.is_absolute()
        or bool(windows.drive)
        or normalized.startswith("~")
        or ".." in posix.parts
    )


def _check_repo_path(
    repo_root: Path,
    value: str,
    label: str,
    report: RegistryValidationReport,
    *,
    required: bool = True,
    must_exist: bool = True,
) -> None:
    if not value:
        if required:
            report.errors.append(f"{label}: repository-relative path is required")
        return
    if _unsafe_relative_path(value):
        report.errors.append(f"{label}: unsafe or absolute path: {value}")
        return
    target = repo_root / value
    if must_exist and not target.exists():
        report.errors.append(f"{label}: referenced path does not exist: {value}")


def _index_rows(
    table_name: str,
    rows: list[dict[str, str]],
    id_field: str,
    report: RegistryValidationReport,
) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    for number, row in enumerate(rows, start=2):
        if row.get("__extra_columns__"):
            report.errors.append(f"{table_name}:{number}: row has surplus columns")
        identifier = row.get(id_field, "")
        if not identifier:
            report.errors.append(f"{table_name}:{number}: empty {id_field}")
            continue
        if not SAFE_ID.fullmatch(identifier):
            report.errors.append(
                f"{table_name}:{number}: unsafe or malformed {id_field} {identifier!r}"
            )
        if identifier in index:
            report.errors.append(f"{table_name}:{number}: duplicate {id_field} {identifier}")
            continue
        index[identifier] = row
    return index


def normalize_doi(value: str) -> str:
    normalized = value.strip().casefold()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix) :]
            break
    return normalized.rstrip("/ .")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256_canonical_text(path: Path) -> str:
    """Hash repository text after normalizing CRLF/CR line endings to LF."""

    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def _python_test_inventory(
    path: Path,
    label: str,
    report: RegistryValidationReport,
) -> tuple[set[str], dict[str, set[str]]]:
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except (OSError, SyntaxError, UnicodeError) as exc:
        report.errors.append(f"{label}: cannot parse Python assertions: {exc}")
        return set(), {}
    tests: set[str] = set()
    classes: dict[str, set[str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith(
            "test_"
        ):
            tests.add(node.name)
    for node in tree.body:
        if not isinstance(node, ast.ClassDef):
            continue
        classes[node.name] = {
            child.name
            for child in node.body
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
            and child.name.startswith("test_")
        }
    return tests, classes


def _validate_passed_unit_command(
    row: dict[str, str],
    assertion_path: str,
    tests: set[str],
    classes: dict[str, set[str]],
    report: RegistryValidationReport,
) -> None:
    case_id = row.get("case_id", "<unknown>")
    try:
        tokens = shlex.split(row.get("command", ""), posix=True)
    except ValueError as exc:
        report.errors.append(f"regression case {case_id}: invalid unit command: {exc}")
        return
    if len(tokens) < 4 or tokens[1:3] != ["-m", "unittest"]:
        report.errors.append(
            f"regression case {case_id}: passed unit command must invoke python -m unittest"
        )
        return
    executable = PurePosixPath(tokens[0].replace("\\", "/")).name.casefold()
    if executable not in {"python", "python.exe", "python3", "python3.exe", "py", "py.exe"}:
        report.errors.append(
            f"regression case {case_id}: unit command has unsupported executable"
        )
    normalized_tokens = [token.replace("\\", "/") for token in tokens]
    if assertion_path not in normalized_tokens:
        dotted_module = assertion_path.removesuffix(".py").replace("/", ".")
        dotted_targets = [
            token for token in tokens[3:] if token.startswith(dotted_module + ".")
        ]
        if dotted_targets:
            for target in dotted_targets:
                selector = target[len(dotted_module) + 1 :].split(".")
                class_name = selector[0] if selector else ""
                test_name = selector[1] if len(selector) > 1 else ""
                if class_name not in classes:
                    report.errors.append(
                        f"regression case {case_id}: unit command names unknown class "
                        f"{class_name!r}"
                    )
                elif test_name and test_name not in classes[class_name]:
                    report.errors.append(
                        f"regression case {case_id}: unit command names unknown test_id "
                        f"{test_name!r}"
                    )
            if any(not part.isidentifier() for part in dotted_module.split(".")):
                report.errors.append(
                    f"regression case {case_id}: dotted unittest module is not importable; "
                    "use the repository-relative test file path"
                )
        else:
            report.errors.append(
                f"regression case {case_id}: unit command does not resolve assertion_path"
            )
    for index, token in enumerate(tokens[:-1]):
        if token != "-k":
            continue
        selector = tokens[index + 1]
        if selector in tests:
            continue
        if "." in selector:
            class_name, test_name = selector.rsplit(".", 1)
            if class_name in classes and test_name in classes[class_name]:
                continue
            if class_name not in classes:
                report.errors.append(
                    f"regression case {case_id}: unit command names unknown class "
                    f"{class_name!r}"
                )
                continue
        report.errors.append(
            f"regression case {case_id}: unit command names unknown test_id {selector!r}"
        )


def _workspace_index(
    repo_root: Path,
    report: RegistryValidationReport | None = None,
) -> dict[str, WorkspaceRecord]:
    root = Path(repo_root)
    index: dict[str, WorkspaceRecord] = {}
    for parent_name, expected_kind in (
        ("pilots", "pilot"),
        ("manuscript-projects", "manuscript_project"),
    ):
        parent = root / parent_name
        if not parent.exists():
            continue
        for workspace_path in sorted(path for path in parent.iterdir() if path.is_dir()):
            manifest_path = workspace_path / "PROJECT.json"
            if not manifest_path.exists():
                continue
            label = workspace_path.relative_to(root).as_posix()
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                if report is not None:
                    report.errors.append(f"{label}/PROJECT.json: cannot parse: {exc}")
                continue
            if not isinstance(manifest, dict):
                if report is not None:
                    report.errors.append(f"{label}/PROJECT.json: root must be an object")
                continue
            project_id = manifest.get("project_id", "")
            kind = manifest.get("workspace_kind", "")
            if not isinstance(project_id, str) or not SAFE_ID.fullmatch(project_id):
                if report is not None:
                    report.errors.append(f"{label}: unsafe or missing project_id")
                continue
            if manifest.get("schema_version") != SCHEMA_VERSION and report is not None:
                report.errors.append(
                    f"{label}: registry linking requires schema {SCHEMA_VERSION}"
                )
            if kind != expected_kind and report is not None:
                report.errors.append(
                    f"{label}: workspace_kind {kind!r} does not match its directory"
                )
            if project_id in index:
                if report is not None:
                    report.errors.append(f"duplicate workspace project_id {project_id}")
                continue
            index[project_id] = WorkspaceRecord(
                project_id=project_id,
                workspace_kind=str(kind),
                path=workspace_path,
                relative_path=label,
                manifest=manifest,
            )
    return index


def load_workspace_index(repo_root: Path) -> dict[str, WorkspaceRecord]:
    return _workspace_index(Path(repo_root))


def _dataset_payload(
    repo_root: Path,
    record: WorkspaceRecord,
    route_id: str,
    input_manifest: str,
) -> dict[str, object]:
    candidates = [
        row
        for row in _workspace_table(record, "evidence/data_candidates.tsv")
        if row.get("route_id") == route_id and row.get("decision") == "use"
    ]
    data_ids = {row.get("data_id", "") for row in candidates}
    resources = [
        row
        for row in _workspace_table(record, "evidence/data_resources.tsv")
        if row.get("data_id") in data_ids
    ]
    usage = [
        row
        for row in _workspace_table(record, "evidence/cohort_usage.tsv")
        if row.get("route_id") == route_id and row.get("data_id") in data_ids
    ]
    candidate_fields = (
        "accession",
        "uri",
        "source",
        "disease",
        "tissue",
        "modality",
        "subjects",
        "biological_unit",
        "checksum",
    )
    resource_fields = (
        "uri",
        "source_version",
        "bytes",
        "checksum",
        "fields_supplied",
        "identifier_field",
    )
    usage_fields = (
        "cohort_key",
        "analysis_step",
        "role",
        "outcome_used",
        "features_influenced",
        "parameters_influenced",
        "cutoff_influenced",
    )
    manifest_sha = ""
    if input_manifest and not _unsafe_relative_path(input_manifest):
        manifest_path = repo_root / input_manifest
        if manifest_path.exists() and manifest_path.is_file():
            # The manifest is repository text, not the downloaded data itself.
            # Canonical line endings keep its identity stable across Windows and
            # Linux checkouts. Resource rows still carry raw-file SHA256 values.
            manifest_sha = _sha256_canonical_text(manifest_path)
    return {
        "candidates": sorted(
            [{field: row.get(field, "") for field in candidate_fields} for row in candidates],
            key=lambda item: json.dumps(item, sort_keys=True),
        ),
        "resources": sorted(
            [{field: row.get(field, "") for field in resource_fields} for row in resources],
            key=lambda item: json.dumps(item, sort_keys=True),
        ),
        "cohort_usage": sorted(
            [{field: row.get(field, "") for field in usage_fields} for row in usage],
            key=lambda item: json.dumps(item, sort_keys=True),
        ),
        "input_manifest_sha256": manifest_sha,
    }


def _derive_dataset_key_for_record(
    repo_root: Path,
    record: WorkspaceRecord,
    route_id: str,
    input_manifest: str,
) -> str:
    payload = _dataset_payload(repo_root, record, route_id, input_manifest)
    if not payload["candidates"] and not payload["resources"] and not payload[
        "cohort_usage"
    ] and not payload["input_manifest_sha256"]:
        return ""
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def derive_dataset_key(
    repo_root: Path,
    project_id: str,
    route_id: str,
    input_manifest: str,
) -> str:
    workspaces = _workspace_index(Path(repo_root))
    record = workspaces.get(project_id)
    if record is None:
        raise ValueError(f"unknown workspace project_id {project_id}")
    key = _derive_dataset_key_for_record(Path(repo_root), record, route_id, input_manifest)
    if not key:
        raise ValueError(f"no dataset evidence for {project_id}:{route_id}")
    return key


def _qualified_license(value: str) -> bool:
    normalized = " ".join(value.split()).casefold()
    disallowed = (
        "not selected",
        "unknown",
        "unlicensed",
        "no license",
        "all rights reserved",
        "internal project use only",
    )
    return bool(normalized) and not any(marker in normalized for marker in disallowed)


def _qualifying_case(row: dict[str, str], *, level: str | None = None) -> bool:
    if row.get("lifecycle") != "active" or row.get("last_outcome") != "passed":
        return False
    if row.get("oracle_source") in NON_QUALIFYING_ORACLES:
        return False
    return level is None or row.get("test_level") == level


def _independent_pair(
    references: Iterable[dict[str, str]],
    transfers: Iterable[dict[str, str]],
) -> bool:
    for reference in references:
        for transfer in transfers:
            if (
                reference.get("pilot_id")
                and transfer.get("pilot_id")
                and reference.get("pilot_id") != transfer.get("pilot_id")
                and reference.get("anchor_key")
                and transfer.get("anchor_key")
                and reference.get("anchor_key") != transfer.get("anchor_key")
                and reference.get("dataset_key")
                and transfer.get("dataset_key")
                and reference.get("dataset_key") != transfer.get("dataset_key")
            ):
                return True
    return False


def _validate_contract(
    repo_root: Path,
    release: dict[str, str],
    report: RegistryValidationReport,
) -> set[str]:
    release_id = release.get("release_id", "<unknown>")
    contract_value = release.get("contract_path", "")
    if not contract_value or _unsafe_relative_path(contract_value):
        return set()
    contract_path = repo_root / contract_value
    if not contract_path.exists():
        return set()
    header, rows = _read_tsv(contract_path)
    if header != list(CONTRACT_HEADER):
        report.errors.append(
            f"module {release_id}: contract header differs from the required schema"
        )
        return set()
    seen: set[str] = set()
    kinds: set[str] = set()
    required_test_ids: set[str] = set()
    for number, row in enumerate(rows, start=2):
        if row.get("__extra_columns__"):
            report.errors.append(
                f"module {release_id} contract:{number}: row has surplus columns"
            )
        item_id = row.get("item_id", "")
        if not item_id:
            report.errors.append(f"module {release_id} contract:{number}: empty item_id")
        elif item_id in seen:
            report.errors.append(
                f"module {release_id} contract:{number}: duplicate item_id {item_id}"
            )
        seen.add(item_id)
        kind = row.get("kind", "")
        if kind not in CONTRACT_KINDS:
            report.errors.append(
                f"module {release_id} contract:{number}: invalid kind {kind!r}"
            )
        kinds.add(kind)
        if row.get("required") not in BOOLEAN_VALUES:
            report.errors.append(
                f"module {release_id} contract:{number}: required must be true or false"
            )
        if row.get("failure_action") not in CONTRACT_ACTIONS:
            report.errors.append(
                f"module {release_id} contract:{number}: invalid failure_action"
            )
        for field_name in ("name", "specification", "failure_action"):
            if not row.get(field_name):
                report.errors.append(
                    f"module {release_id} contract:{number}: empty {field_name}"
                )
        test_ids = _split_ids(row.get("test_id", ""))
        if row.get("required") == "true" and not test_ids:
            report.errors.append(
                f"module {release_id} contract:{number}: required item lacks test_id"
            )
        if row.get("required") == "true":
            required_test_ids.update(test_ids)
    for required_kind in ("input", "output", "invariant"):
        if required_kind not in kinds:
            report.errors.append(
                f"module {release_id}: contract lacks a {required_kind} item"
            )
    return required_test_ids


def _validate_workspace_links(
    repo_root: Path,
    workspaces: dict[str, WorkspaceRecord],
    capabilities: dict[str, dict[str, str]],
    releases: dict[str, dict[str, str]],
    promotions: dict[str, dict[str, str]],
    report: RegistryValidationReport,
) -> dict[str, tuple[WorkspaceRecord, dict[str, str]]]:
    qualified_issues: dict[str, tuple[WorkspaceRecord, dict[str, str]]] = {}
    for project_id, record in workspaces.items():
        routes = {
            row.get("route_id", ""): row
            for row in _workspace_table(record, "evidence/routes.tsv")
            if row.get("route_id")
        }
        requirements = {
            row.get("module_id", ""): row
            for row in _workspace_table(record, "evidence/code_requirements.tsv")
            if row.get("module_id")
        }
        candidates = _workspace_table(record, "evidence/code_candidates.tsv")
        runs = _workspace_table(record, "execution/runs.tsv")
        issues = _workspace_table(record, "evidence/issues.tsv")

        route_capabilities: dict[str, set[str]] = {}
        for route_id, route in routes.items():
            capability_ids = set(_split_ids(route.get("capability_ids", "")))
            route_capabilities[route_id] = capability_ids
            for capability_id in capability_ids:
                capability = capabilities.get(capability_id)
                if capability is None:
                    report.errors.append(
                        f"{project_id}:{route_id}: unknown capability {capability_id}"
                    )
                elif (
                    route.get("decision_status") in {"active", "promoted"}
                    and capability.get("lifecycle") == "deprecated"
                ):
                    report.errors.append(
                        f"{project_id}:{route_id}: active route uses deprecated "
                        f"capability {capability_id}"
                    )

        for module_id, requirement in requirements.items():
            route_id = requirement.get("route_id", "")
            capability_id = requirement.get("capability_id", "")
            if capability_id not in capabilities:
                report.errors.append(
                    f"{project_id}:{module_id}: unknown capability {capability_id!r}"
                )
            if capability_id not in route_capabilities.get(route_id, set()):
                report.errors.append(
                    f"{project_id}:{module_id}: capability {capability_id!r} is not "
                    f"declared by route {route_id}"
                )

        release_candidates_by_route: dict[str, set[str]] = {}
        for candidate in candidates:
            release_id = candidate.get("reusable_release_id", "")
            if not release_id:
                continue
            code_id = candidate.get("code_id", "")
            release = releases.get(release_id)
            if release is None:
                report.errors.append(
                    f"{project_id}:{code_id}: unknown reusable release {release_id}"
                )
                continue
            requirement = requirements.get(candidate.get("module_id", ""), {})
            if requirement and release.get("capability_id") != requirement.get(
                "capability_id"
            ):
                report.errors.append(
                    f"{project_id}:{code_id}: release capability does not match "
                    "the code requirement"
                )
            if candidate.get("decision") == "use" and release.get("lifecycle") != "active":
                report.errors.append(
                    f"{project_id}:{code_id}: selected code uses non-active release"
                )
            if candidate.get("decision") == "use":
                release_candidates_by_route.setdefault(
                    candidate.get("route_id", ""), set()
                ).add(release_id)

        for run in runs:
            run_id = run.get("run_id", "")
            route_id = run.get("route_id", "")
            release_ids = set(_split_ids(run.get("module_release_ids", "")))
            for release_id in release_ids:
                release = releases.get(release_id)
                if release is None:
                    report.errors.append(
                        f"{project_id}:{run_id}: unknown module release {release_id}"
                    )
                elif release.get("capability_id") not in route_capabilities.get(
                    route_id, set()
                ):
                    report.errors.append(
                        f"{project_id}:{run_id}: release capability is not declared "
                        f"by route {route_id}"
                    )
            expected = release_candidates_by_route.get(route_id, set())
            if run.get("status") == "passed" and not expected.issubset(release_ids):
                missing = sorted(expected - release_ids)
                report.errors.append(
                    f"{project_id}:{run_id}: passed run omits selected releases "
                    f"{';'.join(missing)}"
                )

        for issue in issues:
            issue_id = issue.get("issue_id", "")
            if not issue_id:
                continue
            qualified_issues[f"{project_id}#{issue_id}"] = (record, issue)
            affected = set(_split_ids(issue.get("affected_capability_ids", "")))
            for capability_id in affected:
                if capability_id not in capabilities:
                    report.errors.append(
                        f"{project_id}#{issue_id}: unknown affected capability "
                        f"{capability_id}"
                    )
            promotion_id = issue.get("promotion_id", "")
            disposition = issue.get("disposition", "")
            if disposition in {"promote_to_module", "promote_to_core"} and not promotion_id:
                report.errors.append(
                    f"{project_id}#{issue_id}: promotion disposition lacks promotion_id"
                )
            if not promotion_id:
                continue
            promotion = promotions.get(promotion_id)
            if promotion is None:
                report.errors.append(
                    f"{project_id}#{issue_id}: unknown promotion {promotion_id}"
                )
                continue
            if issue.get("candidate_scope") != promotion.get("target_scope"):
                report.errors.append(
                    f"{project_id}#{issue_id}: issue scope differs from promotion target"
                )
            promotion_capabilities = set(
                _split_ids(promotion.get("affected_capability_ids", ""))
            )
            if not affected or not affected.issubset(promotion_capabilities):
                report.errors.append(
                    f"{project_id}#{issue_id}: issue capabilities do not match promotion"
                )
    return qualified_issues


def validate_registries(repo_root: Path) -> RegistryValidationReport:
    """Validate registries and their bidirectional links to schema-v2 workspaces."""

    root = Path(repo_root)
    report = RegistryValidationReport()
    tables: dict[str, list[dict[str, str]]] = {}
    for name, relative in REGISTRY_FILES.items():
        path = root / relative
        if not path.exists():
            report.errors.append(f"missing registry: {relative}")
            tables[name] = []
            continue
        header, rows = _read_tsv(path)
        if header != list(REGISTRY_HEADERS[name]):
            report.errors.append(
                f"{relative}: header differs from the required registry schema"
            )
        tables[name] = rows
        for number, row in enumerate(rows, start=2):
            for (enum_table, field_name), allowed in ENUMS.items():
                if enum_table != name:
                    continue
                value = row.get(field_name, "")
                if value not in allowed:
                    report.errors.append(
                        f"{relative}:{number}: invalid {field_name} {value!r}"
                    )

    capabilities = _index_rows(
        "capabilities", tables.get("capabilities", []), "capability_id", report
    )
    releases = _index_rows(
        "module_releases", tables.get("module_releases", []), "release_id", report
    )
    promotions = _index_rows(
        "promotions", tables.get("promotions", []), "promotion_id", report
    )
    cases = _index_rows(
        "regression_cases", tables.get("regression_cases", []), "case_id", report
    )
    workspaces = _workspace_index(root, report)

    for capability_id, row in capabilities.items():
        for field_name in (
            "family",
            "name",
            "input_modalities",
            "biological_unit",
            "output_kind",
            "claim_ceiling",
            "out_of_scope",
        ):
            if not row.get(field_name):
                report.errors.append(f"capability {capability_id}: empty {field_name}")

    required_tests_by_release: dict[str, set[str]] = {}
    for release_id, row in releases.items():
        capability_id = row.get("capability_id", "")
        if capability_id not in capabilities:
            report.errors.append(
                f"module {release_id}: unknown capability {capability_id!r}"
            )
        for field_name in (
            "module_id",
            "name",
            "version",
            "language",
            "license",
            "unit_test_command",
            "unit_test_path",
            "source_provenance",
        ):
            if not row.get(field_name):
                report.errors.append(f"module {release_id}: empty {field_name}")
        if row.get("version") and not SEMVER.fullmatch(row.get("version", "")):
            report.errors.append(f"module {release_id}: version is not semantic versioning")
        if not _qualified_license(row.get("license", "")):
            report.errors.append(f"module {release_id}: license is not qualified")
        for field_name in (
            "path", "entrypoint", "environment_lock", "contract_path",
            "unit_test_path",
        ):
            _check_repo_path(
                root,
                row.get(field_name, ""),
                f"module {release_id}.{field_name}",
                report,
            )
        for path_field, hash_field in (
            ("entrypoint", "implementation_sha256"),
            ("environment_lock", "environment_sha256"),
            ("contract_path", "contract_sha256"),
            ("unit_test_path", "unit_test_sha256"),
        ):
            expected_hash = row.get(hash_field, "")
            if not SHA256.fullmatch(expected_hash):
                report.errors.append(
                    f"module {release_id}: {hash_field} must be a 64-character SHA256"
                )
                continue
            path_value = row.get(path_field, "")
            if path_value and not _unsafe_relative_path(path_value):
                bound_path = root / path_value
                if bound_path.exists() and bound_path.is_file():
                    actual_hash = _sha256_canonical_text(bound_path)
                    if actual_hash.casefold() != expected_hash.casefold():
                        report.errors.append(
                            f"module {release_id}: {hash_field} is stale for {path_field}; "
                            "release or version must be updated"
                        )
        unit_test_path = row.get("unit_test_path", "").replace("\\", "/")
        try:
            unit_tokens = [
                token.replace("\\", "/")
                for token in shlex.split(row.get("unit_test_command", ""))
            ]
        except ValueError as exc:
            report.errors.append(
                f"module {release_id}: unit_test_command cannot be parsed: {exc}"
            )
            unit_tokens = []
        if unit_test_path and unit_test_path not in unit_tokens:
            report.errors.append(
                f"module {release_id}: unit_test_command does not reference "
                "unit_test_path"
            )
        module_path = row.get("path", "").replace("\\", "/").rstrip("/")
        contract_path = row.get("contract_path", "").replace("\\", "/")
        if module_path and contract_path and not (
            contract_path == module_path or contract_path.startswith(module_path + "/")
        ):
            report.errors.append(
                f"module {release_id}: contract_path must be inside the module path"
            )
        supersedes = row.get("supersedes", "")
        if supersedes:
            if supersedes == release_id:
                report.errors.append(f"module {release_id}: cannot supersede itself")
            elif supersedes not in releases:
                report.errors.append(
                    f"module {release_id}: supersedes unknown release {supersedes}"
                )
        required_tests_by_release[release_id] = _validate_contract(root, row, report)

    workspace_routes: dict[str, dict[str, dict[str, str]]] = {}
    workspace_runs: dict[str, dict[str, dict[str, str]]] = {}
    for project_id, record in workspaces.items():
        workspace_routes[project_id] = {
            row.get("route_id", ""): row
            for row in _workspace_table(record, "evidence/routes.tsv")
            if row.get("route_id")
        }
        workspace_runs[project_id] = {
            row.get("run_id", ""): row
            for row in _workspace_table(record, "execution/runs.tsv")
            if row.get("run_id")
        }

    for case_id, row in cases.items():
        capability_id = row.get("capability_id", "")
        if capability_id not in capabilities:
            report.errors.append(
                f"regression case {case_id}: unknown capability {capability_id!r}"
            )
        target_scope = row.get("target_scope", "")
        target_id = row.get("target_id", "")
        if not target_id:
            report.errors.append(f"regression case {case_id}: empty target_id")
        elif target_scope == "capability":
            if target_id not in capabilities:
                report.errors.append(
                    f"regression case {case_id}: unknown target capability {target_id}"
                )
            elif target_id != capability_id:
                report.errors.append(
                    f"regression case {case_id}: target and recorded capability differ"
                )
        elif target_scope == "module":
            release = releases.get(target_id)
            if release is None:
                report.errors.append(
                    f"regression case {case_id}: unknown target module release {target_id}"
                )
            elif release.get("capability_id") != capability_id:
                report.errors.append(
                    f"regression case {case_id}: module capability does not match case"
                )
        elif target_scope == "core" and not target_id.startswith("CORE-"):
            report.errors.append(
                f"regression case {case_id}: core target must use a CORE- identifier"
            )
        for field_name in ("command", "oracle_ref"):
            if not row.get(field_name):
                report.errors.append(f"regression case {case_id}: empty {field_name}")
        _check_repo_path(
            root,
            row.get("input_manifest", ""),
            f"regression case {case_id}.input_manifest",
            report,
            required=row.get("test_level") == "real_data",
        )
        _check_repo_path(
            root,
            row.get("assertion_path", ""),
            f"regression case {case_id}.assertion_path",
            report,
        )

        assertion_value = row.get("assertion_path", "")
        assertion_file = root / assertion_value if assertion_value else None
        test_inventory: set[str] = set()
        class_inventory: dict[str, set[str]] = {}
        if (
            assertion_file is not None
            and not _unsafe_relative_path(assertion_value)
            and assertion_file.exists()
            and assertion_file.suffix.casefold() == ".py"
        ):
            test_inventory, class_inventory = _python_test_inventory(
                assertion_file, f"regression case {case_id}", report
            )
            for test_id in _split_ids(row.get("covered_test_ids", "")):
                if test_id not in test_inventory:
                    report.errors.append(
                        f"regression case {case_id}: covered_test_ids names unknown "
                        f"test_id {test_id!r}"
                    )
        elif row.get("covered_test_ids"):
            report.errors.append(
                f"regression case {case_id}: covered_test_ids require a parseable "
                "Python assertion_path"
            )

        if row.get("test_level") == "unit" and row.get("last_outcome") == "passed":
            if not test_inventory:
                report.errors.append(
                    f"regression case {case_id}: passed unit case has no discoverable tests"
                )
            else:
                _validate_passed_unit_command(
                    row,
                    assertion_value.replace("\\", "/"),
                    test_inventory,
                    class_inventory,
                    report,
                )

        project_id = row.get("pilot_id", "")
        record = workspaces.get(project_id)
        route_id = row.get("route_id", "")
        if record is None:
            report.errors.append(
                f"regression case {case_id}: unknown workspace project_id {project_id!r}"
            )
        elif route_id not in workspace_routes.get(project_id, {}):
            report.errors.append(
                f"regression case {case_id}: unknown route {project_id}:{route_id}"
            )

        last_outcome = row.get("last_outcome", "")
        last_run_id = row.get("last_run_id", "")
        if last_outcome == "not_run" and last_run_id:
            report.errors.append(
                f"regression case {case_id}: not_run must not record last_run_id"
            )
        if last_outcome != "not_run" and not last_run_id:
            report.errors.append(
                f"regression case {case_id}: recorded outcome lacks last_run_id"
            )
        run = workspace_runs.get(project_id, {}).get(last_run_id) if last_run_id else None
        if last_run_id and run is None:
            report.errors.append(
                f"regression case {case_id}: unknown run {project_id}:{last_run_id}"
            )
        elif run is not None:
            if run.get("route_id") != route_id:
                report.errors.append(
                    f"regression case {case_id}: run belongs to another route"
                )
            expected_run_kind = {
                "real_data": "real_data",
                "unit": "unit",
            }.get(row.get("test_level", ""))
            if expected_run_kind and run.get("run_kind") != expected_run_kind:
                report.errors.append(
                    f"regression case {case_id}: {row.get('test_level')} case "
                    f"must bind a {expected_run_kind} run, not "
                    f"{run.get('run_kind')!r}"
                )
            expected_status = {"passed": "passed", "failed": "failed"}.get(last_outcome)
            if expected_status and run.get("status") != expected_status:
                report.errors.append(
                    f"regression case {case_id}: {last_outcome} case points to "
                    f"run status {run.get('status')!r}"
                )
            if row.get("test_level") == "real_data" and row.get(
                "input_manifest"
            ) != run.get("data_manifest"):
                report.errors.append(
                    f"regression case {case_id}: input_manifest differs from its run"
                )
            if row.get("test_level") == "real_data" and row.get(
                "command"
            ) != run.get("command"):
                report.errors.append(
                    f"regression case {case_id}: command differs from its bound run"
                )

        if row.get("test_level") == "real_data" and record is not None:
            manifest_anchor = record.manifest.get("anchor", {})
            doi = manifest_anchor.get("doi", "") if isinstance(manifest_anchor, dict) else ""
            expected_anchor = normalize_doi(str(doi))
            if not expected_anchor or normalize_doi(row.get("anchor_key", "")) != expected_anchor:
                report.errors.append(
                    f"regression case {case_id}: anchor_key must match PROJECT.json DOI"
                )
            expected_dataset = _derive_dataset_key_for_record(
                root, record, route_id, row.get("input_manifest", "")
            )
            if not expected_dataset:
                report.errors.append(
                    f"regression case {case_id}: no auditable dataset fingerprint"
                )
            elif row.get("dataset_key") != expected_dataset:
                report.errors.append(
                    f"regression case {case_id}: dataset_key does not match workspace "
                    "data/resource/cohort evidence"
                )
        elif row.get("test_level") != "real_data" and (
            row.get("anchor_key") or row.get("dataset_key")
        ):
            report.errors.append(
                f"regression case {case_id}: non-real-data case must not claim "
                "anchor_key or dataset_key"
            )

    active_cases = [row for row in cases.values() if row.get("lifecycle") == "active"]
    for case in active_cases:
        if case.get("case_role") != "transfer":
            continue
        same_target_refs = [
            row
            for row in active_cases
            if row.get("case_role") == "reference"
            and row.get("test_level") == "real_data"
            and row.get("target_scope") == case.get("target_scope")
            and row.get("target_id") == case.get("target_id")
        ]
        if not _independent_pair(same_target_refs, [case]):
            report.errors.append(
                f"regression case {case.get('case_id')}: transfer is not independent "
                "of a reference in workspace, anchor and derived dataset"
            )

    for release_id, release in releases.items():
        maturity = release.get("maturity")
        release_cases = [
            row
            for row in active_cases
            if row.get("target_scope") == "module" and row.get("target_id") == release_id
        ]
        qualifying = [row for row in release_cases if _qualifying_case(row)]
        unit_cases = [row for row in qualifying if row.get("test_level") == "unit"]
        references = [
            row
            for row in qualifying
            if row.get("case_role") == "reference"
            and row.get("test_level") == "real_data"
        ]
        transfers = [
            row
            for row in qualifying
            if row.get("case_role") == "transfer"
            and row.get("test_level") == "real_data"
        ]
        if maturity in {"unit_verified", "reference_verified", "transfer_verified"}:
            if not unit_cases:
                report.errors.append(
                    f"module {release_id}: {maturity} lacks a passed qualifying unit case"
                )
            covered = {
                test_id for row in qualifying for test_id in _split_ids(row.get("covered_test_ids", ""))
            }
            missing_tests = sorted(required_tests_by_release.get(release_id, set()) - covered)
            if missing_tests:
                report.errors.append(
                    f"module {release_id}: required contract tests lack passed case "
                    f"coverage: {';'.join(missing_tests)}"
                )
        if maturity in {"reference_verified", "transfer_verified"} and not references:
            report.errors.append(
                f"module {release_id}: {maturity} lacks a passed real-data reference "
                "with an oracle stronger than reviewed_freeze"
            )
        if maturity == "transfer_verified" and not _independent_pair(references, transfers):
            report.errors.append(
                f"module {release_id}: transfer_verified requires different workspace, "
                "PROJECT DOI and derived dataset"
            )

    qualified_issues = _validate_workspace_links(
        root, workspaces, capabilities, releases, promotions, report
    )

    for promotion_id, promotion in promotions.items():
        for field_name in (
            "source_refs",
            "target_id",
            "affected_capability_ids",
            "change_statement",
            "evidence_refs",
            "counterexample_case_ids",
        ):
            if not promotion.get(field_name):
                report.errors.append(f"promotion {promotion_id}: empty {field_name}")
        target_scope = promotion.get("target_scope")
        target_id = promotion.get("target_id", "")
        if target_scope == "module" and target_id not in releases:
            report.errors.append(
                f"promotion {promotion_id}: unknown target module release {target_id!r}"
            )
        if target_scope == "core" and not target_id.startswith("CORE-"):
            report.errors.append(
                f"promotion {promotion_id}: core target must use a CORE- identifier"
            )
        affected = set(_split_ids(promotion.get("affected_capability_ids", "")))
        for capability_id in affected:
            if capability_id not in capabilities:
                report.errors.append(
                    f"promotion {promotion_id}: unknown affected capability {capability_id}"
                )
        if target_scope == "module" and target_id in releases:
            target_capability = releases[target_id].get("capability_id")
            if affected != {target_capability}:
                report.errors.append(
                    f"promotion {promotion_id}: module promotion must declare exactly "
                    "its release capability"
                )

        source_refs = _split_ids(promotion.get("source_refs", ""))
        for source_ref in source_refs:
            match = SOURCE_REF.fullmatch(source_ref)
            if match is None:
                report.errors.append(
                    f"promotion {promotion_id}: malformed source ref {source_ref!r}"
                )
                continue
            source = qualified_issues.get(source_ref)
            if source is None:
                report.errors.append(
                    f"promotion {promotion_id}: unresolved source issue {source_ref}"
                )
                continue
            _, issue = source
            if issue.get("promotion_id") != promotion_id:
                report.errors.append(
                    f"promotion {promotion_id}: source issue {source_ref} lacks reciprocal link"
                )
            if issue.get("candidate_scope") != target_scope:
                report.errors.append(
                    f"promotion {promotion_id}: source issue {source_ref} has another scope"
                )
            issue_capabilities = set(
                _split_ids(issue.get("affected_capability_ids", ""))
            )
            if not issue_capabilities or not issue_capabilities.issubset(affected):
                report.errors.append(
                    f"promotion {promotion_id}: source issue {source_ref} has unrelated capabilities"
                )

        required_ids = _split_ids(promotion.get("required_case_ids", ""))
        evidence_ids = _split_ids(promotion.get("evidence_refs", ""))
        counter_ids = _split_ids(promotion.get("counterexample_case_ids", ""))
        for label, identifiers in (
            ("required", required_ids),
            ("evidence", evidence_ids),
            ("counterexample", counter_ids),
        ):
            missing = [case_id for case_id in identifiers if case_id not in cases]
            if missing:
                report.errors.append(
                    f"promotion {promotion_id}: unknown {label} cases {', '.join(missing)}"
                )
        selected = [cases[case_id] for case_id in required_ids if case_id in cases]
        counterexamples = [cases[case_id] for case_id in counter_ids if case_id in cases]
        for case in selected + counterexamples:
            if case.get("capability_id") not in affected:
                report.errors.append(
                    f"promotion {promotion_id}: case {case.get('case_id')} has unrelated capability"
                )
            if target_scope == "module" and not (
                case.get("target_scope") == "module" and case.get("target_id") == target_id
            ):
                report.errors.append(
                    f"promotion {promotion_id}: module case {case.get('case_id')} "
                    "does not target its release"
                )
        for case in counterexamples:
            if case.get("case_role") not in {"adversarial", "safety"}:
                report.errors.append(
                    f"promotion {promotion_id}: counterexample {case.get('case_id')} "
                    "is not adversarial or safety"
                )
            if case.get("oracle_source") not in {"failure_expected", "property"}:
                report.errors.append(
                    f"promotion {promotion_id}: counterexample {case.get('case_id')} "
                    "has an unsuitable oracle"
                )

        maturity = promotion.get("maturity")
        if maturity in {"provisional", "confirmed"}:
            qualifying = [row for row in selected if _qualifying_case(row)]
            real_references = [
                row
                for row in qualifying
                if row.get("case_role") == "reference"
                and row.get("test_level") == "real_data"
            ]
            if not real_references:
                report.errors.append(
                    f"promotion {promotion_id}: {maturity} lacks a qualifying "
                    "real-data reference"
                )
            if not counterexamples or not all(_qualifying_case(row) for row in counterexamples):
                report.errors.append(
                    f"promotion {promotion_id}: {maturity} requires passed registered "
                    "adversarial cases"
                )
        if maturity == "confirmed":
            qualifying = [row for row in selected if _qualifying_case(row)]
            for capability_id in affected:
                references = [
                    row
                    for row in qualifying
                    if row.get("capability_id") == capability_id
                    and row.get("case_role") == "reference"
                    and row.get("test_level") == "real_data"
                ]
                transfers = [
                    row
                    for row in qualifying
                    if row.get("capability_id") == capability_id
                    and row.get("case_role") == "transfer"
                    and row.get("test_level") == "real_data"
                ]
                if not _independent_pair(references, transfers):
                    report.errors.append(
                        f"promotion {promotion_id}: confirmed capability {capability_id} "
                        "requires independent reference and transfer cases"
                    )
            if promotion.get("review_status") != "approved" or promotion.get(
                "authority"
            ) not in {"user", "collaborator"}:
                report.errors.append(
                    f"promotion {promotion_id}: confirmed requires approved human authority"
                )
        if maturity in {"provisional", "confirmed"} and not promotion.get("reviewer"):
            report.errors.append(f"promotion {promotion_id}: maturity requires reviewer")
        if maturity == "rejected" and promotion.get("review_status") != "rejected":
            report.errors.append(
                f"promotion {promotion_id}: rejected maturity requires rejected review_status"
            )

    return report


def capability_coverage(repo_root: Path) -> dict[str, str]:
    """Derive coverage only after strict repository and registry validation."""

    report = validate_registries(repo_root)
    if not report.ok:
        raise ValueError("invalid registries: " + "; ".join(report.errors[:8]))
    tables = load_registries(repo_root)
    coverage: dict[str, str] = {
        row.get("capability_id", ""): "uncovered"
        for row in tables.get("capabilities", [])
        if row.get("capability_id")
    }
    qualifying = [
        row
        for row in tables.get("regression_cases", [])
        if _qualifying_case(row, level="real_data")
    ]
    for capability_id in coverage:
        rows = [row for row in qualifying if row.get("capability_id") == capability_id]
        references = [row for row in rows if row.get("case_role") == "reference"]
        transfers = [row for row in rows if row.get("case_role") == "transfer"]
        if _independent_pair(references, transfers):
            coverage[capability_id] = "transfer_exercised"
        elif references:
            coverage[capability_id] = "reference_exercised"
    return coverage
