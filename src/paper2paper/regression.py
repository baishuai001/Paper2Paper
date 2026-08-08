from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath, PureWindowsPath

from .registry import load_registries, load_workspace_index, validate_registries


@dataclass
class RegressionPlanItem:
    case_id: str
    reasons: list[str] = field(default_factory=list)


@dataclass
class RegressionPlan:
    changed_paths: list[str]
    selected: list[RegressionPlanItem]
    unmatched_paths: list[str] = field(default_factory=list)

    @property
    def case_ids(self) -> list[str]:
        return [item.case_id for item in self.selected]


def _normalise_changed_path(value: str | Path) -> str:
    raw = str(value)
    windows = PureWindowsPath(raw)
    normalised = raw.replace("\\", "/")
    if normalised.startswith("./"):
        normalised = normalised[2:]
    posix = PurePosixPath(normalised)
    if (
        not normalised
        or windows.is_absolute()
        or bool(windows.drive)
        or posix.is_absolute()
        or ".." in posix.parts
    ):
        raise ValueError(f"changed path must be repository-relative: {raw}")
    return posix.as_posix()


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [
            {key: (value or "").strip() for key, value in row.items() if key}
            for row in csv.DictReader(handle, delimiter="\t")
        ]


def _split_ids(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def _workspace_bindings(workspace_path: Path) -> tuple[set[str], set[str]]:
    capabilities: set[str] = set()
    releases: set[str] = set()
    for row in _read_rows(workspace_path / "evidence/routes.tsv"):
        capabilities.update(_split_ids(row.get("capability_ids", "")))
    for row in _read_rows(workspace_path / "evidence/code_candidates.tsv"):
        if row.get("decision") == "use":
            releases.update(_split_ids(row.get("reusable_release_id", "")))
    return capabilities, releases


def _documentation_only(path: str) -> bool:
    lower = path.lower()
    name = PurePosixPath(lower).name
    return (
        lower.startswith("docs/")
        or name in {"readme.md", "license", "license.md"}
        or lower.endswith(".md")
    )


def _workspace_semantic_change(relative: str) -> bool:
    lower = relative.lower()
    if lower == "project.json":
        return True
    first = PurePosixPath(lower).parts[0] if PurePosixPath(lower).parts else ""
    if first in {"config", "code", "evidence", "execution"}:
        return True
    if first == "analysis":
        return True
    if first == "anchor" and lower.endswith((".md", ".tsv", ".json")):
        return True
    return False


def build_regression_plan(
    repo_root: Path,
    changed_paths: list[str | Path],
) -> RegressionPlan:
    """Return the smallest justified case set for repository-relative changes.

    Registry validation is a hard precondition: an invalid evidence graph cannot
    be used to manufacture a reassuring regression plan.
    """

    root = Path(repo_root)
    validation = validate_registries(root)
    if not validation.ok:
        details = "; ".join(validation.errors[:8])
        raise ValueError(f"cannot plan regression from invalid registries: {details}")

    paths = [_normalise_changed_path(value) for value in changed_paths]
    tables = load_registries(root)
    active_cases = {
        row.get("case_id", ""): row
        for row in tables.get("regression_cases", [])
        if row.get("case_id") and row.get("lifecycle") == "active"
    }
    releases = {
        row.get("release_id", ""): row
        for row in tables.get("module_releases", [])
        if row.get("release_id") and row.get("lifecycle") == "active"
    }
    promotions = tables.get("promotions", [])
    workspaces = load_workspace_index(root)
    workspace_bindings = {
        project_id: _workspace_bindings(record.path)
        for project_id, record in workspaces.items()
    }

    reasons_by_case: dict[str, list[str]] = {}
    matched_paths: set[str] = set()

    def select(case_id: str, reason: str) -> None:
        if case_id not in active_cases:
            return
        reasons = reasons_by_case.setdefault(case_id, [])
        if reason not in reasons:
            reasons.append(reason)

    def select_all(reason: str) -> None:
        for case_id in active_cases:
            select(case_id, reason)

    def select_release(release_id: str, changed: str) -> None:
        release = releases[release_id]
        capability_id = release.get("capability_id", "")
        for case_id, case in active_cases.items():
            if (
                case.get("target_scope") == "module"
                and case.get("target_id") == release_id
            ) or case.get("capability_id") == capability_id:
                select(
                    case_id,
                    f"{changed} changes release {release_id} or its capability contract",
                )

    for changed in paths:
        lower = changed.lower()

        if lower.startswith(".github/"):
            matched_paths.add(changed)
            select_all(f"{changed} changes shared CI or execution policy")
            continue

        affected_releases: set[str] = set()
        for release_id, release in releases.items():
            module_path = release.get("path", "").replace("\\", "/").rstrip("/")
            assets = {
                release.get("entrypoint", "").replace("\\", "/"),
                release.get("contract_path", "").replace("\\", "/"),
                release.get("environment_lock", "").replace("\\", "/"),
                release.get("unit_test_path", "").replace("\\", "/"),
            }
            if (
                (module_path and (changed == module_path or changed.startswith(module_path + "/")))
                or changed in assets
            ):
                affected_releases.add(release_id)
        if affected_releases:
            matched_paths.add(changed)
            for release_id in sorted(affected_releases):
                select_release(release_id, changed)
            continue

        workspace_match = None
        workspace_relative = ""
        for project_id, record in workspaces.items():
            prefix = record.relative_path.rstrip("/")
            if changed == prefix or changed.startswith(prefix + "/"):
                workspace_match = project_id
                workspace_relative = changed[len(prefix) :].lstrip("/")
                break
        if workspace_match is not None:
            matched_paths.add(changed)
            if _workspace_semantic_change(workspace_relative):
                capabilities, bound_releases = workspace_bindings[workspace_match]
                for case_id, case in active_cases.items():
                    if (
                        case.get("pilot_id") == workspace_match
                        or case.get("capability_id") in capabilities
                        or (
                            case.get("target_scope") == "module"
                            and case.get("target_id") in bound_releases
                        )
                    ):
                        select(
                            case_id,
                            f"{changed} changes semantic evidence or execution for "
                            f"workspace {workspace_match}",
                        )
            continue

        if lower == "registries/capabilities.tsv":
            matched_paths.add(changed)
            select_all(f"{changed} changes shared capability definitions")
            continue
        if lower == "registries/module_releases.tsv":
            matched_paths.add(changed)
            select_all(f"{changed} changes reusable release bindings")
            continue
        if lower == "registries/regression_cases.tsv":
            matched_paths.add(changed)
            select_all(f"{changed} changes the active regression inventory")
            continue
        if lower == "registries/promotions.tsv":
            matched_paths.add(changed)
            for promotion in promotions:
                selected_ids = _split_ids(promotion.get("required_case_ids", ""))
                selected_ids.update(
                    _split_ids(promotion.get("counterexample_case_ids", ""))
                )
                for case_id in selected_ids:
                    select(
                        case_id,
                        f"{changed} changes promotion evidence requiring {case_id}",
                    )
            continue

        if lower.startswith("src/paper2paper/") and lower.endswith(".py"):
            matched_paths.add(changed)
            select_all(f"{changed} changes shared execution or validation logic")
            continue

        if lower.startswith("tests/"):
            matched_paths.add(changed)
            continue
        if _documentation_only(changed):
            matched_paths.add(changed)
            continue

    selected = [
        RegressionPlanItem(case_id=case_id, reasons=reasons_by_case[case_id])
        for case_id in sorted(reasons_by_case)
    ]
    unmatched = [path for path in paths if path not in matched_paths]
    return RegressionPlan(
        changed_paths=paths,
        selected=selected,
        unmatched_paths=unmatched,
    )
