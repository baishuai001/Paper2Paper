from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .registry import (
    SAFE_ID,
    capability_coverage,
    load_workspace_index,
    validate_registries,
)
from .regression import build_regression_plan
from .workspace import (
    init_workspace,
    load_manifest,
    next_actions,
    promote_workspace,
    status_summary,
    validate_workspace,
    write_readiness_report,
)


def _find_repo_root(path: Path) -> Path | None:
    start = Path(path).resolve()
    if start.is_file():
        start = start.parent
    for candidate in (start, *start.parents):
        if (candidate / "pyproject.toml").is_file() and (
            candidate / "registries"
        ).is_dir():
            return candidate
    return None


def _require_registered_pilot(repo_root: Path, pilot_dir: Path) -> None:
    """Reject promotion when the Pilot is not owned by this registry root."""

    root = Path(repo_root).resolve()
    pilot = Path(pilot_dir).resolve()
    pilots_root = (root / "pilots").resolve()
    if not pilot.is_relative_to(pilots_root):
        raise ValueError(
            "promotion Pilot must be inside the selected repository's pilots directory"
        )
    manifest = load_manifest(pilot)
    project_id = manifest.get("project_id", "")
    record = load_workspace_index(root).get(str(project_id))
    if record is None or record.path.resolve() != pilot:
        raise ValueError(
            "promotion Pilot is not registered by project_id and path in this repository"
        )


def _require_manuscript_target(
    repo_root: Path, project_dir: Path, project_id: str
) -> None:
    """Reject manuscript targets that would escape or corrupt the workspace index."""

    root = Path(repo_root).resolve()
    target = Path(project_dir).resolve()
    manuscript_root = (root / "manuscript-projects").resolve()
    if not isinstance(project_id, str) or not SAFE_ID.fullmatch(project_id):
        raise ValueError(
            "promotion project_id must start with an alphanumeric character and "
            "contain only letters, digits, dot, underscore, @ or hyphen"
        )
    if target.parent != manuscript_root:
        raise ValueError(
            "promotion target must be one direct child of this repository's "
            "manuscript-projects directory"
        )
    if target.exists():
        raise ValueError(
            "promotion target must not already exist; choose a new manuscript directory"
        )
    workspace_index = load_workspace_index(root)
    if project_id in workspace_index:
        raise ValueError(
            f"promotion project_id already exists in this repository: {project_id}"
        )
    if any(record.path.resolve() == target for record in workspace_index.values()):
        raise ValueError(
            "promotion target is already registered by another project_id"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="paper2paper",
        description=(
            "Audit real anchor papers, learn only scoped reusable controls, and "
            "promote user-approved routes into manuscript projects."
        ),
    )
    commands = parser.add_subparsers(dest="command", required=True)

    init_parser = commands.add_parser("init", help="create a Pilot workspace")
    init_parser.add_argument("project_dir", type=Path)
    init_parser.add_argument("--project-id", required=True)
    init_parser.add_argument("--title", required=True)
    init_parser.add_argument("--anchor-title", required=True)
    init_parser.add_argument("--doi", default="")

    promote_parser = commands.add_parser(
        "promote", help="create a manuscript project from a user-promoted pilot route"
    )
    promote_parser.add_argument("pilot_dir", type=Path)
    promote_parser.add_argument("route_id")
    promote_parser.add_argument("project_dir", type=Path)
    promote_parser.add_argument("--project-id", required=True)
    promote_parser.add_argument("--title", required=True)
    promote_parser.add_argument(
        "--repo-root",
        type=Path,
        help="Paper2Paper repository root; inferred from the Pilot when omitted",
    )

    validate_parser = commands.add_parser(
        "validate", help="validate structure and stage dependencies"
    )
    validate_parser.add_argument("project_dir", type=Path)
    validate_parser.add_argument("--json", action="store_true")

    status_parser = commands.add_parser(
        "status", help="show route readiness and validation state"
    )
    status_parser.add_argument("project_dir", type=Path)
    status_parser.add_argument("--json", action="store_true")

    next_parser = commands.add_parser(
        "next", help="show evidence-driven next actions"
    )
    next_parser.add_argument("project_dir", type=Path)
    next_parser.add_argument("--json", action="store_true")

    report_parser = commands.add_parser(
        "report", help="write reports/readiness.md"
    )
    report_parser.add_argument("project_dir", type=Path)

    registry_parser = commands.add_parser(
        "validate-registry", help="validate cross-paper capabilities and regression evidence"
    )
    registry_parser.add_argument("repo_root", type=Path)
    registry_parser.add_argument("--json", action="store_true")

    coverage_parser = commands.add_parser(
        "coverage", help="show evidence-derived capability coverage"
    )
    coverage_parser.add_argument("repo_root", type=Path)
    coverage_parser.add_argument("--json", action="store_true")

    regression_parser = commands.add_parser(
        "regression-plan", help="select scoped regression cases for changed paths"
    )
    regression_parser.add_argument("repo_root", type=Path)
    regression_parser.add_argument(
        "--changed-path", action="append", required=True,
        help="repository-relative changed path; repeat for multiple paths",
    )
    regression_parser.add_argument("--json", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "init":
            path = init_workspace(
                args.project_dir,
                args.project_id,
                args.title,
                args.anchor_title,
                args.doi,
            )
            print(f"initialized: {path}")
            return 0

        if args.command == "promote":
            repo_root = args.repo_root or _find_repo_root(args.pilot_dir)
            if repo_root is None:
                raise ValueError(
                    "promotion requires a Paper2Paper repository root with registries"
                )
            _require_registered_pilot(repo_root, args.pilot_dir)
            _require_manuscript_target(
                repo_root, args.project_dir, args.project_id
            )
            registry_report = validate_registries(repo_root)
            if not registry_report.ok:
                raise ValueError(
                    "repository registries and workspace links must validate before "
                    "promotion: " + "; ".join(registry_report.errors)
                )
            path = promote_workspace(
                args.pilot_dir,
                args.route_id,
                args.project_dir,
                args.project_id,
                args.title,
            )
            print(f"promoted manuscript project: {path}")
            return 0

        if args.command == "validate":
            report = validate_workspace(args.project_dir)
            if args.json:
                print(
                    json.dumps(
                        {
                            "ok": report.ok,
                            "errors": report.errors,
                            "warnings": report.warnings,
                        },
                        indent=2,
                        ensure_ascii=False,
                    )
                )
            else:
                for warning in report.warnings:
                    print(f"WARNING: {warning}")
                for error in report.errors:
                    print(f"ERROR: {error}")
                print(
                    "structural and contract validation passed"
                    if report.ok
                    else "structural and contract validation failed"
                )
            return 0 if report.ok else 1

        if args.command == "status":
            summary = status_summary(args.project_dir)
            if args.json:
                print(json.dumps(summary, indent=2, ensure_ascii=False))
            else:
                print(f"project: {summary['project_id']}")
                print(f"workspace kind: {summary['workspace_kind']}")
                print(f"stage: {summary['stage']}")
                print(
                    "selected route: "
                    f"{summary['selected_route_id'] or '(none)'}"
                )
                if not summary["routes"]:
                    print("routes: none")
                for route in summary["routes"]:
                    print(
                        f"{route['route_id']}\texecution_ready="
                        f"{str(route['execution_ready']).lower()}\t"
                        f"promotion_evidence_complete="
                        f"{str(route['promotion_evidence_complete']).lower()}\t"
                        f"execution_gaps={len(route['execution_gaps'])}"
                    )
                print(f"validation errors: {len(summary['validation_errors'])}")
            return 0 if not summary["validation_errors"] else 1

        if args.command == "next":
            actions = next_actions(args.project_dir)
            if args.json:
                print(json.dumps(actions, indent=2, ensure_ascii=False))
            else:
                for index, action in enumerate(actions, start=1):
                    print(f"{index}. {action}")
            return 0

        if args.command == "report":
            path = write_readiness_report(args.project_dir)
            print(f"wrote: {path}")
            return 0


        if args.command == "validate-registry":
            report = validate_registries(args.repo_root)
            payload = {
                "ok": report.ok,
                "errors": report.errors,
                "warnings": report.warnings,
            }
            if args.json:
                print(json.dumps(payload, indent=2, ensure_ascii=False))
            else:
                for warning in report.warnings:
                    print(f"WARNING: {warning}")
                for error in report.errors:
                    print(f"ERROR: {error}")
                print(
                    "cross-paper registry validation passed"
                    if report.ok
                    else "cross-paper registry validation failed"
                )
            return 0 if report.ok else 1

        if args.command == "coverage":
            coverage = capability_coverage(args.repo_root)
            if args.json:
                print(json.dumps(coverage, indent=2, ensure_ascii=False))
            else:
                for capability_id, status in sorted(coverage.items()):
                    print(f"{capability_id}\t{status}")
            return 0

        if args.command == "regression-plan":
            plan = build_regression_plan(args.repo_root, args.changed_path)
            payload = {
                "changed_paths": plan.changed_paths,
                "selected": [
                    {"case_id": item.case_id, "reasons": item.reasons}
                    for item in plan.selected
                ],
                "unmatched_paths": plan.unmatched_paths,
            }
            if args.json:
                print(json.dumps(payload, indent=2, ensure_ascii=False))
            else:
                for item in plan.selected:
                    print(f"{item.case_id}\t{' | '.join(item.reasons)}")
                for path in plan.unmatched_paths:
                    print(f"UNMATCHED\t{path}")
            return 0
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
