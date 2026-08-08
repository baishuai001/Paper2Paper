from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from .workspace import (
    init_workspace,
    load_manifest,
    next_actions,
    promote_workspace,
    status_summary,
    validate_workspace,
    write_readiness_report,
)


SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@-]*$")


def _load_repository_workspaces(repo_root: Path) -> dict[str, Path]:
    """Index the repository's direct Pilot and manuscript workspaces.

    This is deliberately local and small. Paper2Paper needs to prevent an
    accidental duplicate project ID or an out-of-repository promotion; it does
    not need a cross-paper platform to do so.
    """

    root = Path(repo_root).resolve()
    index: dict[str, Path] = {}
    for parent_name in ("pilots", "manuscript-projects"):
        parent = root / parent_name
        if not parent.is_dir():
            continue
        for manifest_path in sorted(parent.glob("*/PROJECT.json")):
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            project_id = manifest.get("project_id")
            if not isinstance(project_id, str) or not SAFE_ID.fullmatch(project_id):
                raise ValueError(
                    f"invalid project_id in {manifest_path.relative_to(root)}"
                )
            workspace = manifest_path.parent.resolve()
            if project_id in index:
                raise ValueError(
                    f"duplicate project_id in repository workspaces: {project_id}"
                )
            index[project_id] = workspace
    return index


def _find_repo_root(path: Path) -> Path | None:
    start = Path(path).resolve()
    if start.is_file():
        start = start.parent
    for candidate in (start, *start.parents):
        if (
            (candidate / "pyproject.toml").is_file()
            and (candidate / "pilots").is_dir()
            and (candidate / "manuscript-projects").is_dir()
        ):
            return candidate
    return None


def _require_repository_pilot(repo_root: Path, pilot_dir: Path) -> None:
    """Reject promotion when the Pilot is not owned by this repository."""

    root = Path(repo_root).resolve()
    pilot = Path(pilot_dir).resolve()
    pilots_root = (root / "pilots").resolve()
    if pilot.parent != pilots_root:
        raise ValueError(
            "promotion Pilot must be one direct child of the selected "
            "repository's pilots directory"
        )
    manifest = load_manifest(pilot)
    project_id = manifest.get("project_id", "")
    record = _load_repository_workspaces(root).get(str(project_id))
    if record is None or record != pilot:
        raise ValueError(
            "promotion Pilot project_id and path do not match this repository"
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
    workspace_index = _load_repository_workspaces(root)
    if project_id in workspace_index:
        raise ValueError(
            f"promotion project_id already exists in this repository: {project_id}"
        )
    if target in workspace_index.values():
        raise ValueError(
            "promotion target is already registered by another project_id"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="paper2paper",
        description=(
            "Turn real anchor-paper frameworks into evidence-backed manuscript "
            "projects while recording workflow improvements discovered in use."
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
                    "promotion requires a Paper2Paper repository root"
                )
            _require_repository_pilot(repo_root, args.pilot_dir)
            _require_manuscript_target(
                repo_root, args.project_dir, args.project_id
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


    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
