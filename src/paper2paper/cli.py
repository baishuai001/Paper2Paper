from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .workspace import (
    init_workspace,
    next_actions,
    status_summary,
    validate_workspace,
    write_readiness_report,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="paper2paper",
        description="Move an anchor-paper adaptation toward an executable manuscript.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    init_parser = commands.add_parser("init", help="create a paper workspace")
    init_parser.add_argument("project_dir", type=Path)
    init_parser.add_argument("--project-id", required=True)
    init_parser.add_argument("--title", required=True)
    init_parser.add_argument("--anchor-title", required=True)
    init_parser.add_argument("--doi", default="")

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
                        f"manuscript_eligible="
                        f"{str(route['manuscript_eligible']).lower()}\t"
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
