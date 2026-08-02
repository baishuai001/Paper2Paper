from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .core import (
    downstream_impact,
    init_project,
    route_priorities,
    status_summary,
    validate_project,
)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="paper2paper",
        description=(
            "Validate execution-first, auditable paper-adaptation projects."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser(
        "init", help="create a new intake-stage project"
    )
    init_parser.add_argument("project_dir", type=Path)
    init_parser.add_argument("--project-id", required=True)
    init_parser.add_argument("--title", required=True)

    validate_parser = subparsers.add_parser(
        "validate", help="validate a project instance"
    )
    validate_parser.add_argument("project_dir", type=Path)
    validate_parser.add_argument("--json", action="store_true")

    status_parser = subparsers.add_parser(
        "status", help="show project status and review queue"
    )
    status_parser.add_argument("project_dir", type=Path)
    status_parser.add_argument("--json", action="store_true")

    priority_parser = subparsers.add_parser(
        "priorities",
        help="show recorded and rule-derived execution priorities",
    )
    priority_parser.add_argument("project_dir", type=Path)
    priority_parser.add_argument("--json", action="store_true")

    impact_parser = subparsers.add_parser(
        "impact", help="show downstream entities affected by a change"
    )
    impact_parser.add_argument("project_dir", type=Path)
    impact_parser.add_argument("entity_id")
    impact_parser.add_argument("--json", action="store_true")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        if args.command == "init":
            target = init_project(args.project_dir, args.project_id, args.title)
            print(f"initialized intake project: {target}")
            return 0

        if args.command == "validate":
            report = validate_project(args.project_dir)
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
                    "validation passed"
                    if report.ok
                    else f"validation failed: {len(report.errors)} error(s)"
                )
            return 0 if report.ok else 1

        if args.command == "status":
            summary = status_summary(args.project_dir)
            if args.json:
                print(json.dumps(summary, indent=2, ensure_ascii=False))
            else:
                print(f"project: {summary['project_id']}")
                print(f"primary output: {summary['primary_output']}")
                print(f"status: {summary['project_status']}")
                print(f"gate: {summary['current_gate']}")
                print(
                    "active route: "
                    f"{summary['active_route_id'] or '(not selected)'}"
                )
                print(
                    "open work items: "
                    + (", ".join(summary["open_work_items"]) or "none")
                )
                print(
                    "pending reviews: "
                    + (", ".join(summary["pending_reviews"]) or "none")
                )
                print(
                    "pending decisions: "
                    + (", ".join(summary["pending_decisions"]) or "none")
                )
                print(
                    "open changes: "
                    + (", ".join(summary["open_change_requests"]) or "none")
                )
                print(f"validation errors: {len(summary['validation_errors'])}")
            return 0 if not summary["validation_errors"] else 1

        if args.command == "priorities":
            priorities = route_priorities(args.project_dir)
            if args.json:
                print(json.dumps(priorities, indent=2, ensure_ascii=False))
            else:
                if not priorities:
                    print("no route assessments")
                for item in priorities:
                    print(
                        f"{item['route_id']}\trecorded={item['recorded']}\t"
                        f"expected={item['expected']}\t"
                        f"match={str(item['matches']).lower()}"
                    )
            return 0

        if args.command == "impact":
            impacted = downstream_impact(args.project_dir, args.entity_id)
            if args.json:
                print(json.dumps(impacted, indent=2, ensure_ascii=False))
            else:
                if not impacted:
                    print(f"no active downstream dependencies: {args.entity_id}")
                for item in impacted:
                    print(
                        f"{item['entity_id']}\tdepth={item['depth']}\t"
                        f"via={item['via_relation']}"
                    )
            return 0

    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    parser.error("unknown command")
    return 2
