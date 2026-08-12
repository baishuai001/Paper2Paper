from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

from .schema import TABLES


SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._@-]*$")

FINDING_COLUMNS = (
    "finding_id",
    "title",
    "scope",
    "capability_area",
    "status",
    "source_refs",
    "generalized_failure",
    "risk_to_paper",
    "core_response",
    "verification_refs",
    "updated_at",
    "notes",
)

FINDING_SCOPES = ("paper_specific", "paper_type", "general")
FINDING_STATUSES = (
    "observed",
    "accepted",
    "implemented",
    "verified",
    "local_only",
    "rejected",
)

WORKFLOW_CANDIDATE_ACTIONS = {"consider", "deferred", "implemented"}


@dataclass
class LearningValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames or []
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    return header, rows


def _split_refs(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def _escape_markdown(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ").strip()


def _load_pilot_issues(
    repo_root: Path,
) -> tuple[list[dict[str, str]], LearningValidationReport, int]:
    root = Path(repo_root).resolve()
    report = LearningValidationReport()
    issues: list[dict[str, str]] = []
    project_ids: set[str] = set()
    pilots_scanned = 0
    pilots_root = root / "pilots"
    if not pilots_root.is_dir():
        report.errors.append("missing pilots directory")
        return issues, report, pilots_scanned

    expected_issue_header = list(TABLES["issues"]["columns"])
    for manifest_path in sorted(pilots_root.glob("*/PROJECT.json")):
        pilots_scanned += 1
        pilot_dir = manifest_path.parent
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            report.errors.append(
                f"cannot read {manifest_path.relative_to(root)}: {exc}"
            )
            continue

        project_id = manifest.get("project_id")
        if not isinstance(project_id, str) or not SAFE_ID.fullmatch(project_id):
            report.errors.append(
                f"invalid project_id in {manifest_path.relative_to(root)}"
            )
            continue
        if project_id in project_ids:
            report.errors.append(f"duplicate Pilot project_id: {project_id}")
            continue
        project_ids.add(project_id)

        issue_path = pilot_dir / TABLES["issues"]["path"]
        if not issue_path.is_file():
            report.errors.append(
                f"missing issue table: {issue_path.relative_to(root)}"
            )
            continue
        header, rows = _read_tsv(issue_path)
        if header != expected_issue_header:
            report.errors.append(
                f"{issue_path.relative_to(root)}: header must match issues contract"
            )
            continue
        seen_issue_ids: set[str] = set()
        for row_number, row in enumerate(rows, start=2):
            issue_id = row.get("issue_id", "")
            if not issue_id or not SAFE_ID.fullmatch(issue_id):
                report.errors.append(
                    f"{issue_path.relative_to(root)}:{row_number}: invalid issue_id"
                )
                continue
            if issue_id in seen_issue_ids:
                report.errors.append(
                    f"{issue_path.relative_to(root)}:{row_number}: duplicate issue_id={issue_id}"
                )
                continue
            seen_issue_ids.add(issue_id)
            issues.append(
                {
                    **row,
                    "project_id": project_id,
                    "pilot_path": pilot_dir.relative_to(root).as_posix(),
                    "source_ref": f"{project_id}#{issue_id}",
                }
            )
    return issues, report, pilots_scanned


def _is_workflow_candidate(issue: dict[str, str]) -> bool:
    # The source scope says where the fact was observed, not whether it is
    # reusable.  An anchor-only or scientific-result issue can still expose a
    # workflow problem.  workflow_action only nominates it for triage; it does
    # not make it a product gap.
    return issue.get("workflow_action") in WORKFLOW_CANDIDATE_ACTIONS


def validate_learning(repo_root: Path) -> LearningValidationReport:
    root = Path(repo_root).resolve()
    issues, report, _ = _load_pilot_issues(root)
    issue_by_ref = {row["source_ref"]: row for row in issues}
    candidate_refs = {
        row["source_ref"] for row in issues if _is_workflow_candidate(row)
    }

    finding_path = root / "learning" / "findings.tsv"
    if not finding_path.is_file():
        report.errors.append("missing learning/findings.tsv")
        return report
    header, findings = _read_tsv(finding_path)
    if header != list(FINDING_COLUMNS):
        report.errors.append(
            "learning/findings.tsv: header must exactly match the learning contract"
        )
        return report

    seen_finding_ids: set[str] = set()
    claimed_source_refs: dict[str, str] = {}
    required_fields = (
        "finding_id",
        "title",
        "scope",
        "capability_area",
        "status",
        "source_refs",
        "generalized_failure",
        "risk_to_paper",
        "core_response",
        "updated_at",
    )
    for row_number, row in enumerate(findings, start=2):
        prefix = f"learning/findings.tsv:{row_number}"
        for field_name in required_fields:
            if not row.get(field_name, ""):
                report.errors.append(f"{prefix}: empty {field_name}")
        finding_id = row.get("finding_id", "")
        if finding_id and not SAFE_ID.fullmatch(finding_id):
            report.errors.append(f"{prefix}: invalid finding_id={finding_id!r}")
        if finding_id in seen_finding_ids:
            report.errors.append(f"{prefix}: duplicate finding_id={finding_id}")
        seen_finding_ids.add(finding_id)

        scope = row.get("scope", "")
        status = row.get("status", "")
        if scope not in FINDING_SCOPES:
            report.errors.append(f"{prefix}: invalid scope={scope!r}")
        if status not in FINDING_STATUSES:
            report.errors.append(f"{prefix}: invalid status={status!r}")
        if scope == "paper_specific" and status not in {"local_only", "rejected"}:
            report.errors.append(
                f"{prefix}: paper_specific finding must be local_only or rejected"
            )
        if status == "local_only" and scope != "paper_specific":
            report.errors.append(
                f"{prefix}: local_only status requires scope=paper_specific"
            )
        if status == "verified" and not row.get("verification_refs", ""):
            report.errors.append(
                f"{prefix}: verified finding requires verification_refs"
            )
        try:
            date.fromisoformat(row.get("updated_at", ""))
        except ValueError:
            report.errors.append(
                f"{prefix}: updated_at must be an ISO date (YYYY-MM-DD)"
            )

        refs = _split_refs(row.get("source_refs", ""))
        if len(refs) != len(set(refs)):
            report.errors.append(f"{prefix}: duplicate source_refs in one finding")
        for source_ref in refs:
            if source_ref not in issue_by_ref:
                report.errors.append(
                    f"{prefix}: unknown source_ref={source_ref}"
                )
                continue
            if source_ref not in candidate_refs:
                report.errors.append(
                    f"{prefix}: source_ref is not marked as a workflow candidate: "
                    f"{source_ref}"
                )
            previous = claimed_source_refs.get(source_ref)
            if previous and previous != finding_id:
                report.errors.append(
                    f"{prefix}: source_ref={source_ref} is already classified by {previous}"
                )
            claimed_source_refs[source_ref] = finding_id

    return report


def learning_summary(repo_root: Path) -> dict[str, Any]:
    root = Path(repo_root).resolve()
    issues, issue_report, pilots_scanned = _load_pilot_issues(root)
    validation = validate_learning(root)
    # Avoid duplicating issue-table errors returned by both passes.
    validation_errors = list(dict.fromkeys(validation.errors))
    validation_warnings = list(dict.fromkeys(validation.warnings))
    if issue_report.errors:
        validation_errors = list(
            dict.fromkeys([*validation_errors, *issue_report.errors])
        )

    finding_path = root / "learning" / "findings.tsv"
    findings: list[dict[str, str]] = []
    if finding_path.is_file():
        header, loaded = _read_tsv(finding_path)
        if header == list(FINDING_COLUMNS):
            findings = loaded

    candidates = [row for row in issues if _is_workflow_candidate(row)]
    classified_refs = {
        source_ref
        for row in findings
        for source_ref in _split_refs(row.get("source_refs", ""))
    }
    untriaged = [
        row for row in candidates if row["source_ref"] not in classified_refs
    ]
    return {
        "repo_root": str(root),
        "ok": not validation_errors,
        "validation_errors": validation_errors,
        "validation_warnings": validation_warnings,
        "pilots_scanned": pilots_scanned,
        "issues_scanned": len(issues),
        "workflow_candidates": len(candidates),
        "findings": findings,
        "classified_source_issues": len(
            {ref for ref in classified_refs if ref in {c["source_ref"] for c in candidates}}
        ),
        "untriaged_candidates": untriaged,
    }


def write_learning_report(repo_root: Path) -> Path:
    root = Path(repo_root).resolve()
    summary = learning_summary(root)
    if not summary["ok"]:
        raise ValueError(
            "learning records are invalid: "
            + "; ".join(summary["validation_errors"])
        )

    lines = [
        "# Paper2Paper 跨 Pilot 学习报告",
        "",
        f"生成日期：{date.today().isoformat()}。",
        "",
        "## 这份报告是什么",
        "",
        "它把不同真实文献 Pilot 中标记为 workflow 候选的问题汇集到一起。"
        "原始问题仍保存在各 Pilot；这里的候选不自动等于 Paper2Paper 产品缺口。",
        "",
        "分类含义：",
        "",
        "- `paper_specific`：只在当前课题解决，不改产品核心；",
        "- `paper_type`：可能适用于同类论文，需用真实实例确认；",
        "- `general`：跨论文类型都可能导致相同误判；",
        "- `verified`：最小修复已实现，并留下来源 Pilot 与复验引用。",
        "- `observed` 只是待归类事实，不表示已确认产品缺口；",
        "",
        "聊天上下文不是持久记忆；本报告、来源 `issues.tsv`、"
        "`learning/findings.tsv` 和 Git 历史才是跨对话依据。",
        "",
        "## 汇总",
        "",
        f"- 扫描 Pilot：{summary['pilots_scanned']}",
        f"- 扫描问题：{summary['issues_scanned']}",
        f"- workflow 原始候选：{summary['workflow_candidates']}",
        f"- 已归类来源问题：{summary['classified_source_issues']}",
        f"- 尚未归类候选：{len(summary['untriaged_candidates'])}",
        "",
        "## 已归类发现",
        "",
        "| ID | 范围 | 状态 | 能力领域 | 归类后的问题 | 来源 |",
        "|---|---|---|---|---|---|",
    ]
    for row in summary["findings"]:
        lines.append(
            "| {finding_id} | {scope} | {status} | {area} | {failure} | {refs} |".format(
                finding_id=_escape_markdown(row["finding_id"]),
                scope=_escape_markdown(row["scope"]),
                status=_escape_markdown(row["status"]),
                area=_escape_markdown(row["capability_area"]),
                failure=_escape_markdown(row["generalized_failure"]),
                refs=_escape_markdown(row["source_refs"]),
            )
        )
    if not summary["findings"]:
        lines.append("| — | — | — | — | 尚无已归类发现 | — |")

    lines.extend(
        [
            "",
            "## 尚未归类的 workflow 候选",
            "",
            "这些记录只说明真实工作中出现过问题并建议考虑 workflow 动作。"
            "AI 仍须判断它是课题特有、论文类型级还是通用问题；不能按数量自动扩建核心。",
            "",
            "| 来源 | Pilot 路径 | 严重度 | 观察 | 对论文的后果 | 当前动作 |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in summary["untriaged_candidates"]:
        lines.append(
            "| {ref} | {path} | {severity} | {observation} | {consequence} | {action} |".format(
                ref=_escape_markdown(row["source_ref"]),
                path=_escape_markdown(row["pilot_path"]),
                severity=_escape_markdown(row.get("severity", "")),
                observation=_escape_markdown(row.get("observation", "")),
                consequence=_escape_markdown(row.get("consequence", "")),
                action=_escape_markdown(row.get("workflow_action", "")),
            )
        )
    if not summary["untriaged_candidates"]:
        lines.append("| — | — | — | 当前没有未归类候选 | — | — |")

    lines.extend(
        [
            "",
            "## 修复闭环",
            "",
            "1. 先在来源 Pilot 解决当前论文问题；",
            "2. 把原始事实写入该 Pilot 的 `issues.tsv`；",
            "3. 将候选归类为课题特有、论文类型级或通用；",
            "4. 只有可减少真实论文错误、时间或审查负担时，才做最小核心修改；",
            "5. 用来源 Pilot 复验；跨类型主张还需另一篇合适 Pilot；",
            "6. 保存测试/运行引用后才标记 `verified`；",
            "7. 立即回到目标论文的数据、Figure、分析与写作。",
            "",
            "本报告不是能力成熟度排行、创新性门槛或路线评分表。",
        ]
    )
    output = root / "reports" / "workflow-learning.md"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output
