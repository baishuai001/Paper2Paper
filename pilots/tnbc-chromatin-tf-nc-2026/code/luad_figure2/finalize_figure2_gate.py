#!/usr/bin/env python3
"""Combine data, anchor and extra-robustness verdicts without moving thresholds."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_outputs(root: Path, receipt: dict, report: str) -> None:
    out_dir = root / "audit/final_gate"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figure2_final_verdict.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    report_dir = root / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "luad-figure2-gate-report.md").write_text(report, encoding="utf-8")


def system_summary(record: dict) -> str:
    if record.get("evaluation_status") == "NOT_FULLY_EVALUATED_AFTER_EARLY_STOP":
        return (
            f"{record['qualified']}/{record['attempted']} 个已完成样本通过硬 QC；"
            f"{record.get('aborted_after_locked_gate', 0)} 个已启动样本在另一系统锁定 FAIL_DATA 后中止；"
            f"其余 {record['not_run']} 个未启动，均未计作 QC 失败"
        )
    if record.get("evaluation_status") == "INCOMPLETE_UNEXPECTEDLY":
        return (
            f"仅完成 {record['attempted']}/{record['expected']} 个样本；"
            "这是非预期的不完整运行"
        )
    return (
        f"{record['qualified']}/{record['observed']} 通过硬 QC；"
        f"失败比例 {100 * float(record['failure_fraction']):.1f}%"
    )


def data_failure(root: Path, raw: dict) -> None:
    atomic_bundle = {
        "status": "NOT_GENERATED_DUE_FAIL_DATA",
        "Figure2_complete": "NOT_GENERATED_DUE_FAIL_DATA",
        "SupplementaryFigure3_complete": "NOT_GENERATED_DUE_FAIL_DATA",
        "SupplementaryFigure4A-C_complete": "NOT_GENERATED_DUE_FAIL_DATA",
        "reason": "Biological-signal steps are forbidden after the pre-registered raw-data gate fails",
    }
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "final_verdict": "FAIL_DATA",
        "eligible_to_plan_Figure3": False,
        "continue_to_Figure3": False,
        "stop_after_this_gate": True,
        "data_gate": {
            "patient_n": 22,
            "PDX": raw["systems"]["PDX"],
            "cell_line": raw["systems"]["cell_line"],
            "early_stop_locked_system": raw.get("early_stop_locked_system"),
            "failure_reasons": raw.get("failure_reasons", ["raw_data_gate"]),
        },
        "anchor_logic": "NOT_EVALUATED_AFTER_DATA_FAILURE",
        "extra_robustness_not_part_of_anchor_method": "NOT_EVALUATED_AFTER_DATA_FAILURE",
        "atomic_bundle": atomic_bundle,
        "interpretation": "预注册数据充分性或硬 QC 失败；未评估生物学信号。",
    }
    reasons = ", ".join(raw.get("failure_reasons", ["raw_data_gate"]))
    report = f"""# LUAD Figure 2 一级闸门报告

## 结论

**最终判定：`FAIL_DATA`。** Figure 2 生物学信号未评估；本次运行在数据闸门停止。

- 原发 LUAD 患者：22 个公开 TCGA-LUAD 固定峰矩阵样本；
- LUAD PDX：{system_summary(raw['systems']['PDX'])}；
- 严格 LUAD 细胞系：{system_summary(raw['systems']['cell_line'])}；
- 失败原因：`{reasons}`。

未运行样本若源于另一系统已不可逆锁定数据闸门，仅标记为“早停后未运行”，不伪装成实验 QC 失败。

Figure 2、Supplementary Figure 3 和 Supplementary Figure 4A–C 均标记为 `NOT_GENERATED_DUE_FAIL_DATA`。这是原子交付规则在数据失败路径上的预期结果，不是漏图；生成其中任何一张生物学信号图都会违反停止规则。

不得通过放宽阈值、替换样本或切换 LUAD 亚型解释本次失败。`continue_to_Figure3 = FALSE`。
"""
    write_outputs(root, receipt, report)
    (root / "audit/final_gate/figure2_atomic_bundle_status.json").write_text(
        json.dumps(atomic_bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()
    raw = load(root / "audit/raw_atac_qc/figure2_raw_atac_qc_receipt.json")
    if raw["raw_data_gate"] != "PASS":
        data_failure(root, raw)
        return

    promoter = load(root / "audit/promoter_gate/promoter_gate_receipt.json")
    promoter_perm = load(root / "audit/promoter_gate/promoter_matched_permutation_receipt.json")
    motif = load(root / "audit/motif_gate/motif_gate_receipt.json")
    motif_perm = load(root / "audit/motif_gate/motif_matched_permutation_receipt.json")

    patient_n = int(promoter["patient_n"])
    data_failures: list[str] = []
    if patient_n < 7:
        data_failures.append("patient:qualified_below_minimum")
    if data_failures:
        verdict = "FAIL_DATA"
    elif motif["anchor_logic_verdict"] != "ANCHOR_PASS":
        verdict = "FAIL_SIGNAL"
    elif not promoter_perm["robustness_pass"] or not motif_perm["robustness_pass"]:
        verdict = "ANCHOR_PASS_ROBUSTNESS_FAIL"
    else:
        verdict = "PASS"

    eligible = verdict == "PASS"
    interpretations = {
        "PASS": "锚点阈值及两项预注册额外稳健性检验均通过；可另行规划 Figure 3，但本次运行仍在此停止。",
        "ANCHOR_PASS_ROBUSTNESS_FAIL": "TNBC 原文式筛选产生了信号，但没有超过匹配零假设；证据不足以进入 Figure 3。",
        "FAIL_SIGNAL": "TNBC 原文式染色质/motif 入口信号不足；停止于 Figure 2。",
        "FAIL_DATA": "一项或多项预注册数据充分性规则失败；不解释生物学信号。",
    }
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "final_verdict": verdict,
        "eligible_to_plan_Figure3": eligible,
        "continue_to_Figure3": False,
        "stop_after_this_gate": True,
        "data_gate": {
            "patient_n": patient_n,
            "PDX": raw["systems"]["PDX"],
            "cell_line": raw["systems"]["cell_line"],
            "failure_reasons": data_failures,
        },
        "anchor_logic": {
            "HC_TFs_after_promoter_activity": motif["HC_TFs_after_promoter_activity"],
            "motif_testable_HC_TFs": motif["motif_testable_HC_TFs"],
            "triple_system_motif_enriched_TFs": motif["triple_system_motif_enriched_TFs"],
            "triple_system_fraction_of_testable": motif["triple_system_fraction_of_testable"],
            "verdict": motif["anchor_logic_verdict"],
        },
        "extra_robustness_not_part_of_anchor_method": {
            "promoter_matched_permutation_p": promoter_perm["empirical_p_ge_observed"],
            "promoter_robustness_pass": promoter_perm["robustness_pass"],
            "motif_matched_permutation_p": motif_perm["empirical_p_ge_observed"],
            "motif_robustness_pass": motif_perm["robustness_pass"],
        },
        "interpretation": interpretations[verdict],
    }

    report = f"""# LUAD Figure 2 一级闸门报告

## 结论

**最终判定：`{verdict}`。** {receipt['interpretation']}

`eligible_to_plan_Figure3 = {str(eligible).upper()}`；`continue_to_Figure3 = FALSE`；`stop_after_this_gate = TRUE`。

## 数据与原子交付

- 原发 LUAD 患者：{patient_n} 个独立肿瘤。TCGA 原始/比对 ATAC 受控，因此使用 GDC 公开的两个技术重复固定峰计数；主定义要求两个重复均 CPM≥1，并同步报告三套阈值敏感性；
- LUAD PDX：{system_summary(raw['systems']['PDX'])}。先去除小鼠宿主 reads，再按 hg38 比对并以 paired-end fragment 模式调用峰；
- 严格 LUAD 细胞系：{system_summary(raw['systems']['cell_line'])}。重复提交、非 LUAD 与污染模型未计入；
- Figure 2、Supplementary Figure 3 与 Supplementary Figure 4A–C 共用同一最终 manifest、97-TF SHA、QC 表及阈值配置。

## TNBC 原文式锚点判据

- 三端启动子开放且三端平均 LUAD NES≥0 的 HC-TF：{motif['HC_TFs_after_promoter_activity']}（要求 ≥10）；
- motif 可检验 HC-TF：{motif['motif_testable_HC_TFs']}；
- 患者、PDX、细胞系三端均在至少一半样本达到 HOMER BH q<1e-5 的 TF：{motif['triple_system_motif_enriched_TFs']}（要求 ≥3）；
- 上述三端 motif TF 占 motif 可检验 HC-TF：{100 * motif['triple_system_fraction_of_testable']:.2f}%（要求 ≥10%）；
- 锚点逻辑判定：`{motif['anchor_logic_verdict']}`。

## 额外稳健性（不是 TNBC 原文判据）

- 三端启动子交集的 GC、启动子长度、TCGA-LUAD 表达匹配置换：经验 p={promoter_perm['empirical_p_ge_observed']:.6g}；pass={str(promoter_perm['robustness_pass']).upper()}；
- 三端 motif 交集按数据库可用性与 motif 数量分层置换：经验 p={motif_perm['empirical_p_ge_observed']:.6g}；pass={str(motif_perm['robustness_pass']).upper()}。

该层仅检验结果是否超过“普遍易开放/易匹配 motif”所能解释的程度，不替换 ARACNe3/VIPER，不重选 97 个 TF，也不修改锚点阈值。

## 停止规则

- `FAIL_DATA`：不解释生物学信号；修复同一比较的数据问题后才可重跑；
- `FAIL_SIGNAL`：当前 LUAD-vs-LUSC TF 程序不进入 Figure 3，不得事后放宽阈值；
- `ANCHOR_PASS_ROBUSTNESS_FAIL`：可报告“按锚点能筛出结果”，但证据不足以进入 Figure 3；
- `PASS`：仅表示可另行规划 Figure 3；本次任务仍在本闸门停止。
"""
    write_outputs(root, receipt, report)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
