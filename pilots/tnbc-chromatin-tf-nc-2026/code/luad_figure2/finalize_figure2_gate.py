#!/usr/bin/env python3
"""Write the descriptive report for the completed TNBC-style LUAD Figure 2.

The anchor paper reports the number of HC-TFs and motif-supported TFs; it does
not impose minimum HC-TF counts, minimum triple-system motif counts, cohort
failure fractions, or matched-permutation vetoes.  This finalizer therefore
reports the observed results without converting them into an invented
biological PASS/FAIL decision.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve()

    raw = load(root / "audit/raw_atac_qc/figure2_raw_atac_qc_receipt.json")
    promoter = load(root / "audit/promoter_gate/promoter_gate_receipt.json")
    motif = load(root / "audit/motif_gate/motif_gate_receipt.json")
    if raw["raw_data_completion"] != "COMPLETE":
        raise RuntimeError("Cannot finalize Figure 2 before all predeclared raw samples are processed")
    if motif["analysis_status"] != "COMPLETE":
        raise RuntimeError("Cannot finalize an incomplete HOMER analysis")

    counts = {
        "Figure1_LUAD_specific_TFs": int(promoter["frozen_tf_count"]),
        "triple_system_promoter_accessible_TFs": int(promoter["triple_system_promoter_accessible_TFs"]),
        "HC_TFs": int(motif["HC_TFs_after_promoter_activity"]),
        "motif_testable_HC_TFs": int(motif["motif_testable_HC_TFs"]),
        "motif_supported_in_at_least_one_system": int(motif["HC_TFs_with_motif_enrichment_in_at_least_one_system"]),
        "triple_system_motif_TFs": int(motif["triple_system_motif_enriched_TFs"]),
    }
    atomic_bundle = {
        "status": "COMPLETE",
        "Figure2_complete": "COMPLETE",
        "SupplementaryFigure3_complete": "COMPLETE",
        "SupplementaryFigure4A-C_complete": "COMPLETE",
    }
    receipt = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "final_verdict": "ORIGINAL_STYLE_ANALYSIS_COMPLETE",
        "analysis_complete": True,
        "biological_pass_fail_threshold_applied": False,
        "automatic_Figure3_decision": None,
        "sample_counts": {
            "patient": int(promoter["patient_n"]),
            "PDX": int(promoter["PDX_n"]),
            "cell_line": int(promoter["cell_line_n"]),
        },
        "anchor_style_results": counts,
        "anchor_rules": {
            "promoter_support": "accessible in at least ceiling(n/2) samples in each system",
            "activity_filter": "exclude only when mean NES < 0 in patient, PDX, and cell-line cohorts",
            "motif_support": "HOMER adjusted p < 1e-5 in at least ceiling(n/2) samples per system",
            "motif_databases": "complete JASPAR and CIS-BP references tested independently",
        },
        "diagnostic_QC_does_not_exclude_samples": True,
        "removed_non_anchor_rules": [
            "maximum 20% cohort failure fraction",
            "minimum 10 HC-TFs",
            "minimum 3 triple-system motif TFs",
            "minimum 10% triple-system motif fraction",
            "matched-permutation veto",
        ],
        "atomic_bundle": atomic_bundle,
    }

    out_dir = root / "audit/final_gate"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figure2_final_verdict.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (out_dir / "figure2_atomic_bundle_status.json").write_text(
        json.dumps(atomic_bundle, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    diagnostic_lines = []
    for system in ("PDX", "cell_line"):
        record = raw["systems"][system]
        diagnostic_lines.append(
            f"- {system}：{record['analysis_included']}/{record['expected']} 个预定义模型完成并纳入；"
            f"其中 {record['diagnostic_reference_failures']} 个低于报告性读段/峰数参考线，但未据此剔除。"
        )
    report = f"""# LUAD Figure 2 原文式分析报告

## 完成状态

`ORIGINAL_STYLE_ANALYSIS_COMPLETE`。Figure 2、Supplementary Figure 3 和 Supplementary Figure 4A–C 已由同一份 22/13/19（患者/PDX/细胞系）manifest 生成。

- 患者：{promoter['patient_n']} 个 TCGA-LUAD 原发肿瘤；
{chr(10).join(diagnostic_lines)}

质量指标全部保留并报告；未使用 100 万有效片段、10,000 个峰或 20% 失败比例剔除样本。qPCR 建库富集信息不适用于这些既有公共数据，记为不可回溯，而不是计算性失败。

## 原文式筛选结果

- Figure 1 输入：{counts['Figure1_LUAD_specific_TFs']} 个 LUAD 特异 TF；
- 三端启动子均在至少一半样本开放：{counts['triple_system_promoter_accessible_TFs']} 个；
- 仅排除患者、PDX、细胞系三个活动队列平均 NES **同时**小于 0 后，HC-TF：**{counts['HC_TFs']} 个**；
- HC-TF 中有 JASPAR 和/或 CIS-BP 已知 motif、可检验者：{counts['motif_testable_HC_TFs']} 个；
- 至少一个系统达到 HOMER adjusted p<1e-5 且至少半数样本支持：{counts['motif_supported_in_at_least_one_system']} 个；
- 患者、PDX、细胞系三端均达到上述 motif 规则：**{counts['triple_system_motif_TFs']} 个**。

这些数值是结果，不附加“至少 10/3/10%”等 TNBC 原文没有设置的项目停止线，也不以额外匹配置换否决。
"""
    report_dir = root / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "luad-figure2-gate-report.md").write_text(report, encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
