# LUAD Figure 2 一级闸门报告

> **历史 v1 收据，已被取代。** 本文件只记录后来废止的额外“样本失败比例不超过 20%”规则及当时的早停路径，不是当前 Figure 2 结论。现行原文式 v2 分析已完整运行并通过 28/28 项原子核验；请参见项目 `reports/luad-figure2-final-report.md`。

## 结论

**最终判定：`FAIL_DATA`。** Figure 2 生物学信号未评估；本次运行在数据闸门停止。

- 原发 LUAD 患者：22 个公开 TCGA-LUAD 固定峰矩阵样本；
- LUAD PDX：0/0 个已完成样本通过硬 QC；1 个已启动样本在另一系统锁定 FAIL_DATA 后中止；其余 12 个未启动，均未计作 QC 失败；
- 严格 LUAD 细胞系：14/19 通过硬 QC；失败比例 26.3%；
- 失败原因：`cell_line:failure_fraction_above_0.20`。

未运行样本若源于另一系统已不可逆锁定数据闸门，仅标记为“早停后未运行”，不伪装成实验 QC 失败。

Figure 2、Supplementary Figure 3 和 Supplementary Figure 4A–C 均标记为 `NOT_GENERATED_DUE_FAIL_DATA`。这是原子交付规则在数据失败路径上的预期结果，不是漏图；生成其中任何一张生物学信号图都会违反停止规则。

不得通过放宽阈值、替换样本或切换 LUAD 亚型解释本次失败。`continue_to_Figure3 = FALSE`。
