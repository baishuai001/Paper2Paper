# LUAD Paper2Paper 科学问题、研究计划与最终图件索引

版本：2026-08-20  
状态：Figure 1–5 及同步补充图已运行至预先规定的证据链终点。本文档是入口索引；详细结果与 Figure Legends 以分图最终报告为准。

## 科学问题

本项目检验：能否从 LUAD 相对 LUSC 的大规模转录组中识别 TF 活性程序，以患者、PDX 和细胞系染色质数据优先化高置信 TF，进而解析 LUAD 内部调控网络异质性、患者结局关联和可重复药物脆弱性。

这不是 LUAD 内部亚型的发现研究。Figure 1 的比较轴是 `LUAD versus LUSC`；Figure 3–5 才在 LUAD 内部检验网络异质性、结局和药敏。

## 完成状态

| 阶段 | 输入与结果 | 判定 |
|---|---|---|
| Figure 1 | TCGA 1,017 例发现；GSE81089 175 例 RNA-seq 验证；GSE41271 263 例 LUAD/LUSC microarray 跨平台验证；PDMR 56 个 PDX 与 DepMap 103 条细胞系投射 | 完成。158 个 TCGA LUAD-specific TF；95/158 在 GSE81089、70/158 在 GSE41271、44/158 被两个外部患者队列共同支持 |
| Figure 2 | 158 个 TF；22 个患者、13 个 PDX、19 条细胞系 ATAC | 完成。31 个 HC-TF；20 个 motif 可测试；FOXA3、NFATC4、XBP1 获严格三系统 motif 支持 |
| Figure 3 | 31 个 HC-TF | 完成。3,011 条 signed edges、2,776 个靶基因、3 个调控模块，并存在患者间 TF 活性异质性 |
| Figure 4 | 31 个 HC-TF；GSE41271 与 TCGA-LUAD 生存终点 | 完成。存在若干名义关联，但跨终点、跨队列和置换支持有限；稳健预后框架未成立 |
| Figure 5 | GDSC2、CTRPv2、PRISM 与 PDXE v2 | 完成。5 个药物–TF 组合在至少两个细胞系药敏库同方向复现；严格 PDX 验证为 0，治疗选择终点未通过 |

## 当前可以支持的主题

> A chromatin-informed TF-activity framework defines a LUAD regulatory network with a selective cross-system motif core, marked intertumour heterogeneity and a small set of reproducible cell-line pharmacogenomic associations.

不能把结果写成“广泛三系统 motif 保守”“稳健跨队列预后分层”或“已能指导 LUAD 治疗选择”。

## 最终报告与 Figure Legends

- [Figure 1、Supplementary Figures 1–2](luad-figure1-completion-audit-and-legends.md)
- [Figure 2、Supplementary Figures 3、4A–C、5](luad-figure2-final-report.md)
- [Figures 3–5、Supplementary Figures 6–9](luad-figures3-to5-final-report-and-legends.md)

## 主图位置

- [Figure 1 PDF](../execution/luad-figure1-v2/results/figures/Figure1_complete.pdf)
- [Figure 2 PDF](../execution/luad-figure2/results/figures/Figure2_complete.pdf)
- [Figure 3 PDF](../execution/luad-figure3/results/figures/Figure3_complete.pdf)
- [Figure 4 PDF](../execution/luad-figure4/results/figures/Figure4_complete.pdf)
- [Figure 5 PDF](../execution/luad-figure5/results/figures/Figure5_complete.pdf)

所有原始数据下载和分析均在指定云服务器完成；本地仅保存代码、协议、可审计结果和图件。
