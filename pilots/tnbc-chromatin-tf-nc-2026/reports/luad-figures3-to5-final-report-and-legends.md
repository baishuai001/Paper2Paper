# LUAD Figures 3–5 结果、判定与 Figure Legends

日期：2026-08-20  
输入：Figure 2 冻结的 31 个 LUAD HC-TF；三系统严格 motif 子集为 FOXA3、NFATC4、XBP1。

## 1. 总判定

| 证据层 | 运行状态 | 科学判定 |
|---|---|---|
| Figure 3 调控网络结构 | 完成并通过独立核验 | 31 个 HC-TF 可形成可解释的有符号 regulon、3 个协同模块和明显的患者间异质性 |
| Figure 4 患者结局 | 完成并通过独立核验 | 有若干名义显著 TF，但跨终点/跨队列和置换检验支持有限；不能宣称稳健独立预后框架 |
| Figure 5A–F 细胞系药敏 | 完成 | 三个药敏库共完成 39,085 个合格关联检验，354 条数据集内显著关联，5 条药物–TF 关联在至少两个数据库同方向重复 |
| Figure 5G–H PDX 验证 | 完成，结果为阴性 | 14 个计算判定的 LUAD-like PDX、16 个样本数合格药物、423 个探索性检验中，没有预先要求的“细胞系重复且 PDX 全局 FDR≤0.05”关联；治疗选择终点未通过 |

因此，分析可以并且已经推进到证据链末端，但结果只支持“染色质优先的 LUAD TF 调控网络及选择性药敏关联”，不支持与 TNBC 原文等强的“可指导治疗选择”结论。

细胞系层面同方向复现的 5 个组合均表现为“TF 活性越高，药物耐受倾向越强”：

| 药物 | TF | 复现数据库 |
|---|---|---|
| 5-Fluorouracil | ZNF254 | CTRPv2；GDSC2 |
| 5-Fluorouracil | ZNF540 | CTRPv2；GDSC2 |
| PHA-793887 | WWC2 | CTRPv2；PRISM |
| Vorinostat | ZNF254 | GDSC2；PRISM |
| Dactolisib | ZNF254 | GDSC2；PRISM |

这些是统计关联，不表示相应 TF 是药物直接靶点；并且没有组合通过预先规定的 PDX 验证。

## 2. 最终图件位置

| 图件 | 文件 |
|---|---|
| Figure 3 | [PDF](../execution/luad-figure3/results/figures/Figure3_complete.pdf) / [PNG](../execution/luad-figure3/results/figures/Figure3_complete.png) |
| Supplementary Figure 6 | [PDF](../execution/luad-figure3/results/figures/SupplementaryFigure6_complete.pdf) / [PNG](../execution/luad-figure3/results/figures/SupplementaryFigure6_complete.png) |
| Figure 4 | [PDF](../execution/luad-figure4/results/figures/Figure4_complete.pdf) / [PNG](../execution/luad-figure4/results/figures/Figure4_complete.png) |
| Supplementary Figure 7 | [PDF](../execution/luad-figure4/results/figures/SupplementaryFigure7_complete.pdf) / [PNG](../execution/luad-figure4/results/figures/SupplementaryFigure7_complete.png) |
| Figure 5 | [PDF](../execution/luad-figure5/results/figures/Figure5_complete.pdf) / [PNG](../execution/luad-figure5/results/figures/Figure5_complete.png) |
| Supplementary Figure 8 | [PDF](../execution/luad-figure5/results/figures/SupplementaryFigure8_complete.pdf) / [PNG](../execution/luad-figure5/results/figures/SupplementaryFigure8_complete.png) |
| Supplementary Figure 9 | [PDF](../execution/luad-figure5/results/figures/SupplementaryFigure9_complete.pdf) / [PNG](../execution/luad-figure5/results/figures/SupplementaryFigure9_complete.png) |

## 3. Figure legends

### Figure 3. Network organization and intertumour heterogeneity of the 31 chromatin-prioritized LUAD HC-TFs

**(A)** Signed regulon composition of the 31 HC-TFs in the TCGA lung-cancer ARACNe3 network. Targets are separated into activated versus repressed and shared versus TF-private components. Stars mark FOXA3, NFATC4 and XBP1, the HC-TFs with strict three-system motif support. The 31 regulons contain 3,011 signed TF–target edges and 2,776 unique target genes.

**(B)** HC-TF communities and associated biological processes. Communities were detected by Louvain clustering of the full positive shared-activated-target graph; activated targets from each of the three modules were tested for GO Biological Process enrichment against the union of HC-TF targets. Up to three top BH-adjusted terms per module are shown.

**(C)** Pairwise Pearson correlations between the activities of the 31 HC-TFs across 516 TCGA-LUAD tumors. Rows and columns were hierarchically clustered using the same correlation matrix.

**(D)** Collaboration network based on shared activated targets. Nodes represent HC-TFs and are colored by Louvain module; node size reflects regulon breadth. The full graph contained 134 TF-pair edges with at least one shared activated target and was used for community detection. For legibility, the panel displays edges in the top 20% of shared-target weights.

**(E)** Relationship between regulon breadth and the number of collaborating HC-TF partners. Each point is one HC-TF, colored by module; the line and band show the linear fit and 95% confidence interval.

**(F)** Intertumour heterogeneity of HC-TF activity in 516 TCGA-LUAD tumors. Skewness was tested using its asymptotic normal statistic and BH correction; strong skew required FDR≤0.05 and absolute skewness >1. ELF3, MAGED2 and PHC2 met this rule. Point shape indicates strict three-system motif support.

**(G)** Representative activity distributions for the two most positively and two most negatively skewed TFs, shown separately in TCGA-LUAD and TCGA-LUSC.

### Supplementary Figure 6. Proliferation association and cross-system heterogeneity of LUAD HC-TFs

**(A)** Pearson correlation between each HC-TF activity and MKI67 expression across TCGA-LUAD tumors. Bars are ordered by correlation; color indicates BH FDR≤0.05.

**(B)** Skewness of each HC-TF activity distribution across TCGA-LUAD primary tumors, LUAD PDMR PDXs and LUAD DepMap cell lines. The panel is descriptive because the three systems have different sample sizes and sampling structures.

### Figure 4. Association of LUAD HC-TF activity with patient outcome

**(A–D)** Multivariable Cox models for **(A)** overall survival (OS) in GSE41271 (n=178; 69 events), **(B)** recurrence-free survival (RFS) in GSE41271 (n=177; 69 events), **(C)** OS in TCGA-LUAD (n=467; 165 events) and **(D)** disease-free survival (DFS) in TCGA-LUAD (n=277; 80 events). TF activity was dichotomized at the cohort/endpoint median and follow-up was censored at 10 years. Models included age, sex, pathologic stage and smoking history. Diamonds and lines show log2 hazard ratios and 95% confidence intervals; blue and red denote nominally favorable and adverse associations, respectively. BH FDR was calculated across 31 TFs within each cohort and endpoint.

**(E)** Overlap of nominally significant favorable and adverse TFs between OS and RFS in GSE41271. ETV1 was the only favorable TF shared by both endpoints; no adverse TF was shared.

**(F)** Kaplan–Meier curves for representative GSE41271 TFs selected from the multivariable cross-endpoint ranking (ETV1 and ZNF444), shown for OS and RFS. Cutpoints were estimated with `surv_cutpoint` and these panels are descriptive rather than independent tests.

**(G)** Overlap of nominally significant favorable and adverse TFs between OS and DFS in TCGA-LUAD. CREBRF was the only favorable TF shared by both endpoints; no adverse TF was shared.

**(H)** Kaplan–Meier curves for representative TCGA-LUAD TFs CREBRF and PHC2, shown for OS and DFS with the same descriptive cutpoint procedure.

### Supplementary Figure 7. Permutation assessment of the number of prognostic HC-TFs

Null distributions of the number of nominally significant HC-TFs obtained from 5,000 survival-outcome permutations for univariable and multivariable models in each cohort and endpoint. The observed multivariable excess was not significant for GSE41271 OS (empirical P=0.3242), GSE41271 RFS (P=0.0700), TCGA OS (P=0.1402) or TCGA DFS (P=0.6782). Only the TCGA OS univariable count exceeded its null distribution (P=0.0404). These results constrain Figure 4 to exploratory prognostic evidence.

### Figure 5. Pharmacogenomic associations of LUAD HC-TF activity and independent PDX validation attempt

**(A–C)** Volcano plots of concordance-index associations between the 31 HC-TF activities and drug sensitivity in **(A)** GDSC2, **(B)** CTRPv2 and **(C)** PRISM. Only drugs with at least 10 paired LUAD cell lines and at least one AAC>0.2 were tested. CI>0.5 denotes higher TF activity associated with greater sensitivity; CI<0.5 denotes resistance. BH correction was performed across all tested drug–TF pairs within each dataset.

**(D)** Membership patterns of significant drug–TF pairs across the three pharmacogenomic datasets. A replicated association required the same canonical compound and TF, FDR≤0.05 in at least two datasets and the same CI direction.

**(E)** Concordance-index effect matrix for the cross-dataset replicated drug–TF pairs. Five pairs met the same-direction replication rule; red and blue denote sensitivity- and resistance-directed effects, respectively.

**(F)** Number and direction of replicated drug associations per HC-TF, annotated where the same TF had a nominal multivariable prognostic association in Figure 4. This integration is hypothesis-generating and does not establish that a TF is the molecular target of a drug.

**(G)** Prespecified in-vivo validation result. PDXE v2 baseline RNA was projected through the frozen TCGA regulon; an externally audited LUAD-versus-LUSC activity classifier identified 14 of 27 RNA-profiled NSCLC PDXs as LUAD-like. No drug–TF pair that replicated in at least two cell-line datasets reached global PDX FDR≤0.05 with at least 10 LUAD-like PDXs; the panel therefore reports the negative validation result instead of a nominal scatterplot.

**(H)** PDX response waterfall was not generated because the Figure 5G validation criterion was not met. The explicit negative panel preserves the frozen stopping rule.

### Supplementary Figure 8. Cell-line mapping, drug coverage and lung-state classifier validation

**(A)** Number of histologically annotated LUAD cell lines with paired TF activity and drug-response identifiers in GDSC2, CTRPv2 and PRISM (48, 59 and 45, respectively).

**(B)** Number of drugs entering association testing after the n≥10 and AAC>0.2 rules (117, 383 and 768, respectively).

**(C)** Distribution of the frozen LUAD-minus-LUSC TF-activity classifier score in TCGA training tumors and independent PDMR and DepMap lung models. The classifier used pre-defined 158-TF LUAD and 193-TF LUSC signatures, label-free row-z normalization within each lung cohort and a threshold fitted once in TCGA. Balanced accuracy was 0.940 in TCGA, 0.904 in PDMR and 0.690 in DepMap. The PDXE assignment is consequently reported as a computational LUAD-like state, not a pathology call.

### Supplementary Figure 9. PDXE drug coverage and exploratory response screen

**(A)** Single-agent coverage among the 14 computationally LUAD-like PDXE models. Bars show the number of evaluable PDXs per drug and whether the compound could be mapped to a cell-line pharmacogenomic identity.

**(B)** Exploratory concordance-index screen across 16 drugs with at least 10 LUAD-like PDXs and 31 HC-TFs (423 finite tests). CI>0.5 denotes higher TF activity associated with lower best average response and therefore greater sensitivity. Points are highlighted when the same drug–TF direction replicated in the cell-line datasets. No highlighted pair passed global BH FDR≤0.05.

## 4. 论文主张边界

当前结果可支持的主题为：

> A chromatin-informed TF-activity framework defines a LUAD regulatory network with a selective cross-system motif core, marked intertumour heterogeneity and a small set of reproducible cell-line pharmacogenomic associations.

当前结果不能支持：

- 31 个 HC-TF 普遍具有三系统 motif 保守；严格支持仅 3 个。
- 31-TF 程序形成跨队列稳健预后分类器；置换和跨队列结果不足。
- TF 活性已经可以指导 LUAD PDX 或患者治疗选择；严格 PDX 终点为阴性。
