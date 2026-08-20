# LUAD Figure 1 完成审计与 Figure Legends

日期：2026-08-20  
状态：**完成并通过独立原子核验**。本报告取代此前将 GSE81089 分支标为暂定的版本。

## 1. 修订结果

- GSE81089 已补入两个官方处理列别名：`L608T_2122 -> L608T`、`L771T_1 -> L771T`。最终为 108 例 LUAD、67 例 LUSC；旧版 106/67 不再使用。
- 修订后，TCGA 发现的 158 个 LUAD-specific TF 中有 95 个在 GSE81089 方向一致复现；旧版 97 个不再使用。
- 已加入 GSE41271 Illumina HumanWG-6 v3 微阵列跨平台验证：183 例 LUAD、80 例 LUSC、12 例其他组织学；70/158 个 LUAD TF 方向一致复现。
- GSE81089 与 GSE41271 共同支持 44/158 个 LUAD TF。
- Figure 2 的冻结输入仍为全部 158 个 TCGA LUAD-specific TF。外部队列复现数的修订不改变 Figure 2 的输入、HC-TF 或 motif 结果。
- 已生成主图和两套补充图的最终合并 PDF/PNG；核验文件为 `results/figure1_complete/figure1_complete_independent_verification.json`。

## 2. 最终图件位置

| 图件 | 最终文件 |
|---|---|
| Figure 1 | [PDF](../execution/luad-figure1-v2/results/figures/Figure1_complete.pdf) / [PNG](../execution/luad-figure1-v2/results/figures/Figure1_complete.png) |
| Supplementary Figure 1 | [PDF](../execution/luad-figure1-v2/results/figures/SupplementaryFigure1_complete.pdf) / [PNG](../execution/luad-figure1-v2/results/figures/SupplementaryFigure1_complete.png) |
| Supplementary Figure 2 | [PDF](../execution/luad-figure1-v2/results/figures/SupplementaryFigure2_complete.pdf) / [page 1 PNG](../execution/luad-figure1-v2/results/figures/SupplementaryFigure2_complete_page1.png) / [page 2 PNG](../execution/luad-figure1-v2/results/figures/SupplementaryFigure2_complete_page2.png) |

## 3. Figure legends

### Figure 1. ARACNe3/VIPER identifies a LUAD-specific TF-activity program that transfers across patients and experimental models

**(A)** Analytical framework. An ARACNe3 interactome was reconstructed from 1,017 independent primary TCGA lung tumors (516 LUAD and 501 LUSC), followed by VIPER/msVIPER analysis to identify 158 LUAD-specific and 193 LUSC-specific TFs. Independent networks were reconstructed in GSE81089 RNA-seq (108 LUAD and 67 LUSC) and GSE41271 microarray data (183 LUAD and 80 LUSC; 12 other lung-cancer histologies retained for network inference but not the LUAD–LUSC contrast). Direction-concordant replication supported 95/158 LUAD TFs in GSE81089, 70/158 in GSE41271 and 44/158 in both cohorts. The frozen TCGA regulon was separately projected onto 56 PDMR PDXs and 103 DepMap cell lines. External replication did not reduce the 158-TF Figure 2 input.

**(B)** Concordance between differential TF activity and TF expression in the TCGA discovery cohort. Each point denotes one regulator and shows its msVIPER normalized enrichment score (NES; LUAD versus LUSC) and limma-voom log2 fold change. Colors distinguish LUAD-specific, LUSC-specific, VIPER-only, expression-only and non-significant regulators. The line and band show the linear fit and 95% confidence interval (Pearson r = 0.7255, P < 2.2 × 10⁻¹⁶).

**(C)** Single-sample TF activity across 1,017 TCGA tumors. Rows are the 351 histology-associated TFs and columns are independent primary tumors. Activity was standardized within TF; red and blue denote relatively higher and lower activity. Tumors were clustered without using histology labels. Annotations show histology, age, sex, vital status, pathologic stage, smoking exposure, published expression subtype and selected driver alterations; gray denotes unavailable values.

**(D)** Projection of the frozen TCGA regulon onto 56 PDMR lung PDX models (22 LUAD and 34 LUSC). The heat map shows row-standardized activity for the evaluable TCGA-derived TFs. Clinical and molecular annotations were retained without imputation.

**(E)** Projection of the frozen TCGA regulon onto 103 DepMap 22Q2 lung-cancer cell lines (76 LUAD and 27 LUSC). The heat map shows row-standardized activity and annotations for histology, reported demographic variables, primary/metastatic origin, collection site, growth mode and selected driver alterations.

### Supplementary Figure 1. Derivation of the patient cohorts used for TF-network discovery and replication

**(A)** TCGA discovery-cohort derivation. From 1,141 expression columns with count and TPM estimates, 1,029 primary-tumor aliquots were retained; removal of 12 duplicate patient aliquots yielded 1,017 independent primary tumors (516 LUAD and 501 LUSC).

**(B)** GSE81089 RNA-seq replication-cohort derivation. Nineteen matched normal tissues and 24 large-cell/not-otherwise-specified tumors were removed from 218 samples. Canonicalization of the two documented processed-column aliases recovered all target tumors, yielding 175 tumors (108 LUAD and 67 LUSC).

**(C)** GSE41271 microarray validation-cohort derivation. The 275 primary lung tumors comprised 183 LUAD, 80 LUSC and 12 other histologies. The whole cohort was available for independent network inference; the exact LUAD and LUSC samples defined the differential contrast. For genes represented by multiple probes, the single-symbol probe with the highest median absolute deviation was retained.

### Supplementary Figure 2. Independent replication and cross-system structure of the LUAD TF program

**(A)** GSE81089 relationship between msVIPER NES and expression log2 fold change in 108 LUAD versus 67 LUSC tumors (Pearson r = 0.72, P = 7.9 × 10⁻³²³). Plot elements and regulator categories are as in Figure 1B.

**(B)** Direction-concordant overlap of LUAD- and LUSC-specific TFs between TCGA and GSE81089. The corrected GSE81089 run supports 95/158 TCGA LUAD-specific TFs.

**(C–E)** Pairwise Pearson correlations between TF-activity profiles in **(C)** TCGA primary tumors, **(D)** PDMR PDX models and **(E)** DepMap cell lines. Rows and columns use the same average-linkage clustering of (1-r); histology is indicated on both axes.

**(F)** Cross-system activity effects for the 95 LUAD TFs replicated in GSE81089. LUAD-minus-LUSC mean activity differences were standardized separately in TCGA, GSE81089, PDMR and DepMap before visualization; this panel describes transferability and does not define the Figure 2 input.

**(G)** Independent GSE41271 microarray relationship between msVIPER NES and LUAD-versus-LUSC expression fold change (Pearson r = 0.7837, P < 2.2 × 10⁻¹⁶).

**(H)** Direction-concordant TF-membership combinations across TCGA, GSE81089 and GSE41271. Among the 158 TCGA LUAD-specific TFs, 70 replicate in GSE41271 and 44 are jointly supported by both external patient cohorts.

## 4. 对后续图的影响

Figure 2–5 使用的是 TCGA 发现的 158 个 TF及其后产生的 31 个 HC-TF，而不是 95 个 GSE81089 复现 TF。因此本次修订只更新 Figure 1、Supplementary Figures 1–2、图例和跨图文字；Figure 2 的 ATAC、promoter gate、motif gate 及其下游 Figure 3–5 无需因这两例样本重跑。
