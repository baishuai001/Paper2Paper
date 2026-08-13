# CRC M-vs-rest 调控程序一级闸门：预先冻结协议

冻结日期：2026-08-13  
冻结状态：任何真实 M-vs-rest TF 活性结果产生前锁定  
真实数据执行位置：仅限云服务器 `/media/desk16/iy13202/projects/Paper2Paper`  
表型清单：`analysis/phenotypes/m-vs-rest.json`

## 1. 闸门回答什么

将 TNBC 锚点论文的入口方法迁移到 CRC，回答：CRC-atlas 作者预先定义的 M 型相对 B、T、desert 三型，是否在原发、未治疗、未富集样本的患者级 Cancer-cell pseudobulk 中对应一个跨数据集可重复的 CRC 特异 TF 活性程序？

本闸门只检验关联和跨数据集可重复性，不声称免疫细胞驱动癌细胞、TF 具有因果作用或已经能够指导治疗。

此前执行的“四种免疫表型能否被独立聚类完整重建”是额外的质量控制，不是 TNBC 原文进入 ARACNe3/VIPER 的必要步骤。其 FAIL 收据永久保留，但不再作为本闸门的硬停止条件。M-vs-rest 标签不在本次分析中重新聚类或优化。

## 2. 与 TNBC 原文一致及新增的部分

主方法保持：大样本肿瘤组织 RNA 表达构建癌种特异 ARACNe3 网络，`viper` 包将网络转为 regulon，VIPER/aREA 计算逐样本 TF 活性，msVIPER 检验组间差异活性。调控因子使用 TNBC Supplementary Data 2 中 PAN-GO “regulators of transcription”清单。论文 Methods 称 2,139 genes，但实际 sheet 为 1 行表头加 2,138 行非空 gene-symbol 数据；其中 79 行为重复 symbol，去重后是 2,059 个唯一符号。ARACNe3 输入采用可审计补充表的 2,059 个唯一符号，并将论文计数、sheet 行数和去重数同时写入收据。

CRC 迁移新增但不替代原方法的部分是：患者级 Cancer-cell pseudobulk、按原始数据集分层的效应量、随机效应 meta 分析和无泄漏 leave-one-dataset-out（LODO）验证。这些步骤用于处理 CRC-atlas 汇总多个研究造成的批次和重复性问题。

CollecTRI、DoRothEA、ULM 或其他先验网络均不是主分析；若以后运行，只能标记为敏感性分析，不能替代 ARACNe3/VIPER 判定。

## 3. 输入与分析单位

### 3.1 CRC 特异网络

- 数据：GDC 的 TCGA-COAD 与 TCGA-READ `STAR - Counts`、`Primary Tumor` 文件；使用其中 `tpm_unstranded`。
- 重复 aliquot：按 TCGA participant 合并，逐基因取 TPM 均值，避免同一患者重复计权。
- 基因符号：同一符号的 Ensembl 行求和；空符号删除。
- 预过滤：TPM >= 1 的患者比例至少 10%，且跨患者方差大于 0。
- ARACNe3：Califano Lab 官方仓库，固定提交 `3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14`。
- 参数：100 个 subnetworks、默认 0.63212 subsampling、每个 subnetwork FDR alpha 0.05、默认 Maximum-Entropy pruning、seed 1729、24 threads。
- 共识边：将 `log.p.values` 还原为 binomial p 值，对全部合并边作 BH 校正，保留 FDR <= 0.05。
- regulon：`viper::aracne2regulon(..., format="3col")`；TF-target mode 由 TCGA CRC `log2(TPM+1)`估计。

### 3.2 CRC-atlas 癌细胞表达

- 输入：已审计 CRC-atlas H5AD；必须匹配 30,875,155,333 bytes 和 SHA256 `718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6`。
- 单位：患者 `donor_id`，不是细胞或样本。
- 范围和分组完全来自表型清单；同一患者范围内标签必须唯一。
- 只汇总作者标注 `Cancer cell` 的 `raw/X` 非负整数 counts；每位患者所有合规原发样本合并。
- 主阈值：至少 50 个 Cancer cells；20、100、200 仅作预设敏感性分析。
- 表达：患者内求和，转为 `log2(CPM+1)`；保留至少 10%患者 CPM >= 1 的基因。

## 4. 数据与网络可判定条件

必须全部满足，否则结论为 `INDETERMINATE`，不能写成生物学阴性：

1. TCGA CRC 独立 participant >= 550，过滤后基因 >= 10,000，表达矩阵无重复样本或重复基因符号；
2. PAN-GO 清单审计明确记录论文称 2,139 genes、sheet 实有 2,138 个非空记录和 2,059 个唯一 gene symbols，且其中至少 1,500 个唯一符号进入 ARACNe3；
3. 共识网络至少 500 个 regulator；在 CRC-atlas 可测基因上，至少 400 个 regulon 具有 >=25 个 targets；
4. 主阈值下 M >=30、non-M >=80；
5. 至少 4 个数据集各含 >=3 M 和 >=3 non-M；
6. 任一数据集贡献的 M 患者不超过全部 M 的 60%；
7. 患者标签、数据集归属、Cancer-cell 数量和 pseudobulk 列顺序均通过一一对应检查。

## 5. VIPER 与 msVIPER

- 单样本活动：对全部合规患者的 gene-by-patient `log2(CPM+1)`矩阵运行 `viper(..., method="scale", minsize=25, nes=TRUE)`，得到逐患者 NES。
- 组间 signature：每个信息性数据集内计算 M-vs-rest 基因 Welch t statistic，再以有效样本量平方根加权 Stouffer 合并。
- null：在每个数据集内部置换 M 标签 500 次，重新计算完整 gene signature；该矩阵传给 `msviper`。
- msVIPER 显著性阈值沿用锚点论文：BH-FDR <= 0.01。

## 6. 跨数据集 TF 程序

对每个 VIPER TF，在每个同时有 >=3 M 和 >=3 non-M 的数据集计算 Hedges g（正值表示 M 活性更高），再作 REML 随机效应 meta 分析。一个“可重复 TF”必须同时满足：

- 至少 4 个信息性数据集；
- meta BH-FDR <= 0.05；
- |meta Hedges g| >= 0.50；
- 至少 75%信息性数据集与 meta 同方向；
- I2 <= 75%。

“可重复 TF 程序”预先定义为至少 10 个可重复 TF，其中至少 5 个同时满足 msVIPER FDR <=0.01 且方向一致。

## 7. 无泄漏 LODO 验证

- 每折留出一个完整数据集；测试集标签在特征选择、标准化和拟合中不可见。
- 仅在训练数据集计算 TF 效应并选择绝对训练 meta 效应最大的 20 个 TF；所有转换仅由训练数据拟合。
- 分类器：class-weight balanced、L2 logistic regression，固定 seed 1729。
- 主性能：把各留出数据集内部的 M-vs-rest concordance 合并为 stratified AUROC，避免数据集构成差异造成虚假 pooled AUROC。
- 区间：按“数据集×类别”分层患者 bootstrap 1,000 次。
- null：数据集内部置换标签 500 次；每次重跑训练内特征选择和整个 LODO。

LODO 必须同时满足 AUROC >=0.65、bootstrap 95% CI 下限 >0.55、置换经验 p <=0.05。

## 8. 最终停止规则

- `PASS`：第4节全部可判定条件、至少10个可重复 TF、至少5个 msVIPER 同向确认 TF、以及第7节三项 LODO 条件全部通过。
- `FAIL`：数据与网络可判定，但任一预设调控程序或 LODO 条件失败。结论仅为“M-vs-rest 不能支持当前 TNBC 式入口”；随后应冻结另一个 CRC 亚型协议，复用 TCGA 网络和患者 VIPER 活性，不在本结果上调阈值。
- `INDETERMINATE`：数据身份、样本量、网络或运行完整性不足。

判定生成后立即停止：本轮不运行 ATAC、药物敏感性、治疗推荐、空间验证或湿实验推断，也不自动尝试其他亚型。

## 9. 可复现性合同

云端保留原始下载、GDC manifest、TCGA 表达、ARACNe3 100 个 subnetworks、合并网络、regulon、患者 pseudobulk、VIPER/msVIPER 矩阵、meta/LODO 明细、日志和 SHA256。Git 仅保留协议、表型清单、代码、小型汇总表、收据和最终报告。任何失败和重试不得覆盖旧日志。
