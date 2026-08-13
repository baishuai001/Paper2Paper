# CRC M-vs-rest 调控程序一级闸门：预先冻结协议

冻结日期：2026-08-13

修订冻结：2026-08-13，在任何正式 ARACNe3 或 TF 活性结果产生前，根据 TNBC Code Ocean v1.0 作者代码审计纠正网络参数

真实数据执行位置：仅限云服务器 `/media/desk16/iy13202/projects/Paper2Paper`

表型清单：`analysis/phenotypes/m-vs-rest.json`

## 1. 闸门问题

CRC-atlas 作者预先定义的 M 型相对 B、T、desert 三型，是否在原发、未治疗、未富集样本的患者级 Cancer-cell pseudobulk 中，对应一个跨独立研究可重复的 CRC 特异 TF 活性程序？

本闸门只检验关联和跨研究可重复性，不声称免疫细胞驱动癌细胞、TF 有因果作用或可指导治疗。此前“四类表型能否由无监督聚类重建”的 FAIL 是额外 QC，不是 TNBC 原文进入 ARACNe3/VIPER 的条件，因此保留为辅助证据但不再硬停止本路线。M-vs-rest 标签不重新聚类、优化或调阈值。

## 2. TNBC 锚点代码审计与主方法选择

主方法保持为：大样本肿瘤 bulk RNA 构建癌种特异 ARACNe3 网络，`viper` 将网络转为 regulon，VIPER/aREA 计算逐患者 TF 活性，msVIPER 检验组间差异活性。

作者公开 Code Ocean 胶囊固定为 tag `v1.0`、commit `edf5314ce5b9ee0e2f88b2310e7c2df5619ad888`。审计结果如下：

- 论文 Methods 说明 ARACNe3、VIPER/aREA 和 TF FDR 0.01，但未规定 ARACNe3 子网络数或 VIPER `minsize`。
- 作者 `Figure1/05B-Run_ARACNe.sh` 不传 `-x`；固定 ARACNe3 版本的默认值是 1 个子网络、subsample 0.63212、子网络内 BH-FDR alpha 0.05、Maximum-Entropy pruning。
- 作者 `Figure1/06-VIPER_TNBC_NonTNBC.R` 明确使用 `subnets/subnet1_defaultid.tsv`，而不是 consolidated network；`viper` 与 `msviper` 均用 `minsize=1`，`ttestNull` 用 1,000 次置换。
- 作者脚本在 `header=TRUE` 后又删除第一行，实际会误删第一条真实边。本实现保留全部数据边，不复制该明显 bug。
- 作者 regulon 构建和 msVIPER 代码被注释，并直接载入预计算 RDS；因此胶囊提供算法意图但不是从原始输入可直接端到端执行的实现。本项目补齐这些步骤。

据此，本闸门的主网络固定为 1 个 ARACNe3 子网络并直接转 regulon；不使用 consolidated network，也不对 consolidated binomial p 值增加第二次 BH。固定 seed 1729 和 24 threads 仅用于确定性与性能，不改变作者算法。CollecTRI、DoRothEA、ULM 等先验网络不得替代本判定。

## 3. 输入与分析单位

### 3.1 CRC 特异网络

- GDC TCGA-COAD/READ `STAR - Counts`、`Primary Tumor`，使用唯一 TPM assay。
- 重复 aliquot 按 TCGA participant 对每个基因取均值；映射到同一 gene symbol 的 Ensembl 行取均值，与作者 helper 一致；删除空符号和非有限行。
- 不增加低表达或方差过滤；作者预处理代码没有这些步骤。
- PAN-GO：Methods 称 2,139 genes；补充表实际为 2,138 条非空 symbol 记录，含 79 条重复，得到 2,059 个唯一符号。使用可审计的 2,059 个唯一符号。
- ARACNe3：官方 commit `3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14`；1 个 subnetwork；subsample 0.63212；子网络内 FDR alpha 0.05；Maximum-Entropy pruning；seed 1729；24 threads。
- regulon：完整 `subnets/subnet1_crc.tsv` 三列边进入 `viper::aracne2regulon(..., format="3col")`；TF-target mode 由 TCGA CRC TPM 估计，与作者对 TCGA 的调用一致。

### 3.2 CRC-atlas 癌细胞表达

- H5AD 必须匹配 30,875,155,333 bytes 和 SHA256 `718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6`。
- 分析单位为患者 `donor_id`；只汇总作者标注 `Cancer cell` 的 `raw/X` 非负整数 counts。
- 范围与 M-vs-rest 分组完全来自冻结 manifest。同一患者的技术数据集可合并，但患者不得跨独立 `study_id`。
- 主阈值为至少 50 个 Cancer cells；20、100、200 为预设敏感性阈值。
- 患者内 counts 求和后转 `log2(CPM+1)`；保留至少 10%患者 CPM >=1 的基因。

## 4. 可判定条件

以下必须全部满足，否则为 `INDETERMINATE`，不能解释为生物学阴性：

1. TCGA CRC participant >=550，过滤后基因 >=10,000，基因与患者标识唯一；
2. PAN-GO 三种计数均记录，且至少 1,500 个唯一符号进入 ARACNe3；
3. H5AD、ARACNe3 commit、TNBC Code Ocean commit 和补充表身份均通过审计；
4. ARACNe3 成功生成恰好 1 个非空子网络，且该子网络至少含 500 个 regulator；
5. 至少 400 个转换后 regulon 在 CRC-atlas 中有至少 1 个可测 target；
6. 50 细胞阈值下 M >=30、non-M >=80；
7. 至少 3 个独立 `study_id` 各含 >=3 M 和 >=3 non-M；任一研究贡献的 M 不超过全部 M 的 60%；
8. 患者标签、研究归属、细胞数和 pseudobulk 列顺序一一对应，无跨研究患者泄漏。

## 5. VIPER 与 msVIPER

- 单患者活动：对全部合规 pseudobulk 运行 `viper(..., minsize=1, nes=TRUE)`；不显式传 `method`，因此使用 `viper` 默认的 `"none"`，与作者调用一致。
- 组间 signature：在每个信息性独立研究内计算 M-vs-rest 基因 Welch t statistic，再按有效样本量平方根加权合并。这是用户明确要求的跨研究迁移层。
- null：每个独立研究内置换 M 标签 1,000 次，每次重算完整 signature，并传给 `msviper(..., minsize=1)`。
- msVIPER 显著阈值：BH-FDR <=0.01，沿用锚点论文。

## 6. 跨研究 TF 程序

每个 VIPER TF 在各信息性研究计算 Hedges g（正值表示 M 更高），再作 REML 随机效应 meta。因仅 3 个信息性研究，使用 modified Knapp-Hartung 标准误和 t 检验。一个“可重复 TF”须同时满足：

- 至少 3 个信息性研究；
- meta BH-FDR <=0.05；
- |meta Hedges g| >=0.50；
- 至少 75%研究同向；3 个研究时即 3/3 同向；
- I² <=75%。

“可重复 TF 程序”定义为至少 10 个可重复 TF，其中至少 5 个同时满足 msVIPER FDR <=0.01 且方向一致。

## 7. 无泄漏 LOSO 验证

- 每折留出一个完整 `study_id`；测试标签不参与特征选择、标准化或拟合。
- 每折仅用训练研究选择绝对训练 meta 效应最大的 20 个 TF；分类器为 class-balanced L2 logistic regression，seed 1729。
- 主性能为按研究内 M-vs-rest concordance 加权的 stratified AUROC。
- 95% CI：按“研究×类别”分层患者 bootstrap 1,000 次。
- null：研究内置换标签 500 次；每次重跑训练内特征选择和完整 LOSO。

AUROC >=0.65、bootstrap 95% CI 下限 >0.55、置换经验 p <=0.05 预先定义为“强 LOSO 外推支持”。LOSO 必须完整运行和报告，但它是额外的预测泛化诊断，不是 TNBC 锚点方法，也不能否决已经由独立研究效应 meta、3/3 同向和 msVIPER 确认建立的可重复 TF 程序。

## 8. 最终停止规则

- `PASS`：全部可判定条件通过，且至少 10 个可重复 TF、其中至少 5 个获 msVIPER 同向确认。LOSO 是否达到“强支持”单独报告，不改变 PASS。
- `FAIL`：数据与网络可判定，但上述可重复 TF 程序条件失败。只表示 M-vs-rest 不支持当前入口；随后应另行冻结其他 CRC 亚型，复用 TCGA 网络和患者 VIPER 矩阵，不能在本结果上调阈值。
- `INDETERMINATE`：数据身份、样本量、网络或执行完整性不足。

判定生成后立即停止：本轮不运行 ATAC、药敏、治疗推荐、空间验证、湿实验推断，也不自动尝试其他亚型。

## 9. 可复现性

云端保留原始下载、GDC manifest、TCGA 表达、单个 ARACNe3 子网络与合并输出、regulon、pseudobulk、VIPER/msVIPER、meta/LOSO 明细、日志和 SHA256。Git 只保留协议、表型 manifest、代码、小型结果、收据和报告。失败与重试日志不得覆盖。
