# 第零阶段执行规范：CRC 基础图谱、恶性 CNA 与 cNMF

冻结日期：2026-08-10  
状态：方法与验收合同已冻结；真实对象 P0.0–P0.6 第一轮诊断已完成，整体科学验收未通过。  
作用：为 CRC 是否能够承接 NPJ HCC 锚点的恶性程序主线提供可信的基础层。

本文件是运行前合同，不是运行结果。真实文件、实际数值、代码失败与修复、分阶段判定见
[`../reports/phase-zero-real-run-report.md`](../reports/phase-zero-real-run-report.md)。合同冻结后
允许修复会导致错误或不可执行的实现，但不得根据结果反向改写科学退出门。

## 一、这一阶段要回答什么

第零阶段不是为了生成一张漂亮 UMAP，而是回答五个会决定后续全部 Figure 是否成立的问题：

1. 公开 CRC 对象与 Supplementary Table S1 能否建立不丢失患者身份的设计表？
2. QC、doublet、整合和注释后，关键状态是否跨患者、跨数据集存在，而不是批次产物？
3. 哪些上皮细胞有足够证据称为恶性，哪些必须保留为 unresolved？
4. 恶性细胞的表达程序是否能在合理 K、随机种子、患者和数据集扰动下复现？
5. 是否存在至少一组冻结程序，可以进入后续生态、空间、临床和扰动分析？

若上述问题不通过，后续分析应缩减或停止，而不是换参数直到出现预期故事。

## 二、冻结的来源和实现边界

### 论文与数据来源

- CRC Cancer Cell 正文、补充方法、Supplementary Table S1；
- CRC-atlas 作者仓库固定提交
  `82c15ecc2ea36baa0c4cffe8818bf58e21dcfc48`；
- CRC-atlas 论文版/Zenodo/CELLxGENE 对象分别记录，不默认同版本；
- NPJ HCC 正文和补充材料只提供锚点比较与诊断 Figure，不作为可执行代码来源。

### 官方与同类实现

- scvi-tools/scANVI：整合和标签传播；
- scib-metrics：批次去除与生物保留的成对评估；
- CopyKAT + SCEVAN：恶性 CNA 主分析与敏感性；
- `dylkot/cNMF`：cNMF 主算法和非负批次校正预处理；
- `tiroshlab/3ca`：跨 rank、跨患者程序复现与合并原则。

所有实际版本、commit、容器 digest 和许可证在首次运行前写入环境收据；这里不以“最新版”
代替冻结版本。

## 三、阶段与退出门

```text
P0.0 资源收据
  → P0.1 患者—样本设计表
  → P0.2 QC、ambient RNA 与 doublet 收据
  → P0.3 整合与批次/生物保留验收
  → P0.4 marker 与标签复核
  → P0.5 逐样本 CNA 恶性判定
  → P0.6 cNMF K、稳定性、患者混合与程序合并
  → P0-GATE 第零阶段科学审查
```

前一阶段未通过时，后一阶段只能做标明目的的诊断运行，不得生成正式结论。

## 四、P0.0：资源收据

### 输入

- 论文版 CRC atlas H5AD 或最小必需对象；
- Supplementary Table S1；
- 作者代码仓库；
- 所选官方方法仓库/软件包。

### 必须记录

- URL/DOI/accession、访问日期、版本或 record ID；
- 文件名、字节数、SHA256、压缩/解压后大小；
- AnnData 维度、`X/raw/layers`、obs/var 字段和数据类型；
- 原始非负计数的真实位置；
- 代码 commit/release、许可证、环境和入口；
- 预计峰值磁盘、内存、CPU/GPU 和临时目录。

### 退出条件

- 真实文件可读取；
- 非负计数层可追溯；
- 患者、样本、研究、组织/状态至少有可连接字段；
- 资源足够完成下一阶段，而不只是网页声称公开。

如果对象版本和论文样本数不一致，建立版本差异表，不静默选择更方便的版本。

## 五、P0.1：患者—样本设计表

以 Supplementary Table S1 为主表，以 H5AD `obs` 为实际计算表，进行双向连接。

### 输出字段

`patient_id, sample_id, study_id, dataset_id, platform, tissue, site, state, treatment_status,
regimen, timepoint, response, recist, msi, stage, paired_group, object_version, cell_count,
join_status, exclusion_reason`。

### 检查

- patient/sample ID 在各自定义层级唯一；
- Table S1→H5AD 和 H5AD→Table S1 的未匹配 ID 数；
- 同一患者跨样本、部位、时间和数据集的重复；
- 状态是否与数据集/平台完全重合；
- 每个决定性比较的独立患者数，而不是细胞数；
- 治疗、response 和 RECIST 缺失及定义是否可比较。

### 退出条件

所有进入后续分析的细胞均能追溯到患者和样本；未知身份的细胞不得进入患者级推断。

## 六、P0.2：QC、环境 RNA 与 doublet

### 运行单位

所有阈值、ambient RNA 和 doublet 均按样本运行；跨样本全局阈值只作敏感性比较。

### 主规则

CRC 作者阈值作为预注册起点：转录本 >400、基因 >100、线粒体 <50%，再按样本 3 MAD。
同时保存连续分布，检查其对组织、平台和稀有状态的差异影响。FASTQ/empty-droplet 信息
可得时使用 scAR 或相应已验证环境 RNA 方法；不可得时明确标记 `not_run`。

SOLO/其他 doublet 方法失败必须使该样本状态为 `unknown` 或 `failed`；禁止将缺失结果填为
singlet。对跨样本 pool、极小样本和平台不兼容情况单独处理。

### 每样本收据

- QC 前后细胞数、基因/UMI/线粒体分布；
- 各规则单独和联合删除数；
- ambient RNA 方法、输入和删除/校正量；
- doublet 方法、阈值、预测比例、失败状态；
- 关键 marker/细胞类型在 QC 前后的变化；
- 参数敏感性和日志。

### 退出条件

不存在静默失败；决定性状态不能由极低质量或 doublet 高风险样本独占；合理阈值扰动不应
消除整个候选状态。未完成 doublet 的样本可保留作敏感性，但不得混入主分析而不披露。

## 七、P0.3：整合与批次/生物保留验收

### 主实现

以 CRC 作者 scVI/scANVI 结构为 adapter，使用可追溯 count layer、冻结的 HVG 规则和样本/
数据集 batch key。训练/验证/留出样本和随机种子写入 manifest。

### 必须同时优化的两个目标

| 目标 | 指标/检查 |
|---|---|
| 去除不期望技术差异 | iLISI、batch ASW、kBET、graph connectivity，逐数据集/患者报告 |
| 保留真实生物差异 | cell-type cLISI/ASW、marker 保留、稀有状态召回、疾病/部位信号、留出标签转移 |

作者 UMAP 只作视觉参照，不是退出门。至少运行两个合理随机种子；如 scANVI 与 Harmony
得出相反生物结论，保留冲突并分析原因，而不是选择最符合故事的方法。

### 退出条件

技术混合得到改善且已知生物结构没有明显崩溃；结果不是由单一数据集或训练 seed 决定；
未通过的细胞群被标为不稳定，不进入主线命名。

## 八、P0.4：marker 与标签复核

作者标签保存在 `author_label`；目标复核标签保存在 `reviewed_label`，不能覆盖原列。

### 每个目标类型的证据

- 正 marker、负 marker、相邻混淆类型；
- marker 在不同患者/数据集中的覆盖，而不只全局均值；
- 训练集与留出数据集的性能；
- doublet/低质量/组织污染比例；
- cluster-level 和 cell-level 不确定性；
- 人工改名的规则、日期和证据。

### 退出条件

标签能够跨多个患者和数据集复核；无法区分的细胞保留 `unresolved`。上皮候选进入 P0.5，
恶性身份不能只由作者的 `tumor cell` 标签决定。

## 九、P0.5：逐样本 CNA 与恶性判定

### 输入

- 经 P0.2 通过的原始非负 UMI；
- patient/sample/tissue；
- reviewed epithelial label；
- 明确的正常参考候选。

### 主流程

1. 每个患者/样本独立运行 CopyKAT；固定主版本、参数、seed 和正常参考策略。
2. 使用 SCEVAN 对同一输入独立运行敏感性。
3. 记录 CopyKAT 的 `aneuploid/diploid/not.defined`、SCEVAN 的 `tumor/normal/filtered` 和连续
   CNA 证据；不把任一工具的类别名直接改写为病理真值。
4. 用冻结规则生成 `two_method_malignancy_support/one_method_support/method_unresolved/
   without_CNA_support`。`without_CNA_support` 不是非恶性证明；只有独立病理、DNA 或明确正常
   组织身份才能另行支持 nonmalignant。
5. 以锚点 CopyKAT 1.1.0 或另一合理参数做版本/参数敏感性，不与主版本结果混写。

### CNA 每样本验收表

`input_cells, input_genes, reference_cells, copykat_defined_fraction, scevan_defined_fraction,
method_concordance, two_method_malignancy_support, one_method_support, method_unresolved,
normal_reference_unexpected_support, cna_burden_summary,
qc_flag, exclusion_reason`。

### 失败与降级

- 正常参考不足、有效判定比例过低、方法严重冲突：该样本不进入主发现集；
- 近二倍体肿瘤风险：不得将 diploid 自动标成正常；
- 某数据集全部被一种方法判成恶性而其他数据集相反：优先调查平台/参考效应；
- CNA 仅作为恶性证据之一，不证明具体驱动事件或克隆演化。

## 十、P0.6：cNMF 程序发现与复现

### 输入冻结

- 主诊断集：作者 Cancer 候选中同时得到 CopyKAT aneuploid 与 SCEVAN tumor 支持的细胞；
- 敏感性集：加入有额外证据的 unresolved；
- 统一基因 ID 和基因宇宙；
- 官方 cNMF `Preprocess` 产生非负批次校正 count matrix；
- 患者、样本、数据集和状态 manifest；
- 每患者最大贡献细胞数，防止大患者支配。

### 运行设计

- pooled donor-balanced discovery；
- 细胞数足够时，逐患者或逐数据集发现；
- K 范围在富集/临床结果前冻结，初始 5–20；
- 每个 K 100–200 次随机初始化；
- 预留患者和数据集作 holdout；
- 所有 seed、抽样和失败记录保留。

### K 选择算法

1. 对每个 K 输出 stability、reconstruction error、冗余和可匹配程序数；
2. 识别 error knee/plateau 和 stability plateau；
3. 在通过误差平台的稳定区选择最小 K；
4. K0±1/2，条件允许时 ±4 做匹配敏感性；
5. 无稳定平台则结论为“不足以冻结 K”，不使用通路解释强行决定。

### 稳定性验收

- seed、平衡抽样、患者留出、数据集留出四类扰动；
- cosine similarity + top-gene Jaccard；
- Hungarian 一对一匹配；
- 与置换程序零分布比较，冻结相似度阈值；
- 报告匹配失败、拆分和合并，不只展示成功程序。

### 患者/数据集混合

按每患者等量细胞或患者级 usage 计算贡献。输出 normalized Shannon entropy、effective donor
number、最低贡献患者数、患者/数据集各自指标、平衡抽样置信区间和分层置换零分布。混合
指标是诊断量；跨患者 holdout 复现才是主要证据。

### 程序合并

1. 先要求同一单元跨相邻 K/rank 重现；
2. 再要求匹配至少另一个独立患者，优先另一个数据集；
3. 统一基因宇宙、患者等权谱；
4. cosine condensed distance + average/complete linkage，或预先冻结的图聚类；
5. cut 由置换零分布和敏感性确定；
6. singleton 不强制合并，标为 context-specific；
7. 输出完整成员关系、相似度、top genes、来源患者/数据集和技术标记；
8. 在查看空间、临床和扰动结果前冻结 MP 定义。

### 退出条件

至少一个非技术程序在多个 seed/抽样、至少两个患者且最好两个数据集中复现，并在合理 K
和 CNA 支持集合敏感性下保持。即使通过，名称也应是“CNA 支持的作者 Cancer 候选程序”，
除非另有 DNA/病理证据。否则完整“恶性程序”路线停止或缩减为描述性 atlas 分析。

## 十一、必须产生的 source tables 和诊断图

### Source tables

1. `P0_resource_manifest.tsv`
2. `P0_patient_sample_design.tsv`
3. `P0_qc_by_sample.tsv`
4. `P0_doublet_receipt.tsv`
5. `P0_integration_metrics.tsv`
6. `P0_annotation_evidence.tsv`
7. `P0_malignancy_calls.tsv`
8. `P0_cna_sample_qc.tsv`
9. `P0_cnmf_k_selection.tsv`
10. `P0_cnmf_stability.tsv`
11. `P0_cnmf_patient_mixing.tsv`
12. `P0_program_merge_membership.tsv`

### 诊断图

- 每样本 QC 与纳入漏斗；
- 数据集×状态×患者覆盖；
- 整合前后批次和生物指标；
- 按患者/数据集分层 marker 图；
- CopyKAT/SCEVAN 一致性与 CNA 热图；
- K stability/error/冗余；
- seed/患者/数据集 holdout 程序匹配；
- 患者和数据集混合及零分布；
- 合并前后程序树/图和完整成员表。

每张图必须由对应 source table 非交互生成；图不能替代 source table。

## 十二、代码可靠性验收

每个步骤至少需要：

- 固定来源、版本/commit、许可证和环境；
- 明确函数/CLI、输入输出、默认参数和随机种子；
- 最小合成 fixture：ID 错配、负值 cNMF 输入、doublet 失败、无正常参考、单患者程序等必须
  失败关闭；
- 小型真实输入 smoke；
- 患者身份、无数据泄漏、矩阵非负、统计单位和结果行数等科学不变量测试；
- 完整命令、退出码、stdout/stderr、运行时间和资源占用；
- 同一冻结输入可重复生成一致 source table；
- 方法/版本/参数敏感性和阴性结果保留。

作者代码、官方包或 AI 重建代码都执行同一验收；来源名气不能跳过测试。

## 十三、第零阶段最终审查单

只有以下全部满足，才写“基础分析已验收”：

- [ ] 真实资源版本、校验和和计数层已确认；
- [ ] 患者—样本双向连接完成，决定性组人数明确；
- [ ] QC/doublet 无静默失败；
- [ ] 批次去除和生物保留同时通过；
- [ ] marker/标签跨患者和数据集复核；
- [ ] CNA 逐样本运行并保留 unresolved；
- [ ] cNMF K 有预先定义的选择证据；
- [ ] seed、患者和数据集稳定性通过；
- [ ] 患者混合不由大样本机械造成；
- [ ] 程序合并可追溯且在下游结果前冻结；
- [ ] source tables、图、日志、环境和失败记录齐全；
- [ ] AI 给出继续、缩减或停止的判断及依据。

这些勾选项不得因文献中已有相似图、软件安装成功或单次命令退出码为 0 而视为通过。当前
P0.0 已通过，P0.1/P0.2 受限通过，P0.3/P0.4 仅完成发布对象诊断；P0.5 的 12 样本面板与
P0.6 的 CNA 支持输入、pooled seed 和 A/B 留出已执行完成，但分别因正常参考异常支持、
K/来源稳定性和敏感性不足而未通过科学验收。`P0_program_merge_membership.tsv` 不应在 K
未冻结时生成；修正后的评估只保留逐 K 的直接留出匹配。最终数值与判断见真实运行报告。
