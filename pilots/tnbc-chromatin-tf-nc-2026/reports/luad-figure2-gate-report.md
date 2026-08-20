# LUAD Figure 2 一级闸门最终报告

> **历史记录，已被取代。** 本报告记录 2026-08-14 采用额外“样本失败比例不超过 20%”规则时的 v1 决策。该规则不属于 TNBC 锚点文献的 Figure 2 入选逻辑，后经用户明确要求废止，并由原文式 v2 协议及完整 Figure 2 运行取代。当前结论请以 [`luad-figure2-final-report.md`](luad-figure2-final-report.md) 为准；下文仅为 Paper2Paper 的决策溯源，不能解释为当前项目状态。

日期：2026-08-14

最终判定：`FAIL_DATA`

后续状态：停止；不得进入 Figure 3

## 1. 执行范围与结论

本轮按冻结协议完成了以下工作：

- 补充 Figure 1 的临床与分子注释；
- 用既有 Figure 1 结果生成 Supplementary Figure 1–2 的核心内容；
- 对 Figure 2 所需的患者、PDX 和细胞系 ATAC 数据执行数据清单审计、下载、统一处理和硬质控；
- 将 Figure 2、Supplementary Figure 3、Supplementary Figure 4A–C 设为不可拆分的原子输出；
- 在数据充分性闸门失败后停止，不进入启动子开放、motif、稳健性或 Figure 3 分析。

最终结果为 `FAIL_DATA`。失败发生在细胞系 ATAC 数据充分性层面，因此本轮**没有检验 LUAD 是否存在 TNBC 式染色质支持的 TF 生物学信号**。

## 2. Figure 1 与 Supplementary Figure 1–2

Figure 1 原有数值结果被冻结，未重新选择 TF、样本或比较组。新增内容仅包括临床/分子注释和补图。

- TCGA-LUAD/LUSC：性别、年龄、分期、吸烟包年、表达亚型及关键驱动突变；
- PDMR PDX：性别、年龄、吸烟状态、取材部位、原发/转移与活检信息；
- DepMap：性别、年龄、原发/转移来源和取材部位；
- Supplementary Figure 1：各队列纳入、排除和最终分析样本流程；
- Supplementary Figure 2A–F：独立队列效应、TF 交集、各系统样本相关性及跨系统重复 TF 效应。

独立验证共 25 项检查，25 项全部通过；冻结的 Figure 1 数值文件哈希未发生变化。

## 3. Figure 2 数据审计

### 3.1 患者肿瘤

审计了 22 例 TCGA-LUAD 患者的公开 ATAC 处理矩阵。该部分达到预定样本清单要求。

### 3.2 LUAD 细胞系

19 个预先冻结的独立 LUAD 细胞系 ATAC 文库全部完成统一处理和硬质控：

- 合格：14/19；
- 不合格：5/19；
- 失败比例：26.32%；
- 冻结的最大允许失败比例：20%。

失败样本及原因如下：

| 细胞系 | 有效分析单位 | 峰数 | FRiP | 失败原因 |
|---|---:|---:|---:|---|
| H1648 | 2,088,866 | 9,413 | 0.0676 | 峰数未达阈值 |
| H1650 | 781,150 | 4,494 | 0.0564 | 有效单位和峰数均未达阈值 |
| H1819 | 855,212 | 6,754 | 0.1067 | 有效单位和峰数均未达阈值 |
| PC-14 | 1,424,761 | 8,270 | 0.0933 | 峰数未达阈值 |
| RERF-LC-KJ | 1,076,326 | 7,253 | 0.0910 | 峰数未达阈值 |

由于 5/19 = 26.32% 超过 20%，数据闸门被锁定为失败。

### 3.3 PDX

13 个 PDX 原始数据均已在云服务器完成下载。细胞系闸门锁定失败时：

- 1 个正在运行的 PDX 任务被有记录地中止；
- 其余 12 个 PDX 未再启动；
- 上述 1 个中止和 12 个未运行样本均**不计为 PDX 质控失败**；
- PDX 系统状态为 `NOT_FULLY_EVALUATED_AFTER_EARLY_STOP`。

这是预注册早停规则的执行结果，不能把未运行样本解释为低质量样本。

## 4. 原子输出状态

因为原始数据闸门在任何生物学信号检验之前失败，以下输出统一标记为：

`NOT_GENERATED_DUE_FAIL_DATA`

- Figure 2；
- Supplementary Figure 3；
- Supplementary Figure 4A–C。

相应的启动子开放、TF 锚点、JASPAR/CIS-BP motif 富集及额外稳健性分析均为 `NOT_EVALUATED_AFTER_DATA_FAILURE`。没有生成可被误解为正式结果的下游信号表或图。

## 5. 统计与代码验证

- Figure 1 注释与 Supplementary Figure 1–2：独立验证 25/25 通过；
- Figure 2 数据失败回执：独立验证 24/24 通过；
- HOMER 解析器使用真实 `knownResults` 结构完成了分母、加权背景计数、效应量和 q 值的烟雾测试；
- 早停汇总明确区分 `COMPLETED`、`ABORTED_AFTER_LOCKED_DATA_GATE` 与 `NOT_RUN_AFTER_LOCKED_DATA_GATE`；
- 数据下载和分析均在云服务器完成；本地仅同步代码、审计回执、小型结果表、图和报告，未同步 FASTQ、BAM 或峰文件。

## 6. 解释与停止规则

本次失败只能解释为：**当前冻结的 LUAD 三系统数据组合未通过预注册的数据充分性标准**。它不等于 LUAD 中不存在可重复的 TF 程序，也不等于 TNBC 方法在 LUAD 中生物学失败，因为生物学层尚未被评估。

按照冻结协议，本轮在此停止：

- 不改变质控阈值；
- 不事后删除失败样本以改变分母；
- 不替换或追加细胞系数据进行当轮补救；
- 不生成 Figure 2 原子图组；
- 不进入 Figure 3 或其余补图。

若以后重新启动，必须作为一个新版本闸门，预先登记新的数据来源、样本清单、阈值和版本号，不能覆盖本次 `FAIL_DATA` 结论。

## 7. 核心审计文件

- `execution/luad-figure2/audit/final_gate/figure2_final_verdict.json`
- `execution/luad-figure2/audit/raw_atac_qc/figure2_raw_atac_qc.tsv`
- `execution/luad-figure2/audit/raw_atac_qc/figure2_raw_atac_qc_receipt.json`
- `execution/luad-figure2/audit/final_gate/figure2_atomic_bundle_status.json`
- `execution/luad-figure2/audit/final_gate/figure2_data_failure_verification_receipt.json`
- `execution/luad-figure1/audit/figure1_annotations/figure1_annotations_independent_verification.json`
