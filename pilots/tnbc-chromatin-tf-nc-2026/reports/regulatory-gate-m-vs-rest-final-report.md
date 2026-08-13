# CRC M-vs-rest 调控程序一级闸门：最终报告

判定日期：2026-08-13

执行位置：云服务器 `/media/desk16/iy13202/projects/Paper2Paper`

冻结表型：`crc_atlas_M_vs_rest`

正式判定：**FAIL，已停止**

## 1. 结论

本次一级闸门可以作出科学判定，不是数据或软件失败：16 项数据、样本量、网络和执行完整性条件全部通过。但在预先冻结的 M-vs-rest 对比中：

- 跨独立研究可重复 TF：0 个，要求至少 10 个；
- 其中获 msVIPER FDR ≤0.01 且方向一致确认的 TF：0 个，要求至少 5 个；
- 20、100、200 个 Cancer cells 三个敏感性阈值下，可重复 TF 仍均为 0。

因此，按照用户指定的唯一硬科学停止规则——“只有没有可重复 TF 程序才判入口失败”——本路线判为 **FAIL**。这表示：当前 CRC-atlas 的作者 M 型相对 non-M 型，未能支持一个满足本闸门标准的、癌细胞层面的 CRC 特异可重复 TF 活性程序。它不等于 M 型没有生物学含义，也不证明任何单个 nominal TF 无作用。

判定后未自动切换其他 CRC 亚型，未运行 ATAC、药物敏感性、治疗选择、空间验证或湿实验推断。

## 2. 一级闸门实际检验的问题

在原发、术前未治疗、未做细胞富集的 CRC 样本中，将作者预定义的 M 型与 B、T、desert 合并的 non-M 型比较：

1. 用 TCGA-COAD/READ bulk RNA 构建 CRC 特异 ARACNe3 网络；
2. 用该网络和 VIPER 计算 CRC-atlas 每位患者 Cancer-cell pseudobulk 的 TF 活性；
3. 在独立研究内估计 M-vs-rest 效应，再进行随机效应 meta；
4. 用 1,000 次研究内标签置换的 msVIPER 作正交确认；
5. 用留一研究验证（LOSO）评价额外的预测外推能力。

此前“四类免疫表型能否由另一个无监督聚类方案重建”的 FAIL 仅保留为辅助证据，不进入本次硬判定。CollecTRI＋ULM 也未用于本次正式结果。

## 3. 与 TNBC 锚点方法的关系

| 环节 | 本次实现 | 性质 |
|---|---|---|
| 调控网络 | ARACNe3，1 个 subnetwork，subsample 0.63212，子网内 FDR alpha 0.05，Maximum-Entropy pruning | 与 TNBC Code Ocean v1.0 作者命令一致 |
| 网络来源 | 完整 `subnets/subnet1_crc.tsv` | 与作者代码明确使用 `subnet1` 一致；未用 consensus 作主网络 |
| regulon | `viper::aracne2regulon(..., format="3col")`，原始 CRC TPM 估计 mode | 与作者 TCGA 调用一致 |
| 单患者活性 | `viper(..., minsize=1, nes=TRUE)`，不显式传 `method`，因此为默认 `none` | 与作者调用一致 |
| 组间确认 | `msviper(..., minsize=1)`，1,000 次 null | 与作者阈值和置换数一致 |
| 癌种数据 | TCGA-COAD/READ GDC STAR TPM；CRC-atlas Cancer-cell pseudobulk | 必要的跨癌种替换；TNBC 锚点使用的是其 TNBC 数据与 Kallisto TPM |
| 跨队列统计 | study_id 分层、REML＋modified Knapp–Hartung、LOSO | 用户明确要求的 CRC 多研究迁移层，不是 TNBC 原文步骤 |

未使用早期试行中的 100 个子网络、第二次 consensus BH、`minsize=25`、低表达过滤、VIPER `scale/auto`、CollecTRI 或 ULM。作者脚本在 `header=TRUE` 后误删第一条真实边的明显错误也未复制。

## 4. 数据与可判定性

### 4.1 输入身份

- CRC-atlas H5AD：30,875,155,333 bytes，3,790,266 × 28,127；完整 SHA256 为 `718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6`。
- TNBC Supplementary Data 2：25,906,587 bytes；SHA256 为 `D43C305BAA7E316AF1A0227708DC5C177E20A9107C3872AFA5E0BF4B379227DD`。
- ARACNe3：官方 commit `3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14`。
- TNBC Code Ocean：tag `v1.0`，commit `edf5314ce5b9ee0e2f88b2310e7c2df5619ad888`。
- 正式输入审计还记录了表型 manifest 与 13 个实际执行源文件的 SHA256。

第一次正式调用在 ARACNe 启动前因跨平台 shell 可执行位缺失而停止，没有生成网络或 TF 结果。错误日志完整保留。修复为显式 `bash run_aracne3.sh` 后重试；重试仅复用刚完成的 H5AD 全哈希收据，并验证 H5AD 的 bytes、mtime、预期 SHA256、旧收据状态及旧收据 SHA256，未跳过身份核验。

### 4.2 PAN-GO 审计

- Methods 声称 2,139 genes；
- 补充表实际为 2,138 条非空 symbol 记录；
- 其中 79 条为重复记录，得到 2,059 个唯一 symbol；
- 2,051 个 symbol 存在于 TCGA CRC 矩阵并进入 ARACNe3，8 个缺失。

三种计数均保留，没有将 2,059 静默写成 2,139。

### 4.3 TCGA CRC 网络输入

- 647 个 GDC primary-tumor 文件：COAD 481、READ 166；
- 624 位 participant，13 位有多个 aliquot；
- 59,427 个唯一 gene symbol；
- 重复 symbol 按锚点代码取均值，重复 aliquot 按 participant 取均值；
- 不增加低表达过滤；0 个非有限 gene row 被删除。

### 4.4 CRC-atlas 患者级表达

- 213 位有作者标签的患者均生成 Cancer-cell pseudobulk；
- 汇总 244,196 个作者标注的 `Cancer cell`；
- 28,088 个唯一 raw gene symbol，过滤后 16,306 个基因；
- raw/X 通过非负整数验证；
- 主阈值 ≥50 个 Cancer cells：206 位，M=46、non-M=160，共来自 16 个研究；
- 满足每组至少 3 人的独立研究只有 3 个，构成实际 meta/msVIPER/LOSO 推断集：

| 独立研究 | M | non-M | LOSO AUROC |
|---|---:|---:|---:|
| Joanito 2022 Nat Genet | 5 | 22 | 0.509 |
| Lee 2020 Nat Genet | 10 | 22 | 0.673 |
| Pelka 2021 Cell | 18 | 42 | 0.665 |

推断集共 119 位患者（M=33、non-M=86）。技术数据集不是独立重复；同一研究内的 fresh/frozen 或 10x v2/v3 不被错误计作多个独立队列。

全部 16 项数据/网络条件均为 true，所以结果是可解释的 `FAIL`，不是 `INDETERMINATE`。

## 5. 网络与 TF 活性结果

- ARACNe3 成功生成恰好 1 个子网络，日志为 `SUCCESS!`；
- `subnet1_crc.tsv` 含 268,758 条边、2,051 个 regulator；
- 转换后 2,045 个 regulon；
- 其中 2,028 个在 CRC-atlas pseudobulk 中至少有 1 个可测 target；
- 213 位患者 × 2,028 个 TF 的 VIPER 活性矩阵完成；
- msVIPER 主分析：206 位主阈值患者、3 个独立研究、1,000 次研究内置换。

网络、可执行文件、输入表达和 regulator 列表均有 SHA256 收据。

## 6. 硬科学结果

### 6.1 跨研究 meta

2,028 个 TF 均在 3 个信息性研究中具有可用效应。独立复核得到：

| 条件 | TF 数 |
|---|---:|
| nominal meta p ≤0.05 | 7 |
| meta BH-FDR ≤0.05 | 0 |
| \|meta Hedges g\| ≥0.50 | 285 |
| 3/3 研究同向 | 823 |
| I² ≤75% | 1,999 |
| 除多重校正外的联合条件全部通过 | 239 |
| 全部可重复条件通过 | **0** |

最小 nominal meta p 来自 PRRX1：meta g=1.103、p=0.0377，但 FDR=0.866。其余 nominal TF 也没有通过多重校正。因此不能把 PRRX1、GLI3、AEBP1、RAX、HOPX、TBX18 或 MEOX2 称为可重复程序。

### 6.2 msVIPER

- nominal p ≤0.05：134 个 TF；
- BH-FDR ≤0.01：**0 个 TF**；
- 最小 nominal p 为 MEF2D：NES=3.463、p=0.000534、FDR=0.438。

因此不存在可用于确认 meta 程序的 FDR 显著 TF。SMAD3 等 nominal 结果也只能作为未确认信号，不能进入核心发现。

### 6.3 癌细胞数阈值敏感性

| 最小 Cancer cells | 患者 | M | non-M | 信息性研究 | 可重复 TF | LOSO AUROC |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 210 | 47 | 163 | 3 | **0** | 0.651 |
| 50（主） | 206 | 46 | 160 | 3 | **0** | 0.651 |
| 100 | 204 | 45 | 159 | 3 | **0** | 0.651 |
| 200 | 196 | 44 | 152 | 3 | **0** | 0.641 |

失败不依赖单一癌细胞数阈值。

## 7. LOSO：有弱外推信号，但未达到“强支持”

- stratified AUROC：0.651；
- 按研究×类别分层 bootstrap 1,000 次：95% CI 0.528–0.762；
- 研究内标签置换 500 次：经验 p=0.01996；
- 预定义强支持要求 AUROC ≥0.65、CI 下限 >0.55、置换 p≤0.05；其中 CI 条件失败。

所以 LOSO 显示一个边界性的、并不稳定的可预测信号；Joanito 留出折接近随机。更重要的是，LOSO 在冻结协议中只是额外泛化诊断，不能建立或否定硬 TF 程序。最终 `FAIL` 来自“0 个可重复 TF、0 个 msVIPER 确认 TF”，不是来自 LOSO。

## 8. 独立统计复核

判定后运行了只读的第二条算术审计路径，不改变原始表或 verdict。18 项检查全部通过：

- meta 与 msVIPER 的 BH-FDR 从原始 p 值重新计算后与保存值一致；最大绝对差分别为 `2.22e-16` 和 `8.88e-15`；
- 重新应用全部联合条件仍得到 0 个可重复 TF；
- 重新计算三研究样本量、每折 AUROC、stratified AUROC、bootstrap 分位数与置换 p，均与主收据一致；
- 三个敏感性阈值均重新确认为 0 个可重复 TF；
- 自动判定确认为全部 data conditions 通过、两个 scientific conditions 失败、`verdict=FAIL`、`stop_now=true`。

## 9. 判定边界

本结果支持的表述是：

> 在当前冻结队列、CRC 特异 ARACNe3 网络、VIPER/msVIPER 和跨研究 mKH meta 标准下，CRC-atlas M-vs-rest 未建立可重复的癌细胞 TF 活性程序，因此不能作为继续迁移 TNBC 染色质—调控—治疗框架的当前入口。

本结果不支持以下表述：

- “M 型 CRC 没有生物学意义”；
- “所有 nominal TF 都是假的”；
- “CRC 不可能使用 TNBC 框架”；
- “可以据此推荐治疗”；
- “四类作者表型错误”。

只有 3 个信息性独立研究，modified Knapp–Hartung 的自由度为 2，再对 2,028 个 TF 做多重校正，功效有限。因此这是一个严格、低假阳性优先的入口失败，不是对 TF 生物学不存在的证明。

## 10. 停止与后续可复用边界

判定生成后已停止。若由人类另行决定更换 CRC 亚型：

- `REGULATORY_GATE_MANIFEST` 可选择新的冻结表型；
- `REGULATORY_GATE_WORK` 隔离新路线结果；
- `REGULATORY_GATE_SHARED_NETWORK_OUTPUTS` 可复用 TCGA、PAN-GO 和 ARACNe3 网络；
- 主阈值、case/control、study 字段和判定 route 均从 manifest 读取；
- 当新表型改变队列时，患者选择、VIPER/msVIPER、meta、LOSO 与判定会重新运行；不会在本次结果上调阈值。

本轮没有选择或运行下一亚型。

## 11. 证据位置

Git 小型证据：`reports/receipts/regulatory-gate-m-vs-rest/`。

云端完整工作目录：`/media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/regulatory-gate-M-vs-rest/`。

正式流水线日志：云端 `logs/14_final_retry.log`。

正式判定 SHA256：`22de5500631b94cb0b47829698350e6d90f2fa996e9f779bba6fe7aa0fb8cb8c`。

独立复核：`post_decision_validation.json`，SHA256 `697b14dd850033c7b3bf2a89e026235948c9408f104bcc8388b2aba8fead9017`。本地小型证据副本与云端原件的 SHA256 已逐项比对一致。

来源入口：[TNBC 论文](https://www.nature.com/articles/s41467-026-76385-8)、[TNBC Code Ocean](https://codeocean.com/capsule/7227095/tree/v1)、[CRC-atlas 代码仓库](https://github.com/icbi-lab/crc-atlas)、[ARACNe3 官方仓库](https://github.com/califano-lab/ARACNe3)、[Bioconductor viper](https://bioconductor.org/packages/viper)。
