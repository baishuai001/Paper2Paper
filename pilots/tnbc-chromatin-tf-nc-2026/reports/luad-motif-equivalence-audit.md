# LUAD motif 技术等价性与统一背景重算报告

日期：2026-08-17  
状态：`COMPLETE`  
运行位置：云服务器；本地仅同步代码、协议、审计回执和小型结果表

## 1. 结论

本次审计得到三个必须同时保留的结论：

1. **TNBC 阳性对照恢复成功。** 从作者 Code Ocean 发布的逐样本 HOMER
   `knownResults.txt` 起步，作者未修改的 Figure2/06B 脚本生成的 JASPAR、
   CIS-BP 合并表与发布表逐字节相同；按作者 07A/08 规则重算后，56 个 motif
   可检验 TF 中 43 个至少获一个系统支持、31 个获三系统支持，TF 交集差异
   为 0，LOR 汇总在 `1e-12` 容差下差异为 0。
2. **严格照作者实际数据库选择规则，LUAD 三系统 motif 仍为 0。** 20 个
   motif 可检验 LUAD HC-TF 中，没有 TF 同时通过患者、PDX、细胞系。NFATC4
   和 XBP1 通过 PDX＋细胞系，但患者不通过。
3. **这个 0 不能解释为干净的生物学缺失。** 严格作者规则在预定义肺腺谱系
   阳性对照上仅有 PDX 的 CEBPB 通过，患者和细胞系均为 0；允许所有已编目的
   候选 motif 后，47-model 候选库和完整 1,944-model 库均稳定恢复三端
   **ELF3**，且两种分析的 20-TF 支持矩阵完全一致。

因此，当前数据**不支持 TNBC 式“广泛的三端染色质支持 TF 程序”**；但也不能
写成“LUAD 完全没有三端共同 motif”。更准确的结论是：严格作者实现得到 0，
而对 motif 模型选择不那么排他的两种预设敏感性分析均得到一个稳健的 ELF3
轴。结果对 JASPAR/CIS-BP 模型选择敏感。

## 2. TNBC 阳性对照的恢复范围

| 项目 | 结果 |
|---|---:|
| 固定系统数 | TCGA 7；PDX 38；细胞系 26 |
| 半数阈值 | 4；19；13 |
| motif 可检验 TF | 56 |
| 至少一个系统支持 | 43 |
| 三系统支持 | 31 |
| 发布交集差异 | 0 |
| LOR 汇总差异（容差 `1e-12`） | 0 |
| 06B 合并表 | JASPAR、CIS-BP 均与发布文件 SHA256 相同 |

这里恢复的是作者发布的**逐样本 HOMER 输出至最终三端交集**这一完整 motif
统计层，不是从 TNBC ATAC 峰重新运行 HOMER。Code Ocean capsule 没有提供
全部 query peak BED，PDX 原始 ATAC 还受控，因此不能诚实声称完成了 TNBC
raw ATAC → peak → HOMER 的全链重算。LUAD 则从本项目的峰文件实际运行了
峰到 HOMER 的链条。

## 3. LUAD 输入统一

| 项目 | 冻结值 |
|---|---|
| 样本 | 患者 22；PDX 13；细胞系 19 |
| 峰数可匹配 | 患者 22；PDX 10；细胞系 19 |
| 每个可匹配样本目标峰 | 3,584 |
| 峰宽 | 精确 200 bp |
| 目标 GC | 所有样本使用相同 5% GC-bin 配额；均值跨度 0.00234 |
| 共同背景 | 28,672 个 TCGA-LUSC 可及峰；200 bp；相同 GC-bin 比例 |
| blacklist | 同一 ENCODE hg38 blacklist |
| 随机种子 | 20260817 |
| HOMER | v5.1 |
| motif | 候选 JASPAR 21＋CIS-BP 26；完整联合库 1,944 |

3 个 PDX 因统一处理后峰数不足而不可进行峰数匹配：MGH1065、MGH1157、
MGH9243。主表仍按原预定义 13 例计算半数阈值 7；同时按 10 个可匹配 PDX
计算阈值 5。两种分母的三端结论完全相同。

所有样本复用同一个背景 BED 和 SHA。HOMER 会自动删除与各目标集重叠的背景
并做内部归一化，所以输出中的有效背景分母为 25,757–27,647，而不是机械等于
28,672；目标有效分母也在 3,581–3,587 间轻微变化。解析器直接读取每个
`knownResults.txt` 表头中的实际分母，没有用输入 BED 行数替代。

## 4. 三层 LUAD motif 结果

| 层级 | 数据库规则 | 三端 HC-TF | 阳性对照支持：患者/PDX/细胞系 |
|---|---|---:|---:|
| 作者等价主分析 | JASPAR、CIS-BP 分开校正；每 TF 有 JASPAR 时优先 JASPAR，否则 CIS-BP | 0/20 | 0/1/0 |
| 候选库敏感性 | 47 个候选模型联合校正；允许任一已编目来源 | 1/20：ELF3 | 1/2/1 |
| 完整库敏感性 | 1,944 个模型联合校正；允许任一已编目来源 | 1/20：ELF3 | 1/2/1 |

候选 47 与完整 1,944 库的 HC-TF 系统支持矩阵逐项相同。共同模式为：

- ELF3：患者、PDX、细胞系均支持；
- FOXA3：仅患者支持；
- NFATC4、XBP1：PDX 和细胞系支持，患者不支持；
- 其余 16 个 motif 可检验 HC-TF 未达到任一系统的半数阈值。

完整库中的 ELF3 由 CIS-BP `M03002_2.00` 驱动：患者 22/22、PDX 10/13
（可匹配样本 10/10）、细胞系 19/19 均为 `q <= 1e-5`；系统平均 log2 odds
ratio 分别为 0.341、0.521、0.462。PDX 的 CEBPB 为 7/13。严格作者规则会
在 ELF3 已有 JASPAR motif 时舍弃 CIS-BP motif；JASPAR ELF3 模型在这些
GC 匹配开放区域中近乎普遍出现，因而没有形成差异富集。这解释了 0 与 1 的
来源：不是峰处理是否能检测 motif，而是每个 TF 采用哪个 PWM 模型。

## 5. 与旧结果的关系

旧的未统一分析对 8 个预定义肺腺谱系阳性对照的系统支持数为：患者 8、PDX 0、
细胞系 0。这种完整的系统分裂本身就是背景、峰宽、峰数或数据库不可比的警报。
统一后，患者不再出现 8/8 的普遍富集，模型端能恢复 ELF3/CEBPB，说明旧的
“三端 0”确有技术构造因素。

但统一并没有恢复 TNBC 量级的结果：TNBC 为 31/56，LUAD 严格作者规则为
0/20，允许所有候选模型后也只有 1/20。两个癌种的 TF 集和数据构成不同，
这些比例不作正式跨癌种假设检验；它们只说明当前 LUAD 证据链明显更窄。

## 6. 判读

- **代码/统计层：通过。** TNBC 发布结果被精确恢复；LUAD 51 个可分析样本
  全部完成；独立结构与算术检查 151/151 通过。
- **峰到 HOMER 技术链：可以工作。** ELF3 在候选库和完整库中跨三端、跨
  全部可匹配样本稳定显著，不支持“PDX/细胞系 motif 流程整体失效”。
- **严格 TNBC 方法迁移：未恢复三端 TF。** 严格作者规则为 0，而且其阳性
  对照校准不完整，因此不能把 0 单独当作生物学否定证据。
- **广泛 LUAD 调控程序：当前不支持。** 最稳健的可继续验证对象只有 ELF3；
  不能据此宣称已复制 TNBC 的广泛 chromatin-informed TF 框架。

本报告不新增“至少几个 TF 才继续”的人为停止线，也不自动进入 Figure 3。
若后续围绕 ELF3 推进，应先用 motif footprinting、ELF3 ChIP/CUT&RUN 或其他
直接结合证据验证，而不是把 motif 富集等同于 TF 实际结合。

## 7. 主要审计文件

- `execution/luad-figure2/motif-equivalence/audit/final_motif_equivalence_receipt.json`
- `execution/luad-figure2/motif-equivalence/audit/motif_equivalence_independent_validation_receipt.json`
- `execution/luad-figure2/motif-equivalence/audit/tnbc_author_motif_replay_receipt.json`
- `execution/luad-figure2/motif-equivalence/audit/harmonized_input_receipt.json`
- `execution/luad-figure2/motif-equivalence/results/compiled/hc_tf_support_comparison.tsv`
- `execution/luad-figure2/motif-equivalence/results/compiled/lineage_positive_control_comparison.tsv`
- `execution/luad-figure2/motif-equivalence/results/anchor_equivalent/harmonized_hc_tf_triple_system_summary.tsv`
- `execution/luad-figure2/motif-equivalence/results/full1944/harmonized_hc_tf_triple_system_summary.tsv`

协议和事后等价性修订的时间顺序分别记录于：

- `analysis/luad-figure2-gate-protocol.md`；
- `analysis/luad-figure2-pre-result-amendments.md` 第 12 节。
