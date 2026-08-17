# LUAD Figure 2 审计修正与版本记录

> **2026-08-15 状态：本文件记录的 v1 规则已被
> `luad-figure2-gate-protocol.md` v2 取代。** v1 曾错误地以 97 个外部复现 TF
> 作为 Figure 2 输入，并加入 20% 样本失败率、最低 HC/motif 数和匹配置换
> 否决。当前原文式运行改用 Figure 1 的 158 个 LUAD 发现 TF，恢复全部
> 13 个 PDX；旧 `FAIL_DATA` 及本文件保留作决策历史，不再控制执行。

第 1–11 节收纳在当时 LUAD Figure 2 启动子或 motif 结果生成之前发现并修正的
问题；第 12 节明确记录候选库首轮结果之后、因逐行核对作者代码而发生的
等价性修订。任何一节均不得被解释为看结果后放宽阈值。

## 1. Figure 1 TF 父表不是 97 行

`externally_replicated_TFs.tsv` 的 SHA256 为 `eace4fdf6509d01e8a9648b34b4dc3fdfb5955627c2c7a526c6ce1495a073844`，但该父表实际同时包含 97 个 LUAD 和 153 个 LUSC 复现 TF，共 250 行。Figure 2 输入现明确固定为：

`externally_replicated == TRUE AND discovery_category == LUAD`

排序后三列规范 TSV 共 97 行，SHA256 为 `64b4aa4e226c6e9142ecd560ec51ec1f49e0c9823d07e332757deeea3aa3c69b`。执行与独立验证均同时检查父表 SHA、派生表 SHA 和 97 行计数，防止把 LUSC TF 混入 LUAD 染色质闸门。

## 2. JASPAR 与 CIS-BP 必须独立运行

锚点 Methods 写明两套数据库独立使用；Supplementary Data 5 也分别提供 `HOMER_JASPAR_TCGA_PDX_CellLine` 和 `HOMER_CISBP_TCGA_PDX_CellLine` 工作表，同一 TF 可同时出现在两表。原先的“JASPAR 优先、CIS-BP 仅补缺”实现已在运行前撤销。当前实现为：

- 每个样本分别运行 JASPAR 2024 与 CIS-BP 2.0 HOMER；
- 每库分别做 BH 校正并使用 `q < 1e-5`；
- TF 层面取任一数据库中显著且 LOR>0 的对应 motif，同时保留两库全量结果；
- 外部 PWM 使用 HOMER 随附 `parseJasparMatrix.pl` 的阈值 0 导入约定。

## 3. 允许范围内的 QC 损耗不能被代码重新禁止

协议允许每个原始数据系统不超过 20% 的硬 QC 失败，并设置患者/PDX/细胞系最小 n 为 7/10/10。原子 manifest、绘图标签、置换与完整性验证已改为读取最终合格 n，而不是重新硬编码 22/13/19。名义起始数仍完整记录。

## 4. MACS2 运行环境必须隔离

服务器预装 MACS2 与全局 NumPy 2.x 存在 ABI 不兼容，A427 冒烟测试在峰调用前明确失败。未修改数据或参数；改为隔离并固定 `MACS2 2.2.9.1 + NumPy 1.26.4`。失败日志保留，A427 从服务器已保存的 FASTQ 重启，不重新下载数据。

## 5. 单端细胞系不能把 ataqv 的 `TSS=NA` 误判为无染色质信号

ataqv 1.3.0 源码将 HQAA 明确定义为正确配对、非重复、MAPQ≥30 的常染色体 reads，并用 HQAA 片段计算 TSS enrichment。因此，DRA 单端细胞系必然得到 `HQAA=0` 和 `TSS=NA`。A427 冒烟运行在任何 TF、启动子或 motif 闸门结果生成前确认了这一工具边界。

当前实现保留 ataqv 原始 JSON；同时对所有 PDX 和细胞系使用同一套 Tn5 插入位点（正链 +4 bp、负链 −5 bp）计算链方向校正的 ±1 kb TSS profile，并按两端各 100 bp 的平均值归一化。精确中心值、中心 101 bp 最大值和均值全部报告，但均不参与硬 QC 或锚点信号判定。这是布局兼容的补充 QC，不是事后新增淘汰标准。

## 6. TSV 行尾必须同时兼容 Python 与 Bash

原始运行清单由 Python `csv.DictWriter` 生成时继承了默认 CRLF 行尾，Bash 将最后一列读成 `TRUE\r`，导致第一次按系统批处理没有启动任何新样本却生成了空的完成标记。该标记在正式批处理前删除，没有产生错误样本结果。当前所有 Python TSV writer 固定为 LF，所有 Bash 清单读取器同时显式去除残留 CR；批处理完成后还必须由 32 行聚合 QC 和样本级产物进行独立核对，不能依赖驱动完成标记本身。

## 7. 自动报告必须与判定 JSON 一致且使用 UTF-8

结果生成前静态检查发现，旧版 `finalize_figure2_gate.py` 的判定 JSON 逻辑可读，但 Markdown 模板中的中文曾被错误编码。当前报告生成器已整体重写为 UTF-8，并把“有资格另行规划 Figure 3”与“本次是否继续运行”拆开：仅 `PASS` 时 `eligible_to_plan_Figure3=TRUE`，但所有结局均固定 `continue_to_Figure3=FALSE`、`stop_after_this_gate=TRUE`。独立验证继续检查显式停止字段。

## 8. 原子交付边界不变

Figure 2、Supplementary Figure 3 与 Supplementary Figure 4A–C 必须共享同一最终 manifest、97-TF SHA、QC 表和阈值配置。任何一组缺失都不得称 Figure 2 完成；最终判定后停止，不进入 Figure 3。

## 9. HOMER 分母与背景计数必须采用程序实际输出

正式 motif 结果生成前，使用 A427/A549 peaks 和 HOMER 自带 E2F1 motif 分别模拟 JASPAR 与 CIS-BP 命名，验证了 `knownResults.txt` 的真实结构。HOMER 会在序列提取时丢弃少量无效序列，并将 GC 归一化后的背景 motif 命中数报告为加权小数；因此不能把输入位置文件行数当作最终分母，也不能把背景命中数四舍五入为整数。

解析器现固定为：从 `# of Target/Background Sequences with Motif(of N)` 表头读取 HOMER 实际分母，保留目标和背景加权小数，再计算带 0.5 连续性校正的 odds ratio。烟雾测试的断言为 target `837.0/1000.0`、background `1662.9/1993.0`，两套数据库均得到 `log2(OR)=0.02582205`、`q=0.4306`。该修正发生在任何正式 HC-TF motif 结果生成之前，不改变 `q<1e-5`、LOR>0 或半数样本阈值。

## 10. 不得把早停后未运行的样本伪装成 QC 失败

细胞系端出现第 4 个硬 QC 失败后，4/19 已不可逆超过 20% 上限，因此全局数据闸门至少为 `FAIL_DATA`。从这一时点起不再启动新的 PDX 计算任务；已下载的 13 个 PDX 原始数据、manifest 和已启动的 PDX 流程验证均保留。

聚合器新增显式早停模式：锁定失败的系统必须全部完成且独立满足失败条件；另一系统未启动的样本记录为 `NOT_RUN_AFTER_LOCKED_DATA_GATE`，已经启动但按停止规则终止的样本记录为 `ABORTED_AFTER_LOCKED_DATA_GATE`，两者均不计入实际 QC 失败率。自动报告必须分别展示“已完成/通过/实际失败”“早停后中止”和“早停后未启动”，不得用缺失输出夸大 PDX 数据质量问题。该改动只影响失败审计的诚实表达，不允许绕过数据闸门进入 promoter、motif 或 Figure 3。

## 11. 候选调控因子启动子定位不得擅自增加 protein-coding 过滤

第一次原文式启动子计算在产生 HC-TF 结果前，因 `MEIS3P1`、`MEIS3P2`、
`NANOGP1` 无法映射到“protein-coding transcript”而停止。核查确认三者并非
拼写错误：它们均来自冻结的 PAN-GO 调控因子输入，GENCODE v47 也存在
exact-symbol transcript，只是基因类型分别为 processed 或 transcribed
pseudogene。TNBC 锚点的启动子规则是候选基因座/TSS 开放规则，并未声明
“仅蛋白编码候选”这一额外过滤。

修订后仍保留全部 158 个 Figure 1 候选：候选启动子按 exact-symbol GENCODE
transcript 定位，三项非蛋白编码候选及类型写入独立 receipt；全基因组峰注释
背景仍限制为 protein-coding transcript。该修订不改变样本、启动子窗口、半数
样本、NES 或 motif 阈值，也没有查看任何 HC-TF 或 motif 正式结果。

## 12. motif 技术等价性审计必须采用作者代码的实际数据库选择规则

2026-08-17 的技术等价性审计先完成了 47 个候选 motif 合并后一次校正的
LUAD 结果（回执时间 `2026-08-17T04:59:08Z`）。随后逐行复核 TNBC 作者
Figure2/06B、07A 与 08 代码以及 `JASPAR_CISBP_comparison_long.csv`，确认
作者的实际实现是：JASPAR 与 CIS-BP 分开运行和校正；每个 TF 有 JASPAR
motif 时使用 JASPAR，只有 `CISBP_only` 时才使用 CIS-BP。TNBC 的 56 个
motif 可检验 TF 中，44 个同时有两库 motif、5 个仅 JASPAR、7 个仅 CIS-BP。

因此，先前合并 47 个模型的一次校正不能再标为“作者等价主分析”，而改列为
更严格的联合校正敏感性结果。新增主分析严格执行：

- JASPAR 21 个候选模型与 CIS-BP 26 个候选模型分别运行 HOMER；
- 两库分别保留其自身的 Benjamini q 值；
- TF 层面优先 JASPAR、无 JASPAR 时才使用 CIS-BP；
- `q <= 1e-5`、半数样本、22/13/19 固定分母、8 个预定义肺腺谱系阳性对照
  均不改变；
- 47-model 联合校正与 1,944-model 完整库联合校正均保留，不得择优替代作者
  等价结果。

这次修订发生在发现联合校正结果之后，故在此明确记录时间顺序。修订方向是
回到作者代码并可能使结果更严格，不是根据结果放宽阈值或挑选数据库。
