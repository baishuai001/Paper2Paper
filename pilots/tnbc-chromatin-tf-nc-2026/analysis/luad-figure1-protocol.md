# LUAD Figure 1 冻结分析协议

## 1. 目标与边界

本轮只回答一个问题：将 TNBC 锚点文献 Figure 1 的调控网络证据链迁移到肺癌后，能否发现并跨系统验证 **LUAD 相对 LUSC 的组织学特异 TF 活性程序**。本轮不做 ATAC、motif、预后、药敏或治疗推荐；这些分别属于后续 Figure 2–5。

比较在运行前冻结为 `LUAD vs LUSC`，不根据结果改成 EGFR、TRU/PIF/PPR 或其他 LUAD 内部分型。

## 2. 数据及独立性

| 层级 | 队列 | 用途 | 独立单位 |
|---|---|---|---|
| 患者发现 | TCGA-LUAD + TCGA-LUSC | 合并建 ARACNe3 网络、msVIPER 与差异表达 | 每位患者一个原发肿瘤 |
| 患者验证 | GSE81089 | 独立建网并验证 TF 方向 | 每位患者一个肿瘤；仅 histology 1/2 |
| PDX | NCI PDMR | 投射冻结的 TCGA regulon | 每位供体一个按冻结名称规则选择的祖先型 PDX；公共 RNA 清单无数值传代字段 |
| 细胞系 | DepMap 22Q2 | 投射冻结的 TCGA regulon | 每个 DepMap ID 一个模型 |

技术重复、同一供体的多个传代和同一患者多个 TCGA aliquot 不作为独立样本。GSE81089 的癌旁和其他组织学不进入比较。

## 3. 锚点方法约束

- 调控因子使用锚点代码所用 PAN-GO 列表；同时记录论文宣称数与实际唯一 symbol 数的差异。
- ARACNe3 固定 commit `3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14`，与作者一样只取默认 `subnet1`，保留默认 subsampling、FDR 和 MaxEnt/DPI pruning；只把随机种子固定为 `2025` 以便复现。
- signature 使用 `viper::rowTtest` 和双侧 p 值的正态分位数转换；null 使用 `ttestNull(per=1000, repos=TRUE, seed=1)`。由于本次 TCGA 比较有 1,017 位患者，作者隐式单核调用的实测预计耗时超过 12 小时；在尚未产生任何 msVIPER/TF 结果时，冻结为使用该函数官方 `cores=32` 参数分发相互独立的置换，并以 `L'Ecuyer-CMRG` 固定并行随机流。这是性能适配，不改变 1,000 次、有放回置换或后续检验定义。
- `msviper` 与逐样本 `viper` 均使用 `minsize=1`。
- 作者脚本在已正确读取表头后又执行 `[-1, ]`，实际额外删除第一条网络边；本实现将这个可执行代码行为原样保留，并在 receipt 中单列删除前后边数。
- TCGA ARACNe3 输入为原始 TPM；计数只用于 limma-voom 表达效应。独立患者验证 GSE81089 使用 `log2(FPKM+1)`，但它仍是 RNA-seq 同平台验证；与 METABRIC 同为 log2 数值尺度不构成 microarray 跨平台等价。跨平台验证须另行加入 GSE41271 等经审计的微阵列队列。PDMR 使用官方 RSEM TPM，DepMap 将固定 22Q2 的 `log2(TPM+1)` 还原成 TPM。
- Figure 1 的主筛选保持作者代码的逻辑结构；独立患者验证要求 msVIPER raw p<=0.05、表达 BH FDR<=0.05 且方向一致。
- 不允许用 CollecTRI、ULM 或其他算法代替未通过的 ARACNe3/VIPER 结果。

## 4. 预先冻结的停止规则

“图已生成”与“生物学闸门通过”分开报告。只有同时满足以下条件，Figure 1 生物学闸门才为 PASS：

1. TCGA 中存在作者规则定义的 LUAD-specific TF；
2. 至少 3 个 LUAD-specific TF 在独立患者队列 GSE81089 中以预定阈值同方向复现；
3. 用这些复现 TF 构建的 LUAD 程序分数，在 PDMR PDX 和 DepMap 细胞系中均为 LUAD 高于 LUSC；
4. 上述两个模型队列均达到 BH FDR<=0.05 且 Hedges' g>=0.5。

任一条件失败：仍交付全部真实结果、失败位置和图，但判为 FAIL，并在 Figure 1 停止；不得事后改变亚型、阈值、训练队列或 TF 集合来“救回”结果。

## 5. 预定输出

- Figure 1A：数据与算法流程；
- Figure 1B：TCGA TF 的 msVIPER NES 与表达 logFC；
- Figure 1C：TCGA 单样本 TF 活性；
- Figure 1D：PDMR PDX TF 活性；
- Figure 1E：DepMap 细胞系 TF 活性；
- C–E 的列均按 TF 活性做无监督层次聚类，不按 LUAD/LUSC 标签预切分；另以 Fisher 精确检验量化两簇与组织学的对应关系；
- 所有 panel 的数值表、样本 manifest、TF 清单、跨系统统计、软件/提交版本、SHA256 和最终 PASS/FAIL receipt。

## 6. 预先冻结后的软件兼容处理

GSE81089 首次正式运行在 `viper 1.38.0::aREA` 内暴露单靶基因 regulon 的矩阵降维错误：对 `1 x n` 权重矩阵逐列 `apply` 后得到向量，随后 `colSums` 因缺少二维维度而停止。该错误发生在全部网络、逐样本 VIPER 活性和 1,000 次 null 置换完成之后、任何 GSE81089 msVIPER/TF 结果产生之前。兼容补丁只把该向量恢复为原本的 `1 x n` 形状；不删除两个单靶基因 regulon，不改变 `minsize=1`、数值、阈值、种子、置换次数或停止规则。首次错误日志、补丁代码、检查点复用和重跑脚本均纳入审计。
