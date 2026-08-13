# LUAD Figure 1 数据审计

## 审计结论

本轮把 TNBC 锚点文献 Figure 1 的比较语法冻结为 **LUAD vs LUSC**，而不是事后选择 EGFR、TRU/PIF/PPR 或其他 LUAD 内部分型。四个分析层级均取得了可区分两种组织学、可按独立患者或模型计数的表达数据；大文件下载和计算均位于云服务器。

| 层级 | 数据源 | LUAD | LUSC | 独立单位 | Figure 1 用途 |
|---|---|---:|---:|---|---|
| 患者发现 | TCGA-LUAD/TCGA-LUSC | 516 | 501 | 每位患者一个原发肿瘤 | 合并建 ARACNe3 网络、差异表达、msVIPER 和逐样本 VIPER |
| 患者验证 | GSE81089 | 106 | 67 | 每位患者一个肿瘤 | 独立建网并验证 TF 方向 |
| PDX 验证 | NCI PDMR | 22 | 34 | 每位供体一个冻结规则选出的 PDX | 投射冻结的 TCGA regulon |
| 细胞系验证 | DepMap 22Q2 | 76 | 27 | 每个 DepMap ID 一个模型 | 投射冻结的 TCGA regulon |

## 数据处理决定

- TCGA：使用 STAR counts 做 limma-voom 表达效应，使用 STAR TPM 恢复后的原始 TPM 做 ARACNe3/VIPER；去除同一患者重复 aliquot。
- GSE81089：依据 GEO 官方编码 `1=squamous`、`2=adenocarcinoma`，排除正常、large cell/NOS 和其余样本；正式 ARACNe3 输入为 `log2(FPKM+1)`。原始 FPKM 触发 ARACNe3 `stof` 解析失败，失败尝试被保留归档，没有静默替换结果。
- PDMR：从 332 个符合条件的肺癌 PDX RNA 测序文件中，按供体去重；排除 originator、organoid、culture 和 repeat mouse，每位供体按预先规定的名称优先级只选一个 PDX。最终使用官方 RSEM TPM。公开清单没有可靠的数值传代字段，因此不宣称“精确最早传代”。
- DepMap：冻结 22Q2；EH7558 提供模型信息，EH7556 提供 `log2(TPM+1)`，计算时还原为 TPM。118 个组织学合格模型中 103 个有对应表达；15 个缺失表达的模型在 receipt 中明列而未插补。
- 调控因子：锚点 PAN-GO 文件名义数量与实际内容不完全一致；去重后为 2,059 个 gene symbol，其中 2,051 个存在于 TCGA 网络输入。未使用 CollecTRI、ULM 或其他替代网络。

## 表达尺度

| 队列 | 正式输入尺度 | 理由 |
|---|---|---|
| TCGA | raw TPM | 与 TNBC 作者 TCGA 代码先将 `log2(TPM+0.001)` 还原后建网的可执行行为一致 |
| GSE81089 | `log2(FPKM+1)` | 对应锚点独立患者验证队列采用的 log2 标准化表达；同时是该公开矩阵可被 ARACNe3 稳定解析的预先记录适配 |
| PDMR | RSEM TPM | 官方逐模型量化输出 |
| DepMap | TPM（由 EH7556 逆变换） | EH7556 的发布尺度为 `log2(TPM+1)` |

## 可追溯资源

- TNBC 锚点代码快照：Code Ocean capsule `7227095/v1`，本次使用的本地快照提交为 `edf5314ce5b9ee0e2f88b2310e7c2df5619ad888`。
- ARACNe3：固定提交 `3d8791a23e3bd8fd0d74f3b8d48f912e81d00f14`，随机种子 `2025`，默认一个 subnet、FDR 与 MaxEnt/DPI 剪枝。
- msVIPER null：仍使用 `viper::ttestNull` 的 1,000 次有放回置换和 seed 1。单核正式尝试在尚未产生任何 TF 结果的 2% 处显示预计运行超过 12 小时，故冻结使用该函数官方 `cores=32` 参数及 `L'Ecuyer-CMRG` 可重复并行流；停止点、理由和完成检查点哈希均单独保留。这是与结果无关的性能适配，未改变置换或检验定义。
- TCGA：UCSC Xena 的 GDC STAR counts、STAR TPM、临床信息和 GENCODE v36 注释。
- GSE81089：NCBI GEO `GSE81089` 的 counts、FPKM 和 SOFT 元数据。
- PDMR：NCI PDMR 公共 RNA-seq/RSEM 文件及其公开模型信息。
- DepMap：Bioconductor ExperimentHub 固定 22Q2 资源 EH7558（metadata）与 EH7556（TPM）。

云端正式目录：`/media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/luad-figure1/`。最终小型 manifests、哈希、软件版本、统计表和判定 receipt 将随 Figure 1 结果回传并纳入版本控制；原始大文件和中间网络保留在云端，不提交 Git。
