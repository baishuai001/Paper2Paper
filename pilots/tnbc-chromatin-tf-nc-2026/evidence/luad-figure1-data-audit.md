# LUAD Figure 1 数据审计

> **2026-08-16 复核更正（优先于下文旧结论）**：GSE81089 官方肿瘤队列包含 108 例腺癌和 67 例鳞癌；现有流程因表达矩阵列名 `L608T_2122`、`L771T_1` 未映射回样本名 `L608T`、`L771T`，静默漏掉 2 例腺癌，故下文 106/67 及由此产生的外部复现统计均为待重跑结果。GSE81089 是 RNA-seq，不是 METABRIC microarray 的跨平台等价替代。既往冻结前未形成可追溯的 LUAD/LUSC 微阵列候选审计；该流程缺口已确认。TCGA 发现队列、158 个 TCGA-LUAD 特异 TF 及 Figure 2 的输入不受此样本映射错误影响。

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

## TNBC 锚点 Figure 1 数据链复核（2026-08-16）

TNBC Figure 1 不是“四套 RNA-seq 并列验证”，而是以下四层证据链：

| 层级 | 真实数据模态与样本数 | 实际用途 |
|---|---|---|
| TCGA-BRCA 患者 | 公共 bulk RNA-seq；初始 1,099 个肿瘤记录，实际活动矩阵 1,001 个原发肿瘤（142 TNBC、859 non-TNBC） | 独立构建 TCGA ARACNe3 网络并发现 150 个 TNBC-specific TF |
| METABRIC 患者 | 公共 Illumina HT-12 v3 microarray；1,980 个有表达数据的肿瘤（233 TNBC、1,747 non-TNBC） | 在不同平台上另建 ARACNe3 网络；复现 76/150 个 TCGA TNBC TF，不与 TCGA 原始表达值拼接 |
| UHN PDX | 作者院内生成的 bulk RNA-seq；73 个表达谱 | 使用 TCGA regulon 做逐样本 VIPER 投射，不以 73 个 PDX 另建网络 |
| 乳腺细胞系 | 公共 GSE73526 bulk RNA-seq；实际 82 条细胞系 | 使用 TCGA regulon 投射；论文正文/Figure 1A 的 83 为计数笔误 |

需保留的锚点勘误：

- Methods 明确称 METABRIC 为 microarray；Data Availability 将其写成 processed RNA-seq，二者矛盾，数据平台与 cBioPortal 元数据支持 Methods。
- TCGA 的 1,099 是未完成最终纳排的起始表达记录；Figure 1C 实际活动矩阵为 1,001。差额由 7 个转移样本和 91 个受体状态缺失/不明确的原发样本构成。
- GSE73526 GEO 官方设计、补充表和活动矩阵均为 82；论文的 83 不能作为真实样本数。
- Figure 2 从 TCGA 发现的全部 150 个 TF 开始，而不是只从 METABRIC 复现的 76 个开始。

## 当前 LUAD Figure 1 与锚点的真实对应关系

| 层级 | 当前 LUAD 数据 | 与锚点关系 |
|---|---|---|
| 患者发现 | TCGA-LUAD/LUSC RNA-seq，516 + 501 = 1,017 | 对应 TCGA-BRCA RNA-seq 发现层；158 个 LUAD-specific TF 仍是 Figure 2 合法入口 |
| 独立患者验证 | GSE81089 RNA-seq，正确目标应为 108 LUAD + 67 LUSC = 175 | 只是独立同平台验证，不等价于 METABRIC microarray 跨平台验证 |
| PDX | PDMR RNA-seq，22 LUAD + 34 LUSC = 56 | 对应 PDX 投射层 |
| 细胞系 | DepMap 22Q2 RNA-seq，76 LUAD + 27 LUSC = 103 | 对应细胞系投射层 |

当前遗漏不是有意排除微阵列，而是冻结前审计不足：候选队列表没有形成实质记录，锚点 Data Availability/Methods 的模态矛盾未先解决，也没有沿 GSE81089 原论文 Methods 追踪其 7 个 Affymetrix 外部队列。GSE81089 本身可独立建网，因此选择它并非无效；错误在于把“结构相似”写成了近似“证据模态等价”。

## 已核验的 LUAD/LUSC 微阵列修复候选

| 优先级 | 队列 | 平台与可用结构 | 决定 |
|---:|---|---|---|
| 1 | GSE41271 | GPL6884 Illumina WG6-v3；275 个肿瘤，183 LUAD、80 LUSC、12 其他 | 最适合单队列独立建 ARACNe3 网络和跨平台复现 |
| 2 | GSE37745 | GPL570；196 个 NSCLC，106 LUAD、66 LUSC、24 large-cell | 可作第二独立网络/投射验证 |
| 3 | GSE50081 | GPL570；181 个早期 NSCLC，约 127 LUAD、42 严格 LUSC | 可作外部验证 |
| 4 | GSE30219 | GPL570；直接 LUAD/LUSC 比较约 85/61，另含多种肺肿瘤 | 异质性较高，优先投射而非单独主网络 |
| 5 | GSE19188 | GPL570；约 45 LUAD、27 LUSC，另有正常/大细胞癌 | 样本偏小，仅建议投射 |

GSE42127 的 176 个样本全部包含于 GSE41271，不能作为独立队列或与其相加。不同研究的已处理芯片矩阵不得直接拼接，也不得与 RNA-seq FPKM/TPM 拼接；若利用多个 GPL570 队列，应分别建网/投射后做效应合并，或从 raw CEL 联合标准化并先审计批次与组织学混杂。

## 修复顺序与影响边界

1. 显式映射 `L608T_2122 -> L608T`、`L771T_1 -> L771T`，将 GSE81089 恢复为 108/67 后重跑该分支和 Supplementary Figure 1–2。
2. 以 GSE41271 建独立 microarray ARACNe3/VIPER 验证层；其余芯片队列作为独立外部复现/效应合并，不盲目拼矩阵。
3. 更新 Figure 1A、协议、样本流图、复现 TF 数和最终报告；旧的 97/158 复现结果在重跑前保持“暂定”。
4. Figure 2 不等待上述修复：它按锚点同样的逻辑使用全部 158 个 TCGA 发现 TF，故 TCGA 入口和当前 ATAC 分析不受 GSE81089 两样本错误或微阵列缺口影响。
