# LUAD Paper2Paper 科学问题、研究计划、创新性与当前 Figure Legends

> 版本：2026-08-17
>
> 锚点文献：*A chromatin-informed transcriptional regulatory framework to stratify patients and guide therapy selection in triple-negative breast cancer*
>
> 当前迁移问题：`LUAD versus LUSC`，即肺腺癌相对肺鳞癌的组织学/谱系特异调控程序。
>
> 重要边界：当前设计不是 LUAD 内部患者亚型研究，也尚不能声称已经实现 LUAD 患者治疗分层。

## 1. 拟解决的科学问题

### 1.1 一句话科学问题

能否从大规模肺癌转录组中发现一个在患者肿瘤、PDX 和细胞系中可重复、并由三系统开放染色质和 DNA motif 共同支持的 LUAD 特异转录因子（TF）程序；进一步解析其调控网络、临床意义和可重复的药物脆弱性？

### 1.2 分层科学问题

1. **发现层**：LUAD 相对 LUSC 是否存在稳定的 TF 活性程序，而不仅是 TF 表达差异？
2. **跨系统验证层**：该程序能否在独立患者队列、PDX 和细胞系中保持方向一致？
3. **染色质层**：候选 TF 的启动子是否在 LUAD 患者、PDX 和细胞系中普遍开放，其 DNA-binding motif 是否在三系统可及染色质中富集？
4. **网络层**：染色质支持的高置信 TF 是否形成可解释的 regulon 社区、协同关系和共享靶基因结构？
5. **临床层**：在仅限 LUAD 患者的分析中，这些 TF 活性或网络模块是否与复发和生存独立相关？
6. **治疗层**：这些 TF 状态是否在多个独立药敏平台中指向相同药物/通路，并能在具有治疗前表达和体内反应数据的 PDX 中获得验证？

### 1.3 中心假设

LUAD 的谱系身份不是由单个标志基因决定，而是由一组具有一致 regulon 活性、启动子可及性和 motif 占位潜力的 TF 网络维持；该网络内部的连续状态可能解释部分 LUAD 患者结局和治疗反应差异。

## 2. 主要研究计划

| 阶段 | 主要问题与分析 | 预定主图 | 当前状态 |
|---|---|---|---|
| 1. TF 程序发现与跨系统投射 | TCGA-LUAD/LUSC 分别按每位患者一个原发肿瘤纳入；ARACNe3 建网，msVIPER 比较 TF 活性，limma-voom 计算表达效应；在独立患者队列另行建网；将冻结 TCGA regulon 投射到 PDMR PDX 和 DepMap 细胞系 | Figure 1A-E | 已生成。TCGA、PDX、DepMap 主体结果可用；GSE81089 漏纳 2 例 LUAD，所有依赖 106/67 和 97 个复现 TF 的面板待重跑 |
| 1b. 跨平台加固 | 修复 GSE81089 别名映射为 108 LUAD/67 LUSC；以 GSE41271 microarray 独立建 ARACNe3/VIPER 网络，其他芯片队列逐队列投射或做效应合并，不拼接不同平台表达值 | Figure 1/Supplementary Figure 1-2 修订版 | 尚未完成；这是补齐 TNBC 中 METABRIC 所承担的跨平台验证 |
| 2. 三系统染色质支持 | 从全部 158 个 TCGA-LUAD 发现 TF 出发；22 个患者、13 个 PDX、19 条细胞系分别按至少半数样本启动子开放；仅排除三系统平均 NES 同时小于 0 的 TF；对 HC-TF 在完整 JASPAR 2024 和 CIS-BP 2.0 motif 库中运行 HOMER，单样本 adjusted p<1e-5，系统支持要求至少半数样本 | Figure 2A-E；Supplementary Figure 3；Supplementary Figure 4A-C | 正在运行。启动子/NES 层已得到 31 个 HC-TF，其中 20 个有可检验 motif；最终三系统 motif TF 数尚未冻结 |
| 3. 调控网络结构 | 使用 TCGA ARACNe3 signed regulon 描述 HC-TF 的 regulon 大小、调控方向、TF-TF 活性协同、共享靶基因和网络社区；在单个 LUAD 患者中量化模块活性分布 | Figure 3 | 仅在 Figure 2 完整结果具有可解释输入后启动；不设置事后“救结果”阈值 |
| 4. 临床关联 | 仅在 LUAD 患者内部检验 HC-TF/模块与 OS、RFS 的关系；TCGA 发现，独立 LUAD 队列验证；Cox 模型预先纳入年龄、性别、分期和吸烟等可获得协变量 | Figure 4 | 数据资源可获得，但尚未运行。必须避免把 LUAD/LUSC 组织学差异误当作 LUAD 内部预后信号 |
| 5. 药物脆弱性与体内验证 | 在 LUAD 细胞系内关联 TF 活性与 PRISM、GDSC2、CTRPv2 药物反应；要求跨平台方向重复，再检查治疗前 RNA 与同药 PDX 反应的样本交集 | Figure 5 | 预计只能先完成部分面板。共同细胞系/药物样本量及同药 PDX 交集仍是主要风险 |

## 3. 预期创新性及其成立条件

### 3.1 可能形成的创新

1. **从表达标志物推进到调控活性**：以 ARACNe3 regulon 和 VIPER 活性识别 TF 程序，避免把 TF mRNA 丰度直接等同于 TF 功能。
2. **患者-模型-细胞系三端染色质交叉验证**：同一批 TF 同时接受启动子开放、平均活性和 motif 富集检查，比单队列表达分析证据更完整。
3. **将谱系程序连接到网络、临床和药物反应**：如果新 HC-TF/网络模块能在 LUAD 内部解释结局并跨药敏平台重复，研究可从组织学分类推进到可转化机制。
4. **可审计的真实数据复现体系**：样本独立性、原始 ATAC 处理、完整 motif 数据库、多重检验分母、阈值和失败记录均留有 receipt，可作为 Paper2Paper 的可复用范式。

### 3.2 目前不能提前宣称的创新

- `NKX2-1/FOXA/HNF/ELF3` 与 `TP63/SOX2` 等 LUAD-LUSC 谱系轴已有充分先验；只重新发现这些 TF，方法迁移成立，但生物学新颖性有限。
- Figure 1 当前回答的是 **LUAD 与 LUSC 的区别**，不是 LUAD 内部患者分层。只有 Figure 3-5 在 LUAD 内部找到可重复状态、结局或治疗差异，才可使用“patient stratification”或“therapy selection”表述。
- motif 富集代表潜在结合位点支持，不等于 ChIP/CUT&RUN 的直接 TF 占位证据。
- 公共 PDX 和药筛数据只有在治疗前表达、明确组织学及同药反应形成足量交集时，才能作为治疗验证；不能用不匹配模型拼接成验证链。

## 4. 当前 Figure 完成度与推荐使用版本

当前本地图件总目录：

`C:\Users\bai\OneDrive\Documents\Codex\2026_07_09 Literature_Reproduction\Paper2Paper\pilots\tnbc-chromatin-tf-nc-2026\execution\luad-figure1\results\figures`

云服务器对应目录：

`/media/desk16/iy13202/projects/Paper2Paper/pilots/tnbc-chromatin-tf-nc-2026/execution/luad-figure1/results/figures`

下表同时提供 PDF（检查文字、矢量排版）与 PNG（快速预览）链接。

### 4.1 主图

| Panel | 当前状态 | 推荐文件 |
|---|---|---|
| Figure 1A | 已生成但待修订；仍显示 GSE81089 106/67 和 97/158 | [PDF](../execution/luad-figure1/results/figures/Figure1A_workflow_stacked.pdf) / [PNG](../execution/luad-figure1/results/figures/Figure1A_workflow_stacked.png) |
| Figure 1B | 当前有效；完全来自 TCGA 发现层 | [PDF](../execution/luad-figure1/results/figures/Figure1B_TCGA_logFC_msviper.pdf) / [PNG](../execution/luad-figure1/results/figures/Figure1B_TCGA_logFC_msviper.png) |
| Figure 1C | 当前有效；推荐使用带临床/分子注释版本 | [PDF](../execution/luad-figure1/results/figures/Figure1C_TCGA_TF_activity_annotated.pdf) / [PNG](../execution/luad-figure1/results/figures/Figure1C_TCGA_TF_activity_annotated.png) |
| Figure 1D | 当前有效；PDX 使用全部可计算 TCGA-specific TF，不依赖 GSE81089 的 97-TF 清单 | [PDF](../execution/luad-figure1/results/figures/Figure1D_PDMR_PDX_TF_activity_annotated.pdf) / [PNG](../execution/luad-figure1/results/figures/Figure1D_PDMR_PDX_TF_activity_annotated.png) |
| Figure 1E | 当前有效；细胞系使用全部可计算 TCGA-specific TF，不依赖 GSE81089 的 97-TF 清单 | [PDF](../execution/luad-figure1/results/figures/Figure1E_DepMap_22Q2_TF_activity_annotated.pdf) / [PNG](../execution/luad-figure1/results/figures/Figure1E_DepMap_22Q2_TF_activity_annotated.png) |
| Figure 2 | 尚未完成原子输出，不应引用半成品图 | 完成后位于 `execution/luad-figure2/results/figures/` |

Figure 1A 的横向备选版位于 [PDF](../execution/luad-figure1/results/figures/Figure1A_workflow.pdf) / [PNG](../execution/luad-figure1/results/figures/Figure1A_workflow.png)。为降低与锚点版式近似的风险，当前推荐上下布局的 `Figure1A_workflow_stacked`。

### 4.2 补充图

| Panel | 当前状态 | 推荐文件 |
|---|---|---|
| Supplementary Figure 1A | 当前有效；TCGA 纳排流程 | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure1A_TCGA_inclusion.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure1A_TCGA_inclusion.png) |
| Supplementary Figure 1B | 待修订；图中 2 个“unmatched”实际是未处理的表达列别名 | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure1B_GSE81089_inclusion.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure1B_GSE81089_inclusion.png) |
| Supplementary Figure 2A | 待 GSE81089 108/67 重跑 | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure2A_GSE81089_logFC_msviper.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure2A_GSE81089_logFC_msviper.png) |
| Supplementary Figure 2B | 待重算 overlap | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure2B_TF_overlap.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure2B_TF_overlap.png) |
| Supplementary Figure 2C | 当前有效；TCGA 样本相关性 | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure2C_TCGA_sample_correlations.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure2C_TCGA_sample_correlations.png) |
| Supplementary Figure 2D | 当前有效；PDMR PDX 样本相关性 | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure2D_PDMR_PDX_sample_correlations.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure2D_PDMR_PDX_sample_correlations.png) |
| Supplementary Figure 2E | 当前有效；DepMap 样本相关性 | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure2E_DepMap_22Q2_sample_correlations.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure2E_DepMap_22Q2_sample_correlations.png) |
| Supplementary Figure 2F | 待重跑；依赖旧的 97 个 GSE81089 复现 TF | [PDF](../execution/luad-figure1/results/figures/SupplementaryFigure2F_replicated_LUAD_TF_effects.pdf) / [PNG](../execution/luad-figure1/results/figures/SupplementaryFigure2F_replicated_LUAD_TF_effects.png) |
| Supplementary Figure 3、4A-C | 尚未完成原子输出 | 完成后位于 `execution/luad-figure2/results/figures/` |

`SupplementaryFigure1_inclusion_flow` 是较早的三队列总览版本，既未完整展示四层数据链，又包含旧 GSE81089 计数，不建议作为正式补图；应以 Supplementary Figure 1A-B 及后续新增的模型纳排面板为准。

## 5. 当前 Figure Legends（中文工作稿）

### Figure 1 | ARACNe3/VIPER 识别并跨实验系统验证 LUAD 特异 TF 活性程序

**(A)** 研究流程。TCGA-LUAD 与 TCGA-LUSC 原发肿瘤 RNA-seq 在合并肺癌表达空间中构建 ARACNe3 interactome，并以 VIPER/msVIPER 计算单样本及组间 TF 活性，得到 158 个 LUAD-specific 和 193 个 LUSC-specific TF。TCGA regulon 随后投射到患者肿瘤、PDMR PDX 和 DepMap 细胞系；GSE81089 使用独立 ARACNe3 网络进行患者队列复现。当前图中的 GSE81089 `106 LUAD/67 LUSC` 和 `97/158` 为待修订数值，正确目标样本数为 108/67。

**(B)** TCGA 中 TF 活性与表达效应的一致性。每个点代表一个进入 ARACNe3/VIPER 分析的调控因子；横轴为 LUAD 相对 LUSC 的 msVIPER normalized enrichment score（NES），纵轴为 limma-voom 估计的 log2 fold change。颜色表示按锚点作者规则得到的 LUAD-specific、LUSC-specific、仅 VIPER 显著、仅表达效应达到条件或均不显著的类别；标注显示绝对 NES 最大的代表性 TF。黑线及阴影为线性拟合和 95% 置信区间。Pearson `r=0.7255`，`p<2.2×10^-16`。

**(C)** 1,017 个 TCGA 原发肺癌（516 LUAD、501 LUSC）的单样本 TF 活性热图。行是 351 个 TCGA-specific TF（158 LUAD、193 LUSC），列是患者；每个 TF 在患者间进行行标准化，红色和蓝色分别表示相对较高和较低的 TF 活性，而不是基因表达量或跨 TF 的绝对活性。TF 行按发现类别固定排列，患者列仅按 TF 活性进行无监督层次聚类。顶部注释依次显示组织学、性别、年龄、生存状态、分期、吸烟包年、已发表表达亚型以及主要癌基因和抑癌通路突变；灰色表示不可获得。两簇与组织学的最佳对应准确率为 93.8%，Fisher 检验 BH FDR=`2.48×10^-207`。

**(D)** 56 个 NCI PDMR 肺癌 PDX（22 LUAD、34 LUSC）中冻结 TCGA regulon 的 VIPER 活性。显示 351 个可计算的 TCGA-specific TF；行和列均基于行标准化后的 TF 活性进行无监督聚类。顶部注释显示组织学、性别、年龄、已知转移病史、吸烟史、活检部位、组织类型和公开分子记录；缺失值不插补。两簇与组织学的最佳对应准确率为 91.1%，BH FDR=`2.76×10^-10`。

**(E)** 103 条 DepMap 22Q2 肺癌细胞系（76 LUAD、27 LUSC）中冻结 TCGA regulon 的 VIPER 活性。一个 TCGA-specific TF 在该表达矩阵中不可计算，故显示 350 个 TF。顶部注释显示组织学、性别、年龄、原发或转移来源、取材部位、生长方式以及主要癌基因和抑癌通路突变。两簇与组织学的最佳对应准确率为 80.6%，BH FDR=`7.03×10^-5`。

**热图读图说明。** Figure 1C 中间大片红色并非绘图错误。颜色是在每个 TF 行内进行标准化后的相对活性；当一组 LUAD-specific TF 在同一批 LUAD 样本中整体升高时，会形成连续红色模块，与对应的 LUSC-specific 蓝色模块相对。它证明谱系程序强，但不代表所有红色 TF 的原始 NES 相同，也不代表患者内每个 TF 都被绝对激活。

### Supplementary Figure 1 | Figure 1 患者队列纳入与排除流程

**(A)** TCGA 肺癌发现队列。下载的 counts/TPM 表达列共 1,141 列；限定 TCGA sample type 01 后保留 1,029 个原发肿瘤列，再按患者去除 12 个重复原发 aliquot，最终得到 1,017 位独立患者，包括 516 例 LUAD 和 501 例 LUSC，用于 ARACNe3、VIPER/msVIPER 和差异表达分析。

**(B)** GSE81089 独立患者队列。218 个 RNA-seq 样本中排除 19 个配对正常组织和 24 个 large-cell/NOS 肿瘤后，元数据定义的 LUAD/LUSC 肿瘤共 175 个。当前旧图因未处理 `L608T_2122→L608T` 和 `L771T_1→L771T` 两个表达列别名，将两例 LUAD 错记为 unmatched，显示 173 例（106 LUAD、67 LUSC）；修订版应为 175 例（108 LUAD、67 LUSC）。

### Supplementary Figure 2 | Figure 1 的独立复现、样本结构与跨系统 TF 效应

**(A)** GSE81089 独立网络中 msVIPER NES 与 limma-voom log2 fold change 的关系。每个点代表一个调控因子，颜色定义同 Figure 1B；黑线及阴影表示线性拟合和 95% 置信区间。当前 173 样本版本的 Pearson `r=0.72`、`p=1.3×10^-322`；该数值须在恢复两例 LUAD 后重新计算。

**(B)** TCGA 与 GSE81089 中方向一致的 LUAD-specific 和 LUSC-specific TF 集合。柱形分别表示仅在 GSE81089、两队列重叠及仅在 TCGA 出现的 TF 数。当前旧版本 LUAD 为 103/97/61，LUSC 为 132/153/40；因依赖 GSE81089 旧样本集，修订前不作为最终计数。

**(C-E)** 分别显示 TCGA 患者、PDMR PDX 和 DepMap 细胞系的样本-样本 Pearson 相关矩阵。相关系数由 Figure 1C-E 的 TF 活性向量计算，行列使用同一个基于 `1-r` 的 average-linkage 层次聚类树；组织学显示于顶部和左侧。TCGA 和 PDMR 使用 351 个 TF，DepMap 使用 350 个可计算 TF。

**(F)** 旧版本中 97 个 GSE81089 同方向复现 LUAD TF 在 TCGA、GSE81089、PDMR PDX 和 DepMap 22Q2 中的 LUAD-LUSC 平均活性差。每个队列内对 TF 效应进行列标准化后绘图，TF 行进行层次聚类，队列顺序固定；红色表示该 TF 在相应队列中相对较强的 LUAD 偏向，蓝色表示相对较弱。该 panel 依赖待修订的 97-TF 清单，不能作为最终图使用。

## 6. 后续 Figure Legends 的冻结方式

Figure 2、Supplementary Figure 3 和 Supplementary Figure 4A-C 必须在同一原子运行完成并通过完整性核验后统一写入图例。最终图例至少要报告：实际 `22/13/19` 样本数、158 个输入 TF、三端启动子开放 TF 数、HC-TF 数、motif 可检验 TF 数、一个/两个/三个系统 motif 支持数、HOMER adjusted p 阈值、半数样本阈值、JASPAR/CIS-BP 数据库版本，以及没有 canonical 峰的 PDX 如何保留在分母中。任何一个数字均不得从运行中间状态提前补写。
