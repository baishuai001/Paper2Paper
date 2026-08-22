# LUAD 主分析分支重算报告：FDR05 → HC45 → shared≥2 → nominal log-rank

生成日期：2026-08-22
分析标识：`luad-fdr05-hc45-primary-v1`（云端实际目录因首次占位重试为 `luad-fdr05-hc45-primary-v2`，两者不是两套科学分析）

## 1. 直接结论

用户指定的三项规则均已完整执行，代码与输出通过技术审计：

1. Figure 1 入口改为 **BH FDR≤0.05、NES>0、表达 logFC>0**；镜像的 LUSC 规则为 FDR≤0.05、NES<0、表达 logFC<0。TCGA-LUAD 得到 **187 个候选 TF**。
2. 187 个 TF 经原文式三系统启动子开放性与平均活性筛选后得到 **45 个 HC-TF**，因此 Figure 2 已随入口变化重算。
3. Figure 3B 按“两个靶基因至少共享 2 个激活型 HC-TF”重绘，得到 **237 个节点、480 条边**。
4. Figure 4 主图按 **未校正 log-rank P≤0.05** 决定显示，45 个 HC-TF 全部进入生存筛选；BH 和 maxstat 校正结果只保留在补充材料中，不用于删掉主图结果。

这里的“通过”指：数据、代码、计数、PDF 和审计记录均按冻结规则成功生成；它不等于所有生物学结论都已获得独立验证。Figure 3B 仍是描述性网络，Figure 4 仍是探索性预后筛选。

## 2. 本次冻结规则到底是什么

| 环节 | 主分析规则 | 没有改变的设置 | 主图如何判定 |
|---|---|---|---|
| Figure 1 | BH FDR≤0.05；LUAD 为 NES>0 且表达 logFC>0；LUSC 为 NES<0 且表达 logFC<0 | 外部队列只做同规则复现，不反向缩减 TCGA 发现集 | TCGA 的 187 个 LUAD TF 进入 Figure 2 |
| Figure 2 | 启动子在患者、PDX、细胞系各至少一半样本开放；只排除三个系统平均 NES 均小于 0 的 TF | motif 仍按 JASPAR 优先、无 JASPAR 时用 CIS-BP；HOMER q≤1×10⁻⁵且各系统至少半数样本 | 得到 45 个 HC-TF；motif 是下游支持证据，不再缩减 Figure 3/4 入口 |
| Figure 3B | 两个靶基因共享的激活型 HC-TF 数量≥2 | 网络构造、Louvain 着色和版式来自 TNBC 官方 Figure3/04、05 代码 | 达阈值的边全部展示 |
| Figure 4 | 数据选择切点后，原始 log-rank P≤0.05 | 仍保留 120 个月行政截尾、最大选择切点、每组至少 20% 样本 | 未校正 P≤0.05 决定主图；校正值仅作敏感性分析 |

因此，本次 Figure 4 的“未校正”只改变**显著性判定**，没有偷偷把生存时间、切点算法或最小组比例改掉。如果要取消最大选择切点，改成固定中位数二分，那是另一套分析，不是本次规则。

## 3. Figure 1：187 个 TF 是怎样得到的

### 3.1 发现集

- TCGA：LUAD 特异 187 个；LUSC 镜像特异 248 个。
- GSE81089 RNA-seq：LUAD 137 个；LUSC 179 个。
- GSE41271 微阵列：LUAD 120 个；LUSC 171 个。

这些数是在每个队列内部独立应用相同的 FDR、NES 和表达方向规则得到的，不是把不同平台的数据直接混合检验。

### 3.2 外部复现

以 TCGA 的 187 个 LUAD TF 为分母：

- 88/187 在 GSE81089 中同方向复现；
- 60/187 在 GSE41271 中同方向复现；
- 33/187 同时在两个外部队列复现。

Figure 2 的入口仍是 **187**，不是 33。原因是外部队列在这里承担“独立复现证据”的角色，不是第二、第三道发现集过滤器；否则会把平台差异和样本量差异误当成生物学排除条件。

## 4. Figure 2：为什么从 187 变成 45 个 HC-TF

在 187 个候选 TF 中：

- 170 个在患者系统满足启动子开放规则；
- 59 个在 PDX 系统满足启动子开放规则；
- 63 个在细胞系系统满足启动子开放规则；
- 46 个同时通过三个系统的启动子开放规则；
- 其中 1 个 TF 在患者、PDX、细胞系三个系统的平均 NES 均小于 0，被原文式活性规则排除；
- 最终得到 **45 个 HC-TF**。

相对旧的 31 个 HC-TF，新增加的 14 个是：ARID4A、CALCOCO1、CREBL2、EPC1、MYBL1、PCGF5、SLC30A9、TADA2B、ZNF34、ZNF438、ZNF652、ZNF672、ZNF799、ZSCAN2。完整名单见 `execution/luad-fdr05-hc45-primary-v1/results/figure2/audit/primary_fdr05_hc45/HC45.tsv`。

### 4.1 motif 结果随之怎样变化

- 45 个 HC-TF 中有 25 个具有可测试的已知 motif；
- 19 个至少在患者、PDX、细胞系中的一个系统获得正式 motif 支持；
- 按作者的 JASPAR 优先规则，三系统共同支持为 **5 个：CREBL2、FOXA3、MYBL1、NFATC4、XBP1**；
- 把 JASPAR 与 CIS-BP 的全部可用 motif 做并集仅作为敏感性分析时，三系统共同支持为 6 个。

后续 Figure 3 和 Figure 4 使用全部 45 个 HC-TF，而不是只使用 5 个三系统 motif TF；这与 TNBC 文章使用全部 HC-TF 进入网络和临床分析的逻辑一致。

## 5. Figure 3B：共享≥2个 TF 后到底画出了什么

严格定义如下：

- 节点：被 HC-TF 激活的靶基因；
- 边：两个靶基因至少被 **2 个相同的 HC-TF** 激活；
- 颜色：官方代码的 Louvain 网络分组，只是按网络拓扑着色。

结果：

- 阈值前有 310 个激活型共享靶基因；
- 阈值后网络含 237 个节点、480 条边；
- 共有 25 个连通分量，最大连通分量含 142 个节点；
- 任意一条边最多共享 3 个 HC-TF。

这说明“共享≥2”能够形成可见网络，但不能直接写成“发现了经过统计验证的高阶调控 community”。原因不是主观否定，而是此前针对同一网络做的保持节点度数置换检验中，没有靶基因对在 BH FDR 后显著。换句话说：当前 Figure 3B 是按用户冻结阈值得到的**描述性共调控投影**；颜色模块可以描述，不能当成独立统计检验的阳性结果。

## 6. Figure 4：未校正 log-rank P≤0.05 的结果

“未校正 log-rank P≤0.05”最直白的含义是：对每个 TF，在选出的高/低活性切点上比较两条生存曲线，只看这个原始 P 值是否≤0.05。本次一共检查 45 个 TF，而且切点也是从数据中选择的，所以这些数是**筛选命中数**，不是已经验证的预后标志物数量。

### 6.1 主分析命中数

| 队列 | 终点 | 45 个 HC-TF 中未校正 P≤0.05 | 仅供敏感性参考：BH FDR≤0.05 | 仅供敏感性参考：Lau92 校正 P≤0.05 |
|---|---:|---:|---:|---:|
| GSE41271-LUAD | OS | 29 | 17 | 9 |
| GSE41271-LUAD | RFS | 20 | 4 | 4 |
| TCGA-LUAD | OS | 34 | 33 | 17 |
| TCGA-LUAD | DFS | 14 | 0 | 2 |

主图 A–D 展示第一列的未校正结果。后两列没有被用来否决主图，只用于告诉读者结果对多重检验和切点选择有多敏感。

### 6.2 同一队列两个终点的重叠

- GSE41271 有利方向：13 个 TF 同时关联 OS 与 RFS——CALCOCO1、CREBL2、CREBRF、E2F5、ETV1、FOXD2、MYBL1、ZBTB18、ZNF141、ZNF33A、ZNF438、ZNF540、ZNF75D。
- GSE41271 不利方向：PHC2 同时关联 OS 与 RFS。
- TCGA 有利方向：8 个 TF 同时关联 OS 与 DFS——CALCOCO1、CREBRF、CRY2、NR3C2、TADA2B、ZNF19、ZNF33B、ZNF540。
- TCGA 不利方向：没有 TF 同时关联 OS 与 DFS。

### 6.3 Kaplan–Meier 示例

- GSE41271 ZNF540：OS P=0.000132；RFS P=0.000554。
- GSE41271 PHC2：OS P=0.0400；RFS P=0.0187。
- TCGA CRY2：OS P=4.20×10⁻⁶；DFS P=0.00395。

这些例子均满足用户指定的未校正 log-rank P≤0.05。它们是否在更严格校正后仍成立，已经在补充表中保留，但不改变本次主图入选规则。

## 7. 绘图代码与计算位置

- 数据下载、TF 活性、染色质、motif、网络与生存分析全部在云服务器完成。
- Figure 1/2 的统计图来自 TNBC Code Ocean 官方绘图代码移植；官方未提供的流程示意图才由项目中的声明式示意图代码生成。
- Figure 3B 使用校验和锁定的 TNBC Figure3/04 与 Figure3/05 代码，只把共享阈值改为 2，并接入 45 个 LUAD HC-TF。
- Figure 4 的 forest、Euler、Kaplan–Meier 和 volcano 构造函数来自校验和锁定的 TNBC 官方代码；LUAD 包装器只提供数据、队列名称与终点名称。为避免 PDF 标题越界，最终版只缩短了显示标题，没有改任何统计量。
- 本地仅下载云端结果、做矢量 PDF 拼接和视觉质检，没有在本地重跑生物信息学分析。

核心代码位于 `code/primary_fdr05_hc45/`；主要机器审计记录位于 `execution/luad-fdr05-hc45-primary-v1/results/`。

## 8. 最终 PDF 与图例

### Figure 1

`Figure1_primary_FDR05.pdf`。A，LUAD TF 发现与外部复现工作流。B，TF 活性 NES 与表达 logFC 的一致性。C，患者肿瘤中的 LUAD/LUSC TF 活性及临床/分子注释。D，PDX 中相同 TF 程序的活性。E，肺癌细胞系中相同 TF 程序的活性。

### Figure 2

`Figure2_primary_HC45.pdf`。A，从 187 个 LUAD TF 经三系统启动子开放和活性规则获得 45 个 HC-TF。B，患者、PDX 与细胞系的启动子开放交集。C，45 个 HC-TF 的 motif 支持矩阵。D，各系统 motif 支持的交集。E，19 个至少获得一个系统 motif 支持的 HC-TF，以及 5 个三系统共同支持 TF 的完整范围和局部放大。

### Figure 3B

`Figure3B_shared_ge2.pdf`。共享至少 2 个激活型 HC-TF 的靶基因网络。节点为靶基因，边表示共享调控，节点颜色为 Louvain 拓扑分组，节点大小表示连接数。

### Figure 4

`Figure4_primary_unadjusted_logrank.pdf`。A–D，GSE41271 OS/RFS 与 TCGA OS/DFS 中未校正 log-rank P≤0.05 的 HC-TF 及二分组 HR。E、G，同一队列两个终点中有利或不利 TF 的重叠。F、H，ZNF540、PHC2 与 CRY2 的 Kaplan–Meier 示例。

### Supplementary Figures

- `SupplementaryFigure1_primary_FDR05.pdf`：三个患者队列的纳入流程与最终 LUAD/LUSC 数量。
- `SupplementaryFigure2_primary_FDR05.pdf`：外部队列复现、队列间 TF 重叠、相关性热图及跨系统 TF 活性分布。
- `SupplementaryFigure4_HC45.pdf`：187→45 的染色质/活性闸门、三系统平均 NES 与 motif 数据库可用性。
- `SupplementaryFigure6_logrank_sensitivity.pdf`：四个终点的单变量与多变量敏感性图；不用于改变主图的未校正入选规则。

## 9. 质量控制

- Figure 1/2 官方代码端：51 个 PDF 原子图，Poppler 渲染检查全部通过。
- Figure 4 最终标题修正版：22 个 PDF 原子图，Poppler 渲染检查 22/22 通过。
- 最终拼接：8 个 PDF 均为单页、矢量嵌入，PDFium 解析与渲染 8/8 通过；逐页联系表已肉眼检查。
- 嵌入来源与 SHA256：`execution/luad-fdr05-hc45-primary-v1/output/primary_delivery_embedding_manifest.tsv` 和 `primary_delivery_sha256.tsv`。
- PDF 质检记录：`execution/luad-fdr05-hc45-primary-v1/output/primary_delivery_pdf_qa.tsv`。
