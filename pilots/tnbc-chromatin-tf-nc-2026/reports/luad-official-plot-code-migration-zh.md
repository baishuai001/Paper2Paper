# LUAD图件迁移到TNBC官方绘图代码的审计说明

> 更新日期：2026-08-21
> 官方来源：CodeOcean capsule `7227095/v1`，commit `edf5314`
> 整体交付状态：`STRICT_OFFICIAL_PORT_REVIEW_ONLY`；`publication_ready=false`。最终7个评审PDF固定输出到`execution/luad-official-code-port/final/`，文件级凭据见[FINAL_REVIEW_SET.md](../execution/luad-official-code-port/FINAL_REVIEW_SET.md)。

## 结论

TNBC作者已经公开了绝大多数统计面板的绘图代码。此前LUAD重建把“数据格式适配”错误扩大为“重新设计图件”，导致字号、配色、布局和部分网络语义偏离锚点文献。此前由`code/publication_rebuild/`生成的自定义主图、补图、流程图、网络图及`anchor_style`拼版统一标记为**withdrawn / deprecated**，只保留历史追溯用途，不再进入论文交付。

严格端口采用以下不可变规则：

1. 有官方构造器的面板直接执行校验锁定的官方脚本或官方原始行块。
2. 只允许路径、LUAD/LUSC标签、数据集显示名、冻结样本数和适配对象名等已声明替换。
3. 不修改官方`theme`、字体、字号、颜色、`geom_*`、`scale_*`、聚类选项、图例布局或绘图几何。
4. 每个官方源文件锁定SHA-256；所有替换均核对出现次数，并输出unified diff和机器可读receipt。
5. capsule没有构造器的面板标记`MANUAL_REQUIRED`，不得用自绘图冒充官方代码复现。
6. LUAD缺少合格统计输入的面板标记`UNSUPPORTED`，不得生成替代结果。
7. 最终拼版只嵌入官方原子PDF；拼版本身不得重画、重着色或改变面板内部信息。

## 外部患者验证队列的冻结角色

- **主要验证：GSE41271微阵列**，183例LUAD、80例LUSC，复现70/158个TCGA发现TF。它占据TNBC代码中METABRIC的主要验证槽位。
- **次级验证：GSE81089 RNA-seq**，108例LUAD、67例LUSC，复现95/158个TF。它在隔离的次级官方代码运行中呈现，不替代主要微阵列验证。
- 两个外部队列共同复现44/158个TF。不同平台不拼接表达强度。

这一设计直接回应了此前“为什么没有像TNBC一样使用微阵列大队列”的问题：正式端口现在以GSE41271作为METABRIC类比验证，GSE81089保留为额外RNA-seq支持。

Figure 1的发现层与跨系统热图层必须分开记数：Figure 1B保留351个显著TF（158 LUAD＋193 LUSC）；Figure 1C–E只包含350个跨TCGA、PDMR PDX和DepMap可估计的TF。ZNF737的regulon虽有35个靶基因，但与DepMap表达矩阵交集为0，因此不做零填充或插补。

## Figure 2的冻结生物学结果

- 发现集合：158个LUAD特异性活性TF。
- 染色质信息支持集合：31个HC-TF。
- 具有可检验已知motif：20/31。
- 按冻结的作者式JASPAR优先规则，正式严格患者、PDX、细胞系三系统共同motif TF：FOXA3、NFATC4、XBP1，共3个。
- 早期candidate motif输入版本的三系统交集为0；它仅作“输入版本边界”敏感性分析，不得替代正式主分析的3个TF。
- q≤0.05且每系统至少50%样本支持的敏感性集合另加入ETV1、ZNF75D；它们不能并入严格主集合。
- Figure 3–5继续使用全部31个HC-TF，与TNBC作者使用全部94个HC-TF的逻辑一致。

因此，LUAD可以表述为“31-TF染色质信息调控网络具有选择性跨系统motif核心”，不能表述为“广泛三系统共同motif程序”。

## 官方代码覆盖和已知边界

| 面板范围 | 官方代码状态 | 边界 |
|---|---|---|
| Figure 1B–E、Supplementary Figure 2A–F | 有官方代码且端口`PASS` | GSE41271为主要验证；GSE81089为隔离次级验证。Figure 1B为351个显著TF，Figure 1C–E为350个可估计TF。Supplementary Figure 2F最终科学版使用一行审计过的`ANCHOR_BUGFIX`，严格错误重放仅作审计 |
| Figure 2B–E、Supplementary Figure 4B–C | 有官方代码且端口`PASS` | 正式三系统motif为FOXA3/NFATC4/XBP1；candidate版交集0仅作版本边界。Figure 2C为`CAPSULE_PRINT_MISMATCH`，且配对RNA/VIPER为21/22；Figure 2B/启动子闸门仍用22例 |
| Supplementary Figure 3A–L | 有官方代码且已严格运行 | 3A输入为`PARTIAL`；3B/3C IDR不适用；3D/3F/3K/3L存在官方固定标签重叠，整组`publication_ready: false` |
| Figure 3A–G、Supplementary Figure 4D、5A–B | 有官方代码且已严格运行 | Figure 3B为6节点/3边；Figure 3D为31节点/134边，无显示性裁边 |
| Figure 4A–H、Supplementary Figure 6A–H | 有官方代码且已严格运行 | 主森林图显示名义显著结果；只有TCGA OS的ZNF75D通过BH校正 |
| Figure 5A–F、Supplementary Figure 7A–C | 有官方代码且已严格运行 | Figure 5A/B/C顺序为PRISM/GDSC2/CTRPv2；Figure 5F高CI分支为`SUPPORTED_ZERO` |
| Figure 5G、Figure 5H、Supplementary Figure 7D | 有官方构造器但LUAD输入不合格 | 0/423个PDX检验进入匹配验证，三面板均为`UNSUPPORTED` |
| Figure 1A、Figure 2A、Supplementary Figure 1、Supplementary Figure 4A | capsule无构造器 | 全部`MANUAL_REQUIRED`；不再使用旧自绘图 |

## 当前执行状态

- **Figure 3及关联补图：`PASS`。** 19个PDF由官方脚本生成；源码、diff、输出和Figure 3B/3D拓扑均通过机器验证。
- **Figure 4–5及关联补图：`PASS`。** 31个官方源码块完成来源与切片校验，实际执行28个官方源码块，共生成39个原子PDF。Figure 5G、Figure 5H和Supplementary Figure 7D按真实输入状态记为`UNSUPPORTED`。
- **Supplementary Figure 3A–L：严格端口通过但不具备最终投稿状态。** 12个PDF均由官方构造器生成，绘图构造器重写数为0；receipt为`PASS_STRICT_PORT_WITH_DECLARED_INPUT_AND_VISUAL_BOUNDARIES`，同时明确`publication_ready: false`。
- **Figure 1–2及Supplementary Figure 2、4B–C：`PASS`。** 官方源码、适配器、主要/次级验证角色、正式motif版本和输出receipt均已冻结并通过。Figure 2C的RNA/VIPER配对为21/22，不对缺失样本插补。
- **最终7个PDF拼版：`STRICT_OFFICIAL_PORT_REVIEW_ONLY`。** 固定输出到`execution/luad-official-code-port/final/`并接受逐文件、逐页机器QA；由于Supplementary Figure 3明确`publication_ready: false`，不得将整体状态升格为投稿就绪。

## Figure 4–5 P0审计修订

2026-08-21的独立审计发现旧版Figure 4–5端口仍有四类实现偏差。修订版`figure45-v2-p0fix`已经完成云端复核，并取代旧`figure45`与v1压缩包：

- Supplementary Figure 5B不再把官方`Helvetica`替换为`Liberation Sans`；现在逐字执行官方热图块，并恢复热图后的FDR星号说明。源码字体替换数为0。
- Supplementary Figure 6D/6H只将图内TCGA终点显示文字从`RFS`改为冻结语义`DFS`；输入文件名和对象名继续保留原始`RFS`命名，避免把显示修正误扩展到数据路径。
- Figure 4A–D标题恢复官方语义`Multivariate {cohort} {endpoint}`。
- Supplementary Figure 7B的四个案例不再由wrapper自定义调用helper；四段官方调用与保存代码均被实际执行，并保留官方`cairo_pdf()`默认7×7英寸画布。

修订版机器验收为`PASS`、0项失败：31个来源块中30个可作为独立R代码严格解析；`fig5h_waterfall`是从官方循环内部截取、仅用于来源追溯且从不执行的`UNSUPPORTED`切片，因此其独立解析明确记为`R_PARSE_NOT_APPLICABLE`，而不是伪报解析通过。39个原子PDF全部具有逐文件SHA-256，原子PDF哈希清单受receipt锁定；任何与清单不一致的同名PDF均不得进入最终拼版。

## 两个必须保留的capsule问题

### Figure 2C印刷版不一致

公开脚本`Figure2/06C-PromoterAccessibility_MotifEnrichment.R`生成四块、TF位于纵轴的组合图，而论文印刷版Figure 2C采用不同的三块横向构图。严格端口执行公开脚本并标记`CAPSULE_PRINT_MISMATCH`；若以后需要贴近印刷版，只能建立独立的人工排版版本，不能称为零修改官方代码输出。

Figure 2C还有一个与排版无关的配对边界：22例ATAC患者中，21例能与RNA/VIPER配对；TCGA-44-A47F不在可用RNA活性对象中。因此Figure 2C配对层是21/22，Figure 2B和启动子闸门仍是22/22的冻结ATAC队列；不使用零填充或其他插补。

### Supplementary Figure 2F对象引用错误

官方源码第475行把PDX长表从TCGA对象生成。为同时满足代码审计和科学正确性，端口保留两套输出：

- `CAPSULE_BUG_RETAINED`：严格重放官方错误，仅供审计；
- `ANCHOR_BUGFIX`：只把第475行的TCGA对象名改为PDX对象名，作为最终科学输出。

除此之外不改变该面板的绘图构造器。

## Supplementary Figure 3的输入和视觉边界

- 3A使用公开固定GDC open-count矩阵，不是重新逐样本调用的MACS2/IDR峰，因此为`PARTIAL`。
- 3B和3C分别使用PDX、细胞系MACS2 q=0.01规范峰；缺少合格真实生物学重复，IDR为不适用。
- 3D–F使用冻结的1,000次置换峰饱和对象；3G–I使用冻结的共识可及性Pearson矩阵。
- 3J–L使用真实hg38 500-bp全基因组背景，共6,062,095个bin。
- 3K的13个冻结PDX中有1个零峰，官方非空长表因此显示n=12。
- 3D、3F、3K和3L这4个面板因官方固定标签坐标发生文字重叠。严格端口不移动标签，所以该组明确`publication_ready: false`；后续若修订，只能明确标记为post-port layout调整。

其他官方固定坐标在LUAD输入上也留下了可复核的排版边界：Figure 1B和主要验证版Supplementary Figure 2A的`stat_cor`文字在绘图区右边界被截；GSE81089次级Supplementary Figure 2A的高密度基因标签发生重叠；Figure 4B的长标题较拥挤；Supplementary Figure 4D的小热图标题与Pearson图例标题叠印。完整数值均保留在结果表中，但这些面板仍需在严格端口归档后另做、且明确披露的post-port layout修订。

## 旧版Figure 3B/3D错误及修正

旧版Figure 3B在官方规则得到3条边、6个端点后，额外加入24个无边孤立靶基因；旧版Figure 3D又用80%边权分位数裁掉101/134条边，并把节点大小映射为regulon靶基因数而非官方网络degree。这些改动同时改变了视觉和科学含义。

严格官方代码结果为：

- Figure 3B：6个节点、3条边，不添加无边节点；
- Figure 3D：31个节点、134条正边、1个连通分量、0个孤立点；节点大小为degree，边宽为共享激活靶基因数。

旧版Figure 3B、3D及其派生拼版均已撤回。

## 交付位置和判定规则

约定最终目录为：

`pilots/tnbc-chromatin-tf-nc-2026/execution/luad-official-code-port/final/`

计划文件为：

1. `Figure1.pdf`
2. `Figure2.pdf`
3. `Figure3.pdf`
4. `Figure4.pdf`
5. `Figure5.pdf`
6. `Supplementary_Figures.pdf`
7. `TNBC_vs_LUAD_comparison_atlas.pdf`

当前已核验的官方端口原子结果覆盖Figure 1–5及对应补图；Figure 3、Figure 4–5 P0修订版和Supplementary Figure 3分别归档于`execution/luad-official-code-port/figure3/`、`figure45-v2-p0fix/`和`supp3/`。旧`figure45/`已标记为**SUPERSEDED**，不得继续用于论文。最终7个评审PDF固定输出到`execution/luad-official-code-port/final/`，页数、字节数和SHA-256见[FINAL_REVIEW_SET.md](../execution/luad-official-code-port/FINAL_REVIEW_SET.md)；固定路径与机器校验通过不改变证据等级，整体仍为`STRICT_OFFICIAL_PORT_REVIEW_ONLY`、`publication_ready=false`。旧目录`execution/luad-publication-rebuild/results/publication/`及其中的自绘PDF也已**SUPERSEDED / WITHDRAWN**，不再是正式输出位置。
