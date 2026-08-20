# LUAD论文重建：科学叙事、逐面板图注与TNBC对照审计（严格官方代码端口版）

> 版本日期：2026-08-21
> 官方来源：TNBC CodeOcean capsule `7227095/v1`，commit `edf5314`
> 重要声明：此前由 `code/publication_rebuild/` 自定义渲染器生成的主图、补图、流程图和 `anchor_style` 拼图已经全部撤回（**withdrawn / deprecated**），不得再作为论文图件、图注依据或“已完成结果”引用。本文只描述严格官方代码端口及其明确边界。
> 当前整体交付等级：`STRICT_OFFICIAL_PORT_REVIEW_ONLY`，`publication_ready=false`。Figure 1–5及对应补图的代码端口状态不等于投稿就绪；最终7个评审PDF固定输出到`execution/luad-official-code-port/final/`，文件级凭据见[FINAL_REVIEW_SET.md](../execution/luad-official-code-port/FINAL_REVIEW_SET.md)。

## 一、当前结论和证据强度

本研究要检验的是：由大规模患者转录组发现的LUAD转录因子（TF）活性状态，能否经染色质可及性优先筛选，形成可在患者肿瘤、患者来源异种移植模型（PDX）和细胞系中继续研究的调控网络，并进一步产生可检验的预后及药物反应线索。

冻结后的证据链为：

1. 在TCGA-LUAD与TCGA-LUSC之间进行ARACNe3/msVIPER分析，得到351个显著发现TF（158个LUAD特异性活性TF＋193个LUSC特异性TF）。Figure 1B保留全351个；Figure 1C–E的跨系统热图仅使用350个可估计TF，因为ZNF737的35个regulon靶基因与DepMap表达矩阵的交集为0，故不进行零填充或插补。
2. **主要独立验证队列为GSE41271微阵列队列**：183例LUAD、80例LUSC，复现70/158个发现TF。该队列在官方代码端口中占据TNBC文献METABRIC的主要验证槽位。
3. **次级验证队列为GSE81089 RNA-seq队列**：108例LUAD、67例LUSC，复现95/158个发现TF。它在隔离的次级官方代码运行中报告，不替代GSE41271的主要验证角色。
4. GSE41271和GSE81089共同复现44/158个TF。两个平台的表达值不拼接，也不将较高的95/158替代微阵列主要验证结果70/158。
5. 采用TNBC式“启动子可及性＋三系统平均TF活性”规则，从158个发现TF中得到31个染色质信息支持TF（HC-TF）。
6. 31个HC-TF中有20个具备可检验的已知motif；按冻结的作者式“JASPAR优先、无JASPAR时才用CIS-BP”输入版本，正式严格三系统motif集合为FOXA3、NFATC4和XBP1，共3个。早期candidate motif输入版本的三系统交集为0，它只作为“输入版本边界”敏感性分析，不得替代正式结果。另外，预先声明的放宽统计规则（HOMER q≤0.05且每个系统至少50%样本支持）再加入ETV1和ZNF75D，但这两个TF不能并入严格主集合。
7. 与TNBC作者在Figure 3–5中使用全部94个HC-TF一致，LUAD的Figure 3–5使用全部31个HC-TF，而不是只使用3个严格三系统motif TF。
8. 临床和治疗证据仍属探索性：冻结分析中只有TCGA总体生存（OS）的ZNF75D通过队列内BH校正；5组“药物－TF”关系在至少两个细胞系药物资源中同方向重复，但0/423个公共PDX检验满足进入匹配体内验证的条件。因此Figure 5G、Figure 5H及Supplementary Figure 7D为`UNSUPPORTED`，不是阴性生物学结论，更不能用自绘漏斗或空瀑布图代替。

当前可以辩护的中心主张是：

> 染色质信息辅助的TF活性框架定义了LUAD的31-TF调控网络，并识别出选择性跨患者、PDX和细胞系保守的motif核心；其预后和药物反应关联提供后续验证假设，而不是已经完成的临床或体内证据。

当前结果**不支持**“LUAD存在类似TNBC的广泛三系统共同motif程序”这一更强表述。

## 二、严格官方代码端口的状态边界

状态词只按机器凭据使用：`PASS`表示官方源码校验、运行及输出验证均已完成；`PENDING`表示尚未完成云端全运行；`MANUAL_REQUIRED`表示capsule没有该面板的绘图构造器；`UNSUPPORTED`表示LUAD没有满足官方构造器输入条件的数据；`CAPSULE_PRINT_MISMATCH`表示公开代码输出与论文印刷版排布不一致。

| 范围 | 当前状态 | 可以说什么 | 不能说什么 |
|---|---|---|---|
| Figure 1B–E、Figure 2B–E、Supplementary Figure 2、Supplementary Figure 4B–C | `PASS` | 官方源码、适配器、全运行输出和receipt已校验；GSE41271为主要、GSE81089为次级；正式motif三系统交集为FOXA3/NFATC4/XBP1 | `PASS`只指代码端口通过，不把整体交付升格为投稿就绪 |
| Figure 1A、Figure 2A、Supplementary Figure 1、Supplementary Figure 4A | `MANUAL_REQUIRED` | 官方capsule没有构造器，最终严格版只能明确标记缺口 | 不得把自绘流程图或筛选图称为“官方代码复现” |
| Figure 2C | `PASS` + `CAPSULE_PRINT_MISMATCH` | 严格端口保留capsule的四块纵向输出；22例ATAC患者中21例有可配对RNA/VIPER，未配对的TCGA-44-A47F不插补 | 不得自行重排成论文印刷版并称为官方输出，也不得把Figure 2C的21例误写成Figure 2B/启动子闸门的分母 |
| Figure 3及Supplementary Figure 4D、5A–B | `PASS` | 19个PDF由校验锁定的官方脚本生成；Figure 3B/3D拓扑另经验证 | 不得再混入旧自绘网络、裁边或人工孤点 |
| Supplementary Figure 3A–L | `PASS_STRICT_PORT_WITH_DECLARED_INPUT_AND_VISUAL_BOUNDARIES` | 12个面板均运行官方构造器，绘图构造器重写数为0 | 整组仍为`publication_ready: false`，不能说已达到最终投稿视觉质量 |
| Figure 4A–H、Figure 5A–F、Supplementary Figure 5–7A–C | `PASS` | P0修订版实际执行28个官方源码块，生成39个原子PDF；全部有逐文件SHA-256 | 不能把`SUPPORTED_ZERO`分支解释为有正结果，也不得回用已superseded的v1输出 |
| Figure 5G、Figure 5H、Supplementary Figure 7D | `UNSUPPORTED` | 没有符合条件的LUAD PDX验证输入 | 不得生成替代图，也不得解释为药物关联在体内已被否定 |
| 最终主图、补图和TNBC-vs-LUAD对照图谱拼版 | `STRICT_OFFICIAL_PORT_REVIEW_ONLY` | 只嵌入经来源锁定的官方原子PDF；固定输出到`execution/luad-official-code-port/final/`并以逐文件SHA和严格QA锁定 | 不得称`publication_ready`或“最终投稿版” |

## 三、主图Figure Legends

以下图注基于已通过的严格官方代码端口。`PASS`描述代码来源、执行和输出验证；整体交付仍只能称`STRICT_OFFICIAL_PORT_REVIEW_ONLY`。

### Figure 1｜LUAD TF活性程序的发现、独立验证和跨模型投射（Figure 1A为`MANUAL_REQUIRED`；B–E为`PASS`）

**A，** `MANUAL_REQUIRED`。研究设计面板拟说明TCGA发现、GSE41271主要微阵列验证、GSE81089次级RNA-seq验证，以及向PDMR PDX和DepMap细胞系的投射关系。TNBC公开capsule没有Figure 1A构造器，因此严格官方代码版不绘制此面板；若以后人工重建，必须明确标记为非官方构造器。

**B，** TCGA-LUAD与TCGA-LUSC比较中，基因表达log2 fold change与msVIPER TF活性NES的关系。横、纵零线、回归线、Pearson相关及TF分类颜色均沿用官方Figure 1B构造器。该图检验TF活性是否提供不同于TF自身mRNA丰度的信息。

**C，** Figure 1B的发现层包含351个显著TF；其中350个在TCGA、PDMR PDX和DepMap细胞系中均有真实活性值，用于TCGA-LUAD/LUSC患者的样本级活性热图。ZNF737由于DepMap中无可估计靶基因而不进入Figure 1C–E，不做0填充或插补。热图以发现类别标记LUAD特异、LUSC特异或非显著TF；158个LUAD特异性TF是研究主集合。官方ComplexHeatmap代码执行逐TF标准化、固定行序和样本聚类，LUAD临床、组织学和分子注释字段只通过数据适配器映射。

**D，** 将同一350-TF跨系统可用集合投射至PDMR肺癌PDX转录组后得到的样本级TF活性热图；其中158个LUAD特异性TF仍由行注释识别。样本聚类、TF行注释、热图和注释条构造均来自官方Figure 1D代码；PDX诊断和状态注释不等同于患者临床亚型。

**E，** 将同一350-TF跨系统可用集合投射至DepMap肺癌细胞系后的TF活性热图。官方Figure 1E构造器保持不变；细胞系谱系和状态标签仅说明模型来源，不能视为患者肿瘤的等价替代。

### Figure 2｜染色质优先筛选和三系统motif证据（Figure 2A为`MANUAL_REQUIRED`；B–E为`PASS`）

**A，** `MANUAL_REQUIRED`。拟说明22例患者、13个LUAD PDX和19条LUAD细胞系的ATAC数据设计。公开capsule没有Figure 2A绘图构造器，严格版不自行重绘。

**B，** 官方ComplexUpset构造器显示158个发现TF在患者、PDX和细胞系中的启动子开放集合及其交集。每个系统的主规则为启动子在至少一半样本中开放，对应患者≥11/22、PDX≥7/13、细胞系≥10/19。

**C，** `CAPSULE_PRINT_MISMATCH`。严格端口执行公开脚本`Figure2/06C-PromoterAccessibility_MotifEnrichment.R`，其输出为TF位于纵轴的四块组合：HC状态、四个转录组系统的平均TF活性、三系统逐样本启动子可及性和三系统逐样本motif证据。主要外部患者槽位为GSE41271，而不是GSE81089。22例ATAC患者中21例可与RNA/VIPER一一配对；TCGA-44-A47F无可用RNA/VIPER，故Figure 2C的配对部分使用21例。Figure 2B和启动子闸门仍使用全22例分母，不插补第22例。公开脚本排布与论文印刷版不一致；严格端口保留capsule输出并标明这一边界。

**D，** 官方ComplexUpset构造器显示31个HC-TF在患者、PDX和细胞系中达到主要motif规则的集合归属。冻结作者优先版本的正式三系统共同集合为FOXA3、NFATC4和XBP1。早期candidate输入版本的交集为0，仅作版本敏感性边界，不进入正式面板。图形面积或交集柱高只表示当前LUAD数据的集合计数，不代表证据强度超过TNBC。

**E，** 官方示例构造器显示预先选定LUAD TF在三个系统中的motif富集和样本支持。数据库选择执行JASPAR优先、无JASPAR时才使用CIS-BP；严格主结果与q≤0.05/50%敏感性结果必须分开标注。

### Figure 3｜31个HC-TF的regulon结构与肿瘤间异质性（`PASS`）

**A，** 31个HC-TF的regulon组成，包括激活/抑制以及共同/特有靶基因类别；附带官方motif和净效应注释方格。Figure 3使用全部31个HC-TF，motif方格不是再次筛掉TF的入口。

**B，** 共同激活靶基因网络。按官方规则，只有共同受到至少3个激活HC-TF调控的靶基因对才连边，且`graph_from_data_frame`不额外接收孤立节点。LUAD严格网络包含6个节点、3条边和3个双节点组分；该结果不支持广泛汇聚的靶基因群落，也不宜对双节点组分强行进行GO命名。

**C，** TCGA-LUAD中31个HC-TF活性的Pearson相关热图，以及官方代码同时生成的最高/最低相关TF摘要条形图。相关性说明共同活动结构，不直接证明TF间物理互作。

**D，** 共享至少1个激活靶基因的HC-TF协作网络。严格官方规则得到31个节点、134条正重叠边、1个连通分量和0个孤立点；节点大小表示网络degree，边宽表示共享激活靶基因数。旧版额外采用80%边权分位数裁边的网络已撤回。

**E，** regulon大小与协作TF数量之间的关系，使用与Figure 3D相同的完整边表。

**F，** TCGA-LUAD患者中31个HC-TF活性分布的偏度。显著性采用官方偏度近似检验并进行BH校正；偏态不能被直接解释成离散分子亚型。

**G，** CRY2、MAGED2、ZNF19和NR3C2的代表性TF活性分布。四个例子由完整Figure 3D网络中degree最高的TF确定性选择，仅替换TNBC脚本中硬编码的示例TF名称；密度、分面和Wilcoxon构造保持官方代码不变。

### Figure 4｜31个HC-TF的探索性预后关联（`PASS`）

**A–D，** 官方显著结果森林图构造器分别显示GSE41271-LUAD OS、GSE41271-LUAD无复发生存（RFS）、TCGA-LUAD OS和TCGA-LUAD无病生存（DFS）的多变量Cox结果。冻结分析对全部31个HC-TF按队列/终点中位数二分并截断10年随访，协变量为年龄、性别、病理分期和吸烟史；主图只绘制名义p≤0.05的TF，分别为3、4、10和2个。点和线表示log2(HR)及95%置信区间。只有TCGA OS中的ZNF75D通过队列内BH校正，因此主图中的其他名义结果不得称为经多重校正的预后标志物。

**E，** GSE41271中名义显著的多变量OS与RFS TF集合，分别按风险和保护方向显示官方Euler图。

**F，** GSE41271中ETV1和ZNF444的OS及RFS Kaplan–Meier曲线。曲线使用冻结的最大选择秩切点，与主要Cox筛查的中位数二分不同，只作为描述性展示。

**G，** TCGA-LUAD中名义显著的多变量OS与DFS TF集合，分别按风险和保护方向显示官方Euler图。

**H，** TCGA-LUAD中CREBRF和PHC2的OS及DFS Kaplan–Meier曲线。曲线使用冻结切点，只用于展示，不构成独立队列复现。

### Figure 5｜细胞系药物关联及缺失的PDX验证环节（A–F为`PASS`；G–H为`UNSUPPORTED`）

**A–C，** 官方concordance-index火山图构造器依次显示PRISM（23,591项检验）、GDSC2（3,627项检验）和CTRPv2（11,867项检验）中31个HC-TF活性与药物反应的关系。每个资源内部独立进行BH校正；不同资源的药物和细胞系覆盖不同，不能仅比较显著点数量。

**D，** 官方UpSet构造器显示资源内FDR≤0.05的“药物－TF”关联在PRISM、GDSC2和CTRPv2之间的集合归属。进入跨资源重复集合要求同一规范化药物－TF组合在至少两个资源出现；方向一致性由冻结统计表另行核实。

**E，** 31个HC-TF的跨队列三状态预后分类计数矩阵。横轴为RFS分类（复发、无复发或非预后），纵轴为OS分类（不良、良好或非预后）；只有在GSE41271和TCGA中均名义显著且方向一致的TF才进入良好/不良预后类别。该面板不是“5组重复药物效应热图”。

**F，** 官方低CI药物－TF点阵将跨资源重复关联按药物通路分组，并并列显示TF的OS/RFS预后类别。当前有低CI主通路输入；高CI主通路分支为`SUPPORTED_ZERO`，不得补画。冻结结果中的5组同方向重复关系为5-fluorouracil–ZNF254、5-fluorouracil–ZNF540、PHA-793887–WWC2、vorinostat–ZNF254和dactolisib–ZNF254。

**G，** `UNSUPPORTED`。423个公共PDX检验中，0个同时对应“细胞系跨资源重复、效应同方向且存在匹配PDX药物反应”的药物－TF组合，官方PDX concordance构造器因此没有合格输入。本面板仅保留明确状态说明，不绘制旧版自定义漏斗。

**H，** `UNSUPPORTED`。没有通过验证的LUAD PDX药物－TF组合可进入官方疗效瀑布图。本缺口表示公共数据覆盖不足，不证明5组细胞系关联在体内无效。

## 四、Supplementary Figure Legends

### Supplementary Figure 1｜样本纳入或研究流程（全部`MANUAL_REQUIRED`）

TNBC公开capsule未提供Supplementary Figure 1的绘图构造器。严格官方代码集合只显示`MANUAL_REQUIRED`说明，不再使用此前自绘的TCGA、GSE41271或GSE81089纳入流程。样本数仍由机器可读manifest报告：GSE41271主要比较为183例LUAD与80例LUSC；GSE81089次级比较为108例LUAD与67例LUSC。

### Supplementary Figure 2｜独立患者验证和跨模型内部一致性（`PASS`）

**A，** 主要面板使用GSE41271微阵列数据，按官方散点构造器展示差异表达log2 fold change与msVIPER NES的关系；隔离的次级页面以同一构造器展示GSE81089 RNA-seq结果。GSE41271复现70/158，GSE81089复现95/158。官方`stat_cor`使用固定坐标，导致主要面板右下角统计文字在LUAD范围内被边界截断；次级页面的若干高密度基因标签发生重叠。统计值保留在本报告和结果表中，图形需标记为待post-port layout修订。

**B，** 官方Euler构造器分别显示TCGA发现集合与GSE41271主要验证集合的重叠，以及隔离次级页面中的TCGA与GSE81089重叠；两个外部队列共同复现44/158个TF。

**C–E，** 官方Purple–Green ComplexHeatmap构造器分别显示TCGA患者、PDMR PDX和DepMap细胞系中的样本间TF活性Pearson相关性及相应样本注释。

**F，** 158个发现TF在TCGA、PDMR PDX、DepMap细胞系和GSE41271中的样本级NES分布。官方源码第475行误把PDX长表从TCGA对象生成；最终科学输出采用仅改正该对象名的一行`ANCHOR_BUGFIX`版本，严格保留错误的重放只用于审计，不进入最终补图。

### Supplementary Figure 3｜ATAC峰数量、饱和度、相关性和基因组分布（严格端口已运行；整体`publication_ready: false`）

**A，** 22例患者的可及区域数量。该面板状态为`PARTIAL`：输入来自公开的固定GDC open-count矩阵，不是本项目重新逐样本调用的MACS2/IDR峰。

**B，** 13个冻结LUAD PDX模型的MACS2 q=0.01规范峰数量；每个模型只有一个公共生物学文库，IDR为不适用。3个模型的规范峰数为0，保留为0而不补造峰。

**C，** 19条LUAD细胞系的MACS2 q=0.01规范峰数量；冻结队列没有合格的真实生物学重复集合，IDR为不适用。

**D–F，** 分别为患者、PDX和细胞系的1,000次置换峰饱和曲线及官方渐近拟合标注。D和F因官方固定文字坐标在LUAD数据范围内发生标注重叠；严格端口没有移动文字，所以这两个面板尚不具备最终投稿视觉质量。

**G–I，** 分别为患者、PDX和细胞系的共识可及性信号Pearson相关热图，使用官方ComplexHeatmap构造器。

**J–L，** 分别为患者、PDX和细胞系峰的基因组注释分布，并使用真实hg38 500-bp全基因组背景（6,062,095个bin）。K中冻结13个PDX模型有1个零峰，官方非空长表因此显示n=12。K和L的官方固定n标签与分布重叠，需在严格端口归档后另行声明版式修订，不能把修订版冒充零改动官方输出。

### Supplementary Figure 4｜HC-TF规则、motif数据库和增殖关联

**A，** `MANUAL_REQUIRED`。公开capsule没有该筛选流程面板的构造器，严格版不自行绘制。

**B，** `PASS`。官方ComplexUpset构造器显示31个HC-TF在患者、PDX和细胞系中的平均活性类别及交集；使用22/13/19的冻结分母。

**C，** `PASS`。官方Euler构造器显示31个HC-TF的JASPAR和CIS-BP motif可用性；JASPAR优先、CIS-BP后备。

**D，** `PASS`。官方volcano和heatmap构造器显示31个HC-TF活性与MKI67表达在TCGA患者、GSE41271患者、PDMR PDX和DepMap细胞系中的Pearson相关性；显著标记使用官方`|r|>0.4`和BH FDR<0.05规则。严格端口保留了官方小热图的固定标题/图例布局，因此热图标题与Pearson图例标题存在叠印；这不改变相关系数或FDR，但使该原子图仍不具备投稿排版质量。

### Supplementary Figure 5｜跨系统偏度和肿瘤微环境混杂（`PASS`）

**A，** 官方构造器比较31个HC-TF在PDMR PDX和DepMap细胞系中的活性偏度，并与Figure 3F患者结果共同解释跨系统异质性。

**B，** 官方volcano及heatmap构造器显示TCGA-LUAD和GSE41271-LUAD中HC-TF活性与ESTIMATE免疫、基质及肿瘤纯度指标的相关性。该分析用于识别患者TF信号是否可能受到非肿瘤细胞含量混杂，不证明细胞来源。

### Supplementary Figure 6｜全部单变量和多变量Cox结果（A–H为`PASS`）

**A–D，** 官方volcano构造器依次显示GSE41271 OS、GSE41271 RFS、TCGA OS和TCGA DFS的31个HC-TF单变量Cox结果。

**E–H，** 以相同顺序显示多变量Cox结果。官方脚本按名义p<0.05着色和标记；BH FDR保留在冻结统计表中，不能把所有带标签点称为多重校正显著。此前自绘的I–P置换分布不是TNBC capsule的Supplementary Figure 6面板，已从严格官方图集合撤回；5,000次置换结果可作为审计表保留，但不得冒充官方补图。

### Supplementary Figure 7｜扩展药物基因组学证据（A–C为`PASS`；D为`UNSUPPORTED`）

**A，** 官方Euler构造器显示PRISM、GDSC2和CTRPv2之间可比较药物的重叠。

**B，** 官方AAC散点构造器展示两组确定性选择的跨资源重复示例：vorinostat–ZNF254和dactolisib–ZNF254，并分别报告贡献数据集中的细胞系活性与药物反应关系。

**C，** 官方扩展低CI药物－TF点阵按补充通路类别显示重复关联，并附OS/RFS预后分类；高CI补充分支为`SUPPORTED_ZERO`，不补画。

**D，** `UNSUPPORTED`。没有通过细胞系重复与方向匹配后还能进入PDX验证的药物－TF组合，官方PDX AAC示例没有合格输入。

旧版“扩展补充图8/9”没有TNBC capsule对应面板，已经从严格官方图集合撤回；其模型分类器和PDX覆盖信息仅可作为审计表或方法附件保存。

## 五、Supplementary Tables与机器凭据

以下表格构成图件的统计依据，但“存在统计输入表”不等于相应面板已完成官方端口：

| 证据包 | 核心内容 | 当前用途 |
|---|---|---|
| Figure 1发现与验证表 | 351个显著发现TF（其中158个LUAD特异）；跨系统可估计350；GSE41271 70/158；GSE81089 95/158；共同44 | Figure 1及Supplementary Figure 2的输入和队列角色审计；ZNF737缺失边界另记 |
| Figure 2 HC-TF/motif表 | 31个HC-TF、20个motif可检验TF、正式严格三系统TF为FOXA3/NFATC4/XBP1；candidate版交集0仅作版本敏感性边界 | Figure 2及Supplementary Figure 4B–C已通过的官方端口输入与凭据 |
| Figure 3网络表 | regulon类别、3B边表、3D 31节点/134边拓扑、偏度和相关性 | 已完成Figure 3及Supplementary Figure 4D、5A–B的复核 |
| Figure 4 Cox/KM表 | 全部31个TF的单/多变量模型、BH FDR、冻结KM切点 | 已完成Figure 4和Supplementary Figure 6；置换结果只作表格审计 |
| Figure 5药物/PDX表 | 三资源全量CI检验、跨资源重复关系、423个PDX检验、0个合格验证输入 | 已完成Figure 5A–F和Supplementary Figure 7A–C；支持G/H及7D的`UNSUPPORTED`判断 |
| 官方端口凭据 | 官方源码SHA-256、unified diff、运行manifest、输出SHA-256、验证receipt | 证明只做了声明允许的数据/标签/路径适配 |

最终表格包计划随PDF写入`pilots/tnbc-chromatin-tf-nc-2026/execution/luad-official-code-port/final/`；在实际打包及校验完成前，本报告不把该目录描述为已交付。

## 六、TNBC与LUAD逐主图子面板比较

“必须改进”分为两类：数据或解释问题应立即修正；纯版式问题若需要改变官方构造器，只能在严格端口归档后建立单独、显式标记的publication-layout分支，不能悄悄改动并仍称“官方代码复现”。

| 子图 | TNBC中的生物学作用 | LUAD中的对应解释 | 解释力、可视化差异与改进边界 |
|---|---|---|---|
| 1A | 定义患者发现、METABRIC验证、PDX和细胞系分支。 | 拟定义TCGA发现、GSE41271主要验证、GSE81089次级验证及两类模型投射。 | `MANUAL_REQUIRED`；当前没有可声称完成的LUAD官方面板。以后人工绘制必须单列来源，不能再用旧自绘版。 |
| 1B | 说明TF活性与TF表达不是同一信号。 | 在LUAD-vs-LUSC中检验同一命题；发现层保留351个显著TF。 | 官方散点构造相同，代码端口已`PASS`；差异来自数据范围和标签密度。 |
| 1C | 以乳腺受体、PAM50、SCMOD2等注释展示患者异质性。 | 以肺癌可获得的临床、组织学和分子字段解释LUAD/LUSC活性。 | 注释功能对应但变量不等价；官方色标和聚类不改。必须避免把肺癌注释写成乳腺亚型的直接替代。 |
| 1D | 检验TNBC程序在PDX中的迁移。 | 检验包含158个LUAD特异TF的350-TF跨系统可用集合在PDMR PDX中的投射。 | LUAD PDX的诊断粒度和规模不同；必须公开模型筛选，不用热图聚类替代可迁移性的统计检验。 |
| 1E | 检验TNBC程序在乳腺细胞系中的迁移。 | 检验LUAD程序在DepMap肺癌细胞系中的投射。 | 官方热图语法相同；细胞系状态标签需如实报告。 |
| 2A | 引入患者、PDX、细胞系三系统ATAC。 | 拟说明22/13/19个LUAD模型。 | `MANUAL_REQUIRED`；不再把自绘流程图当成正式结果。 |
| 2B | 显示三系统启动子开放TF的集合交集。 | 对158个TF执行每系统至少半数样本开放规则。 | 使用官方ComplexUpset，代码端口已`PASS`；患者闸门分母为22，不与Figure 2C的21个RNA/VIPER配对样本混淆。 |
| 2C | 同时显示HC状态、平均活性、逐样本启动子和motif证据。 | 以GSE41271作为外部患者槽位显示158个TF及31个HC-TF的形成；RNA/VIPER可配对ATAC患者为21/22。 | capsule输出与印刷版不一致；严格端口保留四块纵向输出并标记`CAPSULE_PRINT_MISMATCH`，不自行仿画印刷版；Figure 2B/启动子闸门仍用22例。 |
| 2D | 显示各系统motif支持集合及交集。 | 显示LUAD仅3个严格三系统共同TF。 | 官方UpSet语法相同，但LUAD交集更小；必须报告精确成员和分母，不能靠放大图形制造“广泛程序”。 |
| 2E | 展示代表TF的系统别motif富集。 | 展示预先选定LUAD TF的严格和敏感性证据。 | 官方点/分布构造不改；主规则与敏感性规则必须视觉和文字分层。 |
| 3A | 总结HC-TF regulon方向、共享性和motif注释。 | 总结31个HC-TF的相同证据层。 | 已严格运行；LUAD规模较小。需避免把motif方格误写为Figure 3的再次入选门槛。 |
| 3B | 识别共同激活TF连接的靶基因群落。 | 得到3条边、6个端点和3个双节点组分。 | 已严格运行；稀疏不是绘图失败。旧版添加24个孤点和自绘dyad均撤回，不降低阈值美化网络。 |
| 3C | 展示HC-TF共同活性结构及极端相关TF。 | 展示31个TF的相关热图及top/bottom摘要。 | 官方输出包含两个原子图，最终拼版必须保留二者；不只截取热图。 |
| 3D | 展示共享激活靶基因所定义的TF协作网络。 | 31节点、134边、单一连通分量。 | 旧版80%裁边造成的假孤点已纠正；官方网络布局和degree映射保持。 |
| 3E | 检验regulon宽度与协作数量的关系。 | 对同一完整LUAD网络进行对应检验。 | 官方散点构造保持；不为可视化增加无依据的离群点标签。 |
| 3F | 识别仅在部分患者中活跃的偏态TF。 | 识别LUAD患者间TF活性偏态。 | 必须保留BH校正并限制为异质性线索，不能直接命名亚型。 |
| 3G | 展示代表性TF活性分布。 | 展示degree最高的CRY2、MAGED2、ZNF19、NR3C2。 | 选择规则确定且可审计；官方密度/分面构造不改。 |
| 4A | 外部队列OS多变量关联。 | GSE41271 OS中3个名义显著TF的官方森林图。 | 外部队列远小于METABRIC；只称探索性，完整31-TF分母见补图/表。 |
| 4B | 外部队列复发结局多变量关联。 | GSE41271 RFS中4个名义显著TF。 | 事件数有限；不能因官方图只显示显著行而隐藏筛查分母。 |
| 4C | TCGA OS多变量关联。 | TCGA-LUAD OS中10个名义显著TF。 | 仅ZNF75D通过BH；森林图颜色本身不代表FDR通过。 |
| 4D | TCGA复发结局多变量关联。 | TCGA-LUAD DFS中2个名义显著TF。 | 结局定义和可用样本不同于外部RFS，不能直接合并效应。 |
| 4E | 外部队列OS/RFS风险和保护集合重叠。 | GSE41271名义显著TF的方向别Euler图。 | 官方Euler图保留；小集合必须同时报告成员和完整分母。 |
| 4F | 展示外部队列代表TF的生存曲线。 | ETV1、ZNF444的OS/RFS曲线。 | 使用maxstat显示切点，和主要中位数筛查不同；只能描述，不能称独立复现。 |
| 4G | TCGA OS/RFS风险和保护集合重叠。 | TCGA OS/DFS方向别Euler图。 | 同4E；OS与DFS的重叠不是跨队列验证。 |
| 4H | 展示TCGA代表TF的生存曲线。 | CREBRF、PHC2的OS/DFS曲线。 | 必须保留风险表和删失信息；不能把同一队列多终点称为独立验证。 |
| 5A | 在一个药物资源中筛查TF－药物关系。 | PRISM的23,591项CI检验。 | 资源顺序必须按官方端口写为PRISM，而非旧报告中的GDSC2。 |
| 5B | 第二个药物资源筛查。 | GDSC2的3,627项CI检验。 | 覆盖小于PRISM，显著点数不可直接比较。 |
| 5C | 第三个药物资源筛查。 | CTRPv2的11,867项CI检验。 | 保留资源内BH，不跨资源混算FDR。 |
| 5D | 显示显著药物－TF关系的跨资源交集。 | 识别至少两个资源重复的规范化药物－TF组合。 | 已使用官方UpSet；旧版集合条形图撤回。方向一致性需由表格核实。 |
| 5E | 将TF的OS与RFS预后类别交叉汇总。 | 显示31个LUAD HC-TF的3×3预后分类计数。 | 旧报告误写为5组药物效应热图，现已纠正。官方tile构造保持。 |
| 5F | 将重复药物关联、通路与TF预后类别结合。 | 显示低CI主通路药物－TF关系；高CI分支为0。 | 只呈现有输入的官方分支；不得用空矩阵或自绘气泡图补足对称性。 |
| 5G | 在PDX中检验细胞系发现。 | 无合格输入，`UNSUPPORTED`。 | 不绘制漏斗或替代统计图；在拼版中明确显示缺失原因。 |
| 5H | 展示PDX疗效排序。 | 无已验证组合，`UNSUPPORTED`。 | 不能把“没有可进入检验的组合”写成“检验后阴性”。 |

## 七、输出位置与使用规则

约定的最终输出目录为：

`pilots/tnbc-chromatin-tf-nc-2026/execution/luad-official-code-port/final/`

最终7个评审PDF为`Figure1.pdf`、`Figure2.pdf`、`Figure3.pdf`、`Figure4.pdf`、`Figure5.pdf`、`Supplementary_Figures.pdf`和`TNBC_vs_LUAD_comparison_atlas.pdf`，固定输出到`execution/luad-official-code-port/final/`；页数、字节数、SHA-256和边界清单见[FINAL_REVIEW_SET.md](../execution/luad-official-code-port/FINAL_REVIEW_SET.md)。Figure 1–2官方代码端口已经`PASS`，但固定路径与机器校验通过不改变证据等级：整体仍只能称`STRICT_OFFICIAL_PORT_REVIEW_ONLY`、`publication_ready=false`，不得称最终投稿交付。

目前已核验的原子结果位于：

- Figure 3：`execution/luad-official-code-port/figure3/`
- Figure 4–5 P0修订版：`execution/luad-official-code-port/figure45-v2-p0fix/`
- Supplementary Figure 3：`execution/luad-official-code-port/supp3/`

旧目录`execution/luad-publication-rebuild/results/publication/`及其自绘PDF只保留历史追溯用途，不得再作为正式交付路径。
