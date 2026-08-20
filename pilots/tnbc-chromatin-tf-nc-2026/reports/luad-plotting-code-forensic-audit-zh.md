# LUAD图形实现取证审计（历史报告，已被严格官方代码端口取代）

> 原审计日期：2026-08-20
>
> 状态：**SUPERSEDED / WITHDRAWN**
>
> 取代版依据：TNBC CodeOcean capsule `7227095/v1`，commit `edf5314`
>
> 当前整体交付等级：`STRICT_OFFICIAL_PORT_REVIEW_ONLY`

## 如何使用本报告

本报告只保留“为什么旧自绘图被撤回”的取证历史。旧版中建议的自定义wrapper、自定义解释图、统一Liberation Sans、手工面板重排、自定义色板和`tnbc_anchor_theme.R`路线已全部**SUPERSEDED**，不得再执行，也不得作为当前图件的方法依据。

当前唯一允许的主路线是：对有公开构造器的面板，直接执行校验锁定的官方脚本或官方原始行块；只做路径、列名、LUAD/LUSC标签、数据集显示名和冻结样本数等声明过的数据适配。无构造器的面板记为`MANUAL_REQUIRED`，缺少合格LUAD输入的面板记为`UNSUPPORTED`。

## 旧自绘实现为何被撤回

旧实现将“数据格式适配”扩大成“重新设计图件”，因而同时改变了视觉语法和部分科学对象：

1. 没有逐段执行TNBC官方绘图构造器，而是使用通用主题、动态色板和自定义几何重画。
2. 将性质不同的热图、森林图、网络图和药敏图统一塞入报告式模板，造成二次缩小和信息层级丢失。
3. Figure 2B/2D把官方ComplexUpset改成自制UpSet；Figure 2E把motif case examples改成prevalence气泡图。
4. Figure 5E把官方OS×RFS的3×3预后分类计数图误画成药物效应热图。
5. Supplementary Figure 2B、2F和Supplementary Figure 5–7多处用通用柱状图、热图或气泡图替代官方Euler、逐样本NES分布、散点和通路点阵。

因此，旧`code/publication_rebuild/`渲染器、`anchor_style`拼图以及`execution/luad-publication-rebuild/results/publication/`中的PDF只保留历史追溯价值，不得进入严格端口拼版。

## Figure 3B/3D的两个历史性对象错误

### Figure 3B

TNBC官方代码只把满足“靶基因对共享至少3个HC-TF”的边及其端点交给`graph_from_data_frame()`。LUAD在同一规则下只有3条边、6个端点：

| 靶基因对 | 共同HC-TF |
|---|---|
| NDNF–ADGRF5 | CSRNP1、ETV1、NR3C2 |
| ADGRD1–SELENBP1 | CRY2、MXD4、ZBTB18 |
| YPEL3–CACFD1 | MAGED2、MXD4、ZNF444 |

旧自绘代码额外加入24个无边节点，制造了不属于官方统计对象的“散点涂鸦”。当前严格端口只保留6节点/3边，不再另画自定义dyad解释版。

### Figure 3D

LUAD完整TF协作网络为31个节点、134条正边、1个连通分量、0个孤立节点。旧自绘代码用未预注册的80%边权分位数删除101/134条边，并把节点大小误映射为regulon靶基因数。当前严格端口执行官方构造器：使用全134条正边，节点大小映射network degree，边宽映射共享激活靶基因数。

Figure 3当前为`PASS`，共19个官方脚本输出PDF。运行目录中的早期失败验证日志只是历史中间件，已被最终`PASS` receipt取代，不得与最终凭据并列解读。

## 当前严格官方代码端口的冻结事实

### Figure 1–2：`PASS`

- Figure 1B保留351个显著发现TF（158 LUAD＋193 LUSC）；Figure 1C–E使用350个跨系统可估计TF。ZNF737的35个regulon靶基因与DepMap表达矩阵交集为0，因此不零填充、不插补。
- GSE41271微阵列（183 LUAD/80 LUSC）是主要外部验证，复现70/158；GSE81089 RNA-seq（108/67）是隔离的次级验证，复现95/158，不取代主要微阵列队列。
- 31个HC-TF中20个具有可检验已知motif。冻结的作者式JASPAR优先版本得到3个正式三系统共同motif TF：FOXA3、NFATC4、XBP1。
- 早期candidate motif输入版本得到0个三系统交集；这是输入版本敏感性边界，不是正式Figure 2结果，不得用它覆盖FOXA3/NFATC4/XBP1。
- Figure 2C只有21/22例ATAC患者可与RNA/VIPER配对；TCGA-44-A47F无可用RNA/VIPER。Figure 2B和启动子闸门仍使用22例ATAC分母，不对第22例插补。
- Figure 2C公开capsule输出与论文印刷版布局不一致，状态为`CAPSULE_PRINT_MISMATCH`。严格端口保留公开代码输出；本报告旧版“以印刷版为视觉锚点自行重排”的建议已**SUPERSEDED**。

### Figure 4–5：P0修订版`PASS`

独立P0审计后，旧Figure 4–5 v1输出已**SUPERSEDED**。取代版完成了以下修正：

- 取消对官方`Helvetica`的`Liberation Sans`替换，并恢复Supplementary Figure 5B的FDR星号说明；字体替换数为0。
- Supplementary Figure 6D/6H图内终点显示从`RFS`更正为冻结的`DFS`语义，不改数据对象。
- Figure 4A–D标题恢复官方`Multivariate {cohort} {endpoint}`语义。
- Supplementary Figure 7B实际执行4段官方case调用和保存代码，保留7×7英寸官方画布。

机器验收为`PASS`、0项失败；实际执行28个官方源码块，生成39个原子PDF，全部具有逐文件SHA-256并由receipt锁定。Figure 5G、Figure 5H和Supplementary Figure 7D由于0/423个PDX检验能进入匹配验证，依真实输入记为`UNSUPPORTED`，不用自绘图补位。

### Supplementary Figure 3：代码端口通过，但`publication_ready: false`

Supplementary Figure 3A–L的12个PDF均由官方构造器生成，构造器重写数为0，receipt为`PASS_STRICT_PORT_WITH_DECLARED_INPUT_AND_VISUAL_BOUNDARIES`。同时必须保留下列边界：

- 3A为`PARTIAL`；3B/3C没有合格真实生物学重复，IDR为不适用。
- 3J–L使用真实hg38 500-bp全基因组背景，共6,062,095个bin。
- 3K的13个冻结PDX中有1个零峰，官方非空长表显示n=12。
- **3D、3F、3K、3L共4个面板存在官方固定标签重叠。** 严格端口不移动标签，因此整组明确`publication_ready: false`。

最终逐页视觉QA还确认：官方固定坐标使Figure 1B及主要验证版Supplementary Figure 2A的相关性文字在LUAD范围内被边界截断，GSE81089次级散点的部分基因标签重叠，Supplementary Figure 4D的小热图标题与图例标题叠印。它们是严格端口必须如实保留和披露的排版边界，不是数据缺失；若以后修复，只能作为单列的post-port layout版本，不能反写为“官方构造器零改动”。

## 当前交付边界

7个评审PDF为`Figure1.pdf`、`Figure2.pdf`、`Figure3.pdf`、`Figure4.pdf`、`Figure5.pdf`、`Supplementary_Figures.pdf`和`TNBC_vs_LUAD_comparison_atlas.pdf`，固定输出到`execution/luad-official-code-port/final/`；页数、字节数、SHA-256及最终边界见[FINAL_REVIEW_SET.md](../execution/luad-official-code-port/FINAL_REVIEW_SET.md)。

即使固定输出和机器验证已经完成，由于上述Supplementary Figure 3及数据覆盖边界，整体输出仍只能称：

> `STRICT_OFFICIAL_PORT_REVIEW_ONLY`

不得称`publication_ready`、“最终投稿版”或“零边界完整复刻”。
