# LUAD严格官方代码端口：最终评审集

> 生成日期：2026-08-21
> TNBC官方来源：CodeOcean capsule `7227095/v1`，commit `edf5314`
> overall_status：`STRICT_OFFICIAL_PORT_REVIEW_ONLY`
> publication_ready：`false`

本评审集用于检查TNBC官方绘图代码迁移到LUAD冻结结果后的实际输出，不是最终投稿版。固定输出目录为：

`C:\Users\bai\OneDrive\Documents\Codex\2026_07_09 Literature_Reproduction\Paper2Paper\pilots\tnbc-chromatin-tf-nc-2026\execution\luad-official-code-port\final\`

## 1. 文件凭据

| 文件 | 页数 | Bytes | SHA-256 | 状态 |
|---|---:|---:|---|---|
| Figure1.pdf | 1 | 25,825,256 | `07509ca111afdfbe96997013a34f66e915b2ffe67b3f2d61ed81b95983bbd59fb` | REVIEW_ONLY |
| Figure2.pdf | 1 | 857,796 | `0ff961a115cfed312483b9860a96ac6b7b16ad5929fcd9266c023031373d37303` | REVIEW_ONLY |
| Figure3.pdf | 1 | 959,385 | `be5c57d274a3b8f264baa95757bccc571f9b23f5a763959441bcdbb5968377885` | REVIEW_ONLY |
| Figure4.pdf | 1 | 544,402 | `0e5a5f827ba501b88da79c65a9f059718e15aec95881d7b066287c1951bc85bb24` | REVIEW_ONLY |
| Figure5.pdf | 1 | 7,005,919 | `492e656cd9452cbdd8df1018285df928deb3ae26a5914c1e2d26587b86c6fb300` | REVIEW_ONLY |
| Supplementary_Figures.pdf | 8 | 74,249,841 | `33629b8dd2e1b196d34963acbe228141b58f2f796a61eb234b020e28feda99bd8` | REVIEW_ONLY |
| TNBC_vs_LUAD_comparison_atlas.pdf | 73 | 356,937,234 | `9bc31c1c52220bfe5723d5d48f922ca922c5a09d1d9bea288a181d7f9d337470` | REVIEW_ONLY |

总页数：86。

严格机器校验：`STRICT_PASS_7_PDFS_86_PAGES`；失败项：`0`。

最终视觉复核：`REVIEW_ONLY_QA_COMPLETE_WITH_DECLARED_OFFICIAL_LAYOUT_AND_REFERENCE_CROP_BOUNDARIES`。

## 2. Figure 1–5冻结结果

### Figure 1

TCGA-LUAD与TCGA-LUSC比较得到351个显著发现TF，其中158个LUAD特异、193个LUSC特异。Figure 1B保留351个；Figure 1C–E使用350个跨TCGA、PDMR PDX和DepMap可估计TF。ZNF737虽有35个regulon靶基因，但与DepMap表达矩阵交集为0，因此不零填充、不插补。GSE41271微阵列主要验证队列为183例LUAD/80例LUSC，复现70/158；GSE81089次级RNA-seq队列为108/67，复现95/158；两队列共同复现44/158。

### Figure 2

158个LUAD发现TF经原文式闸门得到31个HC-TF，其中20个有可检验已知motif。按冻结的JASPAR优先、无JASPAR才使用CIS-BP规则，严格患者、PDX、细胞系三系统共同motif TF为FOXA3、NFATC4和XBP1。放宽到q≤0.05且每个系统仍要求至少50%样本支持时，另出现ETV1和ZNF75D；二者只属于敏感性结果，不进入严格主集合。Figure 2C为`PASS + CAPSULE_PRINT_MISMATCH`，RNA/VIPER配对为21/22；Figure 2B及启动子闸门仍使用22例ATAC患者，不插补TCGA-44-A47F。

### Figure 3

Figure 3–5沿用全部31个HC-TF，而不是只使用3个三系统motif TF。Figure 3B为6个节点、3条边和3个双节点组分，不添加无边孤立节点。Figure 3D为31个节点、134条正边、1个连通分量和0个孤立点；节点大小映射degree，边宽映射共享激活靶基因数。

### Figure 4

四个主森林图中名义显著TF数依次为3、4、10和2。只有TCGA-LUAD总体生存中的ZNF75D通过队列内BH校正，其余只能表述为探索性预后关联。Kaplan–Meier曲线采用冻结的最大选择秩切点，仅用于描述性展示。

### Figure 5

PRISM、GDSC2和CTRPv2分别包含23,591、3,627和11,867项药物－TF检验。冻结结果包含5组至少在两个细胞系资源中同方向重复的关系：5-fluorouracil–ZNF254、5-fluorouracil–ZNF540、PHA-793887–WWC2、vorinostat–ZNF254和dactolisib–ZNF254。0/423个公共PDX检验进入匹配验证，因此Figure 5G、Figure 5H及Supplementary Figure 7D为`UNSUPPORTED`。

## 3. 显式缺口

### MANUAL_REQUIRED（5）

Figure 1A、Figure 2A、Supplementary Figure 1A、Supplementary Figure 1B和Supplementary Figure 4A没有公开官方绘图构造器；评审集不以旧自绘图冒充官方代码复现。

### UNSUPPORTED（3）

Figure 5G、Figure 5H和Supplementary Figure 7D没有满足冻结匹配规则的LUAD PDX验证输入。这表示数据覆盖不足，不是经过体内检验后的阴性结论。

### SUPPORTED_ZERO（2）

Figure 5F的高CI主通路分支和Supplementary Figure 7C的高CI补充通路分支经过正式计算，但满足条件的记录数为0，因此不生成替代图。

## 4. Supplementary Figure 3边界

- 3A为`PARTIAL`：使用公开固定GDC open-count矩阵，不是本项目重新逐样本调用MACS2/IDR得到的峰。
- 3B和3C缺少合格真实生物学重复，IDR为不适用。
- 3D–F使用冻结的1,000次置换峰饱和对象；3G–I使用冻结的共识可及性Pearson矩阵。
- 3J–L使用真实hg38 500-bp全基因组背景，共6,062,095个bin。
- 3K的13个冻结PDX中有1个零峰，官方非空长表显示n=12。
- 3D、3F、3K和3L存在官方固定标签坐标重叠；严格端口没有移动标签，因此整组仍为`publication_ready=false`。

## 5. TNBC参考图矩形裁切例外

以下边界只涉及`TNBC_vs_LUAD_comparison_atlas.pdf`左侧的TNBC整页参考图裁切，不涉及右侧LUAD官方原子PDF。部分标注跨越印刷版面板边界，无法仅用一个矩形同时做到“保留全部本面板内容”和“完全消除相邻面板碎片”。本次以保留科学内容为优先，并明确记录：

- Figure 3A：为保留横轴标题和Net Effect注释，保留极少量相邻Figure 3B文字。
- Figure 3B：为保留完整的GTPase标签，左下仍可能出现极小相邻字符。
- Figure 3D：继续收紧会裁掉RCOR2节点或网络连线；斜向`ARTICLE IN PRESS`是TNBC论文源图水印，不是重绘污染。
- Figure 5G：参考面板内容与Figure 5H共享边界，单一矩形无法完全分离；最终裁切优先保留Figure 5G的完整横轴标题，因此仍会看见Figure 5H的大写面板字母和少量边缘。
- Supplementary Figure 2C：为移除2A轴题并保留热图，原始C字母顶部被裁；图谱页眉仍标识面板。
- Supplementary Figure 3E：`Model Fit`标签伸入3F列，干净矩形会省略该越界标签的一部分。
- Supplementary Figure 3I与Supplementary Figure 4D：为去除相邻面板残片，原始面板字母被裁；图谱页眉仍标识面板，科学内容完整。
- Supplementary Figure 3K和3L：为保留本面板主体，边界仍可能保留相邻面板的极小线段或字符。
- Supplementary Figure 7A：集合图右缘与7B相撞；去除7B碎片会轻微裁切外轮廓，但集合标签和计数完整。
- Supplementary Figure 7B：为去除7A粉色弧线，原始B字母被裁；图谱页眉仍标识面板。
- Supplementary Figure 7C：底部机制注释靠近7D边界；继续放宽会引入7D。

这些例外不得表述为LUAD官方端口的数据缺失、裁边或重绘；它们只是TNBC整页参考栅格在逐面板比较图谱中的几何限制。

## 6. 官方构造器的固定排版边界

以下问题来自TNBC公开脚本中的固定坐标、固定画布或标签密度；严格端口没有偷偷改动`geom`、字号、坐标或避让规则。机器验证通过并不等于这些面板已达到投稿排版质量：

- Figure 1B和主要验证版Supplementary Figure 2A沿用官方`stat_cor(label.x=2.5, label.y=-3.75)`，在LUAD数值范围内右下角相关性文字会被绘图区边界截断；完整统计量保留在结果表和图例报告中。
- GSE81089次级Supplementary Figure 2A沿用官方标签算法；在351个发现TF的LUAD数据上，顶部和底部若干基因标签发生重叠。
- Figure 4B的长队列/终点标题在四列主图拼版中较拥挤；Supplementary Figure 4D顶部小热图标题与Pearson图例标题发生叠印。
- Supplementary Figure 3D、3F、3K和3L的固定标注重叠已在上一节单列。

这些均属于待明确标记的`post-port layout`工作；若后续修订，必须保留当前严格端口作为不可改写的审计基线，并逐项披露排版调整。

## 7. 最终判定

该评审集证明：有公开构造器且有合格LUAD输入的Figure 1–5及补充面板已经按锁定的TNBC官方代码执行和拼版。三系统motif证据支持的是20个可检验HC-TF中的3个选择性保守核心，不能表述为LUAD存在“广泛三系统共同motif程序”。评审集同时保留5项`MANUAL_REQUIRED`、3项`UNSUPPORTED`、2项`SUPPORTED_ZERO`、TNBC参考裁图例外及官方固定排版边界。因此最终状态固定为：

`STRICT_OFFICIAL_PORT_REVIEW_ONLY`

`publication_ready=false`

在补齐人工面板、独立体内验证输入和post-port排版修订之前，不得改称“最终投稿版”。
