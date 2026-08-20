# LUAD Figure 2 原文式协议（v2.2，2026-08-19）

## 1. 分析问题与输入

本轮只回答：Figure 1 发现的 LUAD 特异 TF，能否像 TNBC 锚点一样，在
LUAD 原发肿瘤、LUAD PDX 和严格鉴定的 LUAD 细胞系中获得一致的启动子
开放和 DNA motif 支持。

Figure 2 的输入是 Figure 1 发现阶段的全部 **158 个 LUAD 特异 TF**：

- 父表：`execution/luad-figure1/results/tables/tcga_specific_TFs.tsv`；
- 父表 SHA256：`6f449276a97b8c66cbf93bd2a4ca6d51095fa426015dc9046a7c084431d5037b`；
- 机械选择条件：`discovery_category == LUAD`；
- 独立患者队列中的 97 个同方向复现 TF 仍作为 Figure 1 外部验证结果，
  但不再错误地替代原文 Figure 2 的发现集输入。

## 2. 与 TNBC 锚点对齐的筛选顺序

锚点实际由 150 个 TNBC 特异 TF 起步：

1. 患者中启动子在至少一半样本开放；
2. PDX 和细胞系中也分别在至少一半样本开放；
3. 仅剔除患者、PDX、细胞系三个活动队列的平均 NES **同时小于 0** 的 TF，
   得到 HC-TF；
4. 对 HC-TF 的已知 DNA motif 以 JASPAR 和 CIS-BP 为独立参考集运行 HOMER；
5. 单样本 HOMER Benjamini `q <= 1e-5`（与作者公开代码的实际边界一致），且至少一半样本显著，定义该系统支持；
6. 报告在一个、两个或三个系统中获得 motif 支持的 TF 数量。

TNBC 原文由 150 → 128（患者启动子）→ 103（三端启动子）→ 94 HC-TF；
其中 9 个是因为三个系统平均 NES 同时小于 0 被剔除。94 个 HC-TF 中
56 个有已知 motif，43 个至少在一个系统获得 motif 支持，31 个在三端均获
支持。这些是锚点的观测结果，不是 LUAD 必须达到的最低数目。

## 3. 样本清单

| 系统 | 数据 | 独立 n | 主分析定义 |
|---|---|---:|---|
| 患者 | GDC TCGA-LUAD ATAC | 22 | 原发实体瘤；每位患者一次；两个技术重复 |
| PDX | GSE269746 / PRJNA1123604 | 13 | `MGH*-adeno` LUAD PDX；排除 tSCLC/de novo SCLC |
| 细胞系 | DRA Dataset 2 | 19 | Cellosaurus/DepMap 支持的 LUAD、DMSO 基线、每条一次 |

预定义样本须全部完成到可审计的 peak/BAM 或患者开放矩阵后才运行主分析。
这是计算完整性要求，不是生物学信号门槛。不能把技术重复、非 LUAD、污染
模型或药物处理库补入样本数。

## 4. ATAC 处理与质量信息

- 所有下载和大规模计算在云服务器完成；
- Nextera adapter 由 cutadapt 去除；PDX 先去除 mouse reads，再比对 hg38；
- 去除 duplicate、chrM、chrY、未比对、secondary/supplementary、MAPQ<30
  和 ENCODE blacklist reads；
- 启动子与 motif 主分析限定常染色体和 chrX；alt/random/unplaced contig 峰保留在
  原始 QC 中，但不伪装成 canonical promoter/motif 支持；
- 细胞系单端 reads 做 Tn5 `+4/-5 bp` 校正，MACS2 参数为
  `--keep-dup all -B --shift -75 --extsize 150 --nomodel --SPMR -q 0.01`；
- PDX 保留真实 paired-end fragment，MACS2 使用 `-f BAMPE ... -q 0.01`；
- 只有真实重复存在时做 IDR 0.05；单一公共文库标为不适用；
- 峰数、有效 reads/fragment、FRiP、ataqv、统一 Tn5-TSS profile 和饱和度
  全量报告。

锚点的“可及/不可及区域 qPCR 至少 10 倍”是建库后、测序前的实验 QC。
公共 LUAD 原始数据不能事后补做该实验，故记为不可回溯。此前添加的“至少
100 万有效单位、至少 10,000 峰、任一系统失败比例不超过 20%”不是锚点规则，
本协议不以这些参考值剔除已完成样本或停止分析。

TCGA 原始/比对文件受控。主分析使用 GDC 开放固定峰计数：同一肿瘤两个技术
重复均 `CPM >= 1` 才将峰记为开放；同步报告 `CPM >= 0.5`、`CPM >= 2` 和
合并重复 `CPM >= 1` 的敏感性结果，不以敏感性结果择优替换主定义。

## 5. 启动子和 HC-TF

- 候选调控因子的启动子按 GENCODE hg38 v47 basic 中 exact-symbol transcript
  定位；全基因组峰注释背景仍使用 protein-coding transcript；
- PAN-GO 冻结输入中的非蛋白编码候选显式写入映射 receipt，不在启动子层暗中
  删除；其已知 DNA-binding motif 可检验性在 motif 目录层单独判定；
- 主定义按 Methods：TSS 上游 2.5 kb、下游 1 kb；
- 图例 `-1 kb/+100 bp` 作为必须报告的敏感性定义；
- 同一 TF 任一转录本启动子与样本峰相交即记为开放；
- 系统阈值是 `ceiling(n/2)`：患者 11/22、PDX 7/13、细胞系 10/19；
- 三端均开放后，仅在患者、PDX、细胞系平均 LUAD NES 都小于 0 时剔除；
  缺失活动值不伪装成负值，也不改变原文的“AND”条件。

## 6. HOMER motif

- JASPAR 2024 CORE vertebrates 与 CIS-BP 2.0 Homo sapiens 独立运行；
- 每个数据库使用其**完整参考 motif 集**运行 HOMER，然后再映射到 HC-TF；
  不先缩成 HC-TF motif 子库，以免改变多重检验分母；
- HOMER Benjamini `q <= 1e-5` 定义单样本显著；LOR 完整报告但不另加原文
  未声明的正值门槛；
- 系统支持仍为 `ceiling(n/2)`；三端 motif TF 是三个系统均获得支持者；
- Figure 2D 显示至少一个系统有支持的 motif 可检验 HC-TF 及其三端组合；
- Figure 2E 显示这些 TF 的平均 LOR、SD、半数样本状态和最高富集 TF 的
  逐样本分布。

肺癌背景冻结为 TCGA LUAD/LUSC 类型峰和全部 13 PDX/19 LUAD 细胞系峰的
合并、去重开放区域。JASPAR/CIS-BP 原始快照、转换文件、SHA 和 HOMER 版本
写入 receipt。

## 7. 同步交付

- Figure 2A–E；
- Supplementary Figure 3A–L：峰数/IDR、1000 次顺序置换饱和度、样本相关性、
  峰基因组注释；
- Supplementary Figure 4A：158 → 患者开放 → 三端开放 → HC → motif 的数量流；
- Supplementary Figure 4B：三端启动子组合，并突出“三个系统 NES 同时<0”的
  唯一活动剔除；
- Supplementary Figure 4C：HC-TF 的 JASPAR only、CIS-BP only、both、none。

所有图和表共用同一 22/13/19 manifest、同一 158-TF 输入和同一阈值配置。

## 8. 判读与停止规则

本轮不再设置 TNBC 原文不存在的生物学 PASS/FAIL 数量线。以下旧规则已撤销：

- 样本失败比例最多 20%；
- HC-TF 至少 10 个；
- 三端 motif TF 至少 3 个或至少占 10%；
- promoter/motif 匹配置换 `p<0.05` 才准继续。

唯一会阻止出图的是执行不完整或完整性错误，例如预定义样本未处理完、文件
损坏、manifest 不一致、基因/TF 无法机械映射或 HOMER 输出缺失。完整运行后
如实报告 HC-TF 和三端 motif TF 的实际数目；这些结果是否足以支持 Figure 3，
另基于 Figure 3 的原文输入需求讨论，不由本协议中的人为最小数量自动裁决。

此前 2026-08-14 的 `FAIL_DATA` 文件保留为被本 v2 协议取代的历史记录，不删除、
不冒充本次生物学结果。

## 9. 三端 motif 为 0 后的技术等价性审计

本节由 2026-08-17 的用户指令追加，角色是解释三端 motif 为 0 是否由技术
不可比造成；它不回写或放宽第 2、5、6 节的原文式生物学规则。

1. **TNBC 阳性对照**：下载作者 Code Ocean 中逐样本 HOMER
   `knownResults.txt`，先用作者未修改的 Figure2/06B 脚本重组 JASPAR 与
   CIS-BP 合并表，再用同一 `q <= 1e-5`、每系统至少半数样本规则重算三端
   交集；逐 TF 与作者发布的交集及 LOR 汇总比较。
2. **LUAD 谱系阳性对照**：预先固定 `NKX2-1、FOXA1、FOXA2、CEBPA、CEBPB、
   GATA6、HNF4A、ELF3`，不得按结果增删；分别报告 PDX 和细胞系能否检出。
3. **统一输入**：所有可比较样本把峰按 summit（无 summit 时用区间中点）
   重设为 200 bp；使用同一 hg38 FASTA、同一 ENCODE hg38 blacklist。motif
   该技术诊断使用冻结的 47-model 候选集合（20 个 motif 可检验 HC-TF 加 8 个
   预设谱系对照在 JASPAR 2024/CIS-BP 2.0 中的全部可映射模型），但严格照
   TNBC 作者代码把 JASPAR 21 个模型与 CIS-BP 26 个模型**分开运行和校正**，
   再按作者规则对每个 TF 优先采用 JASPAR，只有 JASPAR 无该 TF 时才采用
   CIS-BP；每个来源的同一文件在三个系统中复用。把 47 个模型合并后
   一次校正作为更严格的敏感性结果；另以两库合并的完整 1,944-model 快照
   运行最严格的多重检验敏感性分析。后两者不得替代分库主诊断。
4. **峰数与 GC 匹配**：以完整 19 条细胞系中最小的可用峰数定义可匹配下限；
   对达到下限的样本使用完全相同的 5% GC-bin 配额和确定性抽样，使每个样本
   的目标峰数、长度和 GC-bin 计数完全相同。峰数低于下限者标为不可检验，
   不伪装为 motif 阴性。
5. **统一背景**：由 TCGA-LUSC 开放区域构建一份固定 200 bp、同 blacklist、
   同 GC-bin 比例的背景 BED；三个系统和所有样本复用同一文件。HOMER 对每个
   目标集自动去除与目标重叠的背景并做低阶寡核苷酸归一化，因此结果表中的
   有效加权背景分母可随样本略变，必须从 HOMER 表头读取，不能用 BED 行数
   替代。
6. **双分母报告**：主表按原预定义样本数计算半数门槛（患者 11/22、PDX
   7/13、细胞系 10/19）；同时报告仅限峰数可匹配样本的技术敏感性分母。
   两者都不构成新增的项目停止线。

若 TNBC 阳性对照不能恢复作者结果，优先判为代码/解析问题；若 TNBC 能恢复、
但 LUAD 谱系阳性对照在 PDX 或细胞系也失败，则三端为 0 暂不能作生物学缺失
解释；只有阳性对照成立后，LUAD HC-TF motif 缺失才具备较强的生物学解释力。

## 10. motif 主分析与两类敏感性分析的最终层级

2026-08-17 进一步逐行核对作者 `Figure2/05B-HOMER.sh` 后，确认峰到 HOMER
这一层的原文实现是：每个样本使用自身全部查询峰、以 `-size 200` 统一宽度，
并对所有样本复用一个共同乳腺可及性 consensus background；作者脚本没有把
每个样本压缩到最小峰数，也没有用 non-TNBC 可及峰作为竞争背景。

因此，LUAD Figure 2 motif 的最终分析层级冻结如下：

1. **主分析**：22 个患者、13 个 PDX、19 个细胞系均保留全部可用峰；按 summit
   （无 summit 时按峰中点）统一为 200 bp；去除非 canonical、ENCODE blacklist
   和含非 A/C/G/T 的区间。背景由 TCGA-LUAD、TCGA-LUSC 类型峰及全部预定义
   LUAD 患者/PDX/细胞系峰构建共同肺癌可及性 universe，重叠区间合并为
   consensus CRE 后重新居中为 200 bp。所有样本复用该背景，GC 差异交由 HOMER
   对共同背景执行的寡核苷酸归一化处理。
2. **数据库**：每个样本分别用完整 JASPAR 2024 和完整 CIS-BP 2.0 运行 HOMER，
   两库独立计算 Benjamini q；完成后才映射 20 个 motif 可检验 HC-TF 和 8 个
   预定义肺腺谱系阳性对照。TF 有 JASPAR motif 时采用 JASPAR；只有缺少 JASPAR
   时才回退 CIS-BP，严格遵循作者 06B/07A/08 的实际选择逻辑。
3. **峰数匹配敏感性**：每样本 3,584 峰、共同 GC-bin 配额的结果只检验不同
   峰深度是否改变结论，不得替代主分析。
4. **LUSC 背景特异性分析**：共同 TCGA-LUSC 可及峰背景回答“LUAD 相对 LUSC
   是否仍富集”，属于更强的肺癌组织学特异性问题，不得冒充 TNBC 式主分析。

固定样本分母仍为患者 22、PDX 13、细胞系 19。若某个预定义样本在 MACS2
`q=0.01` 和 canonical/blacklist 过滤后没有峰，则保留在原分母中并显式标记
`NOT_TESTABLE_NO_CANONICAL_PEAKS`；不降低 q、不引入 random contig 峰、不用其他
样本替换，也不把“全部峰”曲解为必须人为生成非空峰集。

作者 capsule 只给出其 consensus background 文件名和使用代码，没有提供背景
构建脚本或该 BED 本体。因此上述 LUAD 背景是功能角色同构且完全可复现的实现，
不能宣称与 TNBC 背景具有逐区间或逐字节等价性。此限制必须写入回执和报告。

## 11. 主图、异质性与阈值敏感性的展示层级

Figure 2 主图只展示第 10 节冻结的 B0 主分析，单样本 HOMER 阈值仍为
`q <= 1e-5`，系统支持仍要求至少一半固定样本（患者 11/22、PDX 7/13、
细胞系 10/19）。B1–B3、峰数匹配、LUSC 竞争背景、TRU/PP/PI 状态构成和
TRU-like 诊断重计数均属于解释异质性或技术稳健性的补充分析，不得替换
Figure 2 的 B0 结果。

根据 2026-08-19 的用户指令，新增且只新增一个 motif 阈值敏感性情景：

- 单样本 HOMER `q <= 0.05`；
- 每个系统仍须至少一半固定样本支持，即患者 11/22、PDX 7/13、细胞系
  10/19；
- 仍使用 B0 共同肺癌可及性背景、完整 JASPAR/CIS-BP 分库校正以及
  JASPAR 优先/CIS-BP 回退结果；
- 未检验或缺失的样本/TF 组合保留在固定分母中并计为不支持；
- 不增加 LOR 正值等锚点未声明的筛选条件。

先前讨论过的 `q <= 0.05` 且仅要求 30% 样本不再作为选定敏感性分析，也不
进入图、结论或候选 TF 清单。新的 `q <= 0.05 + 50%` 结果进入补充图，并与
原文阈值结果并列显示；它只能说明结论对统计阈值的敏感程度，不能重命名为
主结果。

同步交付层级最终为：Figure 2A–E；Supplementary Figure 3A–L；
Supplementary Figure 4A–C；以及承载状态异质性、B0–B3 背景稳健性和上述
阈值敏感性的 Supplementary Figure 5。
