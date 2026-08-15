# LUAD Figure 2 原文式协议（v2，2026-08-15）

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
5. 单样本 adjusted `p < 1e-5`，且至少一半样本显著，定义该系统支持；
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

- GENCODE hg38 v47 basic protein-coding transcripts；
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
- HOMER adjusted `p < 1e-5` 定义单样本显著；LOR 完整报告但不另加原文
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
