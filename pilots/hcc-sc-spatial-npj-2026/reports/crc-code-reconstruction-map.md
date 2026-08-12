# NPJ HCC 换 CRC：逐模块代码来源与补建方案

审计日期：2026-08-09。

本文件回答的不是“这些分析在理论上能不能写代码”，而是：锚点作者代码是否可得、CRC 作者代码能复用到哪里、哪些必须从官方实现补建、同类论文代码只能提供哪些处理细节，以及每个模块怎样才算可靠。

主文逐子图对应关系见 [`../anchor/main-figure-panel-map.tsv`](../anchor/main-figure-panel-map.tsv)，补图逐子图对应关系见 [`../anchor/supplementary-figure-panel-map.tsv`](../anchor/supplementary-figure-panel-map.tsv)；数据层解释分别见 [`crc-subfigure-replacement-audit.md`](crc-subfigure-replacement-audit.md) 和 [`crc-supplementary-panel-audit.md`](crc-supplementary-panel-audit.md)。

QC、批次校正、样本覆盖、marker/注释、恶性 CNA，以及 cNMF 的 K、稳定性、患者混合与
程序合并不能只用下表一个“恶性识别与 cNMF”行概括。两篇论文和代码的逐项审计见
[`phase-zero-foundation-audit.md`](phase-zero-foundation-audit.md)，目标实现与验收见
[`../analysis/phase-zero-execution-specification.md`](../analysis/phase-zero-execution-specification.md)。

## 总结判断

1. **锚点作者代码不可公开取得。** 锚点原文只写“custom analysis code ... available from the corresponding author upon reasonable request”。因此不能精确核对作者的预处理、cNMF、空间分析和 Geneformer 细节，只能申请代码并同时独立重建。
2. **CRC-atlas 是首要作者代码来源，但不是开箱即用的整篇流水线。** 已固定审计提交 `82c15ecc2ea36baa0c4cffe8818bf58e21dcfc48`。仓库包含图谱、髓系 cNMF、Xenium、IMC 和 TCGA 脚本，根目录许可证为 BSD-3-Clause-Clear；但存在硬编码作者路径、缺失 README 所称的 downstream workflow/config、notebook 隐含中间文件和无自动测试等问题。
3. **GSE225857 的同癌种作者代码只能作处理参考。** [`jalon9358/LianLab_CRCLM`](https://github.com/jalon9358/LianLab_CRCLM) 提供单细胞和 ST 脚本，能帮助理解文件、患者和原发/肝转移处理；但当前仓库页面未见许可证，脚本含作者本地路径且部分依赖隐含 helper。未获得许可前不得复制进目标代码，只能阅读其分析顺序并用有许可证的官方包重写。
4. **关键算法都有官方底座，但“官方包”不等于目标论文代码。** 官方 cNMF、pySCENIC、CellChat、cell2location、COMMOT、MISTy、decoupler/PROGENy、Geneformer 等只能提供算法。患者切分、队列留出、manifest、统计单位、非交互入口、source table、失败关闭和预期结果仍要由目标项目补建。
5. **当前没有任何模块做过目标服务器真实 smoke。** 下表的“可改造/可补建”是代码审计结论，不是运行合格结论。

## 代码来源分级

| 级别 | 可做什么 | 不可做什么 |
|---|---|---|
| 锚点作者代码 | 若作者后续提供，可用于核对锚点实现和参数 | 当前不存在公开仓库，不能假装已经取得或验证 |
| CRC-atlas 作者代码 | 在固定提交、许可证归属、参数化路径和补测试后改造成 CRC 模块 | 不能把仓库存在等同于整条目标论文可运行；不能直接复用已发表结果充当新结果 |
| 官方方法实现 | 作为新代码的算法底座，固定版本/commit 和环境 | 不提供目标数据设计、患者独立性、Figure 叙事和结果验收 |
| 同类论文作者代码 | 帮助理解同癌种文件结构、合理处理顺序和常见参数 | 无明确许可证或不可移植时不复制；其结果和手工标签不能成为目标真值 |
| 目标项目代码 | 连接数据、算法、统计、图和验收，是最终可复现入口 | 不能用 AI 生成后未经 fixture、负向测试和真实小数据运行就称“可靠” |

## 模块级补建地图

| 模块 | 对应子图 | 作者代码可用部分 | 官方底座 | 同类论文可借细节 | 最终处理 | 晋级前必须验证 |
|---|---|---|---|---|---|---|
| 数据读取、manifest、注释和图谱 | F1A–F1D | CRC-atlas `workflows/build_atlas.nf`、`analyses/01_dataloader`、02–04 系列和 `src` helpers | Scanpy/AnnData 或 Seurat | GSE225857 `data_import_and_filter.R` 仅参考文件映射 | 参数化路径；优先读取处理后对象，不为“模仿”盲目重建 427 万细胞 | H5AD backed 读取；患者/样本/队列/组织 ID 唯一；计数、标签和作者嵌入版本冻结 |
| 组成与相关统计 | F1E–F1F、F4B/D、F5C/I | CRC-atlas 绘图/汇总代码可借输入格式 | scCODA 或合适组成模型；R mixed model/bootstrapping | 无需复制同类论文统计脚本 | 目标项目重写患者级统计和 source table | 细胞不作独立重复；组别人数、队列混杂、配对关系、留一队列结果均输出 |
| 恶性识别与 cNMF | F2A–F2F | CRC-atlas `analyses/20_myeloid_cNMF/01_preprocessing.ipynb`、`02_cNMF.ipynb`、`03_ORA_GSEA.ipynb` 只补流程细节 | [cNMF 官方实现](https://github.com/dylkot/cNMF)；恶性识别可用 CopyKAT/SCEVAN 官方实现 | GSE225857 `tumorcell_analysis.R` 参考恶性细胞处理 | 重新写恶性上皮 adapter、K/seed runner、程序合并和供者/队列留出 | 输入必须是可追溯非负计数；多种子一致；程序跨供者/队列；负向 fixture 能识别批次程序 |
| 功能、泛癌和 CRC 分型 | F2F、F2I–F2K、F4E–N、F5G/L–N | CRC-atlas cNMF 富集 notebook 和 TCGA 脚本可借输出格式 | fgsea/clusterProfiler/decoupler；[CMSclassifier/crcsc](https://github.com/Sage-Bionetworks/crcsc)；公开 CRIS/iCMS 签名 | 同癌种论文只帮助确认分型/基因集来源 | 目标项目统一基因 ID、背景、数据库版本和随机基线 | 基因集版本冻结；多重检验；随机匹配基线；独立队列复现；完整富集表而非只留显著项 |
| 调控与状态坐标 | F2G–H、F3A–C、F5J–K | 锚点无公开代码；CRC-atlas无目标恶性轨迹模块 | pySCENIC、Monocle/Slingshot、GeneSwitches 官方实现 | GSE225857 细胞类型脚本可借常规输入组织 | 独立补建；轨迹输出统一降级为状态连续性 | motif/数据库和种子冻结；多根节点/多算法/供者 bootstrap；证明结果不由患者或队列单独驱动 |
| 临床、生存和肿瘤—正常 | F3D–H、F5D–F、F8D–F、F8N、F8O–P | CRC-atlas `analyses/41_TCGA_survival/CRC_immune_subgroups_TCGA.R` 可作代码 donor | TCGAbiolinks、GSVA/ssGSEA、survival、DESeq2 | 其他 CRC signature 论文只借 endpoint/协变量报告规范 | 参数化数据入口，重写患者 ID、endpoint、连续得分、多变量、PH 与外部验证 | 患者 ID 双向连接；cutoff 不使用结局优化；开发/验证用途冻结；PH 诊断、多重校正和外部队列结果 |
| 药物预测 | F3I–J | 锚点只报告方法，无代码；CRC-atlas无对应模块 | [scTherapy 官方实现](https://github.com/kris-nader/scTherapy) | 相似论文可提供候选比较设计，不能提供目标真值 | 写 adapter、模型版本收据、正负对照和与 DepMap/PRISM/GDSC 的比较 | 已知正负对照；输入域检查；不同方法稳定性；输出明确为计算预测而非疗效 |
| 髓系、内皮和成纤维亚群 | F4A–F5N | CRC-atlas myeloid 与图谱代码可改造 | Scanpy/Seurat、伪 bulk 差异和富集官方包 | GSE225857 `myeloid_cells_analysis.R`、`endothelial_cells_analysis.R`、`fibroblast_cells_analysis.R` 仅参考 | 复用数据读取/绘图，重写跨供者状态验证和目标关系分析 | 亚群跨患者、跨队列；不由单一数据集/治疗决定；已发表 CRC-atlas/GSE225857 结论与目标结果分开 |
| 单细胞通信 | F6A–F6E | CRC-atlas可能提供邻域/下游上下文，但无可直接替代目标模块 | CellChat 官方实现 | GSE225857 同癌种脚本可借输入准备与常见 LR 检查 | 按患者运行、配对汇总，建立发送/接收最低覆盖和全结果表 | 患者级重复；细胞数阈值；数据库/参数敏感性；空间复核；始终使用关联措辞 |
| GSE225857 空间读取和解卷积 | F7A–D | CRC-atlas Xenium代码不能解析该 Visium/ST 资源 | [cell2location 官方实现](https://github.com/BayraktarLab/cell2location)；也可将 SPOTlight 作为方法敏感性 | GSE225857 `ST-seq_data_analysis.R` 对文件和作者处理顺序有价值，但无许可证只参考 | 目标项目重写 GEO archive parser、坐标/图像校验、单细胞参考和全切片批处理 | 实际解包 607 MB archive；单细胞—空间患者关系；所有 4 个原发/2 个肝转移切片；source table 与图像坐标一致 |
| CRC-atlas Xenium/IMC 邻域 | F7C–D、F7N | `analyses/60_Xenium_spatial/00_prepare_data.ipynb`–`05_results.ipynb`、`analyses/70_IMC_spatial/IMC_CN_analysis.py` 和 `src/spatial_helpers` 可改造 | NicheCompass/LIANA 等官方依赖 | 无需用 GSE225857 替代其细胞分辨角色 | 只下载首轮约 0.37 GiB 矩阵/细胞/analysis 对象；为每个程序先出 380 基因覆盖收据 | 患者重叠披露；面板覆盖达到阈值；全合格切片运行；不得声称未测基因的完整程序验证 |
| 空间通信、关系和通路 | F7E–J、F7L–P | 锚点无公开代码；CRC-atlas Xenium notebook可借空间对象结构 | [COMMOT](https://github.com/zcang/COMMOT)、MISTy、decoupler/PROGENy 官方实现 | GSE225857 脚本可借可视化布局，不复制 | 建跨切片 runner、空间置换、距离敏感性、训练/留出或全切片 ROI 分析 | 预先规定通路/ROI；多切片、多患者；空间空模型；所有比较和多重校正输出 |
| Geneformer 虚拟扰动 | F8A | 锚点不公开 fine-tuning、split 和 perturbation 代码；CRC-atlas没有此模块 | [Geneformer 官方模型/代码](https://huggingface.co/ctheodoris/Geneformer)（Apache-2.0） | 最近 CRC Geneformer 论文只用于重复性和设计比较，不复制其问题 | 从零写 target adapter、tokenization receipt、供者隔离 split、多个 seed、label permutation 和简单基线 | 小型 GPU smoke；token 覆盖；供者不泄漏；置换失败关闭；多 seed/split 稳定；优于简单 DE/网络基线 |
| DepMap 依赖与表达 | F8B–C、F8G | 锚点无代码；CRC-atlas无该模块 | DepMap 官方发布文件/API | 无需相关论文脚本 | 写版本化下载、CRC 细胞系筛选、ID join、pan-essential/组织选择性和稳健相关 | 发布季度和 checksum；细胞系清单；ID 唯一；留一细胞系；非 CRC 基线 |
| 非交互执行和产物验收 | 所有子图 | CRC-atlas notebook/Nextflow结构可借，但现有 downstream 入口不完整 | pytest/unittest、Snakemake/Nextflow任选其一，按目标规模决定 | 不复制散乱 shell/R 工作目录 | 每个模块提供 CLI、环境锁、输入/输出合同、日志、source table、Figure 和小 fixture | 路径无作者硬编码；空/错ID/缺基因/错版本时失败；真实小数据 smoke；同一输入可稳定重跑 |

## 哪些代码可以直接改造

“直接改造”仍然要求固定提交、保留许可证/归属、移除硬编码、补入口和测试。当前最值得直接改造的是：

- CRC-atlas 的 H5AD/metadata 读取和图谱组织代码，用于 F1、F4A–C、F5A/H；
- CRC-atlas 髓系 cNMF 的预处理、输出组织和富集 notebook，用作 F2A–F、F4E 的处理细节，但 cNMF 对象必须改成恶性上皮；
- CRC-atlas Xenium/IMC 读取、细胞邻域和绘图代码，用于 F7C/D/N 的有限面板/蛋白空间补充；
- CRC-atlas TCGA 生存脚本，用作 F3D–G、F5D–F、F8O–P 的输入与绘图参考。

不能“直接改造”的地方包括：作者仓库中不存在的恶性上皮 cNMF 稳健性、GSE225857 空间 parser、患者级组成/通信统计、CRC 分型双投影、治疗 estimand、Geneformer 扰动和 DepMap 连接。这些必须在官方实现上新建目标项目模块。

## 哪些只能从同类论文补处理细节

GSE225857 作者仓库提供了 `data_import_and_filter.R`、`tumorcell_analysis.R`、`myeloid_cells_analysis.R`、`endothelial_cells_analysis.R`、`fibroblast_cells_analysis.R` 和 `ST-seq_data_analysis.R`。它能帮助回答：文件怎样组织、哪些患者/组织相配、作者如何准备 ST 对象、哪些分析在原发/肝转移中已做过。

它目前不能成为可靠代码底座，原因是：

- 仓库未见明确许可证；
- 脚本包含作者机器路径和隐含中间对象；
- 部分 helper/手工区域标签不能从仓库完整追溯；
- 没有环境锁、自动测试或非交互端到端入口；
- 其 F3/MCAM 成纤维和 CRCLM 空间结论已经发表，不能被目标论文重新认领。

正确做法是读懂其文件合同和分析顺序，然后使用有许可证的官方包在目标项目重写；若确需复制代码片段，先取得许可证或作者许可并保留归属。

## 可靠代码的最低形态

每个进入真实分析的模块至少需要：

1. 固定数据版本、文件 checksum、代码 commit 和环境锁；
2. 一个非交互命令入口，而不是依赖 notebook 当前内存；
3. 明确的输入列、患者/样本 ID、输出 source table 和 Figure 路径；
4. 小型 fixture：正常路径能产生预期结构；错 ID、负数 cNMF 输入、缺关键基因、错误版本等会失败关闭；
5. 一次目标服务器真实小数据 smoke，保存命令、日志、退出码和产物；
6. 对统计模块额外验证独立单位、配对、队列留出、多重检验和主张边界；
7. 对空间/Geneformer 模块额外验证空间空模型、供者隔离、随机基线和稳定性。

因此，答案不是“作者没代码，所以 AI 把 92 张图重新写一遍”。可行策略是：**能合法改造的 CRC 作者代码保留其成熟数据处理；算法使用官方实现；同类论文只补文件和处理细节；Paper2Paper 补齐患者级设计、可移植入口、测试、真实 smoke 与主张限制。**
