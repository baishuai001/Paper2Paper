# LUAD Figure 1 最终报告

## 结论

按照 TNBC 锚点文献的 Figure 1 代码语法和证据链，以预先冻结的 **LUAD vs LUSC** 比较完成真实数据运行。Figure 1A–E 均已生成；预先定义的生物学闸门判定为 **PASS**。

这里的 PASS 只表示“LUAD 相对 LUSC 的组织学特异 TF 活性程序”完成患者发现、独立患者复现、PDX 投射和细胞系投射。它不等于 Figure 2 的染色质证据已经成立，也不代表已获得 LUAD 内部治疗亚型或治疗靶点。

## 数据与计算

| 层级 | 队列 | LUAD | LUSC | 用途 |
|---|---|---:|---:|---|
| 患者发现 | TCGA | 516 | 501 | ARACNe3、msVIPER、limma、逐样本 VIPER |
| 独立患者验证 | GSE81089 | 106 | 67 | 独立 ARACNe3 网络及同方向 TF 复现 |
| PDX | NCI PDMR | 22 | 34 | 冻结 TCGA regulon 投射 |
| 细胞系 | DepMap 22Q2 | 76 | 27 | 冻结 TCGA regulon 投射 |

全部原始数据下载和计算位于云服务器：
`/media/desk16/iy13202/projects/Paper2Paper/tmp/tnbc-chromatin-tf-nc-2026/luad-figure1/`。

TCGA ARACNe3 网络从 121,882,726 条候选关系，经 FDR 和 MaxEnt/DPI 剪枝得到 201,319 条边；按作者脚本的实际行为再删除首条数据边后，2046 个 regulon 进入 VIPER。GSE81089 独立网络得到 77,615 条边；按同一作者行为删除首条边后，2008 个 regulon 进入 VIPER。

## 核心结果

1. TCGA 发现 158 个 LUAD-specific TF 和 193 个 LUSC-specific TF。msVIPER NES 与表达 log2 fold change 强相关：Pearson `r=0.7255`，`p<2.2e-16`。
2. 独立患者队列 GSE81089 复现 97 个 LUAD-specific TF 和 153 个 LUSC-specific TF。代表性的 LUAD 复现 TF 包括 `NKX2-1`、`FOXA1`、`FOXA2`、`ELF3`、`HNF1B`、`GATA5`、`HOPX`、`SPDEF`、`TOX3` 和 `RORC`；LUSC 侧包括已知鳞状谱系轴 `TP63`、`SOX2`。
3. 由 97 个独立患者复现 LUAD TF 构成的程序，在三个外部系统中均显著升高：

| 队列 | 实际可计算 TF | Hedges' g | BH FDR |
|---|---:|---:|---:|
| GSE81089 | 97 | 3.7886 | 7.58e-46 |
| PDMR PDX | 97 | 3.2026 | 1.10e-15 |
| DepMap 22Q2 | 96 | 1.0772 | 9.71e-05 |

4. 不使用组织学标签切分、只按 TF 活性进行两簇层次聚类时，与真实组织学的一致率为 TCGA 93.8%、PDMR PDX 91.1%、DepMap 80.6%；对应 Fisher 检验 FDR 分别为 2.48e-207、2.76e-10 和 7.03e-05。
5. 预设停止规则要求：至少 3 个独立患者复现 LUAD TF，并且 PDMR、DepMap 均满足方向为正、BH FDR≤0.05、Hedges' g≥0.5。实际为 97 个复现 TF，且两个模型系统全部超过阈值，因此判定 **PASS**。

## 与 TNBC 锚点的对齐及偏差

- 保留：PAN-GO 调控因子、ARACNe3 默认单 subnet/FDR/MaxEnt-DPI、`rowTtest` signature、`ttestNull(per=1000, repos=TRUE, seed=1)`、`msviper(minsize=1)`、逐样本 `viper(minsize=1)`、limma-voom 以及作者代码额外删除首条网络边的实际行为。
- 性能适配：1,017 位 TCGA 患者的单核置换在 2% 时预计超过 12 小时；在尚未产生任何 TF 结果之前，冻结为使用 `ttestNull` 官方 `cores=32` 和 `L'Ecuyer-CMRG` 随机流。置换次数、是否有放回、种子和检验定义均未改变。
- 兼容补丁：`viper 1.38.0::aREA` 对两个单靶基因 regulon 出现 `1 x n` 矩阵降维错误。补丁仅恢复矩阵形状，没有删除 regulon，也未改变 `minsize=1` 或任何数值定义。首次错误日志和补丁均保留。
- GSE81089 建 regulon 时出现一次 mixture component 1000 次未收敛警告，随后另一 component 在 214 次收敛并成功输出网络。所有最终 msVIPER 数值均为有限值；该警告不隐藏，作为局限保留。

## 质量控制与交付

- 独立验证脚本 28 项检查全部通过：样本数、患者/供体去重、TF 交集、有限数值、BH 校正、停止规则重算、兼容补丁披露及 A–E 的 PDF/PNG 完整性。
- 图形经人工目检；仅修正图例/标签裁切和 `p=0` 的显示方式。重绘前后 10 个核心数值文件的解压后内容 SHA256 完全一致。
- 最终云端紧凑包 SHA256：`2d2788860a61a88951e9d7615edad588e7452426de00eec58074edc6e1c84eb1`。
- 原始大文件、完整网络和大型 RDS 保留在云服务器；Git 仅保存代码、协议、manifest、小型结果、图、日志和审计 receipt。

### 临床/分子注释与同步补图

- Figure 1C–E 已在不重算 TF 活性的前提下补充临床和分子注释。TCGA 的性别、年龄、分期、吸烟量、已发表 LUAD 表达亚型和关键驱动/抑癌通路突变覆盖率分别为 100%、97.25%、98.82%、76.20%、40.12% 和 100%；缺失值保留为缺失，不作推断填补。
- GSE81089 的性别、年龄、分期、吸烟和生存状态均为 100% 覆盖。PDMR PDX 的性别、年龄、吸烟史和组织类型均为 100%，取材部位为 98.21%；已知转移病史与公开分子记录分别仅覆盖 30.36% 和 23.21%，图中明确显示缺失。DepMap 的性别、年龄、原发/转移来源和取材部位覆盖率分别为 99.03%、80.58%、91.26% 和 100%。
- Supplementary Figure 1 已补齐四队列纳入/排除流程；Supplementary Figure 2A–F 已同步生成独立队列表达—TF 活性一致性、方向一致 TF 交集、三套样本相关性及冻结 97 TF 的跨系统效应。
- 注释与补图独立验证 25/25 项通过；注释重绘前后冻结数值文件 SHA256 差异为 0，因此原 Figure 1 的 97 TF、效应量和 PASS 判定未改变。

## 科学解释与下一步边界

Figure 1 得到的是一个很强的 **肺腺癌相对肺鳞癌的谱系/组织学 TF 程序**，不是 LUAD 内部患者分层。强信号包含 `NKX2-1/FOXA/HNF/ELF3` 与 `TP63/SOX2` 等已知腺癌—鳞癌谱系轴，因此结果既是方法迁移成功，也提示后续新颖性风险。

只有在另行冻结并通过 Figure 2 的患者原发肿瘤、PDX 和细胞系三端 ATAC/启动子开放/motif 证据后，才能继续称为“染色质支持的高置信 TF”。本报告到 Figure 1 为止，不外推预后、药敏或治疗结论。
