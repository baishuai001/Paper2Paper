# LUAD Figure 2 与同步补充图最终报告

**版本：** 2026-08-19  
**状态：** Figure 2 方法学入口通过；“广泛三系统 motif 保守”不成立，当前结果支持“选择性跨系统 motif 核心”。  
**替代关系：** 本报告以完成后的 B0 主分析为准，取代 2026-08-14 形成的早期 `luad-figure2-gate-report.md` 结论。

## 1. 冻结分析规则

- Figure 1 冻结输入：158 个 LUAD-specific TF。
- ATAC 系统与固定分母：22 例原发 LUAD、13 个 LUAD PDX、19 条严格 LUAD 细胞系。
- 启动子开放：每个系统至少一半样本，即患者至少 11/22、PDX 至少 7/13、细胞系至少 10/19。
- HC-TF 活性排除规则：仅排除在患者、PDX、细胞系三个系统平均 NES 均小于 0 的 TF。
- 主 motif 分析：所有样本峰固定为 200 bp；不把各样本压缩至相同峰数；使用 B0 共同肺癌可及性背景；JASPAR 优先、CIS-BP 回退；HOMER `q <= 1e-5`；每系统至少 50% 固定样本支持。
- B1--B3 背景、LUAD 内在状态、峰数匹配和宽松阈值均只作为补充敏感性分析，不替换 B0。

## 2. 主分析结果

筛选链为：

`158 Figure-1 TF -> 141 患者启动子开放 -> 32 三系统启动子开放 -> 排除 1 个三系统平均 NES 均为负的 TF -> 31 HC-TF -> 20 个具有可测试 motif -> 3 个严格三系统 motif TF。`

严格三系统 motif TF 为：

- **FOXA3：** 患者 22/22、PDX 10/13、细胞系 18/19。
- **NFATC4：** 患者 22/22、PDX 10/13、细胞系 19/19。
- **XBP1：** 患者 22/22、PDX 10/13、细胞系 18/19。

31 个 HC-TF 才是与 TNBC 原文同构的后续 Figure 3/4/5 输入。3 个三系统 motif TF 是其中获得最强远端顺式元件证据的子集，不是再次定义 HC-TF 的过滤条件。

## 3. 用户指定的宽松阈值敏感性分析

仅把 HOMER 阈值放宽至 `q <= 0.05`，仍要求每个系统至少 50% 固定样本支持；未增加 log2 odds ratio 阈值，也未改分母。结果为 **5 个**，而不是此前 30% 样本情景下讨论的 10 个：

| TF | 原发 LUAD | LUAD PDX | LUAD 细胞系 |
|---|---:|---:|---:|
| ETV1 | 22/22 | 10/13 | 16/19 |
| FOXA3 | 22/22 | 10/13 | 19/19 |
| NFATC4 | 22/22 | 10/13 | 19/19 |
| XBP1 | 22/22 | 10/13 | 19/19 |
| ZNF75D | 22/22 | 7/13 | 10/19 |

这个 5-TF 结果只进入 Supplementary Figure 5D。主文 Figure 2 仍报告严格 B0 的 3 个 TF。

## 4. 背景与异质性敏感性结果

### 4.1 共同背景

| 背景 | 定义 | 峰数 | 严格三系统 motif TF |
|---|---|---:|---|
| B0 | TCGA-LUAD + TCGA-LUSC + 三个 LUAD 目标系统 | 292,983 | FOXA3; NFATC4; XBP1 |
| B1 | TCGA-LUAD + 三个 LUAD 目标系统 | 256,914 | FOXA3; NFATC4; XBP1 |
| B2 | TCGA-LUAD only | 138,103 | FOXA3; NFATC4; XBP1 |
| B3 | 三个 LUAD 目标系统 only | 189,542 | E2F5; FOXA3; NFATC4; XBP1; ZNF281 |

FOXA3、NFATC4、XBP1 在四种共同背景中均稳定。B3 新增的 E2F5、ZNF281 不替换 B0 主结果。

### 4.2 LUAD 内在状态

- 原发 LUAD：TRU 6、PP 8、PI 7、未分类 1。
- PDX：TRU 5、PP 3、PI 5。
- 细胞系：TRU 10、PP 5、PI 4。

TRU 限制性分析未得到三系统共同 motif TF，因此它只用于解释全 LUAD 队列的状态异质性，不能作为一个已通过的新主分析替代方案。

## 5. 图件与图注

### Figure 2

**位置：** `execution/luad-figure2/results/figures/Figure2_complete.pdf`

**Legend.** **(A)** Figure 2 的冻结设计：22 例原发 LUAD、13 个 PDX 和 19 条严格 LUAD 细胞系的 ATAC 数据用于检验 Figure 1 冻结的 158 个 LUAD-specific TF。**(B)** 158 个输入 TF 在三系统中的启动子可及性组合。**(C)** Figure 1 各表达系统中的 TF 活性类别，与三类 ATAC 系统逐样本启动子开放状态并列展示；31 个 HC-TF 单独标注。**(D)** 20 个可测试 HC-TF 中，15 个在至少一个系统达到严格 motif 支持阈值的组合；颜色表示 JASPAR/CIS-BP 分配类别。**(E)** B0 共同背景下各 TF 的平均 log2 odds ratio、系统间变异和逐样本分布。主 motif 判定为 `q <= 1e-5` 且每系统至少 50% 固定样本支持。

### Supplementary Figure 3

**位置：** `execution/luad-figure2/results/figures/SupplementaryFigure3_complete.pdf`

**Legend.** **(A--C)** 原发 LUAD、PDX 和细胞系的逐样本峰数。原发肿瘤使用公开固定峰矩阵；PDX 和细胞系公开数据无可重跑的真实生物学重复，故 IDR 记为不适用。**(D--F)** 1,000 次样本顺序置换得到的累计 consensus-peak 饱和曲线，并以 `nls(SSasymp)` 估计渐近峰数和饱和比例。**(G--I)** 三系统样本间全基因组可及性 Pearson 相关矩阵。**(J--L)** 按 promoter、exonic、intronic、distal 的层级规则进行逐样本峰注释。PDX MGH9243 的注释表总峰数为 0，四个类别比例不可定义；该样本仍保留在冻结的 13 个 PDX 主分析分母中，仅从比例箱线图和对应比例检验中排除。

### Supplementary Figure 4A--C

**位置：** `execution/luad-figure2/results/figures/SupplementaryFigure4A-C_complete.pdf`

**Legend.** **(A)** 从 158 个冻结 TF 到 141 个患者启动子开放 TF、32 个三系统启动子开放 TF、31 个 HC-TF、20 个 motif 可测试 HC-TF和 3 个严格三系统 motif TF 的完整证据流。**(B)** 三系统启动子开放组合及唯一一个“三系统平均 NES 均小于 0”排除项。**(C)** 31 个 HC-TF 的 motif 数据库覆盖：12 个同时具有 JASPAR 与 CIS-BP motif，8 个仅有 CIS-BP motif，11 个没有可分配 motif；没有 JASPAR-only TF。

### Supplementary Figure 5

**位置：** `execution/luad-figure2/results/figures/SupplementaryFigure5_complete.pdf`

**Legend.** **(A)** 三个系统中 TRU、PP、PI 和未分类 LUAD 状态的构成。**(B)** 20 个 motif 可测试 HC-TF 在 B0--B3 共同背景中的严格三系统支持稳定性。**(C)** 四种背景的峰数及严格三系统 motif TF 数；B0 为主分析。**(D)** 主阈值 `q <= 1e-5, prevalence >=50%` 与用户指定敏感性阈值 `q <= 0.05, prevalence >=50%` 的三系统 TF 对比；宽松结果不替换 B0。

## 6. 判定与下一步

- **可以继续 Figure 3/4：** 31 个 HC-TF 的定义已按 TNBC 同构规则完成。
- **不能写成“广泛三系统 motif 保守”：** 严格主分析只有 3/20 个可测试 HC-TF获得三系统支持；合理表述是“31-TF 染色质优先调控网络具有一个选择性、背景稳健的三系统 motif 核心”。
- **下一闸门：** Figure 3 应使用全部 31 个 HC-TF构建 regulon 结构、模块、协同与异质性；Figure 4 再检验 31 个 HC-TF 的患者结局关联和独立队列复现。

