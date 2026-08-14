# LUAD Figure 2 染色质闸门冻结协议

## 1. 闸门问题

本闸门只回答：Figure 1 冻结得到的 **97 个、在独立患者队列中同方向复现的 LUAD-vs-LUSC TF**，是否能像 TNBC 锚点一样，在 LUAD 原发肿瘤、LUAD PDX 和严格鉴定的 LUAD 细胞系三端获得一致的启动子开放和 DNA motif 证据。

Figure 1 数值结果不得重算或按 Figure 2 结果改选。唯一 TF 输入由 Figure 1 的合并复现表机械冻结：

- 父表 `externally_replicated_TFs.tsv` 同时含 LUAD 97 个和 LUSC 153 个复现 TF；父表 SHA256 为 `eace4fdf6509d01e8a9648b34b4dc3fdfb5955627c2c7a526c6ce1495a073844`；
- 固定选择条件为 `externally_replicated == TRUE AND discovery_category == LUAD`，不得读取 LUSC 方向的 153 个 TF；
- 排序并写为三列规范 TSV 后的预期 TF 数为 `97`，派生输入 SHA256 为 `64b4aa4e226c6e9142ecd560ec51ec1f49e0c9823d07e332757deeea3aa3c69b`。

本闸门不声称发现 LUAD 内部分型，也不进行预后、药敏或治疗推荐。闸门不通过即停止，不进入 Figure 3–5。

## 2. 与 TNBC 锚点严格对齐的逻辑

锚点 Figure 2 的核心顺序保持不变：

1. 在每一生物系统中生成独立样本的 ATAC 峰；
2. 以“至少一半样本开放”定义系统层面的 TF 编码基因启动子开放；
3. 取患者、PDX、细胞系三端启动子开放的交集；
4. 去除在 TCGA、PDX、细胞系中 LUAD 平均 VIPER NES 任一小于 0 的 TF，得到 LUAD HC-TF；
5. 对开放 CRE 分别使用 JASPAR 与 CIS-BP 独立运行 HOMER motif 富集，并在 TF 层合并两库支持；
6. `BH q < 1e-5` 定义单样本 motif 富集，至少一半样本富集定义系统层面 motif 支持；
7. 同一运行同步生成 Figure 2、Supplementary Figure 3、Supplementary Figure 4A–C。

锚点方法正文/可执行代码的启动子窗口是 TSS 上游 2.5 kb、下游 1 kb；图例写作上游 1 kb、下游 100 bp。主分析冻结为 `-2500/+1000 bp`，并把 `-1000/+100 bp` 作为必须报告的图例窗口敏感性分析，不事后选择更有利的窗口。

## 3. 三端样本与独立性

| 系统 | 主分析资源 | 预期独立 n | 主分析纳入规则 |
|---|---|---:|---|
| 患者 | GDC TCGA-LUAD ATAC | 22 | 原发实体瘤；作者 QC 通过；每位患者/组织块只计一次；两个 ATAC 反应为技术重复 |
| PDX | GSE269746 / PRJNA1123604 | 13 | 标题为 `MGH*-adeno` 的独立 LUAD PDX；转化性和原发 SCLC 均排除 |
| 细胞系 | DRA006908–DRA006930 Dataset 2 | 19 | Cellosaurus 疾病为 lung adenocarcinoma、无误鉴定警告、DMSO 基线；一条细胞系只计一次 |

细胞系作者面板共有 23 条；H1703 为 LUSC，H1299 和 H2126 为肺大细胞癌，RERF-LC-OK 已确认污染并实际为 Marcus/astrocytoma 衍生物。前两条大细胞癌仅进入宽松 NSCLC 敏感性分析；H1703 和 RERF-LC-OK 不得用于 LUAD 主分析。DRA006903–DRA006907 的五条模型在 Dataset 2 再次出现，只能作为批次复核，不能增加独立 n。

## 4. ATAC 预处理

### 4.1 PDX 与细胞系原始读段

- 所有下载和计算只在云服务器进行；
- Nextera adapter：cutadapt；
- PDX 先对 mouse genome 比对并去除宿主读段，再对 hg38 比对；宿主去除是 PDX 公共原始数据所必需的技术处理，不改变生物比较；
- human alignment：Bowtie2；
- 去除未比对、secondary/supplementary、PCR duplicate、chrM、chrY、MAPQ < 30 和 ENCODE hg38 blacklist 重叠读段；
- 单端细胞系按锚点进行 Tn5 `+4/-5 bp` 校正，并用 MACS2：`--keep-dup all -B --shift -75 --extsize 150 --nomodel --SPMR -q 0.01`；
- PDX 为公开数据原有 paired-end 文库，主分析使用片段模式 `MACS2 -f BAMPE --keep-dup all -B --nomodel --SPMR -q 0.01`，不得人为丢弃一端伪装成单端文库；R1 单端锚点参数只作为敏感性分析；
- 只有真实生物/技术重复存在时才运行 IDR 0.05，不为单库模型伪造 IDR。

原始读段硬 QC：每个模型至少 100 万条去重复、MAPQ≥30、非 chrM/chrY、非 blacklist 的人源读段或片段，并至少 10,000 个 MACS2 峰。TSS enrichment 和 FRiP 全量报告；因三资源文库布局和建库流程不同，不用一个事后选择的 FRiP 阈值删样本。

ataqv 原始 JSON 保留，但其 HQAA/TSS 定义要求正确配对 reads，因而对单端细胞系不适用。为避免把工具的 `TSS=NA` 误写成生物学无信号，另对所有原始 PDX/细胞系 BAM 用统一的 Tn5 `+4/-5 bp` 插入位点计算链方向校正的 ±1 kb TSS profile，并按两端各 100 bp 平均值归一化；该指标仅报告，不进入硬 QC 或信号闸门。

### 4.2 TCGA 开放处理矩阵的预注册翻译

TCGA 原始/比对文件受控，但 GDC 公开每个技术重复的 hg38 固定峰原始计数、峰集和完整 QC。主分析合并前要求同一肿瘤的两个技术重复均达到 `CPM >= 1`，以该峰为该患者可及；这是从开放计数矩阵重建逐患者峰的预注册替代，不冒充 MACS2 原始读段重跑。

必须同时报告 `CPM >= 0.5`、`CPM >= 2` 以及“合并技术重复 CPM >= 1”三种敏感性定义。若要求从 FASTQ/BAM 完全复刻 TCGA 患者端，需另行申请 `phs000178`，但不阻断本轮公开处理数据分析。

## 5. 启动子、活动与 motif

### 5.1 启动子开放

- 参考：冻结版本的 GENCODE hg38 gene annotation；下载文件、版本和 SHA256 写入 receipt；
- 同一 gene symbol 的所有 basic/canonical protein-coding TSS 合并为一个基因启动子集合；
- 单样本任一峰与该 TF 的任一启动子相交即记为开放；
- 系统层面阈值使用可执行锚点代码的 `ceiling(n/2)`；无 QC 损耗时为患者 11/22、PDX 7/13、细胞系 10/19，若发生预注册允许范围内的硬 QC 损耗，则分母和阈值自动使用最终合格独立样本数；
- 三端都达标并且 TCGA-LUAD、PDMR-LUAD PDX、DepMap-LUAD 细胞系的平均 NES 均大于等于 0，才称为 LUAD HC-TF。

### 5.2 motif

- 按锚点 Methods 和 Supplementary Data 5 的实际交付，对 JASPAR 2024 CORE vertebrates 与 MEME Suite 发布快照中的 CIS-BP 2.0 Homo sapiens 分别运行 HOMER；同一 HC-TF 在两库均有 motif 时，两套结果均参与检验，不采用“JASPAR 优先、CIS-BP 仅补缺”的简化；
- 两个数据库在单样本内分别进行 Benjamini-Hochberg 校正；TF 层面若任一数据库中的任一对应 motif 达到 `q < 1e-5` 且 LOR>0，即记为该样本支持，同时完整保留两库所有 motif 结果；
- 外部 PWM 转为 HOMER 格式时采用 HOMER 随附 `parseJasparMatrix.pl` 的阈值 0 约定，并在 receipt 中记录该参数；
- 背景在看结果前冻结为：TCGA LUAD+LUSC 类型峰与全部通过 QC 的肺 PDX/细胞系峰的合并、去重可及区域；
- 系统层面仍用 `ceiling(n/2)`；Figure 2E 同时显示每系统平均 log2 odds ratio、样本间 SD 和达到/未达到一半样本的状态。

上述两库独立检验修正发生在任何 LUAD motif 结果生成之前，依据锚点补充表中分列的 `HOMER_JASPAR_TCGA_PDX_CellLine` 与 `HOMER_CISBP_TCGA_PDX_CellLine` 工作表完成，属于原文方法对齐，不是结果驱动的阈值调整。

## 6. Figure 与补图同步交付

- Figure 2A：三端数据与方法示意；
- Figure 2B：97 个冻结 TF 的三端启动子开放交集；
- Figure 2C：HC-TF 启动子开放比例与三端平均 NES；
- Figure 2D：motif 富集 TF 的三端交集；
- Figure 2E：motif LOR/SD/频率和最高平均富集 TF 的逐样本分布；
- Supplementary Figure 3A–C：逐样本峰数及真实重复/IDR 状态；
- Supplementary Figure 3D–F：1,000 次样本顺序置换、`nls(SSasymp)` 饱和曲线；
- Supplementary Figure 3G–I：系统内逐样本 Pearson 相关；
- Supplementary Figure 3J–L：promoter/exonic/intronic/distal peak 注释；
- Supplementary Figure 4A：97 TF → 三端启动子开放 → NES 过滤 → HC-TF 的流程；
- Supplementary Figure 4B：三端启动子开放组合和被 NES<0 排除的 TF；
- Supplementary Figure 4C：HC-TF 在 JASPAR only、CIS-BP only、both、none 的分布。

所有 panel 必须由同一 manifest、同一 QC 表、同一 TF SHA 和同一阈值配置生成。主图通过而补图缺失，不得称 Figure 2 完成。

## 7. 预先冻结的停止规则

### 7.1 数据闸门

任一条件失败即 `FAIL_DATA`：

1. 最终独立样本少于患者 7、PDX 10 或细胞系 10；
2. 任一系统超过 20% 预定义样本因硬 QC 失败；
3. TF 输入不是 97 个或 SHA256 不一致；
4. 三端无法统一到 hg38、GENCODE symbol 或独立模型 ID。

### 7.2 锚点逻辑闸门

同时满足才记为 `ANCHOR_PASS`：

1. 至少 10 个 TF 成为三端启动子开放且三端平均 NES≥0 的 HC-TF；
2. HC-TF 中至少 3 个 motif 可检验 TF 在三端均达到“至少一半样本 q<1e-5”；
3. 上一条件对应的三端 motif 支持 TF 至少占 motif 可检验 HC-TF 的 10%。

任一失败即 `FAIL_SIGNAL`，停止 Figure 3–5。阈值用于保证后续网络/模块分析不是由一两个 TF 支撑，不按结果放宽。

### 7.3 额外统计稳健性（与原文判据分开）

为判断三端交集是否只是高可及启动子或 motif 数据库覆盖造成，另做 10,000 次匹配置换：启动子按 GC、长度和 TCGA 表达分位匹配，motif 按数据库来源和 motif 数匹配。三端启动子交集和三端 motif 交集分别要求经验 `p < 0.05`。

- 锚点逻辑和置换都通过：`PASS`；
- 锚点逻辑通过、置换失败：`ANCHOR_PASS_ROBUSTNESS_FAIL`，明确表示原文式筛选可产生结果，但证据不足以进入 Figure 3，自动停止；
- 锚点逻辑失败：`FAIL_SIGNAL`。

该置换层是额外质量控制，不得写成 TNBC 原文方法，也不得用替代 TF 网络或放宽阈值“救回”结果。

## 8. 不允许的事后操作

- 不得把比较改成 EGFR-mutant vs WT、TRU/PIF/PPR、LUAD vs SCLC 来救回结果；
- 不得把 H1299、H2126、H1703 或误鉴定的 RERF-LC-OK 填入严格 LUAD 主分析；
- 不得把药物处理库、重复提交批次、技术重复当作独立 n；
- 不得用 CollecTRI、ULM 或其他 TF 方法替换 Figure 1 冻结 TF；
- 不得只交付主图而隐去 Supplementary Figure 3、4A–C 的 QC 与失败位置。
