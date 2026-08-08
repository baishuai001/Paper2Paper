# 近邻论文与数据的 donor 资产账

本文件不把近邻论文简单分成“能做/不能做”。每个来源同时回答两件事：它已经发表了什么具体主线，以及它还能为锚点文献到目标论文的模仿提供什么。前者约束实质重复，后者用于组合框架、数据、代码和方法。

截至 2026-08-08，下面的 `code` 状态只表示本轮是否找到并验收了公开实现；“未定位”不等于作者一定没有私有代码。

## 1. 锚点文献

**SPP-DONOR-ANCHOR** — *Targeting SPP1+TAMs associated with liver metastasis reverses immunosuppression and synergizes with immunotherapy in colorectal cancer*（[DOI](https://doi.org/10.1136/jitc-2025-014128)）

- **框架资产**：临床问题 → 跨解剖部位单细胞图谱 → 髓系聚焦 → SPP1/SELENOP 对照状态 → 功能/代谢 → 空间定位 → bulk 临床和免疫生态 → 干预。
- **可借 Figure 角色**：全景到对象、状态到功能、空间定位、患者级临床桥接、效应细胞桥接、干预收束。
- **数据状态**：多个公共来源可做部分分析复现，但精确整合对象、机构 mIF、动物源表和若干作者选择不可得。
- **代码状态**：未定位到作者公开代码；必须按模块从同类论文、方法论文和官方实现补足。
- **已占据的具体主线**：CRC/CRLM 中 SPP1+TAM 富集、CD8 免疫抑制以及药物组合的中心叙事。
- **继续使用方式**：目标论文应尽量保留这条证据阶梯；替换对象、疾病、关系、结局或 donor 框架时，把已登记的统计和因果修复叠加进去。

## 2. Cancer Letters 2026：不是只有一个 TREM2 标签

**SPP-DONOR-TREM2-2026** — *Multi-omic profiling identifies TREM2+ lipid-laden macrophages as inflammatory drivers and therapeutic targets of colorectal cancer liver metastasis*（[PubMed](https://pubmed.ncbi.nlm.nih.gov/42250756/)，[publisher](https://www.sciencedirect.com/science/article/pii/S0304383526004209)）

- **框架资产**：CyTOF、单细胞、空间、bulk 和 lipidomics 联合识别细胞状态；再连接空间生态位、肿瘤来源脂质/凋亡碎片、Kupffer-cell-related 状态假设、ALOX5/ALOX5AP—白三烯通路、中性粒细胞招募、肿瘤干性、预后、TREM2 细胞清除/白三烯抑制和 anti-PD-1。
- **可借 Figure 角色**：跨模态锁定状态 → 证明真实脂质表型 → 定位生态位 → 上游物质来源 → 下游多细胞效应 → 直接干预。这比“把 SPP1 换成 TREM2”丰富得多。
- **数据状态**：本轮从 PubMed 和 publisher record 没有核实到可下载 accession、原始多组学源表或患者映射；目前只能把它列为框架/方法 donor，不能宣称数据可复用。
- **代码状态**：本轮没有定位并验收公开代码仓库；不能把论文中的分析名称当成可运行实现。
- **已占据的具体主线**：在 CRLM 中把 TREM2+ lipid-laden TAM 定义为炎症驱动者和治疗靶点，并以白三烯、中性粒细胞、干性和 PD-1 联合完成机制/干预链。
- **仍可与锚点组合**：比较 SPP1/ECM-remodeling 与 TREM2/APOE/GPNMB lipid-handling 程序是重叠、连续还是不同状态；借用其“上游物质—空间生态位—下游效应细胞”结构研究另一个锚点对象；借用其多组学证据顺序修复锚点的单 marker 代谢分型；其可靠代码若以后找到并验收，可作为对应模块 donor。
- **必须修复/限制**：在公开数据和代码未核实前，只能提出转录层候选；lipid-laden、白三烯生成、细胞来源和治疗效应不能由公共 scRNA 打分代替。

## 3. 胃癌 2024：SPP1+TAM—CD8 论文仍可拆成多个 donor

**SPP-DONOR-GC-SPP1-CD8-2024** — *Potential crosstalk between SPP1+ TAMs and CD8+ exhausted T cells promotes an immunosuppressive environment in gastric metastatic cancer*（[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10870525/)）

- **框架资产**：GSE163558 跨原发/多器官转移单细胞图谱 → CD8 与巨噬细胞分别重聚类 → CellPhoneDB/CellChat → GDF15–TGFBR2 候选 → GDF15 抑制小鼠肝转移模型 → 多个 bulk 队列中的配体受体评分、预后和免疫特征。
- **数据资产**：GSE163558 的 10 个组织样本来自 6 位患者，含 3 个原发、2 个肝转移、2 个淋巴结转移、卵巢/腹膜各 1 个及 1 个邻近组织；处理后 MTX/TSV 和 SRA 原始数据公开。论文另用 TCGA-STAD、GSE15459、GSE57303、GSE62254 和 GSE84437。
- **代码状态**：全文没有给出可核实的代码仓库；Seurat、Monocle、CellPhoneDB、CellChat、GSVA、CIBERSORT 和 LASSO/Cox 名称不是可复用项目代码。
- **已占据的具体主线**：胃癌转移中 SPP1+TAM 与 exhausted CD8 的通信，尤其 GDF15–TGFBR2，以及基于该网络的 LR 评分。
- **仍可与锚点组合**：借用“分别重聚类两类细胞—候选配体受体—功能验证”的证据顺序；把 GSE163558 用作多器官转移数据 donor；把其 ligand–receptor 框架用于另一个发送/接收细胞组合；把其 bulk 桥接作为反面合同，重建无泄漏的锁定评分。
- **必须修复**：CellChat 只能产生假设；细胞不能充当患者重复；异质 bulk 的批次处理、特征筛选和外部验证要重新审计；动物每组 4 只且使用 t 检验，不能支撑宽泛机制结论；没有公开代码时必须另找可靠实现并加反向测试。

## 4. 胃癌 2026：空间多细胞生态位是 donor，不是封死胃癌

**SPP-DONOR-GC-CSF1-SPP1-2026** — *The CSF1+ tumor cell–SPP1+ macrophage axis drives gastric cancer progression and immunotherapy resistance*（[full text](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1817573/full)）

- **框架资产**：GSE167297/GSE163558 单细胞整合 → GSE251950 空间映射 → PRJEB25780 免疫治疗队列去卷积 → 巨噬/成纤维细胞亚群 → CellCall 通信 → CSF1+ malignant cell—SPP1+ macrophage—MFAP5+ CAF/CD8 多细胞生态位 → Transwell/ELISA/增殖侵袭 → mIHC。
- **数据资产**：文章明确列出 GSE167297、GSE163558、GSE251950、PRJEB25780 和 TCGA-STAD，公共数据入口可定位；仍需真实解析患者、组织、治疗和空间切片字段。
- **代码状态**：全文和 data availability 没有给出代码链接；目前不可视为合格代码 donor。
- **已占据的具体主线**：胃癌中 CSF1+ 肿瘤细胞诱导 SPP1+ 巨噬细胞并与免疫治疗抵抗相关。
- **仍可与锚点组合**：借用“上游肿瘤/基质来源 → TAM 状态 → 第三类细胞/治疗结局”的三元生态位；借用公共空间和 ICB 数据作为其他对象的验证数据；把它的 source-attribution 思路用于锚点的 SPP1 来源修复。
- **必须修复**：文章用 `maxstat` 选择连续变量 cutoff，并用非胃癌的 IMvigor210 支持治疗推断；目标论文应锁定模型/阈值并使用疾病匹配队列。空间共定位不等于方向性通信；体外共培养不能替代体内必要性；文章自己也承认单细胞/空间样本量和体外验证边界。

## 5. 胃癌转移数据 donor：不同数据支持不同模仿路线

### SPP-DONOR-GC-LM-NK — GSE246662

- [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246662) 提供 3 位患者的原发胃癌和配对肝转移，以及 3 位健康供者肝组织；197.1 MB 的 CSV TAR 和 SRA 原始数据公开。
- 数据生成论文已经研究“胃癌肝转移的免疫景观和 NK 功能受损”，所以不能把这句话原样当目标论文。
- 它仍可支持锚点式的“髓系状态—NK 功能”候选、配对原发—肝转移比较和健康肝背景控制。3 对患者只适合发现/方向验证，不能被总细胞数放大成大样本。

### SPP-DONOR-GC-MULTISITE — GSE163558

- [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163558) 提供来自 6 位患者的 10 个原发/邻近/肝/淋巴结/卵巢/腹膜样本，处理后矩阵公开。
- 优点是多转移部位；限制是每个部位患者很少、样本并非完整配对，不能把组织标签当作患者独立重复。

### SPP-DONOR-GC-ATLAS — GSE239676

- [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239676) 是 20 位 treatment-naive stage-IV 胃腺癌患者的 68 个样本，含原发灶、9 位肝转移、6 位腹膜转移、1 位两者兼有、4 位腹膜加卵巢转移，以及邻近组织/血液；同时有单细胞转录组和免疫受体组。
- GEO 明确写明 raw data 因隐私不提供；论文声明处理后细胞表达矩阵和 metadata 可下载。它适合作为处理后数据 donor，但原始 FASTQ 复算被阻断。
- 原论文已经比较肝与腹膜转移的肿瘤/TME 共演化并聚焦 ferroptosis；目标路线应利用其配对、多部位和 TCR 结构，而不是重述 ferroptosis 主线。

### SPP-DONOR-GC-ROUTE-BULK — GSE237876

- [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237876) 是 14 位患者的 66 个原发与肝/腹膜/卵巢等转移组织 bulk whole-transcriptome 样本，处理后 TXT 和 SRA 原始数据公开；它不是单细胞队列。
- 一位患者可有多个原发区域、转移灶和重复，因此统计单位必须是患者，不能把 66 个样本都当独立患者。
- 原论文聚焦 route-specific msEMT/CAF。它可作为锁定免疫细胞程序的患者配对 bulk 验证 donor，但不能用同一队列既筛程序又声称外部验证。

## 6. 癌种替换不止胃癌：目前已核实的另一个入口

**SPP-DONOR-PDAC-LM — GSE263733**（[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263733)）

- 公开记录描述 treatment-naive PDAC 原发灶、正常胰腺和部分配对肝转移的单细胞研究，并提供 SRA/处理后资源入口。
- GEO 的 218 个 sample 条目包含测序 lane/拆分记录，不能直接解释为 218 位患者；必须解析患者—组织—library 映射后才判断可用患者数。
- 原论文聚焦恶性导管亚型、克隆演化和 Treg 生态；它既是“换癌种”的数据候选，也可能提供 malignant–immune 框架。是否能形成 SPP1/其他髓系目标路线仍需最近发表和字段级预检。

这只是本轮已核实的第二个癌种入口，不代表癌种搜索已经完成。乳腺癌、黑色素瘤、肺癌、神经内分泌肿瘤等肝转移候选尚未完成同等深度审计，不能因未列出而视为被排除。

## 7. 当前 donor 使用边界

- 论文框架可借，不代表论文的所有统计选择都可照搬；每项硬伤必须在目标论文修复。
- accession 可定位不等于数据合同已通过；本文件不能把 `metadata_checked` 提升为 `sample_parsed`。
- 方法名称可定位不等于代码可复用；没有 commit、license、环境、入口、测试和真实输入 smoke 的实现仍是未验收候选。
- 已发表的字面主线应进入重复风险账，但 donor 资产继续保留。后续路线比较必须同时引用这两张账。
