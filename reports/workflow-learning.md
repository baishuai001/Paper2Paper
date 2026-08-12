# Paper2Paper 跨 Pilot 学习报告

生成日期：2026-08-12。

## 这份报告是什么

它把不同真实文献 Pilot 中标记为 workflow 候选的问题汇集到一起。原始问题仍保存在各 Pilot；这里的候选不自动等于 Paper2Paper 产品缺口。

分类含义：

- `paper_specific`：只在当前课题解决，不改产品核心；
- `paper_type`：可能适用于同类论文，需用真实实例确认；
- `general`：跨论文类型都可能导致相同误判；
- `verified`：最小修复已实现，并留下来源 Pilot 与复验引用。
- `observed` 只是待归类事实，不表示已确认产品缺口；

聊天上下文不是持久记忆；本报告、来源 `issues.tsv`、`learning/findings.tsv` 和 Git 历史才是跨对话依据。

## 汇总

- 扫描 Pilot：2
- 扫描问题：45
- workflow 原始候选：29
- 已归类来源问题：29
- 尚未归类候选：0

## 已归类发现

| ID | 范围 | 状态 | 能力领域 | 归类后的问题 | 来源 |
|---|---|---|---|---|---|
| P2P-L001 | general | verified | learning_loop | 各 Pilot 能记录本地问题，但不同对话此前不能汇集、分类并追踪产品修复与复验。 | P2P-GASTRIC-NRRS#PF-002;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I11 |
| P2P-L002 | general | verified | study_design | workflow 若只检查患者独立性而不记录队列在筛选、拟合、调参和 cutoff 中的角色，会错误标注外部验证。 | P2P-GASTRIC-NRRS#PF-003 |
| P2P-L003 | paper_type | verified | model_contract | 论文给出基因和系数并不等于模型可计算；变换、映射、重复和缺失策略及 cutoff 可能仍缺失。 | P2P-GASTRIC-NRRS#PF-004;P2P-GASTRIC-NRRS#PF-013 |
| P2P-L004 | general | verified | data_provenance | 一个 accession 可能不能提供表达、结局、注释等全部必需字段，若只登记数据集名称会制造假可得性。 | P2P-GASTRIC-NRRS#PF-009 |
| P2P-L005 | general | verified | data_identity | 跨文件连接若不核对标识唯一性、集合覆盖和未匹配数，可能把结局附到错误样本。 | P2P-GASTRIC-NRRS#PF-010 |
| P2P-L006 | general | verified | decision_boundary | 代码能运行和数据能解析不等于路线值得形成目标论文，也不等于用户已经批准中心主张。 | P2P-GASTRIC-NRRS#PF-015 |
| P2P-L007 | paper_type | verified | data_use_fit | 公开、规模大或字段存在不能证明数据能完成拟定 Figure/claim；必须先写明论文具体用途，再核对决定性组、设计混杂、测量空间、运行资源和治疗定义。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I04;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I05;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I07;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I08;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I10 |
| P2P-L013 | general | observed | path_portability | 有效文献材料可能因 Unicode、路径长度或工具路径支持而在科学审计前失败。 | P2P-GASTRIC-NRRS#PF-001 |
| P2P-L014 | general | verified | execution_environment | 依赖在一个用户或解释器环境安装成功，不能证明目标 OS 用户、路径和锁定版本能读取并运行。 | P2P-GASTRIC-NRRS#PF-011;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I25 |
| P2P-L015 | paper_specific | local_only | code_repair | 胃癌 Pilot 的首个绘图调用使用了当前 lifelines 不接受的参数。 | P2P-GASTRIC-NRRS#PF-012 |
| P2P-L016 | general | verified | data_provenance | 稳定 accession 或 URL 的内容仍可能更新，若不固定字节和 checksum，同一命令可读取不同输入。 | P2P-GASTRIC-NRRS#PF-014 |
| P2P-L017 | paper_type | observed | data_provenance | 依赖在线服务的分析若不保存请求、响应、日期和版本，后续无法区分服务变化与分析变化。 | P2P-GASTRIC-NRRS#PF-005 |
| P2P-L018 | general | verified | anchor_mapping | 用一个“可替代/部分可替代”标签概括整幅多面板 Figure，会掩盖不同子图的数据、代码和 claim 条件。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I12 |
| P2P-L019 | general | verified | method_reconstruction | 锚点论文、数据生成论文、作者代码和官方方法分别只提供部分事实；若选择一个来源全盘照搬，会把未报告选择、代码缺陷或不适配输入传播到全部后续 Figure。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I09;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I13;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I17;P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I18 |
| P2P-L020 | paper_type | verified | method_reconstruction | 在 K 尚未冻结时，把跨 seed、相邻 K 和不同输入的一对一匹配边直接做传递并集，会经由中间节点把同一运行、同一 K 的多个不同程序错误合并。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I22 |
| P2P-L021 | paper_specific | local_only | evidence_packaging | 手工证据打包命令把 checksum 清单自身纳入清单，形成不可满足的自引用。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I23 |
| P2P-L022 | paper_type | observed | method_reconstruction | 多工具判定中，正常参考的异常阳性率可以很低，但大量参考仍可能处于 filtered、not.defined 或其他未解决状态。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I27 |
| P2P-L023 | paper_type | observed | method_reconstruction | 官方工具生成的字段名或元数据值可能来自函数默认值，却不代表该参数在当前统计路径中实际生效；若只照抄输出列名，会制造方法冲突。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I28 |
| P2P-L024 | general | verified | data_provenance | 下游步骤若反复读取并重新哈希几十 GB 的上游输入，会把已经完成的分析拖成超时；若完全跳过核验，又会失去输入身份约束。 | P2P-HCC-SC-SPATIAL-NPJ-2026#HCC-I19 |

## 尚未归类的 workflow 候选

这些记录只说明真实工作中出现过问题并建议考虑 workflow 动作。AI 仍须判断它是课题特有、论文类型级还是通用问题；不能按数量自动扩建核心。

| 来源 | Pilot 路径 | 严重度 | 观察 | 对论文的后果 | 当前动作 |
|---|---|---|---|---|---|
| — | — | — | 当前没有未归类候选 | — | — |

## 修复闭环

1. 先在来源 Pilot 解决当前论文问题；
2. 把原始事实写入该 Pilot 的 `issues.tsv`；
3. 将候选归类为课题特有、论文类型级或通用；
4. 只有可减少真实论文错误、时间或审查负担时，才做最小核心修改；
5. 用来源 Pilot 复验；跨类型主张还需另一篇合适 Pilot；
6. 保存测试/运行引用后才标记 `verified`；
7. 立即回到目标论文的数据、Figure、分析与写作。

本报告不是能力成熟度排行、创新性门槛或路线评分表。
