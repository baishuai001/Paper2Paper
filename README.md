# Paper2Paper

Paper2Paper 是一个面向初学者、由 AI 协助执行的论文模仿 workflow。它以真实文献为锚点，把论文中的问题、数据、方法、代码和 Figure 拆开审计，再将其中值得借鉴的结构转化为一条能够真实执行、根据结果调整并最终形成论文的路线。

它的终点是论文，不是方向推荐报告，也不是不断扩建的研究平台。

当前版本、真实 Pilot、已验证范围和未决事项见 [PROJECT_STATE.md](PROJECT_STATE.md)。

## 目标与两个循环

Paper2Paper 同时维护两个循环，但二者有明确主次。

### 主循环：把真实文献推进为论文

```text
锚点审计
  -> 候选路线组合
  -> 数据与代码可得性预检
  -> 最小真实数据运行
  -> 用户选择或停止
  -> 正式论文项目
  -> 冻结分析、完整执行、解释结果和写作
```

主循环回答的是：这篇文献可以怎样借鉴，哪条路线科学上成立、数据可取、代码可运行、Figure 能闭环，而且没有与现有论文实质重复？

### 次循环：让真实文献反过来完善 Paper2Paper

```text
Pilot 中发现问题
  -> 判断问题属于锚点、Pilot 实现、分析模块、核心规则还是科学结果
  -> 先在当前实例处理
  -> 满足晋升门槛后进入 module 或 core
  -> 只重跑受影响的回归案例
  -> 用新规则重新检查相关旧 Pilot
```

次循环只服务主循环。某项平台工作如果不能减少真实论文的错误、时间或审查负担，就应推迟；当一条正式路线已经具备真实数据、可靠代码、Figure 闭环和可写主张时，应停止通用化建设并进入论文生产。

完整机制见 [学习循环](docs/learning-loop.md)。

## 四层架构

| 层次 | 责任 | 不能做什么 |
| --- | --- | --- |
| `core` | 通用 schema、校验器、状态和项目纪律 | 不保存某篇论文的公式、预期结果或生物学结论 |
| `module` | 可复用分析能力及其输入输出合同、实现、许可和测试 | 不因在一个 Pilot 中跑通就声称适用于所有论文 |
| `pilot` | 对一篇真实锚点文献完成审计、路线组合和最小真实运行 | 不冒充正式论文项目，也不直接宣布核心规则已经通用 |
| `manuscript-project` | 用户批准后的正式论文路线、冻结分析、Figure 和稿件 | 不为了平台通用化而无限增加与稿件无关的功能 |

`registries/` 是四层之间的可追溯索引，不是第五个研究层。它记录能力、模块、晋升决定和回归案例。详细边界见 [架构说明](docs/architecture.md)。

## 可以怎样模仿

以下改写都可以进入候选组合：

- 原样复现，作为学习、结论核查或模块基线；
- 替换 marker、单个基因或基因集；
- 替换核心细胞类型；
- 替换癌种；
- 扩展为泛癌；
- 构建或重算转录组 signature；
- 组合多个有明确作用的替换；
- 保留锚点的科学与 Figure 结构，改用其他论文、方法论文或官方实现中的可靠代码；
- 修复锚点论文的数据、统计、计算或验证设计问题。

Paper2Paper 不按“改动有多大”排序，也不设置抽象的创新性准入门槛。简单替换可以接受；但科学设计失效、必需数据不可得、代码无法可靠运行，或与已发表工作实质重复的路线不能进入正式论文项目。

## 每篇真实文献必须产生两类结果

### 文献侧结果

- 锚点论文的科学语法和 Figure—证据映射；
- 可保留、修复、替换、删除和无法复现的部分；
- 完整但不过度扩张的候选路线组合；
- 每条路线的数据、代码、统计单位、重复发表和初学者负担；
- 至少一项能尽早推翻路线的最小真实运行；
- `continue`、`refine`、`reroute` 或 `stop` 的证据化决定；
- 值得进入正式论文项目的候选，或明确的停止理由。

### 产品侧结果

- 新发现的问题及其作用范围；
- 可复用模块、数据解析器、代码 donor 或科学不变量测试；
- 是否晋升到 module/core 的决定与成熟度；
- 哪些既有 Pilot 会受到影响，是否需要重跑。

并不是每篇论文都必须修改 Paper2Paper 核心。只完善该文献而没有发现通用缺陷，也是完整而正常的结果。

## 能力矩阵，而不是单一训练基准

Paper2Paper 不用一篇文献或一个运行作为全局验收基准。每个真实案例只覆盖它实际运行过的能力，例如：

- 固定 bulk signature 重算；
- 微阵列探针映射；
- 患者级生存分析；
- 单细胞髓系重聚类；
- 患者级 pseudobulk 推断；
- 空间共定位；
- 跨尺度证据组合。

能力矩阵分别登记参考案例、独立迁移案例、失败案例和成熟度。胃癌 GN-R01 的数据、公式和预期结果只能回归固定 bulk-signature 示例，不能验证 SPP1+TAM、单细胞、空间或整个 Paper2Paper。

## 数据、代码与执行纪律

每条路线必须留下可以让用户复查的证据：

- 数据需求、检索式、候选队列、真实文件、字段来源、ID 连接和队列用途；
- 代码模块需求、donor、版本、许可证、入口、目标环境、安装和真实输入 smoke test；
- signature 的特征顺序、系数、映射、归一化、缺失策略和 cutoff 合同；
- Figure 与数据需求、代码模块、source table 和验收测试的映射；
- 失败运行、阴性结果、限制和它们对论文路线的影响；
- 最近发表检索、最接近论文和实质重复判断；
- 用户的路线选择、中心主张和高影响变更决定。

结构化记录用于让用户审查 AI，不代替真实下载、运行或科学判断。

## 四种不同的“通过”

- 结构与合同校验：文件和登记关系一致；
- 单元/集成测试：实现满足已写明的技术和科学不变量；
- 真实数据运行：指定资源、环境和参数确实生成了产物；
- 科学审查：设计、统计单位、解释和结论边界合理。

这四层证据不能合并表述为“验证通过”。`paper2paper validate` 只检查结构、登记合同和关键跨阶段依赖；它不会自动重新下载数据、复算 checksum、核对许可证法律含义、执行完整分析或证明科学结论正确。

## 仓库结构

```text
Paper2Paper/
├─ src/paper2paper/       # core：通用 schema、校验器和 CLI
├─ modules/               # 可复用分析模块、合同和测试
├─ registries/            # 能力、模块、晋升和回归索引
├─ pilots/                # 相互隔离的真实文献 Pilot
├─ manuscript-projects/   # 用户批准后的正式论文项目
├─ docs/                  # 操作、架构和治理文档
└─ tests/                 # core 回归测试
```

跨论文复用的详细规则见 [模块复用](docs/module-reuse.md)，变更、PR 和审查纪律见 [项目治理](docs/project-governance.md)。旧 workspace 和字段的迁移边界见 [Schema 2.0 migration](docs/schema-2-migration.md)。

## Quick start

需要 Python 3.10 或更高版本。

```bash
python -m pip install -e .

paper2paper init pilots/my-project \
  --project-id P2P-MY-PROJECT \
  --title "My paper project" \
  --anchor-title "Anchor paper title" \
  --doi "10.xxxx/xxxxx"

paper2paper validate pilots/my-project
paper2paper status pilots/my-project
paper2paper next pilots/my-project
paper2paper report pilots/my-project
paper2paper validate-registry .
paper2paper coverage .
```

`init` 默认只建立 Pilot，不建立稿件。Pilot 完成真实证据审查并记录用户的
`approve_route_evidence` 和 `promote` 决定后，使用：

```bash
paper2paper promote pilots/my-project ROUTE-1 manuscript-projects/my-paper \
  --project-id P2P-MS-MY-PAPER \
  --title "My manuscript project"
```

`promote` 建立独立的 manuscript project，并锁定来源 Pilot、route 和 passed runs；它不会
把 Pilot 运行结果直接当成正式论文结果。`next` 根据当前已登记事实列出下一步，不能替用户选择
路线，也不能生成脱离真实数据和代码证据的方向结论。

## 数据、凭据和许可证边界

Git 只保存自有代码、配置、较小 metadata、检索记录、manifest、经过再分发审查的 source table、汇总结果、报告和稿件。不要提交服务器密码、token、患者可识别信息、受控数据、未经许可的 PDF 或大型测序文件。

逐患者派生表即使只含公开编号，也必须先核对原始数据的再分发和署名条件；未完成时只提交生成代码、manifest、汇总统计和不含逐患者记录的审计产物。

Paper2Paper 自有代码采用 [MIT License](LICENSE)。外部数据、论文附件和第三方代码仍受各自许可证及使用条款约束。

项目的绑定纪律见 [RULES.md](RULES.md)。
