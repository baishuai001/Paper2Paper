# Paper2Paper 四层架构

## 设计目标

Paper2Paper 的架构只为一个目标服务：让初学者能够从真实论文出发，在 AI 协助下形成一条科学上可辩护、数据可获得、代码可运行、结果可审查并能写成论文的路线。

架构必须同时避免两类错误：

- 把某篇论文的特殊事实推广成全局规则；
- 为了追求“通用平台”而延迟正在进行的论文。

因此，仓库分为 `core`、`module`、`pilot` 和 `manuscript-project` 四层。`registries/` 只负责连接和追踪四层，不承担科学分析。

## 四层及其责任

### 1. Core：通用纪律与机械校验

位置：`src/paper2paper/`、`tests/` 以及适用的通用文档。

Core 可以包含：

- workspace 类型、状态和 schema；
- 表格字段、枚举、引用和阶段依赖校验；
- 数据/代码/许可/路径/统计单位等跨论文安全纪律；
- CLI、报告生成和影响范围查询；
- 不依赖某篇论文结果的合成或逻辑回归测试。

Core 不可以包含：

- 某篇论文的疾病、marker、细胞或基因集默认值；
- 某个真实队列、患者列表、signature 系数或 cutoff；
- 某个 Pilot 的冻结预期数值；
- “胃癌案例通过，所以所有 signature 都可靠”一类跨范围结论。

Core 只校验已经能够机械判断的内容。无法机械证明的科学判断必须保留为有依据、可审查但仍带不确定性的决定。

### 2. Module：可复用的分析能力

位置：`modules/`，并在 `registries/module_releases.tsv` 登记。

Module 是一个有边界的能力，而不是从论文仓库复制出来的整套 pipeline。一个模块至少具有：

- 稳定的 `module_id` 和版本；
- 一个或多个 `capability_id`；
- 输入、输出和失败条件；
- 统计单位和科学不变量；
- donor 与方法依据；
- 许可证和依赖；
- 单元、合成、真实参考和迁移测试状态；
- 已验证和未验证的适用范围。

示例能力可以是“从 GPL570 表达矩阵按固定规则重算一个锁定 signature”，而不能笼统写成“复现胃癌论文”。

模块可以由作者代码、官方包、方法论文、同类高质量论文和 Paper2Paper 适配代码共同组成。不同来源各自负责哪一步必须分开记录。

### 3. Pilot：一篇真实锚点文献的审计与最小试验

位置：`pilots/<pilot-id>/`。

Pilot 保存：

- 锚点论文及公开复现边界；
- Figure—数据—代码—claim 映射；
- 候选路线组合；
- 数据和代码检索、实际文件和最小运行；
- 失败、阴性结果、限制和决定；
- 文献侧产出与产品侧发现。

Pilot 的终点是形成明确判断，而不一定是一篇稿件。可能的结束方式包括：

- 找到可晋升路线；
- 保留为训练/参考案例；
- 因数据、代码、科学设计或重复发表停止；
- 等待明确的外部资源或用户决定。

Pilot 不直接修改 Core。它先记录 finding 和当前实例处理，再通过独立 promotion 进入 module/core。

### 4. Manuscript project：正式论文生产

位置：`manuscript-projects/<project-id>/`。

只有用户批准的路线才能建立正式项目。它必须记录来源 Pilot 和路线，但拥有独立的：

- 研究问题、中心 claim 与 claim ceiling；
- 队列和分析规范；
- 冻结参数、模块版本和数据 manifest；
- 完整 runs、results、Figure 和 source table；
- 稿件、限制、数据代码可得性和发布批准。

正式项目可以根据结果 `refine`、`reroute` 或 `stop`。如果 reroute 改变中心问题、主要结局或核心数据，应由用户批准并产生新的分析版本；不能用事后修改掩盖阴性结果。

正式项目必须是 `manuscript-projects/` 的一级子目录，且安全格式的 `project_id` 在整个仓库中
唯一。Pilot 不能预先批准发布；正式项目必须在自己的最新合格真实运行之后，以带时区的完整
时间戳重新取得用户批准。

## Registries：四层之间的索引

`registries/` 不是一个自动替用户做科学决策的数据库。它只保存跨层追踪所需的最小记录。

### `capabilities.tsv`

回答“Paper2Paper 实际验证过什么能力”：

- capability ID 与名称；
- 模态、数据层级和统计单位；
- 输入、输出和 claim 范围；
- 参考、迁移和失败案例；
- 由 regression cases 推导的覆盖状态。

### `module_releases.tsv`

回答“哪个实现提供这项能力”：

- module ID、版本和路径；
- capability ID；
- 合同、科学不变量、donor、许可和依赖；
- 测试层级和适用范围。

### `promotions.tsv`

回答“为什么一个 Pilot 发现进入了 module/core”：

- 来源 finding 和 Pilot；
- finding 类型、目标层和风险等级；
- 来源证据、独立迁移证据和用户/维护者决定；
- 当前成熟度、回滚或废弃原因。

### `regression_cases.tsv`

回答“改变一个规则后应重跑什么”：

- 被验证的 core control、module 或 capability；
- 对应 Pilot、route、run、测试和 artifact；
- reference、transfer 或 failure 案例类型；
- 最近验证版本和状态。

## 依赖和反馈方向

执行依赖从上游到下游：

```text
core contracts
   -> modules
      -> pilots
         -> manuscript-projects
```

知识反馈从下游通过显式晋升返回上游：

```text
真实结果/失败
   -> pilot finding
      -> promotion review
         -> module 或 core 变更
            -> 选择性回归
```

下游不能直接改写上游默认值；上游变更也不能自动改写旧 Pilot 的科学结论。旧结论只有在受影响案例完成新 run 后才能更新。

## 状态必须分轴记录

不要用一个 `ready` 覆盖所有含义。至少区分：

- 科学状态：问题和统计设计是否通过审查；
- 可得性状态：数据和代码是否真实取得并解析；
- 证据阶段：方向审计、预检、最小运行或 Figure 闭环；
- 运行类型：real data/unit/synthetic integration/environment smoke；
- 执行状态：planned/running/passed/failed/invalidated；
- Promotion 成熟度：observed/provisional/confirmed/rejected/deprecated；
- Module 成熟度：draft/unit verified/reference verified/transfer verified；
- 正式稿件状态：draft/frozen/executed/interpreted/release approved。

某项状态通过不能替代另一项。例如，代码测试通过不能证明科学问题成立，也不能满足
`real_data` 运行门槛；真实运行成功不能证明路线非重复；Pilot 完成不能证明正式稿件完成。

## 目标目录结构

```text
Paper2Paper/
├─ README.md
├─ RULES.md
├─ src/paper2paper/
├─ tests/
├─ modules/
│  └─ <module-id>/
│     ├─ README.md
│     ├─ contract.tsv
│     ├─ code/             # 提取出独立实现时才需要
│     └─ tests/            # 提取出独立实现时才需要
├─ registries/
│  ├─ capabilities.tsv
│  ├─ module_releases.tsv
│  ├─ promotions.tsv
│  └─ regression_cases.tsv
├─ pilots/
│  └─ <pilot-id>/
└─ manuscript-projects/
   └─ <project-id>/
```

目录本身不是完成证据。空白 Pilot、空 module 或仅有 registry 行都不能提升能力成熟度。
