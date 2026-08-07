# 从锚点文献到正式稿件

Paper2Paper 使用两种不同 workspace：

- `workspace_kind=pilot`：审计一篇真实锚点、比较路线并完成最小真实试验；
- `workspace_kind=manuscript_project`：承接用户批准的路线，完成冻结分析、Figure 和稿件。

Pilot 不能通过修改 stage 直接变成正式论文项目。晋升时建立新的 manuscript project，并在 `PROJECT.json.provenance` 记录来源 Pilot、route 和 run。

## Pilot workflow

```text
anchor_audit
  -> route_generation
  -> verification
  -> pilot_review
  -> pilot_complete

任何阶段都可以 -> stopped
```

### 1. `anchor_audit`

目标：理解锚点论文真正做了什么，以及哪些证据无法公开复现。

操作：

- 完成 `anchor/audit.md`；
- 按 Figure 拆出科学问题、数据、方法、代码、输出和 claim；
- 还原该论文特有的中心关系、证据顺序、分析模块和 Figure 叙事，把它们视为可利用的框架资产；
- 标记可保留、修复、替换、扩展、删除和无法复现的模块；
- 对每项硬伤同时记录受影响的 claim、仍有效的框架、目标论文中的修复、验证条件和剩余风险；
- 记录私有数据、湿实验和未公开作者选择形成的 claim ceiling；
- 区分患者、样本、细胞、spot 和切片等统计/生物学单位。

退出条件：用户能够理解每幅主图在证据链中的作用、这套框架为什么值得借鉴、每项硬伤如何修复或限制主张，且不可公开复现的边界已写明。

### 2. `route_generation`

目标：形成完整而不过度扩张的候选组合。

操作：

- 从当前锚点实际具有的研究对象、关系、尺度、结局、证据顺序、模块、数据角色和 Figure 结构生成改变轴，不使用固定替换清单限制候选；
- 对每个框架要素比较保留、修复、替换、扩展或删除；原样复现、marker、细胞、癌种、signature 等仅作为非穷举工具箱；
- 把锚点硬伤的修复默认叠加到所有相关目标路线，而不是自动把“修复硬伤”另立为目标论文方向；只有当稳健性或方法学本身被明确选为中心研究问题时，修复才可以成为独立路线；
- 简单路线与复杂路线享有同等候选资格。只能因一条具体路线的科学设计、数据、代码、Figure 或实质重复证据停止它，不能因其属于简单 signature、泛癌或单轴替换而整类排除；
- 每条路线写入 `evidence/routes.tsv`；
- 使用 `route_role` 区分稿件候选、训练复现和辅助分析；
- 使用 `capability_ids` 说明路线需要哪些能力；
- 明确问题、统计单位、比较、主要结局、claim ceiling、falsifier 和停止条件；
- 估计数据、代码、初学者负担和日历时间；
- 将 `evidence_stage` 设为实际达到的阶段，不能提前写成最小真实运行；
- 新生成但尚未完成可比最低预检的路线一律标为 `candidate`。`active` 和 `backup` 表示预检后的执行选择，不是由隐藏评分算法产生的“最佳”和“次佳”。

退出条件：锚点特有框架中的重要可变要素没有被无理由遗漏，且每条 candidate 路线都能说明“从锚点文献保留了什么、修复/改变了什么、怎样形成目标论文”。此阶段只形成候选池，不允许依据未经核查的印象提前指定首选。

### 3. `verification`

目标：按成本漏斗把“可能做”变成有真实证据的判断。

#### 方向审计：`direction_audited`

- 科学问题、统计单位、比较和主要结局清楚；
- 最近发表初查没有立即发现实质重复；
- 最低数据、代码和 Figure 合同已经列出。

#### 可得性预检：`availability_prechecked`

- 必需数据字段已经映射到实际文件；
- 跨来源 ID、患者数和队列用途已核查；
- 必需代码有候选 donor、许可、版本、入口和环境信息；
- signature 或固定模型具有可计算合同；
- 最近发表检索和最接近论文已登记。

#### 最小真实运行：`minimal_real_run`

- 每个关键数据需求至少解析代表性真实文件；
- 最可能推翻路线的代码模块在目标环境运行；
- 至少生成一个代表性 source table 或 Figure 片段；
- `run_kind=real_data` 的 passed run 保留输入来源、退出码、命令、环境、日志和全部 artifact；
- 单元、合成和环境 smoke run 单独登记，不能提高真实执行证据阶段；
- 成功、失败和阴性结果均保留；
- finding 已按 Pilot/module/core/scientific result 判断作用范围。

#### Figure 闭环：`figure_loop_closed`

- 最低主图均有真实数据、模块、source table 和验收条件；
- 结果能够支持的 claim 和不能支持的 claim 已说明；
- 该状态适用于确实已完成 Figure 闭环的路线，不是 Pilot 完成的统一要求。

在 AI 推荐首选和备选之前，进入短名单的路线必须完成同一组最低预检字段：具体科学问题与统计单位、必需数据及字段、代码 donor 或重建方案、最低 Figure 闭环、初学者实施负担，以及最接近论文与实质重复风险。最低预检深度要可比，但不要求所有候选完成同等成本的完整下载或真实运行；高成本验证只对预检后仍值得推进的路线进行。

退出条件：至少一条路线得到可继续、修正、晋升或停止的真实证据；不要求所有候选完成同等成本的运行。

### 4. `pilot_review`

目标：AI 在完成证据搜集和可比最低预检后，先给出明确首选、备选、反对理由和不确定性；用户据此审查高影响路线决定，而不是替 AI 完成原始搜集、逐表核对和初步排序。

AI 提交的审查包必须包括：

- 开头先用普通语言给出 AI 的当前首选、一个备选、为什么这样判断、什么新证据会改变推荐；如果证据不足以推荐，必须说明缺哪项证据并由 AI 继续补齐，不能把“无法选择”包装成让用户自行比较；
- 每条 active/backup 路线的一句话问题和用途；
- 当前 evidence stage；
- 数据、代码、Figure、时间和初学者负担；
- 实际运行结果和未运行部分；
- 最近发表、实质重复风险和可区分点；
- 锚点硬伤、拟采用的修复、修复证据、剩余风险以及保留下来的框架价值；
- 最大风险、停止条件和 AI 的推荐；
- 拟晋升到 module/core 的 finding，以及尚未证明的范围。

详细表格是上述判断的可追溯依据，供用户需要时抽查；它不是要求初学者逐行完成的作业。用户的职责是批准或否决正式路线、中心 claim、重大 reroute 和最终发布，AI 对调查充分性、比较质量和推荐依据负责。

若某个 `manuscript_candidate` 要进入晋升评审，用户需在最后一次合格真实运行之后记录带时区
完整时间戳的 `approve_route_evidence`。这表示用户认可证据包可供决策，不等于批准最终稿件
或保证投稿成功。最终 `promote` 决定必须更晚，旧批准不能被后续新增证据自动消费。

### 5. `pilot_complete`

Pilot 只能以明确决定完成：

- `promote`：用户批准一条 route，随后建立独立 manuscript project；
- `retain_training`：保留为训练、参考或回归案例；
- `close_pilot`：没有晋升路线，但审计和停止判断完整；
- `stop`：存在明确阻断理由。

完成 Pilot 不会自动提升 module/core 成熟度。Promotion 和 regression evidence 仍需在中央 registries 中独立审查。

## Manuscript-project workflow

```text
specification
  -> execution
  -> interpretation
  -> writing
  -> complete

任何阶段都可以 -> stopped
```

### 1. `specification`

目标：在查看完整结果前冻结正式分析。

操作：

- 核对 `PROJECT.json.provenance` 与来源 Pilot/route/run；
- 完成 `analysis/specification.md`；
- 固定队列、排除规则、统计单位、比较、结局、协变量、参数和随机种子；
- 固定数据 manifest 和 module release IDs；
- 区分主要、验证、敏感性和探索性分析；
- 定义 Figure、source table、scientific invariants、falsifier 和 claim ceiling。

退出条件：完整分析不再依赖未记录的临场选择。

### 2. `execution`

目标：从冻结输入稳定生成 source table 和 Figure。

操作：

- 每次运行记录命令、commit、环境、data manifest、module releases、日志和产物；
- 失败运行也必须保留；
- 先验证 source table，再美化图形；
- 输入、关键参数或代码改变时产生新 run，不能覆盖旧记录。

退出条件：最低正文图和关键补充分析可以从已记录输入重新生成。

### 3. `interpretation`

目标：让 claim 服从结果，而不是让结果服从预想故事。

操作：

- 将结果写入 `execution/results.tsv`；
- 标记支持、削弱、矛盾或不确定；
- 记录限制、替代解释和稿件影响；
- 结果后的 `refine`、`reroute` 或 `stop` 进入 decisions；
- 事后新增分析默认为探索性。

退出条件：所有中心 claim 有证据，或已被删除、降低和重新限定。

### 4. `writing`

目标：形成不过度解释且可追溯的稿件。

操作：

- 完成 `manuscript/draft.md`；
- Methods 来自冻结规范和运行记录；
- Results 来自 source table；
- Figure legend 写清统计单位和检验；
- Discussion 区分结果、解释和推测；
- Limitations 包含数据缺失、队列偏倚、代码重建和证据层级；
- 更新最近发表检索。

退出条件：标题、摘要和正文不超过 claim ceiling，稿件与 artifact 能共同审查。

### 5. `complete`

用户记录 `approve_release` 后才能完成。交付至少包括：

- 稿件与 Figure；
- source tables；
- 代码、依赖和配置；
- 数据 manifest 和访问说明；
- 运行与结果记录；
- 限制、失败和不可复现边界。

## 阶段变更原则

阶段不是进度装饰。提高 `PROJECT.json.stage` 前必须满足该 workspace 对应阶段的退出条件。

`paper2paper validate` 只检查结构、登记合同和可机械判断的跨阶段依赖。它不能证明科学内容真实、最近发表检索完整或稿件值得发表，也不能替代用户审查。
