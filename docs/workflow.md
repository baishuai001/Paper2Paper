# 从锚点文献到目标论文

Paper2Paper 使用两种 workspace：Pilot 用于审计和最小真实试验；manuscript project 用于用户批准后的正式论文生产。二者不能通过改一个状态相互替代。

## Pilot workflow

```text
anchor_audit → route_generation → verification → pilot_review → pilot_complete
任何阶段均可 → stopped
```

### 1. `anchor_audit`

目标：理解锚点文献真正做了什么、为什么这样组织证据，以及公开材料能支持到哪里。

完成 `anchor/audit.md`，至少包括：

- 文献身份、科学问题、中心关系、主要结局和 falsifier；
- 每幅主图/补图的数据、metadata、方法、代码、统计单位、输出和 claim；
- 数据和队列在发现、拟合、验证、定位、机制等环节的角色；
- 分析步骤与 Figure 叙事；
- 可保留、修复、替换、扩展、删除和无法复现的部分；
- 每项硬伤的修复设计与剩余风险；
- 私有数据、湿实验、缺失代码和未报告选择造成的 claim ceiling。

退出条件：初学者能够理解这套框架为什么值得借鉴、哪里失效、怎样在目标论文中避免同一问题。

### 2. `route_generation`

目标：从锚点特有框架生成完整但不过度膨胀的目标论文候选。

每条路线写入 `evidence/routes.tsv`，明确：

- 一句话问题、目标疾病/对象、统计单位、比较、主要结局；
- 从锚点保留什么、修复什么、改变什么；
- claim ceiling、falsifier 和停止条件；
- 最低数据、代码和 Figure；
- 时间、计算和初学者负担；
- 最近发表与实质重复风险。

候选轴必须随当前文献扩展。任何常见替换方式都只是例子，不得作为固定清单。新路线先标为 `candidate`；`active` 和 `backup` 只能在可比的最低预检后产生，不来自隐藏评分。

### 3. `verification`

验证按成本逐步升级：

#### `direction_audited`

- 问题、统计单位、比较、结局和 Figure 最低闭环清楚；
- 最近发表初查未发现立即阻断的具体重复；
- 数据与代码需求已写出。

#### `availability_prechecked`

- 必需数据文件已访问并解析代表性内容；
- 必需字段、患者/样本 ID、队列用途和独立性已核对；
- 必需代码有合格实现、可靠 donor 组合或有依据的重建方案；
- 许可、版本、环境、入口和输入输出合同明确；
- 固定 signature/模型具有可计算规范；
- 最近发表和最接近论文已登记。

#### `minimal_real_run`

- 关键模块在目标环境读取真实文件并非交互运行；
- `runs.tsv` 保留输入来源、命令、环境、时间、退出码、日志和产物；
- 至少生成一个代表性 source table 或 Figure 片段；
- `results.tsv` 说明结果支持、削弱、矛盾还是无法判断；
- 失败、阴性结果和限制保留。

单元、合成和环境 smoke 运行不能冒充 `real_data`。

#### `figure_loop_closed`

最低主图都有真实数据、分析代码、source table 和验收条件，且可支持与不可支持的 claim 已分开。该状态只描述 Figure 证据，不自动表示值得发表。

### 4. `pilot_review`

AI 提交普通语言选择包：

- 当前首选与一个备选；
- 每条路线的具体问题和用途；
- 已完成与未完成证据；
- 数据、代码、Figure、时间和初学者负担；
- 最近论文、具体重复风险和可区分点；
- 锚点硬伤、目标论文修复、修复证据和剩余风险；
- 最大风险、停止条件、推荐及会改变推荐的新证据。

若最低证据仍不可比，AI 继续补证，不让用户从未核验清单中猜选。用户只负责批准或否决正式路线、中心 claim 和重大改变。

### 5. `pilot_complete`

Pilot 以明确决定结束：

- `promote`：建立独立目标论文项目；
- `retain_training`：保留为学习案例；
- `close_pilot`：审计完成但无稿件路线；
- `stop`：记录明确阻断。

每个 Pilot 还要记录运行中发现的 workflow 问题、当前实例修复、是否修改 Paper2Paper 以及修改后的复验。没有 workflow 新问题也是正常结果。

## Manuscript-project workflow

```text
specification → execution → interpretation → writing → complete
任何阶段均可 → stopped
```

### 1. `specification`

核对来源 Pilot/route/run，冻结队列、排除规则、统计单位、比较、结局、协变量、预处理、缺失处理、模型、参数、随机种子、数据 manifest、代码版本、Figure、source table、falsifier 和 claim ceiling。

### 2. `execution`

从冻结输入稳定生成 source tables 和 Figures。每次运行保留命令、commit、环境、manifest、日志和产物；失败运行不删除。输入、参数或代码改变时建立新 run，不能覆盖历史。

### 3. `interpretation`

让 claim 服从结果。将支持、削弱、矛盾、阴性或不确定结果写入 `results.tsv`，记录替代解释、限制和 `continue/refine/reroute/stop` 影响。看过结果后新增的分析默认标为探索性。

### 4. `writing`

Methods 来自冻结规范和运行记录，Results 来自 source tables，Figure legend 写清统计单位和检验，Discussion 区分结果、解释和推测，Limitations 说明数据缺失、代码重建和证据边界，并更新最近发表检索。

### 5. `complete`

只有稿件、Figures、source tables、代码/依赖、数据访问说明、运行记录、结果和限制闭环，并由用户在最新合格真实运行后记录 `approve_release`，才可完成。

## 自动校验的边界

`paper2paper validate` 检查结构、枚举、引用、必要文件和部分阶段依赖。它不重新下载数据、不证明 `real_data` 标签真实、不判断文献检索完整，也不替代科学审查或用户决定。
