# 从锚点论文到稿件

Paper2Paper按成本从低到高推进。每个阶段都必须产生具体文件；不能用一段漂亮的AI总结替代。

## 1. `anchor_audit`

目标：理解锚点论文真正做了什么，以及哪些证据无法公开复现。

操作：

- 完成 `anchor/audit.md`；
- 按Figure拆出科学问题、数据、方法、代码、输出和claim；
- 标记可保留、可修复、可替换、可删除和不可复现的模块；
- 记录私有数据、湿实验和未公开作者选择形成的claim ceiling。

退出条件：用户能够解释锚点论文每一幅主图在论文中的作用。

## 2. `route_generation`

目标：形成一个完整而不过度扩张的候选组合。

操作：

- 至少检查原样复现、marker、基因集、细胞、癌种、泛癌、signature和组合替换是否适用；
- 每条路线写入 `evidence/routes.tsv`；
- 明确研究问题、统计单位、比较、结局、claim ceiling、falsifier和停止条件；
- 初步估计数据、代码、初学者负担和日历时间，并在验证后更新；
- 不因路线只是简单替换而排除，也不因生物学叙述吸引人而提前选中。

退出条件：不存在未经说明而遗漏的常见路线模式。

## 3. `verification`

目标：把“可能做”变成“已经用最小证据验证可以做”。

操作：

- 建立数据需求，再搜索数据，而不是先看到数据再编问题；
- 建立代码模块需求，再搜索作者代码、官方实现、相关论文和重建donor；
- 每次数据、代码和文献检索写入 `search_log.tsv`，包括零结果；
- 下载并解析代表性数据文件；
- 核对疾病、组织、模态、统计单位、必需字段和最低受试者数；
- 安装代码并运行真实输入的非交互smoke test；
- 保存smoke input、output和必需测试证据；
- 映射最低主图并至少生成一个代表性输出；
- 检查最近发表和实质重复。

退出条件：至少一条路线没有未解决的科学、数据、代码、Figure或发表重合缺口。

## 4. `selection`

目标：由用户选择一条可执行路线，而不是由AI用模糊总分替用户决定。

操作：

- 比较每条可执行路线的数据规模、代码改造量、Figure闭环、初学者负担和主要风险；
- 记录 `select`、`backup` 或 `reject` 决定；
- `PROJECT.json.selected_route_id`只能指向一条通过验证的路线。

退出条件：一条路线被选择，其余路线保留为备份、训练或停止。

## 5. `specification`

目标：在查看完整结果前冻结分析选择。

操作：

- 完成 `analysis/specification.md`；
- 固定队列、排除规则、统计单位、比较、结局、协变量、参数、随机种子和signature公式；
- 区分discovery、validation、sensitivity和exploratory分析；
- 定义模块输入输出和科学不变量测试。

退出条件：完整分析不再依赖未记录的临场选择。

## 6. `execution`

目标：从已记录的输入稳定生成source table和图件。

操作：

- 每次执行记录命令、commit、环境、数据manifest、日志和产物；
- 失败运行也必须保留；
- 先验证source table，再美化图形；
- 不得覆盖旧运行来隐藏参数变化。

退出条件：最低正文图和关键补充分析可重新生成。

## 7. `interpretation`

目标：让claim服从结果，而不是让结果服从预想故事。

操作：

- 将每个结果写入 `execution/results.tsv`；
- 标记支持、削弱、矛盾或不确定；
- 记录限制、替代解释和稿件影响；
- `refine`、`reroute`或`stop`进入 `decisions.tsv`。

退出条件：所有中心claim有证据或已被删除、降低。

## 8. `writing` 与 `complete`

目标：形成不过度解释且可追溯的论文。

操作：

- 完成 `manuscript/draft.md`；
- 从已验证产物撰写Methods、Results和Limitations；
- 更新最近发表检索；
- 在干净环境重跑投稿配置；
- 用户记录 `approve_release`。

退出条件：稿件、图件、source table、代码、数据manifest和限制能够共同审查。

## 阶段变更原则

阶段不是进度装饰。提高 `PROJECT.json.stage` 前必须先完成该阶段的退出条件。
`paper2paper validate`会检查关键跨阶段依赖，但不能替代用户检查科学内容真实性。
