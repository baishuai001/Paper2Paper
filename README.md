# Paper2Paper

Paper2Paper 是一个面向初学者、由 AI 协助执行的论文模仿 workflow。它从一篇真实的**锚点文献**出发，珍惜并利用该文献自己的问题—证据—分析—Figure 框架，在修复硬伤、取得真实数据、补齐可靠代码并排除实质重复后，形成一篇可执行、可审查、可写作的**目标论文**。

终点是目标论文，不是候选方向报告、通用能力平台或不断增长的测试数量。

## 七条不可改变的要求

1. **面向初学者。** AI 承担资料搜集、技术判断、路线比较、代码构建和结果解释的主要工作，并用普通语言说明依据、反对理由和不确定性；不能把表格当作作业甩给用户。
2. **简单模仿路线可以接受并应被认真考虑。** 原样复现、局部替换、组合替换或其他方式没有先验高低；必须审查具体研究问题，不能按类别淘汰。
3. **数据必须能够取得。** accession 或网页不算可得；必需文件要能真实下载、解析、连接，并包含完成相应 Figure 和 claim 所需的字段。
4. **代码需要可靠、稳定，缺失时要从同类文献补足。** 可以组合锚点文献、同类论文、方法论文和官方实现中的代码，但必须记录来源、版本、许可、输入输出、适配、测试和真实运行证据。
5. **避免实质重复发表。** 只停止已经被充分发表的具体问题—对象—关系—结局—证据链组合；不能据此淘汰某一整类模仿路线，也不能丢掉近邻论文可提供的数据、代码和方法。
6. **最终目的是产出论文。** 当一条获批路线已经具备真实数据、可靠代码、Figure 闭环和可写主张时，停止无关的 workflow 扩建，进入完整分析和写作。
7. **不能重新引入 PaperRoute 的创新性门槛。** Paper2Paper 不按“创新程度”、改动大小或想象中的期刊档次设置准入门槛或隐藏排序。

## 主循环与次循环

Paper2Paper 保留两个互相促进、但有明确主次的循环。

```text
主循环：锚点文献审计
  → 目标论文候选路线
  → 数据/代码/近邻发表预检
  → 最小真实运行
  → 用户批准
  → 正式分析、Figure、解释与稿件

次循环：真实执行中发现 workflow 问题
  → 先修复当前文献
  → 必要时修改 Paper2Paper 文档、schema、代码或防错测试
  → 用同一真实文献复验
  → 回到主循环
```

次循环只服务主循环。一个改动若不能减少真实论文的错误、时间或审查负担，就不进入产品核心；次循环也不能因为“还可继续优化”而延迟已经能够推进的目标论文。

详见 [真实文献驱动的双循环](docs/learning-loop.md)。

## 怎样从锚点文献产生目标论文

每篇文献都有自己的研究对象、中心关系、尺度、数据角色、证据顺序、分析步骤、Figure 叙事和未报告选择。路线生成必须先还原这些文献特有要素，再逐项判断保留、修复、替换、扩展或删除。

marker、基因集、细胞、癌种、泛癌、signature 等只能作为例子，不能成为封闭菜单。遇到其他类型文献时，候选轴必须随该文献扩展。AI 应先形成完整但不过度膨胀的候选组合，再基于具体证据给出首选、备选、反对理由以及会改变判断的新证据；不能使用“代码优先”或其他隐藏加权算法。

### 硬伤必须转化为修复设计

指出硬伤不是目的，更不是嫌弃锚点文献的理由。每项硬伤都要同时回答：

- 哪个 Figure 或 claim 受到影响；
- 哪些框架资产仍然成立和值得保留；
- 目标论文怎样通过数据替代、统计修正、代码重建、验证重排或降低 claim 来避免同一问题；
- 什么最小证据能验证或推翻修复；
- 修复后还剩什么风险。

只有中心问题不可证伪、合理修复和等价替代均失败，或修复后已经不再保留锚点中心逻辑时，才停止整套框架。局部可用部分仍可继续借鉴。

## 两种 workspace

- `pilots/<id>/`：审计一篇真实锚点文献，生成目标论文候选，并完成低成本预检和最小真实运行。
- `manuscript-projects/<id>/`：只存放用户明确批准的目标论文，冻结正式分析、完成 Figure 和稿件。

Pilot 不能通过改一个状态直接冒充目标论文。晋升会建立独立 manuscript project，并记录来源 Pilot、路线和真实运行；正式项目必须重新冻结并运行完整分析。

支持代码位于 `src/paper2paper/`，只负责创建、校验、查看和晋升 workspace。它不是第三种研究层。

## 证据不能混写

必须区分：

- 结构与合同校验；
- 单元/合成测试；
- 真实数据运行；
- 科学审查；
- 用户对目标论文路线和最终发布的批准。

`paper2paper validate` 只能检查可机械判断的结构、引用和阶段依赖。它不能证明数据内容真实、文献检索完整、科学结论正确或论文值得发表；这些仍需结合原始文件、日志、source table、Figure 和人类审查。

## 仓库结构

```text
Paper2Paper/
├─ README.md、RULES.md
├─ src/paper2paper/       # 最小 workspace 工具
├─ tests/                 # 防止关键误判的必要测试
├─ pilots/                # 相互隔离的真实锚点文献 Pilot
├─ manuscript-projects/   # 用户批准后的目标论文
└─ docs/                  # 数据、代码、执行、审查和双循环说明
```

## Quick start

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
```

用户在 Pilot 中审查并批准一条 `manuscript_candidate` 后：

```bash
paper2paper promote pilots/my-project ROUTE-1 manuscript-projects/my-paper \
  --project-id P2P-MS-MY-PAPER \
  --title "My manuscript project"
```

当前真实案例和未决事项见 [PROJECT_STATE.md](PROJECT_STATE.md)。数据寻找见 [docs/data-discovery.md](docs/data-discovery.md)，代码补足见 [docs/code-recovery.md](docs/code-recovery.md)，完整操作见 [docs/workflow.md](docs/workflow.md)，用户审查边界见 [docs/review-guide.md](docs/review-guide.md)。

## 数据与安全边界

Git 只保存自有代码、配置、检索记录、manifest、经许可审查的 source table、汇总结果和稿件。不得提交密码、token、患者可识别信息、受控数据、未经许可的 PDF 或大型测序文件。外部数据、论文附件和第三方代码仍受各自条款约束；MIT License 只覆盖 Paper2Paper 自有代码。
