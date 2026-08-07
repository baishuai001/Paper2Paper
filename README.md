# Paper2Paper

Paper2Paper 是一个面向初学者的论文模仿 workflow。它从一篇锚点论文出发，帮助用户完成：

1. 拆解锚点论文；
2. 生成可比较的改写路线；
3. 系统寻找并验证数据；
4. 系统寻找、组合或重建代码；
5. 用最小真实数据完成执行试验；
6. 选择一条可完成的论文路线；
7. 冻结分析设计并生成结果；
8. 根据结果调整主张；
9. 形成论文稿件和可复现包。

终点是论文，不是方向推荐报告，也不是无限扩建 workflow。

## 可以怎样模仿

以下改写都可以进入候选组合：

- 原样复现；
- 替换 marker 或单个基因；
- 替换基因集；
- 替换核心细胞类型；
- 替换癌种；
- 扩展为泛癌；
- 构建转录组 signature；
- 组合多个明确替换；
- 保留锚点论文的科学和图件结构，采用其他论文或官方方法的可靠代码；
- 修复锚点论文中的数据、统计或计算问题。

Paper2Paper 不按“改动有多大”排序，也不要求候选路线通过创新性门槛。简单替换可以接受；
科学设计失效、数据不可得、代码无法可靠运行或与已发表论文实质重复则不能继续。

`routes.tsv.route_role`明确区分三类用途：`manuscript_candidate`是接受投稿资格审查的
稿件候选，`training`是限定在同类分析模块中的学习或重建实例，`supporting`只支撑其他
稿件路线。执行完成与稿件资格分别报告：训练实例可以完整跑通，但不能因此被当作整个
Paper2Paper的验收基准，也不能被误选为最终稿件主线。

## 新核心：执行和审查各占一半

Paper2Paper 不只保存结论。每条路线必须留下可检查的执行证据：

- 数据需求、检索式、候选数据、字段级资源、代表性文件解析和队列用途；
- 代码模块需求、donor检索、版本、许可、入口、目标路径环境、安装和smoke test；
- signature的特征顺序、系数、映射、归一化、缺失策略和cutoff合同；
- 运行中发现的问题、证据、影响、是否推广为通用约束及处理状态；
- Figure与数据需求、代码模块、source table和验收测试的映射；
- 最近发表检索和可区分点；
- 用户的路线选择和高影响变更决定；
- 运行命令、环境、数据manifest、结果和稿件影响。

结构化记录用于让用户审查AI的工作，不用于代替实际下载、运行和科学判断。
`paper2paper validate`只检查结构和已登记合同之间是否一致；它不会自动重新下载数据、
复算checksum、执行分析、核对许可证含义或证明科学结论正确。

## 工作阶段

```text
anchor_audit
  -> route_generation
  -> verification
  -> selection
  -> specification
  -> execution
  -> interpretation
  -> writing
  -> complete
```

阶段含义见 [完整 workflow](docs/workflow.md)。

## Quick start

需要 Python 3.10 或更高版本，无第三方运行依赖。

```bash
python -m pip install -e .

paper2paper init ../my-project \
  --project-id P2P-MY-PROJECT \
  --title "My paper project" \
  --anchor-title "Anchor paper title" \
  --doi "10.xxxx/xxxxx"

paper2paper validate ../my-project
paper2paper status ../my-project
paper2paper next ../my-project
paper2paper report ../my-project
```

`next` 根据当前事实列出下一步，而不是生成脱离数据和代码证据的方向结论。

## 仓库结构

- `docs/`：从锚点审计到稿件交付的操作手册；
- `src/paper2paper/`：workspace初始化、校验、路线缺口和下一步工具；
- `tests/`：关键科学与执行纪律的回归测试；
- `pilots/`：彼此隔离的真实锚点论文测试实例。

每个论文workspace由少量文件组成：

```text
PROJECT.json
anchor/audit.md
evidence/routes.tsv
evidence/search_log.tsv
evidence/data_requirements.tsv
evidence/data_candidates.tsv
evidence/data_resources.tsv
evidence/cohort_usage.tsv
evidence/code_requirements.tsv
evidence/code_candidates.tsv
evidence/model_specifications.tsv
evidence/figures.tsv
evidence/literature.tsv
evidence/decisions.tsv
evidence/issues.tsv
execution/runs.tsv
execution/results.tsv
analysis/specification.md
manuscript/draft.md
reports/readiness.md
```

## 数据和凭据边界

Git只保存代码、配置、较小的metadata、检索记录、manifest、经过再分发审查的source table、报告和稿件。
不要提交服务器密码、token、患者可识别信息、受控数据、未经许可的PDF或大型测序文件。
逐患者派生表即使只含公开编号，也必须先核对原始数据的再分发和署名条件；未完成时只
提交生成代码、汇总结果和不含逐患者记录的审计产物。

项目的绑定规则见 [RULES.md](RULES.md)。

## License

Paper2Paper自有代码采用 [MIT License](LICENSE)。外部数据、论文附件和第三方代码仍受
各自来源的许可证及使用条款约束，MIT许可证不会改变这些外部材料的权利状态。
