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

## 新核心：执行和审查各占一半

Paper2Paper 不只保存结论。每条路线必须留下可检查的执行证据：

- 数据需求、检索式、候选数据、代表性文件解析和队列独立性；
- 代码模块需求、donor检索、版本、许可、入口、安装和smoke test；
- Figure与数据需求、代码模块、source table和验收测试的映射；
- 最近发表检索和可区分点；
- 用户的路线选择和高影响变更决定；
- 运行命令、环境、数据manifest、结果和稿件影响。

结构化记录用于让用户审查AI的工作，不用于代替实际下载、运行和科学判断。

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
evidence/code_requirements.tsv
evidence/code_candidates.tsv
evidence/figures.tsv
evidence/literature.tsv
evidence/decisions.tsv
execution/runs.tsv
execution/results.tsv
analysis/specification.md
manuscript/draft.md
reports/readiness.md
```

## 数据和凭据边界

Git只保存代码、配置、较小的metadata、检索记录、manifest、source table、报告和稿件。
不要提交服务器密码、token、患者可识别信息、受控数据、未经许可的PDF或大型测序文件。

项目的绑定规则见 [RULES.md](RULES.md)。

## License

尚未选择许可证。在许可证加入前，不要假定代码可以被外部分发或再许可。
