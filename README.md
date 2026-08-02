# Paper2Paper

Paper2Paper 是一个面向初学者的、以执行为优先的文献模仿 workflow。它从一篇
锚点论文出发，允许优先考虑保留论文模板的改写，并把路线落实为可获取的数据、
可审查的代码、可生成的图表和最终稿件。

项目的终点是论文，不是无限扩建 workflow。具有明确科学问题且不构成实质性重复发表
时，以下路线均可接受：

- 原样复现，用于学习、核验或建立可靠基线；
- 替换 marker、单个基因或基因集；
- 替换核心细胞类型；
- 替换癌种；
- 扩展为泛癌；
- 构建转录组 signature；
- 组合上述替换；
- 在保持稿件主线的前提下修复锚点论文的科学或计算问题。

Paper2Paper 不要求候选路线通过“创新性”门槛，也不按创新程度排序。发表重合检查是
另一个问题：`adjacent` 和 `high_overlap_distinguishable` 可以继续，只有
`duplicate` 必须停止。

## 执行优先原则

候选路线先看能否做成，再看是否值得投入。固定优先顺序为：

1. 代码是否达到可复用资格；
2. 数据及关键 metadata 是否已经验证可取得；
3. 是否能覆盖最低充分的正文图和稿件证据链；
4. 研究设计是否科学有效；
5. 能否复用锚点论文的结构与模块；
6. 初学者负担与完成时间；
7. 是否存在实质性重复发表。

路线使用透明分级而非不透明总分：

- `P0`：科学设计合格，数据已验证，代码已达标，图件映射完整，且非重复发表；
- `P1`：数据已验证，代码 donor 可小幅适配，仍需完成一次执行 spike；
- `P2`：数据基本可用，但关键代码需要重建；
- `P3`：关键数据、代码、科学设计、图件闭环或发表重合尚未过关。

只有 `P0` 路线能够被批准为活动稿件路线。

## Quick start

Requires Python 3.10 or later and has no runtime dependencies outside the
standard library.

```bash
python -m pip install -e .

paper2paper validate workspaces/spp1-tam-jitc
paper2paper status workspaces/spp1-tam-jitc
paper2paper priorities workspaces/spp1-tam-jitc
```

Create a new intake workspace:

```bash
paper2paper init ../my-paper-project \
  --project-id P2P-MY-PAPER \
  --title "My paper adaptation"
```

## 仓库结构

- `src/paper2paper/`：schema、验证、状态和依赖影响工具；
- `docs/`：路线选择、数据、代码、发表重合、图件规划和反馈规则；
- `tests/`：产品纪律和通用工程能力的回归测试；
- `workspaces/`：彼此隔离的论文实例；
- `templates/`：常见模仿路线的最小模板。

机器可读的事实和决定是来源真相：`PROJECT.json`、`registry/*.tsv`、运行清单和
source table。Markdown 报告只是它们的可读视图。只存在于聊天中的决定不算项目决定。

## 数据与凭据边界

仓库可以保存 schema、代码、配置、较小的 source table、校验和、数据 manifest、
审查记录和报告。不要提交服务器密码、token、患者可识别信息、受控数据、未经许可的
PDF，或大型 FASTQ/H5AD/RDS/图像文件。

更多约束见 [项目纪律](PROJECT_DISCIPLINE.md)、[workflow](docs/workflow.md) 和
[从 PaperRoute 迁移说明](MIGRATION_FROM_PAPERROUTE.md)。

## License

尚未选择许可证。在许可证加入前，不要假定代码可被外部分发或再许可。
