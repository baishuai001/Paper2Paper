# 跨 Pilot 学习账

这个目录让不同文献、不同对话中真实出现的问题能够汇集，但不把每个课题困难自动升级为 Paper2Paper 产品功能。

## 两层事实

1. `pilots/<id>/evidence/issues.tsv` 是现场记录：具体观察、证据、论文后果、当前课题修复和 workflow 建议。不得为了通用化而改写原始事实。
2. `learning/findings.tsv` 是归类记录：产品级表述、适用范围、处理状态、来源问题和复验引用。

`reports/workflow-learning.md` 是派生报告，可以用下面的命令重建：

```bash
paper2paper learn . --write-report
```

## 范围

- `paper_specific`：只在来源课题解决；
- `paper_type`：可能适用于同类论文，需真实复验；
- `general`：不同论文类型都可能产生相同错误。

## 状态

- `observed`：发现了可能的产品问题，尚未接受；
- `accepted`：已有依据将其作为产品缺口处理；
- `implemented`：最小修复已落地，尚未完成来源复验；
- `verified`：来源 Pilot 已复验并留下引用；
- `local_only`：明确只在当前课题解决；
- `rejected`：证据不足、已有能力已覆盖或补法会造成不必要扩建。

## 每个新对话的使用规则

1. 从最新已合并分支开始，读取 `RULES.md`、`PROJECT_STATE.md`、本文件和跨 Pilot 报告；
2. 在自己的 Pilot 中记录原始问题并先解决当前论文；
3. 运行聚合命令，AI 对尚未归类候选给出范围判断和依据；
4. 产品修复必须引用来源问题，且只做阻止该错误的最小修改；
5. 用来源 Pilot 复验后才能标记 `verified`；
6. 分支未合并前，其他对话看不到其中的记录。并行对话应优先只改各自 Pilot，核心与学习账在同步最新变更后再改；
7. 完成必要修复后立刻回到目标论文主循环。

这张账不是 PaperRoute 式创新性门槛，也不是能力成熟度、模块发布或路线评分系统。
