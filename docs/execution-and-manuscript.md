# 最小真实运行、Pilot 结论与正式稿件

## 最小真实运行不是“看起来能跑”

路线只有在以下证据同时存在时，才可记录为
`evidence_stage=minimal_real_run`：

1. 必需数据字段已经落到真实文件，并完成必要的 ID、样本数和队列用途核查；
2. 必需代码模块在目标环境运行了代表性真实输入；
3. `runs.tsv` 中存在 `run_kind=real_data`、`status=passed`、`exit_code=0`
   的运行，保留输入来源、命令、环境、data manifest、日志和全部声明产物；
4. `results.tsv` 至少有一项结果连接到该 passed run 和 Figure/source table；
5. 最低主图已完成映射，并至少生成一个代表性 source table 或 Figure 片段；
6. 失败、阴性结果和限制没有被后一次成功覆盖。

`unit`、`synthetic_integration` 和 `environment_smoke` 运行可以验证公式、失败条件
和可移植性，但永远不能替代上述真实数据运行。任一声明日志、manifest、artifact 或
result source table 缺失、指向目录或只是零字节占位文件时，`execution_ready` 必须失败。

当前 alpha 版仍有一条明确的人审边界：`run_kind=real_data` 是执行者作出的声明。
校验器会检查 manifest、日志和产物确实是非空文件，并检查运行与结果的引用关系；但它还不能
仅凭这些字段证明文件内容一定来自登记的真实队列，也不能识别故意把合成数据错标为
`real_data`。审查者必须把 data manifest 与 `data_candidates.tsv`、`data_resources.tsv`、
`cohort_usage.tsv` 逐项核对。这个边界不能用文件名关键词检查来假装解决。

这一级证据回答“这个限定路线是否真实执行过”。它不回答“是否值得形成新论文”。

## 两个不同的判断

- `execution_ready`：数据、代码、Figure 映射、真实运行和结果是否满足当前路线的最小执行合同。
- `promotion_evidence_complete`：路线是否为 `manuscript_candidate`，且执行证据、最近发表审查和用户的 `approve_route_evidence` 决定是否齐全。

训练复现和辅助路线可以 `execution_ready=true`，但不能因此晋升为正式论文项目。即使
`promotion_evidence_complete=true`，也只是说明证据包可以进入用户决策；它不等于用户已经批准，
更不等于保证能发表。

## Pilot 结论必须分成两侧

`reports/pilot-outcome.md` 同时交付：

- 文献侧：锚点论文审计、候选路线、真实运行、科学限制和 continue/refine/reroute/stop 判断；
- workflow 侧：发现的问题、当前论文修复、Paper2Paper 是否需要修改，以及修改后怎样用同一 Pilot 复验。

问题记录来源范围（anchor、pilot execution、both 或 scientific result）、证据、后果、
当前修复、workflow action 和状态。科学阴性结果不能被误写成 workflow 故障；一个 Pilot
的局部修复也不能被写成所有论文都已验证。

## 问题状态与 workflow 修改分开

`issues.tsv.status` 说明当前问题实例是否 open、resolved、accepted risk、superseded 或
wont fix；`workflow_action` 只说明 Paper2Paper 修改是 none、consider、implemented 或
deferred。`implemented` 仍只证明这次修改已经落地并用当前案例复验，不证明对其他论文
天然适用。以后在另一文献中使用时，重新核对输入、方法和科学边界即可，不建立中央成熟度体系。

## 从 Pilot 到正式稿件

Pilot 不包含正式 `manuscript/draft.md`。只有下列条件满足后，才可运行
`paper2paper promote` 建立独立 `manuscript-project`：

1. route role 是 `manuscript_candidate`；
2. 最小真实运行和最近发表审查完成；
3. 用户记录 `approve_route_evidence`；
4. 用户另行记录 `promote`；
5. 新项目在 `PROJECT.json.provenance` 锁定来源 Pilot、route 和合格 passed
   real-data runs；单元测试 run 不冒充科学来源运行。

`approve_route_evidence` 和 `promote` 都必须由用户在 `pilot_review` 记录带时区的完整时间戳。
前者不得早于最后一次合格真实数据运行；后者不得早于该运行，也不得早于证据批准。历史批准
不能在 AI 后续补齐数据或结果后被自动消费。

CLI 只允许把新项目创建为当前仓库 `manuscript-projects/` 的一级子目录；目标目录必须尚未
存在，`project_id` 必须在 Pilot 和正式项目中全局唯一。

晋升只复制可审查的方向与证据，不复制 Pilot 的运行和结果冒充正式分析。正式项目必须重新冻结
`analysis/specification.md`、重新验证代码和环境、重新生成 source tables/Figures，之后才进入解释和写作。
Pilot 不允许记录 `approve_release`，晋升代码也会防御性地拒绝复制这种旧记录。正式项目只有
在重新执行和审查后，才能由用户在 `stage=complete` 记录新的 `approve_release`；决定必须使用
带时区的完整时间戳，并且不得早于最后一次合格真实数据运行的完成时间。即使有这项决定，
只要代码 donor、真实运行、必需 Figure、source table 或 verified result 尚未闭环，
`stage=complete` 的校验仍必须失败。

## 正式项目的停止条件

当用户批准的路线已经具备真实数据、可靠代码、Figure 闭环、合理 claim 和可审查稿件时，停止为
“以后也许有用”的通用化继续扩建 Paper2Paper。只有新的真实阻断问题才允许恢复 workflow 修改。
