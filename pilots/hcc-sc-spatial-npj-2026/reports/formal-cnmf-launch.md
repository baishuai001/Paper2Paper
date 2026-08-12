# CRC 正式 cNMF 启动记录

日期：2026-08-12  
状态：服务器后台运行中；尚无 K、程序或 Figure 结论

## 用户批准

用户已审查 `paired-malignancy-human-review.md`，同意按报告建议进入 cNMF。批准范围是 Liu
主发现及其内部稳定性，不包括程序命名、Che 外部投射或下游生物学解释。

## 冻结设计

- 23,766 个 Liu 双方法 CNA 支持细胞作为可用池；
- 11 位配对患者、25 个样本；
- 每位患者每个状态最多 500 个细胞，预计主输入 9,370 个细胞；
- K=5–20，每 K 100 次初始化；
- 两个 NMF seed、一个细胞重抽样、两个患者不重叠留出；
- cNMF 1.7.1，2,000 HVG，density threshold 0.1；
- 五组运行完成后停在程序命名前的 K 人工审查点。

完整合同见 [`../analysis/formal-cnmf-execution-contract.md`](../analysis/formal-cnmf-execution-contract.md)。

## 代码验收

- 新正式设计与运行测试：6/6 通过；
- 既有第零阶段兼容性测试：34/34 通过；
- 代码、合同和服务器监督脚本均绑定 SHA256；
- 并行入口使用 cNMF 官方 `worker_i/total_workers` 任务分片；
- 服务器绝对路径监督脚本作为执行证据保存，不视为跨机器默认入口。

## 服务器位置

- 监督进程首次 PID：`3593125`；
- 输出：`$PAPER2PAPER_SERVER_ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero/outputs/formal_cnmf_v1`；
- 日志：`$PAPER2PAPER_SERVER_ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero/logs/formal_cnmf_v1`；
- 顶层监督日志：`$PAPER2PAPER_SERVER_ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero/logs/formal_cnmf_v1.supervisor.log`；
- 冻结 manifest：`$PAPER2PAPER_SERVER_ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero/manifests/formal_cnmf_v1`；
- 监控心跳：`crc-cnmf-k`，每 15 分钟检查一次。

## 当前不能声称

启动成功只说明任务已进入冻结执行链。两套输入尚需生成并验收，五组 cNMF 尚需全部完成；
在 `READY_FOR_K_REVIEW` 出现前，不能声称 K 稳定、程序存在、原发—肝转移差异成立或 NPJ
HCC 的恶性程序 Figure 已在 CRC 中复现。
