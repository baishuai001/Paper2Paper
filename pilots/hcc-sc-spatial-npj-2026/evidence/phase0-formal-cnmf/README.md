# Phase 0 正式 cNMF 精简证据包

本目录保存可进入 Git 审查的精简证据，不复制服务器上的全部 NMF 分片、表达矩阵、usage 和
spectra 大文件。

- `run-receipts/`：五组正式运行的输入、参数、版本、产物 checksum 和结论边界；
- `stability/`：稳定性运行收据、置换阈值和患者 A/B 留出逐 K 汇总；
- `review/`：供人工审查的逐 K 综合表及其来源收据；
- `code_and_contract.sha256`：启动时冻结的合同和代码 hash；
- `server-bound-wrappers/`：实际服务器监督入口。
- `runtime-snapshots/`：正式 K=5–20 运行实际使用的代码快照；
- `sensitivity/`：两项预定敏感性分析的输入、运行、K 比较与稳定性精简回执。

完整产物位于服务器：
`$PAPER2PAPER_SERVER_ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero/outputs/`。

版本关系不能混写：五组正式 K=5–20 运行实际使用的是
`runtime-snapshots/run_crc_cnmf.formal-E8AC3E755886.py`（SHA256 `E8AC…`）；当前可复用代码及
两项敏感性运行使用的是增强版 `code/phase_zero/run_crc_cnmf.py`（SHA256 `628A…`）。正式运行
回执保留前者，敏感性 `code_and_contract.sha256` 保留后者。

本目录只保存可审查的精简证据，不单独代表当前全部科学状态。K=10 的冻结依据及后续 Che
外部复核分别见 `reports/formal-cnmf-sensitivity-human-review.md`、
`reports/k10-program-che-eligibility-human-review.md` 和
`reports/che-external-projection-human-review.md`。富集、空间、临床和扰动分析仍未因此自动完成。
