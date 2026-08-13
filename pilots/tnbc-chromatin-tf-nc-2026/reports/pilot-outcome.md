# Pilot outcome

当前有效的一级闸门是 2026-08-13 完成的 CRC M-vs-rest 调控程序闸门，判为 **FAIL 并停止**。完整证据见 [`regulatory-gate-m-vs-rest-final-report.md`](regulatory-gate-m-vs-rest-final-report.md)。2026-08-12 的四分类聚类 FAIL 被保留，但仅是辅助 QC，不再作为 TNBC 调控框架的必要入口。

## Paper-side outcome

- 数据、CRC 特异 ARACNe3 网络、VIPER/msVIPER 与跨研究统计均可完整执行；结果不是技术性不可判定。
- 0 个 TF 通过冻结的跨研究可重复标准，0 个 TF 获 msVIPER FDR≤0.01 确认；20/100/200-cell 敏感性均为 0。
- LOSO AUROC 为 0.651、置换 p=0.01996，但 CI 下限只有 0.528；它是边界性辅助信号，不能救回硬 TF 程序。
- 本轮不支持继续做 M-vs-rest 的染色质、治疗或药敏迁移；也不支持宣称 M 型或 nominal TF 没有生物学意义。

## Workflow-side outcome

1. 锚点方法必须在冻结前审计作者代码，而不只看 Methods。作者实际使用 1 个 ARACNe3 子网络、`subnet1`、`minsize=1` 和 1,000 次 null；早期添加的 100 子网、第二次 BH 和较大 minsize 均被撤回。
2. “合理的附加 QC”不能未经授权变成硬停止规则。四类聚类重建和 CollecTRI＋ULM 只能作辅助或敏感性分析，不能替代 ARACNe3＋VIPER。
3. 技术数据集不等于独立研究。fresh/frozen、10x v2/v3 必须归回 `study_id`；本次真实独立证据只有 3 个研究。
4. 失败规则要直接编码。硬 verdict 只由可重复 TF 程序决定；LOSO 始终报告，但不覆盖该规则。
5. 跨平台 shell 入口应显式调用解释器。Windows/Git 同步丢失执行位导致一次无科学结果的失败，随后以 `bash run_aracne3.sh` 修复并保留错误日志。
6. 收据要区分旧结果与本次结果。本次以生成时间、所需字段和源文件 SHA256 识别旧收据，并为非数据重试记录了全哈希收据 SHA256、H5AD mtime 和复用理由。
7. 正式结果后应使用独立算术路径复核 BH、联合条件、AUROC、bootstrap 和 permutation；本次 18 项复核全部通过。

## Reuse boundary

本 Pilot 实际覆盖：CRC-atlas H5AD 数据审计、患者 Cancer-cell pseudobulk、GDC TCGA-COAD/READ TPM、PAN-GO 审计、CRC ARACNe3 网络、VIPER/msVIPER、独立研究 meta、LOSO、敏感性分析、自动停止判定和结果复核。它没有覆盖 ATAC、治疗响应、药物选择、空间验证或湿实验。

表型切换已 manifest 化；TCGA/PAN-GO/ARACNe3 可作为共享网络工件复用。新亚型必须使用新的 work 目录和预先冻结协议，不能在本次失败结果上调阈值。

## Human decision and next boundary

本轮到此关闭，不自动尝试其他 CRC 亚型。若要继续，需由人类明确选择下一个临床或生物学亚型，再冻结新的 manifest 和停止规则。
