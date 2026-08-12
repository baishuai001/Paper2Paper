# CRC 第零阶段真实运行代码

这组脚本把
[`../../analysis/phase-zero-execution-specification.md`](../../analysis/phase-zero-execution-specification.md)
中的 P0.0–P0.6 合同落实到 CRC-atlas 真实 H5AD。实际结果与科学边界见
[`../../reports/phase-zero-real-run-report.md`](../../reports/phase-zero-real-run-report.md)。

代码分三类证据，不能混写：

1. 单元/反向测试：检查字段、身份连接、矩阵、路径和结果合同会失败关闭；
2. 真实数据诊断：确认真实文件、官方软件和 adapter 能运行并暴露限制；
3. 科学验收：还需患者/数据集留出、比较设计和预注册退出门，不能由“30 个测试通过”替代。

## 1. 环境

```bash
python3 -m venv tmp/hcc-sc-spatial-npj-2026/phase-zero/env/bootstrap
tmp/hcc-sc-spatial-npj-2026/phase-zero/env/bootstrap/bin/pip install \
  -r pilots/hcc-sc-spatial-npj-2026/code/phase_zero/requirements-bootstrap.txt

python3.10 -m venv tmp/hcc-sc-spatial-npj-2026/phase-zero/env/cnmf
tmp/hcc-sc-spatial-npj-2026/phase-zero/env/cnmf/bin/pip install \
  -r pilots/hcc-sc-spatial-npj-2026/code/phase_zero/requirements-cnmf.txt
```

CopyKAT 和 SCEVAN 使用独立的 R 安装记录。当前真实运行固定 CopyKAT 1.1.0；SCEVAN
1.0.3 来自提交 `5a49b88ac9445eeffcebb95404e3190992faac04`，其 R library 位于
`tmp/hcc-sc-spatial-npj-2026/phase-zero/env/R-scevan-lib`。安装目录和来源 checksum 不提交 Git。

## 2. P0.0–P0.4

以下命令依次核验 Table S1、H5AD、身份连接、发布对象 QC、发布 embedding 和 marker 证据：

```bash
PY=tmp/hcc-sc-spatial-npj-2026/phase-zero/env/bootstrap/bin/python
ROOT=tmp/hcc-sc-spatial-npj-2026/phase-zero
CODE=pilots/hcc-sc-spatial-npj-2026/code/phase_zero

$PY $CODE/build_crc_table_s1_design.py \
  --workbook $ROOT/raw/crc_atlas_table_s1.xlsx \
  --expected-sha256 EF58CFC7C00ED9826B8E421767A71D45271F4847155A5E5DDD1620089EE7D920 \
  --output-dir $ROOT/outputs/table_s1

$PY $CODE/inspect_crc_h5ad.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --expected-bytes 30875155333 \
  --output-dir $ROOT/outputs/h5ad

$PY $CODE/join_crc_h5ad_table_s1.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --table-samples $ROOT/outputs/table_s1/P0_table_s1_samples.tsv \
  --expected-obs 3790266 \
  --output-dir $ROOT/outputs/h5ad_table_join

$PY $CODE/audit_crc_atlas_qc.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --expected-obs 3790266 \
  --output-dir $ROOT/outputs/qc_audit

$PY $CODE/audit_crc_integration.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --expected-obs 3790266 \
  --output-dir $ROOT/outputs/integration_audit

$PY $CODE/audit_crc_markers.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --marker-table $ROOT/raw/colon_atlas_cell_type_marker_genes.csv \
  --expected-obs 3790266 \
  --output-dir $ROOT/outputs/marker_audit
```

P0.2–P0.4 读取的是发布后对象，所以只能审计可见证据；它们不会伪造已被删除的细胞、
ambient RNA、完整 doublet 流程、整合前 comparator 或重新训练结果。

## 3. P0.5：CNA 面板

### 3.1 冻结样本面板

选择器先统计同时拥有至少 50 个作者 Cancer 候选和 50 个明确正常参考的样本，再在
polyp、tumor、metastasis 中按数据集/患者覆盖确定性选择。它不读取基因、CNA 调用或临床
结局；H5AD checksum 继承自已经通过的 P0.0 收据，同时重新核对文件路径和字节数。

```bash
$PY $CODE/select_crc_cna_panel.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --h5ad-receipt $ROOT/outputs/h5ad/P0_h5ad_receipt.json \
  --output-dir $ROOT/outputs/cna_panel_selection_v2
```

### 3.2 顺序运行 CopyKAT 与 SCEVAN

```bash
$PY $CODE/run_crc_cna_panel.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --panel $ROOT/outputs/cna_panel_selection_v2/P0_cna_selected_panel.tsv \
  --output-root $ROOT/outputs/cna_panel_run_12 \
  --cna-runner $CODE/run_crc_cna.R \
  --scevan-r-lib $ROOT/env/R-scevan-lib \
  --cores 4 --seed 20260810
```

每个样本各自准备 raw-count 输入、运行两种方法并保存日志。一个样本失败不会被另一个样本
覆盖；最终退出码非零且 receipt 标记 `completed_with_failures`。

### 3.3 汇总但不制造“真值”

```bash
$PY $CODE/summarize_crc_cna_panel.py \
  --panel $ROOT/outputs/cna_panel_selection_v2/P0_cna_selected_panel.tsv \
  --run-root $ROOT/outputs/cna_panel_run_12 \
  --output-dir $ROOT/outputs/cna_panel_summary_12
```

CopyKAT=aneuploid 且 SCEVAN=tumor 的交集称为“两方法 CNA 支持候选”，不称 DNA 验证的
恶性真值。diploid/normal 也不自动等于非恶性；`not.defined/filtered` 保留 unresolved。

## 4. P0.6：CNA 支持输入的 cNMF 诊断

先按样本限额构建输入；这里要求三种横断面状态都有合格分析单元：

```bash
$PY $CODE/prepare_crc_cnmf_from_cna.py \
  --h5ad $ROOT/raw/crc_atlas_cellxgene.h5ad \
  --cna-evidence $ROOT/outputs/cna_panel_summary_12/P0_cna_panel_cell_evidence.tsv \
  --output-dir $ROOT/outputs/cnmf_cna_input \
  --min-cells-per-unit 30 --max-cells-per-unit 150 \
  --required-state polyp --required-state tumor --required-state metastasis \
  --seed 20260810
```

随后在独立 cNMF 1.7.1 环境运行两个 pooled seed，并把每个状态的独立样本交替分为患者/
数据集完全不重叠的 A、B 两折。下面仍是诊断强度；每个 K 只有 50 次初始化，不是合同要求
的最终 100–200 次：

```bash
CNMF_PY=$ROOT/env/cnmf/bin/python
for SEED in 20260810 20260811; do
  $CNMF_PY $CODE/run_crc_cnmf.py \
    --counts-h5ad $ROOT/outputs/cnmf_cna_input/crc_two_method_cna_supported_counts.h5ad \
    --output-dir $ROOT/outputs/cnmf_cna_runs/seed_${SEED} \
    --run-name crc_cna_cnmf \
    --components 5 6 7 8 9 10 --n-iter 50 --seed ${SEED}
done

$PY $CODE/split_crc_cnmf_holdouts.py \
  --input-h5ad $ROOT/outputs/cnmf_cna_input/crc_two_method_cna_supported_counts.h5ad \
  --output-dir $ROOT/outputs/cnmf_cna_holdouts --seed 20260810

for FOLD in A B; do
  $CNMF_PY $CODE/run_crc_cnmf.py \
    --counts-h5ad $ROOT/outputs/cnmf_cna_holdouts/crc_two_method_cna_supported_holdout_${FOLD}.h5ad \
    --output-dir $ROOT/outputs/cnmf_cna_runs/holdout_${FOLD} \
    --run-name crc_cna_cnmf \
    --components 5 6 7 8 9 10 --n-iter 50 --seed 20260810
done

$PY $CODE/evaluate_crc_cnmf_stability.py \
  --run seed_20260810=$ROOT/outputs/cnmf_cna_runs/seed_20260810/crc_cna_cnmf \
  --run seed_20260811=$ROOT/outputs/cnmf_cna_runs/seed_20260811/crc_cna_cnmf \
  --run holdout_A=$ROOT/outputs/cnmf_cna_runs/holdout_A/crc_cna_cnmf \
  --run holdout_B=$ROOT/outputs/cnmf_cna_runs/holdout_B/crc_cna_cnmf \
  --output-dir $ROOT/outputs/cnmf_cna_stability \
  --top-n 50 --permutations 100 --seed 20260810
```

程序匹配使用 cosine、top-gene Jaccard、Hungarian 一对一匹配和基因置换 99% 阈值。
即使跨 seed 稳定，若程序被单一数据集支配、患者/数据集留出不复现或 CNA 集合敏感性不
稳定，也不能冻结 K 或命名 meta-program。K 未冻结前，输出只按每个 K 记录 A/B 直接
一对一配对和逐 K 汇总；禁止把跨 seed、相邻 K 和不同输入的匹配边做无约束传递并集来
制造跨 K 程序家族，因为这会把同一运行、同一 K 的不同程序错误合并。

## 5. 回归测试

```bash
cd pilots/hcc-sc-spatial-npj-2026/code/phase_zero
$PAPER2PAPER_SERVER_ROOT/tmp/hcc-sc-spatial-npj-2026/phase-zero/env/bootstrap/bin/python \
  -m unittest test_phase_zero_bootstrap.py
```

测试覆盖资源改变、身份连接不完整、负值/非整数矩阵、重复单元、不完整 CNA 方法连接、
SCEVAN adapter、cNMF spectra 方向、CNA→cNMF 状态约束，以及跨 K 伪家族合并。测试数量
不是项目目标；只保留能防止真实执行中已出现或高风险误判的测试。

## 6. 2026-08-11 配对恶性识别入口

本轮在运行前冻结了 `analysis/paired-cna-review-contract.md`，随后用以下脚本完成 43 个样本的
逐样本 CopyKAT 1.1.0 与 SCEVAN 1.0.3 识别：

- `build_crc_paired_cna_panel.py`：从配对设计表生成冻结面板；
- `run_crc_cna_panel_configurable.py`：按指定分片运行初筛或全部作者 Cancer 候选；
- `build_full_cna_shards.py`：只根据初筛结果和冻结规则生成允许样本及四个平衡分片；
- `combine_crc_paired_cna_summaries.py`：合并分片并生成样本、患者、队列和正常参考审查表；
- `test_paired_cna_helpers.py`：防止面板身份、分片重复、合并重复行和监督入口再次退化。

本轮真实运行命令、绝对路径监督脚本和小型聚合结果保存在
`evidence/phase0-paired-cna-review/`。服务器监督脚本只是执行证据；移植到新机器时应重新配置
路径和环境，不能把其中的绝对路径当作产品默认值。正式 cNMF 在人工审查前没有启动。

## 7. 正式 Liu cNMF

用户批准后，正式运行使用以下新增入口：

- `formal_cnmf_design.py`：纯表格层的队列筛选、患者×状态平衡和患者留出规则；
- `prepare_crc_formal_cnmf_input.py`：合并四个完整 CNA 分片并从 `raw/X` 提取正式计数；
- `split_crc_formal_patient_holdouts.py`：保持同一患者原发与肝转移在同一折；
- `run_crc_cnmf.py --run-role ... --total-workers ...`：使用 cNMF 官方任务分片并区分诊断与正式候选收据；
- `evaluate_crc_cnmf_stability.py`：分别标记 NMF seed、细胞重抽样和患者留出，不能把 Liu 内部留出误写成数据集留出。

完整参数见 `analysis/formal-cnmf-execution-contract.md`。旧 `cnmf_cna_runs_v1` 是 1,561 个细胞、
K=5–10、每 K 50 次初始化的诊断运行，不能与 `formal_cnmf_v1` 混合或冒充正式结果。

## 8. K=10 Che 外部投射

K 与程序资格冻结后，使用以下入口完成 Che 外部复核：

- `freeze_liu_k10_projection_reference.py`：只读 Liu 主运行，冻结 P1–P10 的 top-50 rank-score 投射公式、Liu 方向和 P1/P7/P8 确认性资格；
- `run_che_external_projection.py`：盲法计算 Che 全 10 程序，按患者比较原发与肝转移，并运行两项预定敏感性；
- `test_che_projection_helpers.py`：检查投射基因覆盖、缺失基因位置、精确小样本检验、Hodges–Lehmann、Holm 校正和最终状态标签。

完整参数与结论边界见 `analysis/che-external-projection-contract.md`。真实结果位于
`outputs/phase-zero/che_projection_v1/`，执行日志与运行时副本位于
`evidence/phase0-che-projection/`。本轮没有在 Che 中重做 cNMF、重新选 K、改变资格或进行程序命名。
