# 一级闸门最终报告：CRC 免疫表型能否支撑 TNBC 式 TF 框架

## 结论

**一级闸门不通过（FAIL），并已按预设规则停止。**

CRC-atlas 数据具备足够的患者数、癌细胞数和跨数据集结构，且 M 型相对非 M 型呈现非常强、可跨数据集识别的髓系富集信号；但是，盲于作者标签的四分类重建未达到预先冻结的稳定性标准。因此不能把作者发布的四种免疫表型直接当作足够稳健的暴露变量，继续检验 TNBC 式“肿瘤细胞 TF 程序”。

这不是“CRC 中没有 M 型信号”，也不是“TNBC 框架永远不能用于 CRC”。准确解释是：**本次预先指定的四表型入口变量未通过严格的独立重建标准，当前迁移路线必须停在这里。**

## 冻结问题和停止规则

冻结问题是：作者的 M 型（macrophage/monocyte dominance）相对非 M 型，是否在原发瘤作者标注 Cancer cell 中对应稳定、患者级、跨数据集可复现的 TF 活性程序。

进入 TF 层之前，必须同时满足：

1. 数据身份、唯一标签、患者数和跨数据集最低条件；
2. M 型髓系比例差异 `p < 0.01` 且 Cliff's delta ≥ 0.33；
3. 相关距离四类重建 ARI ≥ 0.50；
4. 50 个固定种子 K-means 的中位 ARI ≥ 0.40；
5. leave-one-dataset-out M-vs-rest AUROC ≥ 0.70。

数据充分但任一科学信号条件失败，结论为 FAIL；随即停止癌细胞 pseudobulk、TF、TCGA-context、ATAC、药物和空间分析。完整冻结协议见 `analysis/gate1-protocol.md`。

## 云端数据审计

所有真实数据读取与计算均在云服务器 `/media/desk16/iy13202/projects/Paper2Paper` 完成。本地只保存代码、精简收据和报告。

- 输入 H5AD：30,875,155,333 bytes；既有 SHA256 收据为 `718774363933B57C6A4661E8AAF0EE7D86B717B047442AFE6457C01986789AC6`。
- 全库：3,790,266 cells、588 patients、1,553 samples、71 datasets。
- 原文限定队列：原发瘤、未治疗、`enrichment_cell_types=naive`；845,306 cells，213 名具作者免疫标签的患者。
- 主癌细胞阈值 ≥50：206 patients；M=46，non-M=160。
- 同时具有 ≥3 M 和 ≥3 non-M 的信息性数据集：4（最低要求 3）。
- 单一数据集贡献的最大 M 比例：21.7%（上限 60%）。
- 合规队列内矛盾标签：0。

初版审计错误地跨正常/转移及治疗前后样本检查标签，在 6 名患者中发现多标签并判失败。逐样本追踪证实这 6 例均为治疗状态不同的样本标签变化；限定到原文分析范围后全部消失。旧输出和失败日志没有删除；修正后的正式审计独立写入 `01_input_audit_v2_retry`。

## 表型统计验证

分析单位始终是患者。使用原发瘤、未治疗、非富集数据的作者细胞类型注释；按原文排除 neutrophil 特征，保留全肿瘤细胞作为分母；用线性模型去除 dataset 效应。聚类过程不读取作者四分类标签，标签只用于事后计算 ARI。

| 预设指标 | 实际结果 | 阈值 | 判定 |
|---|---:|---:|---|
| M vs non-M 髓系比例，中位数 | 20.64% vs 6.48% | M 更高 | 通过 |
| Mann–Whitney 双侧 P | 1.53 × 10^-20 | <0.01 | 通过 |
| Cliff's delta | 0.888 | ≥0.33 | 通过 |
| 相关距离四类 ARI | 0.173 | ≥0.50 | **失败** |
| 50 次 K-means 中位 ARI | 0.232 | ≥0.40 | **失败** |
| 50 次 K-means ARI 范围 | 0.197–0.241 | — | 稳定偏低 |
| LODO M-vs-rest AUROC | 0.905 | ≥0.70 | 通过 |

结果表明二分类 M 信号很强，但四类分区对合理的重建算法不稳健。M-vs-rest 的高 AUROC 不能抵消两个预设四分类失败项，因为协议要求全部条件满足。

## 代码支持、补齐与验证

云端固定了 CRC 作者仓库提交 `82c15ecc2ea36baa0c4cffe8818bf58e21dcfc48`。作者仓库公开了 atlas 处理、myeloid cNMF、空间及 TCGA 生存脚本，但未提供论文方法中患者四分类的可执行实现；TCGA 脚本还保留作者私有绝对路径。因此本次补齐了：

- 范围化数据审计与患者队列合同；
- 独立四表型重建、Cliff's delta、ARI、50-seed 稳健性和 LODO；
- raw/X 癌细胞患者 pseudobulk；
- OmniPath CollecTRI 下载冻结与 ULM TF 活性；
- dataset 固定效应模型、DerSimonian–Laird meta、无泄漏 LODO、bootstrap 和组内置换；
- fail-closed 自动判定器。

云端语法检查和 6 项合成单元测试全部通过。后四项 TF 模块已构建、可在一个未来通过表型闸门的路线中使用，但本次未读取真实表达矩阵运行，因为停止条件已经触发；不能把“未运行”写成“无 TF 信号”。

## 自动判定与复核材料

自动判定器在 Git 提交 `eeac0e70f321fdd17746bcf5dff048cb045de08c` 上返回：

- `verdict=FAIL`
- `stop_now=true`
- 失败项：`hierarchical_ARI_ge_0_50`、`median_kmeans_ARI_ge_0_40`
- 下游状态：`NOT_RUN_DUE_TO_PRE_SPECIFIED_STOP`

云端关键收据：

- `tmp/tnbc-chromatin-tf-nc-2026/outputs/00_decision/G0_gate1_decision.json`
- `tmp/tnbc-chromatin-tf-nc-2026/outputs/01_input_audit_v2_retry/G1_input_audit_receipt.json`
- `tmp/tnbc-chromatin-tf-nc-2026/outputs/02_phenotype/G2_phenotype_gate_receipt.json`
- `tmp/tnbc-chromatin-tf-nc-2026/logs/03_gate_code_tests.log`
- `tmp/tnbc-chromatin-tf-nc-2026/logs/04_phenotype_gate.log`
- `tmp/tnbc-chromatin-tf-nc-2026/logs/06_input_audit_v2_retry.log`
- `tmp/tnbc-chromatin-tf-nc-2026/logs/07_gate_decision.log`

正式审计和表型收据 SHA256 分别为 `DE0867762036929A627B75172F4BFFB203BB97FC78D325425335130655F5DACA` 和 `824C900A47F9237A07B4F4788BDCD13DD6DF1202342AFA1E8E3F5DC9815C099E`。

## 最终边界

本报告支持的最高结论是：CRC-atlas 的作者 M 标签具有强烈、跨数据集可识别的髓系组成含义，但作者四种免疫表型未通过本次预设的独立四类重建标准，因此不能以该四分类为入口继续声称或搜索 TNBC 式肿瘤细胞 TF 程序。

按用户要求，判断后停止；不做结果驱动的阈值放宽、聚类调参或下游探索。

