#!/usr/bin/env python3
"""Build a fail-closed K=10 program eligibility screen for Che validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


ANCHOR_RUN = "seed_20260812"
OTHER_SEED_RUN = "seed_20260813"
RESAMPLE_RUN = "cell_resample_20260813"
HOLDOUT_A_RUN = "patient_holdout_A"
HOLDOUT_B_RUN = "patient_holdout_B"
K = 10
PROGRAMS = list(range(1, K + 1))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_bool(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes"}


def one_to_one_block(
    frame: pd.DataFrame,
    *,
    comparison_type: str,
    left_run: str,
    right_run: str,
    left_k: int = K,
    right_k: int = K,
    expected_rows: int = K,
) -> pd.DataFrame:
    result = frame.loc[
        (frame["comparison_type"] == comparison_type)
        & (frame["left_run"] == left_run)
        & (frame["right_run"] == right_run)
        & (frame["left_k"].astype(int) == left_k)
        & (frame["right_k"].astype(int) == right_k)
    ].copy()
    if len(result) != expected_rows:
        raise ValueError(
            f"Expected {expected_rows} rows for {comparison_type} {left_run} "
            f"K{left_k} -> {right_run} K{right_k}; found {len(result)}"
        )
    if result["left_program"].duplicated().any() or result["right_program"].duplicated().any():
        raise ValueError("Hungarian match is not one-to-one")
    result["left_program"] = result["left_program"].astype(int)
    result["right_program"] = result["right_program"].astype(int)
    result["exceeds_both_nulls"] = result["exceeds_both_nulls"].map(as_bool)
    return result.sort_values("left_program")


def map_by_left(block: pd.DataFrame) -> dict[int, pd.Series]:
    return {int(row.left_program): row for row in block.itertuples(index=False)}


def evidence_row(
    component: str,
    anchor_program: int,
    row: object | None,
    *,
    required: bool,
    source_file: Path,
    expected_right_program: int | None = None,
) -> dict[str, object]:
    if row is None:
        return {
            "component": component,
            "eligibility_required": required,
            "anchor_program": anchor_program,
            "left_run": "",
            "left_k": np.nan,
            "left_program": np.nan,
            "right_run": "",
            "right_k": np.nan,
            "right_program": np.nan,
            "expected_right_program": expected_right_program,
            "connects_expected_program": False,
            "cosine": np.nan,
            "top_gene_jaccard": np.nan,
            "cosine_null_q99": np.nan,
            "jaccard_null_q99": np.nan,
            "exceeds_both_nulls": False,
            "component_pass": False,
            "source_file": str(source_file),
        }
    actual_right = int(row.right_program)
    connects = expected_right_program is None or actual_right == int(expected_right_program)
    return {
        "component": component,
        "eligibility_required": required,
        "anchor_program": anchor_program,
        "left_run": row.left_run,
        "left_k": int(row.left_k),
        "left_program": int(row.left_program),
        "right_run": row.right_run,
        "right_k": int(row.right_k),
        "right_program": actual_right,
        "expected_right_program": expected_right_program,
        "connects_expected_program": connects,
        "cosine": float(row.cosine),
        "top_gene_jaccard": float(row.top_gene_jaccard),
        "cosine_null_q99": float(row.cosine_null_q99),
        "jaccard_null_q99": float(row.jaccard_null_q99),
        "exceeds_both_nulls": bool(row.exceeds_both_nulls),
        "component_pass": bool(row.exceeds_both_nulls) and connects,
        "source_file": str(source_file),
    }


def load_gene_sets(path: Path) -> dict[str, set[str]]:
    table = pd.read_csv(path, sep="\t", dtype=str)
    required = {"category", "gene", "use", "provenance"}
    if set(table.columns) != required:
        raise ValueError(f"Gene-set columns must be exactly {sorted(required)}")
    if table[["category", "gene"]].duplicated().any():
        raise ValueError("Duplicate category/gene entries in frozen gene sets")
    table["gene"] = table["gene"].str.upper()
    return {key: set(group["gene"]) for key, group in table.groupby("category")}


def gene_categories(gene: str, gene_sets: dict[str, set[str]]) -> list[str]:
    symbol = gene.upper()
    categories = [category for category, genes in gene_sets.items() if symbol in genes]
    if symbol.startswith("MT-"):
        categories.append("mitochondrial")
    if re.fullmatch(r"RP[SL](?:\d+|P[0-2])", symbol):
        categories.append("ribosomal")
    if re.fullmatch(r"HB(?:A[12]|B|D|E1|G[12]|M|Q1|Z)", symbol):
        categories.append("hemoglobin")
    return sorted(set(categories))


def technical_decision(counts: dict[str, int], correlations: dict[str, float]) -> tuple[bool, list[str]]:
    combined_mito_ribo = counts.get("mitochondrial", 0) + counts.get("ribosomal", 0)
    technical_top50 = (
        counts.get("mitochondrial", 0)
        + counts.get("ribosomal", 0)
        + counts.get("stress_dissociation", 0)
    )
    reasons: list[str] = []
    checks = [
        (counts.get("mitochondrial", 0) >= 10, "mitochondrial_top50_ge_10"),
        (counts.get("ribosomal", 0) >= 15, "ribosomal_top50_ge_15"),
        (combined_mito_ribo >= 20, "mitochondrial_plus_ribosomal_top50_ge_20"),
        (counts.get("stress_dissociation", 0) >= 10, "stress_dissociation_top50_ge_10"),
        (counts.get("hemoglobin", 0) >= 3, "hemoglobin_top50_ge_3"),
        (
            abs(correlations["rho_mito_fraction"]) >= 0.50
            and counts.get("mitochondrial", 0) >= 5,
            "mitochondrial_usage_correlation",
        ),
        (
            max(abs(correlations["rho_log1p_total_counts"]), abs(correlations["rho_detected_genes"]))
            >= 0.60
            and technical_top50 >= 5,
            "library_complexity_usage_correlation",
        ),
    ]
    reasons.extend(label for passed, label in checks if passed)
    return bool(reasons), reasons


def contribution_row(path: Path, program: int) -> pd.Series:
    table = pd.read_csv(path, sep="\t")
    hit = table.loc[
        (table["k"].astype(int) == K)
        & (table["program"].astype(int) == int(program))
        & (table["unit"] == "analysis_patient_id")
    ]
    if len(hit) != 1:
        raise ValueError(f"Expected one patient contribution row for program {program} in {path}")
    return hit.iloc[0]


def compute_qc_correlations(
    h5ad_path: Path, usages_path: Path
) -> tuple[dict[int, dict[str, float]], dict[str, object]]:
    data = ad.read_h5ad(h5ad_path)
    usage = pd.read_csv(usages_path, sep="\t", index_col=0)
    usage.columns = [int(value) for value in usage.columns]
    if set(data.obs_names) != set(usage.index):
        raise ValueError("Usage cell IDs and H5AD obs_names do not match")
    usage = usage.loc[data.obs_names, PROGRAMS].astype(float)
    row_sum = usage.sum(axis=1)
    if (row_sum <= 0).any():
        raise ValueError("At least one cell has non-positive total cNMF usage")
    usage = usage.div(row_sum, axis=0)

    matrix = data.X
    total = np.asarray(matrix.sum(axis=1)).ravel()
    if sparse.issparse(matrix):
        detected = np.asarray(matrix.getnnz(axis=1)).ravel()
    else:
        detected = np.count_nonzero(np.asarray(matrix), axis=1)
    mito_mask = np.array([str(gene).upper().startswith("MT-") for gene in data.var_names])
    if mito_mask.any():
        mito_counts = np.asarray(matrix[:, mito_mask].sum(axis=1)).ravel()
    else:
        mito_counts = np.zeros(data.n_obs, dtype=float)
    mito_fraction = np.divide(mito_counts, total, out=np.zeros_like(total, dtype=float), where=total > 0)
    log_total = np.log1p(total)

    mito_available = bool(mito_mask.any() and np.nanstd(mito_fraction) > 0)
    result: dict[int, dict[str, float]] = {}
    for program in PROGRAMS:
        values = usage[program].to_numpy()
        result[program] = {
            "rho_log1p_total_counts": float(spearmanr(values, log_total).statistic),
            "rho_detected_genes": float(spearmanr(values, detected).statistic),
            "rho_mito_fraction": (
                float(spearmanr(values, mito_fraction).statistic) if mito_available else np.nan
            ),
        }
    metadata = {
        "mitochondrial_genes_in_input": int(mito_mask.sum()),
        "mito_fraction_metric_available": mito_available,
        "mito_fraction_unavailable_reason": (
            ""
            if mito_available
            else "No varying mitochondrial fraction can be calculated from the frozen cNMF input"
        ),
    }
    return result, metadata


def build_report(table: pd.DataFrame, output_path: Path) -> None:
    eligible = table.loc[table["che_eligibility"] == "eligible_for_confirmatory_che", "program_id"].tolist()
    control = table.loc[table["eligibility_use"] == "eligible_control_or_review", "program_id"].tolist()
    lines = [
        "# CRC K=10 逐程序联合稳定性与技术筛查人工审查报告",
        "",
        "## 一句话结论",
        "",
        f"有资格进入 Che 确认性外部复核的主程序为：{', '.join('P'+str(x) for x in eligible) if eligible else '无'}。",
        "这里的‘有资格’只表示 Liu 内部证据足以接受下一步外部检验，不表示已经被 Che 验证。",
        "",
        "## 本轮做了什么",
        "",
        "每个 Liu 主分析 K=10 程序都接受了跨 seed、细胞重抽样、患者 A/B 独立留出、两项恶性细胞集合敏感性、患者贡献和技术成分检查。",
        "所有阈值在生成资格结论前冻结；没有运行 Che，没有做程序生物学命名。",
        "",
        "## 技术指标可用性说明",
        "",
    ]
    if bool(table["mito_fraction_metric_available"].all()):
        lines.append("冻结输入可计算逐细胞线粒体比例及其与程序 usage 的相关性。")
    else:
        max_mito = int(table["top50_mitochondrial_n"].max())
        lines.extend(
            [
                "冻结 cNMF 输入不能得到有变化的逐细胞线粒体比例，因此该相关性记为不可用，而不是记为通过。",
                f"但所有程序 top-50 中线粒体基因最多为 {max_mito} 个；预注册相关性淘汰规则还要求至少 5 个线粒体 top gene，",
                "所以这个缺失指标不会改变本轮任何程序的资格。总 counts 和检测基因数相关性均已计算。",
            ]
        )
    lines.extend(
        [
            "",
        "## 逐程序结论",
        "",
        ]
    )
    for row in table.itertuples(index=False):
        label = "可进入 Che 确认性复核" if row.che_eligibility == "eligible_for_confirmatory_che" else "不可进入 Che 确认性复核"
        lines.append(f"### P{row.program_id}：{label}")
        lines.append("")
        lines.append(f"- top genes：{row.top20_genes}")
        lines.append(f"- 五项稳定性：{row.required_stability_pass_count}/5 通过。")
        lines.append(
            f"- 患者贡献：{'通过' if row.patient_contribution_pass else '失败'}；"
            f"三套完整输入中最小 effective patient={row.min_effective_patient_number:.2f}，"
            f"最少 ≥5% 贡献患者数={int(row.min_contributors_ge_5pct)}。"
        )
        if row.technical_dominant:
            lines.append(f"- 技术主导：是（{row.technical_failure_reasons}）。")
        else:
            lines.append("- 技术主导：否。")
        if row.review_flags:
            lines.append(f"- 仍需注意的非淘汰标记：{row.review_flags}。")
        lines.append(f"- 最终依据：{row.eligibility_reason_plain}")
        lines.append("")
    lines.extend(
        [
            "## 如何进入下一步",
            "",
            "Che 投射时建议盲法计算全部 10 个程序，以保留阴性对照并避免选择性报告；",
            "但只有本表预先标为 eligible_for_confirmatory_che 的程序，才能承担确认性外部复核结论。",
            "其余程序即使 Che 中偶然出现较高得分，也不能反向抹去 Liu 内部稳定性失败。",
            "",
            "## 仍需人工审查",
            "",
            "- 程序是否具有明确生物学含义尚未判断；本轮没有做命名或富集。",
            "- 细胞周期程序可以是真实肿瘤状态，因此只作控制/复核标记，不自动删除。",
            "- Che 只有 5 位配对患者，后续必须报告效应方向和不确定性，不能把阴性简单写成‘未验证’。",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--main-root", required=True, type=Path)
    parser.add_argument("--sensitivity-root", required=True, type=Path)
    parser.add_argument("--gene-sets", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    main_stability_path = args.main_root / "stability" / "P0_cnmf_stability.tsv"
    sensitivity_path = args.sensitivity_root / "P0_cnmf_sensitivity_same_k_pairs.tsv"
    top_gene_path = args.main_root / "runs" / ANCHOR_RUN / "P0_cnmf_top_genes.tsv"
    h5ad_path = args.main_root / "input_main" / "crc_liu_dual_cna_patient_state_balanced_counts.h5ad"
    usage_path = (
        args.main_root
        / "runs"
        / ANCHOR_RUN
        / "crc_liu_formal_cnmf"
        / "crc_liu_formal_cnmf.usages.k_10.dt_0_1.consensus.txt"
    )
    required_paths = [
        main_stability_path,
        sensitivity_path,
        top_gene_path,
        h5ad_path,
        usage_path,
        args.gene_sets,
        args.contract,
    ]
    for path in required_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty required input: {path}")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    main_stability = pd.read_csv(main_stability_path, sep="\t")
    sensitivity = pd.read_csv(sensitivity_path, sep="\t")
    cross_seed = one_to_one_block(
        main_stability,
        comparison_type="same_k_cross_seed",
        left_run=ANCHOR_RUN,
        right_run=OTHER_SEED_RUN,
    )
    resample = one_to_one_block(
        main_stability,
        comparison_type="same_k_cell_resample",
        left_run=ANCHOR_RUN,
        right_run=RESAMPLE_RUN,
    )
    holdout_a = one_to_one_block(
        main_stability,
        comparison_type="same_k_cross_input",
        left_run=ANCHOR_RUN,
        right_run=HOLDOUT_A_RUN,
    )
    holdout_b = one_to_one_block(
        main_stability,
        comparison_type="same_k_cross_input",
        left_run=ANCHOR_RUN,
        right_run=HOLDOUT_B_RUN,
    )
    bridge = one_to_one_block(
        main_stability,
        comparison_type="same_k_independent_patient_holdout",
        left_run=HOLDOUT_A_RUN,
        right_run=HOLDOUT_B_RUN,
    )
    exclude = one_to_one_block(
        sensitivity,
        comparison_type="same_k_cross_input",
        left_run="primary_seed_20260812",
        right_run="exclude_high_unresolved",
    )
    include = one_to_one_block(
        sensitivity,
        comparison_type="same_k_cross_input",
        left_run="primary_seed_20260812",
        right_run="include_one_method",
    )
    adjacent_k9 = one_to_one_block(
        main_stability,
        comparison_type="adjacent_k_same_seed",
        left_run=ANCHOR_RUN,
        right_run=ANCHOR_RUN,
        left_k=9,
        right_k=10,
        expected_rows=9,
    )
    adjacent_k11 = one_to_one_block(
        main_stability,
        comparison_type="adjacent_k_same_seed",
        left_run=ANCHOR_RUN,
        right_run=ANCHOR_RUN,
        left_k=10,
        right_k=11,
        expected_rows=10,
    )

    maps = {
        "cross_seed": map_by_left(cross_seed),
        "cell_resample": map_by_left(resample),
        "holdout_A": map_by_left(holdout_a),
        "holdout_B": map_by_left(holdout_b),
        "bridge": map_by_left(bridge),
        "exclude_high_unresolved": map_by_left(exclude),
        "include_one_method": map_by_left(include),
        "adjacent_k11": map_by_left(adjacent_k11),
    }
    k9_by_anchor = {int(row.right_program): row for row in adjacent_k9.itertuples(index=False)}

    evidence: list[dict[str, object]] = []
    for program in PROGRAMS:
        for component, source, required in [
            ("cross_seed", main_stability_path, True),
            ("cell_resample", main_stability_path, True),
            ("holdout_A", main_stability_path, True),
            ("holdout_B", main_stability_path, True),
            ("exclude_high_unresolved", sensitivity_path, True),
            ("include_one_method", sensitivity_path, True),
        ]:
            evidence.append(
                evidence_row(component, program, maps[component][program], required=required, source_file=source)
            )
        a_program = int(maps["holdout_A"][program].right_program)
        b_program = int(maps["holdout_B"][program].right_program)
        evidence.append(
            evidence_row(
                "holdout_A_to_B_bridge",
                program,
                maps["bridge"][a_program],
                required=True,
                source_file=main_stability_path,
                expected_right_program=b_program,
            )
        )
        evidence.append(
            evidence_row(
                "adjacent_K9_descriptive",
                program,
                k9_by_anchor.get(program),
                required=False,
                source_file=main_stability_path,
            )
        )
        evidence.append(
            evidence_row(
                "adjacent_K11_descriptive",
                program,
                maps["adjacent_k11"].get(program),
                required=False,
                source_file=main_stability_path,
            )
        )
    evidence_table = pd.DataFrame(evidence)

    top = pd.read_csv(top_gene_path, sep="\t")
    top = top.loc[(top["k"].astype(int) == K) & (top["rank"].astype(int) <= 50)].copy()
    top["program"] = top["program"].astype(int)
    top["rank"] = top["rank"].astype(int)
    if top.groupby("program").size().to_dict() != {program: 50 for program in PROGRAMS}:
        raise ValueError("Each K=10 anchor program must have exactly 50 ranked genes")
    gene_sets = load_gene_sets(args.gene_sets)
    annotations: list[dict[str, object]] = []
    category_counts: dict[int, dict[str, int]] = {}
    for program, group in top.groupby("program"):
        counts: dict[str, int] = {}
        for row in group.sort_values("rank").itertuples(index=False):
            categories = gene_categories(str(row.gene), gene_sets)
            for category in categories:
                counts[category] = counts.get(category, 0) + 1
            annotations.append(
                {
                    "program_id": int(program),
                    "rank": int(row.rank),
                    "gene": str(row.gene),
                    "technical_or_review_categories": ";".join(categories),
                }
            )
        category_counts[int(program)] = counts
    annotation_table = pd.DataFrame(annotations)
    qc_correlations, qc_metadata = compute_qc_correlations(h5ad_path, usage_path)

    mixing_paths = {
        ANCHOR_RUN: args.main_root / "runs" / ANCHOR_RUN / "P0_cnmf_patient_dataset_mixing.tsv",
        OTHER_SEED_RUN: args.main_root / "runs" / OTHER_SEED_RUN / "P0_cnmf_patient_dataset_mixing.tsv",
        RESAMPLE_RUN: args.main_root / "runs" / RESAMPLE_RUN / "P0_cnmf_patient_dataset_mixing.tsv",
    }
    for path in mixing_paths.values():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing patient contribution table: {path}")

    evidence_by_component = {
        component: group.set_index("anchor_program")
        for component, group in evidence_table.groupby("component")
    }
    records: list[dict[str, object]] = []
    for program in PROGRAMS:
        seed_row = maps["cross_seed"][program]
        resample_row = maps["cell_resample"][program]
        contribution = {
            ANCHOR_RUN: contribution_row(mixing_paths[ANCHOR_RUN], program),
            OTHER_SEED_RUN: contribution_row(mixing_paths[OTHER_SEED_RUN], int(seed_row.right_program)),
            RESAMPLE_RUN: contribution_row(mixing_paths[RESAMPLE_RUN], int(resample_row.right_program)),
        }
        min_effective = min(float(row["effective_unit_number"]) for row in contribution.values())
        min_contributors = min(int(row["contributors_ge_5pct"]) for row in contribution.values())
        patient_contribution_pass = min_effective >= 3 and min_contributors >= 3

        component_pass = {
            name: bool(evidence_by_component[name].loc[program, "component_pass"])
            for name in [
                "cross_seed",
                "cell_resample",
                "holdout_A",
                "holdout_B",
                "holdout_A_to_B_bridge",
                "exclude_high_unresolved",
                "include_one_method",
            ]
        }
        holdout_joint_pass = (
            component_pass["holdout_A"]
            and component_pass["holdout_B"]
            and component_pass["holdout_A_to_B_bridge"]
        )
        five_passes = {
            "cross_seed": component_pass["cross_seed"],
            "cell_resample": component_pass["cell_resample"],
            "patient_A_B_joint": holdout_joint_pass,
            "exclude_high_unresolved": component_pass["exclude_high_unresolved"],
            "include_one_method": component_pass["include_one_method"],
        }
        stability_pass = all(five_passes.values())
        counts = category_counts[program]
        correlations = qc_correlations[program]
        technical_dominant, technical_reasons = technical_decision(counts, correlations)
        review_flags: list[str] = []
        if counts.get("cell_cycle", 0) >= 10:
            review_flags.append("cell_cycle_dominant")
        if max(
            counts.get("immune_lineage", 0),
            counts.get("fibroblast_lineage", 0),
            counts.get("endothelial_lineage", 0),
        ) >= 5:
            review_flags.append("lineage_contamination_review")
        if counts.get("liver_plasma_abundance", 0) >= 5:
            review_flags.append("ambient_or_site_review")
        eligible = stability_pass and patient_contribution_pass and not technical_dominant
        eligibility = (
            "eligible_for_confirmatory_che" if eligible else "not_eligible_for_confirmatory_che"
        )
        if eligible and review_flags:
            eligibility_use = "eligible_control_or_review"
        elif eligible:
            eligibility_use = "core_candidate"
        else:
            eligibility_use = "negative_or_exploratory_only"
        failures = [name for name, passed in five_passes.items() if not passed]
        if not patient_contribution_pass:
            failures.append("patient_contribution")
        if technical_dominant:
            failures.append("technical_dominant")
        reason_plain = (
            "五项稳定性、患者贡献及技术筛查均通过。"
            if not failures
            else "未通过：" + "；".join(failures) + "。"
        )
        top20 = top.loc[top["program"] == program].sort_values("rank").head(20)["gene"].astype(str)
        records.append(
            {
                "program_id": program,
                "top20_genes": ",".join(top20),
                "cross_seed_target_program": int(seed_row.right_program),
                "cross_seed_pass": five_passes["cross_seed"],
                "cell_resample_target_program": int(resample_row.right_program),
                "cell_resample_pass": five_passes["cell_resample"],
                "holdout_A_target_program": int(maps["holdout_A"][program].right_program),
                "holdout_A_pass": component_pass["holdout_A"],
                "holdout_B_target_program": int(maps["holdout_B"][program].right_program),
                "holdout_B_pass": component_pass["holdout_B"],
                "holdout_bridge_actual_B_program": int(maps["bridge"][int(maps["holdout_A"][program].right_program)].right_program),
                "holdout_bridge_connects_expected_B": bool(
                    evidence_by_component["holdout_A_to_B_bridge"].loc[program, "connects_expected_program"]
                ),
                "holdout_bridge_pass": component_pass["holdout_A_to_B_bridge"],
                "patient_A_B_joint_pass": five_passes["patient_A_B_joint"],
                "exclude_high_unresolved_target_program": int(maps["exclude_high_unresolved"][program].right_program),
                "exclude_high_unresolved_pass": five_passes["exclude_high_unresolved"],
                "include_one_method_target_program": int(maps["include_one_method"][program].right_program),
                "include_one_method_pass": five_passes["include_one_method"],
                "required_stability_pass_count": int(sum(five_passes.values())),
                "all_required_stability_pass": stability_pass,
                "primary_effective_patient_number": float(contribution[ANCHOR_RUN]["effective_unit_number"]),
                "primary_contributors_ge_5pct": int(contribution[ANCHOR_RUN]["contributors_ge_5pct"]),
                "other_seed_mapped_program": int(seed_row.right_program),
                "other_seed_effective_patient_number": float(contribution[OTHER_SEED_RUN]["effective_unit_number"]),
                "other_seed_contributors_ge_5pct": int(contribution[OTHER_SEED_RUN]["contributors_ge_5pct"]),
                "resample_mapped_program": int(resample_row.right_program),
                "resample_effective_patient_number": float(contribution[RESAMPLE_RUN]["effective_unit_number"]),
                "resample_contributors_ge_5pct": int(contribution[RESAMPLE_RUN]["contributors_ge_5pct"]),
                "min_effective_patient_number": min_effective,
                "min_contributors_ge_5pct": min_contributors,
                "patient_contribution_pass": patient_contribution_pass,
                "top50_mitochondrial_n": counts.get("mitochondrial", 0),
                "top50_ribosomal_n": counts.get("ribosomal", 0),
                "top50_stress_dissociation_n": counts.get("stress_dissociation", 0),
                "top50_hemoglobin_n": counts.get("hemoglobin", 0),
                "top50_cell_cycle_n": counts.get("cell_cycle", 0),
                "top50_immune_lineage_n": counts.get("immune_lineage", 0),
                "top50_fibroblast_lineage_n": counts.get("fibroblast_lineage", 0),
                "top50_endothelial_lineage_n": counts.get("endothelial_lineage", 0),
                "top50_liver_plasma_abundance_n": counts.get("liver_plasma_abundance", 0),
                **correlations,
                "mito_fraction_metric_available": bool(qc_metadata["mito_fraction_metric_available"]),
                "mito_fraction_metric_note": str(qc_metadata["mito_fraction_unavailable_reason"]),
                "technical_dominant": technical_dominant,
                "technical_failure_reasons": ";".join(technical_reasons),
                "review_flags": ";".join(review_flags),
                "adjacent_K9_mapped_program": (
                    int(k9_by_anchor[program].left_program) if program in k9_by_anchor else np.nan
                ),
                "adjacent_K9_pass_descriptive": (
                    bool(k9_by_anchor[program].exceeds_both_nulls) if program in k9_by_anchor else False
                ),
                "adjacent_K11_mapped_program": int(maps["adjacent_k11"][program].right_program),
                "adjacent_K11_pass_descriptive": bool(maps["adjacent_k11"][program].exceeds_both_nulls),
                "che_eligibility": eligibility,
                "eligibility_use": eligibility_use,
                "failed_hard_gates": ";".join(failures),
                "eligibility_reason_plain": reason_plain,
            }
        )

    joint = pd.DataFrame(records).sort_values("program_id")
    if len(joint) != K or set(joint["program_id"]) != set(PROGRAMS):
        raise AssertionError("Joint screen must contain exactly programs 1-10")
    joint_path = args.output_dir / "P0_k10_program_joint_stability_technical_screen.tsv"
    evidence_path = args.output_dir / "P0_k10_program_match_evidence.tsv"
    annotation_path = args.output_dir / "P0_k10_program_top50_technical_annotations.tsv"
    report_path = args.output_dir / "P0_k10_program_che_eligibility_human_review.md"
    joint.to_csv(joint_path, sep="\t", index=False)
    evidence_table.to_csv(evidence_path, sep="\t", index=False)
    annotation_table.to_csv(annotation_path, sep="\t", index=False)
    build_report(joint, report_path)

    eligible_programs = joint.loc[
        joint["che_eligibility"] == "eligible_for_confirmatory_che", "program_id"
    ].astype(int).tolist()
    receipt = {
        "status": "READY_FOR_CHE_ELIGIBILITY_REVIEW",
        "anchor_run": ANCHOR_RUN,
        "k": K,
        "density_threshold": 0.1,
        "programs_screened": PROGRAMS,
        "eligible_for_confirmatory_che": eligible_programs,
        "not_eligible_for_confirmatory_che": [x for x in PROGRAMS if x not in eligible_programs],
        "che_executed": False,
        "biological_naming_executed": False,
        "technical_metric_availability": qc_metadata,
        "input_sha256": {str(path): sha256_file(path) for path in required_paths + list(mixing_paths.values())},
        "outputs": [str(joint_path), str(evidence_path), str(annotation_path), str(report_path)],
        "row_counts": {
            "joint_screen": len(joint),
            "match_evidence": len(evidence_table),
            "top50_annotations": len(annotation_table),
        },
    }
    receipt_path = args.output_dir / "P0_k10_program_che_eligibility_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
