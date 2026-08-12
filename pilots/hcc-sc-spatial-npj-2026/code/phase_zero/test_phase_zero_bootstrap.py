from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


table_module = load_module("build_crc_table_s1_design", "build_crc_table_s1_design.py")
h5ad_module = load_module("inspect_crc_h5ad", "inspect_crc_h5ad.py")
join_module = load_module("join_crc_h5ad_table_s1", "join_crc_h5ad_table_s1.py")
qc_module = load_module("audit_crc_atlas_qc", "audit_crc_atlas_qc.py")
integration_module = load_module("audit_crc_integration", "audit_crc_integration.py")
marker_module = load_module("audit_crc_markers", "audit_crc_markers.py")
cna_input_module = load_module("prepare_crc_cna_inputs", "prepare_crc_cna_inputs.py")
cna_panel_module = load_module("select_crc_cna_panel", "select_crc_cna_panel.py")
cna_panel_run_module = load_module("run_crc_cna_panel", "run_crc_cna_panel.py")
cna_panel_summary_module = load_module("summarize_crc_cna_panel", "summarize_crc_cna_panel.py")
cna_panel_plot_module = load_module("plot_crc_cna_panel", "plot_crc_cna_panel.py")
cnmf_input_module = load_module("prepare_crc_cnmf_input", "prepare_crc_cnmf_input.py")
cnmf_cna_input_module = load_module("prepare_crc_cnmf_from_cna", "prepare_crc_cnmf_from_cna.py")
cnmf_holdout_module = load_module("split_crc_cnmf_holdouts", "split_crc_cnmf_holdouts.py")
cnmf_run_module = load_module("run_crc_cnmf", "run_crc_cnmf.py")
cnmf_stability_module = load_module("evaluate_crc_cnmf_stability", "evaluate_crc_cnmf_stability.py")
copykat_summary_module = load_module("summarize_crc_copykat", "summarize_crc_copykat.py")
cna_concordance_module = load_module(
    "summarize_crc_cna_concordance", "summarize_crc_cna_concordance.py"
)


class PhaseZeroBootstrapTests(unittest.TestCase):
    def write_workbook(self, path: Path, duplicate_sample: bool = False) -> None:
        patients = pd.DataFrame(
            [
                {
                    "patient_id": "P1",
                    "dataset": "D1",
                    "study_id": "S1",
                    "treatment_status_before_resection": "naive",
                    "treatment_response": None,
                    "RECIST": None,
                    "platform": "10x",
                },
                {
                    "patient_id": "P2",
                    "dataset": "D2",
                    "study_id": "S2",
                    "treatment_status_before_resection": "treated",
                    "treatment_response": "PR",
                    "RECIST": "PR",
                    "platform": "10x",
                },
            ]
        )
        samples = pd.DataFrame(
            [
                {
                    "sample_id": "A",
                    "patient_id": "P1",
                    "dataset": "D1",
                    "study_id": "S1",
                    "sample_type": "normal",
                    "sample_tissue": "colon",
                },
                {
                    "sample_id": "A" if duplicate_sample else "B",
                    "patient_id": "P2",
                    "dataset": "D2",
                    "study_id": "S2",
                    "sample_type": "tumor",
                    "sample_tissue": "colon",
                },
            ]
        )
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            patients.to_excel(writer, sheet_name=table_module.PATIENT_SHEET, index=False, startrow=1)
            samples.to_excel(writer, sheet_name=table_module.SAMPLE_SHEET, index=False, startrow=1)

    def write_h5ad(self, path: Path, second_patient: str = "P2") -> None:
        obs = pd.DataFrame(
            {
                "sample_id": ["A", "B"],
                "patient_id": ["P1", second_patient],
                "donor_id": ["P1", second_patient],
                "dataset": ["D1", "D2"],
                "study_id": ["S1", "S2"],
                "sample_type": ["normal", "tumor"],
                "sample_tissue": ["colon", "colon"],
                "tissue": ["colon", "colon"],
                "cell_type_coarse": ["Epithelial", "Immune"],
                "cell_type_fine": ["Cancer Colonocyte-like", "CD8"],
                "cell_type_coarse_crc_atlas": ["Cancer cell", "T cell"],
                "cell_type_fine_crc_atlas": ["Cancer Colonocyte-like", "CD8"],
                "n_counts": [800.0, 900.0],
                "n_genes_by_counts": [300.0, 350.0],
                "pct_counts_mito": [5.0, 6.0],
                "SOLO_doublet_status": ["singlet", "singlet"],
                "SOLO_doublet_prob": [0.1, 0.2],
                "SOLO_singlet_prob": [0.9, 0.8],
            },
            index=["c1", "c2"],
        )
        var = pd.DataFrame(
            {"symbol": ["EPCAM", "PTPRC"], "GeneSymbol": ["EPCAM", "PTPRC"]},
            index=["g1", "g2"],
        )
        matrix = np.array([[2.0, 0.0], [0.0, 3.0]], dtype=np.float32)
        adata = ad.AnnData(X=matrix, obs=obs, var=var)
        adata.layers["counts"] = matrix.copy()
        adata.raw = adata.copy()
        adata.obsm["X_scANVI"] = np.array([[0.0, 0.0], [1.0, 1.0]], dtype=np.float32)
        adata.write_h5ad(path)

    def test_table_s1_builds_design_and_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workbook = root / "table.xlsx"
            self.write_workbook(workbook)
            receipt = table_module.build_design(workbook, root / "out")
            self.assertEqual(receipt["patient_rows"], 2)
            self.assertEqual(receipt["sample_rows"], 2)
            self.assertTrue((root / "out" / "P0_patient_sample_design_table_s1.tsv").is_file())

    def test_table_s1_rejects_duplicate_sample_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workbook = root / "table.xlsx"
            self.write_workbook(workbook, duplicate_sample=True)
            with self.assertRaisesRegex(ValueError, "duplicate identifiers"):
                table_module.build_design(workbook, root / "out")

    def test_h5ad_inspection_records_nonnegative_matrix(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            receipt = h5ad_module.inspect_h5ad(h5ad, root / "out", h5ad.stat().st_size)
            self.assertEqual(receipt["n_obs"], 2)
            self.assertEqual(receipt["n_vars"], 2)
            table = pd.read_csv(root / "out" / "P0_h5ad_matrices.tsv", sep="\t")
            self.assertEqual(bool(table.loc[0, "nonnegative"]), True)

    def test_h5ad_rejects_wrong_expected_size(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            with self.assertRaisesRegex(ValueError, "byte mismatch"):
                h5ad_module.inspect_h5ad(h5ad, root / "out", h5ad.stat().st_size + 1)

    def test_h5ad_table_join_accepts_table_superset(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workbook = root / "table.xlsx"
            self.write_workbook(workbook)
            table_module.build_design(workbook, root / "table")
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            receipt = join_module.build_join(
                h5ad,
                root / "table" / "P0_table_s1_samples.tsv",
                root / "join",
                expected_obs=2,
            )
            self.assertEqual(receipt["status"], "passed")
            self.assertEqual(receipt["samples_in_both"], 2)

    def test_h5ad_table_join_rejects_patient_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workbook = root / "table.xlsx"
            self.write_workbook(workbook)
            table_module.build_design(workbook, root / "table")
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad, second_patient="WRONG")
            with self.assertRaisesRegex(ValueError, "patient_id"):
                join_module.build_join(
                    h5ad,
                    root / "table" / "P0_table_s1_samples.tsv",
                    root / "join",
                    expected_obs=2,
                )

    def test_cellxgene_identity_mapping_is_explicit_and_tested(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workbook = root / "table.xlsx"
            self.write_workbook(workbook)
            table_module.build_design(workbook, root / "table")
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            field_map = {
                "sample_id": "sample_id",
                "patient_id": "donor_id",
                "dataset": "dataset",
                "study_id": "study_id",
                "sample_type": "sample_type",
                "sample_tissue": "tissue",
            }
            join_receipt = join_module.build_join(
                h5ad,
                root / "table" / "P0_table_s1_samples.tsv",
                root / "join",
                expected_obs=2,
                field_map=field_map,
            )
            self.assertEqual(join_receipt["h5ad_field_map"]["patient_id"], "donor_id")
            qc_receipt = qc_module.audit_qc(
                h5ad,
                root / "qc",
                expected_obs=2,
                identity_field_map={
                    "sample_id": "sample_id",
                    "patient_id": "donor_id",
                    "dataset": "dataset",
                    "sample_type": "sample_type",
                    "sample_tissue": "tissue",
                },
            )
            self.assertEqual(qc_receipt["identity_field_map"]["sample_tissue"], "tissue")

    def test_pooled_sample_keeps_donor_units_and_reports_limitation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            workbook = root / "table.xlsx"
            self.write_workbook(workbook)
            table_module.build_design(workbook, root / "table")
            base = root / "base.h5ad"
            self.write_h5ad(base)
            adata = ad.read_h5ad(base)
            extra = adata[[0]].copy()
            extra.obs_names = ["c3"]
            extra.obs.loc["c3", "patient_id"] = "P3"
            extra.obs.loc["c3", "donor_id"] = "P3"
            pooled = ad.concat([adata, extra], merge="same")
            pooled_path = root / "pooled.h5ad"
            pooled.write_h5ad(pooled_path)
            receipt = join_module.build_join(
                pooled_path,
                root / "table" / "P0_table_s1_samples.tsv",
                root / "join",
                expected_obs=3,
            )
            self.assertEqual(receipt["status"], "passed_with_limitations")
            self.assertEqual(receipt["pooled_sample_ids"], ["A"])
            units = pd.read_csv(root / "join" / "P0_h5ad_analysis_units.tsv", sep="\t")
            self.assertEqual(len(units), 3)
            qc_receipt = qc_module.audit_qc(pooled_path, root / "qc", expected_obs=3)
            self.assertEqual(qc_receipt["sample_patient_analysis_units"], 3)

    def test_qc_audit_records_post_qc_and_doublet_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            receipt = qc_module.audit_qc(h5ad, root / "qc", expected_obs=2)
            self.assertEqual(receipt["status"], "passed")
            self.assertEqual(receipt["published_threshold_violations"]["count_below_or_equal_400"], 0)
            self.assertEqual(receipt["samples_without_probability_evidence"], 0)

    def test_integration_audit_is_explicitly_diagnostic_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            receipt = integration_module.audit_integration(
                h5ad,
                root / "integration",
                max_cells=100,
                neighbors=1,
                expected_obs=2,
            )
            self.assertEqual(receipt["status"], "diagnostic_only_no_preintegration_comparator")
            self.assertTrue((root / "integration" / "P0_integration_metrics.tsv").is_file())

    def test_marker_audit_uses_author_markers_without_overwriting_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            marker_table = root / "markers.csv"
            marker_table.write_text("cell_type,symbol\nCancer cell,EPCAM\nCD8,PTPRC\n", encoding="utf-8")
            receipt = marker_module.audit_markers(
                h5ad,
                marker_table,
                root / "markers",
                per_patient_label_cap=2,
                expected_obs=2,
            )
            self.assertEqual(receipt["status"], "diagnostic_only_author_labels_not_overwritten")
            evidence = pd.read_csv(root / "markers" / "P0_annotation_evidence.tsv", sep="\t")
            self.assertEqual(set(evidence["gene"]), {"EPCAM", "PTPRC"})
            gene_map = pd.read_csv(root / "markers" / "P0_annotation_gene_map.tsv", sep="\t")
            self.assertEqual(set(gene_map["selection_rule"]), {"lowest_feature_index_no_var_n_cells_available"})

    def test_duplicate_marker_features_are_resolved_outcome_blind(self):
        adata = ad.AnnData(
            X=np.ones((2, 3), dtype=np.float32),
            var=pd.DataFrame(
                {"GeneSymbol": ["EPCAM", "EPCAM", "PTPRC"], "n_cells": [5, 10, 20]},
                index=["g1", "g2", "g3"],
            ),
        )
        resolved, mapping = marker_module.resolve_marker_features(
            adata,
            marker_module.gene_symbols(adata, "GeneSymbol"),
            ["EPCAM", "PTPRC"],
        )
        self.assertEqual(resolved["EPCAM"], 1)
        self.assertEqual(int(mapping.loc[mapping["gene"] == "EPCAM", "candidate_features"].iloc[0]), 2)

    def test_cna_input_uses_integer_raw_counts_and_known_normals(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            adata = ad.read_h5ad(h5ad)
            adata.obs.loc["c2", ["sample_id", "patient_id", "donor_id"]] = ["A", "P1", "P1"]
            adata.write_h5ad(h5ad)
            receipt = cna_input_module.prepare_cna_input(
                h5ad,
                root / "cna",
                sample_id="A",
                patient_id="P1",
                max_cancer=2,
                max_epithelial=2,
                max_reference=2,
                min_cancer=1,
                min_reference=1,
            )
            self.assertEqual(receipt["raw_matrix"], "raw/X")
            self.assertEqual(receipt["role_counts"], {"author_cancer": 1, "known_normal_reference": 1})
            self.assertTrue((root / "cna" / "counts_genes_by_cells.mtx.gz").is_file())

    def test_cnmf_input_is_patient_balanced_raw_counts_and_diagnostic_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            receipt = cnmf_input_module.prepare_cnmf_input(
                h5ad,
                root / "cnmf",
                max_datasets=1,
                max_patients_per_dataset=1,
                max_cells_per_patient=1,
                min_cells_per_patient=1,
            )
            self.assertEqual(receipt["status"], "diagnostic_author_cancer_input_pending_cna")
            self.assertEqual(receipt["raw_matrix"], "raw/X")
            self.assertEqual(receipt["cells"], 1)
            self.assertEqual(receipt["patients"], 1)
            prepared = ad.read_h5ad(root / "cnmf" / "crc_author_cancer_balanced_counts.h5ad")
            self.assertTrue(np.allclose(prepared.X.data, np.rint(prepared.X.data)))
            self.assertIn("source_cell_id", prepared.obs.columns)

    def test_cnmf_runner_rejects_negative_or_noninteger_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "tiny.h5ad"
            self.write_h5ad(h5ad)
            data = ad.read_h5ad(h5ad)
            data.obs["analysis_patient_id"] = data.obs["donor_id"].astype(str)
            data.write_h5ad(h5ad)
            validated = cnmf_run_module.validate_counts_input(h5ad)
            self.assertEqual(validated["cells"], 2)
            data = ad.read_h5ad(h5ad)
            data.X[0, 0] = -1.0
            data.write_h5ad(h5ad)
            with self.assertRaisesRegex(ValueError, "negative"):
                cnmf_run_module.validate_counts_input(h5ad)

    def test_cnmf_runner_rejects_too_few_consensus_initializations(self):
        with self.assertRaisesRegex(ValueError, "at least 5"):
            cnmf_run_module.run_cnmf(
                Path("unused.h5ad"), Path("unused-output"), "unused", [5], 2, 1, 2.0, 100, 100
            )

    def test_cnmf_program_matching_is_one_to_one_and_gene_aware(self):
        left = pd.DataFrame(
            [[10.0, 5.0, 0.0, 0.0], [0.0, 0.0, 5.0, 10.0]],
            index=["L1", "L2"],
            columns=["A", "B", "C", "D"],
        )
        right = pd.DataFrame(
            [[0.0, 0.0, 6.0, 9.0], [9.0, 6.0, 0.0, 0.0]],
            index=["R1", "R2"],
            columns=["A", "B", "C", "D"],
        )
        matches = cnmf_stability_module.match_programs(left, right, top_n=2)
        observed = set(zip(matches["left_program"], matches["right_program"]))
        self.assertEqual(observed, {("L1", "R2"), ("L2", "R1")})
        self.assertTrue((matches["top_gene_jaccard"] == 1.0).all())
        cosine_q99, jaccard_q99 = cnmf_stability_module.empirical_thresholds(
            left, right, top_n=2, permutations=5, rng=np.random.default_rng(7)
        )
        self.assertTrue(np.isfinite(cosine_q99))
        self.assertGreaterEqual(jaccard_q99, 0.0)

    def test_cnmf_spectra_loader_preserves_program_by_gene_orientation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run_dir = root / "tiny"
            run_dir.mkdir()
            spectra = pd.DataFrame(
                [[1.0, 0.5, 0.0], [0.0, 0.5, 1.0]],
                index=[1, 2],
                columns=["G1", "G2", "G3"],
            )
            spectra.to_csv(
                run_dir / "tiny.gene_spectra_score.k_2.dt_0_1.txt", sep="\t", index=True
            )
            loaded = cnmf_stability_module.discover_spectra(root, "tiny")
            self.assertEqual(loaded[2].shape, (2, 3))
            self.assertEqual(list(loaded[2].columns), ["G1", "G2", "G3"])

    def test_cnmf_spectra_loader_rejects_transposed_or_wrong_k_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run_dir = root / "tiny"
            run_dir.mkdir()
            pd.DataFrame(
                [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]],
                index=["G1", "G2", "G3"],
                columns=[1, 2],
            ).to_csv(run_dir / "tiny.gene_spectra_score.k_2.dt_0_1.txt", sep="\t")
            with self.assertRaisesRegex(ValueError, "orientation/shape mismatch"):
                cnmf_stability_module.discover_spectra(root, "tiny")

    def test_copykat_summary_keeps_diploid_calls_non_exclusionary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = pd.DataFrame(
                {
                    "cell_id": ["c1", "c2", "c3"],
                    "cna_input_role": ["author_cancer", "author_cancer", "known_normal_reference"],
                    "is_known_normal_reference": [False, False, True],
                }
            )
            prediction = pd.DataFrame(
                {"cell.names": ["c1", "c2", "c3"], "copykat.pred": ["aneuploid", "diploid", "diploid"]}
            )
            metadata.to_csv(root / "metadata.tsv", sep="\t", index=False)
            prediction.to_csv(root / "prediction.tsv", sep="\t", index=False)
            receipt = copykat_summary_module.summarize(
                root / "metadata.tsv", root / "prediction.tsv", root / "summary"
            )
            self.assertEqual(receipt["author_cancer_aneuploid"], 1)
            table = pd.read_csv(root / "summary" / "P0_copykat_role_prediction.tsv", sep="\t")
            diploid = table.loc[table["copykat.pred"] == "diploid", "cna_evidence_class"]
            self.assertEqual(set(diploid), {"diploid_does_not_exclude_malignancy"})

    def test_cna_panel_selection_is_dataset_balanced_and_outcome_blind(self):
        table = pd.DataFrame(
            {
                "dataset": ["D1", "D1", "D2", "D3", "D4"],
                "sample_id": ["s1", "s2", "s3", "s4", "s5"],
                "analysis_patient_id": ["p1", "p2", "p3", "p4", "p5"],
                "sample_type": ["tumor", "tumor", "tumor", "polyp", "polyp"],
                "tissue": ["colon"] * 5,
                "cancer_cells": [100, 200, 80, 90, 70],
                "reference_cells": [100, 60, 90, 90, 80],
                "limiting_role_cells": [100, 60, 80, 90, 70],
                "eligible": [True] * 5,
            }
        )
        selected = cna_panel_module.select_panel(table, ["tumor", "polyp"], 2)
        self.assertEqual(len(selected), 4)
        self.assertEqual(
            selected.loc[selected["sample_type"].eq("tumor"), "sample_id"].tolist(),
            ["s1", "s3"],
        )
        self.assertFalse(selected.duplicated(["sample_type", "dataset"]).any())
        self.assertTrue(selected["selection_rule"].str.contains("outcome-blind").all())

    def test_cna_panel_upstream_resource_receipt_fails_closed_on_size_change(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "data.h5ad"
            h5ad.write_bytes(b"abc")
            receipt = root / "receipt.json"
            receipt.write_text(
                json.dumps(
                    {
                        "status": "passed",
                        "h5ad": str(h5ad.resolve()),
                        "bytes": 3,
                        "sha256": "A" * 64,
                    }
                ),
                encoding="utf-8",
            )
            verified = cna_panel_module.verify_upstream_h5ad_receipt(h5ad, receipt)
            self.assertEqual(verified["sha256"], "A" * 64)
            h5ad.write_bytes(b"abcd")
            with self.assertRaisesRegex(ValueError, "byte size changed"):
                cna_panel_module.verify_upstream_h5ad_receipt(h5ad, receipt)

    def test_cna_panel_runner_rejects_duplicate_analysis_units(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "panel.tsv"
            pd.DataFrame(
                {
                    "panel_order": [1, 2],
                    "dataset": ["D1", "D1"],
                    "sample_id": ["S1", "S1"],
                    "analysis_patient_id": ["P1", "P1"],
                    "sample_type": ["tumor", "tumor"],
                    "tissue": ["colon", "colon"],
                    "eligible": [True, True],
                }
            ).to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "duplicate sample-patient"):
                cna_panel_run_module.load_panel(path)

    def test_cna_panel_summary_rejects_incomplete_run_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "P0_cna_panel_run_receipt.json").write_text(
                json.dumps(
                    {
                        "status": "completed_with_failures",
                        "units_failed": 1,
                        "units": [{"panel_order": 1, "status": "scevan_failed"}],
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "not a complete successful run"):
                cna_panel_summary_module.read_completed_state(root)

    def test_cna_panel_cell_classes_preserve_support_and_nonproof(self):
        table = pd.DataFrame(
            {
                "cna_input_role": [
                    "author_cancer",
                    "author_cancer",
                    "author_cancer",
                    "known_normal_reference",
                ],
                "copykat_malignancy_support": [True, True, False, False],
                "scevan_malignancy_support": [True, False, False, False],
                "either_malignancy_support": [True, True, False, False],
                "method_unresolved": [False, False, False, False],
            }
        )
        classes = cna_panel_summary_module.classify_cell_evidence(table)
        self.assertEqual(classes.iloc[0], "author_cancer_two_method_malignancy_support")
        self.assertEqual(classes.iloc[1], "author_cancer_one_method_malignancy_support")
        self.assertEqual(
            classes.iloc[2], "author_cancer_without_cna_support_not_nonmalignancy_proof"
        )
        self.assertEqual(
            classes.iloc[3], "normal_reference_without_malignancy_support_not_truth"
        )

    def test_cna_panel_plot_source_requires_exclusive_classes_to_sum(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "units.tsv"
            row = {
                "panel_order": 1,
                "sample_type": "tumor",
                "author_cancer_cells": 10,
                "author_cancer_two_method_support": 4,
                "author_cancer_one_method_support": 2,
                "author_cancer_unresolved_without_support": 1,
                "author_cancer_defined_without_support_not_nonmalignancy_proof": 3,
                "normal_reference_cells": 20,
                "normal_reference_copykat_aneuploid": 1,
                "normal_reference_scevan_tumor": 2,
                "normal_reference_unexpected_malignancy_support": 2,
            }
            pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
            source = cna_panel_plot_module.build_plot_source(path)
            self.assertAlmostEqual(
                source.loc[0, "author_cancer_two_method_support_fraction"], 0.4
            )
            row["author_cancer_defined_without_support_not_nonmalignancy_proof"] = 4
            pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "do not sum"):
                cna_panel_plot_module.build_plot_source(path)

    def test_cna_derived_cnmf_balancing_is_deterministic_and_state_guarded(self):
        rows = []
        for order, state in ((1, "polyp"), (2, "tumor")):
            for index in range(8):
                rows.append(
                    {
                        "cell_id": f"c{order}_{index}",
                        "panel_order": order,
                        "panel_dataset": f"D{order}",
                        "panel_sample_id": f"S{order}",
                        "panel_patient_id": f"P{order}",
                        "panel_sample_type": state,
                        "panel_tissue": "colon",
                    }
                )
        selected = pd.DataFrame(rows)
        balanced_a, units_a = cnmf_cna_input_module.balance_supported_cells(
            selected, min_cells_per_unit=5, max_cells_per_unit=5, seed=17, required_states=["polyp", "tumor"]
        )
        balanced_b, units_b = cnmf_cna_input_module.balance_supported_cells(
            selected, min_cells_per_unit=5, max_cells_per_unit=5, seed=17, required_states=["polyp", "tumor"]
        )
        pd.testing.assert_frame_equal(balanced_a, balanced_b)
        pd.testing.assert_frame_equal(units_a, units_b)
        self.assertEqual(units_a["selected_cells"].tolist(), [5, 5])
        with self.assertRaisesRegex(ValueError, "required CRC states"):
            cnmf_cna_input_module.balance_supported_cells(
                selected, min_cells_per_unit=5, max_cells_per_unit=5, seed=17, required_states=["metastasis"]
            )

    def test_cnmf_validator_records_two_method_cna_selection_basis(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.h5ad"
            data = ad.AnnData(
                X=np.asarray([[1, 2], [2, 1]], dtype=np.float32),
                obs=pd.DataFrame(
                    {
                        "analysis_patient_id": ["P1", "P2"],
                        "dataset": ["D1", "D2"],
                        "sample_type": ["tumor", "metastasis"],
                        "cna_selection_class": [
                            "author_cancer_two_method_malignancy_support",
                            "author_cancer_two_method_malignancy_support",
                        ],
                    },
                    index=["c1", "c2"],
                ),
                var=pd.DataFrame(index=["g1", "g2"]),
            )
            data.write_h5ad(path)
            receipt = cnmf_run_module.validate_counts_input(path)
            self.assertEqual(
                receipt["selection_basis"],
                "copykat_aneuploid_and_scevan_tumor_intersection",
            )

    def test_cna_derived_cnmf_input_identity_joins_and_sorts_backed_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "atlas.h5ad"
            data = ad.AnnData(
                X=np.asarray([[1, 2], [2, 3], [3, 4]], dtype=np.float32),
                obs=pd.DataFrame(
                    {
                        "sample_id": ["S1", "S1", "S1"],
                        "donor_id": ["P1", "P1", "P1"],
                        "dataset": ["D1", "D1", "D1"],
                        "sample_type": ["tumor", "tumor", "tumor"],
                        "tissue": ["colon", "colon", "colon"],
                    },
                    index=["c1", "c2", "c3"],
                ),
                var=pd.DataFrame(
                    {"GeneSymbol": ["G1", "G2"], "n_cells": [3, 3]},
                    index=["f1", "f2"],
                ),
            )
            data.raw = data.copy()
            data.write_h5ad(h5ad)
            evidence = root / "evidence.tsv"
            pd.DataFrame(
                {
                    "cell_id": ["c3", "c1", "c2", "r1"],
                    "panel_order": [1, 1, 1, 1],
                    "panel_dataset": ["D1", "D1", "D1", "D1"],
                    "panel_sample_id": ["S1", "S1", "S1", "S1"],
                    "panel_patient_id": ["P1", "P1", "P1", "P1"],
                    "panel_sample_type": ["tumor", "tumor", "tumor", "tumor"],
                    "panel_tissue": ["colon", "colon", "colon", "colon"],
                    "cna_input_role": ["author_cancer"] * 3 + ["known_normal_reference"],
                    "copykat.pred": ["aneuploid"] * 3 + ["diploid"],
                    "scevan_class": ["tumor"] * 3 + ["normal"],
                    "cna_selection_class": [
                        "author_cancer_two_method_malignancy_support"
                    ]
                    * 3
                    + ["normal_reference_without_malignancy_support_not_truth"],
                }
            ).to_csv(evidence, sep="\t", index=False)
            receipt = cnmf_cna_input_module.prepare_cnmf_from_cna(
                h5ad,
                evidence,
                root / "output",
                min_cells_per_unit=1,
                max_cells_per_unit=3,
                required_states=["tumor"],
            )
            self.assertEqual(receipt["cells"], 3)
            self.assertEqual(receipt["units_with_unexpected_normal_reference_support"], 0)
            prepared = ad.read_h5ad(
                root / "output" / "crc_two_method_cna_supported_counts.h5ad"
            )
            self.assertEqual(prepared.obs_names.tolist(), ["c1", "c2", "c3"])
            self.assertTrue(np.array_equal(prepared.X.toarray(), data.X))

    def test_cnmf_holdout_folds_are_state_stratified_patient_dataset_disjoint(self):
        rows = []
        index = []
        for state_index, state in enumerate(("polyp", "tumor", "metastasis"), start=1):
            for replicate in (1, 2):
                order = state_index * 10 + replicate
                rows.append(
                    {
                        "panel_order": order,
                        "dataset": f"D{order}",
                        "analysis_patient_id": f"P{order}",
                        "sample_id": f"S{order}",
                        "sample_type": state,
                        "cna_selection_class": "author_cancer_two_method_malignancy_support",
                    }
                )
                index.append(f"c{order}")
        obs = pd.DataFrame(rows, index=index)
        manifest = cnmf_holdout_module.assign_folds(obs, seed=19)
        self.assertEqual(set(manifest["holdout_fold"]), {"A", "B"})
        for _, group in manifest.groupby("sample_type"):
            self.assertEqual(set(group["holdout_fold"]), {"A", "B"})
        fold_a = manifest[manifest["holdout_fold"] == "A"]
        fold_b = manifest[manifest["holdout_fold"] == "B"]
        self.assertFalse(set(fold_a["dataset"]) & set(fold_b["dataset"]))
        self.assertFalse(
            set(fold_a["analysis_patient_id"]) & set(fold_b["analysis_patient_id"])
        )
        duplicated = obs.copy()
        duplicated.loc[duplicated.index[-1], "dataset"] = duplicated.iloc[0]["dataset"]
        with self.assertRaisesRegex(ValueError, "dataset- and patient-disjoint"):
            cnmf_holdout_module.assign_folds(duplicated, seed=19)

    def test_cnmf_stability_distinguishes_seed_and_independent_holdout_evidence(self):
        self.assertEqual(
            cnmf_stability_module.cross_run_comparison_type("seed_1", "seed_2"),
            "same_k_cross_seed",
        )
        self.assertEqual(
            cnmf_stability_module.cross_run_comparison_type("holdout_A", "holdout_B"),
            "same_k_independent_patient_dataset_holdout",
        )
        self.assertEqual(
            cnmf_stability_module.cross_run_comparison_type("seed_1", "holdout_A"),
            "same_k_cross_input",
        )
        self.assertEqual(
            cnmf_stability_module.evaluation_status(False),
            "diagnostic_seed_and_adjacent_k_comparison_only",
        )
        self.assertEqual(
            cnmf_stability_module.evaluation_status(True),
            "diagnostic_seed_adjacent_k_and_independent_holdout_comparison",
        )
        synthetic_matches = pd.DataFrame(
            {
                "comparison_type": [
                    "same_k_independent_patient_dataset_holdout",
                    "same_k_independent_patient_dataset_holdout",
                    "same_k_independent_patient_dataset_holdout",
                    "adjacent_k_same_seed",
                ],
                "left_k": [5, 5, 6, 5],
                "right_k": [5, 5, 6, 6],
                "left_program": ["1", "2", "1", "1"],
                "right_program": ["1", "2", "1", "1"],
                "cosine": [0.9, 0.8, 0.85, 0.95],
                "top_gene_jaccard": [0.5, 0.4, 0.45, 0.6],
                "cosine_null_q99": [0.3, 0.3, 0.3, 0.3],
                "jaccard_null_q99": [0.1, 0.1, 0.1, 0.1],
                "exceeds_both_nulls": [True, False, True, True],
            }
        )
        direct_pairs, by_k = (
            cnmf_stability_module.summarize_independent_holdout_pairs(
                synthetic_matches
            )
        )
        self.assertEqual(direct_pairs["holdout_pair_id"].tolist(), ["K5-M01", "K5-M02", "K6-M01"])
        self.assertEqual(by_k["k"].tolist(), [5, 6])
        self.assertEqual(by_k["matched_pairs"].tolist(), [2, 1])
        self.assertEqual(by_k["pairs_exceeding_both_nulls"].tolist(), [1, 1])

    def test_cna_concordance_preserves_support_conflict_and_unresolved_classes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = pd.DataFrame(
                {
                    "cell_id": ["c1", "c2", "c3", "c4", "c5", "c6"],
                    "cna_input_role": [
                        "author_cancer",
                        "author_cancer",
                        "author_cancer",
                        "author_cancer",
                        "known_normal_reference",
                        "known_normal_reference",
                    ],
                    "is_known_normal_reference": [False, False, False, False, True, True],
                }
            )
            copykat = pd.DataFrame(
                {
                    "cell.names": ["c1", "c2", "c3", "c4", "c5", "c6"],
                    "copykat.pred": [
                        "aneuploid",
                        "aneuploid",
                        "diploid",
                        "not.defined",
                        "diploid",
                        "not.defined",
                    ],
                }
            )
            scevan = pd.DataFrame(
                {
                    "cell.names": ["c1", "c2", "c3", "c4", "c5", "c6"],
                    "class": ["tumor", "normal", "tumor", "filtered", "normal", "normal"],
                }
            )
            metadata.to_csv(root / "metadata.tsv", sep="\t", index=False)
            copykat.to_csv(root / "copykat.tsv", sep="\t", index=False)
            scevan.to_csv(root / "scevan.tsv", sep="\t", index=False)
            receipt = cna_concordance_module.summarize(
                root / "metadata.tsv",
                root / "copykat.tsv",
                root / "scevan.tsv",
                root / "summary",
            )
            self.assertEqual(receipt["author_cancer_both_malignancy_support"], 1)
            self.assertEqual(receipt["author_cancer_either_malignancy_support"], 3)
            self.assertEqual(receipt["author_cancer_method_unresolved"], 1)
            detail = pd.read_csv(root / "summary" / "P0_cna_method_pair_by_role.tsv", sep="\t")
            self.assertIn("discordant_aneuploid_vs_scevan_normal", set(detail["evidence_class"]))
            self.assertIn(
                "scevan_only_support_copykat_diploid_nonexclusive", set(detail["evidence_class"])
            )
            self.assertIn("unresolved_both_methods", set(detail["evidence_class"]))

    def test_cna_concordance_fails_closed_on_incomplete_cell_identity_join(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pd.DataFrame(
                {
                    "cell_id": ["c1", "c2"],
                    "cna_input_role": ["author_cancer", "known_normal_reference"],
                    "is_known_normal_reference": [False, True],
                }
            ).to_csv(root / "metadata.tsv", sep="\t", index=False)
            pd.DataFrame(
                {"cell.names": ["c1", "c2"], "copykat.pred": ["aneuploid", "diploid"]}
            ).to_csv(root / "copykat.tsv", sep="\t", index=False)
            pd.DataFrame({"cell.names": ["c1"], "class": ["tumor"]}).to_csv(
                root / "scevan.tsv", sep="\t", index=False
            )
            with self.assertRaisesRegex(ValueError, "SCEVAN identity join is incomplete"):
                cna_concordance_module.summarize(
                    root / "metadata.tsv",
                    root / "copykat.tsv",
                    root / "scevan.tsv",
                    root / "summary",
                )

    def test_scevan_adapter_is_pinned_to_the_actual_singular_normal_argument(self):
        runner = (HERE / "run_crc_cna.R").read_text(encoding="utf-8")
        self.assertIn('packageVersion("SCEVAN")', runner)
        self.assertIn("norm_cell = known_normal", runner)
        self.assertNotIn("norm_cells = known_normal", runner)
        self.assertIn("output_dir = output_dir", runner)
        self.assertIn('scevan_counts <- as(input$counts, "CsparseMatrix")', runner)
        self.assertIn('inherits(scevan_counts, "dgCMatrix")', runner)
        self.assertIn("SCEVAN sparse-class conversion changed matrix dimensions or dimnames", runner)
        self.assertIn("scevan_counts,", runner)
        self.assertIn('legacy_output_alias <- file.path(output_dir, "output")', runner)
        self.assertIn("alias_is_link <- !is.na(alias_target) && nzchar(alias_target)", runner)
        self.assertIn("SCEVAN legacy ./output alias already exists", runner)
        self.assertIn("file.symlink(output_dir, legacy_output_alias)", runner)
        self.assertIn("on.exit(unlink(legacy_output_alias), add = TRUE)", runner)


if __name__ == "__main__":
    unittest.main()
