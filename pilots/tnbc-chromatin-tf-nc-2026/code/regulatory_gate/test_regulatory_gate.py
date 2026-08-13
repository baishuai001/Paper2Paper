from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from phenotype import PhenotypeManifest, patient_table, scope_mask
from statistics import hedges_g, reml_meta, stratified_auc


def manifest(tmp_path: Path) -> PhenotypeManifest:
    payload = {
        "phenotype_id": "test",
        "label_column": "immune_infiltration_type",
        "case_values": ["M"],
        "control_values": ["B", "T", "desert"],
        "analysis_cell_column": "cell_type_coarse_crc_atlas",
        "analysis_cell_values": ["Cancer cell"],
        "patient_column": "donor_id",
        "dataset_column": "dataset",
        "scope": {
            "sample_type_equals": "Tumor",
            "exclude_tumor_source_contains": ["normal", "metasta"],
            "exclude_medical_condition_contains": ["normal", "polyp"],
            "treatment_status_before_resection_equals": "naive",
            "enrichment_cell_types_equals": "naive",
        },
        "primary_min_cells": 50,
    }
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return PhenotypeManifest.load(path)


class RegulatoryGateTests(unittest.TestCase):
    def test_manifest_scope_and_patient_counts(self) -> None:
        obs = pd.DataFrame(
            {
                "donor_id": ["p1", "p1", "p2", "p3"],
                "sample_id": ["s1", "s1", "s2", "s3"],
                "dataset": ["d1", "d1", "d2", "d3"],
                "study_id": ["a", "a", "b", "c"],
                "sample_type": ["Tumor"] * 4,
                "tumor_source": ["Primary", "Primary", "Primary", "Metastasis"],
                "medical_condition": ["CRC"] * 4,
                "treatment_status_before_resection": ["naive"] * 4,
                "enrichment_cell_types": ["naive"] * 4,
                "immune_infiltration_type": ["M", "M", "B", "T"],
                "cell_type_coarse_crc_atlas": ["Cancer cell", "T cell", "Cancer cell", "Cancer cell"],
            }
        )
        with tempfile.TemporaryDirectory() as directory:
            spec = manifest(Path(directory))
            self.assertEqual(scope_mask(obs, spec).tolist(), [True, True, True, False])
            patients, _ = patient_table(obs, spec)
        self.assertEqual(patients.set_index("donor_id").loc["p1", "analysis_cells"], 1)
        self.assertEqual(set(patients["group"]), {"case", "control"})

    def test_effect_meta_and_stratified_auc(self) -> None:
        effect, variance = hedges_g(np.array([2.0, 2.2, 2.4]), np.array([0.0, 0.2, 0.4]))
        self.assertGreater(effect, 2)
        self.assertGreater(variance, 0)
        meta = reml_meta(np.array([0.7, 0.8, 0.6, 0.9]), np.array([0.1, 0.1, 0.1, 0.1]))
        self.assertTrue(0.6 < meta["effect"] < 0.9)
        self.assertEqual(meta["I2"], 0)
        labels = np.array([False, False, True, True, False, False, True, True])
        predictions = np.array([0.1, 0.2, 0.8, 0.9, 0.2, 0.3, 0.7, 0.8])
        datasets = np.array(["d1"] * 4 + ["d2"] * 4)
        self.assertEqual(stratified_auc(labels, predictions, datasets), 1.0)


if __name__ == "__main__":
    unittest.main()
