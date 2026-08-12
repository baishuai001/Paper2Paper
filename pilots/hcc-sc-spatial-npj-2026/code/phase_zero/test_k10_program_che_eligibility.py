from __future__ import annotations

import unittest
from pathlib import Path
from types import SimpleNamespace

from build_k10_program_che_eligibility import evidence_row, technical_decision


class K10ProgramEligibilityUnitTests(unittest.TestCase):
    def test_stress_program_is_technical_dominant_at_frozen_threshold(self) -> None:
        counts = {"stress_dissociation": 10}
        correlations = {
            "rho_mito_fraction": 0.0,
            "rho_log1p_total_counts": 0.0,
            "rho_detected_genes": 0.0,
        }
        dominant, reasons = technical_decision(counts, correlations)
        self.assertTrue(dominant)
        self.assertIn("stress_dissociation_top50_ge_10", reasons)

    def test_cell_cycle_alone_is_not_automatic_technical_failure(self) -> None:
        counts = {"cell_cycle": 50}
        correlations = {
            "rho_mito_fraction": 0.0,
            "rho_log1p_total_counts": 0.0,
            "rho_detected_genes": 0.0,
        }
        dominant, reasons = technical_decision(counts, correlations)
        self.assertFalse(dominant)
        self.assertEqual(reasons, [])

    def test_library_size_correlation_requires_technical_top_genes(self) -> None:
        counts = {"stress_dissociation": 4}
        correlations = {
            "rho_mito_fraction": 0.0,
            "rho_log1p_total_counts": 0.95,
            "rho_detected_genes": 0.95,
        }
        dominant, reasons = technical_decision(counts, correlations)
        self.assertFalse(dominant)
        self.assertEqual(reasons, [])

    def test_patient_bridge_must_connect_the_two_anchor_mapped_programs(self) -> None:
        row = SimpleNamespace(
            left_run="patient_holdout_A",
            left_k=10,
            left_program=1,
            right_run="patient_holdout_B",
            right_k=10,
            right_program=6,
            cosine=0.9,
            top_gene_jaccard=0.5,
            cosine_null_q99=0.1,
            jaccard_null_q99=0.01,
            exceeds_both_nulls=True,
        )
        result = evidence_row(
            "holdout_A_to_B_bridge",
            2,
            row,
            required=True,
            source_file=Path("stability.tsv"),
            expected_right_program=3,
        )
        self.assertTrue(result["exceeds_both_nulls"])
        self.assertFalse(result["connects_expected_program"])
        self.assertFalse(result["component_pass"])


if __name__ == "__main__":
    unittest.main()
