from __future__ import annotations

import unittest
from pathlib import Path

from evaluate_crc_cnmf_stability import cross_run_comparison_type, evaluation_status
import pandas as pd

from run_crc_cnmf import relabel_k_selection_density_metadata, run_cnmf


class FormalCnmfRuntimeTests(unittest.TestCase):
    def test_k_selection_default_is_not_mislabeled_as_final_consensus_threshold(self) -> None:
        stats = pd.DataFrame(
            {
                "k": [5.0],
                "local_density_threshold": [0.5],
                "silhouette": [0.8],
                "prediction_error": [10.0],
            }
        )
        observed = relabel_k_selection_density_metadata(stats, 0.1)
        self.assertNotIn("local_density_threshold", observed.columns)
        self.assertEqual(observed.loc[0, "cnmf_reported_default_threshold_no_filter"], 0.5)
        self.assertEqual(observed.loc[0, "final_consensus_density_threshold"], 0.1)
        self.assertEqual(
            observed.loc[0, "k_selection_filtering"],
            "disabled_by_cnmf_skip_density_stats_path",
        )

    def test_formal_stability_roles_are_not_mislabeled_as_dataset_holdout(self) -> None:
        self.assertEqual(
            cross_run_comparison_type("patient_holdout_A", "patient_holdout_B"),
            "same_k_independent_patient_holdout",
        )
        self.assertEqual(
            cross_run_comparison_type("seed_20260812", "cell_resample_20260813"),
            "same_k_cell_resample",
        )
        self.assertEqual(
            evaluation_status(True, True),
            "candidate_seed_cell_resample_and_patient_holdout_comparison",
        )

    def test_parallel_worker_count_fails_closed_before_execution(self) -> None:
        with self.assertRaisesRegex(ValueError, "total_workers"):
            run_cnmf(
                Path("missing.h5ad"),
                Path("unused"),
                "unused",
                [5],
                5,
                1,
                0.1,
                2000,
                1000,
                total_workers=0,
            )


if __name__ == "__main__":
    unittest.main()
