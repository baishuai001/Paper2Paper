from __future__ import annotations

import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from freeze_liu_k10_projection_reference import rank_scores
from run_che_external_projection import (
    exact_sign_flip_pvalue,
    finalize_replication_status,
    holm_adjust,
    paired_hodges_lehmann,
)


class RankProjectionTests(unittest.TestCase):
    def test_signature_with_top_ranked_genes_scores_higher(self) -> None:
        genes = [f"G{x}" for x in range(1, 501)]
        matrix = sparse.csr_matrix(
            np.vstack(
                [np.arange(500, 0, -1, dtype=float), np.arange(1, 501, dtype=float)]
            )
        )
        signatures = {program: genes[:50] for program in range(1, 11)}
        scores, coverage = rank_scores(matrix, pd.Index(genes), signatures, max_rank=100)
        self.assertGreater(scores.loc[0, 1], scores.loc[1, 1])
        self.assertTrue((coverage["present_genes"] == 50).all())

    def test_signature_coverage_below_45_fails_closed(self) -> None:
        genes = [f"G{x}" for x in range(1, 45)]
        matrix = sparse.csr_matrix(np.ones((1, len(genes))))
        signatures = {program: [f"G{x}" for x in range(1, 51)] for program in range(1, 11)}
        with self.assertRaisesRegex(ValueError, "below 45/50"):
            rank_scores(matrix, pd.Index(genes), signatures)

    def test_missing_gene_positions_are_compacted_before_scoring(self) -> None:
        genes = [f"G{x}" for x in range(2, 51)]
        matrix = sparse.csr_matrix(np.arange(49, 0, -1, dtype=float).reshape(1, -1))
        signatures = {program: [f"G{x}" for x in range(1, 51)] for program in range(1, 11)}
        scores, coverage = rank_scores(matrix, pd.Index(genes), signatures)
        self.assertEqual(scores.shape, (1, 10))
        self.assertTrue((coverage["present_genes"] == 49).all())

    def test_exact_sign_flip_small_sample_floor(self) -> None:
        self.assertEqual(exact_sign_flip_pvalue(np.ones(5)), 0.0625)

    def test_paired_hodges_lehmann_uses_walsh_averages(self) -> None:
        self.assertEqual(paired_hodges_lehmann(np.array([1.0, 2.0, 3.0])), 2.0)

    def test_holm_adjustment_is_monotone_in_sorted_pvalues(self) -> None:
        adjusted = holm_adjust(pd.Series([0.01, 0.04, 0.03]))
        self.assertTrue(np.allclose(adjusted.to_numpy(), [0.03, 0.06, 0.06]))

    def test_primary_nonreplication_is_not_hidden_by_sensitivity_label(self) -> None:
        self.assertEqual(
            finalize_replication_status(
                primary_status="not_replicated_or_opposite",
                eligible_for_confirmatory_che=True,
                sensitivity_direction_stable=False,
            ),
            "not_replicated_or_opposite",
        )

    def test_primary_directional_result_that_flips_is_sensitivity_unstable(self) -> None:
        self.assertEqual(
            finalize_replication_status(
                primary_status="directionally_consistent_weak",
                eligible_for_confirmatory_che=True,
                sensitivity_direction_stable=False,
            ),
            "sensitivity_unstable",
        )


if __name__ == "__main__":
    unittest.main()
