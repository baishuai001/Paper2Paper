import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from phenotype_gate import cliffs_delta, residualize_by_dataset
from score_ulm import ulm
from tf_statistics import hedges_g, random_effects


class GateStatisticsTests(unittest.TestCase):
    def test_cliffs_delta_direction(self):
        self.assertEqual(cliffs_delta(np.array([3, 4]), np.array([1, 2])), 1.0)

    def test_dataset_residual_means_are_zero(self):
        x = np.array([[1.0], [3.0], [10.0], [14.0]])
        dataset = pd.Series(["a", "a", "b", "b"])
        residual = residualize_by_dataset(x, dataset)
        self.assertTrue(np.allclose([residual[:2].mean(), residual[2:].mean()], 0.0, atol=1e-10))

    def test_ulm_signed_direction(self):
        genes = np.array([f"g{i}" for i in range(10)])
        edges = pd.DataFrame({"source": ["TF"] * 10, "target": genes, "sign": [1] * 5 + [-1] * 5})
        expression = np.array([[3, 4, 5, 4, 3, 0, 1, 0, 1, 0], [0, 1, 0, 1, 0, 3, 4, 5, 4, 3]], dtype=float)
        score, tfs, _ = ulm(expression, genes, edges)
        self.assertEqual(tfs, ["TF"])
        self.assertGreater(score[0, 0], 0)
        self.assertLess(score[1, 0], 0)

    def test_random_effects_direction(self):
        pooled, _, p, i2, _ = random_effects(np.array([0.7, 0.8, 0.9]), np.array([0.04, 0.04, 0.04]))
        self.assertGreater(pooled, 0)
        self.assertLess(p, 0.01)
        self.assertGreaterEqual(i2, 0)

    def test_hedges_g_direction(self):
        g, variance = hedges_g(np.array([2.0, 3.0, 4.0]), np.array([0.0, 1.0, 2.0]))
        self.assertGreater(g, 0)
        self.assertGreater(variance, 0)


if __name__ == "__main__":
    unittest.main()
