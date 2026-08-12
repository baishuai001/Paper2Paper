from __future__ import annotations

import unittest

import numpy as np

from gate1_common import bh_fdr


class GateCommonTests(unittest.TestCase):
    def test_bh_is_monotone_in_rank_order(self) -> None:
        values = np.array([0.04, 0.001, 0.03, 0.20])
        adjusted = bh_fdr(values)
        order = np.argsort(values)
        self.assertTrue(np.all(np.diff(adjusted[order]) >= -1e-12))
        self.assertTrue(np.all(adjusted >= values - 1e-12))


if __name__ == "__main__":
    unittest.main()
