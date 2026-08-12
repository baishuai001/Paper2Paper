from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd


HERE = Path(__file__).resolve().parent
FORMAL_CODE = Path(
    os.environ.get(
        "P2P_PHASE_ZERO_CODE",
        str(HERE),
    )
)
sys.path.insert(0, str(FORMAL_CODE))


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


builder = load("paired_builder", HERE / "build_crc_paired_cna_panel.py")
runner = load("configurable_runner", HERE / "run_crc_cna_panel_configurable.py")
combiner = load("paired_combiner", HERE / "combine_crc_paired_cna_summaries.py")
full_shards = load("full_shards", HERE / "build_full_cna_shards.py")


class PairedCnaHelpersTests(unittest.TestCase):
    def test_full_supervisor_exposes_formal_code_to_configurable_runner(self):
        script = (HERE / "run_full_cna_supervisor.sh").read_text(encoding="utf-8")
        self.assertIn('PYTHONPATH="$CODE" OMP_NUM_THREADS=4', script)

    def test_full_shards_include_only_preapproved_samples_exactly_once(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            panel = root / "panel.tsv"
            review = root / "review.tsv"
            pd.DataFrame(
                {
                    "panel_order": [1, 2, 3],
                    "sample_id": ["A", "B", "C"],
                    "analysis_patient_id": ["P1", "P1", "P2"],
                    "cancer_cells": [100, 300, 200],
                }
            ).to_csv(panel, sep="\t", index=False)
            pd.DataFrame(
                {
                    "panel_order": [1, 2, 3],
                    "sample_id": ["A", "B", "C"],
                    "reference_review": [
                        "pass_le_5pct",
                        "warning_gt_5pct_le_10pct",
                        "sensitivity_required_gt_10pct",
                    ],
                }
            ).to_csv(review, sep="\t", index=False)
            receipt = full_shards.freeze_shards(panel, review, root / "frozen", 2)
            self.assertEqual(receipt["samples_allowed"], 2)
            self.assertEqual(receipt["samples_excluded"], 1)
            orders = [
                order
                for shard in receipt["shards"]
                for order in shard["panel_orders"]
            ]
            self.assertEqual(sorted(orders), [1, 2])
            self.assertEqual(len(orders), len(set(orders)))

    def test_cell_evidence_allows_repeated_sample_rows_but_rejects_cross_shard_overlap(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shard_a = root / "a"
            shard_b = root / "b"
            shard_a.mkdir()
            shard_b.mkdir()
            filename = "cells.tsv"
            pd.DataFrame(
                {"panel_order": [1, 1, 2, 2], "cell_id": ["a", "b", "c", "d"]}
            ).to_csv(shard_a / filename, sep="\t", index=False)
            pd.DataFrame(
                {"panel_order": [3, 3], "cell_id": ["e", "f"]}
            ).to_csv(shard_b / filename, sep="\t", index=False)

            combined = combiner._load_shards(
                [shard_a, shard_b],
                filename,
                {"panel_order", "cell_id"},
                unique_panel_order_per_row=False,
            )
            self.assertEqual(len(combined), 6)

            pd.DataFrame(
                {"panel_order": [2, 3], "cell_id": ["g", "h"]}
            ).to_csv(shard_b / filename, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "overlap across shards"):
                combiner._load_shards(
                    [shard_a, shard_b],
                    filename,
                    {"panel_order", "cell_id"},
                    unique_panel_order_per_row=False,
                )

    def test_builder_keeps_all_samples_only_for_complete_pairs(self):
        rows = []
        for patient, states in {
            "P1": ["tumor", "metastasis", "metastasis"],
            "P2": ["tumor"],
        }.items():
            for i, state in enumerate(states):
                rows.append(
                    {
                        "dataset": "Liu_2024_mixCD45PosCD45Neg",
                        "sample_id": f"{patient}_{state}_{i}",
                        "analysis_patient_id": patient,
                        "sample_type": state,
                        "tissue": "liver" if state == "metastasis" else "colon",
                        "cells": 100,
                        "cancer_cells": 50,
                        "reference_cells": 50,
                        "other_epithelial_cells": 0,
                        "eligible": True,
                        "analysis_unit_id": f"{patient}_{state}_{i}::{patient}",
                    }
                )
        table = pd.DataFrame(rows)
        original = builder.COHORTS
        try:
            builder.COHORTS = {
                "Liu_2024_mixCD45PosCD45Neg": {
                    "analysis_role": "primary_discovery",
                    "contrast": "primary_vs_liver_metastasis",
                    "expected_pairs": 1,
                }
            }
            panel, summary = builder.build_panel(table)
        finally:
            builder.COHORTS = original
        self.assertEqual(set(panel["analysis_patient_id"]), {"P1"})
        self.assertEqual(len(panel), 3)
        self.assertEqual(summary.loc[0, "paired_patients"], 1)

    def test_configurable_runner_forwards_full_cell_limits(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            h5ad = root / "x.h5ad"
            h5ad.write_bytes(b"data")
            cna_runner = root / "run.R"
            cna_runner.write_text("x", encoding="utf-8")
            rlib = root / "rlib"
            rlib.mkdir()
            panel = root / "panel.tsv"
            pd.DataFrame(
                {
                    "panel_order": [1],
                    "dataset": ["D"],
                    "sample_id": ["S"],
                    "analysis_patient_id": ["P"],
                    "sample_type": ["tumor"],
                    "tissue": ["colon"],
                    "eligible": [True],
                }
            ).to_csv(panel, sep="\t", index=False)

            complete = {
                "command": ["fake"],
                "started_at": "x",
                "finished_at": "x",
                "exit_code": 0,
                "log": "fake",
                "log_bytes": 1,
                "log_sha256": "A" * 64,
            }
            with (
                patch.object(runner, "prepare_cna_input", return_value={}) as prepare,
                patch.object(runner, "run_command", return_value=complete),
                patch.object(runner, "summarize_copykat", return_value={}),
                patch.object(runner, "summarize_concordance", return_value={}),
            ):
                result = runner.run_panel(
                    h5ad,
                    panel,
                    root / "out",
                    cna_runner,
                    rlib,
                    "Rscript",
                    1,
                    7,
                    max_cancer=100000,
                    max_epithelial=100,
                    max_reference=300,
                    min_cancer=50,
                    min_reference=50,
                    reference_labels=["T cell", "B cell", "NK"],
                )
            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["cell_limits"]["max_cancer"], 100000)
            self.assertEqual(prepare.call_args.kwargs["max_cancer"], 100000)
            self.assertEqual(
                prepare.call_args.kwargs["reference_labels"], ["T cell", "B cell", "NK"]
            )

    def test_review_contract_marks_reference_sensitivity_without_losing_pair_identity(self):
        units = pd.DataFrame(
            {
                "panel_order": [1, 2],
                "dataset": ["D", "D"],
                "sample_id": ["primary", "met"],
                "analysis_patient_id": ["P", "P"],
                "paired_group": ["D::P", "D::P"],
                "contrast": ["primary_vs_liver_metastasis"] * 2,
                "analysis_role": ["primary_discovery"] * 2,
                "pair_member": ["primary", "liver_metastasis"],
                "author_cancer_cells": [300, 300],
                "author_cancer_two_method_support": [100, 80],
                "author_cancer_method_unresolved": [20, 30],
                "normal_reference_cells": [300, 300],
                "normal_reference_unexpected_malignancy_support": [3, 45],
            }
        )
        reviewed = combiner.classify_samples(units)
        self.assertEqual(reviewed.loc[0, "reference_review"], "pass_le_5pct")
        self.assertEqual(
            reviewed.loc[1, "reference_review"], "sensitivity_required_gt_10pct"
        )
        pairs = combiner.summarize_pairs(reviewed)
        self.assertEqual(pairs.loc[0, "pair_screen_status"], "pending_reference_sensitivity")


if __name__ == "__main__":
    unittest.main()
