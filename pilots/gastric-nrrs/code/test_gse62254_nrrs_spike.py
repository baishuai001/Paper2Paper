from __future__ import annotations

import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).with_name("gse62254_nrrs_spike.py")
SPEC = importlib.util.spec_from_file_location("gse62254_nrrs_spike", MODULE_PATH)
assert SPEC and SPEC.loader
spike = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = spike
SPEC.loader.exec_module(spike)


class GSE62254SpikeTests(unittest.TestCase):
    def test_metadata_path_never_exposes_an_absolute_directory(self) -> None:
        recorded = spike.metadata_path(MODULE_PATH.resolve(), MODULE_PATH.parent)
        self.assertEqual(recorded, MODULE_PATH.name)
        self.assertFalse(Path(recorded).is_absolute())

    def test_resource_verification_fails_closed_on_change(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "resource.bin"
            path.write_bytes(b"audited input")
            resource = spike.Resource(
                resource_id="TEST",
                url="https://example.org/resource.bin",
                local_name="resource.bin",
                expected_bytes=path.stat().st_size,
                sha256=spike.sha256_file(path),
            )
            self.assertTrue(spike.verify_resource(path, resource)["verified"])
            path.write_bytes(b"changed input")
            with self.assertRaises(ValueError):
                spike.verify_resource(path, resource)

    def test_nrrs_uses_sample_sd_and_published_coefficient_order(self) -> None:
        values = pd.DataFrame(
            {
                gene: np.arange(1, 6, dtype=float) * (index + 1)
                for index, gene in enumerate(spike.TARGET_GENES)
            },
            index=[f"GSM{index}" for index in range(5)],
        )
        zscores, score = spike.compute_nrrs(values)
        self.assertTrue(np.allclose(zscores.mean(axis=0), 0.0))
        self.assertTrue(np.allclose(zscores.std(axis=0, ddof=1), 1.0))
        expected = zscores.mul(pd.Series(spike.COEFFICIENTS), axis=1).sum(axis=1)
        self.assertTrue(np.allclose(score, expected))

    def test_probe_aggregation_is_outcome_blind_and_deterministic(self) -> None:
        extra_genes = list(spike.TARGET_GENES[2:])
        expression = pd.DataFrame(
            [
                [1, 2, 3, 4],
                [1, 3, 5, 7],
                [5, 5, 5, 5],
                *[[index, index + 1, index + 2, index + 3] for index in range(6)],
            ],
            index=["p1", "p2", "p3", *[f"extra{index}" for index in range(6)]],
            columns=["s1", "s2", "s3", "s4"],
        )
        mapping = pd.DataFrame(
            {
                "probe_id": ["p1", "p2", "p3", *[f"extra{index}" for index in range(6)]],
                "gene": ["AGT", "AGT", "EPHB3", *extra_genes],
            }
        )
        mean, _ = spike.aggregate_gene_expression(expression, mapping, "mean_probes")
        self.assertTrue(np.allclose(mean["AGT"], [1, 2.5, 4, 5.5]))
        maximum, selected = spike.aggregate_gene_expression(
            expression, mapping, "max_variance_probe"
        )
        self.assertTrue(np.allclose(maximum["AGT"], expression.loc["p2"]))
        self.assertEqual(
            selected.loc[selected["gene"] == "AGT"].query("selected")["probe_id"].tolist(),
            ["p2"],
        )

    def test_geo_metadata_requires_complete_unique_identifiers(self) -> None:
        sample_ids = [f"GSM{index:04d}" for index in range(300)]
        titles = [f"T{index}" for index in range(300)]
        patient_ids = [str(index) for index in range(300)]
        lines = [
            "!Sample_title\t" + "\t".join(f'\"{value}\"' for value in titles),
            "!Sample_geo_accession\t"
            + "\t".join(f'\"{value}\"' for value in sample_ids),
            "!Sample_characteristics_ch1\t"
            + "\t".join(f'\"patient: {value}\"' for value in patient_ids),
            "!series_matrix_table_begin",
        ]
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "matrix.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                handle.write("\n".join(lines) + "\n")
            observed = spike.parse_geo_matrix_metadata(path)
        self.assertEqual(len(observed), 300)
        self.assertEqual(observed.iloc[0].to_dict(), {
            "gsm_id": "GSM0000", "sample_title": "T0", "patient_id": "0"
        })

    def test_median_grouping_is_balanced_with_ties(self) -> None:
        score = pd.Series(
            [0.0, 1.0, 1.0, 2.0], index=["d", "c", "b", "a"], name="nrrs"
        )
        groups, cutoff = spike.median_groups(score)
        self.assertEqual(cutoff, 1.0)
        self.assertEqual((groups == "high").sum(), 2)
        self.assertEqual((groups == "low").sum(), 2)
        self.assertEqual(groups["b"], "high")

    def test_kaplan_meier_plot_is_noninteractive_and_written(self) -> None:
        frame = pd.DataFrame(
            {
                "os_months": [5.0, 8.0, 10.0, 12.0, 15.0, 20.0],
                "death": [1, 1, 0, 1, 0, 0],
                "risk_group": ["high", "high", "high", "low", "low", "low"],
            }
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "km.png"
            spike.plot_kaplan_meier(frame, output, 0.05)
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 1000)

    def test_expected_result_contract_passes_and_fails_closed(self) -> None:
        observed = {
            "main": {"hr": 1.5, "direction": True},
            "sensitivity": {"rho": 0.97},
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "expected.json"
            path.write_text(
                json.dumps(
                    {
                        "main": {"hr": 1.5, "direction": True},
                        "sensitivity": {"rho": 0.97},
                    }
                ),
                encoding="utf-8",
            )
            report = spike.verify_expected_results(observed, path)
            self.assertTrue(report["verified"])
            self.assertEqual(report["metrics_checked"], 3)
            path.write_text(json.dumps({"main": {"hr": 9.0}}), encoding="utf-8")
            with self.assertRaises(ValueError):
                spike.verify_expected_results(observed, path)


if __name__ == "__main__":
    unittest.main()
