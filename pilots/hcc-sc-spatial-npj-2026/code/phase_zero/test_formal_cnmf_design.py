from __future__ import annotations

import unittest

import pandas as pd

from formal_cnmf_design import (
    ONE_METHOD_CLASS,
    SELECTED_CLASS,
    assign_patient_folds,
    balance_patient_state_cells,
    select_cohort_cells,
)


def evidence_fixture() -> pd.DataFrame:
    rows = []
    order = 0
    for patient in ("P1", "P2", "P3", "P4"):
        for state in ("tumor", "metastasis"):
            order += 1
            sample = f"{patient}_{state}"
            for i in range(40):
                rows.append(
                    {
                        "cell_id": f"{sample}_c{i}",
                        "panel_order": order,
                        "panel_dataset": "Liu",
                        "panel_sample_id": sample,
                        "panel_patient_id": patient,
                        "panel_sample_type": state,
                        "panel_tissue": "colon" if state == "tumor" else "liver",
                        "cna_input_role": "author_cancer",
                        "copykat.pred": "aneuploid",
                        "scevan_class": "tumor",
                        "method_unresolved": False,
                        "method_unresolved_bool": False,
                        "cna_selection_class": SELECTED_CLASS,
                    }
                )
            for i in range(2):
                rows.append(
                    {
                        "cell_id": f"{sample}_r{i}",
                        "panel_order": order,
                        "panel_dataset": "Liu",
                        "panel_sample_id": sample,
                        "panel_patient_id": patient,
                        "panel_sample_type": state,
                        "panel_tissue": "colon" if state == "tumor" else "liver",
                        "cna_input_role": "known_normal_reference",
                        "copykat.pred": "diploid",
                        "scevan_class": "normal",
                        "method_unresolved": False,
                        "method_unresolved_bool": False,
                        "cna_selection_class": "known_normal_reference_no_malignancy_support",
                    }
                )
    return pd.DataFrame(rows)


class FormalCnmfDesignTests(unittest.TestCase):
    def test_cohort_selection_requires_paired_states_and_dual_calls(self) -> None:
        selected, sample_qc = select_cohort_cells(
            evidence_fixture(),
            "Liu",
            expected_patients=4,
            expected_samples=8,
            expected_dual_cells=320,
        )
        self.assertEqual(len(selected), 320)
        self.assertEqual(len(sample_qc), 8)
        self.assertEqual(set(selected["cna_selection_class"]), {SELECTED_CLASS})

    def test_two_frozen_sensitivity_selections_are_explicit(self) -> None:
        fixture = evidence_fixture()
        one_method = fixture.loc[
            fixture["cna_selection_class"].eq(SELECTED_CLASS)
        ].head(8).copy()
        one_method["cell_id"] = "one_method_" + one_method["cell_id"]
        one_method["copykat.pred"] = "diploid"
        one_method["scevan_class"] = "tumor"
        one_method["cna_selection_class"] = ONE_METHOD_CLASS
        fixture = pd.concat([fixture, one_method], ignore_index=True)

        expanded, _ = select_cohort_cells(
            fixture,
            "Liu",
            selected_classes=(SELECTED_CLASS, ONE_METHOD_CLASS),
        )
        self.assertEqual(len(expanded), 328)
        self.assertEqual(
            set(expanded["cna_selection_class"]),
            {SELECTED_CLASS, ONE_METHOD_CLASS},
        )

        excluded, _ = select_cohort_cells(
            fixture,
            "Liu",
            excluded_sample_ids=("P1_tumor",),
            require_paired_states=False,
        )
        self.assertEqual(len(excluded), 280)
        self.assertNotIn("P1_tumor", set(excluded["panel_sample_id"]))

    def test_one_method_class_requires_exactly_one_supporting_call(self) -> None:
        fixture = evidence_fixture()
        invalid = fixture.loc[
            fixture["cna_selection_class"].eq(SELECTED_CLASS)
        ].head(1).copy()
        invalid["cell_id"] = "invalid_one_method"
        invalid["cna_selection_class"] = ONE_METHOD_CLASS
        fixture = pd.concat([fixture, invalid], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "recorded CNA support class"):
            select_cohort_cells(
                fixture,
                "Liu",
                selected_classes=(SELECTED_CLASS, ONE_METHOD_CLASS),
            )

    def test_patient_state_cap_is_deterministic_and_not_per_sample(self) -> None:
        fixture = evidence_fixture()
        selected, _ = select_cohort_cells(fixture, "Liu")
        # Add a second technical sample to P1 tumor. The cap must still apply to
        # the combined P1-tumor unit rather than once per sample.
        extra = selected.loc[
            selected["panel_patient_id"].eq("P1")
            & selected["panel_sample_type"].eq("tumor")
        ].copy()
        extra["cell_id"] = "extra_" + extra["cell_id"]
        extra["panel_sample_id"] = "P1_tumor_repeat"
        expanded = pd.concat([selected, extra], ignore_index=True)
        first, units, composition = balance_patient_state_cells(expanded, 30, 12)
        second, _, _ = balance_patient_state_cells(expanded, 30, 12)
        p1_tumor = first.loc[
            first["panel_patient_id"].eq("P1")
            & first["panel_sample_type"].eq("tumor")
        ]
        self.assertEqual(len(p1_tumor), 30)
        self.assertEqual(first["cell_id"].tolist(), second["cell_id"].tolist())
        self.assertEqual(units["selected_cells"].max(), 30)
        self.assertIn("P1_tumor_repeat", set(composition["panel_sample_id"]))

    def test_patient_holdouts_keep_pairs_together(self) -> None:
        selected, _ = select_cohort_cells(evidence_fixture(), "Liu")
        obs = selected.rename(
            columns={"panel_patient_id": "analysis_patient_id", "panel_sample_type": "sample_type"}
        )
        manifest = assign_patient_folds(obs, 20260812)
        self.assertEqual(set(manifest["holdout_fold"]), {"A", "B"})
        self.assertEqual(manifest["analysis_patient_id"].nunique(), 4)
        self.assertTrue((manifest[["tumor", "metastasis"]] > 0).all().all())

    def test_holdout_rejects_unpaired_patient(self) -> None:
        selected, _ = select_cohort_cells(evidence_fixture(), "Liu")
        obs = selected.rename(
            columns={"panel_patient_id": "analysis_patient_id", "panel_sample_type": "sample_type"}
        )
        obs = obs.loc[
            ~(
                obs["analysis_patient_id"].eq("P1")
                & obs["sample_type"].eq("metastasis")
            )
        ]
        with self.assertRaisesRegex(ValueError, "paired states"):
            assign_patient_folds(obs, 20260812)

if __name__ == "__main__":
    unittest.main()
