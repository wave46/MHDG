from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from comparison.adaptive import (  # noqa: E402
    SampledFields,
    compare_adaptive_run,
    compare_sampled_fields,
    reference_sample_points,
)


class AdaptiveComparisonTests(unittest.TestCase):
    @patch("comparison.adaptive.compare_adaptive_files")
    def test_run_comparison_adds_profile_and_convergence(self, compare_files) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            (run / "inputs").mkdir()
            for path in (
                run / "inputs/reference.h5",
                run / "candidate.h5",
                run / "positionFeketeNodesTri2D.h5",
            ):
                path.write_text("synthetic\n", encoding="utf-8")
            (run / "stdout.log").write_text("Error: 1.0E-4\n", encoding="utf-8")
            (run / "run_plan.json").write_text(
                json.dumps(
                    {
                        "case_id": "legacy_case",
                        "workflow_id": "cold_adaptive",
                        "layout_id": "serial_omp1",
                    }
                ),
                encoding="utf-8",
            )
            (run / "run_metadata.json").write_text(
                json.dumps(
                    {
                        "status": "completed",
                        "hdf5_outputs": ["candidate.h5"],
                        "runtime_files": {
                            "positionFeketeNodesTri2D.h5": {
                                "path": str(run / "positionFeketeNodesTri2D.h5")
                            }
                        },
                    }
                ),
                encoding="utf-8",
            )
            compare_files.return_value = {
                "schema_version": 1,
                "status": "passed",
                "failures": [],
                "tolerances": {},
            }

            path, report = compare_adaptive_run(
                run,
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "tolerances.json",
            )

        self.assertEqual(report["status"], "passed")
        self.assertTrue(report["convergence"]["passed"])
        self.assertEqual(report["tolerance_profile"]["id"], "adaptive_reference")
        self.assertEqual(path.name, "comparison.json")
        self.assertEqual(compare_files.call_args.args[3], 4)

    @patch("comparison.adaptive.compare_adaptive_files")
    def test_run_comparison_accepts_stage_file_overrides(self, compare_files) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            candidate = run / "candidate.h5"
            reference = run / "golden.h5"
            fekete = run / "positionFeketeNodesTri2D.h5"
            for path in (candidate, reference, fekete):
                path.write_text("synthetic\n", encoding="utf-8")
            (run / "stdout.log").write_text("Error: 1.0E-4\n", encoding="utf-8")
            (run / "run_plan.json").write_text(
                json.dumps(
                    {
                        "case_id": "legacy_case",
                        "workflow_id": "cold_adaptive",
                        "layout_id": "mpi4_omp4",
                    }
                ),
                encoding="utf-8",
            )
            (run / "run_metadata.json").write_text(
                json.dumps(
                    {
                        "status": "completed",
                        "hdf5_outputs": ["candidate.h5"],
                        "runtime_files": {
                            "positionFeketeNodesTri2D.h5": {"path": str(fekete)}
                        },
                    }
                ),
                encoding="utf-8",
            )
            compare_files.return_value = {
                "schema_version": 1,
                "status": "passed",
                "failures": [],
                "tolerances": {},
            }

            _, report = compare_adaptive_run(
                run,
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "tolerances.json",
                candidate_override=candidate,
                reference_override=reference,
                tolerance_profile_override="adaptive_reference",
            )

        self.assertEqual(report["status"], "passed")
        self.assertEqual(compare_files.call_args.args[:2], (reference, candidate))

    def test_reference_points_are_inside_each_triangle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "reference.h5"
            with h5py.File(path, "w") as handle:
                mesh = handle.create_group("mesh")
                mesh["X"] = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
                mesh["Tlin"] = np.array([[1], [2], [3]])

            points = reference_sample_points(path, 4)

        self.assertEqual(points.shape, (4, 2))
        self.assertTrue((points > 0).all())
        self.assertTrue((points.sum(axis=1) < 1).all())
        np.testing.assert_allclose(points[0], [1 / 3, 1 / 3])

    def test_identical_samples_are_characterized_without_tolerances(self) -> None:
        sampled = self._sampled()

        report = compare_sampled_fields(sampled, sampled, 4)

        self.assertEqual(report["status"], "characterized")
        self.assertEqual(report["sampling"]["common_coverage"], 1.0)
        for dataset in report["datasets"].values():
            for metrics in dataset["equations"].values():
                self.assertEqual(metrics["relative_l2"], 0.0)
                self.assertIsNone(metrics["passed"])

    def test_tolerances_and_coverage_produce_failure(self) -> None:
        reference = self._sampled()
        candidate = self._sampled(offset=0.1, inside=[True, True, False])
        tolerances = {
            "minimum_point_coverage": 1.0,
            "solution": {
                "relative_l2_max": 1.0,
                "normalized_linf_max": 1.0,
            },
            "gradient": {
                "relative_l2_max": 1e-3,
                "normalized_linf_max": 1e-3,
            },
        }

        report = compare_sampled_fields(reference, candidate, 1, tolerances)

        self.assertEqual(report["status"], "failed")
        self.assertAlmostEqual(report["sampling"]["common_coverage"], 2 / 3)
        self.assertIn("common point coverage is below tolerance", report["failures"])
        self.assertTrue(report["datasets"]["solution"]["equations"]["rho"]["passed"])
        self.assertTrue(
            any("gradient_x/rho exceeds tolerance" in item for item in report["failures"])
        )

    @staticmethod
    def _sampled(
        offset: float = 0.0, inside: list[bool] | None = None
    ) -> SampledFields:
        solution = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]) + offset
        gradient = np.stack((solution * 0.1, solution * 0.2), axis=2)
        return SampledFields(
            ["rho", "Gamma"],
            np.array(inside if inside is not None else [True] * 3),
            solution,
            gradient,
        )


if __name__ == "__main__":
    unittest.main()
