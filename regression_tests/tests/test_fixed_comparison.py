from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from comparison.fixed.hdf5 import compare_hdf5_files  # noqa: E402
from tests.fixtures.harness import run_command  # noqa: E402
from tests.fixtures.solutions import write_solution  # noqa: E402


TOLERANCES = {
    "require_finite": True,
    "newton_error_max": 2.0e-4,
    "mesh_coordinate_atol": 1.0e-12,
    "relative_l2_max": 1.0e-10,
    "normalized_linf_max": 1.0e-9,
}


class Hdf5ComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.reference = self.root / "reference.h5"
        self.candidate = self.root / "candidate.h5"
        write_solution(self.reference, grouped=True)
        write_solution(self.candidate, grouped=False)

    def test_grouped_reference_matches_legacy_flat_candidate(self) -> None:
        with h5py.File(self.candidate, "r+") as handle:
            del handle["conservative_variable_names"]

        report = compare_hdf5_files(
            self.reference, self.candidate, TOLERANCES
        )

        self.assertEqual(report["status"], "passed")
        self.assertEqual(
            report["formats"]["candidate"], "flat_solution+flat_mesh"
        )
        self.assertEqual(
            report["solution"]["datasets"]["u"]["equations"]["Gamma"][
                "relative_l2"
            ],
            0.0,
        )
        self.assertEqual(report["solution"]["equation_names"], ["rho", "Gamma"])

    def test_per_equation_difference_fails(self) -> None:
        with h5py.File(self.candidate, "r+") as handle:
            handle["u"][1] += 1.0

        report = compare_hdf5_files(
            self.reference, self.candidate, TOLERANCES
        )

        self.assertEqual(report["status"], "failed")
        equations = report["solution"]["datasets"]["u"]["equations"]
        self.assertTrue(equations["rho"]["passed"])
        self.assertFalse(equations["Gamma"]["passed"])

    def test_nonfinite_candidate_fails(self) -> None:
        with h5py.File(self.candidate, "r+") as handle:
            handle["q"][0] = np.nan

        report = compare_hdf5_files(
            self.reference, self.candidate, TOLERANCES
        )

        self.assertEqual(report["status"], "failed")
        metrics = report["solution"]["datasets"]["q"]["equations"]["rho"]
        self.assertFalse(metrics["finite"])

    def test_connectivity_difference_fails(self) -> None:
        with h5py.File(self.candidate, "r+") as handle:
            handle["T"][0, 0] = 99

        report = compare_hdf5_files(
            self.reference, self.candidate, TOLERANCES
        )

        self.assertEqual(report["status"], "failed")
        self.assertFalse(report["mesh"]["connectivity"]["T"]["passed"])

    def test_magnetic_contract_matches(self) -> None:
        self._write_magnetic(self.reference)
        self._write_magnetic(self.candidate)

        report = compare_hdf5_files(
            self.reference, self.candidate, TOLERANCES
        )

        magnetic = report["magnetic"]
        self.assertEqual(report["status"], "passed")
        self.assertTrue(magnetic["present"])
        self.assertTrue(magnetic["passed"])
        self.assertEqual(magnetic["datasets"]["topology"]["mode"], "exact")
        self.assertEqual(
            magnetic["datasets"]["rho_pol_norm"]["relative_l2"],
            0.0,
        )

    def test_magnetic_region_difference_fails(self) -> None:
        self._write_magnetic(self.reference)
        self._write_magnetic(self.candidate)
        with h5py.File(self.candidate, "r+") as handle:
            handle["magnetic/topology_region"][2] = 3

        report = compare_hdf5_files(
            self.reference, self.candidate, TOLERANCES
        )

        self.assertEqual(report["status"], "failed")
        self.assertFalse(
            report["magnetic"]["datasets"]["topology_region"]["passed"]
        )
        self.assertIn("magnetic/topology_region differs", report["failures"])

    def test_missing_candidate_magnetic_contract_fails(self) -> None:
        self._write_magnetic(self.reference)

        report = compare_hdf5_files(
            self.reference, self.candidate, TOLERANCES
        )

        self.assertEqual(report["status"], "failed")
        self.assertIn("candidate is missing magnetic data", report["failures"])

    @staticmethod
    def _write_magnetic(path: Path) -> None:
        with h5py.File(path, "r+") as handle:
            magnetic = handle.create_group("magnetic")
            magnetic.create_dataset("topology", data=np.asarray([b"limited"]))
            magnetic.create_dataset("topology_id", data=np.asarray([1]))
            magnetic.create_dataset(
                "topology_region", data=np.asarray([1, 1, 2, 2])
            )
            magnetic.create_dataset(
                "rho_pol_norm", data=np.asarray([0.0, 0.5, 1.0, 1.2])
            )
            magnetic.create_dataset(
                "topology_normal",
                data=np.asarray(
                    [[1.0, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.0]]
                ),
            )


class CompareCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.run = self.root / "run"
        (self.run / "inputs").mkdir(parents=True)
        (self.run / "outputs").mkdir()

        self.reference = self.run / "inputs/reference.h5"
        self.checkpoint = self.run / "outputs/result_0000.h5"
        self.final = self.run / "outputs/result.h5"
        write_solution(self.reference, grouped=True)
        write_solution(self.checkpoint, grouped=False, solution_offset=1.0)
        write_solution(self.final, grouped=False)

        (self.run / "run_plan.json").write_text(
            json.dumps(
                {
                    "case_id": "legacy_case",
                    "workflow_id": "warm",
                    "layout_id": "mpi4_omp4",
                }
            ),
            encoding="utf-8",
        )
        (self.run / "run_metadata.json").write_text(
            json.dumps(
                {
                    "status": "completed",
                    "hdf5_outputs": [
                        "outputs/result.h5",
                        "outputs/result_0000.h5",
                    ],
                }
            ),
            encoding="utf-8",
        )
        (self.run / "stdout.log").write_text(
            f"Error: 8.0E-4\nError: 1.0E-5\n"
            f"Output written to file {self.checkpoint}\n"
            f"Output written to file {self.final}\n",
            encoding="utf-8",
        )

    def test_public_command_selects_final_output_and_writes_report(self) -> None:
        report_path = self.root / "comparison.json"
        completed = run_command(
            "compare",
            str(self.run),
            "--report",
            str(report_path),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["candidate"], str(self.final.resolve()))
        self.assertEqual(report["convergence"]["final_newton_error"], 1.0e-5)
        self.assertEqual(report["tolerance_profile"]["id"], "fixed_same_layout")
        self.assertEqual(len(report["files"]["reference"]["sha256"]), 64)
        self.assertIn("Newton error: PASS", completed.stdout)
        self.assertIn("Mesh: PASS", completed.stdout)
        self.assertIn("solution/u: PASS", completed.stdout)
        self.assertIn("transport_1d: PASS", completed.stdout)

    def test_newton_error_above_tolerance_fails(self) -> None:
        (self.run / "stdout.log").write_text(
            f"Error: 3.0E-4\nOutput written to file {self.final}\n",
            encoding="utf-8",
        )
        report_path = self.root / "failed-comparison.json"
        completed = run_command(
            "compare",
            str(self.run),
            "--report",
            str(report_path),
        )

        self.assertEqual(completed.returncode, 1)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "failed")
        self.assertFalse(report["convergence"]["passed"])
        self.assertIn("Newton error: FAIL", completed.stdout)
        self.assertIn(
            "FAIL: final Newton error exceeds tolerance",
            completed.stdout,
        )

    def test_nondefault_layout_uses_cross_layout_profile(self) -> None:
        (self.run / "run_plan.json").write_text(
            json.dumps(
                {
                    "case_id": "legacy_case",
                    "workflow_id": "warm",
                    "layout_id": "serial_omp16",
                }
            ),
            encoding="utf-8",
        )
        report_path = self.root / "cross-layout.json"
        completed = run_command(
            "compare",
            str(self.run),
            "--report",
            str(report_path),
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["tolerance_profile"]["id"], "fixed_cross_layout")

    def test_report_cannot_replace_candidate_hdf5(self) -> None:
        original_size = self.final.stat().st_size
        completed = run_command(
            "compare",
            str(self.run),
            "--report",
            str(self.final),
        )

        self.assertEqual(completed.returncode, 1)
        self.assertEqual(self.final.stat().st_size, original_size)
        with h5py.File(self.final, "r") as handle:
            self.assertIn("u", handle)
