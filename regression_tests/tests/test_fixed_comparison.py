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
    "mesh_connectivity": "exact",
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


class NumberingInvariantComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.reference = self.root / "reference.h5"
        self.candidate = self.root / "candidate.h5"
        write_solution(self.reference, grouped=True)
        write_solution(self.candidate, grouped=True)
        _permute_solution_numbering(self.candidate)
        self.tolerances = {
            **TOLERANCES,
            "mesh_connectivity": "numbering_invariant",
            "mesh_coordinate_atol": 1.0e-8,
        }

    def test_node_element_and_face_renumbering_matches(self) -> None:
        report = compare_hdf5_files(
            self.reference,
            self.candidate,
            self.tolerances,
        )

        self.assertEqual(report["status"], "passed", report["failures"])
        self.assertEqual(report["mesh"]["mode"], "numbering_invariant")
        self.assertTrue(report["mesh"]["coordinates"]["passed"])
        topology = report["mesh"]["connectivity"]["topology"]
        self.assertTrue(topology["passed"])
        self.assertEqual(topology["reversed_face_count"], 1)
        for dataset in ("u", "q", "u_tilde"):
            details = report["solution"]["datasets"][dataset]
            self.assertTrue(details["mesh_alignment_applied"])
            self.assertTrue(details["passed"])

    def test_different_element_topology_still_fails(self) -> None:
        with h5py.File(self.candidate, "r+") as handle:
            handle["mesh/T"][2, 0] = 3

        report = compare_hdf5_files(
            self.reference,
            self.candidate,
            self.tolerances,
        )

        self.assertEqual(report["status"], "failed")
        topology = report["mesh"]["connectivity"]["topology"]
        self.assertFalse(topology["passed"])
        self.assertIn("elements differ", topology["reason"])
        self.assertEqual(
            report["solution"]["reason"], "mesh comparison failed"
        )

    def test_changed_element_local_ordering_fails(self) -> None:
        with h5py.File(self.candidate, "r+") as handle:
            for name in ("T", "Tlin"):
                handle[f"mesh/{name}"][:, 0] = np.roll(
                    handle[f"mesh/{name}"][:, 0], -1
                )

        report = compare_hdf5_files(
            self.reference,
            self.candidate,
            self.tolerances,
        )

        topology = report["mesh"]["connectivity"]["topology"]
        self.assertEqual(report["status"], "failed")
        self.assertIn("local node ordering", topology["reason"])


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


def _permute_solution_numbering(path: Path) -> None:
    reference_node_to_candidate = np.array([3, 1, 4, 2], dtype=np.int32)
    candidate_element_order = np.array([1, 0])
    candidate_boundary_order = np.array([2, 0, 3, 1])

    with h5py.File(path, "r+") as handle:
        mesh = handle["mesh"]
        coordinates = np.asarray(mesh["X"]).T
        candidate_coordinates = np.empty_like(coordinates)
        candidate_coordinates[reference_node_to_candidate - 1] = coordinates
        mesh["X"][...] = candidate_coordinates.T

        for name in ("T", "Tlin"):
            nodes = np.asarray(mesh[name]).T
            candidate_nodes = reference_node_to_candidate[
                nodes[candidate_element_order] - 1
            ]
            mesh[name][...] = candidate_nodes.T

        boundary_nodes = np.asarray(mesh["Tb"]).T
        mesh["Tb"][...] = reference_node_to_candidate[
            boundary_nodes[candidate_boundary_order] - 1
        ].T
        boundary_flags = np.asarray(mesh["boundaryFlag"])
        mesh["boundaryFlag"][...] = boundary_flags[candidate_boundary_order]
        mesh["intfaces"][...] = np.array([[1, 3, 2, 2, 2]]).T
        mesh["extfaces"][...] = np.array(
            [[1, 2], [2, 1], [2, 3], [1, 1]]
        ).T

        solution = handle["solution"]
        equation_count = int(handle["simulation_parameters/Neq"][0])
        for name, shape in {
            "u": (2, 3, equation_count),
            "q": (2, 3, equation_count, 2),
        }.items():
            values = np.asarray(solution[name]).reshape(shape)
            solution[name][...] = values[candidate_element_order].reshape(-1)

        traces = np.asarray(solution["u_tilde"]).reshape(5, 2, equation_count)
        candidate_traces = np.empty_like(traces)
        candidate_traces[0] = traces[0, ::-1]
        candidate_traces[1:] = traces[1 + candidate_boundary_order]
        solution["u_tilde"][...] = candidate_traces.reshape(-1)
