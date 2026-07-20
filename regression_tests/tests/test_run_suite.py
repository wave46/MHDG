from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from create_bundle import create_bundle  # noqa: E402


PARAMETERS = """&INPUT_LST
    transport_model_path = '/old/transport_model.nml'
    field_path = '/old/equilibrium.h5'
    jtor_path = '/old/current_density.h5'
    save_folder = '/old/output/'
/
&ADAPT_LST
    adaptivity = .true.
    rest_adapt = .true.
    geometry_path = '/old/geometry.geo'
/
"""

COLD_PARAMETER_FILES = (
    "param_cold_fixed_time_init.txt",
    "param_cold_fixed_diffusion_reduction.txt",
    *(f"param_cold_fixed_continuation_{index:02d}.txt" for index in range(1, 6)),
)
COLD_TRANSPORT_FILES = (
    "transport_cold_fixed_initial.nml",
    *(
        f"transport_cold_fixed_continuation_{index:02d}.nml"
        for index in range(1, 6)
    ),
)


class SuiteCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.run_root = self.root / "runs"
        self._create_bundle()
        self._create_executables()

        self.settings = self.root / "settings.env"
        self.settings.write_text(
            "MHDG_REGRESSION_SETTINGS_VERSION=1\n"
            f"MHDG_REGRESSION_DATA_ROOT={self.bundle}\n"
            f"MHDG_REGRESSION_RUN_ROOT={self.run_root}\n"
            f"MHDG_SERIAL_EXECUTABLE={self.serial_executable}\n"
            f"MHDG_PARALLEL_EXECUTABLE={self.parallel_executable}\n"
            f"MHDG_MPI_LAUNCHER={self.mpi_launcher}\n",
            encoding="utf-8",
        )

    def test_warm_suite_runs_compares_and_writes_summary(self) -> None:
        completed = self._suite("warm", "warm-test")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = self._summary("warm", "warm-test")
        self.assertEqual(summary["status"], "passed")
        self.assertEqual(len(summary["results"]), 1)
        result = summary["results"][0]
        self.assertEqual(result["layout_id"], "mpi4_omp4")
        self.assertEqual(result["workflow_id"], "warm")
        self.assertEqual(result["run_status"], "completed")
        self.assertEqual(result["comparison_status"], "passed")
        self.assertIn("mpi4_omp4", completed.stdout)

    def test_cold_matrix_run_only_saves_all_workflow_layout_results(self) -> None:
        self._executable(self.serial_executable, WORKFLOW_SOLVER)
        self._executable(self.parallel_executable, WORKFLOW_SOLVER)
        completed = self._suite("cold_matrix", "cold-test", "--run-only")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = self._summary("cold_matrix", "cold-test")
        self.assertEqual(summary["status"], "passed")
        self.assertEqual(summary["comparison_mode"], "deferred")
        self.assertEqual(summary["workflow_ids"], ["cold_fixed", "cold_adaptive"])
        self.assertNotIn("workflow_id", summary)
        self.assertEqual(len(summary["results"]), 10)
        self.assertEqual(
            {
                (result["workflow_id"], result["layout_id"])
                for result in summary["results"]
            },
            {
                (workflow, layout)
                for workflow in ("cold_fixed", "cold_adaptive")
                for layout in (
                    "serial_omp1",
                    "mpi2_omp1",
                    "mpi2_omp4",
                    "mpi4_omp1",
                    "mpi4_omp4",
                )
            },
        )
        self.assertTrue(
            all(result["status"] == "passed" for result in summary["results"])
        )
        self.assertTrue(
            all(
                result["comparison_status"] == "not_run"
                for result in summary["results"]
            )
        )
        self.assertIn("cold_adaptive  mpi4_omp4", completed.stdout)

        self._executable(self.serial_executable, FAILING_SOLVER)
        self._executable(self.parallel_executable, FAILING_SOLVER)
        resumed = self._suite(
            "cold_matrix", "cold-test", "--run-only", "--resume"
        )
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertEqual(
            len(self._summary("cold_matrix", "cold-test")["results"]), 10
        )
        self.assertIn("skipping recorded cold_fixed / serial_omp1", resumed.stdout)
        self.assertIn("suite passed:", completed.stdout)

    def test_warm_parallelism_suite_continues_after_one_layout_fails(self) -> None:
        self._executable(self.serial_executable, FAILING_SOLVER)
        completed = self._suite("warm_parallelism", "parallelism-test")

        self.assertEqual(completed.returncode, 1)
        summary = self._summary("warm_parallelism", "parallelism-test")
        self.assertEqual(summary["status"], "failed")
        self.assertEqual(len(summary["results"]), 5)
        self.assertEqual(summary["results"][0]["status"], "solver_failed")
        self.assertTrue(
            all(result["status"] == "passed" for result in summary["results"][1:])
        )
        self.assertIn("mpi4_omp4", completed.stdout)

    def test_golden_check_uses_default_settings_and_warm_suite(self) -> None:
        self._set_bundle_class("golden")
        completed = subprocess.run(
            [str(REGRESSION_ROOT / "regression.sh"), "golden-check"],
            check=False,
            capture_output=True,
            env={
                **os.environ,
                "PYTHON": sys.executable,
                "MHDG_REGRESSION_GOLDEN_SETTINGS": str(self.settings),
            },
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("suite: warm", completed.stdout)
        self.assertIn("mpi4_omp4", completed.stdout)

    def test_golden_check_rejects_candidate_bundle(self) -> None:
        completed = subprocess.run(
            [
                str(REGRESSION_ROOT / "regression.sh"),
                "--settings",
                str(self.settings),
                "golden-check",
            ],
            check=False,
            capture_output=True,
            env={**os.environ, "PYTHON": sys.executable},
            text=True,
        )

        self.assertEqual(completed.returncode, 1)
        self.assertIn("requires bundle_class=golden", completed.stderr)

    def _create_bundle(self) -> None:
        source = self.root / "source"
        source.mkdir()
        for filename in (
            "mesh.msh",
            "mesh_adaptive_initial.msh",
            "geometry.geo",
            "equilibrium.h5",
            "current_density.h5",
            "transport_model.nml",
            "restart.h5",
        ):
            (source / filename).write_text(f"synthetic {filename}\n", encoding="utf-8")
        (source / "param.txt").write_text(PARAMETERS, encoding="utf-8")
        for filename in COLD_PARAMETER_FILES:
            (source / filename).write_text(PARAMETERS, encoding="utf-8")
        for filename in COLD_TRANSPORT_FILES:
            (source / filename).write_text(
                f"synthetic {filename}\n", encoding="utf-8"
            )
        _write_solution(source / "reference_mpi4_omp4.h5")

        self.bundle = self.root / "bundle"
        create_bundle("legacy_case", source, self.bundle, REGRESSION_ROOT / "cases")

    def _create_executables(self) -> None:
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        self.serial_executable = self._executable(bin_dir / "serial", SOLVER)
        self.parallel_executable = self._executable(bin_dir / "parallel", SOLVER)
        self.mpi_launcher = self._executable(bin_dir / "mpirun", MPI_LAUNCHER)
        (bin_dir / "positionFeketeNodesTri2D.h5").write_text(
            "synthetic Fekete nodes\n", encoding="utf-8"
        )

    def _set_bundle_class(self, bundle_class: str) -> None:
        path = self.bundle / "manifest.json"
        manifest = json.loads(path.read_text(encoding="utf-8"))
        manifest["bundle_class"] = bundle_class
        path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    @staticmethod
    def _executable(path: Path, contents: str) -> Path:
        path.write_text(contents, encoding="utf-8")
        path.chmod(0o755)
        return path

    def _suite(
        self, suite: str, run_id: str, *arguments: str
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                str(REGRESSION_ROOT / "regression.sh"),
                "--settings",
                str(self.settings),
                "suite",
                suite,
                "--run-id",
                run_id,
                *arguments,
            ],
            check=False,
            capture_output=True,
            env={**os.environ, "PYTHON": sys.executable},
            text=True,
        )

    def _summary(self, suite: str, run_id: str) -> dict[str, object]:
        path = self.run_root / "suites" / suite / run_id / "suite_summary.json"
        return json.loads(path.read_text(encoding="utf-8"))


def _write_solution(path: Path) -> None:
    with h5py.File(path, "w") as handle:
        mesh = handle.create_group("mesh")
        mesh.create_dataset("X", data=np.array([[1.0, 2.0], [0.0, 1.0]]))
        mesh.create_dataset("T", data=np.array([[1], [2]], dtype=np.int32))
        mesh.create_dataset("Tlin", data=np.array([[1], [2]], dtype=np.int32))
        mesh.create_dataset("Tb", data=np.array([[1], [2]], dtype=np.int32))
        mesh.create_dataset("Nelems", data=np.array([1], dtype=np.int32))
        mesh.create_dataset("Nnodesperelem", data=np.array([2], dtype=np.int32))

        solution = handle.create_group("solution")
        solution.create_dataset("u", data=np.arange(1.0, 5.0))
        solution.create_dataset("q", data=np.arange(1.0, 9.0))
        solution.create_dataset("u_tilde", data=np.arange(1.0, 7.0))

        parameters = handle.create_group("simulation_parameters")
        parameters.create_dataset("Neq", data=np.array([2], dtype=np.int32))
        physics = parameters.create_group("physics")
        physics.create_dataset(
            "conservative_variable_names", data=np.array([b"rho", b"Gamma"])
        )


SOLVER = """#!/usr/bin/env bash
set -euo pipefail
cp inputs/reference.h5 outputs/result.h5
printf 'Error: 1.0E-5\n'
printf 'Output written to file outputs/result.h5\n'
"""

WORKFLOW_SOLVER = """#!/usr/bin/env bash
set -euo pipefail
if [[ -e inputs/reference.h5 ]]; then
  cp inputs/reference.h5 outputs/result.h5
else
  stage=${PWD##*/}
  if (($# == 1)); then
    history=$stage
  else
    history=$(<"$2.h5")
    history=${history%$'\\n'}">"$stage
  fi
  printf '%s\n' "$history" > outputs/result.h5
fi
printf 'Error: 1.0E-5\n'
printf 'Output written to file outputs/result.h5\n'
"""

FAILING_SOLVER = """#!/usr/bin/env bash
printf 'failed\n' >&2
exit 7
"""

MPI_LAUNCHER = """#!/usr/bin/env bash
set -euo pipefail
test "$1" = '--bind-to'
test "$2" = 'core'
test "$3" = '--map-by'
[[ "$4" == slot:PE=* ]]
shift 4
test "$1" = '-n'
shift 2
exec "$@"
"""


if __name__ == "__main__":
    unittest.main()
