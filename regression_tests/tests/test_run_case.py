from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


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
    geometry_path = '/old/geometry.geo'
/
"""


class RunCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self._create_bundle()

        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        self.serial_executable = self._executable(bin_dir / "serial", SOLVER)
        self.parallel_executable = self._executable(bin_dir / "parallel", SOLVER)
        self.mpi_launcher = self._executable(bin_dir / "mpirun", MPI_LAUNCHER)
        self.runtime_file = bin_dir / "positionFeketeNodesTri2D.h5"
        self.runtime_file.write_text("synthetic Fekete nodes\n", encoding="utf-8")
        self.environment_script = bin_dir / "environment setup.sh"
        self.environment_script.write_text(
            "export MHDG_TEST_ENV=loaded\n", encoding="utf-8"
        )
        self.build_manifest = bin_dir / "build_metadata.json"
        self.build_manifest.write_text("{}\n", encoding="utf-8")
        self.run_root = self.root / "runs"
        self.settings = self.root / "settings.env"
        self.settings.write_text(
            "MHDG_REGRESSION_SETTINGS_VERSION=1\n"
            f"MHDG_REGRESSION_DATA_ROOT={self.bundle}\n"
            f"MHDG_REGRESSION_RUN_ROOT={self.run_root}\n"
            f"MHDG_SERIAL_EXECUTABLE={self.serial_executable}\n"
            f"MHDG_PARALLEL_EXECUTABLE={self.parallel_executable}\n"
            f"MHDG_MPI_LAUNCHER={self.mpi_launcher}\n"
            f"MHDG_ENVIRONMENT_SCRIPT={self.environment_script}\n"
            f"MHDG_SOLVER_REVISION={'a' * 40}\n"
            "MHDG_BUILD_DESCRIPTION=synthetic-test-build\n"
            f"MHDG_BUILD_MANIFEST={self.build_manifest}\n",
            encoding="utf-8",
        )

    def test_parallel_run_records_logs_outputs_and_metadata(self) -> None:
        completed = self._run("mpi4_omp4", "parallel")
        self.assertEqual(completed.returncode, 0, completed.stderr)

        run_dir = self._run_dir("mpi4_omp4", "parallel")
        self.assertEqual(
            (run_dir / "stdout.log").read_text(encoding="utf-8"), "solver stdout\n"
        )
        self.assertEqual(
            (run_dir / "stderr.log").read_text(encoding="utf-8"), "solver stderr\n"
        )
        self.assertEqual(
            (run_dir / "outputs/mpi_ranks.txt").read_text(encoding="utf-8"),
            "4\n",
        )
        self.assertEqual(
            (run_dir / "outputs/omp_threads.txt").read_text(encoding="utf-8"),
            "4\n",
        )
        self.assertEqual(
            (run_dir / "outputs/environment.txt").read_text(encoding="utf-8"),
            "loaded\n",
        )
        self.assertEqual(
            (run_dir / "outputs/omp_places.txt").read_text(encoding="utf-8"),
            "cores\n",
        )
        self.assertEqual(
            (run_dir / "outputs/omp_proc_bind.txt").read_text(encoding="utf-8"),
            "spread\n",
        )

        metadata = json.loads(
            (run_dir / "run_metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["status"], "completed")
        self.assertEqual(metadata["exit_code"], 0)
        self.assertEqual(metadata["environment"]["OMP_NUM_THREADS"], "4")
        self.assertEqual(
            metadata["command"][:7],
            [
                str(self.mpi_launcher),
                "--bind-to",
                "core",
                "--map-by",
                "slot:PE=4",
                "-n",
                "4",
            ],
        )
        self.assertEqual(metadata["environment"]["OMP_PLACES"], "cores")
        self.assertEqual(metadata["environment"]["OMP_PROC_BIND"], "spread")
        self.assertEqual(metadata["executable"]["path"], str(self.parallel_executable))
        self.assertEqual(metadata["solver"]["revision"], "a" * 40)
        self.assertEqual(
            metadata["solver"]["build_manifest"]["path"],
            str(self.build_manifest),
        )
        self.assertEqual(
            metadata["runtime_files"]["positionFeketeNodesTri2D.h5"]["path"],
            str(self.runtime_file),
        )
        self.assertEqual(
            metadata["environment"]["setup_script"]["path"],
            str(self.environment_script),
        )
        self.assertEqual(metadata["hdf5_outputs"], ["outputs/result.h5"])
        self.assertEqual(len(metadata["executable"]["sha256"]), 64)

    def test_serial_run_without_hdf5_is_reported_as_incomplete(self) -> None:
        self._executable(self.serial_executable, NO_OUTPUT_SOLVER)
        completed = self._run("serial_omp1", "missing-output")
        self.assertEqual(completed.returncode, 1)

        run_dir = self._run_dir("serial_omp1", "missing-output")
        metadata = json.loads(
            (run_dir / "run_metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["status"], "missing_hdf5_output")
        self.assertEqual(metadata["exit_code"], 0)
        self.assertEqual(metadata["command"][0], str(self.serial_executable))
        self.assertEqual(metadata["environment"]["OMP_NUM_THREADS"], "1")

    def test_solver_failure_keeps_logs_and_exit_code(self) -> None:
        self._executable(self.serial_executable, FAILING_SOLVER)
        completed = self._run("serial_omp1", "solver-failure")
        self.assertEqual(completed.returncode, 1)

        run_dir = self._run_dir("serial_omp1", "solver-failure")
        metadata = json.loads(
            (run_dir / "run_metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["status"], "solver_failed")
        self.assertEqual(metadata["exit_code"], 7)
        self.assertEqual(
            (run_dir / "stderr.log").read_text(encoding="utf-8"), "failed\n"
        )

    def _create_bundle(self) -> None:
        source = self.root / "source"
        source.mkdir()
        for filename in (
            "mesh.msh",
            "geometry.geo",
            "equilibrium.h5",
            "current_density.h5",
            "transport_model.nml",
            "restart.h5",
            "reference_mpi4_omp4.h5",
        ):
            (source / filename).write_text(
                f"synthetic {filename}\n", encoding="utf-8"
            )
        (source / "param.txt").write_text(PARAMETERS, encoding="utf-8")

        self.bundle = self.root / "bundle"
        create_bundle(
            "legacy_fixed", source, self.bundle, REGRESSION_ROOT / "cases"
        )

    @staticmethod
    def _executable(path: Path, contents: str) -> Path:
        path.write_text(contents, encoding="utf-8")
        path.chmod(0o755)
        return path

    def _run(self, layout: str, run_id: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                str(REGRESSION_ROOT / "regression.sh"),
                "--settings",
                str(self.settings),
                "run",
                "legacy_fixed",
                "warm",
                "--layout",
                layout,
                "--run-id",
                run_id,
            ],
            check=False,
            capture_output=True,
            env={**os.environ, "PYTHON": sys.executable},
            text=True,
        )

    def _run_dir(self, layout: str, run_id: str) -> Path:
        return self.run_root / "legacy_fixed" / "warm" / layout / run_id


SOLVER = """#!/usr/bin/env bash
set -euo pipefail
printf 'solver stdout\n'
printf 'solver stderr\n' >&2
printf '%s\n' "$OMP_NUM_THREADS" > outputs/omp_threads.txt
printf '%s\n' "$OMP_PLACES" > outputs/omp_places.txt
printf '%s\n' "$OMP_PROC_BIND" > outputs/omp_proc_bind.txt
printf '%s\n' "$MHDG_TEST_ENV" > outputs/environment.txt
printf 'synthetic hdf5\n' > outputs/result.h5
"""

NO_OUTPUT_SOLVER = """#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$OMP_NUM_THREADS" > outputs/omp_threads.txt
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
test "$4" = 'slot:PE=4'
shift 4
test "$1" = '-n'
printf '%s\n' "$2" > outputs/mpi_ranks.txt
shift 2
exec "$@"
"""


if __name__ == "__main__":
    unittest.main()
