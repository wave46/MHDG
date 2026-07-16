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

from check_bundle import BundleError  # noqa: E402
from create_bundle import create_bundle  # noqa: E402
from prepare_run import prepare_run, render_parameter_file  # noqa: E402


PARAMETERS = """&INPUT_LST
    transport_model_path = '/old/transport_model.nml'
    field_path = '/old/equilibrium.h5' ! magnetic field
    jtor_path = '/old/current_density.h5'
    save_folder = '/old/output/'
/
&ADAPT_LST
    geometry_path = '/old/geometry.geo'
/
"""


class RunPreparationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.source = self.root / "source"
        self.source.mkdir()

        filenames = (
            "mesh.msh",
            "geometry.geo",
            "equilibrium.h5",
            "current_density.h5",
            "transport_model.nml",
            "restart.h5",
            "reference_mpi4_omp4.h5",
        )
        for filename in filenames:
            (self.source / filename).write_text(
                f"synthetic {filename}\n", encoding="utf-8"
            )
        (self.source / "param.txt").write_text(PARAMETERS, encoding="utf-8")

        self.bundle = self.root / "bundle"
        create_bundle(
            "legacy_fixed", self.source, self.bundle, REGRESSION_ROOT / "cases"
        )

        self.run_root = self.root / "runs"
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        self.serial_executable = self._executable(bin_dir / "serial_solver")
        self.parallel_executable = self._executable(bin_dir / "parallel_solver")
        self.mpi_launcher = self._executable(bin_dir / "mpirun")
        self.runtime_file = bin_dir / "positionFeketeNodesTri2D.h5"
        self.runtime_file.write_text("synthetic Fekete nodes\n", encoding="utf-8")

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

    @staticmethod
    def _executable(path: Path) -> Path:
        path.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
        path.chmod(0o755)
        return path

    def test_prepares_isolated_parallel_run(self) -> None:
        prepared = prepare_run(
            self.settings,
            "legacy_fixed",
            "warm",
            "mpi4_omp4",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "fixture",
        )

        expected = self.run_root / "legacy_fixed" / "warm" / "mpi4_omp4" / "fixture"
        self.assertEqual(prepared.path, expected)
        self.assertEqual(prepared.omp_threads, 4)
        self.assertEqual(
            prepared.command[:4],
            [str(self.mpi_launcher), "-n", "4", str(self.parallel_executable)],
        )
        self.assertEqual(
            prepared.command[-2:],
            [str(expected / "inputs" / "mesh"), str(expected / "inputs" / "restart")],
        )
        self.assertTrue((expected / "outputs").is_dir())
        self.assertTrue((expected / "inputs" / "equilibrium.h5").is_symlink())
        self.assertTrue((expected / "inputs" / "reference.h5").is_symlink())
        self.assertEqual(
            (expected / "positionFeketeNodesTri2D.h5").resolve(),
            self.runtime_file.resolve(),
        )

        parameters = (expected / "param.txt").read_text(encoding="utf-8")
        self.assertNotIn("/old/", parameters)
        rendered_paths = (
            expected / "inputs" / "transport_model.nml",
            expected / "inputs" / "equilibrium.h5",
            expected / "inputs" / "current_density.h5",
            expected / "inputs" / "geometry.geo",
        )
        for path in rendered_paths:
            self.assertIn(str(path), parameters)
        self.assertIn(f"{expected / 'outputs'}/", parameters)

        bundled_parameters = (
            self.bundle / "case_data" / "legacy_fixed" / "param.txt"
        ).read_text(encoding="utf-8")
        self.assertEqual(bundled_parameters, PARAMETERS)

    def test_public_command_prepares_serial_run(self) -> None:
        completed = subprocess.run(
            [
                str(REGRESSION_ROOT / "regression.sh"),
                "--settings",
                str(self.settings),
                "prepare",
                "legacy_fixed",
                "warm",
                "--layout",
                "serial_omp1",
                "--run-id",
                "public",
            ],
            check=False,
            capture_output=True,
            env={**os.environ, "PYTHON": sys.executable},
            text=True,
        )

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("OMP_NUM_THREADS=1", completed.stdout)
        plan_path = (
            self.run_root
            / "legacy_fixed"
            / "warm"
            / "serial_omp1"
            / "public"
            / "run_plan.json"
        )
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        self.assertEqual(plan["command"][0], str(self.serial_executable))
        self.assertNotIn(str(self.mpi_launcher), plan["command"])

    def test_missing_parameter_assignment_fails(self) -> None:
        source = self.root / "incomplete_param.txt"
        source.write_text(
            "&INPUT_LST\n  field_path = '/old/field.h5'\n/\n",
            encoding="utf-8",
        )
        destination = self.root / "rendered_param.txt"

        with self.assertRaisesRegex(BundleError, "must appear once"):
            render_parameter_file(
                source,
                destination,
                {
                    "field_path": "/new/field.h5",
                    "save_folder": "/new/output/",
                },
            )

        self.assertFalse(destination.exists())

    def test_existing_run_directory_is_not_replaced(self) -> None:
        run_dir = (
            self.run_root / "legacy_fixed" / "warm" / "mpi4_omp4" / "existing"
        )
        run_dir.mkdir(parents=True)
        marker = run_dir / "keep.txt"
        marker.write_text("keep\n", encoding="utf-8")

        with self.assertRaisesRegex(BundleError, "run directory already exists"):
            prepare_run(
                self.settings,
                "legacy_fixed",
                "warm",
                "mpi4_omp4",
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "layouts.json",
                "existing",
            )

        self.assertEqual(marker.read_text(encoding="utf-8"), "keep\n")

    def test_missing_runtime_file_fails_before_preparation(self) -> None:
        self.runtime_file.unlink()

        with self.assertRaisesRegex(BundleError, "required runtime file is missing"):
            prepare_run(
                self.settings,
                "legacy_fixed",
                "warm",
                "mpi4_omp4",
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "layouts.json",
                "missing-runtime-file",
            )

        run_dir = (
            self.run_root
            / "legacy_fixed"
            / "warm"
            / "mpi4_omp4"
            / "missing-runtime-file"
        )
        self.assertFalse(run_dir.exists())


if __name__ == "__main__":
    unittest.main()
