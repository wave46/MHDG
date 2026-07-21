from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from support.errors import BundleError  # noqa: E402
from preparation.models import PreparedStagedRun  # noqa: E402
from preparation.parameters import render_parameter_file  # noqa: E402
from prepare_run import prepare_run  # noqa: E402
from tests.fixtures.case_data import PARAMETERS  # noqa: E402
from tests.fixtures.harness import create_harness  # noqa: E402


class RunPreparationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        fixture = create_harness(self.root)
        self.bundle = fixture.bundle
        self.run_root = fixture.run_root
        self.settings = fixture.settings
        self.parallel_executable = fixture.parallel_executable
        self.mpi_launcher = fixture.mpi_launcher
        self.runtime_file = fixture.runtime_file

    def test_prepares_isolated_parallel_run(self) -> None:
        prepared = prepare_run(
            self.settings,
            "legacy_case",
            "warm",
            "mpi4_omp4",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "fixture",
        )

        expected = self.run_root / "legacy_case" / "warm" / "mpi4_omp4" / "fixture"
        self.assertEqual(prepared.path, expected)
        self.assertEqual(prepared.omp_threads, 4)
        self.assertEqual(
            prepared.command[:8],
            [
                str(self.mpi_launcher),
                "--bind-to",
                "core",
                "--map-by",
                "slot:PE=4",
                "-n",
                "4",
                str(self.parallel_executable),
            ],
        )
        self.assertEqual(
            prepared.command[-2:],
            [str(expected / "inputs" / "mesh"), str(expected / "inputs" / "restart")],
        )
        self.assertTrue((expected / "outputs").is_dir())
        self.assertTrue((expected / "res").is_dir())
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
            self.bundle / "inputs" / "param.txt"
        ).read_text(encoding="utf-8")
        self.assertEqual(bundled_parameters, PARAMETERS)

    def test_prepares_seven_stage_fixed_mesh_workflow(self) -> None:
        prepared = prepare_run(
            self.settings,
            "legacy_case",
            "cold_fixed",
            "serial_omp1",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "cold",
        )

        self.assertIsInstance(prepared, PreparedStagedRun)
        self.assertEqual(len(prepared.stages), 7)
        first = prepared.stages[0].run
        second = prepared.stages[1].run
        self.assertEqual(len(first.command), 2)
        self.assertEqual(len(second.command), 3)
        self.assertFalse((first.path / "inputs/restart.h5").exists())
        self.assertFalse((second.path / "inputs/restart.h5").exists())
        self.assertTrue(
            all((stage.run.path / "res").is_dir() for stage in prepared.stages)
        )
        self.assertTrue((prepared.path / "inputs/reference.h5").is_symlink())

        first_parameters = (first.path / "param.txt").read_text(encoding="utf-8")
        self.assertIn(str(first.path / "inputs/transport_model.nml"), first_parameters)
        self.assertIn(f"{first.path / 'outputs'}/", first_parameters)
        self.assertIn("rest_adapt = .false.", first_parameters)

        plan = json.loads(
            (prepared.path / "run_plan.json").read_text(encoding="utf-8")
        )
        self.assertEqual(plan["workflow_kind"], "staged_fixed_mesh")
        self.assertEqual(
            [stage["stage_id"] for stage in plan["stages"]],
            [
                "time_init",
                "diffusion_reduction",
                "continuation_01",
                "continuation_02",
                "continuation_03",
                "continuation_04",
                "continuation_05",
            ],
        )

        adaptive = prepare_run(
            self.settings,
            "legacy_case",
            "cold_adaptive",
            "serial_omp1",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "adaptive",
        )
        adaptive_parameters = [
            (stage.run.path / "param.txt").read_text(encoding="utf-8")
            for stage in adaptive.stages
        ]
        self.assertEqual(
            (adaptive.stages[0].run.path / "inputs/mesh.msh").resolve(),
            (
                self.bundle
                / "inputs/mesh_adaptive_initial.msh"
            ).resolve(),
        )
        self.assertEqual(
            (first.path / "inputs/mesh.msh").resolve(),
            (self.bundle / "inputs/mesh.msh").resolve(),
        )
        self.assertIn("rest_adapt = .true.", adaptive_parameters[0])
        self.assertIn("rest_adapt = .true.", adaptive_parameters[1])
        for index in (2, 3, 4, 5, 6):
            self.assertIn("rest_adapt = .false.", adaptive_parameters[index])

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
            self.run_root / "legacy_case" / "warm" / "mpi4_omp4" / "existing"
        )
        run_dir.mkdir(parents=True)
        marker = run_dir / "keep.txt"
        marker.write_text("keep\n", encoding="utf-8")

        with self.assertRaisesRegex(BundleError, "run directory already exists"):
            prepare_run(
                self.settings,
                "legacy_case",
                "warm",
                "mpi4_omp4",
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "layouts.json",
                "existing",
            )

        self.assertEqual(marker.read_text(encoding="utf-8"), "keep\n")
