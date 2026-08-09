from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from bundle.case_definition import normalize_case_definition  # noqa: E402
from preparation.models import PreparedStagedRun  # noqa: E402
from preparation.parameters import render_parameter_file  # noqa: E402
from preparation.workspace import populate_warm_run  # noqa: E402
from prepare_run import prepare_run  # noqa: E402
from support.errors import BundleError  # noqa: E402
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
            (expected / "inputs/restart.h5").resolve(),
            (self.bundle / "inputs/restart.h5").resolve(),
        )
        self.assertEqual(
            (expected / "inputs/reference.h5").resolve(),
            (self.bundle / "inputs/reference_mpi4_omp4.h5").resolve(),
        )
        self.assertEqual(
            (expected / "positionFeketeNodesTri2D.h5").resolve(),
            self.runtime_file.resolve(),
        )

        parameters = (expected / "param.txt").read_text(encoding="utf-8")
        self.assertNotIn("/old/", parameters)
        rendered_paths = (
            expected / "inputs" / "transport_model.nml",
            expected / "inputs" / "impurity_model.nml",
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

    def test_selects_warm_impurity_configuration(self) -> None:
        prepared = prepare_run(
            self.settings,
            "legacy_case",
            "warm_impurity_n",
            "mpi4_omp4",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "nitrogen",
        )

        parameters = (prepared.path / "param.txt").read_text(encoding="utf-8")
        self.assertIn("impurity_radiation = .true.", parameters)
        impurity_configuration = prepared.path / "inputs/impurity_model.nml"
        self.assertEqual(
            impurity_configuration.resolve(),
            (self.bundle / "inputs/impurity_model_n.nml").resolve(),
        )
        self.assertIn("impurity_names = 'N'", impurity_configuration.read_text())

        plan = json.loads(
            (prepared.path / "run_plan.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            plan["parameter_overrides"],
            {
                "compute_from_flux": True,
                "impurity_radiation": True,
            },
        )

    def test_prepares_pr04_neutral_variants(self) -> None:
        pressure = prepare_run(
            self.settings,
            "legacy_case",
            "warm_neutral_pressure",
            "mpi4_omp4",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "pressure-on",
        )
        pressure_parameters = (pressure.path / "param.txt").read_text(
            encoding="utf-8"
        )
        self.assertIn("neutralp_lambda = 1", pressure_parameters)
        self.assertIn("compute_from_flux = .true.", pressure_parameters)

        neutralgamma = prepare_run(
            self.settings,
            "legacy_case",
            "cold_step_neutralgamma",
            "serial_omp1",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "neutralgamma",
        )
        neutralgamma_parameters = (
            neutralgamma.stages[0].run.path / "param.txt"
        ).read_text(encoding="utf-8")
        for assignment in (
            "impurity_radiation = .false.",
            "nrp = 2",
            "nts = 1",
            "tau(6) = 1.0",
            "tau(7) = 1.0",
        ):
            self.assertIn(assignment, neutralgamma_parameters)

    def test_warm_workflow_selects_restart_and_reference_independently(self) -> None:
        staging = self.root / "selected_staging"
        final = self.root / "selected_final"
        staging.mkdir()
        artifact_roles = {
            "mesh",
            "geometry",
            "equilibrium_magnetic_field",
            "equilibrium_current_density",
            "transport_configuration",
            "warm_parameters",
            "impurity_configuration",
            "selected_restart",
            "selected_reference",
        }
        artifacts = {}
        for role in artifact_roles:
            artifact = self.root / role
            artifact.write_text(
                PARAMETERS if role == "warm_parameters" else role,
                encoding="utf-8",
            )
            artifacts[role] = artifact

        populate_warm_run(
            staging,
            final,
            artifacts,
            {},
            {
                "restart_role": "selected_restart",
                "reference_role": "selected_reference",
                "impurity_configuration_role": "impurity_configuration",
            },
            {},
        )

        self.assertEqual(
            (staging / "inputs/restart.h5").resolve(),
            artifacts["selected_restart"].resolve(),
        )
        self.assertEqual(
            (staging / "inputs/reference.h5").resolve(),
            artifacts["selected_reference"].resolve(),
        )

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

    def test_derived_workflow_overlays_parameter_values(self) -> None:
        case = normalize_case_definition(
            "example",
            {
                "schema_version": 2,
                "description": "Synthetic inheritance contract",
                "reference": {"branch": "develop", "revision": "a" * 40},
                "workflows": {
                    "base": {
                        "type": "warm_same_state",
                        "description": "Base workflow",
                        "parameter_overrides": {
                            "rest_adapt": False,
                            "nrp": 1,
                        },
                    },
                    "variant": {
                        "extends": "base",
                        "description": "Derived parameter variant",
                        "parameter_overrides": {
                            "nrp": 2,
                            "neutral_pressure_option": 1,
                        },
                    },
                },
            },
        )

        self.assertEqual(
            case["workflows"]["variant"]["parameter_overrides"],
            {
                "rest_adapt": False,
                "nrp": 2,
                "neutral_pressure_option": 1,
            },
        )
        self.assertEqual(
            case["workflows"]["base"]["parameter_overrides"]["nrp"],
            1,
        )

    def test_prepares_coarse_one_step_race_workflows(self) -> None:
        fixed = prepare_run(
            self.settings,
            "legacy_case",
            "cold_step_fixed",
            "serial_omp1",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "race-fixed",
        )
        adaptive = prepare_run(
            self.settings,
            "legacy_case",
            "cold_step_adaptive",
            "serial_omp1",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "race-adaptive",
        )

        self.assertEqual(len(fixed.stages), 1)
        self.assertEqual(len(adaptive.stages), 1)
        coarse_mesh = self.bundle / "inputs/mesh_adaptive_initial.msh"
        for prepared in (fixed, adaptive):
            stage = prepared.stages[0].run
            self.assertEqual(
                (stage.path / "inputs/mesh.msh").resolve(),
                coarse_mesh.resolve(),
            )
            self.assertEqual(len(stage.command), 2)
            self.assertEqual(
                (stage.path / "inputs/impurity_model.nml").resolve(),
                (self.bundle / "inputs/impurity_model_w.nml").resolve(),
            )
            parameters = (stage.path / "param.txt").read_text(encoding="utf-8")
            for assignment in (
                "steady = .false.",
                "saveNR = .false.",
                "time_adapt = .false.",
                "NR_adapt = .false.",
                "div_adapt = .false.",
                "osc_adapt = .false.",
                "nts = 1",
            ):
                self.assertIn(assignment, parameters)

        fixed_parameters = (
            fixed.stages[0].run.path / "param.txt"
        ).read_text(encoding="utf-8")
        adaptive_parameters = (
            adaptive.stages[0].run.path / "param.txt"
        ).read_text(encoding="utf-8")
        self.assertIn("adaptivity = .false.", fixed_parameters)
        self.assertIn("rest_adapt = .false.", fixed_parameters)
        self.assertIn("nrp = 2", fixed_parameters)
        self.assertIn("adaptivity = .true.", adaptive_parameters)
        self.assertIn("rest_adapt = .true.", adaptive_parameters)
        self.assertIn("nrp = 2", adaptive_parameters)

        fixed_plan = json.loads(
            (fixed.path / "run_plan.json").read_text(encoding="utf-8")
        )
        adaptive_plan = json.loads(
            (adaptive.path / "run_plan.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            fixed_plan["stages"][0]["parameter_overrides"]["nrp"], 2
        )
        self.assertEqual(
            adaptive_plan["stages"][0]["parameter_overrides"]["nrp"], 2
        )

    def test_prepares_disabled_impurity_scratch_workflow(self) -> None:
        prepared = prepare_run(
            self.settings,
            "legacy_case",
            "cold_step_impurity_off",
            "mpi4_omp4",
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "layouts.json",
            "impurity-off-scratch",
        )

        stage = prepared.stages[0].run
        self.assertFalse((stage.path / "inputs/restart.h5").exists())
        parameters = (stage.path / "param.txt").read_text(encoding="utf-8")
        self.assertIn("impurity_radiation = .false.", parameters)
        self.assertIn("nrp = 2", parameters)
        self.assertIn("nts = 1", parameters)

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
