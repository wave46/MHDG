from __future__ import annotations

import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from bundle.cases import load_case_definition  # noqa: E402
from suite.configuration import load_suite_definition  # noqa: E402
from suite.pairs import compare_generated_meshes  # noqa: E402
from suite.runner import _comparison_mode, run_suite  # noqa: E402
from suite.verification import verify_suite  # noqa: E402
from tests.fixtures.harness import create_harness, run_command  # noqa: E402
from tests.fixtures.solutions import write_solution  # noqa: E402


class SuiteWorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.fixture = create_harness(
            self.root,
            solver=SOLVER,
            reference_writer=write_solution,
        )

    def test_cold_suite_records_deferred_workflows_and_resumes(self) -> None:
        self.fixture.install_solver(WORKFLOW_SOLVER)
        completed = self._run_suite("cold", "cold-test", "--run-only")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = self._summary("cold", "cold-test")
        self.assertEqual(summary["status"], "passed")
        self.assertEqual(summary["comparison_mode"], "deferred")
        self.assertEqual(
            {result["workflow_id"] for result in summary["results"]},
            {"cold_fixed", "cold_adaptive"},
        )
        self.assertTrue(
            all(result["comparison_status"] == "not_run" for result in summary["results"])
        )

        resumed = self._run_suite("cold", "cold-test", "--run-only", "--resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertEqual(len(self._summary("cold", "cold-test")["results"]), 2)
        self.assertIn("skipping recorded cold_fixed / mpi4_omp4", resumed.stdout)

    def test_resume_preserves_and_retries_an_unrecorded_cell(self) -> None:
        run_id = "interrupted-test"
        partial = self.fixture.run_directory("warm", "mpi4_omp4", run_id)
        marker = partial / "partial-output.txt"

        def interrupt_cell(*_args) -> None:
            partial.mkdir(parents=True)
            marker.write_text("keep\n", encoding="utf-8")
            raise KeyboardInterrupt

        with (
            redirect_stdout(StringIO()),
            patch("suite.runner.run_cell", side_effect=interrupt_cell),
            self.assertRaises(KeyboardInterrupt),
        ):
            run_suite(
                self.fixture.settings,
                "warm",
                REGRESSION_ROOT / "cases",
                REGRESSION_ROOT / "layouts.json",
                REGRESSION_ROOT / "suites.json",
                REGRESSION_ROOT / "tolerances.json",
                run_id,
                compare=False,
            )

        resumed = self._run_suite("warm", run_id, "--run-only", "--resume")

        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertEqual(marker.read_text(encoding="utf-8"), "keep\n")
        result = self._summary("warm", run_id)["results"][0]
        self.assertEqual(
            Path(result["run_directory"]).name,
            "interrupted-test-resume-1",
        )
        self.assertIn("preserving unrecorded run", resumed.stdout)

    def test_resume_rejects_changed_execution_inputs_and_new_build(self) -> None:
        run_id = "stable-inputs"
        completed = self._run_suite("warm", run_id, "--run-only")
        self.assertEqual(completed.returncode, 0, completed.stderr)

        self.fixture.install_solver(FAILING_SOLVER, "parallel")
        changed = self._run_suite("warm", run_id, "--run-only", "--resume")
        self.assertEqual(changed.returncode, 1)
        self.assertIn("parallel_executable", changed.stderr)

        rebuilt = self._run_suite(
            "warm",
            run_id,
            "--run-only",
            "--resume",
            "--build",
        )
        self.assertEqual(rebuilt.returncode, 1)
        self.assertIn("--build cannot be used with --resume", rebuilt.stderr)

    def test_suite_check_requires_golden_data_and_runs_warm_defaults(self) -> None:
        rejected = run_command(
            "suite",
            "check",
            "--settings",
            str(self.fixture.settings),
        )
        self.assertEqual(rejected.returncode, 1)
        self.assertIn("requires bundle_class=golden", rejected.stderr)

        self.fixture.set_bundle_class("golden")
        completed = run_command(
            "suite",
            "check",
            environment={
                "MHDG_REGRESSION_GOLDEN_SETTINGS": str(self.fixture.settings)
            },
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("suite: warm", completed.stdout)
        self.assertIn("mpi4_omp4", completed.stdout)

    def test_race_suite_compares_layout_pairs_and_can_recompare(self) -> None:
        write_solution(self.fixture.serial_executable.parent / "race_result.h5")
        self.fixture.install_solver(RACE_SOLVER)
        completed = self._run_suite("race", "race-test")

        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary_path = (
            self.fixture.run_root
            / "suites/race/race-test/suite_summary.json"
        )
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        self.assertEqual(summary["comparison_mode"], "layout_pairs")
        self.assertEqual(
            summary["layout_ids"],
            ["serial_omp1", "serial_omp16"],
        )
        self.assertEqual(len(summary["results"]), 4)
        self.assertTrue(
            all(
                result["comparison_status"] == "not_run"
                for result in summary["results"]
            )
        )
        self.assertEqual(len(summary["comparisons"]), 2)
        self.assertTrue(
            all(item["status"] == "passed" for item in summary["comparisons"])
        )
        self.assertEqual(
            Path(summary["comparisons"][0]["comparison_report"]).name,
            "comparison_from_serial_omp1.json",
        )

        report = json.loads(
            Path(summary["comparisons"][0]["comparison_report"]).read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(report["convergence"]["final_newton_error"], 1.0e5)
        self.assertIsNone(report["convergence"]["maximum"])
        self.assertTrue(report["convergence"]["passed"])

        rechecked = run_command("suite", "compare", str(summary_path))
        self.assertEqual(rechecked.returncode, 0, rechecked.stderr)
        self.assertIn("serial_omp1", rechecked.stdout)
        self.assertIn("serial_omp16", rechecked.stdout)

    def test_cold_matrix_combines_golden_and_all_pair_checks(self) -> None:
        suite = load_suite_definition(
            "cold_matrix",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )

        self.assertEqual(len(suite["layout_comparisons"]), 6)
        self.assertEqual(len(suite["layouts"]), 4)
        self.assertTrue(suite["reference_comparisons"])
        self.assertEqual(suite["layout_comparison_policy"], "fixed_hdf5")
        self.assertEqual(_comparison_mode(True, True, False), "deferred")
        self.assertEqual(
            _comparison_mode(True, False, False),
            "deferred_layout_pairs",
        )
        smoke = load_suite_definition(
            "initialization_smoke",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        self.assertFalse(smoke["reference_comparisons"])
        self.assertEqual(_comparison_mode(False, False, True), "execution_only")

    def test_pr04_neutral_suites_are_minimal(self) -> None:
        pressure = load_suite_definition(
            "neutral_pressure_warm",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        self.assertEqual(pressure["workflow_ids"], ["warm_neutral_pressure"])
        self.assertEqual(pressure["layouts"], ["mpi4_omp4"])
        self.assertFalse(pressure["reference_comparisons"])

        neutralgamma = load_suite_definition(
            "neutralgamma_race",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        self.assertEqual(neutralgamma["workflow_ids"], ["cold_step_neutralgamma"])
        self.assertEqual(
            neutralgamma["layouts"], ["serial_omp1", "serial_omp16"]
        )
        self.assertEqual(
            neutralgamma["layout_comparisons"],
            [{"baseline": "serial_omp1", "candidate": "serial_omp16"}],
        )
        self.assertFalse(neutralgamma["reference_comparisons"])

    def test_pr05_wall_source_suites_are_focused(self) -> None:
        warm = load_suite_definition(
            "neutral_sources_in_elements_warm",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        race = load_suite_definition(
            "neutral_sources_in_elements_race_matrix",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        cold_adaptive = load_suite_definition(
            "neutral_sources_in_elements_cold_adaptive",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        restart_producer = load_suite_definition(
            "neutral_sources_in_elements_restart_producer",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )

        self.assertTrue(warm["reference_comparisons"])
        self.assertEqual(len(race["layouts"]), 4)
        self.assertEqual(len(race["layout_comparisons"]), 6)
        self.assertEqual(cold_adaptive["layouts"], ["mpi4_omp4"])
        self.assertFalse(cold_adaptive["reference_comparisons"])
        self.assertEqual(
            restart_producer["workflow_ids"],
            ["bootstrap_neutral_sources_in_elements"],
        )
        self.assertFalse(restart_producer["reference_comparisons"])

    def test_pr06_neutral_feature_suites_are_focused(self) -> None:
        warm = load_suite_definition(
            "neutral_features_warm",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        producers = load_suite_definition(
            "neutral_feature_references",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        race = load_suite_definition(
            "neutral_feature_race_matrix",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )
        impurity_restarts = load_suite_definition(
            "impurity_restart_producers",
            REGRESSION_ROOT / "suites.json",
            REGRESSION_ROOT / "layouts.json",
            REGRESSION_ROOT / "cases",
        )

        self.assertEqual(
            warm["workflow_ids"],
            [
                "warm_neutral_sources_in_elements",
                "warm_neutral_pressure",
                "warm_neutral_perpendicular",
                "warm_neutral_limiter_fixed",
                "warm_neutral_limiter_ti",
            ],
        )
        self.assertEqual(warm["layouts"], ["mpi4_omp4"])
        self.assertTrue(warm["reference_comparisons"])
        self.assertEqual(
            producers["workflow_ids"],
            [
                "warm_neutral_pressure",
                "warm_neutral_perpendicular",
                "warm_neutral_limiter_fixed",
                "warm_neutral_limiter_ti",
            ],
        )
        self.assertEqual(producers["layouts"], ["mpi4_omp4"])
        self.assertTrue(producers["reference_comparisons"])
        self.assertEqual(
            race["workflow_ids"],
            [
                "race_neutral_sources_in_elements",
                "race_neutral_pressure",
                "race_neutral_perpendicular",
                "race_neutral_limiter_fixed",
                "race_neutral_limiter_ti",
            ],
        )
        self.assertEqual(len(race["layouts"]), 4)
        self.assertEqual(len(race["layout_comparisons"]), 6)
        self.assertEqual(race["tolerance_profile"], "neutral_feature_race_step")
        self.assertFalse(race["reference_comparisons"])
        self.assertEqual(
            impurity_restarts["workflow_ids"],
            [
                "bootstrap_impurity_off",
                "bootstrap_impurity_n",
                "bootstrap_impurity_nw",
            ],
        )
        self.assertFalse(impurity_restarts["reference_comparisons"])

        case = load_case_definition("legacy_case", REGRESSION_ROOT / "cases")
        self.assertEqual(
            case["workflows"]["warm_neutral_sources_in_elements"][
                "tolerance_profile"
            ],
            "fixed_same_layout",
        )
        for workflow in producers["workflow_ids"]:
            self.assertEqual(
                case["workflows"][workflow]["tolerance_profile"],
                "neutral_feature_same_layout",
            )

    def test_generated_adaptive_mesh_comparison_is_byte_exact(self) -> None:
        reference = self.root / "reference/stages/01_single_step/res/temp.msh"
        candidate = self.root / "candidate/stages/01_single_step/res/temp.msh"
        reference.parent.mkdir(parents=True)
        candidate.parent.mkdir(parents=True)
        reference.write_text("same mesh\n", encoding="utf-8")
        candidate.write_text("same mesh\n", encoding="utf-8")

        matching = compare_generated_meshes(
            self.root / "reference", self.root / "candidate"
        )
        self.assertIsNotNone(matching)
        self.assertTrue(matching["passed"])

        candidate.write_text("different mesh\n", encoding="utf-8")
        differing = compare_generated_meshes(
            self.root / "reference", self.root / "candidate"
        )
        self.assertIsNotNone(differing)
        self.assertFalse(differing["passed"])

    @patch("suite.verification.compare_completed_run")
    def test_offline_verification_dispatches_only_completed_runs(self, compare) -> None:
        compare.side_effect = (
            ("fixed_hdf5", self.root / "fixed.json", _passing_report()),
            ("mesh_independent", self.root / "adaptive.json", _passing_report()),
        )
        source = self.root / "suite_summary.json"
        source.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "suite_id": "cold_matrix",
                    "run_id": "offline-test",
                    "case_id": "legacy_case",
                    "results": [
                        _suite_result(self.root, "cold_fixed"),
                        _suite_result(self.root, "cold_adaptive"),
                        _suite_result(self.root, "cold_fixed", "solver_failed"),
                    ],
                }
            ),
            encoding="utf-8",
        )

        path, report = verify_suite(
            source,
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "tolerances.json",
        )

        self.assertEqual(path.name, "verification_summary.json")
        self.assertEqual(compare.call_count, 2)
        self.assertEqual(
            [result["comparison_policy"] for result in report["results"]],
            ["fixed_hdf5", "mesh_independent", None],
        )
        self.assertEqual(
            [result["convergence_status"] for result in report["results"]],
            ["passed", "passed", None],
        )
        self.assertEqual(report["status"], "failed")
        self.assertEqual(
            report["results"][-1]["failures"],
            ["solver run did not complete"],
        )

    @patch("suite.verification.compare_layout_pairs")
    @patch("suite.verification.compare_completed_run")
    def test_offline_verification_combines_references_and_pairs(
        self,
        compare,
        compare_pairs,
    ) -> None:
        compare.return_value = (
            "fixed_hdf5",
            self.root / "reference.json",
            _passing_report(),
        )
        compare_pairs.return_value = [
            {
                "workflow_id": "cold_fixed",
                "baseline_layout_id": "serial_omp1",
                "candidate_layout_id": "serial_omp16",
                "status": "passed",
                "failures": [],
            }
        ]
        source = self.root / "suite_summary.json"
        source.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "suite_id": "cold_matrix",
                    "run_id": "combined-test",
                    "case_id": "legacy_case",
                    "workflow_ids": ["cold_fixed"],
                    "tolerance_profile": "cold_cross_layout",
                    "reference_comparisons": True,
                    "layout_comparisons": [
                        {
                            "baseline": "serial_omp1",
                            "candidate": "serial_omp16",
                        }
                    ],
                    "results": [
                        _suite_result(self.root, "cold_fixed"),
                        _suite_result(
                            self.root,
                            "cold_fixed",
                            layout="serial_omp16",
                        ),
                    ],
                }
            ),
            encoding="utf-8",
        )

        _, report = verify_suite(
            source,
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "tolerances.json",
        )

        self.assertEqual(compare.call_count, 2)
        compare_pairs.assert_called_once()
        self.assertEqual(len(report["results"]), 2)
        self.assertEqual(len(report["comparisons"]), 1)
        self.assertEqual(report["status"], "passed")

        compare_pairs.reset_mock()
        _, references_only = verify_suite(
            source,
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "tolerances.json",
            include_layout_pairs=False,
        )
        compare_pairs.assert_not_called()
        self.assertNotIn("comparisons", references_only)

    @patch("suite.verification.compare_completed_run")
    def test_offline_verification_checks_every_staged_producer_log(
        self,
        compare,
    ) -> None:
        workflow = load_case_definition(
            "legacy_case",
            REGRESSION_ROOT / "cases",
        )["workflows"]["cold_fixed"]
        run_directory = self.root / "cold_fixed/serial_omp1"
        records = []
        for index, stage in enumerate(workflow["stages"], start=1):
            stage_directory = run_directory / f"stage-{index}"
            stage_directory.mkdir(parents=True)
            (stage_directory / "stdout.log").write_text(
                "Error: 1.0E-5\n",
                encoding="utf-8",
            )
            records.append(
                {
                    "stage_id": stage["stage_id"],
                    "status": "completed",
                    "run_directory": str(stage_directory),
                }
            )
        (run_directory / "run_metadata.json").write_text(
            json.dumps({"stages": records}),
            encoding="utf-8",
        )
        compare.return_value = (
            "reference_matrix",
            self.root / "matrix.json",
            {"status": "failed", "failures": ["old fields differ"]},
        )
        source = self.root / "staged_suite_summary.json"
        source.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "suite_id": "cold",
                    "run_id": "staged-convergence",
                    "case_id": "legacy_case",
                    "results": [
                        _suite_result(self.root, "cold_fixed"),
                    ],
                }
            ),
            encoding="utf-8",
        )

        _, passing = verify_suite(
            source,
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "tolerances.json",
        )
        self.assertEqual(
            passing["results"][0]["convergence_status"],
            "passed",
        )

        last_log = Path(records[-1]["run_directory"]) / "stdout.log"
        last_log.write_text("Error: 1.0E-3\n", encoding="utf-8")
        _, failing = verify_suite(
            source,
            REGRESSION_ROOT / "cases",
            REGRESSION_ROOT / "tolerances.json",
        )
        self.assertEqual(
            failing["results"][0]["convergence_status"],
            "failed",
        )

    def _run_suite(self, suite: str, run_id: str, *arguments: str):
        return run_command(
            "suite",
            "run",
            suite,
            "--settings",
            str(self.fixture.settings),
            "--run-id",
            run_id,
            *arguments,
        )

    def _summary(self, suite: str, run_id: str) -> dict:
        path = self.fixture.run_root / "suites" / suite / run_id / "suite_summary.json"
        return json.loads(path.read_text(encoding="utf-8"))


def _suite_result(
    root: Path,
    workflow: str,
    status: str = "completed",
    layout: str = "serial_omp1",
) -> dict:
    return {
        "workflow_id": workflow,
        "layout_id": layout,
        "run_directory": str(root / workflow / layout),
        "run_status": status,
    }


def _passing_report() -> dict:
    return {
        "status": "passed",
        "failures": [],
        "convergence": {"passed": True},
    }


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
elif (($# == 1)); then
  printf '%s\n' "${PWD##*/}" > outputs/result.h5
else
  cp "$2.h5" outputs/result.h5
fi
printf 'Error: 1.0E-5\n'
printf 'Output written to file outputs/result.h5\n'
"""

RACE_SOLVER = """#!/usr/bin/env bash
set -euo pipefail
cp "$(dirname "$0")/race_result.h5" outputs/result.h5
printf 'Error: 1.0E+5\n'
printf 'Output written to file outputs/result.h5\n'
"""

FAILING_SOLVER = """#!/usr/bin/env bash
exit 7
"""
