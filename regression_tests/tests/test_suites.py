from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

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

        self.fixture.install_solver(FAILING_SOLVER)
        resumed = self._run_suite("cold", "cold-test", "--run-only", "--resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.assertEqual(len(self._summary("cold", "cold-test")["results"]), 2)
        self.assertIn("skipping recorded cold_fixed / mpi4_omp4", resumed.stdout)

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
                    "schema_version": 1,
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
        self.assertEqual(report["status"], "failed")
        self.assertEqual(
            report["results"][-1]["failures"],
            ["solver run did not complete"],
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


def _suite_result(root: Path, workflow: str, status: str = "completed") -> dict:
    return {
        "workflow_id": workflow,
        "layout_id": "serial_omp1",
        "run_directory": str(root / workflow),
        "run_status": status,
    }


def _passing_report() -> dict:
    return {"status": "passed", "failures": []}


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

FAILING_SOLVER = """#!/usr/bin/env bash
exit 7
"""
