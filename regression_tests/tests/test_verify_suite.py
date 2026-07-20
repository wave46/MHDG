from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from verify_suite import verify_suite  # noqa: E402


class SuiteVerificationTests(unittest.TestCase):
    def test_shell_command_writes_failed_verification_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = self._result(root, "cold_fixed")
            result["run_status"] = "solver_failed"
            source = self._summary(root, [result])

            completed = subprocess.run(
                [
                    str(REGRESSION_ROOT / "regression.sh"),
                    "suite-verify",
                    str(source),
                ],
                check=False,
                capture_output=True,
                env={**os.environ, "PYTHON": sys.executable},
                text=True,
            )

            self.assertEqual(completed.returncode, 1, completed.stderr)
            self.assertIn("verification failed:", completed.stdout)
            self.assertTrue((root / "verification_summary.json").is_file())

    @patch("verify_suite.compare_completed_run")
    def test_verification_dispatches_without_running_solver(self, compare) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = self._summary(
                root,
                [
                    self._result(root, "cold_fixed"),
                    self._result(root, "cold_adaptive"),
                ],
            )
            compare.side_effect = (
                (
                    "fixed_hdf5",
                    root / "fixed.json",
                    {"status": "passed", "failures": []},
                ),
                (
                    "mesh_independent",
                    root / "adaptive.json",
                    {"status": "passed", "failures": []},
                ),
            )

            path, report = verify_suite(
                source, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json"
            )

            self.assertEqual(report["status"], "passed")
            self.assertEqual(
                [result["comparison_policy"] for result in report["results"]],
                ["fixed_hdf5", "mesh_independent"],
            )
            self.assertEqual(path.name, "verification_summary.json")
            self.assertTrue(path.is_file())
            self.assertEqual(compare.call_count, 2)

    @patch("verify_suite.compare_completed_run")
    def test_incomplete_run_fails_without_comparison(self, compare) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = self._result(root, "cold_fixed")
            result["run_status"] = "solver_failed"
            source = self._summary(root, [result])

            _, report = verify_suite(
                source, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json"
            )

            self.assertEqual(report["status"], "failed")
            self.assertEqual(
                report["results"][0]["failures"], ["solver run did not complete"]
            )
            compare.assert_not_called()

    @staticmethod
    def _summary(root: Path, results: list[dict[str, object]]) -> Path:
        path = root / "suite_summary.json"
        path.write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "suite_id": "cold_matrix",
                    "run_id": "overnight-test",
                    "case_id": "legacy_case",
                    "results": results,
                }
            ),
            encoding="utf-8",
        )
        return path

    @staticmethod
    def _result(root: Path, workflow: str) -> dict[str, object]:
        return {
            "workflow_id": workflow,
            "layout_id": "serial_omp1",
            "run_directory": str(root / workflow),
            "run_status": "completed",
        }


if __name__ == "__main__":
    unittest.main()
