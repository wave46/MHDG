from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from check_neutral_wall_sources import (  # noqa: E402
    check_run,
    check_suite,
)


class NeutralWallSourceCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)

    def test_checks_every_wall_volume_pair(self) -> None:
        run = self._write_run(
            "conservative",
            [
                ("puff", 1.0e20, 1.0e20),
                ("pump", 2.0e20, 2.0e20),
                ("puff", 1.0e20, 1.0e20 * (1.0 + 5.0e-13)),
                ("pump", 3.0e20, 3.0e20 * (1.0 + 2.0e-12)),
            ],
        )

        report = check_run(run)

        self.assertEqual(report["status"], "failed")
        self.assertEqual(report["record_counts"], {"puff": 2, "pump": 2})
        self.assertEqual(len(report["failures"]), 1)
        self.assertIn("pump wall/volume", report["failures"][0])

    def test_suite_compares_final_layout_totals_at_race_tolerance(self) -> None:
        baseline = self._write_run(
            "baseline",
            [("puff", 1.0e20, 1.0e20), ("pump", 2.0e20, 2.0e20)],
        )
        close = self._write_run(
            "close",
            [
                ("puff", 1.0e20 * (1.0 + 4.0e-8), 1.0e20 * (1.0 + 4.0e-8)),
                ("pump", 2.0e20 * (1.0 + 4.0e-8), 2.0e20 * (1.0 + 4.0e-8)),
            ],
        )
        far = self._write_run(
            "far",
            [
                ("puff", 1.0e20 * (1.0 + 6.0e-8), 1.0e20 * (1.0 + 6.0e-8)),
                ("pump", 2.0e20, 2.0e20),
            ],
        )

        passing_summary = self._write_suite("passing", baseline, close)
        self.assertEqual(check_suite(passing_summary)["status"], "passed")

        failing_summary = self._write_suite("failing", baseline, far)
        failed = check_suite(failing_summary)
        self.assertEqual(failed["status"], "failed")
        self.assertTrue(any("final puff" in item for item in failed["failures"]))

    def test_rejects_nonpositive_final_totals(self) -> None:
        run = self._write_run(
            "nonpositive",
            [("puff", 0.0, 0.0), ("pump", -2.0e20, -2.0e20)],
        )

        report = check_run(run)

        self.assertEqual(report["status"], "failed")
        self.assertTrue(
            any("final puff totals must be positive" in item for item in report["failures"])
        )
        self.assertTrue(
            any("final pump totals must be positive" in item for item in report["failures"])
        )

    def _write_run(
        self,
        name: str,
        records: list[tuple[str, float, float]],
    ) -> Path:
        run = self.root / name
        run.mkdir()
        (run / "run_metadata.json").write_text(
            json.dumps({"schema_version": 2, "status": "completed"}),
            encoding="utf-8",
        )
        lines = [
            "NEUTRAL_WALL_SOURCE_CONSERVATION "
            f"{source} wall={wall:.16E} volume={volume:.16E}\n"
            for source, wall, volume in records
        ]
        (run / "stdout.log").write_text("".join(lines), encoding="utf-8")
        return run

    def _write_suite(self, name: str, baseline: Path, candidate: Path) -> Path:
        path = self.root / f"{name}.json"
        path.write_text(
            json.dumps(
                {
                    "schema_version": 2,
                    "workflow_ids": ["cold_step"],
                    "layout_comparisons": [
                        {"baseline": "serial_omp1", "candidate": "mpi4_omp1"}
                    ],
                    "results": [
                        {
                            "workflow_id": "cold_step",
                            "layout_id": "serial_omp1",
                            "run_status": "completed",
                            "run_directory": str(baseline),
                        },
                        {
                            "workflow_id": "cold_step",
                            "layout_id": "mpi4_omp1",
                            "run_status": "completed",
                            "run_directory": str(candidate),
                        },
                    ],
                }
            ),
            encoding="utf-8",
        )
        return path


if __name__ == "__main__":
    unittest.main()
