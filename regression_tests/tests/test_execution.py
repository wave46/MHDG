from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REGRESSION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REGRESSION_ROOT / "tools"))

from tests.fixtures.harness import create_harness, run_command  # noqa: E402


class ExecutionWorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.fixture = create_harness(
            self.root,
            solver=SOLVER,
            include_provenance=True,
        )

    def test_parallel_run_records_logs_outputs_and_metadata(self) -> None:
        completed = self._run("mpi4_omp4", "parallel")
        self.assertEqual(completed.returncode, 0, completed.stderr)

        run_dir = self._run_dir("mpi4_omp4", "parallel")
        self.assertEqual((run_dir / "stdout.log").read_text(), "solver stdout\n")
        self.assertEqual((run_dir / "stderr.log").read_text(), "solver stderr\n")
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
        metadata = json.loads(
            (run_dir / "run_metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["status"], "completed")
        self.assertEqual(metadata["environment"]["OMP_NUM_THREADS"], "4")
        self.assertEqual(metadata["solver"]["revision"], "a" * 40)
        self.assertEqual(
            metadata["solver"]["build_manifest"]["path"],
            str(self.fixture.build_manifest),
        )
        self.assertEqual(metadata["hdf5_outputs"], ["outputs/result.h5"])

    def test_solver_outcomes_have_distinct_statuses(self) -> None:
        cases = (
            (NO_OUTPUT_SOLVER, "warm", "missing-output", "missing_hdf5_output", 0),
            (FAILING_SOLVER, "warm", "solver-failure", "solver_failed", 7),
            (
                FATAL_FILE_ERROR_SOLVER,
                "cold_adaptive",
                "reported-file-error",
                "solver_reported_error",
                0,
            ),
        )
        for solver, workflow, run_id, status, exit_code in cases:
            with self.subTest(status=status):
                self.fixture.install_solver(solver, "serial")
                completed = self._run("serial_omp1", run_id, workflow)
                self.assertEqual(completed.returncode, 1)
                metadata = json.loads(
                    (
                        self._run_dir("serial_omp1", run_id, workflow)
                        / "run_metadata.json"
                    ).read_text()
                )
                self.assertEqual(metadata["status"], status)
                self.assertEqual(metadata["exit_code"], exit_code)

        fatal_stage = self._run_dir(
            "serial_omp1",
            "reported-file-error",
            "cold_adaptive",
        ) / "stages/01_time_init/run_metadata.json"
        stage_metadata = json.loads(fatal_stage.read_text())
        self.assertEqual(len(stage_metadata["fatal_log_messages"]), 2)

    def test_staged_run_passes_each_output_to_the_next_stage(self) -> None:
        self.fixture.install_solver(STAGED_SOLVER, "serial")
        completed = self._run(
            "serial_omp1", "cold-success", workflow="cold_fixed"
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

        run_dir = self._run_dir(
            "serial_omp1", "cold-success", workflow="cold_fixed"
        )
        metadata = json.loads(
            (run_dir / "run_metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["status"], "completed")
        self.assertEqual(
            [stage["status"] for stage in metadata["stages"]],
            ["completed"] * 7,
        )
        self.assertEqual(
            metadata["hdf5_outputs"],
            ["stages/07_continuation_05/outputs/result.h5"],
        )

        stage_dirs = sorted((run_dir / "stages").iterdir())
        expected_history = ">".join(path.name for path in stage_dirs) + "\n"
        self.assertEqual(
            (stage_dirs[-1] / "outputs/result.h5").read_text(encoding="utf-8"),
            expected_history,
        )
        self.assertFalse((stage_dirs[0] / "inputs/restart.h5").exists())
        for previous, current in zip(stage_dirs, stage_dirs[1:]):
            self.assertEqual(
                (current / "inputs/restart.h5").resolve(),
                (previous / "outputs/result.h5").resolve(),
            )
        self.assertEqual(
            (run_dir / "stdout.log").resolve(),
            (stage_dirs[-1] / "stdout.log").resolve(),
        )

    def test_staged_run_stops_after_failed_stage(self) -> None:
        self.fixture.install_solver(FAILING_STAGED_SOLVER, "serial")
        completed = self._run(
            "serial_omp1", "cold-failure", workflow="cold_fixed"
        )
        self.assertEqual(completed.returncode, 1)

        run_dir = self._run_dir(
            "serial_omp1", "cold-failure", workflow="cold_fixed"
        )
        metadata = json.loads(
            (run_dir / "run_metadata.json").read_text(encoding="utf-8")
        )
        self.assertEqual(metadata["status"], "solver_failed")
        self.assertEqual(
            [stage["status"] for stage in metadata["stages"]],
            ["completed", "completed", "solver_failed"] + ["not_run"] * 4,
        )
        self.assertFalse(
            (run_dir / "stages/04_continuation_02/run_metadata.json").exists()
        )

    def _run(
        self, layout: str, run_id: str, workflow: str = "warm"
    ) -> subprocess.CompletedProcess[str]:
        return self._run_workflows(layout, run_id, [workflow])

    def _run_workflows(
        self, layout: str, run_id: str, workflows: list[str]
    ) -> subprocess.CompletedProcess[str]:
        return run_command(
            "run",
            "legacy_case",
            *workflows,
            "--settings",
            str(self.fixture.settings),
            "--layout",
            layout,
            "--run-id",
            run_id,
        )

    def _run_dir(self, layout: str, run_id: str, workflow: str = "warm") -> Path:
        return self.fixture.run_directory(workflow, layout, run_id)


SOLVER = """#!/usr/bin/env bash
set -euo pipefail
test -d res
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

STAGED_SOLVER = """#!/usr/bin/env bash
set -euo pipefail
test -d res
stage=${PWD##*/}
if (($# == 1)); then
  history=$stage
else
  history=$(<"$2.h5")
  history=${history%$'\\n'}">"$stage
fi
printf '%s\n' "$history" > outputs/result.h5
printf 'Error: 1.0E-5\n'
printf 'Output written to file %s\n' "$PWD/outputs/result.h5"
"""

FATAL_FILE_ERROR_SOLVER = """#!/usr/bin/env bash
set -euo pipefail
test -d res
printf 'synthetic hdf5\n' > outputs/result.h5
printf 'Error opening destination file:./res/temp.msh\n'
printf "Error   : Unable to open file './res/temp.msh'\n" >&2
"""

FAILING_STAGED_SOLVER = """#!/usr/bin/env bash
set -euo pipefail
stage=${PWD##*/}
if [[ "$stage" == '03_continuation_01' ]]; then
  printf 'failed stage\n' >&2
  exit 7
fi
printf '%s\n' "$stage" > outputs/result.h5
printf 'Output written to file %s\n' "$PWD/outputs/result.h5"
"""
