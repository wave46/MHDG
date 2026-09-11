"""Process outcomes and restart propagation, using tiny executable fixtures."""

import json
from dataclasses import replace

import pytest

from regression_tests.execute import execute_prepared, final_execution
from regression_tests.prepare import prepare_run
from support.files import file_identity
from tests.fixtures.harness import REGRESSION_ROOT, create_harness, run_command


@pytest.fixture
def harness(tmp_path):
    return create_harness(tmp_path, solver=SOLVER, include_provenance=True)


def run(harness, workflow="warm", layout="serial_omp1"):
    completed = run_command(
        "run", "legacy_case", workflow, "--settings", str(harness.settings),
        "--layout", layout, "--run-id", "execution-test",
    )
    directory = harness.run_directory(workflow, layout, "execution-test")
    metadata = json.loads((directory / "run_metadata.json").read_text())
    return completed, directory, metadata


def test_parallel_logs_environment_and_provenance(harness):
    completed, directory, metadata = run(harness, layout="mpi4_omp4")
    assert completed.returncode == 0, completed.stderr
    assert (directory / "stdout.log").read_text() == "solver stdout\n"
    assert (directory / "stderr.log").read_text() == "solver stderr\n"
    for filename, value in {
        "mpi_ranks": "4", "omp_threads": "4", "omp_places": "cores",
        "omp_proc_bind": "spread", "environment": "loaded",
    }.items():
        assert (directory / f"outputs/{filename}.txt").read_text() == value + "\n"
    assert metadata["status"] == "completed"
    assert metadata["environment"]["OMP_NUM_THREADS"] == "4"
    assert metadata["solver"]["build_manifest"]["path"] == str(harness.build_manifest)
    assert metadata["hdf5_outputs"] == ["outputs/result.h5"]
    for record, path in (
        (metadata["executable"], harness.parallel_executable),
        (metadata["environment"]["setup_script"], harness.environment_script),
        (metadata["solver"]["build_manifest"], harness.build_manifest),
        (metadata["runtime_files"][harness.runtime_file.name], harness.runtime_file),
    ):
        assert record == {"path": str(path), **file_identity(path)}
    for record in metadata["output_files"]:
        assert record == {"path": record["path"], **file_identity(directory / record["path"])}


@pytest.mark.parametrize("solver,workflow,status,exit_code", [
    ("no_output", "warm", "missing_hdf5_output", 0),
    ("failure", "warm", "solver_failed", 7),
    ("file_error", "cold_adaptive", "solver_reported_error", 0),
])
def test_solver_outcomes(harness, solver, workflow, status, exit_code):
    harness.install_solver({
        "no_output": NO_OUTPUT_SOLVER, "failure": FAILING_SOLVER,
        "file_error": FATAL_FILE_ERROR_SOLVER,
    }[solver], "serial")
    completed, directory, metadata = run(harness, workflow)
    assert completed.returncode == 1
    assert metadata["status"] == status
    assert metadata["exit_code"] == exit_code
    if solver == "file_error":
        first = json.loads((directory / "stages/01_time_init/run_metadata.json").read_text())
        assert len(first["fatal_log_messages"]) == 2
        assert all(stage["status"] == "not_run" for stage in metadata["stages"][1:])


def test_stages_pass_selected_output_to_next_restart(harness):
    # Multiple outputs force selection from the solver log, not filename ordering.
    harness.install_solver(STAGED_SOLVER + "printf 'decoy\\n' > outputs/another.h5\n", "serial")
    completed, directory, metadata = run(harness, "cold_fixed")
    assert completed.returncode == 0, completed.stderr
    stages = sorted((directory / "stages").iterdir())
    assert metadata["status"] == "completed"
    assert [stage["status"] for stage in metadata["stages"]] == ["completed"] * len(stages)
    final = stages[-1] / "outputs/result.h5"
    assert metadata["hdf5_outputs"] == [final.relative_to(directory).as_posix()]
    assert final.read_text() == ">".join(path.name for path in stages) + "\n"
    assert not (stages[0] / "inputs/restart.h5").exists()
    for previous, current in zip(stages, stages[1:]):
        assert (current / "inputs/restart.h5").resolve() == previous / "outputs/result.h5"
    assert (directory / "stdout.log").resolve() == stages[-1] / "stdout.log"
    observed_directory, observed = final_execution(directory, metadata)
    assert observed_directory == stages[-1]
    assert observed["executable"] == {"path": str(harness.serial_executable), **file_identity(harness.serial_executable)}


@pytest.mark.parametrize("ambiguous", [False, True])
def test_stages_stop_when_producer_fails_or_output_is_ambiguous(harness, ambiguous):
    solver = FAILING_STAGED_SOLVER
    if ambiguous:
        solver = solver.replace(
            "exit 7", "touch outputs/one.h5 outputs/two.h5\n  exit 0",
        )
    harness.install_solver(solver, "serial")
    completed, directory, metadata = run(harness, "cold_fixed")
    expected = "output_selection_failed" if ambiguous else "solver_failed"
    assert completed.returncode == 1
    assert metadata["status"] == expected
    assert [stage["status"] for stage in metadata["stages"]] == [
        "completed", "completed", expected, *(["not_run"] * 4),
    ]
    assert metadata["hdf5_outputs"] == []
    assert not (directory / "stages/04_continuation_02/run_metadata.json").exists()
    assert not (directory / "stages/04_continuation_02/inputs/restart.h5").exists()


@pytest.mark.parametrize("log_failure", [False, True])
def test_launch_failure_is_recorded(harness, log_failure):
    prepared = prepare_run(
        harness.settings, "legacy_case", "warm", "serial_omp1",
        REGRESSION_ROOT / "cases", REGRESSION_ROOT / "layouts.json", "launch-failure",
    )
    if log_failure:
        (prepared.path / "stdout.log").mkdir()
    else:
        prepared = replace(prepared, command=[str(harness.root / "missing-command")])
    result = execute_prepared(prepared, {})
    metadata = json.loads((prepared.path / "run_metadata.json").read_text())
    assert result.status == metadata["status"] == "launch_failed"
    assert result.exit_code is None
    assert metadata["launch_error"]


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
