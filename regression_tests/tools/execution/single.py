"""Execute and document one prepared solver run."""

from __future__ import annotations

from typing import Any

from execution.artifacts import (
    fatal_log_messages,
    file_record,
    optional_file_record,
    output_records,
)
from execution.models import (
    RunObservations,
    RunResult,
    RuntimeEnvironment,
    SolverProcessResult,
)
from execution.process import launch_solver, runtime_environment
from preparation.models import PreparedRun
from support.documents import write_json_atomic
from support.time import utc_now


def execute_run(prepared: PreparedRun, settings: dict[str, str]) -> RunResult:
    """Execute a prepared command and record logs, outputs, and provenance."""
    runtime = runtime_environment(settings, prepared.omp_threads)
    provenance = _solver_provenance(prepared, settings, runtime)
    process = launch_solver(prepared, runtime.values)
    observations = _observe_run(prepared, process)
    metadata = _run_metadata(
        prepared,
        runtime,
        process,
        observations,
        provenance,
        settings,
    )
    write_json_atomic(
        prepared.path / "run_metadata.json",
        metadata,
        "run metadata",
    )
    return RunResult(
        prepared.path,
        observations.status,
        process.exit_code,
        process.duration_seconds,
        observations.hdf5_outputs,
    )


def _solver_provenance(
    prepared: PreparedRun,
    settings: dict[str, str],
    runtime: RuntimeEnvironment,
) -> dict[str, Any]:
    return {
        "environment_script": (
            file_record(runtime.setup_script, str(runtime.setup_script))
            if runtime.setup_script is not None
            else None
        ),
        "executable": file_record(prepared.executable, str(prepared.executable)),
        "runtime_files": {
            name: file_record(path, str(path))
            for name, path in sorted(prepared.runtime_files.items())
        },
        "build_manifest": optional_file_record(settings, "MHDG_BUILD_MANIFEST"),
    }


def _observe_run(
    prepared: PreparedRun,
    process: SolverProcessResult,
) -> RunObservations:
    outputs = output_records(prepared.path)
    hdf5_outputs = [
        record["path"] for record in outputs if record["path"].endswith(".h5")
    ]
    fatal_messages = fatal_log_messages(
        prepared.path / "stdout.log",
        prepared.path / "stderr.log",
    )
    status = _run_status(
        process.exit_code,
        process.launch_error,
        hdf5_outputs,
        fatal_messages,
    )
    return RunObservations(outputs, hdf5_outputs, fatal_messages, status)


def _run_metadata(
    prepared: PreparedRun,
    runtime: RuntimeEnvironment,
    process: SolverProcessResult,
    observations: RunObservations,
    provenance: dict[str, Any],
    settings: dict[str, str],
) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "status": observations.status,
        "started_utc": process.started_utc,
        "finished_utc": utc_now(),
        "duration_seconds": process.duration_seconds,
        "exit_code": process.exit_code,
        "launch_error": process.launch_error,
        "fatal_log_messages": observations.fatal_log_messages,
        "working_directory": str(prepared.path),
        "environment": {
            **runtime.openmp,
            "setup_script": provenance["environment_script"],
        },
        "command": prepared.command,
        "logs": {"stdout": "stdout.log", "stderr": "stderr.log"},
        "executable": provenance["executable"],
        "runtime_files": provenance["runtime_files"],
        "solver": {
            "revision": settings.get("MHDG_SOLVER_REVISION") or None,
            "build_description": settings.get("MHDG_BUILD_DESCRIPTION") or None,
            "build_manifest": provenance["build_manifest"],
        },
        "output_files": observations.output_files,
        "hdf5_outputs": observations.hdf5_outputs,
    }


def _run_status(
    exit_code: int | None,
    launch_error: str | None,
    hdf5_outputs: list[str],
    fatal_messages: list[str],
) -> str:
    if launch_error is not None:
        return "launch_failed"
    if exit_code != 0:
        return "solver_failed"
    if fatal_messages:
        return "solver_reported_error"
    if not hdf5_outputs:
        return "missing_hdf5_output"
    return "completed"
