"""Execute prepared workflows, recording outputs, provenance and stage failures."""

from __future__ import annotations

import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from comparison.shared.outputs import select_candidate
from support.documents import load_json, write_json_atomic
from support.environments import source_environment
from support.errors import BundleError, ComparisonError
from support.files import file_identity
from support.time import utc_now
from .prepare import PreparedExecution, PreparedRun, PreparedStagedRun, openmp_environment


@dataclass(frozen=True)
class RunResult:
    path: Path
    status: str
    exit_code: int | None
    duration_seconds: float
    hdf5_outputs: list[str]


def execute_prepared(prepared: PreparedExecution, settings: dict[str, str]) -> RunResult:
    """Run one warm command or a sequential cold workflow."""
    if isinstance(prepared, PreparedStagedRun):
        return _execute_staged(prepared, settings)
    return _execute_run(prepared, settings)


def _execute_run(prepared: PreparedRun, settings: dict[str, str]) -> RunResult:
    script = None
    environment = dict(os.environ)
    if configured := settings.get("MHDG_ENVIRONMENT_SCRIPT"):
        script = Path(configured).expanduser()
        if not script.is_absolute():
            raise BundleError("MHDG_ENVIRONMENT_SCRIPT must be an absolute path")
        script, environment = source_environment(script)
    openmp = openmp_environment(prepared.omp_threads)
    environment.update(openmp)

    # Capture input identities before the process can change any files.
    provenance = {
        "environment": {
            **openmp,
            "setup_script": _file_record(script, str(script)) if script else None,
        },
        "executable": _file_record(prepared.executable, str(prepared.executable)),
        "runtime_files": {
            name: _file_record(path, str(path))
            for name, path in sorted(prepared.runtime_files.items())
        },
        "solver": {
            "revision": settings.get("MHDG_SOLVER_REVISION") or None,
            "build_description": settings.get("MHDG_BUILD_DESCRIPTION") or None,
            "build_manifest": _optional_file_record(settings, "MHDG_BUILD_MANIFEST"),
        },
    }
    started_utc = utc_now()
    started_clock = time.monotonic()
    exit_code = None
    launch_error = None
    try:
        with (prepared.path / "stdout.log").open("w", encoding="utf-8") as stdout, (
            prepared.path / "stderr.log"
        ).open("w", encoding="utf-8") as stderr:
            process = subprocess.run(
                prepared.command, cwd=prepared.path, env=environment,
                stdout=stdout, stderr=stderr, check=False,
            )
        exit_code = process.returncode
    except OSError as exc:
        launch_error = str(exc)
    duration = time.monotonic() - started_clock

    outputs = [
        _file_record(path, path.relative_to(prepared.path).as_posix())
        for path in sorted((prepared.path / "outputs").rglob("*"))
        if path.is_file()
    ]
    hdf5_outputs = [record["path"] for record in outputs if record["path"].endswith(".h5")]
    fatal_messages = _fatal_log_messages(prepared.path, missing_ok=launch_error is not None)
    if launch_error is not None:
        status = "launch_failed"
    elif exit_code != 0:
        status = "solver_failed"
    elif fatal_messages:
        status = "solver_reported_error"
    elif not hdf5_outputs:
        status = "missing_hdf5_output"
    else:
        status = "completed"
    metadata = {
        "schema_version": 2,
        "status": status,
        "started_utc": started_utc,
        "finished_utc": utc_now(),
        "duration_seconds": duration,
        "exit_code": exit_code,
        "launch_error": launch_error,
        "fatal_log_messages": fatal_messages,
        "working_directory": str(prepared.path),
        "command": prepared.command,
        "logs": {"stdout": "stdout.log", "stderr": "stderr.log"},
        **provenance,
        "output_files": outputs,
        "hdf5_outputs": hdf5_outputs,
    }
    write_json_atomic(prepared.path / "run_metadata.json", metadata, "run metadata")
    return RunResult(prepared.path, status, exit_code, duration, hdf5_outputs)


def _execute_staged(prepared: PreparedStagedRun, settings: dict[str, str]) -> RunResult:
    started_utc = utc_now()
    started_clock = time.monotonic()
    records = [
        {
            "stage_id": stage.stage_id,
            "restart_from": stage.restart_from,
            "run_directory": str(stage.run.path),
            "status": "not_run",
            "exit_code": None,
            "duration_seconds": None,
            "selected_hdf5": None,
        }
        for stage in prepared.stages
    ]
    selected_output = None
    last_result = None
    status = "completed"
    for index, (stage, record) in enumerate(zip(prepared.stages, records), start=1):
        if stage.restart_from == "previous_stage":
            if selected_output is None:
                raise BundleError(f"stage {stage.stage_id} has no restart source")
            try:
                (stage.run.path / "inputs/restart.h5").symlink_to(selected_output)
            except OSError as exc:
                raise BundleError(f"cannot link restart for {stage.run.path.name}: {exc}") from exc
        print(f"stage {index}/{len(prepared.stages)}: {stage.stage_id}", flush=True)
        last_result = _execute_run(stage.run, settings)
        status = last_result.status
        record.update(
            status=status, exit_code=last_result.exit_code,
            duration_seconds=last_result.duration_seconds,
        )
        if status != "completed":
            break
        try:
            selected_output = select_candidate(
                last_result.path, {"hdf5_outputs": last_result.hdf5_outputs},
            )
        except ComparisonError as exc:
            status = "output_selection_failed"
            record.update(status=status, selection_error=str(exc))
            selected_output = None
            break
        record["selected_hdf5"] = str(selected_output)

    if last_result is not None:
        for filename in ("stdout.log", "stderr.log"):
            link = prepared.path / filename
            try:
                link.symlink_to((last_result.path / filename).relative_to(prepared.path))
            except OSError as exc:
                raise BundleError(f"cannot link workflow log {link}: {exc}") from exc
    hdf5_outputs = (
        [selected_output.relative_to(prepared.path).as_posix()]
        if status == "completed" and selected_output is not None else []
    )
    duration = time.monotonic() - started_clock
    exit_code = last_result.exit_code if last_result else None
    metadata = {
        "schema_version": 2,
        "status": status,
        "started_utc": started_utc,
        "finished_utc": utc_now(),
        "duration_seconds": duration,
        "exit_code": exit_code,
        "working_directory": str(prepared.path),
        "logs": {"stdout": "stdout.log", "stderr": "stderr.log"},
        "stages": records,
        "hdf5_outputs": hdf5_outputs,
    }
    if last_result is not None:
        last_metadata = load_json(last_result.path / "run_metadata.json", "stage run metadata")
        for name in ("environment", "executable", "runtime_files", "solver"):
            metadata[name] = last_metadata[name]
    write_json_atomic(prepared.path / "run_metadata.json", metadata, "run metadata")
    return RunResult(prepared.path, status, exit_code, duration, hdf5_outputs)


def _file_record(path: Path, display_path: str) -> dict[str, Any]:
    return {"path": display_path, **file_identity(path)}


def _optional_file_record(settings: dict[str, str], key: str) -> dict[str, Any] | None:
    if not (value := settings.get(key)):
        return None
    path = Path(value).expanduser()
    if not path.is_absolute() or not path.is_file():
        raise BundleError(f"{key} must be an absolute path to a file")
    return _file_record(path, str(path))


def _fatal_log_messages(run_directory: Path, *, missing_ok: bool = False) -> list[str]:
    """Catch solver file errors even when the process returns zero; bound the report."""
    markers = (
        "Error opening source file:",
        "Error opening destination file:",
        "Error   : Unable to open file",
    )
    messages = []
    for filename in ("stdout.log", "stderr.log"):
        path = run_directory / filename
        if missing_ok and not path.is_file():
            continue  # Opening a log may itself have failed during launch.
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            for line_number, line in enumerate(stream, start=1):
                if any(marker in line for marker in markers):
                    messages.append(f"{filename}:{line_number}: {line.rstrip()}")
                    if len(messages) == 20:
                        return messages
    return messages
