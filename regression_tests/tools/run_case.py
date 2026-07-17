#!/usr/bin/env python3
"""Prepare and execute one isolated MHDG regression run."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from check_bundle import BundleError, read_settings
from compare_hdf5 import ComparisonError
from compare_run import select_candidate
from prepare_run import (
    PreparedExecution,
    PreparedRun,
    PreparedStagedRun,
    openmp_environment,
    prepare_run,
)


FATAL_LOG_MARKERS = (
    "Error opening source file:",
    "Error opening destination file:",
    "Error   : Unable to open file",
)
MAX_FATAL_LOG_MESSAGES = 20


@dataclass(frozen=True)
class RunResult:
    path: Path
    status: str
    exit_code: int | None
    duration_seconds: float
    hdf5_outputs: list[str]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_id", metavar="CASE")
    parser.add_argument("workflow_ids", metavar="WORKFLOW", nargs="+")
    parser.add_argument("--layout", required=True, dest="layout_id")
    parser.add_argument("--run-id")
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--layouts", required=True, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    try:
        settings = read_settings(args.settings)
    except BundleError as exc:
        print(f"run failed: {exc}", file=sys.stderr)
        return 1

    completed = True
    for workflow_id in args.workflow_ids:
        print(f"workflow: {workflow_id}")
        try:
            prepared = prepare_run(
                args.settings,
                args.case_id,
                workflow_id,
                args.layout_id,
                args.cases,
                args.layouts,
                args.run_id,
            )
            _print_prepared(prepared)
            result = execute_prepared(prepared, settings)
        except BundleError as exc:
            print(f"run failed: {exc}", file=sys.stderr)
            completed = False
            continue
        _print_result(result)
        completed = completed and result.status == "completed"
    return 0 if completed else 1


def _print_prepared(prepared: PreparedExecution) -> None:
    print(f"run directory: {prepared.path}")
    if isinstance(prepared, PreparedRun):
        print(f"command: {shlex.join(prepared.command)}")
        print("solver output: stdout.log and stderr.log", flush=True)
    else:
        print(f"stages: {len(prepared.stages)}", flush=True)


def _print_result(result: RunResult) -> None:
    print(f"status: {result.status}")
    print(f"solver exit code: {result.exit_code}")
    print(f"runtime: {result.duration_seconds:.3f} s")
    print(f"HDF5 outputs: {len(result.hdf5_outputs)}")


def execute_prepared(
    prepared: PreparedExecution, settings: dict[str, str]
) -> RunResult:
    """Execute either one warm run or a sequential staged workflow."""
    if isinstance(prepared, PreparedStagedRun):
        return execute_staged_run(prepared, settings)
    return execute_run(prepared, settings)


def execute_run(prepared: PreparedRun, settings: dict[str, str]) -> RunResult:
    """Execute a prepared command and record logs, outputs, and provenance."""
    stdout_path = prepared.path / "stdout.log"
    stderr_path = prepared.path / "stderr.log"
    environment, environment_script = _runtime_environment(settings)
    openmp = openmp_environment(prepared.omp_threads)
    environment.update(openmp)
    executable = _file_record(prepared.executable, str(prepared.executable))
    build_manifest = _optional_file_record(settings, "MHDG_BUILD_MANIFEST")

    started_utc = _utc_now()
    started_clock = time.monotonic()
    exit_code: int | None = None
    launch_error: str | None = None
    try:
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr:
            completed = subprocess.run(
                prepared.command,
                cwd=prepared.path,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        exit_code = completed.returncode
    except OSError as exc:
        launch_error = str(exc)

    duration = time.monotonic() - started_clock
    output_files = _output_records(prepared.path)
    hdf5_outputs = [
        record["path"] for record in output_files if record["path"].endswith(".h5")
    ]
    fatal_log_messages = _fatal_log_messages(stdout_path, stderr_path)
    status = _run_status(
        exit_code, launch_error, hdf5_outputs, fatal_log_messages
    )

    metadata = {
        "schema_version": 1,
        "status": status,
        "started_utc": started_utc,
        "finished_utc": _utc_now(),
        "duration_seconds": duration,
        "exit_code": exit_code,
        "launch_error": launch_error,
        "fatal_log_messages": fatal_log_messages,
        "working_directory": str(prepared.path),
        "environment": {**openmp, "setup_script": environment_script},
        "command": prepared.command,
        "logs": {"stdout": "stdout.log", "stderr": "stderr.log"},
        "executable": executable,
        "runtime_files": {
            name: _file_record(path, str(path))
            for name, path in sorted(prepared.runtime_files.items())
        },
        "solver": {
            "revision": settings.get("MHDG_SOLVER_REVISION") or None,
            "build_description": settings.get("MHDG_BUILD_DESCRIPTION") or None,
            "build_manifest": build_manifest,
        },
        "output_files": output_files,
        "hdf5_outputs": hdf5_outputs,
    }
    _write_json(prepared.path / "run_metadata.json", metadata)

    return RunResult(prepared.path, status, exit_code, duration, hdf5_outputs)


def execute_staged_run(
    prepared: PreparedStagedRun, settings: dict[str, str]
) -> RunResult:
    """Run stages in order, passing each selected HDF5 result to the next."""
    started_utc = _utc_now()
    started_clock = time.monotonic()
    stage_records = []
    selected_output: Path | None = None
    last_result: RunResult | None = None
    status = "completed"

    for index, stage in enumerate(prepared.stages, start=1):
        if stage.restart_from == "previous_stage":
            if selected_output is None:
                raise BundleError(f"stage {stage.stage_id} has no restart source")
            _link_restart(stage.run.path, selected_output)

        print(
            f"stage {index}/{len(prepared.stages)}: {stage.stage_id}", flush=True
        )
        result = execute_run(stage.run, settings)
        last_result = result
        record = {
            "stage_id": stage.stage_id,
            "restart_from": stage.restart_from,
            "run_directory": str(stage.run.path),
            "status": result.status,
            "exit_code": result.exit_code,
            "duration_seconds": result.duration_seconds,
            "selected_hdf5": None,
        }
        stage_records.append(record)
        if result.status != "completed":
            status = result.status
            break

        try:
            selected_output = select_candidate(
                result.path, {"hdf5_outputs": result.hdf5_outputs}
            )
        except ComparisonError as exc:
            record["status"] = "output_selection_failed"
            record["selection_error"] = str(exc)
            status = "output_selection_failed"
            selected_output = None
            break
        record["selected_hdf5"] = str(selected_output)

    for stage in prepared.stages[len(stage_records) :]:
        stage_records.append(
            {
                "stage_id": stage.stage_id,
                "restart_from": stage.restart_from,
                "run_directory": str(stage.run.path),
                "status": "not_run",
                "exit_code": None,
                "duration_seconds": None,
                "selected_hdf5": None,
            }
        )

    if last_result is not None:
        _link_summary_logs(prepared.path, last_result.path)
    hdf5_outputs = []
    if status == "completed" and selected_output is not None:
        hdf5_outputs = [selected_output.relative_to(prepared.path).as_posix()]

    metadata = {
        "schema_version": 1,
        "status": status,
        "started_utc": started_utc,
        "finished_utc": _utc_now(),
        "duration_seconds": time.monotonic() - started_clock,
        "exit_code": last_result.exit_code if last_result is not None else None,
        "working_directory": str(prepared.path),
        "logs": {"stdout": "stdout.log", "stderr": "stderr.log"},
        "stages": stage_records,
        "hdf5_outputs": hdf5_outputs,
    }
    if last_result is not None:
        stage_metadata = _read_json(last_result.path / "run_metadata.json")
        for name in ("environment", "executable", "runtime_files", "solver"):
            metadata[name] = stage_metadata[name]
    _write_json(prepared.path / "run_metadata.json", metadata)

    return RunResult(
        prepared.path,
        status,
        metadata["exit_code"],
        metadata["duration_seconds"],
        hdf5_outputs,
    )


def _link_restart(run_dir: Path, source: Path) -> None:
    try:
        (run_dir / "inputs" / "restart.h5").symlink_to(source)
    except OSError as exc:
        raise BundleError(f"cannot link restart for {run_dir.name}: {exc}") from exc


def _link_summary_logs(workflow_dir: Path, stage_dir: Path) -> None:
    for filename in ("stdout.log", "stderr.log"):
        link = workflow_dir / filename
        try:
            link.symlink_to((stage_dir / filename).relative_to(workflow_dir))
        except OSError as exc:
            raise BundleError(f"cannot link workflow log {link}: {exc}") from exc


def _runtime_environment(
    settings: dict[str, str],
) -> tuple[dict[str, str], dict[str, Any] | None]:
    environment = os.environ.copy()
    configured_path = settings.get("MHDG_ENVIRONMENT_SCRIPT")
    if not configured_path:
        return environment, None

    script = Path(configured_path).expanduser()
    if not script.is_absolute() or not script.is_file():
        raise BundleError(
            "MHDG_ENVIRONMENT_SCRIPT must be an absolute path to a file"
        )

    command = [
        "bash",
        "-c",
        'source "$1" >/dev/null && env -0',
        "mhdg-regression-environment",
        str(script),
    ]
    try:
        completed = subprocess.run(
            command,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except OSError as exc:
        raise BundleError(f"cannot source environment script {script}: {exc}") from exc
    if completed.returncode != 0:
        error = completed.stderr.decode(errors="replace").strip()
        detail = f": {error}" if error else ""
        raise BundleError(f"environment script failed ({completed.returncode}){detail}")

    sourced_environment = {}
    for entry in completed.stdout.split(b"\0"):
        if not entry:
            continue
        key, value = entry.split(b"=", 1)
        sourced_environment[os.fsdecode(key)] = os.fsdecode(value)
    return sourced_environment, _file_record(script, str(script))


def _run_status(
    exit_code: int | None,
    launch_error: str | None,
    hdf5_outputs: list[str],
    fatal_log_messages: list[str],
) -> str:
    if launch_error is not None:
        return "launch_failed"
    if exit_code != 0:
        return "solver_failed"
    if fatal_log_messages:
        return "solver_reported_error"
    if not hdf5_outputs:
        return "missing_hdf5_output"
    return "completed"


def _fatal_log_messages(*paths: Path) -> list[str]:
    messages = []
    for path in paths:
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            for line_number, line in enumerate(stream, start=1):
                if not any(marker in line for marker in FATAL_LOG_MARKERS):
                    continue
                messages.append(f"{path.name}:{line_number}: {line.rstrip()}")
                if len(messages) == MAX_FATAL_LOG_MESSAGES:
                    return messages
    return messages


def _output_records(run_dir: Path) -> list[dict[str, Any]]:
    output_dir = run_dir / "outputs"
    return [
        _file_record(path, path.relative_to(run_dir).as_posix())
        for path in sorted(output_dir.rglob("*"))
        if path.is_file()
    ]


def _optional_file_record(
    settings: dict[str, str], key: str
) -> dict[str, Any] | None:
    value = settings.get(key)
    if not value:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute() or not path.is_file():
        raise BundleError(f"{key} must be an absolute path to a file")
    return _file_record(path, str(path))


def _file_record(path: Path, display_path: str) -> dict[str, Any]:
    try:
        size = path.stat().st_size
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise BundleError(f"cannot inspect file {path}: {exc}") from exc
    return {"path": display_path, "size_bytes": size, "sha256": digest.hexdigest()}


def _read_json(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BundleError(f"cannot read run metadata {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise BundleError(f"run metadata must contain an object: {path}")
    return document


def _write_json(path: Path, document: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        temporary.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)
    except OSError as exc:
        raise BundleError(f"cannot write run metadata {path}: {exc}") from exc


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


if __name__ == "__main__":
    raise SystemExit(main())
