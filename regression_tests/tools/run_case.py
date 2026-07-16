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
from prepare_run import PreparedRun, openmp_environment, prepare_run


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
    parser.add_argument("workflow_id", metavar="WORKFLOW")
    parser.add_argument("--layout", required=True, dest="layout_id")
    parser.add_argument("--run-id")
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--layouts", required=True, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    try:
        prepared = prepare_run(
            args.settings,
            args.case_id,
            args.workflow_id,
            args.layout_id,
            args.cases,
            args.layouts,
            args.run_id,
        )
        print(f"run directory: {prepared.path}")
        print(f"command: {shlex.join(prepared.command)}")
        print("solver output: stdout.log and stderr.log", flush=True)
        result = execute_run(prepared, read_settings(args.settings))
    except BundleError as exc:
        print(f"run failed: {exc}", file=sys.stderr)
        return 1

    print(f"status: {result.status}")
    print(f"solver exit code: {result.exit_code}")
    print(f"runtime: {result.duration_seconds:.3f} s")
    print(f"HDF5 outputs: {len(result.hdf5_outputs)}")
    return 0 if result.status == "completed" else 1


def execute_run(prepared: PreparedRun, settings: dict[str, str]) -> RunResult:
    """Execute a prepared command and record logs, outputs, and provenance."""
    stdout_path = prepared.path / "stdout.log"
    stderr_path = prepared.path / "stderr.log"
    environment, environment_script = _runtime_environment(settings)
    openmp = openmp_environment(prepared.omp_threads)
    environment.update(openmp)
    executable = _file_record(prepared.executable, str(prepared.executable))

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
    status = _run_status(exit_code, launch_error, hdf5_outputs)

    metadata = {
        "schema_version": 1,
        "status": status,
        "started_utc": started_utc,
        "finished_utc": _utc_now(),
        "duration_seconds": duration,
        "exit_code": exit_code,
        "launch_error": launch_error,
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
        },
        "output_files": output_files,
        "hdf5_outputs": hdf5_outputs,
    }
    _write_json(prepared.path / "run_metadata.json", metadata)

    return RunResult(prepared.path, status, exit_code, duration, hdf5_outputs)


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
    exit_code: int | None, launch_error: str | None, hdf5_outputs: list[str]
) -> str:
    if launch_error is not None:
        return "launch_failed"
    if exit_code != 0:
        return "solver_failed"
    if not hdf5_outputs:
        return "missing_hdf5_output"
    return "completed"


def _output_records(run_dir: Path) -> list[dict[str, Any]]:
    output_dir = run_dir / "outputs"
    return [
        _file_record(path, path.relative_to(run_dir).as_posix())
        for path in sorted(output_dir.rglob("*"))
        if path.is_file()
    ]


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
