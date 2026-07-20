"""Execute staged workflows and propagate restart solutions."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from comparison.shared.outputs import select_candidate
from execution.models import RunResult
from execution.single import execute_run
from preparation.models import PreparedStage, PreparedStagedRun
from support.documents import load_json, write_json_atomic
from support.errors import BundleError, ComparisonError
from support.time import utc_now


def execute_staged_run(
    prepared: PreparedStagedRun,
    settings: dict[str, str],
) -> RunResult:
    """Run stages in order, passing each selected HDF5 result to the next."""
    started_utc = utc_now()
    started_clock = time.monotonic()
    stage_records = []
    selected_output: Path | None = None
    last_result: RunResult | None = None
    status = "completed"

    for index, stage in enumerate(prepared.stages, start=1):
        _link_stage_restart(stage, selected_output)
        print(
            f"stage {index}/{len(prepared.stages)}: {stage.stage_id}",
            flush=True,
        )
        result = execute_run(stage.run, settings)
        last_result = result
        record = _stage_record(stage, result)
        stage_records.append(record)
        if result.status != "completed":
            status = result.status
            break

        try:
            selected_output = select_candidate(
                result.path,
                {"hdf5_outputs": result.hdf5_outputs},
            )
        except ComparisonError as exc:
            record["status"] = "output_selection_failed"
            record["selection_error"] = str(exc)
            status = "output_selection_failed"
            selected_output = None
            break
        record["selected_hdf5"] = str(selected_output)

    stage_records.extend(
        _not_run_record(stage)
        for stage in prepared.stages[len(stage_records) :]
    )
    if last_result is not None:
        _link_summary_logs(prepared.path, last_result.path)
    hdf5_outputs = _workflow_outputs(prepared.path, status, selected_output)
    metadata = _staged_metadata(
        prepared,
        status,
        started_utc,
        time.monotonic() - started_clock,
        last_result,
        stage_records,
        hdf5_outputs,
    )
    write_json_atomic(
        prepared.path / "run_metadata.json",
        metadata,
        "run metadata",
    )
    return RunResult(
        prepared.path,
        status,
        metadata["exit_code"],
        metadata["duration_seconds"],
        hdf5_outputs,
    )


def _link_stage_restart(stage: PreparedStage, source: Path | None) -> None:
    if stage.restart_from != "previous_stage":
        return
    if source is None:
        raise BundleError(f"stage {stage.stage_id} has no restart source")
    try:
        (stage.run.path / "inputs" / "restart.h5").symlink_to(source)
    except OSError as exc:
        raise BundleError(
            f"cannot link restart for {stage.run.path.name}: {exc}"
        ) from exc


def _stage_record(stage: PreparedStage, result: RunResult) -> dict[str, Any]:
    return {
        "stage_id": stage.stage_id,
        "restart_from": stage.restart_from,
        "run_directory": str(stage.run.path),
        "status": result.status,
        "exit_code": result.exit_code,
        "duration_seconds": result.duration_seconds,
        "selected_hdf5": None,
    }


def _not_run_record(stage: PreparedStage) -> dict[str, Any]:
    return {
        "stage_id": stage.stage_id,
        "restart_from": stage.restart_from,
        "run_directory": str(stage.run.path),
        "status": "not_run",
        "exit_code": None,
        "duration_seconds": None,
        "selected_hdf5": None,
    }


def _workflow_outputs(
    workflow_directory: Path,
    status: str,
    selected_output: Path | None,
) -> list[str]:
    if status != "completed" or selected_output is None:
        return []
    return [selected_output.relative_to(workflow_directory).as_posix()]


def _staged_metadata(
    prepared: PreparedStagedRun,
    status: str,
    started_utc: str,
    duration_seconds: float,
    last_result: RunResult | None,
    stage_records: list[dict[str, Any]],
    hdf5_outputs: list[str],
) -> dict[str, Any]:
    metadata = {
        "schema_version": 1,
        "status": status,
        "started_utc": started_utc,
        "finished_utc": utc_now(),
        "duration_seconds": duration_seconds,
        "exit_code": last_result.exit_code if last_result is not None else None,
        "working_directory": str(prepared.path),
        "logs": {"stdout": "stdout.log", "stderr": "stderr.log"},
        "stages": stage_records,
        "hdf5_outputs": hdf5_outputs,
    }
    if last_result is not None:
        stage_metadata = load_json(
            last_result.path / "run_metadata.json",
            "stage run metadata",
        )
        for name in ("environment", "executable", "runtime_files", "solver"):
            metadata[name] = stage_metadata[name]
    return metadata


def _link_summary_logs(workflow_directory: Path, stage_directory: Path) -> None:
    for filename in ("stdout.log", "stderr.log"):
        link = workflow_directory / filename
        try:
            link.symlink_to(
                (stage_directory / filename).relative_to(workflow_directory)
            )
        except OSError as exc:
            raise BundleError(f"cannot link workflow log {link}: {exc}") from exc
