"""Construct stable JSON plans for prepared runs and staged workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from preparation.commands import openmp_environment
from preparation.models import PreparationInputs, PreparedStage
from support.documents import write_json_direct
from support.time import utc_now


def write_run_plan(
    staging_directory: Path,
    inputs: PreparationInputs,
    command: list[str],
    stage: dict[str, Any] | None = None,
    applied_overrides: dict[str, Any] | None = None,
) -> None:
    """Write the plan for one warm run or one workflow stage."""
    plan = _base_plan(inputs, inputs.run_directory, command=command)
    if stage is not None:
        plan["stage_id"] = stage["stage_id"]
        plan["restart_from"] = stage["restart_from"]
    if applied_overrides is not None:
        plan["parameter_overrides"] = applied_overrides
    write_json_direct(staging_directory / "run_plan.json", plan)


def write_staged_plan(
    staging_directory: Path,
    inputs: PreparationInputs,
    stages: list[PreparedStage],
) -> None:
    """Write the parent plan describing every prepared workflow stage."""
    plan = _base_plan(
        inputs,
        inputs.run_directory,
        workflow_kind=inputs.workflow["kind"],
    )
    plan["stages"] = [
        {
            "stage_id": stage.stage_id,
            "restart_from": stage.restart_from,
            "working_directory": str(stage.run.path),
            "command": stage.run.command,
            "parameter_overrides": parameter_overrides(
                inputs.workflow,
                inputs.workflow["stages"][index],
                inputs.requested_overrides,
            ),
        }
        for index, stage in enumerate(stages)
    ]
    write_json_direct(staging_directory / "run_plan.json", plan)


def parameter_overrides(
    workflow: dict[str, Any],
    stage: dict[str, Any],
    requested: dict[str, bool | float | int | str] | None = None,
) -> dict[str, Any]:
    """Merge workflow-wide and stage-specific parameter values."""
    return {
        **workflow.get("parameter_overrides", {}),
        **stage.get("parameter_overrides", {}),
        **(requested or {}),
    }


def _base_plan(
    inputs: PreparationInputs,
    run_directory: Path,
    workflow_kind: str | None = None,
    command: list[str] | None = None,
) -> dict[str, Any]:
    plan = {
        "schema_version": 2,
        "created_utc": utc_now(),
        "case_id": inputs.case["case_id"],
        "workflow_id": inputs.workflow_id,
    }
    if workflow_kind is not None:
        plan["workflow_kind"] = workflow_kind
    plan.update(
        {
            "layout_id": inputs.layout_id,
            "layout": inputs.layout,
            "working_directory": str(run_directory),
            "environment": openmp_environment(inputs.layout["omp_threads"]),
        }
    )
    if command is not None:
        plan["command"] = command
    plan.update(
        {
            "bundle": {
                "root": str(inputs.bundle_root),
                "bundle_id": inputs.manifest["bundle_id"],
                "bundle_version": inputs.manifest["bundle_version"],
            },
            "artifacts": {
                role: str(path) for role, path in sorted(inputs.artifacts.items())
            },
            "runtime_files": {
                name: str(path)
                for name, path in sorted(inputs.runtime_files.items())
            },
        }
    )
    return plan
