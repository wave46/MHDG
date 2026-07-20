"""Prepare warm and staged regression workflow directories."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any

from preparation.commands import solver_command
from preparation.models import (
    PreparationInputs,
    PreparedRun,
    PreparedStage,
    PreparedStagedRun,
)
from preparation.plans import logical_overrides, write_run_plan, write_staged_plan
from preparation.workspace import (
    populate_stage_run,
    populate_warm_run,
    temporary_run_directory,
)


def prepare_warm_run(inputs: PreparationInputs) -> PreparedRun:
    """Prepare one warm restart and fixed-reference comparison run."""
    command = solver_command(
        inputs.run_directory,
        inputs.executable,
        inputs.launcher,
        inputs.layout,
        restart=True,
    )
    with temporary_run_directory(inputs.run_directory) as staging:
        populate_warm_run(
            staging,
            inputs.run_directory,
            inputs.artifacts,
            inputs.runtime_files,
        )
        write_run_plan(staging, inputs, command)
    return _prepared_run(inputs, command)


def prepare_staged_run(inputs: PreparationInputs) -> PreparedStagedRun:
    """Prepare every directory in a cold staged workflow."""
    prepared_stages = []
    with temporary_run_directory(inputs.run_directory) as staging:
        (staging / "inputs").mkdir()
        (staging / "stages").mkdir()
        (staging / "inputs" / "reference.h5").symlink_to(
            inputs.artifacts[inputs.workflow["reference_role"]]
        )
        for index, stage in enumerate(inputs.workflow["stages"], start=1):
            prepared_stages.append(_prepare_stage(inputs, staging, index, stage))
        write_staged_plan(staging, inputs, prepared_stages)
    return PreparedStagedRun(inputs.run_directory, prepared_stages)


def _prepare_stage(
    inputs: PreparationInputs,
    staging: Path,
    index: int,
    stage: dict[str, Any],
) -> PreparedStage:
    directory_name = f"{index:02d}_{stage['stage_id']}"
    stage_directory = inputs.run_directory / "stages" / directory_name
    staging_stage = staging / "stages" / directory_name
    command = solver_command(
        stage_directory,
        inputs.executable,
        inputs.launcher,
        inputs.layout,
        restart=stage["restart_from"] == "previous_stage",
    )
    overrides = logical_overrides(inputs.workflow, stage)
    populate_stage_run(
        staging_stage,
        stage_directory,
        inputs.artifacts,
        inputs.runtime_files,
        inputs.workflow,
        stage,
        overrides,
    )
    stage_inputs = replace(inputs, run_directory=stage_directory)
    write_run_plan(
        staging_stage,
        stage_inputs,
        command,
        stage,
        overrides,
    )
    return PreparedStage(
        stage["stage_id"],
        stage["restart_from"],
        _prepared_run(stage_inputs, command),
    )


def _prepared_run(
    inputs: PreparationInputs,
    command: list[str],
) -> PreparedRun:
    return PreparedRun(
        inputs.run_directory,
        command,
        inputs.layout["omp_threads"],
        inputs.executable,
        inputs.runtime_files,
    )
