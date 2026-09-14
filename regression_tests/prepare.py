"""Prepare isolated warm/cold runs from case-owned inputs and shared workflows."""

from __future__ import annotations

import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from .catalog import load_case_definition, workflow_required_roles, load_layout
from .documents import write_json_direct
from .config import (bundle_root_from_settings, absolute_setting, solver_executable,
                     runtime_files, selected_mpi_launcher)
from .bundles import validate_bundle_root, load_manifest, artifact_path
from .support import BundleError, IDENTIFIER_RE, utc_now, utc_run_id
from .parameters import render_parameter_file


@dataclass(frozen=True)
class PreparedRun:
    path: Path
    command: list[str]
    omp_threads: int
    executable: Path
    runtime_files: dict[str, Path]


@dataclass(frozen=True)
class PreparedStage:
    stage_id: str
    restart_from: str
    run: PreparedRun


@dataclass(frozen=True)
class PreparedStagedRun:
    path: Path
    stages: list[PreparedStage]


PreparedExecution = PreparedRun | PreparedStagedRun


@dataclass(frozen=True)
class PreparationInputs:
    run_directory: Path
    executable: Path
    launcher: Path | None
    layout: dict[str, Any]
    artifacts: dict[str, Path]
    runtime_files: dict[str, Path]
    case: dict[str, Any]
    workflow_id: str
    workflow: dict[str, Any]
    layout_id: str
    bundle_root: Path
    manifest: dict[str, Any]
    requested_overrides: dict[str, bool | float | int | str]


def openmp_environment(threads: int) -> dict[str, str]:
    """Return deterministic OpenMP placement for one solver process."""
    return {
        "OMP_NUM_THREADS": str(threads),
        "OMP_PLACES": "cores",
        "OMP_PROC_BIND": "spread",
    }


def solver_command(
    run_directory: Path,
    executable: Path,
    launcher: Path | None,
    layout: dict[str, Any],
    restart: bool,
) -> list[str]:
    """Build the serial or MPI solver command for one prepared run."""
    arguments = [
        str(executable),
        str(run_directory / "inputs" / "mesh"),
    ]
    if restart:
        arguments.append(str(run_directory / "inputs" / "restart"))
    if launcher is None:
        return arguments
    return [
        str(launcher),
        "--bind-to",
        "core",
        "--map-by",
        f"slot:PE={layout['omp_threads']}",
        "-n",
        str(layout["mpi_ranks"]),
        *arguments,
    ]


def load_preparation_inputs(
    settings: dict[str, Any],
    case_id: str,
    workflow_id: str,
    layout_id: str,
    case_directory: Path,
    layouts_path: Path,
    run_id: str | None,
    validate_bundle: bool,
    requested_overrides: dict[str, bool | float | int | str] | None = None,
    artifact_overrides: dict[str, Path] | None = None,
    require_reference: bool = True,
) -> PreparationInputs:
    """Resolve the files and declarations needed to prepare one run."""
    bundle_root = bundle_root_from_settings(settings)
    if validate_bundle:
        validate_bundle_root(bundle_root, case_directory)

    case = load_case_definition(case_id, case_directory)
    workflow = case["workflows"].get(workflow_id)
    if workflow is None:
        raise BundleError(f"case {case_id} has no workflow {workflow_id}")
    if not require_reference:
        workflow = {key: value for key, value in workflow.items() if key != "reference"}
    layout = load_layout(layout_id, layouts_path)
    artifacts, manifest = _case_artifacts(
        bundle_root,
        case,
        workflow_id,
        case_directory,
        workflow,
        artifact_overrides or {},
    )

    run_root = absolute_setting(settings, "MHDG_REGRESSION_RUN_ROOT")
    executable = solver_executable(settings, layout["execution"], workflow["model"])
    runtime_inputs = runtime_files(executable)
    launcher = selected_mpi_launcher(settings) if layout["execution"] == "mpi" else None
    run_directory = _run_directory(
        run_root,
        case_id,
        workflow_id,
        layout_id,
        run_id or utc_run_id(),
    )
    return PreparationInputs(
        run_directory=run_directory,
        executable=executable,
        launcher=launcher,
        layout=layout,
        artifacts=artifacts,
        runtime_files=runtime_inputs,
        case=case,
        workflow_id=workflow_id,
        workflow=workflow,
        layout_id=layout_id,
        bundle_root=bundle_root,
        manifest=manifest,
        requested_overrides=dict(requested_overrides or {}),
    )


def _run_directory(
    run_root: Path,
    case_id: str,
    workflow_id: str,
    layout_id: str,
    run_id: str,
) -> Path:
    if not IDENTIFIER_RE.fullmatch(run_id):
        raise BundleError(f"invalid run identifier: {run_id}")
    run_directory = run_root / case_id / workflow_id / layout_id / run_id
    if run_directory.exists():
        raise BundleError(f"run directory already exists: {run_directory}")
    return run_directory


def _case_artifacts(
    bundle_root: Path,
    case: dict[str, Any],
    workflow_id: str,
    case_directory: Path,
    workflow: dict[str, Any],
    overrides: dict[str, Path],
) -> tuple[dict[str, Path], dict[str, Any]]:
    manifest = load_manifest(bundle_root, case_directory, case_id=case["case_id"])

    paths = {}
    for role in workflow_required_roles(workflow):
        if role in overrides:
            paths[role] = overrides[role].resolve(strict=True)
            if not paths[role].is_file():
                raise BundleError(f"producer input is not a file: {paths[role]}")
            continue
        try:
            artifact_id = manifest["roles"][role]
        except KeyError as exc:
            raise BundleError(
                f"bundle does not provide role {role} for workflow {workflow_id}"
            ) from exc
        relative_path = manifest["artifacts"][artifact_id]["path"]
        paths[role] = artifact_path(bundle_root, relative_path, f"workflow role {role}")
    return paths, manifest




WARM_INPUT_LINKS = {
    "mesh": "mesh.msh",
    "geometry": "geometry.geo",
    "equilibrium_magnetic_field": "equilibrium.h5",
    "equilibrium_current_density": "current_density.h5",
    "transport_configuration": "transport_model.nml",
}


@contextmanager
def temporary_run_directory(final_directory: Path) -> Iterator[Path]:
    """Build a run beside its destination and publish it atomically."""
    try:
        final_directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{final_directory.name}.",
            dir=final_directory.parent,
        ) as workspace:
            staging = Path(workspace) / "run"
            staging.mkdir()
            yield staging
            staging.rename(final_directory)
    except OSError as exc:
        raise BundleError(
            f"cannot prepare run {final_directory}: {exc}"
        ) from exc


def _populate_run(staging: Path, inputs: PreparationInputs, overrides, stage=None):
    """Link immutable inputs and render only the mutable parameter file."""
    directory = staging / "inputs"
    directory.mkdir(parents=True)
    (staging / "outputs").mkdir()
    (staging / "res").mkdir()
    workflow, artifacts = inputs.workflow, inputs.artifacts
    roles = WARM_INPUT_LINKS if stage is None else {
        role: filename for role, filename in WARM_INPUT_LINKS.items()
        if role not in {"mesh", "transport_configuration"}
    }
    sources = {filename: artifacts[role] for role, filename in roles.items()}
    if stage is None:
        parameter_role = "warm_parameters"
        sources["restart.h5"] = artifacts[workflow.get("restart", "warm_restart")]
        if workflow.get("reference"):
            sources["reference.h5"] = artifacts[workflow["reference"]]
    else:
        parameter_role = stage["parameters"]
        sources["mesh.msh"] = artifacts[workflow["mesh"]]
        sources["transport_model.nml"] = artifacts[stage["transport"]]
    if workflow.get("impurity_configuration"):
        sources["impurity_model.nml"] = artifacts[workflow["impurity_configuration"]]
    for filename, source in sources.items():
        (directory / filename).symlink_to(source)
    for filename, source in inputs.runtime_files.items():
        (staging / filename).symlink_to(source)
    render_parameter_file(
        artifacts[parameter_role], staging / "param.txt",
        _parameter_replacements(inputs.run_directory), overrides,
        {**workflow.get("parameter_namelists", {}), **(stage or {}).get("parameter_namelists", {})},
    )


def _parameter_replacements(final_directory: Path) -> dict[str, Path | str]:
    inputs = final_directory / "inputs"
    return {
        "transport_model_path": inputs / "transport_model.nml",
        "impurity_model_path": inputs / "impurity_model.nml",
        "field_path": inputs / "equilibrium.h5",
        "jtor_path": inputs / "current_density.h5",
        "geometry_path": inputs / "geometry.geo",
        "save_folder": f"{final_directory / 'outputs'}/",
    }


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
        plan["stage_id"] = stage["id"]
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
        workflow_kind=inputs.workflow["type"],
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
    return {
        "schema_version": 2, "created_utc": utc_now(),
        "case_id": inputs.case["case_id"], "workflow_id": inputs.workflow_id,
        **({"workflow_kind": workflow_kind} if workflow_kind is not None else {}),
        "layout_id": inputs.layout_id, "layout": inputs.layout, "model": inputs.workflow["model"],
        "working_directory": str(run_directory),
        "environment": openmp_environment(inputs.layout["omp_threads"]),
        **({"command": command} if command is not None else {}),
        "bundle": {
            "root": str(inputs.bundle_root), "bundle_id": inputs.manifest["bundle_id"],
            "bundle_version": inputs.manifest["bundle_version"],
        },
        "artifacts": {role: str(path) for role, path in sorted(inputs.artifacts.items())},
        "runtime_files": {name: str(path) for name, path in sorted(inputs.runtime_files.items())},
    }



def prepare_warm_run(inputs: PreparationInputs) -> PreparedRun:
    """Prepare one warm-restart run, with a reference when declared."""
    overrides = {
        **inputs.workflow.get("parameter_overrides", {}),
        **inputs.requested_overrides,
    }
    command = solver_command(
        inputs.run_directory,
        inputs.executable,
        inputs.launcher,
        inputs.layout,
        restart=True,
    )
    with temporary_run_directory(inputs.run_directory) as staging:
        _populate_run(staging, inputs, overrides)
        write_run_plan(staging, inputs, command, applied_overrides=overrides)
    return _prepared_run(inputs, command)


def prepare_staged_run(inputs: PreparationInputs) -> PreparedStagedRun:
    """Prepare every directory in a cold staged workflow."""
    prepared_stages = []
    with temporary_run_directory(inputs.run_directory) as staging:
        (staging / "inputs").mkdir()
        (staging / "stages").mkdir()
        reference_role = inputs.workflow.get("reference")
        if reference_role is not None:
            (staging / "inputs" / "reference.h5").symlink_to(
                inputs.artifacts[reference_role]
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
    directory_name = f"{index:02d}_{stage['id']}"
    stage_directory = inputs.run_directory / "stages" / directory_name
    staging_stage = staging / "stages" / directory_name
    command = solver_command(
        stage_directory,
        inputs.executable,
        inputs.launcher,
        inputs.layout,
        restart=stage["restart_from"] == "previous_stage",
    )
    overrides = parameter_overrides(
        inputs.workflow,
        stage,
        inputs.requested_overrides,
    )
    stage_inputs = replace(inputs, run_directory=stage_directory)
    _populate_run(staging_stage, stage_inputs, overrides, stage)
    write_run_plan(
        staging_stage,
        stage_inputs,
        command,
        stage,
        applied_overrides=overrides,
    )
    return PreparedStage(
        stage["id"],
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


def prepare_run(
    settings: dict[str, Any],
    case_id: str,
    workflow_id: str,
    layout_id: str,
    case_dir: Path,
    layouts_path: Path,
    run_id: str | None = None,
    validate_bundle: bool = True,
    requested_overrides: dict[str, bool | float | int | str] | None = None,
    *, artifact_overrides: dict[str, Path] | None = None,
    require_reference: bool = True,
) -> PreparedExecution:
    """Create one validated, isolated run or staged workflow directory."""
    inputs = load_preparation_inputs(
        settings,
        case_id,
        workflow_id,
        layout_id,
        case_dir,
        layouts_path,
        run_id,
        validate_bundle,
        requested_overrides,
        artifact_overrides,
        require_reference,
    )
    workflow_kind = inputs.workflow["type"]
    if workflow_kind == "warm_same_state":
        return prepare_warm_run(inputs)
    if workflow_kind in {"staged_fixed_mesh", "staged_adaptive_mesh"}:
        return prepare_staged_run(inputs)
    raise BundleError(
        f"run preparation does not support workflow kind {workflow_kind}"
    )
