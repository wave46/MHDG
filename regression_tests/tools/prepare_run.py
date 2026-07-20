#!/usr/bin/env python3
"""Prepare an isolated MHDG regression run without executing the solver."""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from check_bundle import (
    bundle_root_from_settings,
    load_case_definition,
    load_validated_json,
    read_settings,
    required_case_roles,
    validate_bundle_root,
)
from support.documents import write_json_direct
from support.errors import BundleError, HarnessError
from support.identifiers import IDENTIFIER_RE
from support.time import utc_now, utc_run_id


ASSIGNMENT_RE = re.compile(
    r"^(?P<prefix>\s*(?P<key>[A-Za-z][A-Za-z0-9_]*)\s*=\s*).*$"
)

INPUT_LINKS = {
    "mesh": "mesh.msh",
    "geometry": "geometry.geo",
    "equilibrium_magnetic_field": "equilibrium.h5",
    "equilibrium_current_density": "current_density.h5",
    "transport_configuration": "transport_model.nml",
    "warm_restart": "restart.h5",
    "warm_reference": "reference.h5",
}
RUNTIME_FILENAMES = ("positionFeketeNodesTri2D.h5",)


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


def openmp_environment(threads: int) -> dict[str, str]:
    """Return deterministic OpenMP placement for one solver process."""
    return {
        "OMP_NUM_THREADS": str(threads),
        "OMP_PLACES": "cores",
        "OMP_PROC_BIND": "spread",
    }


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
    except HarnessError as exc:
        print(f"run preparation failed: {exc}", file=sys.stderr)
        return 1

    print(f"run prepared: {prepared.path}")
    if isinstance(prepared, PreparedRun):
        print(f"OMP_NUM_THREADS={prepared.omp_threads}")
        print(f"command: {shlex.join(prepared.command)}")
    else:
        for stage in prepared.stages:
            print(f"stage {stage.stage_id}: {shlex.join(stage.run.command)}")
    return 0


def prepare_run(
    settings_path: Path,
    case_id: str,
    workflow_id: str,
    layout_id: str,
    case_dir: Path,
    layouts_path: Path,
    run_id: str | None = None,
    validate_bundle: bool = True,
) -> PreparedExecution:
    """Create one validated, isolated run or staged workflow directory."""
    settings = read_settings(settings_path)
    bundle_root = bundle_root_from_settings(settings)
    if validate_bundle:
        validate_bundle_root(bundle_root, case_dir)

    case = load_case_definition(case_id, case_dir)
    workflow = case["workflows"].get(workflow_id)
    if workflow is None:
        raise BundleError(f"case {case_id} has no workflow {workflow_id}")
    layout = _load_layout(layout_id, layouts_path)
    artifacts, manifest = _case_artifacts(
        bundle_root, case, workflow_id, case_dir
    )
    run_root = _absolute_setting(settings, "MHDG_REGRESSION_RUN_ROOT")
    executable = _solver_executable(settings, layout["execution"])
    runtime_files = _runtime_files(executable)
    launcher = _mpi_launcher(settings) if layout["execution"] == "mpi" else None

    run_id = run_id or utc_run_id()
    if not IDENTIFIER_RE.fullmatch(run_id):
        raise BundleError(f"invalid run identifier: {run_id}")

    run_dir = run_root / case_id / workflow_id / layout_id / run_id
    if run_dir.exists():
        raise BundleError(f"run directory already exists: {run_dir}")

    if workflow["kind"] == "warm_same_state":
        return _prepare_warm_run(
            run_dir,
            executable,
            launcher,
            layout,
            artifacts,
            runtime_files,
            case,
            workflow_id,
            layout_id,
            bundle_root,
            manifest,
        )
    if workflow["kind"] in {"staged_fixed_mesh", "staged_adaptive_mesh"}:
        return _prepare_staged_run(
            run_dir,
            executable,
            launcher,
            layout,
            artifacts,
            runtime_files,
            case,
            workflow_id,
            workflow,
            layout_id,
            bundle_root,
            manifest,
        )
    raise BundleError(
        f"run preparation does not support workflow kind {workflow['kind']}"
    )


def _prepare_warm_run(
    run_dir: Path,
    executable: Path,
    launcher: Path | None,
    layout: dict[str, Any],
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
    case: dict[str, Any],
    workflow_id: str,
    layout_id: str,
    bundle_root: Path,
    manifest: dict[str, Any],
) -> PreparedRun:
    command = _solver_command(run_dir, executable, launcher, layout, restart=True)
    try:
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{run_dir.name}.", dir=run_dir.parent
        ) as workspace:
            staging = Path(workspace) / "run"
            staging.mkdir()
            _populate_warm_run(staging, run_dir, artifacts, runtime_files)
            _write_plan(
                staging,
                run_dir,
                command,
                case,
                workflow_id,
                layout_id,
                layout,
                bundle_root,
                manifest,
                artifacts,
                runtime_files,
            )
            staging.rename(run_dir)
    except OSError as exc:
        raise BundleError(f"cannot prepare run {run_dir}: {exc}") from exc

    return PreparedRun(
        run_dir, command, layout["omp_threads"], executable, runtime_files
    )


def _prepare_staged_run(
    run_dir: Path,
    executable: Path,
    launcher: Path | None,
    layout: dict[str, Any],
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
    case: dict[str, Any],
    workflow_id: str,
    workflow: dict[str, Any],
    layout_id: str,
    bundle_root: Path,
    manifest: dict[str, Any],
) -> PreparedStagedRun:
    prepared_stages = []
    try:
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{run_dir.name}.", dir=run_dir.parent
        ) as workspace:
            staging = Path(workspace) / "run"
            (staging / "inputs").mkdir(parents=True)
            (staging / "stages").mkdir()
            (staging / "inputs" / "reference.h5").symlink_to(
                artifacts[workflow["reference_role"]]
            )

            for index, stage in enumerate(workflow["stages"], start=1):
                directory_name = f"{index:02d}_{stage['stage_id']}"
                stage_dir = run_dir / "stages" / directory_name
                staging_stage = staging / "stages" / directory_name
                command = _solver_command(
                    stage_dir,
                    executable,
                    launcher,
                    layout,
                    restart=stage["restart_from"] == "previous_stage",
                )
                _populate_stage_run(
                    staging_stage,
                    stage_dir,
                    artifacts,
                    runtime_files,
                    workflow,
                    stage,
                )
                logical_overrides = _logical_overrides(workflow, stage)
                _write_plan(
                    staging_stage,
                    stage_dir,
                    command,
                    case,
                    workflow_id,
                    layout_id,
                    layout,
                    bundle_root,
                    manifest,
                    artifacts,
                    runtime_files,
                    stage,
                    logical_overrides,
                )
                prepared_stages.append(
                    PreparedStage(
                        stage["stage_id"],
                        stage["restart_from"],
                        PreparedRun(
                            stage_dir,
                            command,
                            layout["omp_threads"],
                            executable,
                            runtime_files,
                        ),
                    )
                )

            _write_staged_plan(
                staging,
                run_dir,
                case,
                workflow_id,
                layout_id,
                layout,
                bundle_root,
                manifest,
                artifacts,
                runtime_files,
                prepared_stages,
                workflow,
            )
            staging.rename(run_dir)
    except OSError as exc:
        raise BundleError(f"cannot prepare run {run_dir}: {exc}") from exc

    return PreparedStagedRun(run_dir, prepared_stages)


def render_parameter_file(
    source: Path,
    destination: Path,
    replacements: dict[str, Path | str],
    logical_overrides: dict[str, bool] | None = None,
) -> None:
    """Render selected path and logical assignments in a parameter-file copy."""
    try:
        lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    except OSError as exc:
        raise BundleError(f"cannot read parameter file {source}: {exc}") from exc

    values = {}
    for key, value in replacements.items():
        value = str(value)
        if "'" in value:
            raise BundleError(f"cannot render a path containing a quote: {value}")
        values[key.lower()] = f"'{value}'"
    for key, value in (logical_overrides or {}).items():
        normalized = key.lower()
        if normalized in values:
            raise BundleError(f"duplicate parameter replacement: {key}")
        values[normalized] = ".true." if value else ".false."
    counts = dict.fromkeys(values, 0)
    rendered = []
    for line in lines:
        body = line.rstrip("\r\n")
        ending = line[len(body) :]
        code, marker, comment = body.partition("!")
        match = ASSIGNMENT_RE.match(code)
        key = match.group("key").lower() if match else ""
        if key not in values:
            rendered.append(line)
            continue

        suffix = f" !{comment}" if marker else ""
        rendered.append(f"{match.group('prefix')}{values[key]}{suffix}{ending}")
        counts[key] += 1

    invalid = [key for key, count in counts.items() if count != 1]
    if invalid:
        details = ", ".join(f"{key} ({counts[key]} matches)" for key in invalid)
        raise BundleError(f"parameter assignments must appear once: {details}")

    try:
        destination.write_text("".join(rendered), encoding="utf-8")
    except OSError as exc:
        raise BundleError(f"cannot write parameter file {destination}: {exc}") from exc


def _load_layout(layout_id: str, layouts_path: Path) -> dict[str, Any]:
    schema_path = layouts_path.parent / "schemas" / "layouts.schema.json"
    document = load_validated_json(layouts_path, schema_path, "layout definitions")
    layout = document["layouts"].get(layout_id)
    if layout is None:
        available = ", ".join(sorted(document["layouts"]))
        raise BundleError(f"unknown layout {layout_id}; available: {available}")
    return layout


def _case_artifacts(
    bundle_root: Path,
    case: dict[str, Any],
    workflow_id: str,
    case_dir: Path,
) -> tuple[dict[str, Path], dict[str, Any]]:
    schema_path = case_dir.parent / "schemas" / "bundle-manifest.schema.json"
    manifest = load_validated_json(
        bundle_root / "manifest.json", schema_path, "bundle manifest"
    )
    data_id = case["external_data_id"]
    case_data = manifest["case_data"].get(data_id)
    if case_data is None or case_data["case_id"] not in {
        case["case_id"],
        data_id,
    }:
        raise BundleError(f"bundle does not contain case data for {case['case_id']}")

    paths = {}
    for role in required_case_roles(case, workflow_id):
        try:
            artifact_id = case_data["roles"][role]
        except KeyError as exc:
            raise BundleError(
                f"bundle does not provide role {role} for workflow {workflow_id}"
            ) from exc
        relative_path = manifest["artifacts"][artifact_id]["path"]
        paths[role] = (bundle_root / relative_path).resolve(strict=True)
    return paths, manifest


def _populate_warm_run(
    staging: Path,
    final_run_dir: Path,
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
) -> None:
    staging_inputs = staging / "inputs"
    staging_outputs = staging / "outputs"
    staging_inputs.mkdir()
    staging_outputs.mkdir()
    (staging / "res").mkdir()

    for role, filename in INPUT_LINKS.items():
        (staging_inputs / filename).symlink_to(artifacts[role])
    for filename, source in runtime_files.items():
        (staging / filename).symlink_to(source)

    final_inputs = final_run_dir / "inputs"
    replacements = {
        "transport_model_path": final_inputs / "transport_model.nml",
        "field_path": final_inputs / "equilibrium.h5",
        "jtor_path": final_inputs / "current_density.h5",
        "geometry_path": final_inputs / "geometry.geo",
        "save_folder": f"{final_run_dir / 'outputs'}/",
    }
    render_parameter_file(
        artifacts["warm_parameters"], staging / "param.txt", replacements
    )


def _populate_stage_run(
    staging: Path,
    final_run_dir: Path,
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
    workflow: dict[str, Any],
    stage: dict[str, Any],
) -> None:
    staging_inputs = staging / "inputs"
    staging_inputs.mkdir(parents=True)
    (staging / "outputs").mkdir()
    (staging / "res").mkdir()

    input_sources = {
        "mesh.msh": artifacts[workflow["mesh_role"]],
        "geometry.geo": artifacts["geometry"],
        "equilibrium.h5": artifacts["equilibrium_magnetic_field"],
        "current_density.h5": artifacts["equilibrium_current_density"],
        "transport_model.nml": artifacts[stage["transport_configuration_role"]],
    }
    for filename, source in input_sources.items():
        (staging_inputs / filename).symlink_to(source)
    for filename, source in runtime_files.items():
        (staging / filename).symlink_to(source)

    final_inputs = final_run_dir / "inputs"
    replacements = {
        "transport_model_path": final_inputs / "transport_model.nml",
        "field_path": final_inputs / "equilibrium.h5",
        "jtor_path": final_inputs / "current_density.h5",
        "geometry_path": final_inputs / "geometry.geo",
        "save_folder": f"{final_run_dir / 'outputs'}/",
    }
    render_parameter_file(
        artifacts[stage["parameter_role"]],
        staging / "param.txt",
        replacements,
        _logical_overrides(workflow, stage),
    )


def _logical_overrides(
    workflow: dict[str, Any], stage: dict[str, Any]
) -> dict[str, bool]:
    return {
        **workflow.get("logical_overrides", {}),
        **stage.get("logical_overrides", {}),
    }


def _solver_command(
    run_dir: Path,
    executable: Path,
    launcher: Path | None,
    layout: dict[str, Any],
    restart: bool,
) -> list[str]:
    arguments = [
        str(executable),
        str(run_dir / "inputs" / "mesh"),
    ]
    if restart:
        arguments.append(str(run_dir / "inputs" / "restart"))
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


def _write_plan(
    staging: Path,
    run_dir: Path,
    command: list[str],
    case: dict[str, Any],
    workflow_id: str,
    layout_id: str,
    layout: dict[str, Any],
    bundle_root: Path,
    manifest: dict[str, Any],
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
    stage: dict[str, Any] | None = None,
    logical_overrides: dict[str, bool] | None = None,
) -> None:
    plan = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "case_id": case["case_id"],
        "workflow_id": workflow_id,
        "layout_id": layout_id,
        "layout": layout,
        "working_directory": str(run_dir),
        "environment": openmp_environment(layout["omp_threads"]),
        "command": command,
        "bundle": {
            "root": str(bundle_root),
            "bundle_id": manifest["bundle_id"],
            "bundle_version": manifest["bundle_version"],
        },
        "artifacts": {role: str(path) for role, path in sorted(artifacts.items())},
        "runtime_files": {
            name: str(path) for name, path in sorted(runtime_files.items())
        },
    }
    if stage is not None:
        plan["stage_id"] = stage["stage_id"]
        plan["restart_from"] = stage["restart_from"]
        plan["logical_overrides"] = logical_overrides or {}
    write_json_direct(staging / "run_plan.json", plan)


def _write_staged_plan(
    staging: Path,
    run_dir: Path,
    case: dict[str, Any],
    workflow_id: str,
    layout_id: str,
    layout: dict[str, Any],
    bundle_root: Path,
    manifest: dict[str, Any],
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
    stages: list[PreparedStage],
    workflow: dict[str, Any],
) -> None:
    plan = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "case_id": case["case_id"],
        "workflow_id": workflow_id,
        "workflow_kind": workflow["kind"],
        "layout_id": layout_id,
        "layout": layout,
        "working_directory": str(run_dir),
        "environment": openmp_environment(layout["omp_threads"]),
        "bundle": {
            "root": str(bundle_root),
            "bundle_id": manifest["bundle_id"],
            "bundle_version": manifest["bundle_version"],
        },
        "artifacts": {role: str(path) for role, path in sorted(artifacts.items())},
        "runtime_files": {
            name: str(path) for name, path in sorted(runtime_files.items())
        },
        "stages": [
            {
                "stage_id": stage.stage_id,
                "restart_from": stage.restart_from,
                "working_directory": str(stage.run.path),
                "command": stage.run.command,
                "logical_overrides": _logical_overrides(
                    workflow, workflow["stages"][index]
                ),
            }
            for index, stage in enumerate(stages)
        ],
    }
    write_json_direct(staging / "run_plan.json", plan)


def _absolute_setting(settings: dict[str, str], key: str) -> Path:
    value = settings.get(key)
    if not value:
        raise BundleError(f"settings must define {key}")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise BundleError(f"{key} must be an absolute path")
    return path.resolve()


def _solver_executable(settings: dict[str, str], execution: str) -> Path:
    key = (
        "MHDG_SERIAL_EXECUTABLE"
        if execution == "serial"
        else "MHDG_PARALLEL_EXECUTABLE"
    )
    path = _absolute_setting(settings, key)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise BundleError(f"{key} is not an executable file: {path}")
    return path


def _runtime_files(executable: Path) -> dict[str, Path]:
    files = {}
    for filename in RUNTIME_FILENAMES:
        path = executable.parent / filename
        if not path.is_file():
            raise BundleError(f"required runtime file is missing: {path}")
        files[filename] = path.resolve()
    return files


def _mpi_launcher(settings: dict[str, str]) -> Path:
    value = settings.get("MHDG_MPI_LAUNCHER")
    if not value:
        raise BundleError("settings must define MHDG_MPI_LAUNCHER")
    resolved = shutil.which(value)
    if not resolved:
        raise BundleError(f"MPI launcher is not executable or not found: {value}")
    return Path(resolved).resolve()


if __name__ == "__main__":
    raise SystemExit(main())
