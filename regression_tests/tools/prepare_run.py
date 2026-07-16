#!/usr/bin/env python3
"""Prepare an isolated MHDG regression run without executing the solver."""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from check_bundle import (
    BundleError,
    bundle_root_from_settings,
    load_case_definition,
    load_validated_json,
    read_settings,
    required_case_roles,
    validate_bundle_root,
)


RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
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
    except BundleError as exc:
        print(f"run preparation failed: {exc}", file=sys.stderr)
        return 1

    print(f"run prepared: {prepared.path}")
    print(f"OMP_NUM_THREADS={prepared.omp_threads}")
    print(f"command: {shlex.join(prepared.command)}")
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
) -> PreparedRun:
    """Create one validated, isolated run directory."""
    settings = read_settings(settings_path)
    bundle_root = bundle_root_from_settings(settings)
    if validate_bundle:
        validate_bundle_root(bundle_root, case_dir)

    case = load_case_definition(case_id, case_dir)
    workflow = case["workflows"].get(workflow_id)
    if workflow is None:
        raise BundleError(f"case {case_id} has no workflow {workflow_id}")
    if workflow["kind"] != "warm_same_state":
        raise BundleError("run preparation currently supports warm_same_state only")

    layout = _load_layout(layout_id, layouts_path)
    artifacts, manifest = _case_artifacts(
        bundle_root, case, workflow_id, case_dir
    )
    run_root = _absolute_setting(settings, "MHDG_REGRESSION_RUN_ROOT")
    executable = _solver_executable(settings, layout["execution"])
    runtime_files = _runtime_files(executable)
    launcher = _mpi_launcher(settings) if layout["execution"] == "mpi" else None

    run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if not RUN_ID_RE.fullmatch(run_id):
        raise BundleError(f"invalid run identifier: {run_id}")

    run_dir = run_root / case_id / workflow_id / layout_id / run_id
    if run_dir.exists():
        raise BundleError(f"run directory already exists: {run_dir}")

    command = _solver_command(run_dir, executable, launcher, layout)
    try:
        run_dir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{run_id}.", dir=run_dir.parent
        ) as workspace:
            staging = Path(workspace) / "run"
            staging.mkdir()
            _populate_run(staging, run_dir, artifacts, runtime_files)
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


def render_parameter_file(
    source: Path, destination: Path, replacements: dict[str, Path | str]
) -> None:
    """Replace selected namelist path assignments in a copied parameter file."""
    try:
        lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    except OSError as exc:
        raise BundleError(f"cannot read parameter file {source}: {exc}") from exc

    values = {key.lower(): str(value) for key, value in replacements.items()}
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

        value = values[key]
        if "'" in value:
            raise BundleError(f"cannot render a path containing a quote: {value}")
        suffix = f" !{comment}" if marker else ""
        rendered.append(f"{match.group('prefix')}'{value}'{suffix}{ending}")
        counts[key] += 1

    invalid = [key for key, count in counts.items() if count != 1]
    if invalid:
        details = ", ".join(f"{key} ({counts[key]} matches)" for key in invalid)
        raise BundleError(f"parameter path assignments must appear once: {details}")

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
    if case_data is None or case_data["case_id"] != case["case_id"]:
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


def _populate_run(
    staging: Path,
    final_run_dir: Path,
    artifacts: dict[str, Path],
    runtime_files: dict[str, Path],
) -> None:
    staging_inputs = staging / "inputs"
    staging_outputs = staging / "outputs"
    staging_inputs.mkdir()
    staging_outputs.mkdir()

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


def _solver_command(
    run_dir: Path,
    executable: Path,
    launcher: Path | None,
    layout: dict[str, Any],
) -> list[str]:
    arguments = [
        str(executable),
        str(run_dir / "inputs" / "mesh"),
        str(run_dir / "inputs" / "restart"),
    ]
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
) -> None:
    created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    plan = {
        "schema_version": 1,
        "created_utc": created.replace("+00:00", "Z"),
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
    (staging / "run_plan.json").write_text(
        json.dumps(plan, indent=2) + "\n", encoding="utf-8"
    )


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
