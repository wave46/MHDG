"""Resolve validated inputs for one requested regression run."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition, required_case_roles
from bundle.schemas import load_validated_json
from bundle.settings import bundle_root_from_settings, read_settings
from bundle.validation import validate_bundle_root
from preparation.models import PreparationInputs
from support.errors import BundleError
from support.identifiers import IDENTIFIER_RE
from support.time import utc_run_id


RUNTIME_FILENAMES = ("positionFeketeNodesTri2D.h5",)


def load_preparation_inputs(
    settings_path: Path,
    case_id: str,
    workflow_id: str,
    layout_id: str,
    case_directory: Path,
    layouts_path: Path,
    run_id: str | None,
    validate_bundle: bool,
) -> PreparationInputs:
    """Resolve the files and declarations needed to prepare one run."""
    settings = read_settings(settings_path)
    bundle_root = bundle_root_from_settings(settings)
    if validate_bundle:
        validate_bundle_root(bundle_root, case_directory)

    case = load_case_definition(case_id, case_directory)
    workflow = case["workflows"].get(workflow_id)
    if workflow is None:
        raise BundleError(f"case {case_id} has no workflow {workflow_id}")
    layout = _load_layout(layout_id, layouts_path)
    artifacts, manifest = _case_artifacts(
        bundle_root,
        case,
        workflow_id,
        case_directory,
    )

    run_root = _absolute_setting(settings, "MHDG_REGRESSION_RUN_ROOT")
    executable = _solver_executable(settings, layout["execution"])
    runtime_files = _runtime_files(executable)
    launcher = _mpi_launcher(settings) if layout["execution"] == "mpi" else None
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
        runtime_files=runtime_files,
        case=case,
        workflow_id=workflow_id,
        workflow=workflow,
        layout_id=layout_id,
        bundle_root=bundle_root,
        manifest=manifest,
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
    case_directory: Path,
) -> tuple[dict[str, Path], dict[str, Any]]:
    schema_path = case_directory.parent / "schemas" / "bundle-manifest.schema.json"
    manifest = load_validated_json(
        bundle_root / "manifest.json",
        schema_path,
        "bundle manifest",
    )
    case_id = case["case_id"]
    case_data = manifest["case_data"].get(case_id)
    if case_data is None or case_data["case_id"] != case_id:
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
