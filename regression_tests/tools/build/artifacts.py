"""Publish build artifacts, provenance, and generated settings."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any

from build.models import BuildConfiguration, BuildError
from build.process import toolchain_record
from support.documents import write_json_direct
from support.files import file_identity
from support.time import utc_now


def install_executable(source: Path, destination: Path, label: str) -> Path:
    """Validate and copy one executable into the isolated build directory."""
    if not source.is_file() or not os.access(source, os.X_OK):
        raise BuildError(f"{label} was not produced: {source}")
    shutil.copy2(source, destination)
    return destination


def install_runtime_file(
    configuration: BuildConfiguration,
    filename: str,
) -> Path:
    """Copy one solver runtime file beside the built executables."""
    source = configuration.test_directory / filename
    if not source.is_file():
        raise BuildError(f"generic Fekete-node data does not exist: {source}")
    destination = configuration.binary_directory / filename
    shutil.copy2(source, destination)
    return destination


def write_build_metadata(
    configuration: BuildConfiguration,
    started_utc: str,
    commands: list[list[str]],
    executables: dict[str, Path],
    runtime_file: Path,
    profile: dict[str, Any],
) -> Path:
    """Write the reproducible build record."""
    metadata_path = configuration.build_directory / "build_metadata.json"
    metadata = {
        "schema_version": 2,
        "build_id": configuration.build_id,
        "status": "completed",
        "started_utc": started_utc,
        "finished_utc": utc_now(),
        "repository": {
            "root": str(configuration.repository_root),
            "revision": configuration.revision,
            "dirty": bool(configuration.changes),
            "changes": configuration.changes,
        },
        "profile": profile,
        "commands": commands,
        "environment_script": _file_record(configuration.environment_script),
        "toolchain": toolchain_record(configuration.environment),
        "artifacts": {
            name: _file_record(path) for name, path in executables.items()
        },
        "runtime_files": {runtime_file.name: _file_record(runtime_file)},
    }
    write_json_direct(metadata_path, metadata)
    return metadata_path


def write_generated_settings(
    configuration: BuildConfiguration,
    executables: dict[str, Path],
    metadata_path: Path,
) -> Path:
    """Write settings that select the newly built executables."""
    settings = dict(configuration.settings)
    settings.update(
        {
            "MHDG_REGRESSION_BUILD_ROOT": str(configuration.build_root),
            "MHDG_SERIAL_EXECUTABLE": str(executables["serial"]),
            "MHDG_PARALLEL_EXECUTABLE": str(executables["parallel"]),
            "MHDG_ENVIRONMENT_SCRIPT": str(
                configuration.environment_script.resolve()
            ),
            "MHDG_SOLVER_REVISION": configuration.revision,
            "MHDG_BUILD_DESCRIPTION": (
                f"regression build {configuration.build_id}"
            ),
            "MHDG_BUILD_MANIFEST": str(metadata_path),
        }
    )
    path = configuration.build_directory / "settings.env"
    _write_settings(path, settings)
    return path


def _file_record(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    identity = file_identity(resolved)
    return {
        "path": str(resolved),
        "size": identity["size_bytes"],
        "sha256": identity["sha256"],
    }


def _write_settings(path: Path, settings: dict[str, str]) -> None:
    if any("\n" in value or "\r" in value for value in settings.values()):
        raise BuildError("settings values cannot contain newlines")
    contents = "".join(f"{key}={value}\n" for key, value in sorted(settings.items()))
    path.write_text(contents, encoding="utf-8")
