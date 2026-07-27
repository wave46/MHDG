"""Record immutable files shared by every cell in one suite run."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from support.errors import BundleError
from support.files import file_identity
from support.paths import require_file


def suite_execution_inputs(
    settings_path: Path,
    settings: dict[str, str],
    bundle_root: Path,
    layout_ids: list[str],
    tracked_files: dict[str, Path],
) -> dict[str, Any]:
    """Describe files that must remain stable across suite resumes."""
    records = {
        "settings": _file_record(settings_path, "settings"),
        "bundle_manifest": _file_record(
            bundle_root / "manifest.json",
            "bundle manifest",
        ),
        **{
            name: _file_record(path, name.replace("_", " "))
            for name, path in tracked_files.items()
        },
    }
    if any(layout.startswith("serial_") for layout in layout_ids):
        records["serial_executable"] = _setting_file_record(
            settings,
            "MHDG_SERIAL_EXECUTABLE",
        )
    if any(layout.startswith("mpi") for layout in layout_ids):
        records["parallel_executable"] = _setting_file_record(
            settings,
            "MHDG_PARALLEL_EXECUTABLE",
        )
        launcher = settings.get("MHDG_MPI_LAUNCHER")
        resolved_launcher = shutil.which(launcher) if launcher else None
        if resolved_launcher is None:
            raise BundleError("MHDG_MPI_LAUNCHER is not executable or not found")
        records["mpi_launcher"] = _file_record(
            Path(resolved_launcher),
            "MPI launcher",
        )
    for key, name in (
        ("MHDG_ENVIRONMENT_SCRIPT", "environment_script"),
        ("MHDG_BUILD_MANIFEST", "build_manifest"),
    ):
        if settings.get(key):
            records[name] = _setting_file_record(settings, key)
    return records


def _setting_file_record(
    settings: dict[str, str],
    key: str,
) -> dict[str, Any]:
    value = settings.get(key)
    if not value:
        raise BundleError(f"settings must define {key}")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise BundleError(f"{key} must be an absolute path")
    return _file_record(path, key)


def _file_record(path: Path, label: str) -> dict[str, Any]:
    path = require_file(path, label)
    return {"path": str(path), **file_identity(path)}
