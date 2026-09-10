"""Machine preferences and invocation selections, adapted to existing runners."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from bundle.settings import read_settings
from support.documents import load_json
from support.environments import source_environment
from support.errors import BundleError
from support.files import file_identity

ROOT = Path(__file__).resolve().parent
MACHINE_KEYS = {
    "run_root": "MHDG_REGRESSION_RUN_ROOT",
    "build_root": "MHDG_REGRESSION_BUILD_ROOT",
    "environment_script": "MHDG_ENVIRONMENT_SCRIPT",
    "mpi_launcher": "MHDG_MPI_LAUNCHER",
    "build_jobs": "MHDG_REGRESSION_BUILD_JOBS",
}


def settings(path=None, *, case=None, bundle=None, build_manifest=None, use_build=True):
    """Resolve explicit selections before local defaults; never choose newest data.

    The MHDG_* dictionary is a temporary bridge to the existing runners. No
    intermediate settings files are written; resume records the selected values.
    """
    if path is None:
        configured = os.environ.get("MHDG_REGRESSION_SETTINGS")
        path = Path(configured).expanduser() if configured else ROOT / "settings.local.json"
        optional = not configured
    else:
        path, optional = Path(path).expanduser(), False
    defaults = {}
    if path.suffix == ".env":
        # Explicit migration input for old campaigns/prebuilt environments.
        values = read_settings(path)
    else:
        document = {} if optional and not path.exists() else load_json(path, "machine settings")
        _keys(document, {*MACHINE_KEYS, "defaults"}, "machine settings")
        values = {"MHDG_REGRESSION_SETTINGS_VERSION": "2"}
        for key, value in document.items():
            if key == "defaults":
                continue
            if key == "build_jobs":
                if type(value) is not int or value < 1:
                    raise BundleError("build_jobs must be a positive integer")
                values[MACHINE_KEYS[key]] = str(value)
            elif key == "mpi_launcher" and isinstance(value, str) and "/" not in value:
                values[MACHINE_KEYS[key]] = _text(value, key)
            else:
                values[MACHINE_KEYS[key]] = str(_path(value, path.parent, key))
        defaults = document.get("defaults", {})
        _keys(defaults, {"build", "bundles"}, "defaults")
        bundles = defaults.get("bundles", {})
        if not isinstance(bundles, dict):
            raise BundleError("defaults.bundles must map case identifiers to bundle paths")
        if bundle is None and case in bundles:
            bundle = _path(bundles[case], path.parent, f"defaults.bundles.{case}")
        if build_manifest is None and use_build and "build" in defaults:
            build_manifest = _path(defaults["build"], path.parent, "defaults.build")
    values.setdefault("MHDG_REGRESSION_RUN_ROOT", str(Path.home() / ".cache/mhdg-regression"))
    values.setdefault("MHDG_REGRESSION_BUILD_ROOT", str(Path(values["MHDG_REGRESSION_RUN_ROOT"]) / "builds"))
    if bundle is not None:
        values["MHDG_REGRESSION_DATA_ROOT"] = str(Path(bundle).expanduser().resolve())
    if build_manifest is not None and use_build:
        values.update(build_settings(Path(build_manifest).expanduser()))
    return values


def build_settings(path: Path) -> dict[str, str]:
    """Select a completed existing build and verify its recorded binaries/data."""
    path = path.resolve()
    record = load_json(path, "build manifest")
    if record.get("schema_version") != 2 or record.get("status") != "completed":
        raise BundleError("build manifest must describe a completed version-2 build")
    profile = record.get("profile", {})
    if not isinstance(profile, dict) or (
        profile.get("model"), profile.get("dimension")
    ) != ("NGammaTiTeNeutral", "2D"):
        raise BundleError("this harness currently supports NGammaTiTeNeutral 2D builds")
    build_id = _text(record.get("build_id"), "build_id")
    repository = record.get("repository")
    if not isinstance(repository, dict):
        raise BundleError("build manifest is missing repository provenance")
    values = {
        "MHDG_BUILD_MANIFEST": str(path),
        "MHDG_BUILD_DESCRIPTION": f"regression build {build_id}",
        "MHDG_SOLVER_REVISION": _text(repository.get("revision"), "repository.revision"),
    }
    runtime = _artifact(record.get("runtime_files"), "positionFeketeNodesTri2D.h5", path.parent)
    for variant, key in (("serial", "MHDG_SERIAL_EXECUTABLE"), ("parallel", "MHDG_PARALLEL_EXECUTABLE")):
        executable = _artifact(record.get("artifacts"), variant, path.parent)
        if not os.access(executable, os.X_OK):
            raise BundleError(f"build executable is not executable: {executable}")
        if (executable.parent / "positionFeketeNodesTri2D.h5").resolve() != runtime:
            raise BundleError(f"Fekete data must be beside the executable: {executable}")
        values[key] = str(executable)
    return values


def _artifact(records, name, directory):
    item = records.get(name) if isinstance(records, dict) else None
    if not isinstance(item, dict):
        raise BundleError(f"build manifest is missing artifact {name}")
    path = _path(item.get("path"), directory, f"artifact {name}.path")
    if not path.is_file():
        raise BundleError(f"build artifact is missing: {path}")
    actual = file_identity(path)
    if any(item.get(key) != actual[key] for key in ("size_bytes", "sha256")):
        raise BundleError(f"build artifact checksum/size changed: {path}")
    return path


def execution_environment(values):
    script = values.get("MHDG_ENVIRONMENT_SCRIPT")
    return source_environment(Path(script))[1] if script else dict(os.environ)


def mpi_launcher(values, environment):
    """Prefer Open MPI's explicit executable name; reject another MPI family."""
    configured = values.get("MHDG_MPI_LAUNCHER")
    choices = [configured] if configured else ["mpirun.openmpi", "mpirun"]
    for choice in choices:
        found = shutil.which(choice, path=environment.get("PATH", ""))
        if not found:
            continue
        if not openmpi_version(found, environment):
            continue
        return str(Path(found).resolve())
    raise BundleError("usable Open MPI launcher not found; load its environment or set mpi_launcher")


def openmpi_version(launcher, environment):
    try:
        result = subprocess.run(
            [launcher, "--version"], env=environment, capture_output=True,
            text=True, timeout=5, check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    output = result.stdout + result.stderr
    return output.strip() if result.returncode == 0 and any(
        name in output for name in ("Open MPI", "OpenRTE")
    ) else None


def runtime_settings(values, layouts):
    """Check selected execution prerequisites before creating run directories."""
    from preparation.configuration import _runtime_files, _solver_executable

    for execution in {layout["execution"] for layout in layouts}:
        key = "MHDG_SERIAL_EXECUTABLE" if execution == "serial" else "MHDG_PARALLEL_EXECUTABLE"
        if key not in values:
            raise BundleError("no build selected; pass --build-manifest FILE or use check --build")
        _runtime_files(_solver_executable(values, execution))
    environment = execution_environment(values)
    if any(layout["execution"] == "mpi" for layout in layouts):
        values["MHDG_MPI_LAUNCHER"] = mpi_launcher(values, environment)
    return values


def _keys(document, allowed, label):
    if not isinstance(document, dict):
        raise BundleError(f"{label} must be a JSON object")
    unknown = document.keys() - allowed
    if unknown:
        raise BundleError(f"unknown {label} fields: {', '.join(sorted(unknown))}")


def _text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise BundleError(f"{label} must be a nonempty string")
    return value


def _path(value, directory, label):
    path = Path(_text(value, label)).expanduser()
    return (directory / path).resolve() if not path.is_absolute() else path.resolve()
