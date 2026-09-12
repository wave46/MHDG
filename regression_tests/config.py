"""Machine JSON, generated build selection and runtime environment setup."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


from .documents import load_json
from .support import BundleError
from .files import file_identity, require_file, require_directory

ROOT = Path(__file__).resolve().parent
RUNTIME_FILE = "positionFeketeNodesTri2D.h5"
MACHINE_KEYS = {
    "run_root": "MHDG_REGRESSION_RUN_ROOT",
    "build_root": "MHDG_REGRESSION_BUILD_ROOT",
    "environment_script": "MHDG_ENVIRONMENT_SCRIPT",
    "mpi_launcher": "MHDG_MPI_LAUNCHER",
    "build_jobs": "MHDG_REGRESSION_BUILD_JOBS",
}


def settings(path=None, *, case=None, bundle=None, build_manifest=None, use_build=True):
    """Resolve explicit selections before local defaults; never choose newest data.

    Return one resolved runtime dictionary for preparation/execution; no other
    module reads machine settings or writes intermediate settings files.
    """
    values, defaults = machine_settings(path)
    if bundle is None:
        bundle = defaults.get("bundles", {}).get(case)
    if build_manifest is None and use_build:
        build_manifest = defaults.get("build")
    values.setdefault("MHDG_REGRESSION_RUN_ROOT", str(Path.home() / ".cache/mhdg-regression"))
    values.setdefault("MHDG_REGRESSION_BUILD_ROOT", str(Path(values["MHDG_REGRESSION_RUN_ROOT"]) / "builds"))
    if bundle is not None:
        values["MHDG_REGRESSION_DATA_ROOT"] = str(Path(bundle).expanduser().resolve())
    if build_manifest is not None and use_build:
        values.update(build_settings(Path(build_manifest).expanduser()))
    return values


def machine_settings(path=None):
    """Read machine preferences/default paths without loading a build or bundle."""
    if path is None:
        configured = os.environ.get("MHDG_REGRESSION_SETTINGS")
        path = Path(configured).expanduser() if configured else ROOT / "settings.local.json"
        optional = not configured
    else:
        path, optional = Path(path).expanduser(), False
    document = {} if optional and not path.exists() else load_json(path, "machine settings")
    _keys(document, {*MACHINE_KEYS, "defaults"}, "machine settings")
    values = {}
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
    defaults = {"bundles": {name: _path(value, path.parent, f"defaults.bundles.{name}")
                             for name, value in bundles.items()},
                **({"build": _path(defaults["build"], path.parent, "defaults.build")}
                   if "build" in defaults else {})}
    return values, defaults


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
    _text(record.get("build_id"), "build_id")
    repository = record.get("repository")
    if not isinstance(repository, dict):
        raise BundleError("build manifest is missing repository provenance")
    values = {
        "MHDG_BUILD_MANIFEST": str(path),
        "MHDG_SOLVER_REVISION": _text(repository.get("revision"), "repository.revision"),
    }
    runtime = _artifact(record.get("runtime_files"), RUNTIME_FILE, path.parent)
    artifacts = record.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts or artifacts.keys() - {"serial", "parallel"}:
        raise BundleError("build manifest must contain serial and/or parallel executables")
    if "variants" in profile and profile["variants"] != {
        name: "serial" if name == "serial" else "parall" for name in artifacts
    }:
        raise BundleError("build artifacts do not match the declared variants")
    for variant, key in (("serial", "MHDG_SERIAL_EXECUTABLE"), ("parallel", "MHDG_PARALLEL_EXECUTABLE")):
        if variant not in artifacts:
            continue
        executable = _artifact(artifacts, variant, path.parent)
        if not os.access(executable, os.X_OK):
            raise BundleError(f"build executable is not executable: {executable}")
        if (executable.parent / RUNTIME_FILE).resolve() != runtime:
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
    for execution in {layout["execution"] for layout in layouts}:
        key = "MHDG_SERIAL_EXECUTABLE" if execution == "serial" else "MHDG_PARALLEL_EXECUTABLE"
        if key not in values:
            raise BundleError(f"no build selected for {execution}; pass --build-manifest FILE or use check --build")
        runtime_files(solver_executable(values, execution))
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





def source_environment(
    script: Path,
    base_environment: dict[str, str] | None = None,
) -> tuple[Path, dict[str, str]]:
    """Source one shell script and return its resolved path and environment."""
    script = require_file(script, "environment script")
    command = [
        "bash",
        "-c",
        'source "$1" >/dev/null && env -0',
        "mhdg-regression-environment",
        str(script),
    ]
    try:
        completed = subprocess.run(
            command,
            env=base_environment,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        raise BundleError(f"cannot source environment script {script}: {exc}") from exc
    if completed.returncode != 0:
        error = completed.stderr.decode(errors="replace").strip()
        detail = f": {error}" if error else ""
        raise BundleError(
            f"environment script failed ({completed.returncode}){detail}"
        )
    try:
        return script, {
            os.fsdecode(key): os.fsdecode(value)
            for entry in completed.stdout.split(b"\0")
            if entry
            for key, value in [entry.split(b"=", 1)]
        }
    except ValueError as exc:
        raise BundleError("environment script produced invalid environment data") from exc


def bundle_root_from_settings(values):
    """Resolve the selected external bundle; file contracts are owned by files.py."""
    value = values.get("MHDG_REGRESSION_DATA_ROOT")
    if not value:
        raise BundleError("no bundle selected; pass --bundle DIR or set defaults.bundles")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise BundleError("selected bundle root must be an absolute path")
    return require_directory(path, "bundle root")


def absolute_setting(settings: dict[str, str], key: str) -> Path:
    value = settings.get(key)
    if not value:
        raise BundleError(f"settings must define {key}")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise BundleError(f"{key} must be an absolute path")
    return path.resolve()


def solver_executable(settings: dict[str, str], execution: str) -> Path:
    key = (
        "MHDG_SERIAL_EXECUTABLE"
        if execution == "serial"
        else "MHDG_PARALLEL_EXECUTABLE"
    )
    path = absolute_setting(settings, key)
    if not path.is_file() or not os.access(path, os.X_OK):
        raise BundleError(f"{key} is not an executable file: {path}")
    return path


def runtime_files(executable: Path) -> dict[str, Path]:
    path = executable.parent / RUNTIME_FILE
    return {RUNTIME_FILE: require_file(path, "solver runtime")}


def selected_mpi_launcher(settings: dict[str, str]) -> Path:
    value = settings.get("MHDG_MPI_LAUNCHER")
    if not value:
        raise BundleError("settings must define MHDG_MPI_LAUNCHER")
    resolved = shutil.which(value)
    if not resolved:
        raise BundleError(f"MPI launcher is not executable or not found: {value}")
    return Path(resolved).resolve()
