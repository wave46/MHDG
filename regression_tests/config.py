"""Machine JSON, generated build selection and runtime environment setup."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any


from .documents import load_json
from .support import BundleError
from .files import file_identity, require_file, require_directory

ROOT = Path(__file__).resolve().parent
RUNTIME_FILE = "positionFeketeNodesTri2D.h5"
DEFAULT_MODEL = "NGammaTiTeNeutral"
MODELS = {
    DEFAULT_MODEL: ("N-Gamma-Ti-Te-Neutral", ("rho", "Gamma", "nEi", "nEe", "rhon")),
    "NGammaTiTeNeutralGamma": ("N-Gamma-Ti-Te-NeutralGamma", ("rho", "Gamma", "nEi", "nEe", "rhon", "Gamman")),
}
EXECUTIONS = {"serial": "serial", "mpi": "parall"}


def build_key(model, execution):
    if model not in MODELS or execution not in EXECUTIONS:
        raise BundleError(f"unsupported solver build: {model}/{execution}")
    return f"{model}/{execution}"


def executables_from_manifest(record):
    """Return model/execution executable records from current or existing manifests."""
    if record.get("status") != "completed" or not isinstance(record.get("profile"), dict) or record["profile"].get("dimension") != "2D":
        raise BundleError("build manifest must describe a completed 2D build")
    artifacts = record.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        raise BundleError("build manifest has no executables")
    if record.get("schema_version") == 2:
        model = record["profile"].get("model")
        artifacts = {build_key(model, {"serial": "serial", "parallel": "mpi"}.get(name)): item
                     for name, item in artifacts.items()}
    elif record.get("schema_version") != 3:
        raise BundleError("unsupported build manifest version")
    for key in artifacts:
        parts = key.split("/")
        if len(parts) != 2 or build_key(*parts) != key:
            raise BundleError(f"invalid build key: {key}")
    return artifacts
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


def build_settings(path: Path) -> dict[str, Any]:
    """Select a completed existing build and verify its recorded binaries/data."""
    path = path.resolve()
    record = load_json(path, "build manifest")
    artifacts = executables_from_manifest(record)
    _text(record.get("build_id"), "build_id")
    repository = record.get("repository", {})
    values = {
        "MHDG_BUILD_MANIFEST": str(path),
        "MHDG_SOLVER_REVISION": _text(repository.get("revision"), "repository.revision"),
        "MHDG_EXECUTABLES": {},
    }
    runtime = _artifact(record.get("runtime_files"), RUNTIME_FILE, path.parent)
    for key in artifacts:
        executable = _artifact(artifacts, key, path.parent)
        if not os.access(executable, os.X_OK):
            raise BundleError(f"build executable is not executable: {executable}")
        if (executable.parent / RUNTIME_FILE).resolve() != runtime:
            raise BundleError(f"Fekete data must be beside the executable: {executable}")
        values["MHDG_EXECUTABLES"][key] = str(executable)
    return values


def verify_executable(values, model, execution, executable):
    """Verify the selected artifact again immediately before launching a run."""
    manifest = absolute_setting(values, "MHDG_BUILD_MANIFEST")
    record = load_json(manifest, "build manifest")
    artifacts = executables_from_manifest(record)
    key = build_key(model, execution)
    expected = _artifact(artifacts, key, manifest.parent)
    if expected != executable.resolve():
        raise BundleError(f"selected executable does not match build {key}: {executable}")
    _artifact(record.get("runtime_files"), RUNTIME_FILE, manifest.parent)
    return {"build_id": record["build_id"], "git_commit": record["repository"]["revision"],
            "git_dirty": bool(record["repository"].get("dirty", False))}


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


def runtime_settings(values, requirements):
    """Check the required builds and launcher before creating run directories."""
    for model, execution in sorted(requirements):
        executable = solver_executable(values, execution, model)
        verify_executable(values, model, execution, executable)
    environment = execution_environment(values)
    if any(execution == "mpi" for _, execution in requirements):
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


def absolute_setting(settings: dict[str, Any], key: str) -> Path:
    value = settings.get(key)
    if not value:
        raise BundleError(f"settings must define {key}")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise BundleError(f"{key} must be an absolute path")
    return path.resolve()


def solver_executable(settings, execution, model=DEFAULT_MODEL):
    key = build_key(model, execution)
    selected = settings.get("MHDG_EXECUTABLES", {}).get(key)
    if not selected:
        raise BundleError(f"no build selected for {key}; use check --build or select a compatible --build-manifest")
    path = Path(selected)
    if not path.is_absolute() or not path.is_file() or not os.access(path, os.X_OK):
        raise BundleError(f"{key} is not an executable file: {path}")
    return path.resolve()


def runtime_files(executable: Path) -> dict[str, Path]:
    path = executable.parent / RUNTIME_FILE
    return {RUNTIME_FILE: require_file(path, "solver runtime")}


def selected_mpi_launcher(settings: dict[str, Any]) -> Path:
    value = settings.get("MHDG_MPI_LAUNCHER")
    if not value:
        raise BundleError("settings must define MHDG_MPI_LAUNCHER")
    resolved = shutil.which(value)
    if not resolved:
        raise BundleError(f"MPI launcher is not executable or not found: {value}")
    return Path(resolved).resolve()
