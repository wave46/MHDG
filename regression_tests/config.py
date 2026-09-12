"""Machine preferences and invocation selections, adapted to existing runners."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition
from bundle.schemas import load_validated_json
from catalogs.layouts import layout_pairs, load_layouts

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
    values, defaults = machine_settings(path)
    if bundle is None:
        bundle = defaults.get("bundles", {}).get(case)
    if build_manifest is None and use_build:
        build_manifest = defaults.get("build")
    values.setdefault("MHDG_REGRESSION_RUN_ROOT", str(Path.home() / ".cache/mhdg-regression"))
    values.setdefault("MHDG_REGRESSION_BUILD_ROOT", str(Path(values["MHDG_REGRESSION_RUN_ROOT"]) / "builds"))
    if bundle is not None:
        values["MHDG_REGRESSION_DATA_ROOT"] = str(Path(bundle).expanduser().resolve())
    if not use_build or build_manifest is not None:
        # A partial new build must not inherit an executable from an old selection.
        for key in ("MHDG_SERIAL_EXECUTABLE", "MHDG_PARALLEL_EXECUTABLE", "MHDG_SOLVER_REVISION",
                    "MHDG_BUILD_MANIFEST"):
            values.pop(key, None)
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
    runtime = _artifact(record.get("runtime_files"), "positionFeketeNodesTri2D.h5", path.parent)
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
    from .prepare import _runtime_files, _solver_executable

    for execution in {layout["execution"] for layout in layouts}:
        key = "MHDG_SERIAL_EXECUTABLE" if execution == "serial" else "MHDG_PARALLEL_EXECUTABLE"
        if key not in values:
            raise BundleError(f"no build selected for {execution}; pass --build-manifest FILE or use check --build")
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


def load_suite_definition(
    suite_id: str,
    suites_path: Path,
    layouts_path: Path,
    case_directory: Path,
    *, case_id: str | None = None,
) -> dict[str, Any]:
    """Load one suite and validate its workflow and layout references."""
    schema_path = suites_path.parent / "schemas" / "suites.schema.json"
    document = load_validated_json(suites_path, schema_path, "suite definitions")
    declaration = document["suites"].get(suite_id)
    if declaration is None:
        available = ", ".join(sorted(document["suites"]))
        raise BundleError(f"unknown suite {suite_id}; available: {available}")
    layouts = load_layouts(layouts_path)
    defaults = document["defaults"]
    relations = declaration.get("relations", [])
    selected = declaration.get(
        "layouts", "all" if relations else [defaults["layout"]]
    )
    selected = list(layouts) if selected == "all" else selected
    unknown_layouts = [layout for layout in selected if layout not in layouts]
    if unknown_layouts:
        raise BundleError(
            f"suite {suite_id} has unknown layouts: {', '.join(unknown_layouts)}"
        )

    pairs = layout_pairs({name: layouts[name] for name in selected}, relations)
    suite = {
        "description": declaration["description"],
        "case_id": case_id or declaration.get("case", defaults["case"]),
        "diagnostics": declaration.get("diagnostics", "off"),
        "workflow_ids": declaration["workflows"],
        "layouts": list(dict.fromkeys(
            layout for pair in pairs for layout in pair.values()
        )) if pairs else selected,
        "reference_comparisons": declaration.get("reference_comparisons", not pairs),
    }
    if pairs:
        suite["layout_comparisons"] = pairs
        suite["tolerance_profile"] = declaration["tolerance_profile"]
        if "layout_comparison_policy" in declaration:
            suite["layout_comparison_policy"] = declaration["layout_comparison_policy"]

    case = load_case_definition(suite["case_id"], case_directory)
    unknown_workflows = [
        workflow_id
        for workflow_id in suite["workflow_ids"]
        if workflow_id not in case["workflows"]
    ]
    if unknown_workflows:
        raise BundleError(
            f"suite {suite_id} has unknown workflows: "
            f"{', '.join(unknown_workflows)}"
        )
    return suite


def load_selection(name, suites_path, layouts_path, case_directory, *, case_id=None):
    """Expand a profile into unique case/suite selections; focused suites use the same loader."""
    document = load_validated_json(
        suites_path, suites_path.parent / "schemas/suites.schema.json", "suite definitions",
    )
    profiles = document.get("profiles", {})
    if profiles.keys() & document["suites"].keys():
        raise BundleError("profile and suite names must be distinct")
    is_profile = name in profiles
    if is_profile and case_id:
        raise BundleError("--case applies to a focused suite; profiles declare their cases")

    def expand(profile, ancestors=()):
        if profile in ancestors:
            raise BundleError(f"cyclic profile inclusion: {' -> '.join((*ancestors, profile))}")
        if profile not in profiles:
            raise BundleError(f"unknown included profile: {profile}")
        for parent in profiles[profile].get("include", []):
            yield from expand(parent, (*ancestors, profile))
        yield from profiles[profile]["checks"]

    entries = expand(name) if is_profile else [{"suite": name, "case": case_id}]
    selected = {}
    for entry in entries:
        suite = load_suite_definition(
            entry["suite"], suites_path, layouts_path, case_directory, case_id=entry.get("case"),
        )
        selected[(entry["suite"], suite["case_id"])] = {"suite_id": entry["suite"], **suite}
    if not selected:
        raise BundleError(f"selection {name} contains no checks")
    return is_profile, list(selected.values())


def require_bundle_class(
    bundle_root: Path,
    case_directory: Path,
    required: str,
) -> None:
    """Require the source bundle to declare the requested publication class."""
    schema = case_directory.parent / "schemas" / "bundle-manifest.schema.json"
    manifest = load_validated_json(
        bundle_root / "manifest.json",
        schema,
        "bundle manifest",
    )
    actual = manifest.get("bundle_class", "unspecified")
    if actual != required:
        raise BundleError(
            f"golden-check requires bundle_class={required}; found {actual}"
        )
