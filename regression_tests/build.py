"""Build the required 2-D neutral solver variants and record their provenance."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from support.documents import write_json_atomic
from support.environments import source_environment
from support.errors import BundleError
from support.files import file_identity
from support.time import utc_now

VARIANTS = {"serial": "serial", "parallel": "parall"}
RUNTIME_FILE = "positionFeketeNodesTri2D.h5"


@dataclass(frozen=True)
class BuildResult:
    path: Path
    metadata_path: Path
    executables: dict[str, Path]


def build_solver(
    settings: dict[str, str], repository_root: Path, jobs: int | None = None,
    *, variants=("serial", "parallel"),
) -> BuildResult:
    """Build each selected variant once; return artifacts without writing settings."""
    selected = set(variants)
    if not selected or selected - VARIANTS.keys():
        raise BundleError("build variants must select serial and/or parallel")
    try:
        jobs = int(settings.get("MHDG_REGRESSION_BUILD_JOBS", "8")) if jobs is None else jobs
    except ValueError as exc:
        raise BundleError("build jobs must be a positive integer") from exc
    if type(jobs) is not int or jobs < 1:
        raise BundleError("build jobs must be a positive integer")
    repository_root = repository_root.resolve()
    library = repository_root / "lib"
    if not library.is_dir():
        raise BundleError(f"solver library directory does not exist: {library}")
    runtime = repository_root / "test" / RUNTIME_FILE
    if not runtime.is_file():
        raise BundleError(f"generic Fekete-node data does not exist: {runtime}")
    root = settings.get("MHDG_REGRESSION_BUILD_ROOT")
    if root is None:
        run_root = settings.get("MHDG_REGRESSION_RUN_ROOT")
        if not run_root:
            raise BundleError("settings must define a build root or run root")
        root = str(Path(run_root) / "builds")
    build_root = Path(root).expanduser()
    if not build_root.is_absolute():
        raise BundleError("regression build root must be an absolute path")

    script = None
    environment = dict(os.environ)
    if configured := settings.get("MHDG_ENVIRONMENT_SCRIPT"):
        script, environment = source_environment(Path(configured).expanduser())
    revision = _capture(["git", "rev-parse", "HEAD"], repository_root).strip()
    changes = _capture(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"], repository_root,
    ).splitlines()
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    build_id = f"{timestamp}-{revision[:8]}" + ("-dirty" if changes else "")
    directory = build_root.resolve() / build_id
    try:
        (directory / "bin").mkdir(parents=True)
        (directory / "logs").mkdir()
    except OSError as exc:
        raise BundleError(f"cannot create build directory {directory}: {exc}") from exc

    started = utc_now()
    commands = []
    executables = {}
    # Existing objects/generated sources may come from a manual build with other
    # flags. An empty tree needs no clean; every serial/MPI transition does.
    clean_needed = any(
        any(library.glob(pattern))
        for pattern in ("*.o", "*.mod", "*.smod", "*.kmo", "*.F", "*.f90", "*.F90", "MHDG-*")
    )
    for variant, mode in VARIANTS.items():
        if variant not in selected:
            continue
        target = f"MHDG-NGammaTiTeNeutral-{mode}-2D"
        compile_command = [
            "make", f"-j{jobs}", f"MODE={mode}", "COMPTYPE=opt",
            "MDL=NGammaTiTeNeutral", "DIM=2D",
            f"MHDG_GIT_COMMIT={revision}", f"MHDG_GIT_DIRTY={'true' if changes else 'false'}",
            f"MHDG_BUILD_ID={build_id}", target,
        ]
        steps = [["make", "clean"], compile_command] if clean_needed else [compile_command]
        commands.extend(steps)
        _run_logged(steps, library, environment, directory / "logs" / f"{variant}.log")
        source = library / target
        if not source.is_file() or not os.access(source, os.X_OK):
            raise BundleError(f"built {variant} executable was not produced: {source}")
        destination = directory / "bin" / target
        shutil.copy2(source, destination)
        executables[variant] = destination
        clean_needed = True

    runtime_copy = directory / "bin" / RUNTIME_FILE
    shutil.copy2(runtime, runtime_copy)
    metadata_path = directory / "build_metadata.json"
    metadata = {
        "schema_version": 2,
        "build_id": build_id,
        "status": "completed",
        "started_utc": started,
        "finished_utc": utc_now(),
        "repository": {
            "root": str(repository_root), "revision": revision,
            "dirty": bool(changes), "changes": changes,
        },
        "profile": {
            "model": "NGammaTiTeNeutral", "dimension": "2D", "compile_type": "opt",
            "jobs": jobs, "variants": {name: VARIANTS[name] for name in executables},
        },
        "commands": commands,
        "environment_script": _file_record(script) if script else None,
        "toolchain": {
            "make": _version(["make", "--version"], environment),
            "fortran_compiler": _version(["mpifort", "--version"], environment),
            "compiler_command": _version(["mpifort", "-show"], environment),
        },
        "artifacts": {name: _file_record(path) for name, path in executables.items()},
        "runtime_files": {RUNTIME_FILE: _file_record(runtime_copy)},
    }
    write_json_atomic(metadata_path, metadata, "build metadata")
    return BuildResult(directory, metadata_path, executables)


def _file_record(path):
    return {"path": str(path.resolve()), **file_identity(path)}


def _run_logged(commands, directory, environment, log_path):
    environment = {**environment, "PWD": str(directory)}
    try:
        with log_path.open("w", encoding="utf-8") as log:
            for command in commands:
                rendered = shlex.join(command)
                print(f"build: {rendered}", flush=True)
                print(f"$ {rendered}", file=log, flush=True)
                with subprocess.Popen(
                    command, cwd=directory, env=environment, stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT, text=True, bufsize=1,
                ) as process:
                    for line in process.stdout:
                        print(line, end="")
                        log.write(line)
                    if process.wait() != 0:
                        raise BundleError(f"command failed; see {log_path}: {rendered}")
    except OSError as exc:
        raise BundleError(f"cannot run build; see {log_path}: {exc}") from exc


def _capture(command, directory, environment=None):
    try:
        return subprocess.run(
            command, cwd=directory, env=environment, check=True, capture_output=True, text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        raise BundleError(f"cannot run {shlex.join(command)}") from exc


def _version(command, environment):
    try:
        lines = _capture(command, Path.cwd(), environment).splitlines()
        return lines[0] if lines else None
    except BundleError:
        return None
