#!/usr/bin/env python3
"""Build the serial and parallel regression executables reproducibly."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO

from check_bundle import BundleError, read_settings
from support.documents import write_json_direct
from support.files import sha256_digest
from support.time import utc_now


MODEL = "NGammaTiTeNeutral"
DIMENSION = "2D"
COMPILE_TYPE = "opt"
DEFAULT_JOBS = 8
VARIANTS = {"serial": "serial", "parallel": "parall"}
RUNTIME_FILE = "positionFeketeNodesTri2D.h5"


class BuildError(BundleError):
    """Raised when the solver build cannot be completed or recorded."""


@dataclass(frozen=True)
class BuildResult:
    path: Path
    settings_path: Path
    metadata_path: Path
    executables: dict[str, Path]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--jobs", type=parse_build_jobs)
    parser.add_argument(
        "--repository-root", type=Path, default=Path(__file__).resolve().parents[2],
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    try:
        result = build_solver(args.settings, args.repository_root, args.jobs)
    except BundleError as exc:
        print(f"build failed: {exc}", file=sys.stderr)
        return 1

    print(f"build completed: {result.path}")
    print(f"generated settings: {result.settings_path}")
    return 0


def build_solver(
    settings_path: Path,
    repository_root: Path,
    jobs: int | None = None,
) -> BuildResult:
    """Build clean serial and parallel variants and return their local settings."""
    settings = read_settings(settings_path)
    repository_root = repository_root.resolve()
    lib_dir = repository_root / "lib"
    test_dir = repository_root / "test"
    _require_directory(lib_dir, "solver library directory")

    environment_script = Path(
        settings.get(
            "MHDG_ENVIRONMENT_SCRIPT",
            str(lib_dir / "Make.inc" / "init_vars_libs.sh"),
        )
    ).expanduser()
    _require_file(environment_script, "build environment script")
    environment = _environment_from_script(environment_script)

    jobs = jobs or _configured_jobs(settings)
    revision, changes = _git_state(repository_root)
    build_root = _build_root(settings)
    build_id = _build_id(revision, bool(changes))
    build_dir = build_root / build_id
    bin_dir = build_dir / "bin"
    log_dir = build_dir / "logs"
    try:
        bin_dir.mkdir(parents=True)
        log_dir.mkdir()
    except OSError as exc:
        raise BuildError(f"cannot create build directory {build_dir}: {exc}") from exc

    started = utc_now()
    commands: list[list[str]] = []
    executables: dict[str, Path] = {}
    for variant, mode in VARIANTS.items():
        target = f"MHDG-{MODEL}-{mode}-{DIMENSION}"
        clean = ["make", "clean"]
        compile_command = [
            "make",
            f"-j{jobs}",
            f"MODE={mode}",
            f"COMPTYPE={COMPILE_TYPE}",
            target,
        ]
        commands.extend((clean, compile_command))
        log_path = log_dir / f"{variant}.log"
        _run_logged((clean, compile_command), lib_dir, environment, log_path)

        source = lib_dir / target
        _require_executable(source, f"built {variant} executable")
        destination = bin_dir / target
        shutil.copy2(source, destination)
        executables[variant] = destination

    runtime_source = test_dir / RUNTIME_FILE
    _require_file(runtime_source, "generic Fekete-node data")
    runtime_destination = bin_dir / RUNTIME_FILE
    shutil.copy2(runtime_source, runtime_destination)

    metadata_path = build_dir / "build_metadata.json"
    metadata = {
        "schema_version": 1,
        "build_id": build_id,
        "status": "completed",
        "started_utc": started,
        "finished_utc": utc_now(),
        "repository": {
            "root": str(repository_root),
            "revision": revision,
            "dirty": bool(changes),
            "changes": changes,
        },
        "profile": {
            "model": MODEL,
            "dimension": DIMENSION,
            "compile_type": COMPILE_TYPE,
            "jobs": jobs,
            "variants": VARIANTS,
        },
        "commands": commands,
        "environment_script": _file_record(environment_script),
        "toolchain": _toolchain_record(environment),
        "artifacts": {
            name: _file_record(path) for name, path in executables.items()
        },
        "runtime_files": {RUNTIME_FILE: _file_record(runtime_destination)},
    }
    write_json_direct(metadata_path, metadata)

    generated_settings = dict(settings)
    generated_settings.update(
        {
            "MHDG_REGRESSION_BUILD_ROOT": str(build_root),
            "MHDG_SERIAL_EXECUTABLE": str(executables["serial"]),
            "MHDG_PARALLEL_EXECUTABLE": str(executables["parallel"]),
            "MHDG_ENVIRONMENT_SCRIPT": str(environment_script.resolve()),
            "MHDG_SOLVER_REVISION": revision,
            "MHDG_BUILD_DESCRIPTION": f"regression build {build_id}",
            "MHDG_BUILD_MANIFEST": str(metadata_path),
        }
    )
    generated_settings_path = build_dir / "settings.env"
    _write_settings(generated_settings_path, generated_settings)
    return BuildResult(
        build_dir, generated_settings_path, metadata_path, executables
    )


def _configured_jobs(settings: dict[str, str]) -> int:
    value = settings.get("MHDG_REGRESSION_BUILD_JOBS")
    if value is None:
        return DEFAULT_JOBS
    try:
        return parse_build_jobs(value)
    except argparse.ArgumentTypeError as exc:
        raise BuildError(str(exc)) from exc


def parse_build_jobs(value: str) -> int:
    """Parse a positive build-job count for CLI and settings inputs."""
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "build jobs must be a positive integer"
        ) from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("build jobs must be a positive integer")
    return parsed


def _build_root(settings: dict[str, str]) -> Path:
    configured = settings.get("MHDG_REGRESSION_BUILD_ROOT")
    if configured:
        root = Path(configured).expanduser()
    else:
        run_root = settings.get("MHDG_REGRESSION_RUN_ROOT")
        if not run_root:
            raise BuildError(
                "settings must define MHDG_REGRESSION_BUILD_ROOT or "
                "MHDG_REGRESSION_RUN_ROOT"
            )
        root = Path(run_root).expanduser() / "builds"
    if not root.is_absolute():
        raise BuildError("regression build root must be an absolute path")
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise BuildError(f"cannot create regression build root {root}: {exc}") from exc
    return root.resolve()


def _environment_from_script(script: Path) -> dict[str, str]:
    command = [
        "bash",
        "-c",
        'source "$1" >/dev/null && env -0',
        "mhdg-build-environment",
        str(script.resolve()),
    ]
    try:
        completed = subprocess.run(command, check=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", b"").decode(errors="replace").strip()
        suffix = f": {detail}" if detail else ""
        raise BuildError(f"cannot load build environment {script}{suffix}") from exc
    return {
        key.decode(errors="surrogateescape"): value.decode(errors="surrogateescape")
        for entry in completed.stdout.split(b"\0")
        if entry
        for key, value in [entry.split(b"=", 1)]
    }


def _run_logged(
    commands: tuple[list[str], ...],
    cwd: Path,
    environment: dict[str, str],
    log_path: Path,
) -> None:
    environment = {**environment, "PWD": str(cwd)}
    try:
        with log_path.open("w", encoding="utf-8") as log:
            for command in commands:
                rendered = shlex.join(command)
                print(f"build: {rendered}", flush=True)
                print(f"$ {rendered}", file=log, flush=True)
                _stream_command(command, cwd, environment, log)
    except OSError as exc:
        raise BuildError(f"cannot run build; see {log_path}: {exc}") from exc


def _stream_command(
    command: list[str], cwd: Path, environment: dict[str, str], log: TextIO
) -> None:
    process = subprocess.Popen(
        command,
        cwd=cwd,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert process.stdout is not None
    with process.stdout:
        for line in process.stdout:
            print(line, end="")
            log.write(line)
    if process.wait() != 0:
        raise BuildError(f"command failed; see {log.name}: {shlex.join(command)}")


def _git_state(repository_root: Path) -> tuple[str, list[str]]:
    revision = _capture(["git", "rev-parse", "HEAD"], repository_root).strip()
    status = _capture(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        repository_root,
    )
    return revision, [line for line in status.splitlines() if line]


def _capture(command: list[str], cwd: Path, environment: dict[str, str] | None = None) -> str:
    try:
        completed = subprocess.run(
            command,
            cwd=cwd,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise BuildError(f"cannot run {shlex.join(command)}") from exc
    return completed.stdout


def _toolchain_record(environment: dict[str, str]) -> dict[str, str | None]:
    return {
        "make": _optional_version(["make", "--version"], environment),
        "fortran_compiler": _optional_version(["mpifort", "--version"], environment),
        "compiler_command": _optional_version(["mpifort", "-show"], environment),
    }


def _optional_version(command: list[str], environment: dict[str, str]) -> str | None:
    try:
        output = _capture(command, Path.cwd(), environment)
    except BuildError:
        return None
    return output.splitlines()[0] if output else None


def _file_record(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved),
        "size": resolved.stat().st_size,
        "sha256": sha256_digest(resolved),
    }


def _write_settings(path: Path, settings: dict[str, str]) -> None:
    if any("\n" in value or "\r" in value for value in settings.values()):
        raise BuildError("settings values cannot contain newlines")
    contents = "".join(f"{key}={value}\n" for key, value in sorted(settings.items()))
    path.write_text(contents, encoding="utf-8")


def _build_id(revision: str, dirty: bool) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    suffix = "-dirty" if dirty else ""
    return f"{timestamp}-{revision[:8]}{suffix}"


def _require_directory(path: Path, label: str) -> None:
    if not path.is_dir():
        raise BuildError(f"{label} does not exist: {path}")


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise BuildError(f"{label} does not exist: {path}")


def _require_executable(path: Path, label: str) -> None:
    if not path.is_file() or not os.access(path, os.X_OK):
        raise BuildError(f"{label} was not produced: {path}")


if __name__ == "__main__":
    raise SystemExit(main())
