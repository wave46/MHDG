"""Run and record external build and toolchain commands."""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path
from typing import TextIO

from build.models import BuildError


def run_logged(
    commands: tuple[list[str], ...],
    working_directory: Path,
    environment: dict[str, str],
    log_path: Path,
) -> None:
    """Run commands sequentially while streaming output to one build log."""
    environment = {**environment, "PWD": str(working_directory)}
    try:
        with log_path.open("w", encoding="utf-8") as log:
            for command in commands:
                rendered = shlex.join(command)
                print(f"build: {rendered}", flush=True)
                print(f"$ {rendered}", file=log, flush=True)
                _stream_command(command, working_directory, environment, log)
    except OSError as exc:
        raise BuildError(f"cannot run build; see {log_path}: {exc}") from exc


def git_state(repository_root: Path) -> tuple[str, list[str]]:
    """Return the current revision and porcelain worktree changes."""
    revision = _capture(["git", "rev-parse", "HEAD"], repository_root).strip()
    status = _capture(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        repository_root,
    )
    return revision, [line for line in status.splitlines() if line]


def toolchain_record(environment: dict[str, str]) -> dict[str, str | None]:
    """Capture available build-tool versions without requiring every tool."""
    return {
        "make": _optional_version(["make", "--version"], environment),
        "fortran_compiler": _optional_version(
            ["mpifort", "--version"],
            environment,
        ),
        "compiler_command": _optional_version(["mpifort", "-show"], environment),
    }


def _stream_command(
    command: list[str],
    working_directory: Path,
    environment: dict[str, str],
    log: TextIO,
) -> None:
    process = subprocess.Popen(
        command,
        cwd=working_directory,
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


def _capture(
    command: list[str],
    working_directory: Path,
    environment: dict[str, str] | None = None,
) -> str:
    try:
        completed = subprocess.run(
            command,
            cwd=working_directory,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise BuildError(f"cannot run {shlex.join(command)}") from exc
    return completed.stdout


def _optional_version(
    command: list[str],
    environment: dict[str, str],
) -> str | None:
    try:
        output = _capture(command, Path.cwd(), environment)
    except BuildError:
        return None
    return output.splitlines()[0] if output else None
