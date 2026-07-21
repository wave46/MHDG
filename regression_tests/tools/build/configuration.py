"""Resolve and validate one regression build configuration."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

from build.models import BuildConfiguration, BuildError
from build.process import git_state
from bundle.settings import read_settings
from support.environments import source_environment


DEFAULT_JOBS = 8


def configure_build(
    settings_path: Path,
    repository_root: Path,
    jobs: int | None,
) -> BuildConfiguration:
    """Resolve settings, source the toolchain, and create an isolated build area."""
    settings = read_settings(settings_path)
    repository_root = repository_root.resolve()
    library_directory = repository_root / "lib"
    test_directory = repository_root / "test"
    _require_directory(library_directory, "solver library directory")

    environment_script = Path(
        settings.get(
            "MHDG_ENVIRONMENT_SCRIPT",
            str(library_directory / "Make.inc" / "init_vars_libs.sh"),
        )
    ).expanduser()
    environment_script, environment = source_environment(environment_script)

    selected_jobs = jobs or _configured_jobs(settings)
    revision, changes = git_state(repository_root)
    build_root = _build_root(settings)
    build_id = _build_id(revision, bool(changes))
    build_directory = build_root / build_id
    binary_directory = build_directory / "bin"
    log_directory = build_directory / "logs"
    try:
        binary_directory.mkdir(parents=True)
        log_directory.mkdir()
    except OSError as exc:
        raise BuildError(
            f"cannot create build directory {build_directory}: {exc}"
        ) from exc

    return BuildConfiguration(
        settings=settings,
        repository_root=repository_root,
        library_directory=library_directory,
        test_directory=test_directory,
        build_root=build_root,
        build_id=build_id,
        build_directory=build_directory,
        binary_directory=binary_directory,
        log_directory=log_directory,
        environment_script=environment_script,
        environment=environment,
        jobs=selected_jobs,
        revision=revision,
        changes=changes,
    )


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


def _configured_jobs(settings: dict[str, str]) -> int:
    value = settings.get("MHDG_REGRESSION_BUILD_JOBS")
    if value is None:
        return DEFAULT_JOBS
    try:
        return parse_build_jobs(value)
    except argparse.ArgumentTypeError as exc:
        raise BuildError(str(exc)) from exc


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
        raise BuildError(
            f"cannot create regression build root {root}: {exc}"
        ) from exc
    return root.resolve()


def _build_id(revision: str, dirty: bool) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    suffix = "-dirty" if dirty else ""
    return f"{timestamp}-{revision[:8]}{suffix}"


def _require_directory(path: Path, label: str) -> None:
    if not path.is_dir():
        raise BuildError(f"{label} does not exist: {path}")
