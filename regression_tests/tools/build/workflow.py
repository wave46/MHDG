"""Orchestrate clean serial and parallel regression builds."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from build.artifacts import (
    install_executable,
    install_runtime_file,
    write_build_metadata,
    write_generated_settings,
)
from build.configuration import configure_build
from build.models import BuildConfiguration, BuildResult
from build.process import run_logged
from support.time import utc_now


MODEL = "NGammaTiTeNeutral"
DIMENSION = "2D"
COMPILE_TYPE = "opt"
VARIANTS = {"serial": "serial", "parallel": "parall"}
RUNTIME_FILE = "positionFeketeNodesTri2D.h5"


def build_solver(
    settings_path: Path,
    repository_root: Path,
    jobs: int | None = None,
) -> BuildResult:
    """Build clean serial and parallel variants and return their local settings."""
    configuration = configure_build(settings_path, repository_root, jobs)
    started_utc = utc_now()
    commands, executables = _build_variants(configuration)
    runtime_file = install_runtime_file(configuration, RUNTIME_FILE)
    metadata_path = write_build_metadata(
        configuration,
        started_utc,
        commands,
        executables,
        runtime_file,
        _build_profile(configuration.jobs),
    )
    settings_path = write_generated_settings(
        configuration,
        executables,
        metadata_path,
    )
    return BuildResult(
        configuration.build_directory,
        settings_path,
        metadata_path,
        executables,
    )


def _build_variants(
    configuration: BuildConfiguration,
) -> tuple[list[list[str]], dict[str, Path]]:
    commands: list[list[str]] = []
    executables: dict[str, Path] = {}
    for variant, mode in VARIANTS.items():
        target = f"MHDG-{MODEL}-{mode}-{DIMENSION}"
        clean = ["make", "clean"]
        compile_command = [
            "make",
            f"-j{configuration.jobs}",
            f"MODE={mode}",
            f"COMPTYPE={COMPILE_TYPE}",
            f"MHDG_GIT_COMMIT={configuration.revision}",
            f"MHDG_GIT_DIRTY={'true' if configuration.changes else 'false'}",
            f"MHDG_BUILD_ID={configuration.build_id}",
            target,
        ]
        commands.extend((clean, compile_command))
        run_logged(
            (clean, compile_command),
            configuration.library_directory,
            configuration.environment,
            configuration.log_directory / f"{variant}.log",
        )
        executables[variant] = install_executable(
            configuration.library_directory / target,
            configuration.binary_directory / target,
            f"built {variant} executable",
        )
    return commands, executables


def _build_profile(jobs: int) -> dict[str, Any]:
    return {
        "model": MODEL,
        "dimension": DIMENSION,
        "compile_type": COMPILE_TYPE,
        "jobs": jobs,
        "variants": VARIANTS,
    }
