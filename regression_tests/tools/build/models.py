"""Typed build configuration and results."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from support.errors import BundleError


class BuildError(BundleError):
    """Raised when the solver build cannot be completed or recorded."""


@dataclass(frozen=True)
class BuildConfiguration:
    settings: dict[str, str]
    repository_root: Path
    library_directory: Path
    test_directory: Path
    build_root: Path
    build_id: str
    build_directory: Path
    binary_directory: Path
    log_directory: Path
    environment_script: Path
    environment: dict[str, str]
    jobs: int
    revision: str
    changes: list[str]


@dataclass(frozen=True)
class BuildResult:
    path: Path
    settings_path: Path
    metadata_path: Path
    executables: dict[str, Path]
