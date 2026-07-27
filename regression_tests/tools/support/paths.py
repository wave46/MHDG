"""Explicit path contracts for regression-harness inputs and records."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from support.errors import PathError


def require_file(path: Path, label: str) -> Path:
    """Resolve a direct path that must identify an existing regular file."""
    path = path.expanduser()
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise PathError(f"{label} file does not exist: {path}") from exc
    except OSError as exc:
        raise PathError(f"cannot resolve {label} file {path}: {exc}") from exc
    if not resolved.is_file():
        raise PathError(f"{label} is not a file: {resolved}")
    return resolved


def require_directory(path: Path, label: str) -> Path:
    """Resolve a direct path that must identify an existing directory."""
    path = path.expanduser()
    try:
        resolved = path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise PathError(f"{label} directory does not exist: {path}") from exc
    except OSError as exc:
        raise PathError(f"cannot resolve {label} directory {path}: {exc}") from exc
    if not resolved.is_dir():
        raise PathError(f"{label} is not a directory: {resolved}")
    return resolved


def recorded_file(value: Any, label: str, base: Path | None = None) -> Path:
    """Resolve a file path that must be present as a string in recorded data."""
    if not isinstance(value, str) or not value:
        raise PathError(f"{label} is not recorded")
    path = Path(value).expanduser()
    if base is not None and not path.is_absolute():
        path = base / path
    return require_file(path, label)


def recorded_directory(value: Any, label: str) -> Path:
    """Resolve a directory path that must be present in recorded data."""
    if not isinstance(value, str) or not value:
        raise PathError(f"{label} is not recorded")
    return require_directory(Path(value), label)
