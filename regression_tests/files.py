"""File identities and direct/recorded path contracts."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from .support import PathError

def sha256_digest(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_identity(path: Path) -> dict[str, Any]:
    """Return the stable size and checksum fields shared by harness records."""
    try:
        size = path.stat().st_size
        digest = sha256_digest(path)
    except OSError as exc:
        raise PathError(f"cannot inspect file {path}: {exc}") from exc
    return {"size_bytes": size, "sha256": digest}


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
