"""Shared filesystem operations for regression-harness artifacts."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from support.errors import PathError


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


def is_within(path: Path, directory: Path) -> bool:
    """Return whether a resolved path is contained by a directory."""
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True
