"""Shared filesystem operations for regression-harness artifacts."""

from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_digest(path: Path) -> str:
    """Return the hexadecimal SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_within(path: Path, directory: Path) -> bool:
    """Return whether a resolved path is contained by a directory."""
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True
