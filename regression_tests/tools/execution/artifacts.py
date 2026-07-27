"""Inspect solver logs, outputs, and provenance files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from support.errors import BundleError
from support.files import file_identity


FATAL_LOG_MARKERS = (
    "Error opening source file:",
    "Error opening destination file:",
    "Error   : Unable to open file",
)
MAX_FATAL_LOG_MESSAGES = 20


def fatal_log_messages(*paths: Path) -> list[str]:
    """Return bounded solver file-error messages from captured logs."""
    messages = []
    for path in paths:
        with path.open("r", encoding="utf-8", errors="replace") as stream:
            for line_number, line in enumerate(stream, start=1):
                if not any(marker in line for marker in FATAL_LOG_MARKERS):
                    continue
                messages.append(f"{path.name}:{line_number}: {line.rstrip()}")
                if len(messages) == MAX_FATAL_LOG_MESSAGES:
                    return messages
    return messages


def output_records(run_directory: Path) -> list[dict[str, Any]]:
    """Describe every regular output file produced by one solver run."""
    output_directory = run_directory / "outputs"
    return [
        file_record(path, path.relative_to(run_directory).as_posix())
        for path in sorted(output_directory.rglob("*"))
        if path.is_file()
    ]


def optional_file_record(
    settings: dict[str, str],
    key: str,
) -> dict[str, Any] | None:
    """Describe an optional absolute file configured in settings."""
    value = settings.get(key)
    if not value:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute() or not path.is_file():
        raise BundleError(f"{key} must be an absolute path to a file")
    return file_record(path, str(path))


def file_record(path: Path, display_path: str) -> dict[str, Any]:
    """Describe a file using its chosen display path and stable identity."""
    return {"path": display_path, **file_identity(path)}
