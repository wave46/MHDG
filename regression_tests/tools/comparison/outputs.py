"""Select and resolve files produced or consumed by a completed run."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from support.errors import ComparisonError
from support.paths import recorded_file, require_file


TIME_SAVE_RE = re.compile(r"_\d{4}\.h5$")
OUTPUT_RE = re.compile(r"Output written to file\s+(.+\.h5)\s*$", re.MULTILINE)


def select_candidate(
    run_directory: Path,
    metadata: dict[str, Any],
    override: Path | None = None,
) -> Path:
    """Select the final HDF5 result recorded for a completed run."""
    if override is not None:
        return resolve_run_file(run_directory, override, "", "candidate")

    recorded_outputs = metadata.get("hdf5_outputs")
    if not isinstance(recorded_outputs, list) or not recorded_outputs:
        raise ComparisonError("run metadata contains no HDF5 output")
    declared = [
        recorded_file(value, "HDF5 output", run_directory)
        for value in recorded_outputs
    ]

    logged_output = _last_logged_output(run_directory, declared)
    if logged_output is not None:
        return logged_output

    final_outputs = [path for path in declared if not TIME_SAVE_RE.search(path.name)]
    if len(final_outputs) == 1:
        return final_outputs[0]
    if len(declared) == 1:
        return declared[0]
    raise ComparisonError("cannot select one final HDF5 output; use --candidate")


def resolve_run_file(
    run_directory: Path,
    override: Path | None,
    default: str,
    label: str,
) -> Path:
    """Resolve an override or run-relative default file."""
    path = override if override is not None else Path(default)
    path = path.expanduser()
    if not path.is_absolute():
        path = run_directory / path
    return require_file(path, label)


def _last_logged_output(
    run_directory: Path,
    declared: list[Path],
) -> Path | None:
    stdout_path = run_directory / "stdout.log"
    if not stdout_path.is_file():
        return None

    text = stdout_path.read_text(encoding="utf-8", errors="replace")
    for match in reversed(OUTPUT_RE.findall(text)):
        path = Path(match.strip()).expanduser()
        path = (
            path.resolve()
            if path.is_absolute()
            else (run_directory / path).resolve()
        )
        if path in declared:
            return path
    return None
