"""Shell-environment loading for regression builds and solver runs."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from support.errors import BundleError
from support.paths import require_file


def source_environment(
    script: Path,
    base_environment: dict[str, str] | None = None,
) -> tuple[Path, dict[str, str]]:
    """Source one shell script and return its resolved path and environment."""
    script = require_file(script, "environment script")
    output = _capture_environment(script, base_environment)
    return script, _decode_environment(output)


def _capture_environment(
    script: Path,
    base_environment: dict[str, str] | None,
) -> bytes:
    """Run the shell script and capture the resulting encoded environment."""
    command = [
        "bash",
        "-c",
        'source "$1" >/dev/null && env -0',
        "mhdg-regression-environment",
        str(script),
    ]
    try:
        completed = subprocess.run(
            command,
            env=base_environment,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        raise BundleError(f"cannot source environment script {script}: {exc}") from exc
    if completed.returncode != 0:
        error = completed.stderr.decode(errors="replace").strip()
        detail = f": {error}" if error else ""
        raise BundleError(
            f"environment script failed ({completed.returncode}){detail}"
        )
    return completed.stdout


def _decode_environment(output: bytes) -> dict[str, str]:
    """Decode the NUL-separated output produced by ``env -0``."""
    try:
        return {
            os.fsdecode(key): os.fsdecode(value)
            for entry in output.split(b"\0")
            if entry
            for key, value in [entry.split(b"=", 1)]
        }
    except ValueError as exc:
        raise BundleError("environment script produced invalid environment data") from exc
