"""UTC timestamps used in harness records."""

from __future__ import annotations

from datetime import datetime, timezone


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


def utc_run_id() -> str:
    """Return the timestamp format used for default run identifiers."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
