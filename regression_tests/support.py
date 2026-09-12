"""Shared errors, record identifiers and UTC timestamps."""

from __future__ import annotations

import re
from datetime import datetime, timezone

class HarnessError(ValueError):
    """Base class for user-facing harness failures."""


class BundleError(HarnessError):
    """Raised when settings or bundle data violate the regression contract."""


class MissingArtifactError(BundleError):
    """Raised when a declared bundle artifact does not exist."""


class DocumentError(HarnessError):
    """Raised when a harness document cannot be read or written."""


class ComparisonError(HarnessError):
    """Raised when comparison data violate the regression contract."""


class PathError(HarnessError):
    """Raised when a direct or recorded harness path is invalid."""


IDENTIFIER_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


def utc_run_id() -> str:
    """Return the timestamp format used for default run identifiers."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
