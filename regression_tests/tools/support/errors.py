"""Errors shared across regression harness subsystems."""


class HarnessError(ValueError):
    """Base class for user-facing harness failures."""


class DocumentError(HarnessError):
    """Raised when a harness document cannot be read or written."""


class ComparisonError(HarnessError):
    """Raised when comparison data violate the regression contract."""
