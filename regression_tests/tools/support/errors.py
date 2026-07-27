"""Errors shared across regression harness subsystems."""


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
