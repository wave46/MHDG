"""Parse and evaluate Newton convergence from an MHDG solver log."""

from __future__ import annotations

import re
from math import isfinite
from dataclasses import dataclass
from pathlib import Path

from support.errors import ComparisonError


NEWTON_CONVERGENCE_FAILURE = "final Newton error is missing or exceeds tolerance"
NEWTON_FINITE_FAILURE = "final Newton error is missing or non-finite"
_ERROR_RE = re.compile(r"^\s*Error:\s*([-+0-9.eE]+)\s*$", re.MULTILINE)


@dataclass(frozen=True)
class NewtonConvergence:
    """The final Newton error and its acceptance threshold."""

    final_error: float | None
    maximum: float | None

    @property
    def passed(self) -> bool:
        return self.failure is None

    @property
    def failure(self) -> str | None:
        """Explain why the recorded error does not satisfy this check."""
        if self.final_error is None or not isfinite(self.final_error):
            return NEWTON_FINITE_FAILURE
        if self.maximum is not None and self.final_error > self.maximum:
            return NEWTON_CONVERGENCE_FAILURE
        return None

    def as_report(self) -> dict[str, bool | float | None]:
        """Return the stable convergence section used in comparison reports."""
        return {
            "passed": self.passed,
            "final_newton_error": self.final_error,
            "maximum": self.maximum,
        }


def read_newton_convergence(
    log_path: Path,
    maximum: float | None,
) -> NewtonConvergence:
    """Read the last Newton error and optionally enforce an upper bound."""
    try:
        solver_output = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ComparisonError(f"cannot read solver log {log_path}: {exc}") from exc

    errors = _ERROR_RE.findall(solver_output)
    final_error = float(errors[-1]) if errors else None
    return NewtonConvergence(final_error=final_error, maximum=maximum)
