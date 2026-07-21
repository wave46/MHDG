"""Common construction and publication of comparison reports."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from comparison.inputs import ComparisonInputs
from comparison.shared.convergence import NewtonConvergence
from support.documents import write_json_atomic
from support.errors import ComparisonError
from support.time import utc_now


def merge_comparison_failures(
    failures: list[str],
    convergence: NewtonConvergence,
) -> list[str]:
    """Append a Newton-convergence failure without changing component results."""
    merged = list(failures)
    if convergence.failure is not None:
        merged.append(convergence.failure)
    return merged


def comparison_report_fields(
    inputs: ComparisonInputs,
    tolerance_profile_id: str,
    tolerances: dict[str, Any],
    convergence: NewtonConvergence,
    failures: list[str],
) -> dict[str, Any]:
    """Build fields shared by fixed- and adaptive-mesh reports."""
    return {
        "schema_version": 2,
        "created_utc": utc_now(),
        "status": "passed" if not failures else "failed",
        "run_directory": str(inputs.run_directory),
        "case_id": inputs.plan["case_id"],
        "workflow_id": inputs.plan["workflow_id"],
        "layout_id": inputs.plan["layout_id"],
        "tolerance_profile": {"id": tolerance_profile_id, **tolerances},
        "convergence": convergence.as_report(),
        "failures": failures,
    }


def save_comparison_report(
    run_directory: Path,
    report_path: Path | None,
    report: dict[str, Any],
    protected_paths: tuple[Path, ...],
) -> Path:
    """Atomically write a report without replacing any comparison input."""
    output = (report_path or run_directory / "comparison.json").expanduser().resolve()
    protected = {path.expanduser().resolve() for path in protected_paths}
    if output in protected:
        raise ComparisonError("comparison report cannot replace a comparison input")
    write_json_atomic(output, report, "comparison report")
    return output
