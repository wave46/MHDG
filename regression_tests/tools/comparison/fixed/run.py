"""Compare one completed fixed-mesh run with its bundled reference."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from comparison.shared.convergence import (
    NEWTON_CONVERGENCE_FAILURE,
    NewtonConvergence,
    read_newton_convergence,
)
from comparison.fixed.hdf5 import compare_hdf5_files
from comparison.inputs import ComparisonInputs, ComparisonOverrides
from comparison.shared.outputs import resolve_run_file, select_candidate
from comparison.shared.tolerances import load_fixed_tolerances
from support.documents import write_json_atomic
from support.errors import ComparisonError
from support.files import file_identity
from support.time import utc_now


def compare_fixed_run(
    inputs: ComparisonInputs,
    overrides: ComparisonOverrides,
    report_path: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Compare a completed fixed-mesh run and save its JSON report."""
    profile_id, tolerances = load_fixed_tolerances(
        inputs.tolerances_path,
        inputs.workflow,
        inputs.plan["layout_id"],
        overrides.tolerance_profile,
    )
    candidate = select_candidate(
        inputs.run_directory,
        inputs.metadata,
        overrides.candidate,
    )
    reference = resolve_run_file(
        inputs.run_directory,
        overrides.reference,
        "inputs/reference.h5",
        "reference",
    )

    hdf5_report = compare_hdf5_files(reference, candidate, tolerances)
    convergence = read_newton_convergence(
        inputs.run_directory / "stdout.log", tolerances["newton_error_max"]
    )
    failures = _comparison_failures(hdf5_report, convergence)
    report = _comparison_report(
        inputs,
        candidate,
        reference,
        profile_id,
        tolerances,
        convergence,
        hdf5_report,
        failures,
    )
    report_path = _save_report(
        inputs.run_directory,
        candidate,
        reference,
        report_path,
        report,
    )
    return report_path, report


def _comparison_failures(
    hdf5_report: dict[str, Any],
    convergence: NewtonConvergence,
) -> list[str]:
    failures = list(hdf5_report["failures"])
    if not convergence.passed:
        failures.append(NEWTON_CONVERGENCE_FAILURE)
    return failures


def _comparison_report(
    inputs: ComparisonInputs,
    candidate: Path,
    reference: Path,
    tolerance_profile_id: str,
    tolerances: dict[str, Any],
    convergence: NewtonConvergence,
    hdf5_report: dict[str, Any],
    failures: list[str],
) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "created_utc": utc_now(),
        "status": "passed" if not failures else "failed",
        "run_directory": str(inputs.run_directory),
        "case_id": inputs.plan["case_id"],
        "workflow_id": inputs.plan["workflow_id"],
        "layout_id": inputs.plan["layout_id"],
        "candidate": str(candidate),
        "reference": str(reference),
        "files": {
            "candidate": {"path": str(candidate), **file_identity(candidate)},
            "reference": {"path": str(reference), **file_identity(reference)},
        },
        "tolerance_profile": {"id": tolerance_profile_id, **tolerances},
        "convergence": convergence.as_report(),
        "hdf5": hdf5_report,
        "failures": failures,
    }


def _save_report(
    run_directory: Path,
    candidate: Path,
    reference: Path,
    report_path: Path | None,
    report: dict[str, Any],
) -> Path:
    report_path = report_path or run_directory / "comparison.json"
    if not report_path.is_absolute():
        report_path = report_path.resolve()
    if report_path.resolve() in {candidate, reference}:
        raise ComparisonError("comparison report cannot replace an HDF5 input")
    write_json_atomic(report_path, report, "comparison report")
    return report_path
