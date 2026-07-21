"""Compare one completed fixed-mesh run with its bundled reference."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from comparison.fixed.hdf5 import compare_hdf5_files
from comparison.inputs import ComparisonInputs, ComparisonOverrides
from comparison.shared.convergence import NewtonConvergence, read_newton_convergence
from comparison.shared.outputs import resolve_run_file, select_candidate
from comparison.shared.report import (
    comparison_report_fields,
    merge_comparison_failures,
    save_comparison_report,
)
from comparison.shared.tolerances import load_fixed_tolerances
from support.files import file_identity


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
    failures = merge_comparison_failures(hdf5_report["failures"], convergence)
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
    report_path = save_comparison_report(
        inputs.run_directory,
        report_path,
        report,
        (candidate, reference),
    )
    return report_path, report


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
        **comparison_report_fields(
            inputs,
            tolerance_profile_id,
            tolerances,
            convergence,
            failures,
        ),
        "candidate": str(candidate),
        "reference": str(reference),
        "files": {
            "candidate": {"path": str(candidate), **file_identity(candidate)},
            "reference": {"path": str(reference), **file_identity(reference)},
        },
        "hdf5": hdf5_report,
    }
