#!/usr/bin/env python3
"""Compare one completed fixed-mesh run with its bundled reference."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Any

from check_bundle import load_case_definition
from comparison.convergence import (
    NEWTON_CONVERGENCE_FAILURE,
    read_newton_convergence,
)
from comparison.metrics import format_metric, maximum_metric
from comparison.hdf5 import compare_hdf5_files
from support.documents import load_json, write_json_atomic
from support.errors import ComparisonError, HarnessError
from support.files import file_identity
from support.paths import recorded_file, require_directory, require_file
from support.time import utc_now


TIME_SAVE_RE = re.compile(r"_\d{4}\.h5$")
OUTPUT_RE = re.compile(r"Output written to file\s+(.+\.h5)\s*$", re.MULTILINE)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_directory", metavar="RUN_DIRECTORY", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument(
        "--tolerance-profile",
        dest="tolerance_profile",
    )
    parser.add_argument("--report", type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--tolerances", required=True, type=Path, help=argparse.SUPPRESS
    )
    args = parser.parse_args(argv)

    try:
        report_path, report = compare_run(
            args.run_directory,
            args.cases,
            args.tolerances,
            args.candidate,
            args.reference,
            args.tolerance_profile,
            args.report,
        )
    except HarnessError as exc:
        print(f"comparison failed: {exc}", file=sys.stderr)
        return 1

    print(f"comparison {report['status']}: {report_path}")
    print(f"candidate: {report['candidate']}")
    print(f"reference: {report['reference']}")
    print(f"tolerance profile: {report['tolerance_profile']['id']}")
    print_fixed_summary(report)
    return 0 if report["status"] == "passed" else 1


def compare_run(
    run_directory: Path,
    case_dir: Path,
    tolerances_path: Path,
    candidate_override: Path | None = None,
    reference_override: Path | None = None,
    tolerance_profile_override: str | None = None,
    report_override: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Compare a completed run and atomically write its JSON report."""
    run_directory = require_directory(run_directory, "run directory")
    plan = load_json(run_directory / "run_plan.json", "run plan")
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    if metadata.get("status") != "completed":
        raise ComparisonError("run metadata status is not completed")

    case = load_case_definition(plan["case_id"], case_dir)
    workflow = case["workflows"].get(plan["workflow_id"])
    if workflow is None:
        raise ComparisonError("run plan refers to an unknown workflow")
    profile_id, tolerances = _tolerance_profile(
        tolerances_path,
        workflow,
        plan["layout_id"],
        tolerance_profile_override,
    )

    candidate = select_candidate(run_directory, metadata, candidate_override)
    reference = _input_path(
        run_directory, reference_override, "inputs/reference.h5", "reference"
    )
    hdf5_report = compare_hdf5_files(reference, candidate, tolerances)

    convergence = read_newton_convergence(
        run_directory / "stdout.log", tolerances["newton_error_max"]
    )
    failures = list(hdf5_report["failures"])
    if not convergence.passed:
        failures.append(NEWTON_CONVERGENCE_FAILURE)

    report = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "passed" if not failures else "failed",
        "run_directory": str(run_directory),
        "case_id": plan["case_id"],
        "workflow_id": plan["workflow_id"],
        "layout_id": plan["layout_id"],
        "candidate": str(candidate),
        "reference": str(reference),
        "files": {
            "candidate": {"path": str(candidate), **file_identity(candidate)},
            "reference": {"path": str(reference), **file_identity(reference)},
        },
        "tolerance_profile": {"id": profile_id, **tolerances},
        "convergence": convergence.as_report(),
        "hdf5": hdf5_report,
        "failures": failures,
    }

    report_path = report_override or run_directory / "comparison.json"
    if not report_path.is_absolute():
        report_path = report_path.resolve()
    if report_path.resolve() in {candidate, reference}:
        raise ComparisonError("comparison report cannot replace an HDF5 input")
    write_json_atomic(report_path, report, "comparison report")
    return report_path, report


def print_fixed_summary(report: dict[str, Any]) -> None:
    """Print the human-readable summary for a fixed-mesh report."""
    convergence = report["convergence"]
    print(
        "Newton error: "
        f"{_pass_label(convergence['passed'])} "
        f"{format_metric(convergence['final_newton_error'])} "
        f"<= {format_metric(convergence['maximum'])}"
    )

    mesh = report["hdf5"].get("mesh", {})
    coordinates = mesh.get("coordinates", {})
    connectivity = mesh.get("connectivity", {})
    connectivity_passed = bool(connectivity) and all(
        item.get("passed", False) for item in connectivity.values()
    )
    mesh_passed = connectivity_passed and coordinates.get("passed", False)
    print(
        "Mesh: "
        f"{_pass_label(mesh_passed)} connectivity="
        f"{_pass_label(connectivity_passed)} max|dX|="
        f"{format_metric(coordinates.get('maximum_absolute_error'))}"
    )

    solution = report["hdf5"].get("solution", {}).get("datasets", {})
    for name in ("u", "q", "u_tilde"):
        dataset = solution.get(name, {})
        metrics = list(dataset.get("equations", {}).values())
        print(
            f"solution/{name}: {_pass_label(dataset.get('passed', False))} "
            f"max relL2={format_metric(maximum_metric(metrics, 'relative_l2'))} "
            f"max nLinf={format_metric(maximum_metric(metrics, 'normalized_linf'))}"
        )

    transport = report["hdf5"].get("transport_1d", {})
    transport_metrics = list(transport.get("datasets", {}).values())
    if transport.get("present", False):
        print(
            "transport_1d: "
            f"{_pass_label(transport.get('passed', False))} "
            "max relL2="
            f"{format_metric(maximum_metric(transport_metrics, 'relative_l2'))} "
            "max nLinf="
            f"{format_metric(maximum_metric(transport_metrics, 'normalized_linf'))}"
        )

    failures = report["failures"]
    for failure in failures[:10]:
        print(f"FAIL: {failure}")
    if len(failures) > 10:
        print(f"... {len(failures) - 10} more failures; see the JSON report")


def _pass_label(passed: bool) -> str:
    return "PASS" if passed else "FAIL"


def select_candidate(
    run_directory: Path,
    metadata: dict[str, Any],
    override: Path | None = None,
) -> Path:
    """Select the final HDF5 result recorded for a completed run."""
    if override is not None:
        return _input_path(run_directory, override, "", "candidate")

    recorded_outputs = metadata.get("hdf5_outputs")
    if not isinstance(recorded_outputs, list) or not recorded_outputs:
        raise ComparisonError("run metadata contains no HDF5 output")
    declared = [
        recorded_file(value, "HDF5 output", run_directory)
        for value in recorded_outputs
    ]

    stdout_path = run_directory / "stdout.log"
    if stdout_path.is_file():
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

    final_outputs = [path for path in declared if not TIME_SAVE_RE.search(path.name)]
    if len(final_outputs) == 1:
        return final_outputs[0]
    if len(declared) == 1:
        return declared[0]
    raise ComparisonError("cannot select one final HDF5 output; use --candidate")


def _tolerance_profile(
    path: Path,
    workflow: dict[str, Any],
    layout_id: str,
    override: str | None,
) -> tuple[str, dict[str, Any]]:
    document = load_json(path, "tolerance definitions")
    profile_id = override or workflow.get("tolerance_profile")
    if override is None and layout_id != workflow.get("default_layout"):
        profile_id = workflow.get("cross_layout_tolerance_profile", profile_id)
    profiles = document.get("profiles", {})
    if not profile_id or profile_id not in profiles:
        raise ComparisonError(f"unknown tolerance profile: {profile_id}")

    profile = profiles[profile_id]
    required = {
        "newton_error_max",
        "mesh_coordinate_atol",
        "relative_l2_max",
        "normalized_linf_max",
    }
    missing = sorted(required - profile.keys())
    if missing:
        raise ComparisonError(f"tolerance profile is missing: {', '.join(missing)}")
    return profile_id, profile


def _input_path(
    run_directory: Path,
    override: Path | None,
    default: str,
    label: str,
) -> Path:
    path = override if override is not None else Path(default)
    path = path.expanduser()
    if not path.is_absolute():
        path = run_directory / path
    return require_file(path, label)
if __name__ == "__main__":
    raise SystemExit(main())
