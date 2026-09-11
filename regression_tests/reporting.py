"""Shared terminal conventions; scientific details retain their existing reports."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

import shlex
import sys


def error(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 1


def status(subject, outcome, path) -> None:
    print(f"{subject} {outcome}: {path}", flush=True)


def prepared(run) -> None:
    from .prepare import PreparedRun

    status("run", "prepared", run.path)
    if isinstance(run, PreparedRun):
        print(f"command: {shlex.join(run.command)}", flush=True)
    else:
        print(f"stages: {len(run.stages)}", flush=True)


def run_result(result) -> int:
    print(f"status: {result.status}")
    print(f"solver exit code: {result.exit_code}")
    print(f"runtime: {result.duration_seconds:.3f} s")
    print(f"HDF5 outputs: {len(result.hdf5_outputs)}")
    return 0 if result.status == "completed" else 1


def bundle(summary) -> None:
    print(f"case: {summary.case_id}; version: {summary.bundle_version}")
    print(
        f"verified {summary.verified_artifact_count}/{summary.artifact_count} artifacts "
        f"({summary.verified_bytes} bytes)"
    )
    for warning in summary.warnings:
        print(f"warning: {warning}", file=sys.stderr)


def comparison(policy, path, report) -> None:
    print(f"comparison policy: {policy}")
    if "method_selection" in report:
        print(f"method selection: {report['method_selection']['reason']}")
    status("comparison", report["status"], path)
    if policy in {"fixed_hdf5", "mesh_independent"}:
        files = report if policy == "fixed_hdf5" else report["files"]
        print(f"candidate: {files['candidate']}")
        print(f"reference: {files['reference']}")
        print(f"tolerance profile: {report['tolerance_profile']['id']}")
        if policy == "fixed_hdf5":
            print_fixed_summary(report)
        else:
            print_adaptive_summary(report)
    else:
        print(f"stages checked: {report['checked_stage_count']}/{report['total_stage_count']}")
        for failure in report["failures"]:
            print(f"FAIL: {failure}")


def diagnostics(report, path) -> None:
    status("balance diagnostics", report["status"], path)
    for failure in report["failures"]:
        print(f"FAIL: {failure}")


def maximum_metric(
    metrics: Iterable[Mapping[str, float | None]], key: str
) -> float | None:
    """Return the largest available metric value."""
    values = (metric.get(key) for metric in metrics)
    return max((value for value in values if value is not None), default=None)


def format_metric(value: float | None) -> str:
    """Format an available metric or a readable missing-value marker."""
    return f"{value:.3e}" if value is not None else "n/a"


def print_fixed_summary(report: dict[str, Any]) -> None:
    """Print the human-readable summary for a fixed-mesh report."""
    convergence = report["convergence"]
    maximum = convergence["maximum"]
    acceptance = (
        "finite"
        if maximum is None
        else f"<= {format_metric(maximum)}"
    )
    print(
        "Newton error: "
        f"{_pass_label(convergence['passed'])} "
        f"{format_metric(convergence['final_newton_error'])} "
        f"({acceptance})"
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

    magnetic = report["hdf5"].get("magnetic", {})
    magnetic_metrics = list(magnetic.get("datasets", {}).values())
    if magnetic.get("present", False):
        print(
            "magnetic: "
            f"{_pass_label(magnetic.get('passed', False))} "
            "max relL2="
            f"{format_metric(maximum_metric(magnetic_metrics, 'relative_l2'))} "
            "max nLinf="
            f"{format_metric(maximum_metric(magnetic_metrics, 'normalized_linf'))}"
        )

    failures = report["failures"]
    for failure in failures[:10]:
        print(f"FAIL: {failure}")
    if len(failures) > 10:
        print(f"... {len(failures) - 10} more failures; see the JSON report")


def print_adaptive_summary(report: dict[str, Any]) -> None:
    """Print the human-readable summary for an adaptive-mesh report."""
    sampling = report["sampling"]
    print(
        f"common points: {sampling['common_points']}/{sampling['point_count']} "
        f"({sampling['common_coverage']:.3%})"
    )
    for dataset_name, dataset in report["datasets"].items():
        metrics = list(dataset["equations"].values())
        print(
            f"{dataset_name}: max relL2="
            f"{format_metric(maximum_metric(metrics, 'relative_l2'))} "
            f"max nLinf="
            f"{format_metric(maximum_metric(metrics, 'normalized_linf'))}"
        )


def _pass_label(passed: bool) -> str:
    return "PASS" if passed else "FAIL"
