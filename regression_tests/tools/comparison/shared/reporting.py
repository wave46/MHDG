"""Concise terminal summaries for comparison reports."""

from __future__ import annotations

from typing import Any

from comparison.shared.metrics import format_metric, maximum_metric


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
