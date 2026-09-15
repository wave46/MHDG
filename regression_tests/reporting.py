"""Shared terminal conventions; scientific details retain their existing reports."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any
from pathlib import Path

import shlex
import sys
import textwrap


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
    if result.status != "completed":
        from .execute import execution_failures
        for failure in execution_failures(result.path):
            print(f"FAIL: {failure}")
    return 0 if result.status == "completed" else 1


def bundle(summary) -> None:
    print(f"case: {summary.case_id}; version: {summary.bundle_version}")
    print(
        f"verified {summary.verified_artifact_count}/{summary.artifact_count} artifacts "
        f"({summary.verified_bytes} bytes)"
    )
    for warning in summary.warnings:
        print(f"warning: {warning}", file=sys.stderr)


def readiness(report) -> None:
    status("bundle readiness", report["status"], report["source"])
    print(f"case: {report['case_id']}; workflows: {', '.join(report['workflows']) or 'base bundle'}")
    print("role  |  requirement  |  presence  |  origin  |  file")
    for row in report["artifacts"]:
        print("  |  ".join(str(row[key] or "—") for key in
                           ("role", "requirement", "presence", "origin", "path")))
        if row["presence"] == "missing":
            for workflow, required in row["producers"].items():
                print(f"  producer {workflow}; requires: {', '.join(required)}")
    print("Presence check only; bundle validate checks recorded sizes and checksums.")


def golden_review(path, report) -> None:
    status("golden refresh", report["status"], path)
    print(f"producers: {len(report['producers'])}; validation suites: {len(report['checks'])}")
    for item in report["producers"]:
        previous = item.get("old_reference", {})
        print(f"  {item['workflow_id']} / {item['layout_id']}: old reference {previous.get('status', 'unavailable')}")
    print(f"candidate: {path.parent / 'candidate'}")
    print("Review refresh.json and comparison reports, then use golden publish with a reason and provenance.")


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


def _print_results(summary):
    print("workflow       layout        run                 comparison  runtime     result")
    for result in summary["results"]:
        runtime = result.get("duration_seconds")
        runtime_text = f"{runtime:.3f} s" if runtime is not None else "n/a"
        print(
            f"{result['workflow_id']:<14} "
            f"{result['layout_id']:<13} "
            f"{(result['run_status'] or 'n/a'):<19} "
            f"{(result.get('comparison_status') or 'n/a'):<11} "
            f"{runtime_text:<11} {result['status']}"
        )
        for failure in result["failures"][:3]:
            print(f"  FAIL: {failure}")
        if len(result["failures"]) > 3:
            print(f"  ... {len(result['failures']) - 3} more failures")


def print_run_summary(summary: dict[str, Any], path: Path) -> None:
    """Print execution and comparison status for every suite cell."""
    _print_results(summary)
    _print_layout_comparisons(summary.get("comparisons", []))
    _print_diagnostics(summary, path)
    print(
        f"suite {summary['status']}: {path} "
        f"({summary['duration_seconds']:.3f} s)"
    )


def print_verification_summary(summary: dict[str, Any], path: Path) -> None:
    """Print comparison status for every recorded suite cell."""
    _print_results(summary)
    _print_layout_comparisons(summary.get("comparisons", []))
    _print_diagnostics(summary, path)
    print(f"verification {summary['status']}: {path}")


def suggestions(report):
    scope = f"merge base with {report['base']} plus local changes" if report["base"] else "local changes only"
    print(f"suggestions: {scope}; {len(report['paths'])} changed files (advisory; nothing executed)")
    for check in report["checks"]:
        print(shlex.join(check["command"]))
        grouped = {}
        for path, reason in sorted(check["reasons"].items()):
            grouped.setdefault(reason, []).append(shlex.quote(path))
        for reason, paths in grouped.items():
            print(f"  {reason}:")
            print(textwrap.fill(", ".join(paths), width=100, initial_indent="    ", subsequent_indent="    ",
                                break_long_words=False, break_on_hyphens=False))
    for path in report["unmapped"]:
        print(f"unmapped: {path}; choose verification manually")
    if report["documentation"]:
        print(f"documentation: {len(report['documentation'])} files; review without solver checks")
    if not report["paths"]:
        print("no changes found; use --base REVISION to include committed work")
    if any("--build" in item["command"] for item in report["checks"]):
        print("Solver suggestions build the current tree. To share the first build, replace --build")
        print("on subsequent checks with --build-manifest PATH printed by that build.")
        print("File mapping does not establish complete coverage; assess cold convergence and diagnostic format changes separately.")


def cleanup(report):
    selected = [row for row in report["rows"] if row["selected"]]
    for row in selected or report["rows"]:
        action = report["status"] if row["selected"] else "keep"
        if row["selected"] and row["reasons"]:
            action = "blocked"
        reasons = "; ".join(row["reasons"][:3]) or ("explicit selection" if row["selected"] else "not selected")
        if len(row["reasons"]) > 3:
            reasons += f"; {len(row['reasons']) - 3} more protection records"
        print(f"{action}: {row['path']} ({row['kind']}, {row['bytes']} bytes): {reasons}")
    status("cleanup", report["status"], report["root"])
    print(f"selected: {report['selected_bytes']} bytes")


def _print_diagnostics(summary, path):
    report = summary.get("diagnostics", {})
    if report.get("outputs") or report.get("failures"):
        diagnostics(report, path)
        if report.get("outputs") and all(item["mode"] == "off" for item in report["outputs"]):
            print("  mode off: output absence checked")
    elif summary.get("checks_enabled") is False:
        status("balance diagnostics", "deferred", path)


def _print_layout_comparisons(comparisons: list[dict[str, Any]]) -> None:
    if not comparisons:
        return
    print("workflow             baseline       candidate      policy             result")
    for comparison in comparisons:
        print(
            f"{comparison['workflow_id']:<20} "
            f"{comparison['baseline_layout_id']:<14} "
            f"{comparison['candidate_layout_id']:<14} "
            f"{str(comparison['comparison_policy'] or 'n/a'):<18} "
            f"{comparison['status']}"
        )
        for failure in comparison["failures"][:3]:
            print(f"  FAIL: {failure}")
        if comparison.get("diagnostics"):
            report = comparison["diagnostics"]
            print(f"  diagnostic comparison {report['status']}" + (f": {report['reason']}" if "reason" in report else ""))
