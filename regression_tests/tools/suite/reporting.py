"""Render compact suite command summaries."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def print_run_summary(summary: dict[str, Any], path: Path) -> None:
    """Print execution and comparison status for every suite cell."""
    print("workflow       layout        run                 comparison  runtime     result")
    for result in summary["results"]:
        runtime = result["duration_seconds"]
        runtime_text = f"{runtime:.3f} s" if runtime is not None else "n/a"
        print(
            f"{result['workflow_id']:<14} "
            f"{result['layout_id']:<13} "
            f"{(result['run_status'] or 'n/a'):<19} "
            f"{(result['comparison_status'] or 'n/a'):<11} "
            f"{runtime_text:<11} {result['status']}"
        )
        for failure in result["failures"][:3]:
            print(f"  FAIL: {failure}")
        if len(result["failures"]) > 3:
            print(f"  ... {len(result['failures']) - 3} more failures")
    _print_layout_comparisons(summary.get("comparisons", []))
    print(
        f"suite {summary['status']}: {path} "
        f"({summary['duration_seconds']:.3f} s)"
    )


def print_verification_summary(summary: dict[str, Any], path: Path) -> None:
    """Print comparison status for every recorded suite cell."""
    if summary.get("comparisons") is not None:
        _print_layout_comparisons(summary["comparisons"])
        print(f"verification {summary['status']}: {path}")
        return
    print("workflow       layout        policy             result")
    for result in summary["results"]:
        print(
            f"{str(result['workflow_id']):<14} "
            f"{str(result['layout_id']):<13} "
            f"{str(result['comparison_policy'] or 'n/a'):<18} "
            f"{result['status']}"
        )
        for failure in result["failures"][:3]:
            print(f"  FAIL: {failure}")
    print(f"verification {summary['status']}: {path}")


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
