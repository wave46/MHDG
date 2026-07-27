"""Validate suite summaries before golden-reference publication."""

from __future__ import annotations

from typing import Any, Literal

from references.matrix.collection import validate_matrix_summary
from support.errors import BundleError
from support.identifiers import IDENTIFIER_RE


def validate_promotion_summary(
    summary: dict[str, Any],
) -> Literal["canonical", "matrix"]:
    """Validate accepted suite evidence and identify its reference form."""
    required = {
        "schema_version",
        "status",
        "suite_id",
        "run_id",
        "case_id",
        "workflow_ids",
        "layout_ids",
        "comparison_mode",
        "results",
    }
    missing = sorted(required - summary.keys())
    if missing:
        raise BundleError(f"suite summary is missing: {', '.join(missing)}")
    if summary["schema_version"] != 2 or summary["status"] != "passed":
        raise BundleError("only a passing version-2 suite summary can be promoted")
    for name in ("suite_id", "run_id", "case_id"):
        if not isinstance(summary[name], str) or not IDENTIFIER_RE.fullmatch(
            summary[name]
        ):
            raise BundleError(f"suite summary has invalid {name}")

    if summary["comparison_mode"] == "immediate":
        _validate_canonical_summary(summary)
        return "canonical"
    validate_matrix_summary(summary)
    return "matrix"


def _validate_canonical_summary(summary: dict[str, Any]) -> None:
    workflow_ids = summary["workflow_ids"]
    if (
        not isinstance(workflow_ids, list)
        or len(workflow_ids) != 1
        or not isinstance(workflow_ids[0], str)
        or not IDENTIFIER_RE.fullmatch(workflow_ids[0])
    ):
        raise BundleError("canonical suite must declare exactly one workflow")
    results = summary["results"]
    if not isinstance(results, list) or not results:
        raise BundleError("suite summary contains no layout results")
    layouts = []
    for result in results:
        if not isinstance(result, dict) or any(
            result.get(name) != expected
            for name, expected in (
                ("status", "passed"),
                ("run_status", "completed"),
                ("comparison_status", "passed"),
            )
        ):
            raise BundleError("suite summary contains a non-passing layout")
        layout = result.get("layout_id")
        if not isinstance(layout, str) or not IDENTIFIER_RE.fullmatch(layout):
            raise BundleError("suite summary has an invalid layout identifier")
        layouts.append(layout)
    if len(set(layouts)) != len(layouts):
        raise BundleError("suite summary contains duplicate layouts")
