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
        "results",
    }
    missing = sorted(required - summary.keys())
    if missing:
        raise BundleError(f"suite summary is missing: {', '.join(missing)}")
    if summary["schema_version"] != 1 or summary["status"] != "passed":
        raise BundleError("only a passing version-1 suite summary can be promoted")
    for name in ("suite_id", "run_id", "case_id"):
        if not isinstance(summary[name], str) or not IDENTIFIER_RE.fullmatch(
            summary[name]
        ):
            raise BundleError(f"suite summary has invalid {name}")

    if "workflow_id" in summary:
        _validate_canonical_summary(summary)
        return "canonical"
    validate_matrix_summary(summary)
    return "matrix"


def _validate_canonical_summary(summary: dict[str, Any]) -> None:
    workflow_id = summary["workflow_id"]
    if not isinstance(workflow_id, str) or not IDENTIFIER_RE.fullmatch(workflow_id):
        raise BundleError("suite summary has invalid workflow_id")
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
