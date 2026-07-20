#!/usr/bin/env python3
"""Verify previously saved workflow-by-layout suite results."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from check_bundle import BundleError, load_case_definition
from compare_matrix import compare_completed_run
from support.documents import write_json_atomic
from support.errors import HarnessError
from support.time import utc_now


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_summary", metavar="SUITE_SUMMARY", type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--tolerances", required=True, type=Path, help=argparse.SUPPRESS
    )
    args = parser.parse_args(argv)

    try:
        path, summary = verify_suite(
            args.suite_summary, args.cases, args.tolerances
        )
    except (BundleError, HarnessError) as exc:
        print(f"suite verification failed: {exc}", file=sys.stderr)
        return 1

    _print_summary(summary, path)
    return 0 if summary["status"] == "passed" else 1


def verify_suite(
    suite_summary_path: Path,
    case_dir: Path,
    tolerances_path: Path,
) -> tuple[Path, dict[str, Any]]:
    """Compare all recorded suite runs without executing the solver."""
    suite_summary_path = _existing_file(suite_summary_path, "suite summary")
    source = _load_json(suite_summary_path, "suite summary")
    _validate_source_summary(source)
    case = load_case_definition(source["case_id"], case_dir)

    output_path = suite_summary_path.parent / "verification_summary.json"
    summary = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "running",
        "source_summary": str(suite_summary_path),
        "suite_id": source["suite_id"],
        "run_id": source["run_id"],
        "case_id": source["case_id"],
        "results": [],
    }
    write_json_atomic(output_path, summary, "verification summary")
    for source_result in source["results"]:
        summary["results"].append(
            _verify_result(source_result, case, case_dir, tolerances_path)
        )
        write_json_atomic(output_path, summary, "verification summary")

    summary["status"] = (
        "passed"
        if all(result["status"] == "passed" for result in summary["results"])
        else "failed"
    )
    summary["finished_utc"] = utc_now()
    write_json_atomic(output_path, summary, "verification summary")
    return output_path, summary


def _verify_result(
    source: dict[str, Any],
    case: dict[str, Any],
    case_dir: Path,
    tolerances_path: Path,
) -> dict[str, Any]:
    workflow_id = source.get("workflow_id")
    layout_id = source.get("layout_id")
    result = {
        "workflow_id": workflow_id,
        "layout_id": layout_id,
        "run_directory": source.get("run_directory"),
        "run_status": source.get("run_status"),
        "comparison_policy": None,
        "comparison_report": None,
        "status": "failed",
        "failures": [],
    }
    if source.get("run_status") != "completed":
        result["failures"] = ["solver run did not complete"]
        return result

    workflow = case["workflows"].get(workflow_id)
    if workflow is None:
        result["failures"] = [f"unknown workflow: {workflow_id}"]
        return result
    try:
        run_directory = Path(source["run_directory"])
        policy, report_path, report = compare_completed_run(
            run_directory, case_dir, tolerances_path
        )
        result["comparison_policy"] = policy
        result["comparison_report"] = str(report_path)
        result["failures"] = report["failures"]
        result["status"] = report["status"]
    except (BundleError, HarnessError, TypeError) as exc:
        result["failures"] = [str(exc)]
    return result


def _validate_source_summary(summary: dict[str, Any]) -> None:
    required = {"schema_version", "suite_id", "run_id", "case_id", "results"}
    missing = sorted(required - summary.keys())
    if missing:
        raise BundleError(f"suite summary is missing: {', '.join(missing)}")
    if summary["schema_version"] != 1:
        raise BundleError("suite summary has unsupported schema version")
    if not isinstance(summary["results"], list) or not summary["results"]:
        raise BundleError("suite summary contains no results")
    result_fields = {"workflow_id", "layout_id", "run_directory", "run_status"}
    for result in summary["results"]:
        if not isinstance(result, dict) or not result_fields <= result.keys():
            raise BundleError("suite summary contains an invalid result")


def _print_summary(summary: dict[str, Any], path: Path) -> None:
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


def _existing_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise BundleError(f"{label} does not exist: {path}")
    return path


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BundleError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise BundleError(f"{label} must contain a JSON object")
    return document


if __name__ == "__main__":
    raise SystemExit(main())
