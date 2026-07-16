#!/usr/bin/env python3
"""Run and compare each layout in a tracked regression suite."""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from build_solver import build_solver
from check_bundle import (
    BundleError,
    bundle_root_from_settings,
    load_case_definition,
    load_validated_json,
    read_settings,
    validate_bundle_root,
)
from compare_hdf5 import ComparisonError
from compare_run import compare_run
from prepare_run import prepare_run
from run_case import execute_run


SUITE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_id", metavar="SUITE")
    parser.add_argument("--run-id")
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--layouts", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--suites", required=True, type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--tolerances", required=True, type=Path, help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--require-bundle-class",
        choices=("candidate", "golden"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--build",
        action="store_true",
        help="build clean serial and parallel executables before running",
    )
    parser.add_argument("--build-jobs", type=_positive_integer, metavar="N")
    parser.add_argument(
        "--repository-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)

    try:
        settings_path = args.settings
        if args.build:
            build = build_solver(settings_path, args.repository_root, args.build_jobs)
            settings_path = build.settings_path
            print(f"using build: {build.path}")
        elif args.build_jobs is not None:
            raise BundleError("--build-jobs requires --build")
        summary_path, summary = run_suite(
            settings_path,
            args.suite_id,
            args.cases,
            args.layouts,
            args.suites,
            args.tolerances,
            args.run_id,
            args.require_bundle_class,
        )
    except BundleError as exc:
        print(f"suite failed: {exc}", file=sys.stderr)
        return 1

    _print_summary(summary, summary_path)
    return 0 if summary["status"] == "passed" else 1


def run_suite(
    settings_path: Path,
    suite_id: str,
    case_dir: Path,
    layouts_path: Path,
    suites_path: Path,
    tolerances_path: Path,
    run_id: str | None = None,
    required_bundle_class: str | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Execute all suite layouts, continue after failures, and write a summary."""
    suite = _load_suite(suite_id, suites_path, layouts_path, case_dir)
    settings = read_settings(settings_path)
    bundle_root = bundle_root_from_settings(settings)
    validate_bundle_root(bundle_root, case_dir)
    if required_bundle_class is not None:
        _require_bundle_class(bundle_root, case_dir, required_bundle_class)
    run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if not SUITE_ID_RE.fullmatch(run_id):
        raise BundleError(f"invalid suite run identifier: {run_id}")

    suite_dir = _suite_directory(settings, suite_id, run_id)
    started_utc = _utc_now()
    started_clock = time.monotonic()
    print(f"suite: {suite_id} ({run_id})")
    results = []
    for layout_id in suite["layouts"]:
        print(f"running {layout_id} ...", flush=True)
        results.append(
            _run_layout(
                settings_path,
                settings,
                suite,
                layout_id,
                run_id,
                case_dir,
                layouts_path,
                tolerances_path,
            )
        )

    summary = {
        "schema_version": 1,
        "started_utc": started_utc,
        "finished_utc": _utc_now(),
        "duration_seconds": time.monotonic() - started_clock,
        "status": "passed" if all(_passed(item) for item in results) else "failed",
        "suite_id": suite_id,
        "run_id": run_id,
        "description": suite["description"],
        "case_id": suite["case_id"],
        "workflow_id": suite["workflow_id"],
        "results": results,
    }
    summary_path = suite_dir / "suite_summary.json"
    _write_json(summary_path, summary)
    return summary_path, summary


def _run_layout(
    settings_path: Path,
    settings: dict[str, str],
    suite: dict[str, Any],
    layout_id: str,
    run_id: str,
    case_dir: Path,
    layouts_path: Path,
    tolerances_path: Path,
) -> dict[str, Any]:
    result = {
        "layout_id": layout_id,
        "status": "error",
        "run_status": None,
        "comparison_status": None,
        "duration_seconds": None,
        "run_directory": None,
        "comparison_report": None,
        "failures": [],
    }
    try:
        prepared = prepare_run(
            settings_path,
            suite["case_id"],
            suite["workflow_id"],
            layout_id,
            case_dir,
            layouts_path,
            run_id,
            validate_bundle=False,
        )
        result["run_directory"] = str(prepared.path)
        run = execute_run(prepared, settings)
        result["run_status"] = run.status
        result["duration_seconds"] = run.duration_seconds
        result["status"] = run.status
        if run.status != "completed":
            result["failures"] = [f"solver run status is {run.status}"]
            return result

        report_path, comparison = compare_run(
            run.path, case_dir, tolerances_path
        )
        result["comparison_status"] = comparison["status"]
        result["comparison_report"] = str(report_path)
        result["failures"] = comparison["failures"]
        result["status"] = (
            "passed" if comparison["status"] == "passed" else "comparison_failed"
        )
    except (BundleError, ComparisonError) as exc:
        result["failures"] = [str(exc)]
    return result


def _load_suite(
    suite_id: str, suites_path: Path, layouts_path: Path, case_dir: Path
) -> dict[str, Any]:
    schema_path = suites_path.parent / "schemas" / "suites.schema.json"
    document = load_validated_json(suites_path, schema_path, "suite definitions")
    suite = document["suites"].get(suite_id)
    if suite is None:
        available = ", ".join(sorted(document["suites"]))
        raise BundleError(f"unknown suite {suite_id}; available: {available}")

    layouts_schema = layouts_path.parent / "schemas" / "layouts.schema.json"
    layouts = load_validated_json(
        layouts_path, layouts_schema, "layout definitions"
    )["layouts"]
    unknown = [layout for layout in suite["layouts"] if layout not in layouts]
    if unknown:
        raise BundleError(f"suite {suite_id} has unknown layouts: {', '.join(unknown)}")

    case = load_case_definition(suite["case_id"], case_dir)
    if suite["workflow_id"] not in case["workflows"]:
        raise BundleError(
            f"suite {suite_id} refers to unknown workflow {suite['workflow_id']}"
        )
    return suite


def _require_bundle_class(
    bundle_root: Path, case_dir: Path, required: str
) -> None:
    schema = case_dir.parent / "schemas" / "bundle-manifest.schema.json"
    manifest = load_validated_json(
        bundle_root / "manifest.json", schema, "bundle manifest"
    )
    actual = manifest.get("bundle_class", "unspecified")
    if actual != required:
        raise BundleError(
            f"golden-check requires bundle_class={required}; found {actual}"
        )


def _suite_directory(
    settings: dict[str, str], suite_id: str, run_id: str
) -> Path:
    configured = settings.get("MHDG_REGRESSION_RUN_ROOT")
    if not configured:
        raise BundleError("settings must define MHDG_REGRESSION_RUN_ROOT")
    root = Path(configured).expanduser()
    if not root.is_absolute():
        raise BundleError("MHDG_REGRESSION_RUN_ROOT must be an absolute path")
    path = root.resolve() / "suites" / suite_id / run_id
    try:
        path.mkdir(parents=True)
    except FileExistsError as exc:
        raise BundleError(f"suite summary directory already exists: {path}") from exc
    except OSError as exc:
        raise BundleError(f"cannot create suite summary directory {path}: {exc}") from exc
    return path


def _print_summary(summary: dict[str, Any], path: Path) -> None:
    print("layout        run                 comparison  runtime     result")
    for result in summary["results"]:
        runtime = result["duration_seconds"]
        runtime_text = f"{runtime:.3f} s" if runtime is not None else "n/a"
        print(
            f"{result['layout_id']:<13} "
            f"{(result['run_status'] or 'n/a'):<19} "
            f"{(result['comparison_status'] or 'n/a'):<11} "
            f"{runtime_text:<11} {result['status']}"
        )
        for failure in result["failures"][:3]:
            print(f"  FAIL: {failure}")
        if len(result["failures"]) > 3:
            print(f"  ... {len(result['failures']) - 3} more failures")
    print(
        f"suite {summary['status']}: {path} "
        f"({summary['duration_seconds']:.3f} s)"
    )


def _passed(result: dict[str, Any]) -> bool:
    return result["status"] == "passed"


def _write_json(path: Path, document: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        temporary.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)
    except OSError as exc:
        raise BundleError(f"cannot write suite summary {path}: {exc}") from exc


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


def _positive_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("build jobs must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("build jobs must be a positive integer")
    return parsed


if __name__ == "__main__":
    raise SystemExit(main())
