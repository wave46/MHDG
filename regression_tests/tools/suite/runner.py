"""Orchestrate workflow-by-layout regression suites."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from bundle.settings import bundle_root_from_settings, read_settings
from bundle.validation import validate_bundle_root
from suite.cells import run_cell
from suite.configuration import (
    load_suite_definition,
    require_bundle_class,
    suite_directory,
)
from suite.models import SuiteRunInputs
from suite.summary import (
    completed_cells,
    finalize_summary,
    new_summary,
    resume_summary,
)
from support.documents import write_json_atomic
from support.errors import BundleError
from support.identifiers import IDENTIFIER_RE
from support.time import utc_run_id


def run_suite(
    settings_path: Path,
    suite_id: str,
    case_directory: Path,
    layouts_path: Path,
    suites_path: Path,
    tolerances_path: Path,
    run_id: str | None = None,
    required_bundle_class: str | None = None,
    compare: bool = True,
    resume: bool = False,
) -> tuple[Path, dict[str, Any]]:
    """Execute every workflow-layout cell and write a reusable summary."""
    suite = load_suite_definition(
        suite_id,
        suites_path,
        layouts_path,
        case_directory,
    )
    settings = read_settings(settings_path)
    bundle_root = bundle_root_from_settings(settings)
    validate_bundle_root(bundle_root, case_directory)
    if required_bundle_class is not None:
        require_bundle_class(bundle_root, case_directory, required_bundle_class)

    selected_run_id = run_id or utc_run_id()
    if not IDENTIFIER_RE.fullmatch(selected_run_id):
        raise BundleError(f"invalid suite run identifier: {selected_run_id}")

    output_directory = suite_directory(
        settings,
        suite_id,
        selected_run_id,
        resume,
    )
    summary_path = output_directory / "suite_summary.json"
    selected_workflows = suite["workflow_ids"]
    comparison_mode = "immediate" if compare else "deferred"
    if resume:
        summary = resume_summary(
            summary_path,
            suite_id,
            selected_run_id,
            selected_workflows,
            suite["layouts"],
            comparison_mode,
        )
    else:
        summary = new_summary(
            suite_id,
            selected_run_id,
            suite,
            selected_workflows,
            comparison_mode,
        )
        write_json_atomic(summary_path, summary, "suite summary")

    inputs = SuiteRunInputs(
        settings_path=settings_path,
        settings=settings,
        case_id=suite["case_id"],
        run_id=selected_run_id,
        case_directory=case_directory,
        layouts_path=layouts_path,
        tolerances_path=tolerances_path,
        compare=compare,
    )
    _run_pending_cells(
        inputs,
        suite_id,
        selected_workflows,
        suite["layouts"],
        summary_path,
        summary,
    )
    finalize_summary(summary)
    write_json_atomic(summary_path, summary, "suite summary")
    return summary_path, summary


def _run_pending_cells(
    inputs: SuiteRunInputs,
    suite_id: str,
    workflow_ids: list[str],
    layout_ids: list[str],
    summary_path: Path,
    summary: dict[str, Any],
) -> None:
    print(f"suite: {suite_id} ({inputs.run_id})")
    recorded = completed_cells(summary)
    started_clock = time.monotonic()
    for layout_id in layout_ids:
        for workflow_id in workflow_ids:
            if (workflow_id, layout_id) in recorded:
                print(f"skipping recorded {workflow_id} / {layout_id}", flush=True)
                continue
            print(f"running {workflow_id} / {layout_id} ...", flush=True)
            summary["results"].append(run_cell(inputs, workflow_id, layout_id))
            summary["duration_seconds"] += time.monotonic() - started_clock
            started_clock = time.monotonic()
            write_json_atomic(summary_path, summary, "suite summary")
    summary["duration_seconds"] += time.monotonic() - started_clock
