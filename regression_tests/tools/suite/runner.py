"""Orchestrate workflow-by-layout regression suites."""

from __future__ import annotations

import time
from dataclasses import replace
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
from suite.pairs import compare_layout_pairs
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
    pairs = suite.get("layout_comparisons")
    comparison_mode = _comparison_mode(bool(pairs), compare)
    if resume:
        summary = resume_summary(
            summary_path,
            suite_id,
            selected_run_id,
            suite,
            comparison_mode,
        )
    else:
        summary = new_summary(
            suite_id,
            selected_run_id,
            suite,
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
    cell_inputs = replace(inputs, compare=False) if pairs else inputs
    _run_pending_cells(
        cell_inputs,
        suite_id,
        selected_workflows,
        suite["layouts"],
        summary_path,
        summary,
        resume,
    )
    if pairs and compare:
        summary["comparisons"] = compare_layout_pairs(
            summary,
            case_directory,
            tolerances_path,
        )
    finalize_summary(summary)
    write_json_atomic(summary_path, summary, "suite summary")
    return summary_path, summary


def _comparison_mode(layout_pairs: bool, compare: bool) -> str:
    if layout_pairs:
        return "layout_pairs" if compare else "deferred_layout_pairs"
    return "immediate" if compare else "deferred"


def _run_pending_cells(
    inputs: SuiteRunInputs,
    suite_id: str,
    workflow_ids: list[str],
    layout_ids: list[str],
    summary_path: Path,
    summary: dict[str, Any],
    resume: bool,
) -> None:
    print(f"suite: {suite_id} ({inputs.run_id})")
    recorded = completed_cells(summary)
    started_clock = time.monotonic()
    for layout_id in layout_ids:
        for workflow_id in workflow_ids:
            if (workflow_id, layout_id) in recorded:
                print(f"skipping recorded {workflow_id} / {layout_id}", flush=True)
                continue
            cell_inputs = _resumable_cell_inputs(
                inputs,
                workflow_id,
                layout_id,
                resume,
            )
            print(f"running {workflow_id} / {layout_id} ...", flush=True)
            summary["results"].append(
                run_cell(cell_inputs, workflow_id, layout_id)
            )
            summary["duration_seconds"] += time.monotonic() - started_clock
            started_clock = time.monotonic()
            write_json_atomic(summary_path, summary, "suite summary")
    summary["duration_seconds"] += time.monotonic() - started_clock


def _resumable_cell_inputs(
    inputs: SuiteRunInputs,
    workflow_id: str,
    layout_id: str,
    resume: bool,
) -> SuiteRunInputs:
    """Preserve an unrecorded run and select a new attempt identifier."""
    if not resume:
        return inputs

    run_root = Path(inputs.settings["MHDG_REGRESSION_RUN_ROOT"]).expanduser()
    cell_root = run_root / inputs.case_id / workflow_id / layout_id
    existing = cell_root / inputs.run_id
    if not existing.exists():
        return inputs

    attempt = 1
    while (cell_root / f"{inputs.run_id}-resume-{attempt}").exists():
        attempt += 1
    retry_id = f"{inputs.run_id}-resume-{attempt}"
    print(
        f"preserving unrecorded run {existing}; retrying as {retry_id}",
        flush=True,
    )
    return replace(inputs, run_id=retry_id)
