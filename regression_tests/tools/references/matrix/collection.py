"""Validate matrix suite evidence and collect its stage outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from bundle.artifacts import validate_bundle_identity
from references.matrix.models import MatrixRun, StageReference
from support.documents import load_json
from support.errors import BundleError
from support.files import is_within
from support.identifiers import IDENTIFIER_RE
from support.paths import recorded_directory, recorded_file


def validate_matrix_summary(summary: dict[str, Any]) -> None:
    """Require one completed deferred result for every declared matrix cell."""
    workflow_ids = summary.get("workflow_ids")
    layout_ids = summary.get("layout_ids")
    if (
        summary.get("comparison_mode") != "deferred"
        or not _valid_id_list(workflow_ids)
        or not _valid_id_list(layout_ids)
    ):
        raise BundleError(
            "reference matrix requires deferred workflow_ids and layout_ids"
        )

    results = summary.get("results")
    if not isinstance(results, list) or not results:
        raise BundleError("suite summary contains no matrix results")
    cells = []
    for result in results:
        if not isinstance(result, dict) or any(
            result.get(name) != expected
            for name, expected in (
                ("status", "passed"),
                ("run_status", "completed"),
                ("comparison_status", "not_run"),
            )
        ):
            raise BundleError("reference matrix contains an incomplete run")
        cell = (result.get("workflow_id"), result.get("layout_id"))
        if any(
            not isinstance(value, str) or not IDENTIFIER_RE.fullmatch(value)
            for value in cell
        ):
            raise BundleError("reference matrix contains an invalid cell")
        cells.append(cell)

    expected = {
        (workflow_id, layout_id)
        for workflow_id in workflow_ids
        for layout_id in layout_ids
    }
    if len(cells) != len(set(cells)) or set(cells) != expected:
        raise BundleError("reference matrix is incomplete or contains duplicate cells")


def collect_matrix_runs(
    summary: dict[str, Any],
    case: dict[str, Any],
    source_bundle: Path,
    source_manifest: dict[str, Any],
) -> list[MatrixRun]:
    """Validate staged run provenance and select one output per stage."""
    runs = []
    for result in summary["results"]:
        workflow_id = result["workflow_id"]
        layout_id = result["layout_id"]
        workflow = case["workflows"].get(workflow_id)
        if workflow is None or not workflow.get("stages"):
            raise BundleError(f"matrix refers to unsupported workflow {workflow_id}")

        directory = recorded_directory(result.get("run_directory"), "run")
        plan = load_json(directory / "run_plan.json", "run plan")
        metadata = load_json(directory / "run_metadata.json", "run metadata")
        expected = {
            "case_id": summary["case_id"],
            "workflow_id": workflow_id,
            "layout_id": layout_id,
        }
        if any(plan.get(name) != value for name, value in expected.items()):
            raise BundleError(
                f"matrix run has inconsistent identity: {workflow_id}/{layout_id}"
            )
        if metadata.get("status") != "completed":
            raise BundleError(f"matrix run is incomplete: {workflow_id}/{layout_id}")
        validate_bundle_identity(plan, source_bundle, source_manifest)

        expected_stages = [stage["stage_id"] for stage in workflow["stages"]]
        plan_stages = [stage.get("stage_id") for stage in plan.get("stages", [])]
        metadata_stages = metadata.get("stages", [])
        recorded_stages = [stage.get("stage_id") for stage in metadata_stages]
        if plan_stages != expected_stages or recorded_stages != expected_stages:
            raise BundleError(
                f"matrix run stages differ from {workflow_id} definition"
            )

        stages = tuple(
            _stage_reference(stage, workflow_id, layout_id)
            for stage in metadata_stages
        )
        runs.append(MatrixRun(workflow_id, layout_id, directory, stages))
    return runs


def _stage_reference(
    stage: dict[str, Any],
    workflow_id: str,
    layout_id: str,
) -> StageReference:
    if stage.get("status") != "completed":
        raise BundleError(
            f"matrix stage is incomplete: {workflow_id}/{layout_id}/"
            f"{stage.get('stage_id')}"
        )
    stage_directory = recorded_directory(stage.get("run_directory"), "stage run")
    solution = recorded_file(stage.get("selected_hdf5"), "stage solution")
    if not is_within(solution, stage_directory):
        raise BundleError("selected stage solution is outside its run directory")
    return StageReference(stage["stage_id"], solution)


def _valid_id_list(value: Any) -> bool:
    return (
        isinstance(value, list)
        and bool(value)
        and len(value) == len(set(value))
        and all(
            isinstance(item, str) and IDENTIFIER_RE.fullmatch(item)
            for item in value
        )
    )
