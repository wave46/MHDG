"""Collect staged workflow outputs into a portable golden-reference matrix."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

from check_bundle import load_validated_json
from support.bundles import (
    load_bundle_json,
    register_artifact,
    validate_bundle_identity,
)
from support.documents import write_json_direct
from support.errors import BundleError
from support.files import is_within
from support.identifiers import IDENTIFIER_RE
from support.paths import (
    recorded_directory,
    recorded_file,
    require_directory,
    require_file,
)
from support.time import utc_now


REFERENCE_MATRIX_ROLE = "reference_matrix"
REFERENCE_MATRIX_ID = "golden_matrix_index"


@dataclass(frozen=True)
class StageReference:
    stage_id: str
    solution: Path


@dataclass(frozen=True)
class MatrixRun:
    workflow_id: str
    layout_id: str
    directory: Path
    stages: tuple[StageReference, ...]


@dataclass(frozen=True)
class ReferenceMatrix:
    bundle_root: Path
    bundle_id: str
    bundle_version: str
    references: dict[tuple[str, str, str], Path]

    def reference_for(
        self, workflow_id: str, layout_id: str, stage_id: str
    ) -> Path:
        key = (workflow_id, layout_id, stage_id)
        try:
            return self.references[key]
        except KeyError as exc:
            cell = "/".join(key)
            raise BundleError(f"golden matrix has no reference for {cell}") from exc


def load_reference_matrix(
    bundle_root: Path,
    case: dict[str, Any],
    schema_dir: Path,
) -> ReferenceMatrix | None:
    """Load the optional staged-reference index without rehashing the bundle."""
    bundle_root = require_directory(bundle_root, "golden bundle")
    manifest = load_validated_json(
        bundle_root / "manifest.json",
        schema_dir / "bundle-manifest.schema.json",
        "bundle manifest",
    )
    data_id = case["external_data_id"]
    case_data = manifest["case_data"].get(data_id)
    valid_case_ids = {case["case_id"], data_id}
    if case_data is None or case_data["case_id"] not in valid_case_ids:
        raise BundleError(f"bundle does not contain case data for {case['case_id']}")

    index_id = case_data["roles"].get(REFERENCE_MATRIX_ROLE)
    if index_id is None:
        return None
    index_path = _manifest_artifact_path(
        bundle_root,
        manifest,
        index_id,
        "reference matrix",
        expected_media_type="application/json",
    )
    matrix = load_validated_json(
        index_path,
        schema_dir / "reference-matrix.schema.json",
        "reference matrix",
    )
    if matrix["case_id"] != case_data["case_id"]:
        raise BundleError("reference matrix has the wrong case_id")

    references: dict[tuple[str, str, str], Path] = {}
    for entry in matrix["references"]:
        key = (entry["workflow_id"], entry["layout_id"], entry["stage_id"])
        if key in references:
            cell = "/".join(key)
            raise BundleError(f"reference matrix contains duplicate cell {cell}")
        references[key] = _manifest_artifact_path(
            bundle_root,
            manifest,
            entry["artifact_id"],
            f"reference {'/'.join(key)}",
            expected_media_type="application/x-hdf5",
        )
    return ReferenceMatrix(
        bundle_root,
        manifest["bundle_id"],
        manifest["bundle_version"],
        references,
    )


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
        plan = load_bundle_json(directory / "run_plan.json", "run plan")
        metadata = load_bundle_json(directory / "run_metadata.json", "run metadata")
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


def install_reference_matrix(
    staging: Path,
    manifest: dict[str, Any],
    summary_path: Path,
    summary: dict[str, Any],
    runs: list[MatrixRun],
    case: dict[str, Any],
    source_manifest: dict[str, Any],
) -> None:
    """Copy selected stage outputs and add one generated bundle role."""
    references_dir = staging / "references/matrix"
    provenance_dir = staging / "provenance/golden_matrix"
    for directory in (references_dir, provenance_dir):
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True)
    for artifact_id in list(manifest["artifacts"]):
        if artifact_id.startswith("golden_matrix_"):
            del manifest["artifacts"][artifact_id]

    entries = []
    for run in runs:
        entries.extend(_install_run(staging, manifest, references_dir, run))
        _install_run_provenance(staging, manifest, provenance_dir, run)

    suite_target = provenance_dir / "suite_summary.json"
    shutil.copy2(summary_path, suite_target)
    register_artifact(
        staging,
        manifest,
        "golden_matrix_suite_summary",
        suite_target,
        "application/json",
    )

    index_path = references_dir / "index.json"
    write_json_direct(
        index_path,
        {
            "schema_version": 1,
            "created_utc": utc_now(),
            "case_id": summary["case_id"],
            "suite_id": summary["suite_id"],
            "suite_run_id": summary["run_id"],
            "source_bundle": {
                "bundle_id": source_manifest["bundle_id"],
                "bundle_version": source_manifest["bundle_version"],
            },
            "tracked_reference": {
                "branch": case["reference_branch"],
                "revision": case["reference_revision"],
            },
            "references": entries,
        },
    )
    register_artifact(
        staging, manifest, REFERENCE_MATRIX_ID, index_path, "application/json"
    )
    try:
        case_data = manifest["case_data"][case["external_data_id"]]
    except KeyError as exc:
        raise BundleError("source bundle has no matching case-data entry") from exc
    case_data["case_id"] = case["case_id"]
    roles = case_data["roles"]
    roles[REFERENCE_MATRIX_ROLE] = REFERENCE_MATRIX_ID


def _stage_reference(
    stage: dict[str, Any], workflow_id: str, layout_id: str
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


def _install_run(
    staging: Path,
    manifest: dict[str, Any],
    references_dir: Path,
    run: MatrixRun,
) -> list[dict[str, str]]:
    target_dir = references_dir / run.workflow_id / run.layout_id
    target_dir.mkdir(parents=True)
    entries = []
    for number, stage in enumerate(run.stages, start=1):
        artifact_id = (
            f"golden_matrix_{run.workflow_id}_{run.layout_id}_{stage.stage_id}"
        )
        target = target_dir / f"{number:02d}_{stage.stage_id}.h5"
        shutil.copy2(stage.solution, target)
        register_artifact(
            staging, manifest, artifact_id, target, "application/x-hdf5"
        )
        entries.append(
            {
                "workflow_id": run.workflow_id,
                "layout_id": run.layout_id,
                "stage_id": stage.stage_id,
                "artifact_id": artifact_id,
            }
        )
    return entries


def _install_run_provenance(
    staging: Path,
    manifest: dict[str, Any],
    provenance_dir: Path,
    run: MatrixRun,
) -> None:
    target_dir = provenance_dir / run.workflow_id / run.layout_id
    target_dir.mkdir(parents=True)
    for filename in ("run_plan.json", "run_metadata.json"):
        target = target_dir / filename
        shutil.copy2(require_file(run.directory / filename, filename), target)
        artifact_id = (
            f"golden_matrix_{run.workflow_id}_{run.layout_id}_"
            f"{Path(filename).stem}"
        )
        register_artifact(
            staging, manifest, artifact_id, target, "application/json"
        )


def _manifest_artifact_path(
    bundle_root: Path,
    manifest: dict[str, Any],
    artifact_id: str,
    label: str,
    expected_media_type: str | None = None,
) -> Path:
    artifact = manifest["artifacts"].get(artifact_id)
    if artifact is None:
        raise BundleError(f"{label} refers to unknown artifact {artifact_id}")
    if expected_media_type and artifact["media_type"] != expected_media_type:
        raise BundleError(f"{label} artifact is not {expected_media_type}")

    relative = PurePosixPath(artifact["path"])
    candidate = bundle_root.joinpath(*relative.parts).resolve()
    try:
        candidate.relative_to(bundle_root)
    except ValueError as exc:
        raise BundleError(f"{label} resolves outside the golden bundle") from exc
    return require_file(candidate, label)


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
