"""Install collected matrix references and their provenance."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from bundle.artifacts import register_artifact
from references.matrix.models import (
    REFERENCE_MATRIX_ID,
    REFERENCE_MATRIX_ROLE,
    MatrixRun,
)
from support.documents import write_json_direct
from support.errors import BundleError
from support.files import file_identity
from support.paths import require_file
from support.time import utc_now


def install_reference_matrix(
    staging: Path,
    manifest: dict[str, Any],
    summary_path: Path,
    summary: dict[str, Any],
    runs: list[MatrixRun],
    case: dict[str, Any],
    source_manifest: dict[str, Any],
    canonical_warm_roles: tuple[str, ...] = ("warm_restart", "warm_reference"),
) -> None:
    """Copy selected stage outputs and add one generated bundle role."""
    references_directory = staging / "references/matrix"
    provenance_directory = staging / "provenance/golden_matrix"
    for directory in (references_directory, provenance_directory):
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True)
    for artifact_id in list(manifest["artifacts"]):
        if artifact_id.startswith("golden_matrix_"):
            del manifest["artifacts"][artifact_id]

    entries = []
    for run in runs:
        entries.extend(_install_run(staging, manifest, references_directory, run))
        _install_run_provenance(staging, manifest, provenance_directory, run)

    suite_target = provenance_directory / "suite_summary.json"
    shutil.copy2(summary_path, suite_target)
    register_artifact(
        staging,
        manifest,
        "golden_matrix_suite_summary",
        suite_target,
        "application/json",
    )

    index_path = references_directory / "index.json"
    write_json_direct(
        index_path,
        {
            "schema_version": 2,
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
        staging,
        manifest,
        REFERENCE_MATRIX_ID,
        index_path,
        "application/json",
    )
    if manifest["case_id"] != case["case_id"]:
        raise BundleError("source bundle contains another case")
    manifest["roles"][REFERENCE_MATRIX_ROLE] = REFERENCE_MATRIX_ID
    if summary["suite_id"] == "cold_matrix" and canonical_warm_roles:
        _install_canonical_warm_state(
            staging,
            manifest,
            runs,
            case,
            canonical_warm_roles,
        )


def _install_canonical_warm_state(
    staging: Path,
    manifest: dict[str, Any],
    runs: list[MatrixRun],
    case: dict[str, Any],
    roles: tuple[str, ...],
) -> None:
    """Use the canonical cold final state for warm restart and reference."""
    canonical_layout = case["workflows"]["cold_fixed"]["default_layout"]
    matching = [
        run
        for run in runs
        if run.workflow_id == "cold_fixed" and run.layout_id == canonical_layout
    ]
    if len(matching) != 1:
        raise BundleError(
            "cold matrix has no unique canonical cold_fixed final state"
        )
    solution = matching[0].stages[-1].solution
    for role in roles:
        artifact_id = manifest["roles"].get(role)
        artifact = manifest["artifacts"].get(artifact_id)
        if artifact is None:
            raise BundleError(f"source bundle has no {role} artifact")
        target = staging / artifact["path"]
        shutil.copy2(solution, target)
        artifact.update(file_identity(target))


def _install_run(
    staging: Path,
    manifest: dict[str, Any],
    references_directory: Path,
    run: MatrixRun,
) -> list[dict[str, str]]:
    target_directory = references_directory / run.workflow_id / run.layout_id
    target_directory.mkdir(parents=True)
    entries = []
    for number, stage in enumerate(run.stages, start=1):
        artifact_id = (
            f"golden_matrix_{run.workflow_id}_{run.layout_id}_{stage.stage_id}"
        )
        target = target_directory / f"{number:02d}_{stage.stage_id}.h5"
        shutil.copy2(stage.solution, target)
        register_artifact(
            staging,
            manifest,
            artifact_id,
            target,
            "application/x-hdf5",
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
    provenance_directory: Path,
    run: MatrixRun,
) -> None:
    target_directory = provenance_directory / run.workflow_id / run.layout_id
    target_directory.mkdir(parents=True)
    for filename in ("run_plan.json", "run_metadata.json"):
        target = target_directory / filename
        shutil.copy2(require_file(run.directory / filename, filename), target)
        artifact_id = (
            f"golden_matrix_{run.workflow_id}_{run.layout_id}_"
            f"{Path(filename).stem}"
        )
        register_artifact(
            staging,
            manifest,
            artifact_id,
            target,
            "application/json",
        )
