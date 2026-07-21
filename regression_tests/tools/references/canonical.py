"""Collect and install one canonical warm golden reference."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bundle.artifacts import register_artifact, validate_bundle_identity
from comparison.shared.outputs import select_candidate
from support.documents import load_json, write_json_direct
from support.errors import BundleError
from support.files import file_identity
from support.paths import recorded_directory, recorded_file, require_file
from support.time import utc_now


PROVENANCE_FILES = {
    "suite_summary.json": ("golden_reference_suite_summary", "application/json"),
    "run_plan.json": ("golden_reference_run_plan", "application/json"),
    "run_metadata.json": ("golden_reference_run_metadata", "application/json"),
    "comparison.json": ("golden_reference_comparison", "application/json"),
    "stdout.log": ("golden_reference_stdout", "text/plain"),
    "stderr.log": ("golden_reference_stderr", "text/plain"),
}


@dataclass(frozen=True)
class CanonicalRun:
    directory: Path
    plan: dict[str, Any]
    metadata: dict[str, Any]
    comparison: dict[str, Any]
    solution: Path


@dataclass(frozen=True)
class CanonicalReference:
    run: CanonicalRun
    artifact_id: str


def collect_canonical_reference(
    summary: dict[str, Any],
    case: dict[str, Any],
    source_bundle: Path,
    source_manifest: dict[str, Any],
) -> CanonicalReference:
    """Validate the canonical result and its source-reference identity."""
    workflow = case["workflows"].get(summary["workflow_id"])
    if workflow is None:
        raise BundleError("suite summary refers to an unknown workflow")
    run = _canonical_run(summary, workflow["default_layout"])
    artifact_id, source_reference = _source_reference(
        source_bundle,
        source_manifest,
        case,
    )
    _validate_run_sources(
        run,
        source_bundle,
        source_manifest,
        source_reference,
    )
    return CanonicalReference(run, artifact_id)


def install_canonical_reference(
    staging: Path,
    manifest: dict[str, Any],
    summary_path: Path,
    summary: dict[str, Any],
    reference: CanonicalReference,
    case: dict[str, Any],
    source_manifest: dict[str, Any],
) -> None:
    """Replace the warm reference and record its accepted evidence."""
    _install_solution(
        staging,
        manifest,
        reference.artifact_id,
        reference.run.solution,
    )
    _install_provenance(
        staging,
        manifest,
        summary_path,
        summary,
        reference,
        case,
        source_manifest,
    )


def _canonical_run(summary: dict[str, Any], layout_id: str) -> CanonicalRun:
    matching = [
        result for result in summary["results"] if result["layout_id"] == layout_id
    ]
    if len(matching) != 1:
        raise BundleError(f"suite must contain canonical layout {layout_id} once")

    directory = recorded_directory(matching[0].get("run_directory"), "run")
    plan = load_json(directory / "run_plan.json", "run plan")
    metadata = load_json(directory / "run_metadata.json", "run metadata")
    comparison = load_json(directory / "comparison.json", "comparison")
    expected = {
        "case_id": summary["case_id"],
        "workflow_id": summary["workflow_id"],
        "layout_id": layout_id,
    }
    for name, value in expected.items():
        if plan.get(name) != value or comparison.get(name) != value:
            raise BundleError(f"canonical run has inconsistent {name}")
    if metadata.get("status") != "completed" or comparison.get("status") != "passed":
        raise BundleError("canonical run evidence is not fully passing")

    solution = select_candidate(directory, metadata)
    reported = recorded_file(comparison.get("candidate"), "candidate")
    if solution != reported:
        raise BundleError("selected and compared canonical solutions differ")
    if not _same_file(
        file_identity(solution),
        comparison.get("files", {}).get("candidate"),
    ):
        raise BundleError("canonical solution changed after comparison")
    return CanonicalRun(directory, plan, metadata, comparison, solution)


def _source_reference(
    bundle: Path,
    manifest: dict[str, Any],
    case: dict[str, Any],
) -> tuple[str, Path]:
    data_id = case["external_data_id"]
    try:
        artifact_id = manifest["case_data"][data_id]["roles"]["warm_reference"]
        relative_path = manifest["artifacts"][artifact_id]["path"]
    except KeyError as exc:
        raise BundleError("source bundle has no warm reference artifact") from exc
    return artifact_id, require_file(bundle / relative_path, "warm reference")


def _validate_run_sources(
    run: CanonicalRun,
    source_bundle: Path,
    manifest: dict[str, Any],
    source_reference: Path,
) -> None:
    validate_bundle_identity(run.plan, source_bundle, manifest)
    reported = recorded_file(run.comparison.get("reference"), "reference")
    if reported != source_reference:
        raise BundleError("canonical run compared against another reference")
    if not _same_file(
        file_identity(source_reference),
        run.comparison.get("files", {}).get("reference"),
    ):
        raise BundleError("source reference changed after comparison")


def _install_solution(
    staging: Path,
    manifest: dict[str, Any],
    artifact_id: str,
    solution: Path,
) -> None:
    artifact = manifest["artifacts"][artifact_id]
    target = staging / artifact["path"]
    shutil.copy2(solution, target)
    artifact.update(file_identity(target))
    artifact["path"] = target.relative_to(staging).as_posix()


def _install_provenance(
    staging: Path,
    manifest: dict[str, Any],
    summary_path: Path,
    summary: dict[str, Any],
    reference: CanonicalReference,
    case: dict[str, Any],
    source_manifest: dict[str, Any],
) -> None:
    directory = staging / "provenance/golden_reference"
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    for artifact_id in list(manifest["artifacts"]):
        if artifact_id.startswith("golden_reference_"):
            del manifest["artifacts"][artifact_id]

    run = reference.run
    record = {
        "schema_version": 1,
        "status": "golden",
        "created_utc": utc_now(),
        "case_id": summary["case_id"],
        "workflow_id": summary["workflow_id"],
        "canonical_layout": run.plan["layout_id"],
        "suite_id": summary["suite_id"],
        "suite_run_id": summary["run_id"],
        "reference_artifact": reference.artifact_id,
        "source_bundle": {
            "bundle_id": source_manifest["bundle_id"],
            "bundle_version": source_manifest["bundle_version"],
        },
        "tracked_reference": {
            "branch": case["reference_branch"],
            "revision": case["reference_revision"],
        },
        "solver": run.metadata.get("solver", {}),
        "executable": run.metadata.get("executable", {}),
    }
    record_path = directory / "golden_reference.json"
    write_json_direct(record_path, record)
    register_artifact(
        staging,
        manifest,
        "golden_reference_record",
        record_path,
        "application/json",
    )

    sources = {
        "suite_summary.json": summary_path,
        **{
            filename: run.directory / filename
            for filename in PROVENANCE_FILES
            if filename != "suite_summary.json"
        },
    }
    for filename, source in sources.items():
        artifact_id, media_type = PROVENANCE_FILES[filename]
        target = directory / filename
        shutil.copy2(require_file(source, filename), target)
        register_artifact(staging, manifest, artifact_id, target, media_type)


def _same_file(first: dict[str, Any], second: Any) -> bool:
    return isinstance(second, dict) and (
        first["size_bytes"] == second.get("size_bytes")
        and first["sha256"] == second.get("sha256")
    )
