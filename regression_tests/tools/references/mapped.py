"""Collect suite outputs and install them into declared workflow roles."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bundle.artifacts import validate_bundle_identity
from comparison.shared.outputs import select_candidate
from support.documents import load_json
from support.errors import BundleError
from support.files import file_identity, is_within
from support.paths import recorded_directory


@dataclass(frozen=True)
class MappedReference:
    role: str
    solution: Path


@dataclass(frozen=True)
class MappedReferences:
    items: tuple[MappedReference, ...]


def collect_mapped_references(
    summary: dict[str, Any],
    mappings: list[dict[str, Any]],
    case: dict[str, Any],
    source_bundle: Path,
    source_manifest: dict[str, Any],
) -> MappedReferences:
    """Select completed suite outputs using workflow-declared roles."""
    if (
        summary.get("schema_version") != 2
        or summary.get("status") != "passed"
        or summary.get("case_id") != case["case_id"]
    ):
        raise BundleError("mapped references require a passing version-2 suite")

    references = []
    for mapping in mappings:
        workflow_id = mapping["workflow"]
        workflow = case["workflows"].get(workflow_id)
        if workflow is None:
            raise BundleError(f"mapped reference uses unknown workflow {workflow_id}")
        layout_id = workflow.get("default_layout")
        allowed_roles = {
            role
            for name in ("restart_role", "reference_role")
            if (role := workflow.get(name))
        } | set(workflow.get("output_roles", []))
        unknown = sorted(set(mapping["roles"]) - allowed_roles)
        if unknown:
            raise BundleError(
                f"workflow {workflow_id} does not declare roles: "
                + ", ".join(unknown)
            )
        result = _selected_result(summary, workflow_id, layout_id)
        run = recorded_directory(result.get("run_directory"), "mapped run")
        plan = load_json(run / "run_plan.json", "run plan")
        metadata = load_json(run / "run_metadata.json", "run metadata")
        expected = {
            "case_id": case["case_id"],
            "workflow_id": workflow_id,
            "layout_id": layout_id,
        }
        if any(plan.get(name) != value for name, value in expected.items()):
            raise BundleError(f"mapped run has inconsistent identity: {workflow_id}")
        if metadata.get("status") != "completed":
            raise BundleError(f"mapped run is incomplete: {workflow_id}")
        validate_bundle_identity(plan, source_bundle, source_manifest)
        solution = select_candidate(run, metadata)
        if not is_within(solution, run):
            raise BundleError("mapped solution is outside its run directory")
        references.extend(
            MappedReference(role, solution) for role in mapping["roles"]
        )
    return MappedReferences(tuple(references))


def install_mapped_references(
    staging: Path,
    manifest: dict[str, Any],
    references: MappedReferences,
    case: dict[str, Any],
) -> None:
    """Replace only the declared bundle roles."""
    for reference in references.items:
        specification = case["bundle_files"][reference.role]
        artifact_id = manifest["roles"].get(
            reference.role,
            specification["artifact_id"],
        )
        existing = manifest["artifacts"].get(artifact_id)
        target = (
            staging / existing["path"]
            if existing is not None
            else staging / "inputs" / specification["filename"]
        )
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(reference.solution, target)
        artifact = {
            "path": target.relative_to(staging).as_posix(),
            **file_identity(target),
            "media_type": specification["media_type"],
        }
        if specification["optional"]:
            artifact["optional"] = True
        manifest["roles"][reference.role] = artifact_id
        manifest["artifacts"][artifact_id] = artifact


def _selected_result(
    summary: dict[str, Any],
    workflow_id: str,
    layout_id: str,
) -> dict[str, Any]:
    matching = [
        result
        for result in summary.get("results", [])
        if result.get("workflow_id") == workflow_id
        and result.get("layout_id") == layout_id
        and result.get("status") == "passed"
        and result.get("run_status") == "completed"
    ]
    if len(matching) != 1:
        raise BundleError(
            f"suite has no unique completed result for {workflow_id}/{layout_id}"
        )
    return matching[0]
