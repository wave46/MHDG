"""Compile concise case catalogs into the harness's normalized contract."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

from support.errors import BundleError


MEDIA_TYPES = {
    ".h5": "application/x-hdf5",
    ".json": "application/json",
    ".msh": "application/x-gmsh",
    ".geo": "text/plain",
    ".nml": "text/plain",
    ".txt": "text/plain",
}
COMPARISON_KEYS = {
    "method": "comparison_policy",
    "profile": "tolerance_profile",
    "cross_layout_profile": "cross_layout_tolerance_profile",
    "stage_profile": "stage_tolerance_profile",
}


def normalize_case_definition(
    case_id: str,
    declaration: dict[str, Any],
) -> dict[str, Any]:
    """Derive machine-oriented fields from one readable case declaration."""
    reference = declaration["reference"]
    return {
        "schema_version": declaration["schema_version"],
        "case_id": case_id,
        "description": declaration["description"],
        "reference_branch": reference["branch"],
        "reference_revision": reference["revision"],
        "bundle_files": _bundle_files(case_id, declaration.get("files", {})),
        "workflows": _workflows(declaration["workflows"]),
    }


def _bundle_files(
    case_id: str,
    files: dict[str, dict[str, str]],
) -> dict[str, dict[str, Any]]:
    normalized = {}
    for optional, section in (
        (False, files.get("required", {})),
        (True, files.get("optional", {})),
    ):
        for role, filename in section.items():
            if role in normalized:
                raise BundleError(f"case file role is declared twice: {role}")
            suffix = Path(filename).suffix.lower()
            try:
                media_type = MEDIA_TYPES[suffix]
            except KeyError as exc:
                raise BundleError(
                    f"cannot infer media type for case file: {filename}"
                ) from exc
            normalized[role] = {
                "filename": filename,
                "artifact_id": f"{case_id}_{role}",
                "media_type": media_type,
                "optional": optional,
            }
    return normalized


def _workflows(
    declarations: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    resolved: dict[str, dict[str, Any]] = {}
    active: set[str] = set()

    def resolve(workflow_id: str) -> dict[str, Any]:
        if workflow_id in resolved:
            return resolved[workflow_id]
        if workflow_id in active:
            raise BundleError(f"workflow inheritance cycle at {workflow_id}")
        active.add(workflow_id)
        declaration = declarations[workflow_id]
        parent_id = declaration.get("extends")
        if parent_id:
            if parent_id not in declarations:
                raise BundleError(
                    f"workflow {workflow_id} extends unknown workflow {parent_id}"
                )
            changes = deepcopy(
                {
                    key: value
                    for key, value in declaration.items()
                    if key != "extends"
                }
            )
            merged = _extend_workflow(resolve(parent_id), changes)
        else:
            merged = deepcopy(declaration)
        active.remove(workflow_id)
        resolved[workflow_id] = merged
        return merged

    return {
        workflow_id: _normalize_workflow(resolve(workflow_id))
        for workflow_id in declarations
    }


def _extend_workflow(
    base: dict[str, Any],
    changes: dict[str, Any],
) -> dict[str, Any]:
    """Apply one derived declaration while retaining base parameter values."""
    extended = deepcopy(base)
    extended.update(changes)
    if "parameter_overrides" in changes:
        extended["parameter_overrides"] = {
            **base.get("parameter_overrides", {}),
            **changes["parameter_overrides"],
        }
    return extended


def _normalize_workflow(declaration: dict[str, Any]) -> dict[str, Any]:
    workflow = {
        "kind": declaration["type"],
        "description": declaration["description"],
        "required_artifact_roles": declaration.get("inputs", []),
    }
    _copy_if_present(declaration, workflow, "layout", "default_layout")
    _copy_if_present(declaration, workflow, "mesh", "mesh_role")
    _copy_if_present(declaration, workflow, "reference", "reference_role")
    _copy_if_present(
        declaration,
        workflow,
        "parameter_overrides",
        "parameter_overrides",
    )
    for source, target in COMPARISON_KEYS.items():
        _copy_if_present(declaration.get("comparison", {}), workflow, source, target)

    stages = declaration.get("stages")
    if stages:
        adaptive_stages = set(declaration.get("adaptive_stages", []))
        stage_ids = {stage["id"] for stage in stages}
        unknown = sorted(adaptive_stages - stage_ids)
        if unknown:
            raise BundleError(
                "adaptive_stages contains unknown stages: " + ", ".join(unknown)
            )
        workflow["stages"] = [
            _normalize_stage(stage, index, adaptive_stages)
            for index, stage in enumerate(stages)
        ]
    return workflow


def _normalize_stage(
    declaration: dict[str, Any],
    index: int,
    adaptive_stages: set[str],
) -> dict[str, Any]:
    stage_id = declaration["id"]
    stage = {
        "stage_id": stage_id,
        "parameter_role": declaration["parameters"],
        "transport_configuration_role": declaration["transport"],
        "restart_from": "analytical" if index == 0 else "previous_stage",
        "newton_check": declaration.get("newton_check", "bounded"),
    }
    overrides = declaration.get("parameter_overrides", {}).copy()
    if stage_id in adaptive_stages:
        overrides["rest_adapt"] = True
    if overrides:
        stage["parameter_overrides"] = overrides
    return stage


def _copy_if_present(
    source: dict[str, Any],
    target: dict[str, Any],
    source_key: str,
    target_key: str,
) -> None:
    if source_key in source:
        target[target_key] = source[source_key]
