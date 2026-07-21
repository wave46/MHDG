"""Load case contracts and derive their required artifact roles."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from bundle.case_definition import normalize_case_definition
from bundle.schemas import load_validated_json
from support.errors import BundleError


def load_case_definition(case_id: str, case_directory: Path) -> dict[str, Any]:
    """Load and validate one tracked case definition by its identifier."""
    if not re.fullmatch(r"[a-z][a-z0-9_]*", case_id):
        raise BundleError(f"invalid case identifier: {case_id}")

    path = case_directory / f"{case_id}.json"
    schema_path = case_directory.parent / "schemas" / "case.schema.json"
    declaration = load_validated_json(
        path,
        schema_path,
        f"case definition {path.name}",
    )
    case = normalize_case_definition(case_id, declaration)
    _validate_case_workflows(case)
    return case


def required_case_roles(
    case: dict[str, Any],
    workflow_id: str | None = None,
) -> set[str]:
    """Return base package roles or the roles required by one workflow."""
    if workflow_id is None:
        return {
            role
            for role, spec in case.get("bundle_files", {}).items()
            if not spec.get("optional", False)
        }
    try:
        workflow = case["workflows"][workflow_id]
    except KeyError as exc:
        raise BundleError(
            f"case {case['case_id']} has no workflow {workflow_id}"
        ) from exc
    return workflow_required_roles(workflow)


def workflow_required_roles(workflow: dict[str, Any]) -> set[str]:
    """Return explicit and staged artifact roles for one workflow."""
    roles = set(workflow.get("required_artifact_roles", []))
    for name in ("mesh_role", "reference_role"):
        if workflow.get(name):
            roles.add(workflow[name])
    for stage in workflow.get("stages", []):
        roles.add(stage["parameter_role"])
        roles.add(stage["transport_configuration_role"])
    return roles


def _validate_case_workflows(case: dict[str, Any]) -> None:
    declared_roles = set(case.get("bundle_files", {}))
    for workflow_id, workflow in case["workflows"].items():
        unknown = (
            sorted(workflow_required_roles(workflow) - declared_roles)
            if declared_roles
            else []
        )
        if unknown:
            raise BundleError(
                f"workflow {workflow_id} uses undefined artifact roles: "
                + ", ".join(unknown)
            )
        if workflow["kind"] not in {
            "staged_fixed_mesh",
            "staged_adaptive_mesh",
        }:
            continue

        missing = [
            name
            for name in ("mesh_role", "reference_role", "stages")
            if not workflow.get(name)
        ]
        if missing:
            raise BundleError(
                f"workflow {workflow_id} is missing: {', '.join(missing)}"
            )
        stages = workflow["stages"]
        stage_ids = [stage["stage_id"] for stage in stages]
        if len(stage_ids) != len(set(stage_ids)):
            raise BundleError(f"workflow {workflow_id} has duplicate stage IDs")
        if stages[0]["restart_from"] != "analytical":
            raise BundleError(
                f"workflow {workflow_id} must start from analytical initialization"
            )
        if any(stage["restart_from"] != "previous_stage" for stage in stages[1:]):
            raise BundleError(
                f"workflow {workflow_id} continuation stages must restart "
                "from the previous stage"
            )
