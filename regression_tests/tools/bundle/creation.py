"""Create validated regression bundles from prepared case files."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition, workflow_required_roles
from bundle.models import ValidationSummary
from bundle.validation import validate_bundle_root
from support.documents import write_json_direct
from support.errors import BundleError
from support.files import file_identity, is_within
from support.paths import require_directory
from support.time import utc_now


def create_bundle(
    case_id: str,
    source: Path,
    output: Path,
    case_dir: Path,
    bundle_version: str = "1.0.0",
) -> ValidationSummary:
    """Copy prepared files into a new bundle and validate it before publication."""
    case = load_case_definition(case_id, case_dir)
    file_contract = case.get("bundle_files")
    if not file_contract:
        raise BundleError(f"case {case_id} does not define bundle_files")
    if not bundle_version:
        raise BundleError("bundle version must not be empty")

    source = require_directory(source, "source")
    output = output.expanduser().resolve()
    if output.exists():
        raise BundleError(f"output already exists: {output}")
    if is_within(output, source):
        raise BundleError("output must be outside the prepared source directory")

    workflow_roles = {
        role
        for workflow in case["workflows"].values()
        for role in workflow_required_roles(workflow)
    }
    undefined_roles = sorted(workflow_roles - file_contract.keys())
    if undefined_roles:
        raise BundleError(
            "required roles have no bundle file definition: "
            + ", ".join(undefined_roles)
        )

    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{output.name}.", dir=output.parent
        ) as workspace:
            staging = Path(workspace) / "bundle"
            staging.mkdir()
            manifest = _populate_bundle(staging, source, case, bundle_version)
            write_json_direct(staging / "manifest.json", manifest)
            summary = validate_bundle_root(staging, case_dir)
            staging.rename(output)
    except OSError as exc:
        raise BundleError(f"cannot create bundle {output}: {exc}") from exc
    return summary


def _populate_bundle(
    staging: Path,
    source: Path,
    case: dict[str, Any],
    bundle_version: str,
) -> dict[str, Any]:
    data_id = case["external_data_id"]
    target_directory = staging / "case_data" / data_id
    artifacts: dict[str, dict[str, Any]] = {}
    roles: dict[str, str] = {}

    for role, file_spec in case["bundle_files"].items():
        source_path = source / file_spec["filename"]
        if not source_path.exists():
            if file_spec.get("optional", False):
                continue
            raise BundleError(
                f"required file missing for role {role}: {source_path.name}"
            )
        if not source_path.is_file():
            raise BundleError(f"bundle source is not a file: {source_path}")

        artifact_id = file_spec["artifact_id"]
        if artifact_id in artifacts:
            raise BundleError(f"duplicate bundle artifact_id: {artifact_id}")

        target_path = target_directory / source_path.name
        if target_path.exists():
            raise BundleError(f"duplicate bundle filename: {source_path.name}")
        target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, target_path)
        relative_path = target_path.relative_to(staging).as_posix()
        artifacts[artifact_id] = {
            "path": relative_path,
            **file_identity(target_path),
            "media_type": file_spec["media_type"],
        }
        if "description" in file_spec:
            artifacts[artifact_id]["description"] = file_spec["description"]
        roles[role] = artifact_id

    return {
        "schema_version": 1,
        "bundle_id": f"{data_id}_bundle",
        "bundle_version": bundle_version,
        "bundle_class": "candidate",
        "created_utc": utc_now(),
        "artifacts": artifacts,
        "case_data": {
            data_id: {
                "case_id": case["case_id"],
                "description": case["description"],
                "roles": roles,
            }
        },
    }
