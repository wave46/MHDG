"""Validate external regression bundles and their tracked case contracts."""

from __future__ import annotations

from pathlib import Path, PurePosixPath
from typing import Any

from bundle.cases import load_case_definitions, required_case_roles
from bundle.models import ValidationSummary
from bundle.schemas import load_validated_json
from bundle.settings import bundle_root_from_settings, read_settings
from support.errors import BundleError, MissingArtifactError
from support.files import sha256_digest


def validate_bundle(settings_path: Path, case_dir: Path) -> ValidationSummary:
    """Validate bundle structure, physical artifacts, and tracked case roles."""
    bundle_root = bundle_root_from_settings(read_settings(settings_path))
    return validate_bundle_root(bundle_root, case_dir)


def validate_bundle_root(bundle_root: Path, case_dir: Path) -> ValidationSummary:
    """Validate a bundle root without requiring a local settings file."""
    try:
        bundle_root = bundle_root.resolve(strict=True)
    except FileNotFoundError as exc:
        raise BundleError(f"bundle root does not exist: {bundle_root}") from exc
    if not bundle_root.is_dir():
        raise BundleError(f"bundle root is not a directory: {bundle_root}")

    schema_dir = case_dir.parent / "schemas"

    manifest = load_validated_json(
        bundle_root / "manifest.json",
        schema_dir / "bundle-manifest.schema.json",
        "bundle manifest",
    )
    cases = load_case_definitions(case_dir)

    available, verified_bytes, warnings = _verify_artifacts(
        bundle_root, manifest["artifacts"]
    )
    case_data = _resolve_case_data(manifest)
    _verify_reference_matrices(
        bundle_root, manifest, case_data, available, schema_dir
    )
    checked_cases = _verify_case_requirements(cases, case_data, available)

    return ValidationSummary(
        bundle_id=manifest["bundle_id"],
        bundle_version=manifest["bundle_version"],
        artifact_count=len(manifest["artifacts"]),
        verified_artifact_count=len(available),
        verified_bytes=verified_bytes,
        case_data=sorted(case_data),
        checked_cases=checked_cases,
        warnings=warnings,
    )


def _verify_artifacts(
    bundle_root: Path, artifacts: dict[str, dict[str, Any]]
) -> tuple[set[str], int, list[str]]:
    available: set[str] = set()
    verified_bytes = 0
    warnings: list[str] = []

    for artifact_id, artifact in artifacts.items():
        label = f"manifest.artifacts.{artifact_id}"
        try:
            path = _artifact_path(bundle_root, artifact["path"], label)
        except MissingArtifactError:
            if not artifact.get("optional", False):
                raise
            warnings.append(f"optional artifact unavailable: {artifact_id}")
            continue

        actual_size = path.stat().st_size
        if actual_size != artifact["size_bytes"]:
            raise BundleError(
                f"{label}.size_bytes is {artifact['size_bytes']}, "
                f"but the file contains {actual_size} bytes"
            )
        if sha256_digest(path) != artifact["sha256"]:
            raise BundleError(f"{label}.sha256 does not match the file")
        available.add(artifact_id)
        verified_bytes += actual_size

    return available, verified_bytes, warnings


def _resolve_case_data(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    artifact_ids = set(manifest["artifacts"])
    resolved = {}
    for data_id, case_data in manifest["case_data"].items():
        for role, artifact_id in case_data["roles"].items():
            if artifact_id not in artifact_ids:
                raise BundleError(
                    f"manifest.case_data.{data_id}.roles.{role} refers to "
                    f"unknown artifact {artifact_id}"
                )
        resolved[data_id] = case_data
    return resolved


def _verify_reference_matrices(
    bundle_root: Path,
    manifest: dict[str, Any],
    case_data_by_id: dict[str, dict[str, Any]],
    available: set[str],
    schema_dir: Path,
) -> None:
    artifacts = manifest["artifacts"]
    for data_id, case_data in case_data_by_id.items():
        index_id = case_data["roles"].get("reference_matrix")
        if index_id is None:
            continue
        index_artifact = artifacts[index_id]
        index_path = _artifact_path(
            bundle_root,
            index_artifact["path"],
            f"manifest.artifacts.{index_id}",
        )
        matrix = load_validated_json(
            index_path,
            schema_dir / "reference-matrix.schema.json",
            f"reference matrix {data_id}",
        )
        if matrix["case_id"] != case_data["case_id"]:
            raise BundleError(f"reference matrix {data_id} has the wrong case_id")

        cells = []
        for reference in matrix["references"]:
            artifact_id = reference["artifact_id"]
            artifact = artifacts.get(artifact_id)
            if artifact is None or artifact_id not in available:
                raise BundleError(
                    f"reference matrix {data_id} uses unavailable artifact "
                    f"{artifact_id}"
                )
            if artifact["media_type"] != "application/x-hdf5":
                raise BundleError(
                    f"reference matrix artifact {artifact_id} is not HDF5"
                )
            cells.append(
                (
                    reference["workflow_id"],
                    reference["layout_id"],
                    reference["stage_id"],
                )
            )
        if len(cells) != len(set(cells)):
            raise BundleError(f"reference matrix {data_id} contains duplicate cells")


def _verify_case_requirements(
    cases: list[dict[str, Any]],
    case_data_by_id: dict[str, dict[str, Any]],
    available: set[str],
) -> list[str]:
    checked = []
    for case in cases:
        case_id = case["case_id"]
        if case_id not in case_data_by_id:
            continue

        case_data = case_data_by_id[case_id]
        if case_data["case_id"] != case_id:
            raise BundleError(
                f"manifest.case_data.{case_id}.case_id is "
                f"{case_data['case_id']}, expected {case_id}"
            )

        required = required_case_roles(case)
        mapping = case_data["roles"]
        missing = sorted(required - mapping.keys())
        if missing:
            raise BundleError(
                f"case {case_id} is missing required artifact roles: {', '.join(missing)}"
            )
        unavailable = sorted(role for role in required if mapping[role] not in available)
        if unavailable:
            raise BundleError(
                f"case {case_id} has unavailable required artifacts for roles: "
                f"{', '.join(unavailable)}"
            )
        checked.append(case_id)
    return checked


def _artifact_path(root: Path, relative_path: str, label: str) -> Path:
    if "\\" in relative_path:
        raise BundleError(f"{label}.path must use POSIX separators")
    posix_path = PurePosixPath(relative_path)
    if (
        posix_path.is_absolute()
        or ".." in posix_path.parts
        or posix_path == PurePosixPath(".")
    ):
        raise BundleError(f"{label}.path must stay relative to the bundle root")

    candidate = root.joinpath(*posix_path.parts)
    resolved = candidate.resolve(strict=False)
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise BundleError(f"{label}.path resolves outside the bundle root") from exc
    if not candidate.exists():
        raise MissingArtifactError(f"{label}.path does not exist: {relative_path}")
    if not resolved.is_file():
        raise BundleError(f"{label}.path is not a regular file: {relative_path}")
    return resolved
