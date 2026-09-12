"""Load staged reference indexes used by scientific comparisons."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

from .documents import load_validated_json
from .support import BundleError
from .files import require_directory, require_file

REFERENCE_MATRIX_ROLE = "reference_matrix"


@dataclass(frozen=True)
class ReferenceMatrix:
    bundle_root: Path
    bundle_id: str
    bundle_version: str
    references: dict[tuple[str, str, str], Path]

    def reference_for(
        self,
        workflow_id: str,
        layout_id: str,
        stage_id: str,
    ) -> Path:
        key = (workflow_id, layout_id, stage_id)
        try:
            return self.references[key]
        except KeyError as exc:
            cell = "/".join(key)
            raise BundleError(
                f"golden matrix has no reference for {cell}"
            ) from exc


def load_reference_matrix(
    bundle_root: Path,
    case: dict[str, Any],
    schema_directory: Path,
) -> ReferenceMatrix | None:
    """Load the optional staged-reference index without rehashing the bundle."""
    bundle_root = require_directory(bundle_root, "golden bundle")
    manifest = load_validated_json(
        bundle_root / "manifest.json",
        schema_directory / "bundle-manifest.schema.json",
        "bundle manifest",
    )
    case_id = case["case_id"]
    if manifest["case_id"] != case_id:
        raise BundleError(f"bundle does not contain case data for {case['case_id']}")

    index_id = manifest["roles"].get(REFERENCE_MATRIX_ROLE)
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
        schema_directory / "reference-matrix.schema.json",
        "reference matrix",
    )
    if matrix["case_id"] != case_id:
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
