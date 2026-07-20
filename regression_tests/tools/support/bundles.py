"""Shared operations for creating and publishing regression bundles."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from support.errors import BundleError
from support.files import sha256_digest
from support.paths import recorded_directory


def validate_bundle_identity(
    plan: dict[str, Any], source_bundle: Path, manifest: dict[str, Any]
) -> None:
    """Confirm that a recorded run used the expected source bundle."""
    bundle = plan.get("bundle")
    if not isinstance(bundle, dict) or not isinstance(bundle.get("root"), str):
        raise BundleError("run plan has no bundle identity")
    if recorded_directory(bundle["root"], "run bundle") != source_bundle:
        raise BundleError("run did not use the configured source bundle")
    if (
        bundle.get("bundle_id") != manifest.get("bundle_id")
        or bundle.get("bundle_version") != manifest.get("bundle_version")
    ):
        raise BundleError("run and source bundle identities differ")


def register_artifact(
    staging: Path,
    manifest: dict[str, Any],
    artifact_id: str,
    path: Path,
    media_type: str,
) -> None:
    """Register a staged file in a bundle manifest."""
    manifest["artifacts"][artifact_id] = {
        "path": path.relative_to(staging).as_posix(),
        **file_identity(path),
        "media_type": media_type,
    }


def file_identity(path: Path) -> dict[str, Any]:
    """Return the checksum identity used by bundle manifests."""
    try:
        digest = sha256_digest(path)
        return {"size_bytes": path.stat().st_size, "sha256": digest}
    except OSError as exc:
        raise BundleError(f"cannot checksum {path}: {exc}") from exc


def load_bundle_json(path: Path, label: str) -> dict[str, Any]:
    """Load a JSON object using the existing bundle-error contract."""
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BundleError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise BundleError(f"{label} must contain a JSON object")
    return document
