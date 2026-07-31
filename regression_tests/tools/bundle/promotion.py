"""Promote accepted regression results into complete bundles."""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from bundle.cases import load_case_definition
from bundle.models import ValidationSummary
from bundle.settings import bundle_root_from_settings, read_settings
from bundle.validation import validate_bundle_root
from references.canonical import (
    CanonicalReference,
    collect_canonical_reference,
    install_canonical_reference,
)
from references.matrix.collection import collect_matrix_runs
from references.matrix.installation import install_reference_matrix
from references.matrix.models import MatrixRun
from references.mapped import (
    MappedReferences,
    collect_mapped_references,
    install_mapped_references,
)
from references.summaries import validate_promotion_summary
from support.documents import load_json, write_json_direct
from support.errors import BundleError
from support.files import is_within
from support.paths import require_file
from support.time import utc_now


AcceptedReferences = CanonicalReference | MappedReferences | list[MatrixRun]


def promote_bundle(
    settings_path: Path,
    summary_paths: Path | Sequence[Path],
    output: Path,
    bundle_version: str,
    case_directory: Path,
    *,
    bundle_class: str = "golden",
    matrix_warm_roles: tuple[str, ...] = ("warm_restart", "warm_reference"),
) -> ValidationSummary:
    """Copy a source bundle and install accepted reference results in order."""
    source_bundle = bundle_root_from_settings(read_settings(settings_path))
    validate_bundle_root(source_bundle, case_directory)
    source_manifest = load_json(source_bundle / "manifest.json", "bundle manifest")
    if isinstance(summary_paths, Path):
        summary_paths = [summary_paths]
    if not summary_paths:
        raise BundleError("at least one suite summary is required")
    accepted_summaries = []
    for path in summary_paths:
        path = require_file(path, "suite summary")
        summary = load_json(path, "suite summary")
        promotion_kind = validate_promotion_summary(summary)
        case = load_case_definition(summary["case_id"], case_directory)
        accepted = _collect_references(
            promotion_kind,
            summary,
            case,
            source_bundle,
            source_manifest,
        )
        accepted_summaries.append((path, summary, accepted, case))

    return _publish_bundle(
        source_bundle,
        source_manifest,
        accepted_summaries,
        output,
        bundle_version,
        bundle_class,
        case_directory,
        matrix_warm_roles,
    )


def promote_mapped_bundle(
    settings_path: Path,
    summary_path: Path,
    mappings: list[dict[str, Any]],
    output: Path,
    bundle_version: str,
    case_directory: Path,
    *,
    bundle_class: str = "candidate",
) -> ValidationSummary:
    """Copy a source bundle and install workflow-mapped suite outputs."""
    source_bundle = bundle_root_from_settings(read_settings(settings_path))
    validate_bundle_root(source_bundle, case_directory)
    source_manifest = load_json(source_bundle / "manifest.json", "bundle manifest")
    summary_path = require_file(summary_path, "suite summary")
    summary = load_json(summary_path, "suite summary")
    case = load_case_definition(summary.get("case_id", ""), case_directory)
    accepted = collect_mapped_references(
        summary,
        mappings,
        case,
        source_bundle,
        source_manifest,
    )
    return _publish_bundle(
        source_bundle,
        source_manifest,
        [(summary_path, summary, accepted, case)],
        output,
        bundle_version,
        bundle_class,
        case_directory,
        ("warm_restart", "warm_reference"),
    )


def _publish_bundle(
    source_bundle: Path,
    source_manifest: dict[str, Any],
    accepted_summaries: list[
        tuple[Path, dict[str, Any], AcceptedReferences, dict[str, Any]]
    ],
    output: Path,
    bundle_version: str,
    bundle_class: str,
    case_directory: Path,
    matrix_warm_roles: tuple[str, ...],
) -> ValidationSummary:
    if not bundle_version.strip():
        raise BundleError("bundle version must not be empty")
    if bundle_class not in {"candidate", "golden"}:
        raise BundleError(f"unsupported bundle class: {bundle_class}")

    output = output.expanduser().resolve()
    if output.exists():
        raise BundleError(f"output already exists: {output}")
    if is_within(output, source_bundle):
        raise BundleError("bundle output must be outside the source bundle")

    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{output.name}.",
            dir=output.parent,
        ) as workspace:
            staging = Path(workspace) / "bundle"
            shutil.copytree(source_bundle, staging)
            manifest = load_json(staging / "manifest.json", "bundle manifest")
            for summary_path, summary, accepted, case in accepted_summaries:
                _install_references(
                    staging,
                    manifest,
                    summary_path,
                    summary,
                    accepted,
                    case,
                    source_manifest,
                    matrix_warm_roles,
                )
            manifest["bundle_version"] = bundle_version
            manifest["bundle_class"] = bundle_class
            manifest["created_utc"] = utc_now()
            write_json_direct(staging / "manifest.json", manifest)
            result = validate_bundle_root(staging, case_directory)
            staging.rename(output)
    except OSError as exc:
        raise BundleError(f"cannot create bundle {output}: {exc}") from exc
    return result


def _collect_references(
    promotion_kind: str,
    summary: dict[str, Any],
    case: dict[str, Any],
    source_bundle: Path,
    source_manifest: dict[str, Any],
) -> AcceptedReferences:
    if promotion_kind == "canonical":
        return collect_canonical_reference(
            summary,
            case,
            source_bundle,
            source_manifest,
        )
    return collect_matrix_runs(
        summary,
        case,
        source_bundle,
        source_manifest,
    )


def _install_references(
    staging: Path,
    manifest: dict[str, Any],
    summary_path: Path,
    summary: dict[str, Any],
    accepted: AcceptedReferences,
    case: dict[str, Any],
    source_manifest: dict[str, Any],
    matrix_warm_roles: tuple[str, ...],
) -> None:
    if isinstance(accepted, CanonicalReference):
        install_canonical_reference(
            staging,
            manifest,
            summary_path,
            summary,
            accepted,
            case,
            source_manifest,
        )
        return
    if isinstance(accepted, MappedReferences):
        install_mapped_references(staging, manifest, accepted, case)
        return
    install_reference_matrix(
        staging,
        manifest,
        summary_path,
        summary,
        accepted,
        case,
        source_manifest,
        matrix_warm_roles,
    )
