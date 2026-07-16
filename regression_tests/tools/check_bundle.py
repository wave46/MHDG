#!/usr/bin/env python3
"""Read-only validation for an external MHDG regression data bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

try:
    from jsonschema import Draft202012Validator, FormatChecker
    from jsonschema.exceptions import SchemaError
except ImportError as exc:  # pragma: no cover - depends on the local environment
    raise SystemExit(
        "jsonschema is required; install regression_tests/requirements.txt"
    ) from exc


SETTING_RE = re.compile(r"^[A-Z][A-Z0-9_]*$")


class BundleError(ValueError):
    """Raised when settings or bundle data violate the regression contract."""


class MissingArtifactError(BundleError):
    """Raised when a declared artifact path does not exist."""


@dataclass(frozen=True)
class ValidationSummary:
    bundle_id: str
    bundle_version: str
    artifact_count: int
    verified_artifact_count: int
    verified_bytes: int
    case_data: list[str]
    checked_cases: list[str]
    warnings: list[str]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path)
    args = parser.parse_args(argv)

    try:
        summary = validate_bundle(args.settings, args.cases)
    except BundleError as exc:
        print(f"bundle validation failed: {exc}", file=sys.stderr)
        return 1

    print(f"bundle valid: {summary.bundle_id} version {summary.bundle_version}")
    print(
        f"verified {summary.verified_artifact_count} of "
        f"{summary.artifact_count} artifacts ({summary.verified_bytes} bytes)"
    )
    print(f"case data: {', '.join(summary.case_data)}")
    checked = ", ".join(summary.checked_cases) or "none from this checkout"
    print(f"required roles checked: {checked}")
    for warning in summary.warnings:
        print(f"warning: {warning}", file=sys.stderr)
    return 0


def validate_bundle(settings_path: Path, case_dir: Path) -> ValidationSummary:
    """Validate bundle structure, physical artifacts, and tracked case roles."""
    bundle_root = _bundle_root(_read_settings(settings_path))
    schema_dir = case_dir.parent / "schemas"

    manifest = _load_json(bundle_root / "manifest.json", "bundle manifest")
    _validate_json(
        manifest,
        schema_dir / "bundle-manifest.schema.json",
        "bundle manifest",
    )
    cases = _load_cases(case_dir, schema_dir / "case.schema.json")

    available, verified_bytes, warnings = _verify_artifacts(
        bundle_root, manifest["artifacts"]
    )
    roles_by_case = _resolve_roles(manifest)
    checked_cases = _verify_case_requirements(cases, roles_by_case, available)

    return ValidationSummary(
        bundle_id=manifest["bundle_id"],
        bundle_version=manifest["bundle_version"],
        artifact_count=len(manifest["artifacts"]),
        verified_artifact_count=len(available),
        verified_bytes=verified_bytes,
        case_data=sorted(roles_by_case),
        checked_cases=checked_cases,
        warnings=warnings,
    )


def _load_cases(case_dir: Path, schema_path: Path) -> list[dict[str, Any]]:
    if not case_dir.is_dir():
        raise BundleError(f"tracked case directory does not exist: {case_dir}")
    cases = []
    for path in sorted(case_dir.glob("*.json")):
        case = _load_json(path, f"case definition {path.name}")
        _validate_json(case, schema_path, f"case definition {path.name}")
        cases.append(case)
    return cases


def _validate_json(document: Any, schema_path: Path, label: str) -> None:
    schema = _load_json(schema_path, f"schema {schema_path.name}")
    try:
        Draft202012Validator.check_schema(schema)
    except SchemaError as exc:
        raise BundleError(f"invalid regression schema {schema_path}: {exc.message}") from exc

    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    errors = sorted(
        validator.iter_errors(document),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if not errors:
        return
    error = errors[0]
    location = ".".join(str(part) for part in error.absolute_path)
    where = f".{location}" if location else ""
    raise BundleError(f"{label}{where}: {error.message}")


def _read_settings(settings_path: Path) -> dict[str, str]:
    """Parse KEY=VALUE settings without executing shell code."""
    try:
        lines = settings_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise BundleError(f"cannot read settings file {settings_path}: {exc}") from exc

    settings: dict[str, str] = {}
    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            raise BundleError(
                f"{settings_path}:{line_number}: expected a KEY=VALUE setting"
            )
        key, value = (part.strip() for part in line.split("=", 1))
        if not SETTING_RE.fullmatch(key):
            raise BundleError(f"{settings_path}:{line_number}: invalid setting name")
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
            value = value[1:-1]
        settings[key] = value
    return settings


def _bundle_root(settings: dict[str, str]) -> Path:
    if settings.get("MHDG_REGRESSION_SETTINGS_VERSION") != "1":
        raise BundleError("settings must define MHDG_REGRESSION_SETTINGS_VERSION=1")
    data_root = settings.get("MHDG_REGRESSION_DATA_ROOT")
    if not data_root:
        raise BundleError("settings must define MHDG_REGRESSION_DATA_ROOT")

    root = Path(data_root).expanduser()
    if not root.is_absolute():
        raise BundleError("MHDG_REGRESSION_DATA_ROOT must be an absolute path")
    try:
        root = root.resolve(strict=True)
    except FileNotFoundError as exc:
        raise BundleError(f"bundle root does not exist: {root}") from exc
    if not root.is_dir():
        raise BundleError(f"bundle root is not a directory: {root}")
    return root


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
        if _sha256(path) != artifact["sha256"]:
            raise BundleError(f"{label}.sha256 does not match the file")
        available.add(artifact_id)
        verified_bytes += actual_size

    return available, verified_bytes, warnings


def _resolve_roles(manifest: dict[str, Any]) -> dict[str, dict[str, str]]:
    artifact_ids = set(manifest["artifacts"])
    roles_by_case = {}
    for data_id, case_data in manifest["case_data"].items():
        if case_data["case_id"] != data_id:
            raise BundleError(
                f"manifest.case_data.{data_id}.case_id must equal its case_data identifier"
            )
        for role, artifact_id in case_data["roles"].items():
            if artifact_id not in artifact_ids:
                raise BundleError(
                    f"manifest.case_data.{data_id}.roles.{role} refers to "
                    f"unknown artifact {artifact_id}"
                )
        roles_by_case[data_id] = case_data["roles"]
    return roles_by_case


def _verify_case_requirements(
    cases: list[dict[str, Any]],
    roles_by_case: dict[str, dict[str, str]],
    available: set[str],
) -> list[str]:
    checked = []
    for case in cases:
        case_id = case["case_id"]
        data_id = case["external_data_id"]
        if data_id not in roles_by_case:
            continue

        required = {
            role
            for workflow in case["workflows"].values()
            for role in workflow.get("required_artifact_roles", [])
        }
        mapping = roles_by_case[data_id]
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


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as stream:
            document = json.load(stream)
    except OSError as exc:
        raise BundleError(f"cannot read {label} {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise BundleError(f"invalid JSON in {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise BundleError(f"{label} must be a JSON object")
    return document


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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
