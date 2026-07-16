#!/usr/bin/env python3
"""Create a validated MHDG regression bundle from prepared case files."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from check_bundle import (
    BundleError,
    ValidationSummary,
    load_case_definition,
    required_case_roles,
    validate_bundle_root,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        required=True,
        dest="case_id",
        metavar="CASE",
        help="tracked case identifier",
    )
    parser.add_argument(
        "--source", required=True, type=Path, help="prepared case directory"
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="new bundle path"
    )
    parser.add_argument(
        "--bundle-version",
        default="1.0.0",
        metavar="VERSION",
        help="bundle version label (default: %(default)s)",
    )
    parser.add_argument(
        "--cases", required=True, type=Path, help=argparse.SUPPRESS
    )
    args = parser.parse_args(argv)

    try:
        summary = create_bundle(
            args.case_id,
            args.source,
            args.output,
            args.cases,
            args.bundle_version,
        )
    except BundleError as exc:
        print(f"bundle creation failed: {exc}", file=sys.stderr)
        return 1

    print(f"bundle created: {args.output}")
    print(
        f"recorded {summary.artifact_count} artifacts "
        f"({summary.verified_bytes} bytes)"
    )
    return 0


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

    source = _source_directory(source)
    output = output.expanduser().resolve()
    if output.exists():
        raise BundleError(f"output already exists: {output}")
    if _is_within(output, source):
        raise BundleError("output must be outside the prepared source directory")

    required_roles = required_case_roles(case)
    undefined_roles = sorted(required_roles - file_contract.keys())
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
            _write_manifest(staging, manifest)
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
            "sha256": _sha256(target_path),
            "size_bytes": target_path.stat().st_size,
            "media_type": file_spec["media_type"],
        }
        if "description" in file_spec:
            artifacts[artifact_id]["description"] = file_spec["description"]
        roles[role] = artifact_id

    created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return {
        "schema_version": 1,
        "bundle_id": f"{data_id}_bundle",
        "bundle_version": bundle_version,
        "created_utc": created.replace("+00:00", "Z"),
        "artifacts": artifacts,
        "case_data": {
            data_id: {
                "case_id": case["case_id"],
                "description": case["description"],
                "roles": roles,
            }
        },
    }


def _write_manifest(staging: Path, manifest: dict[str, Any]) -> None:
    (staging / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def _source_directory(source: Path) -> Path:
    try:
        source = source.expanduser().resolve(strict=True)
    except FileNotFoundError as exc:
        raise BundleError(f"source directory does not exist: {source}") from exc
    if not source.is_dir():
        raise BundleError(f"source is not a directory: {source}")
    return source


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
