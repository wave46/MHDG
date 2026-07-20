#!/usr/bin/env python3
"""Promote accepted regression results into a complete golden bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from check_bundle import (
    BundleError,
    ValidationSummary,
    bundle_root_from_settings,
    load_case_definition,
    read_settings,
    validate_bundle_root,
)
from compare_hdf5 import ComparisonError
from compare_run import select_candidate
from reference_matrix import (
    MatrixRun,
    collect_matrix_runs,
    install_reference_matrix,
    validate_matrix_summary,
)


ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
PROVENANCE_FILES = {
    "suite_summary.json": ("golden_reference_suite_summary", "application/json"),
    "run_plan.json": ("golden_reference_run_plan", "application/json"),
    "run_metadata.json": ("golden_reference_run_metadata", "application/json"),
    "comparison.json": ("golden_reference_comparison", "application/json"),
    "stdout.log": ("golden_reference_stdout", "text/plain"),
    "stderr.log": ("golden_reference_stderr", "text/plain"),
}


@dataclass(frozen=True)
class CanonicalRun:
    directory: Path
    plan: dict[str, Any]
    metadata: dict[str, Any]
    comparison: dict[str, Any]
    solution: Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_summary", metavar="SUITE_SUMMARY", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bundle-version", required=True, metavar="VERSION")
    parser.add_argument("--settings", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    try:
        result = promote_bundle(
            args.settings,
            args.suite_summary,
            args.output,
            args.bundle_version,
            args.cases,
        )
    except (BundleError, ComparisonError) as exc:
        print(f"golden bundle promotion failed: {exc}", file=sys.stderr)
        return 1

    print(f"golden bundle created: {args.output.expanduser().resolve()}")
    print(f"bundle: {result.bundle_id} version {result.bundle_version}")
    print(f"verified {result.artifact_count} artifacts ({result.verified_bytes} bytes)")
    return 0


def promote_bundle(
    settings_path: Path,
    summary_path: Path,
    output: Path,
    bundle_version: str,
    case_dir: Path,
) -> ValidationSummary:
    """Copy a source bundle and install accepted reference results."""
    if not bundle_version.strip():
        raise BundleError("golden bundle version must not be empty")

    source_bundle = bundle_root_from_settings(read_settings(settings_path))
    validate_bundle_root(source_bundle, case_dir)
    source_manifest = _load_json(source_bundle / "manifest.json", "bundle manifest")
    summary_path = _existing_file(summary_path, "suite summary")
    summary = _load_json(summary_path, "suite summary")
    promotion_kind = _validate_summary(summary)

    case = load_case_definition(summary["case_id"], case_dir)
    run = None
    matrix_runs: list[MatrixRun] = []
    reference_id = None
    if promotion_kind == "canonical":
        workflow = case["workflows"].get(summary["workflow_id"])
        if workflow is None:
            raise BundleError("suite summary refers to an unknown workflow")
        run = _canonical_run(summary, workflow["default_layout"])
        reference_id, source_reference = _source_reference(
            source_bundle, source_manifest, case
        )
        _validate_run_sources(
            run, source_bundle, source_manifest, source_reference
        )
    else:
        matrix_runs = collect_matrix_runs(
            summary, case, source_bundle, source_manifest
        )

    output = output.expanduser().resolve()
    if output.exists():
        raise BundleError(f"output already exists: {output}")
    if _is_within(output, source_bundle):
        raise BundleError("golden output must be outside the source bundle")

    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=f".{output.name}.", dir=output.parent
        ) as workspace:
            staging = Path(workspace) / "bundle"
            shutil.copytree(source_bundle, staging)
            manifest = _load_json(staging / "manifest.json", "bundle manifest")
            if run is not None and reference_id is not None:
                _install_reference(staging, manifest, reference_id, run.solution)
                _install_provenance(
                    staging,
                    manifest,
                    summary_path,
                    summary,
                    run,
                    case,
                    source_manifest,
                    reference_id,
                )
            else:
                install_reference_matrix(
                    staging,
                    manifest,
                    summary_path,
                    summary,
                    matrix_runs,
                    case,
                    source_manifest,
                )
            manifest["bundle_version"] = bundle_version
            manifest["bundle_class"] = "golden"
            manifest["created_utc"] = _utc_now()
            _write_json(staging / "manifest.json", manifest)
            result = validate_bundle_root(staging, case_dir)
            staging.rename(output)
    except OSError as exc:
        raise BundleError(f"cannot create golden bundle {output}: {exc}") from exc
    return result


def _canonical_run(summary: dict[str, Any], layout_id: str) -> CanonicalRun:
    matching = [
        result for result in summary["results"] if result["layout_id"] == layout_id
    ]
    if len(matching) != 1:
        raise BundleError(f"suite must contain canonical layout {layout_id} once")

    directory = _existing_directory(Path(matching[0]["run_directory"]), "run")
    plan = _load_json(directory / "run_plan.json", "run plan")
    metadata = _load_json(directory / "run_metadata.json", "run metadata")
    comparison = _load_json(directory / "comparison.json", "comparison")
    expected = {
        "case_id": summary["case_id"],
        "workflow_id": summary["workflow_id"],
        "layout_id": layout_id,
    }
    for name, value in expected.items():
        if plan.get(name) != value or comparison.get(name) != value:
            raise BundleError(f"canonical run has inconsistent {name}")
    if metadata.get("status") != "completed" or comparison.get("status") != "passed":
        raise BundleError("canonical run evidence is not fully passing")

    solution = select_candidate(directory, metadata)
    reported = _existing_file(Path(comparison.get("candidate", "")), "candidate")
    if solution != reported:
        raise BundleError("selected and compared canonical solutions differ")
    if not _same_file(_file_identity(solution), comparison.get("files", {}).get("candidate")):
        raise BundleError("canonical solution changed after comparison")
    return CanonicalRun(directory, plan, metadata, comparison, solution)


def _source_reference(
    bundle: Path,
    manifest: dict[str, Any],
    case: dict[str, Any],
) -> tuple[str, Path]:
    data_id = case["external_data_id"]
    try:
        reference_id = manifest["case_data"][data_id]["roles"]["warm_reference"]
        relative_path = manifest["artifacts"][reference_id]["path"]
    except KeyError as exc:
        raise BundleError("source bundle has no warm reference artifact") from exc
    return reference_id, _existing_file(bundle / relative_path, "warm reference")


def _validate_run_sources(
    run: CanonicalRun,
    source_bundle: Path,
    manifest: dict[str, Any],
    source_reference: Path,
) -> None:
    _validate_bundle_identity(run.plan, source_bundle, manifest)

    comparison = run.comparison
    reported = _existing_file(Path(comparison.get("reference", "")), "reference")
    if reported != source_reference:
        raise BundleError("canonical run compared against another reference")
    if not _same_file(
        _file_identity(source_reference),
        comparison.get("files", {}).get("reference"),
    ):
        raise BundleError("source reference changed after comparison")


def _validate_bundle_identity(
    plan: dict[str, Any], source_bundle: Path, manifest: dict[str, Any]
) -> None:
    bundle = plan.get("bundle")
    if not isinstance(bundle, dict) or not isinstance(bundle.get("root"), str):
        raise BundleError("run plan has no bundle identity")
    if _existing_directory(Path(bundle["root"]), "run bundle") != source_bundle:
        raise BundleError("run did not use the configured source bundle")
    if (
        bundle.get("bundle_id") != manifest.get("bundle_id")
        or bundle.get("bundle_version") != manifest.get("bundle_version")
    ):
        raise BundleError("run and source bundle identities differ")


def _install_reference(
    staging: Path,
    manifest: dict[str, Any],
    artifact_id: str,
    solution: Path,
) -> None:
    artifact = manifest["artifacts"][artifact_id]
    target = staging / artifact["path"]
    shutil.copy2(solution, target)
    artifact.update(_file_identity(target))
    artifact["path"] = target.relative_to(staging).as_posix()


def _install_provenance(
    staging: Path,
    manifest: dict[str, Any],
    summary_path: Path,
    summary: dict[str, Any],
    run: CanonicalRun,
    case: dict[str, Any],
    source_manifest: dict[str, Any],
    reference_id: str,
) -> None:
    directory = staging / "provenance/golden_reference"
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    for artifact_id in list(manifest["artifacts"]):
        if artifact_id.startswith("golden_reference_"):
            del manifest["artifacts"][artifact_id]

    record = {
        "schema_version": 1,
        "status": "golden",
        "created_utc": _utc_now(),
        "case_id": summary["case_id"],
        "workflow_id": summary["workflow_id"],
        "canonical_layout": run.plan["layout_id"],
        "suite_id": summary["suite_id"],
        "suite_run_id": summary["run_id"],
        "reference_artifact": reference_id,
        "source_bundle": {
            "bundle_id": source_manifest["bundle_id"],
            "bundle_version": source_manifest["bundle_version"],
        },
        "tracked_reference": {
            "branch": case["reference_branch"],
            "revision": case["reference_revision"],
        },
        "solver": run.metadata.get("solver", {}),
        "executable": run.metadata.get("executable", {}),
    }
    record_path = directory / "golden_reference.json"
    _write_json(record_path, record)
    _register_artifact(
        staging,
        manifest,
        "golden_reference_record",
        record_path,
        "application/json",
    )

    sources = {
        "suite_summary.json": summary_path,
        **{
            filename: run.directory / filename
            for filename in PROVENANCE_FILES
            if filename != "suite_summary.json"
        },
    }
    for filename, source in sources.items():
        artifact_id, media_type = PROVENANCE_FILES[filename]
        target = directory / filename
        shutil.copy2(_existing_file(source, filename), target)
        _register_artifact(staging, manifest, artifact_id, target, media_type)


def _register_artifact(
    staging: Path,
    manifest: dict[str, Any],
    artifact_id: str,
    path: Path,
    media_type: str,
) -> None:
    manifest["artifacts"][artifact_id] = {
        "path": path.relative_to(staging).as_posix(),
        **_file_identity(path),
        "media_type": media_type,
    }


def _validate_summary(summary: dict[str, Any]) -> str:
    required = {
        "schema_version",
        "status",
        "suite_id",
        "run_id",
        "case_id",
        "results",
    }
    missing = sorted(required - summary.keys())
    if missing:
        raise BundleError(f"suite summary is missing: {', '.join(missing)}")
    if summary["schema_version"] != 1 or summary["status"] != "passed":
        raise BundleError("only a passing version-1 suite summary can be promoted")
    for name in ("suite_id", "run_id", "case_id"):
        if not isinstance(summary[name], str) or not ID_RE.fullmatch(summary[name]):
            raise BundleError(f"suite summary has invalid {name}")

    if "workflow_id" in summary:
        _validate_canonical_summary(summary)
        return "canonical"
    validate_matrix_summary(summary)
    return "matrix"


def _validate_canonical_summary(summary: dict[str, Any]) -> None:
    workflow_id = summary["workflow_id"]
    if not isinstance(workflow_id, str) or not ID_RE.fullmatch(workflow_id):
        raise BundleError("suite summary has invalid workflow_id")
    results = summary["results"]
    if not isinstance(results, list) or not results:
        raise BundleError("suite summary contains no layout results")
    layouts = []
    for result in results:
        if not isinstance(result, dict) or any(
            result.get(name) != expected
            for name, expected in (
                ("status", "passed"),
                ("run_status", "completed"),
                ("comparison_status", "passed"),
            )
        ):
            raise BundleError("suite summary contains a non-passing layout")
        layout = result.get("layout_id")
        if not isinstance(layout, str) or not ID_RE.fullmatch(layout):
            raise BundleError("suite summary has an invalid layout identifier")
        layouts.append(layout)
    if len(set(layouts)) != len(layouts):
        raise BundleError("suite summary contains duplicate layouts")


def _same_file(first: dict[str, Any], second: Any) -> bool:
    return isinstance(second, dict) and (
        first["size_bytes"] == second.get("size_bytes")
        and first["sha256"] == second.get("sha256")
    )


def _file_identity(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        return {"size_bytes": path.stat().st_size, "sha256": digest.hexdigest()}
    except OSError as exc:
        raise BundleError(f"cannot checksum {path}: {exc}") from exc


def _existing_file(path: Path, label: str) -> Path:
    try:
        path = path.expanduser().resolve(strict=True)
    except (FileNotFoundError, OSError) as exc:
        raise BundleError(f"{label} does not exist: {path}") from exc
    if not path.is_file():
        raise BundleError(f"{label} is not a file: {path}")
    return path


def _existing_directory(path: Path, label: str) -> Path:
    try:
        path = path.expanduser().resolve(strict=True)
    except (FileNotFoundError, OSError) as exc:
        raise BundleError(f"{label} does not exist: {path}") from exc
    if not path.is_dir():
        raise BundleError(f"{label} is not a directory: {path}")
    return path


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise BundleError(f"cannot read {label} {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise BundleError(f"invalid JSON in {label} {path}: {exc.msg}") from exc
    if not isinstance(document, dict):
        raise BundleError(f"{label} must contain a JSON object")
    return document


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace(
        "+00:00", "Z"
    )


if __name__ == "__main__":
    raise SystemExit(main())
