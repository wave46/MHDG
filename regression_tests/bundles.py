"""Bundle readiness, standalone creation and validation of recorded artifacts."""

from __future__ import annotations

import shutil
import tempfile
import math
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from .catalog import load_case_definition, required_case_roles
from .config import bundle_root_from_settings
from .documents import load_validated_json, write_json_direct
from .support import BundleError, MissingArtifactError, utc_now
from .files import file_identity, sha256_digest, require_directory

CASES = Path(__file__).resolve().parent / "cases"


@dataclass(frozen=True)
class ValidationSummary:
    bundle_id: str
    bundle_version: str
    artifact_count: int
    verified_artifact_count: int
    verified_bytes: int
    case_id: str
    warnings: list[str]


def _required_roles(case, workflows):
    required = required_case_roles(case)
    for workflow in workflows:
        required.update(required_case_roles(case, workflow))
    return required


def load_manifest(root, case_dir=CASES, *, case_id=None):
    manifest = load_validated_json(root / "manifest.json",
                               case_dir.parent / "schemas/bundle-manifest.schema.json",
                               "bundle manifest")
    if case_id is not None and manifest["case_id"] != case_id:
        raise BundleError(f"bundle does not contain case data for {case_id}")
    return manifest


def bundle_readiness(case_id, source, case_dir=CASES, *, workflows=(), catalog=None):
    """Report presence and declared producers; do not run producers or hash files."""
    catalog = {} if catalog is None else catalog
    source = require_directory(source, "source")
    case = load_case_definition(case_id, case_dir, catalog=catalog)
    required = _required_roles(case, workflows)
    indexes = _series_roles(case)
    rows = []
    for role, spec in case["bundle_files"].items():
        filename = spec["filename"]
        path = source / _relative_path(filename, role)
        present = path.is_file()
        if path.exists() and not present:
            raise BundleError(f"bundle source is not a file: {path}")
        if present and role in indexes:
            load_field_series(path)
        producers = {
            name: sorted(required_case_roles(case, name, references=False))
            for name, workflow in case["workflows"].items()
            if role in workflow.get("outputs", [])
        }
        rows.append({"role": role, "path": filename,
                     "origin": "workflow-producible" if producers else "user-supplied",
                     "requirement": "required" if role in required else "optional",
                     "presence": "present" if present else "missing", "producers": producers})
    missing = [row["role"] for row in rows
               if row["requirement"] == "required" and row["presence"] == "missing"]
    return {"case_id": case_id, "source": str(source), "workflows": list(workflows),
            "status": "missing" if missing else "ready", "missing": missing, "artifacts": rows}


def create_bundle(case_id, source, output, case_dir=CASES, bundle_version="1.0.0", *, workflows=(), catalog=None):
    """Copy available declared files, then validate before publishing a candidate."""
    catalog = {} if catalog is None else catalog
    source = require_directory(source, "source")
    output = output.expanduser().absolute()
    if output.exists() or output.is_symlink():
        raise BundleError(f"output already exists: {output}")
    output = output.resolve()
    if output.is_relative_to(source):
        raise BundleError("output must be outside the prepared source directory")
    if not bundle_version:
        raise BundleError("bundle version must not be empty")
    report = bundle_readiness(case_id, source, case_dir, workflows=workflows, catalog=catalog)
    if report["missing"]:
        raise BundleError("required files missing for roles: " + ", ".join(report["missing"]))
    case = load_case_definition(case_id, case_dir, catalog=catalog)
    indexes = _series_roles(case)
    if not case["bundle_files"]:
        raise BundleError(f"case {case_id} does not define bundle_files")
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=f".{output.name}.", dir=output.parent) as workspace:
            staging = Path(workspace) / "bundle"
            (staging / "inputs").mkdir(parents=True)
            manifest = {"schema_version": 2, "bundle_id": f"{case_id}_bundle",
                        "bundle_version": bundle_version, "bundle_class": "candidate",
                        "created_utc": utc_now(), "case_id": case_id, "roles": {}, "artifacts": {}}
            for row in report["artifacts"]:
                if row["presence"] == "missing":
                    continue
                role = row["role"]
                spec = case["bundle_files"][role]
                source_path = source / spec["filename"]
                target = staging / "inputs" / source_path.name
                if role in indexes:
                    target = staging / "inputs" / role / source_path.name
                artifact_id = spec["artifact_id"]
                _copy_artifact(staging, manifest, artifact_id, source_path, target,
                               spec["media_type"], spec["optional"])
                manifest["roles"][role] = artifact_id
                if role in indexes:
                    index = load_field_series(source_path)
                    for number, frame in enumerate(index["frames"], start=1):
                        for kind, relative in frame.items():
                            _copy_artifact(staging, manifest, f"{artifact_id}_{kind}_{number:04d}",
                                           source_path.parent / relative, target.parent / relative,
                                           "application/x-hdf5")
            write_json_direct(staging / "manifest.json", manifest)
            summary = validate_bundle_root(staging, case_dir, workflows=workflows, catalog=catalog)
            staging.rename(output)
    except OSError as exc:
        raise BundleError(f"cannot create bundle {output}: {exc}") from exc
    return summary


def validate_bundle_root(bundle_root, case_dir=CASES, *, workflows=(), required_class=None,
                         required_roles=(), manifest=None, catalog=None):
    """Verify all recorded files and base plus selected workflow requirements."""
    bundle_root = require_directory(bundle_root, "bundle root")
    manifest = manifest if manifest is not None else load_manifest(bundle_root, case_dir)
    if required_class and manifest.get("bundle_class") != required_class:
        actual = manifest.get("bundle_class", "unspecified")
        raise BundleError(f"golden-check requires bundle_class={required_class}; found {actual}")
    case = load_case_definition(manifest["case_id"], case_dir, catalog=catalog)
    available, verified_bytes, warnings = _verify_artifacts(bundle_root, manifest["artifacts"])
    roles = _verify_roles(manifest)
    series_artifact_ids(bundle_root, manifest, case)
    required = _required_roles(case, workflows) | set(required_roles)
    missing = sorted(role for role in required if roles.get(role) not in available)
    if missing:
        raise BundleError(f"case {case['case_id']} has missing required artifact roles: " + ", ".join(missing))
    return ValidationSummary(manifest["bundle_id"], manifest["bundle_version"],
                             len(manifest["artifacts"]), len(available), verified_bytes,
                             case["case_id"], warnings)


def _verify_artifacts(
    bundle_root: Path, artifacts: dict
) -> tuple[set[str], int, list[str]]:
    available: set[str] = set()
    verified_bytes = 0
    warnings: list[str] = []

    for artifact_id, artifact in artifacts.items():
        label = f"manifest.artifacts.{artifact_id}"
        try:
            path = artifact_path(bundle_root, artifact["path"], label)
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


def _verify_roles(manifest: dict) -> dict[str, str]:
    artifact_ids = set(manifest["artifacts"])
    roles = manifest["roles"]
    for role, artifact_id in roles.items():
        if artifact_id not in artifact_ids:
            raise BundleError(
                f"manifest.roles.{role} refers to unknown artifact {artifact_id}"
            )
    return roles


def _copy_artifact(root, manifest, artifact_id, source, target, media_type, optional=False):
    if target.exists() or artifact_id in manifest["artifacts"]:
        raise BundleError(f"duplicate bundle file or artifact: {artifact_id}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    manifest["artifacts"][artifact_id] = {
        "path": target.relative_to(root).as_posix(), **file_identity(target),
        "media_type": media_type, **({"optional": True} if optional else {}),
    }


def _series_roles(case):
    return {workflow["moving_equilibrium"]["series"]
            for workflow in case["workflows"].values() if workflow.get("moving_equilibrium")}


def load_field_series(path):
    """Read ordered path pairs; the bundle manifest owns their file identities."""
    index = load_validated_json(path, CASES.parent / "schemas/field-series.schema.json", "field series")
    if not math.isfinite(index["dt_s"]):
        raise BundleError("field series dt_s must be finite")
    seen = set()
    for frame in index["frames"]:
        for relative in frame.values():
            member = artifact_path(path.parent.resolve(), relative, "field series member")
            if member in seen:
                raise BundleError("field series member paths must be unique")
            seen.add(member)
    return index


def series_artifact_ids(root, manifest, case):
    """Resolve declared series members for validation and compact publication."""
    roles = _series_roles(case)
    if not roles:
        return set()
    by_path = {}
    for artifact_id, spec in manifest["artifacts"].items():
        path = (root / spec["path"]).resolve()
        by_path.setdefault(path, []).append(artifact_id)
    members = set()
    for role in roles:
        artifact_id = manifest["roles"].get(role)
        if artifact_id is None:
            continue
        artifact = manifest["artifacts"][artifact_id]
        if artifact.get("optional") and not (root / artifact["path"]).exists():
            continue
        path = artifact_path(root, artifact["path"], "field series index")
        index = load_field_series(path)
        for frame in index["frames"]:
            for relative in frame.values():
                matches = by_path.get((path.parent / relative).resolve(), [])
                if len(matches) != 1:
                    raise BundleError(f"field series member needs one manifest artifact: {relative}")
                members.add(matches[0])
    return members


def _relative_path(relative_path, label):
    path = PurePosixPath(relative_path)
    if "\\" in relative_path or path.is_absolute() or ".." in path.parts or path == PurePosixPath("."):
        raise BundleError(f"{label}.path must stay relative to the bundle root using POSIX separators")
    return path


def artifact_path(root, relative_path, label):
    candidate = root / _relative_path(relative_path, label)
    resolved = candidate.resolve()
    if not resolved.is_relative_to(root):
        raise BundleError(f"{label}.path resolves outside the bundle root")
    if not candidate.exists():
        raise MissingArtifactError(f"{label}.path does not exist: {relative_path}")
    if not resolved.is_file():
        raise BundleError(f"{label}.path is not a regular file: {relative_path}")
    return resolved


def preflight_bundles(checks, settings_by_case, case_directory, *, required_bundle_class=None, catalog=None):
    """Verify each source bundle once against the union of selected input/reference roles."""
    catalog = {} if catalog is None else catalog
    manifests = {}
    for case_id, values in settings_by_case.items():
        case = load_case_definition(case_id, case_directory, catalog=catalog)
        roles = set().union(*(required_case_roles(case, workflow, references=check["reference_comparisons"])
                              for check in checks if check["case_id"] == case_id
                              for workflow in check["workflow_ids"]))
        bundle = bundle_root_from_settings(values)
        manifest = load_manifest(bundle, case_directory, case_id=case_id)
        validate_bundle_root(bundle, case_directory, required_roles=roles, required_class=required_bundle_class,
                             manifest=manifest, catalog=catalog)
        manifests[case_id] = manifest
    return manifests
