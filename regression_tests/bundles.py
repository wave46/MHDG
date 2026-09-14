"""Bundle readiness, standalone creation and validation of recorded artifacts."""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from .catalog import load_case_definition, required_case_roles
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


def bundle_readiness(case_id, source, case_dir=CASES, *, workflows=()):
    """Report presence and declared producers; do not run producers or hash files."""
    source = require_directory(source, "source")
    case = load_case_definition(case_id, case_dir)
    required = _required_roles(case, workflows)
    rows = []
    for role, spec in case["bundle_files"].items():
        filename = spec["filename"]
        path = source / _relative_path(filename, role)
        present = path.is_file()
        if path.exists() and not present:
            raise BundleError(f"bundle source is not a file: {path}")
        producers = {
            name: sorted(required_case_roles(case, name))
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


def create_bundle(case_id, source, output, case_dir=CASES, bundle_version="1.0.0", *, workflows=()):
    """Copy available declared files, then validate before publishing a candidate."""
    source = require_directory(source, "source")
    output = output.expanduser().absolute()
    if output.exists() or output.is_symlink():
        raise BundleError(f"output already exists: {output}")
    output = output.resolve()
    if output.is_relative_to(source):
        raise BundleError("output must be outside the prepared source directory")
    if not bundle_version:
        raise BundleError("bundle version must not be empty")
    report = bundle_readiness(case_id, source, case_dir, workflows=workflows)
    if report["missing"]:
        raise BundleError("required files missing for roles: " + ", ".join(report["missing"]))
    case = load_case_definition(case_id, case_dir)
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
                artifact_id = spec["artifact_id"]
                if target.exists() or artifact_id in manifest["artifacts"]:
                    raise BundleError(f"duplicate bundle file or artifact: {role}")
                shutil.copy2(source_path, target)
                manifest["roles"][role] = artifact_id
                manifest["artifacts"][artifact_id] = {
                    "path": target.relative_to(staging).as_posix(), **file_identity(target),
                    "media_type": spec["media_type"], **({"optional": True} if spec["optional"] else {}),
                }
            write_json_direct(staging / "manifest.json", manifest)
            summary = validate_bundle_root(staging, case_dir, workflows=workflows)
            staging.rename(output)
    except OSError as exc:
        raise BundleError(f"cannot create bundle {output}: {exc}") from exc
    return summary


def validate_bundle_root(bundle_root, case_dir=CASES, *, workflows=(), required_class=None):
    """Verify all recorded files and base plus selected workflow requirements."""
    bundle_root = require_directory(bundle_root, "bundle root")
    manifest = load_manifest(bundle_root, case_dir)
    if required_class and manifest.get("bundle_class") != required_class:
        actual = manifest.get("bundle_class", "unspecified")
        raise BundleError(f"golden-check requires bundle_class={required_class}; found {actual}")
    case = load_case_definition(manifest["case_id"], case_dir)
    available, verified_bytes, warnings = _verify_artifacts(bundle_root, manifest["artifacts"])
    roles = _verify_roles(manifest)
    load_reference_matrix(bundle_root, case["case_id"], case_dir, manifest=manifest)
    required = _required_roles(case, workflows)
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


@dataclass(frozen=True)
class ReferenceMatrix:
    bundle_id: str
    bundle_version: str
    references: dict[tuple[str, str, str], Path]

    def reference_for(self, workflow_id, layout_id, stage_id):
        key = (workflow_id, layout_id, stage_id)
        try:
            return self.references[key]
        except KeyError as exc:
            raise BundleError(f"golden matrix has no reference for {'/'.join(key)}") from exc


def load_reference_matrix(bundle_root, case_id, case_dir=CASES, *, manifest=None):
    """Read and validate the optional stage index; artifact hashing belongs to bundle validation."""
    bundle_root = require_directory(bundle_root, "golden bundle")
    if manifest is None:
        manifest = load_manifest(bundle_root, case_dir, case_id=case_id)
    index_id = manifest["roles"].get("reference_matrix")
    if index_id is None:
        return None

    def artifact(artifact_id, media_type):
        item = manifest["artifacts"].get(artifact_id)
        if item is None:
            raise BundleError(f"reference matrix refers to unknown artifact {artifact_id}")
        if item["media_type"] != media_type:
            raise BundleError(f"reference matrix artifact {artifact_id} is not {media_type}")
        return artifact_path(bundle_root, item["path"], f"manifest.artifacts.{artifact_id}")

    matrix = load_validated_json(artifact(index_id, "application/json"),
                                 case_dir.parent / "schemas/reference-matrix.schema.json",
                                 "reference matrix")
    if matrix["case_id"] != case_id:
        raise BundleError("reference matrix has the wrong case_id")
    references = {}
    for entry in matrix["references"]:
        key = (entry["workflow_id"], entry["layout_id"], entry["stage_id"])
        if key in references:
            raise BundleError(f"reference matrix contains duplicate cell {'/'.join(key)}")
        references[key] = artifact(entry["artifact_id"], "application/x-hdf5")
    return ReferenceMatrix(manifest["bundle_id"], manifest["bundle_version"], references)


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
