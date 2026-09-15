"""Explicit storage cleanup, using existing records and links; no retention database."""

from pathlib import Path
import os
import shutil

from .config import machine_settings
from .documents import load_json
from .support import BundleError

MARKERS = ("refresh.json", "campaign.json", "manifest.json", "run_plan.json",
           "run_metadata.json", "suite_summary.json", "profile_summary.json", "build_metadata.json")
FINISHED = {"completed", "passed", "deferred", "failed", "ready", "launch_failed", "solver_failed",
            "solver_reported_error", "missing_hdf5_output", "output_selection_failed", "output_contract_failed"}


def _paths(value):
    if isinstance(value, dict):
        for item in value.values():
            yield from _paths(item)
    elif isinstance(value, list):
        for item in value:
            yield from _paths(item)
    elif isinstance(value, str) and value.startswith("/"):
        yield Path(value).resolve()


def _kind(directory, files):
    marker = next((name for name in MARKERS if name in files), None)
    if marker is None:
        return None
    record = load_json(directory / marker, "cleanup record")
    if marker == "manifest.json":
        kind = record.get("bundle_class", "unknown")
        return (kind, [] if kind in {"candidate", "golden"} else ["source or unknown bundle"])
    if marker == "campaign.json":
        return "unknown", ["historical campaign; manual review required"]
    if marker == "run_plan.json":
        metadata = directory / "run_metadata.json"
        record = load_json(metadata, "run metadata") if metadata.is_file() else {}
    kind = {"refresh.json": "refresh", "run_plan.json": "run", "run_metadata.json": "run",
            "suite_summary.json": "suite", "profile_summary.json": "profile",
            "build_metadata.json": "build"}[marker]
    status = record.get("status", "unrecorded")
    return kind, [] if status in FINISHED else [f"active or unfinished: {status}"]


def inventory(root, *, keep=(), settings=None):
    """Inspect one complete storage tree plus optional external retained trees.

    Absolute paths in records and filesystem links protect their targets. Bundle
    provenance is historical, not a runtime dependency of a self-contained bundle.
    Nothing traverses a directory symlink or reads an HDF5 payload.
    """
    root = root.expanduser().resolve(strict=True)
    if not root.is_dir():
        raise BundleError(f"cleanup root is not a directory: {root}")
    rows, links = [], []

    def scan(directory, owner=None, bundled=False):
        entries = sorted(directory.iterdir())
        files = {p.name for p in entries if not p.is_symlink() and p.is_file()}
        kind = _kind(directory, files) if owner is None else None
        row = None
        if kind:
            row = {"path": directory, "kind": kind[0], "bytes": 0, "reasons": kind[1]}
            rows.append(row)
            owner = directory
            bundled = "manifest.json" in files
        start = len(rows)
        size = own_size = 0
        for path in entries:
            if path.is_symlink():
                size += path.lstat().st_size
                own_size += path.lstat().st_size
                links.append((path, path.resolve(), True))
            elif path.is_dir():
                size += scan(path, owner, bundled)
            elif path.is_file():
                size += path.stat().st_size
                own_size += path.stat().st_size
                if path.suffix == ".json" and not bundled:
                    record = load_json(path, "cleanup dependency record")
                    links.extend((path, target, False) for target in _paths(record))
        if row is not None:
            row["bytes"] = size
        elif owner is None and all(item["kind"] == "unknown" for item in rows[start:]):
            del rows[start:]
            rows.append({"path": directory, "kind": "unknown", "bytes": size,
                         "reasons": ["unrecognized data; manual review required"]})
        elif owner is None and own_size:
            rows.append({"path": directory, "kind": "unknown", "bytes": own_size,
                         "reasons": ["loose files in a container; manual review required"]})
        return size

    scan(root)
    protected = [Path(path).expanduser().resolve(strict=True) for path in keep]
    for path in protected:
        if not path.is_dir():
            raise BundleError(f"--keep must name a directory: {path}")
        if not path.is_relative_to(root):
            scan(path)
    _, defaults = machine_settings(settings)
    protected.extend(defaults.get("bundles", {}).values())
    if defaults.get("build"):
        protected.append(defaults["build"])
    for row in rows:
        path = row["path"]
        if path == root or not path.is_relative_to(root):
            row["reasons"].append("scan root or external retained data")
        if any(p.is_relative_to(path) or path.is_relative_to(p) for p in protected):
            row["reasons"].append("configured default or explicitly retained")
    return root, rows, links


def cleanup(root, selected=(), *, delete=False, keep=(), settings=None):
    """Preview exact selections; remove only whole recognized, unneeded directories."""
    root, rows, links = inventory(root, keep=keep, settings=settings)
    chosen = set()
    for value in selected:
        path = Path(value).expanduser()
        path = Path(os.path.abspath(root / path))
        if path.is_symlink() or path.resolve() != path or not path.is_relative_to(root):
            raise BundleError(f"selection must be a real directory within the cleanup root: {path}")
        chosen.add(path)
    indexed = {row["path"]: row for row in rows}
    unknown = chosen - indexed.keys()
    if unknown:
        raise BundleError("select whole inventoried directories: " + ", ".join(map(str, sorted(unknown))))
    if delete and not chosen:
        raise BundleError("--delete requires explicit directory selections")
    # A retained record/link owns its dependency. Selecting both releases that
    # dependency without a graph, recursive deletion policy or cascade option.
    for origin, target, directory_link in links:
        if any(origin.is_relative_to(path) for path in chosen):
            continue
        for path in chosen:
            if target.is_relative_to(path) or (directory_link and path.is_relative_to(target)):
                indexed[path]["reasons"].append(f"needed by {origin}")
    for row in rows:
        row["selected"] = row["path"] in chosen
        row["reasons"] = sorted(set(row["reasons"]))
    blocked = [row for row in rows if row["selected"] and row["reasons"]]
    if delete:
        if blocked:
            raise BundleError("protected cleanup selection: " + "; ".join(
                f"{row['path']}: {row['reasons'][0]}" for row in blocked))
        if not shutil.rmtree.avoids_symlink_attacks:
            raise BundleError("this platform does not provide symlink-safe directory removal")
        for path in sorted(chosen):
            if path.is_symlink() or path.resolve() != path:
                raise BundleError(f"cleanup path changed: {path}")
            shutil.rmtree(path)
    return {"root": root, "rows": [row for row in rows if row["path"].is_relative_to(root)],
            "status": "removed" if delete else "preview", "blocked": bool(blocked),
            "selected_bytes": sum(indexed[path]["bytes"] for path in chosen)}
