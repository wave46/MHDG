"""Build, run ordered producers, collect once, review, then explicitly publish."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import shutil
import tempfile

from .catalog import load_case_definition, workflow_required_roles, load_suite_definition, layout_pairs, load_layouts, required_builds
from .documents import load_json, write_json_atomic
from .support import BundleError, HarnessError, utc_now, utc_run_id
from .files import file_identity
from . import config
from .build import build_solver
from .bundles import artifact_path, validate_bundle_root
from .compare import compare_completed_run, producer_converged, validate_completed_run
from .compare_common import select_candidate
from .execute import execute_prepared, reusable_outputs
from .prepare import prepare_run
from .suites import compare_layout_pairs, run_suite

ROOT = Path(__file__).resolve().parent


def _recipe(case_id, root):
    case = load_case_definition(case_id, root / "cases")
    document = load_json(root / "golden.json", "golden producers")
    declaration = document.get("cases", {}).get(case_id)
    if not declaration:
        raise BundleError(f"no golden producers declared for {case_id}")
    layouts = load_layouts(root / "layouts.json")
    producers = []
    for entry in declaration["producers"]:
        entry = {"workflow": entry} if isinstance(entry, str) else entry
        name = entry["workflow"]
        if name not in case["workflows"]:
            raise BundleError(f"unknown producer workflow: {name}")
        workflow = case["workflows"][name]
        canonical = workflow["layout"]
        selected = list(layouts) if entry.get("matrix") else [canonical]
        if canonical not in layouts:
            raise BundleError(f"unknown producer layout: {canonical}")
        producers.append({"workflow": name, "layouts": selected,
                          "roles": workflow.get("outputs", []), "canonical": canonical,
                          "matrix": document["matrix"] if entry.get("matrix") else None})
    if not producers or len({item["workflow"] for item in producers}) != len(producers):
        raise BundleError("golden producers must be nonempty and unique")
    checks = [load_suite_definition(name, root / "suites.json", root / "layouts.json",
                                          root / "cases", case_id=case_id)
              for name in declaration["checks"]]
    if not checks or len(set(declaration["checks"])) != len(checks):
        raise BundleError("golden checks must be nonempty and unique")
    requirements = set().union(*(required_builds(case, [item["workflow"]], [layouts[name] for name in item["layouts"]])
                                 for item in producers))
    for check in checks:
        requirements.update(required_builds(case, check["workflow_ids"], [layouts[name] for name in check["layouts"]]))
    return case, producers, declaration["checks"], requirements


def _record(path, base=None):
    return {"path": str(path.relative_to(base) if base else path), **file_identity(path)}


def _verify_record(base, record):
    path = artifact_path(base, record["path"], "refresh evidence")
    if file_identity(path) != {key: record[key] for key in ("sha256", "size_bytes")}:
        raise BundleError(f"refresh evidence changed: {record['path']}")
    return path


def refresh(case_id, settings, workspace, jobs=None, *, catalog_root=None):
    """Generate a reviewable candidate; never publish or resume a campaign."""
    root = catalog_root or ROOT
    case, producers, checks, requirements = _recipe(case_id, root)
    source = Path(settings["MHDG_REGRESSION_DATA_ROOT"]).resolve()
    validate_bundle_root(source, root / "cases")
    manifest = load_json(source / "manifest.json", "source bundle")
    if manifest["case_id"] != case_id:
        raise BundleError(f"source bundle does not contain {case_id}")
    available = {role for role, artifact in manifest["roles"].items()
                 if (source / manifest["artifacts"][artifact]["path"]).is_file()}
    for producer in producers:
        workflow = case["workflows"][producer["workflow"]]
        needed = workflow_required_roles({key: value for key, value in workflow.items() if key != "reference"})
        missing = sorted(needed - available)
        if missing:
            raise BundleError(f"producer {producer['workflow']} needs: {', '.join(missing)}")
        available.update(producer["roles"])
    workspace = workspace.expanduser().absolute()
    if workspace.exists() or workspace.is_symlink():
        raise BundleError(f"refresh workspace already exists: {workspace}")
    workspace.mkdir(parents=True)
    workspace = workspace.resolve()
    run_id = utc_run_id()
    report = {"schema_version": 1, "case_id": case_id, "run_id": run_id,
              "created_utc": utc_now(), "status": "running", "producers": [], "checks": [],
              "source": {"bundle_id": manifest["bundle_id"], "bundle_version": manifest["bundle_version"],
                         "manifest": _record(source / "manifest.json")},
              "catalogs": [_record(root / name, root) for name in
                           (f"cases/{case_id}.json", "workflows.json", "layouts.json", "tolerances.json",
                            "suites.json", "golden.json")]}
    path = workspace / "refresh.json"
    write_json_atomic(path, report, "refresh report")
    try:
        build = build_solver(settings, ROOT.parent, jobs, requirements=requirements)
        values = {**settings, **config.build_settings(build.metadata_path),
                  "MHDG_REGRESSION_RUN_ROOT": str(workspace / "runs")}
        config.runtime_settings(values, requirements)
        shutil.copy2(build.metadata_path, workspace / "build.json")
        generated = {}
        for producer in producers:
            results = []
            for layout in producer["layouts"]:
                name = producer["workflow"]
                print(f"producer: {name} / {layout}", flush=True)
                prepared = prepare_run(values, case_id, name, layout, root / "cases", root / "layouts.json",
                                       run_id, validate_bundle=False,
                                       requested_overrides={"balance_diagnostics_mode": "off"},
                                       artifact_overrides=generated, require_reference=False)
                result = {"workflow_id": name, "layout_id": layout, "run_directory": str(prepared.path),
                          "status": "running", "run_status": "running"}
                report["producers"].append(result)
                write_json_atomic(path, report, "refresh report")
                execution = execute_prepared(prepared, values)
                result.update(run_status=execution.status, status=execution.status)
                _validate_producer(result, root)
                result["status"] = "passed"
                result["old_reference"] = _old_reference(result, source, manifest, case, root)
                results.append(result)
                write_json_atomic(path, report, "refresh report")
            if producer["matrix"]:
                matrix = producer["matrix"]
                pairs = layout_pairs(load_layouts(root / "layouts.json"), matrix["relations"])
                comparisons = compare_layout_pairs({**matrix, "workflow_ids": [name], "results": results,
                                                    "layout_comparisons": pairs}, root / "cases", root / "tolerances.json")
                report.setdefault("parallel_checks", []).extend(comparisons)
                if any(item["status"] != "passed" for item in comparisons):
                    raise BundleError(f"producer parallel comparison failed: {name}")
            canonical = next(item for item in results if item["layout_id"] == producer["canonical"])
            for role in producer["roles"]:
                generated[role] = _solution(Path(canonical["run_directory"]))
        candidate = workspace / "candidate"
        _collect(source, candidate, manifest, generated, report, case, root)
        values["MHDG_REGRESSION_DATA_ROOT"] = str(candidate)
        for name in checks:
            summary_path, summary = run_suite(values, name, root / "cases", root / "layouts.json",
                                              root / "suites.json", root / "tolerances.json", f"{run_id}-verify",
                                              case_id=case_id, required_bundle_class="candidate")
            report["checks"].append({"suite": name, "status": summary["status"],
                                     "summary": str(summary_path.relative_to(workspace))})
            if summary["status"] != "passed":
                raise BundleError(f"golden validation failed: {name}")
        validate_bundle_root(candidate, root / "cases")
        report["evidence"] = [_record(item, workspace) for item in sorted(workspace.rglob("*"))
                              if item.is_file() and not item.is_symlink() and item != path
                              and (item.suffix in {".json", ".log"} or item.name == "param.txt")]
        report["status"] = "ready"
    except (HarnessError, OSError, KeyboardInterrupt) as exc:
        if report["producers"] and report["producers"][-1]["status"] != "passed":
            report["producers"][-1]["status"] = "failed"
        report.update(status="failed", failure=str(exc) or "interrupted")
        write_json_atomic(path, report, "refresh report")
        raise
    write_json_atomic(path, report, "refresh report")
    return path, report


def _solution(directory):
    return select_candidate(directory, load_json(directory / "run_metadata.json", "run metadata"))


def _validate_producer(result, root):
    directory = Path(result["run_directory"])
    if result["run_status"] != "completed" or not reusable_outputs(directory):
        raise BundleError(f"producer failed or its outputs changed: {result['workflow_id']}")
    validation = validate_completed_run(directory, root / "cases", root / "tolerances.json")
    result["validation"] = validation
    if validation["status"] != "passed":
        raise BundleError(f"invalid producer output: {result['workflow_id']}: {validation['failures']}")


def _old_reference(result, source, manifest, case, root):
    directory = Path(result["run_directory"])
    workflow = case["workflows"][result["workflow_id"]]
    artifact = manifest["roles"].get(workflow.get("reference"))
    reference = source / manifest["artifacts"][artifact]["path"] if artifact else None
    if reference is None or not reference.is_file():
        return {"status": "unavailable"}
    try:
        use_matrix = bool(workflow.get("stages") and "reference_matrix" in manifest["roles"])
        _, path, report = compare_completed_run(directory, root / "cases", root / "tolerances.json",
                                                reference_override=None if use_matrix else reference,
                                                report_override=directory / "old_reference.json")
        return {"status": report["status"], "report": str(path)}
    except HarnessError as exc:
        return {"status": "unavailable", "reason": str(exc)}


def _install(staging, manifest, artifact_id, source, relative, media_type):
    target = staging / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    manifest["artifacts"][artifact_id] = {"path": relative, **file_identity(target), "media_type": media_type}


def _collect(source, candidate, source_manifest, generated, report, case, root):
    """Collect one candidate; later publication copies this validated bundle."""
    candidate.mkdir()
    manifest = deepcopy(source_manifest)
    matrix_id = manifest["roles"].get("reference_matrix", "golden_matrix_index")
    old_index = manifest["artifacts"].get(matrix_id)
    entries = load_json(source / old_index["path"], "reference matrix")["references"] if old_index else []
    keep = set(manifest["roles"].values()) | {item["artifact_id"] for item in entries}
    manifest["artifacts"] = {key: value for key, value in manifest["artifacts"].items() if key in keep}
    for role, solution in generated.items():
        spec = case["bundle_files"][role]
        artifact_id = manifest["roles"].get(role, spec["artifact_id"])
        _install(candidate, manifest, artifact_id, solution, f"inputs/{spec['filename']}", spec["media_type"])
        manifest["roles"][role] = artifact_id
    # Keep any untouched stage references from the source, replacing produced cells.
    cells = {(item["workflow_id"], item["layout_id"], item["stage_id"]): item for item in entries}
    for result in report["producers"]:
        directory = Path(result["run_directory"])
        metadata = load_json(directory / "run_metadata.json", "producer metadata")
        for stage in metadata.get("stages", []):
            key = (result["workflow_id"], result["layout_id"], stage["stage_id"])
            artifact_id = "golden_matrix_" + "_".join(key)
            previous = cells.get(key, {}).get("artifact_id")
            if previous and previous != artifact_id and previous not in manifest["roles"].values():
                manifest["artifacts"].pop(previous, None)
            relative = "references/matrix/" + "/".join(key) + ".h5"
            _install(candidate, manifest, artifact_id, Path(stage["selected_hdf5"]), relative, "application/x-hdf5")
            cells[key] = dict(zip(("workflow_id", "layout_id", "stage_id"), key), artifact_id=artifact_id)
    if cells:
        index = {"schema_version": 2, "created_utc": utc_now(), "case_id": case["case_id"],
                 "suite_id": "golden_refresh", "suite_run_id": report["run_id"],
                 "source_bundle": {key: source_manifest[key] for key in ("bundle_id", "bundle_version")},
                 "tracked_reference": {"branch": case["reference_branch"], "revision": case["reference_revision"]},
                 "references": list(cells.values())}
        path = candidate / "references/matrix/index.json"
        write_json_atomic(path, index, "reference matrix")
        manifest["artifacts"][matrix_id] = {"path": path.relative_to(candidate).as_posix(),
                                            **file_identity(path), "media_type": "application/json"}
        manifest["roles"]["reference_matrix"] = matrix_id
    for artifact_id, artifact in list(manifest["artifacts"].items()):
        target = candidate / artifact["path"]
        if target.exists():
            continue
        original = source / artifact["path"]
        if not original.is_file() and artifact.get("optional"):
            del manifest["artifacts"][artifact_id]
            manifest["roles"] = {role: value for role, value in manifest["roles"].items() if value != artifact_id}
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, target)
    manifest.update(bundle_class="candidate", bundle_version=f"refresh-{report['run_id']}", created_utc=utc_now())
    write_json_atomic(candidate / "manifest.json", manifest, "candidate manifest")
    validate_bundle_root(candidate, root / "cases")


def publish(workspace, output, bundle_version, reason, provenance, *, catalog_root=None):
    """Revalidate recorded evidence and publish a self-contained golden once."""
    catalog_root = catalog_root or ROOT
    if not all(isinstance(value, str) and value.strip() for value in (bundle_version, reason, provenance)):
        raise BundleError("publication requires a nonempty version, reason and provenance")
    workspace = workspace.expanduser().resolve()
    output = output.expanduser().absolute()
    if output.exists() or output.is_symlink():
        raise BundleError(f"output already exists: {output}")
    output = output.resolve()
    if output.is_relative_to(workspace):
        raise BundleError("golden output must be outside the refresh workspace")
    report = load_json(workspace / "refresh.json", "refresh report")
    if report.get("status") != "ready" or not report.get("evidence"):
        raise BundleError("publication requires a successful, validated refresh")
    for record in report["catalogs"]:
        _verify_record(catalog_root, record)
    case, producers, checks, _ = _recipe(report["case_id"], catalog_root)
    expected = [(item["workflow"], layout) for item in producers for layout in item["layouts"]]
    if [(item["workflow_id"], item["layout_id"]) for item in report["producers"]] != expected:
        raise BundleError("refresh is missing required producers")
    if [item["suite"] for item in report["checks"]] != checks:
        raise BundleError("refresh is missing required checks")
    for record in report["evidence"]:
        _verify_record(workspace, record)
    for item in report["producers"]:
        directory = Path(item["run_directory"])
        if not directory.resolve().is_relative_to(workspace) or item["status"] != "passed" or not reusable_outputs(directory):
            raise BundleError("producer failed or its outputs changed")
        if not producer_converged(item, case["workflows"][item["workflow_id"]], "producer", {}, catalog_root / "tolerances.json"):
            raise BundleError("producer did not converge")
    if any(item["status"] != "passed" for item in report.get("parallel_checks", [])):
        raise BundleError("producer parallel checks failed")
    for check in report["checks"]:
        summary = load_json(artifact_path(workspace, check["summary"], "check summary"), "check summary")
        if check["status"] != "passed" or summary["status"] != "passed":
            raise BundleError("golden validation checks failed")
        if not all(reusable_outputs(Path(item["run_directory"])) for item in summary["results"]):
            raise BundleError("validation run outputs changed")
    candidate = workspace / "candidate"
    validate_bundle_root(candidate, catalog_root / "cases")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output.name}.", dir=output.parent) as temporary:
        staging = Path(temporary) / "golden"
        shutil.copytree(candidate, staging)
        manifest = load_json(staging / "manifest.json", "candidate manifest")
        evidence = [workspace / "refresh.json", *[workspace / item["path"] for item in report["evidence"]
                                                 if not item["path"].startswith("candidate/")]]
        for number, source in enumerate(evidence):
            _install(staging, manifest, f"refresh_evidence_{number}", source,
                     "provenance/refresh/" + source.relative_to(workspace).as_posix(),
                     "application/json" if source.suffix == ".json" else "text/plain")
        publication = {"published_utc": utc_now(), "reason": reason.strip(), "provenance": provenance.strip(),
                       "case_id": case["case_id"], "bundle_version": bundle_version, "refresh_run_id": report["run_id"]}
        receipt = staging / "provenance/publication.json"
        write_json_atomic(receipt, publication, "publication provenance")
        manifest["artifacts"]["publication"] = {"path": "provenance/publication.json", **file_identity(receipt),
                                                "media_type": "application/json"}
        manifest.update(bundle_class="golden", bundle_version=bundle_version, created_utc=publication["published_utc"])
        write_json_atomic(staging / "manifest.json", manifest, "golden manifest")
        summary = validate_bundle_root(staging, catalog_root / "cases")
        staging.rename(output)
    return summary
