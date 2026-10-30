"""Continue candidate validation while retaining successful producer evidence."""

from copy import deepcopy
from pathlib import Path
import shutil

from . import config, golden
from .build import build_solver
from .bundles import artifact_path, load_manifest, preflight_bundles, validate_bundle_root
from .catalog import required_case_roles, newton_probe_builds, load_layouts
from .documents import load_json, write_json_atomic
from .files import file_identity
from .support import BundleError, HarnessError, utc_now, utc_run_id
from .suites import run_suite


def validate_saved_campaign(workspace, root, rerun_workflows):
    """Reject changed recipes, producer outputs or candidate reference artifacts."""
    report = load_json(workspace / "refresh.json", "refresh report")
    if report.get("status") != "failed":
        raise BundleError("continuation requires a failed refresh")
    for record in report["catalogs"]:
        golden._verify_record(root, record)
    catalog = {}
    case, producers, checks, requirements = golden._recipe(report["case_id"], root, catalog)
    selected = {name for check in checks for name in check["workflow_ids"]}
    if not rerun_workflows or set(rerun_workflows)-selected:
        raise BundleError("continuation requires explicit --rerun validation workflows from this campaign")
    expected = [(item["workflow"], item["layout"]) for item in producers]
    if [(item["workflow_id"], item["layout_id"]) for item in report["producers"]] != expected:
        raise BundleError("continuation requires all producers to have finished")
    source_record = report["source"]["manifest"]
    if file_identity(Path(source_record["path"])) != {key: source_record[key] for key in ("sha256", "size_bytes")}:
        raise BundleError("source bundle manifest changed")
    candidate = workspace / "candidate"
    manifest = load_manifest(candidate, root / "cases", case_id=report["case_id"])
    required = set().union(*(required_case_roles(case, name, references=check["reference_comparisons"])
                             for check in checks for name in check["workflow_ids"]))
    validate_bundle_root(candidate, root / "cases", manifest=manifest, catalog=catalog, required_roles=required)
    source = Path(report["source"]["manifest"]["path"]).parent
    source_manifest = load_manifest(source, root / "cases", case_id=report["case_id"])
    generated = {role for producer in producers for role in producer["roles"]}
    for role, artifact in manifest["roles"].items():
        if role not in generated:
            original = source_manifest["roles"].get(role)
            if original is None or file_identity(source / source_manifest["artifacts"][original]["path"]) != file_identity(
                    candidate / manifest["artifacts"][artifact]["path"]):
                raise BundleError(f"candidate input differs from source: {role}")
    for result, producer in zip(report["producers"], producers):
        directory = Path(result["run_directory"])
        if result["status"] != "passed" or not directory.resolve().is_relative_to(workspace):
            raise BundleError("continuation requires successful producers inside the refresh workspace")
        metadata = load_json(directory / "run_metadata.json", "producer metadata")
        config.build_settings(Path(metadata["solver"]["build_manifest"]["path"]))
        for record in result.get("generated_checkpoints", {}).values():
            golden._verify_record(workspace, record)
        for record in result.get("generated_meshes", {}).values():
            golden._verify_record(workspace, record)
        workflow = case["workflows"][producer["workflow"]]
        golden._validate_producer(result, workflow, root, catalog)
        for role in producer["roles"]:
            if role in workflow.get("checkpoint_outputs", {}):
                output = golden._checkpoint_solution(directory, workflow["checkpoint_outputs"][role])
            elif case["bundle_files"][role]["media_type"] == "application/x-gmsh":
                output = golden._generated_mesh(directory)
            else:
                output = golden._solution(directory, result["validation"].get("selected_output"))
            collected = candidate / manifest["artifacts"][manifest["roles"][role]]["path"]
            if file_identity(output) != file_identity(collected):
                raise BundleError(f"candidate reference differs from producer: {role}")
    for check in report["checks"]:
        artifact_path(workspace, check["summary"], "previous check summary")
    return report, checks, requirements, catalog


def continue_refresh(workspace, settings, rerun_workflows, jobs=None, *, catalog_root=None):
    """Revalidate completed producers, then reuse or rerun candidate checks."""
    root = catalog_root or golden.ROOT
    workspace = workspace.expanduser().resolve()
    report, checks, requirements, catalog = validate_saved_campaign(workspace, root, rerun_workflows)
    previous_checks = {}
    for attempt in report.get("continuations", []):
        previous = load_json(golden._verify_record(workspace, attempt["previous_report"]), "previous refresh report")
        previous_checks.update((item["suite"], item) for item in previous["checks"])
    previous_checks.update((item["suite"], item) for item in report["checks"])
    run_id = utc_run_id()
    directory = workspace / "continuations" / run_id
    directory.mkdir(parents=True)
    snapshot = directory / "previous_refresh.json"
    write_json_atomic(snapshot, load_json(workspace / "refresh.json", "refresh report"), "previous refresh report")
    continuation = {"run_id": run_id, "status": "running", "created_utc": utc_now(),
                    "previous_report": golden._record(snapshot, workspace),
                    "rerun_workflows": sorted(set(rerun_workflows)),
                    "reused_producers": [item["workflow_id"] for item in report["producers"]]}
    report = deepcopy(report)
    report.setdefault("continuations", []).append(continuation)
    report.update(status="running", checks=[])
    report.pop("failure", None)
    report.pop("evidence", None)
    path = workspace / "refresh.json"
    write_json_atomic(path, report, "refresh report")
    try:
        build = build_solver(settings, root.parent, jobs, requirements=requirements,
                             newton_probe=newton_probe_builds(checks, root / "cases", load_layouts(root / "layouts.json", catalog=catalog), catalog=catalog))
        shutil.copy2(build.metadata_path, directory / "build.json")
        continuation["build"] = golden._record(directory / "build.json", workspace)
        values = {**settings, **config.build_settings(build.metadata_path),
                  "MHDG_REGRESSION_RUN_ROOT": str(workspace / "runs"),
                  "MHDG_REGRESSION_DATA_ROOT": str(workspace / "candidate")}
        config.runtime_settings(values, requirements)
        candidate = preflight_bundles(checks, {report["case_id"]: values}, root / "cases",
                                      required_bundle_class="candidate", catalog=catalog)[report["case_id"]]
        for suite in checks:
            name = suite["suite_id"]
            previous = previous_checks.get(name)
            summary_path, summary = run_suite(
                values, name, root / "cases", root / "layouts.json", root / "suites.json",
                root / "tolerances.json", f"{run_id}-continue", case_id=report["case_id"],
                required_bundle_class="candidate", suite=suite, catalog=catalog, manifest=candidate,
                reuse_from=workspace / previous["summary"] if previous else None,
                rerun_workflows=rerun_workflows)
            report["checks"].append({"suite": name, "status": summary["status"],
                                     "summary": str(summary_path.relative_to(workspace))})
            write_json_atomic(path, report, "refresh report")
            if summary["status"] != "passed":
                raise BundleError(f"golden validation failed: {name}")
        validate_bundle_root(workspace / "candidate", root / "cases", catalog=catalog)
        report["evidence"] = golden._refresh_evidence(workspace)
        report["status"] = "ready"
        continuation["status"] = "passed"
    except (HarnessError, OSError, KeyboardInterrupt) as exc:
        continuation["status"] = "failed"
        report.update(status="failed", failure=str(exc) or "interrupted")
        write_json_atomic(path, report, "refresh report")
        raise
    write_json_atomic(path, report, "refresh report")
    return path, report
