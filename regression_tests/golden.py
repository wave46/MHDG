"""Build, run ordered producers, collect once, review, then explicitly publish."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import shutil
import tempfile
import h5py

from .catalog import load_case_definition, workflow_required_roles, load_suite_definition, load_layouts, required_builds, required_case_roles, newton_probe_builds
from .documents import load_json, write_json_atomic
from .support import BundleError, HarnessError, utc_now, utc_run_id
from .files import file_identity
from . import config
from .build import build_solver
from .bundles import artifact_path, validate_bundle_root, preflight_bundles, load_manifest, series_artifact_ids
from .compare import compare_completed_run, compare_me_refinement_stop
from .compare_common import select_candidate
from .compare_time import accepted_checkpoints, accepted_checkpoint_files
from .execute import execute_prepared, reusable_outputs
from .prepare import prepare_run
from .compare_fixed import required_scalar
from .suites import run_suite, assess_result, reusable_result

ROOT = Path(__file__).resolve().parent


def _recipe(case_id, root, catalog):
    case = load_case_definition(case_id, root / "cases", catalog=catalog)
    document = load_json(root / "golden.json", "golden producers")
    declaration = document.get("cases", {}).get(case_id)
    if not declaration:
        raise BundleError(f"no golden producers declared for {case_id}")
    layouts = load_layouts(root / "layouts.json", catalog=catalog)
    producers = []
    for name in declaration["producers"]:
        if name not in case["workflows"]:
            raise BundleError(f"unknown producer workflow: {name}")
        workflow = case["workflows"][name]
        canonical = workflow["layout"]
        if canonical not in layouts:
            raise BundleError(f"unknown producer layout: {canonical}")
        roles = list(dict.fromkeys([*workflow.get("outputs", []), *workflow.get("checkpoint_outputs", {})]))
        producers.append({"workflow": name, "layout": canonical, "roles": roles})
    if not producers or len({item["workflow"] for item in producers}) != len(producers):
        raise BundleError("golden producers must be nonempty and unique")
    checks = [{"suite_id": name, **load_suite_definition(name, root / "suites.json", root / "layouts.json",
                                          root / "cases", case_id=case_id, catalog=catalog)}
              for name in declaration["checks"]]
    if not checks or len(set(declaration["checks"])) != len(checks):
        raise BundleError("golden checks must be nonempty and unique")
    requirements = set().union(*(required_builds(case, [item["workflow"]], [layouts[item["layout"]]])
                                 for item in producers))
    for check in checks:
        requirements.update(required_builds(case, check["workflow_ids"], [layouts[name] for name in check["layouts"]]))
    return case, producers, checks, requirements


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
    catalog = {}
    case, producers, checks, requirements = _recipe(case_id, root, catalog)
    source = Path(settings["MHDG_REGRESSION_DATA_ROOT"]).resolve()
    manifest = load_manifest(source, root / "cases", case_id=case_id)
    validate_bundle_root(source, root / "cases", manifest=manifest, catalog=catalog)
    available = {role for role, artifact in manifest["roles"].items()
                 if (source / manifest["artifacts"][artifact]["path"]).is_file()}
    # Outputs scheduled for replacement must come from this refresh, not old data.
    available.difference_update(role for producer in producers for role in producer["roles"])
    for producer in producers:
        workflow = case["workflows"][producer["workflow"]]
        needed = workflow_required_roles(workflow, references=False)
        missing = sorted(needed - available)
        if missing:
            raise BundleError(f"producer {producer['workflow']} needs: {', '.join(missing)}")
        available.update(producer["roles"])
    for check in checks:
        needed = set().union(*(required_case_roles(case, name, references=check["reference_comparisons"])
                               for name in check["workflow_ids"]))
        if needed - available:
            raise BundleError(f"check {check['suite_id']} needs: {', '.join(sorted(needed - available))}")
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
        build = build_solver(settings, ROOT.parent, jobs, requirements=requirements,
                             newton_probe=newton_probe_builds(checks, root / "cases", load_layouts(root / "layouts.json", catalog=catalog), catalog=catalog))
        values = {**settings, **config.build_settings(build.metadata_path),
                  "MHDG_REGRESSION_RUN_ROOT": str(workspace / "runs")}
        config.runtime_settings(values, requirements)
        shutil.copy2(build.metadata_path, workspace / "build.json")
        generated = {}
        mesh_sources = {}
        for producer in producers:
            layout = producer["layout"]
            name = producer["workflow"]
            print(f"producer: {name} / {layout}", flush=True)
            prepared = prepare_run(values, case_id, name, layout, root / "cases", root / "layouts.json",
                                   run_id, validate_bundle=False,
                                   requested_overrides={"balance_diagnostics_mode": "off"},
                                   artifact_overrides=generated, require_reference=False, catalog=catalog, manifest=manifest)
            result = {"workflow_id": name, "layout_id": layout, "run_directory": str(prepared.path),
                      "status": "running", "run_status": "running"}
            report["producers"].append(result)
            write_json_atomic(path, report, "refresh report")
            execution = execute_prepared(prepared, values)
            result.update(run_status=execution.status, status=execution.status)
            _validate_producer(result, case["workflows"][name], root, catalog)
            mesh_source = mesh_sources.get(case["workflows"][name].get("mesh"))
            if mesh_source:
                result["mesh_handoff"] = _check_mesh_handoff(result, mesh_source, root, catalog)
            result["status"] = "passed"
            result["old_reference"] = _old_reference(result, source, manifest, case, root, catalog)
            write_json_atomic(path, report, "refresh report")
            directory = Path(result["run_directory"])
            for role in producer["roles"]:
                checkpoints = case["workflows"][name].get("checkpoint_outputs", {})
                if role in checkpoints:
                    generated[role] = _checkpoint_solution(directory, checkpoints[role])
                    result.setdefault("generated_checkpoints", {})[role] = _record(generated[role], workspace)
                elif case["bundle_files"][role]["media_type"] == "application/x-gmsh":
                    try:
                        mesh = _generated_mesh(directory)
                    except HarnessError:
                        result["status"] = "failed"
                        raise
                    generated[role] = mesh
                    mesh_sources[role] = directory
                    result.setdefault("generated_meshes", {})[role] = _record(mesh, workspace)
                else:
                    generated[role] = _solution(directory, result["validation"].get("selected_output"))
            write_json_atomic(path, report, "refresh report")
        candidate = workspace / "candidate"
        _collect(source, candidate, manifest, generated, report, case)
        values["MHDG_REGRESSION_DATA_ROOT"] = str(candidate)
        candidate_manifest = preflight_bundles(checks, {case_id: values}, root / "cases",
                                               required_bundle_class="candidate", catalog=catalog)[case_id]
        for suite in checks:
            name = suite["suite_id"]
            summary_path, summary = run_suite(values, name, root / "cases", root / "layouts.json",
                                              root / "suites.json", root / "tolerances.json", f"{run_id}-verify",
                                              case_id=case_id, required_bundle_class="candidate",
                                              suite=suite, catalog=catalog, manifest=candidate_manifest)
            report["checks"].append({"suite": name, "status": summary["status"],
                                     "summary": str(summary_path.relative_to(workspace))})
            if summary["status"] != "passed":
                raise BundleError(f"golden validation failed: {name}")
        validate_bundle_root(candidate, root / "cases", catalog=catalog)
        report["evidence"] = _refresh_evidence(workspace)
        report["status"] = "ready"
    except (HarnessError, OSError, KeyboardInterrupt) as exc:
        if report["producers"] and report["producers"][-1]["status"] != "passed":
            report["producers"][-1]["status"] = "failed"
        report.update(status="failed", failure=str(exc) or "interrupted")
        write_json_atomic(path, report, "refresh report")
        raise
    write_json_atomic(path, report, "refresh report")
    return path, report


def _refresh_evidence(workspace):
    return [_record(item, workspace) for item in sorted(workspace.rglob("*"))
            if item.is_file() and not item.is_symlink() and item.name != "refresh.json"
            and (item.suffix in {".json", ".log", ".nml"} or item.name == "param.txt")]


def _solution(directory, selected_output=None):
    return select_candidate(directory, load_json(directory / "run_metadata.json", "run metadata"),
                            Path(selected_output) if selected_output else None)


def _checkpoint_solution(directory, selector):
    """Select a numbered checkpoint or the first accepted interpolated time."""
    import math

    metadata = load_json(directory / "run_metadata.json", "producer metadata")
    if selector != "first_intermediate":
        for path in accepted_checkpoint_files(directory, metadata):
            with h5py.File(path) as handle:
                number = required_scalar(handle, "simulation_parameters/time", "Current_time_step_number")
            if number == selector:
                return path
        raise BundleError(f"producer retained no accepted checkpoint at step {selector}: {directory}")
    plan = load_json(directory / "run_plan.json", "producer plan")
    frames = plan["moving_equilibrium_inputs"]["equilibrium_times_s"]
    initial = plan["moving_equilibrium_inputs"]["initial_step"]
    for path, clock in sorted(accepted_checkpoints(directory, metadata), key=lambda item: item[1][0]):
        number, seconds = clock[:2]
        selected = number > initial and not any(
            math.isclose(seconds, time, rel_tol=1e-10, abs_tol=1e-12) for time in frames)
        if selected:
            return path
    raise BundleError(f"producer retained no accepted checkpoint for {selector}: {directory}")


def _generated_mesh(directory):
    """Select the last retained adaptation mesh in recorded stage order."""
    metadata = load_json(directory / "run_metadata.json", "producer metadata")
    directories = [Path(stage["run_directory"]) for stage in metadata.get("stages", [])] or [directory]
    for stage in reversed(directories):
        mesh = stage / "res/temp.msh"
        if mesh.is_file():
            return mesh
    raise BundleError(f"mesh producer retained no res/temp.msh: {directory}")


def _check_mesh_handoff(result, source, root, catalog):
    """Compare final meshes and fields using the shared fixed-mesh comparator."""
    directory = Path(result["run_directory"])
    reference = _solution(source)
    _, report_path, comparison = compare_completed_run(
        directory, root / "cases", root / "tolerances.json", reference_override=reference,
        comparison_policy_override="fixed_hdf5", report_override=directory / "adaptive_endpoint.json", catalog=catalog,
    )
    if comparison["status"] != "passed":
        raise BundleError(f"fixed/adaptive bootstrap endpoint comparison failed: {report_path}")
    return {"source_run": str(source),
            "comparison_report": str(report_path), "status": comparison["status"]}


def _validate_producer(result, workflow, root, catalog):
    directory = Path(result["run_directory"])
    expectation = workflow.get("moving_equilibrium", {}).get("expected_failure")
    if not expectation and (result["run_status"] != "completed" or not reusable_outputs(directory)):
        raise BundleError(f"producer failed or its outputs changed: {result['workflow_id']}")
    assessed = assess_result(result, root / "cases", root / "tolerances.json", references=False,
                             expected_failure=expectation, catalog=catalog)
    result["validation"] = assessed.get("validation", {})
    if assessed["status"] != "passed":
        raise BundleError(f"invalid producer output: {result['workflow_id']}: {assessed['failures']}")
    if expectation:
        result["expected_failure"] = assessed["expected_failure"]


def _old_reference(result, source, manifest, case, root, catalog):
    directory = Path(result["run_directory"])
    workflow = case["workflows"][result["workflow_id"]]
    artifact = manifest["roles"].get(workflow.get("reference"))
    reference = source / manifest["artifacts"][artifact]["path"] if artifact else None
    if reference is None or not reference.is_file():
        return {"status": "unavailable"}
    try:
        expectation = workflow.get("moving_equilibrium", {}).get("expected_failure")
        compare = compare_me_refinement_stop if expectation else compare_completed_run
        _, path, report = compare(directory, root / "cases", root / "tolerances.json",
                                reference_override=reference, report_override=directory / "old_reference.json",
                                catalog=catalog, **({"completed_steps": expectation["completed_steps"]} if expectation else {}))
        return {"status": report["status"], "report": str(path)}
    except HarnessError as exc:
        return {"status": "unavailable", "reason": str(exc)}


def _install(staging, manifest, artifact_id, source, relative, media_type):
    target = staging / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    manifest["artifacts"][artifact_id] = {"path": relative, **file_identity(target), "media_type": media_type}


def _collect(source, candidate, source_manifest, generated, report, case):
    """Collect one candidate; later publication copies this validated bundle."""
    candidate.mkdir()
    manifest = deepcopy(source_manifest)
    # A complete refresh carries only current case inputs and newly produced
    # references. Retired roles/stage matrices stay in the historical source.
    manifest["roles"] = {role: artifact for role, artifact in manifest["roles"].items()
                         if role in case["bundle_files"]}
    keep = set(manifest["roles"].values()) | series_artifact_ids(source, manifest, case)
    manifest["artifacts"] = {key: value for key, value in manifest["artifacts"].items() if key in keep}
    for role, solution in generated.items():
        spec = case["bundle_files"][role]
        artifact_id = manifest["roles"].get(role, spec["artifact_id"])
        _install(candidate, manifest, artifact_id, solution, f"inputs/{spec['filename']}", spec["media_type"])
        manifest["roles"][role] = artifact_id
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


def _publication_summary(workspace, report):
    """Anchor retained evidence without copying machine paths or campaign logs."""
    build = load_json(workspace / "build.json", "refresh build")
    return {
        "refresh_report": file_identity(workspace / "refresh.json"),
        "candidate_manifest": file_identity(workspace / "candidate/manifest.json"),
        "build": {"build_id": build["build_id"],
                  "repository": {key: build["repository"][key] for key in ("revision", "dirty")
                                 if key in build["repository"]},
                  "profile": build["profile"], "toolchain": build.get("toolchain", {}),
                  "artifacts": {name: {key: item[key] for key in ("sha256", "size_bytes")}
                                for name, item in build["artifacts"].items()},
                  "record": file_identity(workspace / "build.json")},
        "source_bundle": {"bundle_id": report["source"]["bundle_id"],
                          "bundle_version": report["source"]["bundle_version"],
                          "manifest": {key: report["source"]["manifest"][key]
                                       for key in ("sha256", "size_bytes")}},
        "catalogs": report["catalogs"],
        "producers": [{"workflow": item["workflow_id"], "layout": item["layout_id"],
                       "status": item["status"],
                       "stages": [{"stage": stage["stage_id"], "status": stage["status"],
                                   "convergence": stage["convergence"]}
                                  for stage in item["validation"].get("stages", [])],
                       **({"run_status": item["run_status"], "expected_failure": {
                           "reason": item["validation"]["terminal_reason"],
                           "completed_steps": item["validation"]["completed_prefix"]["completed_steps"],
                           "last_attempted_step": item["validation"]["last_attempted_step"]}}
                          if item.get("expected_failure") else {})}
                      for item in report["producers"]],
        "validation_suites": [{"suite": item["suite"], "status": item["status"]}
                              for item in report["checks"]],
        **({"continuations": [_continuation_summary(workspace, item) for item in report["continuations"]]}
           if report.get("continuations") else {}),
    }


def _continuation_summary(workspace, item):
    if "build" not in item:
        return item
    build = load_json(_verify_record(workspace, item["build"]), "continuation build")
    return {**item, "validation_build": {
        "build_id": build["build_id"],
        "repository": {key: build["repository"][key] for key in ("revision", "dirty")
                       if key in build["repository"]},
        "artifacts": {name: {key: value[key] for key in ("sha256", "size_bytes")}
                      for name, value in build["artifacts"].items()}}}


def publish(workspace, output, bundle_version, reason, provenance, *, catalog_root=None):
    """Revalidate full evidence; publish runnable files and acceptance provenance."""
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
    catalog = {}
    case, producers, checks, _ = _recipe(report["case_id"], catalog_root, catalog)
    expected = [(item["workflow"], item["layout"]) for item in producers]
    if [(item["workflow_id"], item["layout_id"]) for item in report["producers"]] != expected:
        raise BundleError("refresh is missing required producers")
    if [item["suite"] for item in report["checks"]] != [check["suite_id"] for check in checks]:
        raise BundleError("refresh is missing required checks")
    for record in report["evidence"]:
        _verify_record(workspace, record)
    for item in report["producers"]:
        directory = Path(item["run_directory"])
        expectation = case["workflows"][item["workflow_id"]].get("moving_equilibrium", {}).get("expected_failure")
        if (not directory.resolve().is_relative_to(workspace) or item["status"] != "passed"
                or not reusable_result(item, catalog_root / "cases", catalog_root / "tolerances.json",
                                        expected_failure=expectation, catalog=catalog)):
            raise BundleError("producer failed or its outputs changed")
        validation = item.get("validation", {})
        if validation.get("status") != "passed" or validation.get("convergence", {}).get("passed") is not True:
            raise BundleError("producer has no successful recorded validation")
        for record in item.get("generated_meshes", {}).values():
            _verify_record(workspace, record)
        for record in item.get("generated_checkpoints", {}).values():
            _verify_record(workspace, record)
    for check in report["checks"]:
        summary = load_json(artifact_path(workspace, check["summary"], "check summary"), "check summary")
        if check["status"] != "passed" or summary["status"] != "passed":
            raise BundleError("golden validation checks failed")
        if not all(reusable_result(item, catalog_root / "cases", catalog_root / "tolerances.json", catalog=catalog,
                                    expected_failure=case["workflows"][item["workflow_id"]].get("moving_equilibrium", {}).get("expected_failure"))
                   for item in summary["results"]):
            raise BundleError("validation run outputs changed")
    candidate = workspace / "candidate"
    validate_bundle_root(candidate, catalog_root / "cases", catalog=catalog)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{output.name}.", dir=output.parent) as temporary:
        staging = Path(temporary) / "golden"
        manifest = load_json(candidate / "manifest.json", "candidate manifest")
        staging.mkdir()
        keep = set(manifest["roles"].values()) | series_artifact_ids(candidate, manifest, case)
        manifest["artifacts"] = {key: value for key, value in manifest["artifacts"].items() if key in keep}
        for artifact in manifest["artifacts"].values():
            source = candidate / artifact["path"]
            if artifact.get("optional") and not source.is_file():
                continue
            target = staging / artifact["path"]
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        publication = {"published_utc": utc_now(), "reason": reason.strip(), "provenance": provenance.strip(),
                       "case_id": case["case_id"], "bundle_version": bundle_version, "refresh_run_id": report["run_id"]}
        publication.update(_publication_summary(workspace, report))
        receipt = staging / "provenance/publication.json"
        write_json_atomic(receipt, publication, "publication provenance")
        manifest["artifacts"]["publication"] = {"path": "provenance/publication.json", **file_identity(receipt),
                                                "media_type": "application/json"}
        manifest.update(bundle_class="golden", bundle_version=bundle_version, created_utc=publication["published_utc"])
        write_json_atomic(staging / "manifest.json", manifest, "golden manifest")
        summary = validate_bundle_root(staging, catalog_root / "cases", catalog=catalog)
        staging.rename(output)
    return summary
