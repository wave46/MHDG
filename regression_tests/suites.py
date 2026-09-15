"""Run and verify scientific suites; resume reuses completed, unchanged outputs."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import time
from typing import Any

from .catalog import load_case_definition, load_suite_definition, load_layouts, required_builds, selection_builds
from .config import (bundle_root_from_settings, absolute_setting, solver_executable,
                     runtime_files, selected_mpi_launcher, runtime_settings, build_key)
from .bundles import preflight_bundles
from .prepare import prepare_run
from .execute import execute_prepared, reusable_outputs, execution_failures
from .compare import compare_completed_run, compare_generated_meshes, validate_completed_run
from .compare_common import select_candidate
from .documents import load_json, write_json_atomic
from .support import BundleError, ComparisonError, HarnessError, IDENTIFIER_RE, utc_now, utc_run_id
from .files import file_identity, require_file, recorded_directory


@dataclass(frozen=True)
class SuiteRunInputs:
    settings: dict[str, Any]
    case_id: str
    run_id: str
    case_directory: Path
    layouts_path: Path
    tolerances_path: Path
    compare: bool
    parameter_overrides: dict[str, bool | float | int | str]
    catalog: dict | None = None
    manifest: dict | None = None
    require_reference: bool = True


def run_cell(
    inputs: SuiteRunInputs,
    workflow_id: str,
    layout_id: str,
) -> dict[str, Any]:
    """Prepare, execute, and optionally compare one suite cell."""
    result = _empty_result(workflow_id, layout_id)
    try:
        prepared = prepare_run(
            inputs.settings,
            inputs.case_id,
            workflow_id,
            layout_id,
            inputs.case_directory,
            inputs.layouts_path,
            inputs.run_id,
            validate_bundle=False,
            requested_overrides=inputs.parameter_overrides,
            catalog=inputs.catalog, manifest=inputs.manifest, require_reference=inputs.require_reference,
        )
        result["run_directory"] = str(prepared.path)
        run = execute_prepared(prepared, inputs.settings)
        result["run_status"] = run.status
        result["duration_seconds"] = run.duration_seconds
        result["status"] = run.status
        return assess_result(result, inputs.case_directory, inputs.tolerances_path,
                             references=inputs.require_reference, enabled=inputs.compare, catalog=inputs.catalog)
    except HarnessError as exc:
        result["failures"] = [str(exc)]
    return result


def _empty_result(workflow_id: str, layout_id: str) -> dict[str, Any]:
    return {
        "workflow_id": workflow_id,
        "layout_id": layout_id,
        "status": "error",
        "run_status": None,
        "comparison_policy": None,
        "comparison_status": "not_run",
        "convergence_status": None,
        "duration_seconds": None,
        "run_directory": None,
        "comparison_report": None,
        "failures": [],
    }


def compare_layout_pairs(
    summary: dict[str, Any],
    case_directory: Path,
    tolerances_path: Path,
    *, catalog=None,
) -> list[dict[str, Any]]:
    """Compare every declared candidate with its same-workflow baseline."""
    catalog = {} if catalog is None else catalog
    return [
        _compare_pair(
            summary,
            workflow_id,
            pair,
            case_directory,
            tolerances_path,
            catalog,
        )
        for workflow_id in summary["workflow_ids"]
        for pair in summary["layout_comparisons"]
    ]


def _compare_pair(
    summary: dict[str, Any],
    workflow_id: str,
    pair: dict[str, str],
    case_directory: Path,
    tolerances_path: Path,
    catalog,
) -> dict[str, Any]:
    baseline_layout = pair["baseline"]
    candidate_layout = pair["candidate"]
    result = {
        "workflow_id": workflow_id,
        "baseline_layout_id": baseline_layout,
        "candidate_layout_id": candidate_layout,
        "baseline_run_directory": None,
        "candidate_run_directory": None,
        "baseline_output": None,
        "comparison_policy": None,
        "comparison_report": None,
        "generated_meshes": None,
        "status": "failed",
        "failures": [],
    }
    try:
        baseline = _completed_run(summary, workflow_id, baseline_layout)
        candidate = _completed_run(summary, workflow_id, candidate_layout)
        result["baseline_run_directory"] = str(baseline)
        result["candidate_run_directory"] = str(candidate)
        reference = _selected_output(baseline)
        result["baseline_output"] = str(reference)
        policy, report_path, report = compare_completed_run(
            candidate,
            case_directory,
            tolerances_path,
            catalog=catalog,
            reference_override=reference,
            tolerance_profile_override=summary["tolerance_profile"],
            comparison_policy_override=summary.get("layout_comparison_policy"),
            report_override=candidate / f"comparison_from_{baseline_layout}.json",
        )
        result["comparison_policy"] = policy
        result["comparison_report"] = str(report_path)
        result["status"] = report["status"]
        result["failures"] = list(report["failures"])
        generated_meshes = compare_generated_meshes(baseline, candidate)
        result["generated_meshes"] = generated_meshes
        if generated_meshes and not generated_meshes["passed"]:
            result["status"] = "failed"
            result["failures"].extend(generated_meshes["failures"])
        if result["status"] == "passed":
            from .diagnostics import compare_outputs

            diagnostic_report = compare_outputs(reference, _selected_output(candidate))
            result["diagnostics"] = diagnostic_report
            if diagnostic_report["status"] == "failed":
                result["status"] = "failed"
                result["failures"].extend(diagnostic_report["failures"])
    except HarnessError as exc:
        result["failures"] = [str(exc)]
    return result


def _completed_run(
    summary: dict[str, Any],
    workflow_id: str,
    layout_id: str,
) -> Path:
    matches = [
        result
        for result in summary["results"]
        if result.get("workflow_id") == workflow_id
        and result.get("layout_id") == layout_id
    ]
    if len(matches) != 1:
        raise ComparisonError(
            f"suite has no unique result for {workflow_id}/{layout_id}"
        )
    result = matches[0]
    if result.get("run_status") != "completed":
        raise ComparisonError(
            f"suite run did not complete: {workflow_id}/{layout_id}"
        )
    return recorded_directory(result.get("run_directory"), "suite run")


def _selected_output(run_directory: Path) -> Path:
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    return select_candidate(run_directory, metadata)


def assess_result(source, case_directory, tolerances_path, *, references=True, enabled=True, catalog=None):
    """Assess fresh, resumed or saved execution with the same checks and result fields."""
    result = {**source, "status": "failed", "failures": [], "comparison_policy": None,
              "comparison_report": None, "comparison_status": "not_run", "convergence_status": None}
    result.pop("validation", None)
    try:
        directory = recorded_directory(source.get("run_directory"), "suite run")
        if source.get("run_status") != "completed":
            result["failures"] = [f"solver run status is {source.get('run_status')}", *execution_failures(directory)]
            return result
        if not enabled:
            result["status"] = "deferred"
            return result
        if references:
            policy, path, report = compare_completed_run(directory, case_directory, tolerances_path, catalog=catalog)
            result.update(comparison_policy=policy, comparison_report=str(path), comparison_status=report["status"])
        else:
            report = validate_completed_run(directory, case_directory, tolerances_path, catalog=catalog)
            result["validation"] = report
        convergence = report["convergence"]["passed"]
        result["convergence_status"] = "passed" if convergence is True else "failed" if convergence is False else "not_checked"
        result["failures"] = list(report["failures"])
        if convergence is False and not result["failures"]:
            result["failures"].append("Newton convergence check failed")
        result["status"] = "passed" if report["status"] == "passed" and convergence is True else "failed"
    except (HarnessError, TypeError) as exc:
        result["failures"] = [str(exc)]
    return result


def _validate_source_summary(summary: dict[str, Any]) -> None:
    required = {"schema_version", "suite_id", "run_id", "case_id", "results"}
    missing = sorted(required - summary.keys())
    if missing:
        raise BundleError(f"suite summary is missing: {', '.join(missing)}")
    if summary["schema_version"] != 2:
        raise BundleError("suite summary has unsupported schema version")
    if not isinstance(summary["results"], list) or not summary["results"]:
        raise BundleError("suite summary contains no results")
    result_fields = {"workflow_id", "layout_id", "run_directory", "run_status"}
    for result in summary["results"]:
        if not isinstance(result, dict) or not result_fields <= result.keys():
            raise BundleError("suite summary contains an invalid result")
    if "layout_comparisons" in summary:
        pair_fields = {"workflow_ids", "tolerance_profile"}
        if not pair_fields <= summary.keys():
            raise BundleError("paired suite summary is incomplete")
        pairs = summary["layout_comparisons"]
        if not isinstance(pairs, list) or not pairs:
            raise BundleError("paired suite summary has no layout comparisons")
        if any(
            not isinstance(pair, dict)
            or set(pair) != {"baseline", "candidate"}
            for pair in pairs
        ):
            raise BundleError("paired suite summary has an invalid comparison")


def suite_execution_inputs(
    settings: dict[str, Any],
    bundle_root: Path,
    requirements,
    tracked_files: dict[str, Path],
) -> dict[str, Any]:
    """Describe files that must remain stable across suite resumes."""
    records = {
        "bundle_manifest": _file_record(
            bundle_root / "manifest.json",
            "bundle manifest",
        ),
        **{
            name: _file_record(path, name.replace("_", " "))
            for name, path in tracked_files.items()
        },
    }
    records["executables"] = {
        build_key(model, execution): _file_record(solver_executable(settings, execution, model), "executable")
        for model, execution in sorted(requirements)
    }
    if any(execution == "mpi" for _, execution in requirements):
        records["mpi_launcher"] = _file_record(selected_mpi_launcher(settings), "MPI launcher")
    records["runtime_files"] = {
        str(path): file_identity(path)
        for item in records["executables"].values()
        for path in runtime_files(Path(item["path"])).values()
    }
    for key, name in (
        ("MHDG_ENVIRONMENT_SCRIPT", "environment_script"),
        ("MHDG_BUILD_MANIFEST", "build_manifest"),
    ):
        if settings.get(key):
            records[name] = _file_record(absolute_setting(settings, key), key)
    return records


def _file_record(path: Path, label: str) -> dict[str, Any]:
    path = require_file(path, label)
    return {"path": str(path), **file_identity(path)}


def _validate_execution_inputs(recorded: Any, current: dict[str, Any]) -> None:
    if not isinstance(recorded, dict):
        raise BundleError("suite summary has no execution input identity")
    changed = sorted(
        name
        for name in recorded.keys() | current.keys()
        if recorded.get(name) != current.get(name)
    )
    if changed:
        raise BundleError(f"suite execution inputs changed: {', '.join(changed)}")


def _validate_recorded_cells(
    results: list[Any],
    workflow_ids: list[str],
    layout_ids: list[str],
) -> None:
    if not isinstance(results, list):
        raise BundleError("suite summary has invalid results")
    planned = {
        (workflow, layout)
        for workflow in workflow_ids
        for layout in layout_ids
    }
    recorded = [
        (result.get("workflow_id"), result.get("layout_id"))
        for result in results
        if isinstance(result, dict)
    ]
    if len(recorded) != len(results) or any(cell not in planned for cell in recorded):
        raise BundleError("suite summary contains an invalid workflow/layout cell")
    if len(set(recorded)) != len(recorded):
        raise BundleError("suite summary contains duplicate workflow/layout cells")


def run_suite(
    settings, suite_id, case_directory, layouts_path, suites_path,
    tolerances_path, run_id=None, required_bundle_class=None, compare=True,
    resume=False, parameter_overrides=None, *, case_id=None, catalog=None, suite=None, manifest=None,
):
    catalog = {} if catalog is None else catalog
    suite = suite or load_suite_definition(suite_id, suites_path, layouts_path, case_directory, case_id=case_id, catalog=catalog)
    bundle_root = bundle_root_from_settings(settings)
    run_id = run_id or utc_run_id()
    if not IDENTIFIER_RE.fullmatch(run_id):
        raise BundleError(f"invalid suite run identifier: {run_id}")
    pairs = suite.get("layout_comparisons")
    references = suite["reference_comparisons"]
    overrides = {"balance_diagnostics_mode": suite["diagnostics"], **(parameter_overrides or {})}
    case = load_case_definition(suite["case_id"], case_directory, catalog=catalog)
    layouts = load_layouts(layouts_path, catalog=catalog)
    requirements = required_builds(case, suite["workflow_ids"], [layouts[name] for name in suite["layouts"]])
    # Profile/CLI callers pass the manifest after preflighting the entire selection.
    if manifest is None:
        manifest = preflight_bundles([suite], {suite["case_id"]: settings}, case_directory,
                                     required_bundle_class=required_bundle_class, catalog=catalog)[suite["case_id"]]
        runtime_settings(settings, requirements)
    identity = suite_execution_inputs(settings, bundle_root, requirements, {
        "case_definition": case_directory / f"{suite['case_id']}.json",
        "workflow_catalog": case_directory.parent / "workflows.json",
        "layout_catalog": layouts_path, "suite_catalog": suites_path,
        "tolerance_catalog": tolerances_path,
    })
    root = absolute_setting(settings, "MHDG_REGRESSION_RUN_ROOT")
    directory = root.resolve() / "suites" / suite_id / suite["case_id"] / run_id
    path = directory / "suite_summary.json"
    expected = {
        "suite_id": suite_id, "run_id": run_id, "case_id": suite["case_id"],
        "workflow_ids": suite["workflow_ids"], "layout_ids": suite["layouts"],
        "checks_enabled": compare, "parameter_overrides": overrides,
        "reference_comparisons": references,
        **{key: suite[key] for key in ("layout_comparisons", "tolerance_profile", "layout_comparison_policy") if key in suite},
    }
    if resume:
        summary = load_json(path, "suite summary")
        changed = [key for key, value in expected.items() if summary.get(key) != value]
        if changed:
            raise BundleError(f"suite summary does not match: {', '.join(changed)}")
        _validate_execution_inputs(summary.get("execution_inputs"), identity)
        _validate_recorded_cells(summary.get("results", []), suite["workflow_ids"], suite["layouts"])
    else:
        try:
            directory.mkdir(parents=True)
        except FileExistsError as exc:
            raise BundleError(f"suite summary directory already exists: {directory}") from exc
        summary = {
            "schema_version": 2, "started_utc": utc_now(), "duration_seconds": 0.,
            "description": suite["description"], "execution_inputs": identity,
            **expected, "results": [],
        }
    summary.update(status="running", finished_utc=None, comparisons=[])
    summary.pop("diagnostics", None)
    write_json_atomic(path, summary, "suite summary")
    inputs = SuiteRunInputs(settings, suite["case_id"], run_id, case_directory,
                            layouts_path, tolerances_path, compare, overrides,
                            catalog=catalog, manifest=manifest, require_reference=references)
    print(f"suite: {suite_id} ({run_id})")
    for layout in suite["layouts"]:
        for workflow in suite["workflow_ids"]:
            started = time.monotonic()
            previous = next((item for item in summary["results"]
                             if (item["workflow_id"], item["layout_id"]) == (workflow, layout)), None)
            if (previous and previous.get("run_status") == "completed"
                    and previous.get("run_directory") and reusable_outputs(Path(previous["run_directory"]))):
                print(f"reusing completed {workflow} / {layout}", flush=True)
                result = assess_result(previous, case_directory, tolerances_path,
                                       references=references, enabled=compare, catalog=catalog)
            else:
                selected = inputs
                existing = root / suite["case_id"] / workflow / layout / run_id
                if existing.exists():
                    if not resume:
                        raise BundleError(f"run directory already exists: {existing}")
                    attempt = 1
                    while (existing.parent / f"{run_id}-resume-{attempt}").exists():
                        attempt += 1
                    selected = replace(inputs, run_id=f"{run_id}-resume-{attempt}")
                    print(f"preserving incomplete or changed run {existing}; restarting workflow as {selected.run_id}")
                print(f"running {workflow} / {layout} ...", flush=True)
                result = run_cell(selected, workflow, layout)
            if previous is not None:
                summary["results"][summary["results"].index(previous)] = result
            else:
                summary["results"].append(result)
            summary["duration_seconds"] += time.monotonic() - started
            write_json_atomic(path, summary, "suite summary")
    if pairs and compare:
        summary["comparisons"] = compare_layout_pairs(summary, case_directory, tolerances_path, catalog=catalog)
    if compare:
        from .diagnostics import check_suite

        summary["diagnostics"] = check_suite(summary, required=False)
    _finish(summary)
    write_json_atomic(path, summary, "suite summary")
    return path, summary


def _finish(summary):
    checks = [*summary["results"], *summary.get("comparisons", [])]
    if summary.get("diagnostics"):
        checks.append(summary["diagnostics"])
    statuses = {item["status"] for item in checks}
    status = "failed" if not statuses or statuses - {"passed", "deferred"} else (
        "deferred" if "deferred" in statuses else "passed")
    summary.update(status=status, finished_utc=utc_now())


def verify_suite(suite_summary_path, case_directory, tolerances_path, *, include_layout_pairs=True):
    catalog = {}
    path = require_file(suite_summary_path, "suite summary")
    source = load_json(path, "suite summary")
    _validate_source_summary(source)
    case = load_case_definition(source["case_id"], case_directory, catalog=catalog)
    if any(item["workflow_id"] not in case["workflows"] for item in source["results"]):
        raise BundleError("suite summary refers to an unknown workflow")
    output = path.parent / "verification_summary.json"
    summary = {
        "schema_version": 2, "created_utc": utc_now(), "status": "running",
        "source_summary": str(path), "suite_id": source["suite_id"],
        "run_id": source["run_id"], "case_id": source["case_id"], "results": [],
    }
    write_json_atomic(output, summary, "verification summary")
    references = source.get("reference_comparisons", not source.get("layout_comparisons"))
    for result in source["results"]:
        summary["results"].append(assess_result(result, case_directory, tolerances_path, references=references, catalog=catalog))
        write_json_atomic(output, summary, "verification summary")
    if include_layout_pairs and source.get("layout_comparisons"):
        summary["comparisons"] = compare_layout_pairs(source, case_directory, tolerances_path, catalog=catalog)
    from .diagnostics import check_suite

    summary["diagnostics"] = check_suite(source, required=False)
    _finish(summary)
    write_json_atomic(output, summary, "verification summary")
    return output, summary


def run_profile(name, checks, settings_by_case, catalog_root, run_id=None, *, resume=False,
                compare=True, required_bundle_class=None, parameter_overrides=None, catalog=None, manifests=None):
    """Run an ordered selection using ordinary suite summaries for resume."""
    from .reporting import print_run_summary

    catalog = {} if catalog is None else catalog
    run_id = run_id or utc_run_id()
    if not IDENTIFIER_RE.fullmatch(run_id):
        raise BundleError(f"invalid profile run identifier: {run_id}")
    roots = {Path(values["MHDG_REGRESSION_RUN_ROOT"]) for values in settings_by_case.values()}
    if len(roots) != 1 or not next(iter(roots)).is_absolute():
        raise BundleError("a profile needs one absolute run root")
    root = next(iter(roots))
    path = root / "profiles" / name / run_id / "profile_summary.json"
    if path.exists() and not resume:
        raise BundleError(f"profile summary already exists: {path}")
    if resume and not path.is_file():
        raise BundleError(f"profile summary not found: {path}")
    if manifests is None:
        manifests = preflight_bundles(checks, settings_by_case, catalog_root / "cases",
                                      required_bundle_class=required_bundle_class, catalog=catalog)
        layouts = load_layouts(catalog_root / "layouts.json", catalog=catalog)
        for case_id, values in settings_by_case.items():
            selected = [check for check in checks if check["case_id"] == case_id]
            runtime_settings(values, selection_builds(selected, catalog_root / "cases", layouts, catalog=catalog))
    summary = {"profile": name, "run_id": run_id, "status": "running", "results": []}
    if resume:
        previous = load_json(path, "profile summary")
        if previous.get("profile") != name or previous.get("run_id") != run_id:
            raise BundleError("profile summary does not match the selection")
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json_atomic(path, summary, "profile summary")
    started = time.monotonic()
    for check in checks:
        suite, case = check["suite_id"], check["case_id"]
        existing = root / "suites" / suite / case / run_id / "suite_summary.json"
        result_path, result = run_suite(
            settings_by_case[case], suite, catalog_root / "cases", catalog_root / "layouts.json",
            catalog_root / "suites.json", catalog_root / "tolerances.json", run_id,
            required_bundle_class, compare, resume and existing.is_file(), parameter_overrides,
            case_id=case, suite=check, catalog=catalog, manifest=manifests[case],
        )
        print_run_summary(result, result_path)
        summary["results"].append({"suite": suite, "case": case, "status": result["status"],
                                   "summary": str(result_path)})
        summary["duration_seconds"] = time.monotonic() - started
        write_json_atomic(path, summary, "profile summary")
        if result["status"] == "failed":
            break
    _finish(summary)
    write_json_atomic(path, summary, "profile summary")
    return path, summary
