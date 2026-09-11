"""Run and verify scientific suites; resume reuses completed, unchanged outputs."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import shutil
import time
from typing import Any

from bundle.cases import load_case_definition
from bundle.settings import bundle_root_from_settings, read_settings
from bundle.validation import validate_bundle_root
from .config import load_suite_definition, require_bundle_class
from .prepare import prepare_run
from .execute import execute_prepared, reusable_outputs
from .compare import compare_completed_run, compare_generated_meshes, producer_converged
from .compare_common import select_candidate
from support.documents import load_json, write_json_atomic
from support.errors import BundleError, ComparisonError, HarnessError
from support.files import file_identity
from support.paths import require_file, recorded_directory
from support.identifiers import IDENTIFIER_RE
from support.time import utc_now, utc_run_id


@dataclass(frozen=True)
class SuiteRunInputs:
    settings: dict[str, str]
    case_id: str
    run_id: str
    case_directory: Path
    layouts_path: Path
    tolerances_path: Path
    compare: bool
    parameter_overrides: dict[str, bool | float | int | str]


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
        )
        result["run_directory"] = str(prepared.path)
        run = execute_prepared(prepared, inputs.settings)
        result["run_status"] = run.status
        result["duration_seconds"] = run.duration_seconds
        result["status"] = run.status
        if run.status != "completed":
            result["failures"] = [f"solver run status is {run.status}"]
            return result
        return _compare_cell(result, inputs)
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
        "comparison_status": None,
        "duration_seconds": None,
        "run_directory": None,
        "comparison_report": None,
        "failures": [],
    }


def compare_layout_pairs(
    summary: dict[str, Any],
    case_directory: Path,
    tolerances_path: Path,
) -> list[dict[str, Any]]:
    """Compare every declared candidate with its same-workflow baseline."""
    return [
        _compare_pair(
            summary,
            workflow_id,
            pair,
            case_directory,
            tolerances_path,
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
            from regression_tests.diagnostics import compare_outputs

            diagnostic_report = compare_outputs(reference, _selected_output(candidate))
            result["diagnostics"] = diagnostic_report
            if diagnostic_report["status"] != "passed":
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


def _verify_result(
    source: dict[str, Any],
    case: dict[str, Any],
    case_directory: Path,
    tolerances_path: Path,
) -> dict[str, Any]:
    workflow_id = source.get("workflow_id")
    result = {
        "workflow_id": workflow_id,
        "layout_id": source.get("layout_id"),
        "run_directory": source.get("run_directory"),
        "run_status": source.get("run_status"),
        "comparison_policy": None,
        "comparison_report": None,
        "convergence_status": None,
        "status": "failed",
        "failures": [],
    }
    if source.get("run_status") != "completed":
        result["failures"] = ["solver run did not complete"]
        return result
    if workflow_id not in case["workflows"]:
        result["failures"] = [f"unknown workflow: {workflow_id}"]
        return result

    try:
        run_directory = Path(source["run_directory"])
        policy, report_path, report = compare_completed_run(
            run_directory,
            case_directory,
            tolerances_path,
        )
        result["comparison_policy"] = policy
        result["comparison_report"] = str(report_path)
        result["convergence_status"] = (
            "passed"
            if producer_converged(
                source,
                case["workflows"][workflow_id],
                policy,
                report,
                tolerances_path,
            )
            else "failed"
        )
        result["failures"] = report["failures"]
        result["status"] = report["status"]
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
    settings: dict[str, str],
    bundle_root: Path,
    layout_ids: list[str],
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
    if any(layout.startswith("serial_") for layout in layout_ids):
        records["serial_executable"] = _setting_file_record(
            settings,
            "MHDG_SERIAL_EXECUTABLE",
        )
    if any(layout.startswith("mpi") for layout in layout_ids):
        records["parallel_executable"] = _setting_file_record(
            settings,
            "MHDG_PARALLEL_EXECUTABLE",
        )
        launcher = settings.get("MHDG_MPI_LAUNCHER")
        resolved_launcher = shutil.which(launcher) if launcher else None
        if resolved_launcher is None:
            raise BundleError("MHDG_MPI_LAUNCHER is not executable or not found")
        records["mpi_launcher"] = _file_record(
            Path(resolved_launcher),
            "MPI launcher",
        )
    # Runtime inputs may also come from a prebuilt executable without a manifest.
    # Record only selected paths, once even when both variants share the file.
    from .prepare import _runtime_files

    records["runtime_files"] = {
        str(path): file_identity(path)
        for name in ("serial_executable", "parallel_executable") if name in records
        for path in _runtime_files(Path(records[name]["path"])).values()
    }
    for key, name in (
        ("MHDG_ENVIRONMENT_SCRIPT", "environment_script"),
        ("MHDG_BUILD_MANIFEST", "build_manifest"),
    ):
        if settings.get(key):
            records[name] = _setting_file_record(settings, key)
    return records


def _setting_file_record(
    settings: dict[str, str],
    key: str,
) -> dict[str, Any]:
    value = settings.get(key)
    if not value:
        raise BundleError(f"settings must define {key}")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise BundleError(f"{key} must be an absolute path")
    return _file_record(path, key)


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
    settings_path, suite_id, case_directory, layouts_path, suites_path,
    tolerances_path, run_id=None, required_bundle_class=None, compare=True,
    resume=False, parameter_overrides=None, *, case_id=None,
):
    suite = load_suite_definition(suite_id, suites_path, layouts_path, case_directory, case_id=case_id)
    settings = read_settings(settings_path)
    bundle_root = bundle_root_from_settings(settings)
    validate_bundle_root(bundle_root, case_directory)
    if required_bundle_class:
        require_bundle_class(bundle_root, case_directory, required_bundle_class)
    run_id = run_id or utc_run_id()
    if not IDENTIFIER_RE.fullmatch(run_id):
        raise BundleError(f"invalid suite run identifier: {run_id}")
    pairs = suite.get("layout_comparisons")
    references = suite["reference_comparisons"]
    if not pairs and not references:
        mode = "execution_only"
    elif not compare:
        mode = "deferred" if references else "deferred_layout_pairs"
    elif pairs:
        mode = "reference_and_layout_pairs" if references else "layout_pairs"
    else:
        mode = "immediate"
    overrides = {"balance_diagnostics_mode": suite["diagnostics"], **(parameter_overrides or {})}
    identity = suite_execution_inputs(settings, bundle_root, suite["layouts"], {
        "case_definition": case_directory / f"{suite['case_id']}.json",
        "workflow_catalog": case_directory.parent / "workflows.json",
        "layout_catalog": layouts_path, "suite_catalog": suites_path,
        "tolerance_catalog": tolerances_path,
    })
    configured = settings.get("MHDG_REGRESSION_RUN_ROOT", "")
    root = Path(configured).expanduser()
    if not root.is_absolute():
        raise BundleError("MHDG_REGRESSION_RUN_ROOT must be an absolute path")
    directory = root.resolve() / "suites" / suite_id / suite["case_id"] / run_id
    path = directory / "suite_summary.json"
    expected = {
        "suite_id": suite_id, "run_id": run_id, "case_id": suite["case_id"],
        "workflow_ids": suite["workflow_ids"], "layout_ids": suite["layouts"],
        "comparison_mode": mode, "parameter_overrides": overrides,
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
                            layouts_path, tolerances_path, compare and references, overrides)
    print(f"suite: {suite_id} ({run_id})")
    for layout in suite["layouts"]:
        for workflow in suite["workflow_ids"]:
            started = time.monotonic()
            previous = next((item for item in summary["results"]
                             if (item["workflow_id"], item["layout_id"]) == (workflow, layout)), None)
            if (previous and previous.get("run_status") == "completed"
                    and previous.get("run_directory") and reusable_outputs(Path(previous["run_directory"]))):
                print(f"reusing completed {workflow} / {layout}", flush=True)
                result = _compare_cell(dict(previous), inputs)
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
        summary["comparisons"] = compare_layout_pairs(summary, case_directory, tolerances_path)
    if compare:
        from .diagnostics import check_suite

        summary["diagnostics"] = check_suite(summary, required=False)
    _finish(summary)
    write_json_atomic(path, summary, "suite summary")
    return path, summary


def _compare_cell(result, inputs):
    if result.get("run_status") != "completed":
        return result
    result.update(status="passed", failures=[], comparison_status="not_run")
    if inputs.compare:
        try:
            policy, path, report = compare_completed_run(
                Path(result["run_directory"]), inputs.case_directory, inputs.tolerances_path,
            )
            result.update(comparison_policy=policy, comparison_report=str(path),
                          comparison_status=report["status"], failures=report["failures"],
                          status="passed" if report["status"] == "passed" else "comparison_failed")
        except HarnessError as exc:
            result.update(status="error", failures=[str(exc)])
    return result


def _finish(summary):
    checks = [*summary["results"], *summary.get("comparisons", [])]
    if summary.get("diagnostics"):
        checks.append(summary["diagnostics"])
    summary.update(status="passed" if checks and all(item["status"] == "passed" for item in checks) else "failed",
                   finished_utc=utc_now())


def verify_suite(suite_summary_path, case_directory, tolerances_path, *, include_layout_pairs=True):
    path = require_file(suite_summary_path, "suite summary")
    source = load_json(path, "suite summary")
    _validate_source_summary(source)
    case = load_case_definition(source["case_id"], case_directory)
    output = path.parent / "verification_summary.json"
    summary = {
        "schema_version": 2, "created_utc": utc_now(), "status": "running",
        "source_summary": str(path), "suite_id": source["suite_id"],
        "run_id": source["run_id"], "case_id": source["case_id"], "results": [],
    }
    write_json_atomic(output, summary, "verification summary")
    if source.get("reference_comparisons", not source.get("layout_comparisons")):
        for result in source["results"]:
            summary["results"].append(_verify_result(result, case, case_directory, tolerances_path))
            write_json_atomic(output, summary, "verification summary")
    if include_layout_pairs and source.get("layout_comparisons"):
        summary["comparisons"] = compare_layout_pairs(source, case_directory, tolerances_path)
    from .diagnostics import check_suite

    summary["diagnostics"] = check_suite(source, required=False)
    _finish(summary)
    write_json_atomic(output, summary, "verification summary")
    return output, summary


def run_profile(name, checks, settings_by_case, catalog_root, run_id=None, *, resume=False,
                compare=True, required_bundle_class=None, parameter_overrides=None):
    """Run an ordered selection using ordinary suite summaries for resume."""
    from .reporting import print_run_summary

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
    # Reject unavailable/wrong-case bundles before running an earlier, costly check.
    for case_id, values in settings_by_case.items():
        bundle = bundle_root_from_settings(values)
        if validate_bundle_root(bundle, catalog_root / "cases").case_id != case_id:
            raise BundleError(f"selected bundle does not contain {case_id}")
        if required_bundle_class:
            require_bundle_class(bundle, catalog_root / "cases", required_bundle_class)
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
            case_id=case,
        )
        print_run_summary(result, result_path)
        summary["results"].append({"suite": suite, "case": case, "status": result["status"],
                                   "summary": str(result_path)})
        summary["duration_seconds"] = time.monotonic() - started
        write_json_atomic(path, summary, "profile summary")
        if result["status"] != "passed":
            break
    summary["status"] = "passed" if (len(summary["results"]) == len(checks)
        and all(result["status"] == "passed" for result in summary["results"])) else "failed"
    write_json_atomic(path, summary, "profile summary")
    return path, summary
