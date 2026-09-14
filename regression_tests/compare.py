"""Compare completed workflows and their recorded cold-stage references."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from .catalog import load_case_definition
from .bundles import ReferenceMatrix, load_reference_matrix
from .documents import load_json, write_json_atomic
from .support import BundleError, ComparisonError, utc_now
from .files import file_identity, recorded_directory, recorded_file, require_directory
from .compare_adaptive import compare_adaptive_files, mesh_differences
from .compare_fixed import compare_hdf5_files, check_output_contract
from .execute import final_execution
from .compare_common import (
    NewtonCheck, effective_newton_maximum, load_adaptive_tolerances,
    load_fixed_tolerances, read_newton_convergence, resolve_run_file, select_candidate,
)


@dataclass(frozen=True)
class ComparisonInputs:
    run_directory: Path
    case_directory: Path
    tolerances_path: Path
    plan: dict[str, Any]
    metadata: dict[str, Any]
    case: dict[str, Any]
    workflow: dict[str, Any]


@dataclass(frozen=True)
class ComparisonOverrides:
    candidate: Path | None = None
    reference: Path | None = None
    tolerance_profile: str | None = None
    newton_check: NewtonCheck = "bounded"
    direct_tolerance_profile: str | None = None


def load_comparison_inputs(
    run_directory: Path,
    case_directory: Path,
    tolerances_path: Path,
) -> ComparisonInputs:
    """Load and validate the documents needed to compare a completed run."""
    run_directory = require_directory(run_directory, "run")
    plan = load_json(run_directory / "run_plan.json", "run plan")
    metadata = _load_completed_metadata(run_directory)
    case = load_case_definition(plan["case_id"], case_directory)
    workflow = case["workflows"].get(plan.get("workflow_id"))
    if workflow is None:
        raise ComparisonError("run plan refers to an unknown workflow")
    return ComparisonInputs(
        run_directory=run_directory,
        case_directory=case_directory,
        tolerances_path=tolerances_path,
        plan=plan,
        metadata=metadata,
        case=case,
        workflow=workflow,
    )


def load_stage_inputs(
    parent: ComparisonInputs,
    run_directory: Path,
) -> ComparisonInputs:
    """Use parent workflow inputs with one completed stage's run data."""
    run_directory = require_directory(run_directory, "stage run")
    return replace(
        parent,
        run_directory=run_directory,
        metadata=_load_completed_metadata(run_directory),
    )


def _load_completed_metadata(run_directory: Path) -> dict[str, Any]:
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    if metadata.get("status") != "completed":
        raise ComparisonError("run metadata status is not completed")
    return metadata


def load_plan_matrix(
    plan: dict[str, Any],
    case: dict[str, Any],
    case_directory: Path,
) -> ReferenceMatrix | None:
    """Load the run's reference matrix and confirm its bundle identity."""
    bundle = plan.get("bundle")
    if not isinstance(bundle, dict) or not isinstance(bundle.get("root"), str):
        raise BundleError("run plan has no bundle identity")
    matrix = load_reference_matrix(Path(bundle["root"]), case["case_id"], case_directory)
    if matrix is not None and (
        bundle.get("bundle_id") != matrix.bundle_id
        or bundle.get("bundle_version") != matrix.bundle_version
    ):
        raise BundleError("run and golden-matrix bundle identities differ")
    return matrix


def _validated_stage_records(
    plan: dict[str, Any],
    workflow: dict[str, Any],
    metadata: dict[str, Any],
) -> list[dict[str, Any]]:
    if metadata.get("status") != "completed":
        raise ComparisonError("run metadata status is not completed")

    expected = [stage["stage_id"] for stage in workflow["stages"]]
    plan_stages = [stage.get("stage_id") for stage in plan.get("stages", [])]
    recorded = metadata.get("stages", [])
    recorded_stages = [stage.get("stage_id") for stage in recorded]
    if plan_stages != expected or recorded_stages != expected:
        raise ComparisonError("recorded stages differ from the tracked workflow")
    return recorded


def compare_completed_run(
    run_directory: Path, case_dir: Path, tolerances_path: Path,
    candidate_override: Path | None = None, reference_override: Path | None = None,
    tolerance_profile_override: str | None = None, comparison_policy_override: str | None = None,
    report_override: Path | None = None,
) -> tuple[str, Path, dict[str, Any]]:
    """Use stage references when available, unless a final-state check is requested."""
    inputs = load_comparison_inputs(run_directory, case_dir, tolerances_path)
    overrides = ComparisonOverrides(candidate_override, reference_override, tolerance_profile_override)
    explicit_final = any(value is not None for value in (
        candidate_override, reference_override, tolerance_profile_override, comparison_policy_override,
    ))
    if inputs.workflow.get("stages") and not explicit_final:
        matrix = load_plan_matrix(inputs.plan, inputs.case, inputs.case_directory)
        if matrix is not None:
            path, report = compare_reference_matrix(inputs, matrix, report_override)
            return "reference_matrix", path, report
    policy = comparison_policy_override or inputs.workflow.get("comparison_policy")
    path, report = compare_run(inputs, overrides, report_override, policy=policy)
    return report["comparison_policy"], path, report


def compare_run(
    inputs: ComparisonInputs, overrides: ComparisonOverrides,
    report_path: Path | None = None, *, policy: str | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Share file selection, Newton acceptance and report writing across methods."""
    policy = policy or inputs.workflow.get("comparison_policy")
    candidate = select_candidate(inputs.run_directory, inputs.metadata, overrides.candidate)
    reference = resolve_run_file(inputs.run_directory, overrides.reference, "inputs/reference.h5", "reference")
    profile_override = overrides.tolerance_profile
    selection = {"requested_policy": policy, "reason": "fixed mesh required"}
    if policy == "mesh_independent":
        coordinate_atol = load_json(inputs.tolerances_path, "tolerance definitions")["fixed_defaults"]["mesh_coordinate_atol"]
        differences = mesh_differences(reference, candidate, coordinate_atol)
        selection.update(mesh_coordinate_atol=coordinate_atol, differences=differences)
        if differences:
            selection["reason"] = "different discrete meshes"
        else:
            policy = "fixed_hdf5"
            selection["reason"] = "matching discrete meshes"
            profile_override = (overrides.direct_tolerance_profile or overrides.tolerance_profile
                                or inputs.workflow.get("direct_tolerance_profile"))
            if not profile_override:
                raise ComparisonError("adaptive workflow must declare a direct comparison profile")
    if policy == "fixed_hdf5":
        profile_id, tolerances = load_fixed_tolerances(
            inputs.tolerances_path, inputs.workflow, inputs.plan["layout_id"], profile_override,
        )
    elif policy == "mesh_independent":
        if inputs.workflow.get("comparison_policy") != "mesh_independent":
            raise ComparisonError("run does not define mesh-independent comparison")
        profile_id, tolerances = load_adaptive_tolerances(
            inputs.tolerances_path, inputs.workflow, overrides.tolerance_profile,
        )
    else:
        raise ComparisonError(f"unsupported comparison policy: {policy}")
    tolerances["newton_error_max"] = effective_newton_maximum(
        tolerances["newton_error_max"], overrides.newton_check,
    )
    protected = [candidate, reference]
    if policy == "fixed_hdf5":
        field_report = compare_hdf5_files(reference, candidate, tolerances)
        details = {
            "candidate": str(candidate), "reference": str(reference),
            "files": {name: {"path": str(path), **file_identity(path)}
                      for name, path in (("candidate", candidate), ("reference", reference))},
            "hdf5": field_report,
        }
    else:
        _, execution = final_execution(inputs.run_directory, inputs.metadata)
        runtime = execution.get("runtime_files", {}).get("positionFeketeNodesTri2D.h5", {})
        fekete = recorded_file(runtime.get("path"), "Fekete-node")
        protected.append(fekete)
        field_report = compare_adaptive_files(
            reference, candidate, fekete, tolerances["samples_per_element"], tolerances,
        )
        details = {key: value for key, value in field_report.items() if key != "tolerances"}
    convergence = read_newton_convergence(inputs.run_directory / "stdout.log", tolerances["newton_error_max"])
    failures = list(field_report["failures"])
    failures.extend(saved_output_contract(inputs)["failures"])
    if convergence.failure is not None:
        failures.append(convergence.failure)
    report = {
        **details,
        "comparison_policy": policy, "method_selection": selection,
        "schema_version": 2, "created_utc": utc_now(),
        "status": "passed" if not failures else "failed",
        "run_directory": str(inputs.run_directory),
        "case_id": inputs.plan["case_id"], "workflow_id": inputs.plan["workflow_id"],
        "layout_id": inputs.plan["layout_id"],
        "tolerance_profile": {"id": profile_id, **tolerances},
        "convergence": convergence.as_report(), "failures": failures,
    }
    output = (report_path or inputs.run_directory / "comparison.json").expanduser().resolve()
    if output in {path.expanduser().resolve() for path in protected}:
        raise ComparisonError("comparison report cannot replace a comparison input")
    write_json_atomic(output, report, "comparison report")
    return output, report


def compare_reference_matrix(
    context: ComparisonInputs, matrix: ReferenceMatrix, report_override: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Check the recorded stages in order, stopping at the first divergent stage."""
    stages = _validated_stage_records(context.plan, context.workflow, context.metadata)
    reports = []
    failures = []
    for stage, definition in zip(stages, context.workflow["stages"]):
        stage_id = stage["stage_id"]
        if stage.get("status") != "completed":
            raise ComparisonError(f"stage {stage_id} did not complete")
        directory = recorded_directory(stage.get("run_directory"), "stage run")
        candidate = recorded_file(stage.get("selected_hdf5"), "stage result")
        reference = matrix.reference_for(context.plan["workflow_id"], context.plan["layout_id"], stage_id)
        inputs = load_stage_inputs(context, directory)
        overrides = ComparisonOverrides(
            candidate, reference, context.workflow.get("stage_tolerance_profile"), definition["newton_check"],
            context.workflow.get("direct_stage_tolerance_profile"),
        )
        path, report = compare_run(
            inputs, overrides, directory / "comparison.json", policy=context.workflow["comparison_policy"],
        )
        reports.append({
            "stage_id": stage_id, "status": report["status"], "comparison_report": str(path),
            "reference": str(reference), "candidate": str(candidate), "failures": report["failures"],
            "comparison_policy": report["comparison_policy"],
        })
        if report["status"] != "passed":
            failures = [f"{stage_id}: {failure}" for failure in report["failures"] or ["comparison failed"]]
            break
    report = {
        "schema_version": 2, "created_utc": utc_now(),
        "status": "passed" if not failures and len(reports) == len(stages) else "failed",
        "run_directory": str(context.run_directory), "case_id": context.case["case_id"],
        "workflow_id": context.plan["workflow_id"], "layout_id": context.plan["layout_id"],
        "comparison_policy": "reference_matrix", "checked_stage_count": len(reports),
        "total_stage_count": len(stages), "first_failed_stage": reports[-1]["stage_id"] if failures else None,
        "stages": reports, "failures": failures,
    }
    output = (report_override or context.run_directory / "matrix_comparison.json").resolve()
    reserved = {Path(stage[name]).resolve() for stage in reports
                for name in ("candidate", "reference", "comparison_report")}
    if output in reserved:
        raise ComparisonError("matrix report cannot replace a stage artifact")
    write_json_atomic(output, report, "matrix report")
    return output, report


def compare_generated_meshes(
    reference_run: Path,
    candidate_run: Path,
) -> dict[str, Any] | None:
    """Require retained Gmsh adaptation outputs to be byte-identical."""
    reference = _generated_meshes(reference_run)
    candidate = _generated_meshes(candidate_run)
    if not reference and not candidate:
        return None

    failures = []
    files = {}
    for relative_path in sorted(reference.keys() | candidate.keys()):
        first = reference.get(relative_path)
        second = candidate.get(relative_path)
        first_identity = file_identity(first) if first else None
        second_identity = file_identity(second) if second else None
        passed = first_identity == second_identity
        files[relative_path] = {
            "passed": passed,
            "reference": first_identity,
            "candidate": second_identity,
        }
        if not passed:
            failures.append(f"generated mesh differs: {relative_path}")

    return {
        "mode": "byte_exact",
        "passed": not failures,
        "files": files,
        "failures": failures,
    }


def _generated_meshes(run_directory: Path) -> dict[str, Path]:
    return {
        str(path.relative_to(run_directory)): path
        for path in run_directory.glob("**/res/temp.msh")
    }



def producer_converged(
    source: dict[str, Any],
    workflow: dict[str, Any],
    policy: str,
    report: dict[str, Any],
    tolerances_path: Path,
) -> bool:
    """Check convergence independently of old-reference field agreement."""
    convergence = report.get("convergence")
    if isinstance(convergence, dict):
        return convergence.get("passed") is True
    if policy not in {"reference_matrix", "producer"}:
        return False

    try:
        return all(check.passed for _, check in _newton_checks(source, workflow, tolerances_path))
    except ComparisonError:
        return False


def _newton_checks(source, workflow, tolerances_path):
    run_directory = Path(source["run_directory"])
    metadata = load_json(run_directory / "run_metadata.json", "run metadata")
    definitions = workflow.get("stages", [{"newton_check": "bounded"}])
    records = metadata.get("stages") if workflow.get("stages") else [
        {"run_directory": str(run_directory), "status": metadata.get("status")}
    ]
    if not isinstance(records, list) or len(records) != len(definitions) or any(
        record.get("stage_id") != definition.get("stage_id")
        or record.get("status") != "completed"
        for record, definition in zip(records, definitions)
    ):
        raise ComparisonError("workflow has missing or incomplete stages")
    maximum = _stage_newton_maximum(workflow, source["layout_id"], tolerances_path)
    return [(record, read_newton_convergence(
        Path(record["run_directory"]) / "stdout.log",
        effective_newton_maximum(maximum, definition["newton_check"]),
    )) for record, definition in zip(records, definitions)]


def validate_completed_run(run_directory, case_directory, tolerances_path):
    """Validate outputs and stage convergence without agreement with a reference."""
    inputs = load_comparison_inputs(run_directory, case_directory, tolerances_path)
    if inputs.workflow.get("stages"):
        _validated_stage_records(inputs.plan, inputs.workflow, inputs.metadata)
    checks = _newton_checks({"run_directory": str(run_directory), "layout_id": inputs.plan["layout_id"]},
                            inputs.workflow, tolerances_path)
    stages, failures = [], []
    for record, convergence in checks:
        directory = Path(record["run_directory"])
        stage_failures = []
        output = None
        try:
            output = select_candidate(directory, _load_completed_metadata(directory))
            # Reuse the mesh, shape and finiteness contracts, not a second set
            # of field validators. Self-agreement is not regression evidence.
            mesh_differences(output, output, 0.)
            fields = compare_hdf5_files(output, output, {
                "mesh_coordinate_atol": 0., "relative_l2_max": 0., "normalized_linf_max": 0.,
            })
            stage_failures.extend(fields["failures"])
            stage_failures.extend(saved_output_contract(load_stage_inputs(inputs, directory))["failures"])
        except (ValueError, TypeError) as exc:
            stage_failures.append(str(exc))
        if convergence.failure:
            stage_failures.append(convergence.failure)
        label = record.get("stage_id", "final")
        failures.extend(f"{label}: {failure}" for failure in stage_failures)
        stages.append({"stage_id": label, "output": str(output) if output else None,
                       "convergence": convergence.as_report(), "failures": stage_failures,
                       "status": "failed" if stage_failures else "passed"})
    return {"status": "failed" if failures else "passed", "failures": failures, "stages": stages,
            "convergence": {"passed": all(check.passed for _, check in checks)}}


def _stage_newton_maximum(
    workflow: dict[str, Any],
    layout_id: str,
    tolerances_path: Path,
) -> float | None:
    profile = workflow.get("stage_tolerance_profile")
    if workflow.get("comparison_policy") == "fixed_hdf5":
        _, tolerances = load_fixed_tolerances(
            tolerances_path,
            workflow,
            layout_id,
            profile,
        )
    else:
        _, tolerances = load_adaptive_tolerances(
            tolerances_path,
            workflow,
            profile,
        )
    return tolerances["newton_error_max"]


def saved_output_contract(inputs):
    """Check the output contract when saved execution evidence is available."""
    from .config import verify_executable
    from .support import HarnessError

    directory, metadata = final_execution(inputs.run_directory, inputs.metadata)
    manifest = metadata.get("solver", {}).get("build_manifest")
    if not manifest:
        return {"status": "unavailable", "failures": []}  # Debug comparisons can lack execution records.
    try:
        plan = load_json(directory / "run_plan.json", "run plan")
        identity = verify_executable({"MHDG_BUILD_MANIFEST": manifest["path"]}, inputs.workflow["model"],
                                     plan["layout"]["execution"], Path(metadata["executable"]["path"]))
        return check_output_contract(select_candidate(directory, metadata), inputs.workflow["model"],
                                     identity, plan.get("parameter_overrides", {}))
    except (HarnessError, OSError) as exc:
        return {"status": "failed", "failures": [str(exc)]}
