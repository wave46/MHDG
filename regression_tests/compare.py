"""Validate workflow stages and compare final states and selected run pairs."""

from __future__ import annotations

from dataclasses import dataclass, replace
import math
from pathlib import Path
import re
from typing import Any

import h5py

from .catalog import load_case_definition
from .documents import load_json, write_json_atomic
from .support import ComparisonError, utc_now
from .files import file_identity, recorded_file, require_directory
from .compare_adaptive import compare_adaptive_files, mesh_differences
from .compare_fixed import compare_hdf5_files, check_output_contract, validate_solution_file, required_scalar
from .execute import final_execution, outputs_unchanged
from .compare_common import (
    select_logged_output, NewtonCheck, effective_newton_maximum, load_adaptive_tolerances,
    load_fixed_tolerances, load_tolerance_catalog, read_newton_convergence, resolve_run_file, select_candidate,
    read_me_completion_errors, NewtonConvergence,
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
    catalog: dict


@dataclass(frozen=True)
class ComparisonOverrides:
    candidate: Path | None = None
    reference: Path | None = None
    tolerance_profile: str | None = None
    newton_check: NewtonCheck = "bounded"


def load_comparison_inputs(
    run_directory: Path,
    case_directory: Path,
    tolerances_path: Path,
    *, catalog=None, require_completed=True,
) -> ComparisonInputs:
    """Load and validate the documents needed to compare a completed run."""
    catalog = {} if catalog is None else catalog
    run_directory = require_directory(run_directory, "run")
    plan = load_json(run_directory / "run_plan.json", "run plan")
    metadata = (_load_completed_metadata(run_directory) if require_completed else
                load_json(run_directory / "run_metadata.json", "run metadata"))
    case = load_case_definition(plan["case_id"], case_directory, catalog=catalog)
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
        workflow=workflow, catalog=catalog,
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


def _validated_stage_records(
    plan: dict[str, Any],
    workflow: dict[str, Any],
    metadata: dict[str, Any],
) -> list[dict[str, Any]]:
    expected = [stage["id"] for stage in workflow["stages"]]
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
    *, catalog=None,
) -> tuple[str, Path, dict[str, Any]]:
    """Validate cold stages and compare the final state with one explicit reference."""
    inputs = load_comparison_inputs(run_directory, case_dir, tolerances_path, catalog=catalog)
    final_stage = inputs.workflow.get("stages", [{}])[-1]
    overrides = ComparisonOverrides(candidate_override, reference_override, tolerance_profile_override,
                                    final_stage.get("newton_check", "bounded"))
    policy = comparison_policy_override or inputs.workflow.get("comparison", {}).get("method")
    path, report = compare_run(inputs, overrides, report_override, policy=policy)
    return report["comparison_policy"], path, report


def compare_run(
    inputs: ComparisonInputs, overrides: ComparisonOverrides,
    report_path: Path | None = None, *, policy: str | None = None, outcome=None,
) -> tuple[Path, dict[str, Any]]:
    """Share file selection, Newton acceptance and report writing across methods."""
    comparison = inputs.workflow.get("comparison", {})
    policy = policy or comparison.get("method")
    candidate = select_candidate(inputs.run_directory, inputs.metadata, overrides.candidate)
    reference = resolve_run_file(inputs.run_directory, overrides.reference, "inputs/reference.h5", "reference")
    profile_override = overrides.tolerance_profile
    selection = {"requested_policy": policy, "reason": "fixed mesh required"}
    if policy == "mesh_independent":
        coordinate_atol = load_tolerance_catalog(inputs.tolerances_path, inputs.catalog)["fixed_defaults"]["mesh_coordinate_atol"]
        differences = mesh_differences(reference, candidate, coordinate_atol)
        selection.update(mesh_coordinate_atol=coordinate_atol, differences=differences)
        if differences:
            selection["reason"] = "different discrete meshes"
        else:
            policy = "fixed_hdf5"
            selection["reason"] = "matching discrete meshes"
            profile_override = (overrides.tolerance_profile
                                or comparison.get("direct_profile"))
            if not profile_override:
                raise ComparisonError("adaptive workflow must declare a direct comparison profile")
    if policy == "fixed_hdf5":
        profile_id, tolerances = load_fixed_tolerances(
            inputs.tolerances_path, inputs.workflow, inputs.plan["layout_id"], profile_override, catalog=inputs.catalog,
        )
    elif policy == "mesh_independent":
        if comparison.get("method") != "mesh_independent":
            raise ComparisonError("run does not define mesh-independent comparison")
        profile_id, tolerances = load_adaptive_tolerances(
            inputs.tolerances_path, inputs.workflow, overrides.tolerance_profile, catalog=inputs.catalog,
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
    if outcome is None:
        convergence = read_newton_convergence(inputs.run_directory / "stdout.log", tolerances["newton_error_max"])
        preceding = _validate_outputs(inputs, _newton_checks(inputs, before_final=True))
        moving = _validate_me_steps(inputs, tolerances["newton_error_max"])
    else:
        moving = outcome["completed_prefix"]
        # A matching-mesh fallback can select a stricter Newton profile.
        for step in moving["steps"][1:]:
            if step["convergence"]["maximum"] != tolerances["newton_error_max"]:
                check = NewtonConvergence(step["convergence"]["final_newton_error"], tolerances["newton_error_max"])
                step["convergence"] = check.as_report()
                if check.failure:
                    step["failures"].append(check.failure)
                    step["status"] = "failed"
                    moving["failures"].append(f"step {step['step']}: {check.failure}")
        moving["convergence"]["passed"] = all(step["convergence"]["passed"] for step in moving["steps"][1:])
        moving["status"] = "failed" if moving["failures"] else "passed"
        convergence = NewtonConvergence(moving["steps"][-1]["convergence"]["final_newton_error"],
                                        tolerances["newton_error_max"])
        preceding = _validate_outputs(inputs, [])
    refinement = _check_refinement(inputs, candidate)
    protected.extend(Path(stage["output"]) for stage in preceding["stages"] if stage["output"])
    if refinement.get("initial_output"):
        protected.append(Path(refinement["initial_output"]))
    failures = list(field_report["failures"])
    failures.extend(saved_output_contract(inputs, candidate)["failures"])
    failures.extend(preceding["failures"])
    failures.extend(refinement.get("failures", []))
    failures.extend(moving.get("failures", []))
    protected.extend(Path(step["output"]) for step in moving.get("steps", []))
    if moving:
        with h5py.File(candidate, "r") as first, h5py.File(reference, "r") as second:
            candidate_time, reference_time = _me_time(first), _me_time(second)
        if (candidate_time[0] != reference_time[0] or any(not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)
                for a, b in zip(candidate_time[1:3], reference_time[1:3]))):
            failures.append("ME reference step/time/dt differs from the completed candidate")
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
        "convergence": {**convergence.as_report(),
                        "passed": convergence.passed and preceding["convergence"]["passed"]
                                  and moving.get("convergence", {}).get("passed", True)},
        **({"moving_equilibrium": moving} if moving else {}),
        **({"expected_failure": {key: outcome[key] for key in
             ("run_status", "terminal_reason", "last_attempted_step")}} if outcome else {}),
        "stages": preceding["stages"], "refinement": refinement, "failures": failures,
    }
    output = (report_path or inputs.run_directory / "comparison.json").expanduser().resolve()
    if output in {path.expanduser().resolve() for path in protected}:
        raise ComparisonError("comparison report cannot replace a comparison input")
    write_json_atomic(output, report, "comparison report")
    return output, report


def compare_me_refinement_stop(directory, case_directory, tolerances_path, *, completed_steps,
                               reference_override=None, report_override=None, catalog=None):
    """Compare a validated completed prefix; raw failure metadata remains intact."""
    inputs = load_comparison_inputs(directory, case_directory, tolerances_path,
                                    catalog=catalog, require_completed=False)
    if not inputs.workflow.get("moving_equilibrium"):
        raise ComparisonError("expected refinement-stop comparison requires a moving-equilibrium workflow")
    maximum = _stage_newton_maximum(inputs.workflow, inputs.plan["layout_id"],
                                    tolerances_path, catalog=inputs.catalog)
    outcome = validate_me_refinement_stop(inputs.run_directory, inputs.plan, inputs.metadata,
                                          maximum, completed_steps=completed_steps)
    if outcome["status"] != "passed":
        raise ComparisonError("invalid expected ME refinement stop: " + "; ".join(outcome["failures"]))
    overrides = ComparisonOverrides(candidate=Path(outcome["selected_output"]), reference=reference_override)
    path, report = compare_run(inputs, overrides, report_override, outcome=outcome)
    return report["comparison_policy"], path, report


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



def _newton_checks(inputs, *, before_final=False):
    workflow = inputs.workflow
    definitions = workflow.get("stages", [{"newton_check": "bounded"}])
    records = _validated_stage_records(inputs.plan, workflow, inputs.metadata) if workflow.get("stages") else [
        {"run_directory": str(inputs.run_directory), "status": "completed"}
    ]
    if any(record.get("status") != "completed" for record in records):
        raise ComparisonError("workflow has incomplete stages")
    if before_final:
        records, definitions = records[:-1], definitions[:-1]
    if not records:
        return []
    maximum = _stage_newton_maximum(workflow, inputs.plan["layout_id"], inputs.tolerances_path, catalog=inputs.catalog)
    return [(record, read_newton_convergence(
        Path(record["run_directory"]) / "stdout.log",
        effective_newton_maximum(maximum, definition["newton_check"]),
    )) for record, definition in zip(records, definitions)]


def validate_completed_run(run_directory, case_directory, tolerances_path, *, catalog=None):
    """Validate outputs and stage convergence without agreement with a reference."""
    inputs = load_comparison_inputs(run_directory, case_directory, tolerances_path, catalog=catalog)
    report = _validate_outputs(inputs, _newton_checks(inputs))
    output = report["stages"][-1]["output"]
    refinement = _check_refinement(inputs, Path(output)) if output else {}
    report["refinement"] = refinement
    report["failures"].extend(refinement.get("failures", []))
    if inputs.workflow.get("moving_equilibrium"):
        moving = _validate_me_steps(inputs, _stage_newton_maximum(
            inputs.workflow, inputs.plan["layout_id"], inputs.tolerances_path, catalog=inputs.catalog))
        report["moving_equilibrium"] = moving
        report["failures"].extend(moving["failures"])
        report["convergence"]["passed"] &= moving["convergence"]["passed"]
    report["status"] = "failed" if report["failures"] else "passed"
    return report


def _validate_me_steps(inputs, maximum):
    if not inputs.workflow.get("moving_equilibrium"):
        return {}
    report = validate_me_completed_steps(inputs.run_directory, inputs.plan, inputs.metadata, maximum)
    if not report["complete"]:
        report["failures"].append("moving equilibrium did not complete all prescribed advances")
    report["status"] = "failed" if report["failures"] else "passed"
    return report


def validate_me_completed_steps(directory, plan, metadata, maximum):
    """Validate a contiguous completed prefix; complete separately marks full coverage."""
    dt = plan["moving_equilibrium_inputs"]["time_step_s"]
    prescribed = plan["parameter_overrides"]["nts"]
    saves = read_me_completion_errors(directory / "stdout.log")
    steps, failures, newton_passed = {}, [], True
    for value in metadata.get("hdf5_outputs", []):
        path = recorded_file(value, "HDF5 output", directory)
        # setSolName distinguishes NR saves; numeric suffixes do not identify time steps.
        if re.search(r"_NR\d+(?:_UNCP)?\.h5$", path.name):
            continue
        problems = validate_solution_file(path)
        with h5py.File(path, "r") as handle:
            number, time, saved_dt, scale = _me_time(handle)
        if not (math.isfinite(number) and number == int(number) and 0 <= number <= prescribed):
            raise ComparisonError(f"invalid completed ME step: {path}")
        number = int(number)
        if not (math.isfinite(scale) and scale > 0 and math.isclose(time, number * dt, rel_tol=1e-10, abs_tol=1e-12)
                and math.isclose(saved_dt, dt, rel_tol=1e-10, abs_tol=1e-12)):
            problems.append("completed step/time differs from the prescribed ME sequence")
        logged_step, error = saves.get(path, (None, None))
        if logged_step != number:
            problems.append("completed output is not logged for its embedded timestep")
        convergence = NewtonConvergence(error, maximum)
        if number:
            newton_passed &= convergence.passed
        if number and convergence.failure:
            problems.append(convergence.failure)
        if number in steps:
            duplicate = compare_hdf5_files(Path(steps[number]["output"]), path,
                {"mesh_coordinate_atol": 0., "relative_l2_max": 0., "normalized_linf_max": 0.})
            if duplicate["status"] != "passed":
                problems.append("conflicting completed copies of the same timestep")
        else:
            steps[number] = {"step": number, "time_s": time if math.isfinite(time) else None, "output": str(path),
                             "convergence": convergence.as_report() if number else None,
                             "failures": [], "status": "passed"}
        steps[number]["failures"].extend(problems)
        if problems:
            steps[number]["status"] = "failed"
        failures.extend(f"step {number}: {problem}" for problem in problems)
    completed = sorted(number for number in steps if number)
    if 0 not in steps:
        failures.append("missing initial ME save")
    if completed != list(range(1, len(completed) + 1)):
        failures.append("completed ME advances are not a contiguous prefix")
    return {"status": "failed" if failures else "passed", "failures": failures,
            "completed_steps": completed, "prescribed_steps": prescribed,
            "complete": completed == list(range(1, prescribed + 1)),
            "convergence": {"passed": newton_passed},
            "steps": [steps[number] for number in sorted(steps)]}


def _me_time(handle):
    number = required_scalar(handle, "simulation_parameters/time", "Current_time_step_number")
    scale = required_scalar(handle, "simulation_parameters/adimensionalization", "time_scale")
    time = required_scalar(handle, "simulation_parameters/time", "Current_time") * scale
    dt = required_scalar(handle, "simulation_parameters/time", "Current_time_step") * scale
    return number, time, dt, scale


def validate_me_refinement_stop(directory, plan, metadata, maximum, *, completed_steps):
    """Validate an expected divergence-budget stop without rewriting raw status."""
    report = {"status": "failed", "run_status": metadata.get("status"), "failures": [],
              "terminal_reason": None, "completed_prefix": None, "selected_output": None}
    failures = report["failures"]
    if not 1 <= completed_steps < plan["parameter_overrides"]["nts"]:
        raise ComparisonError("expected failure requires a nonempty incomplete ME prefix")
    if ((metadata.get("exit_code"), metadata.get("status")) not in
            {(0, "solver_reported_error"), (1, "solver_failed")} or metadata.get("launch_error")):
        failures.append("execution is not a supported numerical stop")
        return report
    if not outputs_unchanged(directory, metadata):
        failures.append("expected-failure outputs are missing or changed")
        return report
    prefix = validate_me_completed_steps(directory, plan, metadata, maximum)
    report["completed_prefix"] = prefix
    failures.extend(prefix["failures"])
    if prefix["completed_steps"] != list(range(1, completed_steps + 1)) or prefix["complete"]:
        failures.append("completed prefix differs from the expected failure")
    log = (directory / "stdout.log").read_text(encoding="utf-8", errors="replace")
    iterations = list(re.finditer(r"Time iteration\s*=\s*(\d+)", log))
    entered = [int(match[1]) for match in iterations]
    report["last_attempted_step"] = entered[-1] if entered else None
    if entered != list(range(1, completed_steps + 2)):
        failures.append("expected failed timestep was not the last attempted advance")
    limit = plan["parameter_overrides"]["max_divergence_refinements"]
    stop = re.compile(rf"NR refinement recovery stopped: per-timestep refinement limit reached:\s*{re.escape(str(limit))}\s*$",
                      re.MULTILINE)
    final_attempt = log[iterations[-1].end():] if iterations else ""
    if not stop.search(final_attempt):
        failures.append("missing expected divergence-refinement-limit terminal reason")
    elif any(not stop.search(message) for message in metadata.get("fatal_log_messages", [])):
        failures.append("execution also reported an unrelated fatal error")
    else:
        report["terminal_reason"] = "divergence_refinement_limit"
    if not failures:
        report.update(status="passed", selected_output=prefix["steps"][-1]["output"])
    return report


def _validate_outputs(inputs, checks):
    """Validate recorded stages once; a final field comparison owns its own output."""
    stages, failures = [], []
    for record, convergence in checks:
        directory = Path(record["run_directory"])
        stage_failures = []
        output = None
        try:
            output = select_candidate(directory, _load_completed_metadata(directory))
            stage_failures.extend(validate_solution_file(output))
            stage_failures.extend(saved_output_contract(load_stage_inputs(inputs, directory), output)["failures"])
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


def _check_refinement(inputs, candidate):
    """Require a valid mesh change and a retained adaptation output."""
    if not inputs.workflow.get("require_refinement"):
        return {}
    directory, metadata = final_execution(inputs.run_directory, inputs.metadata)
    report = {"status": "failed", "failures": []}
    try:
        recorded = [recorded_file(path, "HDF5 output", directory) for path in metadata["hdf5_outputs"]]
        initial = select_logged_output(directory, recorded, first=True)
        if initial is None:
            raise ComparisonError("missing initial output for adaptation check")
        if initial == candidate.resolve():
            raise ComparisonError("adaptation needs a separate recorded initial output")
        with h5py.File(initial, "r") as first, h5py.File(candidate, "r") as last:
            before = int(required_scalar(first, "mesh", "Nelems"))
            after = int(required_scalar(last, "mesh", "Nelems"))
        report.update(initial_output=str(initial), initial_elements=before, final_elements=after)
        coordinate_atol = load_tolerance_catalog(inputs.tolerances_path, inputs.catalog)["fixed_defaults"]["mesh_coordinate_atol"]
        differences = mesh_differences(initial, candidate, coordinate_atol)
        report.update(mesh_coordinate_atol=coordinate_atol, mesh_differences=differences)
        if not differences:
            raise ComparisonError("adaptive probe did not change the mesh")
        generated = recorded_file(str(directory / "res/temp.msh"), "generated adaptation mesh")
        report["generated_mesh"] = {"path": str(generated), **file_identity(generated)}
        report["status"] = "passed"
    except (OSError, KeyError, ValueError, TypeError) as exc:
        report["failures"].append(str(exc))
    return report


def _stage_newton_maximum(
    workflow: dict[str, Any],
    layout_id: str,
    tolerances_path: Path,
    *, catalog=None,
) -> float | None:
    comparison = workflow.get("comparison", {})
    if comparison.get("method") == "fixed_hdf5":
        _, tolerances = load_fixed_tolerances(
            tolerances_path,
            workflow,
            layout_id,
            None, catalog=catalog,
        )
    else:
        _, tolerances = load_adaptive_tolerances(
            tolerances_path,
            workflow,
            None, catalog=catalog,
        )
    return tolerances["newton_error_max"]


def saved_output_contract(inputs, candidate):
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
        return check_output_contract(candidate, directory, inputs.workflow["model"],
                                     identity, plan.get("parameter_overrides", {}))
    except (HarnessError, OSError) as exc:
        return {"status": "failed", "failures": [str(exc)]}
