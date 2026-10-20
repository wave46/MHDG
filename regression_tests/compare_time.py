"""Accepted checkpoint selection, physical clocks and temporal regression checks."""

from __future__ import annotations

import math
from pathlib import Path
import re

import h5py

from .compare_common import read_me_completion_errors, NewtonConvergence
from .compare_fixed import compare_hdf5_files, validate_solution_file, required_scalar, required_array
from .files import recorded_file
from .support import ComparisonError


def accepted_checkpoint_files(directory, metadata):
    """Yield ordinary solution files, excluding NR diagnostic saves."""
    for value in metadata.get("hdf5_outputs", []):
        path = recorded_file(value, "HDF5 output", directory)
        if re.search(r"_NR\d+(?:_UNCP)?\.h5$", path.name):
            continue
        yield path


def accepted_checkpoints(directory, metadata):
    """Pair ordinary solution files with their validated accepted-step clocks."""
    for path in accepted_checkpoint_files(directory, metadata):
        with h5py.File(path) as handle:
            clock = read_checkpoint_time(handle)
        number = clock[0]
        if not math.isfinite(number) or number != int(number) or number < 0:
            raise ComparisonError(f"invalid accepted checkpoint step: {path}")
        yield path, clock


def validate_me_completed_steps(directory, plan, metadata, maximum, *, topology=None, acceptance_checks=None):
    """Validate a contiguous completed prefix; complete separately marks full coverage."""
    preparation = plan["moving_equilibrium_inputs"]
    physical = bool(preparation.get("time_mode"))
    dt = preparation["time_step_s"]
    prescribed = None if physical else plan["parameter_overrides"]["nts"]
    initial = preparation.get("initial_step", 0)
    saves = read_me_completion_errors(directory / "stdout.log")
    steps, failures, newton_passed = {}, [], True
    for path, (number, time, saved_dt, scale) in accepted_checkpoints(directory, metadata):
        problems = validate_solution_file(path)
        proposal = None
        with h5py.File(path, "r") as handle:
            proposal_path = "simulation_parameters/time/Proposed_time_step_seconds"
            if proposal_path in handle:
                proposal = required_scalar(handle, "simulation_parameters/time", "Proposed_time_step_seconds")
            if physical:
                seconds = required_scalar(handle, "simulation_parameters/time", "Current_time_seconds")
                if not math.isclose(time, seconds, rel_tol=1e-10, abs_tol=1e-12):
                    problems.append("physical timestamp differs from normalized checkpoint time")
        if not (math.isfinite(number) and number == int(number) and number >= initial
                and (physical or number <= prescribed)):
            raise ComparisonError(f"invalid completed ME step: {path}")
        number = int(number)
        if physical:
            if not (math.isfinite(scale) and scale > 0 and math.isfinite(time) and math.isfinite(saved_dt)):
                problems.append("nonfinite physical clock or invalid time scale")
        elif not (math.isfinite(scale) and scale > 0 and math.isclose(time, number * dt, rel_tol=1e-10, abs_tol=1e-12)
                  and math.isclose(saved_dt, dt, rel_tol=1e-10, abs_tol=1e-12)):
            problems.append("completed step/time differs from the prescribed ME sequence")
        logged_step, error = saves.get(path, (None, None))
        if logged_step != number and not (number == initial and logged_step == 0):
            problems.append("completed output is not logged for its embedded timestep")
        convergence = NewtonConvergence(error, maximum)
        if number > initial:
            newton_passed &= convergence.passed
        if number > initial and convergence.failure:
            problems.append(convergence.failure)
        if number in steps:
            duplicate = compare_hdf5_files(Path(steps[number]["output"]), path,
                {"mesh_coordinate_atol": 0., "relative_l2_max": 0., "normalized_linf_max": 0.})
            if duplicate["status"] != "passed":
                problems.append("conflicting completed copies of the same timestep")
        else:
            steps[number] = {"step": number, "time_s": time if math.isfinite(time) else None,
                             "dt_s": saved_dt if math.isfinite(saved_dt) else None,
                             "proposed_dt_s": proposal, "output": str(path),
                             "convergence": convergence.as_report() if number > initial else None,
                             "failures": [], "status": "passed"}
        steps[number]["failures"].extend(problems)
        if problems:
            steps[number]["status"] = "failed"
        failures.extend(f"step {number}: {problem}" for problem in problems)
    completed = sorted(number for number in steps if number > initial)
    if initial not in steps:
        failures.append("missing initial ME save")
    if completed != list(range(initial + 1, initial + len(completed) + 1)):
        failures.append("completed ME advances are not a contiguous prefix")
    failures.extend(_check_me_topology(steps, completed, topology))
    if physical:
        clock = _check_physical_sequence(steps, completed, preparation, directory, acceptance_checks or {})
        failures.extend(clock["failures"])
        complete = clock["complete"]
    else:
        clock = None
        complete = completed == list(range(initial + 1, prescribed + 1))
    return {"status": "failed" if failures else "passed", "failures": failures,
            "completed_steps": completed, "prescribed_steps": prescribed, "initial_step": initial,
            "complete": complete, "physical_time": clock,
            "convergence": {"passed": newton_passed},
            "steps": [steps[number] for number in sorted(steps)]}


def _check_physical_sequence(steps, completed, preparation, directory, checks):
    """Check accepted clocks and observable recovery/growth/prediction behavior."""
    initial, final = preparation["initial_time_s"], preparation["final_time_s"]
    failures, previous, intervals = [], initial, []
    initial_save = steps.get(preparation["initial_step"])
    close = lambda a, b: math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)
    if initial_save and (initial_save["time_s"] is None or not close(initial_save["time_s"], initial)):
        failures.append("startup save differs from the accepted restart time")
    for number in completed:
        step = steps[number]
        time, dt = step["time_s"], step["dt_s"]
        if time is None or dt is None or not previous < time <= final + 1e-12 or dt <= 0 or not close(time-previous, dt):
            failures.append(f"step {number}: accepted physical time/dt is inconsistent")
        if time is not None:
            previous = time
        intervals.append(dt)
    complete = bool(completed) and close(previous, final)
    if checks.get("intermediate_restart") and any(close(initial, time) for time in preparation["equilibrium_times_s"]):
        failures.append("restart does not start between supplied equilibrium frames")

    log = (directory / "stdout.log").read_text(encoding="utf-8", errors="replace")
    events = _accepted_step_events(log, completed)
    for name, count in events.items():
        if checks.get(name) and count == 0:
            failures.append(f"required timestep {name} was not exercised")
    if checks.get("rejection"):
        failures.extend(_check_rejection_endpoints(log, events["rejection"], steps, completed))
    if checks.get("growth") and not any(
            a is not None and b is not None and b > a*(1+1e-8) for a, b in zip(intervals, intervals[1:])):
        failures.append("growth did not increase an accepted timestep")
    if checks.get("final_clipping"):
        failures.extend(_check_final_clipping(steps, completed, final))
    return {"complete": complete, "initial_time_s": initial, "final_time_s": final,
            "events": events, "failures": failures}


def _check_final_clipping(steps, completed, final):
    """Require the final interval to be shorter than its starting dt proposal."""
    if not completed:
        return ["no accepted final interval for clipping check"]
    last = steps[completed[-1]]
    previous = steps.get(completed[-1]-1)
    proposal = previous.get("proposed_dt_s") if previous else None
    remaining = final-previous["time_s"] if previous and previous["time_s"] is not None else None
    if (proposal is None or not math.isfinite(proposal) or remaining is None or remaining <= 0
            or proposal <= remaining+1e-12 or last["dt_s"] is None
            or not math.isclose(last["dt_s"], remaining, rel_tol=1e-10, abs_tol=1e-12)):
        return ["final interval was not clipped from its saved timestep proposal to t_final"]
    return []


def _accepted_step_events(log, completed):
    """Count observable temporal events, excluding failed predictor attempts."""
    rejection_count = log.count("TIMESTEP REJECTED (")
    growth_count = log.count("TIMESTEP INCREASED at ")
    # A prediction counts only if that attempt writes an accepted solution before
    # fallback or temporal rejection. Final aliases do not count twice.
    predicted, iteration, predicted_steps = False, None, set()
    for line in log.splitlines():
        if match := re.search(r"Time iteration\s*=\s*(\d+)", line):
            iteration, predicted = int(match[1]), False
        if "NEWTON PREDICTOR weight:" in line:
            predicted = True
        if "NEWTON PREDICTOR FALLBACK" in line or "TIMESTEP REJECTED (" in line:
            predicted = False
        if "Output written to file" in line and predicted and iteration in completed and "_NR" not in line:
            predicted_steps.add(iteration)
    events = {"rejection": rejection_count, "growth": growth_count, "prediction": len(predicted_steps)}
    return events


def _check_rejection_endpoints(log, rejection_count, steps, completed):
    """Require an accepted smaller substep after each logged temporal rejection."""
    pattern = (r"Restored accepted time:\s*([\d.Ee+\-]+)\s*s\s*"
               r"Failed endpoint:\s*([\d.Ee+\-]+)\s*s\s*Retry endpoint:\s*([\d.Ee+\-]+)")
    retries, failures = list(re.finditer(pattern, log)), []
    if len(retries) != rejection_count:
        failures.append("missing physical endpoints in timestep rejection log")
    for retry in retries:
        accepted, failed, midpoint = map(float, retry.groups())
        next_step = next((steps[n] for n in completed if steps[n]["time_s"] is not None
                          and steps[n]["time_s"] > accepted + 1e-12), None)
        if not accepted < midpoint < failed or next_step is None or next_step["time_s"] > midpoint + 1e-10:
            failures.append("rejected trial did not produce an accepted smaller substep")
    return failures


def _check_me_topology(steps, completed, topology):
    """Apply recipe endpoint expectations to completed saves only."""
    if not topology:
        return []
    if not completed:
        return ["no completed ME advance for topology expectations"]
    failures = []
    for label, number in (("first_completed", completed[0]), ("last_completed", completed[-1])):
        if label not in topology:
            continue
        actual, problem = None, None
        try:
            with h5py.File(steps[number]["output"], "r") as handle:
                values = required_array(handle, "magnetic", "topology")
            if values.size != 1:
                raise ComparisonError("magnetic/topology must contain one label")
            value = values.reshape(-1)[0]
            actual = (value.decode() if isinstance(value, bytes) else str(value)).strip(" \x00")
            if actual != topology[label]:
                problem = f"{label} topology: expected {topology[label]}, got {actual}"
        except ComparisonError as exc:
            problem = str(exc)
        steps[number]["topology"] = actual
        if problem:
            steps[number]["failures"].append(problem)
            steps[number]["status"] = "failed"
            failures.append(f"step {number}: {problem}")
    return failures


def read_checkpoint_time(handle):
    number = required_scalar(handle, "simulation_parameters/time", "Current_time_step_number")
    scale = required_scalar(handle, "simulation_parameters/adimensionalization", "time_scale")
    time = required_scalar(handle, "simulation_parameters/time", "Current_time") * scale
    dt = required_scalar(handle, "simulation_parameters/time", "Current_time_step") * scale
    return number, time, dt, scale
