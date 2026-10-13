"""Direct field/mesh contracts and the public comparison command."""

import json
import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest

from regression_tests.compare_fixed import compare_hdf5_files, validate_solution_file
from regression_tests.compare_common import read_newton_convergence
from regression_tests.compare import compare_completed_run, validate_completed_run, validate_me_completed_steps, validate_me_refinement_stop
from regression_tests.files import file_identity
from regression_tests.tests.fixtures.harness import REGRESSION_ROOT, run_command
from regression_tests.tests.fixtures.solutions import write_solution

TOLERANCES = {
    "newton_error_max": 2e-4, "mesh_coordinate_atol": 1e-12,
    "relative_l2_max": 1e-10, "normalized_linf_max": 1e-9,
}


@pytest.fixture
def files(tmp_path):
    reference, candidate = tmp_path / "reference.h5", tmp_path / "candidate.h5"
    write_solution(reference)
    write_solution(candidate)
    return reference, candidate


def test_missing_candidate_equation_names_use_reference_names(files):
    reference, candidate = files
    with h5py.File(candidate, "r+") as handle:
        del handle["simulation_parameters/physics/conservative_variable_names"]
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == "passed"
    assert report["solution"]["equation_names"] == ["rho", "Gamma"]
    assert not validate_solution_file(reference)
    assert not validate_solution_file(candidate)


@pytest.mark.parametrize("dataset,index,value,section", [
    ("solution/u", 1, 3.0, "solution"),
    ("solution/q", 0, np.nan, "solution"),
    ("mesh/T", (0, 0), 99, "mesh"),
    ("transport_1d/coefficients/d_fs", 0, 2.0, "transport_1d"),
])
def test_field_nonfinite_mesh_and_transport_failures(files, dataset, index, value, section):
    reference, candidate = files
    with h5py.File(candidate, "r+") as handle:
        handle[dataset][index] = value
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == "failed"
    assert not report[section]["passed"]
    # A finite changed field is valid data even though it fails regression;
    # bad connectivity and nonfinite fields are invalid on their own.
    assert bool(validate_solution_file(candidate)) == (dataset in {"solution/q", "mesh/T"})
    if section == "transport_1d":
        with h5py.File(candidate, "r+") as handle:
            handle[dataset][index] = np.nan
        assert any("transport_1d" in failure for failure in validate_solution_file(candidate))
    if dataset == "solution/u":
        equations = report["solution"]["datasets"]["u"]["equations"]
        assert equations["rho"]["passed"] and not equations["Gamma"]["passed"]
    elif dataset == "solution/q":
        assert not report["solution"]["datasets"]["q"]["equations"]["rho"]["finite"]
    elif dataset == "mesh/T":
        assert report["solution"]["reason"] == "mesh comparison failed"


@pytest.mark.parametrize("difference", [None, "region", "missing"])
def test_magnetic_contract(files, difference):
    reference, candidate = files
    for path in files:
        if path == candidate and difference == "missing":
            continue
        with h5py.File(path, "r+") as handle:
            magnetic = handle.create_group("magnetic")
            for name, values in {
                "topology": [b"limited"], "topology_id": [1],
                "topology_region": [1, 1, 2, 2], "rho_pol_norm": [0., .5, 1., 1.2],
                "topology_normal": [[1., 1., 1., 1.], [0., 0., 0., 0.]],
            }.items():
                magnetic[name] = np.asarray(values)
            if path == candidate and difference == "region":
                magnetic["topology_region"][2] = 3
    report = compare_hdf5_files(reference, candidate, TOLERANCES)
    assert report["status"] == ("failed" if difference else "passed")
    assert not validate_solution_file(candidate)
    if difference == "missing":
        assert "candidate is missing magnetic data" in report["failures"]
    elif difference == "region":
        assert "magnetic/topology_region differs" in report["failures"]
    else:
        with h5py.File(candidate, "r+") as handle:
            del handle["transport_1d"]  # Optional data need not exist for file validity.
        assert not validate_solution_file(candidate)
        with h5py.File(candidate, "r+") as handle:
            handle["magnetic/rho_pol_norm"][0] = np.inf
        assert any("magnetic/rho_pol_norm" in failure for failure in validate_solution_file(candidate))


@pytest.mark.parametrize("last_value,expected", [
    (None, None), ("", None), ("NaN", None), ("Infinity", None),
    ("********", None), ("1.0E309", None),
    ("1.0D-2", 1e-2), ("1.0d-5", 1e-5), ("1.0E-2", 1e-2),
])
def test_newton_acceptance_uses_last_record(tmp_path, last_value, expected):
    log = tmp_path / "stdout.log"
    log.write_text("no Newton iteration\n" if last_value is None else (
        f" Error: 1.0E-5\n\tError: {last_value}\nOutput written to file result.h5\n"
    ))
    for maximum in (2e-4, None):
        convergence = read_newton_convergence(log, maximum)
        assert convergence.final_error == expected
        assert convergence.passed == (
            expected is not None and (maximum is None or expected <= maximum)
        )
        # Failed parsing/nonfinite values must still produce a valid JSON report.
        json.dumps(convergence.as_report(), allow_nan=False)


@pytest.mark.parametrize("mode", ["same_layout", "cross_layout", "protected_output"])
def test_final_output_selection_tolerance_profile_and_cli(tmp_path, mode):
    run = tmp_path / "run"
    (run / "inputs").mkdir(parents=True)
    (run / "outputs").mkdir()
    reference, checkpoint, final = (run / name for name in (
        "inputs/reference.h5", "outputs/result_0000.h5", "outputs/result.h5",
    ))
    write_solution(reference)
    write_solution(checkpoint, solution_offset=1.0)
    write_solution(final)
    (run / "run_plan.json").write_text(json.dumps({
        "case_id": "legacy_case", "workflow_id": "baseline_warm",
        "layout_id": "serial_omp16" if mode == "cross_layout" else "mpi4_omp4",
    }))
    (run / "run_metadata.json").write_text(json.dumps({
        "status": "completed", "hdf5_outputs": ["outputs/result.h5", "outputs/result_0000.h5"],
    }))
    (run / "stdout.log").write_text(
        f"Error: 1.0E-5\nOutput written to file {checkpoint}\nOutput written to file {final}\n"
    )
    report_path = final if mode == "protected_output" else tmp_path / "comparison.json"
    original = final.read_bytes()
    if mode == "cross_layout":
        compare_completed_run(run, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json",
                              report_override=report_path)
    else:
        completed = run_command("compare", str(run), "--report", str(report_path))
        assert completed.returncode == int(mode == "protected_output"), completed.stderr
    assert final.read_bytes() == original
    if mode == "protected_output":
        assert "cannot replace a comparison input" in completed.stderr
        return
    report = json.loads(report_path.read_text())
    assert report["candidate"] == str(final)
    assert report["status"] == "passed"
    assert report["tolerance_profile"]["id"] == (
        "fixed_cross_layout" if mode == "cross_layout" else "fixed_same_layout"
    )


@pytest.fixture
def me_outputs(tmp_path):
    (tmp_path / 'outputs').mkdir()
    paths = [tmp_path / 'outputs' / name for name in ('result_0000.h5', 'result_0099.h5', 'result_0002.h5')]
    for number, path in enumerate(paths):
        write_solution(path, model=True)
        with h5py.File(path, 'r+') as handle:
            handle['solution/u'][...] = np.tile([2., .2, 4., 3., .5], 6)
            handle['simulation_parameters/adimensionalization/time_scale'] = 2e-7
            for key, value in {'Current_time_step_number': number,
                               'Current_time': number * .08 / 2e-7,
                               'Current_time_step': .08 / 2e-7}.items():
                handle[f'simulation_parameters/time/{key}'] = value
    final = tmp_path / 'outputs/result.h5'
    shutil.copy2(paths[-1], final)
    # Stale/overwritten NR files must never count as completed saves.
    nr = tmp_path / 'outputs/result_NR0001.h5'
    nr.write_bytes(b'not a completed solution')
    (tmp_path / 'stdout.log').write_text(
        f'Output written to file {paths[0]}\nTime iteration = 1\nNR iteration: 1\nError: 1.0D-6\n'
        f'NR iteration: 1\nOutput written to file {nr}\nError: 1.0D-5\nOutput written to file {paths[1]}\n'
        f'Time iteration = 2\nNR iteration: 1\nError: 2.0D-5\nOutput written to file {paths[2]}\n'
        f'Output written to file {final}\n')
    plan = {'case_id': 'legacy_case', 'workflow_id': 'baseline_warm', 'layout_id': 'mpi4_omp4',
            'moving_equilibrium_inputs': {'time_step_s': .08}, 'parameter_overrides': {'nts': 2}}
    metadata = {'status': 'completed', 'hdf5_outputs': [str(path) for path in [*paths, final, nr]]}
    return tmp_path, plan, metadata


def test_me_completion_uses_embedded_steps_and_final_attempts(me_outputs):
    run, plan, metadata = me_outputs
    report = validate_me_completed_steps(run, plan, metadata, 2e-4)
    assert report['status'] == 'passed' and report['complete']
    assert report['completed_steps'] == [1, 2] and len(report['steps']) == 3
    assert report['steps'][1]['convergence']['final_newton_error'] == 1e-5
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize('damage', ['gap', 'time', 'nonfinite_time', 'dt', 'retry_error', 'nonconvergence',
                                  'duplicate'])
def test_me_completion_rejects_invalid_completed_evidence(me_outputs, damage):
    run, plan, metadata = me_outputs
    first = Path(metadata['hdf5_outputs'][1])
    if damage == 'gap':
        metadata['hdf5_outputs'].remove(str(first))
    elif damage in ('retry_error', 'nonconvergence'):
        log = run / 'stdout.log'
        log.write_text(log.read_text().replace('Error: 1.0D-5\n', '' if damage == 'retry_error' else 'Error: 1.0D-2\n'))
    else:
        path = Path(metadata['hdf5_outputs'][3]) if damage == 'duplicate' else first
        with h5py.File(path, 'r+') as handle:
            if damage == 'nonfinite_time':
                handle['simulation_parameters/time/Current_time'][...] = np.nan
            elif damage in ('time', 'dt'):
                key = 'Current_time' if damage == 'time' else 'Current_time_step'
                handle[f'simulation_parameters/time/{key}'][...] *= 2
            else:
                handle['solution/u'][1] = -1.
    report = validate_me_completed_steps(run, plan, metadata, 2e-4)
    assert report['status'] == 'failed' and report['failures']
    if damage in ('retry_error', 'nonconvergence'):
        assert not report['convergence']['passed']
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize('compare_reference', [False, True])
@pytest.mark.parametrize('damage', ['incomplete', 'earlier_error'])
def test_completed_me_run_cannot_hide_missing_or_unconverged_steps(me_outputs, monkeypatch, compare_reference, damage):
    run, plan, metadata = me_outputs
    if damage == 'incomplete':
        metadata['hdf5_outputs'] = metadata['hdf5_outputs'][:2]
        prefix = validate_me_completed_steps(run, plan, metadata, 2e-4)
        assert prefix['status'] == 'passed' and not prefix['complete']
    else:
        log = run / 'stdout.log'
        log.write_text(log.read_text().replace('Error: 1.0D-5', 'Error: 1.0D-2'))
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'equilibrium_series'}, 'layout': 'mpi4_omp4',
                'comparison': {'method': 'fixed_hdf5', 'profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    if compare_reference:
        _, _, report = compare_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                             reference_override=Path(metadata['hdf5_outputs'][-1 if damage == 'incomplete' else 3]))
        assert report['hdf5']['status'] == 'passed'
    else:
        report = validate_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json')
    assert report['status'] == 'failed'
    if damage == 'incomplete':
        assert 'moving equilibrium did not complete all prescribed advances' in report['failures']
    else:
        assert not report['convergence']['passed']


@pytest.fixture
def stopped_me_outputs(me_outputs):
    run, plan, metadata = me_outputs
    for value in metadata['hdf5_outputs'][2:4]:
        Path(value).unlink()
    metadata['hdf5_outputs'] = [*metadata['hdf5_outputs'][:2], metadata['hdf5_outputs'][-1]]
    log = run / 'stdout.log'
    log.write_text(log.read_text().split('Time iteration = 2')[0] +
                   'Time iteration = 2\nNR iteration: 1\nError: 1.0D9\n'
                   'NR refinement recovery stopped: per-timestep refinement limit reached: 0\n')
    plan['parameter_overrides']['max_divergence_refinements'] = 0
    metadata.update(status='solver_failed', exit_code=1, launch_error=None, fatal_log_messages=[])
    metadata['output_files'] = [{'path': value, **file_identity(Path(value))} for value in metadata['hdf5_outputs']]
    return run, plan, metadata


@pytest.mark.parametrize('exit_code,status', [(0, 'solver_reported_error'), (1, 'solver_failed')])
def test_expected_refinement_stop_preserves_raw_status_and_selects_completed_prefix(stopped_me_outputs, exit_code, status):
    run, plan, metadata = stopped_me_outputs
    metadata.update(exit_code=exit_code, status=status)
    original = json.dumps(metadata)
    report = validate_me_refinement_stop(run, plan, metadata, 2e-4, completed_steps=1)
    assert report['status'] == 'passed' and report['run_status'] == status
    assert report['terminal_reason'] == 'divergence_refinement_limit'
    assert report['last_attempted_step'] == 2
    assert report['selected_output'] == metadata['hdf5_outputs'][1]
    assert json.dumps(metadata) == original


@pytest.mark.parametrize('damage', ['launch', 'crash', 'completed', 'missing_reason', 'earlier_reason', 'wrong_budget',
                                  'file_error', 'later_step', 'unconverged_prefix', 'changed_output'])
def test_expected_refinement_stop_rejects_unrelated_or_invalid_failure(stopped_me_outputs, damage):
    run, plan, metadata = stopped_me_outputs
    if damage in ('launch', 'completed'):
        metadata['status'] = 'launch_failed' if damage == 'launch' else 'completed'
    elif damage == 'crash':
        metadata['exit_code'] = 139
    elif damage == 'wrong_budget':
        plan['parameter_overrides']['max_divergence_refinements'] = 1
    elif damage == 'file_error':
        metadata['fatal_log_messages'] = ['stderr.log:1: Error opening source file: missing.h5']
    elif damage == 'changed_output':
        Path(metadata['hdf5_outputs'][1]).write_bytes(b'changed output')
    else:
        log = run / 'stdout.log'
        text = log.read_text()
        if damage == 'missing_reason':
            text = text.replace('NR refinement recovery stopped:', 'Unrelated failure:')
        elif damage == 'earlier_reason':
            stop = 'NR refinement recovery stopped: per-timestep refinement limit reached: 0\n'
            text = text.replace(stop, '').replace('Time iteration = 2\n', stop + 'Time iteration = 2\n')
        elif damage == 'later_step':
            text += 'Time iteration = 3\n'
        else:
            text = text.replace('Error: 1.0D-5', 'Error: 1.0D-2')
        log.write_text(text)
    report = validate_me_refinement_stop(run, plan, metadata, 2e-4, completed_steps=1)
    assert report['status'] == 'failed' and report['failures'] and report['selected_output'] is None
