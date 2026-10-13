"""Direct field/mesh contracts and the public comparison command."""

import json
import shutil
from pathlib import Path
from unittest.mock import Mock

import h5py
import numpy as np
import pytest

from regression_tests.compare_fixed import compare_hdf5_files, validate_solution_file
from regression_tests.compare_common import read_newton_convergence
from regression_tests.compare import compare_completed_run, validate_completed_run, validate_me_completed_steps, validate_me_refinement_stop, compare_me_refinement_stop
from regression_tests.files import file_identity
from regression_tests.support import ComparisonError
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


@pytest.mark.parametrize('reference_check', [False, True])
@pytest.mark.parametrize('damage', [None, 'missing', 'nonfinite', 'nonscalar', 'changed'])
def test_declared_controller_state_is_validated_and_compared(me_outputs, monkeypatch, reference_check, damage):
    run, plan, metadata = me_outputs
    final = Path(metadata['hdf5_outputs'][3])
    with h5py.File(final, 'r+') as handle:
        handle['simulation_parameters/physics/puff'] = 0.
    reference = run / 'reference.h5'
    shutil.copy2(final, reference)
    with h5py.File(final, 'r+') as handle:
        group = handle['simulation_parameters/physics']
        if damage == 'missing':
            del group['puff']
        elif damage == 'nonscalar':
            del group['puff']
            group['puff'] = [0., 1.]
        elif damage:
            group['puff'][...] = np.nan if damage == 'nonfinite' else 2e10
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'fields', 'controller_scales': {'puff': 1e20}},
                'comparison': {'method': 'fixed_hdf5', 'profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    if reference_check:
        _, _, report = compare_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                             reference_override=reference)
        assert report['hdf5']['status'] == 'passed'
    else:
        report = validate_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json')
    passed = damage is None or (damage == 'changed' and not reference_check)
    assert report['status'] == ('passed' if passed else 'failed')
    assert report['controller']['datasets']['puff']['reference_scale'] == 1e20
    assert all('controller' in failure for failure in report['failures'])


@pytest.fixture
def topology_outputs(me_outputs):
    run, plan, metadata = me_outputs
    for value in metadata['hdf5_outputs'][:-1]:
        with h5py.File(value, 'r+') as handle:
            step = int(handle['simulation_parameters/time/Current_time_step_number'][()])
            # The bootstrap label must not substitute for the first completed advance.
            handle['magnetic/topology'] = np.bytes_('limited' if step == 1 else 'lower_single_null')
    return run, plan, metadata


@pytest.mark.parametrize('compare_reference', [False, True])
@pytest.mark.parametrize('damage', [None, 'first', 'last', 'missing'])
def test_me_recipe_checks_completed_topology_endpoints(topology_outputs, monkeypatch, compare_reference, damage):
    run, plan, metadata = topology_outputs
    affected = metadata['hdf5_outputs'][2:4] if damage == 'last' else metadata['hdf5_outputs'][1:2]
    if damage:
        for value in affected:
            with h5py.File(value, 'r+') as handle:
                del handle['magnetic/topology']
                if damage != 'missing':
                    handle['magnetic/topology'] = np.bytes_('limited' if damage == 'last' else 'lower_single_null')
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'fields', 'topology': {
        'first_completed': 'limited', 'last_completed': 'lower_single_null'}},
        'comparison': {'method': 'fixed_hdf5', 'profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    if compare_reference:
        _, _, report = compare_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                             reference_override=Path(metadata['hdf5_outputs'][3]))
        assert report['hdf5']['status'] == 'passed'  # A refreshed reference cannot replace the trajectory assertion.
    else:
        report = validate_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json')
    assert report['status'] == ('passed' if damage is None else 'failed')
    assert report['convergence']['passed']
    assert all('topology' in failure for failure in report['failures'])


@pytest.mark.parametrize('reference_label', ['lower_single_null', 'limited', None, 1.5])
def test_adaptive_me_reference_requires_matching_topology(topology_outputs, monkeypatch, reference_label):
    run, plan, metadata = topology_outputs
    reference, fekete = run / 'reference.h5', run / 'fekete.h5'
    shutil.copy2(metadata['hdf5_outputs'][3], reference)
    with h5py.File(reference, 'r+') as handle:
        handle['mesh/X'][0] *= 1.01
        del handle['magnetic/topology']
        if reference_label:
            handle['magnetic/topology'] = np.bytes_(reference_label) if isinstance(reference_label, str) else reference_label
    fekete.write_bytes(b'not sampled in this metadata check')
    metadata['runtime_files'] = {'positionFeketeNodesTri2D.h5': {'path': str(fekete)}}
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'fields'},
                'comparison': {'method': 'mesh_independent', 'profile': 'adaptive_reference'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    fields = Mock(return_value={'status': 'passed', 'failures': []})
    monkeypatch.setattr('regression_tests.compare.compare_adaptive_files', fields)
    method, _, report = compare_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                              reference_override=reference)
    assert method == 'mesh_independent' and fields.call_count == 1
    assert report['status'] == ('passed' if reference_label == 'lower_single_null' else 'failed')
    assert report['magnetic_topology']['passed'] is (reference_label == 'lower_single_null')


@pytest.mark.parametrize('profile,passed', [('fixed_same_layout', False), ('feature_converged', True)])
def test_adaptive_controller_uses_existing_direct_profile(topology_outputs, monkeypatch, profile, passed):
    run, plan, metadata = topology_outputs
    candidate, reference, fekete = Path(metadata['hdf5_outputs'][3]), run / 'reference.h5', run / 'fekete.h5'
    with h5py.File(candidate, 'r+') as handle:
        handle['simulation_parameters/physics/puff'] = 5e11
    shutil.copy2(candidate, reference)
    with h5py.File(reference, 'r+') as handle:
        handle['mesh/X'][0] *= 1.01
        handle['simulation_parameters/physics/puff'][...] = 0.
    fekete.write_bytes(b'not sampled in this scalar check')
    metadata['runtime_files'] = {'positionFeketeNodesTri2D.h5': {'path': str(fekete)}}
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'fields', 'controller_scales': {'puff': 1e20}},
                'comparison': {'method': 'mesh_independent', 'profile': 'adaptive_reference', 'direct_profile': profile}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    monkeypatch.setattr('regression_tests.compare.compare_adaptive_files',
                        lambda *_, **__: {'status': 'passed', 'failures': []})
    method, _, report = compare_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                              reference_override=reference)
    assert method == 'mesh_independent'
    assert report['status'] == ('passed' if passed else 'failed')
    assert report['controller']['status'] == report['status']


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


@pytest.mark.parametrize('label', ['limited', 'lower_single_null'])
def test_expected_prefix_checks_recipe_topology(stopped_me_outputs, monkeypatch, label):
    from regression_tests.compare import validate_me_stopped_run
    run, plan, metadata = stopped_me_outputs
    with h5py.File(metadata['hdf5_outputs'][1], 'r+') as handle:
        handle['magnetic/topology'] = np.bytes_(label)
    metadata['output_files'] = [{'path': value, **file_identity(Path(value))} for value in metadata['hdf5_outputs']]
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'fields', 'topology': {
        'first_completed': 'limited', 'last_completed': 'limited'}},
        'comparison': {'method': 'fixed_hdf5', 'profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    report = validate_me_stopped_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json', completed_steps=1)
    assert report['status'] == ('passed' if label == 'limited' else 'failed')
    assert report['run_status'] == 'solver_failed'


@pytest.mark.parametrize('difference,passed', [(5e9, True), (2e10, False)])
def test_expected_prefix_compares_controller_scalars_near_zero(stopped_me_outputs, monkeypatch, difference, passed):
    from regression_tests.compare import validate_me_stopped_run
    run, plan, metadata = stopped_me_outputs
    candidate = Path(metadata['hdf5_outputs'][1])
    with h5py.File(candidate, 'r+') as handle:
        handle['simulation_parameters/physics/puff'] = 0.
    reference = run / 'reference.h5'
    shutil.copy2(candidate, reference)
    with h5py.File(candidate, 'r+') as handle:
        handle['simulation_parameters/physics/puff'][...] = difference
    metadata['output_files'] = [{'path': value, **file_identity(Path(value))} for value in metadata['hdf5_outputs']]
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'fields', 'controller_scales': {'puff': 1e20}},
                'comparison': {'method': 'fixed_hdf5', 'profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    validation = validate_me_stopped_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json', completed_steps=1)
    assert validation['status'] == 'passed'
    _, _, report = compare_me_refinement_stop(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                              completed_steps=1, reference_override=reference)
    assert report['status'] == ('passed' if passed else 'failed')
    assert report['candidate'] == str(candidate) and report['expected_failure']['run_status'] == 'solver_failed'


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


@pytest.mark.parametrize('policy', ['fixed_hdf5', 'mesh_independent'])
@pytest.mark.parametrize('reference_damage', [None, 'field', 'step', 'time', 'dt'])
def test_expected_stop_prefix_uses_shared_reference_comparison(stopped_me_outputs, monkeypatch, policy, reference_damage):
    run, plan, metadata = stopped_me_outputs
    (run / 'inputs').mkdir()
    reference = run / 'inputs/reference.h5'
    shutil.copy2(metadata['hdf5_outputs'][1], reference)
    if reference_damage:
        with h5py.File(reference, 'r+') as handle:
            if reference_damage == 'field':
                handle['solution/u'][1] += .1
            else:
                key = {'step': 'Current_time_step_number', 'time': 'Current_time', 'dt': 'Current_time_step'}[reference_damage]
                handle[f'simulation_parameters/time/{key}'][...] *= 2
    (run / 'run_plan.json').write_text(json.dumps(plan))
    metadata_path = run / 'run_metadata.json'
    metadata_path.write_text(json.dumps(metadata))
    original = metadata_path.read_bytes()
    workflow = {'moving_equilibrium': {'series': 'equilibrium_series'}, 'layout': 'mpi4_omp4',
                'comparison': {'method': policy, 'profile': 'fixed_same_layout' if policy == 'fixed_hdf5' else 'adaptive_reference',
                               'direct_profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    method, _, report = compare_me_refinement_stop(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                                  completed_steps=1)
    assert method == 'fixed_hdf5'  # Existing matching-mesh fallback is reused.
    assert report['status'] == ('passed' if reference_damage is None else 'failed')
    assert report['candidate'] == metadata['hdf5_outputs'][1]
    assert report['convergence']['final_newton_error'] == 1e-5  # Never the failed step's 1e9.
    assert report['expected_failure']['run_status'] == 'solver_failed'
    assert metadata_path.read_bytes() == original
    if reference_damage in ('step', 'time', 'dt'):
        assert report['hdf5']['status'] == 'passed' and any('ME reference' in failure for failure in report['failures'])
    with pytest.raises(ComparisonError, match='not completed'):
        compare_completed_run(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json')
    if reference_damage is None and policy == 'fixed_hdf5':
        initial = Path(metadata['hdf5_outputs'][0])
        original_output = initial.read_bytes()
        with pytest.raises(ComparisonError, match='cannot replace a comparison input'):
            compare_me_refinement_stop(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json',
                                       completed_steps=1, report_override=initial)
        assert initial.read_bytes() == original_output


def test_invalid_expected_stop_cannot_enter_reference_comparison(stopped_me_outputs, monkeypatch):
    run, plan, metadata = stopped_me_outputs
    metadata['exit_code'] = 139
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    workflow = {'moving_equilibrium': {'series': 'equilibrium_series'}, 'layout': 'mpi4_omp4',
                'comparison': {'method': 'fixed_hdf5', 'profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    monkeypatch.setattr('regression_tests.compare.compare_run', lambda *_, **__: pytest.fail('invalid outcome reached comparator'))
    with pytest.raises(ComparisonError, match='invalid expected ME refinement stop'):
        compare_me_refinement_stop(run, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json', completed_steps=1)


@pytest.mark.parametrize('damage', [None, 'reference', 'output', 'crash', 'completed', 'log'])
def test_suite_assesses_and_reuses_only_declared_expected_stops(stopped_me_outputs, monkeypatch, damage):
    from regression_tests import suites

    run, plan, metadata = stopped_me_outputs
    (run / 'inputs').mkdir()
    reference = run / 'inputs/reference.h5'
    shutil.copy2(metadata['hdf5_outputs'][1], reference)
    if damage in ('reference', 'output'):
        path = reference if damage == 'reference' else Path(metadata['hdf5_outputs'][1])
        with h5py.File(path, 'r+') as handle:
            handle['solution/u'][1] += .1
    elif damage == 'crash':
        metadata['exit_code'] = 139
    elif damage == 'completed':
        metadata.update(status='completed', exit_code=0)
    elif damage == 'log':
        (run / 'stdout.log').write_text('unrelated failure')
    (run / 'run_plan.json').write_text(json.dumps(plan))
    metadata_path = run / 'run_metadata.json'
    metadata_path.write_text(json.dumps(metadata))
    original = metadata_path.read_bytes()
    expectation = {'reason': 'divergence_refinement_limit', 'completed_steps': 1}
    workflow = {'moving_equilibrium': {'series': 'fields', 'expected_failure': expectation},
                'comparison': {'method': 'fixed_hdf5', 'profile': 'fixed_same_layout'}}
    case = {'workflows': {'baseline_warm': workflow}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition', lambda *_, **__: case)
    monkeypatch.setattr(suites, 'load_case_definition', lambda *_, **__: case)
    source = {'workflow_id': 'baseline_warm', 'layout_id': 'mpi4_omp4',
              'run_status': metadata['status'], 'run_directory': str(run)}
    arguments = (source, REGRESSION_ROOT / 'cases', REGRESSION_ROOT / 'tolerances.json')
    result = suites.assess_result(*arguments, expected_failure=expectation)
    assert result['status'] == ('passed' if damage is None else 'failed')
    assert result['run_status'] == metadata['status']
    reusable = damage in (None, 'reference')  # Reference disagreement never requires a solver replay.
    assert suites._reusable_result(*arguments, expected_failure=expectation) is reusable
    validation = suites.assess_result(*arguments, references=False, expected_failure=expectation)
    assert validation['status'] == ('passed' if reusable else 'failed')
    if damage is None:
        assert result['expected_failure']['last_attempted_step'] == 2
        assert suites.assess_result(*arguments)['status'] == 'failed'  # Undeclared failures stay failures.
    summary_path = run / 'suite_summary.json'
    summary_path.write_text(json.dumps({'schema_version': 2, 'suite_id': 'example', 'run_id': 'example',
                                        'case_id': 'legacy_case', 'reference_comparisons': True, 'results': [source]}))
    from regression_tests import diagnostics
    check_output = Mock(wraps=diagnostics.check_output)
    monkeypatch.setattr(diagnostics, 'check_output', check_output)
    verified = suites.verify_suite(summary_path, arguments[1], arguments[2], include_layout_pairs=False)[1]
    assert verified['results'][0]['status'] == result['status']
    assert check_output.call_count == int(reusable or damage == 'completed')
    if reusable:
        assert check_output.call_args.args[0] == Path(metadata['hdf5_outputs'][1])
    assert metadata_path.read_bytes() == original


def test_prefix_comparison_applies_stricter_fallback_newton_profile_to_earlier_steps(me_outputs, monkeypatch):
    run, plan, metadata = me_outputs
    plan['parameter_overrides'].update(nts=3, max_divergence_refinements=0)
    log = run / 'stdout.log'
    log.write_text(log.read_text().replace('Error: 2.0D-5', 'Error: 1.0D-7') +
                   'Time iteration = 3\nNR iteration: 1\nError: 1.0D9\n'
                   'NR refinement recovery stopped: per-timestep refinement limit reached: 0\n')
    metadata.update(status='solver_failed', exit_code=1)
    metadata['output_files'] = [{'path': value, **file_identity(Path(value))} for value in metadata['hdf5_outputs']]
    (run / 'run_plan.json').write_text(json.dumps(plan))
    (run / 'run_metadata.json').write_text(json.dumps(metadata))
    tolerances = json.loads((REGRESSION_ROOT / 'tolerances.json').read_text())
    tolerances['profiles']['fixed_same_layout']['newton_error_max'] = 1e-6
    path = run / 'tolerances.json'
    path.write_text(json.dumps(tolerances))
    workflow = {'moving_equilibrium': {'series': 'equilibrium_series'}, 'layout': 'mpi4_omp4',
                'comparison': {'method': 'mesh_independent', 'profile': 'adaptive_reference', 'direct_profile': 'fixed_same_layout'}}
    monkeypatch.setattr('regression_tests.compare.load_case_definition',
                        lambda *_, **__: {'workflows': {'baseline_warm': workflow}})
    _, _, report = compare_me_refinement_stop(run, REGRESSION_ROOT / 'cases', path,
                                             completed_steps=2, reference_override=Path(metadata['hdf5_outputs'][2]))
    assert report['hdf5']['status'] == 'passed' and report['convergence']['final_newton_error'] == 1e-7
    assert report['status'] == 'failed' and not report['convergence']['passed']
    assert report['moving_equilibrium']['steps'][1]['status'] == 'failed'
