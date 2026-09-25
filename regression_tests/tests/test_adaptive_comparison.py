"""Adaptive sampling, field tolerances and run-level convergence."""

import json

import h5py
import numpy as np
import pytest

from regression_tests import compare
from regression_tests.compare_adaptive import SampledFields, compare_sampled_fields, reference_sample_points
from regression_tests.tests.fixtures.harness import REGRESSION_ROOT
from regression_tests.tests.fixtures.solutions import write_solution
from regression_tests.support import ComparisonError


def test_run_overrides_and_finite_only_convergence(tmp_path, monkeypatch):
    candidate, reference, fekete = (tmp_path / name for name in (
        "candidate.h5", "golden.h5", "positionFeketeNodesTri2D.h5",
    ))
    for path in (candidate, reference):
        write_solution(path)
    with h5py.File(candidate, "r+") as handle:
        handle["mesh/X"][0] *= 1.01
    fekete.write_text("fixture\n")
    (tmp_path / "stdout.log").write_text("Error: 1.0\n")
    (tmp_path / "run_plan.json").write_text(json.dumps({
        "case_id": "legacy_case", "workflow_id": "bootstrap_adaptive", "layout_id": "mpi4_omp4",
        "stages": [{"stage_id": "final"}],
    }))
    stage = tmp_path / "final-stage"
    stage.mkdir()
    (stage / "run_metadata.json").write_text(json.dumps({
        "runtime_files": {fekete.name: {"path": str(fekete)}},
    }))
    (tmp_path / "run_metadata.json").write_text(json.dumps({
        "status": "completed", "hdf5_outputs": ["candidate.h5"],
        "stages": [{"stage_id": "final", "run_directory": str(stage), "status": "completed"}],
    }))
    calls = []

    def fields(*args):
        calls.append(args)
        return {"status": "passed", "failures": [], "tolerances": {}}

    monkeypatch.setattr(compare, "compare_adaptive_files", fields)
    inputs = compare.load_comparison_inputs(tmp_path, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json")
    inputs.workflow["stages"] = [{"id": "final", "newton_check": "finite_only"}]
    monkeypatch.setattr(compare, "load_comparison_inputs", lambda *args, **kwargs: inputs)
    _, _, report = compare.compare_completed_run(
        tmp_path, inputs.case_directory, inputs.tolerances_path,
        candidate_override=candidate, reference_override=reference,
    )
    assert report["status"] == "passed"
    assert report["convergence"] == {"passed": True, "final_newton_error": 1., "maximum": None}
    assert report["tolerance_profile"]["id"] == "adaptive_reference"
    assert calls[0][:4] == (reference, candidate, fekete, 4)
    assert report["comparison_policy"] == "mesh_independent"
    assert report["method_selection"]["differences"] == ["X"]
    inputs.workflow["stages"][-1]["newton_check"] = "bounded"
    _, _, report = compare.compare_completed_run(
        tmp_path, inputs.case_directory, inputs.tolerances_path,
        candidate_override=candidate, reference_override=reference,
    )
    assert report["status"] == "failed"
    assert "final Newton error exceeds tolerance" in report["failures"]


def test_reference_points_are_deterministic_and_interior(tmp_path):
    path = tmp_path / "reference.h5"
    with h5py.File(path, "w") as handle:
        mesh = handle.create_group("mesh")
        mesh["X"] = np.array([[0., 1., 0.], [0., 0., 1.]])
        mesh["Tlin"] = np.array([[1], [2], [3]])
    points = reference_sample_points(path, 4)
    assert points.shape == (4, 2)
    assert (points > 0).all() and (points.sum(axis=1) < 1).all()
    np.testing.assert_allclose(points[0], [1/3, 1/3])


def test_gradient_tolerances_and_point_coverage():
    def fields(offset, inside):
        solution = np.array([[1., 2.], [2., 4.], [3., 6.]]) + offset
        return SampledFields(["rho", "Gamma"], np.array(inside), solution,
                             np.stack((solution * .1, solution * .2), axis=2))

    report = compare_sampled_fields(
        fields(0., [True] * 3), fields(.1, [True, True, False]), 1,
        {"minimum_point_coverage": 1.,
         "solution": {"relative_l2_max": 1., "normalized_linf_max": 1.},
         "gradient": {"relative_l2_max": 1e-3, "normalized_linf_max": 1e-3}},
    )
    assert report["status"] == "failed"
    assert report["sampling"]["common_coverage"] == 2/3
    assert "common point coverage is below tolerance" in report["failures"]
    assert report["datasets"]["solution"]["equations"]["rho"]["passed"]
    assert any("gradient_x/rho exceeds tolerance" in item for item in report["failures"])


@pytest.fixture
def adaptive_run(tmp_path, monkeypatch):
    reference, candidate = tmp_path / "reference.h5", tmp_path / "candidate.h5"
    write_solution(reference)
    write_solution(candidate)
    (tmp_path / "stdout.log").write_text("Error: 1e-5\n")
    (tmp_path / "run_plan.json").write_text(json.dumps({
        "case_id": "legacy_case", "workflow_id": "bootstrap_adaptive", "layout_id": "mpi4_omp4",
    }))
    (tmp_path / "run_metadata.json").write_text(json.dumps({
        "status": "completed", "hdf5_outputs": [candidate.name],
    }))
    inputs = compare.load_comparison_inputs(tmp_path, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json")
    # These file-comparison checks are independent of cold stage sequencing.
    inputs.workflow.pop("stages")
    monkeypatch.setattr(compare, "load_comparison_inputs", lambda *args, **kwargs: inputs)
    return inputs, reference, candidate


def test_refinement_requires_recorded_initial_output_and_element_growth(adaptive_run):
    inputs, initial, candidate = adaptive_run
    inputs.workflow["require_refinement"] = True
    inputs.metadata["hdf5_outputs"].append(initial.name)
    (inputs.run_directory / "stdout.log").write_text(
        "Mesh converted from gmsh to hdf5. Output written to file inputs/mesh_1_4.h5\n"
        f"Output written to file {initial}\nError: 1e-5\n"
        "Mesh converted from gmsh to hdf5. Output written to file ./res/new_mesh_n1.h5\n"
        f"Output written to file {candidate}\n")
    assert compare.select_candidate(inputs.run_directory, inputs.metadata) == candidate
    assert compare._check_refinement(inputs, candidate)["status"] == "failed"  # Identical meshes.
    with h5py.File(candidate, "r+") as handle:
        handle["mesh/Nelems"][...] = 3
    report = compare._check_refinement(inputs, candidate)
    assert report["status"] == "passed"
    assert (report["initial_elements"], report["final_elements"]) == (2, 3)
    inputs.metadata["hdf5_outputs"].remove(initial.name)
    assert compare._check_refinement(inputs, candidate)["status"] == "failed"


@pytest.mark.parametrize("perturbed", [False, True])
def test_matching_mesh_uses_direct_check_without_fallback(adaptive_run, monkeypatch, perturbed):
    inputs, reference, candidate = adaptive_run
    with h5py.File(candidate, "r+") as handle:
        handle["mesh/X"][0, 0] += 1e-13  # Within the declared mesh tolerance.
        if perturbed:
            handle["solution/u"][1] += .001
    def no_interpolation(*args):
        pytest.fail("matching-mesh comparisons must not interpolate, even after a failed direct check")
    monkeypatch.setattr(compare, "compare_adaptive_files", no_interpolation)
    policy, path, report = compare.compare_completed_run(
        inputs.run_directory, inputs.case_directory, inputs.tolerances_path, reference_override=reference,
    )
    assert policy == "fixed_hdf5"
    assert report["status"] == ("failed" if perturbed else "passed")
    assert report["method_selection"]["reason"] == "matching discrete meshes"
    assert report["tolerance_profile"]["id"] == "cold_fixed_reference"
    assert "hdf5" in report and "sampling" not in report
    assert json.loads(path.read_text())["comparison_policy"] == policy


@pytest.mark.parametrize("corruption", ["nonfinite", "index", "fractional_index", "degenerate", "order", "field_size"])
def test_malformed_mesh_or_storage_never_interpolates(adaptive_run, monkeypatch, corruption):
    inputs, reference, candidate = adaptive_run
    with h5py.File(candidate, "r+") as handle:
        if corruption == "nonfinite":
            handle["mesh/X"][0, 0] = np.nan
        elif corruption == "index":
            handle["mesh/Tlin"][0, 0] = 99
        elif corruption == "fractional_index":
            data = handle["mesh/Tlin"][()].astype(float)
            del handle["mesh/Tlin"]
            handle["mesh/Tlin"] = data + .1
        elif corruption == "degenerate":
            handle["mesh/X"][1] = 0.
        elif corruption == "order":
            handle["mesh/Nnodesperface"][0] = 3
        else:
            del handle["solution/u"]
            handle["solution/u"] = [1.]
    monkeypatch.setattr(compare, "compare_adaptive_files", lambda *args: pytest.fail("invalid mesh fell back to interpolation"))
    with pytest.raises(ComparisonError):
        compare.compare_run(inputs, compare.ComparisonOverrides(reference=reference))


def test_face_order_is_part_of_mesh_identity_and_parallel_checks_stay_strict(adaptive_run, monkeypatch):
    inputs, reference, candidate = adaptive_run
    with h5py.File(candidate, "r+") as handle:
        faces = handle["mesh/F"][()]
        # Permute two boundary-face IDs and their records consistently.
        faces[faces == 2] = 99
        faces[faces == 3] = 2
        faces[faces == 99] = 3
        handle["mesh/F"][...] = faces
        for name in ("Tb", "extfaces"):
            data = handle[f"mesh/{name}"][()]
            data[:, [0, 1]] = data[:, [1, 0]]
            handle[f"mesh/{name}"][...] = data
    differences = compare.mesh_differences(reference, candidate, 1e-12)
    assert "F" in differences and "extfaces" in differences
    monkeypatch.setattr(compare, "compare_adaptive_files", lambda *args: pytest.fail("parallel check interpolated"))
    _, report = compare.compare_run(inputs, compare.ComparisonOverrides(
        reference=reference, tolerance_profile="cold_fixed_reference",
    ), policy="fixed_hdf5")
    assert report["status"] == "failed"
    assert report["tolerance_profile"]["id"] == "cold_fixed_reference"
