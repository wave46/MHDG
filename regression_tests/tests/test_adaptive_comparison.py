"""Adaptive sampling, field tolerances and run-level convergence."""

import json

import h5py
import numpy as np

from regression_tests import compare
from regression_tests.compare_adaptive import SampledFields, compare_sampled_fields, reference_sample_points
from tests.fixtures.harness import REGRESSION_ROOT


def test_run_overrides_and_finite_only_convergence(tmp_path, monkeypatch):
    candidate, reference, fekete = (tmp_path / name for name in (
        "candidate.h5", "golden.h5", "positionFeketeNodesTri2D.h5",
    ))
    for path in (candidate, reference, fekete):
        path.write_text("fixture\n")
    (tmp_path / "stdout.log").write_text("Error: 1.0\n")
    (tmp_path / "run_plan.json").write_text(json.dumps({
        "case_id": "legacy_case", "workflow_id": "cold_adaptive", "layout_id": "mpi4_omp4",
    }))
    (tmp_path / "run_metadata.json").write_text(json.dumps({
        "status": "completed", "hdf5_outputs": ["candidate.h5"],
        "runtime_files": {fekete.name: {"path": str(fekete)}},
    }))
    calls = []

    def fields(*args):
        calls.append(args)
        return {"status": "passed", "failures": [], "tolerances": {}}

    monkeypatch.setattr(compare, "compare_adaptive_files", fields)
    inputs = compare.load_comparison_inputs(tmp_path, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json")
    path, report = compare.compare_run(inputs, compare.ComparisonOverrides(
        candidate, reference, "adaptive_reference", "finite_only",
    ))
    assert report["status"] == "passed"
    assert report["convergence"] == {"passed": True, "final_newton_error": 1., "maximum": None}
    assert report["tolerance_profile"]["id"] == "adaptive_reference"
    assert path.name == "comparison.json"
    assert calls[0][:4] == (reference, candidate, fekete, 4)


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
