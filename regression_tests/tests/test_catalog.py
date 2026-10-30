"""Behavior of shared workflow composition, independent of production cases."""

import json
import shutil
from pathlib import Path

import pytest


from regression_tests.catalog import load_case_definition
from regression_tests.support import BundleError


REGRESSION_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("case_id", ["legacy_case", "diverted_case", "moving_equilibrium"])
def test_production_workflows_inherit_solver_defaults(case_id):
    case = load_case_definition(case_id, REGRESSION_ROOT / "cases")
    policy = {"stab": 2, "newton_formulation": "absolute", "newton_convergence": "update",
              "residual_diagnostics_mode": "off", "bohm_momentum_always_include_diffusion": True,
              "flux_limiter": False}
    for name, workflow in case["workflows"].items():
        expected = dict(policy)
        if name == "heat_flux_limiter":
            expected["flux_limiter"] = True
        if name.startswith("newton_"):
            expected["residual_diagnostics_mode"] = "equations"
            expected["newton_formulation"] = name.removeprefix("newton_")
        overrides = workflow["parameter_overrides"]
        assert {key: overrides[key] for key in expected} == expected, name
        assert all(not set(stage.get("parameter_overrides", {})) & set(policy)
                   for stage in workflow.get("stages", [])), name


@pytest.mark.parametrize("case_id", ["legacy_case", "diverted_case"])
def test_newton_pair_uses_same_unlimited_restart_and_controls(case_id):
    workflows = load_case_definition(case_id, REGRESSION_ROOT / "cases")["workflows"]
    absolute, incremental = (workflows[name] for name in ("newton_absolute", "newton_incremental"))
    assert absolute["restart"] == incremental["restart"] == workflows["baseline_warm"]["restart"]
    controls = [dict(workflow["parameter_overrides"]) for workflow in (absolute, incremental)]
    assert controls[0].pop("newton_formulation") == "absolute"
    assert controls[1].pop("newton_formulation") == "incremental"
    assert controls[0] == controls[1]
    assert controls[0]["flux_limiter"] is False
    warm_controls = dict(workflows["baseline_warm"]["parameter_overrides"])
    warm_controls.pop("newton_formulation")
    assert controls[0] == {**warm_controls, "residual_diagnostics_mode": "equations", "saveNR": True}


@pytest.fixture
def catalog(tmp_path):
    (tmp_path / "schemas").mkdir()
    shutil.copy2(
        REGRESSION_ROOT / "schemas/case.schema.json",
        tmp_path / "schemas/case.schema.json",
    )
    cases = tmp_path / "cases"
    cases.mkdir()
    shared = {
        "schema_version": 2,
        "sequences": {
            "bootstrap": [
                {"id": "initial", "parameters": "parameters", "transport": "transport"},
                {"id": "settle", "parameters": "parameters", "transport": "transport",
                 "parameter_overrides": {"tNR": 0.001}},
            ],
        },
        "workflows": {
            "cold_fixed": {
                "description": "Fixed bootstrap",
                "type": "staged_fixed_mesh",
                "mesh": "mesh",
                "stages": ["bootstrap"],
                "parameter_overrides": {"nrp": 40, "rest_adapt": False},
            },
            "cold_adaptive": {
                "extends": "cold_fixed",
                "description": "Adaptive bootstrap",
                "type": "staged_adaptive_mesh",
                "adaptive_stages": ["initial"],
            },
        },
    }
    (tmp_path / "workflows.json").write_text(json.dumps(shared))
    declaration = {
        "schema_version": 2,
        "description": "An independently supplied case",
        "files": {"required": {
            "mesh": "third_case.msh", "parameters": "third_case.txt",
            "transport": "third_case.nml",
        }},
        "workflows": {"cold_fixed": {}, "cold_adaptive": {}},
    }
    return cases, declaration


def load_case(catalog, name="third_case"):
    cases, declaration = catalog
    (cases / f"{name}.json").write_text(json.dumps(declaration))
    return load_case_definition(name, cases)


def test_case_overrides_reach_derived_workflows_without_changing_shared_inputs(catalog):
    cases, declaration = catalog
    original_shared = (cases.parent / "workflows.json").read_bytes()
    declaration["workflows"]["cold_fixed"] = {
        "parameter_overrides": {"nrp": 12},
        "parameter_namelists": {"nrp": "numerics"},
        "stage_overrides": {"settle": {
            "parameter_overrides": {"tNR": 0.0001},
            "parameter_namelists": {"tNR": "convergence"},
        }},
    }
    declaration["parameter_namelists"] = {"new_knob": "case_physics"}
    case = load_case(catalog)
    adaptive = case["workflows"]["cold_adaptive"]
    assert adaptive["parameter_overrides"] == {"nrp": 12, "rest_adapt": False}
    assert adaptive["parameter_namelists"] == {"new_knob": "case_physics", "nrp": "numerics"}
    assert adaptive["stages"][1]["parameter_namelists"] == {"tNR": "convergence"}
    assert adaptive["stages"][0]["parameter_overrides"] == {"rest_adapt": True}
    assert adaptive["stages"][1]["parameter_overrides"] == {"tNR": 0.0001}

    declaration["files"]["required"]["mesh"] = "another_case.msh"
    declaration["workflows"]["cold_fixed"] = {}
    other = load_case(catalog, "another_case")
    assert other["workflows"]["cold_adaptive"]["stages"][1]["parameter_overrides"] == {
        "tNR": 0.001,
    }
    assert (cases.parent / "workflows.json").read_bytes() == original_shared


def test_custom_workflow_can_extend_an_unselected_shared_parent(catalog):
    catalog[1]["workflows"] = {
        "custom": {"extends": "cold_adaptive", "parameter_overrides": {"nrp": 2}},
    }
    case = load_case(catalog)
    assert list(case["workflows"]) == ["custom"]
    assert case["workflows"]["custom"]["parameter_overrides"]["nrp"] == 2


@pytest.mark.parametrize("workflows, message", [
    ({"a": {"extends": "b"}, "b": {"extends": "a"}}, "inheritance cycle"),
    ({"a": {"extends": "missing"}}, "extends unknown workflow"),
    ({"a": {}}, "resolved workflow is missing"),
    ({"cold_fixed": {"stages": ["missing"]}}, "unknown stage sequence"),
    ({"cold_fixed": {"stage_overrides": {"missing": {}}}}, "unknown stages"),
    ({"cold_fixed": {"mesh": "missing"}}, "undefined artifact roles"),
])
def test_invalid_composition_fails_before_preparation(catalog, workflows, message):
    catalog[1]["workflows"] = workflows
    with pytest.raises(BundleError, match=message):
        load_case(catalog)
