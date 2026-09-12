"""Behavior of shared workflow composition, independent of production cases."""

import json
import shutil
from pathlib import Path

import pytest


from regression_tests.catalog import load_case_definition
from regression_tests.support import BundleError


REGRESSION_ROOT = Path(__file__).resolve().parents[1]


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
        "reference": {"branch": "develop", "revision": "a" * 40},
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
    assert [stage["restart_from"] for stage in adaptive["stages"]] == [
        "analytical", "previous_stage",
    ]
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
