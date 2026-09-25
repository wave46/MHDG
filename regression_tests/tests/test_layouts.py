"""Generated relations and suite selection, independent of production suites."""

import json
import shutil
from pathlib import Path

import pytest


from regression_tests.catalog import layout_pairs, load_layouts
from regression_tests.catalog import load_suite_definition
from regression_tests.support import BundleError


REGRESSION_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def catalog(tmp_path):
    (tmp_path / "schemas").mkdir()
    for name in ("layouts", "suites"):
        shutil.copy2(
            REGRESSION_ROOT / f"schemas/{name}.schema.json",
            tmp_path / f"schemas/{name}.schema.json",
        )
    # Different counts and ordering from production: relations follow execution
    # properties, not hardcoded names or the order of their baseline/candidate.
    (tmp_path / "layouts.json").write_text(json.dumps({
        "schema_version": 2,
        "layouts": ["mpi2_omp3", "serial_omp3", "mpi2_omp1", "serial_omp1"],
    }))
    return tmp_path, load_layouts(tmp_path / "layouts.json")


def test_relations_isolate_execution_changes_and_deduplicate_pairs(catalog):
    _, layouts = catalog
    pairs = layout_pairs(layouts, ["openmp", "mpi", "hybrid"])
    assert [(p["baseline"], p["candidate"]) for p in pairs] == [
        ("serial_omp1", "serial_omp3"),
        ("serial_omp1", "mpi2_omp1"),
        ("mpi2_omp1", "mpi2_omp3"),
    ]
    all_pairs = layout_pairs(layouts, ["all_pairs", "openmp", "mpi", "hybrid"])
    identities = {frozenset(p.values()) for p in all_pairs}
    assert len(all_pairs) == len(identities) == 6
    assert all(len(pair) == 2 for pair in identities)


@pytest.mark.parametrize("relation, selected", [
    ("openmp", ["serial_omp3"]),
    ("mpi", ["mpi2_omp1", "mpi2_omp3"]),
    ("hybrid", ["mpi2_omp1", "mpi3_omp3"]),
])
def test_requested_relation_cannot_silently_lose_its_comparison(catalog, relation, selected):
    _, layouts = catalog
    layouts["mpi3_omp3"] = {**layouts["mpi2_omp3"], "mpi_ranks": 3}
    with pytest.raises(BundleError, match="has no matching layouts"):
        layout_pairs({name: layouts[name] for name in selected}, [relation])


def test_suite_defaults_overrides_and_relations_select_runs_once(catalog, monkeypatch):
    root, layouts = catalog
    monkeypatch.setattr(
        "regression_tests.catalog.load_case_definition",
        lambda *_, **__: {"workflows": {"probe": {}}},
    )
    declaration = {"description": "Probe", "workflows": ["probe"]}

    def load():
        (root / "suites.json").write_text(json.dumps({
            "schema_version": 2,
            "defaults": {"case": "example", "layout": "mpi2_omp3"},
            "suites": {"probe": declaration},
        }))
        return load_suite_definition(
            "probe", root / "suites.json", root / "layouts.json", root / "cases",
        )

    suite = load()
    assert suite["case_id"] == "example"
    assert suite["layouts"] == ["mpi2_omp3"]
    assert suite["reference_comparisons"]
    declaration["layouts"] = "all"
    assert load()["layouts"] == list(layouts)

    declaration.update({
        "case": "another_case", "relations": ["mpi", "hybrid"],
        "tolerance_profile": "probe_tolerance",
    })
    suite = load()
    assert suite["case_id"] == "another_case"
    assert suite["layouts"] == ["serial_omp1", "mpi2_omp1", "mpi2_omp3"]
    assert not suite["reference_comparisons"]
    declaration["reference_comparisons"] = True
    assert load()["reference_comparisons"]

    declaration["layouts"] = ["serial_omp1", "mpi2_omp1"]
    declaration["relations"] = ["mpi"]
    assert load()["layouts"] == ["serial_omp1", "mpi2_omp1"]
    declaration["layouts"] = ["unknown_layout"]
    with pytest.raises(BundleError, match="unknown layouts"):
        load()


def test_profiles_compose_cases_without_duplicate_runs(catalog, monkeypatch):
    from regression_tests.catalog import load_selection

    root, _ = catalog
    monkeypatch.setattr("regression_tests.catalog.load_case_definition",
                        lambda *_, **__: {"workflows": {"probe": {}}})
    document = {
        "schema_version": 2, "defaults": {"case": "first", "layout": "mpi2_omp3"},
        "suites": {"probe": {"description": "Probe", "workflows": ["probe"]}},
        "profiles": {
            "daily": {"description": "Daily", "checks": [{"suite": "probe"}]},
            "all": {"description": "All", "include": ["daily"],
                    "checks": [{"suite": "probe"}, {"suite": "probe", "case": "second"}]},
        },
    }
    def select():
        (root / "suites.json").write_text(json.dumps(document))
        return load_selection("all", root / "suites.json", root / "layouts.json", root / "cases")
    is_profile, checks = select()
    assert is_profile
    assert [(item["suite_id"], item["case_id"]) for item in checks] == [("probe", "first"), ("probe", "second")]
    document["profiles"]["daily"]["include"] = ["all"]
    with pytest.raises(BundleError, match="cyclic profile inclusion"):
        select()
