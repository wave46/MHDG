"""Final-state overrides and cold-stage reference selection."""

import json

import pytest

from regression_tests.catalog import load_case_definition
from regression_tests import compare
from regression_tests.tests.fixtures.harness import REGRESSION_ROOT


@pytest.mark.parametrize("workflow,profile,policy_override", [
    ("cold_fixed", "race_step", None),
    ("cold_adaptive", "cold_cross_layout", "fixed_hdf5"),
])
def test_explicit_final_check_bypasses_stage_matrix(tmp_path, monkeypatch, workflow, profile, policy_override):
    run, references = fixture(tmp_path, workflow)
    calls = []
    report_path = tmp_path / "direct.json"

    def run_comparison(inputs, overrides, path, **kwargs):
        calls.append((overrides, kwargs["policy"]))
        return report_path, {"status": "passed", "failures": [], "comparison_policy": kwargs["policy"]}

    monkeypatch.setattr(compare, "compare_run", run_comparison)
    policy, path, _ = compare.compare_completed_run(
        run, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json",
        reference_override=references["time_init"], tolerance_profile_override=profile,
        comparison_policy_override=policy_override,
    )
    assert policy == calls[0][1] == "fixed_hdf5"
    assert path == report_path and len(calls) == 1
    assert calls[0][0].reference == references["time_init"]
    assert calls[0][0].tolerance_profile == profile


@pytest.mark.parametrize("workflow,fail_second", [("cold_fixed", True), ("cold_adaptive", False)])
def test_stage_references_and_stop_at_first_divergence(tmp_path, monkeypatch, workflow, fail_second):
    run, references = fixture(tmp_path, workflow)
    calls = []

    def stage_comparison(inputs, overrides, path, **kwargs):
        calls.append(overrides)
        failures = ["solution/u exceeds tolerance"] if fail_second and len(calls) == 2 else []
        return path, {"status": "failed" if failures else "passed", "failures": failures,
                      "comparison_policy": kwargs["policy"]}

    monkeypatch.setattr(compare, "compare_run", stage_comparison)
    policy, path, report = compare.compare_completed_run(
        run, REGRESSION_ROOT / "cases", REGRESSION_ROOT / "tolerances.json",
    )
    expected = 2 if fail_second else len(references)
    assert policy == "reference_matrix" and path.name == "matrix_comparison.json"
    assert report["status"] == ("failed" if fail_second else "passed")
    assert report["checked_stage_count"] == len(calls) == expected
    assert report["first_failed_stage"] == ("diffusion_reduction" if fail_second else None)
    assert [overrides.reference for overrides in calls] == list(references.values())[:expected]
    assert calls[0].newton_check == "finite_only" and calls[-1].newton_check == "bounded"
    assert calls[-1].tolerance_profile == ("fixed_stage_reference" if fail_second else "adaptive_reference")
    if workflow == "cold_adaptive":
        assert calls[-1].direct_tolerance_profile == "fixed_stage_reference"


def fixture(root, workflow_id):
    case = load_case_definition("legacy_case", REGRESSION_ROOT / "cases")
    stage_ids = [stage["stage_id"] for stage in case["workflows"][workflow_id]["stages"]]
    bundle = root / "golden"
    (bundle / "references").mkdir(parents=True)
    references, artifacts, entries = {}, {}, []

    def artifact(path, media_type):
        return {"path": path.relative_to(bundle).as_posix(), "sha256": "0" * 64,
                "size_bytes": 1, "media_type": media_type}

    for stage_id in stage_ids:
        reference = bundle / "references" / f"{stage_id}.h5"
        reference.write_text(f"golden {stage_id}\n")
        artifact_id = f"golden_{stage_id}"
        references[stage_id] = reference
        artifacts[artifact_id] = artifact(reference, "application/x-hdf5")
        entries.append({"workflow_id": workflow_id, "layout_id": "mpi4_omp4",
                        "stage_id": stage_id, "artifact_id": artifact_id})
    index = bundle / "references/index.json"
    index.write_text(json.dumps({
        "schema_version": 2, "created_utc": "2026-07-20T00:00:00Z", "case_id": "legacy_case",
        "suite_id": "cold_matrix", "suite_run_id": "golden-1",
        "source_bundle": {"bundle_id": "legacy_bundle", "bundle_version": "candidate-1"},
        "tracked_reference": {"branch": "develop", "revision": "0" * 40}, "references": entries,
    }))
    artifacts["golden_index"] = artifact(index, "application/json")
    (bundle / "manifest.json").write_text(json.dumps({
        "schema_version": 2, "bundle_id": "golden_bundle", "bundle_version": "golden-1",
        "bundle_class": "golden", "created_utc": "2026-07-20T00:00:00Z", "case_id": "legacy_case",
        "roles": {"reference_matrix": "golden_index"}, "artifacts": artifacts,
    }))
    run = root / "run"
    records = []
    for number, stage_id in enumerate(stage_ids, start=1):
        directory = run / "stages" / f"{number:02d}_{stage_id}"
        directory.mkdir(parents=True)
        candidate = directory / "candidate.h5"
        candidate.write_text(f"candidate {stage_id}\n")
        (directory / "run_metadata.json").write_text(json.dumps({"status": "completed"}))
        records.append({"stage_id": stage_id, "run_directory": str(directory),
                        "selected_hdf5": str(candidate), "status": "completed"})
    (run / "run_plan.json").write_text(json.dumps({
        "case_id": "legacy_case", "workflow_id": workflow_id, "layout_id": "mpi4_omp4",
        "bundle": {"root": str(bundle), "bundle_id": "golden_bundle", "bundle_version": "golden-1"},
        "stages": [{"stage_id": stage_id} for stage_id in stage_ids],
    }))
    (run / "run_metadata.json").write_text(json.dumps({"status": "completed", "stages": records}))
    return run, references
