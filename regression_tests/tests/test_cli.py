"""Smoke checks for the public package, without running the scientific solver."""

import json
import subprocess
import sys
from pathlib import Path


from regression_tests.tests.fixtures.harness import create_harness, run_command


ROOT = Path(__file__).resolve().parents[1]


def test_help_and_discovery_need_no_machine_setup():
    overview = run_command()
    assert overview.returncode == 0, overview.stderr
    assert "usage:" in overview.stdout
    # -S excludes installed scientific packages: even a fresh Python can show help.
    help_result = subprocess.run(
        [sys.executable, "-B", "-S", "-m", "regression_tests", "check", "--help"],
        cwd=ROOT.parent, capture_output=True, text=True,
    )
    assert help_result.returncode == 0, help_result.stderr
    assert "usage:" in help_result.stdout
    doctor = subprocess.run(
        [sys.executable, "-B", "-S", "-m", "regression_tests", "doctor"],
        cwd=ROOT.parent, capture_output=True, text=True,
    )
    assert doctor.returncode == 1 and "FAIL Python jsonschema" in doctor.stdout
    assert "Traceback" not in doctor.stderr
    environment = {"MHDG_REGRESSION_SETTINGS": "/nonexistent/settings.json"}
    for listing in ("cases", "suites", "layouts", "workflows"):
        result = run_command("list", listing, environment=environment)
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip()
        if listing == "cases":
            case = result.stdout.split(":", 1)[0]
    result = run_command("list", "workflows", case, environment=environment)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip()


def test_usage_and_runtime_failures_have_distinct_exit_codes(tmp_path):
    usage = run_command("check", "--unknown-option")
    assert usage.returncode == 2
    assert "unrecognized arguments" in usage.stderr
    missing = run_command("check", "--settings", str(tmp_path / "missing.json"))
    assert missing.returncode == 1
    assert missing.stderr.startswith("error: cannot read machine settings")
    assert "Traceback" not in missing.stderr


def test_debug_prepare_uses_discovered_catalogs(tmp_path):
    harness = create_harness(tmp_path)
    result = run_command(
        "prepare", "legacy_case", "cold_step_adaptive",
        "--run-id", "cli-smoke", "--settings", str(harness.settings), catalog=harness.catalog,
    )
    assert result.returncode == 0, result.stderr
    run = harness.run_directory("cold_step_adaptive", "serial_omp1", "cli-smoke")
    plan = json.loads((run / "run_plan.json").read_text())
    assert plan["stages"]
    assert "run prepared:" in result.stdout
    assert not (run / "run_metadata.json").exists()

    result = run_command(
        "prepare", "legacy_case", "cold_step_adaptive", "--layout", "mpi4_omp4",
        "--run-id", "override", "--settings", str(harness.settings), catalog=harness.catalog,
    )
    assert result.returncode == 0, result.stderr
    assert harness.run_directory("cold_step_adaptive", "mpi4_omp4", "override").is_dir()
